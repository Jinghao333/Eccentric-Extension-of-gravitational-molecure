"""Spectral measurements of the exported conservative interior operator."""
import numpy as np
from scipy import linalg, sparse, special
from scipy.sparse.linalg import eigsh, LinearOperator


def load_operator(case, size):
    records = np.memmap(case / "interior.coo", dtype=[
        ("row", "<i4"), ("col", "<i4"), ("real", "<f8"), ("imag", "<f8")], mode="r")
    data = records["real"]
    if np.any(records["imag"]):
        data = data + 1j * records["imag"]
    h = sparse.csr_matrix((data, (records["row"], records["col"])), shape=(size, size))
    difference = h - h.conj().T
    error = np.max(np.abs(difference.data), initial=0.)
    scale = np.max(np.abs(h.data))
    if error > 100 * np.finfo(float).eps * scale:
        raise ValueError(f"Interior Hamiltonian fails Hermiticity check: {error:g}")
    return h, float(error)


def analytical_target(meta, order):
    """Direct COM angular integral of the translated Cartesian +m wavefunction."""
    modes = meta["modes"]
    r = meta["epsilon"] + np.arange(meta["N"]) * meta["dr"]
    a = 1 / meta["mu"]**2
    x = meta["x"]
    out = np.zeros((len(modes), len(r)), complex)
    if x == 0:
        out[modes.index([1, 1])] = r**2 * np.exp(-r/(2*a)) / (2*np.sqrt(6)*a**2.5)
        return out
    z, weights = np.polynomial.legendre.leggauss(order)
    phi = np.arange(2*order) * np.pi/order
    xy = (np.sqrt(1-z[:, None]**2) * np.exp(1j*phi)).ravel()
    weights = np.repeat(weights, len(phi)) * np.pi/order
    harmonics = []
    for l, m in modes:
        if hasattr(special, "sph_harm_y"):
            y = special.sph_harm_y(l, m, np.arccos(z)[:, None], phi)
        else:
            y = special.sph_harm(m, l, phi, np.arccos(z)[:, None])
        harmonics.append(y.ravel().conj() * weights)
    projection = np.array(harmonics)
    for start in range(0, len(r), 64):
        radius = r[None, start:start+64]
        distance = np.sqrt(np.maximum(radius**2 + x*x - 2*x*xy.real[:, None]*radius, 0))
        field = -(xy[:, None]*radius-x)*np.exp(-distance/(2*a))/(8*np.sqrt(np.pi)*a**2.5)
        out[:, start:start+64] = (projection @ field) * radius
    return out


def lanczos_weights(h, u, iterations, history=None):
    """Gauss quadrature for the spectral measure of u, with threshold brackets.

    At either Gaussian node adjacent to zero the Chebyshev-Markov inequalities
    put the CDF between its exclusive and inclusive quadrature sums. Brackets
    are numerical estimates in finite precision; check at multiple orders.
    """
    q = u / np.linalg.norm(u)
    previous = np.zeros_like(q)
    diagonal, off = [], []
    for i in range(min(iterations, len(u))):
        residual = h @ q
        if i:
            residual -= off[-1]*previous
        alpha = np.vdot(q, residual).real
        residual -= alpha*q
        beta = np.linalg.norm(residual)
        diagonal.append(alpha)
        if history is not None and len(diagonal) in (32,64):
            history.append(threshold_bracket(diagonal,off))
        if beta < np.finfo(float).eps * max(1., abs(alpha)) or i+1 == min(iterations, len(u)):
            break
        previous, q = q, residual/beta
        off.append(beta)
    return threshold_bracket(diagonal,off)


def threshold_bracket(diagonal, off):
    energies, vectors = linalg.eigh_tridiagonal(diagonal, off)
    weights = np.abs(vectors[0])**2
    k = np.searchsorted(energies, 0., side="right")
    return dict(order=len(diagonal), estimate=float(weights[k:].sum()),
                lower=float(weights[k+1:].sum()), upper=float(weights[max(0,k-1):].sum()))


def sign_transform(h, nr, e211):
    """Reorder the *exported* operator by radius and factor its narrow band.

    Once H+sI is positive definite, E/(E+s) has exactly the sign of E.
    This resolves zero without spending most of the Lanczos order on the
    high positive kinetic/centrifugal spectrum. No time propagation uses it.
    """
    nmode=h.shape[0]//nr
    band=np.zeros((nmode+1,h.shape[0]),dtype=h.dtype,order="F")
    for offset in range(nmode):
        block=h.diagonal(-offset*nr).reshape(nmode-offset,nr)
        band[offset].reshape(nr,nmode)[:,:nmode-offset]=block.T
    links=np.append(h.diagonal(-1),0).reshape(nmode,nr).T[:-1]
    band[nmode,:-nmode]=links.ravel()
    shift=8*abs(e211)
    for _ in range(12):
        band[0]+=shift
        try:
            factor=linalg.cholesky_banded(band,lower=True,check_finite=False)
            break
        except linalg.LinAlgError:
            band[0]-=shift
            shift*=2
    else:
        raise ValueError("Could not certify a positive shift of the interior operator")

    def transform(v):
        if np.isrealobj(factor) and np.iscomplexobj(v):
            rhs=np.asfortranarray(np.column_stack((v.real,v.imag)))
            sol=linalg.cho_solve_banded((factor,True),rhs,overwrite_b=True,check_finite=False)
            solved=sol[:,0]+1j*sol[:,1]
        else:
            solved=linalg.cho_solve_banded((factor,True),v,check_finite=False)
        return v-shift*solved
    # Verify the band layout by applying the original sparse H to a solved probe.
    rng=np.random.default_rng(211)
    probe=rng.normal(size=h.shape[0])
    solved=linalg.cho_solve_banded((factor,True),probe,check_finite=False)
    unstack=lambda v:v.reshape(nr,nmode).T.flatten()
    residual=h@unstack(solved)+shift*unstack(solved)-unstack(probe)
    error=float(np.linalg.norm(residual)/np.linalg.norm(probe))
    if error>1e-9:
        raise ValueError(f"Shifted solve disagrees with exported H: {error:g}")
    return LinearOperator(h.shape,matvec=transform,dtype=h.dtype),shift,error


def bound_cluster(h, e211, dense=None):
    # The smaller adjacent-level gap is E3 - E2 = (5/9)|E2|.
    half_width = abs(e211)*5/18
    if dense is not None:
        values, basis = dense
    else:
        k = 8
        while True:
            values, basis = eigsh(h, k=k, sigma=e211*(1+1e-5), which="LM",
                                  tol=1e-11, v0=np.random.default_rng(211).normal(size=h.shape[0]))
            if np.max(abs(values-e211)) > half_width:
                break
            k *= 2
            if k >= h.shape[0]-1:
                raise ValueError("Cannot isolate the 211 cluster")
    keep = abs(values-e211) < half_width
    return values[keep], basis[:, keep]


def centered_spectrum(h, u, meta, e211):
    """C0 has independent tridiagonal channels. Resolve *all* negative states
    in the populated channel and use the norm of the orthogonal remainder.
    """
    nr = meta["N"]-2
    j = meta["modes"].index([1, 1])
    section = slice(j*nr, (j+1)*nr)
    radial = h[section, section]
    v = u[section]
    diagonal, side = radial.diagonal().real, radial.diagonal(1).real
    lower = diagonal.min()-2*np.max(abs(side))-1
    levels, bound = linalg.eigh_tridiagonal(diagonal, side, select="v",
                                            select_range=(lower, 0.), tol=0.)
    amplitudes = bound.T @ v
    positive = float(np.linalg.norm(v-bound @ amplitudes)**2)
    window = abs(levels-e211) < abs(e211)*5/18
    cluster = float(np.linalg.norm(amplitudes[window])**2)
    # Record multiplicities in the retained basis, including unpopulated partners.
    cluster_levels = []
    for ell in sorted({l for l,m in meta["modes"]}):
        indices = [i for i,(l,m) in enumerate(meta["modes"]) if l==ell]
        start = indices[0]*nr
        block = h[start:start+nr,start:start+nr]
        nearby = linalg.eigh_tridiagonal(block.diagonal().real, block.diagonal(1).real,
            eigvals_only=True, select="v", select_range=(e211-abs(e211)*5/18,e211+abs(e211)*5/18))
        cluster_levels.extend(nearby.tolist()*len(indices))
    return dict(spectral_status="all negative eigenstates of populated centered radial channel",
                cluster_energies=sorted(cluster_levels), cluster_weight=cluster,
                positive_estimate=positive, positive_interval=[positive,positive],
                other_bound_interval=[float(np.linalg.norm(amplitudes[~window])**2)]*2,
                quadrature_checks=[], unresolved_threshold_weight=0.,
                target_polarization_weight=cluster, other_polarization_weight_in_cluster=0.)


def characterize(case, meta, state, e211):
    full = state.ravel()
    norm = np.linalg.norm(full)
    endpoint_fraction = np.linalg.norm(state[:, [0, -1]])**2/norm**2
    target = analytical_target(meta, max(16, 2*max(l for l, m in meta["modes"])+4))
    refined = analytical_target(meta, max(32, 4*max(l for l, m in meta["modes"])+8))
    reference = refined / np.linalg.norm(refined)
    unit = state / norm
    infidelity = float(np.linalg.norm(unit - reference*np.vdot(reference,unit))**2)
    result = dict(endpoint_fraction=float(endpoint_fraction),
                  reference_quadrature_change=float(np.linalg.norm(refined-target)/np.linalg.norm(refined)),
                  fidelity=max(0., 1.-infidelity), infidelity=infidelity,
                  mode_populations=(meta["dr"]*np.sum(abs(state)**2, axis=1)).tolist())
    if endpoint_fraction > 1e-14:
        result["spectral_status"] = "material frozen endpoints; interior fractions withheld"
        return result, None
    u = state[:, 1:-1].flatten()
    u /= np.linalg.norm(u)
    h, result["hermiticity_max_abs"] = load_operator(case, len(u))
    if meta["case"]=="C0" and len(u)>700:
        result.update(centered_spectrum(h,u,meta,e211))
        return result, None
    dense = linalg.eigh(h.toarray()) if len(u) <= 700 else None
    energies, basis = bound_cluster(h, e211, dense)
    coeff = basis.conj().T @ u
    cluster = float(np.linalg.norm(coeff)**2)
    if dense is None:
        # Deflate the explicitly measured negative-energy cluster first. This
        # keeps its large spectral mass out of the near-zero quadrature error.
        residual = u-basis @ coeff
        weight = float(np.linalg.norm(residual)**2)
        transformed,shift,solve_error=sign_transform(h,meta["N"]-2,e211)
        result.update(spectral_transform_shift=shift,shifted_solve_relative_residual=solve_error)
        reordered=residual.reshape(state.shape[0],-1).T.flatten()
        history=[]
        item=lanczos_weights(transformed,reordered,128,history) if weight else dict(order=0,estimate=0.,lower=0.,upper=0.)
        history.append(item)
        checks=[dict(order=item["order"], **{k:weight*item[k] for k in ("estimate","lower","upper")})
                for item in history]
        positive = checks[-1]["estimate"]
        low = min(v["lower"] for v in checks[-2:])
        high = max(v["upper"] for v in checks[-2:])
        if len(energies):
            # Bound the positive-energy amplitude in approximate negative Ritz vectors.
            leakage = np.sum(abs(coeff)*np.linalg.norm(h@basis-basis*energies,axis=0)/abs(energies))
            low = max(0.,np.sqrt(low)-leakage)**2
            high = min(1.,(np.sqrt(high)+leakage)**2)
    else:
        values, vectors = dense
        positive = float(np.linalg.norm(vectors[:, values>0].conj().T @ u)**2)
        low = high = positive
        checks = []
    result.update(spectral_status="finite-box interior spectrum", cluster_energies=energies.tolist(),
                  cluster_weight=cluster, positive_estimate=positive, positive_interval=[low, high],
                  quadrature_checks=checks, unresolved_threshold_weight=high-low,
                  other_bound_interval=[max(0., 1-cluster-high), max(0., 1-cluster-low)])
    target_interior = refined[:, 1:-1].flatten()
    projected = basis @ (basis.conj().T @ target_interior)
    if not len(energies) or np.linalg.norm(projected) < 0.5*np.linalg.norm(target_interior):
        result["candidate_status"] = "target cluster unresolved"
        return result, None
    projected /= np.linalg.norm(projected)
    result["target_polarization_weight"] = float(abs(np.vdot(projected, u))**2)
    result["other_polarization_weight_in_cluster"] = max(0., cluster-result["target_polarization_weight"])
    corrected = np.zeros_like(state)
    corrected[:, 1:-1] = projected.reshape(state.shape[0], -1)
    corrected *= np.sqrt(meta["input_norm"]/meta["dr"])
    return result, corrected
