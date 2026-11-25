#!/usr/bin/env python3
"""Faster, auditable computation of ONE Newtonian Xu Appendix-A Fig. 8 ordinate.

USAGE: python one_point_fast.py --config configs/fig8_point.json --output-root runs

Local 01--09 modules in physics/ implement the frozen Xu Appendix-A model.
This driver changes the computation graph and integration scheduling, NOT the
physics model, angular projection, Coulomb branches, amplitude conversion, or
flux convention. It has an original adaptive method for each expensive kernel.

Model: prescribed Newtonian eccentric Kepler orbit, linear O(q) scalar-source
response, leading-order hydrogenic cloud, Coulomb/Whittaker NR radial solver,
Xu KG asymptotic momentum for the flux. Not an eccentric finite-q molecule.

This is an optimization implementation, not a convergence certificate.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
import time
import traceback
from typing import Any

import numpy as np
import scipy
import mpmath
import sympy

from batched_quadrature import fourier_at_radius, radial_overlap_batch
from primitive_store import write_physics, write_config_first, radial_grid, SCHEMA

MODULE_NAMES = (
    '01_atom', '02_kepler_orbit', '03_tidal_field', '04_fourier_forcing',
    '05_angular_coupling', '06_radial_source', '07_green_response',
    '08_flux', '09_response_channel_catalogue',
)
FORCING_MODES = ('adaptive_cached', 'gauss_batched_fallback')
RADIAL_MODES = ('reference_per_channel', 'grouped_quad_vec')


def load_modules(module_dir: Path) -> dict[str, Any]:
    module_dir = (module_dir or Path(__file__).resolve().parent/'physics').resolve()
    if not all((module_dir / (name + '.py')).is_file() for name in MODULE_NAMES):
        missing = [name + '.py' for name in MODULE_NAMES if not (module_dir / (name + '.py')).is_file()]
        raise FileNotFoundError(f'Missing modules at {module_dir}: {missing}')
    if str(module_dir) in sys.path:
        sys.path.remove(str(module_dir))
    sys.path.insert(0, str(module_dir))
    # Prevent silently reusing modules imported from another checkout during
    # programmatic runs. Tests can still explicitly monkeypatch modules.
    for name in MODULE_NAMES:
        cached = sys.modules.get(name)
        existing = getattr(cached, '__file__', None)
        if existing and Path(existing).resolve().parent != module_dir:
            sys.modules.pop(name, None)
    return {name: importlib.import_module(name) for name in MODULE_NAMES}


def read_config(path: Path) -> dict[str, Any]:
    values = json.loads(path.read_text(encoding='utf-8'))
    fields = ('M','alpha','M_b_over_M','q','N_i','ell_i','m_i','p','e',
              'ell_f_max','g_values','r_min','r_max','log_points',
              'orbit_interval_points','forcing_abs_tol','forcing_rel_tol',
              'forcing_limit','gauss_initial_order','gauss_max_order',
              'response_abs_tol','response_rel_tol','response_limit','mpmath_dps',
              'forcing_engine','radial_engine')
    missing = [key for key in fields if key not in values]
    if missing:
        raise ValueError(f'Missing config entries: {missing}')
    if 'save_source_arrays' in values:
        raise ValueError('Obsolete config key save_source_arrays: v3 always retains essential Fourier primitives')
    if values['forcing_engine'] not in FORCING_MODES:
        raise ValueError(f'forcing_engine must be one of {FORCING_MODES}')
    if values['radial_engine'] not in RADIAL_MODES:
        raise ValueError(f'radial_engine must be one of {RADIAL_MODES}')
    g = values['g_values']
    if not g or not all(isinstance(v,int) and not isinstance(v,bool) for v in g) or len(set(g)) != len(g):
        raise ValueError('g_values must be a nonempty set of distinct integer labels')
    if values['N_i'] <= 0 or values['ell_i'] < 0 or abs(values['m_i']) > values['ell_i']:
        raise ValueError('invalid hydrogenic quantum numbers')
    if values['ell_f_max'] < 0 or values['r_min'] <= 0 or values['r_max'] <= values['r_min']:
        raise ValueError('invalid truncation/radial interval')
    if values['log_points'] < 4 or values['orbit_interval_points'] < 0:
        raise ValueError('invalid radial sampling')
    if values['M'] <= 0 or values['alpha'] <= 0 or values['M_b_over_M'] <= 0 or values['q'] <= 0:
        raise ValueError('physical scales must be positive')
    if not 0 <= values['e'] < 1 or values['p'] <= 0:
        raise ValueError('must use a bound Kepler orbit (0<=e<1, p>0)')
    return values


def build_grid(c: dict[str, Any], orbit: Any) -> np.ndarray:
    """Compatibility wrapper: the config alone determines the grid."""
    return radial_grid(c)


@dataclass(frozen=True, slots=True)
class CachedForcingAdapterBase:
    """Only an explanation type; the real subclass is built after imports.

    Its sample method is the exact StageIIIForcingAdapter interface. We create
    the subclass dynamically because the authoritative class lives in the repo.
    """


def make_cached_provider(base_class: type, forcing_mod, tidal_mod):
    @dataclass(frozen=True, slots=True)
    class CachedForcingAdapter(base_class):
        g_values: tuple[int, ...] = (1,)
        engine: str = 'gauss_batched_fallback'
        gauss_initial_order: int = 16
        gauss_max_order: int = 128
        _cached: dict[tuple[int,int,int], Any] = field(default_factory=dict, init=False, repr=False, compare=False)
        _grid: np.ndarray | None = field(default=None, init=False, repr=False, compare=False)
        _stats: dict[str, int] = field(default_factory=lambda: {
            'distinct_lm': 0, 'cache_hits': 0, 'cache_misses': 0,
            'batched_radius_integrals': 0, 'gauss_node_evaluations': 0,
            'reference_fallback_radius_integrals': 0,
            'adaptive_unconverged_coefficients': 0,
        }, init=False, repr=False, compare=False)

        def sample(self, l: int, m: int, g: int, radial_points) -> Any:
            radii = np.asarray(radial_points, dtype=float)
            if g not in self.g_values:
                raise ValueError(f'g={g} not in configured g_values')
            if self.available_degrees is not None and l not in self.available_degrees:
                raise KeyError(f'forcing l={l} absent from coverage')
            if self._grid is None:
                object.__setattr__(self, '_grid', radii.copy())
            elif not np.array_equal(radii, self._grid):
                raise ValueError('A cached provider is valid only on one fixed radial grid')
            key = (l, m, g)
            if key in self._cached:
                self._stats['cache_hits'] += 1
                return self._cached[key]
            self._stats['cache_misses'] += 1
            if self.engine == 'adaptive_cached':
                result = base_class.sample(self, l, m, g, radii)
                self._cached[key] = result
                self._stats['adaptive_unconverged_coefficients'] += int(
                    np.count_nonzero(~np.asarray(result.quadrature_converged, dtype=bool))
                )
                return result
            self._calculate_lm_batch(l, m, radii)
            return self._cached[key]

        def _calculate_lm_batch(self, l: int, m: int, radii: np.ndarray) -> None:
            self._stats['distinct_lm'] += 1
            size, n_g = len(radii), len(self.g_values)
            coeff = np.zeros((size, n_g), dtype=np.complex128)
            error = np.zeros((size, n_g), dtype=float)
            evaluations = np.zeros((size, n_g), dtype=np.int64)
            success = np.ones((size, n_g), dtype=bool)
            structures = []
            # No interpolation between orbital radii: each original Fourier
            # coefficient is computed at the SAME physical radial-grid point.
            for i, radius in enumerate(radii):
                structure = tidal_mod.radial_branch_structure(self.orbit, float(radius))
                structures.append(structure)
                if self.orbit.eccentricity == 0.0:
                    u = complex(tidal_mod.tidal_multipole(
                        self.orbit.state_from_eccentric_anomaly(0.0),
                        float(radius), l, m, self.mass_ratio,
                    ))
                    coeff[i, :] = np.array([u if g == m else 0j for g in self.g_values])
                    continue
                batch = fourier_at_radius(
                    self.orbit, tidal_mod.tidal_multipole,
                    tidal_mod.radial_branch_structure, l, m, self.g_values,
                    float(radius), self.mass_ratio,
                    initial_order=self.gauss_initial_order,
                    max_order=self.gauss_max_order,
                    epsabs=self.settings.epsabs,
                    epsrel=self.settings.epsrel,
                )
                self._stats['batched_radius_integrals'] += 1
                self._stats['gauss_node_evaluations'] += batch.evaluations
                if batch.accepted:
                    coeff[i, :] = batch.values
                    error[i, :] = batch.errors
                    evaluations[i, :] = batch.evaluations
                else:
                    # Fail closed: the original branch-aware adaptive
                    # integrals replace an inconclusive finite-order Gauss rule.
                    self._stats['reference_fallback_radius_integrals'] += 1
                    for j, g in enumerate(self.g_values):
                        original = forcing_mod.forcing_coefficient(
                            l, m, g, float(radius), self.orbit,
                            self.mass_ratio, self.settings,
                        )
                        coeff[i,j] = original.value
                        error[i,j] = original.estimated_absolute_error
                        evaluations[i,j] = original.evaluations
                        success[i,j] = original.converged
                        if not original.converged:
                            self._stats['adaptive_unconverged_coefficients'] += 1
            for j, g in enumerate(self.g_values):
                phase = complex(np.exp(1j * g * self.orbit.omega_orb * self.orbit.periapsis_time))
                self._cached[(l,m,g)] = forcing_mod.ForcingSample(
                    l, m, g, radii.copy(), self.mass_ratio,
                    coeff[:, j] * phase, error[:, j], evaluations[:, j],
                    ('mean_anomaly' if self.orbit.eccentricity == 0.0 else 'eccentric_anomaly'), self.output_basis,
                    forcing_mod.FourierBasis.MEAN_ANOMALY, phase,
                    tuple(structures), success[:, j],
                )

    return CachedForcingAdapter


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def make_initial_row(meta) -> dict[str,Any]:
    return dict(
        ell_f=meta.label.ell_f, m_f=meta.label.m_f, g=meta.label.g,
        forcing_m=meta.forcing_m, delta_g=meta.delta_g,
        spectral_sector=meta.spectral_sector,
        parity_allowed=meta.equatorial_parity_allowed,
        allowed_l=meta.allowed_l_values,
        equatorial_l=meta.equatorial_l_values,
        gaunt_values=meta.gaunt_values,
        status='pending' if meta.eligible_for_flux_computation else 'excluded',
        rejection_reasons=meta.rejection_reasons,
        orbital_energy_flux=0.0,
        normalized_flux=0.0,
        abs_amplitude_sq=0.0,
        field_energy_flux=0.0,charge_flux=0.0,
        k_nr=None,k_kg=None,
        kinematic_factor_orbital=None,
        I_re=None, I_im=None, W_re=None, W_im=None,
        A_out_re=None, A_out_im=None,
        radial_error_norm=None, radial_evaluations=None,
        radial_success=None, source_max_error=None,
        source_max_abs=None,source_integral_abs_sq=None,
        forcing_quadrature_all_converged=None,
    )


def run(c: dict[str,Any], module_dir: Path | None, output_dir: Path, *, catalogue_only: bool=False) -> dict[str,Any]:
    module_dir = (module_dir or Path(__file__).resolve().parent/'physics').resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite an existing run folder: {output_dir}")
    if not catalogue_only:
        # Record effective inputs BEFORE loading code or beginning any physics.
        write_config_first(output_dir, c)
    time_start = time.perf_counter()
    marks: dict[str,float] = {}
    modules = load_modules(module_dir)
    atom, kepler, tidal, forcing, source_mod, response_mod, flux_mod, catalogue_mod = (
        modules[n] for n in ('01_atom','02_kepler_orbit','03_tidal_field','04_fourier_forcing',
                             '06_radial_source','07_green_response','08_flux',
                             '09_response_channel_catalogue')
    )
    state = atom.HydrogenicState(
        atom.AtomicParameters(c['M'],c['alpha']/c['M'],c['M']*c['M_b_over_M']),
        atom.QuantumNumbers.from_principal(c['N_i'],c['ell_i'],c['m_i']),
    )
    orbit = kepler.KeplerOrbit(c['M'],c['p'],c['e'],c.get('t_p',0.0),c.get('varpi',0.0))
    radii = radial_grid(c)
    grid_bounds=(float(radii[0]),float(radii[-1]))
    c_g=tuple(int(g) for g in c['g_values'])
    forcing_l_max=c['ell_f_max']+c['ell_i']
    catalogue = catalogue_mod.build_response_channel_catalogue(
        point_name=c.get('point_name','unnamed'), initial_state=state,
        orbit=orbit,ell_f_max=c['ell_f_max'],g_values=c_g,
        forcing_l_max=forcing_l_max,equatorial=True,
    )
    marks['setup_and_catalogue_seconds'] = time.perf_counter()-time_start
    rows_by_label = {(m.label.ell_f,m.label.m_f,m.label.g):make_initial_row(m)
                     for m in catalogue.records}
    all_rows = list(rows_by_label.values())
    base_summary = {
        'point':c.get('point_name','unnamed'),
        'alpha':c['alpha'],'e':c['e'],'p':c['p'],
        'M':c['M'],'mu':c['alpha']/c['M'],
        'M_b':c['M']*c['M_b_over_M'],'q':c['q'],
        'atom_energy':state.energy,'omega_orb':orbit.omega_orb,
        'periapsis':orbit.periapsis_radius,'apoapsis':orbit.apoapsis_radius,
        'semi_major_axis':orbit.semi_major_axis,
        'radii_count':len(radii),
        'catalogue_counts':catalogue.summary(),
        'model':'Xu small-q, Newtonian Kepler, NR Coulomb response, KG flux kinematics',
        'warning':'NOT a convergence-tested Figure 8 reproduction',
    }
    if catalogue_only:
        print('Catalogue only:',base_summary['catalogue_counts'])
        return base_summary

    CachedProvider=make_cached_provider(forcing.StageIIIForcingAdapter,forcing,tidal)
    provider=CachedProvider(
        orbit=orbit,mass_ratio=c['q'],
        output_basis=forcing.FourierBasis.COORDINATE_TIME,
        settings=forcing.QuadratureSettings(c['forcing_abs_tol'],c['forcing_rel_tol'],c['forcing_limit']),
        available_degrees=frozenset(range(forcing_l_max+1)),
        g_values=c_g,engine=c['forcing_engine'],
        gauss_initial_order=c['gauss_initial_order'],gauss_max_order=c['gauss_max_order'],
    )
    backend = response_mod.MpmathWhittakerBackend(c['mpmath_dps'])
    response_settings = response_mod.ResponseIntegrationSettings(
        *grid_bounds, c['response_abs_tol'],c['response_rel_tol'],c['response_limit'],
        breakpoints=tuple(x for x in (orbit.periapsis_radius,orbit.apoapsis_radius)
                          if grid_bounds[0] < x < grid_bounds[1]),
    )
    partial_rows=[]
    group_rows=[]
    source_cache: dict[tuple[int,int,int],np.ndarray] = {}
    amplitude_cache: dict[tuple[int,int,int],complex] = {}
    sampled_zero_keys: list[tuple[int,int,int]] = []
    groups=defaultdict(list)
    for meta in catalogue.flux_records():
        groups[(meta.label.ell_f,meta.label.g)].append(meta)
    t_compute=time.perf_counter()
    complete_groups=0
    print(f"Candidate={len(all_rows)}, physics-eligible={sum(len(v) for v in groups.values())}, "
          f"(ell_f,g) groups={len(groups)}",flush=True)
    try:
        for (ell_f,g),metas in sorted(groups.items()):
            group_t=time.perf_counter()
            active=[]
            for meta in metas:
                source=source_mod.scalar_source_channel(
                    radii,initial_state=state,forcing_provider=provider,
                    forcing_basis=forcing.FourierBasis.COORDINATE_TIME,
                    **meta.source_kwargs(),
                )
                if source.total_source is None or source.missing_l_values.size:
                    raise RuntimeError(f'Incomplete forcing source: {meta.label}')
                row=rows_by_label[(ell_f,meta.label.m_f,g)]
                row['forcing_quadrature_all_converged']=all(
                    np.all(np.asarray(f.quadrature_converged))
                    for f in source.inherited_forcing_metadata
                )
                if not row['forcing_quadrature_all_converged']:
                    raise ArithmeticError(f'Unconverged Fourier coefficient: {meta.label}')
                if not np.any(np.asarray(source.total_source) != 0j):
                    row['status']='sampled_zero'
                    sampled_zero_keys.append((ell_f,meta.label.m_f,g))
                    continue
                active.append((meta,source))
            if not active:
                group_rows.append(dict(ell_f=ell_f,g=g,n_active=0,radial_neval=0,
                                       radial_error_norm=0.0,elapsed_seconds=time.perf_counter()-group_t))
                continue

            if c['radial_engine']=='reference_per_channel':
                for meta,source in active:
                    sampled=response_mod.SampledRadialSource.from_stage_iv(source)
                    res=response_mod.solve_asymptotic_response(
                        operator=source.radial_operator,source=sampled,
                        settings=response_settings,backend=backend,
                    )
                    if res.status is not response_mod.StageVResponseStatus.COMPLETE:
                        raise RuntimeError(f'Stage-V response incomplete: {res.status} for {meta.label}')
                    fl=flux_mod.flux_channel(res,kinematics_mode=flux_mod.FluxKinematicsMode.XU_KG_KINEMATIC)
                    if fl.status is not flux_mod.FluxChannelStatus.OPEN:
                        raise RuntimeError(f'Non-open flux status: {fl.status} for {meta.label}')
                    overlaps=np.asarray(res.overlap.partial_overlaps,dtype=complex)
                    W=complex(res.wronskian)
                    I=complex(res.overlap.total_overlap)
                    A=complex(fl.amplitude)
                    error_norm=float(res.overlap.quadrature_error_bound)
                    evals=-1  # source implementation has several separate adaptive integrals
                    row=rows_by_label[(ell_f,meta.label.m_f,g)]
                    update_completed_row(row,I,W,A,fl,source,error_norm,evals,True,c)
                    key=(ell_f,meta.label.m_f,g)
                    source_cache[key]=np.asarray(source.total_source,dtype=np.complex128).copy()
                    amplitude_cache[key]=A
                    append_partial_rows(partial_rows,meta,source,overlaps)
            else:
                # One basis / Wronskian per (ell_f,g), independent of m_f.
                parameters=response_mod.WhittakerParameters.from_operator(active[0][1].radial_operator)
                basis=response_mod.CoulombHomogeneousBasis(parameters,backend)
                basis.require_invertible()
                W=complex(basis.analytic_wronskian())
                if not (np.isfinite(W.real) and np.isfinite(W.imag)) or W == 0j:
                    raise ArithmeticError(f'Invalid Wronskian at ell_f={ell_f},g={g}')
                flat_parts=[]
                spans=[]
                for meta,source in active:
                    parts=np.asarray(source.partial_sources,dtype=np.complex128)
                    start=len(flat_parts)
                    flat_parts.extend(parts)
                    spans.append((start,len(flat_parts)))
                if not flat_parts:
                    raise RuntimeError('active source unexpectedly has no partials')
                integral=radial_overlap_batch(
                    radii,np.stack(flat_parts,axis=0),basis.inner,
                    epsabs=response_settings.epsabs,epsrel=response_settings.epsrel,
                    limit=response_settings.limit,
                    breakpoints=response_settings.breakpoints,
                )
                if not integral.success:
                    raise ArithmeticError(f'Batch radial integral failed: {ell_f},{g}: {integral.message}')
                eta=flux_mod.coulomb_eta(c['M'],c['alpha']/c['M'],parameters.wave_number)
                factor=math.exp(math.pi*eta/2.0)
                for (meta,source),(a,b) in zip(active,spans):
                    overlaps=integral.partial_integrals[a:b]
                    I=complex(np.sum(overlaps))  # coherent l-sum before |.|^2
                    A=factor*I/W               # Xu outgoing-unit-Coulomb-phase convention
                    fl=flux_mod.channel_flux_from_amplitude(
                        A,scalar_mass=c['alpha']/c['M'],initial_energy=state.energy,
                        g=g,omega_orb=orbit.omega_orb,ell_f=ell_f,m_f=meta.label.m_f,
                        ell_i=c['ell_i'],m_i=c['m_i'],
                        kinematics_mode=flux_mod.FluxKinematicsMode.XU_KG_KINEMATIC,
                    )
                    if fl.status is not flux_mod.FluxChannelStatus.OPEN:
                        raise RuntimeError(f'Non-open flux status {fl.status}: {meta.label}')
                    row=rows_by_label[(ell_f,meta.label.m_f,g)]
                    update_completed_row(row,I,W,A,fl,source,integral.error_norm,integral.evaluations,True,c)
                    key=(ell_f,meta.label.m_f,g)
                    source_cache[key]=np.asarray(source.total_source,dtype=np.complex128).copy()
                    amplitude_cache[key]=A
                    append_partial_rows(partial_rows,meta,source,overlaps)
                evals=integral.evaluations
                error_norm=integral.error_norm
            complete_groups+=1
            group_rows.append(dict(
                ell_f=ell_f,g=g,n_active=len(active),radial_neval=evals,
                radial_error_norm=error_norm,
                elapsed_seconds=time.perf_counter()-group_t,
            ))
            if complete_groups % c.get('progress_every_groups',10)==0 or complete_groups==len(groups):
                print(f"completed group {complete_groups}/{len(groups)} "
                      f"(ell_f={ell_f},g={g}); Fourier cached={len(provider._cached)}",flush=True)
    except BaseException:
        # A failure is never represented by a completed HDF5 record.
        (output_dir/'FAILED.txt').write_text(traceback.format_exc(),encoding='utf-8')
        raise

    for row in all_rows:
        if row['status']=='pending':
            raise RuntimeError(f'eligible channel was not completed: {row}')
    marks['physics_seconds']=time.perf_counter()-t_compute
    total=math.fsum(row['orbital_energy_flux'] for row in all_rows)
    norm=total/(c['q']**2*c['M_b_over_M'])
    summary=dict(base_summary)
    summary.update(
        raw_total_orbital_energy_flux=total,
        normalized_fig8_flux=norm,
        number_computed_open=sum(row['status']=='open' for row in all_rows),
        number_sampled_zero=sum(row['status']=='sampled_zero' for row in all_rows),
        number_structurally_excluded=sum(row['status']=='excluded' for row in all_rows),
        forcing_cache_statistics=provider._stats,
        forcing_cache_entries=len(provider._cached),
    )
    # Store only unique open-channel numerical data. This intentionally retains
    # the complex S_c and Z_infinity=A_out requested for direct physical analysis,
    # even though parts of them are reconstructible from u_lmg and I_l.
    source_hashes={name+'.py':file_sha256(module_dir/(name+'.py')) for name in MODULE_NAMES}
    for name in ('one_point_fast.py','batched_quadrature.py','primitive_store.py',
                 'derive_quantities.py','radial_diagnostics.py','compare_runs.py'):
        path=Path(__file__).parent/name
        if path.is_file():
            source_hashes[name]=file_sha256(path)
    provenance={
        'code_sha256':source_hashes,
        'python':sys.version.split()[0], 'numpy':np.__version__,
        'scipy':scipy.__version__, 'mpmath':mpmath.__version__,
        'sympy':sympy.__version__,
    }
    # These few scalars are non-derivable numerical diagnostics, NOT error bars.
    forcing_error_max=max((float(np.max(np.asarray(s.estimated_absolute_error)))
                           for s in provider._cached.values()),default=0.0)
    quality={
        'fourier_all_accepted':True,
        'max_reported_fourier_abs_error':forcing_error_max,
        'max_reported_radial_group_error_norm':max(
            (float(row['radial_error_norm']) for row in group_rows),default=0.0),
        'fourier_fallback_radius_integrals':provider._stats.get(
            'reference_fallback_radius_integrals',0),
        'elapsed_seconds':time.perf_counter()-time_start,
    }
    write_physics(output_dir/'physics.h5',config=c,
        forcing_cache=provider._cached,source_cache=source_cache,
        amplitude_cache=amplitude_cache,sampled_zero_keys=sampled_zero_keys,
        partial_records=partial_rows,provenance=provenance,quality=quality)
    print(f'DONE: normalized Fig8 ordinate={norm:.16e} | raw={total:.16e}',flush=True)
    print(f'Completed run: config.json + physics.h5 in {output_dir}',flush=True)
    return summary


def update_completed_row(row: dict[str,Any],I:complex,W:complex,A:complex,flux,source,
                         error_norm:float,evals:int,success:bool,c:dict[str,Any]) -> None:
    raw=float(flux.orbital_energy_loss)
    if not math.isfinite(raw):
        raise ArithmeticError(f'non-finite flux for {source.ell_f,source.m_f,source.g}')
    row.update(
        status='open',I_re=I.real,I_im=I.imag,W_re=W.real,W_im=W.imag,
        A_out_re=A.real,A_out_im=A.imag,abs_amplitude_sq=abs(A)**2,
        orbital_energy_flux=raw,normalized_flux=raw/(c['q']**2*c['M_b_over_M']),
        radial_error_norm=error_norm,radial_evaluations=evals,radial_success=success,
        field_energy_flux=flux.field_energy_flux,charge_flux=flux.charge_flux,
        k_nr=flux.k_nr,k_kg=flux.k_kg,
        kinematic_factor_orbital=flux.kinematic_factor_orbital,
    )


def append_partial_rows(out: list[dict[str,Any]],meta,source,values: np.ndarray) -> None:
    for l,I_l in zip(source.l_values,values):
        out.append(dict(
            ell_f=meta.label.ell_f,m_f=meta.label.m_f,g=meta.label.g,
            forcing_l=int(l),
            I_l_re=float(I_l.real),I_l_im=float(I_l.imag),
        ))


def parse_args():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config',type=Path,default=root/'configs'/'fig8_point.json')
    p.add_argument('--module-dir',type=Path,default=root/'physics',
                   help='Optional override; defaults to bundled physics/ modules')
    p.add_argument('--output',type=Path,default=None,help='Explicit new run folder (must not exist)')
    p.add_argument('--output-root',type=Path,default=root/'runs',help='Parent folder for automatically named run directories')
    p.add_argument('--catalogue-only',action='store_true',help='Build metadata without Fourier/Green solves')
    p.add_argument('--profile',action='store_true',help='Save cProfile results to profile.pstats')
    p.add_argument('--forcing',choices=FORCING_MODES,help='Override forcing engine in JSON')
    p.add_argument('--radial',choices=RADIAL_MODES,help='Override radial engine in JSON')
    return p.parse_args()


def main():
    args=parse_args()
    cfg=read_config(args.config)
    if args.forcing:
        cfg['forcing_engine']=args.forcing
    if args.radial:
        cfg['radial_engine']=args.radial
    if args.output is not None:
        destination=args.output
    elif args.catalogue_only:
        destination=args.output_root/'catalogue_only_no_files'
    else:
        import re
        name=re.sub(r'[^A-Za-z0-9_-]+','_',str(cfg.get('point_name','fig8_point'))).strip('_')[:64] or 'fig8_point'
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
        destination=args.output_root/f'{name}_{stamp}'
    if args.profile:
        import cProfile
        prof=cProfile.Profile()
        prof.enable()
        try:
            run(cfg,args.module_dir,destination,catalogue_only=args.catalogue_only)
        finally:
            prof.disable()
            destination.mkdir(parents=True,exist_ok=True)
            prof.dump_stats(str(destination/'profile.pstats'))
    else:
        run(cfg,args.module_dir,destination,catalogue_only=args.catalogue_only)


if __name__=='__main__':
    main()
