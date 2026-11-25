"""Scientific analysis API: no production solver calls; no writes to run folders.

The v4 primitives are the authority.  Reconstruct cheap analytic quantities
on demand, and calculate Whittaker M only for a user-selected open channel.

Important: the *production* hydrogenic field R(r) has
    integral r**2 |R|**2 dr = M_b / (2*mu**2),
so the unperturbed mass density is 2*mu**2*|R(r)Y_lm|**2.
The normalized particle/charge convention is independently labeled below.
"""
from __future__ import annotations
from functools import lru_cache
from collections import defaultdict
import math
import os
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.special import eval_genlaguerre, gammaln, sph_harm_y, gammaincinv
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
from sympy.physics.wigner import wigner_3j

from derive_quantities import derive, _eval_wronskian


def _js_float(x):
    z=float(x)
    return z if math.isfinite(z) else None


def _complex(z):
    z=complex(z)
    return {'re':float(z.real),'im':float(z.imag),'abs':float(abs(z)),
            'phase':float(math.atan2(z.imag,z.real)) if z else None}


def radial_field(r, cfg):
    """Exact functional formula of production 01_atom.unit_radial, including 1/mu.

    This is a field profile, NOT a mass-normalized Schr wavefunction.
    """
    r=np.asarray(r,dtype=float)
    M=float(cfg['M']); mu=float(cfg['alpha'])/M
    N=int(cfg['N_i']); ell=int(cfg['ell_i']); nr=N-ell-1
    if nr<0:raise ValueError('invalid hydrogenic quantum numbers')
    beta=2*mu*float(cfg['alpha'])/N
    rho=beta*r
    log_pref=(-.5*math.log(2)-math.log(mu) + .5*(3*math.log(beta)
            +gammaln(nr+1)-math.log(2*N)-gammaln(N+ell+1)))
    return math.sqrt(float(cfg['M_b_over_M'])*M)*np.exp(log_pref)*rho**ell*eval_genlaguerre(nr,2*ell+1,rho)*np.exp(-rho/2)


def orbital_parameters(cfg):
    M=float(cfg['M']); a0=float(cfg['alpha']);p=float(cfg['p']);e=float(cfg['e'])
    mu=a0/M; N=int(cfg['N_i']); q=float(cfg['q'])
    omega=(1-e**2)**1.5/(M*p**1.5)
    return dict(M=M,alpha=a0,mu=mu,q=q,m_p=q*M,M_b=float(cfg['M_b_over_M'])*M,
                N_i=N,ell_i=int(cfg['ell_i']),m_i=int(cfg['m_i']),
                E_i=-mu*a0**2/(2*N**2),a=p*M/(1-e**2),p=p,e=e,
                periapsis=p*M/(1+e),apoapsis=p*M/(1-e),
                omega=omega,period=2*math.pi/omega,
                t_p=float(cfg.get('t_p',0)),varpi=float(cfg.get('varpi',0)),
                a_B=M/a0**2,
                model='Prescribed Newtonian Kepler orbit; small-q linear scalar response; nonrelativistic Coulomb radial operator with KG momentum in outgoing flux')


def kepler_curve(cfg,n=460):
    x=orbital_parameters(cfg);e=x['e']; E=np.linspace(0,2*math.pi,int(n))
    mean=E-e*np.sin(E)
    nu=np.arctan2(np.sqrt(1-e**2)*np.sin(E),np.cos(E)-e)
    phi=x['varpi']+nu
    r=x['a']*(1-e*np.cos(E))
    return dict(x=(r*np.cos(phi)).tolist(),y=(r*np.sin(phi)).tolist(),
                r=r.tolist(),M_a=mean.tolist(),E=E.tolist(),nu=nu.tolist(),
                phi=phi.tolist(),period=x['period'])


def shell_and_scales(cfg,n=2800):
    """Mass profile and exact/derived characteristic radii.

    For n=2,l=1 the Gamma shape-5 law gives exact quantiles and peaks;
    other states use a finely resolved analytic background quadrature.
    This is *not* a response-convergence sweep.
    """
    x=orbital_parameters(cfg); aB=x['a_B']; N=x['N_i'];ell=x['ell_i'];mu=x['mu'];Mb=x['M_b']
    radius_end=max(x['apoapsis']*1.45, aB*12*N**2)
    # Hybrid near-origin/linear radial plot; sufficient for analytic cloud measures.
    r=np.linspace(0,radius_end,int(n))
    radial=radial_field(r,cfg)
    shell=2*mu**2*r**2*radial**2
    integral=cumulative_trapezoid(shell,r,initial=0)
    norm=float(integral[-1])
    if not norm>0:raise ArithmeticError('nonpositive analytical cloud norm')
    cdf=np.clip(integral/norm,0,1)
    r50=float(np.interp(.5,cdf,r));r90=float(np.interp(.9,cdf,r))
    shell_peak=float(r[int(np.argmax(shell))]); mean=float(np.trapezoid(r*shell,r)/norm)
    # equatorial angular factor constant in phi, so the radial volume-density peak
    # is the maximum of R(r)**2 for a |lm> cloud.
    dens_peak=float(r[int(np.argmax(radial**2))])
    if N==2 and ell==1:
        dens_peak=2*aB; shell_peak=4*aB;mean=5*aB
        r50=float(gammaincinv(5,.5)*aB);r90=float(gammaincinv(5,.9)*aB)
    equatorial_factor=abs(complex(sph_harm_y(ell,x['m_i'],math.pi/2,0)))**2
    rho_eq=2*mu**2*radial**2*equatorial_factor
    return dict(radius=r.tolist(),rho_equatorial=rho_eq.tolist(),
                dM_dr=shell.tolist(),enclosed_fraction=cdf.tolist(),
                mass_integral_from_analytic_profile=norm,
                M_b=Mb,bohr_radius=aB,density_peak=dens_peak,
                shell_peak=shell_peak,mean_radius=mean,r50=r50,r90=r90,
                density_definition='rho_M=2*mu^2*|R_i Y_lm|^2; unperturbed field')


def density_slice(cfg,n=170,extent=None):
    """Equatorial unperturbed physical mass density, not a time snapshot."""
    scales=shell_and_scales(cfg,n=1200)
    params=orbital_parameters(cfg)
    default_extent=1.15*max(scales['r90'],params['apoapsis'],scales['shell_peak'])
    extent=float(default_extent if extent is None else extent)
    if extent<=0:raise ValueError('extent must be positive')
    coords=np.linspace(-extent,extent,int(n))
    xx,yy=np.meshgrid(coords,coords)
    r=np.hypot(xx,yy)
    phi=np.mod(np.arctan2(yy,xx),2*math.pi)
    Y=sph_harm_y(params['ell_i'],params['m_i'],math.pi/2,phi)
    rho=2*params['mu']**2*np.abs(radial_field(r,cfg)*Y)**2
    return dict(coords=coords.tolist(),rho=rho.tolist(),extent=extent,
                max_density=float(rho.max()))


@lru_cache(maxsize=8)
def _cached_derive(path_str,config_mtime,h5_mtime):
    return derive(Path(path_str))


def load_run(path):
    path=Path(path).resolve()
    return _cached_derive(str(path),(path/'config.json').stat().st_mtime_ns,
                          (path/'physics.h5').stat().st_mtime_ns)


def channel_atlas(result):
    rows=result['channels'];cfg=result['run']['config']
    open_rows=[r for r in rows if r['status']=='open']
    total=float(result['summary']['normalized_flux'])
    ranked=sorted(open_rows,key=lambda r:r['normalized_flux'],reverse=True)
    fig11=[]
    group=defaultdict(float)
    group_raw=defaultdict(float)
    group_harmonics=defaultdict(int)
    for row in rows:
        key=(int(row['ell_f']),int(row['m_f']))
        group[key]+=float(row['normalized_flux'])
        group_raw[key]+=float(row['orbital_flux'])
        if row['status']=='open':
            group_harmonics[key]+=1
    angular_modes=[]
    for (ell,m),power in group.items():
        if power<=0:
            continue
        angular_modes.append(dict(
            ell=ell,m=m,power=float(power),raw_power=float(group_raw[(ell,m)]),
            share=float(power/total) if total else None,
            harmonic_count=int(group_harmonics[(ell,m)])))
    angular_modes.sort(key=lambda row:row['power'],reverse=True)
    L=int(cfg['ell_f_max'])
    for k in range(4):
        pnts=[]
        for ell in range(L+1):
            m=ell-2*k
            if abs(m)>ell:continue
            pnts.append(dict(ell=ell,m=m,power=group[(ell,m)]))
        fig11.append(dict(k=k,name=f'm_f = ℓ_f − {2*k}' if k else 'm_f = ℓ_f',points=pnts))
    top=[];cumulative=0
    for idx,row in enumerate(ranked):
        cumulative+=row['normalized_flux']
        if idx<30:
            top.append(dict(ell=int(row['ell_f']),m=int(row['m_f']),g=int(row['g']),
                            power=float(row['normalized_flux']),raw_power=float(row['orbital_flux']),
                            fraction=float(row['normalized_flux']/total) if total else None,
                            cumulative_fraction=float(cumulative/total) if total else None,
                            delta_g=float(row['delta_g']),k_kg=row['k_kg']))
    mvals=list(range(-L,L+1));ells=list(range(L+1))
    Z=[];STATUS=[]
    for m in mvals:
        zrow=[];statusrow=[]
        for ell in ells:
            if abs(m)>ell or (ell+m)%2 != (int(cfg['ell_i'])+int(cfg['m_i']))%2:
                zrow.append(None);statusrow.append('forbidden')
            else:
                flux=group[(ell,m)]
                zrow.append(float(flux) if flux>0 else None)
                statusrow.append('positive' if flux>0 else 'zero_or_closed')
        Z.append(zrow);STATUS.append(statusrow)
    harmonic=[]
    for row in result['harmonics']:
        g=int(row['g']);harmonic.append(dict(g=g,power=float(row['normalized_flux']),
                                        cumulative_fraction=row.get('cumulative_fraction'),
                                        delta_g=float(result['summary']['E_i']+g*result['summary']['omega_orb'])))
    return dict(fig11=fig11,top=top,harmonics=harmonic,angular_modes=angular_modes,
                angular_heatmap=dict(ell=ells,m=mvals,power=Z,status=STATUS),
                channels=[dict(ell=int(r['ell_f']),m=int(r['m_f']),g=int(r['g']),
                               status=r['status'],power=float(r['normalized_flux']),
                               fraction=r['relative_flux_fraction'],delta_g=float(r['delta_g']))
                          for r in rows],
                total_normalized_flux=total,
                top_10_fraction=sum(r['normalized_flux'] for r in ranked[:10])/total if total else None)


def run_overview(result):
    cfg=result['run']['config']; pars=orbital_parameters(cfg)
    scales=shell_and_scales(cfg)
    s=result['summary']; omega=pars['omega']; mu=pars['mu']; Mb=pars['M_b']
    radiating=[r for r in result['channels'] if r['status']=='open']
    # In the existing Stage-VI convention J=2*k_KG*|A|**2.
    J=math.fsum(2*float(r['k_kg'])*abs(complex(r['A_out_re'],r['A_out_im']))**2
                for r in radiating)
    gamma=(mu/Mb)*J if Mb>0 else None
    # Verifying expected scalar-field mass norm is not a certified convergence test.
    return dict(parameters=pars,scales={k:v for k,v in scales.items() if k not in ('radius','rho_equatorial','dM_dr','enclosed_fraction')},
                radial={k:scales[k] for k in ('radius','rho_equatorial','dM_dr','enclosed_fraction')},
                orbit=kepler_curve(cfg),
                totals=dict(orbital_power=float(s['raw_orbital_flux']),
                            normalized_power=float(s['normalized_flux']),
                            outgoing_number_current=J,
                            fractional_cloud_depletion_rate=gamma,
                            approximate_depletion_timescale=(1/gamma if gamma and gamma>0 else None),
                            period=pars['period'],
                            power_per_orbit=float(s['raw_orbital_flux'])*pars['period'],
                            n_open=int(s['n_open']),n_sampled_zero=int(s['n_sampled_zero']),
                            n_candidates=int(s['n_candidates'])),
                warnings=['Single fixed orbit: no secular orbital evolution computed.',
                          'Cloud density is the unperturbed mass density (not the driven cloud).',
                          'Gamma = mu*sum J/M_b assumes the stated leading-order KG charge/cloud-mass convention.',
                          'Finite numerical result; convergence and agreement with paper not yet certified.'],
                quality=result['run'].get('quality',{}))


@lru_cache(maxsize=10000)
def gaunt(l,m,li,mi,lf,mf):
    from math import sqrt,pi
    return ((-1)**mf*sqrt((2*l+1)*(2*li+1)*(2*lf+1)/(4*pi))
            *float(wigner_3j(l,li,lf,0,0,0))*float(wigner_3j(l,li,lf,m,mi,-mf)))


@lru_cache(maxsize=16)
def _whittaker_samples(ell,g,M,mu,Ei,Omega,precision,radius_tuple):
    r=np.asarray(radius_tuple,dtype=float)
    delta=Ei+g*Omega
    if not delta>0:raise ValueError('Whittaker radial profile requires open channel')
    k=math.sqrt(2*mu*delta)
    lam=-2j*k
    kappa=2*mu**2*M/lam
    result=np.empty(len(r),dtype=np.complex128)
    with mp.workdps(int(precision)):
        aa=ell+1-mp.mpc(kappa); bb=2*ell+2
        for i,rad in enumerate(r):
            z=mp.mpc(lam*rad)
            result[i]=complex(mp.exp(-z/2)*mp.power(z,ell+1)*mp.hyp1f1(aa,bb,z))
    return result


def channel_detail(result,ell,mf,g,include_whittaker=True):
    key=(int(ell),int(mf),int(g))
    raw=result['raw'];cfg=result['run']['config'];
    rows={(int(r['ell_f']),int(r['m_f']),int(r['g'])):r for r in result['channels']}
    if key not in rows:raise KeyError('unknown channel')
    row=rows[key]
    if row['status']!='open':raise ValueError(f'channel {key} has status {row["status"]}; no outgoing source was saved')
    r=np.asarray(raw['radius']);mu=float(cfg['alpha'])/float(cfg['M']); M=float(cfg['M']);Ei=float(result['summary']['E_i']);Omega=float(result['summary']['omega_orb']);
    R=radial_field(r,cfg)
    fdata=result['forcing']; fi=(int(cfg['ell_i']),int(cfg['m_i']));forcing_m=mf-fi[1]
    partials=sorted([p for p in result['partials'] if (int(p['ell_f']),int(p['m_f']),int(p['g']))==key],key=lambda p:p['forcing_l'])
    ls=[int(p['forcing_l']) for p in partials]
    u=[];source_part=[];overlap=[]
    for l in ls:
        u_l=np.asarray(fdata[(l,forcing_m,g)],dtype=np.complex128)
        source=2*mu**2*r*R*gaunt(l,forcing_m,fi[0],fi[1],ell,mf)*u_l
        u.append(u_l);source_part.append(source)
    partial_array=np.asarray(source_part,dtype=np.complex128)
    total=np.sum(partial_array,axis=0)
    stored=np.asarray(result['source'][key])
    source_mismatch=float(np.max(np.abs(stored-total)))
    source_scale=max(float(np.max(np.abs(stored))),float(np.max(np.abs(total))),1e-300)
    if source_mismatch/source_scale>1e-10:
        raise ValueError(f'Source formula does not match stored source: {key} relative error={source_mismatch/source_scale:.3e}')
    Il={int(p['forcing_l']):complex(float(p['I_l_re']),float(p['I_l_im'])) for p in partials}
    I=sum(Il.values(),0j)
    coherence=abs(I)/sum(map(abs,Il.values())) if any(Il.values()) else None
    incoh=sum(abs(v)**2 for v in Il.values()); enhancement=abs(I)**2/incoh if incoh else None
    interference=[{'l1':l,'l2':m,'value':2*(Il[l]*Il[m].conjugate()).real}
                  for i,l in enumerate(ls) for m in ls[i+1:]]
    W=_eval_wronskian(mu,M,ell,float(row['delta_g']),int(cfg['mpmath_dps']))
    zminus=None;cumulative=None;integrands=None;cumulative_partial=None
    if include_whittaker:
        zminus=_whittaker_samples(ell,g,M,mu,Ei,Omega,int(cfg['mpmath_dps']),tuple(map(float,r)))
        integrands=zminus[None,:]*partial_array
        cumulative_partial=cumulative_trapezoid(integrands,r,axis=1,initial=0)
        cumulative=np.sum(cumulative_partial,axis=0)
    A=complex(float(row['A_out_re']),float(row['A_out_im']))
    return dict(label=dict(ell=ell,m=mf,g=g,forcing_m=forcing_m),
                spectral=dict(Delta_g=float(row['delta_g']),omega_g=float(row['omega_g']),
                              k_nr=float(row['k_nr']),k_kg=float(row['k_kg']),
                              eta=float(row['eta']),g_omega=g*Omega),
                response=dict(I=_complex(I),W=_complex(W),A_out=_complex(A),
                              normalized_power=float(row['normalized_flux']),
                              orbital_power=float(row['orbital_flux']),
                              fraction=row['relative_flux_fraction'],
                              multipole_coherence=coherence,multipole_enhancement=enhancement,
                              source_relative_difference=source_mismatch/source_scale,
                              grid_overlap=_complex(cumulative[-1]) if cumulative is not None else None,
                              grid_minus_production_I=_complex(cumulative[-1]-I) if cumulative is not None else None),
                radial=dict(r=r.tolist(),
                            forcing=[{'l':l,'abs':np.abs(t).tolist(), 'phase':np.angle(t).tolist()} for l,t in zip(ls,u)],
                            source_partials=[{'l':l,'abs':np.abs(s).tolist(),'phase':np.angle(s).tolist()} for l,s in zip(ls,partial_array)],
                            source_total_abs=np.abs(stored).tolist(),source_total_phase=np.angle(stored).tolist(),
                            overlap_integrands=[{'l':l,'abs':np.abs(v).tolist()} for l,v in zip(ls,integrands)] if integrands is not None else [],
                            overlap_total_abs=np.abs(integrands.sum(axis=0)).tolist() if integrands is not None else None,
                            cumulative_real=cumulative.real.tolist() if cumulative is not None else None,
                            cumulative_imag=cumulative.imag.tolist() if cumulative is not None else None,
                            Z_minus_abs=np.abs(zminus).tolist() if zminus is not None else None),
                partial_overlaps=[{'l':l,**_complex(Il[l])} for l in ls],
                interference=interference,
                caveats=['Cumulative overlap uses the saved sampling grid and trapezoidal integration; compare only diagnostically with saved adaptive I.',
                         'A_out is the physical asymptotic coefficient; Z_minus is the regular Coulomb solution.',
                         'A large amplitude alone is not proof of a resonance.'])
