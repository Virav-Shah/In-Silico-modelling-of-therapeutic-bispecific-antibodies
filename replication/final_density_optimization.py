"""Our reference-informed calibrated model; not recovered author code."""
import argparse,csv,json,os,random,itertools
from pathlib import Path
from dataclasses import replace,asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import least_squares,minimize
from platypus import Problem,Real,Solution,InjectedPopulation,NSGAII
from model import Parameters,SCENARIOS,simulate
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'output/final_density_optimization'
RATES=('kon1','koff1','kon2','koff2')
BOUNDS=[(2.,7.),(-5.,2.),(2.,7.),(-5.,2.)]
LIMIT=.0095
AREA=4*np.pi*5**2*1e-10
VOLUME=4*np.pi*5**3/3*1e-15
CASES=[('very_low',1e-13,1e-13,4e-6),('low',1e-12,1e-12,4e-6),
 ('medium',1e-11,1e-11,4e-6),('equal_reference',1e-10,1e-10,4e-6),
 ('high',1e-9,1e-9,4e-6),('unequal_10',1e-9,1e-10,4e-6),
 ('unequal_20',1e-9,5e-11,4e-6),('unequal_100',1e-9,1e-11,4e-6),
 ('very_low_potent_payload',1e-13,1e-13,1e-8)]
def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False))
def csvwrite(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def score(p,tight=False,method='BDF'):
    r=simulate(p,points=3,rtol=1e-9 if tight else 1e-7,atol=1e-12 if tight else 1e-10,method=method)
    return float(r['tumor'][-1]),float(r['healthy'][-1])
def calibrate():
    base=replace(Parameters(),area_dm2=AREA,cell_volume_l=VOLUME,
                 cl_l_h=Parameters().vc*np.log(2)/(21*24),surface_off_rule='symmetric')
    anchors={'figure4_bispecific':56.,'figure6':280.}
    def residual(z):
        n,f=10**np.asarray(z)
        return [np.log(simulate(replace(base,**{k:getattr(SCENARIOS[name],k) for k in RATES},total_cells=n,cell_volume_l=VOLUME*f),points=3)['intracellular'][-1,0]*1e6/target) for name,target in anchors.items()]
    fit=least_squares(residual,[12.9,np.log10(.5)],bounds=([10,-2],[14,0]),xtol=1e-9,ftol=1e-9,gtol=1e-9)
    if not fit.success:raise RuntimeError(fit.message)
    n,f=10**fit.x;base=replace(base,total_cells=float(n),cell_volume_l=float(VOLUME*f))
    config=dict(status='OUR_REFERENCE_INFORMED_CALIBRATED_MODEL',parameters=asdict(base),geometry_diameter_um=10.,whole_cell_volume_l=VOLUME,effective_intracellular_fraction=float(f),calibration_anchors_uM=anchors,calibration_residual_log=fit.fun.tolist(),calibration_source='Approximate visual 24h uptake readings, not raw data. Fitting is not validation.',design_log10_bounds=BOUNDS,healthy_design_fraction=LIMIT,paper_healthy_limit=.01,
      decisions={'capacity':'N_j = fitted accessible cell count times literal paper population frequency; total surface = N_j * geometric area',
      'clearance':'Explicit binding/internalisation handles target-mediated removal. Background CL=Vc*ln(2)/(21*24) L/h is our assumption, not a Mould fitted parameter or terminal half-life of the two-compartment system.',
      'no_double_counting':'No empirical Mould target-related Michaelis-Menten elimination on top of explicit target removal.',
      'volume':'Effective intracellular fraction of a 10um sphere, fitted once and frozen across densities.',
      'dissociation':'q1=koff1*sqrt(r), q2=koff2/sqrt(r), r=kon2*k2d1/(kon1*k2d2). Symmetric detailed balance is our extension.',
      'search':'Two seeded NSGAII runs then constrained SLSQP. Best verified feasible candidate, not certified global optimum.',
      'domain':'Broad engineering domain, not author-reported bounds or demonstrated engineering reachability.',
      'toxicity':'Weighted Emax at 24h. Nominal limit 0.95% gives numerical headroom, not robustness assurance.'})
    dump(OUT/'configuration.json',config);print('CALIBRATED',n,f,flush=True);return base

def run_case(args):
    name,rho1,rho2,ec50,base_dict,nfe,seeds=args
    base=replace(Parameters(**base_dict),rho1=rho1,rho2=rho2,ec50=ec50)
    out=OUT/name;out.mkdir(exist_ok=True);runs=[]
    for seed in seeds:
        random.seed(seed);np.random.seed(seed);cache={};rows=[];stage='global'
        def ev(z):
            key=tuple(float(x) for x in z)
            if key not in cache:
                rates=10**np.asarray(z);e,h=score(replace(base,**dict(zip(RATES,rates))))
                rows.append(dict(seed=seed,stage=stage,**dict(zip(RATES,map(float,rates))),tumor_kill_pct=100*e,healthy_kill_pct=100*h))
                cache[key]=(100*e,100*h)
            return cache[key]
        problem=Problem(4,2,1);problem.types[:]=[Real(*b) for b in BOUNDS];problem.constraints[:]='<=0'
        def fn(z):
            e,h=ev(z);return [-e,h],[h-100*LIMIT]
        problem.function=fn;initial=[]
        for on,off in itertools.product([1e2,1e3,1e4],[.001,.01,.1,1.,10.,100.]):
            s=Solution(problem);s.variables[:]=np.log10([on,off,on,off]);initial.append(s)
        for vals in ([1e4,1,5e4,1],[1e2,1e2,1e7,1e-5],[1e7,1e-5,1e2,1e2]):
            s=Solution(problem);s.variables[:]=np.log10(vals);initial.append(s)
        alg=NSGAII(problem,population_size=32,generator=InjectedPopulation(initial));alg.run(nfe)
        feasible=sorted([x for x in rows if x['healthy_kill_pct']<=100*LIMIT],key=lambda x:x['tumor_kill_pct'],reverse=True)
        local=[];stage='local'
        if feasible:
            z0=np.log10([feasible[0][k] for k in RATES])
            fit=minimize(lambda z:-ev(z)[0],z0,method='SLSQP',bounds=BOUNDS,
                         constraints=[{'type':'ineq','fun':lambda z:100*LIMIT-ev(z)[1]}],
                         options={'maxiter':60,'ftol':1e-6,'eps':1e-3})
            ev(fit.x);local.append(dict(success=bool(fit.success),message=str(fit.message),iterations=int(fit.nit)))
        best=None
        for row in sorted([x for x in rows if x['healthy_kill_pct']<=100*LIMIT+1e-5],key=lambda x:x['tumor_kill_pct'],reverse=True):
            p=replace(base,**{k:row[k] for k in RATES});e,h=score(p,tight=True)
            if h<=LIMIT+1e-7:
                best=dict(row,tumor_kill_pct=100*e,healthy_kill_pct=100*h,parameters=asdict(p));break
        csvwrite(out/f'evaluations_seed{seed}.csv',rows)
        runs.append(dict(seed=seed,algorithm_nfe=alg.nfe,ode_evaluations=len(rows),local=local,selected=best))
        print('SEARCH',name,seed,len(rows),'best',None if best is None else [best['tumor_kill_pct'],best['healthy_kill_pct']],flush=True)
    selected=max([x['selected'] for x in runs if x['selected']],key=lambda x:x['tumor_kill_pct'],default=None)
    result=dict(name=name,rho1=rho1,rho2=rho2,ec50=ec50,runs=runs,selected=selected)
    dump(out/'result.json',result);return result

def trajectory(name,p):
    r=simulate(p,points=241,rtol=1e-9,atol=1e-12);rows=[]
    for i,t in enumerate(r['t']):
        rows.append(dict(time_h=float(t/3600),tumor_kill_pct=float(100*r['tumor'][i]),healthy_kill_pct=float(100*r['healthy'][i]),tumor_uptake_uM=float(1e6*r['intracellular'][i,0]),tumor_T1_RO_pct=float(100*r['occupancy'][i,0,0]),tumor_T2_RO_pct=float(100*r['occupancy'][i,0,1]),central_uM=float(r['y'][i,0])))
    csvwrite(OUT/(name+'_trajectory.csv'),rows);return r

def finalise(base,results):
    selected_ref=Parameters(**results['equal_reference']['selected']['parameters'])
    rows=[];checks=[];sens=[];traces={}
    for name,rho1,rho2,ec50 in CASES:
        result=results[name];sel=result['selected'];p0=replace(base,rho1=rho1,rho2=rho2,ec50=ec50,kon1=1e4,koff1=1,kon2=1e4,koff2=1)
        before=score(p0,tight=True);fixed=score(replace(selected_ref,rho1=rho1,rho2=rho2,ec50=ec50),tight=True)
        if sel is None:raise RuntimeError('No feasible design for '+name)
        p=Parameters(**sel['parameters']);r=trajectory(name,p);traces[name]=r
        rad=score(p,tight=True,method='Radau');bdf=(r['tumor'][-1],r['healthy'][-1])
        err=max(abs(rad[0]-bdf[0]),abs(rad[1]-bdf[1]));mass=float(max(abs(r['mass_error']))/p.dose)
        assert err<1e-6 and mass<1e-6 and max(r['healthy'])<.01
        exchanged=replace(p,rho1=p.rho2,rho2=p.rho1,kon1=p.kon2,kon2=p.kon1,koff1=p.koff2,koff2=p.koff1,k2d1=p.k2d2,k2d2=p.k2d1)
        swap=score(exchanged,tight=True);swaperr=max(abs(swap[0]-bdf[0]),abs(swap[1]-bdf[1]));assert swaperr<1e-6
        checks.append(dict(scenario=name,solver_error_fraction=float(err),mass_relative_error=mass,target_exchange_error_fraction=float(swaperr)))
        extrema=[]
        variants={'N_80pct':replace(p,total_cells=.8*p.total_cells),'N_120pct':replace(p,total_cells=1.2*p.total_cells),
          'volume_80pct':replace(p,cell_volume_l=.8*p.cell_volume_l),'volume_120pct':replace(p,cell_volume_l=1.2*p.cell_volume_l),
          'turnover_80pct':replace(p,turnover=.8*p.turnover),'turnover_120pct':replace(p,turnover=1.2*p.turnover),
          'original_sengers_off_rule':replace(p,surface_off_rule='sengers'),
          'background_CL_zero':replace(p,cl_l_h=0),'background_CL_double':replace(p,cl_l_h=2*p.cl_l_h)}
        for label,pp in variants.items():
            e,h=score(pp,tight=True);sens.append(dict(scenario=name,perturbation=label,tumor_kill_pct=100*e,healthy_kill_pct=100*h,below_1pct=h<.01));extrema.append(100*h)
        rates={k:getattr(p,k) for k in RATES};candidates=[x['selected']['tumor_kill_pct'] for x in result['runs'] if x['selected']]
        rows.append(dict(scenario=name,rho1=rho1,rho2=rho2,ec50_uM=ec50*1e6,**rates,KD1_uM=p.koff1/p.kon1*1e6,KD2_uM=p.koff2/p.kon2*1e6,
          paper_equal_rates_tumor_pct=100*before[0],paper_equal_rates_healthy_pct=100*before[1],
          transferred_equal_design_tumor_pct=100*fixed[0],transferred_equal_design_healthy_pct=100*fixed[1],
          optimized_tumor_pct=float(100*bdf[0]),optimized_healthy_pct=float(100*bdf[1]),seed_tumor_spread_pp=max(candidates)-min(candidates),sensitivity_max_healthy_pct=max(extrema),
          near_search_boundary=any(min(abs(np.log10(rates[k])-lo),abs(np.log10(rates[k])-hi))<.03 for k,(lo,hi) in zip(RATES,BOUNDS))))
    csvwrite(OUT/'summary.csv',rows);csvwrite(OUT/'sensitivity.csv',sens);dump(OUT/'verification.json',checks)
    baselines=[]
    for name,old in SCENARIOS.items():
        p=replace(base,**{k:getattr(old,k) for k in RATES},rho1=old.rho1,rho2=old.rho2,mode=old.mode)
        r=trajectory('published_'+name,p)
        baselines.append(dict(scenario=name,tumor_pct=float(100*r['tumor'][-1]),healthy_pct=float(100*r['healthy'][-1]),tumor_uptake_uM=float(1e6*r['intracellular'][-1,0]),tumor_T2_RO_pct=float(100*r['occupancy'][-1,0,1])))
    csvwrite(OUT/'published_candidates.csv',baselines)
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,3,figsize=(12,9),layout='constrained')
    for ax,(name,rho1,rho2,ec50) in zip(axes.flat,CASES):
        r=traces[name];ax.plot(r['t']/3600,100*r['tumor'],label='Tumour');ax.plot(r['t']/3600,100*r['healthy'],label='Healthy');ax.axhline(1,c='grey',lw=.6,ls=':');ax.set(title=f'{rho1:.0e} / {rho2:.0e}; EC50 {ec50*1e6:g} µM',xlabel='Time (h)',ylabel='Cell kill (%)',ylim=(0,100),xlim=(0,24))
    axes.flat[0].legend();fig.suptitle('Our calibrated model: independently optimised density scenarios\nNominal healthy kill ≤0.95%; model predictions, not paper data');fig.savefig(OUT/'optimized_timecourses.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(2,1,figsize=(11,7),layout='constrained');x=np.arange(len(rows))
    for a,key,label in [(ax[0],'tumor','Tumour kill (%)'),(ax[1],'healthy','Healthy kill (%)')]:
        a.bar(x-.2,[r[f'transferred_equal_design_{key}_pct'] for r in rows],width=.4,label='Equal-density design transferred')
        a.bar(x+.2,[r[f'optimized_{key}_pct'] for r in rows],width=.4,label='Density-specific optimisation');a.set_ylabel(label);a.set_xticks(x,[r['scenario'].replace('_','\n') for r in rows],fontsize=8)
    ax[0].legend();ax[1].axhline(1,c='red',ls='--');fig.savefig(OUT/'density_reoptimization_comparison.png',dpi=170);plt.close(fig)
    dump(OUT/'results.json',results);print('COMPLETE',json.dumps(rows),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--nfe',type=int,default=384);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seeds',type=int,nargs='+',default=[8407,2026]);ap.add_argument('--resume',action='store_true');a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    base=Parameters(**json.loads((OUT/'configuration.json').read_text())['parameters']) if a.resume and (OUT/'configuration.json').exists() else calibrate()
    results={};jobs=[]
    for name,r1,r2,ec in CASES:
        path=OUT/name/'result.json'
        if a.resume and path.exists():results[name]=json.loads(path.read_text())
        else:jobs.append((name,r1,r2,ec,asdict(base),a.nfe,a.seeds))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for fut in as_completed([pool.submit(run_case,j) for j in jobs]):
            result=fut.result();results[result['name']]=result
    finalise(base,results)
if __name__=='__main__':main()
