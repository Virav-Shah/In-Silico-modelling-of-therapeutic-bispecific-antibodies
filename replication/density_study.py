"""Section 3.4 study: fixed-candidate density sweeps, then unequal-density search.
Reported settings are preserved. Unreported physiology is frozen, not re-fitted.
Search bounds are inferred from Fig.5 axes, NOT reported optimiser settings.
"""
import argparse,csv,json,itertools,random,os
from pathlib import Path
from dataclasses import replace,asdict
import numpy as np
from model import Parameters,SCENARIOS,simulate,densities
from platypus import Problem,Real,Solution,InjectedPopulation,NSGAII
ROOT=Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp'/'matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RATES=('kon1','koff1','kon2','koff2')
BOUNDS=[(4,7),(-3,0),(4,7),(-3,0)]
LIMIT=.01-1e-10

def write(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def endpoint(p,tight=False):
    r=simulate(p,points=3,rtol=1e-9 if tight else 1e-7,atol=1e-12 if tight else 1e-10)
    return float(100*r['tumor'][-1]),float(100*r['healthy'][-1])

def grid(out,name,base,n,lo,hi):
    axis=np.logspace(lo,hi,n);e=np.zeros((n,n));h=e.copy();rows=[]
    for j,y in enumerate(axis):
        for i,x in enumerate(axis):
            p=replace(base,rho1=float(x),rho2=float(y));a,b=endpoint(p);e[j,i]=a;h[j,i]=b
            rows.append(dict(tumor_T1_mol_dm2=float(x),tumor_T2_mol_dm2=float(y),healthy_T1_mol_dm2=.25*x,healthy_T2_mol_dm2=.25*y,tumor_kill_pct=a,healthy_kill_pct=b))
        if j%5==0:print(name,'row',j+1,'/',n,flush=True)
    write(out/(name+'.csv'),rows)
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    for ax,z,cmap,label,levels,maximum in [(axes[0],e,'Greens','Tumor kill (%)',[15,30,45,60,75,90],100),(axes[1],h,'Reds','Healthy kill (%)',[1.5,3,4.5,6,7.5],max(8,float(h.max())))]:
        m=ax.pcolormesh(axis,axis,z,shading='nearest',cmap=cmap,vmin=0,vmax=maximum)
        cs=ax.contour(axis,axis,z,levels=[v for v in levels if z.min()<v<z.max()],colors='black',linewidths=.7)
        ax.clabel(cs,inline=True,fontsize=8,fmt='%g');fig.colorbar(m,ax=ax,label=label)
        ax.set(xscale='log',yscale='log',xlabel='Tumor T1 density (mol/dm²)',ylabel='Tumor T2 density (mol/dm²)',title=label)
        if lo<=-10<=hi:ax.plot(1e-10,1e-10,'wo',ms=5,mec='black')
        ax.plot(1e-9,5e-11,'w^',ms=6,mec='black')
    fig.suptitle(f'{name}: fixed paper Figure 6 rates, {n} × {n} computed grid\nReconstruction; no rate reoptimisation at individual pixels',fontsize=11)
    fig.savefig(out/(name+'.png'),dpi=180);plt.close(fig)
    sym=max(float(np.max(abs(e-e.T))),float(np.max(abs(h-h.T))))
    assert sym<1e-4,('exchange symmetry',sym)
    return dict(points=n*n,max_exchange_asymmetry_percentage_points=sym,tumor_range_pct=[float(e.min()),float(e.max())],healthy_range_pct=[float(h.min()),float(h.max())])

def trajectory(out,name,p):
    r=simulate(p,points=241,rtol=1e-9,atol=1e-12)
    rows=[]
    for i,t in enumerate(r['t']):
        row=dict(time_h=t/3600,tumor_kill_pct=100*r['tumor'][i],healthy_kill_pct=100*r['healthy'][i],tumor_uptake_uM=1e6*r['intracellular'][i,0],healthy_all_uptake_uM=1e6*float(r['intracellular'][i,1:]@np.array([.1,.1,.01])))
        row['healthy_positive_uptake_uM']=row['healthy_all_uptake_uM']/.21
        for j,label in enumerate(['tumor','healthy_T1','healthy_T2','healthy_T12']):
            for k in range(2):row[f'{label}_T{k+1}_occupancy_pct']=100*r['occupancy'][i,j,k]
        rows.append(row)
    write(out/(name+'_trajectory.csv'),rows)
    fig,ax=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
    for j,label in enumerate(['Tumor (T2)','Healthy T1 (T1)','Healthy T2 (T2)','Healthy T12 (T2)']):
        target=0 if j==1 else 1
        ax[0].plot(r['t']/3600,100*r['occupancy'][:,j,target],label=label)
    ax[1].plot(r['t']/3600,1e6*r['intracellular'][:,0],label='Tumor')
    ax[1].plot([v['time_h'] for v in rows],[v['healthy_positive_uptake_uM'] for v in rows],label='Healthy target-positive mean')
    ax[2].plot(r['t']/3600,100*r['tumor'],label='Tumor');ax[2].plot(r['t']/3600,100*r['healthy'],label='All healthy')
    for a,label in zip(ax,['Target occupancy (%)','Internalised ADC (µM)','Cell kill (%)']):
        a.set(xlabel='Time (h)',ylabel=label,xlim=(0,24));a.legend(fontsize=7)
    ax[0].set_ylim(0,100);ax[2].set_ylim(0,100)
    fig.suptitle(name+' | frozen physiology; paper dose, potency and population fractions',fontsize=10)
    fig.savefig(out/(name+'.png'),dpi=180);plt.close(fig)
    return r,dict(name=name,**{k:float(getattr(p,k)) for k in RATES},rho1=p.rho1,rho2=p.rho2,tumor_kill_pct=float(100*r['tumor'][-1]),healthy_kill_pct=float(100*r['healthy'][-1]),tumor_uptake_uM=float(1e6*r['intracellular'][-1,0]),tumor_T2_occupancy_pct=float(100*r['occupancy'][-1,0,1]))

def optimize(out,base,nfe,seed):
    random.seed(seed);np.random.seed(seed)
    problem=Problem(4,2,1);problem.types[:]=[Real(*b) for b in BOUNDS];problem.directions[:]=[Problem.MINIMIZE]*2;problem.constraints[:]='<=0'
    rows=[]
    def fn(z):
        rates=10**np.asarray(z);e,h=endpoint(replace(base,**dict(zip(RATES,rates))))
        rows.append(dict(seed=seed,**{k:float(v) for k,v in zip(RATES,rates)},tumor_kill_pct=e,healthy_kill_pct=h))
        if len(rows)%200==0:print('search',seed,len(rows),'min healthy',min(r['healthy_kill_pct'] for r in rows),flush=True)
        return [-e/100,h/100],[h/100-LIMIT]
    problem.function=fn;seeds=[]
    for values in itertools.product((4.,7.),(-3.,0.),(4.,7.),(-3.,0.)):
        s=Solution(problem);s.variables[:]=values;seeds.append(s)
    for values in ([4,0,np.log10(5e4),0],[4,0,4,0]):
        s=Solution(problem);s.variables[:]=values;seeds.append(s)
    algorithm=NSGAII(problem,population_size=48,generator=InjectedPopulation(seeds));algorithm.run(nfe)
    write(out/f'figure8_search_seed{seed}.csv',rows)
    best=None
    for row in sorted((r for r in rows if r['healthy_kill_pct']<100*LIMIT),key=lambda r:r['tumor_kill_pct'],reverse=True):
        pp=replace(base,**{k:row[k] for k in RATES});e,h=endpoint(pp,tight=True)
        if h<100*LIMIT:best=dict(row,tumor_kill_pct=e,healthy_kill_pct=h,parameters=asdict(pp));break
    return dict(seed=seed,algorithm_nfe=algorithm.nfe,evaluator_calls=len(rows),best_feasible=best,minimum_healthy_observed=min(rows,key=lambda r:r['healthy_kill_pct']))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--profile',choices=['calibrated','original'],default='calibrated');ap.add_argument('--grid',type=int,default=31);ap.add_argument('--nfe',type=int,default=960);ap.add_argument('--seeds',nargs='+',type=int,default=[8407,2026]);a=ap.parse_args()
    if a.grid<3 or a.nfe<48:ap.error('grid >= 3 and nfe >= 48 required')
    out=ROOT/'output'/'density_study'/a.profile;out.mkdir(parents=True,exist_ok=True)
    config={} if a.profile=='original' else json.loads((ROOT/'output/prototype/configuration.json').read_text())['calibration_inputs']
    fixed=replace(SCENARIOS['figure6'],**config)
    unequal=replace(fixed,rho1=1e-9,rho2=5e-11)
    paper8=replace(unequal,kon2=5e4)
    manifest=dict(profile=a.profile,frozen_parameters=asdict(fixed),density_rule='Each healthy expressed target remains 25% of corresponding tumor target; population frequencies unchanged',scoring_time_h=24,figure7_log10_density_range=[-12,-8],S4_log10_density_range=[-14,-9],grid_per_axis=a.grid,optimizer_log10_bounds=BOUNDS,bounds_provenance='Inferred from Figure 5 axes; original optimiser bounds unreported',optimizer='NSGAII, population 48; original algorithm unreported',seeds=a.seeds,nfe_per_seed=a.nfe,healthy_fraction_constraint='strictly < 0.01',source='Stoddard et al. 2026 DOI 10.20935/AcadDrug8407, sections 2.3/2.4/3.4, Table 1/2, Figs 7/8/S4')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    summaries=[];traces={}
    for name,p in [('equal_reference',fixed),('unequal_before_reoptimization',unequal),('unequal_paper_table2',paper8)]:
        traces[name],row=trajectory(out,name,p);summaries.append(row)
    write(out/'candidate_comparison.csv',summaries)
    g7=grid(out,'figure7',fixed,a.grid,-12,-8);s4=grid(out,'figureS4',fixed,a.grid,-14,-9)
    # Explicit factor-of-ten cross-sections distinguish absolute density and ratio.
    rows=[]
    for rho1 in [1e-11,1e-10,1e-9]:
        for ratio in [.1,1,10]:
            pp=replace(fixed,rho1=rho1,rho2=rho1*ratio);e,h=endpoint(pp,tight=True)
            rows.append(dict(rho1=rho1,rho2=rho1*ratio,T2_T1_ratio=ratio,tumor_kill_pct=e,healthy_kill_pct=h))
    write(out/'tenfold_density_checks.csv',rows)
    searches=[optimize(out,unequal,a.nfe,s) for s in a.seeds]
    feasible=[s['best_feasible'] for s in searches if s['best_feasible']]
    selected=max(feasible,key=lambda r:r['tumor_kill_pct']) if feasible else None
    if selected:
        traces['unequal_search_selected'],row=trajectory(out,'unequal_search_selected',Parameters(**selected['parameters']));summaries.append(row);write(out/'candidate_comparison.csv',summaries)
    # Cross-solver verification at 3 fixed candidates and 4 density corners.
    checks=[]
    validation=[('equal_reference',fixed),('unequal_before',unequal),('paper_table2',paper8)]+[(f'corner_{x}_{y}',replace(fixed,rho1=x,rho2=y)) for x,y in itertools.product([1e-12,1e-8],repeat=2)]
    for name,p in validation:
        b=simulate(p,points=3,rtol=1e-9,atol=1e-12);r=simulate(p,points=3,rtol=1e-9,atol=1e-12,method='Radau')
        err=max(abs(b['tumor'][-1]-r['tumor'][-1]),abs(b['healthy'][-1]-r['healthy'][-1]));mass=max(abs(b['mass_error']))/p.dose
        assert err<1e-6 and mass<1e-6,(name,err,mass)
        checks.append(dict(case=name,endpoint_solver_error_fraction=float(err),mass_relative_error=float(mass)))
    result=dict(figure7=g7,S4=s4,searches=searches,selected=selected,candidates=summaries,checks=checks)
    (out/'results.json').write_text(json.dumps(result,indent=2))
    print('COMPLETE',a.profile,'selected',selected,flush=True)
if __name__=='__main__':main()
