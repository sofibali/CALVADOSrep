#!/usr/bin/env python3
"""Convergence verdict + dissociation kinetics + bound-window contact maps.

Since no run stays associated, a plain time-averaged contact map is dominated by
frames where the chains are metres apart. Dissociated frames contribute exactly
zero contacts, so the conditional (bound-window) map is recovered exactly as
   map_bound = map_all / frac_frames_bound
without re-reading any trajectory.
"""
import itertools
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

OUT=Path(__file__).resolve().parent/'analysis'
FIG=OUT/'figures'; FIG.mkdir(exist_ok=True)
NS=0.5                                  # ns per analysed sample
BOUND_NM=1.0                            # min CA-CA distance defining "bound"
plt.rcParams.update({'font.size':8,'figure.dpi':150,'savefig.bbox':'tight',
                     'axes.spines.top':False,'axes.spines.right':False})
store=np.load(OUT/'data'/'raw.npy',allow_pickle=True)[0]

def surv_time(mind):
    """ns until the pair leaves contact for good (never returns below BOUND_NM)."""
    b=mind<BOUND_NM
    if not b.any(): return 0.0
    return float((np.max(np.where(b)[0])+1)*NS)

rows=[]; maps={}
for key,reps in store.items():
    system,variant=key.split('|')
    chains=reps[0][1]['chains']
    for p in reps[0][1]['ncon'].keys():
        i,j=p; pair=f'{chains[i]}-{chains[j]}'
        md=[r[1]['mind'][p] for r in reps]
        cm=[r[1]['cmaps'][p] for r in reps]
        fb=np.array([(m<BOUND_NM).mean() for m in md])          # bound fraction per rep
        st=np.array([surv_time(m) for m in md])
        nfin=int((fb>0.98).sum())                                # reps still bound throughout
        # conditional bound-window map, pooled over reps that spend any time bound
        sel=[(c,f) for c,f in zip(cm,fb) if f>0.02]
        bmap=(np.sum([c/f for c,f in sel],axis=0)/len(sel)) if sel else None
        rows.append(dict(system=system,variant=variant,pair=pair,n_rep=len(reps),
            bound_frac_mean=round(float(fb.mean()),3),
            bound_frac_range=f'{fb.min():.2f}-{fb.max():.2f}',
            survival_ns_median=round(float(np.median(st)),1),
            survival_ns_range=f'{st.min():.0f}-{st.max():.0f}',
            reps_bound_throughout=nfin,
            final_min_dist_nm=round(float(np.mean([m[-1] for m in md])),1)))
        if bmap is not None:
            maps[(system,variant,pair)]=dict(map=bmap,ci=chains[i],cj=chains[j],
                bi=reps[0][1]['blocks'][chains[i]],bj=reps[0][1]['blocks'][chains[j]],
                nsel=len(sel))
df=pd.DataFrame(rows).sort_values(['pair','system','variant'])
df.to_csv(OUT/'dissociation_table.csv',index=False)
print(df.to_string(index=False))
print(f"\nruns where ALL replicates stayed bound: {int((df.reps_bound_throughout==df.n_rep).sum())} / {len(df)}")

# ---- Fig A: min inter-chain distance traces ----
systems=sorted(df.system.unique())
fig,axes=plt.subplots(len(systems),1,figsize=(8,2.3*len(systems)),squeeze=False)
for ax,syst in zip(axes[:,0],systems):
    for key,reps in store.items():
        s,v=key.split('|')
        if s!=syst: continue
        chains=reps[0][1]['chains']
        for p in reps[0][1]['mind'].keys():
            arr=np.array([r[1]['mind'][p] for r in reps],float)
            t=np.arange(arr.shape[1])*NS
            ax.plot(t,np.median(arr,0),lw=0.8,
                    label=f"{v} · {chains[p[0]]}-{chains[p[1]]}")
    ax.axhline(BOUND_NM,color='k',ls='--',lw=0.8)
    ax.set_yscale('log'); ax.set_ylabel('min CA-CA (nm)'); ax.set_xlabel('time (ns)')
    ax.set_title(f'{syst} — median inter-chain distance; dashed = 1 nm contact cutoff',loc='left')
    ax.legend(fontsize=4.5,ncol=3,frameon=False)
fig.tight_layout(); fig.savefig(FIG/'01_dissociation_traces.png'); plt.close(fig)

# ---- Fig B: bound fraction + survival scorecard ----
d=df.copy(); d['label']=d.system+' · '+d.variant+' · '+d['pair']
fig,ax=plt.subplots(figsize=(7,0.25*len(d)+1.4))
y=np.arange(len(d))
ax.barh(y,d.bound_frac_mean,color=np.where(d.bound_frac_mean>0.5,'#2e7d32','#c62828'))
ax.set_yticks(y); ax.set_yticklabels(d.label,fontsize=5); ax.invert_yaxis()
ax.axvline(0.5,color='k',ls='--',lw=0.8)
ax.set_xlabel('fraction of the run spent in contact (min CA-CA < 1 nm)')
ax.set_title('Bound fraction per set. Dashed line = the 50% bar for a usable contact map',
             loc='left',fontsize=8)
for i,(f,s) in enumerate(zip(d.bound_frac_mean,d.survival_ns_median)):
    ax.text(f+0.005,i,f'  {f:.2f} · lost by {s:.0f} ns',va='center',fontsize=4.5)
fig.tight_layout(); fig.savefig(FIG/'02_bound_fraction.png'); plt.close(fig)
np.save(OUT/'data'/'bound_maps.npy',np.array([maps],dtype=object),allow_pickle=True)
print('saved bound maps:',len(maps))
