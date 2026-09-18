from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
OUT=Path('/tmp/claude-64170/-home-vle-Vincent/a5f0a8dd-a59a-412d-b58f-35aebe7f4fab/scratchpad/conv')
FIG=OUT/'figures'
maps=np.load(OUT/'data'/'bound_maps.npy',allow_pickle=True)[0]
plt.rcParams.update({'font.size':8,'figure.dpi':150,'savefig.bbox':'tight'})

def heat(ax,M,bi,bj,title,cmap='magma_r',norm=None,lab='contacts / frame while bound'):
    im=ax.imshow(M,cmap=cmap,aspect='auto',norm=norm)
    ax.set_xticks(range(len(bj))); ax.set_xticklabels(bj,rotation=90,fontsize=5)
    ax.set_yticks(range(len(bi))); ax.set_yticklabels(bi,fontsize=5)
    ax.set_title(title,fontsize=6.5,loc='left')
    cb=plt.colorbar(im,ax=ax,fraction=0.046); cb.ax.tick_params(labelsize=5); cb.set_label(lab,fontsize=5)

V='prod_dom-all_kgo4'
panels=[('parp14_parp9_dtx3l',V,'parp14-parp9'),('parp14_parp9_dtx3l',V,'parp14-dtx3l'),
        ('parp14_parp9_dtx3l',V,'parp9-dtx3l'),('parp9_dtx3l',V,'parp9-dtx3l')]
panels=[p for p in panels if p in maps]
fig,axes=plt.subplots(2,2,figsize=(9.5,8))
for ax,k in zip(axes.ravel(),panels):
    v=maps[k]
    heat(ax,v['map'],v['bi'],v['bj'],f"{k[0]}\n{k[1]} · {k[2]}  (n={v['nsel']} reps)")
    ax.set_ylabel(v['ci'],fontsize=6); ax.set_xlabel(v['cj'],fontsize=6)
for ax in axes.ravel()[len(panels):]: ax.axis('off')
fig.suptitle('Inter-chain domain contacts, conditioned on the chains being bound\n'
             '(no run stays bound for its full 100 ns — these are the interfaces as they exist before dissociation)',
             fontsize=8.5,x=0.01,ha='left')
fig.tight_layout(); fig.savefig(FIG/'03_bound_contact_maps.png'); plt.close(fig)

comps=[(('parp9_dtx3l',V,'parp9-dtx3l'),('parp14_parp9_dtx3l',V,'parp9-dtx3l'),
        'PARP9–DTX3L: inside the trimer  −  on its own'),
       (('parp14_parp9_dtx3l','output_0.002dt_no-restraint','parp14-dtx3l'),
        ('parp14_parp9_dtx3l','output_0.002dt_go-kh','parp14-dtx3l'),
        'PARP14–DTX3L: Gō-KH restrained  −  unrestrained'),
       (('dtx3l_dtx3l','dtx3l_dtx3l_box-90','dtx3l-dtx3l'),
        ('dtx3l_dtx3l','dtx3l_dtx3l_box-30','dtx3l-dtx3l'),
        'DTX3L homodimer: 30 nm box  −  90 nm box')]
comps=[c for c in comps if c[0] in maps and c[1] in maps]
fig,axes=plt.subplots(1,len(comps),figsize=(5.0*len(comps),4.2),squeeze=False)
for ax,(ka,kb,title) in zip(axes[0],comps):
    A,B=maps[ka],maps[kb]
    if A['map'].shape!=B['map'].shape: ax.axis('off'); continue
    D=B['map']-A['map']; lim=np.abs(D).max() or 1.0
    heat(ax,D,B['bi'],B['bj'],title,cmap='RdBu_r',
         norm=TwoSlopeNorm(vcenter=0,vmin=-lim,vmax=lim),lab='Δ contacts / frame')
    ax.set_ylabel(B['ci'],fontsize=6); ax.set_xlabel(B['cj'],fontsize=6)
fig.suptitle('Difference maps — red = more contact in the first-named condition',fontsize=8.5,x=0.01,ha='left')
fig.tight_layout(); fig.savefig(FIG/'04_difference_maps.png'); plt.close(fig)

# top contacts text summary
print("TOP BOUND-WINDOW CONTACTS")
for k in panels:
    v=maps[k]; M=v['map']
    idx=np.dstack(np.unravel_index(np.argsort(-M,axis=None),M.shape))[0][:5]
    print(f"\n{k[0]} / {k[2]}:")
    for a,b in idx:
        if M[a,b]<=0: break
        print(f"   {v['ci']}:{v['bi'][a]:<10s} – {v['cj']}:{v['bj'][b]:<10s} {M[a,b]:7.2f}")
print(f"\nfigures -> {FIG}")
