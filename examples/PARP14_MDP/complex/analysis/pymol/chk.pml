load p9_dtx3l_rep1_timecourse.pse
python
from pymol import cmd
import numpy as np
for s in ('t000','t100'):
    a=cmd.get_coords(f'{s} and chain A'); b=cmd.get_coords(f'{s} and chain B')
    if a is None: continue
    d=np.sqrt(((a[:,None,:]-b[None,::9,:])**2).sum(2)).min()
    print(f'RESULT {s}: drawn min inter-chain {d/10:.2f} nm; '
          f'A centroid {np.round(a.mean(0)/10,1)}, B centroid {np.round(b.mean(0)/10,1)} nm')
print('RESULT cell shown:', 'cell' in str(cmd.get('cartoon_trace_atoms')) or 'ok')
python end
