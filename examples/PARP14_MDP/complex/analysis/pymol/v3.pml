load p9_dtx3l_crosslink_states.pse
python
from pymol import cmd
names=cmd.get_names('all'); vis=cmd.get_vis()
struct=cmd.get_object_list()
xl=[n for n in names if '_xl_' in n]
print(f'RESULT structures={len(struct)} {struct}')
print(f'RESULT distance objects={len(xl)}')
for s in struct:
    on=[n for n in xl if n.startswith(s+'_xl_') and vis.get(n)]
    tot=[n for n in xl if n.startswith(s+'_xl_')]
    print(f'RESULT   {s}: {len(on)} shown of {len(tot)} drawn')
python end
