load p9_dtx3l_xl_h20_crosslink_states.pse
python
from pymol import cmd
objs=cmd.get_object_list()
struct=[o for o in objs if not o.startswith(('best_xl','mode1_xl','mode2_xl','mode3_xl','mode4_xl'))]
xl=[o for o in objs if '_xl_' in o]
print('STRUCTURES:', struct)
print('crosslink distance objects:', len(xl))
print('enabled on open:', [o for o in struct if cmd.get_vis().get(o)])
for s in ('best','mode1'):
    on=[o for o in xl if o.startswith(s+'_xl_') and cmd.get_vis().get(o)]
    print(f'  {s}: {len(on)} satisfied crosslinks shown')
python end
