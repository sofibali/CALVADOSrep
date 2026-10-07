load p9_dtx3l_xl_h20_crosslink_states.pse
python
from pymol import cmd
print('all names:', cmd.get_names('all')[:12])
print('group children of best_xlinks:', cmd.get_object_list('best_xlinks') if 'best_xlinks' in cmd.get_names('all') else 'no such group')
python end
