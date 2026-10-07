# Publication renders of the key states. Domain colours match the figures.
python
from pymol import cmd

DOM = {'A': [('1/2KH1a',60,96,'0x16A07B'),('MD1',107,296,'0x6A3FA0'),
             ('MD2',311,493,'0x9169C7'),('1/2KH1b',516,539,'0x16A07B'),
             ('KH2',543,622,'0x1CC497'),('ART',644,823,'0x9C5C45')],
       'B': [('RRM',11,87,'0xB88C13'),('KH1',138,201,'0x117C60'),
             ('KH2',233,302,'0x16A07B'),('KH3',305,362,'0x1CC497'),
             ('KH4',378,447,'0x16A07B'),('KH5',449,507,'0x117C60'),
             ('RING',556,607,'0xC8481A'),('DTC',608,737,'0x9C4090')]}

def style(ob, cell=False):
    cmd.hide('everything', ob); cmd.show('cartoon', ob); cmd.cartoon('tube', ob)
    cmd.color('grey75', ob)
    for ch, doms in DOM.items():
        for nm, lo, hi, col in doms:
            cmd.color(col, f'{ob} and chain {ch} and resi {lo}-{hi}')
    if cell: cmd.show('cell', ob)

def xlinks(ob, csv_row):
    n = 0
    for pair in (csv_row or '').split(';'):
        if not pair: continue
        a, b = pair.split('_') if '_' in pair else pair.replace('K','').split('-')
        d = f'{ob}_d{a}_{b}'
        cmd.distance(d, f'{ob} and chain A and resi {a} and name CA',
                        f'{ob} and chain B and resi {b} and name CA')
        cmd.color('limon', d); cmd.set('dash_radius', 0.4, d)
        cmd.hide('labels', d)
        cmd.show('spheres', f'{ob} and chain A and resi {a} and name CA')
        cmd.show('spheres', f'{ob} and chain B and resi {b} and name CA')
        n += 1
    return n

import csv
cmd.bg_color('white'); cmd.set('ray_opaque_background', 1)
cmd.set('cartoon_trace_atoms', 1); cmd.set('cartoon_tube_radius', 1.1)
cmd.set('sphere_scale', 0.55); cmd.set('dash_gap', 0.35); cmd.set('dash_width', 3)
cmd.set('cell_color', 'grey70'); cmd.set('ray_shadows', 0)

JOBS = [
  ('p9_dtx3l_xl_h20_best.pdb',  'p9_dtx3l_xl_h20_states.csv', 'best',
   'restrained_best', 'CALVADOS +3 XL restrained - 18 of 36 crosslinks satisfied', False),
  ('p9_dtx3l_best.pdb',         'p9_dtx3l_states.csv',        'best',
   'unrestrained_best', 'CALVADOS unrestrained - best frame, 6 of 36', False),
  ('docked_rep1_t050.pdb',      'docked_rep1_timecourse.csv', 't050',
   'hyres_docked_250ns', 'HyRes matched docked, 250 ns - 17 of 36', False),
  ('p9_dtx3l_rep1_t000.pdb',    'p9_dtx3l_rep1_timecourse.csv','t000',
   'setup_box', 'Setup: chains 17 nm apart in the 40 nm box', True),
]
for pdb, csvf, key, out, title, cell in JOBS:
    cmd.delete('all')
    cmd.load(pdb, 'S'); style('S', cell=cell)
    sat = ''
    try:
        for r in csv.DictReader(open(csvf)):
            if r.get('state') == key:
                sat = r.get('crosslinks') or r.get('sat') or ''
    except Exception: pass
    n = xlinks('S', sat)
    cmd.orient('S')
    if cell: cmd.zoom('S', 12)
    cmd.ray(1400, 1000)
    cmd.png(f'render_{out}.png', dpi=150)
    print(f'RENDER {out}: {n} crosslinks drawn')
python end
