# Chain- and disorder-resolved renders: PARP9 in blues, DTX3L in magentas,
# with the intrinsically disordered regions saturated against pale folded cores.
# Companion to render_states.pml / render_box.pml, which colour by domain identity.
python
from pymol import cmd
import csv

# Four steps, all at distinct lightness so no cross-chain pair rests on hue
# alone -- blue vs magenta at matched lightness is a red-green confusion pair.
# Validated with the dataviz validator (--mode light --pairs all):
# worst CVD dE 14.3 (protan), worst normal-vision dE 20.2.
P9_FOLD,  P9_DIS  = '0x63A2FE', '0x0052B4'      # PARP9  pale / deep blue
DTX_FOLD, DTX_DIS = '0xCC629C', '0x8F1563'      # DTX3L  pale / deep magenta
XL_DASH = '0x1A1A1A'                            # near-black: reads as a measurement, not an entity

# Restrained blocks from binding/p9_dtx3l/input/domains.yaml. Everything outside
# them is unrestrained in CALVADOS, i.e. the disordered regions.
FOLDED = {'A': [(60,96),(107,296),(311,493),(516,539),(543,622),(644,823)],
          'B': [(11,87),(138,201),(233,302),(305,362),(378,447),(449,507),
                (556,607),(608,737)]}
CHAIN = {'A': ('PARP9', P9_FOLD, P9_DIS), 'B': ('DTX3L', DTX_FOLD, DTX_DIS)}

def style(ob, cell=False):
    cmd.hide('everything', ob); cmd.show('cartoon', ob); cmd.cartoon('tube', ob)
    for ch, (_, fold, dis) in CHAIN.items():
        # paint the whole chain disordered, then over-paint the folded blocks
        cmd.color(dis, f'{ob} and chain {ch}')
        for lo, hi in FOLDED[ch]:
            cmd.color(fold, f'{ob} and chain {ch} and resi {lo}-{hi}')
    if cell: draw_box()

BOX = 400.0                                     # A; CRYST1 of every box PDB (40 nm)
BOX_RGB = (0.45, 0.45, 0.45)

def draw_box(L=BOX):
    """Draw the periodic cell as an explicit CGO.

    `show cell` ignores cell_color in this PyMOL build -- it came out yellow for
    one object and salmon for another -- so the edges are built by hand.
    """
    from pymol.cgo import CYLINDER
    c = [(x, y, z) for x in (0, L) for y in (0, L) for z in (0, L)]
    obj = []
    for i in range(8):
        for j in range(i + 1, 8):
            # an edge differs in exactly one coordinate
            if sum(a != b for a, b in zip(c[i], c[j])) == 1:
                obj += [CYLINDER, *c[i], *c[j], 0.9, *BOX_RGB, *BOX_RGB]
    cmd.load_cgo(obj, 'box')


def xlinks(ob, csv_row):
    n = 0
    for pair in (csv_row or '').split(';'):
        if not pair: continue
        a, b = pair.split('_') if '_' in pair else pair.replace('K','').split('-')
        d = f'{ob}_d{a}_{b}'
        cmd.distance(d, f'{ob} and chain A and resi {a} and name CA',
                        f'{ob} and chain B and resi {b} and name CA')
        cmd.color(XL_DASH, d); cmd.set('dash_radius', 0.4, d); cmd.hide('labels', d)
        cmd.show('spheres', f'{ob} and chain A and resi {a} and name CA')
        cmd.show('spheres', f'{ob} and chain B and resi {b} and name CA')
        n += 1
    return n

cmd.bg_color('white'); cmd.set('ray_opaque_background', 1)
cmd.set('cartoon_trace_atoms', 1); cmd.set('sphere_scale', 0.55)
cmd.set('dash_gap', 0.35); cmd.set('dash_width', 3)
cmd.set('ray_shadows', 0)

# (pdb, states csv, state key, output stem, draws cell?, tube radius)
JOBS = [
  ('p9_dtx3l_xl_h20_best.pdb', 'p9_dtx3l_xl_h20_states.csv',  'best', 'restrained_best',    False, 1.1),
  ('p9_dtx3l_best.pdb',        'p9_dtx3l_states.csv',         'best', 'unrestrained_best',  False, 1.1),
  ('docked_rep1_t050.pdb',     'docked_rep1_timecourse.csv',  't050', 'hyres_docked_250ns', False, 1.1),
  ('p9_dtx3l_rep1_t000.pdb',   'p9_dtx3l_rep1_timecourse.csv','t000', 'setup_separated',    True,  1.4),
  ('raw_dock_t0.pdb',          None,                          None,   'setup_docked',       True,  1.4),
]
for pdb, csvf, key, out, cell, rad in JOBS:
    cmd.delete('all')
    cmd.set('cartoon_tube_radius', rad)
    cmd.load(pdb, 'S'); style('S', cell=cell)
    sat = ''
    if csvf:
        try:
            for r in csv.DictReader(open(csvf)):
                if r.get('state') == key:
                    sat = r.get('crosslinks') or r.get('sat') or ''
        except Exception: pass
    n = xlinks('S', sat)
    if cell:
        # frame the CELL, not the molecules: 8 corner pseudoatoms fix the extent,
        # otherwise orient/zoom crops the box away
        for x in (0, 400):
            for y in (0, 400):
                for z in (0, 400):
                    cmd.pseudoatom('corners', pos=[x, y, z])
        cmd.hide('everything', 'corners')
        cmd.reset(); cmd.zoom('corners', 95)
        cmd.turn('x', -14); cmd.turn('y', 22)
        cmd.ray(1300, 1100)
    else:
        cmd.orient('S')
        cmd.ray(1400, 1000)
    cmd.png(f'chain_{out}.png', dpi=150)
    print(f'RENDER chain_{out}: {n} crosslinks drawn')
python end
