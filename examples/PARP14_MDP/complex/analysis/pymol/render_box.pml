python
from pymol import cmd
DOM = {'A': [('MD1',107,296,'0x6A3FA0'),('MD2',311,493,'0x9169C7'),
             ('KH2',543,622,'0x1CC497'),('ART',644,823,'0x9C5C45'),
             ('1/2KH1a',60,96,'0x16A07B'),('1/2KH1b',516,539,'0x16A07B')],
       'B': [('RRM',11,87,'0xB88C13'),('KH1',138,201,'0x117C60'),
             ('KH2',233,302,'0x16A07B'),('KH3',305,362,'0x1CC497'),
             ('KH4',378,447,'0x16A07B'),('KH5',449,507,'0x117C60'),
             ('RING',556,607,'0xC8481A'),('DTC',608,737,'0x9C4090')]}
cmd.bg_color('white'); cmd.set('ray_opaque_background',1)
cmd.set('cartoon_trace_atoms',1); cmd.set('cartoon_tube_radius',1.4)
cmd.set('ray_shadows',0); cmd.set('cell_color','grey60')

for pdb, out, lab in (('p9_dtx3l_rep1_t000.pdb','setup_separated','t = 0 ns, chains 17 nm apart'),
                      ('raw_dock_t0.pdb','setup_docked','t = 0 ns, AF3 docked pose')):
    cmd.delete('all')
    cmd.load(pdb,'S')
    cmd.hide('everything','S'); cmd.show('cartoon','S'); cmd.cartoon('tube','S')
    cmd.color('grey75','S')
    for ch,doms in DOM.items():
        for nm,lo,hi,col in doms:
            cmd.color(col, f'S and chain {ch} and resi {lo}-{hi}')
    cmd.show('cell','S')
    # zoom on the CELL, not the molecules: 8 corner pseudoatoms define the extent,
    # otherwise orient/zoom frames the protein and crops the box away
    for i,(x,y,z) in enumerate([(a,b,c) for a in (0,400) for b in (0,400) for c in (0,400)]):
        cmd.pseudoatom('corners', pos=[x,y,z])
    cmd.hide('everything','corners')
    cmd.reset()
    cmd.zoom('corners', 95)
    cmd.turn('x',-14); cmd.turn('y',22)
    cmd.ray(1300,1100); cmd.png(f'render_{out}.png', dpi=150)
    print(f'RENDER {out}')
python end
