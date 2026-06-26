# Load KH domain comparison PDBs
# Run: pymol load_kh_comparison.pml


# KH1
load kh_pdbs/KH1_AF2_reference.pdb, KH1_AF2
color yellow, KH1_AF2
load kh_pdbs/KH1_baseline_current.pdb, KH1_baseline_current
color gray60, KH1_baseline_current
align KH1_baseline_current, KH1_AF2
load kh_pdbs/KH1_tight_uniform_700.pdb, KH1_tight_u700
color red, KH1_tight_u700
align KH1_tight_u700, KH1_AF2
load kh_pdbs/KH1_tight_uniform_350.pdb, KH1_tight_u350
color salmon, KH1_tight_u350
align KH1_tight_u350, KH1_AF2
load kh_pdbs/KH1_medium_uniform_700.pdb, KH1_medium_u700
color blue, KH1_medium_u700
align KH1_medium_u700, KH1_AF2
load kh_pdbs/KH1_medium_uniform_350.pdb, KH1_medium_u350
color lightblue, KH1_medium_u350
align KH1_medium_u350, KH1_AF2
load kh_pdbs/KH1_loose_uniform_700.pdb, KH1_loose_u700
color green, KH1_loose_u700
align KH1_loose_u700, KH1_AF2
load kh_pdbs/KH1_loose_uniform_350.pdb, KH1_loose_u350
color palegreen, KH1_loose_u350
align KH1_loose_u350, KH1_AF2

# Group KH1
group KH1_group, KH1_AF2 KH1_baseline_current KH1_tight_u700 KH1_tight_u350 KH1_medium_u700 KH1_medium_u350 KH1_loose_u700 KH1_loose_u350

# KH2
load kh_pdbs/KH2_AF2_reference.pdb, KH2_AF2
color yellow, KH2_AF2
load kh_pdbs/KH2_baseline_current.pdb, KH2_baseline_current
color gray60, KH2_baseline_current
align KH2_baseline_current, KH2_AF2
load kh_pdbs/KH2_tight_uniform_700.pdb, KH2_tight_u700
color red, KH2_tight_u700
align KH2_tight_u700, KH2_AF2
load kh_pdbs/KH2_tight_uniform_350.pdb, KH2_tight_u350
color salmon, KH2_tight_u350
align KH2_tight_u350, KH2_AF2
load kh_pdbs/KH2_medium_uniform_700.pdb, KH2_medium_u700
color blue, KH2_medium_u700
align KH2_medium_u700, KH2_AF2
load kh_pdbs/KH2_medium_uniform_350.pdb, KH2_medium_u350
color lightblue, KH2_medium_u350
align KH2_medium_u350, KH2_AF2
load kh_pdbs/KH2_loose_uniform_700.pdb, KH2_loose_u700
color green, KH2_loose_u700
align KH2_loose_u700, KH2_AF2
load kh_pdbs/KH2_loose_uniform_350.pdb, KH2_loose_u350
color palegreen, KH2_loose_u350
align KH2_loose_u350, KH2_AF2

# Group KH2
group KH2_group, KH2_AF2 KH2_baseline_current KH2_tight_u700 KH2_tight_u350 KH2_medium_u700 KH2_medium_u350 KH2_loose_u700 KH2_loose_u350

# KH3
load kh_pdbs/KH3_AF2_reference.pdb, KH3_AF2
color yellow, KH3_AF2
load kh_pdbs/KH3_baseline_current.pdb, KH3_baseline_current
color gray60, KH3_baseline_current
align KH3_baseline_current, KH3_AF2
load kh_pdbs/KH3_tight_uniform_700.pdb, KH3_tight_u700
color red, KH3_tight_u700
align KH3_tight_u700, KH3_AF2
load kh_pdbs/KH3_tight_uniform_350.pdb, KH3_tight_u350
color salmon, KH3_tight_u350
align KH3_tight_u350, KH3_AF2
load kh_pdbs/KH3_medium_uniform_700.pdb, KH3_medium_u700
color blue, KH3_medium_u700
align KH3_medium_u700, KH3_AF2
load kh_pdbs/KH3_medium_uniform_350.pdb, KH3_medium_u350
color lightblue, KH3_medium_u350
align KH3_medium_u350, KH3_AF2
load kh_pdbs/KH3_loose_uniform_700.pdb, KH3_loose_u700
color green, KH3_loose_u700
align KH3_loose_u700, KH3_AF2
load kh_pdbs/KH3_loose_uniform_350.pdb, KH3_loose_u350
color palegreen, KH3_loose_u350
align KH3_loose_u350, KH3_AF2

# Group KH3
group KH3_group, KH3_AF2 KH3_baseline_current KH3_tight_u700 KH3_tight_u350 KH3_medium_u700 KH3_medium_u350 KH3_loose_u700 KH3_loose_u350

# KH4
load kh_pdbs/KH4_AF2_reference.pdb, KH4_AF2
color yellow, KH4_AF2
load kh_pdbs/KH4_baseline_current.pdb, KH4_baseline_current
color gray60, KH4_baseline_current
align KH4_baseline_current, KH4_AF2
load kh_pdbs/KH4_tight_uniform_700.pdb, KH4_tight_u700
color red, KH4_tight_u700
align KH4_tight_u700, KH4_AF2
load kh_pdbs/KH4_tight_uniform_350.pdb, KH4_tight_u350
color salmon, KH4_tight_u350
align KH4_tight_u350, KH4_AF2
load kh_pdbs/KH4_medium_uniform_700.pdb, KH4_medium_u700
color blue, KH4_medium_u700
align KH4_medium_u700, KH4_AF2
load kh_pdbs/KH4_medium_uniform_350.pdb, KH4_medium_u350
color lightblue, KH4_medium_u350
align KH4_medium_u350, KH4_AF2
load kh_pdbs/KH4_loose_uniform_700.pdb, KH4_loose_u700
color green, KH4_loose_u700
align KH4_loose_u700, KH4_AF2
load kh_pdbs/KH4_loose_uniform_350.pdb, KH4_loose_u350
color palegreen, KH4_loose_u350
align KH4_loose_u350, KH4_AF2

# Group KH4
group KH4_group, KH4_AF2 KH4_baseline_current KH4_tight_u700 KH4_tight_u350 KH4_medium_u700 KH4_medium_u350 KH4_loose_u700 KH4_loose_u350

# KH5
load kh_pdbs/KH5_AF2_reference.pdb, KH5_AF2
color yellow, KH5_AF2
load kh_pdbs/KH5_baseline_current.pdb, KH5_baseline_current
color gray60, KH5_baseline_current
align KH5_baseline_current, KH5_AF2
load kh_pdbs/KH5_tight_uniform_700.pdb, KH5_tight_u700
color red, KH5_tight_u700
align KH5_tight_u700, KH5_AF2
load kh_pdbs/KH5_tight_uniform_350.pdb, KH5_tight_u350
color salmon, KH5_tight_u350
align KH5_tight_u350, KH5_AF2
load kh_pdbs/KH5_medium_uniform_700.pdb, KH5_medium_u700
color blue, KH5_medium_u700
align KH5_medium_u700, KH5_AF2
load kh_pdbs/KH5_medium_uniform_350.pdb, KH5_medium_u350
color lightblue, KH5_medium_u350
align KH5_medium_u350, KH5_AF2
load kh_pdbs/KH5_loose_uniform_700.pdb, KH5_loose_u700
color green, KH5_loose_u700
align KH5_loose_u700, KH5_AF2
load kh_pdbs/KH5_loose_uniform_350.pdb, KH5_loose_u350
color palegreen, KH5_loose_u350
align KH5_loose_u350, KH5_AF2

# Group KH5
group KH5_group, KH5_AF2 KH5_baseline_current KH5_tight_u700 KH5_tight_u350 KH5_medium_u700 KH5_medium_u350 KH5_loose_u700 KH5_loose_u350

# KH6
load kh_pdbs/KH6_AF2_reference.pdb, KH6_AF2
color yellow, KH6_AF2
load kh_pdbs/KH6_baseline_current.pdb, KH6_baseline_current
color gray60, KH6_baseline_current
align KH6_baseline_current, KH6_AF2
load kh_pdbs/KH6_tight_uniform_700.pdb, KH6_tight_u700
color red, KH6_tight_u700
align KH6_tight_u700, KH6_AF2
load kh_pdbs/KH6_tight_uniform_350.pdb, KH6_tight_u350
color salmon, KH6_tight_u350
align KH6_tight_u350, KH6_AF2
load kh_pdbs/KH6_medium_uniform_700.pdb, KH6_medium_u700
color blue, KH6_medium_u700
align KH6_medium_u700, KH6_AF2
load kh_pdbs/KH6_medium_uniform_350.pdb, KH6_medium_u350
color lightblue, KH6_medium_u350
align KH6_medium_u350, KH6_AF2
load kh_pdbs/KH6_loose_uniform_700.pdb, KH6_loose_u700
color green, KH6_loose_u700
align KH6_loose_u700, KH6_AF2
load kh_pdbs/KH6_loose_uniform_350.pdb, KH6_loose_u350
color palegreen, KH6_loose_u350
align KH6_loose_u350, KH6_AF2

# Group KH6
group KH6_group, KH6_AF2 KH6_baseline_current KH6_tight_u700 KH6_tight_u350 KH6_medium_u700 KH6_medium_u350 KH6_loose_u700 KH6_loose_u350

# KH7a
load kh_pdbs/KH7a_AF2_reference.pdb, KH7a_AF2
color yellow, KH7a_AF2
load kh_pdbs/KH7a_baseline_current.pdb, KH7a_baseline_current
color gray60, KH7a_baseline_current
align KH7a_baseline_current, KH7a_AF2
load kh_pdbs/KH7a_tight_uniform_700.pdb, KH7a_tight_u700
color red, KH7a_tight_u700
align KH7a_tight_u700, KH7a_AF2
load kh_pdbs/KH7a_tight_uniform_350.pdb, KH7a_tight_u350
color salmon, KH7a_tight_u350
align KH7a_tight_u350, KH7a_AF2
load kh_pdbs/KH7a_medium_uniform_700.pdb, KH7a_medium_u700
color blue, KH7a_medium_u700
align KH7a_medium_u700, KH7a_AF2
load kh_pdbs/KH7a_medium_uniform_350.pdb, KH7a_medium_u350
color lightblue, KH7a_medium_u350
align KH7a_medium_u350, KH7a_AF2
load kh_pdbs/KH7a_loose_uniform_700.pdb, KH7a_loose_u700
color green, KH7a_loose_u700
align KH7a_loose_u700, KH7a_AF2
load kh_pdbs/KH7a_loose_uniform_350.pdb, KH7a_loose_u350
color palegreen, KH7a_loose_u350
align KH7a_loose_u350, KH7a_AF2

# Group KH7a
group KH7a_group, KH7a_AF2 KH7a_baseline_current KH7a_tight_u700 KH7a_tight_u350 KH7a_medium_u700 KH7a_medium_u350 KH7a_loose_u700 KH7a_loose_u350

# KHb
load kh_pdbs/KHb_AF2_reference.pdb, KHb_AF2
color yellow, KHb_AF2
load kh_pdbs/KHb_baseline_current.pdb, KHb_baseline_current
color gray60, KHb_baseline_current
align KHb_baseline_current, KHb_AF2
load kh_pdbs/KHb_tight_uniform_700.pdb, KHb_tight_u700
color red, KHb_tight_u700
align KHb_tight_u700, KHb_AF2
load kh_pdbs/KHb_tight_uniform_350.pdb, KHb_tight_u350
color salmon, KHb_tight_u350
align KHb_tight_u350, KHb_AF2
load kh_pdbs/KHb_medium_uniform_700.pdb, KHb_medium_u700
color blue, KHb_medium_u700
align KHb_medium_u700, KHb_AF2
load kh_pdbs/KHb_medium_uniform_350.pdb, KHb_medium_u350
color lightblue, KHb_medium_u350
align KHb_medium_u350, KHb_AF2
load kh_pdbs/KHb_loose_uniform_700.pdb, KHb_loose_u700
color green, KHb_loose_u700
align KHb_loose_u700, KHb_AF2
load kh_pdbs/KHb_loose_uniform_350.pdb, KHb_loose_u350
color palegreen, KHb_loose_u350
align KHb_loose_u350, KHb_AF2

# Group KHb
group KHb_group, KHb_AF2 KHb_baseline_current KHb_tight_u700 KHb_tight_u350 KHb_medium_u700 KHb_medium_u350 KHb_loose_u700 KHb_loose_u350

# KH8
load kh_pdbs/KH8_AF2_reference.pdb, KH8_AF2
color yellow, KH8_AF2
load kh_pdbs/KH8_baseline_current.pdb, KH8_baseline_current
color gray60, KH8_baseline_current
align KH8_baseline_current, KH8_AF2
load kh_pdbs/KH8_tight_uniform_700.pdb, KH8_tight_u700
color red, KH8_tight_u700
align KH8_tight_u700, KH8_AF2
load kh_pdbs/KH8_tight_uniform_350.pdb, KH8_tight_u350
color salmon, KH8_tight_u350
align KH8_tight_u350, KH8_AF2
load kh_pdbs/KH8_medium_uniform_700.pdb, KH8_medium_u700
color blue, KH8_medium_u700
align KH8_medium_u700, KH8_AF2
load kh_pdbs/KH8_medium_uniform_350.pdb, KH8_medium_u350
color lightblue, KH8_medium_u350
align KH8_medium_u350, KH8_AF2
load kh_pdbs/KH8_loose_uniform_700.pdb, KH8_loose_u700
color green, KH8_loose_u700
align KH8_loose_u700, KH8_AF2
load kh_pdbs/KH8_loose_uniform_350.pdb, KH8_loose_u350
color palegreen, KH8_loose_u350
align KH8_loose_u350, KH8_AF2

# Group KH8
group KH8_group, KH8_AF2 KH8_baseline_current KH8_tight_u700 KH8_tight_u350 KH8_medium_u700 KH8_medium_u350 KH8_loose_u700 KH8_loose_u350

# ════════════════════════════════════════════════════════════
# Full-length AF2 structure with different colorings
# ════════════════════════════════════════════════════════════

set fetch_path, /tmp

# Load FL structure 4 times for different coloring schemes
load /home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb, FL_pLDDT
load /home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb, FL_tight
load /home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb, FL_medium
load /home/sbali/CALVADOS/examples/PARP14_MDP/input/parp14.pdb, FL_loose

# ── FL_pLDDT: color by B-factor (pLDDT) ──
show cartoon, FL_pLDDT
spectrum b, red_white_blue, FL_pLDDT, 30, 95
# AlphaFold pLDDT color scheme: blue=high, red=low

# ── FL_tight: gray base, restrained cores in red ──
show cartoon, FL_tight
color gray80, FL_tight
# Tight boundaries (pLDDT >= 80 cores)
color firebrick, FL_tight and resi 22-47       # RRM1
color firebrick, FL_tight and resi 150-224      # RRM2
color firebrick, FL_tight and resi 267-310      # RRM3
color red, FL_tight and resi 329-350            # KH1
color red, FL_tight and resi 432-450            # KH2
color red, FL_tight and resi 468-494            # KH3
color red, FL_tight and resi 521-570            # KH4
color red, FL_tight and resi 594-620            # KH5
# KH6 too small at tight
color red, FL_tight and resi 749-775            # KH7a
color firebrick, FL_tight and resi 801-946      # MD1L1
color firebrick, FL_tight and resi 1013-1098    # MD2
color firebrick, FL_tight and resi 1219-1279    # MD3
color red, FL_tight and resi 1428-1460          # KHb
color red, FL_tight and resi 1493-1519          # KH8
color firebrick, FL_tight and resi 1545-1571    # WWE
color firebrick, FL_tight and resi 1617-1699    # ART

# ── FL_medium: gray base, restrained cores in blue ──
show cartoon, FL_medium
color gray80, FL_medium
# Medium boundaries (pLDDT >= 70 cores)
color marine, FL_medium and resi 7-90           # RRM1
color marine, FL_medium and resi 149-224        # RRM2
color marine, FL_medium and resi 227-310        # RRM3
color blue, FL_medium and resi 320-374          # KH1
color blue, FL_medium and resi 385-454          # KH2
color blue, FL_medium and resi 467-520          # KH3
color blue, FL_medium and resi 521-593          # KH4
color blue, FL_medium and resi 594-626          # KH5
color blue, FL_medium and resi 714-737          # KH6
color blue, FL_medium and resi 738-775          # KH7a
color marine, FL_medium and resi 791-978        # MD1L1
color marine, FL_medium and resi 1005-1190      # MD2
color marine, FL_medium and resi 1216-1387      # MD3
color blue, FL_medium and resi 1427-1461        # KHb
color blue, FL_medium and resi 1462-1522        # KH8
color marine, FL_medium and resi 1535-1601      # WWE
color marine, FL_medium and resi 1604-1801      # ART

# ── FL_loose: gray base, restrained cores in green ──
show cartoon, FL_loose
color gray80, FL_loose
# Loose boundaries (pLDDT >= 60 cores)
color forest, FL_loose and resi 6-91            # RRM1
color forest, FL_loose and resi 148-224         # RRM2
color forest, FL_loose and resi 225-312         # RRM3
color green, FL_loose and resi 317-384          # KH1
color green, FL_loose and resi 385-454          # KH2
color green, FL_loose and resi 466-520          # KH3
color green, FL_loose and resi 521-593          # KH4
color green, FL_loose and resi 594-640          # KH5
color green, FL_loose and resi 666-707          # KH6
color green, FL_loose and resi 738-776          # KH7a
color forest, FL_loose and resi 791-978         # MD1L1
color forest, FL_loose and resi 1005-1191       # MD2
color forest, FL_loose and resi 1208-1388       # MD3
color green, FL_loose and resi 1425-1461        # KHb
color green, FL_loose and resi 1462-1533        # KH8
color forest, FL_loose and resi 1534-1601       # WWE
color forest, FL_loose and resi 1603-1801       # ART

# Group the FL views
group FL_views, FL_pLDDT FL_tight FL_medium FL_loose

# ════════════════════════════════════════════════════════════
# Display settings
# ════════════════════════════════════════════════════════════
show cartoon, FL_*
hide cartoon, KH*_group
set cartoon_trace_atoms, 1, KH*
set cartoon_tube_radius, 0.2, KH*
bg_color white
set ray_opaque_background, 1

# ════════════════════════════════════════════════════════════
# Scenes
# ════════════════════════════════════════════════════════════

# Scene 1: pLDDT colored full-length
hide everything, FL_tight FL_medium FL_loose
hide everything, KH*_group
show cartoon, FL_pLDDT
orient FL_pLDDT
scene pLDDT_view, store, Full-length colored by pLDDT (blue=high red=low)

# Scene 2: Tight boundaries
hide everything, FL_pLDDT FL_medium FL_loose
hide everything, KH*_group
show cartoon, FL_tight
orient FL_tight
scene Tight_boundaries, store, Tight restrained cores (pLDDT>=80) in red

# Scene 3: Medium boundaries
hide everything, FL_pLDDT FL_tight FL_loose
hide everything, KH*_group
show cartoon, FL_medium
orient FL_medium
scene Medium_boundaries, store, Medium restrained cores (pLDDT>=70) in blue

# Scene 4: Loose boundaries
hide everything, FL_pLDDT FL_tight FL_medium
hide everything, KH*_group
show cartoon, FL_loose
orient FL_loose
scene Loose_boundaries, store, Loose restrained cores (pLDDT>=60) in green

# Scene 5: All KH domain comparisons
hide everything, FL_*
show cartoon, KH*
orient KH*
scene KH_all_tests, store, All KH domains - test comparison (yellow=AF2 ref)

# Scene 6: KH7a + KHb together (split KH pair)
hide everything, all
show cartoon, KH7a_group KHb_group
orient KH7a_group KHb_group
scene KH7a_KHb, store, KH7a + KHb split pair comparison

# Scene 7: KH1-KH4 (best-restrained KH domains)
hide everything, all
show cartoon, KH1_group KH2_group KH3_group KH4_group
orient KH1_group KH2_group KH3_group KH4_group
scene KH1_to_KH4, store, KH1-4 comparison (best-behaved)

# Scene 8: KH5-KH6 (problematic KH domains)
hide everything, all
show cartoon, KH5_group KH6_group
orient KH5_group KH6_group
scene KH5_KH6, store, KH5 + KH6 comparison (hardest to restrain)

# Start with pLDDT view
scene pLDDT_view, recall

# Save session
save kh_restraint_comparison.pse