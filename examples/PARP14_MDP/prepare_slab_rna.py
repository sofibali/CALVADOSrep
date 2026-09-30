#!/usr/bin/env python
"""
Prepare the narrow RNA-effect slab panel, with every lesson from the first
campaign applied.

WHAT THIS PANEL ASKS
--------------------
Does RNA promote self-association, and does the effect scale with RNA-binding
valence? One axis (valence), everything else pinned, each construct run with and
without polyU40.

WHAT CHANGED FROM THE FIRST CAMPAIGN
------------------------------------
1. REPLICATES. The first campaign was n=1 -- all 20 runs shared
   random_number_seed 12345 -- so the +/-RNA and +/-ART differences rested on
   single trajectories whose width wanders by +/-30-45 nm. Now n seeds per run.

2. PAIRED SEEDS. Within a construct, the homotypic and +RNA arms of replicate i
   get the SAME seed, so the protein starting configuration is identical and RNA
   is the only difference. The seed varies across replicates, not across arms.

3. CONSTANT CONCENTRATION. The first panel spanned 26.7x in chain concentration
   because chain count was traded against construct size to equalise COST, not
   concentration -- so construct identity and concentration were confounded.
   Every run here is solved to the same residue molarity.

4. OPTIMISED RESTRAINTS. No run in the first campaign used them
   (custom_restraints: false, untrimmed domain cores). This applies the
   per-domain trims from the restraint sweep AND the 141-pair KH7a-KHb bridge,
   both mapped into each construct's local numbering.

5. EQUILIBRATION DENSITY IS RECORDED. md_full at 406 mM vitrified because
   slab_eq compressed it to ~1020 mg/mL, 3x a physiological condensate. Each run
   records the density its slab would have at the target thickness so an
   at-risk setup is visible before it burns 11 h.

6. CHECKPOINT CADENCE. checkpoint_interval = wfreq, so a preemption costs ~16 s
   rather than ~156 s and trim_dcd.py's boundary stays exact.

Written runs are classified by slab/analyze_slab.py as dispersing / arrested /
liquid BEFORE any c_sat is quoted -- a slab that never disperses can be a jammed
solid, which is how md_full nearly entered the record as a condensate.

    python prepare_slab_rna.py --report
    python prepare_slab_rna.py --seeds 3 --conc 100
"""
import argparse
import os
import shutil
import sys

import numpy as np
import yaml

CWD = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CWD)
from calvados.cfg import Config, Components          # noqa: E402

NA = 6.02214076e23
TEMP, IONIC, PH = 293.15, 0.15, 7.5
STEPS, STEPS_EQ, WFREQ = 200_000_000, 5_000_000, 100_000
RNA_NAME, RNA_NT, RNA_PER_PROTEIN = 'polyU40', 40, 0.30

# The CALVADOS RNA force field, as used in run 1 (prepare_slab.py:94). These are
# NOT the package defaults in calvados/data/default_component.yaml -- those are
# much stiffer (kb 8033/8033, ka 7.24) and use a 0.6 nm nonbonded cutoff instead
# of 2.0. Omitting them silently substitutes a different polymer, which would
# make the +/-RNA arms here incomparable with the first campaign.
RNA_OPTS = dict(rna_kb1=1400.0, rna_kb2=2200.0, rna_ka=4.20, rna_pa=3.14,
                rna_nb_sigma=0.4, rna_nb_scale=15, rna_nb_cutoff=2.0)
TARGET_BEADS = 50_000
MIN_CHAINS = 30

# Per-domain trims from the restraint sweep (prepare_fl_optimized.py).
# KH7a / KHb / KH8 stay at full extent because the KH7a-KHb bridge is defined
# on their full boundaries.
TRIMS = {'RRM1': 10, 'RRM2': 10, 'RRM3': 10, 'KH1': 5, 'KH2': 5, 'KH3': 5,
         'KH4': 10, 'KH5': 5, 'KH6': 10, 'KH7a': 0, 'MD1L1': 10, 'MD2': 10,
         'MD3': 10, 'KHb': 0, 'KH8': 0, 'WWE': 15, 'ART': 10}

# Structured extents these trims are applied to, in FL numbering.
EXTENTS = {'RRM1': (6, 88), 'RRM2': (150, 223), 'RRM3': (227, 301),
           'KH1': (315, 386), 'KH2': (387, 458), 'KH3': (459, 530),
           'KH4': (531, 596), 'KH5': (597, 666), 'KH6': (667, 737),
           'KH7a': (738, 789), 'MD1L1': (791, 978), 'MD2': (1003, 1190),
           'MD3': (1216, 1387), 'KHb': (1389, 1461), 'KH8': (1462, 1533),
           'WWE': (1523, 1601), 'ART': (1605, 1801)}

# The panel. `span` is the FL residue range; RNA-binding units are the RRMs and
# KH groups fully inside it.
PANEL = {
    'fl':            dict(span=(1, 1801),    note='valence 6, +ART — full length'),
    'kh1_art_full':  dict(span=(315, 1801),  note='valence 3, +ART'),
    'core_full_go':  dict(span=(738, 1801),  note='valence 2, +ART — KH7a→ART'),
    'md_full':       dict(span=(790, 1388),  note='valence 0, −ART — the one construct that condenses'),
    # −WWE−ART slices: cut at 1533, after KHb-KH8. Optional tier.
    'fl_khb':        dict(span=(1, 1533),    note='valence 6, −WWE −ART', slice_of='fl'),
    'kh1_khb':       dict(span=(315, 1533),  note='valence 3, −WWE −ART', slice_of='fl'),
    'core_khb':      dict(span=(738, 1533),  note='valence 2, −WWE −ART', slice_of='fl'),
}
CORE_TIER = ['fl', 'kh1_art_full', 'core_full_go', 'md_full']
NOWWEART_TIER = ['fl_khb', 'kh1_khb', 'core_khb']

RNA_UNITS = {'RRM1': (1, 145), 'RRM2': (146, 224), 'RRM3': (225, 314),
             'KH1-KH6': (315, 737), 'KH7a': (738, 789), 'KHb-KH8': (1389, 1533)}

RUN_PY = ("import os\nfrom calvados import sim\n\n"
          "path = os.path.dirname(os.path.abspath(__file__))\n"
          "sim.run(path=path, fconfig='config.yaml', fcomponents='components.yaml')\n")


def valence(span):
    return sum(1 for a, b in RNA_UNITS.values() if a >= span[0] and b <= span[1])


def local(span, fl):
    """FL residue -> 1-based index inside a contiguous slice, or None."""
    return fl - span[0] + 1 if span[0] <= fl <= span[1] else None


def trimmed_domains(span):
    """Trimmed restraint domains, in the construct's own numbering."""
    out = []
    for name, (s, e) in EXTENTS.items():
        t = TRIMS[name]
        s2, e2 = s + t, e - t
        if e2 <= s2:
            continue
        s2, e2 = max(s2, span[0]), min(e2, span[1])
        if e2 - s2 < 5:                      # too short to restrain meaningfully
            continue
        out.append([local(span, s2), local(span, e2)])
    return sorted(out)


def custom_restraints(span, sysname, nmol):
    """
    The KH7a-KHb bridge, remapped to local numbering and written out for EVERY
    chain copy. Only applies where both ends of the bridge are present.

    The per-copy expansion is not optional. CALVADOS's custom-restraint format
    is `name copy bead | name copy bead | r k` and `Sim.map_custom_restraints`
    resolves each line to one absolute bead pair via
    `start_bead + (copy-1)*nbeads + (bead-1)` -- it does NOT broadcast a line
    across copies. Emitting only `copy 1`, which is all the monomer runs ever
    needed, would bridge chain 1 of 47 and leave the other 46 unrestrained.

    Returns the file's lines, or None if this construct does not span both.
    """
    src = os.path.join(CWD, 'input', 'custom_restraints.txt')
    if not os.path.isfile(src):
        return None
    pairs = []
    for raw in open(src):
        raw = raw.strip()
        if not raw or raw.startswith('#'):
            continue
        a, b, c = [p.strip() for p in raw.split('|')]
        i, j = int(a.split()[2]), int(b.split()[2])
        li, lj = local(span, i), local(span, j)
        if li is None or lj is None:
            return None                      # bridge incomplete -> do not apply
        d, k = c.split()[0], c.split()[1]
        pairs.append((li, lj, d, k))
    if not pairs:
        return None
    return [f'{sysname} {c} {li} | {sysname} {c} {lj} | {d} {k}\n'
            for c in range(1, nmol + 1) for li, lj, d, k in pairs]


def fl_reference():
    """The full-length model, as raw PDB lines. Sub-constructs are exact slices
    of it -- verified: fl_wwe_full_go, core_full_go and kh1_art_full all sit at
    RMSD 0.000 nm against the corresponding slice."""
    p = os.path.join(CWD, 'slab', 'homotypic', 'fl', 'input', 'parp14.pdb')
    return [l for l in open(p) if l.startswith(('ATOM', 'HETATM'))]


def write_sliced_pdb(lines, span, out):
    """
    Residues span[0]..span[1], renumbered 1..N.

    Done at the text level rather than through MDAnalysis: slicing an
    AtomGroup and renumbering its residues in place does not survive a copy,
    and silently produced an empty selection.
    """
    n = 0
    prev = None
    with open(out, 'w') as fh:
        for l in lines:
            r = int(l[22:26])
            if not (span[0] <= r <= span[1]):
                continue
            if r != prev:
                n += 1
                prev = r
            fh.write(f'{l[:22]}{n:>4}{l[26:]}')
        fh.write('END\n')
    return n


def geometry(nres, conc_mM, xy_hint=None):
    """chains near TARGET_BEADS, then a box solved to the target concentration."""
    n = max(MIN_CHAINS, int(round(TARGET_BEADS / nres)))
    beads = n * nres
    dens = conc_mM * 1e-3 * NA * 1e-24        # beads per nm^3
    V = beads / dens
    xy = xy_hint or float(int(round((V / 13.0) ** (1 / 3.0))))
    lz = round(V / (xy * xy))
    return n, [xy, xy, float(lz)], beads


def dense_thickness(beads, xy, mg_ml=300.0):
    """Slab thickness if it condensed at mg_ml, nm."""
    mass_g = beads * 110 * 1.66053907e-24
    vol = mass_g / (mg_ml / 1000.0 / 1e21)
    return vol / (xy * xy)


def prepare(key, arm, seed, conc, outroot, fl_u):
    spec = PANEL[key]
    span = spec['span']
    nres = span[1] - span[0] + 1
    sysname = f'parp14_{key}'
    nchain, box, beads_p = geometry(nres, conc)
    n_rna = int(round(nchain * RNA_PER_PROTEIN)) if arm == 'rna' else 0

    d = os.path.join(outroot, f'{key}__{arm}__seed{seed}')
    inp = os.path.join(d, 'input')
    os.makedirs(inp, exist_ok=True)

    # structure
    pdb = os.path.join(inp, f'{sysname}.pdb')
    if spec.get('slice_of') or key != 'md_full':
        write_sliced_pdb(fl_u, span, pdb)
    else:
        src = None
        for c in ('parp14_md1md3_full.pdb', f'{sysname}.pdb'):
            p = os.path.join(CWD, 'slab', 'homotypic', 'md_full', 'input', c)
            if os.path.isfile(p):
                src = p
                break
        shutil.copy(src, pdb)

    # restraints
    with open(os.path.join(inp, 'domains.yaml'), 'w') as fh:
        yaml.safe_dump({sysname: trimmed_domains(span)}, fh)
    cres = custom_restraints(span, sysname, nchain)
    if cres:
        with open(os.path.join(inp, 'custom_restraints.txt'), 'w') as fh:
            fh.writelines(cres)

    # residue table (+ RNA rows when needed)
    base = os.path.join(CWD, 'slab', 'homotypic', 'fl', 'input',
                        'residues_CALVADOS3.csv')
    if arm == 'rna':
        rna_src = os.path.join(CWD, 'slab', 'rna', 'fl', 'input',
                               'residues_CALVADOS3_RNA.csv')
        fres = os.path.join(inp, 'residues_CALVADOS3_RNA.csv')
        shutil.copy(rna_src, fres)
        ffasta = os.path.join(inp, 'rna.fasta')
        with open(ffasta, 'w') as fh:
            fh.write(f'>{RNA_NAME}\n{"r" * RNA_NT}\n')
    else:
        fres = os.path.join(inp, 'residues_CALVADOS3.csv')
        shutil.copy(base, fres)
        ffasta = None

    cfg = Config(
        sysname=sysname, box=box, temp=TEMP, ionic=IONIC, pH=PH,
        topol='slab', slab_eq=True, k_eq=0.02, steps_eq=STEPS_EQ,
        slab_width=max(20, int(box[0])), wfreq=WFREQ, logfreq=WFREQ,
        steps=STEPS, runtime=0, platform='CUDA',
        restart='checkpoint', frestart='restart.chk', verbose=True,
        random_number_seed=seed,
    )
    cfg.write(d, name='config.yaml')

    # Config() has no keyword for these, so add them after the fact.
    extra = yaml.safe_load(open(os.path.join(d, 'config.yaml')))
    extra['checkpoint_interval'] = WFREQ
    if cres:
        extra['custom_restraints'] = True
        extra['custom_restraint_type'] = 'harmonic'
        extra['fcustom_restraints'] = 'input/custom_restraints.txt'
    with open(os.path.join(d, 'config.yaml'), 'w') as fh:
        yaml.safe_dump(extra, fh)

    ck = dict(fresidues=fres, pdb_folder=inp,
              fdomains=os.path.join(inp, 'domains.yaml'))
    if arm == 'rna':
        ck['ffasta'] = ffasta
        ck.update(RNA_OPTS)
    comps = Components(**ck)
    comps.add(name=sysname, molecule_type='protein', nmol=nchain,
              restraint=True, charge_termini='both')
    if arm == 'rna':
        comps.add(name=RNA_NAME, molecule_type='rna', nmol=n_rna)
    comps.write(d, name='components.yaml')

    with open(os.path.join(d, 'run.py'), 'w') as fh:
        fh.write(RUN_PY)

    val = valence(span)
    thick = dense_thickness(beads_p, box[0])
    meta = dict(construct=key, arm=arm, seed=seed, sysname=sysname, nres=int(nres),
                nchain=int(nchain), n_rna=int(n_rna), box=[float(x) for x in box], steps=STEPS,
                beads=int(beads_p + n_rna * RNA_NT * 2),
                conc_mM_residues=float(round(beads_p / NA / (box[0]*box[1]*box[2] * 1e-24) * 1e3, 1)),
                rna_valence=int(val),
                rna_per_chain=float(round(n_rna / nchain, 3)) if nchain else 0.0,
                rna_per_binding_unit=float(round(n_rna / nchain / val, 4)) if val else None,
                custom_restraints=bool(cres),
                n_custom_pairs=len(cres) if cres else 0,
                n_restraint_domains=len(trimmed_domains(span)),
                dense_slab_nm_at_300mgml=float(round(thick, 1)),
                dilute_nm=float(round(box[2] - thick, 1)),
                note=spec['note'])
    with open(os.path.join(d, 'slab_meta.yaml'), 'w') as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--conc', type=float, default=100.0, help='mM residues')
    ap.add_argument('--seeds', type=int, default=3)
    ap.add_argument('--tier', choices=('core', 'nowweart', 'both'), default='both')
    ap.add_argument('--out', default=os.path.join(CWD, 'slab', 'rnapanel'))
    ap.add_argument('--report', action='store_true')
    a = ap.parse_args()

    keys = {'core': CORE_TIER, 'nowweart': NOWWEART_TIER,
            'both': CORE_TIER + NOWWEART_TIER}[a.tier]
    seeds = [12345 + 11111 * i for i in range(a.seeds)]

    if a.report:
        print(f'target {a.conc:.0f} mM residues · {a.seeds} seed(s) · '
              f'{len(keys)} constructs × 2 arms = {len(keys) * 2 * a.seeds} runs\n')
        print(f"{'construct':<15}{'span':<12}{'val':>4}{'ch':>5}{'box (nm)':>17}"
              f"{'slab':>6}{'dil':>7}{'KHbr':>6}{'dom':>5}")
        print('-' * 80)
        for k in keys:
            sp = PANEL[k]['span']
            nres = sp[1] - sp[0] + 1
            n, box, beads = geometry(nres, a.conc)
            th = dense_thickness(beads, box[0])
            br = custom_restraints(sp, 'x', 1)
            print(f"{k:<15}{f'{sp[0]}-{sp[1]}':<12}{valence(sp):>4}{n:>5}"
                  f"{f'{box[0]:.0f}x{box[1]:.0f}x{box[2]:.0f}':>17}{th:>5.0f}n"
                  f"{box[2] - th:>6.0f}n{(len(br) if br else 0):>6}"
                  f"{len(trimmed_domains(sp)):>5}")
        hrs = len(keys) * 2 * a.seeds * 2.05e8 / 4970 / 3600
        print(f"\n~{hrs:.0f} GPU-hours = {hrs / 24:.1f} GPU-days "
              f"(~{hrs / 24 / 2:.1f} days on two free cards)")
        return

    fl_u = fl_reference()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for k in keys:
        for arm in ('homotypic', 'rna'):
            for sd in seeds:
                rows.append(prepare(k, arm, sd, a.conc, a.out, fl_u))
                print(f"  {rows[-1]['construct']:<15}{arm:<11}seed {sd}  "
                      f"{rows[-1]['nchain']:>3}ch  {rows[-1]['conc_mM_residues']:>6.1f} mM  "
                      f"{'KHbridge' if rows[-1]['custom_restraints'] else '—':>9}")
    print(f"\n{len(rows)} runs written to {a.out}")


if __name__ == '__main__':
    main()
