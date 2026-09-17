#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Puncta / self-association feature table for the PARP14 combinatorial library.

CONTEXT
-------
There is a pooled plasmid library of ~2042 PARP14 domain-deletion constructs.
Cells are sorted into DIFFUSE vs PUNCTATE pools and each pool is sequenced, so
every construct gets a read frequency in each pool. The punctate-vs-diffuse
log-ratio is a quantitative, per-construct self-association score. Single
constructs (FL, plus ART / MD1 / MD2 catalytic-dead mutants) serve as controls.

This script produces the MODEL-SIDE half of that comparison: an objective
feature table for all 2^11-1 = 2047 domain combinations, keyed so it can be
joined to the sequencing counts as soon as they exist.

WHY THESE FEATURES
------------------
Puncta here are expected to be a multivalency/percolation phenomenon: PARP14
cross-links RNA through its RNA-binding modules and cross-links ADP-ribosylated
substrate through its ADPr-reader modules, while its own catalytic domains
change how many ADPr marks exist to be read.

    RNA-binding valence      RRM1-3 + the KH domains
    ADPr-reader valence      MD2, MD3, WWE   (mark readers)
    writer                   ART             (creates marks -> creates crosslinks)
    eraser                   MD1             (removes marks -> removes crosslinks)

That decomposition is exactly what the catalytic-dead controls interrogate:
ART-dead removes mark creation, MD1-dead removes mark removal, MD2-dead removes
mark reading. Each predicts a different direction of change in puncta score.

Sequence-level LLPS descriptors (FCR, NCPR, kappa, SCD, aromatic/Arg content)
are included because they capture the classic sticker-and-spacer determinants
that domain counts alone miss.

IMPORTANT -- WEIGHTS ARE NOT INVENTED HERE
------------------------------------------
This script deliberately does NOT hand-tune a propensity score. It emits
features. `--prior` optionally adds ONE clearly-labelled, unfitted hypothesis
score (the multivalency product) for eyeballing only. The real weights should
be FIT against the sort-seq enrichment via `--fit`, which is how you find out
whether the multivalency picture is right rather than assuming it.

Usage:
    python predict_puncta_propensity.py                     # features -> CSV
    python predict_puncta_propensity.py --prior             # + hypothesis score
    python predict_puncta_propensity.py --fit counts.csv    # join + fit + validate
"""

import os
import sys
import csv
import itertools
import numpy as np
from argparse import ArgumentParser

CWD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(CWD, 'data')
sys.path.insert(0, CWD)
from _fig_layout import get_fig_dir as _get_fig_dir

PARP14_ROOT = '/home/sbali/CALVADOS/parp14'
FASTA = os.path.join(PARP14_ROOT, 'input', 'PARP14.fasta')
BOUNDS = os.path.join(PARP14_ROOT, 'input', 'domain_boundaries.csv')

# The 11 grouped units, in N->C order. This is the vocabulary of the AF3 library
# and of the plasmid library -- NOT the finer 17-domain naming used by
# parp14/domain_combinations/, which is a different, older enumeration.
UNITS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'MD1L1',
         'MD2', 'MD3', 'KHb-KH8', 'WWE', 'ART']

# How many INDIVIDUAL folded domains each grouped unit actually contains.
# Valence is about the number of independent binding modules, so the grouped
# KH units must not be counted as one.
N_DOMAINS_IN_UNIT = {
    'RRM1': 1, 'RRM2': 1, 'RRM3': 1,
    'KH1-KH6': 6,     # KH1..KH6
    'KH7a': 1,
    'MD1L1': 1, 'MD2': 1, 'MD3': 1,
    'KHb-KH8': 2,     # KHb + KH8
    'WWE': 1, 'ART': 1,
}

RNA_BINDING_UNITS = ['RRM1', 'RRM2', 'RRM3', 'KH1-KH6', 'KH7a', 'KHb-KH8']
# MD1 is the hydrolase (eraser) and a weak ADPr binder; the READERS that would
# cross-link ADP-ribosylated substrate are MD2, MD3 and WWE.
ADPR_READER_UNITS = ['MD2', 'MD3', 'WWE']
WRITER_UNIT = 'ART'
ERASER_UNIT = 'MD1L1'

AROMATIC = set('FWY')
POSITIVE = set('KR')
NEGATIVE = set('DE')


# ============================================================
# Sequence construction
# ============================================================

def load_reference():
    seq = ''.join(l.strip() for l in open(FASTA).readlines()[1:])
    bounds = {}
    with open(BOUNDS) as fh:
        for row in csv.DictReader(fh):
            bounds[row['Domain'].strip()] = (int(row['Start']), int(row['End']))
    missing = [u for u in UNITS if u not in bounds]
    if missing:
        sys.exit(f"ERROR: {BOUNDS} is missing unit(s): {missing}")
    return seq, bounds


def combo_sequence(combo, seq, bounds):
    """Concatenate the units' residues in N->C order.

    Uses a residue SET rather than string concatenation because the boundary
    table is not a clean partition: MD1L1 ends at 1004 and MD2 starts at 1004,
    so naive concatenation would duplicate residue 1004 whenever both units are
    present. (The 1194-1206 gap between MD2 and MD3 is intentional and stays
    excluded.)
    """
    resids = set()
    for u in combo:
        s, e = bounds[u]
        resids.update(range(s, e + 1))
    return ''.join(seq[i - 1] for i in sorted(resids)), len(resids)


# ============================================================
# Sequence descriptors
# ============================================================

def charge_vector(s):
    return np.array([1.0 if c in POSITIVE else -1.0 if c in NEGATIVE else 0.0
                     for c in s])


def compute_scd(s):
    """Sequence Charge Decoration (Sawle & Ghosh 2015, J Chem Phys 143:085101).

        SCD = (1/N) * sum_{i>j} q_i q_j * sqrt(i-j)

    More negative = charges better segregated into blocks = more
    self-association-prone. Vectorized over charged positions only, since
    uncharged residues contribute nothing and PARP14 constructs run to 1801 aa.
    """
    q = charge_vector(s)
    idx = np.nonzero(q)[0]
    if len(idx) < 2:
        return 0.0
    qc = q[idx]
    # pairwise over charged positions: i > j
    di = idx[:, None] - idx[None, :]
    qq = qc[:, None] * qc[None, :]
    lower = np.tril(np.ones_like(di, dtype=bool), k=-1)
    return float(np.sum(qq[lower] * np.sqrt(di[lower])) / len(s))


def compute_kappa(s):
    """localCIDER kappa: charge patterning, 0 (mixed) -> 1 (segregated).

    Returns NaN when undefined (no charges, or all-charged), which localCIDER
    signals with -1.
    """
    try:
        from localcider.sequenceParameters import SequenceParameters
        k = SequenceParameters(s).get_kappa()
        return float(k) if k is not None and k >= 0 else np.nan
    except Exception:
        return np.nan


def sequence_features(s, want_kappa=True):
    n = len(s)
    q = charge_vector(s)
    npos = int((q > 0).sum())
    nneg = int((q < 0).sum())
    return {
        'length': n,
        'fcr': (npos + nneg) / n,
        'ncpr': (npos - nneg) / n,
        'frac_pos': npos / n,
        'frac_neg': nneg / n,
        'frac_aromatic': sum(c in AROMATIC for c in s) / n,
        'frac_Y': s.count('Y') / n,
        'frac_W': s.count('W') / n,
        'frac_F': s.count('F') / n,
        'frac_R': s.count('R') / n,
        'frac_K': s.count('K') / n,
        'frac_G': s.count('G') / n,
        'frac_S': s.count('S') / n,
        'frac_P': s.count('P') / n,
        'kappa': compute_kappa(s) if want_kappa else np.nan,
        'scd': compute_scd(s),
    }


def architecture_features(combo):
    cs = set(combo)
    rna_val = sum(N_DOMAINS_IN_UNIT[u] for u in RNA_BINDING_UNITS if u in cs)
    adpr_val = sum(N_DOMAINS_IN_UNIT[u] for u in ADPR_READER_UNITS if u in cs)
    feats = {
        'n_units': len(combo),
        'n_domains': sum(N_DOMAINS_IN_UNIT[u] for u in combo),
        'rna_valence': rna_val,
        'adpr_reader_valence': adpr_val,
        'has_writer_ART': int(WRITER_UNIT in cs),
        'has_eraser_MD1': int(ERASER_UNIT in cs),
        # net mark-generating capacity: writes marks but also erases them
        'net_mark_capacity': int(WRITER_UNIT in cs) - int(ERASER_UNIT in cs),
    }
    for u in UNITS:
        feats[f'has_{u.replace("-", "_")}'] = int(u in cs)
    return feats


def prior_score(f):
    """UNFITTED hypothesis score -- for eyeballing only, never for conclusions.

    Multivalent cross-linking needs BOTH arms, so the product (not the sum) of
    the two valences is the natural form: a construct with ten RNA-binding
    domains and zero ADPr readers cross-links nothing. log1p keeps the huge
    KH1-KH6 valence from swamping everything. The writer term is a modest
    bonus for being able to create marks at all.
    """
    return (np.log1p(f['rna_valence']) * np.log1p(f['adpr_reader_valence'])
            * (1.0 + 0.5 * f['has_writer_ART']))


# ============================================================
# Build
# ============================================================

_WORKER_REF = {}


def _init_worker():
    _WORKER_REF['ref'] = load_reference()


def _build_row(args):
    combo, with_prior, want_kappa = args
    seq, bounds = _WORKER_REF['ref']
    s, nres = combo_sequence(combo, seq, bounds)
    name = '_'.join(combo)
    row = {
        'combo': name,
        'combo_key': name.lower(),
        # 11-bit presence barcode in canonical N->C unit order: the most
        # robust join key, immune to naming/ordering drift between the
        # AF3 pipeline, the plasmid library and this script.
        'barcode': ''.join('1' if u in combo else '0' for u in UNITS),
        'n_residues': nres,
    }
    row.update(architecture_features(combo))
    row.update(sequence_features(s, want_kappa=want_kappa))
    if with_prior:
        row['prior_multivalency_score'] = prior_score(row)
    return row


def build_table(with_prior=False, want_kappa=True, workers=None):
    """Build the full 2047-row feature table.

    Parallelized because localCIDER's kappa costs ~2.5 s on a 1801-residue
    sequence (everything else is microseconds), so a serial build of all 2047
    combinations runs ~85 min. Pass --no-kappa for a near-instant table when
    charge patterning isn't needed.
    """
    combos = [c for r in range(1, len(UNITS) + 1)
              for c in itertools.combinations(UNITS, r)]
    tasks = [(c, with_prior, want_kappa) for c in combos]

    if workers is None:
        workers = max(1, (os.cpu_count() or 4) - 2)
    workers = max(1, min(workers, len(tasks)))

    if workers == 1 or not want_kappa:
        # Without kappa the whole build is ~1 s; process startup would dominate.
        _init_worker()
        return [_build_row(t) for t in tasks]

    from concurrent.futures import ProcessPoolExecutor
    rows = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
        for i, row in enumerate(ex.map(_build_row, tasks, chunksize=8), 1):
            rows.append(row)
            if i % 250 == 0:
                print(f"    {i}/{len(tasks)} combinations", flush=True)
    # ProcessPoolExecutor.map preserves input order, but sort defensively so the
    # output CSV is stable regardless of executor behaviour.
    order = {'_'.join(c): i for i, c in enumerate(combos)}
    rows.sort(key=lambda r: order[r['combo']])
    return rows


def write_table(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = list(rows[0].keys())
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return path


# ============================================================
# Experimental join + fit (used once sort-seq counts exist)
# ============================================================

def load_counts(path):
    """Read the sort-seq table.

    Expected columns (case-insensitive), one row per construct:
        combo | combo_key | barcode      -- any one is enough to join on
        punctate_reads, diffuse_reads    -- raw counts
      or
        enrichment                       -- a precomputed log2 ratio

    Returns {join_key_kind: {key: enrichment}}.
    """
    with open(path) as fh:
        rdr = csv.DictReader(fh)
        fields = {f.lower().strip(): f for f in (rdr.fieldnames or [])}
        keycol = next((fields[k] for k in ('barcode', 'combo_key', 'combo')
                       if k in fields), None)
        if keycol is None:
            sys.exit(f"ERROR: {path} needs one of: barcode / combo_key / combo")
        has_counts = 'punctate_reads' in fields and 'diffuse_reads' in fields
        has_enr = 'enrichment' in fields
        if not (has_counts or has_enr):
            sys.exit(f"ERROR: {path} needs punctate_reads+diffuse_reads, "
                     f"or a precomputed enrichment column")

        out = {}
        tot_p = tot_d = 0
        raw = []
        for row in rdr:
            key = (row[keycol] or '').strip()
            if not key:
                continue
            if has_counts:
                p = float(row[fields['punctate_reads']] or 0)
                d = float(row[fields['diffuse_reads']] or 0)
                tot_p += p
                tot_d += d
                raw.append((key, p, d))
            else:
                out[key.lower()] = float(row[fields['enrichment']])
        if has_counts:
            # log2 enrichment of pool frequencies, with a +1 pseudocount so
            # constructs absent from one pool stay finite. Normalizing by pool
            # totals matters: the two sorted pools are not equally deep.
            for key, p, d in raw:
                fp = (p + 1) / (tot_p + len(raw))
                fd = (d + 1) / (tot_d + len(raw))
                out[key.lower()] = float(np.log2(fp / fd))
        return keycol.lower(), out


def fit_model(rows, counts_path):
    kind, obs = load_counts(counts_path)
    keyfield = {'barcode': 'barcode', 'combo_key': 'combo_key',
                'combo': 'combo_key'}[kind]

    X_cols = ['rna_valence', 'adpr_reader_valence', 'has_writer_ART',
              'has_eraser_MD1', 'n_domains', 'length', 'fcr', 'ncpr',
              'frac_aromatic', 'frac_R', 'scd']
    X, y, names = [], [], []
    for r in rows:
        key = str(r[keyfield]).lower()
        if key in obs:
            vals = [r[c] for c in X_cols]
            if any(isinstance(v, float) and np.isnan(v) for v in vals):
                continue
            X.append(vals)
            y.append(obs[key])
            names.append(r['combo'])
    if len(X) < 20:
        sys.exit(f"ERROR: only {len(X)} constructs joined between the feature "
                 f"table and {counts_path} -- check the key column.")

    X = np.asarray(X, float)
    y = np.asarray(y, float)
    print(f"  Joined {len(y)} constructs on '{kind}'")

    # z-score so coefficients are comparable across features on wildly
    # different scales (length ~1800 vs fcr ~0.3)
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Xz = (X - mu) / sd
    Xd = np.column_stack([np.ones(len(Xz)), Xz])

    # Held-out evaluation: an in-sample R^2 on 11 features would flatter any
    # model. 5-fold CV on shuffled rows is the honest number.
    rng = np.random.default_rng(0)
    order = rng.permutation(len(y))
    folds = np.array_split(order, 5)
    preds = np.zeros_like(y)
    for i in range(5):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(5) if j != i])
        beta, *_ = np.linalg.lstsq(Xd[tr], y[tr], rcond=None)
        preds[te] = Xd[te] @ beta
    ss_res = float(np.sum((y - preds) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2_cv = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rho = float(np.corrcoef(y, preds)[0, 1])

    beta_full, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    print(f"\n  5-fold CV R^2 = {r2_cv:.3f}   Pearson r = {rho:.3f}")
    print(f"\n  {'feature':<24}{'std. coefficient':>18}")
    print('  ' + '-' * 42)
    for c, b in sorted(zip(X_cols, beta_full[1:]), key=lambda t: -abs(t[1])):
        print(f"  {c:<24}{b:>+18.4f}")
    print(f"\n  Coefficients are in units of log2 enrichment per 1 SD of feature.")
    print(f"  A positive coefficient = that feature pushes constructs toward the")
    print(f"  PUNCTATE pool.")
    return {'r2_cv': r2_cv, 'pearson': rho,
            'coefficients': dict(zip(X_cols, beta_full[1:].tolist())),
            'n': len(y)}


def main():
    ap = ArgumentParser(description='Puncta/self-association feature table for '
                                    'the PARP14 combinatorial library.')
    ap.add_argument('--prior', action='store_true',
                    help='Add the unfitted multivalency hypothesis score '
                         '(for eyeballing only -- not a prediction)')
    ap.add_argument('--fit', metavar='COUNTS_CSV', default=None,
                    help='Sort-seq table to join and fit against. Needs a '
                         'barcode/combo_key/combo column plus either '
                         'punctate_reads+diffuse_reads or enrichment.')
    ap.add_argument('--no-kappa', action='store_true',
                    help='Skip localCIDER kappa (~2.5 s/sequence, the only slow '
                         'feature). Everything else builds in about a second.')
    ap.add_argument('--workers', type=int, default=None,
                    help='Parallel workers for the feature build '
                         '(default: CPU count - 2)')
    ap.add_argument('--out', default=os.path.join(DATA_PATH, 'puncta_features.csv'))
    args = ap.parse_args()

    rows = build_table(with_prior=args.prior, want_kappa=not args.no_kappa,
                       workers=args.workers)
    path = write_table(rows, args.out)
    print(f"  Built {len(rows)} domain combinations "
          f"({len(UNITS)} units -> 2^{len(UNITS)}-1)")
    print(f"  Saved: {os.path.relpath(path, CWD)}")

    if args.no_kappa:
        print(f"  note: kappa skipped (--no-kappa); the column is present but all NaN")
    else:
        nan_k = sum(1 for r in rows
                    if isinstance(r['kappa'], float) and np.isnan(r['kappa']))
        if nan_k:
            print(f"  note: kappa undefined for {nan_k} construct(s) "
                  f"(localCIDER returns -1 for no-charge / all-charge sequences)")

    if args.prior:
        top = sorted(rows, key=lambda r: -r['prior_multivalency_score'])[:10]
        print(f"\n  Top 10 by UNFITTED multivalency prior "
              f"(hypothesis only -- not a prediction):")
        for r in top:
            print(f"    {r['prior_multivalency_score']:6.2f}  "
                  f"rna={r['rna_valence']:2d} adpr={r['adpr_reader_valence']:2d} "
                  f"ART={r['has_writer_ART']}  {r['combo']}")

    if args.fit:
        print()
        fit_model(rows, args.fit)
    else:
        print(f"\n  No sort-seq counts supplied. Once the diffuse/punctate")
        print(f"  sequencing table exists, fit real weights with:")
        print(f"      python {os.path.basename(__file__)} --fit <counts.csv>")


if __name__ == '__main__':
    main()
