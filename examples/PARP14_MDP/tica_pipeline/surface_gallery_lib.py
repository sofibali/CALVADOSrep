#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared config/helpers for the PARP14 surface-property gallery
(make_surface_gallery.py). Reuses the existing domain-color and
active-site-remapping infrastructure instead of re-deriving it:

  - FL_DOMAIN_COLORS, CATALYTIC_FL, DOMAIN_UNITS, CONSTRUCT_UNITS,
    compute_fl_to_construct  <- tica_pipeline/make_cluster_pse.py
  - ACTIVE_SITES_FL, get_active_sites_for_set, get_units,
    UNIT_TO_SITE             <- sim_registry.py

Adds only what didn't already exist anywhere in the project: the
hydrophobicity scale, discrete residue-chemistry classes, and the
H-bond donor/acceptor atom table used for the active-site pocket
close-ups.
"""
import os
import sys
from pathlib import Path

TICA_DIR = Path(__file__).resolve().parent
ROOT_DIR = TICA_DIR.parent
sys.path.insert(0, str(TICA_DIR))
sys.path.insert(0, str(ROOT_DIR))

from make_cluster_pse import (               # noqa: E402
    FL_DOMAIN_COLORS, CATALYTIC_FL, DOMAIN_UNITS, CONSTRUCT_UNITS,
    compute_fl_to_construct,
)
import sim_registry as sr                     # noqa: E402

# ============================================================
# Chothia Gly-X-Gly max ASA (A^2), Wu et al. 2017 BioData Mining.
# Same table as sim_analysis/figure_sasa_faces.py, duplicated here
# (it's a short fixed reference table, not derived project logic).
# ============================================================
MAX_ASA_CHOTHIA = {
    'PHE': 210, 'ILE': 175, 'LEU': 170, 'VAL': 155, 'PRO': 145,
    'ALA': 115, 'GLY':  75, 'MET': 185, 'CYS': 135, 'TRP': 255,
    'TYR': 230, 'THR': 140, 'SER': 115, 'GLN': 180, 'ASN': 160,
    'GLU': 190, 'ASP': 150, 'HIS': 195, 'LYS': 200, 'ARG': 225,
}
RSA_BURIED_MAX = 0.05
RSA_PARTIAL_MAX = 0.20

# ============================================================
# Kyte-Doolittle hydrophobicity scale
# ============================================================
KYTE_DOOLITTLE = {
    'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
    'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
    'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
    'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5,
}

# ============================================================
# Discrete residue chemistry classification (for pocket close-ups)
# ============================================================
RESIDUE_CLASS = {}
for _r in ('ASP', 'GLU'):
    RESIDUE_CLASS[_r] = ('acidic', '#e6194b')
for _r in ('LYS', 'ARG', 'HIS'):
    RESIDUE_CLASS[_r] = ('basic', '#4363d8')
for _r in ('PHE', 'TRP', 'TYR'):
    RESIDUE_CLASS[_r] = ('aromatic', '#f58231')
for _r in ('SER', 'THR', 'ASN', 'GLN', 'CYS'):
    RESIDUE_CLASS[_r] = ('polar', '#3cb44b')
for _r in ('ALA', 'VAL', 'LEU', 'ILE', 'MET', 'GLY', 'PRO'):
    RESIDUE_CLASS[_r] = ('nonpolar', '#dcdcdc')

RESIDUE_CLASS_LEGEND = {
    'acidic':   '#e6194b',
    'basic':    '#4363d8',
    'aromatic': '#f58231',
    'polar':    '#3cb44b',
    'nonpolar': '#dcdcdc',
}

# ============================================================
# H-bond donor / acceptor atoms (heavy-atom based; these all-atom
# models carry no hydrogens, so classification is by heavy atom
# identity, matching standard PyMOL/Chimera H-bond heuristics).
# role in {'donor', 'acceptor', 'both'}
# Backbone N (donor, except Pro) / O (acceptor) added generically below.
# ============================================================
SIDECHAIN_HBOND_ATOMS = {
    'ASP': [('OD1', 'acceptor'), ('OD2', 'acceptor')],
    'GLU': [('OE1', 'acceptor'), ('OE2', 'acceptor')],
    'ASN': [('OD1', 'acceptor'), ('ND2', 'donor')],
    'GLN': [('OE1', 'acceptor'), ('NE2', 'donor')],
    'LYS': [('NZ', 'donor')],
    'ARG': [('NE', 'donor'), ('NH1', 'donor'), ('NH2', 'donor')],
    'HIS': [('ND1', 'both'), ('NE2', 'both')],
    'SER': [('OG', 'both')],
    'THR': [('OG1', 'both')],
    'TYR': [('OH', 'both')],
    'TRP': [('NE1', 'donor')],
    'CYS': [('SG', 'both')],
    'MET': [('SD', 'acceptor')],
}

DONOR_COLOR = '#1f78ff'
ACCEPTOR_COLOR = '#ff3b1f'
BOTH_COLOR = '#a020f0'
HBOND_ROLE_COLOR = {'donor': DONOR_COLOR, 'acceptor': ACCEPTOR_COLOR, 'both': BOTH_COLOR}


def get_pocket_for_set(set_key, site_name):
    """(catalytic_resids, pocket_resids) in the SET's construct numbering,
    via sim_registry's FL->construct remapping. site_name in sr.SITE_NAMES."""
    sites = sr.get_active_sites_for_set(set_key)
    if site_name not in sites:
        return [], []
    return sites[site_name]['catalytic'], sites[site_name]['pocket']


def get_domain_ranges_for_set(set_key):
    """{domain_name: (color_name, rgb, (start,end))} in construct numbering,
    for every FL_DOMAIN_COLORS domain actually present in this set."""
    fl_to_c = compute_fl_to_construct(set_key)
    out = {}
    for dname, (color_name, (fs, fe), rgb) in FL_DOMAIN_COLORS.items():
        in_construct = [fl_to_c[r] for r in range(fs, fe + 1) if r in fl_to_c]
        if in_construct:
            out[dname] = (color_name, rgb, (min(in_construct), max(in_construct)))
    return out


def classify_rsa(rsa):
    if rsa is None:
        return 'unknown'
    if rsa < RSA_BURIED_MAX:
        return 'buried'
    elif rsa < RSA_PARTIAL_MAX:
        return 'partial'
    return 'surface'
