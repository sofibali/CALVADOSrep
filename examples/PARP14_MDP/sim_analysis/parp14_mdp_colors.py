"""
PARP14_MDP color palette -- adapts the project-wide style guide
(sofibali/parp14-structure-analysis, docs/parp14_colors.py + STYLE_GUIDE.md,
commit 671cd70) to THIS repo's own domain-unit/construct naming.

Two things are reused directly, one is deliberately NOT:

1. DOMAIN identity colors for RRM1/RRM2/RRM3/MD1(->md1l1)/MD2/MD3/KH_b(->
   khb-kh8)/WWE/ART are reused verbatim -- these domains are named and
   bounded the same way in both projects.
2. Sequential/diverging colormaps and status colors are reused verbatim.
3. KH_N / KH_Ca are NOT reused for this project's kh1-kh6 / kh7a, because
   the two projects split the KH domains at DIFFERENT residues:

     Style guide (experimental):  KH_N = KH1-3 (FL 315-520)
                                   KH_Ca = KH4-7 (FL 521-789)
     This project (simulation):   kh1-kh6 = KH1-6 (FL 315-737)
                                   kh7a = KH7 (FL 738-789)

   (KH1-7 sub-boundaries confirmed from this repo's own
   prepare_all_fragments.py DOMAIN_EXTENTS: KH1 315-384, KH2 385-454,
   KH3 455-520, KH4 521-593, KH5 594-665, KH6 666-737, KH7a 738-789.)

   kh1-kh6 spans ALL of KH_N plus part of KH_Ca (KH4-6); kh7a is the tail
   end of KH_Ca only. Borrowing KH_N's/KH_Ca's hex values for kh1-kh6/kh7a
   would silently claim a residue-boundary equivalence that isn't true, so
   this module generates its OWN 2-step teal-green ramp instead -- same
   accessory-KH hue family (same "kind of thing" signal the style guide
   uses elsewhere), independently stepped, with no claim about where the
   experimental scheme's KH_N/KH_Ca boundary falls.
"""
import colorsys

# ---------------------------------------------------------------------------
# Reused verbatim from parp14-structure-analysis/docs/parp14_colors.py --
# same domain, same boundaries in both projects.
# ---------------------------------------------------------------------------
DOMAIN = {
    "rrm1": "#93700F",
    "rrm2": "#B88C13",
    "rrm3": "#DDA816",
    "md1l1": "#0072B2",   # style guide's "MD1" (eraser, blue)
    "md2": "#6A3FA0",
    "md3": "#9169C7",
    "khb-kh8": "#1CC497",  # style guide's "KH_b" -- boundaries match exactly
    "wwe": "#9C4090",
    "art": "#C8481A",
}


def _extend_ramp(base_hex, n_steps, l_min=0.28, l_max=0.72):
    """Same generator as parp14_colors.py's extend_ramp(): n_steps same-hue
    swatches from base_hex, light -> dark, within the validated ordinal-ramp
    lightness band."""
    r = int(base_hex[1:3], 16) / 255
    g = int(base_hex[3:5], 16) / 255
    b = int(base_hex[5:7], 16) / 255
    h, _, s = colorsys.rgb_to_hls(r, g, b)
    steps = []
    for i in range(n_steps):
        L = l_max - (l_max - l_min) * i / max(n_steps - 1, 1)
        r2, g2, b2 = colorsys.hls_to_rgb(h, L, s)
        steps.append("#{:02x}{:02x}{:02x}".format(int(r2 * 255), int(g2 * 255), int(b2 * 255)))
    return steps

# Independently-generated 2-step teal-green ramp (base hue = style guide's
# KH_N, #117C60, reused only for the HUE -- "this is a KH accessory domain"
# -- not the boundary it was validated for). Style guide convention
# (confirmed from its own hexes: RRM1 L=0.318 darkest -> RRM3 L=0.476
# lightest, same for KH_N->KH_b) is N-terminal-most = darkest, C-terminal-
# most = lightest; _extend_ramp emits index 0 = l_max (lightest) first, so
# the N-terminal-most domain takes the LAST step, not the first.
# kh1-kh6 (N-terminal-most of this project's split) -> darker step;
# kh7a (C-terminal-most of this project's split) -> lighter step.
#
# l_max is capped at 0.45 (not the function's 0.72 default) specifically for
# kh7a: it's the one ramp member that sits directly adjacent to md1l1 in the
# real sequence (a "family handoff" boundary, validated as an adjacent
# categorical pair per the style guide's own method) rather than only ever
# appearing next to other ramp members -- HLS lightness and the validator's
# OKLCH lightness don't map linearly, so 0.72 HLS produced an OKLCH L of
# 0.873 (outside the 0.43-0.77 categorical band); 0.45 was the first value
# that cleared validate_palette.js on the kh7a<->md1l1 pair (checked
# directly, see git history of this file).
_kh_ramp = _extend_ramp("#117C60", 2, l_max=0.45)
DOMAIN["kh1-kh6"] = _kh_ramp[1]  # darker (N-terminal-most)
DOMAIN["kh7a"] = _kh_ramp[0]     # lighter (C-terminal-most)

DOMAIN_ORDER = ["rrm1", "rrm2", "rrm3", "kh1-kh6", "kh7a", "md1l1", "md2", "md3",
                "khb-kh8", "wwe", "art"]

# ---------------------------------------------------------------------------
# Reused verbatim -- magnitude / polarity / status, not identity-specific.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# CONSTRUCT / simulation-set identity -- a second, independent categorical
# channel from DOMAIN above (different question: "which named CALVADOS run is
# this," not "which domain"). The style guide's own CONSTRUCT dict doesn't
# cover this project's set names (EV/MD1/MD2-ART/etc. are a different,
# simpler experimental library), so this is built fresh with the same method
# (Okabe-Ito seed, validated all-pairs, snap-to-passing on failures) rather
# than adapted.
#
# 13 named sets split into two groups that are NEVER drawn in the same figure
# (the Go-restraint and non-Go-restraint comparison plots are already
# separate, per how the analysis scripts are run) -- so each group only needs
# to be internally all-pairs distinct, and hues may repeat across groups.
#
# Group 1 -- Go-restraint sets (6): standard Okabe-Ito qualitative palette
# with yellow (#F0E442) dropped. Confirmed by direct validation (not just
# citing the style guide's note) that Okabe-Ito's yellow fails outright here:
# OKLCH L=0.902 (band is 0.43-0.77) and contrast 1.29:1 vs the light surface.
# Tried repeatedly to rescue it (re-stepped darker along its own hue: still
# collides with the orange already in the set once it's dark enough to clear
# the lightness band; shifted hue toward green: collides with teal instead) --
# every fix that keeps Okabe-Ito's other 6 anchors fixed and adds a 7th/lets
# yellow back in fails all-pairs somewhere. Dropping yellow leaves exactly 6
# hues, which is exactly this group's size, and passes clean (ALL CHECKS
# PASS, light mode, --pairs all).
CONSTRUCT_GO = {
    "core_full_go": "#E69F00",
    "kh1_art_full": "#56B4E9",
    "kh1_wwe_full": "#009E73",
    "core_wwe_full_go": "#0072B2",
    "fl_wwe_full_go": "#D55E00",
    # sim_registry.py's canonical key for this run is 'fl_go' (the original
    # 25-rep FL Go-restraint set); 'fl_go_5rep1us' is the local alias
    # compare_full_sims_grid.py registers via --sim-folder-as for the
    # length/replicate-matched 5x1000ns subset. Same color either name.
    "fl_go": "#CC79A7",
    "fl_go_5rep1us": "#CC79A7",
}
CONSTRUCT_GO_ORDER = ["core_full_go", "kh1_art_full", "kh1_wwe_full",
                      "core_wwe_full_go", "fl_wwe_full_go", "fl_go"]

# Group 2 -- non-Go-restraint sets (7). Okabe-Ito minus yellow only has 6
# hues, one short of this group's size, and every 7th-hue candidate tried
# against those same 6 anchors (olive/mustard region between orange and
# vermillion, green region between orange and teal) failed CVD or
# normal-vision separation against a neighbor -- Okabe-Ito's 6 remaining
# anchors already tile the hue wheel too densely to fit a 7th safely. Rather
# than force a collision, this group uses an independently-derived 7-hue set
# (greedy max-min-Delta-E search over an OKLCH-band-and-chroma-filtered pool,
# then verified with the same validator): ALL CHECKS PASS, light mode,
# --pairs all (worst-case CVD dE 10.7, worst-case normal-vision dE 17.6, both
# above target).
CONSTRUCT_NOGO = {
    "md_full": "#c8831c",
    "mka_full": "#2845bd",
    "md2_art_full": "#eb79eb",
    "md2_wwe_full": "#a82366",
    "md3_art_full": "#1cc895",
    "md3_wwe_full": "#2387a8",
    "mka_wwe_full": "#8a62e8",
}
CONSTRUCT_NOGO_ORDER = ["md_full", "mka_full", "md2_art_full", "md2_wwe_full",
                        "md3_art_full", "md3_wwe_full", "mka_wwe_full"]

# Isolation-test pair (md1_md2, md2_md3) -- always its own small figure, only
# needs to be distinct from itself. Reuses two already-validated CONSTRUCT_NOGO
# hues (dE 33-39, far above any floor) rather than deriving new ones.
CONSTRUCT_ISOLATION = {
    "md1_md2": "#c8831c",
    "md2_md3": "#2845bd",
}

# Merged dict for callers that just want "the color for set X" regardless of
# group (e.g. a single shared lookup in sim_registry.py). Group membership is
# still what determines which sets are safe to co-plot -- see the CONTRAST
# note above: this dict does NOT guarantee all 13 keys are mutually
# distinguishable if drawn together in one figure (they were never validated
# that way, since the analysis scripts don't do that).
CONSTRUCT = {**CONSTRUCT_GO, **CONSTRUCT_NOGO, **CONSTRUCT_ISOLATION}

# NOTE on contrast WARNs (both groups): #E69F00/#56B4E9/#CC79A7 (Go group) and
# #eb79eb/#1cc895 (non-Go group) sit below the 3.0:1 contrast-vs-surface
# floor. Per the dataviz skill this is a WARN, not a FAIL, but it obligates
# visible direct labels or a table view wherever these appear -- never
# color-alone identification. All current comparison figures in this repo
# already label bars/lines directly, so this is satisfied in practice, not a
# TODO.

# ---------------------------------------------------------------------------
SEQUENTIAL_STOPS = ["#CDE2FB", "#86B6EF", "#3987E5", "#2A78D6", "#256ABF", "#184F95", "#0D366B"]
DIVERGING_STOPS = ["#0D366B", "#2A78D6", "#86B6EF", "#F0EFEC", "#EC835A", "#E34948", "#8A1F1F"]
DIVERGING_MIDPOINT_LIGHT = "#F0EFEC"
DIVERGING_MIDPOINT_DARK = "#383835"
STATUS = {"good": "#0CA30C", "warning": "#FAB219", "serious": "#EC835A", "critical": "#D03B3B"}
INK = {"primary": "#0B0B0B", "secondary": "#52514E", "muted": "#898781",
       "gridline": "#E1E0D9", "baseline": "#C3C2B7"}
