#!/usr/bin/env python3
"""
Aggregate final_design_stats.csv across all PARP14 BindCraft targets and plot
cross-target comparisons. Safe to re-run at any time -- targets with no
results yet (job not run / still running / 0 accepted designs) are reported
as such and skipped in the plots, not treated as errors.

Usage (from bindcraft_md/, in any env with pandas+matplotlib -- the BindCraft
conda env already has both):
    python analyze_campaign.py
"""
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DESIGNS_DIR = os.path.join(HERE, "designs")
SUMMARY_CSV = os.path.join(HERE, "targets_summary.csv")
OUT_DIR = os.path.join(HERE, "figures")

# validated categorical palette (dataviz skill) -- fixed order, 2 series here
COLOR_BLOCK = "#2a78d6"   # slot 1 blue
COLOR_CLAMP = "#eb6834"   # slot 2 orange
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"


def goal_category(target_name):
    return "clamp" if target_name.startswith("clamp_") else "block"


def load_targets():
    summary = pd.read_csv(SUMMARY_CSV)
    return summary["target"].tolist(), dict(zip(summary["target"], summary["goal"]))


def load_all_final_stats(targets):
    rows = []
    status = []
    for t in targets:
        csv_path = os.path.join(DESIGNS_DIR, t, "final_design_stats.csv")
        if not os.path.exists(csv_path):
            status.append((t, "not run yet"))
            continue
        df = pd.read_csv(csv_path)
        if df.empty:
            status.append((t, "ran, 0 accepted designs"))
            continue
        df = df.copy()
        df["Target"] = t
        df["Goal"] = goal_category(t)
        rows.append(df)
        status.append((t, f"{len(df)} accepted design(s)"))
    combined = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return combined, status


def print_status(status):
    print("\n=== Campaign status ===")
    for t, s in status:
        print(f"  {t:24s} {s}")
    print()


def plot_accepted_counts(combined, targets, out_path):
    counts = {t: 0 for t in targets}
    if not combined.empty:
        for t, n in combined.groupby("Target").size().items():
            counts[t] = n
    goals = [goal_category(t) for t in targets]
    colors = [COLOR_BLOCK if g == "block" else COLOR_CLAMP for g in goals]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(range(len(targets)), [counts[t] for t in targets], color=colors, width=0.6)
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets, rotation=45, ha="right", fontsize=8, color=INK)
    ax.set_ylabel("Accepted designs", color=INK)
    ax.set_title("Accepted BindCraft designs per target", color=INK)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLOR_BLOCK, label="block"),
        plt.Rectangle((0, 0), 1, 1, color=COLOR_CLAMP, label="clamp"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_best_iptm(combined, targets, out_path):
    if combined.empty or "Average_i_pTM" not in combined.columns:
        print("Skipping best-i_pTM plot -- no accepted designs with i_pTM yet.")
        return
    best = combined.groupby("Target")["Average_i_pTM"].max()
    ranked = best.reindex([t for t in targets if t in best.index]).sort_values(ascending=False)
    colors = [COLOR_BLOCK if goal_category(t) == "block" else COLOR_CLAMP for t in ranked.index]

    fig, ax = plt.subplots(figsize=(8, max(3, 0.4 * len(ranked))))
    ax.barh(range(len(ranked)), ranked.values, color=colors, height=0.6)
    ax.set_yticks(range(len(ranked)))
    ax.set_yticklabels(ranked.index, fontsize=8, color=INK)
    ax.invert_yaxis()
    ax.set_xlabel("Best Average_i_pTM (top accepted design)", color=INK)
    ax.set_title("Best interface pTM per target (higher = more confident interface)", color=INK)
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLOR_BLOCK, label="block"),
        plt.Rectangle((0, 0), 1, 1, color=COLOR_CLAMP, label="clamp"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_iptm_vs_ipae(combined, out_path):
    needed = {"Average_i_pTM", "Average_i_pAE"}
    if combined.empty or not needed.issubset(combined.columns):
        print("Skipping i_pTM-vs-i_pAE scatter -- no accepted designs yet.")
        return
    fig, ax = plt.subplots(figsize=(6, 6))
    for goal, color in (("block", COLOR_BLOCK), ("clamp", COLOR_CLAMP)):
        sub = combined[combined["Goal"] == goal]
        if sub.empty:
            continue
        ax.scatter(
            sub["Average_i_pTM"], sub["Average_i_pAE"],
            s=36, color=color, alpha=0.85, edgecolors="white", linewidths=0.6,
            label=goal, zorder=3,
        )
    ax.set_xlabel("Average_i_pTM (interface confidence, higher better)", color=INK)
    ax.set_ylabel("Average_i_pAE (interface error, lower better)", color=INK)
    ax.set_title("Accepted designs: interface confidence vs error", color=INK)
    ax.grid(color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    targets, goals_from_summary = load_targets()
    combined, status = load_all_final_stats(targets)
    print_status(status)

    campaign_csv = os.path.join(OUT_DIR, "campaign_summary.csv")
    if not combined.empty:
        combined.to_csv(campaign_csv, index=False)
        print(f"Wrote {campaign_csv} ({len(combined)} rows across "
              f"{combined['Target'].nunique()} target(s))")
    else:
        print("No accepted designs anywhere yet -- nothing to write to campaign_summary.csv")

    plot_accepted_counts(combined, targets, os.path.join(OUT_DIR, "accepted_counts.png"))
    plot_best_iptm(combined, targets, os.path.join(OUT_DIR, "best_iptm_per_target.png"))
    plot_iptm_vs_ipae(combined, os.path.join(OUT_DIR, "iptm_vs_ipae.png"))


if __name__ == "__main__":
    sys.exit(main())
