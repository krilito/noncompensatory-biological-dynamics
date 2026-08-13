"""Figure 4 quantitative panels from source_data (Python authority)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_style as S

SD = ROOT / "source_data"
OUTDIR = ROOT / "figures" / "reproduced"


def main() -> int:
    S.apply()
    STEM = "Figure4"
    share = pd.read_csv(SD / "figure4.csv")
    lin = pd.read_csv(SD / "figure4_lineage.csv")
    auc = pd.read_csv(SD / "figure4_auc.csv")

    AXES = ["T", "E", "X", "C"]
    t = share[share["axis"] == "T"].sort_values("axis_signal_share", ascending=True).reset_index(drop=True)

    fig = plt.figure(figsize=(S.W_DOUBLE, 4.15))
    gs = GridSpec(1, 4, figure=fig, left=0.105, right=0.985, top=0.885, bottom=0.315,
                  wspace=0.80, width_ratios=[1.0, 0.55, 0.75, 1.0])

    # =====================================================================================
    # a — which lineages carry the T coordinate, within CD45+
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 0])
    ys = np.arange(len(t))
    ax.barh(ys, t["axis_signal_share"] * 100, color=S.C_MID, height=0.66, zorder=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([str(x).replace("_", " ") for x in t["lineage"]], fontsize=S.FS_MIN)
    ax.set_xlabel("Share of T-axis signal (%)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "a", dx=-0.52)
    S.panel_title(ax, "Carriers of T, within CD45+")

    # =====================================================================================
    # b — how many cells each share rests on
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(ys, t["n_cells"], color="white", edgecolor=S.C_MID, linewidth=0.9, height=0.66, zorder=2)
    ax.set_xscale("log")
    ax.set_yticks(ys)
    ax.set_yticklabels([])
    ax.set_xlabel("Cells (log)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "b", dx=-0.12)
    S.panel_title(ax, "Cells behind each share")

    # =====================================================================================
    # c — descriptive lineage x axis response contrast (effect size, no significance stars)
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 2])
    lineages = list(t["lineage"])          # same bottom-to-top order as panels a and b
    piv = (lin.pivot(index="lineage", columns="axis", values="cohen_d")
              .reindex(index=lineages, columns=AXES))
    vmax = float(np.nanmax(np.abs(piv.to_numpy(float))))
    cmap = plt.get_cmap("PuOr_r").copy()
    cmap.set_bad(S.C_PANEL_BG)             # lineage not present in the response table
    im = ax.imshow(np.ma.masked_invalid(piv.to_numpy(float)), cmap=cmap,
                   vmin=-vmax, vmax=vmax, aspect="auto", origin="lower")
    ax.set_xticks(range(len(AXES)))
    ax.set_xticklabels(AXES, fontsize=S.FS_TICK)
    ax.set_yticks(range(len(lineages)))
    ax.set_yticklabels([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cax = ax.inset_axes([0.0, -0.165, 1.0, 0.038])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Cohen's $d$, R − NR", fontsize=S.FS_MIN, labelpad=1.5)
    cb.ax.tick_params(labelsize=S.FS_MIN, length=1.5, pad=1)
    cb.outline.set_visible(False)
    S.panel_label(ax, "c", dx=-0.10)
    S.panel_title(ax, "Contrast per lineage")

    # =====================================================================================
    # d — does the CD45+ projection rank response?
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 3])
    # Accept both historical producer labels and short export labels.
    label_map = {
        "all_samples": "all samples",
        "pretreatment_only": "pretreatment",
        "on_treatment_only": "on-treatment",
        "pretreatment": "pretreatment",
        "on-treatment": "on-treatment",
        "all samples": "all samples",
    }
    keep = [a for a in ("all_samples", "pretreatment_only", "on_treatment_only",
                        "pretreatment", "on-treatment", "all samples")
            if a in set(auc["analysis"].astype(str))]
    # Prefer a stable order: pretreatment then on-treatment when only two rows.
    preferred = ["pretreatment", "pretreatment_only", "on-treatment", "on_treatment_only",
                 "all samples", "all_samples"]
    ordered = [a for a in preferred if a in keep]
    for a in keep:
        if a not in ordered:
            ordered.append(a)
    sub = auc.set_index("analysis").loc[ordered].reset_index()
    sub = sub.iloc[::-1].reset_index(drop=True)
    nice = label_map
    for i, r in sub.iterrows():
        ax.plot([r["CI_low"], r["CI_high"]], [i, i], color=S.C_MID, lw=S.LW_CI,
                solid_capstyle="butt", zorder=2)
        ax.plot([r["CI_low"], r["CI_high"]], [i, i], marker="|", ms=4.5, color=S.C_MID,
                linestyle="none", zorder=2)
        ax.plot(r["AUC"], i, marker="h", ms=S.MS_POINT + 0.6, color=S.C_INK,
                markerfacecolor="white", markeredgewidth=1.2, linestyle="none", zorder=3)

    S.refline(ax, 0.5, axis="x")
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([
        f"{nice.get(str(r['analysis']), str(r['analysis']))}\n$n$ = {int(r['n'])}"
        for _, r in sub.iterrows()
    ], fontsize=S.FS_MIN, linespacing=1.3)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.6, len(sub) - 0.3)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Response AUC (95% CI)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "d", dx=-0.48)
    S.panel_title(ax, "Every interval spans 0.5")

    # Long substrate qualifications live in the figure caption, not in the panel.

    paths = S.export(fig, OUTDIR, STEM)
    print(json.dumps({k: str(v.relative_to(ROOT).as_posix()) for k, v in paths.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
