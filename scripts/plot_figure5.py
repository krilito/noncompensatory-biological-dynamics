"""Figure 5 quantitative panels from source_data (Python authority).\n\nOwner-template composition is separate artwork."""
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
    CR, CNR = S.C_RESPONDER, S.C_NONRESPONDER
    STEM = "Figure5"
    env_all = pd.read_csv(SD / "figure5.csv")
    env = env_all.copy()
    axis_dir = pd.read_csv(SD / "figure5_axis_direction.csv")
    phen = pd.read_csv(SD / "figure5_phenotype.csv")
    imv_samples = pd.read_csv(SD / "figure5_imvigor_samples.csv")

    AXIS_MIN, AXIS_MAX = 0.05, 1.00      # frozen; verified to contain every CI endpoint

    fig = plt.figure(figsize=(S.W_DOUBLE, 5.2))
    # two independent gridspecs so the forest gutters do not distort the bottom row
    gs_top = GridSpec(1, 1, figure=fig, left=0.205, right=0.775, top=0.925, bottom=0.545)
    # Panel d carries four two-line tick labels ("X" + its ANOVA P) and needs the extra width;
    # the wider gutter keeps panel c's title clear of panel d's label.
    gs_bot = GridSpec(1, 3, figure=fig, left=0.085, right=0.995, top=0.395, bottom=0.135,
                      wspace=0.66, width_ratios=[0.98, 0.86, 1.40])

    # =====================================================================================
    # a — operating-envelope AUC on ONE common axis with COMPLETE intervals
    # =====================================================================================
    ax = fig.add_subplot(gs_top[0, 0])
    env = env.sort_values("row_order", ascending=False).reset_index(drop=True)

    for i, r in env.iterrows():
        crosses = bool(r["ci_crosses_0.5"])
        # interval: unbroken line, no detached caps
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=S.C_MID, lw=S.LW_CI,
                solid_capstyle="butt", zorder=2)
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], marker="|", ms=4.5,
                color=S.C_MID, linestyle="none", zorder=2)
        # point estimate: shape encodes analysis unit; fill encodes whether the CI resolves
        mk = "h" if "pseudo-bulk" in str(r["substrate"]) else "o"
        ax.plot(r["auc"], i, marker=mk, ms=S.MS_POINT + 0.8, color=S.C_INK,
                markerfacecolor=("white" if crosses else S.C_INK),
                markeredgewidth=1.2, linestyle="none", zorder=3)

    S.refline(ax, 0.5, axis="x")
    ax.annotate("0.5", xy=(0.5, len(env) - 0.32), fontsize=S.FS_ANNOT, color=S.C_REF,
                ha="center", va="bottom")

    labels = [f"{r['dataset']}\n{r['timing']} · {r['substrate']} · $n$ = {int(r['n'])}"
              for _, r in env.iterrows()]
    ax.set_yticks(range(len(env)))
    ax.set_yticklabels(labels, fontsize=S.FS_MIN, linespacing=1.35)
    ax.set_xlim(AXIS_MIN, AXIS_MAX)
    ax.set_ylim(-0.65, len(env) - 0.15)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.set_xlabel("Response AUC of the frozen boundary (complete 95% CI)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "a", dx=-0.235)
    S.panel_title(ax, "GSE78220 primary is near null; comparison rows are bounded")

    ROLE_TEXT = {
        "PRIMARY_NEAR_NULL": "primary near null",
        "COMPARISON_ONLY": "comparison only",
        "NEAR_NULL": "near null",
        "CARRIER_CONTEXT_UNRESOLVED_AUC": "carrier context, unresolved",
        "RESOLVED_PREDICTIVE_INVERSION_SEMANTICS_UNRESOLVED":
            "resolved inversion, semantics unresolved",
    }
    for i, r in env.iterrows():
        ax.annotate(ROLE_TEXT.get(str(r["evidence_role"]), str(r["evidence_role"]).replace("_", " ").lower()), xy=(1.025, i),
                    xycoords=("axes fraction", "data"), fontsize=S.FS_ANNOT,
                    color=(S.C_INK if not r["ci_crosses_0.5"] else S.C_MID),
                    ha="left", va="center", linespacing=1.25, annotation_clip=False)

    # =====================================================================================
    # b — GSE78220 baseline: the axis separation that supports the boundary collapses
    # =====================================================================================
    ax = fig.add_subplot(gs_bot[0, 0])
    order = ["T axis", "E axis", "X axis", "C axis", "Frozen boundary score"]
    short = ["T", "E", "X", "C", "score"]
    on = axis_dir[(axis_dir["cohort"] == "GSE91061") & (axis_dir["timepoint"] == "on-treatment")]
    pre78 = axis_dir[axis_dir["cohort"] == "GSE78220"]

    xs = np.arange(len(order))
    v_on = [float(on[on["metric"] == m]["PRCR_minus_PD"].iloc[0]) for m in order]
    v_78 = [float(pre78[pre78["metric"] == m]["PRCR_minus_PD"].iloc[0]) for m in order]

    ax.bar(xs - 0.20, v_on, width=0.38, color=S.C_COHORT["GSE91061"], label="GSE91061 on-treatment")
    ax.bar(xs + 0.20, v_78, width=0.38, color="white", edgecolor=S.C_COHORT["GSE78220"],
           linewidth=1.1, hatch="////", label="GSE78220 pretreatment")
    S.refline(ax, 0.0, axis="y")
    ax.set_xticks(xs)
    ax.set_xticklabels(short, fontsize=S.FS_TICK)
    ax.set_ylabel("Responder − non-responder\nmean (axis z units)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "b", dx=-0.30)
    S.panel_title(ax, "At baseline the separation collapses")
    ax.legend(fontsize=S.FS_MIN, loc="lower left", bbox_to_anchor=(-0.02, -0.46),
              frameon=False, handlelength=1.1)

    # =====================================================================================
    # c — IMvigor210: the composite reversed
    # =====================================================================================
    ax = fig.add_subplot(gs_bot[0, 1])
    sc = imv_samples.dropna(subset=["response_strict"])
    groups = [("R", sc[sc["response_strict"] == "R"]["boundary_score"].to_numpy(float)),
              ("NR", sc[sc["response_strict"] == "NR"]["boundary_score"].to_numpy(float))]
    rng = np.random.default_rng(5)
    for j, (lab, vals) in enumerate(groups):
        c = CR if lab == "R" else CNR
        jit = rng.uniform(-0.17, 0.17, size=len(vals))
        ax.plot(np.full(len(vals), j) + jit, vals, marker="o" if lab == "R" else "D",
                ms=1.9, linestyle="none", color=c, alpha=0.5, markeredgewidth=0, zorder=2)
        ax.plot([j - 0.30, j + 0.30], [np.median(vals)] * 2, color=S.C_INK, lw=1.5, zorder=4)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"CR/PR\n$n$ = {len(groups[0][1])}", f"PD\n$n$ = {len(groups[1][1])}"],
                       fontsize=S.FS_TICK, linespacing=1.3)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Frozen boundary score", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "c", dx=-0.28)
    S.panel_title(ax, "IMvigor210 comparison context")
    imv_row = env_all[env_all["dataset"].astype(str) == "IMvigor210"]
    if not imv_row.empty:
        r0 = imv_row.iloc[0]
        auc_txt = f"AUC {float(r0['auc']):.3f} ({float(r0['ci_low']):.3f}–{float(r0['ci_high']):.3f})"
    else:
        auc_txt = "AUC unavailable in source_data/figure5.csv"
    ax.annotate(auc_txt, xy=(0.5, -0.235), xycoords="axes fraction",
                fontsize=S.FS_ANNOT, color=S.C_MID, ha="center", va="top")

    # =====================================================================================
    # d — but the axes did not fail together
    # =====================================================================================
    ax = fig.add_subplot(gs_bot[0, 2])
    phenos = ["inflamed", "excluded", "desert"]
    hatches = ["", "///", "..."]
    axes_order = ["T", "E", "X", "C"]
    xs = np.arange(len(axes_order))
    w = 0.26
    for k, ph in enumerate(phenos):
        sub = phen[phen["phenotype"] == ph].set_index("axis").loc[axes_order]
        means = sub["mean_score"].to_numpy(float)
        se = (sub["std"] / np.sqrt(sub["n"])).to_numpy(float)
        ax.bar(xs + (k - 1) * w, means, width=w * 0.92, yerr=se, capsize=1.6,
               color="white", edgecolor=S.C_INK, linewidth=0.8, hatch=hatches[k],
               error_kw=dict(lw=0.8, ecolor=S.C_MID), label=ph, zorder=2)

    S.refline(ax, 0.0, axis="y")
    pvals = [float(phen[phen["axis"] == a]["anova_p"].iloc[0]) for a in axes_order]
    ax.set_xticks(xs)
    # "P =" is factored into the axis label: repeating it four times leaves the values touching.
    ax.set_xticklabels([f"{a}\n{pv:.3f}" for a, pv in zip(axes_order, pvals)],
                       fontsize=S.FS_MIN, linespacing=1.35)
    ax.set_xlabel("MELD axis · phenotype ANOVA $P$", fontsize=S.FS_AXIS_LABEL)
    for tick, pv in zip(ax.get_xticklabels(), pvals):
        tick.set_color(S.C_INK if pv < 0.05 else S.C_MID)
    ax.set_ylim(-0.42, 0.34)
    ax.set_ylabel("Mean axis score ± SE", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "d", dx=-0.22)
    S.panel_title(ax, "The axes did not fail together")
    ax.legend(fontsize=S.FS_MIN, loc="lower left", bbox_to_anchor=(-0.06, -0.46),
              frameon=False, ncol=3, handlelength=1.1, columnspacing=0.8)

    handles = [
        Line2D([], [], marker="o", ls="none", color=S.C_INK, markerfacecolor="white",
               markeredgewidth=1.2, label="interval spans 0.5"),
        Line2D([], [], marker="o", ls="none", color=S.C_INK, markerfacecolor=S.C_INK,
               markeredgewidth=1.2, label="interval excludes 0.5"),
        Line2D([], [], marker="h", ls="none", color=S.C_INK, markerfacecolor="white",
               markeredgewidth=1.2, label="CD45+ pseudo-bulk substrate"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.53, 1.008),
               ncol=3, fontsize=S.FS_LEGEND, frameon=False)

    paths = S.export(fig, OUTDIR, STEM)
    print(json.dumps({k: str(v.relative_to(ROOT).as_posix()) for k, v in paths.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
