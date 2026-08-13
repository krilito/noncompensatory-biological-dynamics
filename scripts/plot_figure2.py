"""Figure plotter for QDA.

Reads repository-relative source_data tables and writes PDF/PNG under figures/reproduced/.
Scientific panel authority: Python (canonical_interface quantitative producers).
Formal multi-panel layouts for Figure 2 and Figure 5 remain owner artwork;
these scripts regenerate the quantitative scientific content from frozen source tables.
"""
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
    STEM = "Figure2"
    scores = pd.read_csv(SD / "figure2.csv")
    auc = pd.read_csv(SD / "figure2_auc.csv")

    # Map display names used in panel data
    name_map = {
        "MGH": "MGH / GSE115821 + GSE168204",
        "PRJEB23709": "Gide / PRJEB23709",
        "Gide / PRJEB23709": "Gide / PRJEB23709",
        "MORRISON-1-public": "MORRISON-1-public",
    }
    if scores["cohort"].isin(["MGH", "PRJEB23709"]).any():
        scores = scores.copy()
        scores["cohort"] = scores["cohort"].map(lambda x: name_map.get(x, x))

    COHORT_ORDER = ["Gide / PRJEB23709", "MORRISON-1-public", "MGH / GSE115821 + GSE168204"]
    # accept short MGH labels in scores
    if not scores["cohort"].isin(COHORT_ORDER).any():
        # try alternate short names from export
        alt = {"MGH": "MGH / GSE115821 + GSE168204", "PRJEB23709": "Gide / PRJEB23709"}
        scores["cohort"] = scores["cohort"].replace(alt)
    COHORT_SHORT = {"Gide / PRJEB23709": "PRJEB23709", "MORRISON-1-public": "MORRISON-1",
                    "MGH / GSE115821 + GSE168204": "MGH", "MGH": "MGH"}

    # normalize scores cohorts for plotting
    present = [c for c in COHORT_ORDER if (scores["cohort"] == c).any()]
    if not present:
        # panel export uses short names
        COHORT_ORDER = list(scores["cohort"].dropna().unique())
        present = COHORT_ORDER
        COHORT_SHORT = {c: str(c).split()[0][:12] for c in present}

    scores = scores[scores["cohort"].isin(present)].copy()

    fig = plt.figure(figsize=(S.W_DOUBLE, 3.1))
    gs = GridSpec(1, 3, figure=fig, left=0.065, right=0.995, top=0.845, bottom=0.215,
                  wspace=0.62, width_ratios=[1.25, 1.15, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(3)
    xt, xl = [], []
    for i, coh in enumerate(present):
        sub = scores[scores["cohort"] == coh]
        for k, (lab, sel, c, mk) in enumerate([
                ("R", sub["response_binary_responder1"] == 1, CR, "o"),
                ("NR", sub["response_binary_responder1"] == 0, CNR, "D")]):
            v = sub.loc[sel, "boundary_score"].to_numpy(float)
            x = i * 2.4 + k * 0.85
            if len(v):
                ax.plot(x + rng.uniform(-0.19, 0.19, len(v)), v, marker=mk, ms=2.6,
                        linestyle="none", color=c, alpha=0.75, markeredgewidth=0, zorder=2)
                ax.plot([x - 0.30, x + 0.30], [np.median(v)] * 2, color=S.C_INK, lw=1.4, zorder=4)
            xt.append(x)
            xl.append(f"{lab}\n{len(v)}")
        ax.annotate(COHORT_SHORT.get(coh, str(coh)), xy=(i * 2.4 + 0.42, -0.245),
                    xycoords=("data", "axes fraction"), fontsize=S.FS_ANNOT,
                    color=S.C_INK, ha="center", va="top")
    thr = float(scores["frozen_threshold"].iloc[0]) if "frozen_threshold" in scores.columns else 0.735
    ax.axhline(thr, color=S.C_REF, lw=S.LW_REF, ls=(0, (2.5, 2.0)), zorder=1)
    ax.annotate(f"- - -  frozen threshold {thr:g}", xy=(0.995, 0.985), xycoords="axes fraction",
                fontsize=S.FS_ANNOT, color=S.C_REF, ha="right", va="top")
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=S.FS_MIN, linespacing=1.25)
    ax.set_ylim(-4.15, 3.70)
    ax.set_yticks([-4, -3, -2, -1, 0, 1, 2])
    ax.set_ylabel("Frozen boundary score", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "a", dx=-0.145)
    S.panel_title(ax, f"On-treatment state, {len(scores)} biopsies")

    ax = fig.add_subplot(gs[0, 1])
    rows = auc.iloc[::-1].reset_index(drop=True)
    for i, r in rows.iterrows():
        is_patient = str(r.get("analysis_unit", "")) == "patient"
        if pd.notna(r.get("ci_low")) and str(r.get("ci_low")) != "":
            ax.plot([float(r["ci_low"]), float(r["ci_high"])], [i, i], color=S.C_MID,
                    lw=S.LW_CI, solid_capstyle="butt", zorder=2)
            ax.plot([float(r["ci_low"]), float(r["ci_high"])], [i, i], marker="|",
                    ms=4.5, color=S.C_MID, linestyle="none", zorder=2)
        ax.plot(float(r["auc_canonical"]), i,
                marker="s" if is_patient else "o", ms=S.MS_POINT + (0.7 if is_patient else 0),
                color=S.C_INK, markerfacecolor=S.C_INK if is_patient else "white",
                markeredgewidth=1.2, linestyle="none", zorder=3)
    S.refline(ax, 0.5, axis="x")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([
        f"{r['cohort']}\nn = {int(r['n'])} "
        f"{'patients' if str(r.get('analysis_unit')) == 'patient' else 'samples'}"
        for _, r in rows.iterrows()
    ], fontsize=S.FS_MIN, linespacing=1.3)
    ax.set_xlim(0.3, 1.02)
    ax.set_ylim(-0.65, len(rows) - 0.25)
    ax.set_xlabel("Canonical response AUC (95% CI)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "b", dx=-0.52)
    S.panel_title(ax, "Two core cohorts; MGH is supportive")

    ax = fig.add_subplot(gs[0, 2])
    drift = auc[auc["analysis_unit"].astype(str).str.contains("biopsy|sample", case=False, na=False)].iloc[::-1].reset_index(drop=True)
    if drift.empty:
        drift = auc.iloc[::-1].reset_index(drop=True)
    plotted = 0
    for i, r in drift.iterrows():
        if pd.isna(r.get("auc_historical")) or str(r.get("auc_historical")).strip() == "":
            continue
        h, c = float(r["auc_historical"]), float(r["auc_canonical"])
        ax.plot([h, c], [plotted, plotted], color=S.C_MID, lw=1.0, zorder=2)
        ax.plot(h, plotted, marker="o", ms=S.MS_SMALL + 0.6, color=S.C_GREY,
                markerfacecolor="white", markeredgewidth=1.0, linestyle="none", zorder=3)
        ax.plot(c, plotted, marker="o", ms=S.MS_SMALL + 1.2, color=S.C_INK,
                markerfacecolor=S.C_INK, markeredgewidth=0, linestyle="none", zorder=4)
        ax.annotate(f"{c - h:+.3f}", xy=(0.985, plotted), xycoords=("axes fraction", "data"),
                    fontsize=S.FS_MIN, color=S.C_MID, ha="right", va="center")
        plotted += 1
    ax.set_yticks(range(plotted))
    labels = []
    for _, r in drift.iterrows():
        if pd.isna(r.get("auc_historical")) or str(r.get("auc_historical")).strip() == "":
            continue
        labels.append(str(r["cohort"]).replace(" sample-level", ""))
    ax.set_yticklabels(labels, fontsize=S.FS_MIN)
    ax.set_xlim(0.60, 0.95)
    ax.set_ylim(-0.6, max(plotted - 0.15, 0.5))
    ax.set_xticks([0.65, 0.75, 0.85])
    ax.set_xlabel("Response AUC", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "c", dx=-0.40)
    S.panel_title(ax, "Rebuild moved each cohort")

    handles = [
        Line2D([], [], marker="o", ls="none", color=CR, markerfacecolor=CR, label="Responder"),
        Line2D([], [], marker="D", ls="none", color=CNR, markerfacecolor=CNR, label="Non-responder"),
        Line2D([], [], marker="o", ls="none", color=S.C_INK, markerfacecolor="white",
               markeredgewidth=1.2, label="sample-level"),
        Line2D([], [], marker="s", ls="none", color=S.C_INK, markerfacecolor=S.C_INK,
               markeredgewidth=1.2, label="patient-level"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.015),
               ncol=4, fontsize=S.FS_LEGEND, frameon=False)
    paths = S.export(fig, OUTDIR, STEM)
    print(json.dumps({k: str(v.relative_to(ROOT).as_posix()) for k, v in paths.items()}, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
