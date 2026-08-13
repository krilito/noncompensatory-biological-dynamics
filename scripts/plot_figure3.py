"""Figure 3 quantitative panels from source_data (Python authority)."""
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
    CHELP, CHURT = S.C_HELPED, S.C_HURT
    STEM = "Figure3"
    paired = pd.read_csv(SD / "figure3.csv")
    b2 = pd.read_csv(SD / "figure3_b2_summary.csv").iloc[0]
    b3 = pd.read_csv(SD / "figure3_b3.csv")
    loss = pd.read_csv(SD / "figure3_b4_patients.csv")
    b4 = pd.read_csv(SD / "figure3_b4_summary.csv").iloc[0]
    signflip = pd.read_csv(SD / "figure3_signflip.csv")
    influence = pd.read_csv(SD / "figure3_influence.csv")

    fig = plt.figure(figsize=(S.W_DOUBLE, 5.5))
    gs = GridSpec(2, 3, figure=fig, hspace=0.62, wspace=0.42,
                  left=0.075, right=0.985, top=0.93, bottom=0.085,
                  width_ratios=[1.0, 1.05, 1.15])

    # =====================================================================================
    # a — data structure: 16 therapy records from 15 patients
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 0])
    therapies = ["PD1", "ipiPD1"]
    patients = sorted(paired["patient_id"].unique())
    ypos = {p: i for i, p in enumerate(patients)}

    dual = [p for p in patients if (paired["patient_id"] == p).sum() > 1]
    for _, r in paired.iterrows():
        x = therapies.index(r["therapy_type"])
        isR = r["y_true"] == 1
        c = CR if isR else CNR
        ax.plot(x, ypos[r["patient_id"]],
                marker="o" if isR else "D", ms=S.MS_SMALL + 0.8,
                color=c, markerfacecolor="white", markeredgewidth=1.0,
                linestyle="none", zorder=3)

    # connect the dual-record patient's two records so the shared patient is visible
    for p in dual:
        ax.plot([0, 1], [ypos[p], ypos[p]], color=S.C_INK, lw=1.0, zorder=2)

    ax.set_xticks(range(len(therapies)))
    ax.set_xticklabels(therapies, fontsize=S.FS_TICK)
    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels([(f"P{p} *" if p in dual else f"P{p}") for p in patients],
                       fontsize=S.FS_MIN)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-1.9, len(patients) - 0.3)
    ax.set_xlabel("Therapy", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "a", dx=-0.30)
    S.panel_title(ax, "16 therapy records, 15 patients")
    ax.annotate("* two therapies, discordant labels",
                xy=(0.0, -0.235), xycoords="axes fraction", fontsize=S.FS_ANNOT,
                color=S.C_INK, ha="left", va="top")
    ax.annotate("B2, B3: record unit (16)\nB4: patient unit (15)",
                xy=(0.0, -0.325), xycoords="axes fraction", fontsize=S.FS_ANNOT,
                color=S.C_MID, ha="left", va="top", linespacing=1.3)

    # =====================================================================================
    # b — B2 movement: paired PRE -> on-treatment boundary score
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 1])
    for _, r in paired.iterrows():
        isR = r["y_true"] == 1
        c = CR if isR else CNR
        ax.plot([0, 1], [r["PRE_boundary"], r["EDT_boundary"]],
                color=c, lw=S.LW_LINK, alpha=0.75,
                ls="-" if isR else (0, (3, 1.4)), zorder=2)
        ax.plot(0, r["PRE_boundary"], marker="o" if isR else "D", ms=S.MS_SMALL,
                color=c, markerfacecolor="white", markeredgewidth=0.9, zorder=3)
        ax.plot(1, r["EDT_boundary"], marker="o" if isR else "D", ms=S.MS_SMALL,
                color=c, markerfacecolor=c, markeredgewidth=0.0, zorder=3)

    ax.axhline(0.735, color=S.C_REF, lw=S.LW_REF, ls=(0, (2.5, 2.0)), zorder=1)
    ax.annotate("frozen\nthreshold\n0.735", xy=(1.015, 0.735),
                xycoords=("axes fraction", "data"), fontsize=S.FS_ANNOT,
                color=S.C_REF, va="center", ha="left", linespacing=1.2,
                annotation_clip=False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["pre", "on-treatment"], fontsize=S.FS_TICK)
    ax.set_xlim(-0.22, 1.22)
    ax.set_ylabel("Frozen boundary score", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "b")
    S.panel_title(ax, f"Displacement in {int(b2['n_delta_gt0'])} of {int(b2['n_records'])} records")
    ax.annotate(f"paired Wilcoxon $P$ = {b2['wilcoxon_p']:.4f}",
                xy=(0.5, -0.19), xycoords="axes fraction", fontsize=S.FS_ANNOT,
                color=S.C_MID, ha="center", va="top")

    # =====================================================================================
    # c — B2 direction: record-level transition cosine
    # =====================================================================================
    ax = fig.add_subplot(gs[0, 2])
    cos = paired["cosine_to_reference"].to_numpy(float)
    isR = paired["y_true"].to_numpy(int) == 1
    rng = np.random.default_rng(11)          # deterministic jitter
    jit = rng.uniform(-0.16, 0.16, size=len(cos))
    for m, sel, c in [("o", isR, CR), ("D", ~isR, CNR)]:
        ax.plot(cos[sel], jit[sel], marker=m, ms=S.MS_SMALL + 0.6, linestyle="none",
                color=c, markerfacecolor="white", markeredgewidth=0.95, zorder=3)

    med = float(b2["median_record_level_cosine"])
    ax.axvline(med, color=S.C_INK, lw=1.2, zorder=2)
    ax.annotate(f"median {med:.3f}", xy=(med, 0.30), fontsize=S.FS_ANNOT,
                color=S.C_INK, ha="center", va="bottom")
    ax.plot([b2["cohort_transition_cosine"]], [-0.30], marker="^", ms=S.MS_POINT,
            color=S.C_INK, linestyle="none", zorder=4)
    ax.annotate(f"cohort-mean vector {b2['cohort_transition_cosine']:.3f}",
                xy=(b2["cohort_transition_cosine"], -0.36), fontsize=S.FS_ANNOT,
                color=S.C_MID, ha="right", va="top")
    S.refline(ax, 0.0, axis="x")
    ax.set_xlim(-0.35, 1.12)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xlabel("Transition cosine to frozen reference\n(dimensionless)", fontsize=S.FS_AXIS_LABEL)
    ax.spines["left"].set_visible(False)
    S.subtle(ax)
    S.panel_label(ax, "c")
    S.panel_title(ax, "Direction, per record")

    # =====================================================================================
    # d — B3 response specificity (NOT SUPPORTED)
    # =====================================================================================
    ax = fig.add_subplot(gs[1, 0])
    rows = b3.iloc[::-1].reset_index(drop=True)
    for i, r in rows.iterrows():
        if pd.notna(r.get("ci_low")) and str(r.get("ci_low")) != "":
            ax.plot([float(r["ci_low"]), float(r["ci_high"])], [i, i],
                    color=S.C_MID, lw=S.LW_CI, solid_capstyle="butt", zorder=2)
            ax.plot([float(r["ci_low"]), float(r["ci_high"])], [i, i], marker="|",
                    ms=4, color=S.C_MID, linestyle="none", zorder=2)
        ax.plot(r["auc"], i, marker="v", ms=S.MS_POINT, color=S.C_INK,
                markerfacecolor="white", markeredgewidth=1.1, linestyle="none", zorder=3)

    S.refline(ax, 0.5, axis="x")
    ax.annotate("no discrimination", xy=(0.5, 1.72), fontsize=S.FS_ANNOT,
                color=S.C_REF, ha="center", va="bottom")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(x) for x in rows["readout"]], fontsize=S.FS_TICK)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.6, 1.9)
    ax.set_xlabel("Response AUC (record level)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "d", dx=-0.30)
    S.panel_title(ax, "Response specificity not supported")
    ax.annotate(f"11 responder / 5 non-responder records\nMann–Whitney $P$ = "
                f"{float(b3.iloc[0]['mannwhitney_p']):.3f}",
                xy=(0.5, -0.30), xycoords="axes fraction", fontsize=S.FS_ANNOT,
                color=S.C_MID, ha="center", va="top", linespacing=1.25)

    # =====================================================================================
    # e — B4 patient-level delta_L
    # =====================================================================================
    ax = fig.add_subplot(gs[1, 1])
    dl = loss.sort_values("patient_delta_L_X0_minus_X1").reset_index(drop=True)
    vals = dl["patient_delta_L_X0_minus_X1"].to_numpy(float)
    cols = [CHELP if v > 0 else CHURT for v in vals]
    ax.hlines(range(len(vals)), 0, vals, color=cols, lw=1.6, zorder=2)
    ax.plot(vals, range(len(vals)), marker="s", ms=S.MS_SMALL + 0.4, linestyle="none",
            color=S.C_INK, markerfacecolor="white", markeredgewidth=0.9, zorder=3)

    S.refline(ax, 0.0, axis="x")
    mean_dl, med_dl = float(b4["mean_delta_L"]), float(b4["median_delta_L"])
    ax.axvline(mean_dl, color=S.C_INK, lw=1.1, ls="-", zorder=1)
    # Anchored to the RIGHT of the mean rule: the sign-flip inset sits to its left and the two
    # labels collide if either is centred.
    ax.annotate(f"mean {mean_dl:.3f}", xy=(mean_dl + 0.012, 15.9), fontsize=S.FS_ANNOT,
                color=S.C_INK, ha="left", va="bottom")
    # Two lines, anchored just right of the median rule. The single-line form is wider than the
    # 0.235-unit right-hand margin, so right-aligning it ran the text back across the mean rule
    # and over the P33/P13 bars. Every row is empty to the right of zero, so this height is safe.
    ax.annotate(f"median\n{med_dl:+.5f}", xy=(0.018, 1.6),
                fontsize=S.FS_ANNOT, color=S.C_MID, ha="left", va="center", linespacing=1.25)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels([f"P{int(p)}" for p in dl["patient_id"]], fontsize=S.FS_MIN)
    ax.set_xlim(-0.63, 0.235)
    ax.set_ylim(-0.8, len(vals) + 2.1)
    ax.set_xlabel("Δ loss (baseline − augmented), log-loss units", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "e")
    S.panel_title(ax, "Heterogeneous, with a negative mean")
    # The sign split is stated where the sign is *defined*, so neither count can be read as a
    # verdict. n_favor_augmented is the count of delta_L > 0; the rest favour the baseline.
    n_help = int(b4["n_favor_augmented"])
    n_hurt = int(b4["n_favor_baseline"])
    # One centred line: two anchored labels collide at this panel width.
    ax.annotate(f"◀ adding change hurt ({n_hurt})", xy=(0.485, -0.20), xycoords="axes fraction",
                fontsize=S.FS_ANNOT, color=CHURT, ha="right", va="top")
    ax.annotate(f"helped ({n_help}) ▶", xy=(0.515, -0.20), xycoords="axes fraction",
                fontsize=S.FS_ANNOT, color=CHELP, ha="left", va="top")

    # inset: exhaustive sign-flip reference distribution (top-left, clear of the lollipops)
    axi = ax.inset_axes([0.045, 0.755, 0.38, 0.20])
    null = signflip["null_mean_delta_L"].to_numpy(float)
    axi.hist(null, bins=90, color=S.C_LIGHT, edgecolor="none")
    axi.axvline(mean_dl, color=S.C_INK, lw=1.0)
    axi.set_yticks([])
    axi.set_xticks([-0.08, 0, 0.08])
    axi.tick_params(labelsize=S.FS_MIN, pad=1)
    for sp in ("top", "right", "left"):
        axi.spines[sp].set_visible(False)
    axi.spines["bottom"].set_color(S.C_MID)
    axi.annotate(f"$2^{{15}}$ sign flips\ntwo-sided $P$ = {float(b4['P_two_sided']):.3f}",
                 xy=(0.06, 0.97), xycoords="axes fraction", fontsize=S.FS_MIN,
                 color=S.C_MID, ha="left", va="top", linespacing=1.25)

    # =====================================================================================
    # f — influence concentration
    # =====================================================================================
    ax = fig.add_subplot(gs[1, 2])
    neg = dl[dl["patient_delta_L_X0_minus_X1"] < 0].copy()
    neg["mass"] = neg["patient_delta_L_X0_minus_X1"].abs()
    neg = neg.sort_values("mass", ascending=False).reset_index(drop=True)
    share = neg["mass"] / float(b4["total_negative_loss_mass"])
    cum = share.cumsum()

    xs = np.arange(len(neg))
    ax.bar(xs, share * 100, color=CHURT, width=0.62, zorder=2)
    ax.plot(xs, cum * 100, marker="o", ms=S.MS_SMALL, color=S.C_INK,
            lw=1.0, markerfacecolor="white", markeredgewidth=0.9, zorder=3)

    top2 = float(b4["top2_share_of_negative_mass"]) * 100
    ax.axhline(top2, color=S.C_INK, lw=S.LW_REF, ls=(0, (2.5, 2.0)), zorder=1)
    ax.annotate(f"{top2:.1f}% from 2 patients", xy=(len(neg) - 0.4, top2 + 2),
                fontsize=S.FS_ANNOT, color=S.C_INK, ha="right", va="bottom")

    ax.set_xticks(xs)
    ax.set_xticklabels([f"P{int(p)}" for p in neg["patient_id"]], fontsize=S.FS_MIN)
    ax.set_ylim(0, 108)
    ax.set_ylabel("Share of total negative loss mass (%)", fontsize=S.FS_AXIS_LABEL)
    S.subtle(ax)
    S.panel_label(ax, "f")
    S.panel_title(ax, "The mean is carried by two patients")
    # One annotation only. The leave-one-patient rerun answers the question this panel raises, so
    # it replaces the (weaker) largest-opposite-gain note. All 15 patients stay in every analysis;
    # the range is stated as a range precisely so it cannot be read as an exclusion result.
    ax.annotate(f"leave-one-patient rerun: mean stays negative, "
                f"{influence['mean_delta_L'].min():.3f} to {influence['mean_delta_L'].max():.3f}",
                xy=(0.5, -0.22), xycoords="axes fraction", fontsize=S.FS_ANNOT,
                color=S.C_MID, ha="center", va="top")

    # =====================================================================================
    # shared legend — response encoding is colour AND shape
    # =====================================================================================
    handles = [
        Line2D([], [], marker="o", ls="none", color=CR, markerfacecolor="white",
               markeredgewidth=1.0, label="Responder"),
        Line2D([], [], marker="D", ls="none", color=CNR, markerfacecolor="white",
               markeredgewidth=1.0, label="Non-responder"),
        Line2D([], [], marker="v", ls="none", color=S.C_INK, markerfacecolor="white",
               markeredgewidth=1.0, label="record-level estimate"),
        Line2D([], [], marker="s", ls="none", color=S.C_INK, markerfacecolor="white",
               markeredgewidth=1.0, label="patient-level estimate"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.005),
               ncol=4, fontsize=S.FS_LEGEND, frameon=False)

    paths = S.export(fig, OUTDIR, STEM)
    print(json.dumps({k: str(v.relative_to(ROOT).as_posix()) for k, v in paths.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
