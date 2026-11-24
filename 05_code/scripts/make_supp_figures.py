"""Supplementary figures S1-S10 for UroSig-Reeval (BMC Medicine, Additional file). Python/matplotlib only.

Usage: python make_supp_figures.py <project_root> <out_dir> [figure numbers...]
Same inputs, style and colour semantics as make_figures_v2.py (imported for the shared helpers and tables).
Every plotted table is written to <out_dir>/source_data/.
"""
import sys

import make_figures_v2 as M
from make_figures_v2 import (CLIN, COH, DASH, EV, GRADE, INK, MUTE, ORDER, PROLIF, RAND, RAND_L, RAND_LL, ROOT, SD, SIG, SIG_L, SIG_LL, SIG_M,
                             axmm, elig, fig_mm, half_violin, locked, minus, np, pairs, pd, plt, save, sig, stats, tag)

WHICH = set(sys.argv[3:]) or {str(i) for i in range(1, 11)}
ALLC = ORDER + ["GSE116918"]
NAME = dict(COH, GSE116918="GSE116918")


def strip(ax, s, **kw):
    """Short statement of the panel's numbers, set above the axes."""
    kw.setdefault("fontsize", 5.8); kw.setdefault("linespacing", 1.2); kw.setdefault("color", INK)
    return ax.text(0, 1.035, s, transform=ax.transAxes, va="bottom", ha="left", **kw)


def ranked(ax, d, c="harrell_c", lo="harrell_c_lo", hi="harrell_c_hi", ylim=(0.2, 0.95), grade=None):
    d = d.sort_values(c).reset_index(drop=True); x = np.arange(1, len(d) + 1); up = d[lo] > 0.5
    ax.vlines(x[~up], d[lo].clip(lower=ylim[0])[~up], d[hi].clip(upper=ylim[1])[~up], color=RAND_LL, lw=0.5, zorder=1)
    ax.vlines(x[up], d[lo][up], d[hi].clip(upper=ylim[1])[up], color=SIG_L, lw=0.5, zorder=1)
    ax.scatter(x[~up], d[c][~up], s=3.2, color=RAND, lw=0, zorder=3); ax.scatter(x[up], d[c][up], s=4.2, color=SIG, lw=0, zorder=4)
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=2)
    if grade is not None:
        ax.axhline(grade, color=GRADE, lw=0.9, zorder=2)
    ax.set_xlim(0.3, len(d) + 0.7); ax.set_ylim(*ylim)
    return d, int(up.sum())


def scatter_id(ax, x, y, lims, col=SIG_M, s=3.5):
    ax.plot(lims, lims, color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(x, y, s=s, color=col, alpha=0.5, lw=0, zorder=2)
    ax.set_xlim(*lims); ax.set_ylim(*lims)


# =============================================================== S1: gene coverage and evaluability
def s1():
    fig = fig_mm(72)
    e = elig[elig["分析角色"] == "main"].copy(); nmain = e["sig_id"].nunique()
    # ---- a
    ax = axmm(fig, 22, 12, 62, 46); tag(fig, 0, 0, "a")
    cats = [("Cohort used by the article", "#B9B9B9"), ("Gene coverage < 80%", "#EDEDED"), ("80% to < 100%", SIG_M), ("100%", SIG)]
    rows = []
    for c in ALLC:
        d = e[e["cohort"] == c]; used = int(d["原文用过"].sum()); full = int(d["可评估_1.0"].sum()); part = int(d["可评估_0.8"].sum()) - full
        rows.append({"cohort": NAME[c], cats[0][0]: used, cats[1][0]: len(d) - used - part - full, cats[2][0]: part, cats[3][0]: full})
    t = pd.DataFrame(rows).set_index("cohort"); t.to_csv(SD / "FigS1a_evaluability_by_cohort.csv")
    y = np.arange(len(t))[::-1]; left = np.zeros(len(t))
    for name, col in cats:
        ax.barh(y, t[name], left=left, height=0.62, color=col, ec="white", lw=0.5)
        for yy, v, l in zip(y, t[name], left):
            if v >= 14:
                ax.text(l + v / 2, yy, "%d" % v, ha="center", va="center", fontsize=5.4, color="white" if col in (SIG, SIG_M) else INK)
        left += t[name].to_numpy()
    ax.set_yticks(y); ax.set_yticklabels(t.index); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, nmain); ax.set_ylim(-0.5, len(t) - 0.5); ax.set_xlabel("Signatures with published coefficients (n = %d)" % nmain)
    k = axmm(fig, 22, 3.5, 62, 6); k.axis("off"); k.set_xlim(0, 62); k.set_ylim(0, 6)
    for (name, col), (x0, y0) in zip(cats, ((0, 4.3), (33, 4.3), (0, 1.3), (33, 1.3))):
        k.add_patch(plt.Rectangle((x0, y0 - 0.9), 1.8, 1.8, fc=col, ec=RAND, lw=0.3)); k.text(x0 + 2.7, y0, name, va="center", fontsize=5.6)
    # ---- b
    ax = axmm(fig, 104, 12, 26, 46); tag(fig, 91, 0, "b")
    kk = e[e["cohort"].isin(ORDER)].groupby("sig_id")["可评估_0.8"].sum().value_counts().reindex(range(5), fill_value=0)
    kk.rename("n_signatures").to_csv(SD / "FigS1b_cohorts_per_signature.csv")
    ax.bar(kk.index, kk.values, color=[RAND_L] + [SIG] * 4, width=0.7)
    for i, v in zip(kk.index, kk.values):
        ax.text(i, v + 1.5, "%d" % v, ha="center", va="bottom", fontsize=5.6)
    ax.set_xticks(range(5)); ax.set_xlabel("Independent prostatectomy\ncohorts available"); ax.set_ylabel("Signatures"); ax.set_ylim(0, kk.max() * 1.15)
    strip(ax, "%d of %d entered the\nprimary analysis" % (kk[1:].sum(), nmain))
    # ---- c
    ax = axmm(fig, 148, 12, 21, 46); tag(fig, 135, 0, "c")
    d = e[e["cohort"].isin(ORDER) & ~e["原文用过"].astype(bool) & e["缺失基因"].notna()]
    g = d.assign(g=d["缺失基因"].astype(str).str.split(";")).explode("g"); g = g[~g["g"].str.contains(r"\|", regex=True) & (g["g"] != "")]
    vc = g.groupby("g")["sig_id"].nunique().rename_axis("gene").reset_index(name="n").sort_values(["n", "gene"], ascending=[False, True]).set_index("gene")["n"]
    vc.rename("n_signatures").to_csv(SD / "FigS1c_missing_genes.csv")
    top = vc.head(12)[::-1]
    ax.barh(range(len(top)), top.values, color=RAND_L, height=0.68); ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontstyle="italic", fontsize=5.6)
    ax.tick_params(axis="y", length=0); ax.set_ylim(-0.6, len(top) - 0.4); ax.set_xlabel("Signatures in which\nthe gene was missing")
    ax.xaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
    save(fig, "FigureS1_gene_coverage")


# =============================================================== S2: every signature in every cohort
def s2():
    fig = fig_mm(112)
    pairs[["sig_id", "cohort", "n", "events", "harrell_c", "harrell_c_lo", "harrell_c_hi", "grade_c"]].to_csv(SD / "FigS2_c_by_cohort.csv", index=False)
    for i, c in enumerate(ORDER):
        x0, y0 = (15, 100)[i % 2], (9, 64)[i // 2]
        ax = axmm(fig, x0, y0, 69, 38); tag(fig, x0 - 15, y0 - 9, "abcd"[i])
        d = pairs[pairs.cohort == c]; g = d["grade_c"].median()
        d, nup = ranked(ax, d, grade=g)
        ax.text(2, g + 0.012, "Gleason grade group, C = %.2f" % g, color=GRADE, fontsize=5.8, ha="left", va="bottom")
        strip(ax, "%s: %d BCR events in %d patients\n%d signatures, 95%% CI above 0.5 in %d (%.0f%%)" % (COH[c], d["events"].max(), d["n"].max(), len(d), nup, 100 * nup / len(d)))
        ax.set_yticks([0.2, 0.35, 0.5, 0.65, 0.8, 0.95]); ax.set_ylabel("C-index (95% CI)" if i % 2 == 0 else "")
        ax.set_xlabel("Signatures ranked within the cohort" if i // 2 == 1 else "")
    save(fig, "FigureS2_every_signature_by_cohort")


# =============================================================== S3: agreement between discrimination measures
def s3():
    fig = fig_mm(50)
    p = pairs[["sig_id", "cohort", "harrell_c", "uno_c", "auc36", "auc60"]]; p.to_csv(SD / "FigS3_metrics_by_evaluation.csv", index=False)
    s = sig[["sig_id", "pooled_c", "median_uno_c", "median_auc36"]]; s.to_csv(SD / "FigS3_metrics_by_signature.csv", index=False)
    specs = [(p, "harrell_c", "uno_c", "Harrell C-index", "Uno C-index", "evaluations"), (p, "harrell_c", "auc60", "Harrell C-index", "AUC at 60 months", "evaluations"),
             (p, "auc36", "auc60", "AUC at 36 months", "AUC at 60 months", "evaluations"), (s, "pooled_c", "median_uno_c", "Pooled Harrell C-index", "Median Uno C-index", "signatures")]
    for i, (d, a, b, la, lb, unit) in enumerate(specs):
        x0 = 13 + i * 42.5; ax = axmm(fig, x0, 11, 27, 27); tag(fig, x0 - 13, 0, "abcd"[i])
        d = d.dropna(subset=[a, b]); scatter_id(ax, d[a], d[b], (0.22, 0.92) if unit == "evaluations" else (0.35, 0.8), s=3 if unit == "evaluations" else 5)
        ax.set_xlabel(la); ax.set_ylabel(lb); ax.set_xticks([0.4, 0.6, 0.8]); ax.set_yticks([0.4, 0.6, 0.8])
        strip(ax, minus("Spearman ρ = %.2f\nn = %d %s" % (stats.spearmanr(d[a], d[b]).statistic, len(d), unit)))
    save(fig, "FigureS3_agreement_between_measures")


# =============================================================== S4: hazard ratio per SD
def s4():
    fig = fig_mm(64)
    p = pairs[["sig_id", "cohort", "hr_sd", "hr_sd_lo", "hr_sd_hi", "hr_sd_p", "harrell_c"]].copy(); p.to_csv(SD / "FigS4_hazard_ratio_per_sd.csv", index=False)
    p["l"] = np.log2(p["hr_sd"])
    ax = axmm(fig, 14, 11, 84, 42); tag(fig, 0, 0, "a")
    rng = np.random.default_rng(3)
    for i, c in enumerate(ORDER):
        d = p[p.cohort == c]
        half_violin(ax, d["l"], i + 0.06, 0.36, SIG_LL)
        ax.scatter(i - 0.08 - rng.uniform(0, 0.24, len(d)), d["l"], s=1.6, color=SIG_M, alpha=0.7, lw=0, zorder=2)
        q = d["l"].quantile([.25, .5, .75]).to_numpy()
        ax.plot([i + 0.06] * 2, [q[0], q[2]], color=SIG, lw=1.6, solid_capstyle="butt", zorder=3); ax.scatter([i + 0.06], [q[1]], s=9, color="white", ec=SIG, lw=0.7, zorder=4)
        ax.text(i, np.log2(0.365), "median %.2f\n> 1 in %.0f%%, < 1 in %.0f%%" % (2 ** q[1], 100 * (d["hr_sd_lo"] > 1).mean(), 100 * (d["hr_sd_hi"] < 1).mean()), ha="center", va="bottom", fontsize=5.2, color=MUTE, linespacing=1.15)
    ax.axhline(0, color=INK, lw=0.5, ls=DASH, zorder=1)
    tk = [0.5, 0.75, 1, 1.5, 2, 3]; ax.set_yticks(np.log2(tk)); ax.set_yticklabels(["%g" % v for v in tk]); ax.set_ylim(np.log2(0.36), np.log2(3.3))
    ax.set_xticks(range(4)); ax.set_xticklabels([COH[c] for c in ORDER]); ax.set_xlim(-0.55, 3.6); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("Hazard ratio per SD of the score"); strip(ax, "Percentages: evaluations whose 95% CI lies entirely above or below 1")
    ax = axmm(fig, 118, 11, 42, 42); tag(fig, 104, 0, "b")
    ax.axhline(0, color=INK, lw=0.5, ls=DASH, zorder=1); ax.axvline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(p["harrell_c"], p["l"], s=3, color=SIG_M, alpha=0.5, lw=0, zorder=2)
    ax.set_yticks(np.log2(tk)); ax.set_yticklabels(["%g" % v for v in tk]); ax.set_ylim(np.log2(0.36), np.log2(3.3)); ax.set_xlim(0.26, 0.84)
    ax.set_xlabel("Harrell C-index"); ax.set_ylabel("Hazard ratio per SD of the score")
    strip(ax, minus("Spearman ρ = %.2f (n = %d evaluations)" % (stats.spearmanr(p["harrell_c"], p["l"]).statistic, len(p))))
    save(fig, "FigureS4_hazard_ratio_per_sd")


# =============================================================== S5: each alternative analysis against the primary analysis
def s5():
    fig = fig_mm(104)
    specs = [("S1_cov100", "Complete gene coverage only"), ("S2a_raw_log2", "Unstandardised log2 expression"), ("S2b_tcga_sd", "Scaled to TCGA standard deviations"),
             ("S3_plus_GSE116918", "Plus radiotherapy cohort"), ("S4_without_MSKCC", "Without MSKCC"), ("S5_GSE54460_only", "GSE54460 only"),
             ("S10_GSE54460_no_repeats", "GSE54460 without repeat libraries"), ("E3_met_GSE116918", "Metastasis (GSE116918)")]
    base = sig[["sig_id", "pooled_c"]].rename(columns={"pooled_c": "primary"}); out = []
    for i, (tg, name) in enumerate(specs):
        x0, y0 = 13 + (i % 4) * 42.5, (11, 61)[i // 4]
        ax = axmm(fig, x0, y0, 27, 27); tag(fig, x0 - 13, y0 - 11, "abcdefgh"[i])
        d = pd.read_csv(EV / tg / "signatures.csv", dtype={"sig_id": str})[["sig_id", "pooled_c"]].merge(base, on="sig_id"); out.append(d.assign(analysis=name))
        scatter_id(ax, d["primary"], d["pooled_c"], (0.33, 0.82), col=SIG, s=4)
        diff = d["pooled_c"] - d["primary"]
        strip(ax, minus("%s\nn = %d; change %+.3f; ρ = %.2f" % (name, len(d), diff.median(), stats.spearmanr(d["primary"], d["pooled_c"]).statistic)), fontsize=5.5)
        ax.set_xticks([0.4, 0.6, 0.8]); ax.set_yticks([0.4, 0.6, 0.8])
        ax.set_xlabel("Primary analysis" if i // 4 == 1 else ""); ax.set_ylabel("Alternative analysis" if i % 4 == 0 else "")
    pd.concat(out).to_csv(SD / "FigS5_alternative_vs_primary.csv", index=False)
    save(fig, "FigureS5_alternative_analyses")


# =============================================================== S6: signatures outside the primary analysis
def s6():
    fig = fig_mm(69)
    s7 = pd.read_csv(EV / "S7_plus_low_precision/signatures.csv", dtype={"sig_id": str}); s8 = pd.read_csv(EV / "S8_equal_weight/signatures.csv", dtype={"sig_id": str})
    lowp = s7[~s7["sig_id"].isin(sig["sig_id"])]
    groups = [("Published\ncoefficients", sig, SIG), ("Read from figures\nor hazard ratios", lowp, SIG), ("Equal\nweights", s8, SIG)]
    pd.concat([g.assign(group=n.replace("\n", " "))[["group", "sig_id", "pooled_c", "pooled_c_lo", "pooled_c_hi", "pooled_grade_c", "k_cohorts"]] for n, g, _ in groups]).to_csv(SD / "FigS6_signature_groups.csv", index=False)
    ax = axmm(fig, 14, 14, 50, 40); tag(fig, 0, 0, "a"); rng = np.random.default_rng(4)
    for i, (n, g, col) in enumerate(groups):
        half_violin(ax, g["pooled_c"], i + 0.06, 0.34, SIG_LL)
        ax.scatter(i - 0.08 - rng.uniform(0, 0.24, len(g)), g["pooled_c"], s=2.2, color=SIG_M, alpha=0.7, lw=0, zorder=2)
        q = g["pooled_c"].quantile([.25, .5, .75]).to_numpy()
        ax.plot([i + 0.06] * 2, [q[0], q[2]], color=col, lw=1.6, solid_capstyle="butt", zorder=3); ax.scatter([i + 0.06], [q[1]], s=9, color="white", ec=col, lw=0.7, zorder=4)
        ax.text(i, 0.335, "n = %d\nmedian %.3f" % (len(g), q[1]), ha="center", va="bottom", fontsize=5.2, color=MUTE, linespacing=1.15)
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=5.6, linespacing=1.05); ax.set_xlim(-0.55, 2.6); ax.set_ylim(0.33, 0.8); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("Pooled external C-index")
    for (x0, w, lab, g, name) in ((82, 55, "b", s8, "Equal-weight signatures (n = %d)\n" % len(s8)), (146, 22, "c", lowp, "From figures\nor HRs (n = %d)\n" % len(lowp))):
        ax = axmm(fig, x0, 14, w, 40); tag(fig, x0 - 13 if lab == "b" else x0 - 6, 0, lab)
        gm = g["pooled_grade_c"].median(); d, nup = ranked(ax, g, "pooled_c", "pooled_c_lo", "pooled_c_hi", ylim=(0.3, 0.9), grade=gm)
        strip(ax, "%s95%% CI above 0.5 in %d" % (name, nup))
        ax.set_yticks([0.3, 0.5, 0.7, 0.9]); ax.set_xlabel("Signatures ranked" if lab == "b" else "Ranked")
        if lab == "b":
            ax.set_ylabel("Pooled external C-index (95% CI)"); ax.text(1, gm + 0.012, "Gleason grade group, median C = %.2f" % gm, color=GRADE, fontsize=5.8, va="bottom")
        else:
            ax.tick_params(labelleft=False)
    save(fig, "FigureS6_signatures_outside_primary")


# =============================================================== S7: calibration slope by cohort
def s7():
    fig = fig_mm(88)
    v = pd.read_csv(EV / "posthoc/V5_calibration_slope_pairs.csv", dtype={"sig_id": str}); v.to_csv(SD / "FigS7_calibration_slope_by_cohort.csv", index=False)
    for i, c in enumerate(ORDER):
        x0, y0 = (13, 58)[i % 2], (11, 52)[i // 2]
        ax = axmm(fig, x0, y0, 34, 24); tag(fig, x0 - 13, y0 - 11, "abcd"[i])
        d = v[v.cohort == c]; med = d["slope"].median()
        ax.hist(d["slope"].clip(-1.5, 3.0), bins=np.arange(-1.5, 3.01, 0.25), color=SIG_L, ec="white", lw=0.5)
        top = ax.get_ylim()[1]
        ax.plot([0, 0], [0, top], color=INK, lw=0.5, ls=DASH); ax.plot([1, 1], [0, top], color=INK, lw=0.8); ax.plot([med, med], [0, top], color=SIG, lw=1.0)
        ax.set_xticks([-1, 0, 1, 2, 3]); ax.set_xticklabels([minus("-1"), "0", "1", "2", "≥ 3"]); ax.set_ylabel("Signatures" if i % 2 == 0 else "")
        ax.set_xlabel("Calibration slope" if i // 2 == 1 else "")
        strip(ax, "%s (n = %d): median %.2f\n95%% CI below 1 in %.0f%%" % (COH[c], len(d), med, 100 * (d["slope_hi"] < 1).mean()))
    ax = axmm(fig, 116, 11, 52, 65); tag(fig, 100, 0, "e")
    ax.axhline(1, color=INK, lw=0.8, zorder=1); ax.axhline(0, color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(v["lp_sd"], v["slope"].clip(-1.5, 3.0), s=3.5, color=SIG_M, alpha=0.5, lw=0, zorder=2)
    ax.set_xscale("log"); ax.set_xlim(0.008, 25); ax.set_ylim(-1.6, 3.1); ax.set_yticks([-1, 0, 1, 2, 3]); ax.set_yticklabels([minus("-1"), "0", "1", "2", "≥ 3"])
    ax.set_xticks([0.01, 0.1, 1, 10]); ax.set_xticklabels(["0.01", "0.1", "1", "10"])
    ax.set_xlabel("SD of the published linear predictor\nin the evaluation cohort"); ax.set_ylabel("Calibration slope")
    ax.text(24, 1.06, "Ideal", fontsize=5.8, ha="right", va="bottom")
    strip(ax, minus("Spearman ρ = %.2f (n = %d evaluations)" % (stats.spearmanr(v["lp_sd"], v["slope"]).statistic, len(v))))
    save(fig, "FigureS7_calibration_by_cohort")


# =============================================================== S8: radiotherapy cohort, metastasis and biochemical recurrence
def s8():
    fig = fig_mm(112)
    a = pd.read_csv(EV / "run_met_116918/pairs.csv", dtype={"sig_id": str}); b = pd.read_csv(EV / "run_116918/pairs.csv", dtype={"sig_id": str})
    cols = ["sig_id", "n", "events", "harrell_c", "harrell_c_lo", "harrell_c_hi", "grade_c", "rand_pct"]
    a[cols].to_csv(SD / "FigS8_metastasis_GSE116918.csv", index=False); b[cols].to_csv(SD / "FigS8_bcr_GSE116918.csv", index=False)
    for i, (d, name) in enumerate(((a, "Metastasis"), (b, "Biochemical recurrence"))):
        y0 = (9, 64)[i]
        ax = axmm(fig, 15, y0, 92, 38); tag(fig, 0, y0 - 9, "ac"[i])
        g = d["grade_c"].iloc[0]; dd, nup = ranked(ax, d, grade=g)
        ax.text(2, g + 0.012, "Gleason grade group, C = %.2f" % g, color=GRADE, fontsize=5.8, ha="left", va="bottom")
        strip(ax, "%s: %d events in %d patients\n%d signatures, median C = %.3f, 95%% CI above 0.5 in %d (%.0f%%)" % (name, d["events"].max(), d["n"].max(), len(dd), dd["harrell_c"].median(), nup, 100 * nup / len(dd)))
        ax.set_yticks([0.2, 0.35, 0.5, 0.65, 0.8, 0.95]); ax.set_ylabel("C-index (95% CI)"); ax.set_xlabel("Signatures ranked" if i == 1 else "")
        ax = axmm(fig, 125, y0, 42, 38); tag(fig, 112, y0 - 9, "bd"[i])
        r = d["rand_pct"].dropna(); nn, _, patches = ax.hist(r, bins=np.arange(0, 101, 5), color=RAND_L, ec="white", lw=0.5); patches[-1].set_facecolor(SIG)
        ax.axhline(len(r) / 20, color=INK, lw=0.5, ls=DASH); ax.set_xlim(0, 100); ax.set_ylabel("Signatures")
        ax.set_xlabel("Percentile among 1000 random gene sets" if i == 1 else "")
        strip(ax, "%.0f%% above the 95th percentile" % (100 * (r > 95).mean()), color=SIG)
    save(fig, "FigureS8_radiotherapy_cohort")


# =============================================================== S9: heterogeneity and precision
def s9():
    fig = fig_mm(54)
    s = sig[["sig_id", "k_cohorts", "pooled_c", "pooled_c_lo", "pooled_c_hi", "i2", "tau2"]].copy(); s["half"] = (s["pooled_c_hi"] - s["pooled_c_lo"]) / 2
    s.to_csv(SD / "FigS9_heterogeneity_and_precision.csv", index=False)
    ax = axmm(fig, 13, 11, 40, 32); tag(fig, 0, 0, "a")
    h = s.dropna(subset=["i2"]); ax.hist(h["i2"] * 100, bins=np.arange(0, 101, 10), color=SIG_L, ec="white", lw=0.5)
    ax.set_xlabel("I² (%)"); ax.set_ylabel("Signatures"); ax.set_xlim(0, 100)
    strip(ax, "%d signatures with two or more cohorts\nI² = 0 in %.0f%%, I² > 50%% in %.0f%%" % (len(h), 100 * (h["i2"] == 0).mean(), 100 * (h["i2"] > 0.5).mean()))
    ax = axmm(fig, 72, 11, 40, 32); tag(fig, 59, 0, "b"); rng = np.random.default_rng(5)
    for k_ in (1, 2, 3, 4):
        d = s[s["k_cohorts"] == k_]
        ax.scatter(k_ + rng.uniform(-0.22, 0.22, len(d)), d["half"].clip(upper=0.5), s=3.5, color=SIG_M, alpha=0.6, lw=0, zorder=2)
        ax.plot([k_ - 0.32, k_ + 0.32], [d["half"].clip(upper=0.5).median()] * 2, color=SIG, lw=1.4, zorder=3)
        ax.text(k_, 0.535, "n = %d" % len(d), ha="center", va="center", fontsize=5.2, color=MUTE)
    ax.set_xticks([1, 2, 3, 4]); ax.set_xlim(0.4, 4.6); ax.set_ylim(0, 0.56); ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5]); ax.set_yticklabels(["0", "0.1", "0.2", "0.3", "0.4", "≥ 0.5"])
    ax.set_xlabel("Cohorts contributing to the pooled estimate"); ax.set_ylabel("Half-width of the 95% CI")
    strip(ax, "Lines show medians; two cohorts give one degree\nof freedom for the Hartung–Knapp interval")
    ax = axmm(fig, 131, 11, 38, 32); tag(fig, 118, 0, "c")
    ax.axvline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    for k_, col in zip((1, 2, 3, 4), (RAND_L, SIG_L, SIG_M, SIG)):
        d = s[s["k_cohorts"] == k_]; ax.scatter(d["pooled_c"], d["half"].clip(upper=0.5), s=5, color=col, lw=0, zorder=2 + k_)
        ax.text((0, 0.5)[(k_ - 1) % 2], 1.035 + 0.078 * (1 - (k_ - 1) // 2), "● %d cohort%s" % (k_, "" if k_ == 1 else "s"), transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8, color=col if k_ > 1 else RAND)
    ax.set_xlim(0.36, 0.74); ax.set_ylim(0, 0.56); ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5]); ax.set_yticklabels(["0", "0.1", "0.2", "0.3", "0.4", "≥ 0.5"])
    ax.set_xlabel("Pooled external C-index"); ax.set_ylabel("Half-width of the 95% CI")
    save(fig, "FigureS9_heterogeneity_and_precision")


# =============================================================== S10: route by which the article was identified
def s10():
    fig = fig_mm(58)
    s = sig.merge(locked[["sig_id", "来源版本"]], on="sig_id", how="left"); s["route"] = np.where(s["来源版本"].str.startswith("v1"), "Keyword-prioritised\nrecords", "Dual-screened\nrecords")
    routes = ["Keyword-prioritised\nrecords", "Dual-screened\nrecords"]
    s[["sig_id", "route", "pooled_c", "ci_excludes_0.5", "beats_random_95", "beats_grade", "lr_q_bh_combined"]].assign(route=s["route"].str.replace("\n", " ")).to_csv(SD / "FigS10_by_identification_route.csv", index=False)
    ax = axmm(fig, 14, 11, 46, 36); tag(fig, 0, 0, "a"); rng = np.random.default_rng(6)
    for i, r in enumerate(routes):
        g = s[s["route"] == r]
        half_violin(ax, g["pooled_c"], i + 0.06, 0.34, SIG_LL)
        ax.scatter(i - 0.08 - rng.uniform(0, 0.24, len(g)), g["pooled_c"], s=2.2, color=SIG_M, alpha=0.7, lw=0, zorder=2)
        q = g["pooled_c"].quantile([.25, .5, .75]).to_numpy()
        ax.plot([i + 0.06] * 2, [q[0], q[2]], color=SIG, lw=1.6, solid_capstyle="butt", zorder=3); ax.scatter([i + 0.06], [q[1]], s=9, color="white", ec=SIG, lw=0.7, zorder=4)
        ax.text(i, 0.335, "n = %d\nmedian %.3f" % (len(g), q[1]), ha="center", va="bottom", fontsize=5.2, color=MUTE, linespacing=1.15)
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    mw = stats.mannwhitneyu(s.loc[s["route"] == routes[0], "pooled_c"], s.loc[s["route"] == routes[1], "pooled_c"])
    ax.set_xticks(range(2)); ax.set_xticklabels(routes, linespacing=1.05); ax.set_xlim(-0.55, 1.6); ax.set_ylim(0.33, 0.8); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("Pooled external C-index"); strip(ax, "Mann–Whitney P = %.2f" % mw.pvalue)
    ax = axmm(fig, 104, 11, 64, 36); tag(fig, 68, 0, "b")
    items = [("95% CI of pooled C above 0.5", s["ci_excludes_0.5"].astype(bool) & (s["pooled_c"] > 0.5)), ("Better than expression-matched\nrandom gene sets", s["beats_random_95"].astype(bool)),
             ("Adds to the clinical model\n(FDR < 0.05)", s["lr_q_bh_combined"] < 0.05), ("Better than Gleason\ngrade group", s["beats_grade"].astype(bool))]
    y = np.arange(len(items))[::-1]
    for j, (r, col) in enumerate(zip(routes, (SIG, SIG_L))):
        m = s["route"] == r; vals = [100 * v[m].mean() for _, v in items]
        ax.barh(y + (0.19 if j == 0 else -0.19), vals, height=0.36, color=col)
        for yy, v in zip(y, vals):
            ax.text(v + 1.2, yy + (0.19 if j == 0 else -0.19), "%.0f" % v, va="center", fontsize=5.4)
        ax.text(0.99, 0.30 - 0.09 * j, r.replace("\n", " "), transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color=col if j == 0 else SIG_M)
    ax.set_yticks(y); ax.set_yticklabels([n for n, _ in items], linespacing=1.05); ax.tick_params(axis="y", length=0); ax.set_xlim(0, 80); ax.set_xlabel("Signatures (%)")
    save(fig, "FigureS10_by_identification_route")


if __name__ == "__main__":
    for k, f in (("1", s1), ("2", s2), ("3", s3), ("4", s4), ("5", s5), ("6", s6), ("7", s7), ("8", s8), ("9", s9), ("10", s10)):
        if k in WHICH:
            f(); print("figure S" + k, "done")
