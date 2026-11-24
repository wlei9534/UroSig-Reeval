"""Main figures 1-6 for UroSig-Reeval (BMC Medicine), second design. Python/matplotlib only.

Usage: python make_figures_v2.py <project_root> <out_dir> [figure numbers...]
Inputs: 06_results/evaluation_v2 (incl. figure_inputs/ from figure_inputs.py), 06_results/eligibility, 03_extraction/locked.
BMC Medicine: full width 170 mm, height <= 225 mm incl. legend, keys inside the graphic, no title/legend in the file,
fonts embedded, lines > 0.25 pt, <= 10 MB per file. Every plotted table is written to <out_dir>/source_data/.
Colour semantics are fixed across all figures:
  blue = published signatures; vermilion = Gleason grade group; grey = random gene sets; ochre = proliferation index;
  charcoal = clinical model; purple = values reported by the original articles.
"""
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Patch
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

ROOT, OUT = Path(sys.argv[1]), Path(sys.argv[2])
WHICH = set(sys.argv[3:]) or set("123456")
EV = ROOT / "06_results" / "evaluation_v2"
FI = EV / "figure_inputs"
SD = OUT / "source_data"
SD.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6.5,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.right": False, "axes.spines.top": False,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 2.5, "legend.frameon": False, "lines.linewidth": 0.8,
    "patch.linewidth": 0.4, "hatch.linewidth": 0.3, "figure.dpi": 150, "savefig.facecolor": "white",
})
MM = 1 / 25.4
W = 170.0  # mm
SIG, SIG_M, SIG_L, SIG_LL = "#24588F", "#6D94BF", "#B5CAE1", "#E6EDF5"
GRADE, GRADE_L = "#C4472B", "#EBC2B6"
RAND, RAND_L, RAND_LL = "#8E8E8E", "#CFCFCF", "#EDEDED"
PROLIF, CLIN, REPORT, REPORT_L, TEAL = "#C38E21", "#2B2B2B", "#7B4F93", "#D9CBE2", "#2C8C84"
INK, MUTE = "#1E1E1E", "#6B6B6B"
DASH = (0, (3.5, 2))
C_CMAP = LinearSegmentedColormap.from_list("c", ["#8A5A22", "#D9B98C", "#FBFBFB", "#8FB0D3", "#173E69"])
R_CMAP = LinearSegmentedColormap.from_list("r", ["#8A5A22", "#E6D2B5", "#FFFFFF", "#A9C2DE", "#173E69"])
COH = {"GSE54460": "GSE54460", "GSE70769": "GSE70769", "GSE70768": "GSE70768", "MSKCC_GSE21032": "MSKCC"}
ORDER = ["GSE54460", "GSE70769", "GSE70768", "MSKCC_GSE21032"]


def fig_mm(h):
    return plt.figure(figsize=(W * MM, h * MM))


def axmm(fig, x, y, w, h, **kw):
    """Axes placed in mm; x from the left edge, y from the TOP edge of the figure."""
    H = fig.get_figheight() / MM
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def tag(fig, x, y, s):
    H = fig.get_figheight() / MM
    fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight="bold", va="top", ha="left", color=INK)


def note(ax, x, y, s, **kw):
    kw.setdefault("fontsize", 6); kw.setdefault("color", INK); kw.setdefault("transform", ax.transAxes)
    return ax.text(x, y, s, **kw)


def save(fig, name):
    for ext, kw in (("pdf", {}), ("svg", {}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}), ("png", {"dpi": 300})):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", pad_inches=0.01, **kw)
    plt.close(fig)


def minus(s):
    return s.replace("-", "−")


def half_violin(ax, v, pos, width, color, side=1, lw=0):
    v = np.asarray(v, float)
    k = stats.gaussian_kde(v)
    yy = np.linspace(v.min(), v.max(), 200)
    d = k(yy); d = d / d.max() * width
    ax.fill_betweenx(yy, pos, pos + side * d, color=color, lw=lw, zorder=1)


sig = pd.read_csv(EV / "primary/signatures.csv", dtype={"sig_id": str})
pairs = pd.read_csv(EV / "run_primary/pairs.csv", dtype={"sig_id": str})
e12 = pd.read_csv(EV / "exploratory/signatures_E1E2.csv", dtype={"sig_id": str})
feat = pd.read_csv(EV / "posthoc/V3_features.csv", dtype={"sig_id": str})
locked = pd.read_csv(ROOT / "03_extraction/locked/signatures_locked_v2_2026-11-14.csv", dtype=str).fillna("")
terms = pd.read_csv(ROOT / "03_extraction/locked/signature_terms_locked_v2_2026-11-14.csv", dtype=str).fillna("")
endp = pd.read_csv(ROOT / "06_results/eligibility/original_endpoint_class_v2.csv", dtype=str)
leak = pd.read_csv(ROOT / "06_results/eligibility/leakage_matrix_v2.csv", dtype=str).fillna("")
elig = pd.read_csv(ROOT / "06_results/eligibility/eligibility_signature_by_cohort_v2.csv", dtype={"sig_id": str})
meta = sig[["sig_id", "pooled_c", "k_cohorts"]].merge(feat[["sig_id", "tcga_only_training", "ext_reported", "log2_genes", "method"]], on="sig_id", how="left") \
    .merge(endp[["sig_id", "终点类别"]].rename(columns={"终点类别": "endpoint"}), on="sig_id", how="left") \
    .merge(locked[["sig_id", "年份", "基因数"]].rename(columns={"年份": "year", "基因数": "n_genes"}), on="sig_id", how="left")
meta["year"] = meta["year"].astype(int); meta["n_genes"] = meta["n_genes"].astype(int)
cm = pairs.groupby("cohort")[["n", "events"]].max()
NPT, NEV = int(cm["n"].sum()), int(cm["events"].sum())


# =============================================================== Figure 1: design and evidence base
def fig1():
    fig = fig_mm(121)
    # ---- a: flow
    ax = axmm(fig, 0.5, 0.5, 169, 48.5); ax.set_xlim(-0.6, 170.6); ax.set_ylim(0, 48.5); ax.axis("off")
    tag(fig, 0, 0, "a")

    def box(x, yc, w, h, head, n, fc="white", ec=INK, sub=None):
        ax.add_patch(FancyBboxPatch((x, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.8", fc=fc, ec=ec, lw=0.5, zorder=2))
        if sub:
            ax.text(x + w / 2, yc + h / 2 - 1.6, head, ha="center", va="top", fontsize=6, linespacing=1.15, zorder=3)
            ax.text(x + w / 2, yc - h / 2 + 1.5, n, ha="center", va="bottom", fontsize=6.5, fontweight="bold", zorder=3)
        else:
            ax.text(x + w / 2, yc + 0.2, head + "\n" + n, ha="center", va="center", fontsize=6, linespacing=1.25, zorder=3)

    def line(xs, ys, arrow=True, color=INK):
        ax.plot(xs, ys, color=color, lw=0.5, solid_capstyle="butt", zorder=1)
        if arrow:
            ax.annotate("", xy=(xs[-1], ys[-1]), xytext=(xs[-2], ys[-2]), arrowprops=dict(arrowstyle="-|>", lw=0.5, color=color, mutation_scale=5, shrinkA=0, shrinkB=0), zorder=1)

    up, lo, mid, h = 41, 17, 29, 12
    box(0, mid, 24, h, "Records identified\nin PubMed", "n = 1011", fc=SIG_LL, sub=True)
    box(34, up, 38, h, "Keyword-prioritised records,\nfull text assessed", "n = 292", sub=True)
    box(34, lo, 38, h, "Title and abstract screened\nindependently in duplicate", "n = 719  (κ = 0.97)", sub=True)
    box(82, lo, 34, h, "Full text assessed,\ndual independent extraction", "n = 300", sub=True)
    box(124, mid, 20, h + 3, "Analysable\nsignatures", "n = 375", fc=SIG_LL, sub=True)
    box(150, mid, 20, h + 3, "Primary\nanalysis", "211 signatures", fc=SIG_L, sub=True)
    line([24, 29, 29, 34], [mid, mid, up, up]); line([24, 29, 29, 34], [mid, mid, lo, lo])
    line([72, 82], [lo, lo]); line([72, 134, 134], [up, up, mid + (h + 3) / 2]); line([116, 134, 134], [lo, lo, mid - (h + 3) / 2])
    line([144, 150], [mid, mid])
    g = dict(fontsize=5.6, color=MUTE, linespacing=1.25)
    ax.text(77, lo - 7.5, "Excluded, 396\nNot retrievable, 19\nDuplicate preprint or\nnon-English, 4", ha="center", va="top", **g)
    line([77, 77], [lo, lo - 7], arrow=False, color=MUTE)
    ax.text(103, up + 1.2, "Excluded, 70", ha="center", va="bottom", **g)
    ax.text(125, lo - 7.5, "Excluded, 148: not a new signature, 58; non-linear model, 29;\nnot prognostic, 21; not reproducible, 14; other, 26", ha="center", va="top", **g)
    line([125, 125], [lo, lo - 7], arrow=False, color=MUTE)
    ax.text(169.5, up + 1.2, "Weights not published as coefficients, 128 (sensitivity analyses)\nNo independent cohort with ≥ 80% gene coverage, 36", ha="right", va="bottom", **g)
    line([147, 147], [mid, up + 1], arrow=False, color=MUTE)
    ax.text(169.6, mid - (h + 3) / 2 - 1.5, "569 evaluations in 4 radical\nprostatectomy cohorts\n(%d patients, %d BCR events)" % (NPT, NEV), ha="right", va="top", fontsize=5.6, color=INK, linespacing=1.25)

    # ---- b: signatures per year
    axb = axmm(fig, 9, 61, 46, 50); tag(fig, 0, 55, "b")
    an = locked[locked["分析角色"] != "excluded"].copy(); an["year"] = an["年份"].astype(int)
    an["primary"] = an["sig_id"].isin(sig["sig_id"])
    yr = an.groupby(["year", "primary"]).size().unstack(fill_value=0).reindex(range(2010, 2027), fill_value=0)
    yr.to_csv(SD / "Fig1b_signatures_per_year.csv")
    axb.bar(yr.index, yr[True], color=SIG, width=0.78, label="In primary analysis (n = %d)" % yr[True].sum())
    axb.bar(yr.index, yr[False], bottom=yr[True], color=SIG_L, width=0.78, label="Other analysable (n = %d)" % yr[False].sum())
    axb.set_xlim(2009.3, 2026.7); axb.set_xticks([2010, 2014, 2018, 2022, 2026]); axb.set_ylabel("Signatures published"); axb.set_xlabel("Year of publication")
    axb.legend(loc="upper left", handlelength=0.9, handleheight=0.8, borderaxespad=0.1, labelspacing=0.3)

    # ---- c: characteristics of the 211 signatures
    axc = axmm(fig, 84, 61, 44, 50); tag(fig, 62, 55, "c")
    g_bins = pd.cut(meta["n_genes"], [1, 5, 10, 20, 1e4], labels=["2–5", "6–10", "11–20", "> 20"])
    rows = [("Training data", [("TCGA only", int(meta["tcga_only_training"].sum())), ("Other", int((~meta["tcga_only_training"].astype(bool)).sum()))]),
            ("Final model", [("Multivariable Cox", int((meta["method"] == "多因素 Cox").sum())), ("LASSO-type", int((meta["method"] == "LASSO 类").sum())), ("", int((meta["method"] == "单因素或其他").sum()))]),
            ("Original endpoint", [("BCR or RFS", int((meta["endpoint"] == "BCR/RFS").sum())), ("Other", int((meta["endpoint"] != "BCR/RFS").sum()))]),
            ("External validation\nby the authors", [("Reported", int(meta["ext_reported"].sum())), ("Not reported", int((~meta["ext_reported"].astype(bool)).sum()))]),
            ("Genes in the model", [(k, int(v)) for k, v in g_bins.value_counts().sort_index().items()]),
            ("Independent cohorts\navailable here", [({4: "Four", 3: "Three", 2: "Two", 1: "One"}[int(k)], int(v)) for k, v in meta["k_cohorts"].value_counts().sort_index(ascending=False).items()])]
    pd.DataFrame([(r, k, v) for r, items in rows for k, v in items], columns=["feature", "level", "n"]).to_csv(SD / "Fig1c_characteristics.csv", index=False)
    shades = [SIG, SIG_M, SIG_L, SIG_LL]
    for i, (name, items) in enumerate(rows):
        y, left, tot = len(rows) - 1 - i, 0, sum(v for _, v in items)
        for j, (k, v) in enumerate(items):
            axc.barh(y, v / tot * 100, left=left, height=0.62, color=shades[j], ec="white", lw=0.5)
            if v / tot > 0.13 and k:
                axc.text(left + v / tot * 50, y, "%s\n%d" % (k, v), ha="center", va="center", fontsize=5.3, color="white" if j < 2 else INK, linespacing=1.05)
            left += v / tot * 100
    axc.set_yticks(range(len(rows))); axc.set_yticklabels([r[0] for r in rows][::-1], linespacing=1.05)
    axc.set_xlim(0, 100); axc.set_ylim(-0.5, len(rows) - 0.5); axc.set_xlabel("Signatures in the primary analysis (%)")
    axc.tick_params(axis="y", length=0); axc.spines["left"].set_visible(False)

    # ---- d: prior use of the cohorts
    axd = axmm(fig, 147, 61, 20, 50); tag(fig, 131, 55, "d")
    lk = leak[leak["sig_id"].isin(sig["sig_id"])]
    use = pd.Series({k: (lk[c] == "用过").mean() * 100 for k, c in (("TCGA-PRAD", "TCGA"), ("GSE70769", "GSE70769"), ("MSKCC", "MSKCC"), ("GSE70768", "GSE70768"), ("GSE116918", "GSE116918"), ("GSE54460", "GSE54460"))})
    use.to_csv(SD / "Fig1d_prior_use_of_cohorts.csv")
    y = np.arange(len(use))[::-1]
    axd.barh(y, use.values, height=0.62, color=[GRADE if k == "TCGA-PRAD" else RAND_L for k in use.index])
    for yy, v in zip(y, use.values):
        axd.text(v + 2 if v < 70 else v - 2, yy, "%.0f" % v, ha="left" if v < 70 else "right", va="center", fontsize=5.6, color=INK if v < 70 else "white")
    axd.set_yticks(y); axd.set_yticklabels(use.index); axd.set_xlim(0, 100); axd.set_xticks([0, 50, 100]); axd.set_ylim(-0.5, len(use) - 0.5)
    axd.set_xlabel("Signatures whose article\nhad used the cohort (%)"); axd.tick_params(axis="y", length=0)
    save(fig, "Figure1_design_and_evidence_base")


# =============================================================== Figure 2: landscape
def fig2():
    s = sig.merge(meta.drop(columns=["pooled_c", "k_cohorts"]), on="sig_id").sort_values("pooled_c").reset_index(drop=True)
    n = len(s); x = np.arange(1, n + 1)
    gmed = s["pooled_grade_c"].median()
    fig = fig_mm(160)
    L, WD = 28, 129
    # ---- a1: ranked pooled C
    ax = axmm(fig, L, 3, WD, 40); tag(fig, 0, 0, "a")
    up = s["pooled_c_lo"] > 0.5
    ax.vlines(x[~up], s["pooled_c_lo"].clip(lower=0.3)[~up], s["pooled_c_hi"].clip(upper=0.9)[~up], color=RAND_LL, lw=0.5, zorder=1)
    ax.vlines(x[up], s["pooled_c_lo"][up], s["pooled_c_hi"].clip(upper=0.9)[up], color=SIG_L, lw=0.5, zorder=1)
    ax.scatter(x[~up], s["pooled_c"][~up], s=3.2, color=RAND, lw=0, zorder=3)
    ax.scatter(x[up], s["pooled_c"][up], s=4.2, color=SIG, lw=0, zorder=4)
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=2); ax.axhline(gmed, color=GRADE, lw=0.9, zorder=2)
    ax.text(2, gmed + 0.012, "Gleason grade group in the same patients, median C = %.2f" % gmed, color=GRADE, fontsize=6, va="bottom")
    ax.text(2, 0.875, "95%% CI above 0.5 (n = %d)" % up.sum(), color=SIG, fontsize=6, va="top")
    ax.text(2, 0.835, "95%% CI includes 0.5 (n = %d)" % (~up).sum(), color=MUTE, fontsize=6, va="top")
    q1, med, q3 = s["pooled_c"].quantile([.25, .5, .75])
    ax.text(n - 1, 0.33, "Median %.3f (IQR %.3f–%.3f)" % (med, q1, q3), color=INK, fontsize=6, va="bottom", ha="right")
    ax.set_xlim(0.5, n + 0.5); ax.set_ylim(0.3, 0.9); ax.set_yticks([0.3, 0.5, 0.7, 0.9]); ax.set_ylabel("Pooled external\nC-index (95% CI)")
    ax.tick_params(axis="x", bottom=False, labelbottom=False); ax.spines["bottom"].set_visible(False)
    # ---- a2: cohort heatmap
    M = pairs.pivot(index="cohort", columns="sig_id", values="harrell_c").reindex(index=ORDER, columns=s["sig_id"])
    E = elig[elig["cohort"].isin(ORDER)].pivot(index="cohort", columns="sig_id", values="原文用过").reindex(index=ORDER, columns=s["sig_id"])
    M.to_csv(SD / "Fig2a_c_by_cohort.csv"); s.to_csv(SD / "Fig2a_signatures_ranked.csv", index=False)
    axh = axmm(fig, L, 45, WD, 13)
    bg = np.where(M.isna().to_numpy(), np.where(E.to_numpy(dtype=bool), 1.0, 0.0), np.nan)
    axh.imshow(bg, aspect="auto", cmap=ListedColormap(["#F4F4F4", "#B9B9B9"]), vmin=0, vmax=1, extent=(0.5, n + 0.5, 4, 0), interpolation="none")
    im = axh.imshow(np.ma.masked_invalid(M.to_numpy()), aspect="auto", cmap=C_CMAP, norm=TwoSlopeNorm(0.5, 0.3, 0.8), extent=(0.5, n + 0.5, 4, 0), interpolation="none")
    axh.set_yticks(np.arange(4) + 0.5); axh.set_yticklabels([COH[c] for c in ORDER]); axh.hlines([1, 2, 3], 0.5, n + 0.5, color="white", lw=0.6)
    axh.tick_params(axis="both", length=0, labelbottom=False)
    for sp in axh.spines.values():
        sp.set_visible(False)
    cax = axmm(fig, L + WD + 2.5, 45, 1.8, 13); cb = fig.colorbar(im, cax=cax, ticks=[0.3, 0.5, 0.8]); cb.outline.set_linewidth(0.4); cb.ax.tick_params(length=1.5, labelsize=5.5)
    cax.text(0, 1.16, "C-index", transform=cax.transAxes, fontsize=5.5, ha="left", va="bottom")
    # ---- a3: annotation tracks
    axt = axmm(fig, L, 60, WD, 15.5)
    yr = s["year"].clip(lower=2018).to_numpy(); ng = np.log2(s["n_genes"].clip(upper=64)).to_numpy()
    tracks = [("Trained on TCGA only", s["tcga_only_training"].astype(float).to_numpy(), ListedColormap(["#F1F1F1", "#3A3A3A"]), 0, 1),
              ("External validation reported", s["ext_reported"].astype(float).to_numpy(), ListedColormap(["#F1F1F1", "#3A3A3A"]), 0, 1),
              ("Developed for BCR or RFS", (s["endpoint"] == "BCR/RFS").astype(float).to_numpy(), ListedColormap(["#F1F1F1", "#3A3A3A"]), 0, 1),
              ("Number of genes", ng, LinearSegmentedColormap.from_list("g", ["#F3EEF6", REPORT]), 1, 6),
              ("Year of publication", yr, LinearSegmentedColormap.from_list("y", ["#EEF4F3", TEAL]), 2018, 2026)]
    for i, (name, v, cmap, lo, hi) in enumerate(tracks):
        axt.imshow(v[None, :], aspect="auto", cmap=cmap, vmin=lo, vmax=hi, extent=(0.5, n + 0.5, i + 0.86, i), interpolation="none")
    axt.set_ylim(len(tracks), 0); axt.set_xlim(0.5, n + 0.5); axt.set_yticks(np.arange(len(tracks)) + 0.43); axt.set_yticklabels([t[0] for t in tracks])
    axt.tick_params(axis="y", length=0); axt.set_xlabel("Signatures ranked by pooled external C-index (n = %d)" % n)
    for sp in ("left",):
        axt.spines[sp].set_visible(False)
    axt.set_xticks([1, 50, 100, 150, 200])
    # keys for the heatmap and tracks: one strip under the axis title
    axk = axmm(fig, L, 84, WD, 3.2); axk.axis("off")
    g = np.linspace(0, 1, 64)[None, :]
    axk.imshow(g, aspect="auto", cmap=tracks[3][2], extent=(88, 96, 0.7, 2.5), interpolation="none")
    axk.imshow(g, aspect="auto", cmap=tracks[4][2], extent=(121, 129, 0.7, 2.5), interpolation="none")
    axk.set_xlim(0, WD); axk.set_ylim(0, 3.2)
    for x0, col, lab in ((0, "#3A3A3A", "Yes"), (9, "#F1F1F1", "No"), (20, "#B9B9B9", "Cohort used by the article"), (52, "#F4F4F4", "Gene coverage < 80%")):
        axk.add_patch(plt.Rectangle((x0, 0.7), 1.8, 1.8, fc=col, ec=RAND, lw=0.3)); axk.text(x0 + 2.7, 1.55, lab, va="center", fontsize=5.4)
    for x0, x1, a_, b_ in ((87, 97, "Genes  2", "≥ 64"), (120, 130, "Year  ≤ 2018", "2026")):
        axk.text(x0, 1.55, a_, va="center", ha="right", fontsize=5.4); axk.text(x1, 1.55, b_, va="center", ha="left", fontsize=5.4)

    # ---- b: per-cohort distributions
    axb = axmm(fig, 12, 97, 68, 52); tag(fig, 0, 92, "b")
    rng = np.random.default_rng(1)
    pc = pairs[["sig_id", "cohort", "harrell_c", "grade_c", "c_clin", "n", "events"]]; pc.to_csv(SD / "Fig2b_per_cohort_c.csv", index=False)
    for i, c in enumerate(ORDER):
        d = pc[pc.cohort == c]
        half_violin(axb, d["harrell_c"], i + 0.06, 0.36, SIG_LL)
        axb.scatter(i - 0.08 - rng.uniform(0, 0.24, len(d)), d["harrell_c"], s=1.6, color=SIG_M, alpha=0.7, lw=0, zorder=2)
        q = d["harrell_c"].quantile([.25, .5, .75]).to_numpy()
        axb.plot([i + 0.06, i + 0.06], [q[0], q[2]], color=SIG, lw=1.6, solid_capstyle="butt", zorder=3); axb.scatter([i + 0.06], [q[1]], s=9, color="white", ec=SIG, lw=0.7, zorder=4)
        axb.plot([i - 0.36, i + 0.44], [d["grade_c"].median()] * 2, color=GRADE, lw=1.0, zorder=3)
        axb.text(i, 0.225, "%d signatures\n%d events / %d" % (len(d), d["events"].max(), d["n"].max()), ha="center", va="bottom", fontsize=5.2, color=MUTE, linespacing=1.1)
    axb.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    axb.text(3.44, pc[pc.cohort == ORDER[3]]["grade_c"].median() + 0.012, "Gleason grade group", color=GRADE, fontsize=5.8, va="bottom", ha="right")
    axb.set_xticks(range(4)); axb.set_xticklabels([COH[c] for c in ORDER]); axb.set_xlim(-0.55, 3.6); axb.set_ylim(0.22, 0.9); axb.set_yticks([0.3, 0.5, 0.7, 0.9])
    axb.set_ylabel("C-index in each cohort"); axb.tick_params(axis="x", length=0)

    # ---- c-e: Kaplan-Meier
    km = pd.read_csv(FI / "km_curves.csv", dtype={"group": str}); kmeta = pd.read_csv(FI / "km_meta.csv", dtype={"group": str, "sig_id": str})
    km.to_csv(SD / "Fig2c-e_km_curves.csv", index=False); kmeta.to_csv(SD / "Fig2c-e_km_meta.csv", index=False)
    ids = {r.panel: str(r.sig_id).split(".")[0] for r in kmeta.itertuples()}
    cval = {p: float(sig.loc[sig["sig_id"] == ids[p], "pooled_c"].iloc[0]) for p in ("top", "median")}
    tert = ["#8DB0D6", "#3F74AD", "#122F52"]
    panels = [("grade", "Gleason grade group", ["1", "2", "3", "4-5"], ["#EFA58E", "#D9694A", "#A8361C", "#5E180B"], "c", "C = %.2f" % gmed),
              ("median", "Median signature", ["Low", "Middle", "High"], tert, "d", "C = %.2f" % cval["median"]),
              ("top", "Top signature", ["Low", "Middle", "High"], tert, "e", "C = %.2f" % cval["top"])]
    for j, (p, title, groups, cols, lab, ctext) in enumerate(panels):
        x0 = 93 + j * 28
        a = axmm(fig, x0, 97, 20, 33); tag(fig, x0 - 10 if j == 0 else x0 - 5.5, 92, lab)
        for g, col in zip(groups, cols):
            d = km[(km.panel == p) & (km.group == g)]
            a.step(d["time"], d["surv"], where="post", color=col, lw=1.0)
        pv = kmeta.loc[kmeta.panel == p, "p_stratified_lr"].iloc[0]
        a.set_xlim(0, 120); a.set_ylim(0, 1.02); a.set_xticks([0, 48, 96]); a.set_yticks([0, 0.5, 1.0])
        a.text(0, 1.035, title, transform=a.transAxes, ha="left", va="bottom", fontsize=6)
        a.text(0.05, 0.04, ctext + "\nP " + ("< 0.001" if pv < 0.001 else "= %.2f" % pv), transform=a.transAxes, ha="left", va="bottom", fontsize=5.6, linespacing=1.15)
        if j == 0:
            a.set_ylabel("BCR-free survival")
        else:
            a.tick_params(labelleft=False)
        if j == 1:
            a.set_xlabel("Months after radical prostatectomy")
        # numbers at risk, below the axis title
        t = axmm(fig, x0, 138.5, 20, 10.5); t.set_xlim(0, 120); t.set_ylim(len(groups) - 0.4, -1.1); t.axis("off")
        t.text(-46 if j == 0 else -16, -0.95, "At risk, by grade group" if j == 0 else "by score tertile", fontsize=5.2, color=MUTE, ha="left", va="center")
        for k, (g, col) in enumerate(zip(groups, cols)):
            r = kmeta[(kmeta.panel == p) & (kmeta.group == g)].iloc[0]
            t.text(-16, k, g.replace("-", "–"), fontsize=5.2, color=col, ha="right", va="center", fontweight="bold")
            for tt in (0, 48, 96):
                t.text(tt, k, "%d" % r["at_risk_%d" % tt], fontsize=5.2, ha="center", va="center", color=INK)
    save(fig, "Figure2_external_discrimination")


# =============================================================== Figure 3: benchmarks
def fig3():
    fig = fig_mm(140)
    rh = pd.read_csv(FI / "random_hist.csv"); mp = pd.read_csv(EV / "exploratory/metapcna_cohorts.csv").set_index("cohort")
    # ---- a: ridgelines per cohort
    tag(fig, 0, 0, "a")
    xs = np.linspace(0.28, 0.92, 400); rows = []
    for i, c in enumerate(ORDER):
        a = axmm(fig, 21, 4 + i * 19.5, 62, 17.5)
        r = rh[rh.cohort == c]; dens = r["count"] / r["count"].sum() / 0.01
        a.fill_between(r["bin_left"] + 0.005, dens, step="mid", color=RAND_L, lw=0, zorder=1)
        d = pairs[pairs.cohort == c]
        k = stats.gaussian_kde(d["harrell_c"]); a.plot(xs, k(xs), color=SIG, lw=1.0, zorder=3); a.fill_between(xs, k(xs), color=SIG, alpha=0.12, lw=0, zorder=2)
        top = max(dens.max(), k(xs).max()) * 1.18
        refs = [("Proliferation index", mp.loc[c, "c_metapcna"], PROLIF), ("Gleason grade group", d["grade_c"].median(), GRADE), ("Clinical model", d["c_clin"].median(), CLIN)]
        for name, v, col in refs:
            a.plot([v, v], [0, top * 0.80], color=col, lw=1.0, zorder=4)
            rows.append({"cohort": c, "benchmark": name, "c": v})
        a.axvline(0.5, color=INK, lw=0.5, ls=DASH, zorder=0)
        a.set_xlim(0.28, 0.92); a.set_ylim(0, top); a.set_yticks([]); a.spines["left"].set_visible(False)
        a.text(-0.02, 0.35, "%s\n%d events / %d" % (COH[c], d["events"].max(), d["n"].max()), transform=a.transAxes, ha="right", va="center", fontsize=5.8, linespacing=1.15)
        if i == 0:
            a.text(0.488, top * 0.93, "Random gene sets", color=MUTE, fontsize=5.8, ha="right", va="center")
            a.text(0.60, top * 0.93, "Published signatures", color=SIG, fontsize=5.8, ha="left", va="center")
        if i < 3:
            a.tick_params(labelbottom=False)
        else:
            a.set_xlabel("C-index for biochemical recurrence")
    # shared key for the three reference lines
    k = axmm(fig, 21, 4 + 4 * 19.5 + 6.5, 62, 4); k.axis("off"); k.set_xlim(0, 62); k.set_ylim(0, 4)
    for x0, (name, col) in zip((0, 22, 46), (("Proliferation index", PROLIF), ("Gleason grade group", GRADE), ("Clinical model", CLIN))):
        k.plot([x0, x0 + 3], [2, 2], color=col, lw=1.0); k.text(x0 + 4, 2, name, va="center", fontsize=5.8, color=col)
    pd.DataFrame(rows).to_csv(SD / "Fig3a_benchmarks_by_cohort.csv", index=False); rh.to_csv(SD / "Fig3a_random_gene_set_histograms.csv", index=False)

    # ---- b: percentile among random gene sets
    ax = axmm(fig, 103, 4, 64, 32); tag(fig, 90, 0, "b")
    r = pairs["rand_pct"].dropna(); pairs[["sig_id", "cohort", "harrell_c", "rand_pct", "beats_random_95"]].to_csv(SD / "Fig3b_random_percentile.csv", index=False)
    nn, _, patches = ax.hist(r, bins=np.arange(0, 101, 5), color=RAND_L, ec="white", lw=0.5)
    patches[-1].set_facecolor(SIG)
    ax.axhline(len(r) / 20, color=INK, lw=0.5, ls=DASH); ax.text(1, len(r) / 20 + 5, "Expected for random gene sets", fontsize=5.8, va="bottom")
    ax.text(93, nn[-1] * 0.93, "%.0f%% above the\n95th percentile" % (100 * pairs.loc[pairs["rand_pct"].notna(), "beats_random_95"].astype(bool).mean()), ha="right", va="top", color=SIG, fontsize=6, linespacing=1.15)
    ax.set_xlim(0, 100); ax.set_xlabel("Percentile of the signature among 1000 expression-matched\nrandom gene sets"); ax.set_ylabel("Evaluations")

    # ---- c: incremental value
    ax = axmm(fig, 103, 53, 64, 34); tag(fig, 90, 49, "c")
    inc = sig[["sig_id", "mean_delta_c_corrected", "lr_q_bh_combined"]].dropna(subset=["mean_delta_c_corrected"]); inc.to_csv(SD / "Fig3c_incremental_value.csv", index=False)
    q = inc["lr_q_bh_combined"] < 0.05
    ax.hist([inc.loc[~q, "mean_delta_c_corrected"], inc.loc[q, "mean_delta_c_corrected"]], bins=np.arange(-0.016, 0.0441, 0.004), stacked=True, color=[RAND_L, SIG], ec="white", lw=0.5)
    ax.axvline(0, color=INK, lw=0.5, ls=DASH)
    ax.text(0.99, 0.95, "Likelihood-ratio test, FDR < 0.05 (n = %d)" % q.sum(), transform=ax.transAxes, ha="right", va="top", color=SIG, fontsize=5.8)
    ax.text(0.99, 0.84, "FDR ≥ 0.05 (n = %d)" % (~q).sum(), transform=ax.transAxes, ha="right", va="top", color=MUTE, fontsize=5.8)
    ax.text(0.99, 0.66, minus("Median gain %+.4f" % inc["mean_delta_c_corrected"].median()) + "\nClinical model alone, C = %.2f to %.2f" % (pairs.groupby("cohort")["c_clin"].median().min(), pairs.groupby("cohort")["c_clin"].median().max()), transform=ax.transAxes, ha="right", va="top", fontsize=5.8, linespacing=1.2)
    ax.set_xlabel("Optimism-corrected gain in C-index over the clinical model"); ax.set_ylabel("Signatures")

    # ---- d, e: paired differences
    for (col, lo, hi, ylab, key, lab, x0, lims, ccol) in (
            ("d_grade", "d_grade_lo", "d_grade_hi", "Difference in C-index,\nsignature − Gleason grade group", "d", "Fig3d_delta_vs_grade", 21, (-0.55, 0.25), GRADE),
            ("d_mp", "d_mp_lo", "d_mp_hi", "Difference in C-index,\nsignature − proliferation index", "e", "Fig3e_delta_vs_proliferation", 108, (-0.45, 0.45), PROLIF)):
        ax = axmm(fig, x0, 100, 59, 32); tag(fig, x0 - 21 if x0 == 21 else x0 - 18, 96, key)
        d = e12.sort_values(col).reset_index(drop=True); d["rank"] = np.arange(1, len(d) + 1); d.to_csv(SD / (lab + ".csv"), index=False)
        dn, upp = d[hi] < 0, d[lo] > 0; mid = ~(dn | upp)
        ax.vlines(d["rank"], d[lo].clip(lower=lims[0]), d[hi].clip(upper=lims[1]), color=RAND_LL, lw=0.5, zorder=1)
        ax.scatter(d["rank"][mid], d[col][mid], s=3, color=RAND, lw=0, zorder=3)
        ax.scatter(d["rank"][dn], d[col][dn], s=4, color=ccol, lw=0, zorder=4); ax.scatter(d["rank"][upp], d[col][upp], s=4, color=SIG, lw=0, zorder=4)
        ax.axhline(0, color=INK, lw=0.5, ls=DASH, zorder=2)
        ax.set_xlim(0, len(d) + 1); ax.set_ylim(*lims); ax.set_xlabel("Signatures ranked by difference"); ax.set_ylabel(ylab)
        ax.text(0.02, 0.97, "Signature better (n = %d)" % upp.sum(), transform=ax.transAxes, color=SIG, fontsize=5.8, va="top")
        ax.text(0.02, 0.88, "Signature worse (n = %d)" % dn.sum(), transform=ax.transAxes, color=ccol, fontsize=5.8, va="top")
        ax.text(0.02, 0.79, "95%% CI includes 0 (n = %d)" % mid.sum(), transform=ax.transAxes, color=MUTE, fontsize=5.8, va="top")
        ax.text(0.98, 0.05, minus("Median %+.3f" % d[col].median()), transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8)
    save(fig, "Figure3_benchmarks")


# =============================================================== Figure 4: reported vs external
def fig4():
    ir = pd.read_csv(EV / "exploratory/internal_reported.csv", dtype={"sig_id": str})
    three = pd.read_csv(EV / "posthoc/V1_three_step.csv", dtype={"sig_id": str})
    used = pd.read_csv(EV / "posthoc/V1_used_cohort_recomputation.csv", dtype={"sig_id": str})
    v5 = pd.read_csv(EV / "posthoc/V5_calibration_slope_signatures.csv", dtype={"sig_id": str})
    fig = fig_mm(132)
    # ---- a: scatter with marginals
    tag(fig, 0, 0, "a")
    d = ir.dropna(subset=["reported_training", "pooled_c"]); d.to_csv(SD / "Fig4a_reported_vs_external.csv", index=False)
    ax = axmm(fig, 14, 12, 52, 52); axt = axmm(fig, 14, 3.5, 52, 7.5, sharex=ax); axr = axmm(fig, 67.5, 12, 7.5, 52, sharey=ax)
    isc = d["training_metric"] == "C-index"
    ax.plot([0.35, 1], [0.35, 1], color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(d.loc[~isc, "reported_training"], d.loc[~isc, "pooled_c"], s=7, color=REPORT, alpha=0.55, lw=0, zorder=3)
    ax.scatter(d.loc[isc, "reported_training"], d.loc[isc, "pooled_c"], s=8, facecolor="white", edgecolor=REPORT, lw=0.7, zorder=4)
    rho = stats.spearmanr(d["reported_training"], d["pooled_c"]).statistic
    ax.set_xlim(0.52, 1.0); ax.set_ylim(0.38, 1.0); ax.set_xlabel("Discrimination reported in the training set"); ax.set_ylabel("Pooled external C-index (this study)")
    ax.text(0.03, 0.97, "Reported AUC (n = %d)" % (~isc).sum(), transform=ax.transAxes, color=REPORT, fontsize=5.8, va="top")
    ax.text(0.03, 0.90, "Reported C-index (n = %d), open" % isc.sum(), transform=ax.transAxes, color=REPORT, fontsize=5.8, va="top")
    ax.text(0.03, 0.80, minus("%.0f%% below identity\nSpearman ρ = %.2f" % (100 * (d["gap_training"] > 0).mean(), rho)), transform=ax.transAxes, ha="left", va="top", fontsize=5.8, linespacing=1.2)
    ax.text(0.79, 0.83, "Identity", rotation=45, rotation_mode="anchor", fontsize=5.4, color=MUTE, ha="center", va="bottom")
    axt.hist(d["reported_training"], bins=np.arange(0.52, 1.001, 0.02), color=REPORT_L, ec="white", lw=0.4); axt.axis("off")
    axr.hist(d["pooled_c"], bins=np.arange(0.38, 1.001, 0.02), orientation="horizontal", color=SIG_L, ec="white", lw=0.4); axr.axis("off")
    axt.axvline(d["reported_training"].median(), color=REPORT, lw=0.9); axr.axhline(d["pooled_c"].median(), color=SIG, lw=0.9)

    # ---- b: three-step decline
    ax = axmm(fig, 96, 6, 73, 58); tag(fig, 82, 0, "b")
    t = three[["sig_id", "reported", "ours", "pooled_c"]].dropna(); t.to_csv(SD / "Fig4b_three_step.csv", index=False)
    cols, xs, cc = ["reported", "ours", "pooled_c"], [0, 1, 2], [REPORT, RAND, SIG]
    for _, r in t.iterrows():
        ax.plot([0.13, 0.87, 1.13, 1.87], [r["reported"], r["ours"], r["ours"], r["pooled_c"]], color=RAND_L, lw=0.25, alpha=0.55, zorder=1)
    for xx, c, col in zip(xs, cols, cc):
        q = t[c].quantile([.25, .5, .75]).to_numpy()
        ax.add_patch(plt.Rectangle((xx - 0.11, q[0]), 0.22, q[2] - q[0], fc="white", ec=col, lw=0.9, zorder=3)); ax.plot([xx - 0.11, xx + 0.11], [q[1], q[1]], color=col, lw=1.5, zorder=4)
        ax.text(xx, 0.975, "%.2f" % q[1], ha="center", va="top", color=col, fontsize=7.5, fontweight="bold")
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=0)
    ax.set_xticks(xs); ax.set_xticklabels(["Reported by the article\n(training set)", "Recomputed here\n(same training cohort)", "Recomputed here\n(independent cohorts)"])
    ax.set_xlim(-0.4, 2.4); ax.set_ylim(0.38, 1.0); ax.set_ylabel("C-index or AUC"); ax.tick_params(axis="x", length=0)
    ax.set_xlabel("n = %d signatures with all three values; boxes show the median and IQR" % len(t), color=MUTE, fontsize=5.8)

    # ---- c: used vs unused cohorts
    ax = axmm(fig, 14, 85, 38, 38); tag(fig, 0, 77, "c")
    geo = used[used["cohort"].isin(list(COH) + ["GSE116918"])].groupby("sig_id")["harrell_c"].median().rename("used")
    j = pd.concat([geo, pairs.groupby("sig_id")["harrell_c"].median().rename("unused")], axis=1).dropna(); j.to_csv(SD / "Fig4c_used_vs_unused.csv")
    diff = j["used"] - j["unused"]
    ax.plot([0.3, 0.9], [0.3, 0.9], color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(j["unused"], j["used"], s=6, color=SIG, alpha=0.55, lw=0, zorder=3)
    ax.set_xlim(0.33, 0.85); ax.set_ylim(0.33, 0.85); ax.set_xticks([0.4, 0.6, 0.8]); ax.set_yticks([0.4, 0.6, 0.8])
    ax.set_xlabel("C-index in cohorts not used by the article"); ax.set_ylabel("C-index in cohorts\nused by the article")
    ax.text(0, 1.035, minus("%.0f%% above identity (n = %d)\nMedian difference %+.3f, P = %.4f" % (100 * (diff > 0).mean(), len(j), diff.median(), stats.wilcoxon(diff).pvalue)), transform=ax.transAxes, va="bottom", fontsize=5.8, linespacing=1.2)

    # ---- d: external validation reported by the authors vs ours
    ax = axmm(fig, 72, 85, 38, 38); tag(fig, 58, 77, "d")
    e = ir.dropna(subset=["reported_external", "pooled_c"]); e.to_csv(SD / "Fig4d_reported_external_vs_ours.csv", index=False)
    ax.plot([0.3, 1], [0.3, 1], color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(e["reported_external"], e["pooled_c"], s=6, color=REPORT, alpha=0.55, lw=0, zorder=3)
    ax.set_xlim(0.45, 1.0); ax.set_ylim(0.38, 0.93); ax.set_xlabel("External validation reported by the article"); ax.set_ylabel("Pooled external C-index\n(this study)")
    ax.text(0, 1.035, minus("%.0f%% below identity (n = %d)\nMedians %.2f and %.2f" % (100 * (e["gap_external"] > 0).mean(), len(e), e["reported_external"].median(), e["pooled_c"].median())), transform=ax.transAxes, va="bottom", fontsize=5.8, linespacing=1.2)

    # ---- e: calibration slope, by cohort (the pooled value is dominated by GSE54460, see FigS7)
    ax = axmm(fig, 130, 85, 39, 38); tag(fig, 117, 77, "e")
    v5p = pd.read_csv(EV / "posthoc/V5_calibration_slope_pairs.csv", dtype={"sig_id": str}); v5p.to_csv(SD / "Fig4e_calibration_slope_by_cohort.csv", index=False)
    v5.to_csv(SD / "Fig4e_calibration_slope_pooled.csv", index=False)
    rng = np.random.default_rng(7); below = []
    ax.axhline(1, color=INK, lw=0.8, zorder=1); ax.axhline(0, color=INK, lw=0.5, ls=DASH, zorder=1)
    for i, c in enumerate(ORDER):
        d = v5p[v5p.cohort == c]; med = d["slope"].median(); below.append(100 * (d["slope_hi"] < 1).mean())
        ax.scatter(i + rng.uniform(-0.26, 0.26, len(d)), d["slope"].clip(-1.5, 3.0), s=2.2, color=SIG_M, alpha=0.55, lw=0, zorder=2)
        ax.plot([i - 0.36, i + 0.36], [med, med], color=SIG, lw=1.5, zorder=3)
        ax.text(i, 3.28, "%.2f" % med, ha="center", va="bottom", fontsize=5.8, color=SIG)
    ax.set_xticks(range(4)); ax.set_xticklabels([COH[c] for c in ORDER], fontsize=5.2); ax.tick_params(axis="x", length=0); ax.set_xlim(-0.55, 3.55)
    ax.set_ylim(-1.65, 3.75); ax.set_yticks([-1, 0, 1, 2, 3]); ax.set_yticklabels([minus("-1"), "0", "1", "2", "≥ 3"]); ax.set_ylabel("Calibration slope of the\npublished linear predictor")
    ax.text(0, 1.035, "Median slope by cohort (ideal = 1); 95%% CI\nbelow 1 in %.0f%%, %.0f%%, %.0f%% and %.0f%%" % tuple(below), transform=ax.transAxes, va="bottom", fontsize=5.8, linespacing=1.2)
    save(fig, "Figure4_reported_vs_external")


# =============================================================== Figure 5: redundancy and correlates
def fig5():
    R = pd.read_csv(EV / "posthoc/V4_spearman_TCGA.csv", index_col=0); R.index = R.index.astype(str); R.columns = R.columns.astype(str)
    cl = pd.read_csv(EV / "posthoc/V4_clusters_TCGA.csv", dtype={"sig_id": str}).set_index("sig_id")
    fig = fig_mm(155)
    # ---- a: clustered correlation heatmap
    tag(fig, 0, 0, "a")
    Z = linkage(squareform(1 - R.to_numpy(), checks=False), "average"); o = leaves_list(Z)
    Rs = R.iloc[o, o]; ids = Rs.index; n = len(ids)
    S = 78
    ax = axmm(fig, 9, 11, S, S)
    im = ax.imshow(Rs.to_numpy(), cmap=R_CMAP, vmin=-1, vmax=1, interpolation="none", extent=(0, n, n, 0))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.4)
    ax.set_xlabel("%d signatures scored in TCGA-PRAD (n = 497 tumours),\nordered by hierarchical clustering" % n)
    P_CMAP = LinearSegmentedColormap.from_list("p", ["#3E6E8E", "#FFFFFF", PROLIF])
    axp = axmm(fig, 9, 3.5, S, 2.6); imp = axp.imshow(cl.loc[ids, "rho_metapcna"].to_numpy()[None, :], aspect="auto", cmap=P_CMAP, vmin=-0.9, vmax=0.9, interpolation="none")
    axc = axmm(fig, 9, 7, S, 2.6); pc = sig.set_index("sig_id")["pooled_c"].reindex(ids).to_numpy()
    imc = axc.imshow(np.ma.masked_invalid(pc)[None, :], aspect="auto", cmap=C_CMAP, norm=TwoSlopeNorm(0.5, 0.3, 0.8), interpolation="none"); axc.set_facecolor("#E4E4E4")
    for a_, labl in ((axp, "ρ with proliferation index"), (axc, "Pooled external C-index")):
        a_.set_xticks([]); a_.set_yticks([]); a_.text(1.012, 0.5, labl, transform=a_.transAxes, fontsize=5.4, va="center", ha="left")
        for sp in a_.spines.values():
            sp.set_visible(False)
    big = cl["cluster"].value_counts().index[0]; pos = np.where(cl.loc[ids, "cluster"].to_numpy() == big)[0]
    if len(pos):
        ax.plot([pos.min(), pos.max() + 1, pos.max() + 1, pos.min(), pos.min()], [pos.min(), pos.min(), pos.max() + 1, pos.max() + 1, pos.min()], color=INK, lw=0.6)
        ax.text(pos.max() + 3 if pos.max() < n * 0.6 else pos.min() - 3, pos.mean(), "Largest cluster\n%d signatures\nmedian ρ with\nproliferation %.2f" % (len(pos), cl.loc[cl["cluster"] == big, "rho_metapcna"].median()),
                fontsize=5.4, va="center", ha="left" if pos.max() < n * 0.6 else "right", linespacing=1.15, bbox=dict(fc="white", ec="none", alpha=0.82, pad=1.2))
    for x0, im_, ticks, labs, ttl in ((9, im, [-1, 0, 1], [minus("-1"), "0", "1"], "ρ between signature scores"),
                                      (39, imp, [-0.9, 0, 0.9], [minus("-0.9"), "0", "0.9"], "ρ with proliferation index"),
                                      (69, imc, [0.3, 0.5, 0.8], ["0.3", "0.5", "0.8"], "Pooled external C-index")):
        cax = axmm(fig, x0, 11 + S + 13, 18, 1.6); cb = fig.colorbar(im_, cax=cax, orientation="horizontal", ticks=ticks); cb.outline.set_linewidth(0.4); cb.ax.tick_params(length=1.5, labelsize=5.4)
        cb.ax.set_xticklabels(labs); cax.text(0, 1.5, ttl, transform=cax.transAxes, fontsize=5.4, va="bottom", ha="left")
    Rs.to_csv(SD / "Fig5a_signature_correlation_TCGA.csv")

    # ---- b: recurrent genes
    ax = axmm(fig, 127, 4, 20, 86); tag(fig, 113, 0, "b")
    tt = terms[terms["sig_id"].isin(sig["sig_id"])]
    gl = tt.assign(g=tt["HGNC 现行符号"].str.split("|")).explode("g").drop_duplicates(["sig_id", "g"])
    vc = gl.loc[~gl["g"].isin(["", "未映射"]), "g"].value_counts()
    vc = vc.rename_axis("gene").reset_index(name="n").sort_values(["n", "gene"], ascending=[False, True]).set_index("gene")["n"]  # ties in alphabetical order
    vc.rename("n_signatures").to_csv(SD / "Fig5b_gene_recurrence.csv")
    mpg = set(pd.read_csv(EV / "exploratory/metapcna_genes.csv")["gene"])
    top = vc[vc >= 5][::-1]  # every gene used by at least five signatures
    ax.barh(range(len(top)), top.values, color=[PROLIF if g in mpg else SIG_M for g in top.index], height=0.68)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontstyle="italic", fontsize=5.6); ax.set_ylim(-0.6, len(top) - 0.4)
    ax.set_xlabel("Signatures containing\nthe gene (of %d)" % len(sig)); ax.tick_params(axis="y", length=0); ax.set_xticks([0, 3, 6, 9]); ax.set_xlim(0, 9.6)
    ax.add_patch(plt.Rectangle((5.5, 8.66), 0.75, 0.68, fc=PROLIF, ec="none", clip_on=False)); ax.text(6.6, 9.0, "Gene in the proliferation index", fontsize=5.4, va="center", clip_on=False)
    ins = axmm(fig, 153, 57, 16, 22)
    cnt = vc.clip(upper=4).value_counts().sort_index()
    ins.bar(cnt.index, cnt.values / len(vc) * 100, color=RAND_L, width=0.7); ins.set_xticks([1, 2, 3, 4]); ins.set_xticklabels(["1", "2", "3", "≥ 4"])
    ins.set_ylabel("Genes (%)"); ins.set_xlabel("Signatures\nper gene"); ins.set_ylim(0, 100); ins.set_yticks([0, 50, 100])
    ins.text(1, cnt.iloc[0] / len(vc) * 100 + 3, "%.0f%%" % (cnt.iloc[0] / len(vc) * 100), ha="center", fontsize=5.6)
    ins.text(0.98, 0.98, "%d distinct\ngenes" % len(vc), transform=ins.transAxes, ha="right", va="top", fontsize=5.4, color=MUTE, linespacing=1.1)

    # ---- c: meta-regression
    ax = axmm(fig, 44, 117, 36, 30); tag(fig, 0, 111, "c")
    mr = pd.read_csv(EV / "posthoc/V3_metaregression.csv"); mr = mr[mr["model"] == "多变量"].set_index("term")
    names = [("tcga_only_training", "Trained on TCGA only"), ("ext_reported", "External validation reported\nin the article"), ("log2_genes", "Number of genes\n(per doubling)"),
             ("median_rho_metapcna", "Correlation with\nproliferation index (per 1.0)")]
    m = mr.loc[[k for k, _ in names]].assign(label=[v for _, v in names]); m.to_csv(SD / "Fig5c_metaregression.csv")
    yy = np.arange(len(m))[::-1]
    col = [GRADE if hi < 0 else (SIG if lo > 0 else RAND) for lo, hi in zip(m["lo"], m["hi"])]
    ax.hlines(yy, m["lo"], m["hi"], color=col, lw=1.0); ax.scatter(m["coef_logitC"], yy, s=13, color=col, zorder=3, lw=0)
    ax.axvline(0, color=INK, lw=0.5, ls=DASH); ax.set_yticks(yy); ax.set_yticklabels(m["label"], linespacing=1.05); ax.tick_params(axis="y", length=0)
    ax.text(1.05, len(m) - 0.5, "P", transform=ax.get_yaxis_transform(), fontsize=6, fontweight="bold", va="bottom", clip_on=False)
    for y_, r in zip(yy, m.itertuples()):
        ax.text(1.05, y_, "%.3f" % r.p, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, clip_on=False)
    ax.set_xlim(-0.3, 0.22); ax.set_xticks([-0.2, 0, 0.2]); ax.set_xticklabels([minus("-0.2"), "0", "0.2"]); ax.set_ylim(-0.6, len(m) - 0.4)
    ax.set_xlabel("Difference in logit C-index (95% CI)"); ax.spines["left"].set_visible(False)

    # ---- d: external C by year
    ax = axmm(fig, 108, 117, 61, 30); tag(fig, 94, 111, "d")
    d = meta[["sig_id", "year", "pooled_c"]].copy(); d["yb"] = d["year"].clip(lower=2019); d.to_csv(SD / "Fig5d_c_by_year.csv", index=False)
    rng = np.random.default_rng(2); yrs = sorted(d["yb"].unique())
    for y_ in yrs:
        v = d.loc[d["yb"] == y_, "pooled_c"]
        ax.scatter(y_ + rng.uniform(-0.22, 0.22, len(v)), v, s=4, color=SIG_M, alpha=0.6, lw=0, zorder=2)
        ax.plot([y_ - 0.32, y_ + 0.32], [v.median()] * 2, color=SIG, lw=1.4, zorder=3)
        ax.text(y_, 0.775, "n = %d" % len(v), ha="center", va="top", fontsize=5.2, color=MUTE)
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    sp = stats.spearmanr(meta["year"], meta["pooled_c"])
    ax.set_xticks(yrs); ax.set_xticklabels(["≤ 2019"] + [str(y_) for y_ in yrs[1:]]); ax.set_ylim(0.38, 0.78); ax.set_xlim(2018.4, 2026.6); ax.set_yticks([0.4, 0.5, 0.6, 0.7])
    ax.set_ylabel("Pooled external C-index"); ax.set_xlabel("Year of publication")
    ax.text(0, 1.035, minus("Spearman ρ = %.2f, P = %.2f; lines show medians" % (sp.statistic, sp.pvalue)), transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8)
    save(fig, "Figure5_redundancy_and_correlates")


# =============================================================== Figure 6: robustness
def fig6():
    rows = [("primary", "Primary analysis"), ("S1_cov100", "Complete gene coverage only"), ("S2a_raw_log2", "Unstandardised log2 expression"),
            ("S2b_tcga_sd", "Scaled to TCGA standard deviations"), ("S3_plus_GSE116918", "Plus radiotherapy cohort (GSE116918)"),
            ("S4_without_MSKCC", "Without MSKCC"), ("S5_GSE54460_only", "GSE54460 only"), ("S6_endpoint_BCR_RFS", "Developed for BCR or RFS"),
            ("S7_plus_low_precision", "Plus coefficients from figures or HRs"), ("S8_equal_weight", "Equal-weight signatures"),
            ("S10_GSE54460_no_repeats", "GSE54460 without repeat libraries"), ("E3_met_GSE116918", "Metastasis endpoint (GSE116918)")]
    tab = []
    for tg, name in rows:
        s = pd.read_csv(EV / tg / "signatures.csv", dtype={"sig_id": str}); j = json.loads((EV / tg / "summary.json").read_text(encoding="utf-8"))
        q = s["pooled_c"].quantile([.1, .25, .5, .75, .9]).to_numpy()
        tab.append({"analysis": name, "n": len(s), "p10": q[0], "q1": q[1], "median": q[2], "q3": q[3], "p90": q[4], "beats_grade": 100 * j["share_beats_grade_papers"], "beats_random": 100 * j["share_beats_random_papers"]})
    tab = pd.DataFrame(tab); tab.to_csv(SD / "Fig6a_sensitivity.csv", index=False)
    fig = fig_mm(125)
    # ---- a
    ax = axmm(fig, 50, 8, 62, 52); tag(fig, 0, 0, "a")
    y = np.arange(len(tab))[::-1]
    ax.axvspan(tab["q1"][0], tab["q3"][0], color=SIG_LL, lw=0, zorder=0)
    ax.hlines(y, tab["p10"], tab["p90"], color=RAND_L, lw=0.7, zorder=1); ax.hlines(y, tab["q1"], tab["q3"], color=[SIG] + [INK] * (len(tab) - 1), lw=2.0, zorder=2)
    ax.scatter(tab["median"], y, s=11, color="white", ec=[SIG] + [INK] * (len(tab) - 1), lw=0.8, zorder=3)
    ax.axvline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(tab["analysis"]); ax.set_xlim(0.44, 0.70); ax.set_ylim(-0.7, len(tab) - 0.3); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xlabel("Pooled external C-index across signatures: median, IQR and 10th–90th percentiles")
    tr = ax.get_yaxis_transform()
    for xx, h in ((1.10, "Signatures"), (1.33, "Better than\ngrade (%)"), (1.60, "Better than random\ngene sets (%)")):
        ax.text(xx, len(tab) - 0.25, h, transform=tr, ha="center", va="bottom", fontsize=6, fontweight="bold", clip_on=False, linespacing=1.1)
    for yy, r in zip(y, tab.itertuples()):
        ax.text(1.10, yy, "%d" % r.n, transform=tr, ha="center", va="center", clip_on=False)
        ax.text(1.33, yy, "%.1f" % r.beats_grade, transform=tr, ha="center", va="center", clip_on=False)
        ax.text(1.60, yy, "%.0f" % r.beats_random, transform=tr, ha="center", va="center", clip_on=False)

    # ---- b: Harrell vs Uno / AUC
    ax = axmm(fig, 14, 83, 34, 34); tag(fig, 0, 75, "b")
    p = pairs[["sig_id", "cohort", "harrell_c", "uno_c", "auc36", "auc60"]]; p.to_csv(SD / "Fig6b_metrics.csv", index=False)
    ax.plot([0.25, 0.9], [0.25, 0.9], color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.scatter(p["harrell_c"], p["auc36"], s=3.5, color=SIG_M, alpha=0.5, lw=0, zorder=2)
    ax.set_xlim(0.28, 0.88); ax.set_ylim(0.28, 0.88); ax.set_xticks([0.4, 0.6, 0.8]); ax.set_yticks([0.4, 0.6, 0.8])
    ax.set_xlabel("Harrell C-index"); ax.set_ylabel("Time-dependent AUC at 36 months")
    ax.text(0, 1.035, minus("Spearman ρ = %.2f (n = %d evaluations)\nMedians %.2f and %.2f" % (stats.spearmanr(p["harrell_c"], p["auc36"], nan_policy="omit").statistic, p["auc36"].notna().sum(), p["harrell_c"].median(), p["auc36"].median())), transform=ax.transAxes, va="bottom", fontsize=5.8, linespacing=1.2)

    # ---- c: BCR vs metastasis in GSE116918
    ax = axmm(fig, 68, 83, 34, 34); tag(fig, 54, 75, "c")
    a = pd.read_csv(EV / "run_116918/pairs.csv", dtype={"sig_id": str}); b = pd.read_csv(EV / "run_met_116918/pairs.csv", dtype={"sig_id": str})
    j = a.merge(b, on="sig_id", suffixes=("_bcr", "_met")); j[["sig_id", "harrell_c_bcr", "harrell_c_met", "grade_c_bcr", "grade_c_met", "events_bcr", "events_met"]].to_csv(SD / "Fig6c_bcr_vs_metastasis.csv", index=False)
    ax.plot([0.25, 0.9], [0.25, 0.9], color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.axvline(j["grade_c_bcr"].iloc[0], color=GRADE, lw=0.8, zorder=1); ax.axhline(j["grade_c_met"].iloc[0], color=GRADE, lw=0.8, zorder=1)
    ax.scatter(j["harrell_c_bcr"], j["harrell_c_met"], s=5, color=SIG, alpha=0.55, lw=0, zorder=3)
    ax.set_xlim(0.3, 0.8); ax.set_ylim(0.3, 0.8); ax.set_xticks([0.4, 0.6, 0.8]); ax.set_yticks([0.4, 0.6, 0.8])
    ax.set_xlabel("C-index for biochemical recurrence\n(%d events)" % j["events_bcr"].iloc[0]); ax.set_ylabel("C-index for metastasis\n(%d events)" % j["events_met"].iloc[0])
    ax.text(0, 1.035, minus("GSE116918 (n = %d signatures)\nSpearman ρ = %.2f" % (len(j), stats.spearmanr(j["harrell_c_bcr"], j["harrell_c_met"]).statistic)), transform=ax.transAxes, va="bottom", fontsize=5.8, linespacing=1.2)
    ax.text(j["grade_c_bcr"].iloc[0] + 0.012, 0.31, "Gleason\ngrade group", color=GRADE, fontsize=5.6, ha="left", va="bottom", linespacing=1.1)

    # ---- d: precision available in each cohort
    ax = axmm(fig, 123, 83, 46, 34); tag(fig, 109, 75, "d")
    # smallest C-index distinguishable from 0.5 with 80% power (two-sided 5%): 2.8 x the median bootstrap SE of logit C.
    # (Until 20 Nov 2026 the bar was 0.5 + 2.8 x SE of the signature-minus-grade difference, which is a different quantity
    #  and overstated the threshold; V2_mdd_by_cohort.csv still holds that difference and is reported as such.)
    md = pd.read_csv(EV / "posthoc/V2_mdd_by_cohort.csv").set_index("cohort").reindex(ORDER).rename(columns={"median_se": "median_se_diff_vs_grade", "MDD_80pct": "min_detectable_diff_vs_grade"})
    md["median_se_logit_c"] = [pairs.loc[pairs.cohort == c, "logit_c_se"].median() for c in ORDER]
    md["min_c_distinguishable_from_0.5"] = 1 / (1 + np.exp(-2.8 * md["median_se_logit_c"]))
    md["share_of_evaluations_at_or_above"] = [float((pairs.loc[pairs.cohort == c, "harrell_c"] >= md.loc[c, "min_c_distinguishable_from_0.5"]).mean()) for c in ORDER]
    md.to_csv(SD / "Fig6d_precision_by_cohort.csv")
    for i, c in enumerate(ORDER):
        d = pairs[pairs.cohort == c]
        thr = md.loc[c, "min_c_distinguishable_from_0.5"]
        ax.scatter(np.full(len(d), i) + np.random.default_rng(i).uniform(-0.2, 0.2, len(d)), d["harrell_c"], s=2.5, color=SIG_M, alpha=0.5, lw=0, zorder=2)
        ax.plot([i - 0.34, i + 0.34], [thr, thr], color=INK, lw=1.0, zorder=3)
        ax.text(i, thr + 0.014, "%.2f" % thr, ha="center", va="bottom", fontsize=5.6, zorder=5, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.5))
        ax.text(i, 0.225, "%.0f%%" % (100 * (d["harrell_c"] >= thr).mean()), ha="center", va="bottom", fontsize=5.6, color=SIG)
    ax.axhline(0.5, color=INK, lw=0.5, ls=DASH, zorder=1)
    ax.set_xticks(range(4)); ax.set_xticklabels([COH[c] for c in ORDER], fontsize=5.6); ax.set_xlim(-0.55, 3.55); ax.set_ylim(0.22, 0.9); ax.set_yticks([0.3, 0.5, 0.7, 0.9])
    ax.set_ylabel("C-index"); ax.tick_params(axis="x", length=0)
    ax.text(0, 1.107, "Bars: smallest C-index distinguishable from 0.5", transform=ax.transAxes, va="bottom", fontsize=5.8)
    ax.text(0, 1.035, "Percentages: evaluations at or above the bar", transform=ax.transAxes, va="bottom", fontsize=5.8, color=SIG)
    save(fig, "Figure6_robustness")


if __name__ == "__main__":
    for k, f in (("1", fig1), ("2", fig2), ("3", fig3), ("4", fig4), ("5", fig5), ("6", fig6)):
        if k in WHICH:
            f(); print("figure", k, "done")
