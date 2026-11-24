
External validation of 211 published transcriptomic signatures for biochemical recurrence of prostate cancer: comparison with Gleason grade and random gene sets.

The source data of every panel are in the public repository, in the folder `source_data` beside the figure files.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS1_gene_coverage.pdf}\end{center}

\textbf{Fig. S1 Gene coverage and evaluability of the signatures with published coefficients.} \textbf{a} For each cohort, the number of the 247 signatures with published coefficients that could not be evaluated because the article had used the cohort, that had less than 80\% of their gene terms measured, or that were evaluable with 80\% to less than 100\% or with 100\% coverage. \textbf{b} Number of independent prostatectomy cohorts available for each signature; the 211 signatures with at least one cohort entered the primary analysis. \textbf{c} The twelve genes most often missing from a prostatectomy cohort that the article had not used.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS2_every_signature_by_cohort.pdf}\end{center}

\textbf{Fig. S2 C-index of every signature in each prostatectomy cohort.} \textbf{a} GSE54460. \textbf{b} GSE70769. \textbf{c} GSE70768. \textbf{d} MSKCC. Signatures are ranked within each cohort by Harrell's C-index for biochemical recurrence; vertical bars are 95\% bootstrap confidence intervals. The dashed line marks 0.5 and the labelled line the C-index of the Gleason grade group in the same cohort. BCR: biochemical recurrence; CI: confidence interval.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS3_agreement_between_measures.pdf}\end{center}

\textbf{Fig. S3 Agreement between measures of discrimination.} \textbf{a} Harrell's C-index against Uno's C-index truncated at 60 months, for each signature-by-cohort evaluation. \textbf{b} Harrell's C-index against the time-dependent AUC at 60 months. \textbf{c} Time-dependent AUC at 36 months against 60 months. \textbf{d} Pooled Harrell's C-index against the median Uno's C-index of each signature. The line is the line of identity. The AUC was not computed at a time point with fewer than 10 patients at risk or fewer than 10 events, so panels b and c hold fewer evaluations. AUC: area under the curve.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS4_hazard_ratio_per_sd.pdf}\end{center}

\textbf{Fig. S4 Hazard ratio per standard deviation of the score.} \textbf{a} Hazard ratio for biochemical recurrence per standard deviation of each signature score, by cohort, on a logarithmic scale. Points are evaluations; the white circle is the median and the bar the interquartile range. The percentages are the evaluations whose 95\% confidence interval lies entirely above or entirely below 1. \textbf{b} Hazard ratio per standard deviation against Harrell's C-index for the 569 evaluations of the primary analysis. CI: confidence interval; SD: standard deviation.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS5_alternative_analyses.pdf}\end{center}

\textbf{Fig. S5 Each alternative analysis against the primary analysis.} Pooled C-index of each signature in the alternative analysis (vertical axis) against the primary analysis (horizontal axis). \textbf{a} Signatures with complete gene coverage only. \textbf{b} Scores computed from unstandardised log2 expression. \textbf{c} Genes rescaled to their standard deviation in TCGA-PRAD. \textbf{d} With the radiotherapy cohort GSE116918 added. \textbf{e} Without the MSKCC cohort. \textbf{f} GSE54460 only. \textbf{g} GSE54460 without the six patients who had repeat libraries. \textbf{h} Metastasis as the outcome in GSE116918. Each panel gives the number of signatures, the median change from the primary analysis and the Spearman correlation. The line is the line of identity.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS6_signatures_outside_primary.pdf}\end{center}

\textbf{Fig. S6 Signatures analysed outside the primary analysis.} \textbf{a} Pooled external C-index of the 211 signatures with published coefficients (primary analysis), of 25 signatures whose weights were read from figures, given as hazard ratios or taken from the LASSO step, and of 65 signatures published without coefficients and scored with equal weights. White circles are medians and bars interquartile ranges. \textbf{b} The equal-weight signatures ranked by pooled C-index with 95\% confidence intervals. \textbf{c} The 25 signatures with lower-precision weights, ranked in the same way. The labelled line is the median C-index of the Gleason grade group. CI: confidence interval; HR: hazard ratio.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS7_calibration_by_cohort.pdf}\end{center}

\textbf{Fig. S7 Calibration slope of the published linear predictors, by cohort.} \textbf{a} GSE54460. \textbf{b} GSE70769. \textbf{c} GSE70768. \textbf{d} MSKCC. Distribution of the calibration slope, the Cox coefficient of the published linear predictor computed on unstandardised log2 expression. The solid black line marks the ideal slope of 1, the dashed line 0 and the coloured line the cohort median; values are truncated at $-1.5$ and 3. \textbf{e} Calibration slope against the standard deviation of the linear predictor in the evaluation cohort. The slope depends on the expression scale of each platform, which is why it is reported by cohort (post hoc analysis). CI: confidence interval; SD: standard deviation.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS8_radiotherapy_cohort.pdf}\end{center}

\textbf{Fig. S8 Signatures in the radiotherapy cohort GSE116918.} \textbf{a} Harrell's C-index for metastasis (22 events in 248 patients) of the 112 signatures evaluable in this cohort, ranked, with 95\% bootstrap confidence intervals. \textbf{b} Percentile of each signature among 1000 expression-matched random gene sets for metastasis. \textbf{c, d} The same for biochemical recurrence (56 events). In a and c the labelled line is the C-index of the Gleason score group in this cohort. In b and d the dashed line is the count expected in each bin if signatures behaved like random gene sets, and the coloured bar holds signatures above the 95th percentile. The metastasis analysis is exploratory. CI: confidence interval.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS9_heterogeneity_and_precision.pdf}\end{center}

\textbf{Fig. S9 Heterogeneity between cohorts and precision of the pooled estimates.} \textbf{a} Distribution of $I^2$ among the 179 signatures evaluated in two or more cohorts. \textbf{b} Half-width of the 95\% confidence interval of the pooled C-index by the number of cohorts contributing; lines are medians and half-widths are truncated at 0.5. With two cohorts the Hartung--Knapp interval has one degree of freedom and is wide. \textbf{c} Half-width against the pooled C-index, by number of cohorts. CI: confidence interval.

\newpage

\begin{center}\includegraphics[width=\textwidth]{FigureS10_by_identification_route.pdf}\end{center}

\textbf{Fig. S10 Results by the route through which the article was identified.} \textbf{a} Pooled external C-index of the signatures from records taken to full text by the keyword rule (first stage, 142 signatures) and from records screened in duplicate by two runs of an AI model (second stage, 69 signatures). White circles are medians and bars interquartile ranges; the $P$ value is from a Mann--Whitney test. \textbf{b} Percentage of signatures in each group that met each criterion. The figure labels the second-stage records as dual-screened. AI: artificial intelligence; CI: confidence interval; FDR: false discovery rate.

