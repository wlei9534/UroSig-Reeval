# UroSig-Reeval

External re-evaluation of published transcriptomic prognostic signatures for prostate cancer.

Published gene-expression signatures for biochemical recurrence (BCR) after radical prostatectomy were collected from the literature and scored, with the weights printed in each article, in independent public cohorts that the article had not used. Each signature is compared with the Gleason grade group, with a clinical model, and with random gene sets of the same size.

**Status:** public since the manuscript was submitted (24 November 2026). The state of the repository at submission is archived on Zenodo: https://doi.org/10.5281/zenodo.24505570.

## At a glance

| | |
|---|---|
| PubMed records screened (search of 26 October 2026) | 1,011 |
| Analysable signatures (second lock, 14 November 2026) | 375 |
| Signatures in the primary analysis | 211 |
| Signature × cohort evaluations | 569 |
| Cohorts in the primary analysis | GSE54460, GSE70768, GSE70769, MSKCC (GSE21032) |
| Patients / BCR events | 434 / 140 |

## Statistical analysis plan

The plan (SAP v2, finalised 11 November 2026), its English translation and the list of amendments and post hoc analyses are archived on Zenodo: https://doi.org/10.5281/zenodo.24385587. The plan was not registered prospectively. Copies are in `00_protocol/`.

## Layout

The repository keeps the folder names of the working project so that the code runs unchanged.

| Folder | Content |
|---|---|
| `00_protocol/` | SAP (original Chinese PDF, English translation), amendments, PubMed search strategy |
| `02_screening/` | Screening decisions for all 1,011 records (abstracts removed), instructions given to the two screening runs |
| `03_extraction/locked/` | Locked signature tables: one row per signature and one row per gene term, with checksums (`MANIFEST*.txt`). First lock 1 November 2026, second lock 14 November 2026 |
| `03_extraction/v2/` | Comparison, reconciliation and adjudication tables of the duplicate extraction of the second stage |
| `05_code/` | Python package `urosig`, analysis scripts, tests |
| `06_results/eligibility/` | Cohort-use (leakage) matrix, gene coverage and evaluability of every signature in every cohort |
| `06_results/evaluation/` | Results on the first lock |
| `06_results/evaluation_v2/` | Results on the second lock: primary analysis, sensitivity analyses S1–S10, exploratory analyses E1–E4, post hoc analyses V1–V5 |
| `07_submission/` | Figures (PDF, SVG, PNG) with the source data of every panel; tables and additional files of the manuscript |
| `logs/decision_log.md` | Dated log of the decisions and deviations that concern the design, the data, the analysis and its reporting (in Chinese) |
| `docs/data_dictionary.md` | English names for the Chinese column headers and category labels |

Not included: full texts and abstracts of the articles, and expression or clinical data of the cohorts (see below).

## Cohort data

The cohort data are public and are not redistributed here. The scripts that build the processed matrices in `04_cohorts/processed/` are in `05_code/scripts/` (`build_prad_cohorts.py`, `qc_cohorts.py`).

| Cohort | Source |
|---|---|
| GSE54460, GSE70768, GSE70769, GSE116918 | NCBI Gene Expression Omnibus |
| MSKCC | GEO GSE21032; outcome from cBioPortal |
| TCGA-PRAD (reference only, never used for evaluation) | TCGA |

## Reproducing the analysis

Python 3.12 with the packages in `05_code/requirements.txt`, plus matplotlib for the figures.

```
python 05_code/scripts/build_prad_cohorts.py        # build the processed cohort matrices
python 05_code/scripts/run_all.py                   # primary and sensitivity analyses
python 05_code/scripts/pool_results.py              # pooling per signature
python 05_code/scripts/exploratory.py               # E1-E4
python 05_code/scripts/posthoc_v1_v4.py             # V1-V4
python 05_code/scripts/calibration_slope.py         # V5
python 05_code/scripts/figure_inputs.py
python 05_code/scripts/make_figures_v2.py . 07_submission/figures
python 05_code/scripts/make_supp_figures.py . 07_submission/supplementary_figures
```

Set `UROSIG_LOCK=v1` or `UROSIG_LOCK=v2` (default) to choose the extraction lock. Seeds for the bootstrap and the random gene sets are fixed. Tests: `pytest 05_code/tests`.

## Use of artificial intelligence

An AI model (Claude, Anthropic) was used for screening, duplicate extraction, writing and running the code, and drafting documents. Screening decisions and extracted data were checked against the source articles by three authors. Details are in the amendments file and the decision log.

## Authors

Jin Zhao¹*, Song Xue²*, Ruiping Su², Wei Lei³, Tubin Yuan⁴, Xiang Wan¹ (corresponding author)

\* contributed equally

1. Department of Nephrology, Xi'an People's Hospital (Xi'an Fourth Hospital), Xi'an, China
2. Department of Urology, Xi'an People's Hospital (Xi'an Fourth Hospital), Xi'an, China
3. Department of Urology, Pingshan District Central Hospital of Shenzhen, Shenzhen, China
4. Department of Urology, Kweichow Moutai Hospital, Renhuai, Guizhou, China

## Funding

Pingshan District Health System Research Project, Shenzhen (grant number 2025201). The funder had no role in the design of the study, in the collection, analysis and interpretation of data, or in writing.

## Licence

Code: MIT (`LICENSE`). Tables, results and figures produced by this project: CC BY 4.0 (`LICENSE-DATA.md`). The cohort data remain under the terms of their original providers.
