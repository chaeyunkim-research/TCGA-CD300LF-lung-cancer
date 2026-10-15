# CD300LF in lung adenocarcinoma: external validation pipeline

Analysis code for the study of CD300LF expression and prognosis in non-small cell lung cancer. This repository reproduces the external validation of CD300LF in lung adenocarcinoma (LUAD) using public Gene Expression Omnibus (GEO) cohorts, and prepares the TCGA inputs used for gene set enrichment analysis (GSEA).

## Scope

| Manuscript item | Produced by |
| --- | --- |
| Figure 3: forest plot of CD300LF hazard ratios in five GEO cohorts | `analysis_specific/cox_regression_specific.py`, `analysis_specific/forest_plot_specific.py` |
| Figure 4: leave-one-out sensitivity analysis | `analysis_specific/leave_one_out_specific.py` |
| GEO Kaplan–Meier curves (supplementary) | `analysis_specific/km_plot_specific.py` |
| GEO co-expression analyses (supplementary) | `analysis_specific/correlation_table_specific.py`, `analysis_specific/coexpression_clustermap_specific.py` |
| Table 4 and Figure 2: GSEA input files (primary tumors only) for TCGA LUAD/LUSC | `tcga/make_gsea_tumor_only.py` |

TCGA clinicopathological comparisons and multivariable Cox models (Tables 1–2) are not included in this repository.

## Requirements

- Python 3.11
- pandas, numpy, scipy, matplotlib, seaborn, lifelines

```bash
pip install -r requirements.txt
```

GSEA was run separately in the GSEA desktop application (Broad Institute) using the files prepared by `tcga/`.

## Data

Raw data are not included. Download the files below from GEO and place them in `raw_data/<accession>/`. Compressed (`.gz`) files can be used as is.

| Cohort | Platform | Files | Patients analyzed |
| --- | --- | --- | --- |
| GSE31210 | GPL570 (Affymetrix HG-U133 Plus 2.0) | `GSE31210_series_matrix.txt`, `GPL570.txt` | 204 |
| GSE72094 | GPL15048 (Rosetta/Merck custom) | `GSE72094_series_matrix.txt`, `GPL15048.txt` | 398 |
| GSE36471 | GPL9053 (Agilent) | `GSE36471-GPL9053_series_matrix.txt`, `GPL9053.txt` | 115 |
| GSE13213 | GPL6480 (Agilent Whole Human Genome) | `GSE13213_series_matrix.txt.gz`, `GPL6480.txt` | 117 |
| GSE42127 | GPL6884 (Illumina HumanWG-6 v3) | `GSE42127_series_matrix.txt.gz`, `GPL6884.txt` | 133 (adenocarcinoma only) |

TCGA LUAD and LUSC RNA-seq data (Level 3, RSEM normalized, Broad GDAC Firehose run 2016-01-28) go in `raw_data/TCGA_LUAD/` or `raw_data/TCGA_LUSC/`.

### Cohort eligibility

GEO cohorts were included if they contained LUAD patients with overall survival data and were profiled on a platform with a probe targeting CD300LF (NM_139018). GSE68465 and GSE14814 were excluded because their platform (Affymetrix HG-U133A, GPL96) has no CD300LF probe. `cohort_runners/GSE68465_runner_com.py` is kept for reference only and is not used in the analysis.

## How to run

Scripts can be run from any working directory; all paths are resolved by `paths_config.py`.

```bash
# 1. Process each GEO cohort (raw series matrix -> clinical, expression and merged tables)
python cohort_runners/GSE31210_runner_com.py
python cohort_runners/GSE72094_runner_com.py
python cohort_runners/GSE36471_runner_com.py
python cohort_runners/GSE13213_runner_com.py
python cohort_runners/GSE42127_runner_com.py

# 2. Survival and co-expression analyses
python analysis_specific/cox_regression_specific.py
python analysis_specific/forest_plot_specific.py
python analysis_specific/leave_one_out_specific.py
python analysis_specific/km_plot_specific.py
python analysis_specific/correlation_table_specific.py
python analysis_specific/coexpression_clustermap_specific.py

# 3. TCGA inputs for GSEA (primary tumors only, one sample per patient)
#    Place the Firehose RNASeqV2 folders next to the script, then run:
python tcga/make_gsea_tumor_only.py
#    -> CD300LF_GSEA_tumor_only/{LUAD,LUSC}_tumor_only.gct and *_CD300LF_tumor_only.cls
```

GSEA settings used in the manuscript: continuous phenotype (CD300LF), metric Pearson, MSigDB Hallmark v2026.1, 10,000 permutations of type gene_set, gene set size 15–500, no collapse. Tied CD300LF values are offset by 1e-6 so that GSEA accepts the continuous phenotype.

`tcga/build_gsea_cls_specific.py` and `tcga/check_gct_cls_match_com.py` are earlier helpers kept for reference; the manuscript results use `make_gsea_tumor_only.py`.

Outputs are written to `data_processed/<cohort>/` (step 1) and `results/CD300LF/` (step 2).

## Methods implemented

- **Probe mapping:** probes are matched to gene symbols from the platform annotation file. When several probes map to one gene, their values are averaged.
- **Survival time:** converted to months (days / 30.4375) where reported in days.
- **Cox regression:** CD300LF expression is standardized to a z-score within each cohort; univariable Cox models estimate the hazard ratio per 1 SD increase.
- **Meta-analysis:** DerSimonian–Laird random-effects pooling of log hazard ratios, with Cochran's Q and I²; leave-one-out analysis repeats the pooling with each cohort removed.
- **Kaplan–Meier:** median split, log-rank test.
- **Correlation:** Pearson correlation between CD300LF and 14 immune, proliferation and other marker genes defined in `gene_config/CD300LF_panel.py`.
- **GSEA:** CD300LF expression as a continuous phenotype (Pearson metric), MSigDB Hallmark v2026.1, 10,000 gene-set permutations, gene set size 15–500.

## Repository structure

```
paths_config.py          shared paths (raw_data, data_processed, results)
pipeline_common.py       GEO parsing, probe mapping, expression extraction, merging
gene_config/             target gene, marker panel and gene aliases
cohort_runners/          one script per GEO cohort
analysis_specific/       Cox, forest plot, leave-one-out, Kaplan–Meier, correlation
tcga/                    TCGA parsing and GSEA input preparation
requirements.txt         Python dependencies
```

To analyze a different gene, create a new file in `gene_config/` and change `GENE_CONFIG_MODULE` at the top of each script.
