# Pan-Glioma IDH and 1p/19q Codeletion Survival Stratification (TCGA)
## Objective
To assess how molecular subtype (IDH status and 1p/19q codeletion) stratifies overall survival across adult diffuse glioma, combining low-grade glioma (LGG) and glioblastoma (GBM) cohorts from The Cancer Genome Atlas (TCGA).

## Data
Source: cBioPortal, TCGA PanCancer Atlas 2018 clinical patient files for LGG (lgg_tcga_pan_can_atlas_2018) and GBM (gbm_tcga_pan_can_atlas_2018).
- Only patients in the PanCancer Pathways freeze (IN_PANCANPATHWAYS_FREEZE == Yes) were kept.
- Patients missing overall survival time or status were excluded.
- The data are not included in this repository. Download the two cohorts from cBioPortal and set the file paths at the top of Code.py.

## Methodology
- Cleaning: removed administrative and unused columns; converted survival status fields (OS, DSS, PFS) to binary event indicators (1 = event, 0 = censored).
- Subtype assignment: extracted molecular subtype from the TCGA SUBTYPE field into three groups: IDH-wildtype (IDHwt), IDH-mutant with 1p/19q codeletion (IDHmut-codel), and IDH-mutant without codeletion (IDHmut-non-codel). GBM samples labelled GBM were assigned IDHwt.
- Cohort: LGG and GBM cohorts were concatenated into one pan-glioma dataset.
- Survival analysis (lifelines):
Kaplan-Meier estimation of overall survival
Pairwise log-rank tests between molecular subtype groups
Cox proportional hazards regression of overall survival on age and IDH/1p19q group (IDHwt as reference)
