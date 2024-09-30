# Lymph node macrophages drive immune tolerance and resistance to cancer therapy by expression of the immune-regulatory cytokine IL-33.
## Bioinformatic Analysis

### 1. scRNA Analysis
#### A) Differential and pathway analysis
OVERVIEW:
This section performs all the analysis on the ST2-sorted and CD45-sorted single-cell RNAseq data. The annotations for each of theese datasets can be found in the data directory:
- CD45_Tumor_cd8.csv
- ST2_LN_Tregs.csv
- ST2_LN_allcells.csv
- ST2_Tumor_Tregs.csv
- ST2_Tumor_allcells.csv

TReg reference dataset from Gomes et al., 2019 (DOI: 10.1016/j.immuni.2019.01.001) was downloaded from https://figshare.com/projects/Treg_scRNA-seq/38864.

#### B) Velocity/trajectory analysis
Differential analysis comparing local velocities in the scVelo analysis can be performed using the `scvelo_analysis.py` script using the `scvelo_helper` python package.


### 2. RNA-seq Analysis
#### A) MSM-SSM Analysis
OVERVIEW:
This section performs all the analysis on the TDLN and LN MSM/SSM samples that were treated with Cisplating (Cis) or PBS. The main analysis is performed in the `msm_ssm_analysis.R` code. The relevant count/tpm and metadata needed is found in the `data/msm_ssm.{counts|samples|tpm}.tsv` files. Everything from preprocessing to the differential expression, GSEA, ssGSEA, and regulon analysis can be found in this section.

#### B) rDC-mDC analysis
OVERVIEW:
This section performs all the analysis on the regulatory and migratory DCs in the Tumor/WT samples that were treated with Cisplating (Cis) or PBS. The main analysis is performed in the `mdc_rdc_analysis.R` code. The relevant count/tpm and metadata needed is found in the `data/mdc_rdc.{counts|samples|tpm}.tsv` files. Everything from preprocessing to the differential expression, GSEA, ssGSEA, and regulon analysis can be found in this section.

### 3. Meta-Analysis
OVERVIEW:
This section performs the survival curve analysis on the TCGA and GTEx dataset, as well as the expanded analysis on the 3 RNA-seq dataasets from Liu, Gide, and Riaz.

REQUIREMENTS:
The TReg signature derived in the earlier section and Supp Table 2 will be used here, stored as `treg_signature_paper.human.csv`

Liu_PMID31792460:
- TPM matrix from Supp. Data 2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6898788/bin/41591_2019_654_MOESM3_ESM.txt
- Metadata from Supp. Table 1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6898788/bin/41591_2019_654_MOESM4_ESM.xlsx

Gide_PMID30753825:
- Count matrix from Git repo: https://github.com/miabioinformatics/Gide_Quek_CancerCell2019
- Metadata from Supp. Table 1 and 2: https://www.sciencedirect.com/science/article/pii/S1535610819300376?via%3Dihub#app2

Riaz_PMID29033130:
- FPKM matrix from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
- Metadata from Supp. Table 2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5685550/bin/NIHMS907788-supplement-11.xlsx

TCGA and GTEX:
- Gene expression FPKM matrix downloaded from Xena: https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx
- Phenotype metadata downloaded from Xena: https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx
