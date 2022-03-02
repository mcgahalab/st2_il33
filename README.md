# ST2 and Il-33 Repo
## Description
ST2 Gene Name: IL1RL1

_IL1RL1_ should be higher expressed in the bulk tumor samples, as it is not expressed in the tumor cells themselves, while IL-33 should be highly expressed in the lymph nodes.

Main focus of this project is to:
  * Look for associations between high ST2 and tumor progression/size/tumor subtypes (ST2 not expressed in tumor cells so it should be fine to look at bulkRNA)
  * Look for correlations between ST2 and 
    *   PD1 expression  - ST2 promotes PD1 expression
    *   CTLA4 expression  - marker similar to PD1
    *   TOX
    *   KLRG1

If possible, find a database on lymph nodes and look for IL-33 expression in senital lymph nodes.
If possible, investigate single cell more and look at levels of ST2 in terms of monocyte activation and deactivation

## Dataset Details
* **GDC TCGA dataset**: https://portal.gdc.cancer.gov/
  * 33 cancer types found across the TCGA project, with the following sample sizes per cancer type:
```TCGA-ACC    79
TCGA-BLCA  433
TCGA-BRCA 1222
TCGA-CESC  309
TCGA-CHOL   45
TCGA-COAD  521
TCGA-DLBC   48
TCGA-ESCA  173
TCGA-GBM   174
TCGA-HNSC  546
TCGA-KICH   89
TCGA-KIRC  611
TCGA-KIRP  321
TCGA-LAML  151
TCGA-LGG   529
TCGA-LIHC  424
TCGA-LUAD  594
TCGA-LUSC  551
TCGA-MESO   86
TCGA-OV    379
TCGA-PAAD  182
TCGA-PCPG  186
TCGA-PRAD  551
TCGA-READ  177
TCGA-SARC  265
TCGA-SKCM  472
TCGA-STAD  407
TCGA-TGCT  156
TCGA-THCA  568
TCGA-THYM  121
TCGA-UCEC  587
TCGA-UCS    56
TCGA-UVM    80
```   

## Workflow
Reproducibility: The conda environment used for this analysis is provided in `envs/environment.yaml`. Packages that had to be installed outside of conda are listed in the `envs/r-env.txt` file.

The dataset directory is set up as followed:
```
data
├── clinical
│   ├── TCGA-ACC.clinical.rds
|   ...
│   └── TCGA-UVM.clinical.rds
├── counts
│   ├── exprs
│   │   ├── TCGA-ACC.HTSeq-Counts.assay.rds
│   │   ...
│   │   └── TCGA-UVM.HTSeq-Counts.assay.rds
│   └── obj
│       ├── TCGA-ACC.HTSeq-Counts.obj.rds
│       ...
│       └── TCGA-UVM.HTSeq-Counts.obj.rds
└── fpkm
    ├── exprs
    │   ├── TCGA-ACC.HTSeq-FPKM-UQ.assay.rds
    |   ...
    │   └── TCGA-UVM.HTSeq-FPKM-UQ.assay.rds
    └── obj
        ├── TCGA-ACC.HTSeq-FPKM-UQ.obj.rds
        ...
        └── TCGA-UVM.HTSeq-FPKM-UQ.obj.rds

```
All these files were generated using the `code/downloadTCGA.R` script in a web-accessible node, then relocated for downstream analysis.


A preliminary analysis of ST2 gene expression across all TCGA projects, as well as a survival analysis and simple ANOVA analysis on the metadata can be found in the `code/st2_expr-survcurv-meta.R` script


 

