renv::load("/cluster/home/quever/downloads/renvs/")
## sara Tumor/LN KO and WT samples
# Core
library(org.Mm.eg.db)
library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(monocle3)
library(BiocParallel)
library(GenomicFeatures)
library(tidyr)
# library(Seurat)
library(dplyr)
library(plyr)
library(reshape2)
library(stringr)
# library(SeuratWrappers)
# library(STACAS)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
# Annotation
library(SingleR)
# library(ProjecTILs)
# QC 
# library(scater)
# library(DoubletFinder)

###############
#### Setup ####

visualize <- FALSE
seed <- 1234
set.seed(seed)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/cd8_ref'
ds1 <- file.path(PDIR, 'data/GSE157072/GSE157072_Processed_D55_LCMV_RNAseq.txt')
ds2 <- file.path(PDIR, 'preprocess_GSE181785/results/counts/all.tsv')
outdir <- file.path(PDIR, "results")
dir.create(outdir, recursive = T, showWarnings = F)
setwd(PDIR)


# Create ENZ -> SYMBOL mapping key
# genome_gse <- org.Mm.eg.db
# txby <- keys(genome_gse, 'SYMBOL')

# Read in gtf file to map ENS->Biotype
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf_file)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

###########################
#### Global Parameters ####
doublet_quantile_cutoff <- 0.95
visualize_qc <- FALSE

mt_miqc <- FALSE # whether to use miQC for percent.mt filter, else use isOutlier
call_doublets <- FALSE # whether to call doublets

# Preprocessed datasets:
#     Must contain a '*counts.csv.gz' file and '*metadata.csv.gz' file
#     in the /cluster/projects/mcgahalab/ref/scrna/projectils/gse_dat/GSEXXXX/
#     directory
pp_datasets <- c('GSE161345') #, 'GSE184423')

# Annotation
atlas_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils'
scgate_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/scgate_model'
projectils_refdir <- file.path(atlas_refdir, "ProjecTILs")

###################
#### Functions ####
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/geneMap.R")
gm <- geneMap("Mus musculus")

###########################
#### 0. Global objects ####
# Reads in a metadata file containing the GSE, GSM, and dataset label
# to be used when annotating the datasets downstream
df1 <- read.table(ds1, sep="\t", header = T, stringsAsFactors = F) %>%
  dplyr::select(c("NAME", grep("^(IEL|Sp)", colnames(.), value=T))) %>%
  magrittr::set_colnames(gsub(".bam$", "", colnames(.))) %>%
  mutate(NAME=make.unique(gsub(" ", "", NAME)))
df2 <- read.table(ds2, sep="\t", header = T, stringsAsFactors = F) %>%
  mutate("NAME"=make.unique(gm$ENSEMBL$SYMBOL[gene])) %>%
  select(-gene) %>% #tibble::column_to_rownames("gene") %>%
  magrittr::set_colnames(gsub("^GSM.*?_", "", colnames(.)))

df <- full_join(df1, df2, by='NAME')
saveRDS(df, file=file.path(outdir, "counts.rds"))

################################
#### 1.a Make DESeq2 object ####
# This section of the code reads in all the 10x cellranger
# files and loads them into seurat objects. It saves an
# RDS object which contains a list of seurat objects, one
# for each dataset
df <- readRDS(file=file.path(outdir, "counts.rds")) 
df2 <- df %>% 
  filter(!grepl("^NA\\.", NAME)) %>%
  filter(!grepl("^Gm[0-9]*", NAME)) %>%
  filter(!grepl("[0-9]*Rik$", NAME)) %>%
  filter(!grepl("^Rpl[0-9]*", NAME)) %>%
  filter(!grepl("^Rps[0-9]*", NAME)) %>%
  filter(!grepl("^Mir[0-9]*", NAME)) %>%
  filter(!grepl("^Or[0-9]*", NAME)) %>%
  filter(!grepl("^AT[0-9]*", NAME)) %>%
  filter(!grepl("^CT[0-9]*", NAME)) %>%
  filter(!grepl("^AC[0-9]*", NAME)) %>%
  filter(!is.na(NAME)) %>%
  tibble::column_to_rownames("NAME")
df2[is.na(df2)] <- 0
df2 <- as.matrix(df2)
storage.mode(df2) <- 'integer'

coldata <- data.frame(batch=c(rep('GSE157072', '9'), rep('GSE18178', 3)),
                      condition=gsub("(\\.[12]|_30)$", "", colnames(df2)) %>%
                        gsub("(IEL.|Sp.|CD8_)", "", .) %>% tolower) 
rownames(coldata) <- colnames(df2)
dds <- DESeqDataSetFromMatrix(countData = df2,
                              colData = coldata,
                              design= ~ batch + condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
write.table(assay(vsd), file=file.path(outdir, "vst_counts_merged.csv"),
            sep=",", col.names = T, row.names = T, quote = F)

# resultsNames(dds) # lists the coefficients
# res <- results(dds, name="condition_trm_vs_conv.tem")

