## sara Tumor/LN KO and WT samples
# Core
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(monocle3)
library(BiocParallel)
library(GenomicFeatures)
library(tidyr)
library(Seurat)
library(dplyr)
library(plyr)
library(reshape2)
library(stringr)
library(SeuratWrappers)
library(STACAS)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
# Annotation
library(SingleR)
library(ProjecTILs)
# QC 
library(scater)
library(DoubletFinder)
#PacMAP
library(reticulate)

################
#### Python ####
# Python functions for reticulate
use_python("/cluster/home/quever/miniconda3/bin/python3.9")
python_pacmap <- import("pacmap")
python_pandas <- import("pandas")
python_numpy <- import("numpy")
###############
#### Setup ####

visualize <- FALSE
seed <- 1234
set.seed(seed)

ref_seu_dir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcell_ref/GSE124691'
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcell_ref/GSE124691'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
dir.create(datadir, recursive = T, showWarnings = F)
dir.create(outdir, recursive = T, showWarnings = F)
setwd(PDIR)


# Create ENZ -> SYMBOL mapping key
# genome_gse <- org.Mm.eg.db
# txby <- keys(genome_gse, 'SYMBOL')

# # Read in gtf file to map ENS->Biotype
# gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
# GTF <- rtracklayer::import(gtf_file)
# ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

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
runPacmap <- function(seu, reduction){
  # PacMAP clustering
  dat_pd <- reticulate::r_to_py(as.data.frame(seu@reductions[[reduction]]@cell.embeddings))
  nparray <- dat_pd$values
  nparray <- nparray$astype(python_numpy$float)
  tfmr <- python_pacmap$PaCMAP(n_neighbors=10L, MN_ratio=0.5, FP_ratio=2.0) 
  X_transformed <- tfmr$fit_transform(nparray)#, init="pca")
  seu[["pacmap"]] <- CreateDimReducObject(embeddings = data.frame(X_transformed) %>%
                                            magrittr::set_rownames(., Cells(seu)) %>%
                                            as.matrix, 
                                          key = "PacMAP_", 
                                          assay = DefaultAssay(seu))
  return(seu)
}
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/STACAS_fix.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")

###########################
#### 0. Global objects ####
# Reads in a metadata file containing the GSE, GSM, and dataset label
# to be used when annotating the datasets downstream


#######################################
#### 1.a Make ProjecTILs reference ####
# This section of the code reads in all the 10x cellranger
# files and loads them into seurat objects. It saves an
# RDS object which contains a list of seurat objects, one
# for each dataset
ref_groups <- list.files(ref_seu_dir, pattern = "^GSM")
outrds <- file.path(datadir, "seurat_obj", "ref1_unfiltered.rds")
relabel_id <- gsub("^GSM[0-9]*_", "", ref_groups) %>% 
  gsub("Experiment_", "rep", .) %>%
  setNames(., ref_groups)

if(file.exists(outrds)){
  cat("Reading in existing 'seus' object...\n")
  seus <- readRDS(file=outrds)
} else {
  seus <- NULL
}

seus <- lapply(ref_groups, function(grp){
  if(relabel_id[grp] %in% names(seus)){
    cat(paste0("Existing ID: ", grp, "\n"))
    seu <- seus[[relabel_id[grp]]]
  } else {
    gsm_id <- gsub("_.*", "", grp)
    print(paste0(grp, "..."))
    seu <- tryCatch({
      mtx <- Read10X(data.dir = file.path(ref_seu_dir, grp), 
                     strip.suffix=TRUE)
      if(class(mtx)=='list') mtx <- mtx[['Gene Expression']]
      CreateSeuratObject(counts = mtx, project = gsm_id)
    }, error=function(e){cat("Failed\n"); NULL})
  }
  return(seu)
})
names(seus) <- relabel_id
seus <- seus[which(!sapply(seus, is.null))]

# For each sampel, drop all cells where there are fewer than 50 genes detected
sapply(seus,dim) %>% as.data.frame.matrix %>% t
seus <- lapply(seus, function(seu_i){
  subset(seu_i, subset = nFeature_RNA > 50)
})
sapply(seus,dim) %>% as.data.frame.matrix %>% t

dir.create(dirname(outrds), recursive = T, showWarnings = F)
saveRDS(seus, file=outrds)

###########################################
#### 1. Preprocess to remove doublets  ####
if(!exists("seus")) seus <- readRDS(file.path(datadir, "seurat_obj", "ref1_unfiltered.rds"))
outrds_int <- file.path(datadir, "seurat_obj", "ref2_sct.rds")

### Preprocess and Doublet Calling
names(seus) <- as.character(sapply(seus, function(i) unique(i$orig.ident)))
seus <- lapply(seus, function(seu_i){
  PercentageFeatureSet(seu_i, pattern = "^mt-", col.name = "percent.mt")
})

#.... 1) Identify mitochondrial/nCount/nFeature outliers: ----
outliers <- lapply(seus, function(seu_i){
  # 3-MAD/Mindiff from Median as outliers
  # mitochondrial outlier: scater::isOutlier
  mt_outlier_scuttle <- isOutlier(seu_i$percent.mt, nmads=3, 
                                  min.diff=(10-median(seu_i$percent.mt)), 
                                  batch=seu_i$orig.ident, share.mads=T, 
                                  share.medians=T, type="higher")
  # nCount outlier: scater::isOutlier
  ncnt_outlier_scuttle <- isOutlier(seu_i$nCount_RNA, nmads=4, 
                                    # min.diff=median(seu_i$nCount_RNA), 
                                    batch=seu_i$orig.ident, share.mads=F, 
                                    share.medians=F, type="higher")
  # nFeature outlier: scater::isOutlier
  nft_outlier_scuttle <- isOutlier(seu_i$nFeature_RNA, nmads=3, 
                                   # min.diff=median(seu$nFeature_RNA), 
                                   batch=seu_i$orig.ident, share.mads=F, 
                                   share.medians=F, type="higher")
  
  nmt.max <- min(seu_i$percent.mt[mt_outlier_scuttle]) # 6530
  nfeature.max <- min(seu_i$nFeature_RNA[nft_outlier_scuttle]) # 6530
  ncnt.max <- min(seu_i$nCount_RNA[ncnt_outlier_scuttle]) # 37252
  c('mt'=max(c(nmt.max, 10)), 
    'nft'=max(c(nfeature.max, 5000)), 
    'ncnt'=max(c(ncnt.max, 17000)))
})

#.... 2) Flag Potential Doublets ----
#a) Preprocess for Doublet Calling
s.features <- stringr::str_to_title(cc.genes$s.genes)
g2m.features <- stringr::str_to_title(cc.genes$g2m.genes)

seu_ids <- names(seus)
seus_post <- lapply(seu_ids, function(seu_id){
  cat(paste0("Preprocessing and adding on ", seu_id, "\n"))
  seu_i <- seus[[seu_id]]
  
  gsm_id <- as.character(unique(seu_i$orig.ident))
  outlier_i <- outliers[[gsm_id]]
  seu_i <- preprocessSeu(seu_i, ncount_min=1000, ncount_max=Inf, 
                         nfeature_min=50, nfeature_max=outlier_i['nft'], 
                         mt_max=outlier_i['mt'], org='mouse', numpcs=30, getPCs=FALSE)
  return(seu_i)
})
names(seus_post) <- seu_ids

saveRDS(seus_post, file=outrds_int)

######################################
#### 2. Identify Treg populations ####
if(!exists("seus_post")) seus_post <- readRDS(file.path(datadir, "seurat_obj", "ref2_sct.rds"))
seu.mnn <- RunFastMNN(object.list =seus_post, 
                      features = 3000, assay='SCT') %>%
  RunUMAP(., reduction = "mnn", dims = 1:30, n.neighbors=30L,
          min.dist = 0.1, return.model=TRUE) %>%
  FindNeighbors(., reduction = "mnn", dims = 1:30) %>%
  FindClusters(., resolution = 0.9, graph.name='SCT_snn') %>%
seu.mnn <- seu.mnn %>%   runPacmap(., reduction='mnn')

seu.stacas <- lapply(seus_post, function(seu_i){
  DefaultAssay(seu_i) <- 'RNA'
  return(seu_i)
}) %>%
  Run.STACAS(., future.maxSize=25) %>%
  RunUMAP(., dims = 1:30, n.neighbors=30L, min.dist = 0.1) %>%
  FindNeighbors(., dims = 1:30) %>%
  FindClusters(., resolution = 1.2)
seu.stacas <- seu.stacas %>%   runPacmap(., reduction='pca')

if(do_visualize){
  DefaultAssay(seu.mnn) <- 'RNA'
  VariableFeatures(seu.mnn) <- rownames(seu.mnn)
  seu.mnn <- ScaleData(seu.mnn)
  DefaultAssay(seu.stacas) <- 'RNA'
  VariableFeatures(seu.stacas) <- rownames(seu.stacas)
  seu.stacas <- ScaleData(seu.stacas)
  genes <- c('Foxp3', 'Cxcr3', 'Cxcr6', 'Ccr7', 'Tcf7', 'Cxcr5')
  
  pdf("~/xfer/x.pdf", width = 18)
  dp1 <- DimPlot(seu.mnn, group.by='seurat_clusters', raster=T, label=T)
  fp1 <- FeaturePlot(seu.mnn, features=genes, raster=T, pt.size=4)
  cowplot::plot_grid(dp1, fp1, rel_widths = c(1,3))
  
  dp2 <- DimPlot(seu.stacas, group.by='seurat_clusters', raster=T, label=T)
  fp2 <- FeaturePlot(seu.stacas, features=genes, raster=T, pt.size=4)
  cowplot::plot_grid(dp2, fp2, rel_widths = c(1,3))
  
  dp2 <- DimPlot(seu.stacas, group.by='seurat_clusters', raster=T, label=T, reduction='pacmap')
  fp2 <- FeaturePlot(seu.stacas, features=genes, raster=T, pt.size=4, reduction='pacmap')
  cowplot::plot_grid(dp2, fp2, rel_widths = c(1,3))
  dev.off()
}


recode_map <- function(x, method=NULL){
  x %>% 
    recode("0"=ifelse(method=='stacas', "V.Tfh", "II"),
           "1"=ifelse(method=='stacas', "II.Isc_nRes", "IV"),
           "2"=ifelse(method=='stacas', "III.Treg", "V"),
           "3"=ifelse(method=='stacas', "II.Isc_nRes", 'III'),
           "4"=ifelse(method=='stacas', "IV.Ccr7neg", "V"),
           "5"=ifelse(method=='stacas', "IV.Ccr7neg", "IV"),
           "6"=ifelse(method=='stacas', "I.Th1", "I"),
           "7"=ifelse(method=='stacas', "Unk", 'III'),
           "8"=ifelse(method=='stacas', "V.Tfh", "NA"),
           "9"=ifelse(method=='stacas', "Unk", "NA"),
           "10"=ifelse(method=='stacas', "III.Treg", "NA"))
}


seu.stacas$anno <- seu.stacas$seurat_clusters %>% recode_map(., method='stacas')
seu <- seu.stacas
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 1500, verbose = FALSE)
ndim = 20
seed = 1234
seu <- ScaleData(seu, verbose = TRUE)

seu <- RunPCA(seu, features = VariableFeatures(seu), ndims.print = 1:5,
                 nfeatures.print = 5, npcs = ndim)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:ndim, seed.use = seed)

ref.tcell <- ProjecTILs::make.reference(ref = seu, ndim = ndim, seed = seed, 
                                       recalculate.umap = TRUE, annotation.column = "anno")


saveRDS(ref.tcell, file=file.path(projectils_refdir, "seu.sara_tcells.rds"))
