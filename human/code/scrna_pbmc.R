library(BiocParallel)
library(GenomicFeatures)
library(harmony)
library(org.Hs.eg.db)
library(chromVAR)
library(ChIPseeker)
library(SingleR)
library(reshape2)
library(enrichplot)
library(Signac)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(plyr)
library(cowplot)
library(ggrastr)
# library("GOstats")
library(ggplot2)
library(GenomicRanges)
library("clusterProfiler")
library("enrichplot")
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
# library(ArchR)

seed <- 1234
set.seed(seed)
args = commandArgs(trailingOnly=TRUE)

sidx <- args[1] # 1  # Where to start preprocessing individual samples from
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_pbmc_giselle'
seurat_dir <- file.path(PDIR, "data", "seurat_obj")
dataset <- 'pgmc'
macs2_path <- '/cluster/home/quever/miniconda3/envs/r4/bin/macs2'
pbmc_annot <- '/cluster/projects/mcgahalab/ref/scrna/pbmc_dataset/pbmc_multimodal.h5seurat'
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/GTF/genome.gtf'
panglao_dbf <- '/cluster/projects/mcgahalab/ref/scrna/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv'

count_blacklist_frags <- FALSE # counting balcklist from 10x fragments is extremely lengthy process

outdir <- file.path(PDIR, "results")
dir.create(outdir, recursive = T, showWarnings = F)
setwd(PDIR)

# Make TxDB object
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

# doublet_barcodes <- file.path(outdir, "doubletfinder", "doublet_barcodes.rds")
# doublet_quantile_cutoff <- 0.95

###################
#### Functions ####
checkObj <- function(seu){
  seu$test <- rep(FALSE, ncol(seu))
  seu$test[1:5] <- TRUE
  subset(seu, subset=test)
}

appendDatasetId <- function(seus, id, datatype){
  i <- seus[[id]]
  seu_p <- i[[datatype]]
  seu_p$dataset <- id
  return(seu_p)
}
 
#####################################################
#### 5.a Cell type annotation - ENCODE/Blueprint ####
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))

dir.create(file.path(outdir, "annotation"))
getBlueprint <- TRUE
if(getBlueprint){
  #ref <- BlueprintEncodeData()
  #saveRDS(ref, file="~/downloads/BlueprintEncodeData.rds")
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/BlueprintEncodeData.rds") 
  for(lbl in c('label.fine', 'label.main')){
    for(clus in c(TRUE, FALSE)){
      rds_file <- paste0("celldex_blueprint.", gsub("label.", "", lbl),
                         ".", if(clus) 'cluster' else 'cell', ".rds")
      id <- gsub("^.*blueprint(.*).rds", "bp\\1", rds_file)
      if(!file.exists(file.path(outdir, "annotation", rds_file))){
        print(rds_file)
        blueprint_anno <- SingleR(test=GetAssayData(seu_integ), ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            clusters=if(clus) seu_integ$seurat_clusters else NULL)
        saveRDS(blueprint_anno, file=file.path(outdir, "annotation", rds_file))
      } else {
        blueprint_anno <- readRDS(file.path(outdir, "annotation", rds_file))
      }
      
      if(clus){
        cluster_ids <- setNames(blueprint_anno$labels, 
                                as.character(rownames(blueprint_anno)))
        seu_integ@meta.data[,id] <- cluster_ids[as.character(seu_integ$seurat_clusters)]
      } else {
        seu_integ@meta.data[,id] <- blueprint_anno$labels
      }
    }
  }
} 


saveRDS(seu_integ, file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))

###########################################
#### 5.b Cell type annotation - PBMC3k ####
# load PBMC reference
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))
checkObj(seu_integ)

reference <- SeuratDisk::LoadH5Seurat(pbmc_annot)
DefaultAssay(seu_integ) <- "SCT"

# transfer cell type labels from reference to query
tr_anchors <- FindTransferAnchors(reference = reference,
                                  query = seu_integ,
                                  normalization.method = "SCT",
                                  reference.reduction = "spca",
                                  recompute.residuals = FALSE,
                                  dims = 1:50)

preds <- TransferData(anchorset = tr_anchors, 
                      refdata = reference$celltype.l2,
                      weight.reduction = seu_integ[['pca']],
                      dims = 1:50)

seu_integ <- AddMetaData(object = seu_integ,
                         metadata = preds)
saveRDS(seu_integ, file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))

################################
#### 6.a Multimodal Dimplot ####
if(!exists("seu"))  seu <- readRDS(file = file.path(seurat_dir, "scIntegrated_anno.rds"))
Idents(seu) <- 'group'
seu <- subset(seu, idents=c('healthy', 'lowISS', 'highISS'))
seu <- seu_grp; rm(seu_grp)

dir.create(file.path(outdir, "dimplot"), showWarnings = F)

DefaultAssay(seu) <- 'SCT'
# clusters
dp_clus <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='seurat_clusters', pt.size=0.5, shuffle=TRUE)

# cell-type annotations
dp_anno <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
               group.by='bp.fine.cell', pt.size=0.5, shuffle=TRUE)
dp_annocl <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='bp.fine.cluster', pt.size=0.5, shuffle=TRUE)
dp_pred <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='predicted.id', pt.size=0.5, shuffle=TRUE)

# samples/groups
dp_lbl <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
               group.by='orig.ident', pt.size=0.5, shuffle=TRUE)
dp_grp <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by='group', pt.size=0.5, shuffle=TRUE)

# split by group
dp_anno_lbl <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by='bp.fine.cluster', split.by='orig.ident', pt.size=0.5, shuffle=TRUE)
dp_anno_grp <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by='bp.fine.cluster', split.by='group', pt.size=0.5, shuffle=TRUE)



pdf(file.path(outdir, "dimplot", "dimplot.pdf"), width=17)
dp_clus + dp_annocl
dp_anno + dp_annocl
dp_clus + dp_lbl
dp_clus + dp_grp

plot_grid(dp_anno_lbl, ncol=1)
plot_grid(dp_anno_grp, ncol=1)
dev.off()

###################################
#### 6.b ST2 dimplot per group ####
DefaultAssay(seu) <- 'SCT'

# ST2, IL33 plot
features <- c('IL1RL1', 'IL33')
featp <- FeaturePlot(seu, features = features)
pdf(file.path(outdir, "dimplot", "featplots.pdf"), width=15)
featp
dev.off()

# Counting ST2/IL33 positive cells (reads > 0)
DefaultAssay(seu) <- 'RNA'
cnts <- GetAssayData(seu)
cell_types <- unique(seu$bp.fine.cluster)
names(cell_types) <- cell_types

cnts_l <- lapply(cell_types, function(ct){
  idx <- which(seu$bp.fine.cluster == ct)
  apply(cnts[features,idx]>0, 1, function(i){
    table(factor(i, levels=c(TRUE,FALSE)))
  })
})

cnts_df <- lapply(setNames(features,features), function(f){
  do.call(rbind, lapply(cnts_l, function(i) i[,f]))
})

# $IL1RL1
#                               TRUE FALSE
# Monocytes                        0 27894
# CD4+ T-cells                     1 19815
# CD8+ Tem                         1  9154
# naive B-cells                    0  7673
# CD4+ Tcm                         2  2935
# CD8+ Tcm                         1  1758
# NK cells                         9  7448
# Class-switched memory B-cells    0   231
# Macrophages                      0   128
# MPP                              3    88

# $IL33
#                               TRUE FALSE
# Monocytes                       16 27878
# CD4+ T-cells                    10 19806
# CD8+ Tem                        15  9140
# naive B-cells                    5  7668
# CD4+ Tcm                         3  2934
# CD8+ Tcm                         4  1755
# NK cells                         9  7448
# Class-switched memory B-cells    0   231
# Macrophages                      0   128
# MPP                              0    91

