renv::load("/cluster/home/quever/downloads/renvs/")
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
# Preprocess
library(STACAS)
library(SeuratWrappers)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
library(scCustomize)
# Annotation
library(SingleR)
library(ProjecTILs)
# library(cellpypes)
# Enrichment Analysis
library(msigdbr)
library(SCPA)
library(AUCell)
library(GSEABase)
library(enrichplot)
library(pagoda2)
# Regulon Analysis
library(SCopeLoomR)
library(SCENIC)
# QC 
library(scater)
library(DoubletFinder)
# Trajcectory/Velocity
library(parsnip)
# library(velocyto.R)
# library(slingshot)
# library(tradeSeq)
library(reticulate)
# library(harmony)
# # library("GOstats")
# library("clusterProfiler")

# BiocManager::install("AUCell")
# BiocManager::install("tradeSeq")
# install.packages("cellpypes")
# devtools::install_github("dynverse/dyno")
#  conda install -c conda-forge r-uwot
# install.packages(c("Matrix", "sctransform", "SeuratObject"))
# devtools::install_github("satijalab/seurat", ref = "develop")
# install.packages("scSorter")

visualize <- FALSE
seed <- 1234
set.seed(seed)
args <- commandArgs(trailingOnly = TRUE)
group <- args[1]

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d'
setwd(PDIR)

gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
groups <- list.files(datadir, pattern="[^seurat_obj|velocyto|ids|metrics_summary.csv|annotation_map]")

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))
sym2biotype_ids <- with(GTF, setNames(gene_biotype, gene_name))


doublet_quantile_cutoff <- 0.95

###################
#### Functions ####
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
# source("~/git/mini_projects/mini_functions/singlecell/publicationDimPlot.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")
gm <- geneMap(species='Mus musculus')

##-- Recode functions ----
## recode_clusters.R
# recode_list()
# recode_map()
# cd45_recode_map()
# treg_recode_map()
# cd45_treg_recode_map()
source(file.path(PDIR, "ref", "recode_clusters.R"))

# Python functions for reticulate
use_python("/cluster/home/quever/miniconda3/bin/python3.9")
python_pacmap <- import("pacmap")
python_pandas <- import("pandas")
python_numpy <- import("numpy")

#### Other Functions ####
aucellFun <- function(msig_ds, expr_mat, gm, mapfrom='SYMBOL', mapto='ENTREZID'){
  mapped_id <- if(mapfrom==mapto){
    msig_ds$entrez_gene
  } else {
    gm[[mapfrom]][[mapto]][as.character(msig_ds$entrez_gene)]
  }
  
  msig_l <- lapply(split(mapped_id, msig_ds$gs_name), function(i){
    i[!is.na(i)]
  })
  auc <- tryCatch({
    AUCell_run(expr_mat, msig_l)
  }, error=function(e){NULL})
  return(auc)
}

# PacMAP clustering
runPacmap <- function(seu, reduction='pca'){
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

.monocle3Cluster <- function(seu_i, ...){
  set.seed(1234)
  data <- as(as.matrix(seu_i@assays$mnn.reconstructed@data), 'sparseMatrix')
  pd <- seu_i@meta.data
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  cds <- new_cell_data_set(expression_data=data,
                           cell_metadata  = seu_i@meta.data,
                           gene_metadata = fData)
  cds <- preprocess_cds(cds, num_dim = 30, norm_method='size_only', pseudo_count=0)
  reducedDims(cds)$UMAP <- seu_i@reductions$umap@cell.embeddings
  # reducedDims(cds)$PCA <- seu_i@reductions$pca@cell.embeddings[,1:30]
  cds <- cluster_cells(cds, ...)
  seu_i$monocle3_clusters <- cds@clusters$UMAP$clusters
  return(seu_i)
}

# Remove TReg cells that are scattered in other clusters
cleanLnTumorTreg <- function(seul){
  for(seuid in c('LN', 'Tumor')){
    seu_i = seul[[seuid]]
    seu_i <- .monocle3Cluster(seu_i)
    
    seu_i$treg <- grepl("TReg", seu_i$manual_anno)
    pdf("~/xfer/x.pdf")
    DimPlot(seu_i, reduction='umap', label=T, group.by='monocle3_clusters', raster=T)
    dev.off()
    if(seuid == 'Tumor'){
      Idents(seu_i) <- 'monocle3_clusters'
      seu_i <- subset(seu_i, ident=setdiff(unique(Idents(seu_i)), '16'))
      idx <- grepl("^TReg", seu_i$manual_anno) & (seu_i$monocle3_clusters %in% as.character(c(13, 5, 19, 20, 14, 10, 9, 15)))
    } else if(seuid == 'LN'){
      idx <- grepl("^TReg", seu_i$manual_anno) & (seu_i$monocle3_clusters %in% as.character(c(8,9,4,17,12,15,16,11,2,3,14)))
    }
    seu_i$treg_mispositioned <- idx
    seu_i <- subset(seu_i, cells=setdiff(Cells(seu_i), Cells(seu_i)[which(idx)]))
    seul[[seuid]] <- seu_i
  }
  return(seul)
}

####################################
#### 0.a Create Seurat objects  ####
dir.create(file.path(datadir, "seurat_obj"), 
           recursive = TRUE, showWarnings = FALSE)

seus <- lapply(groups, function(grp){
  print(paste0(grp, "..."))
  mtx <- Read10X(data.dir = file.path(datadir, grp), strip.suffix=TRUE)
  seu <- CreateSeuratObject(counts = mtx, project = grp)
  return(seu)
})
names(seus) <- groups

## Merge the scRNA data and save
seu <- merge(x=seus[[1]], y = seus[-1], 
             add.cell.ids = names(seus), 
             project = 'st2_il33_7d')
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu.rds"))
rm(seus)

###########################
#### 0.b QC - RNA data ####
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu.rds"))
dir.create(file.path(outdir, "qc"), recursive = T, showWarnings = F)

DefaultAssay(seu) <- "RNA"
seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mt")
# Plot percent.mt per sample
if(visualize){
  pdf(file.path(outdir, "qc", "percent_mt.pdf"))
  VlnPlot(object = seu, features = 'percent.mt', split.by = 'orig.ident')
  dev.off()
  
  # Remove high mt cells
  mt_fracs <- sapply(split(seu$percent.mt, f=seu$orig.ident), quantile, 
                     probs=seq(0.9, 0.99, by=0.01)) %>%
    round(., 2) %>% 
    as.data.frame() %>%
    mutate(median=rowMedians(as.matrix(.)))
}

# Remove high mitochondrial read samples
mt_cutoff <- 10
seu_mt <- seu[,which(seu$percent.mt <= mt_cutoff)]
seu_qcsumm <- cbind(table(seu$orig.ident),
                    table(seu_mt$orig.ident)) %>%
  as.data.frame() %>%
  rename_with(., ~ c("original", "low.mt")) %>%
  mutate("low.mt.frac"=round(low.mt / original, 2))


if(visualize){
  pdf(file.path(outdir, "qc", "featureScatter.pdf"), width = 14, height = 7)
  plot1 <- FeatureScatter(seu_mt, feature1 = "nCount_RNA", 
                          feature2 = "percent.mt", shuffle = TRUE, raster=T)
  plot2 <- FeatureScatter(seu_mt, feature1 = "nFeature_RNA", 
                          feature2 = "nCount_RNA", shuffle = TRUE, raster=T)
  plot_grid(plot1, plot2, nrow=1)
  dev.off()
  
  # Identify count/feature outliers cells
  feat_fracs <- sapply(split(seu$nFeature_RNA, f=seu$orig.ident), quantile, 
                       probs=c(seq(0, 0.1, by=0.01),
                               seq(0.9, 1, by=0.01))) %>%
    round(., 0) %>% 
    as.data.frame() %>%
    mutate(median=rowMedians(as.matrix(.)))
  mfeat_fracs <- melt(t(feat_fracs)) %>% 
    mutate(Var2=as.integer(gsub("%", "", Var2)))
  count_fracs <- sapply(split(seu$nCount_RNA, f=seu$orig.ident), quantile, 
                        probs=c(seq(0, 0.1, by=0.01),
                                seq(0.9, 1, by=0.01))) %>%
    round(., 0) %>% 
    as.data.frame() %>%
    mutate(median=rowMedians(as.matrix(.)))
  mcount_fracs <- melt(t(count_fracs)) %>% 
    mutate(Var2=as.integer(gsub("%", "", Var2)))
  
  pdf(file.path(outdir, "qc", "count_feat_extremes.pdf"), width = 14, height = 7)
  gp1 <- ggplot(mfeat_fracs, aes(x=Var2, y=value, fill=Var1, col=Var1))+
    geom_point(alpha=0.8) + #scale_color_manual(values=c("median"="black")) +
    theme_bw()
  gp2 <- ggplot(mcount_fracs, aes(x=Var2, y=value, fill=Var1, col=Var1))+
    geom_point(alpha=0.8) + #scale_color_manual(values=c("median"="black")) +
    theme_bw()  
  plot_grid(gp1, gp2, nrow=1)
  dev.off()
  
}

## https://www.nature.com/articles/s41591-021-01323-8#Sec15
# All cells expressing <200 or >6,000 genes were removed
# cells that contained <400 unique molecular identifiers (UMIs) 
# >15% mitochondrial counts

# Remove outlier samples for count and features
seu_qc <- subset(seu_mt,
                 subset =nCount_RNA > 1000 &
                   nFeature_RNA < 6686)
seu_qcsumm2 <- cbind(seu_qcsumm,
                     table(seu_qc$orig.ident) %>%
                       as.integer) %>%
  as.data.frame() %>%
  rename_with(., ~ c("original", "low.mt", "low.mt.frac", "count.filt")) %>%
  mutate("count.frac"=round(count.filt / original, 2))

seu <- seu_qc
saveRDS(seu, file = file.path(datadir, "seurat_obj", "seu_filt.rds"))
rm(seu_qc, seu_mt)

#############################################################
#### 0.c Doublet Finding on a per-sample basis: RUN ONCE ####
doublet_dir <- file.path(outdir, "doublets")
dir.create(doublet_dir, recursive = T, showWarnings = F)
doublet_barcodes <- file.path(doublet_dir, "seu_doublets.rds")
if(file.exists(doublet_barcodes)) warning("Doublet files already identified")

if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu_filt.rds"))

## Apply DoubletFinder to every scRNA sample independently
samples <- unique(seu$orig.ident)
Idents(seu) <- 'orig.ident'
gen_sct <- TRUE
get_pk <- FALSE
lbl <- 'label.main' #label.main, 'label.fine
rds_dir <- if(gen_sct) 'sct_rds' else 'scale_rds'
dir.create( file.path(doublet_dir, rds_dir), showWarnings = F)
seu_doublet_list <- lapply(setNames(samples,samples), function(sample_i){
  seu_i <- subset(seu, idents=sample_i)
  
  doublet_rds_f <- file.path(doublet_dir, rds_dir, paste0(sample_i, ".rds"))
  
  # Process doublets or read in existing files
  if(file.exists(doublet_rds_f)){
    print(paste0("Loading existing analysis: ", sample_i))
    if(file.info(doublet_rds_f)$size>0){
      doublet <- readRDS(doublet_rds_f)
    } else {
      print(paste0(" ... ", sample_i, " being processed in other node"))
      doublet <- NULL
    }
    return(doublet)
  } else {
    print(paste0("DoubletFinding for sample: ", sample_i))
    write.table(data.frame(), file=doublet_rds_f, col.names=FALSE)
  }
  
  # pre-process
  if(gen_sct){
    seu_i <- SCTransform(seu_i)
  } else{ 
    seu_i <- NormalizeData(seu_i, normalization.method = "LogNormalize")
    seu_i <- FindVariableFeatures(seu_i, selection.method = "vst", nfeatures = 2000)
    seu_i <- ScaleData(object = seu_i, vars.to.regress = c("percent.mt"))
  }
  seu_i <- RunPCA(seu_i, npcs=200)
  
  dvar <- seu_i@reductions$pca@stdev^2
  PCNum <- c(min(which(cumsum(dvar)/sum(dvar) >= 0.99)) + 1,
             min(which(cumsum(dvar)/sum(dvar) >= 0.9)))
  PCNum <- min(PCNum)
  
  seu_i <- FindNeighbors(object = seu_i, dims = 1:PCNum)
  res <- getResForStableSmallClus(seu_i)
  if(length(res) == 0) res <- 1.2
  seu_i <- FindClusters(object = seu_i, resolution = res)
  seu_i <- RunUMAP(seu_i, dims = 1:PCNum)
  
  # cell-type identification
  print("Identifying cell type")
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
  rownames(bed.se) <- str_to_title(rownames(bed.se))
  
  if(gen_sct) DefaultAssay(seu_i) <- 'RNA'
  seu_i <- NormalizeData(seu_i, normalization.method = "LogNormalize")
  blueprint_anno <- SingleR(test=GetAssayData(seu_i, assay = "RNA", slot = "data"), 
                            ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            de.method="wilcox", genes='sd',
                            clusters=seu_i$seurat_clusters)
  cluster_ids <- setNames(blueprint_anno$pruned.labels, 
                          as.character(rownames(blueprint_anno)))
  cluster_ids[is.na(cluster_ids)] <- 'Unknown'
  seu_i@meta.data[,'immgen_fine_cell'] <- cluster_ids[as.character(seu_i$seurat_clusters)]
  
  # pK identification
  if(gen_sct) DefaultAssay(seu_i) <- 'SCT'
  sweep_res_list <- NULL
  bcmvn <- NULL
  k=100
  pk <- round(k/ncol(seu_i),5)
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(seu_i$immgen_fine_cell)
  nExp_poi <- getExpMultiplet(ncol(seu_i))  # multiplet chart from https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
  nExp_poi.adj <- nExp_poi*(1-homotypic.prop)
  
  seu_i <- doubletFinder_v4(seu_i, PCs = 1:PCNum, pN = 0.25, pK = pk,  # first bump under 0.1 from bcmvn_seu
                            nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = gen_sct,
                            annotations=seu_i$immgen_fine_cell, grp.size=50)
  seu_i$doublets <- seu_i@meta.data[,'DF.classifications']
  bcds <- runBcds(seu_i, n=25)
  seu_i$bcds <- bcds$bcds[,1]
  cxds <- runCxds(seu_i, n=25)  
  seu_i$cxds <- cxds$cxds$cxds
  
  seul <- list("seu"=seu_i, "sweep"=sweep_res_list,  "res"=res,
               "homotypic.prop"=homotypic.prop, "nExp_poi.adj"=nExp_poi.adj,
               "pk"=pk, "nExp_poi"=nExp_poi, "bcmvn"=bcmvn,
               "bcds"=bcds$top_pairs, "cxds"=cxds$top_pairs)
  print("Saving...")
  saveRDS(seul, file=doublet_rds_f)
  return(seul)
})



{
  # Set celltype colors
  cell_cols <- getCellCols()
  cells <- lapply(seu_doublet_list, function(i) unique(i$seu$immgen_fine_cell)) %>%
    unlist %>% unique %>% sort
  cell_cols_lbl <- lapply(names(cell_cols), function(cl){
    cls <- grep(paste0("^", cl), cells, value=T)
    setNames(cell_cols[[cl]](length(cls)), cls)
  }) %>% unlist
  
  
  # Contour maps
  df_contours <- lapply(seu_doublet_list, function(seu_i){
    doubletFinder_contour(seu_i$seu) %>% #, rm.homotypic=FALSE, max.break=0.5) %>%
      cowplot::plot_grid(plotlist=., ncol=4)
  })
  
  # Heatmap of ES per cluster
  df_heatmap <-  lapply(seu_doublet_list, function(seu_i){
    doubletFinder_heatmap(seu_i$seu, cell_cols_lbl)
  })
  
  # Manual FeaturePlot review
  Bcells <- c('CD79a', 'CD79b', 'CD19', 'EBF1' ) # cd45r
  Tcells <- c('CD3D', 'CD3G', 'CD3e', 'CD4')
  p <- stringr::str_to_title(unlist(c(Bcells, Tcells)))
  df_featplots <- lapply(seu_doublet_list, function(seu_i){
    FeaturePlot(seu_i$seu, features=p, pt.size=0.5, order=T, ncol=4)
  })
  
  
  # SCDS gene-coexpression and binary doublet classifier
  scds_scores <- c("bcds", "cxds")
  df_scdsplots <- lapply(seu_doublet_list, function(seu_i){
    FeaturePlot(seu_i$seu, features=scds_scores, 
                pt.size=0.5, order=T, blend=T)
  })
  
  # Plot everything
  pdf("~/xfer/test13_feat.pdf", width = 15, height = 7)
  df_featplots
  dev.off()
  
  pdf("~/xfer/test13_con.pdf", width = 12, height = 10)
  df_contours
  dev.off()
  
  pdf("~/xfer/test13_hm.pdf", width = 14)
  df_heatmap
  dev.off()
  
  pdf("~/xfer/test13_scds.pdf", height = 5)
  df_scdsplots
  dev.off()
  
}

doublets <- sapply(seu_doublet_list, function(seu_i){
  seu_i <- seu_i$seu
  ndoublets <- sum(seu_i$doublets == 'Doublet')
  dub_idx <- seu_i@meta.data[,'doublets'] == 'Doublet'
  colnames(seu_i)[which(dub_idx)]
}) %>% 
  unlist() %>% as.character()

seu$doublets <- 'Singlet'
seu@meta.data[which(colnames(seu) %in% doublets),'doublets'] <- 'Doublet'
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_filt_dubrm.rds"))


#########################################
#### 1.a Integrate the scRNA samples ####
# Preprocessing of the RNA counts
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu_filt_dubrm.rds"))

# split the dataset into a list of seurat objects by Sample ID
seu.list <- SplitObject(seu, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
preproc_f <- file.path(datadir, "seurat_obj", "seu_separate.sct.rds")
if(!file.exists(preproc_f)){
  seu.list <- lapply(seu.list, function(seu_i){
    # https://github.com/satijalab/seurat/issues/1739
    seu_i <- SCTransform(seu_i, assay = 'RNA', new.assay.name = 'SCT', vst.flavor='v2',
                         vars.to.regress = c('percent.mt'), conserve.memory = T) %>%
      RunPCA(npcs = 30, verbose = FALSE) %>%
      RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
      FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
      FindClusters(resolution = 0.7, verbose = FALSE)
    # seu_i <- CellCycleScoring(seu_i, s.features = str_to_title(cc.genes$s.genes), 
    #                           g2m.features = str_to_title(cc.genes$g2m.genes),
    #                           assay='SCT', set.ident = TRUE)
    # seu_i$CC.Difference <- seu_i$S.Score - seu_i$G2M.Score
    # seu_i <- SCTransform(seu_i, assay = 'RNA', new.assay.name = 'SCT',
    #                      vars.to.regress = c('percent.mt', 'CC.Difference'),
    #                      conserve.memory = T)
    return(seu_i)
  })
  saveRDS(seu.list, file=preproc_f)
} else {
  seu.list <- readRDS(file=preproc_f)
}


seu.mnn <- RunFastMNN(object.list = seu.list, 
                        features = 2000, assay='SCT')
seu.mnn <- RunUMAP(seu.mnn, reduction = "mnn", dims = 1:30, n.neighbors=30L,
                     min.dist = 0.1, return.model=TRUE)
seu.mnn <- FindNeighbors(seu.mnn, reduction = "mnn", dims = 1:30)
seu.mnn <- FindClusters(seu.mnn, resolution = 0.9, graph.name='SCT_snn')
DefaultAssay(seu.mnn) <- 'RNA'
seu.mnn <- NormalizeData(seu.mnn,
                           normalization.method = "LogNormalize") %>%
  FindVariableFeatures(., selection.method = "vst",
                       nfeatures = 3000, verbose = FALSE) %>% 
  ScaleData(.)
seu.mnn <- RunPCA(seu.mnn, features=VariableFeatures(seu.mnn), ndims.print=1:30,
                    nfeatures.print=5, npcs=30)
saveRDS(seu.mnn, file=file.path(datadir, "seurat_obj", "seu_integ.mnn.rds"))



# Select integration features (genes)
integration_features <- SelectIntegrationFeatures(object.list = seu.list, 
                                                  nfeatures = 3000)
seu.list <- PrepSCTIntegration(object.list = seu.list, 
                               anchor.features = integration_features)

# Find integration anchors across all samples based on variable features
integration_anchors <- FindIntegrationAnchors(object.list = seu.list, 
                                              normalization.method = "SCT",
                                              anchor.features = integration_features)
seu <- IntegrateData(anchorset = integration_anchors, 
                     normalization.method = "SCT")
# VariableFeatures(seu) <- integration_features
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50, n.neighbors=15L,
               min.dist = 0.1)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:50)
seu <- FindClusters(seu, resolution = 1.1)


clusters <- scSHC::scSHC(GetAssayData(seu, slot='data'),
                         parallel=F, cores=1)
clusters_validation <- scSHC::testClusters(GetAssayData(seu, assay='RNA', slot='data'), 
                                           as.character(seu$seurat_clusters))

saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_integ.anchors.rds"))

###########################################
#### 1.b Cell type annotation - ImmGen ####
if(!exists("seu"))seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu_integ.mnn.rds"))

dir.create(file.path(outdir, "annotation"))
getImmgen <- TRUE
overwrite <- TRUE

# RNA preprocess and cell-cycle score
DefaultAssay(seu) <- 'RNA'

# 
# x <- split(as.data.frame(t(assay(bed.se))), colData(bed.se)$label.fine)
# tcells <- c(paste0("T cells (", 
#                    c('T.8MEM', 'T.8Mem', 'T.8MEMKLRG1-CD127+.D8.LISOVA',
#                      paste0("T.CD8.", c('1H', '24H', '48H', '5H', '96H', 'CTR'))), 
#                    ")"))
# split(gsub("_.*", "", rownames(colData(bed.se))), colData(bed.se)$label.fine)[tcells] %>%
#   lapply(., function(i) as.data.frame(t(i))) %>%
#   plyr::rbind.fill(.) %>% 
#   magrittr::set_rownames(tcells) %>%
#   write.table(., sep=",", row.names = T, col.names = T, quote = F)
# xdf <- do.call(cbind, lapply(x[tcells], t)) %>%
#   as.matrix
# storage.mode(xdf) <- 'numeric'
# idx <- list('T8MEM'=c(1:3), 'T8Mem'=c(4:6), 'T8MEMKLRG1'=c(7:9))
# x_pval <- lapply(idx, function(idx_i){
#   sapply(idx, function(idx_j){
#     apply(xdf, 1, function(i){
#       t.test(i[idx_i], i[idx_j])$p.value
#     })
#   })
# })
# x_stat <- lapply(idx, function(idx_i){
#   sapply(idx, function(idx_j){
#     apply(xdf, 1, function(i){
#       t.test(i[idx_i], i[idx_j])$statistic
#     })
#   })
# })
# ggs <- lapply(names(x_pval), function(ids){
#   i <- x_pval[[ids]]
#   i <- i[,-which(colSums(i) == nrow(i)), drop=F]
#   ikp <- i[order(rowSums(i)),]
#   
#   x_stat_sel <- x_stat[[ids]][rownames(ikp[1:100,]), colnames(ikp[1:100,])]
#   b <- c(-40, -20, 0, 20, 40)
#   ggplot(melt(t(x_stat_sel)), aes(x=Var1, y=Var2, fill=value)) + 
#     geom_tile() +
#     theme_cowplot() + 
#     scale_fill_gradientn(limits = c(min(b),max(b)),
#                          colours=c("black", "blue", "white", "yellow", "red"),
#                          breaks=b, labels=format(b))+
#     ggtitle(ids) +
#     theme(axis.text.x=element_text(angle=90))
# })
# pdf("~/xfer/x.pdf", height = 18, width = 16)
# cowplot::plot_grid(plotlist=ggs, nrow=1)
# dev.off()


if(getImmgen){
  #ref <- BlueprintEncodeData()
  #saveRDS(ref, file="~/downloads/BlueprintEncodeData.rds")
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
  
  for(lbl in c('label.fine', 'label.main')){
    for(clus in c(TRUE, FALSE)){
      rds_file <- paste0("celldex_immgen.", gsub("label.", "", lbl),
                         ".", if(clus) 'cluster' else 'cell', "X.rds")
      id <- gsub("^.*immgen(.*).rds", "immgen\\1", rds_file)
      if(!file.exists(file.path(outdir, "annotation", rds_file)) | overwrite){
        print(paste0(" Annotating: ", rds_file))
        singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
                                ref=bed.se, 
                                assay.type.test=1, labels=bed.se[[lbl]],
                                de.method="wilcox", genes='sd',
                                clusters=if(clus) seu$seurat_clusters else NULL)
        df <- lapply(c(1:nrow(singler_anno$scores)), function(i){
          sort(singler_anno$scores[i,], decreasing = T) %>% 
                            head(., 50) %>% 
                            as.data.frame %>%
                            tibble::rownames_to_column("celltype")
        }) %>%
          purrr::reduce(., full_join, by='celltype') %>%
          tibble::column_to_rownames('celltype') %>%
          magrittr::set_colnames(as.character(rownames(singler_anno)))
        write.table(df, "~/xfer/singler.csv", sep=",", quote = F, col.names = T, row.names = )
        
        saveRDS(singler_anno, file=file.path(outdir, "annotation", rds_file))
      } else {
        print(paste0(" Reading in annotation: ", rds_file))
        singler_anno <- readRDS(file.path(outdir, "annotation", rds_file))
        # pdf("~/xfer/test3.pdf", width=15)
        # plotScoreHeatmap(singler_anno, show_colnames=TRUE)
        # dev.off()
      }
      
      if(clus){
        cluster_ids <- setNames(singler_anno$pruned.labels, 
                                as.character(rownames(singler_anno)))
        seu@meta.data[,id] <- cluster_ids[as.character(seu$seurat_clusters)]
      } else {
        seu@meta.data[,id] <- singler_anno$pruned.labels
      }
    }
  }
} 
DefaultAssay(seu) <- 'mnn.reconstructed'
colnames(seu@meta.data) <- gsub("X$", "", colnames(seu@meta.data))

# Makes annotations simpler (e.g. "B Cells (B.Fo) -> B.Fo)
id_map <- unique(seu$immgen.fine.cluster) %>%
  gsub("^.*\\(", "", .) %>%
  gsub("\\)", "", .) %>%
  gsub("(^.*?\\..*)\\..*", "\\1", .) %>%
  setNames(., unique(seu$immgen.fine.cluster))
seu$immgen.fine.cluster.simple <- id_map[seu$immgen.fine.cluster]
seu$immgen.fine.cell.simple <- id_map[seu$immgen.fine.cell]

# # Relabel clusters based on manually curated markers
# relabel_map <- unique(seu@meta.data[,c('immgen.fine.cluster.simple', 'seurat_clusters')]) %>%
#   as.data.frame %>% 
#   tibble::remove_rownames(.) %>% 
#   tibble::column_to_rownames(., 'seurat_clusters') %>%
#   rename_with(., ~"celltype")
# rename_key <- c('8'='T.Tregs.Klf2+', '10'='T.Tregs.Klf2+', '19'='T.Tregs', 
#                 '20'='Unknown', '7'='T.4MEM', '17'='T.8MEM')
# relabel_map[names(rename_key), 'celltype'] <- rename_key
# seu$`immgen.fine.cluster.simple` <- relabel_map[as.character(seu$seurat_clusters),'celltype']
if(visualize){
  dpcl <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="seurat_clusters", pt.size=0.5, shuffle=TRUE)
  dp1 <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="immgen.main.cell", pt.size=0.5, shuffle=TRUE)
  dp1b <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="immgen.main.cluster", pt.size=0.5, shuffle=TRUE)
  dp1c <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="immgen.fine.cluster.simple", pt.size=0.5, shuffle=TRUE)
  dp1d <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="immgen.fine.cell.simple", pt.size=0.5, shuffle=TRUE)
  dp2 <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="orig.ident", pt.size=0.5, shuffle=TRUE)
  fp <- FeaturePlot(seu, features=c('nCount_RNA', 'G2M.Score', 'percent.mt'), 
                    pt.size=0.5, order=T, ncol=3)
  pdf(file.path("~/xfer", "dimplot_immgen_anno.pdf"), width=15)
  plot_grid(dpcl, dp1c, ncol=2)
  plot_grid(dp2, dp1d, ncol=2)
  fp
  dev.off()
  pdf(file.path(outdir, "annotation", "dimplot_immgen_anno.pdf"), width=15)
  pdf(file.path("~/xfer", "dimplot_immgen_anno.pdf"), width=15)
  plot_grid(dpcl, dp1, ncol=2)
  plot_grid(dpcl, dp1b, ncol=2)
  plot_grid(dpcl, dp1c, ncol=2)
  plot_grid(dpcl, dp1d, ncol=2)
  plot_grid(dp1c, dp2, ncol=2)
  dev.off()
}


saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_integ.anno.rds"))
##########################################
#### 1.c Cell type annotation - Other ####
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.rds"))
# if(!exists("seul"))  seul <- readRDS(file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.rds"))
# if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))
sctype_dir <- '/cluster/projects/mcgahalab/bin/sc-type/R'
sctype_db <- 'ScTypeDB_full.xlsx'
sctype_tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 


#--- i) outdated - ScType -------------------------------------------------------
for(f in list.files(sctype_dir)) source(file.path(sctype_dir, f))

# prepare gene sets
gs_list = gene_sets_prepare(sctype_db, sctype_tissue)

# get cell-type by cell matrix
es_max = sctype_score(scRNAseqData = seu[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seu@meta.data$seurat_clusters), function(cl){
  es_max_cl = sort(rowSums(es_max[ ,rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl, ])]), 
                   decreasing = !0)
  head(data.frame(cluster = cl, 
                  type = names(es_max_cl), 
                  scores = es_max_cl, 
                  ncells = sum(seu@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
seu$sctype_anno <- setNames(sctype_scores$type, sctype_scores$cluster)[seu$seurat_clusters]


pdf("~/xfer/sara_sctype.pdf")
DimPlot(seu, group.by="immgen.fine.cluster.simple")
DimPlot(seu, group.by="sctype_anno", label = TRUE)
dev.off

#--- ii) outdated - scSorter ----------------------------------------------------
library(scSorter)

# prepare gene sets
source(file.path(sctype_dir, 'gene_sets_prepare.R'))
gs_list = gene_sets_prepare(sctype_db, sctype_tissue)
anno <- lapply(names(gs_list), function(gs_id){
  lapply(names(gs_list[[gs_id]]), function(ct_id){
    if(length(gs_list[[gs_id]][[ct_id]]) == 0) return(NULL)
    data.frame("Celltype"=ct_id,
               "gs"=gs_id,
               "Type"=paste0(ct_id, ".", gs_id),
               "Marker"=as.character(gs_list[[gs_id]][[ct_id]]),
               "Weight"=2)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
# anno$Marker <- stringr::str_to_title(anno$Marker)

DefaultAssay(seu) <- 'RNA'
topgenes <- head(VariableFeatures(seu), 2000)
topgenes <- toupper(topgenes)

expr = GetAssayData(seu)
rownames(expr) <- toupper(rownames(expr))
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

annol <- split(anno, f=anno$Type)
idx <- sapply(annol, function(a){sum(a$Marker %in% rownames(expr))>5})
anno <- do.call(rbind, annol[idx])

picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]
anno = anno[anno$Marker %in% picked_genes, ]
rts <- scSorter(expr, anno)

seu$scsort_anno <- gsub(".gs_positive", "", rts$Pred_Type)


pdf("~/xfer/sara_sctype.pdf")
DimPlot(seu, group.by="immgen.fine.cluster.simple")
DimPlot(seu, group.by="sctype_anno", label = TRUE)
dev.off()
pdf("~/xfer/sara_sctype.pdf", width = 14)
DimPlot(seu, group.by="scsort_anno", label = TRUE)
dev.off()

saveRDS(seu, file = file.path(datadir, "seurat_obj", "4_seu_degPrepX.rds"))
#--- iii) ProjecTILs ----
dir.create(file.path(outdir, "annotation", "projectils"), showWarnings = F)
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/spica'

# refobj_dc <- readRDS(file.path(projectils_dir, 'DC_human_ref_v1.rds'))
# refobj_dc <- orthogene.scgate(refobj_dc, from='human', 'mouse')
# refobj_cd4 <- readRDS(file.path(projectils_dir, 'CD4T_human_ref_v1.rds'))
# refobj_cd8 <- readRDS(file.path(projectils_dir, 'CD8T_human_ref_v1.rds'))
# refobj_momac <- readRDS(file.path(projectils_dir, "APC_atlas_v1_SPICA.rds"))
# refobj_tit <- readRDS(file.path(projectils_dir, "ref_TILAtlas_mouse_v1.rds"))
# refobj_custom2 <- readRDS(file.path(projectils_dir, 'custom', 'dc.cd4.cd8.rds'))
# refobj_custom <- readRDS(file.path(projectils_dir, 'custom', 'custom_ilc2.TILatlas.rds'))
# refobj_custom <- readRDS(file.path(projectils_dir, 'custom', 'custom_ilc2.mouse_atlas.rds'))
refobj_custom <- readRDS(file.path(projectils_dir, 'custom', 'custom_ilc2.mouse_atlasB.rds'))
refobj_custom <- readRDS(file.path(projectils_dir, 'custom', 'custom_ilc2.treg.rds'))
# refobj_custom <- readRDS(file.path(projectils_dir,'custom',  "seu.mnn.small.sara.postSCT.rds")) #seu.mnn.small.sara.rds"))
DefaultAssay(refobj_custom) <- 'RNA'

DefaultAssay(seu) <- 'RNA'
subset_size <- 3
set.seed(1234)
cell_idx <- split(Cells(seu), f = sample(c(1:subset_size), size=ncol(seu), replace = T))
batch <- 1
ptils_map <- lapply(cell_idx, function(cells_i){
  batchf <- file.path(outdir, "annotation", "projectils", 
                      paste0("batch", batch, ".rds"))
  
  if(!file.exists(batchf)){
    print(paste0("batch: ", batch, "..."))
    seu_i <- subset(seu, cells=cells_i)
    seu_i <- ProjecTILs.classifier(seu_i, refobj_custom, ncores = 1, 
                                   split.by = "orig.ident", filter.cells=F)
    seu_i$functional.cluster.custom <- seu_i$functional.cluster
    seu_lbl <- seu_i$functional.cluster.custom
    saveRDS(seu_lbl, file=batchf)
    
    print("plotting...")
    pdf("~/xfer/test.pdf", width = 16, height = 12)
    dp1 <- DimPlot(seu_i, group.by='seurat_clusters', reduction='umap', raster=T, label = T)
    dp2 <- DimPlot(seu_i, group.by='functional.cluster.custom', reduction='umap', raster=T, label = T)
    plot(dp1 + dp2 )
    
    Idents(seu_i) <- 'functional.cluster.custom'
    dp_clusters <- lapply(unique(seu_i$functional.cluster.custom), function(clid){
      scCustomize::Cluster_Highlight_Plot(seu_i, cluster_name = clid, highlight_color = "red",
                                          reduction='umap', background_color = "lightgray", raster=T) + 
        NoLegend() +
        ggtitle(clid)
    })
    plot(cowplot::plot_grid(plotlist = dp_clusters, ncol=7))
    dev.off()
  } else {
    seu_lbl <- readRDS(file=batchf)
  }
  batch <<- batch+1
  return(seu_lbl)
})
ptils_map <- do.call(c, ptils_map)
ptils_map <- ptils_map[which(ptils_map >= (length(ptils_map) * 0.0005))]
names(ptils_map) <- gsub("^.*?\\.", "", names(ptils_map))

seu$functional.cluster.custom <- NA
seu@meta.data[names(ptils_map),'functional.cluster.custom'] <- ptils_map

## Visualize
dp1 <- DimPlot(seu, group.by='seurat_clusters', reduction='umap', raster=T, label = T) + NoLegend()
dp2 <- DimPlot(seu, group.by='functional.cluster.custom', reduction='umap', raster=T, label = T) + NoLegend()
dp3 <- DimPlot(seu, group.by='immgen.fine.cluster.simple', reduction='umap', raster=T, label = T) + NoLegend()

Idents(seu) <- 'functional.cluster.custom'
dp_clusters_immgen <- lapply(na.omit(unique(seu$functional.cluster.custom)), function(clid){
  scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = adjustcolor('red', alpha.f=0.4),
                                      reduction='umap', background_color = "lightgray", raster=T) + 
    NoLegend() +
    ggtitle(clid)
})

Idents(seu) <- 'immgen.fine.cluster.simple'
dp_clusters_ptils <- lapply(na.omit(unique(seu$immgen.fine.cluster.simple)), function(clid){
  scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = adjustcolor(c('red', 'red'), alpha.f=0.4),
                                      reduction='umap', background_color = "lightgray", raster=T) + 
    NoLegend() +
    ggtitle(clid)
})

Idents(seu) <- 'orig.ident'
dp_clusters_id <- lapply(na.omit(unique(seu$orig.ident)), function(clid){
  scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = adjustcolor('red', alpha.f=0.4),
                                      reduction='umap', background_color = "lightgray", raster=T) + 
    NoLegend() +
    ggtitle(clid)
})

.barplotfun <- function(x, colid='functional.cluster.custom', 
                        percentage=F, ncol=8){
  meltdat <- table(seu@meta.data[,colid], seu$orig.ident) 
  if(percentage) meltdat <- apply(meltdat, 2, function(i) i/sum(i))
  topids <- sort(rowSums(meltdat), decreasing = T)
  meltdat <- melt(meltdat) %>%
    mutate(Var1 = factor(Var1, levels=names(topids)))
  
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
  cols=setNames(c(c25[1:ncol], rep("grey", (length(topids)-(ncol)))),
                names(topids))
  
  ggplot(meltdat, 
         aes(x=Var2, fill=Var1, y=value, color=Var1)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values=cols) +
    theme_bw() +
    xlab("Samples") + ylab("Count") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}
bp1a <- .barplotfun(seu, 'functional.cluster.custom', TRUE)
bp2a <- .barplotfun(seu, 'functional.cluster.custom', FALSE)
bp1b <- .barplotfun(seu, 'immgen.fine.cluster.simple', TRUE)
bp2b <- .barplotfun(seu, 'immgen.fine.cluster.simple', FALSE)


pdf(file.path(outdir, "annotation", "dimplot_annotation.pdf"), width = 18, height = 12)
# pdf(file.path("~/xfer", "dimplot_annotation.pdf"), width = 18, height = 12)
dp1 + dp2  + dp3
cowplot::plot_grid(plotlist = dp_clusters_immgen, ncol=6)
cowplot::plot_grid(plotlist = dp_clusters_ptils, ncol=6)
cowplot::plot_grid(plotlist = dp_clusters_id, ncol=7)
cowplot::plot_grid(plotlist=list(bp1a, bp2a), ncol=1)
cowplot::plot_grid(plotlist=list(bp1b, bp2b), ncol=1)
dev.off()

saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_integ.anno.rds"))

#--- iv.a) TReg - ST2 condition ----
dir.create(file.path(outdir, "annotation", "projectils"), showWarnings = F)
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/ProjecTILs'
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)
do_visualize <- FALSE
do_deg <- FALSE

if(FALSE){
  seu <- seul[[1]]
  x <- seu@reductions$pacmap@cell.embeddings
  xidx <- which(x[,1] > quantile(x[,1], 0.99) | x[,1] < quantile(x[,1], 0.01))
  Cells(seu)[xidx]
  seurm <- subset(seu, cells=Cells(seu)[xidx])
  seu <- subset(seu, cells=setdiff(Cells(seu), Cells(seu)[xidx]))
  
  if(TRUE){
    library(reticulate)
    use_python("/cluster/home/quever/miniconda3/bin/python3.9")
    python_pacmap <- import("pacmap")
    python_pandas <- import("pandas")
    python_numpy <- import("numpy")
    
    dat_pd <- reticulate::r_to_py(as.data.frame(t(GetAssayData(seu, slot='data'))))
    # dat_pd <- reticulate::r_to_py(as.data.frame(seu@reductions$mnn@cell.embeddings))
    nparray <- dat_pd$values
    nparray <- nparray$astype(python_numpy$float)
    tfmr <- python_pacmap$PaCMAP(n_neighbors=10L, MN_ratio=0.5, FP_ratio=2.0) 
    X_transformed <- tfmr$fit_transform(nparray)#, init="pca")
    seu[["pacmap"]] <- CreateDimReducObject(embeddings = data.frame(X_transformed) %>%
                                              magrittr::set_rownames(., Cells(seu)) %>%
                                              as.matrix, 
                                            key = "PacMAP_", 
                                            assay = DefaultAssay(seu))
  }
  
  pdf("~/xfer/y2.pdf", height = 5, width = 20)
  # TMP
  dp1 <- DimPlot(seu, group.by="manual_anno", raster=T, reduction='umap', label=T)
  dp2 <- DimPlot(seu, group.by="manual_anno", raster=T, reduction='mnn', label=T)
  dp3 <- DimPlot(seu, group.by="manual_anno", raster=T, reduction='pacmap', label=T)
  cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
  
  dp1 <- DimPlot(seu, group.by="seurat_clusters", raster=T, reduction='umap', label=T)
  dp2 <- DimPlot(seu, group.by="seurat_clusters", raster=T, reduction='mnn', label=T)
  dp3 <- DimPlot(seu, group.by="seurat_clusters", raster=T, reduction='pacmap', label=T)
  cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
  
  dp1 <- DimPlot(seurm, group.by="seurat_clusters", raster=T, reduction='umap', label=T)
  dp2 <- DimPlot(seurm, group.by="manual_anno", raster=T, reduction='umap', label=T)
  dp3 <- DimPlot(seurm, group.by="treg_anno", raster=T, reduction='umap', label=T)
  cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
  dev.off()
}

### Set up the Dataset data
# Generated from Dataset2 in stacas_treg_ref.R
refobj_ds2 <- readRDS(file.path(projectils_dir, 'seu.treg.rds'))
DefaultAssay(refobj_ds2) <- 'RNA'
# Generated from Dataset2 in stacas_treg_ref.R
refobj_ds5 <- readRDS(file.path(projectils_dir, 'seu.treg_NLT_LT.rds'))
DefaultAssay(refobj_ds5) <- 'RNA'
# Dataset3 is an expression of all the genes for the different cluster
dataset3_expr <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/treg_ref/genesets/dataset3_expr.csv'
d3_expr <- read.csv(dataset3_expr) %>%
  as.data.frame %>%
  filter(!duplicated(SYMBOL)) %>%
  tibble::column_to_rownames(., "SYMBOL")

seul_treg <- lapply(seul, function(seu){
  Idents(seu) <- 'manual_anno'
  DefaultAssay(seu) <- 'RNA'
  subset(seu, ident=intersect(c('Treg', 'CD4_CCR7hi'),
                              unique(Idents(seu))))
})
rm(seul); gc()

seul_tregs <- lapply(seul_treg, function(treg_i){
  # ----- Dataset 2: ProjecTILs projection
  print("Dataset 2...")
  treg_i <- ProjecTILs.classifier(treg_i, refobj_ds2, ncores = 1, 
                                  split.by = "orig.ident", filter.cells=F)
  treg_i$treg_dataset2 <- treg_i$functional.cluster
  
  # ----- Dataset 5: ProjecTILs projection
  print("Dataset 5...")
  treg_i <- ProjecTILs.classifier(treg_i, refobj_ds5, ncores = 1, 
                                  split.by = "orig.ident", filter.cells=F)
  treg_i$treg_dataset5 <- treg_i$functional.cluster
  
  # ----- Dataset 3: Spearman Correlation with clusters
  print("Dataset 3...")
  expr_i <- GetAssayData(treg_i, slot='data')
  genes_i <- rownames(expr_i)
  refgenes <- rownames(d3_expr)
  genes_idx <- which(genes_i %in% refgenes)
  d3_expr_ord <- d3_expr[na.omit(match(genes_i, refgenes)),]
  expr_i_ord <- expr_i[genes_idx,]
  
  cormat <- apply(expr_i_ord, 2, function(i){
    apply(d3_expr_ord, 2, function(ref){
      cor(i, ref, method='spearman', use='complete.obs')
    })
  })
  
  max_clusts <- rownames(cormat)[apply(cormat, 2, which.max)]
  treg_i$treg_dataset3 <- max_clusts
  
  return(treg_i)
})
saveRDS(seul_tregs, file=file.path(datadir, "seurat_obj", "tregs.rds"))

seul_tregs <- readRDS(file=file.path(datadir, "seurat_obj", "tregs.rds"))

anno_list <- list()
treg_ds <- c(paste0('treg_dataset', c('2', '3', '3_agg', '5')), "treg_manual_anno")
process_combined <- TRUE
remove_cells <- FALSE
seul_tregs_final <- lapply(setNames(names(seul_tregs),names(seul_tregs)), function(seu_id){
  seu_i <- seul_tregs[[seu_id]]
  seul_days <- lapply(c("Day3", "Day7"), function(batchid){
    seu_i$old_clusters <- seu_i$seurat_clusters
    # VariableFeatures(seu_i) <- vargenes
    batch <- c('Day3'='^B1', 'Day7'='^B2')[batchid]
    
    dims_use <- 1:20
    print("Treating samples Day 3 and 7 separate")
    # seu_j <- SCTransform(seu_j, vst.flavor = "v2", 
    #                      verbose = FALSE, vars.to.regress='CC.Difference')
    seu_j <- subset(seu_i, cells=Cells(seu_i)[grep(batch, seu_i$orig.ident)])
    seu_j <- seu_j %>%
      FindVariableFeatures %>%
      ScaleData %>%
      RunPCA %>%
      FindNeighbors(.,dims = dims_use) %>%
      FindClusters(., resolution = 1.4) %>%
      RunUMAP(., dims=dims_use) 
    seu_j <- CellCycleScoring(seu_j, 
                              s.features = stringr::str_to_title(cc.genes$s.genes), 
                              g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                              set.ident = TRUE)
    seu_j$CC.Difference <- seu_j$S.Score - seu_j$G2M.Score
    d3_map <- setNames(paste0('cl', c('123', '123', '123', '4', '56', '56')), 
                       paste0("spleen_tr_cluster", c(1:6)))
    seu_j$treg_dataset3_agg <- d3_map[seu_j$treg_dataset3]
    
    # Manually annotate the clusters (res 1.4) and remove the problematic clusters
    seu_j$treg_manual_anno <- as.character(seu_j$seurat_clusters) %>% 
      treg_recode_map(., grp=seu_id, day=batchid, anno=TRUE)
    
    Idents(seu_j) <- 'treg_manual_anno'
    if(remove_cells){
      seu_j <- subset(seu_j, ident=setdiff(unique(Idents(seu_j)), 'Remove'))
      seu_j <- RunUMAP(seu_j, dims=dims_use, reduction='mnn')
      seu_j <- runPacmap(seu_j, reduction='mnn')
    }
    
    return(seu_j)
  })
  names(seul_days) <- c('Day3', 'Day7')
  
  # Get the day-specific treg annotations
  anno_list <- lapply(seul_days, function(seu_j){
    lapply(setNames(treg_ds, treg_ds), function(i) setNames(seu_j@meta.data[,i], Cells(seu_j)))
  })
  
  if(process_combined){
    print("Treating samples as Day 3 and 7 combined")
    seu_i$old_clusters <- seu_i$seurat_clusters
    dims_use <- 1:20
    
    seu_j <- seu_i  %>% 
      FindVariableFeatures %>%
      ScaleData %>%
      RunPCA %>%
      RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
      RunTSNE(., dims=dims_use, reduction='mnn') %>%
      FindNeighbors(.,dims = dims_use, reduction='mnn') %>%
      FindClusters(., resolution = 1.4, graph.name='SCT_snn')
    
    seu_j <- CellCycleScoring(seu_j, 
                              s.features = stringr::str_to_title(cc.genes$s.genes), 
                              g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                              set.ident = TRUE)
    seu_j$CC.Difference <- seu_j$S.Score - seu_j$G2M.Score
    d3_map <- setNames(paste0('cl', c('123', '123', '123', '4', '56', '56')), 
                       paste0("spleen_tr_cluster", c(1:6)))
    seu_j$treg_dataset3_agg <- d3_map[seu_j$treg_dataset3]
    
    # Carry over the manual annotations from Day3 and 7
    for(ds_i in treg_ds){
      for(day_j in names(anno_list)){
        colid <- paste0(day_j, "_", ds_i)
        seu_j@meta.data[,colid] <- anno_list[[day_j]][[ds_i]][as.character(Cells(seu_j))]
      }
    }
    seu_j$treg_manual_anno <- c(anno_list$Day3$treg_manual_anno,
                                anno_list$Day3$treg_manual_anno)[as.character(Cells(seu_j))]
    
    cells <- sapply(anno_list, function(i) names(i$treg_manual_anno)) %>% unlist
    
    if(remove_cells){
      seu_j <- subset(seu_j, cells=cells)
      seu_j <- RunUMAP(seu_j, dims=dims_use, reduction='mnn')
      seu_j <- runPacmap(seu_j, reduction='mnn')
    }
    
    seul_days[['combined']] <- seu_j
  }
  return(seul_days)
})
saveRDS(seul_tregs_final, file=file.path(datadir, "seurat_obj", "tregs_final.rds"))
seul_tregs_final <- readRDS(file=file.path(datadir, "seurat_obj", "tregs_final.rds"))

seul_tregs_anno <- lapply(names(seul_tregs_final), function(seu_id){
  seu_i <- seul_tregs_final[[seu_id]]
  dayids <- c('Day3', 'Day7')
  seul_tregs_days <- lapply(setNames(dayids, dayids), function(day_id){
    seu_j <- seu_i[[day_id]]
    # Manually annotate the clusters (res 1.4) and remove the problematic clusters
    seu_j$treg_manual_anno <- as.character(seu_j$seurat_clusters) %>% 
      treg_recode_map(., grp=seu_id, day=day_id, anno=TRUE)
    return(seu_j)
  })
  return(seul_tregs_days)
}) %>%
  setNames(., names(seul_tregs_final))


seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)
st2_treg_map <- lapply(seul_tregs_anno, function(i) {
  lapply(i, function(j) j$treg_manual_anno) %>%
    do.call(c, .)
  }) %>% 
  do.call(c, .)
day_st2_treg_map <- lapply(seul_tregs_anno, function(i) {
  lapply(names(i[1:2]), function(dayid) {
    setNames(paste0(dayid, "_", i[[dayid]]$treg_manual_anno),
             Cells(i[[dayid]]))
  }) %>%
    do.call(c, .)
}) %>% 
  do.call(c, .)
names(st2_treg_map) <- gsub("^.*\\.", "", names(st2_treg_map))
saveRDS(st2_treg_map, file=file.path(datadir, "annotation_map", "st2_tregs.rds"))
st2_tcell_map <- readRDS(file=file.path(datadir, "annotation_map", "st2_tcells.rds"))
st2_treg_map <- readRDS(file=file.path(datadir, "annotation_map", "st2_tregs.rds"))

names(day_st2_treg_map) <- gsub("^.*\\.", "", names(day_st2_treg_map))
for(lt in c('LN', 'Tumor')){
  seul_tregs_final[[lt]]$combined$treg_manual_anno <- st2_treg_map[Cells(seul_tregs_final[[lt]]$combined)]
  seul[[lt]]$treg_manual_anno <- st2_treg_map[Cells(seul[[lt]])]
  seul[[lt]]$tcell_manual_anno <- st2_tcell_map[Cells(seul[[lt]])]
  cell_idx <- which(Cells(seul[[lt]]) %in% names(st2_treg_map))
  seul[[lt]]$manual_anno <- as.character(seul[[lt]]$manual_anno)
  seul[[lt]]$manual_anno[cell_idx] <- st2_treg_map[Cells(seul[[lt]])[cell_idx]]
  cell_idx <- which(Cells(seul[[lt]]) %in% names(st2_tcell_map))
  seul[[lt]]$manual_anno <- as.character(seul[[lt]]$manual_anno)
  seul[[lt]]$manual_anno[cell_idx] <- st2_tcell_map[Cells(seul[[lt]])[cell_idx]]
  
  day3map <- grep("Day3_", day_st2_treg_map, value=T)
  seul_tregs_final[[lt]]$combined$Day3_treg_manual_anno <- day3map[Cells(seul_tregs_final[[lt]]$combined)]
  day7map <- grep("Day7_", day_st2_treg_map, value=T)
  seul_tregs_final[[lt]]$combined$Day7_treg_manual_anno <- day7map[Cells(seul_tregs_final[[lt]]$combined)]
}

if(FALSE){
  if(FALSE){
    if(do_visualize){
      pdf("~/xfer/Treg_LN_Day7.res09.pdf", width = 10)
      dp1 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='umap', label=lbl)
      dp2 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='pacmap', label=lbl)
      cowplot::plot_grid(dp1,dp2)
      dev.off()
      
      if(batchid != 'combined'){
        for(treg_ds  in c(paste0('treg_dataset', c('3_agg', '5', '2')), 'seurat_clusters')){
          anno_list[[treg_ds]][[batchid]] <- setNames(as.character(seu_j@meta.data[,treg_ds]),
                                                      Cells(seu_j))
        }
      } else {
        seu_j$Day3_seurat_clusters <- anno_list[['seurat_clusters']][['Day3']][as.character(Cells(seu_j))]
        seu_j$Day7_seurat_clusters <- anno_list[['seurat_clusters']][['Day7']][as.character(Cells(seu_j))]
        
        pdf("~/xfer/day37_combined.tumor.pdf", height = 5, width = 15)
        lbl=T
        dp1 <- DimPlot(seu_j, group.by="Day3_seurat_clusters", raster=T, reduction='umap', label=lbl)
        dp2 <- DimPlot(seu_j, group.by="Day7_seurat_clusters", raster=T, reduction='umap', label=lbl)
        dp3 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='umap', label=lbl)
        cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
        
        dp1 <- DimPlot(seu_j, group.by="Day3_seurat_clusters", raster=T, reduction='pacmap', label=lbl)
        dp2 <- DimPlot(seu_j, group.by="Day7_seurat_clusters", raster=T, reduction='pacmap', label=lbl)
        dp3 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='pacmap', label=lbl)
        cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
        dev.off()
      }
      
      pdf("~/xfer/y12.pdf", height = 5, width = 20)
      lbl=FALSE
      # TMP
      dp1 <- DimPlot(seu_j, group.by="manual_anno", raster=T, reduction='umap', label=lbl)
      dp2 <- DimPlot(seu_j, group.by="manual_anno", raster=T, reduction='mnn', label=lbl)
      dp3 <- DimPlot(seu_j, group.by="manual_anno", raster=T, reduction='pacmap', label=lbl)
      # dp3 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='tsne', label=lbl)
      fps <- FeaturePlot(seu_j, features ="S.Score", raster=T, reduction='umap')
      fpg2m <- FeaturePlot(seu_j, features ="G2M.Score", raster=T, reduction='umap')
      fp <- FeaturePlot(seu_j, features ="CC.Difference", raster=T, reduction='umap')
      cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
      cowplot::plot_grid(plotlist=list(dp1, fps, fpg2m, fp), nrow=1)
      
      dp1 <- DimPlot(seu_j, group.by="manual_clusters", raster=T, reduction='umap', label=lbl)
      dp2 <- DimPlot(seu_j, group.by="manual_clusters", raster=T, reduction='mnn', label=lbl)
      dp3 <- DimPlot(seu_j, group.by="manual_clusters", raster=T, reduction='pacmap', label=lbl)
      cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
      
      dp2 <- DimPlot(seu_j, group.by="treg_dataset2", raster=T, reduction='umap', label=lbl)
      dp5 <- DimPlot(seu_j, group.by="treg_dataset5", raster=T, reduction='umap', label=lbl)
      dp3 <- DimPlot(seu_j, group.by="treg_dataset3_agg", raster=T, reduction='umap', label=lbl)
      dpsc <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='umap', label=T)
      cowplot::plot_grid(plotlist=list(dp2, dp5, dp3, dpsc), nrow=1)
      
      dp2 <- DimPlot(seu_j, group.by="treg_dataset2", raster=T, reduction='pacmap', label=lbl)
      dp5 <- DimPlot(seu_j, group.by="treg_dataset5", raster=T, reduction='pacmap', label=lbl)
      dp3 <- DimPlot(seu_j, group.by="treg_dataset3_agg", raster=T, reduction='pacmap', label=lbl)
      dpsc <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='pacmap', label=T)
      cowplot::plot_grid(plotlist=list(dp2, dp5, dp3, dpsc), nrow=1)
      dev.off()
      table(seu_j$treg_dataset5, seu_j$seurat_clusters)
      
      pdf("~/xfer/y10a.pdf")
      DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='umap', label=T)
      dev.off()
    }
    
    if(do_deg){
      Idents(seu_j) <- 'seurat_clusters'
      seu_cl <- subset(seu_j, ident='3')
      x2 <- FindMarkers(seu_j, ident.1='B1_LN_KO_3d', ident.2='B1_LN_WT_3d', 
                        group.by='orig.ident') %>% 
        # latent.vars='CC.Difference', test.use='LR') %>%
        tibble::rownames_to_column(., "gene") %>%
        mutate(ens=gm$SYMBOL$ENSEMBL[gene],
               biotype=ens2biotype_ids[ens]) %>%
        relocate(., c(ens, biotype))
      write.table(x, file="~/xfer/deg_cluster11_vs_cluster3.csv",
                  sep=",", col.names = T, row.names = F, quote = F)
    }
  }
}

if(do_visualize){
  total_ggs <- lapply(names(seul_tregs_final), function(seu_id){
    seul_i <- seul_tregs_final[[seu_id]]
    day_ggs <- lapply(setNames(names(seul_i),names(seul_i)), function(day_id){
      seu_i <- seul_i[[day_id]]
      
      if(day_id != 'combined'){
        #Day3 and Day7 plots
        gg1 <- lapply(treg_ds, function(feature){
          DimPlot(seu_i, group.by=feature, raster=T, label=F, reduction='umap')
        }) %>%
          cowplot::plot_grid(plotlist=., nrow=1)
        gg2 <- lapply(treg_ds, function(feature){
          DimPlot(seu_i, group.by=feature, raster=T, label=F, reduction='pacmap')
        }) %>%
          cowplot::plot_grid(plotlist=., nrow=1)
        gg <- cowplot::plot_grid(plotlist=list(gg1, gg2), nrow=2)
      } else {
        #Combined plots compaed to Day3 and Day7 labels
        gg1 <- lapply(treg_ds, function(feature){
          lapply(paste0(c("Day3_", "Day7_", ""), feature), function(feat_i){
            print(feat_i)
            tryCatch({
              DimPlot(seu_i, group.by=feat_i, raster=T, label=F, reduction='umap')
            }, error=function(e){NULL})
          }) %>%
            cowplot::plot_grid(plotlist=., nrow=1)
        }) %>%
          cowplot::plot_grid(plotlist=., ncol=1)
        
        gg2 <- lapply(treg_ds, function(feature){
          lapply(paste0(c("Day3_", "Day7_", ""), feature), function(feat_i){
            tryCatch({
              DimPlot(seu_i, group.by=feat_i, raster=T, label=F, reduction='pacmap')
            }, error=function(e){NULL})
          }) %>%
            cowplot::plot_grid(plotlist=., nrow=1)
        }) %>%
          cowplot::plot_grid(plotlist=., ncol=1)
        gg <- list(gg1, gg2)
      }
      return(gg)
    })
    return(day_ggs)
  })
  names(total_ggs) <- names(seul_tregs_final)
  pdf("~/xfer/treg_annotations.days.pdf", width = 20)
  total_ggs$LN[c('Day3', 'Day7')]
  dev.off()
  pdf("~/xfer/treg_annotations.combined.pdf", height = 20, width = 13)
  total_ggs$LN['combined']
  dev.off()
}

if(do_deg){
  seu_i <- seul_tregs_final$LN$Day3
  Idents(seu_i) <- 'treg_manual_anno'
  celltypes <- c('NLTlike.Effector', 'NLTlike.Central')
  degs <- lapply(setNames(celltypes,celltypes), function(id){
    seu_ij <- subset(seu_i, ident=id)
    Idents(seu_ij) <- 'orig.ident'
    table(seu_ij$orig.ident)
    deg <- FindMarkers(seu_ij, group.by='orig.ident',
                       ident.1='B1_LN_KO_3d', ident.2='B1_LN_WT_3d') %>%
      tibble::rownames_to_column(., "gene") %>%
      mutate(celltype=id, 
             comparison='LNday3_KOtrx_vs_WTtrx',
             ens=gm$SYMBOL$ENSEMBL[gene],
             biotype=ens2biotype_ids[ens]) %>%
      relocate(., c(ens, biotype))
    return(deg)
  })
  for(deg_ctid in names(degs)){
    write.table(degs[[deg_ctid]], file=file.path("~/xfer", paste0("deg.LNday3.", deg_ctid, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  }
}










# Use ssGSEA to validate the same gene expression as in the treg paper
treg_markers <- list("stat1"=c('Ifit1', 'Ifit3', 'Rsad2', 'Usp18', 'Isg15',
                               'Rtp4', 'Irf7', 'Bst2', 'Stat1', 'Stat2',
                               'Ms4a4b', 'Ms4a6b', 'Samhd1'),
                     "Central"=c('Sell', 'Ccr7', 'Bcl2', 'Ly6c1', 'S1pr1'),
                     "Effector"=c('Izumo1r', 'Tnfrsf9', 'Tnfrsf18', 'Tnfrsf4',
                                  'Pdcd1', 'Nrp1', 'Ctla4', 'Icos', 'Cd83'),
                     'NLT_like'=c('S100a4', 'S100a6', 'Ccr2', 'Itgb1', 'Itgb7',
                                  'Cxcr3', 'Capg', 'Ifngr1', 'Ly6a', 'Il7r'))
seul_tregs <- lapply(seul_tregs, function(treg_i){
  DefaultAssay(treg_i) <- 'RNA'
  expr <- GetAssayData(treg_i, slot='counts')
  idx <- which(rowSums(expr==0) < (ncol(treg_i)*0.95))
  #--- ssGSEA from GSVA package
  ssgsea_score = GSVA::gsva(expr = expr[idx,], treg_markers, 
                            method = "ssgsea", 
                            parallel.sz = 1, verbose = T)
  treg_i$ssgsea_max <- NA
  treg_i$ssgsea_max <-  rownames(ssgsea_score)[apply(ssgsea_score, 2, which.max)]
  return(treg_i)
})
saveRDS(seul_tregs, file=file.path(datadir, "seurat_obj", "tregs.rds"))


#--- iv.b) TReg - CD45 cells ----
dir.create(file.path(outdir, "annotation", "projectils"), showWarnings = F)
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/ProjecTILs'
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)
do_visualize <- FALSE
do_deg <- FALSE

### Set up the Dataset data
# Generated from Dataset2 in stacas_treg_ref.R
refobj_ds2 <- readRDS(file.path(projectils_dir, 'seu.treg.rds'))
DefaultAssay(refobj_ds2) <- 'RNA'
# Generated from Dataset2 in stacas_treg_ref.R
refobj_ds5 <- readRDS(file.path(projectils_dir, 'seu.treg_NLT_LT.rds'))
DefaultAssay(refobj_ds5) <- 'RNA'
# Dataset3 is an expression of all the genes for the different cluster
dataset3_expr <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/treg_ref/genesets/dataset3_expr.csv'
d3_expr <- read.csv(dataset3_expr) %>%
  as.data.frame %>%
  filter(!duplicated(SYMBOL)) %>%
  tibble::column_to_rownames(., "SYMBOL")

seul_tregs <- lapply(seul, function(treg_i){
  # ----- Dataset 2: ProjecTILs projection
  print("Dataset 2...")
  treg_i <- ProjecTILs.classifier(treg_i, refobj_ds2, ncores = 1, 
                                  split.by = "orig.ident", filter.cells=F)
  treg_i$treg_dataset2 <- treg_i$functional.cluster
  
  # ----- Dataset 5: ProjecTILs projection
  print("Dataset 5...")
  treg_i <- ProjecTILs.classifier(treg_i, refobj_ds5, ncores = 1, 
                                  split.by = "orig.ident", filter.cells=F)
  treg_i$treg_dataset5 <- treg_i$functional.cluster
  
  # ----- Dataset 3: Spearman Correlation with clusters
  print("Dataset 3...")
  expr_i <- GetAssayData(treg_i, slot='data')
  genes_i <- rownames(expr_i)
  refgenes <- rownames(d3_expr)
  genes_idx <- which(genes_i %in% refgenes)
  d3_expr_ord <- d3_expr[na.omit(match(genes_i, refgenes)),]
  expr_i_ord <- expr_i[genes_idx,]
  
  cormat <- apply(expr_i_ord, 2, function(i){
    apply(d3_expr_ord, 2, function(ref){
      cor(i, ref, method='spearman', use='complete.obs')
    })
  })
  
  max_clusts <- rownames(cormat)[apply(cormat, 2, which.max)]
  treg_i$treg_dataset3 <- max_clusts
  
  return(treg_i)
})
saveRDS(seul_tregs, file=file.path(datadir, "seurat_obj", "cd45_tregs.rds"))

seul_tregs <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tregs.rds"))
treg_ds <- c(paste0('treg_dataset', c('2', '3', '3_agg', '5')), "treg_manual_anno")
seul_tregs_firstpass <- lapply(setNames(names(seul_tregs),names(seul_tregs)), function(seu_id){
  seu_i <- seul_tregs[[seu_id]]
  seu_i$old_clusters <- seu_i$seurat_clusters
  seu_i$temp_id <- paste0(seu_i$immgen.simple, seu_i$projectils.simple)
  Idents(seu_i) <- 'seurat_clusters'
  
  dims_use <- 1:20
  seu_i[['umap_orig']] <- seu_i@reductions$umap
  # seu_j <- subset(seu_i, ident=grep('Treg', as.character(unique(Idents(seu_i))), ignore.case=T, value=T))
  idents <- if(seu_id == 'LN') c('8', '19', '18') else c('19', '4', '20', '17')
  seu_j <- subset(seu_i, ident=idents)
  seu_j <- seu_j %>%
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>%
    FindNeighbors(.,dims = dims_use) %>%
    FindClusters(., resolution = 1.4) %>%
    RunUMAP(., dims=dims_use) 
  seu_j <- CellCycleScoring(seu_j, 
                            s.features = stringr::str_to_title(cc.genes$s.genes), 
                            g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                            set.ident = TRUE)
  seu_j$CC.Difference <- seu_j$S.Score - seu_j$G2M.Score
  seu_j$CC.Sum <- seu_j$S.Score + seu_j$G2M.Score
  d3_map <- setNames(paste0('cl', c('123', '123', '123', '4', '56', '56')), 
                     paste0("spleen_tr_cluster", c(1:6)))
  seu_j$treg_dataset3_agg <- d3_map[seu_j$treg_dataset3]
  
  # PacMAP clustering
  dat_pd <- reticulate::r_to_py(as.data.frame(seu_j@reductions$mnn@cell.embeddings))
  nparray <- dat_pd$values
  nparray <- nparray$astype(python_numpy$float)
  tfmr <- python_pacmap$PaCMAP(n_neighbors=10L, MN_ratio=0.5, FP_ratio=2.0) 
  X_transformed <- tfmr$fit_transform(nparray)#, init="pca")
  seu_j[["pacmap"]] <- CreateDimReducObject(embeddings = data.frame(X_transformed) %>%
                                              magrittr::set_rownames(., Cells(seu_j)) %>%
                                              as.matrix, 
                                            key = "PacMAP_", 
                                            assay = DefaultAssay(seu_j))
  return(seu_j)
})

if(do_visualize){
  pdf("~/xfer/CD45_Treg_sub.pdf", width = 13, height=13)
  gg1 <- lapply(c("umap", "pacmap"), function(reduction){
    lapply(c('old_clusters', 'seurat_clusters'), function(feat_i){
      DimPlot(seul_tregs_firstpass[['LN']], group.by=feat_i, raster=T, 
              label=T, reduction=reduction)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=2)
  gg1
  dev.off()
  
  pdf("~/xfer/CD45_Treg.pdf", width = 26)
  #Combined plots compaed to Day3 and Day7 labels
  gg0 <- lapply(c("LN", "Tumor"), function(feature){
    lapply(c('projectils.simple', 'immgen.simple', 'seurat_clusters'), function(feat_i){
      DimPlot(seul_tregs[[feature]], group.by=feat_i, raster=T, label=T, reduction='umap')
    }) %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=1)
  gg1 <- lapply(c("LN", "Tumor"), function(feature){
    lapply(c(treg_ds[-c(2, length(treg_ds))],  'old_clusters', 'seurat_clusters'), function(feat_i){
      DimPlot(seul_tregs_firstpass[[feature]], group.by=feat_i, raster=T, label=F, reduction='umap_orig')
    }) %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=1)
  gg2 <- lapply(c("LN", "Tumor"), function(feature){
    lapply(c(treg_ds[-c(2, length(treg_ds))],  'old_clusters', 'seurat_clusters'), function(feat_i){
      DimPlot(seul_tregs_firstpass[[feature]], group.by=feat_i, raster=T, label=F, reduction='umap')
    }) %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=1)
  gg3 <- lapply(c("LN", "Tumor"), function(feature){
    lapply(c(treg_ds[-c(2, length(treg_ds))],  'old_clusters', 'seurat_clusters'), function(feat_i){
      DimPlot(seul_tregs_firstpass[[feature]], group.by=feat_i, raster=T, label=F, reduction='pacmap')
    }) %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=1)
  gg4 <- lapply(c("LN", "Tumor"), function(feature){
    FeaturePlot(seul_tregs_firstpass[[feature]], features=c('CC.Sum', 'Il2ra', 'Foxp3', 'Cd3e', 'Cd4', 'Cd8a'), 
                raster=T, reduction='pacmap', combine=F) %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=1)
  gg5 <- lapply(c("LN", "Tumor"), function(feature){
    lapply(c('old_clusters', 'seurat_clusters'), function(feat_i){
      DimPlot(seul_tregs_firstpass[[feature]], group.by=feat_i, raster=T, 
              label=T, reduction='pacmap')
    }) %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., ncol=1)
  gg0
  gg1
  gg2
  gg3
  gg4
  gg5
  dev.off()
}


anno_list <- list()
treg_ds <- c(paste0('treg_dataset', c('2', '3', '3_agg', '5')), "treg_manual_anno")
process_combined <- TRUE
seul_tregs_final <- lapply(setNames(names(seul_tregs_firstpass),names(seul_tregs_firstpass)), function(seu_id){
  seu_j <- seul_tregs_firstpass[[seu_id]]
  seu_j$old_clusters <- seu_j$seurat_clusters
  dims_use <- 1:20
  
  # Manually annotate the clusters (res 1.4) and remove the problematic clusters
  seu_j$treg_manual_anno <- as.character(seu_j$seurat_clusters) %>% 
    cd45_treg_recode_map(., grp=seu_id, anno=TRUE)
  
  Idents(seu_j) <- 'treg_manual_anno'
  seu_j <- subset(seu_j, ident=setdiff(unique(Idents(seu_j)), 'Remove'))
  seu_j <- RunUMAP(seu_j, dims=dims_use, reduction='mnn')
  seu_j <- runPacmap(seu_j, reduction='mnn')
})

# Get the cell-specific treg annotations
anno_list <- lapply(seul_tregs_final, function(seu_j){
  lapply(setNames(treg_ds, treg_ds), function(i) setNames(seu_j@meta.data[,i], Cells(seu_j)))
})



saveRDS(seul_tregs_final, file=file.path(datadir, "seurat_obj", "cd45_tregs_final.rds"))
seul_tregs_final <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tregs_final.rds"))
cd45_treg_map <- lapply(seul_tregs_final, function(i) i$treg_manual_anno) %>% 
  do.call(c, .)
names(cd45_treg_map) <- gsub("^.*\\.", "", names(cd45_treg_map))
dir.create(file.path(datadir, "annotation_map"), recursive = T, showWarnings = F)
saveRDS(cd45_treg_map, file=file.path(datadir, "annotation_map", "cd45_tregs.rds"))










seul_debug2 <- lapply(setNames(names(seul_tregs),names(seul_tregs)), function(seu_id){
  seu_i <- seul_tregs[[seu_id]]
  seu_i$old_clusters <- seu_i$seurat_clusters
  seu_i$temp_id <- paste0(seu_i$immgen.simple, seu_i$projectils.simple)
  Idents(seu_i) <- 'seurat_clusters'
  
  dims_use <- 1:30
  seu_i[['umap_orig']] <- seu_i@reductions$umap
  # seu_j <- subset(seu_i, ident=grep('Treg', as.character(unique(Idents(seu_i))), ignore.case=T, value=T))
  # idents <- if(seu_id == 'LN') c('6') else c('5', '6')
  # seu_j <- subset(seu_i, ident=idents)
  seu_j <- seu_i
  seu_j <- seu_j %>%
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>%
    RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
    # RunTSNE(., dims=dims_use, reduction='mnn') %>%
    FindNeighbors(.,dims = dims_use, reduction='mnn') %>%
    FindClusters(., resolution = 1.4, graph.name='SCT_snn')
  seu_j <- CellCycleScoring(seu_j, 
                            s.features = stringr::str_to_title(cc.genes$s.genes), 
                            g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                            set.ident = TRUE)
  seu_j$CC.Difference <- seu_j$S.Score - seu_j$G2M.Score
  seu_j$CC.Sum <- seu_j$S.Score + seu_j$G2M.Score
  d3_map <- setNames(paste0('cl', c('123', '123', '123', '4', '56', '56')), 
                     paste0("spleen_tr_cluster", c(1:6)))
  seu_j$treg_dataset3_agg <- d3_map[seu_j$treg_dataset3]
  
  # PacMAP clustering
  dat_pd <- reticulate::r_to_py(as.data.frame(seu_j@reductions$mnn@cell.embeddings))
  nparray <- dat_pd$values
  nparray <- nparray$astype(python_numpy$float)
  tfmr <- python_pacmap$PaCMAP(n_neighbors=10L, MN_ratio=0.5, FP_ratio=2.0)
  X_transformed <- tfmr$fit_transform(nparray)#, init="pca")
  seu_j[["pacmap"]] <- CreateDimReducObject(embeddings = data.frame(X_transformed) %>%
                                              magrittr::set_rownames(., Cells(seu_j)) %>%
                                              as.matrix,
                                            key = "PacMAP_",
                                            assay = DefaultAssay(seu_j))
  return(seu_j)
})
if(FALSE){
  if(do_visualize){
    pdf("~/xfer/CD45_debug3.pdf", width = 12)
    gg0 <- lapply(c("LN", "Tumor"), function(feature){
      lapply(c('old_clusters', 'seurat_clusters'), function(feat_i){
        DimPlot(seul_debug2[[feature]], group.by=feat_i, raster=T, label=T, reduction='pacmap')
      }) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg1 <- lapply(c("LN", "Tumor"), function(feature){
      lapply(c('old_clusters', 'seurat_clusters'), function(feat_i){
        seu_x <- seul_debug2[[feature]]
        Idents(seu_x) <- 'old_clusters'
        seu_x <- subset(seu_x, ident=c('8', '18', '19'))
        DimPlot(seu_x, group.by=feat_i, raster=T, label=T, reduction='pacmap')
      }) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg2 <- lapply(c("LN", "Tumor"), function(feature){
      seu_x <- seul_debug2[[feature]]
      Idents(seu_x) <- 'old_clusters'
      seu_x <- subset(seu_x, ident=c('8', '18', '19'))
      
      FeaturePlot(seu_x, features=c('Il2ra', 'Foxp3', 'Cd4', 'Cd3e', 'Cd8a'), 
                  raster=T, reduction='pacmap', combine=F, pt.size = 4) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg0
    gg1
    gg2
    dev.off()
    
    pdf("~/xfer/CD45_debug.pdf", width = 20)
    #Combined plots compaed to Day3 and Day7 labels
    gg0 <- lapply(c("LN", "Tumor"), function(feature){
      lapply(c('projectils.simple', 'immgen.simple', 'seurat_clusters'), function(feat_i){
        DimPlot(seul_tregs[[feature]], group.by=feat_i, raster=T, label=T, reduction='umap')
      }) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg1 <- lapply(c("LN", "Tumor"), function(feature){
      lapply(c('old_clusters', 'seurat_clusters'), function(feat_i){
        DimPlot(seul_debug[[feature]], group.by=feat_i, raster=T, label=F, reduction='umap_orig')
      }) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg2 <- lapply(c("LN", "Tumor"), function(feature){
      lapply(c('old_clusters', 'seurat_clusters'), function(feat_i){
        DimPlot(seul_debug[[feature]], group.by=feat_i, raster=T, label=F, reduction='umap')
      }) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg4 <- lapply(c("LN", "Tumor"), function(feature){
      FeaturePlot(seul_debug[[feature]], features=c('CC.Sum', 'Cd19', 'Cd4', 'Cd3e', 'Cd8a'), 
                  raster=T, reduction='umap', combine=F, pt.size = 4) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg5 <- lapply(c("LN", "Tumor"), function(feature){
      FeaturePlot(seul_debug[[feature]], features=c('Adgre1', 'Foxp3', 'Itgax', 'Cd4', 'Cd3', 'Itgam'), 
                  raster=T, reduction='umap', combine=F, pt.size = 4) %>%
        cowplot::plot_grid(plotlist=., nrow=1)
    }) %>%
      cowplot::plot_grid(plotlist=., ncol=1)
    gg0
    gg1
    gg2
    gg4
    gg5
    dev.off()
  }
}

if(do_visualize){
  total_ggs <- lapply(names(seul_tregs_final), function(seu_id){
    seul_i <- seul_tregs_final[[seu_id]]
    day_ggs <- lapply(setNames(names(seul_i),names(seul_i)), function(day_id){
      seu_i <- seul_i[[day_id]]
      
      if(day_id != 'combined'){
        #Day3 and Day7 plots
        gg1 <- lapply(treg_ds, function(feature){
          DimPlot(seu_i, group.by=feature, raster=T, label=F, reduction='umap')
        }) %>%
          cowplot::plot_grid(plotlist=., nrow=1)
        gg2 <- lapply(treg_ds, function(feature){
          DimPlot(seu_i, group.by=feature, raster=T, label=F, reduction='pacmap')
        }) %>%
          cowplot::plot_grid(plotlist=., nrow=1)
        gg <- cowplot::plot_grid(plotlist=list(gg1, gg2), nrow=2)
      } else {
        #Combined plots compaed to Day3 and Day7 labels
        gg1 <- lapply(treg_ds, function(feature){
          lapply(paste0(c("Day3_", "Day7_", ""), feature), function(feat_i){
            print(feat_i)
            tryCatch({
              DimPlot(seu_i, group.by=feat_i, raster=T, label=F, reduction='umap')
            }, error=function(e){NULL})
          }) %>%
            cowplot::plot_grid(plotlist=., nrow=1)
        }) %>%
          cowplot::plot_grid(plotlist=., ncol=1)
        
        gg2 <- lapply(treg_ds, function(feature){
          lapply(paste0(c("Day3_", "Day7_", ""), feature), function(feat_i){
            tryCatch({
              DimPlot(seu_i, group.by=feat_i, raster=T, label=F, reduction='pacmap')
            }, error=function(e){NULL})
          }) %>%
            cowplot::plot_grid(plotlist=., nrow=1)
        }) %>%
          cowplot::plot_grid(plotlist=., ncol=1)
        gg <- list(gg1, gg2)
      }
      return(gg)
    })
    return(day_ggs)
  })
  names(total_ggs) <- names(seul_tregs_final)
  pdf("~/xfer/treg_annotations.days.pdf", width = 20)
  total_ggs$LN[c('Day3', 'Day7')]
  dev.off()
  pdf("~/xfer/treg_annotations.combined.pdf", height = 20, width = 13)
  total_ggs$LN['combined']
  dev.off()
}

if(do_deg){
  seu_i <- seul_tregs_final$LN$Day3
  Idents(seu_i) <- 'treg_manual_anno'
  celltypes <- c('NLTlike.Effector', 'NLTlike.Central')
  degs <- lapply(setNames(celltypes,celltypes), function(id){
    seu_ij <- subset(seu_i, ident=id)
    Idents(seu_ij) <- 'orig.ident'
    table(seu_ij$orig.ident)
    deg <- FindMarkers(seu_ij, group.by='orig.ident',
                       ident.1='B1_LN_KO_3d', ident.2='B1_LN_WT_3d') %>%
      tibble::rownames_to_column(., "gene") %>%
      mutate(celltype=id, 
             comparison='LNday3_KOtrx_vs_WTtrx',
             ens=gm$SYMBOL$ENSEMBL[gene],
             biotype=ens2biotype_ids[ens]) %>%
      relocate(., c(ens, biotype))
    return(deg)
  })
  for(deg_ctid in names(degs)){
    write.table(degs[[deg_ctid]], file=file.path("~/xfer", paste0("deg.LNday3.", deg_ctid, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  }
}

#--- iv.c) Correlation - Treg annotation [Dataset 3] ----
seul_tregs <- readRDS(file=file.path(datadir, "seurat_obj", "tregs.rds"))
seul_tregs <- lapply(seul_tregs, function(seu){
  seu$orig.ident <- .relabelid(seu$orig.ident)
  return(seu)
})

dataset3_expr <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/treg_ref/genesets/dataset3_expr.csv'
d3_expr <- read.csv(dataset3_expr) %>%
  as.data.frame %>%
  filter(!duplicated(SYMBOL)) %>%
  tibble::column_to_rownames(., "SYMBOL")

seul_tregs <- lapply(seul_tregs, function(seu_i){
  expr_i <- GetAssayData(seu_i, slot='data')
  genes_i <- rownames(expr_i)
  refgenes <- rownames(d3_expr)
  genes_idx <- which(genes_i %in% refgenes)
  d3_expr_ord <- d3_expr[na.omit(match(genes_i, refgenes)),]
  expr_i_ord <- expr_i[genes_idx,]
  
  cormat <- apply(expr_i_ord, 2, function(i){
    apply(d3_expr_ord, 2, function(ref){
      cor(i, ref, method='spearman', use='complete.obs')
    })
  })
  
  max_clusts <- rownames(cormat)[apply(cormat, 2, which.max)]
  seu_i$d3_anno <- max_clusts
  return(seu_i)
})
saveRDS(seul_tregs, file=file.path(datadir, "seurat_obj", "tregs.rds"))
d3_map <- setNames(paste0('cl', c('123', '123', '123', '4', '56', '56')), 
                   paste0("spleen_tr_cluster", c(1:6)))
seul_tregs[[1]]$d3_anno_agg <- d3_map[seul_tregs[[1]]$d3_anno]
seul_tregs[[2]]$d3_anno_agg <- d3_map[seul_tregs[[2]]$d3_anno]

pdf("~/xfer/dataset3_annotations.pdf", width = 20, height = 12)
lapply(seul_tregs, DimPlot, group.by='functional.cluster', 
       split.by='d3_anno', 
       raster=T, label=T, reduction='umap') %>%
  cowplot::plot_grid(plotlist=., nrow=2)

# genes <- read.csv("~/genes.csv", header = F, col.names = 'gene')
activated_markers <- list('Naive'=c('Il2ra', 'Bach2', 'Sell', 'Bcl2',
                                    'Ccr7', 'Satb1'),
                          'Effector'=c('Lag3', 'Maf', 'Klrg1', 'Pdcd1',
                                       'Tigit', 'Icos', 'Irf4', 'Gata3',
                                       'Itgae', 'Fgl2', 'Cxcr3', 'Tbx21',
                                       'Nrp1', 'Tnfrsf18', 'Ctla4', 'Lgals1',
                                       'Prdm1', 'Il10', 'Tgfb1'))
# seul_tregs <- lapply(seul_tregs, function(treg_i){
#   DefaultAssay(treg_i) <- 'RNA'
#   expr <- GetAssayData(treg_i, slot='counts')
#   idx <- which(rowSums(expr==0) < (ncol(treg_i)*0.95))
#   #--- ssGSEA from GSVA package
#   ssgsea_score = GSVA::gsva(expr = expr[idx,], activated_markers, 
#                             method = "ssgsea", 
#                             parallel.sz = 1, verbose = T)
#   treg_i$ssgsea_max <- NA
#   treg_i$ssgsea_max <-  rownames(ssgsea_score)[apply(ssgsea_score, 2, which.max)]
#   return(treg_i)
# })
# vargenes <- as.character(unlist(treg_markers))
vargenes <- as.character(unlist(activated_markers))
# vargenes <- genes$gene

lapply(c("^B1", "^B2"), function(batch){
  lapply(seul_tregs, function(seu_i){
    seu_i$old_clusters <- seu_i$seurat_clusters
    # VariableFeatures(seu_i) <- vargenes
    
    dims_use <- 1:20
    seu_j <- subset(seu_i, cells=Cells(seu_i)[grep(batch, seu_i$orig.ident)])
    seu_j <- seu_j %>%
      FindVariableFeatures %>%
      ScaleData %>%
      RunPCA %>%
      FindNeighbors(.,dims = dims_use) %>%
      FindClusters(., resolution = 0.9) %>%
      RunUMAP(., dims=dims_use) 
    markers <- FindAllMarkers(seu_j, assay = 'RNA')
    # markers2 <- FindMarkers(seu_j, assay = 'RNA', ident.1='10', ident.2='11')
    
    VariableFeatures(seu_j) <- rownames(seu_j)
    seu_j <- ScaleData(seu_j)
    multi_dimplot <- function(seu, features, group.by){
      lapply(features, function(id){
        DoHeatmap(seu, features=id, group.by=group.by)
      }) %>%
        cowplot::plot_grid(plotlist=., ncol=1)
    }
    
    seu_j$d3_anno <- gsub("spleen_tr_cluster", "cl", seu_j$d3_anno)
    pdf("~/xfer/x.pdf")
    DimPlot(seu_i, group.by='seurat_clusters', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='old_clusters', label=T, raster=T)
    multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'old_clusters')
    multi_dimplot(seu_j, activated_markers, 'old_clusters')
    multi_dimplot(seu_j, treg_markers, 'old_clusters')
    cc.genes_v <- unlist(lapply(cc.genes, str_to_title), recursive = T)
    top_markers <- markers %>% 
      filter(!gene %in% cc.genes_v) %>% 
      head(., 50) %>% pull(gene)
    
    DoHeatmap(seu_j, features=top_markers, group.by='old_clusters')
    seu_i$new_cluster <- seu_j$seurat_clusters[Cells(seu_i)]
    DimPlot(seu_i, group.by='new_cluster', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='functional.cluster', label=T, raster=T)
    DimPlot(seu_j, group.by='d3_anno', label=T, raster=T)
    DimPlot(seu_j, group.by='seurat_clusters', label=T, raster=T)
    multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'seurat_clusters')
    multi_dimplot(seu_j, activated_markers, 'seurat_clusters')
    multi_dimplot(seu_j, treg_markers, 'seurat_clusters')
    DoHeatmap(seu_j, features=top_markers, group.by='seurat_clusters')
    dev.off()
   
    
    plotl <- list("d3"=DimPlot(seu_j, group.by='d3_anno', 
                               label=T, raster=T, reduction='umap'),
                  "d3agg"=DimPlot(seu_j, group.by='ssgsea_max', 
                               label=T, raster=T, reduction='umap'),
                  "fc"=DimPlot(seu_j, group.by='functional.cluster', 
                               label=T, raster=T, reduction='umap'),
                  "clus"=DimPlot(seu_j, group.by='seurat_clusters', 
                                 label=T, raster=T, reduction='umap'),
                  "clusnew"=DimPlot(seu_j, group.by='old_clusters', 
                                    label=T, raster=T, reduction='umap'))
    cowplot::plot_grid(plotlist=plotl, nrow=1)
  }) %>% 
    cowplot::plot_grid(plotlist=., nrow=2)
})
dev.off()

for(seu_i in seul_tregs){
  for(batchid in c("Day3", "Day7")){
    seu_i$old_clusters <- seu_i$seurat_clusters
    # VariableFeatures(seu_i) <- vargenes
    batch <- c('Day3'='^B1', 'Day7'='^B2')[batchid]
    
    dims_use <- 1:20
    seu_j <- subset(seu_i, cells=Cells(seu_i)[grep(batch, seu_i$orig.ident)])
    seu_j <- seu_j %>%
      FindVariableFeatures %>%
      ScaleData %>%
      RunPCA %>%
      FindNeighbors(.,dims = dims_use) %>%
      FindClusters(., resolution = 0.9) %>%
      RunUMAP(., dims=dims_use) 
    
    markers <- FindAllMarkers(seu_j, assay = 'RNA')
    # markers2 <- FindMarkers(seu_j, assay = 'RNA', ident.1='10', ident.2='11')
    
    VariableFeatures(seu_j) <- rownames(seu_j)
    seu_j <- ScaleData(seu_j)
    multi_dimplot <- function(seu, features, group.by){
      lapply(features, function(id){
        DoHeatmap(seu, features=id, group.by=group.by)
      }) %>%
        cowplot::plot_grid(plotlist=., ncol=1)
    }
    
    seu_j$d3_anno <- gsub("spleen_tr_cluster", "cl", seu_j$d3_anno)
    pdf(paste0("~/xfer/Tregs_", batchid, ".pdf"))
    DimPlot(seu_i, group.by='seurat_clusters', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='old_clusters', label=T, raster=T)
    multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'old_clusters')
    multi_dimplot(seu_j, activated_markers, 'old_clusters')
    multi_dimplot(seu_j, treg_markers, 'old_clusters')
    cc.genes_v <- unlist(lapply(cc.genes, str_to_title), recursive = T)
    top_markers <- markers %>% 
      filter(!gene %in% cc.genes_v) %>% 
      head(., 50) %>% pull(gene)
    
    DoHeatmap(seu_j, features=top_markers, group.by='old_clusters')
    seu_i$new_cluster <- seu_j$seurat_clusters[Cells(seu_i)]
    DimPlot(seu_i, group.by='new_cluster', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='functional.cluster', label=T, raster=T)
    DimPlot(seu_j, group.by='d3_anno', label=T, raster=T)
    DimPlot(seu_j, group.by='seurat_clusters', label=T, raster=T)
    multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'seurat_clusters')
    multi_dimplot(seu_j, activated_markers, 'seurat_clusters')
    multi_dimplot(seu_j, treg_markers, 'seurat_clusters')
    DoHeatmap(seu_j, features=top_markers, group.by='seurat_clusters')
    pdf(paste0("~/xfer/Tregs_", batchid, ".2.pdf"))
    DoHeatmap(seu_j, features=c('Cd4', 'Il2ra', 'Foxp3', 'Ifng', 'Il17a'), 
              group.by='seurat_clusters')
    dev.off()
  }
}
# saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_integ.anno.rds"))

#--- v) ProjecTILs - Treg LT, NLT, NLT-like annotation ----
## Dataset 5: TReg Lymphoid, LT-like, NLT, NLT-like
# Fig 2.d: PMID 30737144
# Three populations are LT (), NLT (), NLT-like
# Have a scoring across all the cells and clusters for the NLT-like 
seul_tregs <- readRDS(file=file.path(datadir, "seurat_obj", "tregs.rds"))
seul_tregs <- lapply(seul_tregs, function(seu){
  seu$orig.ident2 <- .relabelid(seu$orig.ident)
  return(seu)
})

tcelldat <- file.path(PDIR, "..", "tcell_ref/data", '10X_all_data_preQC.rdata')
tcellmeta <- file.path(PDIR, "..", "tcell_ref/data", 'all_data_cl_metadata.csv')
load(tcelldat)
tcelldat <- UpdateSeuratObject(all_data); rm(all_data); gc()
tcellmeta <- read.csv(tcellmeta, header=T) %>%
  tibble::column_to_rownames(., "X")

cnts <- GetAssayData(tcelldat)
rownames(cnts) <- gm$ENSEMBL$SYMBOL[rownames(cnts)]
cnts <- cnts[which(!is.na(rownames(cnts))),]
tcelldat <- CreateSeuratObject(counts=cnts,
                               assay='RNA', 
                               meta.data=tcellmeta)



tcelldat <- CellCycleScoring(tcelldat, 
                             s.features = stringr::str_to_title(cc.genes$s.genes), 
                             g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                             set.ident = TRUE)
tcelldat$CC.Difference <- tcelldat$S.Score - tcelldat$G2M.Score

tcelldat <- SCTransform(tcelldat, vst.flavor = "v2", 
                        verbose = FALSE, vars.to.regress='CC.Difference') %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
saveRDS(tcelldat, file=file.path(PDIR, "..", "tcell_ref/data", "treg.rds"))



tcelldat <- readRDS(file=file.path(PDIR, "..", "tcell_ref/data", "treg.rds"))
Idents(tcelldat) <- 'cl_annot'
tcellsubdat <- subset(tcelldat, ident=unique(grep("Treg", Idents(tcelldat), value=T)))
tcellsubdat <- SCTransform(tcellsubdat, vst.flavor = "v2", 
                           verbose = FALSE, vars.to.regress='CC.Difference') %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
ref.tregsub <- ProjecTILs::make.reference(ref = tcellsubdat, ndim = 30, seed = 1234, 
                                          recalculate.umap = TRUE, annotation.column = "cl_annot")
saveRDS(ref.tregsub, file=file.path(projectils_dir, 'seu.treg_NLT_LT.rds'))
####


pdf("~/xfer/treg_lt_nlt.pdf")
DimPlot(tcellsubdat, group.by='cl_annot', reduction='umap', label=T, raster=T)
DimPlot(tcellsubdat, group.by='cl_allCells', reduction='umap', label=T, raster=T)
DimPlot(tcellsubdat, group.by='tissue_gen', reduction='umap', label=T, raster=T)
dev.off()

seul_tregs <- lapply(seul_tregs, function(treg_i){
  treg_i$ptils_treg_old2 <- treg_i$functional.cluster
  treg_i <- ProjecTILs.classifier(treg_i, ref.tregsub, ncores = 1, 
                                  split.by = "orig.ident", filter.cells=F)
  treg_i$ptils_tregnlt <- treg_i$functional.cluster
  return(treg_i)
})
pdf("~/xfer/x6.pdf")
DimPlot(seul_tregs[[1]], group.by='ptils_tregnlt', reduction='umap', label=T, raster=T)
DimPlot(seul_tregs[[2]], group.by='ptils_tregnlt', reduction='umap', label=T, raster=T)
dev.off()

for(seu_i in seul_tregs){
  for(batchid in c("Day3", "Day7")){
    seu_i$old_clusters <- seu_i$seurat_clusters
    # VariableFeatures(seu_i) <- vargenes
    batch <- c('Day3'='^B1', 'Day7'='^B2')[batchid]

    dims_use <- 1:20
    seu_j <- subset(seu_i, cells=Cells(seu_i)[grep(batch, seu_i$orig.ident)])
    # seu <- seul[[1]]
    # Idents(seu) <- 'manual_anno'
    # seu_j2 <- merge(seu_j, subset(seu, ident='CD4_CCR7hi'))
   
    
    seu_j <- seu_j %>%
      FindVariableFeatures %>%
      ScaleData %>%
      RunPCA %>%
      FindNeighbors(.,dims = dims_use) %>%
      FindClusters(., resolution =0.9) %>%
      RunUMAP(., dims=dims_use) 
    d3_map <- setNames(paste0('cl', c('123', '123', '123', '4', '56', '56')), 
                       paste0("spleen_tr_cluster", c(1:6)))
    seu_j$dataset5 <- seu_j$ptils_tregnlt
    seu_j$dataset2 <- seu_j$ptils_treg_old
    seu_j$dataset3 <- d3_map[seu_j$d3_anno]
    seu_j$dataset5[which(seu_j$dataset5 %in% paste0("Treg_", c("stress")))] <- 'NA'
    
    seu_j <- CellCycleScoring(seu_j, 
                              s.features = stringr::str_to_title(cc.genes$s.genes), 
                              g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                              set.ident = TRUE)
    seu_j$CC.Difference <- seu_j$S.Score - seu_j$G2M.Score
   
    
    if(FALSE){
      # TMP
      pdf("~/xfer/sara.day7.pdf", height = 5, width = 20)
      dp1 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='umap', label=T)
      fps <- FeaturePlot(seu_j, features ="S.Score", raster=T, reduction='umap')
      fpg2m <- FeaturePlot(seu_j, features ="G2M.Score", raster=T, reduction='umap')
      fp <- FeaturePlot(seu_j, features ="CC.Difference", raster=T, reduction='umap')
      cowplot::plot_grid(plotlist=list(dp1, fps, fpg2m, fp), nrow=1)
      
      dp1 <- DimPlot(seu_j, group.by="seurat_clusters", raster=T, reduction='umap', label=T)
      dp2 <- DimPlot(seu_j, group.by="dataset2", raster=T, reduction='umap', label=F)
      dp3 <- DimPlot(seu_j, group.by="dataset3", raster=T, reduction='umap', label=F)
      dp5 <- DimPlot(seu_j, group.by="dataset5", raster=T, reduction='umap', label=F)
      cowplot::plot_grid(plotlist=list(dp1, dp2,dp3, dp5), nrow=1)
      
      DimPlot(seu_j, group.by="dataset5", split.by='dataset5', raster=T, reduction='umap', label=F)
      
      dev.off()
    }
    
    
    markers <- FindAllMarkers(seu_j, assay = 'RNA')
    # markers2 <- FindMarkers(seu_j, assay = 'RNA', ident.1='10', ident.2='11')
    
    VariableFeatures(seu_j) <- rownames(seu_j)
    seu_j <- ScaleData(seu_j)
    multi_dimplot <- function(seu, features, group.by){
      lapply(features, function(id){
        DoHeatmap(seu, features=id, group.by=group.by)
      }) %>%
        cowplot::plot_grid(plotlist=., ncol=1)
    }
    
    seu_j$d3_anno <- gsub("spleen_tr_cluster", "cl", seu_j$d3_anno)
    pdf(paste0("~/xfer/Tregs_", batchid, ".pdf"))
    DimPlot(seu_i, group.by='seurat_clusters', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='old_clusters', label=T, raster=T)
    multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'old_clusters')
    multi_dimplot(seu_j, activated_markers, 'old_clusters')
    multi_dimplot(seu_j, treg_markers, 'old_clusters')
    cc.genes_v <- unlist(lapply(cc.genes, str_to_title), recursive = T)
    top_markers <- markers %>% 
      filter(!gene %in% cc.genes_v) %>% 
      head(., 50) %>% pull(gene)
    
    DoHeatmap(seu_j, features=top_markers, group.by='old_clusters')
    seu_i$new_cluster <- seu_j$seurat_clusters[Cells(seu_i)]
    DimPlot(seu_i, group.by='new_cluster', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='functional.cluster', label=T, raster=T)
    DimPlot(seu_j, group.by='d3_anno', label=T, raster=T)
    DimPlot(seu_j, group.by='seurat_clusters', label=T, raster=T)
    multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'seurat_clusters')
    multi_dimplot(seu_j, activated_markers, 'seurat_clusters')
    multi_dimplot(seu_j, treg_markers, 'seurat_clusters')
    DoHeatmap(seu_j, features=top_markers, group.by='seurat_clusters')
    pdf(paste0("~/xfer/Tregs_", batchid, ".2.pdf"))
    DimPlot(seu_i, group.by='new_cluster', label=T, raster=T, reduction='umap')
    DimPlot(seu_j, group.by='seurat_clusters', label=T, raster=T)
    DimPlot(seu_j, group.by='ptils_tregnlt', label=T, raster=T)
    DoHeatmap(seu_j, features=c('Cd4', 'Il2ra', 'Foxp3', 'Ifng', 'Il17a'), 
              group.by='seurat_clusters')
    dev.off()
  }
}

saveRDS(seul_tregs, file=file.path(datadir, "seurat_obj", "tregs.rds"))


#--- vi) Monocyte subset ----
dir.create(file.path(outdir, "annotation", "projectils"), showWarnings = F)
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/ProjecTILs'
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)

refobj_custom <- readRDS(file.path(projectils_dir, 'seu.treg.rds'))
DefaultAssay(refobj_custom) <- 'RNA'

seul_mono <- lapply(seul, function(seu){
  Idents(seu) <- 'manual_anno'
  DefaultAssay(seu) <- 'RNA'
  markers2 <- FindMarkers(seu, assay = 'RNA', ident.1='CD4_CCR7hi', ident.2='CD4_Tmem')
  markers2 %>%
    tibble::rownames_to_column(., "gene") %>%
    head(., 100) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    mutate(dir=avg_log2FC>0) %>%
    group_split(dir)
  
  Idents(seul[[1]]) <- 'manual_anno'
  seu <- subset(seul[[1]], ident=intersect(c('DC_CCR7hi', 'CD4_CCR7hi', 'CD4_Tmem', 'Treg'),
                              unique(seul[[1]]$manual_anno)))
  
  dims_use <- 1:30
  seu2 <- seu %>%
    FindVariableFeatures(., nfeatures=4000) %>%
    ScaleData %>%
    RunPCA %>%
    FindNeighbors(.,dims = dims_use) %>%
    FindClusters(., resolution = 0.9) %>%
    RunUMAP(., dims=dims_use) 
  # VariableFeatures(seu2) <- rownames(seu2)
  # seu2 <- ScaleData(seu2)
  pdf("~/xfer/x4.pdf")
  DimPlot(seu2, group.by='manual_anno', raster=T, label=T)
  DimPlot(seu2, group.by='seurat_clusters', raster=T, label=T)
  dev.off()
  
  markers <- FindAllMarkers(seu, assay = 'RNA')
  cc.genes_v <- unlist(lapply(cc.genes, str_to_title), recursive = T)
  markers <- markers %>%
    arrange(desc(abs(avg_log2FC)))
  top_markers <- lapply(split(markers, f=markers$cluster), head, n=20) %>%
    do.call(rbind, .) %>% 
    filter(!gene %in% cc.genes_v) %>% 
     pull(gene)
  pdf("~/xfer/x.pdf")
  DoHeatmap(seu, features = top_markers, group.by='manual_anno')
  DoHeatmap(seu, features = y$CD4_CCR7hi %>% head(50) %>% pull(gene), group.by='manual_anno')
  dev.off()
  
})

for(seu_i in seul_mono){
  for(batchid in c("Day3", "Day7")){
    seu_i$old_clusters <- seu_i$seurat_clusters
    # VariableFeatures(seu_i) <- vargenes
    batch <- c('Day3'='^B1', 'Day7'='^B2')[batchid]
    
    dims_use <- 1:20
    seu_j <- subset(seu_i, cells=Cells(seu_i)[grep(batch, seu_i$orig.ident)])
    seu_j <- seu_j %>%
      FindVariableFeatures %>%
      ScaleData %>%
      RunPCA %>%
      FindNeighbors(.,dims = dims_use) %>%
      FindClusters(., resolution = 0.9) %>%
      RunUMAP(., dims=dims_use) 
    markers <- FindAllMarkers(seu_j, assay = 'RNA')
    Idents(seu_j) <- c('old_clusters')
    markers2 <- FindMarkers(seu_j, assay = 'RNA', ident.1='11', ident.2='3')
    markers2 %>%
      tibble::rownames_to_column(., "gene") %>%
      head(., 50) %>%
      arrange(desc(abs(avg_log2FC))) %>%
      mutate(dir=avg_log2FC>0) %>%
      group_split(dir)
    
    VariableFeatures(seu_j) <- rownames(seu_j)
    seu_j <- ScaleData(seu_j)
    multi_dimplot <- function(seu, features, group.by){
      lapply(features, function(id){
        DoHeatmap(seu, features=id, group.by=group.by)
      }) %>%
        cowplot::plot_grid(plotlist=., ncol=1)
    }
    
    pdf(paste0("~/xfer/Mono_", batchid, ".pdf"))
    print(DimPlot(seu_i, group.by='seurat_clusters', label=T, raster=T, reduction='umap'))
    print(DimPlot(seu_j, group.by='old_clusters', label=T, raster=T))
    print(multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'old_clusters'))
    print(multi_dimplot(seu_j, activated_markers, 'old_clusters'))
    print(multi_dimplot(seu_j, treg_markers, 'old_clusters'))
    cc.genes_v <- unlist(lapply(cc.genes, str_to_title), recursive = T)
    top_markers <- markers %>% 
      filter(!gene %in% cc.genes_v) %>% 
      head(., 50) %>% pull(gene)
    
    print(DoHeatmap(seu_j, features=top_markers, group.by='old_clusters'))
    seu_i$new_cluster <- seu_j$seurat_clusters[Cells(seu_i)]
    print(DimPlot(seu_i, group.by='new_cluster', label=T, raster=T, reduction='umap'))
    print(DimPlot(seu_j, group.by='seurat_clusters', label=T, raster=T))
    print(multi_dimplot(seu_j, lapply(cc.genes, str_to_title), 'seurat_clusters'))
    print(multi_dimplot(seu_j, activated_markers, 'seurat_clusters'))
    print(multi_dimplot(seu_j, treg_markers, 'seurat_clusters'))
    print(DoHeatmap(seu_j, features=top_markers, group.by='seurat_clusters'))
    dev.off()
  }
}
#--- vii.a) CD4 CD8 subset - ST2 condition----
th1_ref <- c('high'=file.path(PDIR, "ref", 'Th1_ref_PMID30737144_STable3.high.csv'), 
             'low'=file.path(PDIR, "ref", 'Th1_ref_PMID30737144_STable3.low.csv'))
th1_ref <- lapply(th1_ref, read.csv, header=F)
cd48_ref <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/ref/cd4_cd8_genes.csv'
cd48_ref <- read.csv(cd48_ref, header=T, sep=",")

reftcell <- '/cluster/home/quever/downloads/gse/GSE116390/filtered'
# mtx <- Read10X(data.dir = reftcell, strip.suffix=TRUE)
# seu <- CreateSeuratObject(counts = mtx, project = grp)
spica_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/spica'
dir.create(file.path(outdir, "annotation", "projectils"), showWarnings = F)
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/ProjecTILs'
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)

# refobj_custom <- readRDS(file.path(projectils_dir, 'seu.treg.rds'))
# DefaultAssay(refobj_custom) <- 'RNA'

seul_tcell <- lapply(seul, function(seu_j){
  Idents(seu_j) <- 'manual_anno'
  DefaultAssay(seu_j) <- 'RNA'
  cd48_idents <- setdiff(grep("CD[48]", unique(seu_j$manual_anno), value=T),
                         'CD4_CCR7hi')
  seu_j <- subset(seu_j, ident=cd48_idents)
  seu_j$old_clusters <- seu_j$seurat_clusters
  
  # Preprocess and merge CD4 and CD8 cells
  dims_use <- 1:20
  seu_j <- seu_j  %>% 
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>%
    RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
    RunTSNE(., dims=dims_use, reduction='mnn') %>%
    FindNeighbors(.,dims = dims_use, reduction='mnn') %>%
    FindClusters(., resolution = 1.4, graph.name='SCT_snn')
  seu_j <- runPacmap(seu_j, reduction='mnn')
  VariableFeatures(seu_j) <- rownames(seu_j)
  seu_j <- ScaleData(seu_j)
  
  # pdf("~/xfer/th1.pdf", width = 20)
  # DoHeatmapSplit(seu_j, features=th1_ref$high[,1], split.by='manual_anno', group.by='orig.ident')
  # DoHeatmapSplit(seu_j, features=th1_ref$low[,1], split.by='manual_anno', group.by='orig.ident')
  # dev.off()
  
  # ProjecTILs annotation
  ptils_datasets <- list("titc"="../spica/ref_TILAtlas_mouse_v1.rds",
                         "vscd4"="../spica/ref_LCMV_CD4_mouse_release_v1.rds",
                         "vscd8"="../spica/ref_LCMV_Atlas_mouse_v1.rds",
                         'tcells'='seu.sara_tcells.rds')
  for(ptils_id in names(ptils_datasets)){
    refobj <- readRDS(file.path(projectils_dir, "../spica", ptils_datasets[[ptils_id]]))
    seu_j <- ProjecTILs.classifier(seu_j, refobj, ncores = 1, 
                                   split.by = "orig.ident", filter.cells=F)
    seu_j@meta.data[ptils_id] <- seu_j$functional.cluster
  }
  
  return(seu_j)
})
saveRDS(seul_tcell, file=file.path(datadir, "seurat_obj", "tcells.rds"))


seul_tcell <- readRDS(file=file.path(datadir, "seurat_obj", "tcells.rds"))
anno_list <- list()
days <- c("Day3", "Day7", 'combined')
tcell_anno <- c('titc', 'vscd4', 'vscd8', 'tcells')
seul_tcell_all <- lapply(setNames(names(seul_tcell),names(seul_tcell)), function(seu_id){
  seu_i <- seul_tcell[[seu_id]]
  lapply(setNames(days, days), function(batchid){
    seu_i$old_clusters <- seu_i$seurat_clusters
    # VariableFeatures(seu_i) <- vargenes
    batch <- c('Day3'='^B1', 'Day7'='^B2', 'combined'='combined')[batchid]
    
    dims_use <- 1:20
    if(batchid == 'combined'){
      print("Treating samples as Day 3 and 7 combined")
      seu_j <- seu_i
    } else {
      print("Treating samples Day 3 and 7 separate")
      seu_j <- subset(seu_i, cells=Cells(seu_i)[grep(batch, seu_i$orig.ident)])
      seu_j <- seu_j %>%
        FindNeighbors(.,dims = dims_use, reduction='mnn') %>%
        FindClusters(., resolution = 1.4, graph.name='SCT_snn') %>%
        RunUMAP(., dims=dims_use, reduction='mnn') 
      seu_j <- runPacmap(seu_j, reduction='mnn')
    }
    seu_j <- CellCycleScoring(seu_j, 
                              s.features = stringr::str_to_title(cc.genes$s.genes), 
                              g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
                              set.ident = TRUE)
    seu_j$CC.Sum <- seu_j$S.Score + seu_j$G2M.Score
    return(seu_j)
  })
})
seul_tcell_comb <- lapply(setNames(names(seul_tcell), names(seul_tcell)), function(seuid){
  seu_i <- seul_tcell[[seuid]]
  for(tid in c(tcell_anno, 'seurat_clusters')){
    for(day in c('Day3', 'Day7')){
      seu_i@meta.data[,paste0(day, "_", tid)] <- NA
      tid_map <- setNames(seul_tcell_all[[seuid]][[day]]@meta.data[,tid], Cells(seul_tcell_all[[seuid]][[day]]))
      seu_i@meta.data[,paste0(day, "_", tid)] <- tid_map[Cells(seu_i)]
    }
  }
  return(seu_i)
})
seul_tcell_all$LN[['combined']] <- seul_tcell_comb$LN
seul_tcell_all$Tumor[['combined']] <- seul_tcell_comb$Tumor

process_combined <- TRUE
subset_cells <- FALSE
seul_tcell_final <- lapply(setNames(names(seul_tcell_all), names(seul_tcell_all)), function(seu_id){
  dims_use <- 1:20
  seul_days <- lapply(c('Day3', 'Day7'), function(batchid){
    seu_j <- seul_tcell_all[[seu_id]][[batchid]]
    # Manually annotate the clusters (res 1.4) and remove the problematic clusters
    seu_j$tcell_manual_anno <- as.character(seu_j$seurat_clusters) %>% 
      tcell_recode_map(., grp=seu_id, day=batchid, anno=TRUE)
    
    if(subset_cells){
      Idents(seu_j) <- 'tcell_manual_anno'
      seu_j <- subset(seu_j, ident=setdiff(unique(Idents(seu_j)), 'Remove'))
      seu_j <- RunUMAP(seu_j, dims=dims_use, reduction='mnn')
      seu_j <- runPacmap(seu_j, reduction='mnn')
    }
    
    return(seu_j)
  })
  names(seul_days) <- c('Day3', 'Day7')
  
  # Get the day-specific treg annotations
  anno_list <- lapply(seul_days, function(seu_j){
    lapply(setNames('tcell_manual_anno', 'tcell_manual_anno'), function(i) setNames(seu_j@meta.data[,i], Cells(seu_j)))
  }) %>% unlist
  names(anno_list) <- gsub("^.*tcell_manual_anno.", "", names(anno_list))
  
  if(process_combined){
    print("Treating samples as Day 3 and 7 combined")
    seu_i <- seul_tcell_all[[seu_id]]$combined
    
    # Carry over the manual annotations from Day3 and 7
    seu_i$tcell_manual_anno <- 'Remove'
    cell_idx <- which(Cells(seu_i) %in% names(anno_list))
    seu_i$tcell_manual_anno[cell_idx] <- anno_list[Cells(seu_i)[cell_idx]]
    
    if(subset_cells){
      if(any(seu_i$tcell_manual_anno == 'Remove')){
        Idents(seu_i) <- 'tcell_manual_anno'
        seu_i <- subset(seu_i, ident=setdiff(unique(seu_i$tcell_manual_anno), 'Remove'))
      }
      
      seu_i <- RunUMAP(seu_i, dims=dims_use, reduction='mnn')
      seu_i <- runPacmap(seu_i, reduction='mnn')
    }
    seul_days[['combined']] <- seu_i
  }
  return(seul_days)
})
saveRDS(seul_tcell_final, file.path(datadir, "seurat_obj", "tcells_final.rds"))
st2_tcell_map <- lapply(seul_tcell_final, function(i) {
  lapply(i[c('Day3', 'Day7')], function(j) j$tcell_manual_anno) %>%
    do.call(c, .)
}) %>% 
  do.call(c, .)
names(st2_tcell_map) <- gsub("^.*\\.", "", names(st2_tcell_map))
saveRDS(st2_tcell_map, file=file.path(datadir, "annotation_map", "st2_tcells.rds"))


seul_tcell_all <- readRDS(file=file.path(datadir, "seurat_obj", "tcells_final.rds"))

tcell_dplot <- function(seu, title){
  dp1 <- DimPlot(seu, group.by='seurat_clusters', raster=T, label=T, reduction='umap') + 
    ggtitle(paste0("UMAP: " , title))
  dp2 <- DimPlot(seu, group.by='manual_anno', raster=T, label=T, reduction='umap')
  row1 <- cowplot::plot_grid(plotlist=list(dp1, dp2), nrow=1)
  row2 <- lapply(tcell_anno, function(grp){
    DimPlot(seu, group.by=grp, raster=T, label=F, reduction='umap')
  }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  umap_plot <- cowplot::plot_grid(plotlist=list(row1, row2), ncol=1)
  
  dp1 <- DimPlot(seu, group.by='seurat_clusters', raster=T, label=T, reduction='pacmap')  + 
    ggtitle(paste0("PacMAP: " , title))
  dp2 <- DimPlot(seu, group.by='manual_anno', raster=T, label=T, reduction='pacmap')
  row1 <- cowplot::plot_grid(plotlist=list(dp1, dp2), nrow=1)
  row2 <- lapply(tcell_anno, function(grp){
    DimPlot(seu, group.by=grp, raster=T, label=F, reduction='pacmap')
  }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  pacmap_plot <- cowplot::plot_grid(plotlist=list(row1, row2), ncol=1)
  
  pacmap_fplot <- FeaturePlot(seu, features = c('Cd4', 'Cd8a', 'Tbx21', 
                                                'Stat4', 'Ccr7', 'CC.sum'), 
                              raster=T, reduction='pacmap', combine=F) %>%
    cowplot::plot_grid(plotlist=., nrow = 2)
  return(list('umap'=umap_plot, 'pacmap'=pacmap_plot, "pacmap_feature"=pacmap_fplot))
}
tcell_comb_dplot <- function(seu, title){
  gg_umap <- lapply(c('Day3', 'Day7'), function(day){
    lapply(c(tcell_anno, 'seurat_clusters'), function(tid){
      DimPlot(seu, group.by=paste0(day, "_", tid), raster=T, 
              label=if(tid=='seurat_clusters') T else F, reduction='umap') + 
        ggtitle(paste0("UMAP[", title, "]: " , day, ".", tid))
    }) %>% 
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., nrow=2)
  
  gg_pacmap <- lapply(c('Day3', 'Day7'), function(day){
    lapply(c(tcell_anno, 'seurat_clusters'), function(tid){
      DimPlot(seu, group.by=paste0(day, "_", tid), raster=T, 
              label=if(tid=='seurat_clusters') T else F, reduction='pacmap') + 
        ggtitle(paste0("PacMAP[", title, "]: " , day, ".", tid))
    }) %>% 
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>%
    cowplot::plot_grid(plotlist=., nrow=2)
  
  pacmap_fplot <- FeaturePlot(seu, features = c('Cd4', 'Cd8a', 'Tbx21', 
                                                'Ccr7', 'CC.sum'), 
                              raster=T, reduction='pacmap', combine = F) %>%
    cowplot::plot_grid(plotlist=., nrow=2)
  
  return(list('umap'=gg_umap, 'pacmap'=gg_pacmap, "pacmap_feature"=pacmap_fplot))
}
pdf("~/xfer/Tcell_LN.pdf", width = 23)
tcell_dplot(seul_tcell_all$LN$Day3, title='LN_Day3')
tcell_dplot(seul_tcell_all$LN$Day7, title='LN_Day7')
tcell_dplot(seul_tcell_all$Tumor$Day3, title='Tumor_Day3')
tcell_dplot(seul_tcell_all$Tumor$Day7, title='Tumor_Day7')
tcell_comb_dplot(seul_tcell_all$LN$combined, title='LN')
tcell_comb_dplot(seul_tcell_all$Tumor$combined, title='Tumor')
dev.off()

for(seuid in names(seul_tcell_all)){
  for(day in c('Day3', 'Day7')){
    pdf(paste0("~/xfer/Tcell_heatmap.", seuid, "_", day, ".pdf"), height = 10, width = 15)
    seu <- seul_tcell_all[[seuid]][[day]]
    VariableFeatures(seu) <- rownames(seu)
    seu <- ScaleData(seu)
    seu$all <- 1
    for(gs_id in colnames(cd48_ref)){
      hm1 <- DoHeatmapSplit(seu, features = cd48_ref[,gs_id], 
                raster=T, group.by='seurat_clusters', 
                split.by='seurat_clusters', assay = 'RNA', 
                main_title=gs_id) +
        ggtitle(paste0(seuid, "_", day, ": ", gs_id))
      print(hm1)
    }
    dev.off()
    cat(paste0("xfer Tcell_heatmap.", seuid, "_", day, ".pdf\n"))
  }
}


{
  pdf("~/xfer/sara_tcells.pdf", height=5, width=17)
  dp1 <- DimPlot(seu_j, reduction='umap', group.by='functional.cluster', raster=T, label=T)
  dp2 <- DimPlot(seu_j, reduction='umap', group.by='old_clusters', raster=T, label=T)
  cowplot::plot_grid(plotlist=list(dp1, dp2), nrow=1)
  dp1 <- DimPlot(seu, reduction='umap', group.by='manual_anno', raster=T)
  dp2 <- DimPlot(seu_j, reduction='umap', group.by='seurat_clusters', raster=T, label=T)
  cowplot::plot_grid(plotlist=list(dp1, dp2), nrow=1)
  dp1 <- DimPlot(seu_j, reduction='umap', group.by='titc', raster=T, label=T)
  dp2 <- DimPlot(seu_j, reduction='umap', group.by='vscd4', raster=T, label=T)
  dp3 <- DimPlot(seu_j, reduction='umap', group.by='vscd8', raster=T, label=T)
  cowplot::plot_grid(plotlist=list(dp1, dp2, dp3), nrow=1)
  dev.off()
}

#--- vii.b) CD4 CD8 subset - CD45 condition ----
cd48_ref <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/ref/cd4_cd8_genes.csv'
cd48_ref <- read.csv(cd48_ref, header=T, sep=",")

cd48_ds1_ref <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/ref/GSE157072_Processed_D55_LCMV_RNAseq.txt.gz'
cd48_ds1_ref <- read.table(cd48_ds1_ref, header=T, sep="\t", check.names = F) %>%
  mutate(NAME = make.unique(gsub(" ", "", NAME))) %>%
  tibble::column_to_rownames(., 'NAME')

spica_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/spica'
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/ProjecTILs'
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)
do_visualize <- FALSE
do_deg <- FALSE

### Set up the Dataset data
seul_tcell <- lapply(seul, function(seu_i){
  seu_i$old_clusters <- seu_i$seurat_clusters
  seu_i$temp_id <- paste0(seu_i$immgen.simple, seu_i$projectils.simple)
  Idents(seu_i) <- 'immgen.simple'
  # Idents(seu_i) <- 'projectils.simple'
  
  dims_use <- 1:20
  seu_i[['umap_orig']] <- seu_i@reductions$umap
  
  # Preprocess and merge CD4 and CD8 cells
  seu_j <- subset(seu_i, ident=grep('^(T|CD)', as.character(unique(Idents(seu_i))), ignore.case=T, value=T))  %>% 
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>%
    RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
    RunTSNE(., dims=dims_use, reduction='mnn') %>%
    FindNeighbors(.,dims = dims_use, reduction='mnn') %>%
    FindClusters(., resolution = 1.4, graph.name='SCT_snn')
  rm(seu_i); gc()
  seu_j <- runPacmap(seu_j, reduction='mnn')
  # VariableFeatures(seu_j) <- rownames(seu_j)
  # seu_j <- ScaleData(seu_j)
  
  # ProjecTILs annotation
  ptils_datasets <- list('tcells'='seu.sara_tcells.rds')
  ptils_id = 'tcells'
  for(ptils_id in names(ptils_datasets)){
    print(ptils_id)
    refobj <- readRDS(file.path(projectils_dir, ptils_datasets[[ptils_id]]))
    seu_j <- ProjecTILs.classifier(seu_j, refobj, ncores = 1,  skip.normalize=T,
                                   split.by = "orig.ident", filter.cells=F)
    print(class(seu_j))
    seu_j@meta.data[ptils_id] <- seu_j$functional.cluster
  }
  return(seu_j)
})
rm(seul); gc()
saveRDS(seul_tcell, file=file.path(datadir, "seurat_obj", "cd45_tcells.rds"))


### Do last minute annotations on the LN, Tumor, combined datatset
seul_tcell <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tcells.rds"))
seul_tcell <- lapply(seul_tcell, function(seu_i){
  print(unique(as.character(seu_i$orig.ident)))
  expr_i <- GetAssayData(seu_i, slot='data')
  genes_i <- rownames(expr_i)
  refgenes <- rownames(cd48_ds1_ref)
  genes_idx <- which(genes_i %in% refgenes)
  d1_expr_ord <- cd48_ds1_ref[na.omit(match(genes_i, refgenes)),]
  expr_i_ord <- expr_i[genes_idx,]
  
  cormat <- apply(expr_i_ord, 2, function(i){
    apply(d1_expr_ord, 2, function(ref){
      cor(i, ref, method='spearman', use='complete.obs')
    })
  })
  
  max_clusts <- rownames(cormat)[apply(cormat, 2, which.max)]
  seu_i$d1_anno <- max_clusts
  return(seu_i)
})

saveRDS(seul_tcell, file=file.path(datadir, "seurat_obj", "cd45_tcells.rds"))

### Read in the preprocessedd Tcell and combine the LN and Tumor together
seul_tcell <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tcells.rds"))
cd45_treg_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tregs.rds"))
dims_use <- 1:30

seu_comb <- lapply(seul_tcell, RunPCA, dims=dims_use) %>%
  RunFastMNN(object.list = ., features = 2000) %>%
  RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
  FindNeighbors(.,dims = dims_use, reduction='mnn') %>%
  FindClusters(., resolution = 1.4, graph.name='RNA_snn')# %>%
  # runPacmap(., reduction='mnn')
for(seuid in names(seul_tcell)[1:2]){
  for(tid in c('tcells', 'immgen.simple', 'seurat_clusters', 'd1_anno')){
    seu_comb@meta.data[,paste0(seuid, "_", tid)] <- NA
    tid_map <- setNames(seul_tcell[[seuid]]@meta.data[,tid], Cells(seul_tcell[[seuid]]))
    seu_comb@meta.data[,paste0(seuid, "_", tid)] <- tid_map[Cells(seu_comb)]
  }
}
seul_tcell$combined <- seu_comb
seul_tcell <- lapply(seul_tcell, function(seu_i){
  seu_i <- CellCycleScoring(seu_i, s.features = str_to_title(cc.genes$s.genes), 
                                                      g2m.features = str_to_title(cc.genes$g2m.genes),
                                                      assay='SCT', set.ident = TRUE)
  seu_i$CC.Sum <- seu_i$S.Score + seu_i$G2M.Score
  cell_ids <- Cells(seu_i)[which(Cells(seu_i) %in% names(cd45_treg_map))]
  seu_i$tcells[cell_ids] <- cd45_treg_map[cell_ids]
  seu_i$tregs <- NA
  seu_i$tregs <- cd45_treg_map[Cells(seu_i)]
  return(seu_i)
})

# Rectify an anonymouse cluster
seu_i <- .monocle3Cluster(seul_tcell$combined)
Idents(seu_i) <- 'monocle3_clusters'
seu_j <- subset(seu_i, ident=c('5', '16', '13'))
seu_j$lnTumor <- c('TRUE'='LN', 'FALSE'='Tumor')[as.character(grepl("_LN_", as.character(seu_j$orig.ident)))]
seu_j2 <- seu_j %>% 
  NormalizeData(.) %>%
  FindVariableFeatures(.)  
seu_j2 <- RunFastMNN(object.list = SplitObject(seu_j2, split.by = "lnTumor")) %>%
  RunUMAP(., reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(., reduction = "mnn", dims = 1:30) %>%
  FindClusters(.)
pdf("~/xfer/anonymous_tcell_clusters.pdf", width = 15)
dp0 <- DimPlot(seu_i, group.by='monocle3_clusters', raster=T, label=T)
dp1 <- DimPlot(seul_tcell$combined, group.by='Tumor_seurat_clusters', raster=T, label=T)
dp2 <- DimPlot(seul_tcell$combined, group.by='LN_seurat_clusters', raster=T, label=T)
cowplot::plot_grid(dp0, dp1, dp2, nrow=1)

dp0a <- DimPlot(seu_j2, group.by='LN_seurat_clusters', raster=T, label=T)
dp0b <- DimPlot(seu_j2, group.by='Tumor_seurat_clusters', raster=T, label=T)
dp0c <- DimPlot(seu_j2, group.by='seurat_clusters', raster=T, label=T)
cowplot::plot_grid(dp0a, dp0b, dp0c, nrow=1)
FeaturePlot(seu_j2, features=c('Cd4', 'Cd8a', 'Foxp3', 'Cd3e'), raster=T)
dev.off()

cells_rm <- Cells(seu_j2)[seu_j2$seurat_clusters %in% c('10', '4')]
seul_tcell <- lapply(seul_tcell, function(seu_i){
  subset(seu_i, cells=setdiff(Cells(seu_i), cells_rm))
})

seu_i <- seul_tcell$combined
seu_i <- .monocle3Cluster(seu_i, partition_qval=0.1) 
cl1_cd4 <-Cells(seu_i)[which((seu_i$seurat_clusters == '4') & grepl('_LN_', seu_i$orig.ident))]
table(seu_i$seurat_clusters, seu_i$LN_tcells)['4',]
table(seu_i$seurat_clusters, seu_i$Tumor_tcells)['4',]
pdf("~/xfer/z3.pdf")
# .makeHighlightPlot(seu_i, ident = 'LN_seurat_clusters')
cd45_tcell_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tcells.rds"))
cd45_treg_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tregs.rds"))
seu_i$tcell <- cd45_tcell_map[Cells(seu_i)]
seu_i$treg <- cd45_treg_map[Cells(seu_i)]
# DimPlot(seu_i, reduction='umap', raster=T, group.by='tcell_manual_anno')
DimPlot(seu_i, reduction='umap', raster=T, group.by='tcell')
DimPlot(seu_i, reduction='umap', raster=T, group.by='treg')
DimPlot(seu_i, reduction='umap', raster=T, group.by='monocle3_clusters', label=T)
DimPlot(seu_i, reduction='umap', raster=T, group.by='seurat_clusters', label=T)
DimPlot(seu_i, reduction='umap', raster=T, group.by='old_clusters', label=T)
DimPlot(seu_i, reduction='umap', raster=T, group.by='Tumor_seurat_clusters', label=T)
DimPlot(seu_i, reduction='umap', raster=T, group.by='LN_seurat_clusters', label=T)
dev.off()
Idents(seu_i) <- 'tcell'
seu_itreg <- subset(seu_i, ident='TReg')
DimPlot(seu_itreg, reduction='umap', raster=T, group.by='tcell')
DimPlot(seu_itreg, reduction='umap', raster=T, group.by='treg')
FeaturePlot(seu_i, features=c('Cd3e', 'Foxp3', 'Cd4'), reduction='umap', raster=T)
dev.off()

### Add the final annotations to the Tcells and create a mapping
seul_tcell <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tcells.rds"))

anno_list <- list()
remove_cells <- FALSE
process_combined <- TRUE
seul_tcell_final <- lapply(c('LN', 'Tumor'), function(seu_id){
  seu_j <- seul_tcell[[seu_id]]
  seu_j$old_clusters <- seu_j$seurat_clusters
  dims_use <- 1:20
  
  # Manually annotate the clusters (res 1.4) and remove the problematic clusters
  seu_j$tcell_manual_anno <- as.character(seu_j$seurat_clusters) %>% 
    cd45_tcell_recode_map(., grp=seu_id, anno=TRUE)
  
  # pdf("~/xfer/tcell.pdf")
  # DimPlot(seu_j, group.by='seurat_clusters', label=T, reduction='umap', raster=T)
  # DimPlot(seu_j, group.by='tcell_manual_anno', label=T, reduction='umap', raster=T)
  # FeaturePlot(seu_j, features='Foxp3', reduction='umap', raster=T, slot='counts')
  # dev.off()
  
  Idents(seu_j) <- 'tcell_manual_anno'
  
  if(remove_cells){
    print("Removing cells")
    seu_j <- subset(seu_j, ident=setdiff(unique(Idents(seu_j)), 'Remove'))
    seu_j <- RunUMAP(seu_j, dims=dims_use, reduction='mnn')
    seu_j <- runPacmap(seu_j, reduction='mnn')
  }
  return(seu_j)
})
names(seul_tcell_final) <- c('LN', 'Tumor')

saveRDS(seul_tcell_final, file=file.path(datadir, "seurat_obj", "cd45_tcells_final.rds"))
seul_tcell_final <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tcells_final.rds"))
cd45_tcell_map <- lapply(seul_tcell_final, function(i) i$tcell_manual_anno) %>% 
  do.call(c, .)
names(cd45_tcell_map) <- gsub("^.*\\.", "", names(cd45_tcell_map))
dir.create(file.path(datadir, "annotation_map"), recursive = T, showWarnings = F)
saveRDS(cd45_tcell_map, file=file.path(datadir, "annotation_map", "cd45_tcells.rds"))
seul_tcell <- seul_tcell_final


if(getImmgen){
  #ref <- BlueprintEncodeData()
  #saveRDS(ref, file="~/downloads/BlueprintEncodeData.rds")
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
  bed.se <- bed.se[,setdiff(seq_along(colnames(bed.se)), grep("T.CD8", colData(bed.se)$label.fine))]
  lbl <- 'label.fine'
  clus=FALSE
  
  seul_tcell <- lapply(seul_tcell, function(seu){
    singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
                            ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            de.method="wilcox", genes='sd',
                            clusters=if(clus) seu$seurat_clusters else NULL)
    seu@meta.data[,'immgen.fine'] <- singler_anno$pruned.labels
    return(seu)
  })
} 

if(do_visualize){
  pdf("~/xfer/immgen_cd45_tcells.pdf", height = 14, width = 20)
  lapply(names(seul_tcell), function(ctype){
    seu <- seul_tcell[[ctype]]
    Idents(seu) <- 'immgen.simple'
    seu$immgen.fine <- as.character(seu$immgen.fine)
    seu_j <- tryCatch({
      seu_j <- subset(seu, ident='T.CD8')
      keep_id <- names(which(table(seu_j$immgen.fine) > 25))
      Idents(seu_j) <- 'immgen.fine'
      subset(seu_j, ident=keep_id)
    }, error=function(e){NULL})
    
    dp1 <- DimPlot(seu, group.by='immgen.simple', raster=T, reduction='umap', 
            combine = T, label=F)
    dp2 <- tryCatch({
      DimPlot(seu_j, group.by='immgen.fine', raster=T, reduction='umap', 
            combine = T, label=F)
    }, error=function(e){NULL})
    cowplot::plot_grid(dp1, dp2, nrow=1)
  }) %>% cowplot::plot_grid(plotlist=., nrow=3)
  dev.off()
}

if(do_visualize){
  pdf("~/xfer/CD45_Tcell.pdf", width = 23, height = 10)
  lapply(c('LN', 'Tumor'), function(ctype){
    lapply(c("tcells", "immgen.simple", "seurat_clusters", 'd1_anno'), function(grp){
      DimPlot(seul_tcell[[ctype]], group.by=grp, raster=T, reduction='umap', 
              combine = T, label=ifelse(grp=='seurat_clusters', T, F))
    }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  }) %>% cowplot::plot_grid(plotlist=., nrow=2)
  
  # lapply(c('LN', 'Tumor'), function(ctype){
  #   lapply(c("tcells", "immgen.simple", "seurat_clusters", 'd1_anno'), function(grp){
  #     DimPlot(seul_tcell[[ctype]], group.by=grp, raster=T, reduction='pacmap', 
  #             combine = T, label=ifelse(grp=='seurat_clusters', T, F))
  #   }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  # }) %>% cowplot::plot_grid(plotlist=., nrow=2)
  
  # lapply(c('LN', 'Tumor'), function(ctype){
  #   lapply(c("tcells", 'immgen.simple', "seurat_clusters", 'd1_anno'), function(grp){
  #     DimPlot(seul_tcell$combined, group.by=ifelse(grp=='tregs', grp, paste0(ctype, "_", grp)), 
  #             raster=T, reduction='pacmap', combine = T,
  #             label=ifelse(grp=='seurat_clusters', T, F))
  #   }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  # }) %>% cowplot::plot_grid(plotlist=., nrow=2)
  
  lapply(c('LN', 'Tumor'), function(ctype){
    lapply(c("tcells", 'immgen.simple', "seurat_clusters", 'd1_anno'), function(grp){
      DimPlot(seul_tcell$combined, group.by=ifelse(grp=='tregs', grp, paste0(ctype, "_", grp)), 
              raster=T, reduction='umap', combine = T,
              label=ifelse(grp=='seurat_clusters', T, F))
    }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  }) %>% cowplot::plot_grid(plotlist=., nrow=2)
  
  lapply(c('LN', 'Tumor'), function(ctype){
    lp <- lapply(c("tcells", 'immgen.simple', "seurat_clusters"), function(grp){
      DimPlot(seul_tcell$combined, group.by=ifelse(grp=='tregs', grp, paste0(ctype, "_", grp)), 
              raster=T, reduction='umap', combine = T,
              label=ifelse(grp=='seurat_clusters', T, F))
    }) 
    lp <- c(lp, list(DimPlot(seul_tcell$combined, group.by='seurat_clusters', 
                             raster=T, reduction='umap', combine = T,
                             label=T)))
    cowplot::plot_grid(plotlist=., nrow=1)
  }) %>% cowplot::plot_grid(plotlist=., nrow=2)
  # ggp <- lapply(seul_tcell, function(seu){
  #   FeaturePlot(seu, features = c('CC.Sum', 'Cd3e', 'Cd4', 'Cd8a', 'Foxp3'),
  #               raster=T, pt.size = 4, combine = F, reduction='pacmap') %>%
  #     cowplot::plot_grid(plotlist=., nrow=1)
  # }) %>% cowplot::plot_grid(plotlist=., nrow=3)
  # ggtitle <-  lapply(names(seul_tcell), function(main_title) ggdraw() + draw_label(label = main_title)) %>%
  #   cowplot::plot_grid(plotlist=., nrow=3)
  # cowplot::plot_grid(ggtitle, ggp, nrow=1, rel_widths = c(0.1, 1))
  
  ggp <- lapply(seul_tcell, function(seu){
    FeaturePlot(seu, features = c('CC.Sum', 'Cd3e', 'Cd4', 'Cd8a', 'Foxp3'),
                raster=T, pt.size = 4, combine = F, reduction='umap') %>%
      cowplot::plot_grid(plotlist=., nrow=1)
  }) %>% cowplot::plot_grid(plotlist=., nrow=3)
  ggtitle <-  lapply(names(seul_tcell), function(main_title) ggdraw() + draw_label(label = main_title)) %>%
    cowplot::plot_grid(plotlist=., nrow=3)
  cowplot::plot_grid(ggtitle, ggp, nrow=1, rel_widths = c(0.1, 1))
  dev.off()
}

if(do_deg){
  seu_i <- seul_tcell$LN
  temp_map <- function(x){
    x %>% 
      recode("0"='0.1.5.6',
             "1"='0.1.5.6',
             "5"='0.1.5.6',
             "6"='0.1.5.6')
  }
  seu_i$temp_clus <- temp_map(seu_i$seurat_clusters)
  Idents(seu_i) <- 'temp_clus'
  deg_all <- FindMarkers(seu_i, group.by='temp_clus',
                         ident.1='0.1.5.6') %>%
    tibble::rownames_to_column(., "gene") %>%
    mutate(celltype='Tumor_LN', 
           comparison='0.1.5.6_vs_all',
           ens=gm$SYMBOL$ENSEMBL[gene],
           biotype=ens2biotype_ids[ens]) %>%
    relocate(., c(ens, biotype))
  write.table(deg_all, file=file.path("~/xfer", "cd45_ln_tcells_deg.csv"),
              sep=",", col.names = T, row.names = F, quote = F)
  
  seu_i <- seul_tcell$Tumor
  temp_map <- function(x){
    x %>% 
      recode("1"='1.22.5',
             "22"='1.22.5',
             "5"='1.22.5',
             "10"='10.7.14',
             "7"='10.7.14',
             "14"='10.7.14',
             .default=NA_character_)
  }
  seu_i$temp_clus <- temp_map(seu_i$seurat_clusters)
  deg <- FindMarkers(seu_i, group.by='temp_clus',
                     ident.1='1.22.5', ident.2='10.7.14') %>%
    tibble::rownames_to_column(., "gene") %>%
    mutate(celltype='Tumor_tcells', 
           comparison='1.22.5_vs_10.7.14',
           ens=gm$SYMBOL$ENSEMBL[gene],
           biotype=ens2biotype_ids[ens]) %>%
    relocate(., c(ens, biotype))
  write.table(deg, file=file.path("~/xfer", "cd45_tumor_tcells_deg.csv"),
              sep=",", col.names = T, row.names = F, quote = F)
  
  
  seu_i <- seul_tcell$combined
  temp_map <- function(x, mode=NULL){
    if(grepl("tumor", mode, ignore.case = T)){
      print("Relabelling tumor")
      x %>% 
        recode("1"='1.22.5',
               "22"='1.22.5',
               "5"='1.22.5',
               "10"='10.7.14',
               "7"='10.7.14',
               "14"='10.7.14',
               .default=NA_character_)
    } else {
      print("Relabelling LN")
      x %>% 
        recode("0"='0.1.5.6',
               "1"='0.1.5.6',
               "5"='0.1.5.6',
               "6"='0.1.5.6')
    }
  }
  tumor_id <- as.character(temp_map(seu_i$Tumor_seurat_clusters, mode='Tumor'))
  ln_id <- as.character(temp_map(seu_i$LN_seurat_clusters, mode='LN'))
  ln_id[which(is.na(ln_id))] <- tumor_id[which(is.na(ln_id))]
  seu_i$temp_clus <- ln_id
  deg <- FindMarkers(seu_i, group.by='temp_clus',
                     ident.1='0.1.5.6', ident.2='1.22.5') %>%
    tibble::rownames_to_column(., "gene") %>%
    mutate(celltype='Tumor_LN_tcells', 
           comparison='0.1.5.6_vs_1.22.5',
           ens=gm$SYMBOL$ENSEMBL[gene],
           biotype=ens2biotype_ids[ens]) %>%
    relocate(., c(ens, biotype))
  write.table(deg, file=file.path("~/xfer", "cd45_tumor_ln_tcells_deg.csv"),
              sep=",", col.names = T, row.names = F, quote = F)
  
}


pdf("~/xfer/CD45_Tcell.pdf", width = 20, height = 10)
ggp <- lapply(seul_tcell, function(seu){
  # cowplot::plot_grid(ggtitle, ggp, nrow=1, rel_widths = c(0.1, 1))
  clusters <- unique(as.character(sort(as.integer(seu$seurat_clusters))))
  lapply(clusters, function(cl){
    Idents(seu) <- 'seurat_clusters'
    clhigh <- Cluster_Highlight_Plot(seurat_object = seu, cluster_name = as.character(cl), 
                           highlight_color = "red",
                           background_color = "lightgray",
                           reduction='pacmap') + ggtitle(cl) + NoLegend()
    fp <- FeaturePlot(seu, features = c('CC.Sum', 'Cd3e', 'Cd4', 'Cd8a', 'Foxp3'),
                raster=T, pt.size = 4, combine = F, reduction='pacmap')
    fpl <- lapply(fp, function(ggf){
      ggf + 
        geom_point(data=clhigh$data, 
                   aes(x=PacMAP_1, y=PacMAP_2, fill=highlight), alpha=0.2)
    })
    
    pdf("~/xfer/x.pdf", width = 20)
    # cowplot::plot_grid(plotlist=c(list(clhigh), fp), nrow=1)
    cowplot::plot_grid(plotlist=fpl, nrow=1)
    # fpl
    dev.off()
  })
}) 


dev.off()

for(seuid in rev(names(seul_tcell))){
  pdf(paste0("~/xfer/CD45_Tcell_heatmap.", seuid, ".pdf"), height = 10, width = 15)
  seu <- seul_tcell[[seuid]]
  VariableFeatures(seu) <- rownames(seu)
  seu <- ScaleData(seu)
  for(gs_id in colnames(cd48_ref)){
    hm1 <- DoHeatmapSplit(seu, features = cd48_ref[,gs_id],
                          raster=T, group.by='seurat_clusters',
                          split.by='seurat_clusters', assay = 'RNA', 
                          main_title=gs_id) +
      ggtitle(paste0("CD45_", seuid, ": ", gs_id))
    # hm1 <- DoHeatmap(seu, features = cd48_ref[,gs_id], 
    #                       raster=T, group.by='seurat_clusters',assay = 'RNA') +
    #   ggtitle(paste0("CD45_", seuid, ": ", gs_id))
    print(hm1)
  }
  dev.off()
  cat(paste0("xfer CD45_Tcell_heatmap.", seuid, ".pdf\n"))
}



#--- viii.a) SingleR CD8 annotations ----
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
cd8refdir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/cd8_ref/results'
cd8_vsd <- read.table(file.path(cd8refdir, "vst_counts_merged.csv"), sep=",", header=T, stringsAsFactors = F)
# bed.se2 <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
bed.se <- SummarizedExperiment(log(cd8_vsd+1), colData=data.frame("label"=colnames(cd8_vsd)))
assayNames(bed.se) <- 'logcounts'
variable_genes <- TRUE
singler_annos <- lapply(seul, function(seu){
  # lapply(setNames(c(TRUE, FALSE), c("clusters", "cells")), function(clus){
    Idents(seu) <- 'manual_anno'
    cd8s <- unique(grep("cd8", seu$manual_anno, ignore.case = T, value=T))
    seu_cd8 <- subset(seu, ident=cd8s)
    
    ## Run a simpler spearman correlation between the normalized CD8 matrix and seurat data
    expr_i <- GetAssayData(seu_cd8, slot='data')
    if(variable_genes){
      seu_cd8 <- FindVariableFeatures(seu_cd8)
      genes <- seu_cd8@assays$RNA@var.features
      expr_i <- expr_i[genes,]
    }
    genes_i <- rownames(expr_i)
    refgenes <- rownames(cd8_vsd)
    genes_idx <- which(genes_i %in% refgenes)
    cd8_vsd_ord <- cd8_vsd[na.omit(match(genes_i, refgenes)),]
    expr_i_ord <- expr_i[genes_idx,]
    
    cormat <- apply(expr_i_ord, 2, function(i){
      apply(cd8_vsd_ord, 2, function(ref){
        cor(i, ref, method='spearman', use='complete.obs')
      })
    })
    
    max_clusts <- rownames(cormat)[apply(cormat, 2, which.max)]
    seu_cd8$cd8_anno <- max_clusts
    id_map <- setNames(seu_cd8$cd8_anno, colnames(seu_cd8))
    seu$cd8_anno <- NA
    seu@meta.data[,'cd8_anno'] <- id_map[colnames(seu)]
    return(seu)
    
    # SingleR works if you have multiple samples per celltype
    # singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
    #                         ref=bed.se, 
    #                         assay.type.test=1, labels=bed.se$label,
    #                         de.method="wilcox", genes='sd',
    #                         clusters=if(clus) seu$seurat_clusters else NULL)
    # return(singler_anno)
  # })
})

pdf(paste0("~/xfer/cd8_anno.", ifelse(variable_genes, "2kvargenes", "allGenes"), ".pdf"), 
    height = 12, width=17)
lapply(singler_annos, function(seu){
  Idents(seu) <- 'manual_anno'
  ids <- unique(grep("cd8", seu$manual_anno, ignore.case = T, value=T))
  seusub <- subset(seu, ident=ids)
  # seusub <- seusub %>%
  #   FindVariableFeatures(.)  %>%
  #   ScaleData(.)  %>%
  #   RunPCA(.)  %>%
  #   RunUMAP(., dims = 1:20, reduction = "pca")
  seusub2 <- seusub
  seusub2$cd8_anno <- gsub("^(CD8_|IEL.|Sp.|Sp.Conv.)", "", seusub$cd8_anno) %>%
    gsub("(_30|.1|.2)$", "", .) %>% tolower
  
  lapply(c('manual_anno', 'cd8_anno'), function(clid){
    dp1 <- scCustomize::DimPlot_scCustom(
      seu, group.by=clid, figure_plot=T, reduction='umap', raster=F, repel=T, label=T
    )
    dp2 <- scCustomize::DimPlot_scCustom(
      seusub, group.by=clid, figure_plot=T, reduction='umap', raster=F, repel=T, label=T
    )
    dp3 <- scCustomize::DimPlot_scCustom(
      seusub2, group.by=clid, figure_plot=T, reduction='umap', raster=F, repel=T, label=T
    )
    cowplot::plot_grid(dp1, dp2, dp3, nrow=1)
  }) %>% 
    cowplot::plot_grid(plotlist=., ncol=1)
  # print(.makeHighlightPlot(seusub, ident='cd8_anno'))
})
dev.off()
cat(paste0("~/xfer/cd8_anno.", ifelse(variable_genes, "2kvargenes", "allGenes"), ".pdf\n"))
  
lapply(singler_annos, function(seu){
  tbl <- table(seu$manual_anno, seu$cd8_anno)
  tbl[grep("cd8", rownames(tbl), ignore.case = T, value=T),] %>%
    write.table(., sep=",", quote = F, col.names = T, row.names = T)
})


## DEG of each CD8 subpopulation
deg_cd8s <- lapply(names(seul), function(seu_id){
  seu <- seul[[seu_id]]
  # lapply(setNames(c(TRUE, FALSE), c("clusters", "cells")), function(clus){
  Idents(seu) <- 'manual_anno'
  cd8s <- unique(grep("^cd8", seu$manual_anno, ignore.case = T, value=T))
  seu_cd8 <- subset(seu, ident=cd8s)
  subset_of_seu <- T
  
  lapply(setNames(cd8s,cd8s), function(cd8_i){
    FindMarkers(seu_cd8, assay = "RNA", test.use='wilcox',
                ident.1= cd8_i, 
                verbose = FALSE,
                logfc.threshold = 0,
                recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
      tibble::rownames_to_column(., "symbol") %>% 
      mutate(biotype = sym2biotype_ids[symbol],
             ensemble = gm$SYMBOL$ENSEMBL[symbol],
             ident.1 = cd8_i, 
             tissue=seu_id) %>% 
      relocate(., c('ensemble', 'biotype'))
  })
})
names(deg_cd8s) <- names(seul)

lapply(names(deg_cd8s), function(deg_id){
  do.call(rbind, deg_cd8s[[deg_id]]) %>% 
    write.table(., file=paste0("~/xfer/cd8.", deg_id, ".csv"),
                sep=",", col.names = T, row.names = F, quote = F)
})
saveRDS(deg_cd8s, file=file.path(outdir, "deg_cd8s.rds"))
  
#--- viii.c) CD45 Cells - Fixing TReg and TCells ----

seul_tregs_final <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tregs_final.rds"))
seul_tcell_final <- readRDS(file=file.path(datadir, "seurat_obj", "cd45_tcells_final.rds"))

seul <- readRDS(file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.rds"))
seu_all_cd8t <- list()

lt <- 'Tumor'
seul[[lt]]$manual_anno <- cd45_recode_map(seul[[lt]]$seurat_clusters, grp=lt)
cd8cells <- Cells(seul[[lt]])[seul[[lt]]$manual_anno == 'CD8']
seu_all_cd8t[[lt]] <- subset(seul[[lt]], cells=cd8cells)
lt <- 'LN'
seul[[lt]]$manual_anno <- cd45_recode_map(seul[[lt]]$seurat_clusters, grp=lt)
cd8cells <- Cells(seul[[lt]])[seul[[lt]]$manual_anno == '22']
seu_all_cd8t[[lt]] <- subset(seul[[lt]], cells=cd8cells)
rm(seul); gc()

# Merge the TReg and TCell cells, then find UMAP and clusters
seulx <- lapply(c('LN', 'Tumor'), function(i){
  cells_tcell <- Cells(seul_tcell_final[[i]])
  cells_treg <- Cells(seul_tregs_final[[i]])
  tcell_only <- tryCatch({seul_tcell_final[[i]][,setdiff(cells_tcell, cells_treg)]}, error=function(e){NULL})
  treg_only <- tryCatch({seul_tregs_final[[i]][,setdiff(cells_treg, cells_tcell)]}, error=function(e){NULL})
  extra_cd8 <- tryCatch({seu_all_cd8t[[i]][,setdiff(Cells(seu_all_cd8t[[i]]), c(cells_treg, cells_tcell))]}, error=function(e){NULL})
  tcell_intersect <- seul_tcell_final[[i]][,intersect(cells_tcell, cells_treg)]
  tregs_intersect <- seul_tregs_final[[i]][,intersect(cells_tcell, cells_treg)]
  tcell_intersect$treg_manual_anno <- tregs_intersect$treg_manual_anno
  cells <- list(tcell_only, treg_only, tcell_intersect, extra_cd8)
  if(any(sapply(cells, is.null))) cells <- cells[-which(sapply(cells, is.null))]
  
  merge(cells[[1]], cells[-1]) %>% 
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>%
    RunUMAP(., dims=1:30, n.neighbors = 20L) %>%
    FindNeighbors(.,dims = 1:30) %>%
    FindClusters(., resolution = 1.2)
})
names(seulx) <-c('LN', 'Tumor')

# cd22_map_immgen <- setNames(seu$immgen, Cells(seu))
# cd22_map_projectils <- setNames(seu$functional.cluster, Cells(seu))
# idx <- which(seulx$LN$manual_anno == '22')
# seulx$LN$tcells[idx] <- cd22_map_projectils[Cells(seulx$LN)[idx]]
# seulx$LN$immgen.simple[idx] <- cd22_map_immgen[Cells(seulx$LN)[idx]]
# newid <- gsub("^.*\\(", "", cd22_map_immgen) %>% gsub ("\\)$", "", .) %>% toupper(.)
# newid[newid %in% names(which(table(newid) <= 10))] <- NA
# seulx$LN$immgen.simple[idx] <- newid[Cells(seulx$LN)[idx]]
# 
# pdf("~/xfer/cluster22.pdf")
# Idents(seulx$LN) <- 'tcell_manual_anno'
# seu_x <- subset(seulx$LN, ident='TReg')
# head(seu_x)
# DimPlot(seulx$LN, group.by='tcells', raster=T, reduction='umap', label=T) + ylim(-10, 15) + xlim(-10,10)
# DimPlot(seulx$LN, group.by='immgen.simple', raster=T, reduction='umap', label=T) + ylim(-10, 15) + xlim(-10,10)
# DimPlot(seu_x, group.by='tcells', raster=T, reduction='umap', label=T) + ylim(-10, 15) + xlim(-10,10)
# DimPlot(seu_x, group.by='immgen.simple', raster=T, reduction='umap', label=F) + ylim(-10, 15) + xlim(-10,10)
# DimPlot(seu_x, group.by='immgen.simple', raster=T, reduction='umap', label=F, pt.size=5)
# dev.off()


# Identity clusters where the TCell annotation states there's a TReg, but the 
# TReg annotation is not annotated.  In the LN, there's one cluster of TCells
# labelled as "TReg", but is not Foxp3+, while the other "TReg" cells are scattered
# across multiple TReg clusters
seulx2 <- lapply(seulx, function(seu_i){
  min_size = 10
  set.seed(1234)
  seu_i <- .monocle3Cluster(seu_i, resolution=0.015)
  # pdf("~/xfer/y.pdf", width=15)
  # DimPlot(seu_i, group.by='monocle3_clusters', reduction='umap', label=T)
  # DimPlot(seu_i, group.by='manual_anno', reduction='umap', label=T)
  # DimPlot(seu_i, group.by='tcell_manual_anno', reduction='umap', label=T)
  # DimPlot(seu_i, group.by='treg_manual_anno', reduction='umap', label=T)
  # dev.off()
  clusters <- 'monocle3_clusters'
  # identify seurat_clusters that are populated by the TCell "TRegs"
  treg_clusters <- table(seu_i$tcell_manual_anno, seu_i@meta.data[,clusters])['TReg',]
  treg_cls <- names(which(treg_clusters > min_size))
  # identify which TReg populations populate those selected seurat_clusters
  seurat_clusters <- table(seu_i$treg_manual_anno, seu_i@meta.data[,clusters])[,treg_cls]
  keep_idx <- colSums(seurat_clusters) > min_size
  seurat_labels <- unlist(apply(seurat_clusters[,keep_idx], 2, function(i) names(which.max(i))))
  # Use the max TReg populations per cluster to label the TRegs and TCells appropriately
  TReg_idx <- which((seu_i$tcell_manual_anno == 'TReg') & 
                      (seu_i@meta.data[,clusters] %in% names(seurat_labels)) & 
                      (is.na(seu_i$treg_manual_anno)))
  seu_i$treg_manual_anno[TReg_idx] <- seurat_labels[as.character(seu_i@meta.data[,clusters][TReg_idx])]
  seu_i$tcell_manual_anno[TReg_idx] <- NA
  
  seu_i$tcell_manual_anno[which(seu_i$tcell_manual_anno == 'TReg')] <- 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_001'
  seu_i$tcell_manual_anno[which((seu_i$manual_anno == 'CD8') & (is.na(seu_i$tcell_manual_anno)))] <- 'CD8_memory'
  seu_i$tcell_manual_anno[which(seu_i$tcell_manual_anno == 'CD4_CCR7pos_TCF7pos')] <- 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002'
  
  # harcoded segment
  seu_i$treg_manual_anno[which(grepl("_LN_", seu_i$orig.ident) & 
                                  seu_i$monocle3_clusters %in% c('165', '150', '121', '85', '173', '184'))] <- 'NLTlike_cycling'
  seu_i$tcell_manual_anno[which(grepl("_LN_", seu_i$orig.ident) & 
                                 seu_i$monocle3_clusters %in% c('165', '150', '121', '85', '173', '184'))] <- 'NLTlike_cycling'
  seu_i$tcell_manual_anno[which((seu_i$monocle3_clusters == '53') & (is.na(seu_i$tcell_manual_anno)))] <- 'CD8_naive2'
  seu_i$tcell_manual_anno[which((seu_i$monocle3_clusters == '151') & (is.na(seu_i$tcell_manual_anno)))] <- 'CD4_naive2'
  return(seu_i)
})

## Validation check
run_validation_check <- FALSE
if(run_validation_check){
  seulx_treg <- lapply(seulx, function(seu_i){
    cells <- Cells(seu_i)[which((seu_i$tcell_manual_anno == 'TReg') & is.na(seu_i$treg_manual_anno))]
    subset(seu_i, cells=cells)
  })
  
  
  pdf("~/xfer/tcell_treg.pdf", width = 13)
  lapply(seulx_treg, function(seu_i){
    FeaturePlot(seu_i, features=c('Foxp3', 'Cd3e', 'Cd4', 'Cd8a'), reduction='umap', raster=T, slot='counts', combine=T)
  })
  
  
  lapply(seulx, function(seu_i){
    l <- list(DimPlot(seu_i, group.by='treg_manual_anno', reduction='umap', raster=T),
              DimPlot(seu_i, group.by='tcell_manual_anno', reduction='umap', raster=T),
              DimPlot(seu_i, group.by='seurat_clusters', reduction='umap', raster=T, label=T))
    cowplot::plot_grid(plotlist=l, nrow=1)
  }) %>% 
    cowplot::plot_grid(plotlist=., nrow=2)
  lapply(seulx, function(seu_i){
    FeaturePlot(seu_i, features=c('Foxp3', 'Cd3e', 'Cd4', 'Cd8a'), reduction='umap', raster=T, slot='counts')
  })
  dev.off()
}

## Save the updated maps
cd45_treg_map <- lapply(seulx, function(i) i$treg_manual_anno) %>% 
  do.call(c, .)
names(cd45_treg_map) <- gsub("^.*\\.", "", names(cd45_treg_map)) %>% 
  gsub("_[12]$", "", .)

cd45_tcell_map <- lapply(seulx, function(i) i$tcell_manual_anno) %>% 
  do.call(c, .)
names(cd45_tcell_map) <- gsub("^.*\\.", "", names(cd45_tcell_map))

saveRDS(cd45_treg_map, file=file.path(datadir, "annotation_map", "cd45_tregs.rds"))
saveRDS(cd45_tcell_map, file=file.path(datadir, "annotation_map", "cd45_tcells.rds"))
saveRDS(seulx, file=file.path(datadir, "seurat_obj", "cd45_tregs_tcells_final.rds"))



#--- viii.a) ST2 - total cell plot ----
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.rds"))
st2_tcell_map <- readRDS(file=file.path(datadir, "annotation_map", "st2_tcells.rds"))
st2_treg_map <- readRDS(file=file.path(datadir, "annotation_map", "st2_tregs.rds"))
overwrite <- FALSE

for(lt in c('LN', 'Tumor')){
  seul[[lt]]$treg_manual_anno <- paste0("TReg_", st2_treg_map[Cells(seul[[lt]])])
  seul[[lt]]$treg_manual_anno[grep("_NA$", seul[[lt]]$treg_manual_anno)] <- NA
  seul[[lt]]$tcell_manual_anno <- paste0("T_", st2_tcell_map[Cells(seul[[lt]])])
  seul[[lt]]$tcell_manual_anno[grep("_NA$", seul[[lt]]$tcell_manual_anno)] <- NA
  
  cell_idx <- which(Cells(seul[[lt]]) %in% names(st2_treg_map))
  seul[[lt]]$manual_anno <- as.character(seul[[lt]]$manual_anno)
  seul[[lt]]$manual_anno[cell_idx] <- paste0("TReg_", st2_treg_map[Cells(seul[[lt]])[cell_idx]])
  
  cell_idx <- which(Cells(seul[[lt]]) %in% names(st2_tcell_map))
  seul[[lt]]$manual_anno <- as.character(seul[[lt]]$manual_anno)
  seul[[lt]]$manual_anno[cell_idx] <- paste0("T_", st2_tcell_map[Cells(seul[[lt]])[cell_idx]])
  seul[[lt]]$manual_anno[grep("_NA$", seul[[lt]]$manual_anno)] <- NA
  
  seul[[lt]]  <- seul[[lt]] %>%
    .monocle3Cluster(.) 
}

# pdf("~/xfer/st2.pdf", width=20)
# DimPlot(seul[['LN']], group.by='monocle3_clusters', reduction='umap', raster=T, label=T) +
# DimPlot(seul[['LN']], group.by='manual_anno', reduction='umap', raster=T, label=T)
# DimPlot(seul[['Tumor']], group.by='monocle3_clusters', reduction='umap', raster=T, label=T) +
# DimPlot(seul[['Tumor']], group.by='manual_anno', reduction='umap', raster=T, label=T)
# dev.off()

pdf("~/xfer/st2.pdf", width=15)
DimPlot(seul[['LN']], group.by='manual_anno', reduction='umap', raster=T, label=T)
DimPlot(seul[['Tumor']], group.by='manual_anno', reduction='umap', raster=T, label=T)
dev.off()

table(seul$Tumor$manual_anno, seul$Tumor$monocle3_clusters)
seul$LN$manual_anno <- gsub("TReg_NLTlike.Effector.Central", "TReg_NLTlike.Effector.Central_001", seul$LN$manual_anno)
idx <- which((seul$LN$manual_anno == 'TReg_NLTlike.Effector.Central_001') & 
               (seul$LN$monocle3_clusters == '10'))
seul$LN$manual_anno[idx] <- 'TReg_NLTlike.Effector.Central_002'

table(seul$Tumor$manual_anno, seul$Tumor$monocle3_clusters)
seul$Tumor$manual_anno <- gsub("TReg_NLTlike.Effector.Central", "TReg_NLTlike.Effector.Central_001", seul$Tumor$manual_anno)
idx <- which((seul$Tumor$manual_anno == 'TReg_NLTlike.Effector.Central_001') & 
               (seul$Tumor$monocle3_clusters %in% c('6', '18')))
seul$Tumor$manual_anno[idx] <- 'TReg_NLTlike.Effector.Central_002'

seul$Tumor$manual_anno <- gsub("TReg_NLTlike.STAT1", "TReg_NLTlike.STAT1_001", seul$Tumor$manual_anno)
idx <- which((seul$Tumor$manual_anno == 'TReg_NLTlike.STAT1_001') & 
               (seul$Tumor$monocle3_clusters %in% c('6', '17')))
seul$Tumor$manual_anno[idx] <- 'TReg_NLTlike.STAT1_002'


# if(overwrite) saveRDS(seul, file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
# seul <- readRDS(file=file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))

seul_clean <- cleanLnTumorTreg(seul)
seul_clean <- lapply(names(seul_clean), function(seuid){
  seu_i = seul_clean[[seuid]]
  cell_idx <- (seu_i$manual_anno %in% c(unique(grep("Remove", seu_i$manual_anno, value=T)), 
                                        c('CD4_CCR7hi', 'Patrolling_monocyte')))
  return(subset(seu_i, cells=Cells(seu_i)[which(!cell_idx)]))
}) %>%
  setNames(., names(seul))
if(overwrite) saveRDS(seul_clean, file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
seul_clean <- readRDS(file=file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))

if(do_visualize){
  pdf("~/xfer/ST2_ln.pdf", width = 7,height = 7)
  print(.makeHighlightPlot(seul_clean$LN, ident='manual_anno'))
  dev.off()
  pdf("~/xfer/ST2_tumor.pdf", width = 7,height = 7)
  print(.makeHighlightPlot(seul_clean$Tumor, ident='manual_anno'))
  dev.off()
  
  pdf("~/xfer/ST2_total_cells_with_tregs_tcells.2.pdf", width = 32,height = 15)
  ggdim_seul <- lapply(names(seul), function(seuid){
    ggdim0 <- lapply(c('treg_manual_anno', 'tcell_manual_anno', 'manual_anno', 'seurat_clusters'), function(grp){
      # Relabel the clusters from Treg to 1_Treg to simplify the overlay
      seu_i = seul_clean[[seuid]]
      
      grp_ids <- na.omit(unique(seu_i@meta.data[,grp]))
      grp_map <- setNames(seq_along(grp_ids), grp_ids)
      legend_map <- setNames(paste(grp_map, names(grp_map), sep="_"), grp_map)
      
      seu_i@meta.data[,grp] <- factor(grp_map[seu_i@meta.data[,grp]])
      labels <- levels(seu_i@meta.data[,grp])
      
      DimPlot(seu_i, group.by=grp, raster=T, label=T, reduction='umap') + 
        scale_color_discrete(labels=legend_map[labels])
    }) %>% 
      cowplot::plot_grid(plotlist=., nrow=1)
    ggseuid <- ggdraw() + draw_label(label = seuid)
    cowplot::plot_grid(ggseuid, ggdim0, nrow=1, rel_widths = c(0.2, 1))
  }) %>%
    cowplot::plot_grid(plotlist=., nrow=2)
  ggdim_seul
  dev.off()
  
  
  
  grp <- 'manual_anno'
  pdf("~/xfer/st2_publication.pdf", width=7, height=5)
  grps <- sapply(seul_clean, function(seu_i) unique(seu_i@meta.data[,grp])) %>%
    unlist %>% 
    unique
  colors_use <- scCustomize::scCustomize_Palette(num_groups = length(grps), 
                                                 ggplot_default_colors = FALSE, 
                                                 color_seed = 231) %>%
    setNames(., grps)
  lapply(seul_clean, function(seu_i){
    colors_use_i <- as.character(colors_use[unique(seu_i@meta.data[,grp])])
    publicationDimPlot(seu_i, grp=grp, simplify_labels=TRUE, 
                       colors_use=colors_use_i, pt.size=0.2, aspect.ratio=1,
                       legend.text = element_text(size=10))
  })
  dev.off()
  
  pdf("~/xfer/st2_treg_publication.pdf", width=7, height=5)
  lapply(seul_clean, function(seu_i){
    cells <- Cells(seu_i)[grep("^TReg", seu_i@meta.data[,grp])]
    seu_i <- subset(seu_i, cells=cells)
    colors_use_i <- as.character(colors_use[unique(seu_i@meta.data[,grp])])
    publicationDimPlot(seu_i, grp=grp, simplify_labels=TRUE, colors_use=colors_use_i, 
                       pt.size=0.2, aspect.ratio=1,
                       legend.text = element_text(size=10))
  })
  dev.off()
 
}

#--- viii.a) CD45 - total cell plot ----
if(!exists("seul"))  seul <- readRDS(file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.rds"))
cd45_tcell_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tcells.rds"))
cd45_treg_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tregs.rds"))
overwrite <- FALSE


rm_removecells <- T
for(lt in c('LN', 'Tumor')){
  seul[[lt]]$manual_anno <- NA
  seul[[lt]]$manual_anno <- cd45_recode_map(seul[[lt]]$seurat_clusters, grp=lt)
  
  if(any(is.na(cd45_treg_map))) cd45_treg_map <- cd45_treg_map[which(!is.na(cd45_treg_map))]
  if(any(is.na(cd45_tcell_map))) cd45_tcell_map <- cd45_tcell_map[which(!is.na(cd45_tcell_map))]
  seul[[lt]]$treg_manual_anno <- paste0("TReg_", cd45_treg_map[Cells(seul[[lt]])])
  seul[[lt]]$treg_manual_anno[grep("_NA$", seul[[lt]]$treg_manual_anno)] <- NA
  seul[[lt]]$tcell_manual_anno <- paste0("T_", cd45_tcell_map[Cells(seul[[lt]])])
  seul[[lt]]$tcell_manual_anno[grep("_NA$", seul[[lt]]$tcell_manual_anno)] <- NA
  
  cell_idx <- which(Cells(seul[[lt]]) %in% names(cd45_tcell_map))
  seul[[lt]]$manual_anno <- as.character(seul[[lt]]$manual_anno)
  seul[[lt]]$manual_anno[cell_idx] <- paste0("T_", cd45_tcell_map[Cells(seul[[lt]])[cell_idx]])
  seul[[lt]]$manual_anno[grep("_NA$", seul[[lt]]$manual_anno)] <- NA
  
  cell_idx <- which(Cells(seul[[lt]]) %in% names(cd45_treg_map))
  seul[[lt]]$manual_anno <- as.character(seul[[lt]]$manual_anno)
  seul[[lt]]$manual_anno[cell_idx] <- paste0("TReg_", cd45_treg_map[Cells(seul[[lt]])[cell_idx]])
  
  if(rm_removecells){
    Idents(seul[[lt]]) <- 'manual_anno'
    if(any(grepl("Remove", as.character(Idents(seul[[lt]])), ignore.case = T))){
      seul[[lt]] <- subset(seul[[lt]], ident=setdiff(unique(as.character(Idents(seul[[lt]]))), 'Remove'))
    }
  }
}

remove_cd3_cd19 <- TRUE
visualize=FALSE
if(remove_cd3_cd19){
  
  # Removes certain cells from clusters they're not suppose to belong to
  .rmCellsUsingMonocle <- function(seu, cellid, ident, mclusters, gene_expr=NULL, rmfrom=FALSE, keepat=FALSE){
    if(!rmfrom & !keepat){
      stop("Either rmfrom or keepat must be TRUE")
    }
    cells_idx <- Cells(seu)[which(seu$monocle3_clusters %in% mclusters)]
    if(mclusters[1] == -1) cells_idx <- Cells(seu)[1]
    seu_subset <-  subset(seu, cells=cells_idx)
    if(!is.null(gene_expr)){
      print(paste0("Removing cells expressing: ", gene_expr))
      if(keepat){
        cell_map <-  GetAssayData(seu, slot='counts')[gene_expr,] > 0
        cell_map[Cells(seu_subset)] <- FALSE
      } else {
        cell_map <- GetAssayData(seu_subset, slot='counts')[gene_expr,] > 0
      }
    } else {
      if(keepat){
        cell_map <- setNames(seu@meta.data[,ident] == cellid, Cells(seu))
        cell_map[Cells(seu_subset)] <- FALSE
      } else {
        cell_map <- setNames(seu_subset@meta.data[,ident] == cellid, Cells(seu_subset))
      }
    }
    cells_rm <- names(which(cell_map))
    rm(seu_subset);
    
    print(paste0(cellid, ": Removing: ", length(cells_rm), " cells..."))
    return(subset(seu, cells=setdiff(Cells(seu), cells_rm)))
  }
  
  # Iterates through each unique identity of a given cluster ID and plots the highlight plot for each
  .makeHighlightPlot <- function(seu, ident){
    Idents(seu) <- ident
    highlight_plot <- lapply(na.omit(unique(seu@meta.data[,ident])), function(clid){
      scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = adjustcolor(c('red', 'red'), alpha.f=0.4),
                                          reduction='umap', background_color = "lightgray", raster=T, pt.size=0.5) + 
        NoLegend() +
        ggtitle(clid)
    })
    return(highlight_plot)
  }
  
  set.seed(1234)
  seul[['LN']]  <- seul[['LN']] %>%
    .monocle3Cluster(.) 
  # pdf("~/xfer/ln2.pdf", width=16)
  # DimPlot(seul[['LN']], group.by='monocle3_clusters', reduction='umap', raster=T, label=T)
  # DimPlot(seul[['LN']], group.by='manual_anno', reduction='umap', raster=T, label=T)
  # dev.off()
  seu_ln <- seul[['LN']] %>%
    .rmCellsUsingMonocle(., cellid='B', ident='manual_anno', 
                                mclusters=c('20', '17', '16', '21'), rmfrom=TRUE) %>%
    .rmCellsUsingMonocle(., cellid='T_CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', ident='manual_anno', 
                         mclusters=c('4', '9', '15', '16', '22'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_001', ident='manual_anno', 
                         mclusters=c('11', '18'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='TReg_NLTlike_cycling', ident='manual_anno', 
                         mclusters=c('18', '11'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='TReg_Effector', ident='manual_anno', 
                         mclusters=c('11', '18', '15', '16'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD8_naive_memory', ident='manual_anno', 
                         mclusters=c('2', '13', '7'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_TReg', ident='manual_anno', 
                         mclusters=c('11', '18'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='TReg_Central', ident='manual_anno', 
                         mclusters=c('4', '15', '11', '18', '9'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD8_memory', ident='manual_anno', 
                         mclusters=c('2', '13', '7'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='TReg_NLTlike.Central', ident='manual_anno', 
                         mclusters=c('18', '11', '16', '21'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD4_naive', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD8_naive', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_Remove', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='22', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='18', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='16', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_Remove', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='30', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='26', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='14', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='25', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T)
  # Idents(seu_ln) <- 'manual_anno'
  # seu_lnsub <- subset(seu_ln, ident='T_CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_001')
  # pdf("~/xfer/y.pdf")
  # FeaturePlot(seu_lnsub, features=c('Foxp3', 'Il2ra', 'Cd4', 'Cd3e'), reduction='umap', raster=T)
  # dev.off()
  # Last minute fix for cluster IDs
  idx <- which((seu_ln$manual_anno == 'T_CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002') & 
                 (seu_ln$monocle3_clusters == '22'))
  seu_ln$manual_anno[idx] <- 'T_CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_003'
  idx <- which((seu_ln$manual_anno == 'TReg_Central') & 
                 (seu_ln$monocle3_clusters %in% c('15', '4', '9')))
  seu_ln$manual_anno[idx] <- 'TReg_Central_002'
  idx <- which((seu_ln$manual_anno == 'TReg_Central') & 
                 (seu_ln$monocle3_clusters %in% c('11')))
  seu_ln$manual_anno[idx] <- 'TReg_Central_001'
  idx <- which((seu_ln$manual_anno == 'TReg_Effector') & 
                 (seu_ln$monocle3_clusters %in% c('11', '18')))
  seu_ln$manual_anno[idx] <- 'TReg_Effector_001'
  idx <- which((seu_ln$manual_anno == 'TReg_Effector') & 
                 (seu_ln$monocle3_clusters %in% c('15', '4', '9', '10', '16', '22')))
  seu_ln$manual_anno[idx] <- 'TReg_Effector_002'
  seu_ln$manual_anno <- gsub("naive2$", "naive", seu_ln$manual_anno) %>%
    gsub("^T_", "", .)
  seu_ln$manual_anno <- gsub("^B$", "B_cells", seu_ln$manual_anno)
  seu_ln$manual_anno <- gsub("TReg_Central_001", "TReg_Effector_Central", seu_ln$manual_anno) %>%
    gsub("TReg_Effector_001", "TReg_Effector_Central", .)
  seu_ln$manual_anno <- gsub("TReg_Central_002", "TReg_Central", seu_ln$manual_anno) %>%
    gsub("TReg_Effector_002", "TReg_Effector", .)
  
  
  
  
  if(visualize){
    # pdf("~/xfer/LNfix_vln.pdf", width = 12, height = 13)
    # Idents(seu_ln) <- 'manual_anno'
    # VlnPlot(seu_ln, features=c('Cd8a', 'Cd4', 'Cd3e', 'Itgam', 'Itgax',
    #                            'Adgre1', 'Foxp3', 'Cd19', 'Cd14'))
    # dev.off()
    pdf("~/xfer/lnfix.pdf", width = 7,height = 7)
    DimPlot(seu_ln, group.by='manual_anno', raster=T, reduction='umap')
    print(.makeHighlightPlot(seu_ln, ident='manual_anno'))
    dev.off()
  }
  
  
  # Removes CD3 positive cells from non-Tcell clusters
  seul[['Tumor']] <- seul[['Tumor']] %>%
    .monocle3Cluster(.)
  # pdf("~/xfer/tumorfix.pdf", width = 7,height = 7)
  # DimPlot(seul[['Tumor']], group.by='monocle3_clusters', reduction='umap', raster=T, label=T)
  # FeaturePlot(seul$Tumor, features=c('Cd3e', 'Cd4', 'Foxp3', 'Cd8a'), raster=T, reduction='umap')
  # print(.makeHighlightPlot(seul[['Tumor']], ident='manual_anno'))
  # dev.off()
  seu_tumor <- seul[['Tumor']] %>%
    .rmCellsUsingMonocle(., gene_expr = 'Cd3e', cellid="NA", ident='manual_anno', 
                       mclusters=c('2', '1', '6', '18', '12'), rmfrom=T) %>%
    .rmCellsUsingMonocle(., cellid='NK', ident='manual_anno', 
                         mclusters=c('2'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD8_memory', ident='manual_anno', 
                         mclusters=as.character(c(3, 7, 10, 4, 9, 11)), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='CD8', ident='manual_anno', 
                         mclusters=as.character(c(3, 7, 10, 4, 9, 11, 8)), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD8_effector_cycling', ident='manual_anno', 
                         mclusters=as.character(c(3)), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='DC', ident='manual_anno', 
                         mclusters=c('6'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='CD3_DN_Naive', ident='manual_anno', 
                         mclusters=c('16'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD4_CCR7pos_TCF7pos', ident='manual_anno', 
                         mclusters=as.character(c(3, 7, 10, 4, 9, 11)), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_001', ident='manual_anno', 
                         mclusters=c('-1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD3_DN_activated', ident='manual_anno', 
                         mclusters=c('11', '9'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='Monocyte_derived_macrophage', ident='manual_anno', 
                         mclusters=c('1'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='T_CD8_memory2', ident='manual_anno', 
                         mclusters=c('8'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='TReg_NLTlike', ident='manual_anno', 
                         mclusters=c('5', '14'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='Monocyte', ident='manual_anno', 
                         mclusters=c('19'), keepat=T) %>%
    .rmCellsUsingMonocle(., cellid='ILC3', ident='manual_anno', 
                         mclusters=c('13'), keepat=T)
  
  idx <- which(seu_tumor$manual_anno == 'T_CD3_DN')
  seu_tumor$manual_anno[idx] <- 'T_CD3_DN_activated'
  seu_tumor$manual_anno <- gsub("naive2$", "naive", seu_tumor$manual_anno) %>%
    gsub("^T_", "", .)
  seu_tumor$manual_anno <- gsub("^B$", "B_cells", seu_tumor$manual_anno)
  
  if(visualize){
    pdf("~/xfer/tumorfix.pdf", width = 14,height = 7)
    DimPlot(seu_tumor, group.by='manual_anno', raster=T, reduction='umap')
    # print(.makeHighlightPlot(seu_tumor, ident='manual_anno'))
    dev.off()
  }
  seul[['Tumor']] <- seu_tumor
  seul[['LN']] <- seu_ln
  
}

#--- viii.a) CD45 - total cell plot - April 2024 ----
if(!exists("seul"))  seul <- readRDS(file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
outdir <- file.path(PDIR, "results")
seu <- seul$Tumor 
grouping <- 'cd45_all_cells'
overwrite=F

## Remove and relabel cells
seu@meta.data$manual_anno[grep("macrophage", seu@meta.data$manual_anno, ignore.case=T)] <- 'Macrophage'
seu <- subset(seu, cells=Cells(seu)[grep("CD3_DN_Naive", seu$manual_anno, invert = T)])
seu$newid <- gsub("(CD8|TReg|CD4|NK)_.*", "\\1\\2", seu$manual_anno)
seu$newid[seu$manual_anno == 'CD8_memory2'] <- 'Gamma_Delta_Tcells'
seu$newidcycling <- gsub("(CD8|TReg|CD4|NK)_.*((_cycling)?)?", "\\1\\2", seu$manual_anno)
seu$newidcycling[seu$manual_anno == 'CD8_memory2'] <- 'Gamma_Delta_Tcells'
# table(seu$manual_anno, gsub("(CD8|TReg|CD4|NK)_.*", "\\1", seu$manual_anno))
# table(seu$manual_anno, gsub("(CD8|TReg|CD4|NK)_.*((_cycling)?)?", "\\1\\2", seu$manual_anno))
# gsub("(CD8|TReg|CD4|NK)_.*((_cycling)?)?", "\\1\\2", c('CD8_effector', 'CD8_effector_cycling', 'TReg_NLTlike'))
seul$Tumor <- seu
pdf("~/xfer/y.pdf", width = 13)
publicationDimPlot(seul$Tumor, grp='manual_anno', reduction='umap')
publicationDimPlot(seul$Tumor, grp='newid', reduction='umap')
publicationDimPlot(seul$Tumor, grp='newidcycling', reduction='umap')
dev.off()
seul$LN$newidcycling <- seul$LN$newid <- seul$LN$manual_anno


## X
# 
# id <- 'Tumor'
# seux <- seul[[id]]@meta.data %>% 
#   dplyr::select(c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, 
#                   manual_anno)) %>%
#   rename_with(., ~gsub("manual_anno", "anno_allCells", .))
# seux <- seul[[id]]@reductions$umap@cell.embeddings  %>%
#   as.data.frame %>% 
#   rename_with(., ~c('umap_allCells_1', 'umap_allCells_2')) %>%
#   cbind(seux, .) %>%
#   tibble::rownames_to_column(., "barcode")
# 
# seuy <- seu@meta.data %>% 
#   dplyr::select(manual_anno) %>%
#   cbind(.,  seu@reductions$umap@cell.embeddings %>% 
#           as.data.frame %>% 
#           rename_with(., ~c('umap_Cd8_1', 'umap_Cd8_2')))  %>%
#   tibble::rownames_to_column(., "barcode")
# 
# seuxy <- left_join(seux, seuy, by='barcode')
# dir='/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/geo_submission/scrna/processed'
# write.table(seuxy, file=file.path(dir, paste0("CD45_", id, "_meta.csv")),
#             sep=",", col.names = T, row.names = F, quote = F)
# cnts <- GetAssayData(seul[[id]], assay='RNA', slot='counts')
# as.data.frame(cnts) %>% tibble::rownames_to_column(., "feature") %>%
#   write.table(., file=file.path(dir, paste0("CD45_", id, "_cnts.csv")),
#               sep=",", col.names = T, row.names = F, quote = F)
## X

## FindAllMarkers on everything and run DotPlot
id = 'newidcycling'
allmarkers <- lapply('Tumor', function(seuid){
  seu <- seul[[seuid]]
  Idents(seu) <- id
  seuid_markers <- FindAllMarkers(object = seu) %>%
    scCustomize::Add_Pct_Diff() 
  return(seuid_markers)
})
names(allmarkers) <-'Tumor'

for(seuid in names(allmarkers)){
  amarkers <- allmarkers[[seuid]] %>%
    dplyr::filter(pct_diff > 0.6)
  top_markers <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 20, 
                                     named_vector = FALSE, make_unique = TRUE)
  top_markers_tbl <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 50, 
                                     named_vector = FALSE, make_unique = TRUE)
  seu <- seul[[seuid]]
  Idents(seu) <- id
  genes <- read.csv(file.path(PDIR, "ref", "dotplots", paste0("CD45_", seuid, ".csv")),
                    header = T, sep=",")
  
  write.table(amarkers %>% dplyr::filter(gene %in% top_markers_tbl),
              file=file.path("~/xfer", paste0("CD45_", seuid, ".csv")), 
              sep=",", col.names = T, row.names = F, quote = F)
  
  pdf(file.path("~/xfer", paste0("CD45_", seuid, ".pdf")), height = 20)
  scCustomize::DotPlot_scCustom(seurat_object = seu, 
                                features = unique(genes$gene),
                                flip_axes = T,
                                x_lab_rotate = TRUE)
  scCustomize::Clustered_DotPlot(seurat_object = seu, 
                                 features = unique(genes$gene),
                                 k=length(unique(Idents(seu))),
                                 x_lab_rotate = TRUE, print_exp_quantiles = T,
                                 exp_color_min = -1, exp_color_max = 3)
  
  scCustomize::DotPlot_scCustom(seurat_object = seu, 
                                features = top_markers,
                                flip_axes = T,
                                x_lab_rotate = TRUE)
  scCustomize::Clustered_DotPlot(seurat_object = seu, 
                                 features = top_markers,
                                 k=length(unique(Idents(seu))),
                                 x_lab_rotate = TRUE, print_exp_quantiles = T,
                                 exp_color_min = -1, exp_color_max = 3)
  dev.off()
  
  cat(paste0("xfer CD45_", seuid, ".csv\n"))
  cat(paste0("xfer CD45_", seuid, ".pdf\n"))
  
}


### DEG UMAP highlight plot
method <- 'wilcox'
sigdig <- 3
deg_umapgg <- lapply(names(seul)[2], function(seuid){
  seu <- seul[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$manual_anno)
  
  # Create testing groups
  uids <- unique(seu$group)
  x <- split(uids, f=gsub("^.*\\.", "", uids))
  xt_all <- lapply(x, function(i){
    data.frame("V1"=grep("KO_7d", i, value=T),
               "V2"=grep("WT_7d", i, value=T),
               "baselvl"=grep("KO_7d", i, value=T),
               'id1'='KO_7d',
               'id2'='WT_7d')
  }) %>% do.call(rbind, .)
 
  outf = file.path(outdir, grouping, paste0(seuid, ".wilcox_deg.rds"))
  if(!file.exists(outf) | overwrite){
    markers_all <- apply(xt_all, 1, function(ids){
      print(paste0("Testing ", ids[4], " - ", ids[5]))
      ids_test <- factor(c(ids[1], ids[2]))
      ids_test <- relevel(ids_test, ids[3])
      ord <- order(ids_test)
      
      Idents(seu) <- 'group'
      tryCatch({
        FindMarkers(seu, 
                    ident.1=ids[c(1,2)[ord[1]]],
                    ident.2=ids[c(1,2)[ord[2]]], 
                    logfc.threshold=0, 
                    test.use=method) %>%
          tibble::rownames_to_column(., "gene") %>%
          mutate(ens=gm$SYMBOL$ENSEMBL[gene],
                 biotype=gm$ENSEMBL$gene_biotype[ens],
                 id1=ids[c(4,5)[ord[1]]],
                 id2=ids[c(4,5)[ord[2]]],
                 fullid1=ids[c(1,2)[ord[1]]],
                 fullid2=ids[c(1,2)[ord[2]]],
                 p_val = round(p_val, sigdig),
                 avg_log2FC=round(avg_log2FC, sigdig),
                 p_val_adj = ifelse(p_val_adj < 10^(-1*sigdig), 
                                    p_val_adj, 
                                    round(p_val_adj, sigdig))) %>%
          relocate(., c(gene, ens, biotype, id1, id2)) %>%
          rename_with(., ~gsub("p_val_adj", paste0(method, ".p_val_adj"), .))
      }, error=function(e){NULL})
      
    })
    saveRDS(markers_all, file=outf)
  } else {
    markers_all <- readRDS(outf)
  }
  
  

  markers_sub <- lapply(markers_all, function(i){
    i %>% 
      dplyr::filter(wilcox.p_val_adj <= 0.05,
                    abs(avg_log2FC) >= 0.5,
                    biotype == 'protein_coding',
                    !grepl("Rik$", gene),
                    !grepl("^Gm", gene))
  })
  names(markers_sub) <- gsub("\\.[0-9]*$", "", names(markers_sub))
  tbl <- table(seu@meta.data$manual_anno, seu@meta.data$newid) %>%
    as.data.frame %>%  dplyr::filter(Freq > 0)
  idmap.v <- with(tbl, setNames(as.character(Var2), as.character(Var1)))
  ctcnts <- sapply(split(markers_sub, idmap.v[names(markers_sub)]), function(i) sum(sapply(i, nrow)))
  
  
  
  maxval <- if(seuid == 'LN') 150 else 101
  cols <- setNames(viridis_plasma_dark_high[ceiling(ctcnts * 1/(maxval/250))+1],
                   names(ctcnts))
  Idents(seu) <- 'newid'
  pdp <- publicationDimPlot(seu, grp='newid', colors_use=as.character(cols[levels(Idents(seu))]), 
                            reduction='umap', return='list')
  
  
  ggleg <- ggplot(data=data.frame(id=seq_len(maxval)), 
                  aes(x=id, y=id, color=id)) +
    geom_point() +
    scale_color_gradientn(colors = scCustomize::viridis_plasma_dark_high) +
    labs(color='# of DEGs')
  leg <- ggpubr::get_legend(ggleg)
  
  pdp$plot <- pdp$plot + ggtitle(paste0(seuid))
  plot_figure <- pdp$plot + pdp$axis_plot + ggpubr::as_ggplot(leg) +
    patchwork::plot_layout(design = c(
      patchwork::area(t = 1, l = 2, b = 11, r = 11),
      patchwork::area(t = 10, l = 1, b = 12, r = 2),
      patchwork::area(t = 1, l = 11, b = 4, r = 12))) & 
    theme(aspect.ratio=pdp$aspect.ratio)
  return(plot_figure)
  pdf("~/xfer/CD45_Tumor.degUMAP.pdf"); plot_figure; dev.off()
})









pdf("~/xfer/CD45_total_cells.pdf", width = 15,height = 15)
ggdim_seul <- lapply(names(seul), function(seuid){
  ggdim0 <- lapply(c('manual_anno'), function(grp){
    # Relabel the clusters from Treg to 1_Treg to simplify the overlay
    seu_i = seul[[seuid]]
    
    grp_ids <- na.omit(unique(seu_i@meta.data[,grp]))
    grp_map <- setNames(seq_along(grp_ids), grp_ids)
    legend_map <- setNames(paste(grp_map, names(grp_map), sep="_"), grp_map)
    
    seu_i@meta.data[,grp] <- factor(grp_map[seu_i@meta.data[,grp]])
    labels <- levels(seu_i@meta.data[,grp])
    
    DimPlot(seu_i, group.by=grp, raster=T, label=T, reduction='umap') + 
      scale_color_discrete(labels=legend_map[labels]) +
      guides(color = guide_legend(ncol = 1))
  }) %>% 
    cowplot::plot_grid(plotlist=., nrow=1)
  ggseuid <- ggdraw() + draw_label(label = seuid)
  cowplot::plot_grid(ggseuid, ggdim0, nrow=1, rel_widths = c(0.2, 1))
}) %>%
  cowplot::plot_grid(plotlist=., nrow=2, align = "v")
ggdim_seul
dev.off()

grp <- 'manual_anno'
pdf("~/xfer/cd45_publication.pdf", width=7, height=5)
grps <- sapply(seul, function(seu_i) unique(seu_i@meta.data[,grp])) %>%
  unlist %>% 
  unique
colors_use <- scCustomize::scCustomize_Palette(num_groups = length(grps), 
                                               ggplot_default_colors = FALSE, 
                                               color_seed = 231) %>%
  setNames(., grps)
lapply(seul, function(seu_i){
  colors_use_i <- as.character(colors_use[unique(seu_i@meta.data[,grp])])
  publicationDimPlot(seu_i, grp=grp, simplify_labels=TRUE, colors_use=colors_use_i, 
                     pt.size=0.2, aspect.ratio=1,
                     legend.text = element_text(size=10))
})
dev.off()



if(overwrite) saveRDS(seul, file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
seul <- readRDS(file=file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))


write_loupe <- F
if(write_loupe){
  seul <- readRDS(file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
  for(id in names(seul)){
    seu <- seul[[id]]
    seu@meta.data <- seu@meta.data[,c('orig.ident', 'seurat_clusters', 'manual_anno', 'monocle3_clusters')]
    seu[['RNA2']] <- CreateAssayObject(counts=seu@assays$RNA$counts)
    DefaultAssay(seu) <- 'RNA2'
    
    dir.create(file.path(PDIR, "results", "cloupe"), showWarnings = F)
    loupeR::create_loupe_from_seurat(seu, output_dir=file.path(PDIR, "results", "cloupe"),
                                     output_name=paste0("CD45_", id), force=T, tmpdir=PDIR)
    file.copy(file.path(PDIR, "results", "cloupe", paste0("CD45_", id, ".cloupe")), to = "~/xfer", overwrite = T)
    cat(paste0("xfer ", paste0("CD45_", id, ".cloupe\n")))
    
  }
  
}

if(do_visualize){
  pdf("~/xfer/xTfixed.pdf", width = 7,height = 7)
  .makeHighlightPlot(seul$Tumor, ident='manual_anno')
  dev.off()
  pdf("~/xfer/LN.pdf", width = 7,height = 7)
  .makeHighlightPlot(seul$LN, ident='manual_anno')
  dev.off()
  
  pdf("~/xfer/y.pdf", width = 7,height = 7)
  DimPlot(seul$LN, group.by='seurat_clusters', reduction='umap', raster=T, label=T)
  dev.off()
  
  pdf("~/xfer/z.pdf", width = 14,height = 14)
  seu <- seul$LN
  seu <- CellCycleScoring(seu, s.features = str_to_title(cc.genes$s.genes), 
                                                      g2m.features = str_to_title(cc.genes$g2m.genes),
                                                      assay='SCT', set.ident = TRUE)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  FeaturePlot(seu, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mt',  'CC.Difference',
                              'S.Score', 'G2M.Score'), reduction='umap', raster=T)
  dev.off()
  
  
  pdf("~/xfer/CD45_total_cells_with_tregs_tcells.pdf", width = 45,height = 15)
  # seul <- cleanLnTumorTreg(seul)
  ggdim_seul <- lapply(names(seul), function(seuid){
    ggdim0 <- lapply(c('immgen.simple', 'treg_manual_anno', 'tcell_manual_anno', 'manual_anno'), function(grp){
      # Relabel the clusters from Treg to 1_Treg to simplify the overlay
      seu_i = seul[[seuid]]
      cell_idx <- (seu_i$manual_anno %in% c(unique(grep("Remove", seu_i$manual_anno, value=T)), 
                                            c('CD4_CCR7hi', 'Patrolling_monocyte')))
      seu_i <- subset(seu_i, cells=Cells(seu_i)[which(!cell_idx)])
      
      grp_ids <- na.omit(unique(seu_i@meta.data[,grp]))
      grp_map <- setNames(seq_along(grp_ids), grp_ids)
      legend_map <- setNames(paste(grp_map, names(grp_map), sep="_"), grp_map)
      
      seu_i@meta.data[,grp] <- factor(grp_map[seu_i@meta.data[,grp]])
      labels <- levels(seu_i@meta.data[,grp])
      
      DimPlot(seu_i, group.by=grp, raster=T, label=T, reduction='umap') + 
        scale_color_discrete(labels=legend_map[labels])
    }) %>% 
      cowplot::plot_grid(plotlist=., nrow=1)
    ggseuid <- ggdraw() + draw_label(label = seuid)
    cowplot::plot_grid(ggseuid, ggdim0, nrow=1, rel_widths = c(0.02, 1))
  }) %>%
    cowplot::plot_grid(plotlist=., nrow=2)
  ggdim_seul
  dev.off()
}

#--- Subset CD8s and DC's in CD45-Tumor ----
seul <- readRDS(file=file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
seu <- seul$Tumor

Idents(seu) <- 'manual_anno'
seuT_cd8 <- subset(seu, ident=grep('^cd8', unique(Idents(seu)), value=T, ignore.case = T))
seuT_dc <- subset(seu, ident=grep('dc', unique(Idents(seu)), value=T, ignore.case = T))
seuT_l <- list("cd8"=seuT_cd8, "dc"=seuT_dc)
rm(seuT_cd8, seuT_dc); gc()

seuT_l <- lapply(seuT_l, function(seuT_i){
  DefaultAssay(seuT_i) <- 'RNA'
  seuT_i <- seuT_i %>%
    FindVariableFeatures(., selection.method = "vst",
                         nfeatures = 3000, verbose = FALSE) %>% 
    ScaleData(.) %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(., resolution=0.6)
  
  return(seuT_i)
  
})

## Section to remove CD4's
{
  seu_i <- seuT_l$cd8
  
  cd4pos_cells <- Cells(seu_i)[which(GetAssayData(seu_i, slot='counts')['Cd4',] > 0)]
  clus_cells <- Cells(seu_i)[seu_i$seurat_clusters %in% c('5', '10')]
  rm_cells <- intersect(cd4pos_cells, clus_cells)
  seuT_l$cd8 <- subset(seu_i, cells=setdiff(Cells(seu_i), rm_cells)) %>%
    FindVariableFeatures(., selection.method = "vst",
                         nfeatures = 3000, verbose = FALSE) %>% 
    ScaleData(.) %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(., resolution=0.3)
}

## Section to remove the doublet cluster
{
  seu_i <- seuT_l$cd8
  seu_i$old_clusters <- seu_i$seurat_clusters
  Idents(seu_i) <- 'seurat_clusters'

  seu_i2 <- subset(seu_i, ident=setdiff(unique(Idents(seu_i)), '5')) %>%
    FindVariableFeatures(., selection.method = "vst",
                         nfeatures = 3000, verbose = FALSE) %>% 
    ScaleData(.) %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(., resolution=0.3)
  # pdf("~/xfer/cd8_minusCd4.minusDubClus.pdf", width = 14, height = 16)
  # p1 <- DimPlot(seu_i, group.by='seurat_clusters', raster=T, label=T)
  # p2 <- DimPlot(seu_i2, group.by='old_clusters', raster=T, label=T)
  # p3 <- DimPlot(seu_i2, group.by='seurat_clusters', raster=T, label=T)
  # p4 <- DimPlot(seu_i2, group.by='manual_anno', raster=T, label=T)
  # cowplot::plot_grid(p1, p2, p3, p4, ncol=2)
  # feats <- c('Pdcd1', 'Havcr2', 'Lag3', 'Cd69', 'Tcf7', 'Gzmb', 'Ly108', 'Slamf6')
  # VlnPlot(seu_i2, group.by='manual_anno', features=feats, raster=T)
  # feats <- c('Cd4','Cd3e','Cd3g','Cd3d','Cd8a','Cd8b1','Klrd1','Klrb1c','Ncam1','Itga2','Il2rb','Klra7')
  # FeaturePlot(seu_i2, features=feats, raster=T, pt.size=3)
  # dev.off()
  seuT_l$cd8 <- seu_i2
}

## Section to remove the Ncr1, Klrb1c expressing cells
{
  seu_i <- seuT_l$cd8
  markerpos_cells <- Cells(seu_i)[which((GetAssayData(seu_i, slot='counts')['Ncr1',] > 0) |
                                          (GetAssayData(seu_i, slot='counts')['Klrb1c',] > 0))]
  clus_cells <- Cells(seu_i)[seu_i$seurat_clusters %in% c('0')]
  rm_cells <- intersect(markerpos_cells, clus_cells)
  seuT_l$cd8 <- subset(seu_i, cells=setdiff(Cells(seu_i), rm_cells)) %>%
    FindVariableFeatures(., selection.method = "vst",
                         nfeatures = 3000, verbose = FALSE) %>% 
    ScaleData(.) %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(., resolution=0.3)
  # pdf("~/xfer/x.pdf", width = 12, height = 9)
  # p1 <- publicationDimPlot(seuT_l$cd8, grp='seurat_clusters', pt.size=2)
  # p2 <- publicationDimPlot(seuT_l$cd8, grp='old_clusters', pt.size=2)
  # p4 <- publicationDimPlot(seuT_l$cd8, grp='manual_anno', pt.size=2)
  # cowplot::plot_grid(p1, p2,p4, ncol=2, align='hv')
  # feats <- c('Klrb1c', 'Ncr1')
  # FeaturePlot(seuT_l$cd8, features=feats, raster=T, pt.size=3)
  # dev.off()
  
}

seu <- DietSeurat(seuT_l$cd8, assays='RNA', dimreducs=c('umap', 'mnn'))
seu@meta.data <- seu@meta.data[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'seurat_clusters', 'manual_anno')]
saveRDS(seu, file=file.path(datadir, "seurat_obj", "CD45_Tumor_CD8.seuratobj.rds"))

write.cloupe=F
if(write.cloupe){
  seu = readRDS(file=file.path(datadir, "seurat_obj", "CD45_Tumor_CD8.seuratobj.rds"))
  seu@meta.data <- seu@meta.data[,c('orig.ident', 'seurat_clusters', 'manual_anno')]
  seu[['RNA2']] <- CreateAssayObject(counts=seu@assays$RNA$counts)
  DefaultAssay(seu) <- 'RNA2'
  
  dir.create(file.path(PDIR, "results", "cloupe"), showWarnings = F)
  loupeR::create_loupe_from_seurat(seu, output_dir=file.path(PDIR, "results", "cloupe"),
                                   output_name="cd45_cd8_tumor", force=T)
  file.copy(file.path(PDIR, "results", "cloupe", "tregs_all_samples.selgenes.cloupe"), to = "~/xfer", overwrite = T)
  cat(paste0("xfer tregs_all_samples.selgenes.cloupe\n"))
  
  
}

# FindAllMarkers deg plot
{
  seu = readRDS(file=file.path(datadir, "seurat_obj", "CD45_Tumor_CD8.seuratobj.rds"))
  Idents(seu) <- 'manual_anno'
  seu <- subset(seu, ident=grep("memory2", unique(Idents(seu)), value=T, invert = T))
  seus <- list('Tumor'=seu)
  grp <- 'manual_anno'
  
  allmarkers <- lapply(names(seus), function(seuid){
    seu <- seus[[seuid]]
    Idents(seu) <- grp
    seuid_markers <- FindAllMarkers(object = seu) %>%
      scCustomize::Add_Pct_Diff() 
    return(seuid_markers)
  })
  names(allmarkers) <- names(seus)
  
  for(seuid in names(seus)){
    amarkers <- allmarkers[[seuid]] %>%
      dplyr::filter(pct_diff > 0.2)
    top_markers <- amarkers %>%
      scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 20, 
                                       named_vector = FALSE, make_unique = TRUE)
    top_markers_tbl <- amarkers %>%
      scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 50, 
                                       named_vector = FALSE, make_unique = TRUE)
    seu <- seus[[seuid]]
    Idents(seu) <- grp
    seu <- subset(seu, idents=grep(".+", unique(Idents(seu)), value=T))
    seu <- seus[[seuid]]
    Idents(seu) <- grp
    genes <- read.csv(file.path(PDIR, "ref", "dotplots", paste0("CD45_", seuid, "_CD8.csv")),
                      header = T, sep=",")
    
    write.table(amarkers %>% dplyr::filter(gene %in% top_markers_tbl),
                file=file.path("~/xfer", paste0("CD45.CD8_", seuid, ".csv")), 
                sep=",", col.names = T, row.names = F, quote = F)
    
    pdf(file.path("~/xfer", paste0("CD45.CD8_", seuid, ".pdf")), height = 12)
    scCustomize::DotPlot_scCustom(seurat_object = seu, 
                                  features = unique(genes$gene),
                                  flip_axes = T,
                                  x_lab_rotate = TRUE)
    scCustomize::Clustered_DotPlot(seurat_object = seu, 
                                   features = unique(genes$gene),
                                   k=length(unique(Idents(seu))),
                                   x_lab_rotate = TRUE, print_exp_quantiles = T,
                                   exp_color_min = -1, exp_color_max = 3)
    
    scCustomize::DotPlot_scCustom(seurat_object = seu, 
                                  features = top_markers,
                                  flip_axes = T,
                                  x_lab_rotate = TRUE)
    scCustomize::Clustered_DotPlot(seurat_object = seu, 
                                   features = top_markers,
                                   k=length(unique(Idents(seu))),
                                   x_lab_rotate = TRUE, print_exp_quantiles = T,
                                   exp_color_min = -1, exp_color_max = 3)
    dev.off()
    
    pdf(file.path("~/xfer", paste0("CD45.CD8_", seuid, "_umap.pdf")), width = 7, height = 7)
    publicationDimPlot(seu, grp=grp, reduction='umap', 
                       simplify_labels=if(grp=='manual_anno') T else F) %>% plot
    dev.off()
    
    cat(paste0("xfer CD45.CD8_", seuid, ".csv\n"))
    cat(paste0("xfer CD45.CD8_", seuid, "_umap.pdf\n"))
    cat(paste0("xfer CD45.CD8_", seuid, ".pdf\n"))
    
    
    seu@meta.data %>% 
      dplyr::select(c(orig.ident, nCount_RNA, nFeature_RNA, manual_anno)) %>%
      cbind(., seu@reductions$umap@cell.embeddings) %>%
      write.table(., file="~/xfer/CD45_Tumor_cd8.csv", sep=",", 
                  quote=F, col.names = T, row.names = T)
  }
}

# ViolinPlot of AUCell score for CD8 signature
{
  expr <- GetAssayData(seu, slot='counts')
  expr <- expr[rowSums(expr)>=50, ]
  f <- list.files(file.path(PDIR, "ref", "cd8_signature"))
  gs <- lapply(f, function(fi) read.table(file.path(PDIR, "ref", "cd8_signature", fi),
                                    sep=",", header = F)$V1) %>% 
    setNames(., gsub(".csv$", "", f))
  
  score <- AUCell::AUCell_run(expr, gs)
  score_df <- cbind(t(assay(score)),
                    seu@meta.data[,c('orig.ident', 'manual_anno')]) %>% 
    as.data.frame %>%
    tidyr::pivot_longer(., cols=!c(manual_anno, orig.ident))
  
  pdf("~/xfer/cd8_signature_scores.pdf")
  X <- split(score_df, f=score_df$name)
  lapply(X, function(x){
    model <- aov(value~manual_anno, data=as.data.frame(x))
    TukeyHSD(model, conf.level=.95)
  })
  ggplot(score_df, aes(x=manual_anno, y=value, fill=manual_anno))+
    facet_grid(name~., space='free', scales='free') +
    geom_violin() + 
    cowplot::theme_cowplot() +
    ylim(0, 0.125) +
    ylab("AUCell score") +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1),
          axis.title.x = element_blank())
  dev.off()
}

seu <- ScaleData(seu, features=Features(seu))
pdf("~/xfer/x.pdf", height = 6, width = 8)
features <- c('Mtor', 'Pten', 'Il33')
scCustomize::Stacked_VlnPlot(
  seu,
  features=features,
  split.by='orig.ident',
  group.by='manual_anno',
  x_lab_rotate = TRUE
)
lapply(features, function(f){
  scCustomize::VlnPlot_scCustom(
    seu,
    features=f,
    split.by='orig.ident',
    group.by='manual_anno',
    raster=T
  )
})
dev.off()
## Output mappings for loupe browser
{
  map_obj <- mapSmergeToCRaggr(seuT_l$cd8, 
                               file.path(PDIR, "data", 'cellranger_aggr', 'cd45', 'outs', 'count', 'filtered_feature_bc_matrix'),
                               return_metadata=T, return_reduction=T)
  
  metadata <- map_obj$metadata %>%
    select(c(Barcode, manual_anno, seurat_clusters, projectils.simple)) 
  umap <- map_obj$umap
  write.table(metadata, file="~/xfer/metadata.csv", sep = ",", quote = F, row.names = F, col.names = T)
  write.table(umap, file="~/xfer/umap.csv", sep = ",", quote = F, row.names = F, col.names = T)
}

## Output mappings for scvelo
{
  
  ### SeuratV5 SeuratDisk WORKAROUND ###
  # assigning the previous version of the `[[` function for the Assay class to the SeuratDisk package environment
  "[[.Assay" <- function(x, i, ..., drop = FALSE) {
    if (missing(x = i)) {
      i <- colnames(x = slot(object = x, name = 'meta.features'))
    }
    data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
    if (drop) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
    }
    return(data.return)
  }
  environment(`[[.Assay`) <- asNamespace("SeuratObject")
  rlang::env_unlock(asNamespace("SeuratDisk"))
  assign("[[.Assay", `[[.Assay`, asNamespace("SeuratDisk"))
  lockEnvironment(asNamespace("SeuratDisk"), bindings = TRUE)
  rm(`[[.Assay`)
  ### WORKAROUND ###
  
  ## Read in the velocyto loom files
  velocyto_dir <- file.path(datadir, "velocyto")
  sample_ids <- grep("^CD45_Tumor", list.files(velocyto_dir), value=T)
  seu_velo_raw <- lapply(setNames(sample_ids,sample_ids), function(sid){
    print(sid)
    f <- list.files(file.path(velocyto_dir, sid))
    loom.data <- ReadVelocity(file = file.path(velocyto_dir, sid, f)) %>%
      as.Seurat(.)
    loom.data$orig.ident <- sid
    loom.data <- RenameCells(loom.data, 
                             new.names=gsub("^M[ca]Gaha_.*__(.*)\\:([ACGT]*).*", "\\1_\\2", Cells(loom.data)) %>%
                               gsub("^T", "Tumor", .))
    return(loom.data)
  })
  
  
  
  ## Preprocess the seurat velocyto file
  # Extract relevant seu_velo data
  # seu_velo <- seu_velo_raw[grep("^CD45_Tumor", names(seu_velo_raw), value=T)]
  if(class(seu_velo_raw) == 'list'){
    seu_velo <- merge(seu_velo_raw[[1]], seu_velo_raw[-1], 
                      project = 'CD45_TDLN_LN.Tumor')
  }
  rm(seu_velo_raw); gc()
  
  

  DefaultAssay(object = seu_velo) <- "spliced"
  gsub("_[ACGT]*$", "", Cells(seu_velo)) %>% table
  gsub("_[ACGT]*$", "", Cells(seuT_l$cd8)) %>% table
  seu_velo <- subset(seu_velo, cells=gsub("_7d", "_B", Cells(seuT_l$cd8)) %>%
                       gsub("_Un", "_U", .))
  
  # Normalize and cluster
  set.seed(seed)
  seu_velo[["RNA"]] <- seu_velo[["spliced"]]
  seu_velo <- seu_velo %>%
    SCTransform(.) %>%
    FindVariableFeatures(.) %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(.) 

  DefaultAssay(seu_velo) <- "spliced"
  seu_velo <- FindVariableFeatures(seu_velo)
  umap <- Embeddings(seuT_l$cd8, reduction='umap')
  rownames(umap) <- gsub("_7d", "_B", rownames(umap)) %>% gsub("_Un", "_U", .)
  seu_velo[['umap_orig']] <- CreateDimReducObject(embeddings = umap[Cells(seu_velo),],
                                                        key = 'UMAP_', assay = 'RNA')
  anno <- seuT_l$cd8$manual_anno
  names(anno) <- gsub("_7d", "_B", names(anno)) %>% gsub("_Un", "_U", .)
  seu_velo$functional_cluster <- anno[Cells(seu_velo)]
  
  print("saving...")
  DefaultAssay(seu_velo) <- "RNA"
  SeuratDisk::SaveH5Seurat(seu_velo, 
                           filename = file.path(outdir, "scVelo", 
                                                paste0("seu_velo.CD45_Tumor_cd8.h5Seurat")),
                           overwrite=T)
  cwd <- getwd()
  setwd(file.path(outdir, "scVelo"))
  SeuratDisk::Convert(source = paste0("seu_velo.CD45_Tumor_cd8.h5Seurat"), 
                      dest = "h5ad", overwrite=T)
  setwd(cwd)
  
}


## Run DEG
{
  comparisons <- list(
    "WTtx_vs_WTun" = c("CD45_Tumor_WT_7d", "CD45_Tumor_WT_Un"),
    "WTtx_vs_KOtx" = c("CD45_Tumor_WT_7d", "CD45_Tumor_KO_7d"),
    "KOtx_vs_KOun" = c("CD45_Tumor_KO_7d", "CD45_Tumor_KO_Un"),
    "WTun_vs_KOun" = c("CD45_Tumor_WT_Un", "CD45_Tumor_KO_Un"))
  Idents(seuT_l$cd8) <- 'manual_anno'
  uniqclusters <- unique(as.character(Idents(seuT_l$cd8)))
  
  for(clustid in uniqclusters){
    seusub <- subset(seuT_l$cd8, ident=clustid)
    Idents(seusub) <- 'orig.ident'
    
    for(compid in names(comparisons)){
      print(compid)
      idents <- comparisons[[compid]]
      deg_markers <- FindMarkers(seusub, 
                                 assay.use='RNA', test.use='wilcox',
                                 ident.1=idents[1], ident.2=idents[2],
                                 verbose = FALSE,
                                 logfc.threshold = 0,
                                 recorrect_umi =  FALSE) %>% 
        tibble::rownames_to_column(., "symbol") %>% 
        mutate(biotype = sym2biotype_ids[symbol],
               ensemble = gm$SYMBOL$ENSEMBL[symbol],
               ident.1 = idents[1],
               ident.2 = idents[2], 
               tissue=compid) %>% 
        relocate(., c('ensemble', 'biotype'))
      write.table(deg_markers, 
                  file=paste0("~/xfer/cd8_minusCd4DubsNcriKlrb1.", clustid, ".", compid, ".csv"), 
                  sep=",", col.names = T, row.names =F, quote = F)
      
    }
  }
}

{
  Idents(seuT_l$dc) <- 'manual_anno'
  markers <- FindMarkers(seuT_l$dc, 
                         assay.use='RNA', test.use='wilcox',
                         ident.1='DC', ident.2='DC_Cd4pos_Cd8pos',
                         verbose = FALSE,
                         logfc.threshold = 0,
                         recorrect_umi =  FALSE)
  markers %<>% 
    tibble::rownames_to_column(., "symbol") %>% 
    mutate(biotype = sym2biotype_ids[symbol],
           ensemble = gm$SYMBOL$ENSEMBL[symbol],
           ident.1 = 'DC',
           ident.2 = 'DC_Cd4pos_Cd8pos', 
           tissue='DC_CD45_Tumor') %>% 
    relocate(., c('ensemble', 'biotype'))
  
  Idents(seuT_l$cd8) <- 'manual_anno'
  cd8_markers_memory <- FindMarkers(seuT_l$cd8, 
                         assay.use='RNA', test.use='wilcox',
                         ident.1='CD8_memory', ident.2='CD8_memory2',
                         verbose = FALSE,
                         logfc.threshold = 0,
                         recorrect_umi =  FALSE) %>% 
    tibble::rownames_to_column(., "symbol") %>% 
    mutate(biotype = sym2biotype_ids[symbol],
           ensemble = gm$SYMBOL$ENSEMBL[symbol],
           ident.1 = 'CD8_memory',
           ident.2 = 'CD8_memory2', 
           tissue='cd8memory_CD45_Tumor') %>% 
    relocate(., c('ensemble', 'biotype'))
  write.table(cd8_markers_memory, file="~/xfer/cd8memory_minuscd4_deg.csv", sep=",",
              col.names = T, row.names =F, quote = F)
  
  Idents(seuT_l$cd8) <- 'seurat_clusters'
  cd8_markers <- FindAllMarkers(seuT_l$cd8, 
                         assay.use='RNA', test.use='wilcox',
                         verbose = FALSE,
                         logfc.threshold = 0,
                         recorrect_umi =  FALSE) %>% 
    tibble::rownames_to_column(., "symbol") %>% 
    mutate(biotype = sym2biotype_ids[symbol],
           ensemble = gm$SYMBOL$ENSEMBL[symbol],
           tissue='cd8_CD45_Tumor') %>% 
    relocate(., c('ensemble', 'biotype'))
  write.table(cd8_markers, file="~/xfer/cd8_minuscd4_minusdoublet_deg.clusters.csv", sep=",",
              col.names = T, row.names =F, quote = F)
  pdf("~/xfer/cd8_minuscd4_minusdoublet.clusters.pdf")
  DimPlot(seuT_l$cd8, group.by='seurat_clusters', label=T, reduction='umap', raster=T)
  FeaturePlot(seuT_l$cd8, feature=c('Cd4', 'Cd8a'), label=T, reduction='umap', raster=T, pt.size=2)
  dev.off()
  
  
  
  uniqclusters <- unique(as.character(seuT_l$cd8$seurat_clusters))
  for(clusid in uniqclusters){
    Idents(seuT_l$cd8) <- 'seurat_clusters'
    seusub <- subset(seuT_l$cd8, ident=clusid)
    
    Idents(seusub) <- 'orig.ident'
    uniqid <- unique(as.character(Idents(seusub)))
    
    for(uid1 in uniqid){
      for(uid2 in uniqid[grep(uid1, uniqid, invert = T)]){
        if(!file.exists(paste0("~/xfer/cd8_cd4neg_dubneg.cl", clusid, ".comp", uid1, "_", uid2, ".csv"))){
          cd8_markers <- FindMarkers(seusub, 
                                     assay.use='RNA', test.use='wilcox',
                                     verbose = FALSE,
                                     ident.1=uid1, ident.2=uid2,
                                     logfc.threshold = 0,
                                     recorrect_umi =  FALSE) %>% 
            tibble::rownames_to_column(., "symbol") %>% 
            mutate(biotype = sym2biotype_ids[symbol],
                   ensemble = gm$SYMBOL$ENSEMBL[symbol],
                   ident.1=uid1,
                   ident.2=uid2,
                   tissue='cd8_CD45_Tumor') %>% 
            relocate(., c('ensemble', 'biotype'))
          write.table(cd8_markers, 
                      file=paste0("~/xfer/cd8_cd4neg_dubneg.cl", clusid, ".comp", uid1, "_", uid2, ".csv"), 
                      sep=",", col.names = T, row.names =F, quote = F)
          cat(paste0("xfer cd8_cd4neg_dubneg.cl", clusid, ".comp", uid1, "_", uid2, ".csv\n"))
        }
      }
    }
  }
  DimPlot(seuT_l$cd8, group.by='seurat_clusters', label=T, reduction='umap', raster=T)
  FeaturePlot(seuT_l$cd8, feature=c('Cd4', 'Cd8a'), label=T, reduction='umap', raster=T, pt.size=2)
  dev.off()
  
  itgaxneg_cells <- Cells(seuT_l$cd8)[GetAssayData(seuT_l$cd8, slot='counts')['Itgax',] == 0]
  cd8_itgaxneg_markers <- FindAllMarkers(subset(seuT_l$cd8, cells=itgaxneg_cells), 
                                assay.use='RNA', test.use='wilcox',
                                verbose = FALSE,
                                logfc.threshold = 0,
                                recorrect_umi =  FALSE) %>% 
    tibble::rownames_to_column(., "symbol") %>% 
    mutate(biotype = sym2biotype_ids[symbol],
           ensemble = gm$SYMBOL$ENSEMBL[symbol],
           tissue='cd8_itgaxNeg_CD45_Tumor') %>% 
    relocate(., c('ensemble', 'biotype'))
  write.table(cd8_itgaxneg_markers, file="~/xfer/cd8_minuscd4_iminusitgax_deg.csv", sep=",",
              col.names = T, row.names =F, quote = F)
}

features <- c('Cd8a', 'Cd4', 'Cd3d', 'Itgax', 'Itgam', 'Cd19', 'Klrb1')
pdf("~/xfer/cd45_tumor_cd8_and_dc.pdf", width = 17)
lapply(seuT_l, function(seuT_i){
  publicationDimPlot(seuT_i, grp='seurat_clusters', simplify_labels=TRUE, 
                     pt.size=1.2, aspect.ratio=1,
                     legend.text = element_text(size=10))
}) %>%
  cowplot::plot_grid(plotlist=., nrow=1)

lapply(seuT_l, function(seuT_i){
  publicationDimPlot(seuT_i, grp='manual_anno', simplify_labels=TRUE, 
                     pt.size=1.2, aspect.ratio=1,
                     legend.text = element_text(size=10))
}) %>%
  cowplot::plot_grid(plotlist=., nrow=1)

lapply(seuT_l, function(seuT_i){
  FeaturePlot_scCustom(seuT_i, features=features, raster=T, pt.size=3)
}) %>%
  cowplot::plot_grid(plotlist=., nrow=1)

dev.off()

#--- Identify clusters in CD45 ----
if(!exists("seul"))  seul <- readRDS(file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.rds"))
cd45_tcell_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tcells.rds"))
cd45_treg_map <- readRDS(file=file.path(datadir, "annotation_map", "cd45_tregs.rds"))
overwrite <- FALSE

seu_id <- 'LN'
seu_i <- seul[[seu_id]]
subset_of_seu <- FALSE
ident_id <- 'seurat_clusters'
Idents(seu_i) <- ident_id
idents=c('20', '27')
markers_all <- FindMarkers(seu_i, assay = "RNA", test.use='wilcox',
                           ident.1= idents[1], ident.2= idents[2],
                           verbose = FALSE,
                           logfc.threshold = 0,
                           recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
  tibble::rownames_to_column(., "symbol") %>% 
  mutate(biotype = sym2biotype_ids[symbol],
         ensemble = gm$SYMBOL$ENSEMBL[symbol],
         ident.1 = idents[1],
         ident.2 = idents[2], 
         comparison=ident_id,
         tissue=seu_id) %>% 
  relocate(., c('ensemble', 'biotype'))
write.table(markers_all, file=file.path("~/xfer", "LN_20vs27.csv"),
            sep=",", col.names = T, row.names = F, quote = F)

idents=c('26')
max_cells <- 5000
cells_sel <- sample(Cells(seu_i)[Idents(seu_i) != idents[1]], size=max_cells)
cells_sel <- c(cells_sel, Cells(seu_i)[Idents(seu_i) == idents[1]])
seu_j <- subset(seu_i, cells=cells_sel)
markers_all <- FindMarkers(seu_j, assay = "RNA", test.use='wilcox',
                           ident.1= idents[1],
                           verbose = FALSE,
                           logfc.threshold = 0,
                           recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
  tibble::rownames_to_column(., "symbol") %>% 
  mutate(biotype = sym2biotype_ids[symbol],
         ensemble = gm$SYMBOL$ENSEMBL[symbol],
         ident.1 = idents[1],
         ident.2 = 'All', 
         comparison=ident_id,
         tissue=seu_id) %>% 
  relocate(., c('ensemble', 'biotype'))
write.table(markers_all, file=file.path("~/xfer", "LN_26vsAll.csv"),
            sep=",", col.names = T, row.names = F, quote = F)


pdf("~/xfer/x.pdf", width = 15, height = 15)
DimPlot(seul$LN, group.by='seurat_clusters', reduction='umap', raster=T, label=T)
DimPlot(seul$Tumor, group.by='seurat_clusters', reduction='umap', raster=T, label=T)
FeaturePlot(seul$LN, features=c('Cd8a', 'Cd4', 'Cd3e', 'Itgam', 'Adgre1', 'Itgax', 'Cd19', 'Foxp3'), reduction='umap', raster=T)
FeaturePlot(seul$Tumor, features=c('Cd8a', 'Cd4', 'Cd3e', 'Itgam', 'Adgre1', 'Itgax', 'Cd19', 'Foxp3'), reduction='umap', raster=T)
dev.off()
####################################################
#### 2. Integrate with the 3d treatment samples ####
preproc_7d_f <- file.path(datadir, "seurat_obj", "seu_separate.sct.rds")
preproc_3d_f <- gsub("scrna_tdln_tumor_7d", "scrna_tdln_tumor", preproc_7d_f)
seu_b3d7d_f <- file.path(datadir, "seurat_obj", "seu_integ_B3D7D.mnn.rds")
seu_b3d7d_split_f <- file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.rds")
seu_3d7d_f <- file.path(datadir, "seurat_obj", "seu_integ_3D7D.mnn.rds")
seu_3d7d_split_f <- file.path(datadir, "seurat_obj", "seu_integ_3D7D.split.rds")
seu_cd45_7d_f <- file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.mnn.rds")
seu_cd45_7d_split_f <- file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.rds")
rm_untreated=FALSE

seu7d.list <- readRDS(preproc_7d_f)
seu3d.list <- readRDS(preproc_3d_f)
if(rm_untreated){
  seu_integ_f <- seu_3d7d_f
  seu_split_f <- seu_3d7d_split_f
  uid <- 'only3d7d'
} else {
  seu_integ_f <- seu_b3d7d_f
  seu_split_f <- seu_b3d7d_split_f
  uid <- 'Un3d7d'
}
if(rm_untreated){
  day7idx <- grep("ST2", names(seu7d.list), value = T)
  day7idx <- grep("Un$|PBS$", day7idx, value=T, invert = T)
  day7cd45_idx <- grep("CD45", names(seu7d.list), value = T)
  day3idx <- grep("Un$|PBS$", names(seu3d.list), value = T,invert = T)
} else {
  day7idx <- grep("ST2", names(seu7d.list), value = T)
  day7cd45_idx <- grep("CD45", names(seu7d.list), value = T)
  day3idx <- names(seu3d.list)
}

#--- a) Integrate 7d and 3d ST2 samples ----
seu.mnn <- RunFastMNN(object.list = c(seu7d.list[day7idx],
                                      seu3d.list[day3idx]), 
                      features = 2000, assay='SCT')
seu.cd45.mnn <- RunFastMNN(object.list = c(seu7d.list[day7cd45_idx]), 
                           features = 2000, assay='SCT')
rm(seu7d.list, seu3d.list); gc()

seu_st2_cd45 <- lapply(list("st2"=seu.mnn, "cd45"=seu.cd45.mnn), function(seu_i){
  seu_i <- RunUMAP(seu_i, reduction = "mnn", dims = 1:30, n.neighbors=30L,
                   min.dist = 0.1, return.model=TRUE)
  seu_i <- FindNeighbors(seu_i, reduction = "mnn", dims = 1:30)
  seu_i <- FindClusters(seu_i, resolution = 0.9, graph.name='SCT_snn')
  DefaultAssay(seu_i) <- 'RNA'
  seu_i <- NormalizeData(seu_i,
                         normalization.method = "LogNormalize") %>%
    FindVariableFeatures(., selection.method = "vst",
                         nfeatures = 3000, verbose = FALSE) %>% 
    ScaleData(.)
  seu_i <- RunPCA(seu_i, features=VariableFeatures(seu_i), ndims.print=1:30,
                  nfeatures.print=5, npcs=30)
  return(seu_i)
})
seu.mnn <- seu_st2_cd45$st2
seu.cd45.mnn <- seu_st2_cd45$cd45

saveRDS(seu.mnn, file=seu_integ_f)
saveRDS(seu.cd45.mnn, file=seu_cd45_7d_f)


#--- b) Split into Tumor and LN cohorts ----
seu7d.list <- readRDS(preproc_7d_f)
seu3d.list <- readRDS(preproc_3d_f)

ids <- names(c(seu7d.list[day7idx], seu3d.list[day3idx]))
cd45_ids <- names(seu7d.list[day7cd45_idx])

splitid <- grepl("LN", ids)
cd45splitid <- grepl("LN", cd45_ids)

seul.ln <- c(seu7d.list[day7idx], seu3d.list[day3idx])[which(splitid)]
seul.tumor <- c(seu7d.list[day7idx], seu3d.list[day3idx])[which(!splitid)]
seul.cd45.ln <- seu7d.list[day7cd45_idx][which(cd45splitid)]
seul.cd45.tumor <- seu7d.list[day7cd45_idx][which(!cd45splitid)]
rm(seu7d.list, seu3d.list); gc()

grps <- list("st2"=list("LN"=seul.ln, "Tumor"=seul.tumor),
             "cd45"=list("LN"=seul.cd45.ln, "Tumor"=seul.cd45.tumor))
seul_grps <- lapply(grps[2], function(seul_j){
  seul <- lapply(seul_j, function(seul_i){
    seul_i <- lapply(seul_i, DietSeurat, counts=T, data=T, scale.data=F)
    seu_i.mnn <- RunFastMNN(object.list = seul_i,
                            features = 2000, assay='SCT')
    
    seu_i.mnn <- RunUMAP(seu_i.mnn, reduction = "mnn", dims = 1:30, n.neighbors=30L,
                         min.dist = 0.1, return.model=TRUE)
    seu_i.mnn <- FindNeighbors(seu_i.mnn, reduction = "mnn", dims = 1:30)
    seu_i.mnn <- FindClusters(seu_i.mnn, resolution = 0.9, graph.name='SCT_snn')
    DefaultAssay(seu_i.mnn) <- 'RNA'
    seu_i.mnn <- NormalizeData(seu_i.mnn,
                               normalization.method = "LogNormalize") %>%
      FindVariableFeatures(., selection.method = "vst",
                           nfeatures = 3000, verbose = FALSE) %>% 
      ScaleData(.)
    seu_i.mnn <- RunPCA(seu_i.mnn, features=VariableFeatures(seu_i.mnn), ndims.print=1:30,
                        nfeatures.print=5, npcs=30)
    return(seu_i.mnn)
  })
  return(seul)
})
rm(seul.ln, seul.tumor); gc()
seul <- seul_grps$st2
seul.cd45 <- seul_grps$cd45
saveRDS(seul, file=seu_split_f)
saveRDS(seul.cd45, file=seu_cd45_7d_split_f)

#--- c) Assess the stability of the clusters ----
seul <- readRDS(file=seu_split_f)
subsetto <- function(obj, n=10000){
  ncells <- length(Cells(obj))
  n <- ceiling(ncells / n)
  grps <- sample(c(1:n), size=ncells, replace=T)
  
  obj_l <- lapply(c(1:n), function(i){
    print(paste0(i, "..."))
    cells_i <- Cells(obj)[which(grps==i)]
    subset(obj, cells=cells_i)
  })
  return(obj_l)
}
set.seed(1234)
seu_subs <- subsetto(seul[[1]], n=10000)
new_clusters_l <- list()
new_clusters_new_l <- list()

idx <-1
seu_i <- seu_subs[[idx]]

new_clusters_l[[idx]] <- scSHC::testClusters(GetAssayData(seu_i, assay='RNA', slot='counts'),
                                             as.character(seu_i$seurat_clusters),
                                             batch=as.character(seu_i$orig.ident),
                                             var.genes=as.character(seu_i@assays$RNA@var.features), 
                                             alpha=1,
                                             num_features=length(seu_i@assays$RNA@var.features),
                                             parallel=FALSE)
sapply(new_clusters_l, function(i) table(i[[1]]))

new_clusters_new_l[[idx]] <- scSHC::scSHC(GetAssayData(seu_i, assay='RNA', slot='counts'),
                              batch=as.character(seu_i$orig.ident),
                              alpha=0.25,
                              parallel=F)
sapply(new_clusters_new_l, function(i) table(i[[1]]))
saveRDS(new_clusters_new_l, file="new_clusters_new_l.rds")

table(new_clusters[[1]])
seu <- seul[[1]]
seu$scshc_clusters_1 <- seu$scshc_clusters_2 <- seu$scshc_clusters_3 <- seu$scshc_clusters_4 <- NA
seu@meta.data[names(new_clusters_new_l[[1]][[1]]),]$scshc_clusters_1 <- new_clusters_new_l[[1]][[1]]
seu@meta.data[names(new_clusters_new_l[[2]][[1]]),]$scshc_clusters_2 <- new_clusters_new_l[[2]][[1]]
seu@meta.data[names(new_clusters_new_l[[3]][[1]]),]$scshc_clusters_3 <- new_clusters_new_l[[3]][[1]]
seu@meta.data[names(new_clusters_new_l[[4]][[1]]),]$scshc_clusters_4 <- new_clusters_new_l[[4]][[1]]
pdf("~/xfer/x3.pdf", width = 15)
table(seu$seurat_clusters, seu$scshc_clusters_1)
table(seu$seurat_clusters, seu$scshc_clusters_2)
table(seu$seurat_clusters, seu$scshc_clusters_3)
table(seu$seurat_clusters, seu$scshc_clusters_4)
dp0 <- DimPlot(seu, reduction='umap', group.by='seurat_clusters', raster=T, label=T) + NoLegend()
dp1 <- DimPlot(seu, reduction='umap', group.by='scshc_clusters_1', raster=T, label=T) + NoLegend()
dp2 <- DimPlot(seu, reduction='umap', group.by='scshc_clusters_2', raster=T, label=T) + NoLegend()
dp3 <- DimPlot(seu, reduction='umap', group.by='scshc_clusters_3', raster=T, label=T) + NoLegend()
dp4 <- DimPlot(seu, reduction='umap', group.by='scshc_clusters_4', raster=T, label=T) + NoLegend()
plot(dp0 + dp1 + dp2)
plot(dp0 + dp3 + dp4)
plot(1)
dev.off()

#--- c) Annotate the split and integrated cohorts ----
projectils_dir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/spica'
rm(seu7d.list, seu3d.list); gc()
seul <- readRDS(file=seu_split_f)
seul.cd45 <- readRDS(file=seu_cd45_7d_split_f)
db <- 'dice'  # immgen or dice
bed.se <- readRDS(paste0("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/", db, ".rds"))
refobj_custom <- readRDS(file.path(projectils_dir, 'custom', 'custom_ilc2.mouse_atlasB.rds'))
DefaultAssay(refobj_custom) <- 'RNA'
visualize=TRUE

## ImmGen annotation using SingleR
seul_grps <- list("st2"=seul, "cd45"=seul.cd45)
seul_grps <- lapply(setNames(names(seul_grps), names(seul_grps)), function(grp_id){
  seul <- seul_grps[[grp_id]]
  seul <- lapply(setNames(names(seul),names(seul)), function(seu_id){
    seu <- seul[[seu_id]]
    DefaultAssay(seu) <- 'RNA'
    lbl <- 'label.fine'
    clus <- TRUE
    
    rds_file <- paste0(grp_id, ".", seu_id, "_", uid, "_", db, ".fine.cluster.rds")
    id <- gsub(paste0("^.*", db, "(.*).rds"), 
               paste0(seu_id, "_", db, "\\1"), rds_file)
    
    if(!file.exists(file.path(outdir, "annotation", rds_file))){
      print(paste0(" Annotating: ", rds_file))
      singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
                              ref=bed.se, 
                              assay.type.test=1, labels=bed.se[[lbl]],
                              de.method="wilcox", genes='sd',
                              clusters=if(clus) seu$seurat_clusters else NULL)
      saveRDS(singler_anno, file=file.path(outdir, "annotation", rds_file))
    } else {
      print(paste0(" Reading in annotation: ", rds_file))
      singler_anno <- readRDS(file.path(outdir, "annotation", rds_file))
    }
    
    if(clus){
      cluster_ids <- setNames(singler_anno$pruned.labels, 
                              as.character(rownames(singler_anno)))
      seu@meta.data[,id] <- cluster_ids[as.character(seu$seurat_clusters)]
    } else {
      seu@meta.data[,id] <- singler_anno$pruned.labels
    }
    return(seu)
  })
  
  seul <- lapply(seul, function(seu){
    id <- grep(paste0("", db, ".[fine|main]"), colnames(seu@meta.data), value=T)
    # Makes annotations simpler (e.g. "B Cells (B.Fo) -> B.Fo)
    id_map <- unique(seu@meta.data[,id]) %>%
      gsub("^.*\\(", "", .) %>%
      gsub("\\)", "", .) %>%
      gsub("(^.*?\\..*)\\..*", "\\1", .) %>%
      setNames(., unique(seu@meta.data[,id]))
    seu@meta.data[,paste0(db, ".simple")] <- id_map[seu@meta.data[,id]]
    seu@meta.data[,id] <- NULL
    return(seu)
  })
  
  return(seul)
})

## ProjecTILs annotation
seul_grps <- lapply(setNames(names(seul_grps), names(seul_grps)), function(grp_id){
  seul <- seul_grps[[grp_id]]
  seul <- lapply(setNames(names(seul), names(seul)), function(seuid){
    seu <- seul[[seuid]]
    DefaultAssay(seu) <- 'RNA'
    subset_size <- 3
    set.seed(1234)
    cell_idx <- split(Cells(seu), f = sample(c(1:subset_size), size=ncol(seu), replace = T))
    batch <- 1
    
    ptils_map <- lapply(cell_idx[batch:length(cell_idx)], function(cells_i){
      batchf <- file.path(outdir, "annotation", "projectils", 
                          paste0("batch", batch, ".", grp_id, ".", uid, ".", seuid, ".rds"))
      print(batchf)
      if(!file.exists(batchf)){
        print(paste0("batch: ", batch, "..."))
        seu_i <- subset(seu, cells=cells_i)
        seu_i <- ProjecTILs.classifier(seu_i, refobj_custom, ncores = 1, 
                                       split.by = "orig.ident", filter.cells=F)
        seu_lbl <- seu_i$functional.cluster
        saveRDS(seu_lbl, file=batchf)
      } else {
        print("reading in...")
        seu_lbl <- readRDS(file=batchf)
      }
      batch <<- batch+1
      return(seu_lbl)
    })
    
    ptils_map <- do.call(c, ptils_map)
    ptils_map <- ptils_map[which(ptils_map >= (length(ptils_map) * 0.0005))]
    names(ptils_map) <- gsub("^.*?\\.", "", names(ptils_map))
    
    seu$projectils.simple <- NA
    lowidents <- (table(ptils_map) < 50)
    if(any(lowidents)){
      ptils_map <- ptils_map[-which(ptils_map %in% names(which(lowidents)))]
    }
    
    seu@meta.data[names(ptils_map),'projectils.simple'] <- ptils_map
    
    return(seu)
  })
  return(seul)
})

seul <- seul_grps$st2
seul.cd45 <- seul_grps$cd45
saveRDS(seul, file=seu_split_f)
saveRDS(seul.cd45, file=seu_cd45_7d_split_f)

#--- d) Integrate the split cohorts ----
seul <- readRDS(file=seu_b3d7d_split_f)
seul <- lapply(seul, function(i){
  DefaultAssay(i) <- 'RNA'
  i
})
seul_reproc <- lapply(seul, function(seu_sub){
  seux <- preprocessSeu(seu_sub, ncount_min=0, ncount_max=Inf,
                        nfeature_min=0, nfeature_max=Inf,
                        mt_max=100, org='mouse', numpcs=30, getPCs=FALSE,
                        variable_features=2000, res=0.8)
  seux[['integrated']] <- NULL
  seux@assays$RNA@scale.data <- as.matrix(0)
  seux@assays$SCT@scale.data <- as.matrix(0)
  return(seux)
})

seu_i.mnn <- RunFastMNN(object.list = seul_reproc,
                        features = 2000, assay='SCT')

seu_i.mnn <- RunUMAP(seu_i.mnn, reduction = "mnn", dims = 1:30, n.neighbors=30L,
                     min.dist = 0.1, return.model=TRUE)
seu_i.mnn <- FindNeighbors(seu_i.mnn, reduction = "mnn", dims = 1:30)
seu_i.mnn <- FindClusters(seu_i.mnn, resolution = 0.9, graph.name='SCT_snn')
DefaultAssay(seu_i.mnn) <- 'RNA'
seu_i.mnn <- NormalizeData(seu_i.mnn,
                           normalization.method = "LogNormalize") %>%
  FindVariableFeatures(., selection.method = "vst",
                       nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(.)
seu_i.mnn <- RunPCA(seu_i.mnn, features=VariableFeatures(seu_i.mnn), ndims.print=1:30,
                    nfeatures.print=5, npcs=30)
return(seu_i.mnn)




#--- e) Visualize ----
.makeHighlightPlot <- function(seu, ident, n_rm=50){
  Idents(seu) <- ident
  highlight_plot <- lapply(na.omit(unique(seu@meta.data[,ident])), function(clid){
    scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = adjustcolor(c('red', 'red'), alpha.f=0.4),
                                        reduction='umap', background_color = "lightgray", raster=T) + 
      NoLegend() +
      ggtitle(clid)
  })
  return(highlight_plot)
}

.barplotfun <- function(x, colid='functional.cluster.custom', 
                        percentage=F, ncol=8){
  meltdat <- table(x@meta.data[,colid], x$orig.ident) 
  if(percentage) meltdat <- apply(meltdat, 2, function(i) i/sum(i))
  topids <- sort(rowSums(meltdat), decreasing = T)
  meltdat <- melt(meltdat) %>%
    mutate(Var1 = factor(Var1, levels=names(topids)))
  
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
  cols=setNames(c(c25[1:ncol], rep("grey", (length(topids)-(ncol)))),
                names(topids))
  
  ggplot(meltdat, 
         aes(x=Var2, fill=Var1, y=value, color=Var1)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) +
    theme_bw() +
    xlab("Samples") + ylab("Count") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}
## Visualize
pdf(file.path("~/xfer", "dimplot_annotation.split.sep.cd45.pdf"), width = 20, height = 12)
dps <- lapply(names(seul.cd45), function(seuid){
  seu <- seul.cd45[[seuid]]
  seu$orig.ident <- .relabelid(seu$orig.ident)
  lvls <- unique(seu$orig.ident)[order(gsub("B[12]_", "", unique(seu$orig.ident)))]
  seu$orig.ident <- factor(seu$orig.ident, levels=lvls)

  dp1 <- DimPlot(seu, group.by='seurat_clusters', reduction='umap', 
                 raster=T, label = T) + ggtitle(seuid) + NoLegend()
  dp2 <- DimPlot(seu, group.by='immgen.simple', reduction='umap', 
                 raster=T, label = T) + NoLegend()
  dp3 <- DimPlot(seu, group.by='dice.simple', reduction='umap', 
                 raster=T, label = T) + NoLegend()
  dp4 <- DimPlot(seu, group.by='projectils.simple', reduction='umap',
                 raster=T, label = T) + NoLegend()

  immgen_hp <- .makeHighlightPlot(seu, 'immgen.simple')
  ptils_hp <- .makeHighlightPlot(seu, 'projectils.simple')
  ident_hp <- .makeHighlightPlot(seu, 'orig.ident')

  bp1a <- .barplotfun(seu, 'immgen.simple', TRUE)
  bp1b <- .barplotfun(seu, 'immgen.simple', FALSE)
  bp2a <- .barplotfun(seu, 'projectils.simple', TRUE)
  bp2b <- .barplotfun(seu, 'projectils.simple', FALSE)
  
  plot(dp1 + dp2 + dp3)
  plot(dp1 + dp2 + dp4)
  plot(cowplot::plot_grid(plotlist = ident_hp, ncol=6))
  plot(cowplot::plot_grid(plotlist = immgen_hp, ncol=6))
  plot(cowplot::plot_grid(plotlist = ptils_hp, ncol=6))
  plot(cowplot::plot_grid(plotlist=list(bp1a, bp1b), nrow=2))
  plot(cowplot::plot_grid(plotlist=list(bp2a, bp2b), nrow=2))
})
dev.off()


################################################
#### 3. ScVelo trajectory/velocity analysis ####
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))

seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()
outrds <- file.path(datadir, "velocyto", "seu_velocyto.rds")
dir.create(file.path(outdir, "scVelo"), showWarnings = F)
DefaultAssay(seu) <- 'RNA'
make.small <- FALSE

## This analysis assumes that velocyto/0.17.17 was run on your 10x data already:
# gtf='/cluster/tools/data/commondata/cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf'
# mask_gtf='/cluster/projects/mcgahalab/ref/ucsc/repeat_masker/mm10_rmsk.gtf'
# 
# velocyto run10x \
# -m ${mask_gtf} \
# --samtools-threads 1 \
# --samtools-memory 20000 \
# ${sampledir} \
# ${gtf}
# 
# velo_obj <- readRDS(file=outrds)


## Read in the velocyto loom files
velocyto_dir <- file.path(datadir, "velocyto")
sample_ids <- list.files(velocyto_dir)
seu_velo_raw <- lapply(setNames(sample_ids,sample_ids), function(sid){
  f <- list.files(file.path(velocyto_dir, sid))
  loom.data <- ReadVelocity(file = file.path(velocyto_dir, sid, f)) %>%
    as.Seurat(.)
  loom.data$orig.ident <- sid
  loom.data <- RenameCells(loom.data, 
                           new.names=gsub("^M[ca]Gaha_.*__(.*)\\:([ACGT]*).*", "\\1_\\2", Cells(loom.data)) %>%
                             gsub("^T", "Tumor", .))
  return(loom.data)
})

## Preprocess the seurat velocyto files
seu_velo <- merge(seu_velo_raw[[1]], seu_velo_raw[-1], 
                  project = 'ST2_TDLN_LN')
seu_velo_l <- list("Merged"=seu_velo,
                   "Tumor"=seu_velo_raw[grep("^Tumor", names(seu_velo_raw), value=T)],
                   "LN"=seu_velo_raw[grep("^LN", names(seu_velo_raw), value=T)])
rm(seu_velo, seu_velo_raw); gc()
# seu_velo <- seu_velo_raw[[1]]
paths <- list("Treg"=c('1'),
              "Tall"=c('1', '3', '5'),
              "All"=c('All'))

Idents(seu) <- 'monocle3_partitions'
x <- lapply(names(seu_velo_l)[-1], function(seu_velo_id){
  print(seu_velo_id)
  seu_velo <- seu_velo_l[[seu_velo_id]]
  if(class(seu_velo) == 'list'){
    seu_velo <- merge(seu_velo[[1]], seu_velo[-1], 
                      project = paste0('ST2_TDLN_LN.', seu_velo_id))
  }
  DefaultAssay(object = seu_velo) <- "spliced"
  if(make.small){
    seu_velo <- subsetSeu(seu_velo,colid='orig.ident', n=3000) %>%
      subset(., nFeature_spliced>100)
  }
  
  lapply(names(paths), function(path_j){
    print(paste0(seu_velo_id, " - ", path_j))
    seu_sub <- if(path_j=='All'){
      seu
    } else { 
      subset(seu, ident=paths[[path_j]])
    }
    seu_small <- subset(seu_sub,
                        cells=Cells(seu_velo))
    seu_velo_small <- subset(seu_velo, cells=Cells(seu_small))
    
    # Normalize and cluster
    set.seed(seed)
    seu_velo_small[["RNA"]] <- seu_velo_small[["spliced"]]
    seu_velo_small <- seu_velo_small %>%
      SCTransform(.) %>%
      RunPCA %>%
      RunUMAP(., dims=1:30) %>%
      FindNeighbors(., dims = 1:30) %>%
      FindClusters(.)
    DefaultAssay(seu_velo_small) <- "RNA"
    seu_velo_small[['umap_orig']] <- CreateDimReducObject(embeddings = Embeddings(seu_small, reduction='umap'),
                                                          key = 'UMAP_', assay = 'RNA')
    seu_velo_small$functional_cluster <- as.character(seu_small$manual_anno)
    
    print("saving...")
    SeuratDisk::SaveH5Seurat(seu_velo_small, 
                             filename = file.path(outdir, "scVelo", 
                                                  paste0("seu_velo.", seu_velo_id, ".", path_j, ".h5Seurat")),
                             overwrite=T)
    cwd <- getwd()
    setwd(file.path(outdir, "scVelo"))
    SeuratDisk::Convert(source = paste0("seu_velo.", seu_velo_id, ".", path_j, ".h5Seurat"), 
                        dest = "h5ad", overwrite=T)
    setwd(cwd)
    return(NULL)
  })
})



#--- a) TRegs or CD8 only ----
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
celltype <- c('cd8'="^T_CD8")
celltype <- c('tregs'="^TReg")

seul_celltype <- lapply(seul, function(seu){
  Idents(seu) <- 'manual_anno'
  seu <- subset(seu, ident=setdiff(unique(grep(celltype, Idents(seu), value=T, ignore.case=T)),
                                   unique(grep("remove", Idents(seu), value=T, ignore.case = T))))
  return(seu)
})

## Read in the velocyto loom files
velocyto_dir <- file.path(datadir, "velocyto")
dir.create(file.path(velocyto_dir, "merge"), showWarnings = F)
sample_ids <- list.files(velocyto_dir) %>%
  grep("^CD45", ., value=T, invert = T) %>%
  grep("^merge", ., value=T, invert=T)
new_ids <- .relabelid(sample_ids)

if(!file.exists(file.path(velocyto_dir, "merge", "seu_velo.rds"))){
  print("Loading individual velocyto files and merging into one seurat object...")
  seu_velo_raw <- lapply(setNames(sample_ids,new_ids), function(sid){
    print(sid)
    f <- list.files(file.path(velocyto_dir, sid))
    loom.data <- ReadVelocity(file = file.path(velocyto_dir, sid, f)) %>%
      as.Seurat(.)
    loom.data$orig.ident <- sid
    loom.data <- RenameCells(loom.data, 
                             new.names=gsub("^M[ca]Gaha_.*__(.*)\\:([ACGT]*).*", paste0(sid, "_\\2"), Cells(loom.data)) %>%
                               gsub("^T", "Tumor", .) %>% 
                               gsub("Tumorumor", "Tumor", .))  
    
    return(loom.data)
  })
  
  ## Preprocess the seurat velocyto files
  seu_velo <- merge(seu_velo_raw[[1]], seu_velo_raw[-1], 
                    project = 'ST2_TDLN_LN')
  saveRDS(seu_velo, file=file.path(velocyto_dir, "merge", "seu_velo.rds"))
} else {
  print("Reading in merged velocyto seurat object...")
  seu_velo <- readRDS(file=file.path(velocyto_dir, "merge", "seu_velo.rds"))
}


unionize_untreated <- FALSE
for(seu_id in names(seul_celltype)){
  seu_j <- seul_celltype[[seu_id]]
  seu_j$batch <- gsub("_.*", "", seu_j$orig.ident)
  # Idents(seu_j) <- 'treg_manual_anno'
  # seu_j <- subset(seu_j, ident=setdiff(na.omit(as.character(unique(Idents(seu_j)))),
  #                                      unique(grep("remove", Idents(seu_j), value=T, ignore.case = T))))
  
  
  for(batchid in c("B1", "B2")){
    dims_use <- 1:20
    Idents(seu_j) <- 'batch'
    seu_jbatch <- subset(seu_j, ident=batchid) %>%
      RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
      RunTSNE(., dims=dims_use, reduction='mnn') #%>%
      # runPacmap(., reduction='mnn')
    
    # pdf("~/xfer/x.pdf")
    # DimPlot(seu_jbatch, group.by='orig.ident', label=T, raster=T, reduction='pacmap')
    # DimPlot(seu_jbatch, group.by='manual_anno', label=T, raster=T, reduction='pacmap')
    # dev.off()
    
    #### Create velocyto seurat object
    DefaultAssay(object = seu_velo) <- "spliced"
    cells_sample <- sample(Cells(seu_jbatch), 
                           size=min(ncol(seu_jbatch), 10000), replace=F)
    # setdiff(unique(gsub("_[ACGT]*$", "", Cells(seu_velo))), unique(gsub("_[ACGT]*$", "", cells_sample))) %>% table
    # setdiff(unique(gsub("_[ACGT]*$", "", cells_sample)), unique(gsub("_[ACGT]*$", "", Cells(seu_velo)))) %>% table
    seu_j_subset <- subset(seu_jbatch, cells=cells_sample)
    seu_velo_j <- subset(seu_velo, cells=cells_sample)
    
    # Normalize and cluster
    set.seed(1234)
    dims <- 1:20
    seu_velo_j[["RNA"]] <- seu_velo_j[["spliced"]]
    seu_velo_j <- seu_velo_j %>%
      SCTransform(., vst.flavor = "v2") %>%
      RunPCA %>%
      FindNeighbors(., dims = dims) %>%
      FindClusters(.) %>%
      RunUMAP(., dims=dims)
    DefaultAssay(seu_velo_j) <- "RNA"
    seu_velo_j[['umap_orig']] <- CreateDimReducObject(embeddings = Embeddings(seu_j_subset, reduction='umap')[Cells(seu_velo_j),],
                                                      key = 'UMAP_', assay = 'RNA')
    seu_velo_j[['pacmap_orig']] <- CreateDimReducObject(embeddings = Embeddings(seu_j_subset, reduction='pacmap')[Cells(seu_velo_j),],
                                                        key = 'PacMAP_', assay = 'RNA')
    seu_velo_j$manual_anno <- as.character(seu_j_subset$manual_anno[Cells(seu_velo_j)]) # $manual_anno
    seu_velo_j$orig.ident <- as.character(seu_j_subset$orig.ident[Cells(seu_velo_j)]) # $manual_anno
    
    if(unionize_untreated){
      wt_un_idx <- grep("WT_Un$", as.character(seu_velo_j$orig.ident))
      ko_un_idx <- grep("KO_Un$", as.character(seu_velo_j$orig.ident))
      seu_velo_j$orig.ident[wt_un_idx] <- 'WT_Un'
      seu_velo_j$orig.ident[ko_un_idx] <- 'KO_Un'
    }
    
    for(subset_type in c('WT', 'KO', '.*')){
      print(paste0(seu_id, "_", subset_type, "..."))
      Idents(seu_velo_j) <- 'orig.ident'
      cells_sub <- grep(paste0("_?", subset_type, "_"), as.character(seu_velo_j$orig.ident), 
                        ignore.case = T)
      seu_velo_j_sub <- subset(seu_velo_j, cells=Cells(seu_velo_j)[cells_sub])
      subset_type <- ifelse(subset_type == '.*', 'all', subset_type)
      
      print("saving...")
      dir.create(file.path(outdir, "scVelo", names(celltype)), showWarnings = F)
      SeuratDisk::SaveH5Seurat(seu_velo_j_sub, 
                               filename = file.path(outdir, "scVelo", names(celltype),
                                                    paste0("seu_velo.", seu_id, ".", batchid, ".", subset_type, ".h5Seurat")),
                               overwrite=T)
      cwd <- getwd()
      setwd(file.path(outdir, "scVelo", names(celltype)))
      SeuratDisk::Convert(source = paste0("seu_velo.", seu_id, ".", batchid, ".", subset_type, ".h5Seurat"), 
                          dest = "h5ad", overwrite=T)
      print(paste0("seu_velo.", seu_id, ".", batchid, ".", subset_type, ".h5Seurat"))
      setwd(cwd)
    }
  }
}


#--- b) TRegs BGPLVM analysis ----
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
celltype <- c('cd8'="^T_CD8")
celltype <- c('tregs'="^TReg")

prepSubData <- function(seurat_obj, genes = NULL){
  sub_data = seurat_obj
  DefaultAssay(sub_data) <- 'RNA'
  sub_data <- NormalizeData(object = sub_data, 
                            scale.factor = 10000, display.progress = F)
  # sub_data = Seurat::FindVariableFeatures(sub_data, display.progress = F, 
  #                                         num.bin = 100, 
  #                                         binning.method = "equal_frequency")
  # 
  exp_mat = if(is.null(genes)){
    t(as.matrix(sub_data@assays$RNA@data[sub_data@assays$RNA@var.features,]))
  } else{
    t(as.matrix(sub_data@assays$RNA@data[genes,]))
  }
  save_mat = merge(sub_data@meta.data[,c("manual_anno", "seurat_clusters")],
                   as.matrix(exp_mat), by = 0)
  rownames(save_mat) = save_mat[,1]
  save_mat = save_mat[,-1]
  
  return(list("seurat" = sub_data, "mat" = save_mat))
}


treg_ids <- sapply(seul, function(seu_i){
  grep(celltype, unique(seu_i$manual_anno), value=T, ignore.case=T) %>%
    grep("remove|cycling", ., value=T, invert = T)
}) %>% as.character %>% unique
# treg_ids <- setNames(scales::hue_pal()(length(treg_ids)), treg_ids)
treg_ids <- setNames(c("#0083da", "#ff8281", "#f22300", "#4d6b9a",
                       "#be7fdc",  "#006d2c", "#00441b", '#c7e9c0', '#f7fcf5'),
                     c("TReg_Central", "TReg_NLT.NLTlike", "TReg_NLT", "TReg_NLTlike.STAT1",
                       "TReg_Effector", "TReg_NLTlike.Effector.Central_001",
                       "TReg_NLTlike.Effector.Central_002", 'TReg_NLTlike.STAT1_001',
                       'TReg_NLTlike.STAT1_002'))

seul_celltype <- lapply(seul, function(seu){
  sampleids <- unique(seu$orig.ident)
  lapply(sampleids, function(sampleid){
    print(sampleid)
    # Subset for non-cycling TRegs
    Idents(seu) <- 'manual_anno'
    seu <- subset(seu, ident=setdiff(unique(grep(celltype, Idents(seu), value=T, ignore.case=T)),
                                     unique(grep("remove|cycling", Idents(seu), value=T, ignore.case = T))))
    
    # Subset for a specific sample
    Idents(seu) <- 'orig.ident'
    seu <- prepSubData(subset(seu, ident=sampleid))
    write.csv(seu$mat, file = paste0("./results/bgplvm/", sampleid, ".varfeat.csv"), 
              row.names = T, quote = F, col.names = T)
    return(seu)
  })
})


{
  # run GPy
}

allids <- as.character(unique(seul$LN$orig.ident))
lv1_idx <- c('B1_LN_WT_3d')
lv0_idx <- setdiff(allids, lv1_idx)
# allids <- c('KO_3d', 'KO_7d', 'WT_3d', 'WT_7d')
lvdf_l <- lapply(setNames(allids, allids), function(fileid){
  lv <- read.csv(paste0("./results/bgplvm/", fileid, ".latent_var.csv"), header = T, row.names = 1)
  ard <- read.csv(paste0("./results/bgplvm/", fileid, ".ARD.csv"), header = T, row.names = 1)
  lv_i <- ifelse(any(fileid %in% lv1_idx), 'LV1', 'LV0')
  lv_2 <- setdiff(c('LV0', 'LV1'), lv_i)
  
  seu_sub <- prepSubData(subset(seul$LN, cells=lv$cell))

  lvdf <- seu_sub$seurat@meta.data[,'manual_anno', drop=F] %>% 
    tibble::rownames_to_column(., "cell") %>% 
    left_join(., lv, by='cell')
  logexp <- as.matrix(seu_sub$seurat@assays$RNA@data[,lv$cell])
  
  avg_grp_lv0 <- sapply(split(lvdf[,lv_i], lvdf$manual_anno), mean)
  if(all(avg_grp_lv0['TReg_Effector'] >= avg_grp_lv0)){
    for(id in grep("^LV", colnames(lvdf), value=T)){
      lvdf[,id] <- -1 * lvdf[,id]
    }
  }
  eff_ids <- grep("Effector", names(avg_grp_lv0), ignore.case = T, value=T)
  nlt_ids <- grep("NLT", names(avg_grp_lv0), ignore.case = T, value=T)
  if(!mean(avg_grp_lv0[eff_ids]) < mean(avg_grp_lv0[nlt_ids])){
    lvdf[,lv_i] <- -1*lvdf[,lv_i]
  }
  
  
  ggp <- ggplot(lvdf, aes_string(x = lv_i, y = lv_2, colour = 'manual_anno'))+
    scale_color_manual(values=treg_ids) +
    geom_point(shape = 19, size = 0.17)+
    theme_classic()+
    theme(legend.position = "none", plot.margin = unit(c(0,0.1,0,0), "cm")) + 
    ggtitle(fileid)
  
  ggl <- ggplot(lvdf, aes_string(x = lv_i, fill = 'manual_anno'))+
    scale_fill_manual(values=treg_ids) +
    geom_density(alpha = 0.5, colour = "grey25")+
    theme_cowplot()+
    theme(legend.position = "bottom", 
          legend.key.width = unit(0.45, "lines"),
          legend.key.height = unit(0.45, "lines"),
          axis.text = element_text(size = 4.5, colour = "black"),
          axis.title = element_text(size = 5.5, colour = "black"),
          legend.text = element_text(size = 5.5),
          legend.title = element_text(size = 5.5),
          plot.margin = unit(c(0,0.1,0,0), "cm"), 
          legend.background = element_blank()) +
    ylim(0, 1.5)
  pdf(paste0("~/xfer/", fileid, ".lv.pdf"))
  print(cowplot::plot_grid(ggp, ggl, ncol=1))
  dev.off()
  
  lvdf$LV_sel <- lvdf[,lv_i]
  return(lvdf)
})
lvdf_l <- c(lvdf_l,
            list("BX_LN_KO_7d"=lvdf_l[['B2_LN_WT_7d']],
                 "BX_LN_WT_7d"=lvdf_l[['B2_LN_KO_7d']],
                 "BX_LN_KO_Un"=lvdf_l[['B2_LN_KO_Un']],
                 "BX_LN_WT_Un"=lvdf_l[['B2_LN_WT_Un']]))
  
sapply(allids, function(fileid) cat(paste0("xfer ", fileid, ".lv.pdf\n")))

comparisons <- list("KO_7d.vs.3d"=c('B2_LN_KO_7d', 'B1_LN_KO_3d'),
                    "WT_7d.vs.3d"=c('B2_LN_WT_7d', 'B1_LN_WT_3d'),
                    "7d_WT.vs.KO"=c('B2_LN_WT_7d', 'B2_LN_KO_7d'),
                    "3d_WT.vs.KO"=c('B1_LN_WT_3d', 'B1_LN_KO_3d'),
                    "swapKO_7d.vs.3d"=c('BX_LN_KO_7d', 'B1_LN_KO_3d'),
                    "swapWT_7d.vs.3d"=c('BX_LN_WT_7d', 'B1_LN_WT_3d'),
                    "swap7d_WT.vs.KO"=c('BX_LN_WT_7d', 'BX_LN_KO_7d'))
delta <- lapply(comparisons, function(comp_i){
  simp_i <- gsub("^.*LN_", "", comp_i)
  comp_j <- c(comp_i, gsub("7d|3d", "Un", comp_i)) %>%
    setNames(., .)
  
  ## 1. Get the LV dat stratifying cell types for each sample
  .getLvDat <- function(dat, lv_lvl){
    dat %>%
      select(manual_anno, lv_lvl)  %>%
      magrittr::set_colnames(c("celltype", "LV")) %>%
      mutate(LV = scales::rescale(LV, to=c(0,1)))
  }
  lv_idx <- sapply(comp_j, function(i) ifelse(any(i %in% lv1_idx), 'LV1', 'LV0'))
  lv_grps <- lapply(comp_j, function(i) .getLvDat(lvdf_l[[i]], lv_idx[i]))
  
  ## 2. Get the LV dat stratifying cell types for each sample, invert if effector is low on LV
  lv_grps_spl <- lapply(lv_grps, function(i) split(i, f=i$celltype))
  grps_mean <- lapply(lv_grps_spl, function(samplei) sapply(samplei, function(ct) mean(ct$LV)))
  grps_dir <- sapply(grps_mean, function(mean_cts){
    # TRUE if effector group is on the low end of the LV spectrum and needs to inverted
    if_else(mean(mean_cts[grep("effector$", names(mean_cts), ignore.case=T)])<0.5,
            TRUE, FALSE)
  })
  lv_grps_spl <- lapply(names(lv_grps), function(id){
    i <- lv_grps[[id]]
    if(grps_dir[id]) i$LV <- 1-i$LV
    split(i, f=i$celltype)
  }) %>% setNames(., names(lv_grps))
  
  ## 3. Get the null distance between the celltypes for two samples
  .getNullKS <- function(lv_sample1, lv_sample2, nsims=10000, size=0.4){
    null_deltaD <- sapply(c(1:nsims), function(x){
      ks1 <- ks.test(sample(lv_sample1$LV, size=(nrow(lv_sample1)*size)),
                     sample(lv_sample1$LV, size=(nrow(lv_sample1)*size)))$statistic
      ks2 <- ks.test(sample(lv_sample2$LV, size=(nrow(lv_sample2)*size)),
                     sample(lv_sample2$LV, size=(nrow(lv_sample2)*size)))$statistic
      ks1-ks2
    })
    return(null_deltaD)
  }
  null_1 <- .getNullKS(lv_grps[[comp_i[1]]], lv_grps[[gsub("7d|3d", "Un", comp_i[1])]])
  null_2 <- .getNullKS(lv_grps[[comp_i[2]]], lv_grps[[gsub("7d|3d", "Un", comp_i[2])]])
  null_d12 <- null_1 - null_2
 
  ## 4. Calculate differential
  .getDiff <- function(lv_grp1_spl, lv_grp2_spl, null,
                       id1='T', id2='U'){
    uniq_combos <- intersect(names(lv_grp1_spl), names(lv_grp2_spl)) %>%
      combn(., m=2)
    
    res <- sapply(list(lv_grp1_spl, lv_grp2_spl), function(lv_grpX_spl){
      apply(uniq_combos, 2, function(uniq_combo_j){
        ksval <- ks.test(lv_grpX_spl[[uniq_combo_j[1]]]$LV,
                         lv_grpX_spl[[uniq_combo_j[2]]]$LV)
        c('D'=round(ksval$statistic,3))
      }) %>% t
    }) %>% 
      as.data.frame %>%
      mutate(delta=V1 - V2,
             celltype1=uniq_combos[1,],
             celltype2=uniq_combos[2,]) %>%
      dplyr::rename_with(., ~paste0("D.",c(id1, id2)), .cols=c('V1', 'V2')) %>%
      dplyr::relocate(., c('celltype1', 'celltype2'))
    res$p <- round((sapply(res$delta, function(d) {
      min(c(sum(d > null),
            sum(d < null)))
    }) / length(null)) * 2, 5)
    res$p.adj <- round(p.adjust(res$p), 5)
    return(res)
  }
  res_1 <- .getDiff(lv_grps_spl[[comp_i[1]]], 
                    lv_grps_spl[[gsub("7d|3d", "Un", comp_i[1])]],
                    null_1)
  res_2 <- .getDiff(lv_grps_spl[[comp_i[2]]], 
                    lv_grps_spl[[gsub("7d|3d", "Un", comp_i[2])]],
                    null_2)
  res_full <- full_join(res_1, res_2, by=c('celltype1', 'celltype2'),
                        suffix=paste0(".", simp_i)) 
  res_full$delta <- res_full[,paste0("delta.", simp_i[1])] - res_full[,paste0("delta.", simp_i[2])]
  res_full$p <- round((sapply(res_full$delta, function(d) {
    min(c(sum(d > null_d12),
          sum(d < null_d12)))
  }) / length(null_d12)) * 2, 5)
  res_full$p.adj <- round(p.adjust(res_full$p), 5)
  return(res_full)
})
  # res <- res %>% mutate(
  #   mean_ct1.1=mean_ranks[[1]][celltype1,'mean'],
  #   mean_ct2.1=mean_ranks[[1]][celltype2,'mean'],
  #   rank_ct1.1=mean_ranks[[1]][celltype1,'rank'],
  #   rank_ct2.1=mean_ranks[[1]][celltype2,'rank'],
  #   mean_ct1.2=mean_ranks[[2]][celltype1,'mean'],
  #   mean_ct2.2=mean_ranks[[2]][celltype2,'mean'],
  #   rank_ct1.2=mean_ranks[[2]][celltype1,'rank'],
  #   rank_ct2.2=mean_ranks[[2]][celltype2,'rank']
  # ) %>% 
  #   dplyr::rename_with(., ~gsub("\\.1", paste0(".", comp_i[1]), .)) %>%
  #   dplyr::rename_with(., ~gsub("\\.2", paste0(".", comp_i[2]), .))
  # 
  # return(res)
# })
lapply(names(delta), function(id){
  write.table(delta[[id]], file=file.path("~/xfer", paste0("deltaLV.", id, ".csv")),
            col.names=T, row.names=F, sep=",", quote=F)
  cat(paste0("xfer deltaLV.", id, ".csv\n"))
})


x <- lapply(lvdf_l, function(lv1){
  lv1_spl <- with(lv1, split(LV_sel, manual_anno))
  
  lapply(lvdf_l, function(lv2){
    lv2_spl <- with(lv2, split(LV_sel, manual_anno))
    sapply(intersect(names(lv1_spl), names(lv2_spl)), function(iid){
      data.frame("p"=wilcox.test(lv1_spl[[iid]], lv2_spl[[iid]])$p.value,
                 "fc"=mean(lv1_spl[[iid]])/mean(lv2_spl[[iid]]),
                 "mean_lv1"=mean(lv1_spl[[iid]]),
                 "mean_lv2"=mean(lv2_spl[[iid]]))
                        
    })
  })
})
sapply(lvdf_l, function(lvi){
  unique(lvi$manual_anno)
})

sampleid <- c('B2_.*_KO_7d') #B1_LN_WT_3d
fileid <- gsub("^.*_(WT|KO)", "\\1", sampleid)

lv <- read.csv(paste0("./results/bgplvm/", fileid, ".LN_latent_var.csv"), header = T, row.names = 1)
ard <- read.csv(paste0("./results/bgplvm/", fileid, ".LN_ARD.csv"), header = T, row.names = 1)
lv_i <- 'LV0'

lvdf <- seul_celltype$LN$seurat@meta.data[,'manual_anno', drop=F] %>% 
  tibble::rownames_to_column(., "cell") %>% 
  left_join(., lv, by='cell')
logexp <- as.matrix(seul_celltype$LN$seurat@assays$RNA@data[,lv$cell])

avg_grp_lv0 <- sapply(split(lvdf$LV0, lvdf$manual_anno), mean)
if(all(avg_grp_lv0['TReg_Effector'] >= avg_grp_lv0)){
  for(id in grep("^LV", colnames(lvdf), value=T)){
    lvdf[,id] <- -1 * lvdf[,id]
  }
}

if(visualize){
  ggp <- ggplot(lvdf, aes(x = LV0, y = LV1, colour = manual_anno))+
    scale_color_manual(values=treg_ids) +
    geom_point(shape = 19, size = 0.17)+
    theme_classic()+
    theme(legend.position = "none", plot.margin = unit(c(0,0.1,0,0), "cm")) + 
    ggtitle(fileid)
  
  ggl <- ggplot(lvdf, aes_string(x = lv_i, fill = 'manual_anno'))+
    scale_fill_manual(values=treg_ids) +
    geom_density(alpha = 0.5, colour = "grey25")+
    theme_cowplot()+
    theme(legend.position = "bottom", 
          legend.key.width = unit(0.45, "lines"),
          legend.key.height = unit(0.45, "lines"),
          axis.text = element_text(size = 4.5, colour = "black"),
          axis.title = element_text(size = 5.5, colour = "black"),
          legend.text = element_text(size = 5.5),
          legend.title = element_text(size = 5.5),
          plot.margin = unit(c(0,0.1,0,0), "cm"), 
          legend.background = element_blank()) +
    ylim(0, 1.5)
  pdf(paste0("~/xfer/", fileid, ".lv.pdf"))
  print(cowplot::plot_grid(ggp, ggl, ncol=1))
  dev.off()
  print(paste0( fileid, ".lv.pdf"))
}


## Run switchde
de <- switchde::switchde(logexp[seul_celltype$LN$seurat@assays$RNA@var.features,], lv[,lv_i]) 
# de <- readRDS("/Users/rquevedo/Projects/mcgaha/st2_il33/results/mouse/scrna_tumor_ln_7d/deg/tmp/de.rds")
# de <- readRDS("~/xfer/de.rds")

sig_g = de$gene[de$qval < 0.05]

if(!exists("msig_ds")){
  gprofiler_dir <- '/cluster/projects/mcgahalab/ref/gprofiler'
  gmt <- GSA::GSA.read.gmt(file.path(gprofiler_dir, 'gprofiler_full_mmusculus.ENSG.gmt'))
  gprof_ds <-setNames(gmt$genesets, 
                      paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
  gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]
  
  
  msig_ds <- lapply(names(gprof_ds), function(sublvl){
    data.frame("gs_name"=sublvl,
               "entrez_gene"=gprof_ds[[sublvl]])
  }) %>% do.call(rbind,.)
  msig_ds$classification <- gsub(":.*", "", msig_ds$gs_name)
}
ora <- clusterProfiler::enricher(gene = na.omit(gm$SYMBOL$ENSEMBL[sig_g]), TERM2GENE = msig_ds, maxGSSize=5000)@result

de_df <- as.data.frame(de) %>% 
  tibble::column_to_rownames('gene')
resl <- sapply(ora$geneID, function(genes_i){
  genes_i <- gm$ENSEMBL$SYMBOL[strsplit(genes_i, split="/")[[1]]]
  res <- sapply(c(mean, sd, median, mad), function(fun) 
    apply(de_df[genes_i,], 2, function(i) fun(na.omit(i)))) %>%
    round(., 3) %>% 
    magrittr::set_colnames(., c('mean', 'sd', 'median', 'mad'))
  resv <- as.data.frame(res[c('k', 't0'),]) %>% unlist
  names(resv) <- gsub("1$", ".k", names(resv)) %>%
    gsub("2$", ".t0", .)
  return(resv)
})
ora_res <- t(resl) %>% 
  as.data.frame %>%
  cbind(ora, .) %>% 
  select(-c(geneID, Description)) %>%
  tibble::remove_rownames(.)


# ora_l <- ora_resl <- list()
ora_l[[fileid]] <- ora
ora_resl[[fileid]] <- ora_res

for(id in names(ora_resl)){
  write.table(ora_resl[[id]], 
              file=file.path("~/xfer", paste0(id, ".ora.csv")),
              sep=",", col.names = T, row.names = F)
  print(paste0(id, ".ora.csv"))
}




idmerge <- lapply(ora_resl, function(i){
  i %>% 
    select(ID, median.t0, p.adjust)
}) %>%
  purrr::reduce(full_join, by='ID')
ids <- lapply(ora_resl, function(i){
  i %>% 
    filter(p.adjust < 0.001) %>%
    pull(ID)
})
setdiff(ids[[1]], ids[[2]]) %>% 
  head(., 20)



test = gprofiler2::gost(query=as.character(sig_g), 
                        organism = "mmusculus", 
                        correction_method = "fdr")
gostplot(test, capped = FALSE, interactive=FALSE)

downsample_list = list()
downsize_n <- min(table(seul_celltype$LN$seurat$manual_anno)) * 0.8
for(i in 1:50){
  print(paste0("s", i))
  downsample_cells = unlist(tapply(Cells(seul_celltype$LN$seurat), 
                                   seul_celltype$LN$seurat$manual_anno, 
                                   sample, size = downsize_n, replace = F))
  
  subexp <- logexp[seul_celltype$LN$seurat@assays$RNA@var.features, downsample_cells]
  subexp = subexp[rowSums(subexp)>0,]
  downsample_list[[paste0("s", i)]] = switchde::switchde(subexp, lv[downsample_cells,"LV0"])
}
gene_names <- data.frame('gene'=rownames(seul_celltype$LN$seurat))
switchde_downsample = lapply(downsample_list, function(x) merge(x, gene_names, by = 'gene'))

sig_g = de$gene[de$qval < 0.05]
test = gProfileR::gprofiler(as.character(sig_g), 
                            organism = "mmusculus", 
                            correction_method = "fdr",
                            hier_filtering = "moderate", 
                            max_p_value = 0.05)

getTabSummary = function(var = "t0"){
  vals = data.frame(row.names = genes_use,
                    "skin" = de$skin_proj_10xGenes_proj[match(genes_use,switchde_list_10x_notscaled$skin_proj_10xGenes_proj$gene),var],
                    "colon" = de$colon_10xGenes_proj[match(genes_use,switchde_list_10x_notscaled$colon_10xGenes_proj$gene),var],
                    "skin_mean" = NA,
                    "colon_mean" = NA,
                    "colon_mean_more" = NA)
  
  for(g in genes_use){
    vals$skin_mean[rownames(vals)==g] = median(unlist(sapply(switchde_skindownsample, function(x) x[x$gene==g,var])), 
                                               na.rm = T)
    vals$colon_mean[rownames(vals)==g] = median(unlist(sapply(switchde_colondownsample, function(x) x[x$gene==g,var])), 
                                                na.rm = T)
    vals$colon_mean_more[rownames(vals)==g] = median(unlist(sapply(switchde_colondownsample_more, 
                                                                   function(x) x[x$gene==g,var])), na.rm = T)
  }
  
  return(vals)
}

#################################################################
#### 4. Trajectory analysis - MONOCLE3 ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "2_seu_annoX.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

DefaultAssay(seu) <- 'integrated'
seu <- PrepSCTFindMarkers(seu)
seu_sub <- seu 

# data <- as(as.matrix(seu_sub@assays$RNA@data), 'sparseMatrix')
data <- as(as.matrix(seu@assays$integrated@scale.data), 'sparseMatrix')
pd <- seu@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(expression_data=data,
                         cell_metadata  = seu@meta.data,
                         gene_metadata = fData)

## Step 1: Normalize and pre-process the data
# cds <- preprocess_cds(cds, num_dim = 50, norm_method='log')
cds <- preprocess_cds(cds, num_dim = 100, norm_method='size_only', pseudo_count=0)
## Step 2: Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "orig.ident")
## Step 3: Reduce the dimensions using UMAP
DefaultAssay(seu) <- 'integrated'
reducedDims(cds)$UMAP <- seu@reductions$umap@cell.embeddings
reducedDims(cds)$PCA <- seu@reductions$pca@cell.embeddings[,1:50]
# cds <- reduce_dimension(cds)
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)
## Step 6: Order cells
# cds <- order_cells(cds, reduction_method='UMAP', root_pr_nodes='1')

seu$monocle3_clusters <- cds@clusters$UMAP$clusters
seu$monocle3_partitions <- cds@clusters$UMAP$partitions

pdf("~/xfer/c1i.pdf", width=12); 
DefaultAssay(seu) <- 'SCT'
dp_mc <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="monocle3_partitions", pt.size=0.5, shuffle=TRUE)
dp_mp <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="monocle3_clusters", pt.size=0.5, shuffle=TRUE)
dp_sc <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="seurat_clusters", pt.size=0.5, shuffle=TRUE)
cowplot::plot_grid(dp_mc, dp_sc, ncol=2)
cowplot::plot_grid(dp_mp, dp_sc, ncol=2)
FeaturePlot(seu, features=c("Cd8a", "Marco", "Cd4", "Cd3e", "Foxp3", "Ly6g"), 
            pt.size=0.5, order=T, ncol=3)
DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
        group.by="immgen.fine.cluster", pt.size=0.5, shuffle=TRUE)
plot_cells(cds, color_cells_by='immgen.fine.cluster',show_trajectory_graph=T)
dev.off()


# Prep SCT for DEG testing
saveRDS(seu, file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
saveRDS(cds, file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))


#--- a) TRegs or CD8 only ----
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
celltype <- c('cd8'="^T_CD8")
celltype <- c('tregs'="^TReg")

seul_celltype <- lapply(seul, function(seu){
  Idents(seu) <- 'manual_anno'
  seu <- subset(seu, ident=setdiff(unique(grep(celltype, Idents(seu), value=T, ignore.case=T)),
                                   unique(grep("remove", Idents(seu), value=T, ignore.case = T))))
  return(seu)
})


unionize_untreated <- FALSE
redo_reductions <- FALSE
seu_traj_plots <- lapply(names(seul_celltype), function(seu_id){
  seu_j <- seul_celltype[[seu_id]]
  seu_j$batch <- gsub("_.*", "", seu_j$orig.ident)
  
  traj_plots <- lapply(c("B1", "B2"), function(batchid){
    dims_use <- 1:20
    Idents(seu_j) <- 'batch'
    seu_jbatch <- subset(seu_j, ident=batchid)
    if(redo_reductions){
      seu_jbatch <- seu_jbatch %>%
        RunUMAP(., dims=dims_use, n.neighbors = 20L, reduction='mnn')  %>%
        RunTSNE(., dims=dims_use, reduction='mnn') %>%
        runPacmap(., reduction='mnn')
    }
    
    
    # DefaultAssay(seu_jbatch) <- 'RNA'
    # seu_jbatch <- PrepSCTFindMarkers(seu_jbatch)
    # data <- as(as.matrix(seu_jbatch@assays$RNA@data), 'sparseMatrix')
    data <- as(as.matrix(seu_jbatch@assays$mnn.reconstructed@data), 'sparseMatrix')
    pd <- seu_jbatch@meta.data
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    cds <- new_cell_data_set(expression_data=data,
                             cell_metadata  = seu_jbatch@meta.data,
                             gene_metadata = fData)
    
    # https://github.com/satijalab/seurat/issues/1658
    cds <- preprocess_cds(cds, num_dim = 30, norm_method='size_only', pseudo_count=0)
    DefaultAssay(seu_jbatch) <- 'mnn.reconstructed'
    reducedDims(cds)$UMAP <- seu_jbatch@reductions$umap@cell.embeddings
    reducedDims(cds)$PCA <- seu_jbatch@reductions$mnn@cell.embeddings[,1:30]
    cds <- cluster_cells(cds)
    cds@clusters$UMAP$clusters <- cds$manual_anno
    cds <- learn_graph(cds, use_partition=FALSE, learn_graph_control=list('rann.k'=20))
    
    seu_jbatch$monocle3_clusters <- cds@clusters$UMAP$clusters
    seu_jbatch$monocle3_partitions <- cds@clusters$UMAP$partitions
    dp <- plot_cells(cds, color_cells_by='manual_anno',show_trajectory_graph=T) + 
      ggtitle(paste0(seu_id, "_", batchid))
    return(dp)
  }) %>% cowplot::plot_grid(plotlist=., nrow=1)
  return(traj_plots)
})  %>% cowplot::plot_grid(plotlist=., nrow=2)

pdf("~/xfer/st2_trajectory_plots.2.pdf", width=12); 
seu_traj_plots
dev.off()



#--- a0) Visualize trajectory ----
cds <- readRDS(file=file.path(datadir, "seurat_obj", "3_cds_anno.rds"))
png("~/xfer/trj.png", units="cm", width=7, height=7, res=400)
cds$singlecol <- 1
plot_cells(cds, color_cells_by='seurat_clusters',
           label_groups_by_cluster=FALSE,  label_cell_groups=FALSE,
           label_principal_points = FALSE, 
           label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE,
           show_trajectory_graph=T) +
  NoLegend()
dev.off()
# cds <- order_cells(cds, root_pr_nodes='Y_83')


#--- a) Branch-level DEGs ----
cds <- readRDS(file=file.path(datadir, "seurat_obj", "3_cds_anno.rds"))
# plot_cells(cds, color_cells_by='manual_anno', 
#            label_principal_points = TRUE, show_trajectory_graph=T)
# cds <- order_cells(cds, root_pr_nodes='Y_83')


ncores <- 1 # 4
trajectories <- list("6to4"=c('Y_83', "Y_90"),
                     "6to9"=c('Y_83', 'Y_163'),
                     "6to1"=c('Y_83', 'Y_11'),
                     "6to2"=c('Y_83', 'Y_64'),
                     "6to0"=c('Y_83', 'Y_96'),
                     '17to0'=c('Y_46', 'Y_96'),
                     '17to2'=c('Y_46', 'Y_64'))
groupings <- c('All'='All', 'DAB'='DAB', 'DMSO'='DMSO')
grp_cds_branch_dat <- lapply(groupings, function(grp){
  print(paste0("Groupings: ", grp))
  cds_branch_dat <- lapply(trajectories, function(traj_i){
    print(traj_i)
    # Extract cells along a branch manually as choose_graph_segments() gives error
    # https://github.com/cole-trapnell-lab/monocle3/issues/644
    # cds_i <- choose_graph_segments(cds, 
    #                                starting_pr_node=traj_i[1],
    #                                ending_pr_nodes=traj_i[2],
    #                                return_list=F,
    #                                clear_cds=F)
    dp_mst <- principal_graph(cds)[["UMAP"]]
    path <- igraph::shortest_paths(
      dp_mst, from <- traj_i[1], to = traj_i[2], mode = "all",
      algorithm = "unweighted")
    path.nodes <- path$vpath[[1]]
    
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    cells.branch <- closest_vertex[closest_vertex %in% path.nodes,]
    print(paste0("Number of cells: ", length(cells.branch)))
    if(grp != 'All'){
      cells.branch <- cells.branch[grep(grp, names(cells.branch))]
      print(paste0("Adjusted number of cells: ", length(cells.branch)))
    }
    cds.branch <- cds[,rownames(cds@colData) %in% names(cells.branch)]
    
    # identify genes with interesting patterns of expression that fall 
    # only within the region of the trajectory you selected
    # "test whether cells at similar positions on the trjaectory have correlated expression"
    subset_pr_test_res <- graph_test(cds.branch, 
                                     neighbor_graph="principal_graph", 
                                     cores=ncores)
    pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
    
    # Grouping these genes into modules 
    gene_module_df <- find_gene_modules(cds.branch[pr_deg_ids,], resolution=0.001)
    
    # organize the modules by their similarity (using hclust) over the trajectory 
    agg_mat <- aggregate_gene_expression(cds.branch, gene_module_df)
    module_dendro <- hclust(dist(agg_mat))
    gene_module_df$module <- factor(gene_module_df$module, 
                                    levels = row.names(agg_mat)[module_dendro$order])
    
    return(list("cds_branch"=cds.branch,
                "degs"=subset_pr_test_res,
                "modules"=gene_module_df))
  })
})

saveRDS(grp_cds_branch_dat, file=file.path(outdir, "markers", "grp_cds_branch_dat.rds"))

#--- c) Post-branch analysis ----
cds <- readRDS(file=file.path(datadir, "seurat_obj", "3_cds_anno.rds"))
grp_cds_branch_dat <- readRDS(file=file.path(outdir, "markers", "grp_cds_branch_dat.rds"))
# cds_branch_dat <- grp_cds_branch_dat$All # Legacy code: DAB+DMSO in one dataset
seu_regulon <- readRDS(file.path(outdir, "regulons", 'seu_regulon.rds'))

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets
mlvl <- c('H', 'C2', 'C5')
n = 5
root_node = 'Y_83'
regoi <- c('Mxd1(+)', 'Fosl1(+)', 'Tcf7l2(+)', 
           'Gatad1(+)', 'Atf4(+)', 'Cebpg(+)',
           'Ets2(+)', 'Spi1(+)', 'Nfkb2(+)',
           'Trp53(+)', 'Drap1(+)', 'Mlx(+)',  'Atf5(+)' )
all_regulons <- c(seu_regulon$regulons)
regulons_ds <- lapply(names(all_regulons), function(sublvl){
  msig_ds <- data.frame("gs_name"=sublvl,
                        "entrez_gene"=all_regulons[[sublvl]])
}) %>% do.call(rbind, .)
genesets_regulons <- list("custom"=c(all_regulons))
regulons <- c(seu_regulon$regulons[regoi])
genesets_c <- list("custom"=c(genesets, genesets_public, regulons))

#Iterate through each branch processed
branchid <- 1
grpids <- names(grp_cds_branch_dat)
grp_cds_branch_res <- lapply(setNames(grpids,grpids), function(grpid){
  cds_branch_dat <- grp_cds_branch_dat[[grpid]]
  print(sprintf("Group: %s...", grpid))
  cds_branch_res <- lapply(cds_branch_dat, function(cds_branch_i){
    print(sprintf("Analyzing branch %i...", branchid))
    branchid <<- branchid + 1
    degs_i <- cds_branch_i$degs       # Genes expressed differentially across the clusters
    cds_i <- cds_branch_i$cds_branch  # Monocle3 SCT object
    module_i <- cds_branch_i$modules  # DEGs broken into co-expressed gene modules
    
    # Top significant genes differentiating across pseudotime
    sig_genes <- degs_i %>% arrange(q_value, desc(morans_test_statistic))
    sig_genes <- head(sig_genes$gene_short_name, 9)
    
    # aggregate expression of all genes in each module across all the clusters
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds_i)), 
                                    cell_group=colData(cds_i)$seurat_clusters)
    agg_mat <- aggregate_gene_expression(cds_i, module_i, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    
    # Calculate the average expression per cluster
    seu_branch <- subset(seu, cells=cell_group_df$cell)
    rm_clus <- which(table(seu_branch$seurat_clusters) < 50)
    Idents(seu_branch) <- 'seurat_clusters'
    seu_branch <- subset(seu_branch, ident=setdiff(unique(seu_branch$seurat_clusters), names(rm_clus)))
    expr <- AverageExpression(seu_branch, assay='RNA', slot='data', group.by='seurat_clusters')
    expr <- expr$RNA[,intersect(colnames(agg_mat), colnames(expr$RNA))]
    
    # Calculate the pseudotime for each cell within each cluster used and order
    cds_i <- order_cells(cds_i, reduction_method='UMAP', root_pr_nodes=root_node)
    pseudotime_i <- cds_i@principal_graph_aux$UMAP$pseudotime
    pseudotime_grp_i <- split(pseudotime_i[cell_group_df$cell], 
                              f=cell_group_df$cell_group) %>%
      sapply(., mean) %>%
      sort
    pseudotime_grp_i <- pseudotime_grp_i[intersect(names(pseudotime_grp_i), colnames(expr))]
    
    sce_pseudo <- lapply(names(pseudotime_grp_i), function(i){
      sce_extract(cds_i, assay_name='counts', meta1='seurat_clusters', value_meta1 = i)
    })
    pathways <- "/cluster/projects/mcgahalab/ref/msigdb/aggregate/h_c2_c5.mus_musculus.msigdb.csv"
    sce_scpa <- compare_pathways(samples = sce_pseudo, 
                                 pathways = pathways)
    
    # Perform an over-representation analysis on each Module across MsigDB & Regulons
    top_oras <- lapply(split(module_i, f=module_i$module), function(module_ij){
      #pathway ora
      entrez_genes <- gm$SYMBOL$ENTREZID[module_ij$id]
      ora_ij <- iterateMsigdb(species='Mus musculus', msig_lvls=msig_lvls[mlvl], 
                              fun=oraFun, entrez_genes=entrez_genes) %>% 
        unlist(., recursive=F)
      
      #regulon ora
      ora_ij$regulon <-  enricher(gene = na.omit(names(entrez_genes)), 
                                  TERM2GENE = regulons_ds)@result
      
      ora_ij <- lapply(ora_ij, function(ora_ijk){
        ora_ijk <- ora_ijk %>% 
          select(-c(Description, geneID, Count)) %>%
          arrange(qvalue) %>%
          mutate(top5="NA",
                 Dataset=gsub("_.*", "", ID))
        if(any(grepl("\\(\\+\\)", ora_ijk$Dataset))){
          ora_ijk$Dataset <- 'regulon'
        }
        if(n > nrow(ora_ijk)) n <- nrow(ora_ijk)
        ora_ijk[1:n, 'top'] <- 'Top'
        return(ora_ijk)
      }) %>% 
        do.call(rbind, .) %>%
        as.data.frame
      return(ora_ij)
    })
    top_oras_merge <- lapply(top_oras, function(i){ 
      i %>% 
        filter(!is.na(ID)) %>% 
        select(Dataset, ID, qvalue, top)
    }) %>% 
      purrr::reduce(., full_join, by=c('ID')) %>% 
      tibble::column_to_rownames("ID")
    top_idx <- (rowSums(top_oras_merge=='Top', na.rm=T) > 0) &
      (rowSums(top_oras_merge=='regulon', na.rm=T) == 0)
    top_gs <- rownames(top_oras_merge[which(top_idx),])
    
    # Calculate the gene-set AUC activity for the average expression per cluster
    auc_l <- iterateMsigdb(species='Mus musculus', msig_lvls=msig_lvls[mlvl], 
                           fun=aucellFun, expr_mat=expr, gm=gm,
                           mapfrom='ENTREZID', mapto='SYMBOL') %>% 
      unlist(., recursive=F)
    auc_ij_all <- lapply(auc_l, function(auc_i){
      as.data.frame(assay(auc_i))
    }) %>% 
      do.call(rbind, .) %>%
      as.data.frame
    rownames(auc_ij_all) <- gsub("^.*\\.", "", rownames(auc_ij_all))
    auc_ij <- auc_ij_all[top_gs,]
    
    # Calculate the Regulon GENE-SETS AUC activity for the average expression per cluster
    regulon_auc_l <- iterateMsigdb(species='Mus musculus', msig_lvls=genesets_regulons, 
                                   fun=aucellFun, expr_mat=expr, gm=gm, 
                                   mapfrom='SYMBOL', mapto='SYMBOL') %>% 
      unlist(., recursive=F)
    regulon_auc_ij_all <- lapply(regulon_auc_l, function(auc_i){
      as.data.frame(assay(auc_i))
    }) %>% 
      do.call(rbind, .) %>%
      as.data.frame
    rownames(regulon_auc_ij_all) <- gsub("^.*\\.", "", rownames(regulon_auc_ij_all))
    (rowSums(top_oras_merge=='Top', na.rm=T) > 0) &
      (rowSums(top_oras_merge=='regulon', na.rm=T) == 0)
    top_regulons  <- top_oras_merge %>% 
      filter((rowSums(.=='regulon', na.rm=T) > 0)) %>% 
      select(grep("qvalue", colnames(.), value=T)) %>%
      filter(rowSums(. < 0.01) > 0) %>% 
      rownames(.)
    regulon_auc_ij <- regulon_auc_ij_all[top_regulons,]
    
    # Calculate the CUSTOM GENE-SETS AUC activity for the average expression per cluster
    custom_auc_l <- iterateMsigdb(species='Mus musculus', msig_lvls=genesets_c, 
                                  fun=aucellFun, expr_mat=expr, gm=gm, 
                                  mapfrom='SYMBOL', mapto='SYMBOL') %>% 
      unlist(., recursive=F)
    custom_auc_ij_all <- lapply(custom_auc_l, function(auc_i){
      as.data.frame(assay(auc_i))
    }) %>% 
      do.call(rbind, .) %>%
      as.data.frame
    rownames(custom_auc_ij_all) <- gsub("^.*\\.", "", rownames(custom_auc_ij_all))
    
    return(list("auc_ij"=auc_ij,                         # Top genesets AUCell score for cluster average expression
                "auc_ij_all"=auc_ij_all,                 # All genesets AUCell score for cluster average expression
                "custom_auc_ij_all"=custom_auc_ij_all,   # Custom geneset AUCell score
                "regulon_auc_ij_all"=regulon_auc_ij_all, # Similar to custom_auc_ij_all, but for regulons instead of pathways
                "regulon_auc_ij"=regulon_auc_ij,         # Top regulon AUCell scores for cluster average expression
                "pseudotime_grp_i"=pseudotime_grp_i,     # Average pseudotime per cluster
                "top_oras"=top_oras_merge,               # ORA of module-level genes within each pathway and regulon
                "expr"=expr,                             # RNA count matrix: gene by cell
                "cds_i"=cds_i,                           # Monocle3 SCT object
                "module_i"=module_i,                     # DEG split into co-expressed modules
                "agg_mat"=agg_mat,                       # Expression of modules per cluster
                "sig_genes"=sig_genes,                   # Significant DEG differentiating branch-level pseudotime
                "degs_i"=degs_i))                        # DEGS differentiating branch-level pseudotime
  })
})

saveRDS(grp_cds_branch_res, file=file.path(outdir, "markers", "grp_cds_branch_res.rds"))

grp_cds_branch_res <- readRDS(file=file.path(outdir, "markers", "grp_cds_branch_res.rds"))

geneset_patterns <- list(c("keep", "HALLMARK|REACTOME|KEGG"),
                         c("keep", "GOMF|GOCC|GOBP"),
                         c("remove", "HALLMARK|REACTOME|KEGG|GOMF|GOCC|GOBP"),
                         c("keep", "\\(\\+\\)|\\(-\\)"))# "*" # "HALLMARK|REACTOME|KEGG" "GOMF|GOCC|GOBP"
for(grp_id in names(grp_cds_branch_res)){
  cds_branch_res <- grp_cds_branch_res[[grp_id]]
  for(bri_id in names(cds_branch_res)){
    bri <- cds_branch_res[[bri_id]]
    pdf(paste0("~/xfer/traj.",grp_id, "_", bri_id, ".pdf"), width = 13, height = 10)
    
    
    # Pseudotime of cells along the principal path for the selected branch
    gp <- plot_cells(bri$cds_i,
                     color_cells_by = "pseudotime",
                     label_cell_groups=T,
                     label_leaves=T,
                     label_branch_points=T,
                     graph_label_size=1.5) +
      ggtitle(bri_id)
    dp <- Cluster_Highlight_Plot(seurat_object = seu, reduction='umap',
                                 cluster_name = names(bri$pseudotime_grp_i),
                                 raster=T, label=T, highlight_color='cyan')
    plot(gp + dp)
    
    # Top signficant DEGs alogn the principal path
    plot_cells(bri$cds_i, genes=bri$sig_genes,
               reduction_method='UMAP',
               show_trajectory_graph=TRUE,
               label_cell_groups=FALSE,
               label_leaves=FALSE) %>% plot
    
    # DEG Module louvain clustering
    ggp <- ggplot(as.data.frame(bri$module_i), aes(x=dim_1, y=dim_2, colour=module))+
      geom_point() +
      theme_bw()
    plot(ggp)
    
    # DEG Module activity plotted on UMAP along principal path
    plot_cells(bri$cds_i, 
               genes=bri$module_i,
               group_cells_by="partition",
               color_cells_by="partition",
               show_trajectory_graph=T) %>% plot
    
    # Heatmap of DEG Module activity across all clusters
    # pheatmap::pheatmap(bri$agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
    #                    scale="column", clustering_method="ward.D2",
    #                    fontsize=6)
    pheatmap::pheatmap(bri$agg_mat[,names(bri$pseudotime_grp_i)], cluster_rows=TRUE, 
                       cluster_cols=FALSE, scale="column", 
                       clustering_method="ward.D2", fontsize=6)
    
    # Heatmap of gene-sets per cluster ordered by pseudotime
    for(geneset_pattern in geneset_patterns){
      auc_ij_comb <- rbind(bri$custom_auc_ij_all, bri$auc_ij)
      
      if(geneset_pattern[1]=='keep'){
        custom_idx <- grepl(geneset_pattern[2], rownames(auc_ij_comb))
      } else {
        custom_idx <- !grepl(geneset_pattern[2], rownames(auc_ij_comb))
      }
      auc_ij_comb <- auc_ij_comb[custom_idx,]
      
      auc_pseudotime <- auc_ij_comb[,names(bri$pseudotime_grp_i)]
      auc_pseudotime[is.na(auc_pseudotime)] <- 0
      zero_idx <- (rowSums(auc_pseudotime) == 0)
      if(any(zero_idx)) auc_pseudotime <- auc_pseudotime[-which(zero_idx),]
      pheatmap::pheatmap(auc_pseudotime, cluster_rows = T, cluster_cols = F, scale='row')
    }
    
    dev.off()
  }
}

for(grp_id in names(grp_cds_branch_res)){
  for(bri_id in names(cds_branch_res)){
    cat(paste0("\nxfer traj.", grp_id, "_", bri_id, ".pdf\n "))
  }
}

#--- d) Add on to post-branch analysis ----
cds <- readRDS(file=file.path(datadir, "seurat_obj", "3_cds_anno.rds"))
cds_branch_dat <- readRDS(file=file.path(outdir, "markers", "cds_branch_dat.rds"))
outdirregulon <- file.path(outdir, "regulons")
seu_regulon <- readRDS(file = file.path(outdirregulon, 'seu_regulon.rds'))
cds_branch_res <- readRDS(file=file.path(outdir, "markers", "cds_branch_markers.rds"))


msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets
mlvl <- c('H', 'C2', 'C5')
n = 5
root_node = 'Y_83'
regoi <- c('Mxd1(+)', 'Fosl1(+)', 'Tcf7l2(+)', 'Drap1(+)', 'Mlx(+)',  'Atf5(+)' )
regulons <- c(seu_regulon$regulons[regoi])
genesets_c <- list("custom"=c(genesets, genesets_public, regulons))

#Iterate through each branch processed
cds_branch_res_custom <- lapply(cds_branch_dat, function(cds_branch_i){
  degs_i <- cds_branch_i$degs       # Genes expressed differentially across the clusters
  cds_i <- cds_branch_i$cds_branch  # Monocle3 SCT object
  module_i <- cds_branch_i$modules  # DEGs broken into co-expressed gene modules
  
  # aggregate expression of all genes in each module across all the clusters
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds_i)), 
                                  cell_group=colData(cds_i)$seurat_clusters)
  agg_mat <- aggregate_gene_expression(cds_i, module_i, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  
  # Calculate the average expression per cluster
  seu_branch <- subset(seu, cells=cell_group_df$cell)
  rm_clus <- which(table(seu_branch$seurat_clusters) < 50)
  Idents(seu_branch) <- 'seurat_clusters'
  seu_branch <- subset(seu_branch, ident=setdiff(unique(seu_branch$seurat_clusters), names(rm_clus)))
  expr <- AverageExpression(seu_branch, assay='RNA', slot='data', group.by='seurat_clusters')
  expr <- expr$RNA[,intersect(colnames(agg_mat), colnames(expr$RNA))]
  
  # Calculate the CUSTOM GENE-SETS AUC activity for the average expression per cluster
  custom_auc_l <- iterateMsigdb(species='Mus musculus', msig_lvls=genesets_c, 
                                fun=aucellFun, expr_mat=expr, gm=gm, mapfrom='SYMBOL', mapto='SYMBOL') %>% 
    unlist(., recursive=F)
  custom_auc_ij_all <- lapply(custom_auc_l, function(auc_i){
    as.data.frame(assay(auc_i))
  }) %>% 
    do.call(rbind, .) %>%
    as.data.frame
  rownames(custom_auc_ij_all) <- gsub("^.*\\.", "", rownames(custom_auc_ij_all))
  
  return(list("custom_auc_ij_all"=custom_auc_ij_all))
})

cds_branch_res[[1]]$custom_auc_ij_all <- cds_branch_res_custom[[1]]$custom_auc_ij_all
cds_branch_res[[2]]$custom_auc_ij_all <- cds_branch_res_custom[[2]]$custom_auc_ij_all
cds_branch_res[[3]]$custom_auc_ij_all <- cds_branch_res_custom[[3]]$custom_auc_ij_all
cds_branch_res[[4]]$custom_auc_ij_all <- cds_branch_res_custom[[4]]$custom_auc_ij_all
cds_branch_res[[5]]$custom_auc_ij_all <- cds_branch_res_custom[[5]]$custom_auc_ij_all
# saveRDS(cds_branch_res, file=file.path(outdir, "markers", "cds_branch_markers.rds"))



######################################
#### XXX. Testing ST2 sample swap ####
seulst2 <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
seulcd45 <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
c(sapply(seulst2, function(i) unique(i$manual_anno)),
  sapply(seulcd45, function(i) unique(i$manual_anno))) %>%
  unlist %>% as.character %>% table

seul <- lapply(c(seulst2, seulcd45), function(seu_i){
  Idents(seu_i) <- 'orig.ident'
  seu_i <- subset(seu_i, ident=unique(grep("^B1", Idents(seu_i), ignore.case = T, value=T)))
  Idents(seu_i) <- 'manual_anno'
  subset(seu_i, ident=unique(grep("Treg.*NLT", Idents(seu_i), ignore.case = T, value=T)))
})
# rm(seulst2); rm(seulcd45); gc()

seu <- merge(seul[[1]], seul[-1])
seu <- seul$LN

options(Seurat.object.assay.version = "v5")
seu <- UpdateSeuratObject(object = seu)

seu <- tryCatch({
  JoinLayers(seu)
}, error=function(e){seu})

find.var.feats.kowt <- FALSE
if(find.var.feats.kowt){
  seu$kowt <- gsub("^.*_(KO|WT)_.*", "\\1", seu$orig.ident)
  seu$id <- gsub("_(KO|WT)_", "_", seu$orig.ident)
  
  kowt_markers <- SplitObject(seu, split.by='id') %>%
    lapply(., function(seu_i){
      FindMarkers(seu_i, 
                  ident.1='KO', ident.2='WT',
                  group.by='kowt')
    })
  kowt_markers_v <- lapply(kowt_markers, function(marker_i){
    marker_i %>% 
      filter(p_val_adj < 0.05) %>% 
      rownames(.)
  }) %>% unlist %>%
    table %>% sort
  
} else {
  kowt_markers_v <- c('Gata3','Foxp3','Itgae','Ahr','Rara','Tgfb1','Nfkbib',
                      'Il1rl1','MyD88','Ccr1','Kit','Adam8','Ikzf2','Ccr4',
                      'Ccr5','Ccr9','Cxcr3','Itgav','Pdcd1')
  kowt_markers_v <- setNames(rep(2, length(kowt_markers_v)), kowt_markers_v)
}


DefaultAssay(seu) <- 'RNA'
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

VariableFeatures(seu) <- names(which(kowt_markers_v > 1))
ndim <- 15
seu <- NormalizeData(seu) %>%
  # FindVariableFeatures(.)  %>%
  ScaleData(.)  %>%
  RunPCA(.) %>% 
  FindNeighbors(., dims = 1:ndim, reduction = "pca")  %>%
  FindClusters(., resolution = 0.7, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(., dims = 1:ndim, reduction = "pca", reduction.name = "umap.unintegrated")

# seu <- SCTransform(seu, vst.flavor = "v2") %>%
#   RunPCA(npcs = 30, verbose = FALSE) %>%
#   RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
#   FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
#   FindClusters(resolution = 0.7, verbose = FALSE)

DefaultAssay(seu) <- 'RNA'
seu_integ <- IntegrateLayers(
  object = seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = T,
  group.by='orig.ident',
  features=VariableFeatures(seu)
)
seu_integ <- seu_integ %>% 
  FindNeighbors(., reduction = "integrated.mnn", dims = 1:ndim) %>%
  FindClusters(., resolution = 0.6, cluster.name = "mnn_clusters") %>% 
  RunUMAP(., reduction = "integrated.mnn", dims = 1:ndim, reduction.name = "umap.mnn")

seu_integ <- JoinLayers(seu_integ)

seu_integ@meta.data$ST2wt_ko <- gsub(".*(KO|WT).*", "\\1", seu_integ$orig.ident)
seu_integ@meta.data$CD45_ST2 <- ifelse(grepl("CD45", seu_integ$orig.ident), "CD45", "ST2")
seu_integ@meta.data$Day <- gsub(".*(3d|7d|Un)", "\\1", seu_integ$orig.ident)
seu_integ.bkup <- seu_integ
seu_integ@meta.data <- seu_integ@meta.data[,c('ST2wt_ko', 'CD45_ST2', 'Day', 'orig.ident')]
RNA <- seu_integ@assays$RNA$counts
seu_integ[['RNA2']] <- CreateAssayObject(counts=RNA[which(rowSums(RNA) > 100),])
DefaultAssay(seu_integ) <- 'RNA2'

for(pc_i in c(1:ndim)){
  pc_j <- pc_i + 1
  pckey <- paste0("PC", paste(c(pc_i, pc_j), collapse="."))
  seu_integ[[pckey]] <- CreateDimReducObject(
    embeddings = Embeddings(seu_integ)[,paste0('PC_', c(pc_i, pc_j))], 
    key = 'PCA_', assay = 'RNA'
  )
}


dir.create(file.path(PDIR, "results", "cloupe"), showWarnings = F)
loupeR::create_loupe_from_seurat(seu_integ, output_dir=file.path(PDIR, "results", "cloupe"),
                         output_name="tregs_all_samples.selgenes", force=T)
file.copy(file.path(PDIR, "results", "cloupe", "tregs_all_samples.selgenes.cloupe"), to = "~/xfer", overwrite = T)
cat(paste0("xfer tregs_all_samples.selgenes.cloupe\n"))



                                                      
##############################################
#### 4.a Differential expression analysis ####
# if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.mnn.rds"))
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))

if(rm_tmp_aug2023){
  ln <- c("Tox", "Runx1", "Ikzf2", "Ikzf1", "Rora", "Akt3", "Ncor1", "Elf1", "Prkcq", "Nfatc3", "Fyn", "Fli1", "Foxo1", "Sik3", "Tank", "Mif", "Btf3", "Hopx", "Traf3", "Itgae", "Tnik", "Nfat5", "Nt5e", "Nfkb1")
  tumor <- c("Ms4a4b", "Ms4a4c", "Cxcl10", "Dgat1", "Areg", "Isg15", "Gzmb", "Nfkbia", "Gadd45b", "Cd274", "Nfkbiz", "Tnfaip3", "Cd69")
  day3deg <- c("Tox", "Tnfrsf18", "Mki67", "Maf", "Lgals1", "Klrg1", "Itgae", "Il4i1", "Il2rb", "Il2ra", "Ikzf4", "Ikzf2", "Icos", "Gzmb", "Gata3", "Foxp3", "Cxcr3", "Ctla4", "Cd69", "Ccr7", "Ccr2")
  pdf("~/xfer/dotplots.pdf")
  Idents(seul$LN) <- 'manual_anno'
  Idents(seul$Tumor) <- 'manual_anno'
  DotPlot_scCustom(seul$LN, features=ln, flip_axes=T, x_lab_rotate=T) + ggtitle("LN")
  DotPlot_scCustom(seul$Tumor, features=tumor, flip_axes=T, x_lab_rotate=T) + ggtitle("Tumor")
  dev.off()
  
  pdf("~/xfer/vlnplots.pdf", width = 20, height = 10)
  cells <- Cells(seul$LN)[which(grepl("B1.*(3d|Un)$", seul$LN$orig.ident) &
                            grepl("^TReg", seul$LN$manual_anno))]
  seuln <- seul$LN %>%
    subset(., cells=cells) %>%
    ScaleData(., day3deg)
  
  Idents(seuln) <- 'manual_anno'
  pdf("~/xfer/vlnplots.pdf", width = 15, height = 25)
  tryCatch({print(Stacked_VlnPlot(seuln, features=day3deg, 
                   split.by='orig.ident', raster=T, x_lab_rotate=T,
                   plot_legend=T))},
           error=function(e){NULL})
  dev.off()
  
  pdf("~/xfer/cd45_cd8_vlnplots.pdf")
  lapply(seul, function(seu_i){
    idents <- grep("cd8", seu_i$manual_anno, ignore.case=T, value=T)
    Idents(seu_i) <- 'manual_anno'
    
    seu_subi <- subset(seu_i, idents=idents) %>%
      ScaleData(., c('Il7r', 'Klrg1'))
    Idents(seu_subi) <- 'manual_anno'
    print(Stacked_VlnPlot(seu_subi, features=c('Il7r', 'Klrg1'), 
                          split.by='orig.ident', raster=T, x_lab_rotate=T,
                          plot_legend=T))
  })
  dev.off()
  
  
  lapply(list(seulcd45, seulst2), function(seul_i){
    seu_i <- seul_i$LN
    cells <- Cells(seu_i)[grep("WT_7d", seu_i$orig.ident)]
    seu_isub <- subset(seu_i, cells=cells)
    cnts <- GetAssayData(seu_isub, slot='counts')
    seu_isub$st2pos <- cnts['Il1rl1',] > 0
    t(table(seu_isub$st2pos, seu_isub$manual_anno))
    # t(apply(mat, 1, function(i) round(i/sum(i),2)))
  })
}


outdirdeg <- file.path(outdir, "degs")
dir.create(outdirdeg, recursive = F)
anno_ident <- 'seurat_clusters' #'manual_anno', 'seurat_clusters'
cols=c("KO.down_WT.down" = "#fdb462", 
       "KO.down_WT.up" = "#80b1d3", 
       "KO.up_WT.up" = "#fdb462", 
       "KO.up_WT.down" = "#fb8072")

anno_ident <- 'manual_anno' #'manual_anno', 'manual_clusters'
pval_sig <- 0.05
fc_sig <- 0.5
pseudop <- 1*10^-250
  
#--- a) ST2: Differential expression per cluster (seurat_clusters) or celltype (manual_anno) ----
# Need to compare treatment-vs-untreated for
#  -  [LN/Tumor]  |  [WT/KO]  |  [Day3/Day7] vs WT_Untreated
#  -  [LN/Tumor]  |  [Day3/Day7]  |  KO_Treated vs WT_Treated
#  -  [LN/Tumor]  |  [Day3/Day7]  |  KO_Untreated vs WT_Untreated
comparisons <- list('LN'=list(
    'WT_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='WT', "order"='7d'),
    'WT_day3_vs_un'=list("day"=c('3d', 'Un'), "batch"='B1', "kowt"='WT', "order"='3d'),
    'KO_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='KO', "order"='7d'),
    'KO_day3_vs_un'=list("day"=c('3d', 'Un'), "batch"='B1', "kowt"='KO', "order"='3d'),
    'DAY3_ko_vs_wt'=list("day"='3d', "batch"='B1', "kowt"=c('KO', 'WT'), "order"='KO'),
    'DAY7_ko_vs_wt'=list("day"='7d', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO'),
    'B1UN_ko_vs_wt'=list("day"='Un', "batch"='B1', "kowt"=c('KO', 'WT'), "order"='KO'),
    'B2UN_ko_vs_wt'=list("day"='Un', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO')),
  'Tumor'=list(
    'WT_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='WT', "order"='7d'),
    'WT_day3_vs_un'=list("day"=c('3d', 'Un'), "batch"='B1', "kowt"='WT', "order"='3d'),
    'KO_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='KO', "order"='7d'),
    'KO_day3_vs_un'=list("day"=c('3d', 'Un'), "batch"='B1', "kowt"='KO', "order"='3d'),
    'DAY3_ko_vs_wt'=list("day"='3d', "batch"='B1', "kowt"=c('KO', 'WT'), "order"='KO'),
    'DAY7_ko_vs_wt'=list("day"='7d', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO'),
    'B1UN_ko_vs_wt'=list("day"='Un', "batch"='B1', "kowt"=c('KO', 'WT'), "order"='KO'),
    'B2UN_ko_vs_wt'=list("day"='Un', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO')))

compare_grp <- list("LN_day3_7"=c("day3_kowt", "day7_kowt", "_.*$"),
                    "LN_day3_un"=c("day3_kowt", "unb1_kowt", "_.*$"),
                    "LN_day7_un"=c("day7_kowt", "unb2_kowt", "_.*$"),
                    "Tumor_day3_7"=c("day3_kowt", "day7_kowt", "_.*$"),
                    "Tumor_day3_un"=c("day3_kowt", "unb1_kowt", "_.*$"),
                    "Tumor_day7_un"=c("day7_kowt", "unb2_kowt", "_.*$"))
for(seuid in c('LN', 'Tumor')){
  grep("Treg.*cycling", seul[[seuid]]$manual_anno, value=T, ignore.case = T) %>% table
  seul[[seuid]]$custom_comp <- ifelse(grepl("Treg.*cycling", seul[[seuid]]$manual_anno, ignore.case = T), 'TReg_Cycling', 'Other')
  oidx <- which(seul[[seuid]]$custom_comp == 'Other')
  grep("Treg.*(Central|Effector)", seul[[seuid]]$manual_anno[oidx], ignore.case = T, value=T) %>% table
  seul[[seuid]]$custom_comp[oidx] <- ifelse(grepl("Treg.*(Central|Effector)", seul[[seuid]]$manual_anno[oidx], ignore.case = T), 'TReg_Effector.Central', 'Other')
}
anno_ident <- 'manual_anno' #'manual_anno', 'manual_clusters'

overwrite <- FALSE
degf <- file.path(outdirdeg, paste0(anno_ident, "_degs.rds")) #"_deg2s.rds"))
degf <- file.path(outdirdeg, paste0(anno_ident, "_degs_tregCustom.rds")) #"_deg2s.rds"))
if(any(!file.exists(degf) | overwrite)){
  print(paste0("Creating new DEG file for: ", anno_ident))
  options(future.globals.maxSize= 891289600)
  ct_markers_all <- lapply(names(seul), function(seu_id){
    print(paste0("> ", seu_id))
    seu_i <- seul[[seu_id]]
    seu_i$day <- gsub("^.*_", "", seu_i$orig.ident)
    seu_i$kowt <- gsub("^.*_(.*)_[a-zA-Z0-9]*$", "\\1", seu_i$orig.ident)
    seu_i$batch <- gsub("_.*", "", seu_i$orig.ident)
    
    Idents(seu_i) <- 'orig.ident'
    ct_markers <- lapply(names(comparisons[[seu_id]]), function(comp_id){
      comp_i <- comparisons[[seu_id]][[comp_id]]
      print(paste0(">> ", comp_id))
      print(paste0("  -- ", paste(unlist(comp_i), collapse="-")))
      
      subset_of_seu <- TRUE
      idents <- with(seu_i@meta.data, orig.ident[(day %in% comp_i$day) & 
                                                   (kowt %in% comp_i$kowt) & 
                                                   (batch == comp_i$batch)]) %>%
        unique
      idents <- idents[order(grepl(comp_i$order, idents))]
      
      # DEG across all cells (not cluster specific)
      markers_all <- FindMarkers(seu_i, assay = "RNA", test.use='wilcox',
                                 ident.1= idents[1], ident.2= idents[2],
                                 verbose = FALSE,
                                 logfc.threshold = 0,
                                 recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
        tibble::rownames_to_column(., "symbol") %>% 
        mutate(biotype = sym2biotype_ids[symbol],
               ensemble = gm$SYMBOL$ENSEMBL[symbol],
               ident.1 = idents[1],
               ident.2 = idents[2], 
               cluster='all_cells',
               comparison=comp_id,
               tissue=seu_id)
      
      # DEG For each cluster
      Idents(seu_i) <- anno_ident
      cluster_ids <- as.character(unique(Idents(seu_i)))
      markers_clusters <- lapply(setNames(cluster_ids,cluster_ids), function(cl_j){
        print(paste0("  -- ", cl_j))
        seu_j <- subset(seu_i, ident=cl_j)
        Idents(seu_j) <- 'orig.ident'
        tryCatch({
          FindMarkers(seu_j, assay = "RNA", test.use='wilcox',
                      ident.1= idents[1], ident.2= idents[2],
                      verbose = FALSE,
                      logfc.threshold = 0,
                      recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
            tibble::rownames_to_column(., "symbol") %>% 
            mutate(biotype = sym2biotype_ids[symbol],
                   ensemble = gm$SYMBOL$ENSEMBL[symbol],
                   ident.1 = idents[1],
                   ident.2 = idents[2], 
                   cluster=cl_j,
                   comparison=comp_id,
                   tissue=seu_id)
        }, error=function(e){
          print("     - No differential genes")
          NULL
        })
      })
      
      # Aggregate and return
      markers <- c(list("All"=markers_all),
                   markers_clusters)
      return(markers)
    })
    
    return(ct_markers)
  })
  names(ct_markers_all) <- names(seul)
  names(ct_markers_all$LN) <- names(comparisons$LN)
  names(ct_markers_all$Tumor) <- names(comparisons$Tumor)
  saveRDS(ct_markers_all, file=degf)
} else {
  print(paste0("Reading in existing DEG file for: ", anno_ident))
  ct_markers_all <- readRDS(degf)
}

# Find all cluster names (e.g. 1,2,3 or Tregs, B, NKT)
all_cts <- sapply(unlist(ct_markers_all, recursive = F), names) %>%
  unlist %>% as.character %>% unique
# all_cts <- all_cts[grep("TReg", all_cts)]

## Write out the DEG files
# ct_markers[[LN/Tumor]][[names(comparisons)]][[celltype]]

lapply(names(ct_markers_all), function(seu_id){
  deg <- ct_markers_all[[seu_id]]
  lapply(names(deg), function(comp_id){
    deg_comp_i <- deg[[comp_id]]
    ct_ids <- names(deg_comp_i)
    merge_ids <- c('symbol', 'ensemble', 'biotype', 'ident.1', 'ident.2', 'comparison', 'tissue')
    deg_merge <- lapply(names(deg_comp_i), function(deg_id){
      if(!is.null(deg_comp_i[[deg_id]])){
        deg_comp_i[[deg_id]] %>%
          select(merge_ids, avg_log2FC, p_val, p_val_adj) %>%
          rename_with(., ~paste0(deg_id, ".", .), .cols=c(avg_log2FC, p_val, p_val_adj))
      } else {
        as.data.frame(matrix(rep("NA", length(merge_ids)), nrow=1)) %>%
          magrittr::set_colnames(merge_ids)
      }
    }) %>%
      purrr::reduce(., left_join, by=merge_ids)
    write.table(deg_merge, 
                file=file.path(outdirdeg, paste0("deg.", seu_id, ".", comp_id, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  })
})


#--- b) CD45: Differential expression per cluster (seurat_clusters) or celltype (manual_anno) ----
# Need to compare treatment-vs-untreated for
#  -  [LN/Tumor]  |  [WT/KO]  |  [Day3/Day7] vs WT_Untreated
#  -  [LN/Tumor]  |  [Day3/Day7]  |  KO_Treated vs WT_Treated
#  -  [LN/Tumor]  |  [Day3/Day7]  |  KO_Untreated vs WT_Untreated
comparisons <- list('LN'=list(
  'WT_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='CD45', "kowt"='WT', "order"='7d'),
  'KO_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='CD45', "kowt"='KO', "order"='7d'),
  'DAY7_ko_vs_wt'=list("day"='7d', "batch"='CD45', "kowt"=c('KO', 'WT'), "order"='KO'),
  'UN_ko_vs_wt'=list("day"='Un', "batch"='CD45', "kowt"=c('KO', 'WT'), "order"='KO')),
  'Tumor'=list(
    'WT_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='CD45', "kowt"='WT', "order"='7d'),
    'KO_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='CD45', "kowt"='KO', "order"='7d'),
    'DAY7_ko_vs_wt'=list("day"='7d', "batch"='CD45', "kowt"=c('KO', 'WT'), "order"='KO'),
    'UN_ko_vs_wt'=list("day"='Un', "batch"='CD45', "kowt"=c('KO', 'WT'), "order"='KO')))

compare_grp <- list("LN_day3_7"=c("day3_kowt", "day7_kowt", "_.*$"),
                    "LN_day3_un"=c("day3_kowt", "unb1_kowt", "_.*$"),
                    "LN_day7_un"=c("day7_kowt", "unb2_kowt", "_.*$"),
                    "Tumor_day3_7"=c("day3_kowt", "day7_kowt", "_.*$"),
                    "Tumor_day3_un"=c("day3_kowt", "unb1_kowt", "_.*$"),
                    "Tumor_day7_un"=c("day7_kowt", "unb2_kowt", "_.*$"))
for(seuid in c('LN', 'Tumor')){
  # grep("Treg.*cycling", seul[[seuid]]$manual_anno, value=T, ignore.case = T) %>% table
  seul[[seuid]]$custom_comp <- ifelse(grepl("Treg.*cycling", seul[[seuid]]$manual_anno, ignore.case = T), 'TReg_Cycling', 'Other')
  oidx <- which(seul[[seuid]]$custom_comp == 'Other')
  # grep("Treg.*", seul[[seuid]]$manual_anno[oidx], ignore.case = T, value=T) %>% table
  seul[[seuid]]$custom_comp[oidx] <- ifelse(grepl("Treg.*", seul[[seuid]]$manual_anno[oidx], ignore.case = T), 'TReg', 'Other')
}
anno_ident <- 'manual_anno' #'manual_anno', 'manual_clusters'

overwrite <- FALSE
degf <- file.path(outdirdeg, paste0("CD45_", anno_ident, "_degs.rds")) #"_deg2s.rds"))
degf <- file.path(outdirdeg, paste0("CD45_", anno_ident, "_degs_tregCustom.rds")) #"_deg2s.rds"))
if(any(!file.exists(degf) | overwrite)){
  print(paste0("Creating new DEG file for: ", anno_ident))
  options(future.globals.maxSize= 891289600)
  ct_markers_all <- lapply(names(seul), function(seu_id){
    print(paste0("> ", seu_id))
    seu_i <- seul[[seu_id]]
    seu_i$day <- gsub("^.*_", "", seu_i$orig.ident)
    seu_i$kowt <- gsub("^.*_(.*)_[a-zA-Z0-9]*$", "\\1", seu_i$orig.ident)
    seu_i$batch <- gsub("_.*", "", seu_i$orig.ident)
    
    Idents(seu_i) <- 'orig.ident'
    ct_markers <- lapply(names(comparisons[[seu_id]]), function(comp_id){
      comp_i <- comparisons[[seu_id]][[comp_id]]
      print(paste0(">> ", comp_id))
      print(paste0("  -- ", paste(unlist(comp_i), collapse="-")))
      
      subset_of_seu <- TRUE
      idents <- with(seu_i@meta.data, orig.ident[(day %in% comp_i$day) & 
                                                   (kowt %in% comp_i$kowt) & 
                                                   (batch == comp_i$batch)]) %>%
        unique
      idents <- idents[order(grepl(comp_i$order, idents))]
      
      # idents <- c('CD8_memory2', 'CD8_memory')
      # Idents(seu_i) <- 'manual_anno'
      # # DEG across all cells (not cluster specific)
      # markers_all <- FindMarkers(seu_i, assay = "RNA", test.use='wilcox',
      #                            ident.1= idents[2],
      #                            verbose = FALSE,
      #                            logfc.threshold = 0,
      #                            recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>%
      #   tibble::rownames_to_column(., "symbol") %>%
      #   mutate(biotype = sym2biotype_ids[symbol],
      #          ensemble = gm$SYMBOL$ENSEMBL[symbol],
      #          ident.1 = idents[1],
      #          cluster='all_cells',
      #          comparison=comp_id,
      #          tissue=seu_id)
      # markers_l <- list()
      # markers_l[[idents[2]]] <- markers_all
      
      # DEG across all cells (not cluster specific)
      markers_all <- FindMarkers(seu_i, assay = "RNA", test.use='wilcox',
                                 ident.1= idents[1], ident.2= idents[2],
                                 verbose = FALSE,
                                 logfc.threshold = 0,
                                 recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
        tibble::rownames_to_column(., "symbol") %>% 
        mutate(biotype = sym2biotype_ids[symbol],
               ensemble = gm$SYMBOL$ENSEMBL[symbol],
               ident.1 = idents[1],
               ident.2 = idents[2], 
               cluster='all_cells',
               comparison=comp_id,
               tissue=seu_id)
      
      # DEG For each cluster
      Idents(seu_i) <- anno_ident
      cluster_ids <- as.character(unique(Idents(seu_i)))
      seu_i <- NormalizeData(seu_i) %>%
        FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
      markers_clusters <- lapply(setNames(cluster_ids,cluster_ids), function(cl_j){
        print(paste0("  -- ", cl_j))
        seu_j <- subset(seu_i, ident=cl_j)
        Idents(seu_j) <- 'orig.ident'
        tryCatch({
          x <- FindMarkers(seu_j, assay = "RNA", slot='data', test.use='wilcox',
                      ident.1= idents[1], ident.2= idents[2],
                      verbose = FALSE,
                      logfc.threshold = 0,
                      recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
            tibble::rownames_to_column(., "symbol") %>% 
            mutate(biotype = sym2biotype_ids[symbol],
                   ensemble = gm$SYMBOL$ENSEMBL[symbol],
                   ident.1 = idents[1],
                   ident.2 = idents[2], 
                   cluster=cl_j,
                   comparison=comp_id,
                   tissue=seu_id)
        }, error=function(e){
          print("     - No differential genes")
          NULL
        })
      })
      
      # Aggregate and return
      markers <- c(list("All"=markers_all),
                   markers_clusters)
      return(markers)
    })
    
    return(ct_markers)
  })
  names(ct_markers_all) <- names(seul)
  names(ct_markers_all$LN) <- names(comparisons$LN)
  names(ct_markers_all$Tumor) <- names(comparisons$Tumor)
  saveRDS(ct_markers_all, file=degf)
} else {
  print(paste0("Reading in existing DEG file for: ", anno_ident))
  ct_markers_all <- readRDS(degf)
}

# Find all cluster names (e.g. 1,2,3 or Tregs, B, NKT)
all_cts <- sapply(unlist(ct_markers_all, recursive = F), names) %>%
  unlist %>% as.character %>% unique
# all_cts <- all_cts[grep("TReg", all_cts)]

## Write out the DEG files
# ct_markers[[LN/Tumor]][[names(comparisons)]][[celltype]]

lapply(names(ct_markers_all), function(seu_id){
  deg <- ct_markers_all[[seu_id]]
  lapply(names(deg), function(comp_id){
    deg_comp_i <- deg[[comp_id]]
    ct_ids <- names(deg_comp_i)
    merge_ids <- c('symbol', 'ensemble', 'biotype', 'ident.1', 'ident.2', 'comparison', 'tissue')
    deg_merge <- lapply(names(deg_comp_i), function(deg_id){
      if(!is.null(deg_comp_i[[deg_id]])){
        deg_comp_i[[deg_id]] %>%
          select(merge_ids, avg_log2FC, p_val, p_val_adj) %>%
          rename_with(., ~paste0(deg_id, ".", .), .cols=c(avg_log2FC, p_val, p_val_adj))
      } else {
        as.data.frame(matrix(rep("NA", length(merge_ids)), nrow=1)) %>%
          magrittr::set_colnames(merge_ids)
      }
    }) %>%
      purrr::reduce(., left_join, by=merge_ids)
    write.table(deg_merge, 
                file=file.path(outdirdeg, paste0("CD45_deg.", seu_id, ".", comp_id, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  })
})


#--- Current term point ----



## Aggregate the KO and WT DEGs into the same dataframe and flag significant changes
volc_comps <- lapply(all_cts, function(ct_i){
  # Iterate through each celltype/cluster
  print(ct_i)
  
  grp_degs <- lapply(names(compare_grp), function(grp_id){
    # Iterate through each comparison group (e.g. LN_day3: KO and WT)
    grp_i <- compare_grp[[grp_id]][1:2]  # where grp_id = 'LN_day3'
    grp_pattern <- compare_grp[[grp_id]][3]  # where pattern = "_.*$", or "^.*_"
    celltype <- gsub("_.*", "", grp_id)
    
    kowt_deg <- lapply(setNames(grp_i,grp_i), function(comp_id){
      # comp_id = 'day3_KO'
      # celltype = 'LN'
      print(comp_id)
      i <- ct_markers_all[[celltype]][[comp_id]]
      if(is.null(i[[ct_i]])){
        colids <- colnames(i[[1]])
        x <- as.data.frame(matrix(ncol=7, nrow=0)) %>%
          magrittr::set_colnames(., c('symbol', 'ensembl', 'biotype', 
                                      'avg_log2FC', 'sig', 'log10p', 'n')) %>%
          rename_with(., ~paste0(., ".", gsub(grp_pattern, "", comp_id)), 
                      .cols = -c('symbol', 'ensembl', 'biotype'))
        x$symbol <- as.character(x$symbol)
        x$ensembl <- as.character(x$ensembl)
        x$biotype <- as.character(x$biotype)
        return(x)
      }
      
      # Find the IDs associated with the DEG comparison
      ids <- c(unique(i[[ct_i]]$ident.1),
               unique(i[[ct_i]]$ident.2))
      
      # Find the number of cells that were compared
      clus_idx <- as.character(seul[[celltype]]@meta.data[,anno_ident]) == as.character(unique(i[[ct_i]]$cluster))
      clus_idx <- if(!any(clus_idx)) (!clus_idx) else (clus_idx)
      n <- setNames(c(sum((seul[[celltype]]$orig.ident == ids[1]) & clus_idx),
                      sum((seul[[celltype]]$orig.ident == ids[2]) & clus_idx)), 
                    ids) %>%
        paste(., collapse="/")
      
      # Relabel and subset key columns for the DEGs
      i[[ct_i]] %>%
        mutate(sig=case_when((avg_log2FC >= fc_sig) & (p_val_adj <= pval_sig) ~ "up",
                             (avg_log2FC <= (-1*fc_sig)) & (p_val_adj <= pval_sig) ~ "down",
                             TRUE ~ "ns"),
               biotype=sym2biotype_ids[symbol],
               ensembl=gm$SYMBOL$ENSEMBL[symbol],
               log10p=(-1*log10(p_val_adj + pseudop)),
               group=paste0(celltype, "_", gsub("_.*", "", comp_id)),
               n=paste(paste0("n_", names(n), "=", n), collapse=";")) %>%
        select(symbol,ensembl, biotype,  avg_log2FC, sig, p_val_adj, log10p, n) %>%
        rename_with(., ~paste0(., ".", gsub(grp_pattern, "", comp_id)), 
                    .cols = -c('symbol', 'ensembl', 'biotype'))
    })
    return(kowt_deg)
  })
  names(grp_degs) <- names(compare_grp)
  
  #  if(any(sapply(volc_comp[c("LN_KO", "LN_WT", 
  #                           "Tumor_KO", "Tumor_WT")], is.null))){
  #   print("null")
  #   return(NULL)
  # }
  
  reduced_degs <- lapply(names(grp_degs), function(grp_id){
    print(grp_id)
    grp_i <- compare_grp[[grp_id]][1:2]  # where grp_id = 'LN_day3'
    grp_pattern <- compare_grp[[grp_id]][3]  # where pattern = "_.*$", or "^.*_"
    ids <- gsub(grp_pattern, "", grp_i)
    
    grp_deg_x <- grp_degs[[grp_id]]
    celltype <- gsub("_.*", "", grp_id)  # e.g. LN_day3 -> LN
    
    # grp_deg_x is a list of two dataframes, containing:
    # KO in first slot and WT in second slot
    if(all(sapply(grp_deg_x, nrow)==0)) return(NULL)
    ln_lfc <- grp_deg_x %>% 
      purrr::reduce(full_join, by=c('symbol', 'ensembl', 'biotype'))
    id1 <- paste0("sig.", ids[1])
    id2 <- paste0("sig.", ids[2])
    sig <- paste0(paste0(ids[1], ".", ln_lfc[,id1]), "_",
                  paste0(ids[2], ".", ln_lfc[,id2]))
    ln_lfc <- ln_lfc %>% 
      mutate(sig=sig,
             celltype=celltype) %>% 
      select(-c(id1, id2))
    ln_lfc$sig <- gsub("NA", "ns", ln_lfc$sig)
    return(ln_lfc)
  })
  names(reduced_degs) <- names(grp_degs)
  return(reduced_degs)
})
names(volc_comps) <- all_cts

## Create a delta-FC metric between KO and WT
volc_comps_dfc <- lapply(volc_comps, function(ct_dat){
  # Iterating through each celltype/cluster
  ## Get the mean and sd delta-FC across all days and samples
  all_dfc <- sapply(ct_dat, function(ct_grp_i){
    tryCatch({
      lfc_id <- grep("avg_log2FC", colnames(ct_grp_i), value=T)
      (2^ct_grp_i[,lfc_id[1]]) - (2^ct_grp_i[,lfc_id[2]])
    }, error=function(e){NULL})
  })
  sd_dfc <- sd(unlist(all_dfc), na.rm=T)
  mean_dfc <- mean(unlist(all_dfc), na.rm=T)
  
  ### Get the mean and sd delta-FC across all days, split by Tumor and LN
  # cts <- c('Tumor', 'LN')
  # meansd_dfc_cts <- lapply(setNames(cts,cts), function(ct){
  #   all_dfc <- sapply(ct_dat[grep(ct, names(ct_dat))], function(ct_grp_i){
  #     tryCatch({
  #       (2^ct_grp_i$avg_log2FC.KO) - (2^ct_grp_i$avg_log2FC.WT)
  #     }, error=function(e){NULL})
  #   })
  #   sd_dfc <- sd(unlist(all_dfc), na.rm=T)
  #   mean_dfc <- mean(unlist(all_dfc), na.rm=T)
  #   c('mean'=mean_dfc, "sd"=sd_dfc)
  # })
  
  lapply(ct_dat, function(ct_grp_i){
    # Iterating through each group (e.g. LN_Day3)
    if(is.null(ct_grp_i)) return(NULL)
    # mean_dfc <- meansd_dfc_cts[[unique(ct_grp_i$celltype)]]['mean']
    # sd_dfc <- meansd_dfc_cts[[unique(ct_grp_i$celltype)]]['sd']
    
    # Calculate the euclidean distance between FC (KO vs WT)
    lfc_id <- grep("avg_log2FC", colnames(ct_grp_i), value=T)
    ct_grp_i$deltafc <- (2^ct_grp_i[,lfc_id[1]]) - (2^ct_grp_i[,lfc_id[2]])
    
    # Assume normally distributed, use parameters to find p-values
    lowp <- with(ct_grp_i, pnorm(deltafc, mean=mean_dfc, sd=sd_dfc, lower.tail=T)*2)
    highp <- with(ct_grp_i, pnorm(deltafc, mean=mean_dfc, sd=sd_dfc, lower.tail=F)*2)
    idx <- (ct_grp_i$deltafc < mean_dfc)
    deltafc_p <- lowp 
    deltafc_p[which(!idx)] <- highp[which(!idx)]
    
    # Assign significance-labels based on FC and padj
    ct_grp_i <- ct_grp_i %>% 
      mutate(deltafc_p=deltafc_p,
             deltafc_padj=p.adjust(deltafc_p, method='BH'),
             deltafc_log10padj=-log10(deltafc_padj + pseudop))
    sig <- with(ct_grp_i, 
                ifelse(deltafc_padj > pval_sig, "n.s.",
                       ifelse(deltafc >= fc_sig, "sig.up", 
                              ifelse(deltafc <= fc_sig, "sig.down", "sig"))))
    ct_grp_i$deltafc_sig <- sig
    # ct_grp_i[which(ct_grp_i$deltafc_padj < 0.05),]
    return(ct_grp_i)
  })
})

## Create a delta-dayFC metric between 7d and 3d (e.g. Day7 KOvsWT: Day3 KOvsWT)
# Initial step to merge together day7 and day3 for tumor and ln independently
day73_dfc_dat <- lapply(volc_comps_dfc, function(ct_dat){
  ## Merge together the dFC for Day7 and Day3, putting Day7 first
  tln <- c('Tumor', 'LN')
  lapply(setNames(tln,tln), function(tissueid){
    day_ids <- grep(tissueid, names(ct_dat), value=T)
    day_ids <- day_ids[order(grepl("day3", day_ids))]
    if(any(sapply(ct_dat[day_ids], is.null))) return(NULL)
    day_tumor_ln <- lapply(day_ids, function(day_id){
      ct_dat[[day_id]] %>% 
        dplyr::select(symbol, ensembl, biotype, deltafc, deltafc_p, deltafc_padj, deltafc_sig) %>% 
        dplyr::rename_with(., .fn = ~ paste0(gsub("^.*_", "", day_id), ".", .), 
                           .cols = -c('symbol', 'ensembl', 'biotype'))
    }) %>% 
      purrr::reduce(., full_join, by=c('symbol', 'ensembl', 'biotype'))
    return(day_tumor_ln)
  } )
})

day73_dfc <- lapply(day73_dfc_dat, function(day73_i){
  ## Get the mean and sd delta-FC across all days and samples
  all_dfc <- sapply(day73_i, function(day73_ij){
    tryCatch({
      (day73_ij$day7.deltafc) - (day73_ij$day3.deltafc)
    }, error=function(e){NULL})
  })
  sd_dfc <- sd(unlist(all_dfc), na.rm=T)
  mean_dfc <- mean(unlist(all_dfc), na.rm=T)
  
  # Iterating through each tissuetype
  day73_i <- lapply(day73_i, function(day73_ij){
    day73_ij$day_deltafc <- (day73_ij$day7.deltafc) - (day73_ij$day3.deltafc)
    
    # Assume normally distributed, use parameters to find p-values
    lowp <- with(day73_ij, pnorm(day_deltafc, mean=mean_dfc, sd=sd_dfc, lower.tail=T)*2)
    highp <- with(day73_ij, pnorm(day_deltafc, mean=mean_dfc, sd=sd_dfc, lower.tail=F)*2)
    idx <- (day73_ij$day_deltafc < mean_dfc)
    deltafc_p <- lowp 
    deltafc_p[which(!idx)] <- highp[which(!idx)]
    
    # Assign significance-labels based on FC and padj
    day73_ij <- day73_ij %>% 
      mutate(day_deltafc_p=deltafc_p,
             day_deltafc_padj=p.adjust(deltafc_p, method='BH'),
             day_deltafc_log10padj=-log10(day_deltafc_padj + pseudop))
    sig <- with(day73_ij, 
                ifelse(day_deltafc_padj > pval_sig, "n.s.",
                       ifelse(day_deltafc >= fc_sig, "sig.up", 
                              ifelse(day_deltafc <= fc_sig, "sig.down", "sig"))))
    day73_ij$day_deltafc_sig <- sig
    return(day73_ij)
  })
  return(day73_i)
})

### Writing out DEG tables:
wtbl <- function(dat, file){
  write.table(dat, file=file, sep=",", quote = F, col.names = T, row.names = F)
}
for(each_ct in all_cts){
  print(each_ct)
  # wtbl(day73_dfc[[each_ct]]$Tumor, file.path("~/xfer", paste0("degtable_Day7vs3Tumor.", each_ct,".csv")))
  # wtbl(day73_dfc[[each_ct]]$LN, file.path("~/xfer", paste0("degtable_Day7vs3LN.", each_ct, ".csv")))
  for(id in names(volc_comps_dfc[[each_ct]])){
    wtbl(volc_comps_dfc[[each_ct]][[id]], 
         file.path("~/xfer", paste0("degtable_", id, ".", each_ct, ".csv")))
  }
}



## Create the correlation barplots
grps_corr_ggs <- lapply(names(compare_grp), function(grp_id){
  # Create the correlation dataframe
  corr_df <- lapply(volc_comps, function(i){
    df <- i[[grp_id]]
    res <- tryCatch({
      corres <- with(df, cor.test(avg_log2FC.KO, avg_log2FC.WT, 
                                  method='pearson', use='complete.obs'))
      data.frame('r'=corres$estimate, 'p'=corres$p.value, 'celltype'=unique(df$celltype))
    }, error=function(e){
      data.frame("r"=NA, "p"=NA, "celltype"=NA)
    })
  }) %>% 
    do.call(rbind, .) %>%
    as.data.frame %>%
    mutate(padj=p.adjust(p, method='BH'),
           log10p=-log10(padj+pseudop),
           sig=ifelse(padj < pval_sig, 'sig', 'nonSig')) %>%
    tibble::rownames_to_column(., "Cluster")
  
  # Make the ggplot barplot
  corr_df$sig[is.na(corr_df$sig)] <- 'nonSig'
  corr_df[is.na(corr_df)] <- 0
  gg <- ggplot(corr_df, aes(y=Cluster, x=r, fill=sig)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('nonSig'='black', 'sig'='darkgoldenrod')) +
    theme_cowplot() + 
    xlim(-1,1) +
    ggtitle(grp_id) + 
    geom_vline(xintercept = 0, linetype='dashed')
  list("df"=corr_df, "gg"=gg)
})
names(grps_corr_ggs) <- names(compare_grp)

# Combine all the correlation dataframes into one large dataframe
all_df <- lapply(names(grps_corr_ggs), function(i){
  grps_corr_ggs[[i]]$df %>%
    mutate(group=i)
}) %>% do.call(rbind, .)


write.table(all_df, file=file.path("~/xfer", paste0("corr_dat.", anno_ident, ".csv")),
            quote = F, sep=",", col.names = T, row.names = F)
pdf(paste0("~/xfer/corr_bars.", anno_ident, ".pdf"))
lapply(grps_corr_ggs, function(i) i$gg)
dev.off()


## Create volcano plots for the delta-FC between KO and WT 
# where KO and WT had their log2FC between Treated and Untreated calculated
ct_i = 'All'
grps <- list("Day7"=c("LN_day7", "Tumor_day7"),
             "Day3"=c("LN_day3", "Tumor_day3"))
make_gg_volcano_plot <- function(ct_df, title=NULL,  xvar='deltafc', 
                                 yvar='deltafc_log10padj', fillvar='deltafc_sig'){
  ggplot(ct_df, aes_string(x=xvar, y=yvar,fill=fillvar)) +
    facet_grid(.~celltype) +
    geom_point(shape = 21, colour = "black") +
    theme_classic() +
    xlim(-4,4) + xlab("Delta FC") + ylab("-log10(padj)") +
    ggtitle(paste0(ct_i, ": ")) + #,suffix))#comp_id, suffix)) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed") +
    scale_fill_manual(values = c("n.s." = "grey",
                                 "sig" = "white",
                                 "sig.down" = "#80b1d3",
                                 "sig.up" = "#fb8072")) +
    # scale_fill_manual(values = c("KO.down_WT.down" = "#fdb462", 
    #                              "KO.up_WT.up" = "#fdb462",
    #                              "KO.down_WT.up" = "#80b1d3", 
    #                              "KO.up_WT.down" = "#fb8072")) +
    ylab("-log10p") + xlab("log2FC") +
    ggtitle(title)
}
filter_for_wtsig <- FALSE
volc_params <- list("KO"=c('avg_log2FC.KO', 'p_val_adj.KO', 'log10p.KO', 'KO'),
                    "WT"=c('avg_log2FC.WT', 'p_val_adj.WT', 'log10p.WT', 'WT'),
                    "Delta"=c('deltafc', 'deltafc_padj', 'deltafc_log10padj', 'KO-WT'))
gg_volcanos <- lapply(grps, function(grp_i){
  ct_df <- rbind(volc_comps_dfc[[ct_i]][[grp_i[1]]],
                 volc_comps_dfc[[ct_i]][[grp_i[2]]]) 
  if(filter_for_wtsig){
    ct_df <- ct_df %>%
      filter(!grepl("WT.ns$", sig))
  }
  
  all_volcano <- lapply(volc_params, function(params_i){
    # params_i <- c('avg_log2FC.KO', 'p_val_adj.KO', 'log10p.KO', 'KO')
    ct_df$sig <- ifelse(ct_df[,params_i[2]] < pval_sig, 
                        ifelse(ct_df[,params_i[1]] > fc_sig, "sig.up",
                               ifelse(ct_df[,params_i[1]] < (-1*fc_sig), "sig.down", "sig")),
                        "n.s.")
    day <- stringr::str_to_title(gsub("^.*_", "", grp_i[1]))  # e.g. LN_day7 -> Day7
    # title <- paste0(day, ", Cells: ", ct_i, " | delta FC: [KO:", day, "vsUn] - [WT:", day, "vsUn]")
    title <- paste0(day, ", ST2_", params_i[4], " Cells: ", ct_i, " | Treated vs Untreated")
    
    
    # Volcano plot (ggplot)
    gg_volcano <- make_gg_volcano_plot(ct_df, title=title,
                                       xvar=params_i[1], yvar=params_i[3], fillvar='sig')
    return(gg_volcano)
  })
  
  # Get the number of DEGs per celltype
  n_df <- lapply(setNames(all_cts,all_cts), function(ct_i){
    # print(ct_i)
    ndf_i <- lapply(grp_i, function(i){
      # print(i)
      df1 <- volc_comps_dfc[[ct_i]][[i]] 
      if(filter_for_wtsig){
        if(!is.null(df1)) df1 <- df1 %>% filter(!grepl("WT.ns$", sig))
      }
      n <- c("cluster"=ct_i,
             "sample"=gsub("_.*", "", i), "day"=gsub("^.*_", "", i),
             'up'=sum(df1$deltafc_sig == 'sig.up', na.rm=T), 
             'down'=sum(df1$deltafc_sig == 'sig.down', na.rm=T))
      return(n)
    }) %>% 
      do.call(rbind, .)
    return(ndf_i)
  }) %>%
    do.call(rbind,.) %>%
    as.data.frame %>%
    mutate(up=as.numeric(up),
           down=as.numeric(down))
  
  return(list("all_volcano"=all_volcano, "volcano_csv"=ct_df, "n_df"=n_df))
})
pdf(file.path("~/xfer", paste0("volcanos.allCells.pdf")), 
    height = 4, width = 7)
lapply(gg_volcanos, function(i) i$all_volcano$Delta)
dev.off()
write.table(gg_volcanos$volcano_csv, file=file.path("~/xfer", "volcanos.allCells.csv"), 
            sep=",", col.names = T, row.names = F, quote=F)

pdf(file.path("~/xfer", paste0("volcanos.Tregs.KO_WT.pdf")), 
    height = 8, width = 7)
lapply(gg_volcanos, function(i) cowplot::plot_grid(plotlist=i$all_volcano[c('KO', 'WT')], nrow=2))
dev.off()
day3_treg_deg <- gg_volcanos$Day3$volcano_csv %>%
  filter(grepl("up|down", sig)) %>%
  select(-grep("delta", colnames(.), value=T))
write.table(day3_treg_deg, file=file.path("~/xfer", "volcanos.Treg.Day3.csv"), 
            sep=",", col.names = T, row.names = F, quote=F)
day7_treg_deg <- gg_volcanos$Day7$volcano_csv %>%
  filter(grepl("up|down", sig)) %>%
  select(-grep("delta", colnames(.), value=T))
write.table(day7_treg_deg, file=file.path("~/xfer", "volcanos.Treg.Day7.csv"), 
            sep=",", col.names = T, row.names = F, quote=F)






grps <- list("Day73"=c("LN_day3_7", "Tumor_day3_7"))
ct_i = 'All'
ct_df <- rbind(volc_comps_dfc[[ct_i]][[grps[[1]][1]]],
               volc_comps_dfc[[ct_i]][[grps[[1]][2]]]) 
ct_sig <- strsplit(ct_df$sig, split="_") %>% do.call(rbind, .) %>%
  as.data.frame %>%
  magrittr::set_colnames(., c('day3_sig', 'day7_sig'))
ct_df <- as.data.frame(cbind(ct_df, ct_sig)) %>%
  mutate(day3_sig=ifelse(p_val_adj.day3 < pval_sig, 
                         ifelse(avg_log2FC.day3 > fc_sig, "sig.up",
                                ifelse(avg_log2FC.day3 < (-1*fc_sig), "sig.down", "sig")),
                         "n.s."),
         day7_sig=ifelse(p_val_adj.day7 < pval_sig, 
                         ifelse(avg_log2FC.day7 > fc_sig, "sig.up",
                                ifelse(avg_log2FC.day7 < (-1*fc_sig), "sig.down", "sig")),
                         "n.s."))
title <- paste0("Cells: ", ct_i, " | delta FC: [Day3:KOvsWT] - [Day7:KOvsWT]")
gg_volcano1 <- make_gg_volcano_plot(ct_df, title=title)
title <- paste0("Cells: ", ct_i, " | Day3: KO-vs-WT")
gg_volcano2 <- make_gg_volcano_plot(ct_df, title=title, xvar='avg_log2FC.day3', 
                                    yvar='log10p.day3', fillvar='day3_sig')
title <- paste0("Cells: ", ct_i, " | Day7: KO-vs-WT")
gg_volcano3 <- make_gg_volcano_plot(ct_df, title=title, xvar='avg_log2FC.day7', 
                                    yvar='log10p.day7', fillvar='day7_sig')
pdf("~/xfer/x3.pdf", height = 4, width = 7)
gg_volcano1
gg_volcano2
gg_volcano3
dev.off()

# Parse the volcano analysis for the number of DEG of delta-FC and barplot it
all_ndf <- lapply(gg_volcanos, function(i) i$n_df) %>% 
  do.call(rbind, .) %>%
  as.data.frame
pdf(file.path("~/xfer", paste0("n_df.", anno_ident, ".pdf")))
ggplot(all_ndf, aes(y=cluster, x=(up+down))) +
  facet_grid(sample~day, scales='fixed', space='free') +
  theme_cowplot() +
  geom_bar(stat='identity', position='dodge')
dev.off()
write.table(all_ndf, file=file.path("~/xfer", paste0("n_df.", anno_ident, ".csv")),
            sep=",", quote=F, row.names = F, col.names = T)

table(seul$LN$manual_anno, gsub("_.*_", "", seul$LN$orig.ident)) %>%
  as.matrix %>%
  apply(., 1, function(i) i[which.min(i)])
table(seul$Tumor$manual_anno, gsub("_.*_", "", seul$Tumor$orig.ident)) %>%
  as.matrix %>%
  apply(., 1, function(i) i[which.min(i)])

pdf(paste0("~/xfer/corr_dimplots.pdf"), width = 12)
lndp1 <- DimPlot(seul$LN, group.by='manual_clusters', raster=T, label=T, reduction='umap') + ggtitle("LN")
lndp2 <- DimPlot(seul$LN, group.by='manual_anno', raster=T, label=T, reduction='umap') + ggtitle("LN")
tdp1 <- DimPlot(seul$Tumor, group.by='manual_clusters', raster=T, label=T, reduction='umap') + ggtitle("Tumor")
tdp2 <- DimPlot(seul$Tumor, group.by='manual_anno', raster=T, label=T, reduction='umap') + ggtitle("Tumor")
lndp1 + lndp2
tdp1 + tdp2
dev.off()

## OUTDATED SCATTERPLOT: Replaced with Correlation plot
{
  ggscatter <- function(df, cols=NULL, title=NULL){
    if(is.null(cols)) {
      cols=c("KO.down_WT.down" = "#fdb462", 
             "KO.down_WT.up" = "#80b1d3", 
             "KO.up_WT.up" = "#fdb462", 
             "KO.up_WT.down" = "#fb8072")
    }
    ggplot(df, aes(x=avg_log2FC.KO, y=avg_log2FC.WT, fill=sig)) +
      facet_grid(.~celltype) +
      geom_point(shape = 21, colour = "black", alpha=0.9) +
      theme_classic() +
      xlim(-5,5) + ylim(-5,5) +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = cols) +
      ylab(paste0("ST2_wt log2FC\n", unique(df$n.KO))) +
      xlab(paste0("ST2_KO log2FC\n", unique(df$n.WT))) +
      # ylab("ST2_wt log2FC [72h/PBS]") + xlab("ST2_KO log2FC [72h/PBS]") +
      ggtitle(if(is.null(title)) ct_i else paste0(title, ": ", ct_i)) + 
      NoLegend()
  }
  gg_scatter_plots <- lapply(names(reduced_degs), function(id){
    ggscatter(reduced_degs[[id]], title=id)
  })
  pdf("~/xfer/x.pdf")
  gg_scatter_plots
  dev.off()
  
  
  
  
  gg_volcano <- ggplot(reduced_degs[[1]] %>% 
                         filter(!grepl("WT.ns", sig)), 
                       aes(x=avg_log2FC.KO, y=log10p.KO, fill=sig)) +
    facet_grid(.~celltype) +
    geom_point(shape = 21, colour = "black") +
    theme_classic() +
    xlim(-4,4) +
    ggtitle(paste0(ct_i, ": ")) + #,suffix))#comp_id, suffix)) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed") +
    scale_fill_manual(values = c("KO.down_WT.down" = "#fdb462", 
                                 "KO.down_WT.up" = "#80b1d3", 
                                 "KO.up_WT.up" = "#fdb462", 
                                 "KO.up_WT.down" = "#fb8072")) +
    ylab("-log10p") + xlab("log2FC") +
    ggtitle(paste0("Cells: ", ct_i, " | ST2_KO [72h/PBS], subset for sig.ST2_WT"))
  
  return(list("gg_scatter"=gg_scatter, 
              "gg_scatter_ln"=gg_scatter_ln,
              "gg_scatter_tumor"=gg_scatter_tumor,
              "gg_volcano"=gg_volcano,
              "ln_lfc"=ln_lfc,
              "tumor_lfc"=tumor_lfc))
  
  
  # 
  names(volc_comps) <- all_cts
  pdf("~/xfer/LN_Tumor_lfc.pdf", height = 4, width = 7)
  lapply(volc_comps, function(i) i$gg_scatter)
  dev.off()
  pdf("~/xfer/LN_lfc.pdf", height = 10, width = 8)
  cowplot::plot_grid(plotlist=lapply(volc_comps[-1], function(i) i$gg_scatter_ln), ncol=3, nrow=4)
  dev.off()
  pdf("~/xfer/Tumor_lfc.pdf", height = 10, width = 8)
  cowplot::plot_grid(plotlist=lapply(volc_comps[-1], function(i) i$gg_scatter_tumor), ncol=3, nrow=4)
  dev.off()
  
  pdf("~/xfer/LN_subset_volcano.pdf", height = 4, width = 7)
  lapply(volc_comps, function(i) i$gg_volcano)
  dev.off()
  
}

## Print the significant gene list to run in enrichr
.printSymbolOv <- function(df){
  # Select all non-significant genes
  x <- df %>% 
    filter(!grepl("ns|NA", sig)) %>%
    split(., .$sig) 
  
  # Print the gene list for enrichr
  for(id in names(x)){
    print(id)
    x[[id]] %>% select(symbol) %>% print(., row.names = FALSE)
  }
}
.printSymbolOv(volc_comps[['Tregs; Intermediate']]$ln_lfc)
.printSymbolOv(volc_comps[['Tregs; Cxcl10_hi']]$ln_lfc)

#--- b) Differential expression between select clusters ----
grp_comparisons <- function(group){
  switch(group,
         LN=list(c("5", "6", "7", "4", "0", "9", "13"),
                 c('15', '20', '19')),
         Tumor=list(c('4', '1', '7', '0', '5', '11', '3'),
                    c('15', '10'),
                    c('6', '18')))
}

seul_degs <- lapply(names(seul)[1], function(seu_type){
  seu <- seul[[seu_type]]
  grp_comps <- grp_comparisons(seu_type)
  
  ids <- sapply(grp_comps, paste, collapse="_")
  grp_degs <- lapply(grp_comps, function(grp_comp_i){
    id <- paste(grp_comp_i, collapse="_")
    dir.create(file.path(outdirdeg, "identify_clusterid"), showWarnings = F)
    fileid <- paste0("deg_", seu_type, "_", anno_ident,".", id, ".rds")
    outf <- file.path(outdirdeg, "identify_clusterid",fileid)
    print(fileid)
    
    if(file.exists(outf)){
      markers_grp <- readRDS(file=outf)
    } else {
      Idents(seu) <- anno_ident
      seu_sub <- subset(seu, ident=grp_comp_i)
      markers_grp <- .betweenClustersFindMarkers(seu_sub, anno_ident)
      saveRDS(markers_grp, file=outf)
    }
    cl_markers <- .identifyClusterSpecificMarkers(markers_grp, slot='score', n=20)
    
    return(cl_markers)
  })
  names(grp_degs) <- ids
  return(grp_degs)
})
names(seul_degs) <- names(seul)[1]

#--- c) Differential expression between groups across all clusters ----
# Iterate through Tumor and LN separately to calculate DEGs between treatment times
seul_degs <- lapply(names(seul), function(seu_type){
  seu <- seul[[seu_type]]
  outf=file.path(outdirdeg, paste0("degs.", seu_type, ".", anno_ident, ".rds"))
  if(file.exists(outf)){
    ct_markers_all <- readRDS(outf)
    return(ct_markers_all)
  }
  
  # Create the groupings
  seu$tissue_cond <- gsub("_[a-zA-Z0-9]*$" , "", seu$orig.ident) %>% gsub("^B[12]_", "", .)
  seu$kowt_comp <- gsub("WT_|KO_", "", seu$orig.ident) %>% gsub("^B[12]_", "", .)
  seu.list <- SplitObject(seu, split.by = "tissue_cond") # Compares 72h timepoint to BL
  seu.list2 <- SplitObject(seu, split.by = "kowt_comp")  # Compares KO to WT at each timepoint
  seu.lists <- c(seu.list2, seu.list)
  seu.lists <- seu.list
  rm(seu.list, seu.list2); gc()
  
  # Iterate through all group comparisons (e.g. WT vs KO, or 7d/3d vs Un)
  ct_markers_all <- lapply(names(seu.lists), function(seu_id){
    print(seu_id)
    # seu_id <- names(seu.lists)[4]
    seu_i <- seu.lists[[seu_id]]
    DefaultAssay(seu_i) <- 'RNA'
    all_ids <- unique(seu_i$orig.ident)
    all_ids_comb <- combn(all_ids, m=2)
    all_ids_comb
    
    colid <- c('tissue_cond', 'kowt_comp')
    colid <- colid[grep(seu_id, seu_i@meta.data[,colid])]
    # Create a grouping of all comparisons to be made for DEG testing
    idents <- switch(colid,
                     kowt_comp={
                       comb <- apply(all_ids_comb, 2, function(i){
                         fac <- factor(gsub("^.*_(KO|WT)_.*$", "\\1", i),
                                       levels=c('KO', 'WT'))
                         i[order(fac)]
                       })
                       # remove non-informative comparisons
                       cs1 <- colSums(gsub("^.*_(KO|WT)_.*$", "\\1", comb) == 'KO')
                       cs2 <- colSums(gsub("_.*", "", comb) == 'B1')
                       comb[,which(cs1==1 & cs2!=1), drop=F]
                     },
                     tissue_cond={
                       comb <- apply(all_ids_comb, 2, function(i){
                         fac <- factor(gsub("^.*_", "", i),
                                       levels=c('7d', '3d', 'Un'))
                         i[order(fac)]
                       })
                       # remove non-informative comparisons
                       cs1 <- colSums(matrix(grepl("7d|3d", comb), nrow=nrow(comb)))
                       cs2 <- colSums(gsub("_.*", "", comb) == 'B1')
                       comb[,-which(cs2==1 & cs1!=2), drop=F]
                     })
    
    # Isolate individual clusters and celltypes
    clusters <- unique(as.character(seu_i@meta.data[,anno_ident]))
    cl_marker_partitions <- lapply(clusters, function(cl_j){
      cells_j <- Cells(seu_i)[which(seu_i@meta.data[,anno_ident] == cl_j)]
      seu_ij <- subset(seu_i, cells=cells_j)
      
      # DEG across all cells within that cluster
      Idents(seu_ij) <- 'orig.ident'
      marker_partitions <- apply(idents, 2, function(ident_i){
        subset_of_seu <- TRUE
        marker_partition <- tryCatch({
          FindMarkers(seu_ij, assay = "RNA", test.use='wilcox',
                      ident.1= ident_i[1], ident.2= ident_i[2],
                      verbose = FALSE,
                      logfc.threshold = 0,
                      recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
            tibble::rownames_to_column(., "symbol") %>%
            mutate(ident1=ident_i[1],
                   ident2=ident_i[2],
                   cluster=cl_j,
                   clustertype=anno_ident)
        }, error=function(e){ NULL })
        return(marker_partition)
      })
      return(marker_partitions)
    })
    setNames(cl_marker_partitions) <- clusters
    return(cl_marker_partitions)
  })
  names(ct_markers_all) <- names(seu.lists)
  
  saveRDS(ct_markers_all,  file=outf)
})
names(seul_degs) <- names(seul)


#--- d) Checking stability of DEGs across clusters ----
ident_to_check <- 'Treg'
genes_to_check <- c('Tox', 'Tnfrsf18', 'Mki67', 'Maf', 'Lgals1', 'Klrg1', 
                    'Itgae', 'Il4i1', 'Il2rb', 'Il2ra', 'Ikzf4', 'Ikzf2',
                    'Icos', 'Gzmb', 'Gata3', 'Foxp3', 'Cxcr3', 'Ctla4',
                    'Cd69', 'Ccr7', 'Ccr2')
fps <- lapply(names(seul), function(seu_id){
  seu <- seul[[seu_id]]
  Idents(seu) <- 'manual_anno'
  cnts <- GetAssayData(seu, slot='counts')
  seu$genes_cnt <- colSums(cnts[genes_to_check,] >0)
  seu$genes_sum <- log2(colSums(cnts[genes_to_check,]) + 1)
  
  
  fp1=FeaturePlot(seu, features='genes_cnt', 
          split.by='orig.ident', reduction='umap', raster=T) +
    ggtitle(seu_id)
  fp2=FeaturePlot(seu, features='genes_sum', 
              split.by='orig.ident', reduction='umap', raster=T)
  fp12 <- plot_grid(fp1, fp2, nrow=2)
  return(fp12)
})
pdf("~/xfer/number_of_genes_expr.pdf", width = 17, height = 7)
fps
dev.off()

summcnts <- lapply(names(seul), function(seu_id){
  seu <- seul[[seu_id]]
  Idents(seu) <- 'manual_anno'
  seusub <- subset(seu, ident='Treg')
  cnts <- GetAssayData(seusub, slot='counts')
  hi_idx <- which((rowSums(cnts>0)/nrow(seusub)) > 0.20)
  keep_idx <- which(rownames(seusub) %in% genes_to_check)
  cnts_hi <- cnts[unique(c(hi_idx, keep_idx)),]
  
  spl_g2c <- split(as.data.frame(t(cnts_hi[genes_to_check,])), f=seusub$seurat_clusters)
  
  spl_rand <- split(as.data.frame(t(cnts_hi)), f=seusub$seurat_clusters)
  summcnts <- lapply(names(spl_rand), function(j){
    print(j)
    rand_clus_j <- spl_rand[[as.character(j)]]
    if(nrow(rand_clus_j)==0) return(NULL)
    
    rand_cnts <- sapply(1:100, function(i){
      set.seed(i)
      ridx <- sample(1:nrow(cnts_hi), size=length(genes_to_check))
      rowSums(rand_clus_j[,ridx] > 0)
    }) %>% as.numeric
    
    g2c_cnts <- as.numeric(rowSums(spl_g2c[[as.character(j)]] > 0))
    tt <- t.test(g2c_cnts, rand_cnts)
    c(tt$estimate, "median_x"=median(g2c_cnts),
      "median_y"=median(rand_cnts), 'cluster'=j)
  }) %>% do.call(rbind, .)
  storage.mode(summcnts) <- 'numeric'
  
  return(summcnts)
})
  
# seusub <- ScaleData(seusub, features=genes_to_check)
# pdf("~/xfer/x.pdf")
# DoHeatmap(seusub, features=genes_to_check, group.by='orig.ident', raster=T)
# dev.off()


#--- e) Comparing DEGs between 3d analysis and current ----
seul_tregs <- readRDS(file=file.path(datadir, "seurat_obj", "tregs.rds"))
old_dir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor'
seu3d <- readRDS(file = file.path(old_dir, "data", "seurat_obj", "3_seu_degPrepX.rds"))
recode_map3d <- function(x){
  x %>% 
    recode("0"='Tregs',
           "1"="B",
           "2"='Tregs',
           "3"="B",
           "4"='Tregs',
           "5"='Tregs; Cycling',
           "6"='Tregs; Cxcl10_hi',
           "7"='Tregs; Cycling',
           "8"="DC",
           "9"='Tregs; Intermediate',
           "10"='CD4 Th2',
           "11"='CD8 Effector',
           "12"='Monocytes',
           "13"='CD8 Naive',
           "14"="B",
           "15"='Tregs; highMt',
           "16"='Tregs; Intermediate',
           "17"='CD4 Naive',
           "18"='Tregs; Unknown',
           "19"='B.GC',
           "20"='Eosinophils',
           "21"='CD4 Th2',
           "22"='NKcells')
}
seu3d$manual_anno <- seu3d$seurat_clusters %>% recode_map3d

.mapCellsFromOldSeuratToNew <- function(seu.old, seu.new, old.col, new.col){
  cells.old <- split(Cells(seu.old), seu.old@meta.data[,old.col])
  cells.old <- lapply(names(cells.old), function(i){
    data.frame("celltype"=i,
               "cells"=cells.old[[i]]) %>%
      tibble::column_to_rownames(., 'cells')
  }) %>% do.call(rbind, .)
 
  seu.new@meta.data[,new.col] <- NA
  intersect_cells <- intersect(Cells(seu.new), rownames(cells.old))
  seu.new@meta.data[intersect_cells,new.col] <- cells.old[intersect_cells,,drop=F]$celltype
  return(seu.new)
}

seul <- lapply(seul, function(seu.new){
  .mapCellsFromOldSeuratToNew(seu.old=seu3d, seu.new=seu.new, 
                              old.col='manual_anno', new.col='anno_3d')
})
seul <- lapply(seul, function(seu.new){
  .mapCellsFromOldSeuratToNew(seu.old=merge(seul_tregs[[1]], seul_tregs[-1]), seu.new=seu.new, 
                              old.col='functional.cluster', new.col='treg_anno')
})
seul <- lapply(seul, function(seu.new){
  .mapCellsFromOldSeuratToNew(seu.old=merge(seul_tregs[[1]], seul_tregs[-1]), seu.new=seu.new, 
                              old.col='ssgsea_max', new.col='treg_anno_ssgsea')
})



tbls <- lapply(seul, function(seu_i){
  # t(table(seu_i$treg_anno_ssgsea, seu_i$treg_anno)) %>%
  t(table(seu_i$seurat_clusters, seu_i$treg_anno)) %>%
    as.data.frame.matrix %>%
    # select(grep("Tregs", colnames(.), value=T)) %>%
    filter(rowSums(.) > 5)
})


pdf("~/xfer/x.pdf", width = 12)
lapply(seul, DimPlot, group.by='anno_3d', label=T, raster=T, reduction='umap') %>%
  cowplot::plot_grid(plotlist=., nrow=1)
lapply(seul, DimPlot, group.by='treg_anno', label=T, raster=T, reduction='umap') %>%
  cowplot::plot_grid(plotlist=., nrow=1)
lapply(seul, DimPlot, group.by='treg_anno_ssgsea', label=T, raster=T, reduction='umap') %>%
  cowplot::plot_grid(plotlist=., nrow=1)
lapply(seul, DimPlot, group.by='seurat_clusters', label=T, raster=T, reduction='umap') %>%
  cowplot::plot_grid(plotlist=., nrow=1)
lapply(seul, function(seu_i){
  Idents(seu_i) <- 'manual_anno'
  subset(seu_i, ident='Treg') %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    DimPlot(., group.by='anno_3d', label=T, raster=T, reduction='umap')
}) %>% cowplot::plot_grid(plotlist=., nrow=1)
lapply(seul, function(seu_i){
  Idents(seu_i) <- 'manual_anno'
  subset(seu_i, ident='Treg') %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    DimPlot(., group.by='treg_anno', label=T, raster=T, reduction='umap')
}) %>% cowplot::plot_grid(plotlist=., nrow=1)
lapply(seul, function(seu_i){
  Idents(seu_i) <- 'manual_anno'
  subset(seu_i, ident='Treg') %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    DimPlot(., group.by='seurat_clusters', label=T, raster=T, reduction='umap')
}) %>% cowplot::plot_grid(plotlist=., nrow=1)
dev.off()

comps <- list('LN'=list(c('B1_LN_KO_3d', 'B1_LN_KO_Un'),
                        c('B1_LN_WT_3d', 'B1_LN_WT_Un')),
              'Tumor'=list(c('B1_Tumor_KO_3d', 'B1_Tumor_KO_Un'),
                           c('B1_Tumor_WT_3d', 'B1_Tumor_WT_Un')))
subset_id <- 'seurat_clusters'# 'seurat_clusters', 'anno_3d'
cnts <- list()
seul <- lapply(names(seul), function(seu_id){
  seu_i <- seul[[seu_id]]
  comp_i <- comps[[seu_id]]
  clus <- unique(seu_i$seurat_clusters[grep("Treg", seu_i$manual_anno)])
  # clus <- grep("Treg", na.omit(unique(seu_i$anno_3d)), value=T)
  
  degs <- lapply(comp_i, function(idents){
    lapply(clus, function(clus_i){
      Idents(seu_i) <- subset_id
      seu_ij <- subset(seu_i, ident=clus_i)
      Idents(seu_ij) <- 'orig.ident'
      
      marker_partition <- tryCatch({
        FindMarkers(seu_ij, assay = "RNA", test.use='wilcox',
                                      ident.1= idents[1], ident.2= idents[2],
                                      verbose = FALSE) %>% 
          tibble::rownames_to_column(., "symbol") %>%
          filter((p_val_adj < pval_sig) & (abs(avg_log2FC) > fc_sig))
        }, error=function(e){
          NULL
        })
      return(marker_partition)
    }) %>% 
      setNames(., clus)
  }) %>% setNames(., sapply(comp_i, paste, collapse="_"))
  cnts[[subset_id]] <- sapply(degs, function(i) sapply(i, nrow))
})
lapply(degs, function(i) sapply(i, nrow))


#############################
#### 5. Pathway analysis ####
cellgrp <- 'st2'
if(cellgrp=='st2'){
  if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
} else {
  if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
}

outdir_pway <- file.path(outdir, "pathway")
dir.create(outdir_pway, recursive = F)

# seul <- lapply(seul, recode_list, newcol='manual_clusters', grp='orig.ident', anno=F)
# seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)

anno_ident <- 'manual_anno' #'manual_anno', 'manual_clusters'
pval_sig <- 0.05
fc_sig <- 0.5

gprofiler_dir <- '/cluster/projects/mcgahalab/ref/gprofiler'
gmt <- GSA::GSA.read.gmt(file.path(gprofiler_dir, 'gprofiler_full_mmusculus.ENSG.gmt'))
gprof_ds <-setNames(gmt$genesets, 
                    paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]


msig_ds <- lapply(names(gprof_ds), function(sublvl){
  data.frame("gs_name"=sublvl,
             "entrez_gene"=gprof_ds[[sublvl]])
}) %>% do.call(rbind,.)

msig_ds$classification <- gsub(":.*", "", msig_ds$gs_name)
msig_l <- lapply(split(gm$ENSEMBL$SYMBOL[msig_ds$entrez_gene], msig_ds$gs_name), function(i){
  i[!is.na(i)]
})
msig_l <- split(msig_l, f=gsub(":.*", "", names(msig_l)))

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets
mlvl <- c('H', 'C2', 'C5')


msig_l <- msig_l['GO']
goids <- c('GO:0038172', 'GO:0004908', 'GO:0007165','GO:0002113', 
           'GO:0002114', 'GO:0038172', 'GO:0038172', 'GO:0050729',
           'GO:0042102', 'GO:0042130', 'GO:0046006', 'GO:0042129')
msig_l$GO <- msig_l[['GO']][sapply(goids, function(i) grep(i, names(msig_l$GO)))]
msig_l$GO[['amigo_il33']] <- c('Myd88', 'Il1rap', 'Traf6', 'Il33', 'Irak1', 
                               'Irak4', 'Il1rl1', 'Map3k7')
msig_l$GO <- msig_l$GO['amigo_il33']


#--- a) CD45: Pagoda2 activity score ####
cal_pagoda2 = function(counts,
                       gSets,
                       trim = 5,
                       n_cores=1,
                       min_gset_size=5){
  
  
  ### must be counts matrix !!!!!
  
  ### other parameters for knn.error.models
  # min.nonfailed = 5, min.count.threshold = 1,
  # max.model.plots = 50,
  # min.size.entries = 2000, min.fpm = 0, cor.method = "pearson",
  # verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE,
  # local.theta.fit = linear.fit, theta.fit.range = c(0.01, 100),
  # alpha.weight.power = 1/2
  
  ### other parameters for pagoda.varnorm
  # batch = NULL, prior = NULL,
  # fit.genes = NULL, minimize.underdispersion = FALSE,
  # n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9,
  # verbose = 0, weight.df.power = 1, smooth.df = -1,
  # theta.range = c(0.01, 100), gene.length = NULL
  
  nPcs = min(round(ncol(counts)/5),5)
  #counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
  tryCatch({
    p2 = Pagoda2$new(as(counts, "dgCMatrix"), n.cores = n_cores,log.scale=F, modelType="raw")
    p2$adjustVariance(plot=F)
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    
    proper.gene.names <- rownames(counts)
    sum_cnt <- sapply(gSets, function(sn) sum(proper.gene.names %in% sn))
    rm_idx <- which(sum_cnt <= min_gset_size)
    if(length(rm_idx)>0) {
      warning(paste("Removing genesets for being too small: ", names(gSets)[rm_idx], sep="\n"))
      gSets <- gSets[-rm_idx]
    }
    
    path_names = c()
    env <- list2env(gSets)
    
    # p2$calculatePcaReduction(nPcs=50,n.odgenes=3e3)
    # p2$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
    # p2$getKnnClusters(method=infomap.community,type='PCA')
    # p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='community',z.threshold=3)
    
    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = F,
                                 min.pathway.size =1)
    
    path_names = names(p2$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(counts))
    rownames(score) = path_names
    colnames(score) = colnames(counts)
    for(i in 1:length(p2$misc$pwpca)){
      if(!is.null(p2$misc$pwpca[[i]]$xp$score)){
        score[i,] = as.numeric(p2$misc$pwpca[[i]]$xp$scores)
      }
    }
    
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

ln_or_tumor <- 'LN'
DefaultAssay(seul[[ln_or_tumor]]) <- 'RNA'
expr <- GetAssayData(seul[[ln_or_tumor]], slot='counts')
expr <- expr[rowSums(expr)>=50, ]
for(id in names(msig_l)){
  print(paste0(id, "..."))
  dir.create(file.path(outdir_pway, "pathways", "pagoda2"), recursive = T, showWarnings = F)
  pagoda2_scores <- cal_pagoda2(expr, msig_l[[id]])
  saveRDS(pagoda2_scores, file=file.path(outdir_pway, "pathways", "pagoda2", 
                                         paste0("pagoda2.", id, ".", cellgrp, ".rds")))
}




## 
pagoda2_scores_l <- lapply(names(msig_l), function(id){
  readRDS(file=file.path(outdir_pway, "pathways", "pagoda2", 
                         paste0("pagoda2.", id, ".rds")))
})
names(pagoda2_scores_l) <- names(msig_l)

seu <- seul[[ln_or_tumor]]
ids <- c('nf-k', 'cell proliferation')
sample_ids <- as.character(unique(seu$orig.ident))
comparisons <- list('KO.D7.trx_vs_un'=c('B2_LN_KO_7d', 'B2_LN_KO_Un'),
                    'WT.D7.trx_vs_un'=c('B2_LN_WT_7d', 'B2_LN_WT_Un'),
                    'KO.D3.trx_vs_un'=c('B1_LN_KO_3d', 'B1_LN_KO_Un'),
                    'WT.D3.trx_vs_un'=c('B1_LN_WT_3d', 'B1_LN_WT_Un'),
                    'D7.trx.wt_vs_ko'=c('B2_LN_WT_7d', 'B2_LN_KO_7d'),
                    'D3.trx.wt_vs_ko'=c('B1_LN_WT_3d', 'B1_LN_KO_3d'))
assay <- 'GO'
per_annotation <- TRUE
per_treg  <- TRUE

# Append the pagoda2 assays
Idents(seu) <- 'orig.ident'
assay_scores <- pagoda2_scores_l[[assay]]
assay_scores[is.na(assay_scores)] <- 0
table(colnames(assay_scores) %in% Cells(seu))
seu[[assay]] <- CreateAssayObject(assay_scores)
DefaultAssay(seu) <- assay

Idents(seu) <- 'manual_anno'
seu$treg <- ifelse(grepl("treg", seu$manual_anno, ignore.case=T), 'Treg', 'Non-treg')
seutreg <- SplitObject(seu, split.by='treg')[['Treg']]

Idents(seu) <- 'manual_anno'
annoids <- unique(as.character(seu$manual_anno))
seusub <- SplitObject(seu, split.by='manual_anno')
Idents(seu) <- 'orig.ident'

de_pas <- lapply(names(comparisons), function(compid){
  idents <- comparisons[[compid]]
  
  .runFindMarkers <- function(seu, idents, id_i, compid, base.ident='orig.ident'){
    feature_i <- grep(id_i, rownames(seu), value=T, ignore.case=T)
    tryCatch({
      Idents(seu) <- base.ident
      FindMarkers(seu, test.use='wilcox', features=feature_i, slot='counts',
                  logfc.threshold=-1, min.pct=0, 
                  ident.1= idents[1], ident.2= idents[2],
                  verbose = T) %>% 
        tibble::rownames_to_column(.) %>%
        mutate(ident.1=idents[1],
               ident.2=idents[2],
               comparison=compid,
               search_pattern=id_i)
    }, error=function(e){
      as.data.frame(matrix(nrow=0, ncol = 10)) %>%
        magrittr::set_colnames(c('rowname', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2',
                                 'p_val_adj', 'ident.1', 'ident.2', 'comparison', 'search_pattern'))
    })
  }
  iterateFindMarkers <- function(seu, ids, compid){
    de_pa <- lapply(ids, function(id_i){
      print(paste0(compid, ": ", id_i, "..."))
      .runFindMarkers(seu, idents, id_i, compid)
    }) %>% 
      do.call(rbind, .) %>% as.data.frame %>%
      arrange(p_val_adj)
    return(de_pa)
  }
  
  de_pa_l <- list()
  de_pa_l[['all']] <- iterateFindMarkers(seu, ids, compid)
  
  if(per_treg){
    de_pa_l[['treg']] <- iterateFindMarkers(seutreg, ids, compid)
  }
  
  if(per_annotation){
   de_pa_annos_l <- lapply(seusub, function(seu_i){
      print(paste0("Running: ", unique(seu_i$manual_anno), "..."))
      Idents(seu_i) <- 'orig.ident'
      tryCatch({
        iterateFindMarkers(seu_i, ids, compid)
      }, error=function(e){NA})
    })
    de_pa_l[['clusters']] <- de_pa_annos_l
  }
  
  return(de_pa_l)
})
names(de_pas) <- names(comparisons)

saveRDS(de_pas, file=file.path(outdir_pway, "pathways", "pagoda2", 
                               paste0("de_wilcoxon.", assay, ".rds")))

de_pas <- readRDS(file=file.path(outdir_pway, "pathways", "pagoda2", 
                       paste0("de_wilcoxon.", assay, ".rds")))

for(id in names(de_pas)){
  write.table(de_pas[[id]]$treg, file=paste0("~/xfer/nfkb_cellprol.", id, ".csv"),
              sep=",", col.names = T, row.names = F, quote = F)
  cat(paste0("xfer nfkb_cellprol.", id, ".csv\n"))
}
# de_pas_sig <- lapply(de_pas, function(de){
#   de %>% 
#     filter(p_val_adj< 0.05)
# })

#--- a) AUCell pathway activity score ####
# Clean up the labels in the seurat object
seul <- lapply(seul, function(seu){
  seu$orig.ident <- .relabelid(seu$orig.ident)
  # seu$manual_anno <- seu$seurat_clusters %>% 
  #   recode_map()
  return(seu)
})


seul_auc <- lapply(names(seul), function(group){
  print(paste0("AUCell analysis of: ", group, "..."))
  seu <- seul[[group]]
  DefaultAssay(seu) <- 'RNA'
  expr <- GetAssayData(seu, slot='counts')
  idx <- which(rowSums(expr==0) < (ncol(seu)*0.95))
  # idx <- which(rowSums(expr==0) < 20000)
  
  auc_ls <- lapply(names(msig_lvls), function(mlvl){
    print(paste0("  - Msigdb lvl: ", mlvl, "..."))
    outf <- file.path(outdir, "aucell", paste0("msigdb_", mlvl, ".", group, ".rds"))
    if(!file.exists(outf)){
      auc_l <- iterateMsigdb(species='Mus musculus', msig_lvls=msig_lvls[mlvl], 
                             fun=aucellFun, expr_mat=expr[idx,], gm=gm,
                             mapfrom='ENTREZID', mapto='SYMBOL') %>% 
        unlist(., recursive=F)
      saveRDS(auc_l, file=outf)
    } else {
      auc_l <- readRDS(outf)
    }
    return(auc_l)
  })
  return(auc_ls)
})




auc_ij_all <- lapply(auc_l, function(auc_i){
  as.data.frame(assay(auc_i))
}) %>% 
  do.call(rbind, .) %>%
  as.data.frame
rownames(auc_ij_all) <- gsub("^.*\\.", "", rownames(auc_ij_all))
auc_ij <- auc_ij_all[top_gs,]


#--- b) SCPA differential pathways ####
# Clean up the labels in the seurat object
seul <- lapply(seul, function(seu){
  seu$orig.ident <- .relabelid(seu$orig.ident)
  # seu$manual_anno <- seu$seurat_clusters %>% 
  #   recode_map()
  return(seu)
})


species <- 'Mus musculus'
pathways <- lapply(names(msig_lvls), function(mlvl){
  lapply(msig_lvls[[mlvl]], function(sublvl){
    msigdbr(species = species, category = mlvl, subcategory=sublvl) %>%
      format_pathways()
  })
}) %>%  unlist(., recursive=F) %>% do.call(c, .)

seu <- seul[[group]]
grp_comparisons <- switch(group,
                         LN=list(c("5", "6", "7", "4", "0", "9", "13"),
                                 c('15', '20', '19')),
                         Tumor=list(c('4', '1', '7', '0', '5', '11', '3'),
                                    c('15', '10'),
                                    c('6', '18')))
scpa_grps_all <- lapply(grp_comparisons, function(grp_comparison){
  id <- paste(grp_comparison, collapse="_")
  print(paste0("Group: ", group, " - ", id))
  outf <- file.path(outdir, "scpa", paste0("SCPA_", id, ".", group, ".rds"))
  if(file.exists(outf)){
    scpa_grps <- readRDS(file=outf)
  } else {
    grp_combn <- combn(grp_comparison, m=2)
    scpa_grps <- apply(grp_combn, 2, function(grp_i){
      print(paste(grp_i, collapse="_"))
      scpa_out <- compare_seurat(seu,
                                 group1 = "manual_clusters", 
                                 group1_population = grp_i,
                                 pathways = pathways, parallel=T, cores=4)
      return(scpa_out)
    })
    saveRDS(scpa_grps, file=outf)
  }
  return(scpa_grps)
})
# grp_comparison <- c('13', '9')
# grp_i <- grp_comparison



scpa_out <- scpa_out %>%
  mutate(Dataset=gsub("_.*", "", Pathway))
scpa_out %>% filter(Dataset=='REACTOME') %>% filter(qval > quantile(qval, 0.95))
saveRDS(scpa_out, file=file.path(outdir, "scpa", paste0("scpa_tregs.",group, ".rds")))



auc_ij_all <- lapply(auc_l, function(auc_i){
  as.data.frame(assay(auc_i))
}) %>% 
  do.call(rbind, .) %>%
  as.data.frame
rownames(auc_ij_all) <- gsub("^.*\\.", "", rownames(auc_ij_all))
auc_ij <- auc_ij_all[top_gs,]

###########################################################
#### 6. Re-analyze Day7 ST2 independently from scratch ####
#--- a) Preprocess ----
groups_7d <- grep("ST2*", groups, value=T)
seus <- lapply(groups_7d, function(grp){
  print(paste0(grp, "..."))
  mtx <- Read10X(data.dir = file.path(datadir, grp), strip.suffix=TRUE)
  seu <- CreateSeuratObject(counts = mtx, project = grp)
  return(seu)
})
names(seus) <- groups_7d

## Merge the scRNA data and save
seu <- merge(x=seus[[1]], y = seus[-1], 
             add.cell.ids = names(seus), 
             project = 'st2_il33_7d')
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu.7dst2.rds"))
seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu.7dst2.rds"))
rm(seus)

## QC removal
seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mt")

mt_cutoff <- 10
# table(seu$percent.mt <= mt_cutoff, seu$orig.ident)
seu_mt <- seu[,which(seu$percent.mt <= mt_cutoff)]
# Remove outlier samples for count and features
seu_qc <- subset(seu_mt,
                 subset =nCount_RNA > 1000 &
                   nFeature_RNA < 6686)

seu_qcsumm <- cbind(table(seu$orig.ident),
                    table(seu_mt$orig.ident),
                    table(seu_qc$orig.ident)) %>%
  as.data.frame() %>%
  rename_with(., ~ c("original", "low.mt", "count.feature")) %>%
  mutate("low.mt.frac"=round(low.mt / original, 2),
         "count.frac"=round(count.feature / original, 2))
seu <- seu_qc
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_filt.7dst2.rds"))
rm(seu)

SeuratData:::InstalledData()$Dataset
azimuthpath <- function(x){
  path <- "/cluster/home/quever/downloads/renvs/renv/library/R-4.2/x86_64-pc-linux-gnu/XXXX.SeuratData/azimuth"
  gsub("XXXX", x, path)
}
seuazi <- RunAzimuth(seu, reference = azimuthpath("pbmcref"))
# Azimuth:::RunAzimuth.Seurat
# SeuratData:::AvailableData()$Dataset
# SeuratData:::InstallData('pbmc3k')

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
DefaultAssay(seu) <- 'RNA'

seu <- NormalizeData(seu) %>%
  FindVariableFeatures(.)  %>%
  ScaleData(.)  %>%
  RunPCA(.)  %>%
  FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
  FindClusters(., resolution = 2, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
seu <- SCTransform(seu, vst.flavor = "v2") %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_filt.7dst2.rds"))
seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu_filt.7dst2.rds"))

DefaultAssay(seu) <- 'RNA'
seu_integ <- IntegrateLayers(
  object = seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
seu_integ <- JoinLayers(seu_integ)

seu_integ <- FindNeighbors(seu_integ, reduction = "integrated.mnn", dims = 1:30) %>%
  FindClusters(., resolution = 2, cluster.name = "mnn_clusters") %>%
  RunUMAP(., reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
saveRDS(seu_integ, file=file.path(datadir, "seurat_obj", "seu_integ.7dst2.rds"))

p1 <- DimPlot(
  seu,
  reduction = "umap.unintegrated",
  group.by = c("orig.ident", "unintegrated_clusters"),
  combine = FALSE
)
p2 <- DimPlot(
  seu_integ,
  reduction = "umap.mnn",
  group.by = c("orig.ident", "mnn_clusters"),
  combine = FALSE
)

pdf("~/xfer/dimplot.st2_d7.pdf", width = 9, height = 9)
p1
p2
dev.off()
#--- b) Compare DEGs ----
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(datadir, "seurat_obj", "seu_integ.7dst2.rds"))
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))

ids <- c(setNames(seul$LN$manual_anno, colnames(seul$LN)),
         setNames(seul$Tumor$manual_anno, colnames(seul$Tumor)))

seu_integ$ids <- paste(seu_integ$orig.ident, gsub("_[0-9]*$", "", colnames(seu_integ)), sep="_")
seu_integ@meta.data[,'manual_anno'] <- ids[seu_integ$ids]
seu_integ$ln_tumor <- gsub("^.*(LN|Tumor).*$", "\\1", seu_integ$orig.ident)
seu_integ$orig.ident2 <- .relabelid(seu_integ$orig.ident)
rm(seul); gc()
seul <- SplitObject(seu_integ, split.by='ln_tumor')

if(visualize){
  
  pdf("~/xfer/dimplot.st2_d7.2.pdf", width = 20, height = 9)
  DimPlot(
    seu_integ,
    reduction = "umap.mnn",
    group.by = c("orig.ident", "mnn_clusters", "manual_anno"),
    combine = TRUE,
    raster=T
  )
  dev.off()
}

comparisons <- list('LN'=list(
  'WT_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='WT', "order"='7d'),
  'KO_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='KO', "order"='7d'),
  'DAY7_ko_vs_wt'=list("day"='7d', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO'),
  'B2UN_ko_vs_wt'=list("day"='Un', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO')),
  'Tumor'=list(
    'WT_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='WT', "order"='7d'),
    'KO_day7_vs_un'=list("day"=c('7d', 'Un'), "batch"='B2', "kowt"='KO', "order"='7d'),
    'DAY7_ko_vs_wt'=list("day"='7d', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO'),
    'B2UN_ko_vs_wt'=list("day"='Un', "batch"='B2', "kowt"=c('KO', 'WT'), "order"='KO')))

compare_grp <- list("LN_day7_un"=c("day7_kowt", "unb2_kowt", "_.*$"),
                    "Tumor_day7_un"=c("day7_kowt", "unb2_kowt", "_.*$"))
anno_ident <- 'manual_anno' #'manual_anno', 'manual_clusters'

overwrite <- FALSE
outdirdeg <- file.path(outdir, "degs")
degf <- file.path(outdirdeg, paste0(anno_ident, "_degs.7dreprocess.rds")) #"_deg2s.rds"))
if(any(!file.exists(degf) | overwrite)){
  print(paste0("Creating new DEG file for: ", anno_ident))
  options(future.globals.maxSize= 891289600)
  ct_markers_all <- lapply(names(seul)[1], function(seu_id){
    print(paste0("> ", seu_id))
    seu_i <- seul[[seu_id]]
    seu_i$orig.ident <- seu_i$orig.ident2
    seu_i$day <- gsub("^.*_", "", seu_i$orig.ident)
    seu_i$kowt <- gsub("^.*_(.*)_[a-zA-Z0-9]*$", "\\1", seu_i$orig.ident)
    seu_i$batch <- gsub("_.*", "", seu_i$orig.ident)
    
    Idents(seu_i) <- 'orig.ident'
    ct_markers <- lapply(names(comparisons[[seu_id]]), function(comp_id){
      comp_i <- comparisons[[seu_id]][[comp_id]]
      print(paste0(">> ", comp_id))
      print(paste0("  -- ", paste(unlist(comp_i), collapse="-")))
      
      subset_of_seu <- TRUE
      idents <- with(seu_i@meta.data, orig.ident[(day %in% comp_i$day) & 
                                                   (kowt %in% comp_i$kowt) & 
                                                   (batch == comp_i$batch)]) %>%
        unique
      idents <- idents[order(grepl(comp_i$order, idents))]
      
      # DEG across all cells (not cluster specific)
      markers_all <- FindMarkers(seu_i, assay = "RNA", test.use='wilcox',
                                 ident.1= idents[1], ident.2= idents[2],
                                 verbose = FALSE,
                                 logfc.threshold = 0,
                                 recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
        tibble::rownames_to_column(., "symbol") %>% 
        mutate(biotype = sym2biotype_ids[symbol],
               ensemble = gm$SYMBOL$ENSEMBL[symbol],
               ident.1 = idents[1],
               ident.2 = idents[2], 
               cluster='all_cells',
               comparison=comp_id,
               tissue=seu_id)
      
      # DEG For each cluster
      Idents(seu_i) <- anno_ident
      cluster_ids <- as.character(na.omit(unique(Idents(seu_i))))
      markers_clusters <- lapply(setNames(cluster_ids,cluster_ids), function(cl_j){
        print(paste0("  -- ", cl_j))
        seu_j <- subset(seu_i, ident=cl_j)
        Idents(seu_j) <- 'orig.ident'
        tryCatch({
          FindMarkers(seu_j, assay = "RNA", test.use='wilcox',
                      ident.1= idents[1], ident.2= idents[2],
                      verbose = FALSE,
                      logfc.threshold = 0,
                      recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
            tibble::rownames_to_column(., "symbol") %>% 
            mutate(biotype = sym2biotype_ids[symbol],
                   ensemble = gm$SYMBOL$ENSEMBL[symbol],
                   ident.1 = idents[1],
                   ident.2 = idents[2], 
                   cluster=cl_j,
                   comparison=comp_id,
                   tissue=seu_id)
        }, error=function(e){
          print("     - No differential genes")
          NULL
        })
      })
      
      # Aggregate and return
      markers <- c(list("All"=markers_all),
                   markers_clusters)
      return(markers)
    })
    
    return(ct_markers)
  })
  ct_markers_all <- readRDS(file=degf)
  
  
  names(ct_markers_all) <- names(seul)[1]
  names(ct_markers_all$LN) <- names(comparisons$LN)
  # names(ct_markers_all$Tumor) <- names(comparisons$Tumor)
  saveRDS(ct_markers_all, file=degf)
} else {
  print(paste0("Reading in existing DEG file for: ", anno_ident))
  ct_markers_all <- readRDS(degf)
}

# Find all cluster names (e.g. 1,2,3 or Tregs, B, NKT)
all_cts <- sapply(unlist(ct_markers_all, recursive = F), names) %>%
  unlist %>% as.character %>% unique
# all_cts <- all_cts[grep("TReg", all_cts)]

## Write out the DEG files
# ct_markers[[LN/Tumor]][[names(comparisons)]][[celltype]]

lapply(names(ct_markers_all), function(seu_id){
  deg <- ct_markers_all[[seu_id]]
  lapply(names(deg), function(comp_id){
    deg_comp_i <- deg[[comp_id]]
    ct_ids <- names(deg_comp_i)
    merge_ids <- c('symbol', 'ensemble', 'biotype', 'ident.1', 'ident.2', 'comparison', 'tissue')
    deg_merge <- lapply(names(deg_comp_i), function(deg_id){
      if(!is.null(deg_comp_i[[deg_id]])){
        deg_comp_i[[deg_id]] %>%
          select(merge_ids, avg_log2FC, p_val, p_val_adj) %>%
          rename_with(., ~paste0(deg_id, ".", .), .cols=c(avg_log2FC, p_val, p_val_adj))
      } else {
        as.data.frame(matrix(rep("NA", length(merge_ids)), nrow=1)) %>%
          magrittr::set_colnames(merge_ids)
      }
    }) %>%
      purrr::reduce(., left_join, by=merge_ids)
    write.table(deg_merge, 
                file=file.path(outdirdeg, paste0("deg_st2scratch.", seu_id, ".", comp_id, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  })
})
###############################################
#### 6. Change of KO-treatment across UMAP ####
# The task of this section is to visualize local clusters (in umap) where
# there is a KO-specific change in treatment. To disambiguate this from
# WT-changes in treatment, we are first evaluating what are the changes
# in treatment in WT condition. Using those genes, I can score the 
# activity of the up-regulated and down-regulated genes in the WT condition
# and compare it to the KO-condition. The final will be a log-ratio between
# the KO and WT scores

if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.rds"))
pathoutdir <- file.path(outdir, "pathway")
dir.create(pathoutdir, recursive = F)

seul <- lapply(seul, recode_list, newcol='manual_clusters', grp='orig.ident', anno=F)
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)
# Clean up the labels in the seurat object
seul <- lapply(seul, function(seu){
  seu$orig.ident <- .relabelid(seu$orig.ident)
  # seu$manual_anno <- seu$seurat_clusters %>% 
  #   recode_map()
  return(seu)
})

anno_ident <- 'manual_clusters' #'manual_anno', 'manual_clusters'
pval_sig <- 0.05
fc_sig <- 0.5

anno_ident <- 'manual_clusters' #'manual_anno', 'manual_clusters'
degf <- file.path(outdir, "degs", paste0(anno_ident, "_degs.rds")) #"_deg2s.rds"))
ct_markers_all <- readRDS(degf)



mode <- 'All' # 'All', 'clusters', 'local'
reduction='umap'
day7_3_list <- list("Day3"=c("ko"="day3_KO", "wt"="day3_WT"),
                    "Day7"=c("ko"="day7_KO", "wt"="day7_WT"))
k <- 100
cow_fps <- lapply(names(seul), function(seuid){
  print(seuid)
  seu_i <- seul[[seuid]]
  expr <- GetAssayData(seu_i, slot='counts')
  idx <- which(rowSums(expr==0) < (ncol(seu_i)*0.95))
  
  fps <- lapply(day7_3_list, function(day37){
    if(mode=='clusters'){
      grps <- as.character(unique(Idents(seu_i)))
      cellidx <- lapply(grps, function(grp_i){
        which(Idents(seu_i) == grp_i)
      }) %>% setNames(., grps)
    } else if(mode=='All'){
      grps <- list("All")
      cellidx <- list("All"=seq_along(Cells(seu_i)))
    } else if(mode=='local'){
      #XXXX
    }
    
    scores_l <- lapply(grps, function(grp){
      ko_ids <- ct_markers_all[[seuid]][[day37['ko']]][[grp]] %>%
        select(ident.1, ident.2) %>% unique %>% unlist
      wt_ids <- ct_markers_all[[seuid]][[day37['wt']]][[grp]] %>%
        select(ident.1, ident.2) %>% unique %>% unlist
      all_ids <- setNames(c(ko_ids, wt_ids), c("KO_Trx", "KO_Un", "WT_Trx", "WT_Un"))
      
      top_markers <- ct_markers_all[[seuid]][[day37['wt']]][[grp]] %>% 
        filter((abs(avg_log2FC) > fc_sig) & (p_val_adj < pval_sig))
      msig_l <- list("up"=top_markers %>% filter(avg_log2FC > 0) %>% pull(symbol),
                     "down"=top_markers %>% filter(avg_log2FC < 0) %>% pull(symbol))
      
      #--- ssGSEA from GSVA package
      cellidx_grp <- cellidx[[grp]]
      ssgsea_score = GSVA::gsva(expr = expr[idx,cellidx_grp], msig_l, 
                                method = "ssgsea", 
                                parallel.sz = 1, verbose = T)
      
      if(any(sapply(msig_l, length)==0)){
        null_idx <- which(sapply(msig_l, length)==0)
        null_mat <- matrix(ncol=ncol(ssgsea_score), 
                           nrow=length(null_idx), data=0)
        rownames(null_mat) <- names(null_idx)
        ssgsea_score <- rbind(ssgsea_score, 
                              null_mat)
        ssgsea_score <- ssgsea_score[names(msig_l),]
      }
      
      #--- nearest neighbours for each cell
      ## Calculate the nearest neighbour for k for the latent-space reductions
      labels <- seu_i@meta.data[cellidx_grp,'orig.ident']
      lspace_obj <- seu_i@reductions[[reduction]]@cell.embeddings[cellidx_grp,]
      nn_obj <- RANN::nn2(data = lspace_obj, query = lspace_obj, 
                          k=k, treetype='kd')
      # Index of every cell of interest
      ko_treated_idx <- which(labels == all_ids['KO_Trx']) # e.g. All KO_day3trx cells
      ko_trx_cells <- nn_obj$nn.idx[ko_treated_idx,]
      all_comp_scores <- apply(ko_trx_cells, 1, function(ko_trx_cell_i){
        # Extract the local ssGSEA scores for each cell comparison
        # should be one local cells that are: KO_Un, WT_Trx, and WT_Un
        local_ssgsea_scores <- sapply(all_ids[-1], function(id_x){
          # Find the KO/WT Trx/Un nearby cells and extract their ssgsea score
          xcell_idx <- which(labels[ko_trx_cell_i] == id_x)
          xcell_cells <- Cells(seu_i)[cellidx_grp][ko_trx_cell_i[xcell_idx]]
          xcell_score <- apply(ssgsea_score[,xcell_cells,drop=F], 1, mean)
          return(xcell_score)
        })
        
        return(as.data.frame(local_ssgsea_scores))
      })
      
      # For every KO_treated cell, find its comparative ssGSEA score comparing
      # it tot he KO_un cell, and compare that to the local WT-Trx/Un score
      ko_wt_ratios <- sapply(seq_along(ko_treated_idx), function(ko_idx_i){
        # Get the KO_Trx ssgsea score
        ko_trx_cell <- Cells(seu_i)[cellidx_grp][ko_treated_idx[ko_idx_i]]
        ko_trx_score <- ssgsea_score[,ko_trx_cell]
        
        # Compare the KO_Trx to the KO_Un;  compare that FC to WT_Trx/Un FC
        comp_score <- all_comp_scores[[ko_idx_i]]
        ko_trx_vs_un <- ko_trx_score / comp_score[,'KO_Un']
        wt_trx_vs_un <- comp_score[,'WT_Trx'] / comp_score[,'WT_Un']
        ko_wt_ratio <- (ko_trx_vs_un - wt_trx_vs_un)
        return(ko_wt_ratio)
      }) %>% t 
      
      # Flip direction of down ssGSEA score and add together up- and down-regulated genes
      ko_wt_ratios[,'down'] <- -1*ko_wt_ratios[,'down']
      ko_wt_scores <- rowSums(ko_wt_ratios, na.rm=T)
      
      quantiles <- quantile(ko_wt_scores, c(0.015, 0.99), na.rm=T)
      ko_wt_scores[ko_wt_scores>quantiles[2]] <- quantiles[2]
      ko_wt_scores[ko_wt_scores<quantiles[1]] <- quantiles[1]
      
      ko_wt_score <- rep(NA, ncol(seu_i))
      ko_wt_score[cellidx_grp][ko_treated_idx] <- ko_wt_scores
      return(ko_wt_score)
    })
    
    seu_i$ko_wt_score <-  do.call(cbind, scores_l) %>% rowSums(., na.rm=T)
    Idents(seu_i) <- 'orig.ident'
    kotrx_ids <- grep("KO_[3d|7d]", seu_i$orig.ident, value=T) %>% unique
    seu_j <- subset(seu_i, ident=kotrx_ids)
    
    max_score <- max(abs(seu_i$ko_wt_score))
    b <- c(-1*max_score, 0, max_score)
    fp <- FeaturePlot(seu_j, features = 'ko_wt_score', reduction='umap', raster = F) + 
      ggtitle(paste0(seuid, ": ", day37[1])) +
      scale_color_gradientn(limits = c(min(b),max(b)),
                            colours=c("black", "blue", "grey", "red", "yellow"),
                            breaks=b, labels=format(b))
    return(fp)
  })
  
  cow_fp <- cowplot::plot_grid(plotlist=fps, nrow=1)
  return(cow_fp)
})
saveRDS(cow_fps, file="~/xfer/cow_fps.rds")
pdf("~/xfer/KOtrx_vs_WTtrx_signature.pdf", width = 9, height = 6)
cow_fps
dev.off()

# ---- a) Change in TRegs ----
seul_tregs <- readRDS(file=file.path(datadir, "seurat_obj", "tregs_final.rds"))
pathoutdir <- file.path(outdir, "ko_change")
dir.create(pathoutdir, recursive = F)

anno_ident <- 'treg_manual_anno' #'manual_anno', 'manual_clusters'
pval_sig <- 0.05
fc_sig <- 0.5

degf <- file.path(outdir, "degs", paste0(anno_ident, "_degs.rds")) #"_deg2s.rds"))
ct_markers_all <- readRDS(degf)



mode <- 'All' # 'All', 'clusters', 'local'
reduction='umap'
day7_3_list <- list("Day3"=c("ko"="day3_KO", "wt"="day3_WT"),
                    "Day7"=c("ko"="day7_KO", "wt"="day7_WT"))
k <- 100
dayid <- 'combined'
cow_fps <- lapply(names(seul), function(seuid){
  print(seuid)
  seu_i <- seul[[seuid]][[dayid]]
  expr <- GetAssayData(seu_i, slot='counts')
  idx <- which(rowSums(expr==0) < (ncol(seu_i)*0.95))
  
  fps <- lapply(day7_3_list, function(day37){
    if(mode=='clusters'){
      grps <- as.character(unique(Idents(seu_i)))
      cellidx <- lapply(grps, function(grp_i){
        which(Idents(seu_i) == grp_i)
      }) %>% setNames(., grps)
    } else if(mode=='All'){
      grps <- list("All")
      cellidx <- list("All"=seq_along(Cells(seu_i)))
    } else if(mode=='local'){
      #XXXX
    }
    
    scores_l <- lapply(grps, function(grp){
      Idents(seu_i) <- 'orig.ident'
      all_ids <- unique(as.character(Idents(seu_i)))
      
      ko_ids <- grep("KO_", all_ids, value=T) %>% sort(.)
      wt_ids <- grep("WT_", all_ids, value=T) %>% sort(.)
      all_ids <- c(ko_ids,wt_ids)
      
      top_markers <- FindMarkers(seu_i, ident.1=wt_ids[1], ident.2=wt_ids[2]) %>%
        as.data.frame %>% 
        filter((p_val_adj < pval_sig) & (abs(avg_log2FC) > fc_sig)) %>%
        tibble::rownames_to_column(., "symbol")
      msig_l <- list("up"=top_markers %>% filter(avg_log2FC > 0) %>% pull(symbol),
                     "down"=top_markers %>% filter(avg_log2FC < 0) %>% pull(symbol))
      
      #--- ssGSEA from GSVA package
      cellidx_grp <- cellidx[[grp]]
      ssgsea_score = GSVA::gsva(expr = expr[idx,cellidx_grp], msig_l, 
                                method = "ssgsea", 
                                parallel.sz = 1, verbose = T)
      
      if(any(sapply(msig_l, length)==0)){
        null_idx <- which(sapply(msig_l, length)==0)
        null_mat <- matrix(ncol=ncol(ssgsea_score), 
                           nrow=length(null_idx), data=0)
        rownames(null_mat) <- names(null_idx)
        ssgsea_score <- rbind(ssgsea_score, 
                              null_mat)
        ssgsea_score <- ssgsea_score[names(msig_l),]
      }
      
      #--- nearest neighbours for each cell
      ## Calculate the nearest neighbour for k for the latent-space reductions
      labels <- seu_i@meta.data[cellidx_grp,'orig.ident']
      lspace_obj <- seu_i@reductions[[reduction]]@cell.embeddings[cellidx_grp,]
      nn_obj <- RANN::nn2(data = lspace_obj, query = lspace_obj, 
                          k=k, treetype='kd')
      # Index of every cell of interest
      KO_Un_id <- grep("KO_Un", all_ids, value=T)
      KO_Trx_id <- grep("KO_[37]d", all_ids, value=T)
      WT_Un_id <- grep("WT_Un", all_ids, value=T)
      WT_Trx_id <- grep("WT_[37]d", all_ids, value=T)
      
      ko_treated_idx <- which(labels == KO_Trx_id) # e.g. All KO_day3trx cells
      ko_trx_cells <- nn_obj$nn.idx[ko_treated_idx,]
      all_comp_scores <- apply(ko_trx_cells, 1, function(ko_trx_cell_i){
        # Extract the local ssGSEA scores for each cell comparison
        # should be one local cells that are: KO_Un, WT_Trx, and WT_Un
        local_ssgsea_scores <- sapply(c(KO_Un_id, WT_Un_id, WT_Trx_id), function(id_x){
          # Find the KO/WT Trx/Un nearby cells and extract their ssgsea score
          xcell_idx <- which(labels[ko_trx_cell_i] == id_x)
          xcell_cells <- Cells(seu_i)[cellidx_grp][ko_trx_cell_i[xcell_idx]]
          xcell_score <- apply(ssgsea_score[,xcell_cells,drop=F], 1, mean)
          return(xcell_score)
        })
        
        return(as.data.frame(local_ssgsea_scores))
      })
      
      # For every KO_treated cell, find its comparative ssGSEA score comparing
      # it tot he KO_un cell, and compare that to the local WT-Trx/Un score
      ko_wt_ratios <- sapply(seq_along(ko_treated_idx), function(ko_idx_i){
        # Get the KO_Trx ssgsea score
        ko_trx_cell <- Cells(seu_i)[cellidx_grp][ko_treated_idx[ko_idx_i]]
        ko_trx_score <- ssgsea_score[,ko_trx_cell]
        
        # Compare the KO_Trx to the KO_Un;  compare that FC to WT_Trx/Un FC
        comp_score <- all_comp_scores[[ko_idx_i]]
        
        ko_trx_vs_un <- ko_trx_score / comp_score[,KO_Un_id]
        wt_trx_vs_un <- comp_score[,WT_Trx_id] / comp_score[,WT_Un_id]
        ko_wt_ratio <- (ko_trx_vs_un - wt_trx_vs_un)
        return(ko_wt_ratio)
      }) %>% t 
      
      # Flip direction of down ssGSEA score and add together up- and down-regulated genes
      ko_wt_ratios[,'down'] <- -1*ko_wt_ratios[,'down']
      ko_wt_scores <- rowSums(ko_wt_ratios, na.rm=T)
      
      quantiles <- quantile(ko_wt_scores, c(0.015, 0.99), na.rm=T)
      ko_wt_scores[ko_wt_scores>quantiles[2]] <- quantiles[2]
      ko_wt_scores[ko_wt_scores<quantiles[1]] <- quantiles[1]
      
      ko_wt_score <- rep(NA, ncol(seu_i))
      ko_wt_score[cellidx_grp][ko_treated_idx] <- ko_wt_scores
      return(ko_wt_score)
    })
    
    seu_i$ko_wt_score <-  do.call(cbind, scores_l) %>% rowSums(., na.rm=T)
    Idents(seu_i) <- 'orig.ident'
    kotrx_ids <- grep("KO_[3d|7d]", seu_i$orig.ident, value=T) %>% unique
    seu_j <- subset(seu_i, ident=kotrx_ids)
    
    max_score <- max(abs(seu_i$ko_wt_score))
    b <- c(-1*max_score, 0, max_score)
    fp <- FeaturePlot(seu_j, features = 'ko_wt_score', reduction='pacmap', raster = F) + 
      ggtitle(paste0(seuid, ": ", day37[1])) +
      scale_color_gradientn(limits = c(min(b),max(b)),
                            colours=c("black", "blue", "grey", "red", "yellow"),
                            breaks=b, labels=format(b))
    return(fp)
  })
  
  cow_fp <- cowplot::plot_grid(plotlist=fps, nrow=1)
  return(cow_fp)
})
saveRDS(cow_fps, file="~/xfer/cow_fps.rds")
pdf(paste0("~/xfer/KOtrx_vs_WTtrx_signature.", dayid, ".pdf"), width = 9, height = 6)
cow_fps
dev.off()
########################################################
#### 7. Mapping the change in celltype proportions ####
seurds <- file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds")
if(!exists("seul"))  seul <- readRDS(file = seurds)

ident_meta <- lapply(seul, function(seu){
  data.frame(sampleID=unique(seu$orig.ident)) %>% 
    mutate(Condition=gsub("^.*_", "", sampleID))
})

sampleid <- 'orig.ident'
clusid <- 'manual_anno'
partid <- 'manual_anno'
annoid <- 'manual_anno'

#---- a) Cell type proportions per cluster ----
lapply(names(seul), function(seuid){
  seu_i <- seul[[seuid]]
  sample_cnts <- table(seu_i@meta.data[,sampleid])
  sample_cnts <- (1/(sample_cnts / min(sample_cnts)))
  
  table(seu_i@meta.data[,annoid], seu_i@meta.data[,sampleid]) %>%
    apply(., 1, function(i) i * sample_cnts) %>% 
    t %>% as.data.frame %>% 
    write.table(., #file=file.path("~/xfer", paste0(seuid, ".cell_anno_counts.csv")),
                sep=",", quote = F, col.names = F, row.names = F)
})
#---- b) Cluster based approach -----
conditions <- c('Condition')
gg_deltas <- lapply(setNames(conditions,conditions), function(delta){
  delta_title <- switch(delta,
                        Condition="Condition [DAB-DMSO]")
  dtbl_melt <- getDiffProportionsPerClusters(seu, clusid=clusid, 
                                             sampleid=sampleid, delta=delta,
                                             ident_meta, annoid=annoid,
                                             partid=partid) %>%
    select(-c(partitions, annos))
  
  gg <- ggplot(dtbl_melt, aes(x=Delta, y=Cluster, fill=Direction)) +
    geom_bar(stat='identity', position='dodge') +
    # facet_grid(partitions+annos~., scales = 'free_y', 
    #            space = 'free_y') +
    scale_fill_manual(values=c('Pos'='#d6604d', 'Neg'='#4393c3')) +
    theme_classic() +
    ggtitle(delta_title) + xlab("Standardized residuals") +
    geom_vline(xintercept=c(-3,3), linetype='dashed', col='grey') +
    theme(strip.text.y.right = element_text(angle = 0))
  return(list("gg"=gg, 'dtbl'=dtbl_melt))
})

delta_map <- with(gg_deltas$Condition$dtbl, setNames(Delta, Cluster))
seu$delta.DAB_DMSO <- delta_map[seu$seurat_clusters]

pdf("~/xfer/anno_prop_cluster.pdf", width = 10, height = 10)
sapply(gg_deltas, function(i) print(i$gg))
DimPlot(seu, group.by='seurat_clusters', 
        split.by='functional.cluster', ncol = 4)
FeaturePlot(seu, feature='delta.DAB_DMSO') +
  scale_color_gradientn( colours = c('blue', 'grey', 'red'),  limits = c(-26,26))
dev.off()
write.table(gg_deltas$Condition$dtbl, file="~/xfer/anno_prop_cluster.tsv",
            sep=",", col.names = T, row.names = F, quote = F)

saveRDS(seu, file=outrds)


#---- c) Nearest neighbor based approach -----
reduction <- 'mnn'
category_label <- 'orig.ident'

seu <- seul$LN
## Calculate the nearest neighbour for k for the latent-space
# reductions
labels <- factor(seu@meta.data[,category_label])
lspace_obj <- seu@reductions[[reduction]]@cell.embeddings
ks <- seq(5, 50, by=5)
k_wproportions <- lapply(setNames(ks, ks), function(k){
  nn_obj <- RANN::nn2(data = lspace_obj, query = lspace_obj, 
                      k=k, treetype='kd')
  wdist <- round(1-(nn_obj$nn.dist/max(nn_obj$nn.dist)), 2) * 100
  ulbl <- unique(labels)
  wproportions <- sapply(c(1:nrow(lspace_obj)), function(idx){
    wcnt <- split(wdist[idx,-1], labels[nn_obj$nn.idx[idx,-1]]) %>%
      sapply(., sum)
    wcnt / sum(wcnt)
  })
  t(wproportions)
})
w <- (ks/max(ks))
wproportion <- sapply(1:ncol(k_wproportions[[1]]), function(colidx){
  sapply(k_wproportions, function(i){
    i[,colidx]
  }) %>%
    apply(., 1, weighted.mean, w=w)
}) %>%
  magrittr::set_colnames(paste0("lprop_", colnames(k_wproportions[[1]])))

## weigh the 
wproportions
seu$localProportion <- wproportion
seu@meta.data <- cbind(seu@meta.data, wproportion)

pdf("~/xfer/test.pdf", width = 20, height = 15)
ids <- colnames(wproportion)
lapply(split(ids, f=grepl("_B2_", ids)), function(lprop_set_i){
  lapply(lprop_set_i, function(lprop_id){
    FeaturePlot(seu, features=lprop_id, reduction='umap', order=F) +
      scale_color_gradientn( colours = alpha(c("white", "grey","red", "darkred"), 0.5),  limits = c(0,1)) +
      ggtitle(lprop_id)
  }) %>%
    cowplot::plot_grid(plotlist=., nrow=1)
}) %>%
  cowplot::plot_grid(plotlist=., nrow=2)
dev.off()


saveRDS(seu, file = outrds)
####################################################
#### 5. GSEA between 72hr and PBS for WT and KO ####
## Take the seurat object and split them into a tumor and TDLN subset,
# then reprocess them all indepedently, recluster, revisualize, and then
# integrate them together using MNN
outdirdeg <- file.path(outdir, "degs")
outdirgsea <- file.path(outdir, "gseas")
# if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
# seu$manual_anno <- seu$seurat_clusters %>% 
#   recode_map()
# DefaultAssay(seu) <- 'RNA'

# Order is important here, deltaNES will be [2] - [1]
group_comp <- list("LN"=c('LN_WT', 'LN_KO'),
                   "Tumor"=c('Tumor_WT', 'Tumor_KO'))
groups <- unlist(group_comp)
groups <- setNames(groups, groups)

# ---- a) Call GSEA on log2Ratios ----
## Differential expression per cluster (seurat_clusters) or celltype (manual_anno)
if(!exists("ct_markers_all2")) ct_markers_all2 <- readRDS(file.path(outdirdeg, "clusters_degs.rds"))
if(!exists("ct_markers_all")) ct_markers_all <- readRDS(file.path(outdirdeg, "celltypes_degs.rds"))
ct_markers <- list("clusters"=ct_markers_all2,
                   "celltypes"=ct_markers_all)
rm(ct_markers_all2, ct_markers_all); gc()

gsea_markers <- lapply(names(ct_markers), function(ct_marker_id){
  ct_marker_i <- ct_markers[[ct_marker_id]]
  # Iterate through clustering type (e.g. seurat_clusters or celltype)
  lapply(groups, function(group_id){
    ct_grp_ij <- ct_marker_i[[group_id]]
    print(paste0(ct_marker_id, " - ", group_id))
    # Iterate through grouping (e.g. LN_KO 72hr vs )
    
    outf <- file.path(outdirgsea, paste0("gsea.", ct_marker_id, ".", group_id, ".rds"))
    #if((!file.exists(outf)) | (file.info(outf)$size < 100)){
    if(!file.exists(outf)){
      saveRDS(NULL,file=outf)
      gsea_grp_ij <- lapply(ct_grp_ij, function(ct_grp_ijk){
        # Iterate through each celltype/cluster (e.g. seurat_cluster 1)
        msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                          'C2'=list('CP:REACTOME'),             # curated gene sets
                          'C5'=list('GO:BP', 'GO:CC', 'GO:MF'))
        obj <- iterateMsigdb(species='Mus musculus', gseaFun, msig_lvls=msig_lvls, 
                             lfc_v=setNames(ct_grp_ijk$avg_log2FC,
                                            gm$SYMBOL$ENTREZ[ct_grp_ijk$symbol]))
        unlist(obj, recursive=F)
      })
      saveRDS(gsea_grp_ij, file=outf)
    } else {
      print("Reading in existing")
      gsea_grp_ij <- readRDS(file=outf)
    }
    return(gsea_grp_ij)
  })
})
names(gsea_markers) <- names(ct_markers)
saveRDS(gsea_markers, file=file.path(outdirgsea, "gsea.clusters_celltypes.rds"))

# ---- b) Compare NES between KO and WT conditions ----
if(!exists("gsea_markers")) gsea_markers <- readRDS(file=file.path(outdirgsea, "gsea.clusters_celltypes.rds"))

# Calculate the delta between KO and WT NES scores
dgsea <- lapply(gsea_markers, function(gsea_i){
  lapply(group_comp, function(group_j){
    gsea1_ij <- gsea_i[group_j][[1]]
    gsea2_ij <- gsea_i[group_j][[2]]
    clusters <- intersect(names(gsea1_ij), names(gsea2_ij))
    clusters <- setNames(clusters,clusters)
    
    lapply(clusters, function(clid){
      .subset <- function(x){
        do.call(rbind, lapply(x, as.data.frame)) %>%
          as.data.frame %>%
          select(ID, NES)
      }
      gsea12 <- dplyr::left_join(.subset(gsea1_ij[[clid]]), 
                                 .subset(gsea2_ij[[clid]]),   
                                 by='ID', suffix=paste0(".", group_j)) %>%
        as.data.frame %>% 
        mutate(Database=gsub("_.*", "", ID),
               dNES=!!rlang::sym(paste0("NES.", group_j[2])) - !!rlang::sym(paste0("NES.", group_j[1]))) %>%
        relocate(Database, .after=ID)
      # arrange(desc(abs(dNES)))
      return(gsea12)
    })
  })
})
saveRDS(dgsea, file=file.path(outdirgsea, "delta_gsea.rds"))

dgsea$celltypes$Tumor$Tregs %>% arrange(desc(abs(dNES)))

# ---- c) Create Cytoscape plots for the comparison GSEAs ----
if(!exists("gsea_markers")) gsea_markers <- readRDS(file=file.path(outdirgsea, "gsea.clusters_celltypes.rds"))
if(!exists("dgsea")) dgsea <- readRDS(file=file.path(outdirgsea, "delta_gsea.rds"))
dir.create(file.path(outdirgsea, "cytoscape"), showWarnings = F)

group <- 'LN_72h'
gs_map <- .mapGsToExactSource(species="Mus musculus")

databases <- names(gsea_markers$celltypes$LN_KO$all)
celltypeid='LN_KO' # 'Tumor_KO'
db='H.base' # "H.base", "C2.CP:REACTOME", "C5.GO:BP", "C5.GO:CC", "C5.GO:MF"      
celltype='Tregs; Cycling' # names(gsea_markers$celltypes$LN_KO)
x <- lapply(setNames(databases,databases), function(db){
  dat <- as.data.frame(gsea_markers$celltypes[[celltypeid]][[celltype]][[db]]) %>% 
    select(-c(NES, pvalue, p.adjust,qvalues))
  databaseid <- switch(db,
                       "H.base"='HALLMARK',
                       "C2.CP:REACTOME"='REACTOME',
                       "C5.GO:BP"='GOBP',
                       "C5.GO:CC"='GOCC',
                       "C5.GO:MF"='GOMF')
  dnes <- dgsea$celltypes$LN[[celltype]] %>%
    filter(Database == databaseid)
  set.seed(1234)
  ridx1 <- sample(1:nrow(dnes), size=10000, replace = T)
  ridx2 <- sample(1:nrow(dnes), size=10000, replace = T)
  null_dnes <- na.omit(dnes[ridx1,4]-dnes[ridx2,3])
  dat_dnes <- left_join(dat, dnes %>% select(ID, dNES), by='ID') %>% 
    rename_with(., ~'NES', .cols='dNES')
  pval <- sapply(dat_dnes$NES, function(i){
    num <- sum(abs(i) >= abs(null_dnes))
    denom <- length(null_dnes)
    1-((num/denom)/2)
  })
  data.frame(dat_dnes$NES, pval)
  dat_dnes <- dat_dnes %>% 
    mutate(pvalue=pval,
           p.adjust=p.adjust(pval, method='BH'),
           qvalues=p.adjust)
  return(dat_dnes)
})

xdf <- do.call(rbind, x) %>% 
  as.data.frame %>% 
  arrange(p.adjust)



gsea_dat_i <-  gsea2CytoscapeFormat(dat, gs_map)
write.table(gsea_dat_i$up, 
            file=file.path(outdir, "gse", "cytoscape", 
                           paste0("gsea_", group, ".", gsub(" ", "_", celltype), ".up.xls")),
            sep="\t", col.names = T, row.names = F, quote = F)
write.table(gsea_dat_i$down, 
            file=file.path(outdir, "gse", "cytoscape", 
                           paste0("gsea_", group, ".", gsub(" ", "_", celltype), ".down.xls")),
            sep="\t", col.names = T, row.names = F, quote = F)


#########################################################################
#### 6. SCENIC Regulon analysis between 72hr-vs-PBS for Tumor and LN ####
## Take the seurat object and split them into a tumor and TDLN subset,
# then reprocess them all indepedently, recluster, revisualize, and then
# integrate them together using MNN
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.rds"))
outdirregulon <- file.path(outdir, "regulons")
dir.create(outdirregulon, showWarnings = F)
setwd(outdirregulon)

anno_ident <- 'seurat_clusters' #'manual_anno', 'seurat_clusters'
cols=c("KO.down_WT.down" = "#fdb462", 
       "KO.down_WT.up" = "#80b1d3", 
       "KO.up_WT.up" = "#fdb462", 
       "KO.up_WT.down" = "#fb8072")

anno_ident <- 'manual_anno' #'manual_anno', 'manual_clusters'
pval_sig <- 0.05
fc_sig <- 0.5
pseudop <- 1*10^-250

seul <- lapply(seul, recode_list, newcol='manual_clusters', grp='orig.ident', anno=F)
seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)

# Clean up the labels in the seurat object
seul <- lapply(seul, function(seu){
  seu$orig.ident <- .relabelid(seu$orig.ident)
  # seu$manual_anno <- seu$seurat_clusters %>% 
  #   recode_map()
  return(seu)
})


#--- a) Pre-pySCENIC ----
# Pre-filtering of the gene expression matrix
min_pct_expr <- 0.01
overwrite_loom <- FALSE

add_cell_annotation <- function(loom, cellAnnotation){
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}


lapply(names(seul), function(seu_id){
  print(seu_id)
  seu <- seul[[seu_id]]
  loomf <- file.path(outdirregulon, paste0("seu.", seu_id, ".loom"))
  
  # Get expression data and filter for expressed in min_pct_expr cutoff
  exprMat <- GetAssayData(seu, slot="data")
  loci1 <- which(rowSums(exprMat) > (min_pct_expr*ncol(exprMat)))
  exprMat_filt <- as.matrix(exprMat[loci1, ])
  cellInfo <- seu@meta.data
  if(!file.exists(loomf) | overwrite_loom){
    loom <- build_loom(loomf, dgem=exprMat_filt)
    loom <- add_cell_annotation(loom, cellInfo)
    close_loom(loom)
  } else {
    print(paste0("Loom file already exists: ", basename(loomf)))
  }
})




#--- b) Post-pySCENIC ----
loom <- open_loom(file.path(outdirregulon, 'seu_scenic.loom'))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC), 
        file = file.path(outdirregulon, 'seu_regulon.rds'))
#--- c) Comparing regulons to localProportions ----
if(!exists("seu_regulon")) seu_regulon <- readRDS(file = file.path(outdirregulon, 'seu_regulon.rds'))
regauc <- assay(seu_regulon$regulonAUC)
stopifnot(all(colnames(regauc) == Cells(seu))) # check to ensure cells are in same order

# Perform a permutation test to estimate the null dstribution of each regulon
localprop <- seu$localProportion
regulon_cor_null <- lapply(rownames(regauc), function(xid){
  x <- regauc[xid,]
  set.seed(1234)
  rs <- sapply(1:1000, function(i){
    localprop_i <- sample(localprop, size=length(localprop), replace=F)
    cor(x, localprop_i, method='pearson')
  })
  data.frame("r"=rs, "regulon"=xid)
}) 
regulon_cor_null_df <- as.data.frame(do.call(rbind, regulon_cor_null))

# Calculate the correlation and p-value for each regulon and localProportion
regulon_cor <- apply(regauc, 1, function(x){
  corres <- cor.test(x, seu$localProportion, method = 'pearson')
  c('cor'=corres$estimate, 'p'=corres$p.value)
}) %>% t %>% as.data.frame %>%
  rename_with(., ~c('correlation', 'p')) %>% 
  mutate(q=p.adjust(p, method='bonferroni'),
         log10padj=(-1*log10(q)), 
         sig=ifelse((q < 0.001) & (abs(correlation)>0.2), 'sig', 'nonSig')) %>% 
  tibble::rownames_to_column(., "regulon") %>% 
  arrange(., q)


# VOLCANO PLOT
winsor <- quantile(regulon_cor$log10padj, 0.95)
regulon_cor$log10padj[regulon_cor$log10padj>winsor] <-  winsor
pdf("~/xfer/regulon_localProp_cor.pdf", width = 7, height = 5)
ggplot(regulon_cor, aes(x=correlation, y=log10padj, colour=sig)) +
  geom_point(alpha=0.7) +
  ggrepel::geom_text_repel(data=regulon_cor[regulon_cor$sig=='sig',], 
                           aes(label=regulon), size=3, max.overlaps=40,
                           min.segment.length = unit(1, "lines"), max.iter = Inf,) +
  scale_color_manual(values=c('sig'='darkred', 'nonSig'='darkgrey')) +
  theme_bw() +
  xlim(-0.5, 0.5) + ylim(0, 200) +
  NoLegend()
dev.off()

# Dimension Plot of localProportion and significant regulons
fp1 <- FeaturePlot(seu, features='localProportion', reduction='umap', order=F, raster=T) +
  scale_color_gradientn( colours = alpha(c("blue", "grey", "white", "grey","red"), 0.5),  limits = c(0,1)) +
  ggtitle("nearest neighbor-based")
sig_regulons <- regulon_cor %>% 
  filter(sig=='sig') %>% arrange(correlation)
regulon_fps <- apply(sig_regulons, 1, function(reg_i){
  seu$regulon <- regauc[reg_i['regulon'],]
  FeaturePlot(seu, features='regulon', reduction='umap', order=F, raster=T) +
    scale_color_gradientn( colours = alpha(c("white", "red","black"), 1),  
                           limits = c(0,0.15)) +
    ggtitle('') + NoLegend() +
    xlab(reg_i['regulon']) + ylab("") + theme(axis.text = element_blank())
})
fp2 <- cowplot::plot_grid(plotlist=regulon_fps, ncol=3)
pdf("~/xfer/top_regulon_dimplot.pdf", width = 13)
plot_grid(fp1, fp2)
dev.off()

# Plot significant regulons in context of their permuted NULL distribution
null_df <- regulon_cor_null_df %>% 
  filter(regulon %in% sig_regulons$regulon) %>%
  rename_with(., ~c('correlation', 'regulon'))
pdf("~/xfer/regulon_null_corr.pdf")
ggplot(null_df, aes(x=correlation)) + 
  # geom_histogram(aes(y=..density..), bins = 200, alpha=0.5, colour='blue', fill='blue') +
  geom_density(alpha=0.4, fill='cyan') +
  geom_point(data=sig_regulons, aes(x=correlation), y=0, color='red') +
  facet_wrap(vars(regulon), scales = 'free') +
  xlim(-0.4, 0.4) +
  theme_bw()
dev.off()


#--- d) Differential regulon per cluster ----
if(!exists("seu_regulon")) seu_regulon <- readRDS(file = file.path(outdirregulon, 'seu_regulon.rds'))
AUCmat <- AUCell::getAUC(seu_regulon$regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))
stopifnot(all(colnames(AUCmat) == Cells(seu))) # check to ensure cells are in same order

# Add the Regulon AUC data to an assay and scale
seu[['AUC']] <- CreateAssayObject(data = AUCmat)
seu <- ScaleData(seu, assay = 'AUC', features = rownames(AUCmat))

# Params
DefaultAssay(seu) <- 'AUC'
Idents(seu) <- 'orig.ident'
splitids <- list("LN_WT"=c("LN_WT_72h", "LN_WT_PBS"),
                 "LN_KO"=c("LN_KO_72h", "LN_KO_PBS"),
                 "Tumor_WT"=c("Tumor_WT_72h", "Tumor_WT_PBS"),
                 "Tumor_KO"=c("Tumor_KO_72h", "Tumor_KO_PBS"))
comp_grps <- list("LN"=c('LN_KO', 'LN_WT'),
                  "Tumor"=c('Tumor_KO', 'Tumor_WT'))

# Go through each 72hr vs PBS grouping and calculate differential regulon  
de_reg_id <- lapply(splitids, function(grp_idents){
  # fileid <- paste(c("regulon", grp_idents, "pdf"), collapse=".")
  seusub <- subset(seu, ident=grp_idents)
  
  # Create an order of regulons for better heatmap visualization
  all.markers <- FindMarkers(object = seusub, slot = 'data', 
                             ident.1=grp_idents[1], ident.2 = grp_idents[2], 
                             logfc.threshold = 0)
  ord_regulons <- rownames(all.markers[hclust(dist(all.markers$avg_log2FC))$order,])
  
  # Calculate the differential regulons between grp_ident[1] and grp_ident[2] per cluster
  clusters <- c('manual_anno', 'seurat_clusters')
  diffexp_reg <- lapply(setNames(clusters,clusters), function(cl_x){
    all_clid <- as.character(unique(seusub@meta.data[,cl_x]))
    
    lapply(setNames(all_clid,all_clid), function(cl_xi){
      print(paste0(cl_x, " - ", cl_xi, "... "))
      cells <- Cells(seusub)[seusub@meta.data[,cl_x] == cl_xi]
      seu_sub <- subset(seusub, cells=cells)
      if((length(unique(seu_sub$orig.ident)) ==2) & all(table(seu_sub$orig.ident) > 50)){
        deg <- FindMarkers(seu_sub, ident.1 = grp_idents[1], ident.2 = grp_idents[2], 
                           logfc.threshold = 0) %>% 
          tibble::rownames_to_column(., "gene") %>%
          mutate(biotype=sym2biotype_ids[gene],
                 dir=ifelse(avg_log2FC>0, 'up', 'dn'))
      } else {
        deg <- NULL
      }
      return(deg)
    })
  })
  
  # Identify regulons specific to a certain cluster
  regulon_specific_hms <- lapply(names(diffexp_reg), function(cl_x){
    print(cl_x)
    grp_diffexp_reg <- diffexp_reg[[cl_x]]
    null_idx <- sapply(grp_diffexp_reg, is.null)
    if(any(null_idx)) grp_diffexp_reg <- grp_diffexp_reg[-which(null_idx)]
    
    # Create an identity matrix of regulons by celltypes/clusters, where a +1 
    # for significnat up-regulated, -1 for down and 0 for non-significant regulons
    id_mat <- lapply(grp_diffexp_reg, function(i){
      i %>% 
        mutate(sig=ifelse(p_val_adj < 0.001 &
                            abs(avg_log2FC) > 0.01, 
                          (sign(avg_log2FC)*1), 
                          0)) %>% 
        select(gene, sig)
    }) %>% 
      purrr::reduce(., full_join, by='gene') %>% 
      tibble::column_to_rownames(., "gene") %>%
      rename_with(., ~ names(grp_diffexp_reg))
    id_mat[is.na(id_mat)] <- 0
    
    # Remove non-informative regulons
    nonsig_idx <- which(rowSums(abs(id_mat)) ==0)
    clus_hc <- hclust(dist(t(id_mat[-nonsig_idx,]), method='binary'))
    reg_hc <- hclust(dist(id_mat[-nonsig_idx,], method='binary'))
    
    phm <- pheatmap::pheatmap(id_mat[-nonsig_idx,][reg_hc$order, clus_hc$order],
                              cluster_rows=F, cluster_cols=F)
    return(list("heatmap"=phm, "id_mat"=id_mat, "clus_hc"=clus_hc, "reg_hc"=reg_hc))
  })
  names(regulon_specific_hms) <- names(diffexp_reg)
  
  return(list("diffexp_reg"=diffexp_reg,
              "regulon_specific_hms"=regulon_specific_hms))
})

# Iterate through LN and Tumor and compare the KO to WT after treatment conditions
delta_fcs <- lapply(comp_grps, function(grp_idents){
  de_grp <- de_reg_id[grp_idents]
  ko_grp <- de_grp[[1]]$diffexp_reg
  wt_grp <- de_grp[[2]]$diffexp_reg
  
  lapply(setNames(names(ko_grp),names(ko_grp)), function(clust_category){
    ko_grp_cat <- ko_grp[[clust_category]]
    wt_grp_cat <- wt_grp[[clust_category]]
    
    lapply(setNames(names(ko_grp_cat),names(ko_grp_cat)), function(clustid){
      print(clustid)
      ko_clust <- ko_grp_cat[[clustid]]
      wt_clust <- wt_grp_cat[[clustid]]
      if(is.null(ko_clust) | is.null(wt_clust)) return(NULL)
      
      l2fc_mat <- dplyr::full_join(ko_clust %>% 
                                     select(gene, avg_log2FC),
                                   wt_clust %>%
                                     select(gene, avg_log2FC),
                                   by='gene') %>%
        tibble::column_to_rownames(., "gene")
      padj_mat <- dplyr::full_join(ko_clust %>% 
                                     select(gene, p_val_adj),
                                   wt_clust %>%
                                     select(gene, p_val_adj),
                                   by='gene') %>%
        tibble::column_to_rownames(., "gene")
      
      # Get difference between KO and WT FC, calculate exact significance
      deltal2fc=2^l2fc_mat
      set.seed(1234)
      x=sample(deltal2fc[,1], size=100000, replace = T)
      y=sample(deltal2fc[,2], size=100000, replace = T)
      nulld <- na.omit(log2(1+(x-y)))
      deltal2fc <- data.frame("gene"=rownames(l2fc_mat),
                              "l2fc"=log2(1+(deltal2fc[,1]-deltal2fc[,2])))
      pospvals <- sapply(deltal2fc$l2fc, function(i){1-(sum(i>nulld,na.rm=T)/length(nulld))})
      negpvals <- sapply(deltal2fc$l2fc, function(i){1-(sum(i<nulld,na.rm=T)/length(nulld))}) 
      pvalsmat <- matrix(c(pospvals, negpvals), nrow=2, byrow=T)
      pvals <- as.numeric(apply(pvalsmat, 2, min))
      deltal2fc <- deltal2fc %>% 
        mutate(pval=(pvals/2),
               padj=p.adjust(pval, method='BH')) %>%
        arrange(padj) %>% 
        tibble::column_to_rownames(., "gene")
      return(list(l2fc=l2fc_mat, padj=padj_mat, deltal2fc=deltal2fc))
    })
  })
})

fc_mat <- lapply(delta_fcs$LN$manual_anno, function(i) i$deltal2fc)
padj_mat <- lapply(delta_fcs$LN$manual_anno, function(i){
  i$deltal2fc %>%
    as.data.frame %>% 
    tibble::rownames_to_column(., "gene") %>%
    select(gene, padj)
})  %>%
  purrr::reduce(., .f=full_join, by='gene')


## Visualization: Cluster/Celltype specific regulons
pdf(file.path("~/xfer", paste0("celltype.", fileid)), width = 7, height = 10)
regulon_specific_hms$manual_anno$heatmap
plot(regulon_specific_hms$manual_anno$clus_hc)
dev.off()
pdf(file.path("~/xfer", paste0("cluster.", fileid)), width = 8, height = 10)
regulon_specific_hms$seurat_clusters$heatmap
plot(regulon_specific_hms$seurat_clusters$clus_hc)
dev.off()
cat(paste0("\nxfer cluster.", fileid, "\n"))
cat(paste0("xfer celltype.", fileid, "\n"))






## Visualization: Heatmap of cell AUC scores per celltype/cluster
library(viridis)
grp_hms <- lapply(names(diffexp_reg), function(grpid){
  print(grpid)
  
  ## Isolate for only features that are significantly different within clusters 
  # or celltypes
  grp_diffexp_reg <- diffexp_reg[[grpid]]
  features <- lapply(grp_diffexp_reg, function(i){
    isig <- i %>% 
      filter(p_val_adj < 0.001,
             abs(avg_log2FC) > 0.01)
    isig$gene
  }) %>% unlist %>% unique
  # Order based on the clustering above
  ord_features <- ord_regulons[ord_regulons %in% features]
  
  
  # Iterate through each celltype/cluster and plot a heatmap to be
  # merged later on
  seu$plotid <- gsub("MDSC_", "", seu$orig.ident)
  Idents(seu) <- 'plotid'
  celltypes <- unique(seu@meta.data[,grpid])
  hms <- lapply(celltypes, function(celltype_j){
    cells_j <- Cells(seu)[which(seu@meta.data[,grpid] == celltype_j)]
    addlbl <- ifelse(celltype_j == celltypes[1], FALSE, TRUE)
    hm <- DoHeatmap(seu, cells = cells_j, features = ord_features,
                    slot = 'scale.data', raster = T) + 
      scale_fill_viridis() +
      NoLegend() +
      xlab(str_wrap(celltype_j, width = 10)) +
      theme(strip.text = if(addlbl) element_blank(),
            axis.text.y = if(addlbl) element_blank())
    return(hm)
  })
  return(hms)
})
names(grp_hms) <- names(diffexp_reg)

pdf("~/xfer/celltype_regulon_deg.pdf", width = 16, height = 12)
widths <- sapply(grp_hms$manual_anno, function(x) x$data$Cell %>% unique %>% length)
widths <- ceiling(scales::rescale((widths / min(widths)), to=c(2,10)))
cowplot::plot_grid(plotlist=grp_hms$manual_anno, nrow=1, align = "h", 
                   rel_widths = widths)
dev.off()
pdf("~/xfer/cluster_regulon_deg.pdf", width = 24, height = 12)
widths <- sapply(grp_hms$seurat_clusters, function(x) x$data$Cell %>% unique %>% length)
widths <- ceiling(scales::rescale((widths / min(widths)), to=c(2,10)))
cowplot::plot_grid(plotlist=grp_hms$seurat_clusters, nrow=1, align = "h", 
                   rel_widths = widths)
dev.off()



##############################################
#### 2 T-cell subset to identify velocity ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))
seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()
seu$tcells <- ifelse(grepl("^T", seu$manual_anno), "T", "nonT")

Idents(seu) <- 'tcells'
seu_tcell <- subset(seu, ident='T')

seu_tcell <- SCTransform(seu_tcell, assay = 'RNA', new.assay.name = 'SCT',
                         vars.to.regress = c('percent.mt'), conserve.memory = T) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

Idents(seu_tcell) <- 'seurat_clusters'
seu_tcell <- subset(seu_tcell, ident=c(7,14),  invert = TRUE)
seu_tcell <- SCTransform(seu_tcell, assay = 'RNA', new.assay.name = 'SCT',
                         vars.to.regress = c('percent.mt'), conserve.memory = T) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

pdf("~/xfer/sara_tcellseu.pdf")
DimPlot(seu_tcell, group.by = "seurat_clusters", label = T)
DimPlot(seu_tcell, group.by = "manual_anno")
dev.off()

#--- a) Refining cell type annotations ----
celltype_feats <- list(
  "Treg_effector_up"=c('Neb', 'Prdm1','Ccr8','Gm20053','Arl5a','Plxnd1',
                       'Ncmap','N4bp1','Snx9','Glcci1','Capg','Ky','Gcnt1',
                       'Gsg2','Fut8','Bcl2l1'),
  "Treg_effector_down"=c('Nod1', 'Scml4','Ccr7','Ahcyl2','Cd96','Mfhas1',
                         'Kbtbd11','St8sia6','Actn1','Ms4a6b','Qser1','Dzip1',
                         'St8sia1','Sell','Fsd1l','Klf3','Gm4956','Ift80',
                         'Lrrc32','Nsg2','Tmtc4','Klhl6','Gprc5b'),
  "mLN_Treg_centralMem"=c('Ikzf2','Capg','Igfbp4','Satb1','Ms4a4b',
                          'Rasgrp2','S1pr1','Sell','Crip1','Ly6c1', 'Klf2'),
  "mLN_Treg_effector"=c('Klf2', 'Crip1','Ly6c1','Rasgrp2','Ly6a','S1pr1','Sell','Lgals1'),
  "mLN_Treg_NLTlike"=c('Ccr2','Ly6a','AW112010','S100a4','S100a11',
                       'Glrx','Maf','Lag3','Icos','S100a6'),
  "naive_markers"=c('Irf4', 'Tcf7', 'Cxcf3', 'Ccr7', 'Ccl5',  'Sell', 'Il7r')
)
DefaultAssay(seu_tcell) <- 'RNA'
object.cc <- AddModuleScore(object = seu_tcell, features = celltype_feats, 
                            name = names(celltype_feats), ctrl = 5)
colnames(object.cc@meta.data) <- gsub("[0-9]$", "", colnames(object.cc@meta.data))

dat <- GetAssayData(seu_tcell, slot='data', assay='RNA')
sum_by_gene <- rowSums(dat)
dat <- dat[which(sum_by_gene>median(sum_by_gene)),]
cell_rankings <- AUCell::AUCell_buildRankings(dat, plotStats=F)

cells_AUC <- AUCell::AUCell_calcAUC(
  geneSets=celltype_feats,
  rankings=cell_rankings, 
  aucMaxRank=nrow(cell_rankings)*0.05)

ids <- rownames(assay(cells_AUC))
for(id in ids){
  object.cc@meta.data[,paste0(id, "_aucell")] <- assay(cells_AUC)[id,]
}

pdf("~/xfer/sara_tcellanno.pdf", height = 12, width = 12)
DimPlot(object.cc, group.by="manual_anno")
DimPlot(object.cc, group.by="manual_anno", reduction = 'pca')
FeaturePlot(object.cc, features=ids)
FeaturePlot(object.cc, feature=paste0(ids, "_aucell"))
FeaturePlot(object.cc, feature=paste0(ids, "_aucell"), reduction='pca')
DefaultAssay(object.cc) <- 'RNA'
gois <- c('Ccl5', 'Cxcr3', 'Sell')
object.cc.goi <- ScaleData(object = object.cc, features=gois, 
                           object.cc.goivars.to.regress = c("percent.mt"))
FeaturePlot(object.cc.goi, slot = 'scale.data', feature=gois)
dev.off()

cell_assignment <- AUCell::AUCell_exploreThresholds(
  cells_AUC, plotHist = F, 
  nCores = 1, assignCells = TRUE)
selectedThresholds <- getThresholdSelected(cell_assignment)
for(geneSetName in names(selectedThresholds)){
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )kji9{
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
    # Plot
    plot(cellsTsne, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(cellsTsne)], pch=16) 
  }
}


cells_assigned <- lapply(cell_assignment, function(x) x$assignment)
assignment_tbl <- reshape2::melt(cells_assigned, value.name="cell")
assignment_mat <- with(assignment_tbl, table(cell, L1))
colnames(assignmentTable)[2] <- "geneSet"


#--- b) Slingshot velocity inference ----
DefaultAssay(seu_tcell) <- 'SCT'
sce_tcell <- as.SingleCellExperiment(seu_tcell)

dataset <- wrap_expression(
  counts = GetAssayData(seu_tcell, slot="counts"),
  expression = GetAssayData(seu_tcell, slot="data")
)
model <- infer_trajectory(dataset, ti_slingshot(), pull_if_needed=FALSE)
ti_slingshot()

container_id = "dynverse/ti_slingshot:v1.0.3"
config <- babelwhale::get_default_config()
test_singularity_installation()
babelwhale::pull_container
processx::run("singularity", c("exec", paste0("docker://", container_id), 
                               "echo", "hi"), 
              echo_cmd = TRUE, echo = FALSE)

singularity exec docker://dynverse/ti_slingshot:v1.0.3

sce_tcell <- slingshot(sce_tcell, clusterLabels = 'seurat_clusters', 
                       reducedDim = "PCA", allow.breaks = FALSE,
                       start.clus="2")
lnes <- getLineages(reducedDim(sce_tcell,"PCA"),
                    sce_tcell$seurat_clusters, start.clus = "2")
crvs <- getCurves(lnes)

# this define the cluster color. You can change it with different color scheme.
slingshot_df <- colData(sce_tcell)
slingshotsub_df <- data.frame(ident=slingshot_df@listData$ident,
                              slingshot=slingshot_df@listData$slingPseudotime_1)
pdf("~/xfer/sara_slingshot.pdf")
DimPlot(seu_tcell, group.by = "seurat_clusters", label = T, cols = my_color)

ggplot(slingshotsub_df, aes(x = slingshot, y = ident, 
                            colour = ident)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime") +
  scale_colour_manual(values = my_color)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

my_color <- setNames(gg_color_hue(length(levels(sce_tcell$ident))),
                     unique(as.character(sce_tcell$ident)))
plot(reducedDims(sce_tcell)$PCA, 
     col = my_color[as.character(sce_tcell$ident)], 
     pch=16, 
     asp = 1)
legend("bottomleft",
       legend = names(my_color[levels(sce_tcell$ident)]),  
       fill = my_color[levels(sce_tcell$ident)])
lines(SlingshotDataSet(crvs), lwd=2, type = 'lineages')
dev.off()


# data <- as(as.matrix(seu_sub@assays$RNA@data), 'sparseMatrix')
data <- as(as.matrix(seu@assays$integrated@scale.data), 'sparseMatrix')
pd <- seu@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(expression_data=data,
                         cell_metadata  = seu@meta.data,
                         gene_metadata = fData)

## Step 1: Normalize and pre-process the data
# cds <- preprocess_cds(cds, num_dim = 50, norm_method='log')
cds <- preprocess_cds(cds, num_dim = 100, norm_method='size_only', pseudo_count=0)
## Step 2: Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "orig.ident")
## Step 3: Reduce the dimensions using UMAP
DefaultAssay(seu) <- 'integrated'
reducedDims(cds)$UMAP <- seu@reductions$umap@cell.embeddings
reducedDims(cds)$PCA <- seu@reductions$pca@cell.embeddings[,1:50]
# cds <- reduce_dimension(cds)
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)
## Step 6: Order cells
# cds <- order_cells(cds, reduction_method='UMAP', root_pr_nodes='1')

seu$monocle3_clusters <- cds@clusters$UMAP$clusters
seu$monocle3_partitions <- cds@clusters$UMAP$partitions

pdf("~/xfer/c1i.pdf", width=12); 
DefaultAssay(seu) <- 'SCT'
dp_mc <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="monocle3_partitions", pt.size=0.5, shuffle=TRUE)
dp_mp <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="monocle3_clusters", pt.size=0.5, shuffle=TRUE)
dp_sc <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="seurat_clusters", pt.size=0.5, shuffle=TRUE)
cowplot::plot_grid(dp_mc, dp_sc, ncol=2)
cowplot::plot_grid(dp_mp, dp_sc, ncol=2)
FeaturePlot(seu, features=c("Cd8a", "Marco", "Cd4", "Cd3e", "Foxp3", "Ly6g"), 
            pt.size=0.5, order=T, ncol=3)
DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
        group.by="immgen.fine.cluster", pt.size=0.5, shuffle=TRUE)
plot_cells(cds, color_cells_by='immgen.fine.cluster',show_trajectory_graph=T)
dev.off()


# Prep SCT for DEG testing
saveRDS(seu, file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
saveRDS(cds, file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))


#######################################################################
#### 1.d Differential between groups/clusters to identify subtypes ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))

# seu <- PrepSCTFindMarkers(seu)
subset_of_seu <- FALSE
markers <- list()
min_clus_size <- 100

# Identify cluster-specific genes
marker_f <- file.path(datadir, "seurat_obj", "3_degClustersX.rds")
if(!file.exists(marker_f)){
  ident_cls <- c("monocle3_clusters", "seurat_clusters")
  markers <- lapply(setNames(ident_cls,ident_cls), function(ident_cl){
    Idents(seu) <- ident_cl
    marker_clus_res <- FindAllMarkers(seu, assay = "RNA", slot='data',
                                      test.use='wilcox', max.cells.per.ident=2000, 
                                      verbose = FALSE,
                                      recorrect_umi =  if(subset_of_seu) FALSE else TRUE)
    
    
    clust_ids <- as.character(unique(Idents(seu)))
    cl.markers <- lapply(setNames(clust_ids,clust_ids), function(ident){
      FindConservedMarkers(seu, assay = "SCT", ident.1 = ident, 
                           grouping.var = "orig.ident", verbose = FALSE)
    })
    return(list("findAllMrk"=marker_clus_res, "consMrk"=cl.markers))
  })
  saveRDS(markers, file=marker_f)
} else {
  if(!exists("markers")) markers <- readRDS(file=marker_f)
}


# Identify genes that differentiate clusters within the same partition
seu.prts <- SplitObject(seu, split.by = "monocle3_partitions")
partition_markers <- lapply(seu.prts, function(seu_prt_i){
  print(unique(seu_prt_i$monocle3_partitions))
  subset_of_seu <- TRUE
  markers <- lapply(setNames(ident_cls,ident_cls), function(ident_cl){
    Idents(seu_prt_i) <- ident_cl
    
    marker_clus_res <- FindAllMarkers(seu_prt_i, assay = "SCT", test.use='wilcox',
                                      max.cells.per.ident=2000, verbose = FALSE,
                                      min.cells.group = min_clus_size,
                                      recorrect_umi =  if(subset_of_seu) FALSE else TRUE)
    if(nrow(marker_clus_res) > 0){
      sig_marker_clus_res <- marker_clus_res %>% 
        filter((p_val_adj < 0.001) &
                 (avg_log2FC >= 2)) %>%
        group_by(cluster) %>%
        top_n(10, avg_log2FC) %>% 
        arrange(., cluster, avg_log2FC) %>% as.data.frame
    } else {
      sig_marker_clus_res <- NULL
    }
    
    
    clust_ids <- as.character(unique(Idents(seu_prt_i)))
    cl_markers <- lapply(setNames(clust_ids,clust_ids), function(ident){
      print(ident)
      tryCatch({
        FindConservedMarkers(seu_prt_i, assay = "SCT", ident.1 = ident, 
                             grouping.var = "orig.ident", verbose = FALSE,
                             min.cells.group = min_clus_size/length(unique(seu_prt_i$orig.ident)),
                             recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>%
          select(minimump_p_val, max_pval) %>% 
          tibble::rownames_to_column(., "genes") %>%
          mutate(cluster=ident)
      }, error=function(e){NULL})
    })
    
    if(all(is.null(unlist(cl_markers)))){
      sig_cl_markers <- NULL
    } else {
      sig_cl_markers <- cl_markers %>%
        plyr::rbind.fill(.) %>%
        filter(max_pval < 0.05) %>%
        group_by(cluster) %>%
        top_n(5, desc(max_pval)) %>%
        arrange(cluster, minimump_p_val)
    }
    
    
    return(list("findAllMrk"=marker_clus_res, "consMrk"=cl_markers,
                "sigFindAllMrk"=sig_marker_clus_res,
                "sigConsMrk"=sig_cl_markers))
  })
  return(markers)
})
saveRDS(partition_markers, file=file.path(datadir, "seurat_obj", "3_degPartitionClustersX.rds"))


## Visualize the cluster and partition markers
markers <- readRDS(file.path(datadir, "seurat_obj", "3_degClustersX.rds"))
partition_markers <- readRDS(file.path(datadir, "seurat_obj", "3_degPartitionClustersX.rds"))

# INTERACTIVE
all_partition_markers <- lapply(names(partition_markers), function(partition){
  print(partition)
  cons <- tryCatch({
    partition_markers[[partition]]$seurat_clusters$sigConsMrk %>% 
      mutate(partition=partition,
             marker_type='ConservedMarkers')
  }, error=function(e){NULL})
  fam <- tryCatch({
    partition_markers[[partition]]$seurat_clusters$sigFindAllMrk %>% 
      mutate(partition=partition,
             marker_type='FindAllMarkers') %>%
      rename('genes'='gene')
  }, error=function(e){NULL})
  tryCatch({
    full_join(as.data.frame(cons), as.data.frame(fam), 
              by=c('cluster', 'partition', 'genes')) %>%
      arrange(cluster)
  }, error=function(e){NULL})
}) %>% 
  do.call(rbind, .)
write.table(all_partition_markers, sep=",", quote=F, col.names = T, row.names = F)
# /INTERACTIVE

# Cluster-level markers
clus_markers <- lapply(names(markers$seurat_clusters$consMrk), function(clus_i){
  markers$seurat_clusters$consMrk[[clus_i]] %>% 
    select(minimump_p_val, max_pval) %>% 
    tibble::rownames_to_column(., "genes") %>%
    mutate(cluster=clus_i) %>%
    group_by(cluster) %>%
    top_n(5, desc(max_pval)) %>%
    arrange(cluster, minimump_p_val)
}) %>% 
  plyr::rbind.fill(.)

clus_markers <- markers$seurat_clusters$findAllMrk %>% 
  filter((p_val_adj < 0.001) &
           (avg_log2FC >= 2)) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>% 
  arrange(., cluster, avg_log2FC) %>% 
  as.data.frame



pdf("~/xfer/dotplot.pdf", width=14)
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- 'seurat_clusters'
DotPlot(seu, features = unique(clus_markers$gene), 
        cols = c("blue", "red"), dot.scale = 4,
        group.by='monocle3_partitions') + RotatedAxis()
DefaultAssay(seu) <- 'RNA'
DoHeatmap(seu, features = unique(clus_markers$gene),
          group.by='monocle3_partitions')
dev.off()

########################################################
#### 1.d Mapping the change in celltype proportions ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))
seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()

ident_meta <- unique(seu$orig.ident) %>%
  strsplit(., split="_") %>%
  do.call(rbind, .) %>%
  as.data.frame %>% 
  rename_with(., ~c("Tissue", "Condition", "Time")) %>%
  mutate(sampleID=unique(seu$orig.ident)) 

.assignClusterToPartition <- function(clusid, partid, annoid=NULL){
  # clusid <- 'seurat_clusters'
  # partid <- 'monocle3_partitions'
  tbl <- table(seu@meta.data[,clusid], seu@meta.data[,partid])
  partitions <- apply(tbl / rowSums(tbl), 1, which.max)
  partitions <- setNames(colnames(tbl)[partitions], names(partitions))
  
  prt_clus_df <- as.data.frame(partitions) %>% 
    tibble::rownames_to_column(., "clusters") %>%
    arrange(partitions)
  
  if(!is.null(annoid)){
    tbl <- table(seu@meta.data[,clusid], seu@meta.data[,annoid])
    annos <- apply(tbl / rowSums(tbl), 1, which.max)
    annos <- setNames(colnames(tbl)[annos], names(annos))
    
    prt_clus_df <- prt_clus_df %>% 
      left_join(as.data.frame(annos) %>% 
                  tibble::rownames_to_column(., "clusters"),
                by='clusters')
  }
  return(prt_clus_df)
} 

.getDeltaTbl <- function(clusid, sampleid, dvar, method='rr', winsorize=TRUE, base_wt=FALSE){
  sort_ord <- unique(c(dvar, colnames(ident_meta), 'sampleID'))
  ident_spl <- ident_meta %>% 
    dplyr::arrange(!!! rlang::syms(sort_ord)) %>% 
    group_split(!! rlang::sym(dvar))
  
  sub_ident_meta <- ident_spl[[1]][,sort_ord[c(2,3,1)]]
  sub_samples <- sapply(ident_spl, function(i) i%>% pull('sampleID'))
  
  # Count the number of cells in each cluster for each sample, normalize to total number of cells
  tbl <- table(list(seu@meta.data[,clusid], seu@meta.data[,sampleid])) %>% 
    apply(., 2, function(i) i/sum(i)) %>%
    round(., 3)
  cluster_hc <- hclust(dist(tbl))$order
  cluster_hc <- rownames(tbl)[cluster_hc]
  dtbl <- switch(method,
                 rr=log2(tbl[,sub_samples[,1]] / tbl[,sub_samples[,2]]),
                 dist=tbl[,sub_samples[,1]] - tbl[,sub_samples[,2]],
                 everything=tbl)
  if(base_wt){
    dtbl <- dtbl[,grep("KO", colnames(dtbl))] - dtbl[,grep("WT", colnames(dtbl))]
  }
  if(winsorize){
    qv <- quantile(dtbl[!is.infinite(dtbl)], c(0.005, 0.995))
    dtbl[dtbl <= qv[1]] <- qv[1]
    dtbl[dtbl >= qv[2]] <- qv[2]
  }
  list("dtbl"=dtbl, "meta"=sub_ident_meta, "cluster_ord"=cluster_hc)
}

clusid <- 'seurat_clusters'
partid <- 'monocle3_partitions'
annoid <- 'manual_anno'
cluster_meta <- .assignClusterToPartition(clusid,partid,annoid)

## Calculate the number and proportion of cells per cluster/partition/annotation
prop_tbl <- table(seu$orig.ident, seu$seurat_clusters) %>% 
  as.data.frame.matrix %>% 
  apply(., 1, function(i) round(i/sum(i),3)) %>%
  #t %>% 
  as.data.frame %>%
  tibble::rownames_to_column("clusters") %>%
  mutate(partitions=setNames(cluster_meta$partitions, cluster_meta$clusters)[clusters],
         anno=setNames(cluster_meta$annos, cluster_meta$clusters)[clusters]) %>%
  relocate(c(partitions, anno), .after=clusters)
write.table(prop_tbl, file=file.path("~/xfer", "celltype.proportions.csv"),
            sep=",", col.names = T, row.names = F, quote = F)


## Calculate the change in proportions between tissue samples
delta <- 'Tissue'  # Tissue [LN/Tumor], Condition [KO/WT], Time [72h/PBS]
sampleid <- 'orig.ident'
gg_deltas <- lapply(c('Tissue', 'Condition', 'Time'), function(delta){
  delta_title <- switch(delta,
                        Condition="Condition [KO-WT]",
                        Tissue="Tissue [LN-Tumor]",
                        Time="Time [72hr-PBS]")
  
  dtbl_l <- .getDeltaTbl(clusid, sampleid, dvar=delta, method='rr')
  # dtbl_l <- .getDeltaTbl(clusid, sampleid, dvar=delta, method='rr', base_wt = TRUE)
  
  # Annotate the delta-proportions based on the paired Tissue (LN/Tumor) 72hr-KO
  dtbl_label <- split(as.data.frame(t(dtbl_l$dtbl)), dtbl_l$meta$Tissue)
  dtbl_labels <- lapply(names(dtbl_label), function(dtbl_id){
    as.data.frame(t(dtbl_label[[dtbl_id]])) %>%
      rename_with(., ~gsub("^.*?_", "", .)) %>%
      tibble::rownames_to_column("Cluster") %>% 
      left_join(., cluster_meta, by=c('Cluster'='clusters')) %>%
      mutate(Group=dtbl_id)
  }) %>% do.call(rbind, .)
  pdf("~/xfer/cluster_prop_diff.rr.scatter.pdf", width = 10, height = 10)
  ggplot(dtbl_labels, aes(x=WT_72h, y=KO_72h, shape=Group, 
                          color=annos, group=annos)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype='dashed', color='grey') +
    geom_segment(aes(x=WT_72h, xend=WT_72h, 
                     y=WT_72h, yend=KO_72h, color=Group))
  dev.off()
  
  # Annotate based on the melted delta-proportions
  dtbl_melt <- melt(dtbl_l$dtbl) %>%
    rename_with(., ~c('Cluster', 'Sample','Delta')) %>%
    mutate(Cluster=factor(as.character(Cluster)), 
           Direction=(Delta>0) %>%
             if_else(., "Pos","Neg")) %>%
    left_join(., cluster_meta, by=c('Cluster'='clusters'))
  dtbl_melt$Cluster <- factor(dtbl_melt$Cluster, levels=dtbl_l$cluster_ord)
  
  
  dtbl_melt <- dtbl_melt %>% 
    mutate(SampleType=factor(gsub(".*_(.*)_.*", "\\1", Sample), c("WT", "KO")),
           TumorType=gsub("_.*", "", Sample),
           SampleID=gsub("_.*_", "_", Sample))
  dtbl_melt_l <- split(dtbl_melt, dtbl_melt$TumorType)
  pdf("~/xfer/cluster_prop_diff.rr.time.pdf", width = 12, height = 5)
  dtbl_ggs <- lapply(dtbl_melt_l, function(dtbl_melt_i){
    ggplot(dtbl_melt_i, aes(y=Delta, x=Cluster, fill=SampleType)) +
      geom_bar(stat='identity', position='dodge') +
      facet_grid(.~annos, scales = 'free',
                 space = 'free') +
      scale_fill_manual(values=c('WT'='#4d4d4d', 'KO'='#bf812d')) +
      theme_classic() +
      geom_hline(yintercept=0, linetype='solid', color = "black", size=0.75) +
      geom_hline(yintercept=c(-5,5), linetype="dashed",  color = "grey", size=0.25) +
      ggtitle(unique(dtbl_melt_i$TumorType)) +
      ylim(-10,10) +
      theme(strip.text.x = element_text(angle = 90),
            axis.text.x = element_text(angle=90),
            legend.position="bottom")
  })
  cowplot::plot_grid(plotlist=dtbl_ggs, nrow = 1)
  dev.off()
  
  ggplot(dtbl_melt, aes(x=Delta, y=Cluster, fill=Direction)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(partitions+annos~Sample, scales = 'free_y', 
               space = 'free_y') +
    scale_fill_manual(values=c('Pos'='#d6604d', 'Neg'='#4393c3')) +
    theme_classic() +
    ggtitle(delta_title) +
    theme(strip.text.y.right = element_text(angle = 0))
})

pdf("~/xfer/cluster_prop_diff.relativerisk.pdf", width = 9)
gg_deltas
dev.off()

pdf("~/xfer/dimplot_clus_anno.pdf", width=12); 
DefaultAssay(seu) <- 'SCT'
dp_mp <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="monocle3_partitions", pt.size=0.5, shuffle=TRUE)
dp_mc <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="monocle3_clusters", pt.size=0.5, shuffle=TRUE)
dp_sc <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="seurat_clusters", pt.size=0.5, shuffle=TRUE)
dp_anno <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by="manual_anno", pt.size=0.5, shuffle=TRUE)
cowplot::plot_grid(dp_mp, dp_mc, ncol=2)
cowplot::plot_grid(dp_mp, dp_sc, ncol=2)
cowplot::plot_grid(dp_anno, dp_sc, ncol=2)
FeaturePlot(seu, features=c("Cd8a", "Il2rb", 'Klra1', 'Prf1', "Cd4", "Cd3e", "Foxp3"), 
            pt.size=0.5, order=T, ncol=3)
Idents(seu) <- 'manual_anno'
VlnPlot(object = seu, features = c('Il2rb', 'Klra1', 'Prf1', 'Cd28', 'Foxp3', 'Il2rb'), stack=T); 
DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
        group.by="immgen.fine.cluster", pt.size=0.5, shuffle=TRUE)
# plot_cells(cds, color_cells_by='manual_anno',show_trajectory_graph=T)
dev.off()







# 
# 
# 
# lapply( markers$seurat_clusters$consMrk, function(i){
#   i %>%
#     filter(max_pval < 0.05) %>%
#     top_n(5, desc(max_pval)) %>%
#     select(max_pval,minimump_p_val)
# })
# 
# 
# top_markers <- marker_clus_res %>%
#   filter((pct.1 >= 0.10) & (pct.2 >= 0.10) & ((avg_log2FC)>=1.5)) %>%
#   group_by(cluster) %>%
#   top_n(5, desc(p_val_adj)) %>%
#   arrange(., cluster, desc(p_val_adj))
# 
# marker_clus_res <- top_markers(cds, group_cells_by="cluster", 
#                                reference_cells=1000, cores=8)
# # marker_test_res <- top_markers(cds, group_cells_by="partition", 
# #                                reference_cells=1000, cores=8)
# 
# 
# DefaultAssay(seu_sub) <- 'RNA'
# Idents(seu_sub) <- 'monocle3_clusters'
# marker_clus_res <- FindAllMarkers(seu_sub, test.use = "wilcox", slot='data',
#                                   assay='RNA', logfc.threshold = 0.25, 
#                                   max.cells.per.ident=2000)
# top_markers <- marker_clus_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(5, pseudo_R2) %>%
#   arrange(., cell_group, marker_test_q_value)
# 
# pdf("~/xfer/c2h.pdf"); 
# plot_genes_by_group(cds,
#                     unique(top_markers %>% pull(gene_id)),
#                     group_cells_by="cluster",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
# dev.off()

###############################
#### 2. Velocity Analysis #####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))
subset_of_seu <- TRUE
outdirdeg <- file.path(outdir, "degs")
dir.create(outdirdeg, recursive = F)

seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()

##... b.3) TransferAnchors ----
ref.mnn <- readRDS(file.path(projectils_refdir, "seu.mnn.small.sara.postSCT.rds")) #seu.mnn.small.sara.rds"))
ref.stacas <- readRDS(file.path(projectils_refdir, "seu.mnn.small.sara.rds"))

# ref.mnn <- SCTransform(ref.mnn, assay = 'RNA', new.assay.name = 'SCT', vst.flavor='v2',
#                        vars.to.regress = c('percent.mt', "CC.Difference"),
#                        conserve.memory = T)
# saveRDS(ref.mnn, file=file.path(projectils_refdir, "seu.mnn.small.sara.postSCT.rds"))

# Manual annotation: mapping seurat clusters to annotations
if(TRUE){
  manual_func_anno <- c("0"="ILC2",
                        "1"="Other",
                        "2"="Other",
                        "3"="Other",
                        "4"="Other",
                        "5"="Other",
                        "6"="Other",
                        "7"="Other",
                        "8"="Other",
                        "9"="ILC2",
                        "10"="Other",
                        "11"="Other",
                        "12"="Other",
                        "13"="Other",
                        "14"="Other",
                        "15"="Other",
                        "16"="Other",
                        "17"="Other",
                        "18"="Other",
                        "19"="Other",
                        "20"="Other",
                        "21"="Other",
                        "22"="Other",
                        "23"="Other")
}
ref.mnn$functional.cluster <- manual_func_anno[ref.mnn$seurat_clusters]

pdf("~/xfer/x.pdf", width = 15)
DimPlot(ref.mnn, group.by='functional.cluster', raster=T, label=T)
DimPlot(ref.mnn, group.by='dataset_celltype', raster=T, label=T)
dev.off()

DefaultAssay(seu) <- 'RNA'
if(file.exists(mnn_projected_anchors_rds)) {
  print("reading existing rds file for mnn transferAnchors...")
  seu <- readRDS(mnn_projected_anchors_rds)
} else {
  ref.mnn <- RunUMAP(ref.mnn, reduction = "mnn", dims = 1:50, n.neighbors=30L,
                     min.dist = 0.1, return.model=TRUE)
  anchors <- FindTransferAnchors(
    reference = ref.mnn,
    query = seu,
    normalization.method = "SCT",
    reference.reduction = "mnn",
    reference.assay='SCT',
    dims = 1:50
  )
  seu <- MapQuery(
    anchorset = anchors,
    query = seu,
    reference = ref.mnn,
    refdata = list(
      cluster = "seurat_clusters",
      anno = "functional.cluster",
      dataset='dataset_celltype'
    ),
    reference.reduction = "mnn", 
    reduction.model = "umap"
  )
  
  pdf("~/xfer/mnn_transferAnchors.pdf", width=15);
  DimPlot(ref.mnn, group.by='functional.cluster', reduction='umap', label=TRUE, 
          label.size = 3, repel = TRUE, raster=T) + ggtitle("Reference datasets")
  DimPlot(seu, group.by='predicted.cluster', reduction='ref.umap', label=TRUE, 
          label.size = 3, repel = TRUE, raster=T) + ggtitle("Query-projection [cluster-map]")
  DimPlot(seu, group.by='predicted.anno', reduction='ref.umap', label=TRUE, 
          label.size = 3, repel = TRUE, raster=T) + ggtitle("Query-projection [annotation-map]")
  DimPlot(seu, group.by='predicted.anno', reduction='umap', label=TRUE, 
          label.size = 3, repel = TRUE, raster=T) + ggtitle("Query [annotation-map]")
  DimPlot(ref.mnn, group.by='functional.cluster', reduction='umap', 
          label=TRUE, label.size = 3, repel = TRUE) + 
    # geom_point(data.frame(query@reductions$umap@cell.embeddings), 
    #                                          mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.6, 
    #                                          size = pointsize, shape = 17, color = "gray10") + 
    geom_density_2d(data = data.frame(seu@reductions$ref.umap@cell.embeddings), 
                    mapping = aes(x = refUMAP_1, y = refUMAP_2), color = "black", 
                    n = 200, h = 2, size = 1) + 
    ggtitle("Density projection of query on reference map") + 
    theme(aspect.ratio = 1)
  
  ggplot(as.data.frame(table(seu$predicted.anno)),
         aes(x=Var1, y=Freq)) +
    geom_bar(stat = 'identity', position = 'dodge') + 
    theme_minimal() +
    theme(axis.text.x =element_text(angle=90)) +
    xlab("") + ylab("Frequency")
  ggplot(as.data.frame(table(seu$predicted.dataset)),
         aes(x=Var1, y=Freq)) +
    geom_bar(stat = 'identity', position = 'dodge') + 
    theme_minimal() +
    theme(axis.text.x =element_text(angle=90)) +
    xlab("") + ylab("Frequency")
  dev.off()
  
  
  saveRDS(seu, file=mnn_projected_anchors_rds)
}


#### Stop! ####

lapply(all_cts, function(ct_i){
  volc_comp <- lapply(names(ct_markers_X), function(comp_id){
    if(comp_id %in% names(seu.list)){
      suffix <- c(' [72h vs PBS]')        # tissue_cond condition
    } else {
      suffix <- c(' [KO vs WT]')        # KO_WT condition
    }
    
    i <- tryCatch({
      ct_markers_X[[comp_id]][[ct_i]] %>%
        mutate(sig=case_when(avg_log2FC >= 0.5 & p_val_adj <= 0.05 ~ "up",
                             avg_log2FC <= -0.5 & p_val_adj <= 0.05 ~ "down",
                             TRUE ~ "ns"),
               log10p=(-1*log10(p_val_adj)))
    }, error=function(e){NULL})
    if(is.null(i)) {
      return(ggplot(data.frame()) + geom_point())
    }
    
    ggplot(i, aes(x=avg_log2FC, y=log10p, fill=sig)) +
      ggrastr::rasterise(geom_point(shape = 21,
                                    colour = "black")) +
      theme_classic() +
      xlim(-4,4) +
      ggtitle(paste0(ct_i, ": ", comp_id, suffix)) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed") + 
      geom_vline(xintercept = c(-0.5, 0.5),
                 linetype = "dashed") +
      scale_fill_manual(values = c("up" = "#ffad73", "down" = "#26b3ff", 
                                   "ns" = "grey"))
  })
  cowplot::plot_grid(plotlist=volc_comp,
                     ncol=4)
})
dev.off()

##  Merge all the DEGs across the different groups and the clusters/celltypes 
# together into one matrix
merge_grp_deg <- lapply(ct_markers_X, function(grp_deg){
  merge_deg <- lapply(names(grp_deg), function(ct_i){
    if(is.null(grp_deg[[ct_i]])) return(data.frame("symbol"=NA))
    grp_deg[[ct_i]] %>% 
      rename_with(., ~paste0(ct_i, ".", .),
                  .cols=grep("avg_log2FC|p_val", colnames(.))) %>%
      select(c(symbol, grep("avg_log2FC|p_val", colnames(.), value=T)))
  }) %>% 
    purrr::reduce(full_join, by='symbol')
  return(merge_deg)
})

for(id in names(merge_grp_deg)){
  idx <- merge_grp_deg[[id]] %>%
    mutate(biotype=ens2biotype_ids[gm$SYMBOL$ENSEMBL[symbol]]) %>%
    relocate(biotype, .after=symbol)
  
  write.table(idx, 
              file=file.path("~/xfer", paste0("degs_celltypes.", id, ".csv")),
              sep=",", row.names = F, col.names = T, quote = F)
  cat(paste0("xfer degs_celltypes.", id, ".csv\n"))
}

## Plot the LFC between IL33 KO and WT for the 72hr-vs-PBS effect
cluster_ids <- sapply(ct_markers_X[c('LN_WT', 'LN_KO')], names) %>% 
  unlist %>% as.character %>% unique

id_list <- list("Tumor"=c('Tumor_KO', 'Tumor_WT'),
                "LN"=c('LN_KO', 'LN_WT'))
grp_cor_ggs <- lapply(names(id_list), function(id_i){
  id <- id_list[[id_i]]
  cor_gg <- lapply(cluster_ids, function(cl_i){
    if(is.null(ct_markers_X[[id[1]]][[cl_i]]) | 
       is.null(ct_markers_X[[id[2]]][[cl_i]])) return(NULL)
    
    ct_df <- ct_markers_X[[id[1]]][[cl_i]] %>%
      select(symbol, avg_log2FC) %>%
      full_join(., ct_markers_X[[id[2]]][[cl_i]] %>%
                  select(symbol, avg_log2FC), by='symbol')
    ct_df[,-1] %>% cor(., use='complete.obs')
    
    gg <- ggplot(ct_df, aes(x=avg_log2FC.x, y=avg_log2FC.y)) +
      geom_point(alpha=0.7) +
      xlim(-5,5) + ylim(-5,5) +
      geom_hline(yintercept=0) + geom_vline(xintercept = 0) +
      xlab(paste0("Log2FC [", id[1], "]")) + 
      ylab(paste0("Log2FC [", id[2], "]")) +
      theme_classic() +
      ggtitle(cl_i)
    
    ct_cor <- ct_df[,-1] %>% cor(., use='complete.obs') #pearson correlation
    ct_cor <- rowDiffs(as.matrix(ct_df[,-1])) %>% abs %>% sum(., na.rm=T)
    ct_cor <- matrix(rep(ct_cor,4), ncol = 2)
    return(list("gg"=gg, "df"=data.frame("ID"=cl_i, "cor"=ct_cor[1,2])))
  }) 
  
  # Merge the cor and cowplot grid the ggs
  cor_df <- lapply(cor_gg, function(i) i$df) %>% 
    do.call(rbind,.) %>% 
    mutate("Sample"=id_i)
  ggs <- lapply(cor_gg, function(i) i$gg) %>%
    cowplot::plot_grid(plotlist=., ncol=3)
  
  return(list("gg"=ggs, "cor"=cor_df))
})

# aggregate the cor into a barplot
maxcor <- lapply(grp_cor_ggs, function(i) i$cor$cor) %>% 
  unlist %>% max
cor_gg <- lapply(grp_cor_ggs, function(i) i$cor) %>%
  do.call(rbind,.) %>%
  mutate(dir=ifelse(cor<0, 'neg', 'pos'),
         Sample=factor(Sample, levels=unique(Sample)),
         ID=factor(ID, levels=unique(ID)),
         cor=cor/maxcor) %>%
  ggplot(., aes(x=ID, y=cor, fill=Sample, group=Sample)) +
  geom_bar(stat='identity', position='dodge') +
  theme_classic() + 
  ylim(0,1) +
  scale_fill_manual(values=c('#cab2d6', '#fb9a99')) +
  xlab("Celltype") + ylab("Pearson_Correlation") +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))
pdf("~/xfer/sara_lm.pdf")
grp_cor_ggs
dev.off()
pdf("~/xfer/sara_lm_cor.pdf", height = 3.5, width = 9)
cor_gg
dev.off()


pdf("~/xfer/sara_dimplot.pdf")
DimPlot(seu, group.by='manual_anno', label = T)
DimPlot(seu, group.by='seurat_clusters', label = T)
DimPlot(seu, group.by='monocle3_partitions', label = T)
dev.off()



###########################################################
#### 3.b Gene-Set Analysis between celltypes and groups ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
if(!exists("ct_degs")) ct_degs <- readRDS(file.path(outdir, "degs", "celltypes_degs.rds"))
ident_meta <- data.frame(sampleID=unique(seu$orig.ident)) %>%
  mutate(Site=gsub("_.*", "", sampleID),
         ST2=gsub("^.*_(.*)_.*$", "\\1", sampleID),
         Time=gsub("^.*_", "", sampleID))
seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()
seu$Site <- with(ident_meta, setNames(Site, sampleID))[seu$orig.ident]
seu$ST2 <- with(ident_meta, setNames(ST2, sampleID))[seu$orig.ident]
seu$Time <- with(ident_meta, setNames(Time, sampleID))[seu$orig.ident]
dir.create(file.path(outdir, "gse"))

##-- a) Using GSEA on the DEGs ----
gsea_dat <- lapply(names(ct_degs)[c(1,3)], function(grp_id){
  grp_i <- ct_degs[[grp_id]]
  ct_gsea <- lapply(names(grp_i)[2], function(ct_id){
    print(paste0(grp_id, ": ", ct_id, "..."))
    grp_ct_i <- grp_i[[ct_id]]
    msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                      'C2'=list('CP:REACTOME'),             # curated gene sets
                      'C5'=list('GO:BP', 'GO:CC', 'GO:MF'))
    obj <- iterateMsigdb(species='Mus musculus', gseaFun, msig_lvls=msig_lvls, 
                         lfc_v=setNames(grp_ct_i$avg_log2FC,
                                        gm$SYMBOL$ENTREZ[grp_ct_i$symbol]))
    
    
    obj_df <- tryCatch({
      obj <- unlist(obj, recursive=F)
      obj_df <- lapply(names(obj), function(obj_id){
        obj[[obj_id]]@result %>% 
          # dplyr::select(ID, enrichmentScore, NES, pvalue, p.adjust) %>% 
          mutate(Group=obj_id,
                 Comparison=grp_id,
                 Celltype=ct_id)
      }) %>% 
        do.call(rbind, .)
      obj_df
    }, error=function(e){NULL})
    return(obj_df)
  })
  names(ct_gsea) <- names(grp_i)[2]
  return(ct_gsea)
})
names(gsea_dat) <- names(ct_degs)[c(1,3)]
saveRDS(gsea_dat, file=file.path(outdir, "gse", "celltypes_gsea.merge.rds"))
gsea_dat <- readRDS(file.path(outdir, "gse", "celltypes_gsea.merge.rds"))


id <- 'NFKB$'
gois <- sapply(gsea_dat, function(i){
  i$Tregs %>% filter(grepl(id, ID)) %>% pull(core_enrichment)
}) %>% as.character %>% 
  strsplit(., split="\\/") %>% 
  unlist %>% gm$ENTREZID$SYMBOL[.]
Idents(seu) <- 'manual_anno'
seusub <- subset(seu, ident='Tregs')
DefaultAssay(seusub) <- 'RNA'
seusub2 <- ScaleData(object = seusub, features=gois, 
                     vars.to.regress = c("percent.mt"))

pdf("~/xfer/sara_hm.pdf")
DoHeatmap(seusub2, assay='RNA', slot='scale.data',features=gois, 
          group.by='orig.ident', size=3.5) +
  theme(text = element_text(size = 10))
dev.off()




dir.create(file.path(outdir, "gse", "cytoscape"))
group <- 'LN_72h'
gs_map <- .mapGsToExactSource(species="Mus musculus")
for(group in names(gsea_dat)){
  for(celltype in names(gsea_dat[[group]])){
    gsea_dat_i <- gsea2CytoscapeFormat(gsea_dat[[group]][[celltype]], gs_map)
    write.table(gsea_dat_i$up, 
                file=file.path(outdir, "gse", "cytoscape", 
                               paste0("gsea_", group, ".", gsub(" ", "_", celltype), ".up.xls")),
                sep="\t", col.names = T, row.names = F, quote = F)
    write.table(gsea_dat_i$down, 
                file=file.path(outdir, "gse", "cytoscape", 
                               paste0("gsea_", group, ".", gsub(" ", "_", celltype), ".down.xls")),
                sep="\t", col.names = T, row.names = F, quote = F)
  }
}


.mapGsToExactSource <- function(msig_lvls, species){
  if(missing(msig_lvls)){
    warning("'msig_lvls' undefined: Defaulting to H, C2, C5, and C8")
    msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                      'C2'=list('CP:REACTOME'),             # curated gene sets
                      'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                      'C8'=list(NULL)) 
  }
  if(missing(species)){
    warning("'species' undefined: Defaulting to Homo sapiens")
    species <- 'Homo sapiens'
  }
  
  main_obj <- lapply(names(msig_lvls), function(mlvl){
    sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
      print(paste0(">", mlvl, ":", sublvl, "..."))
      msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>% 
        as.data.frame %>%
        select(gs_name, gs_exact_source) %>% unique
      with(msig_ds, setNames(gs_exact_source, gs_name))
    })
  }) %>% unlist %>% 
    gsub("[ :-?]", "_", .)
}

gsea2CytoscapeFormat <- function(dat, gs_map, ...){
  if(missing(gs_map)){
    gs_map <- .mapGsToExactSource(...)
  }
  
  dat <- dat %>% 
    mutate(Details=gs_map[ID]) %>%
    select(ID, Description, Details, setSize, enrichmentScore, NES, pvalue,
           p.adjust, qvalues, rank, leading_edge) %>%
    rename_with(., ~c('NAME', 'GS.br..follow.link.to.MSigDB', 'GS.DETAILS', 
                      'SIZE', 'ES', 'NES', 'NOM.p.val', 'FDR.q.val', 
                      'FWER.p.val', 'RANK.AT.MAX', 'LEADING.EDGE'))
  datl <- split(dat, (dat$NES > 0))
  return(list("up"=datl[['TRUE']],
              "down"=datl[['FALSE']]))
}

for(i in names(gsea_dat$HFD_ND)){
  gsea_dat$HFD_ND[[i]] %>% 
    write.table(., file=file.path(outdir, "gse", 
                                  paste0("gsea.", gsub(" ", "_", i), ".csv")),
                sep=",", row.names = F, col.names = T, quote = F)
}
lapply(gsea_dat$HFD_ND, function(i) {
  i %>% 
    arrange(p.adjust) %>% 
    dplyr::group_by(Group) %>% 
    slice_head(n=3)
  # top_n(., n=-2, wt=p.adjust)
})

##-- b) Using AUCell on a per-cell basis [INCOMPLETE] ----
species <- 'Homo sapiens'
mlvl <- 'H'
sublvl <- NULL
msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame()

auccellFun <- function(msig_ds, expr_mat){
  msig_list <- with(msig_ds, split(entrez_gene, gs_name)) %>%
    lapply(., function(i){
      gm$ENTREZID$SYMBOL[i]
    })
  dat <- GetAssayData(seu, slot='data', assay='RNA2')
  sum_by_gene <- rowSums(dat)
  dat <- dat[which(sum_by_gene>median(sum_by_gene)),]
  cell_rankings <- AUCell::AUCell_buildRankings(dat, plotStats=F)
  
  gsea <- tryCatch({
    cells_AUC <- AUCell::AUCell_calcAUC(
      geneSets=msig_list,
      rankings=cell_rankings, 
      aucMaxRank=nrow(cell_rankings)*0.05)
    # seu@meta.data$mTORC <- assay(cells_AUC)['HALLMARK_MTORC1_SIGNALING',]
    # seu@meta.data$h_g2m <- assay(cells_AUC)['HALLMARK_G2M_CHECKPOINT',]
    # pdf("~/xfer/aucell2.pdf")
    # FeaturePlot(seu, feature=c('mTORC', 'h_g2m'))
    # dev.off()
    
    cell_assignment <- AUCell::AUCell_exploreThresholds(
      cells_AUC, plotHist = F, 
      nCores = 1, assignCells = TRUE
    )
  }, error=function(e){NULL})
  return(gsea)
}



