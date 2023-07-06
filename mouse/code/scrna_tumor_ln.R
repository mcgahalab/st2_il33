library(monocle3)
library(parsnip)
library(slingshot)
library(tradeSeq)
library(tidyr)
library(msigdbr)
library(BiocParallel)
library(GenomicFeatures)
library(harmony)
library(org.Mm.eg.db)
library(SingleR)
library(slingshot)
library(reshape2)
library(enrichplot)
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)
library(cellpypes)

visualize <- FALSE
seed <- 1234
set.seed(seed)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor'
setwd(PDIR)

gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
groups <- list.files(datadir)

# Create ENZ -> SYMBOL mapping key
genome_gse <- org.Mm.eg.db
txby <- keys(genome_gse, 'SYMBOL')
ens2sym_ids <- mapIds(genome_gse, keys=txby, column='ENSEMBL',
                  keytype='SYMBOL', multiVals="first")
sym2ens_ids <- ens2sym_ids
sym2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                         keytype='SYMBOL', multiVals="first")
entrez2sym_ids <- setNames(names(sym2entrez_ids), sym2entrez_ids)

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

doublet_quantile_cutoff <- 0.95

###################
#### Functions ####
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")

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
             project = 'st2_il33')
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu.rds"))
rm(seus)

###########################
#### 0.b QC - RNA data ####
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu.rds"))
dir.create(file.path(outdir, "qc"), recursive = F, showWarnings = F)

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
                          feature2 = "percent.mt", shuffle = TRUE)
  plot2 <- FeatureScatter(seu_mt, feature1 = "nFeature_RNA", 
                          feature2 = "nCount_RNA", shuffle = TRUE)
  plot_grid(plot1, plot2, nrow=1)
  dev.off()
  
  # Identify count/feature outliers cells
  feat_fracs <- sapply(split(seu$nFeature_RNA, f=seu$orig.ident), quantile, 
                       probs=seq(0, 0.1, by=0.01)) %>%
    round(., 0) %>% 
    as.data.frame() %>%
    mutate(median=rowMedians(as.matrix(.)))
  count_fracs <- sapply(split(seu$nCount_RNA, f=seu$orig.ident), quantile, 
                        probs=seq(0, 0.1, by=0.01)) %>%
    round(., 0) %>% 
    as.data.frame() %>%
    mutate(median=rowMedians(as.matrix(.)))
}

## https://www.nature.com/articles/s41591-021-01323-8#Sec15
# All cells expressing <200 or >6,000 genes were removed
# cells that contained <400 unique molecular identifiers (UMIs) 
# >15% mitochondrial counts

# Remove outlier samples for count and features
seu_qc <- subset(seu_mt,
                 subset =nCount_RNA > 1000 &
                   nFeature_RNA < 6000)
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
    doubletFinder_contour(seu_i$seu, rm.homotypic=FALSE, max.break=0.5) %>%
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


saveRDS(seu, file=file.path(datadir, "seurat_obj", "2_seu_integratedX.rds"))
seu <- readRDS(file=file.path(datadir, "seurat_obj", "2_seu_integratedX.rds"))

############################################
#### 1.b Preprocessing the counts - RNA ####
# Preprocessing of the RNA counts
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "2_seu_integrated.rds"))

## Find Neighbours and Cluster 
DefaultAssay(seu) <- 'integrated'
seu <- RunPCA(seu, verbose = FALSE)
# # Plot PCA
# pdf("~/xfer/c4.pdf", width=20)
# PCAPlot(seu,
#         split.by = "orig.ident")
# dev.off()

#Confirm #PC's determined explain > 95% of variance
var <- seu@reductions$pca@stdev^2
dvar <- abs(diff(var))
PCNum <- c(min(which(cumsum(dvar)/sum(dvar) >= 0.99)) + 1,
           min(which(cumsum(var)/sum(var) >= 0.8)))
PCNum <- 50 #80% of variance

# Visualizing and selecting PCs
visualize_pcs <- FALSE
if(visualize_pcs){
  pdf(file.path(outdir, "clusters", "pca_dimplot.pdf"))
  pdf(file.path("~/xfer", "pca_dimplot.pdf"))
  .getVar <- function(x) paste0(x, " (", round(max(cumsum(var)[1:x], na.rm=T)/sum(var) * 100,1), "%)")
  for(dim in c(5, 10, 30, 60, 100, PCNum)){
    seu <- RunUMAP(object = seu, dims = 1:dim, n.epochs=200L)
    DimPlot(seu, label = FALSE, reduction = "umap", group.by="orig.ident") + ggtitle(.getVar(dim)) 
  }
  dev.off()
}
DefaultAssay(seu) <- 'integrated'
seu <- RunUMAP(object = seu, dims = 1:PCNum[1], n.neighbors=15L,
               min.dist = 0.1)
# pdf("~/xfer/c5.pdf"); DimPlot(seu, label = FALSE, reduction = "umap",
#                               group.by="seurat_clusters"); dev.off()
seu <- FindNeighbors(seu, dims = 1:PCNum)
seu <- FindClusters(seu, resolution = 1.5)
seu <- BuildClusterTree(seu, reduction='pca')

saveRDS(seu, file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))

#-----------------------------------------------------------------------------
DefaultAssay(seu) <- "integrated"
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 200, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:PCNum)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:PCNum)
seu <- FindClusters(seu, resolution = 1.5)
# Annotate cells using SingleR
bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
overwrite <- TRUE
scell_anno <- lapply(c('label.fine', 'label.main'), function(lbl){
  rds_file <- paste0("celldex_immgen.", gsub("label.", "", lbl), ".cell.rds")
  id <- gsub("^.*immgen(.*).rds", "immgen\\1", rds_file)
  if(!file.exists(file.path(outdir, "annotation", rds_file)) | overwrite){
    print(paste0(" Annotating: ", rds_file))
    singler_anno <- SingleR(test=GetAssayData(seu), ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            clusters=NULL)
    saveRDS(singler_anno, file=rds_file)
  } else {
    singler_anno <- readRDS(file=rds_file)
  }
  return(setNames(list(singler_anno), id))
}) %>% unlist(., recursive=F)
for(id in names(scell_anno)){
  seu@meta.data[,id] <- scell_anno[[id]]$labels
}

# Visualizing and selecting resolution
resolutions <- c(0.01, 0.05, seq(0.1, 3.0, by=0.1))
visualize_res <- FALSE
if(visualize_res){
  dps <- lapply(setNames(resolutions,resolutions), function(res){
    seu <- FindClusters(object = seu, resolution = res)
    dp <- DimPlot(seu, label = FALSE, reduction = "umap", group.by="seurat_clusters") + ggtitle(res) 
    Jl <- lapply(names(scell_anno), function(id){
      tbl <- table(seu$seurat_clusters, seu@meta.data[,id])
      Jclus <- apply(tbl, 1, vegan::diversity) / log(ncol(tbl))
      Jcelltype <- apply(tbl, 2, vegan::diversity) / log(nrow(tbl))
      list("clus"=Jclus,
           "celltype"=Jcelltype)
    })
    names(Jl) <- names(scell_anno)
    Jl_dp <- list("Jl"=sapply(unlist(Jl, recursive=F), mean),
                  "Jl_sd"=sapply(unlist(Jl, recursive=F), sd),
                  "dp"=dp)
    return(Jl_dp)
  })
  Jres <- sapply(dps, function(i) i$Jl)
  Jres_sd <- sapply(dps, function(i) i$Jl_sd) 
  
    

  pdf(file.path(outdir, "clusters", "res_dimplot.pdf"))
  .meltAndLbl <- function(x){
    melt(t(x)) %>%
      mutate(dim=gsub("^.*\\.", "", Var2),
             label=gsub("\\.[a-zA-Z]*$", "", Var2))
  }
  # full_join(.meltAndLbl(Jres),
  #           .meltAndLbl(Jres_sd),
  #           by=c('dim', 'label'), 
  #           suffix=c('.mean', '.sd')) %>% 
  .meltAndLbl(Jres) %>% 
    ggplot(., aes(x=Var1, y=value, color=label)) +
    facet_grid(dim ~ ., space = 'free', scales='free') +
    geom_line() +
    theme_classic() +
    xlab('res') + ylab('J') + ylim(0,1) +
    scale_x_continuous(breaks = seq(0,2,by=0.2)) +
    theme(axis.text.x = element_text(angle=90))
  
  lapply(dps, function(i) i$dp)
  dev.off()
}
vegan::diversity(c(100,1,1,1,1,1))/6
do.call(rbind, dps)

seu <- BuildClusterTree(seu, reduction='pca')

saveRDS(seu, file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))
seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))


#-----------------------------------------------------------------------------
DefaultAssay(seu) <- "integrated"
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 200, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:PCNum)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:PCNum)
seu <- FindClusters(seu, resolution = 1.5)

# https://github.com/satijalab/seurat/issues/1739
seu <- SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT',
                   vars.to.regress = c('percent.mt'), conserve.memory = T) #nCount_RNA regressed in vst
seu <- CellCycleScoring(seu, s.features = str_to_title(cc.genes$s.genes), 
                        g2m.features = str_to_title(cc.genes$g2m.genes),
                        assay='SCT', set.ident = TRUE)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT',
                   vars.to.regress = c('percent.mt', 'CC.Difference'),
                   conserve.memory = T)

## Find Neighbours and Cluster 
seu <- RunPCA(seu, npcs=200, verbose = FALSE)

#Confirm #PC's determined explain > 95% of variance
stdev <- seu@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))

seu <- FindNeighbors(object = seu, dims = 1:PCNum)
seu <- RunUMAP(object = seu, dims = 1:PCNum, n.epochs=200L)
seu <- FindClusters(object = seu, resolution = 1.5)
seu <- BuildClusterTree(seu, reduction='pca')
# clusters <- sapply(dps, function(i) i$cluster) %>% 
#   as.matrix() %>% as.data.frame() %>%
#   rename_with(., ~ resolutions)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu, normalization.method = "LogNormalize", 
                     scale.factor = median(seu$nCount_RNA))
seu <- ScaleData(object = seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(seu, file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))
seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))






resolutions <- seq(0.1, 2.0, by=0.1)
dps <- lapply(resolutions, function(res){
  seu <- FindClusters(object = seu, resolution = res)
  # markersnode <- FindAllMarkers(seu)
  dp <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                group.by="seurat_clusters", pt.size=0.5, shuffle=TRUE) +
    ggtitle(res) + theme(legend.position='none')
  list("gg"=dp, "cluster"=seu$seurat_clusters)
})
seu <- FindClusters(object = seu, resolution = 1.5)
seu <- BuildClusterTree(seu, reduction='pca')

clusters <- sapply(dps, function(i) i$cluster) %>% 
  as.matrix() %>% as.data.frame() %>%
  rename_with(., ~ resolutions)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu, normalization.method = "LogNormalize", 
                      scale.factor = median(seu$nCount_RNA))
seu <- ScaleData(object = seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

markers <- FindAllMarkers(object = seu, only.pos = TRUE, 
                          min.pct = 0.25, thresh.use = 0.25, 
                          test.use='wilcox', slot='data')
top.markers <- do.call(rbind, lapply(split(markers, markers$cluster), head))

# markers <- FindAllMarkers(object = seu, only.pos = TRUE, 
#                           min.pct = 0.25, thresh.use = 0.25, 
#                           test.use='DESeq2', slot='counts')

## manually curated cluster_mapping
clus_key <- setNames(as.character(c(0,1,2,3,3,4,4,4,4,5,6,7,8,9,9,10,11,12,12,
                                    12,12,12,13,13,14,14,14,15,16,17,18)),
                     as.character(c(23,29,24,19,22,25,1,0,17,14,5,8,6,21,26,
                                    18,27,4,2,12,10,20,15,16,9,3,11,28,30,7,13)))
seu$seurat_clusters2 <- clus_key[as.character(seu$seurat_clusters)]

if(visualize){
  pdf(file.path(outdir, "clusters", "resolutions_dimplot.pdf"))
  PlotClusterTree(seu)
  lapply(dps, function(i) i$gg)
  dev.off()
  
  pdf(file.path(outdir, "clusters", "markers_deg.pdf"), width = 15)
  DoHeatmap(seu, features = top.markers$gene, raster=F)
  dev.off()
  
  pdf(file.path(outdir, "clusters", "clusters_dimplot.pdf"), width = 10)
  pdf(file.path("~/xfer", "clusters_dimplot.pdf"), width = 10)
  dp_raw <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                group.by="seurat_clusters", pt.size=0.5, shuffle=TRUE) +
    ggtitle("default_clusters") + theme(legend.position='none')
  dp_manual <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                group.by="seurat_clusters2", pt.size=0.5, shuffle=TRUE) +
    ggtitle("manual_reviewed_clusters") + theme(legend.position='none')
  plot_grid(dp_raw, dp_manual, ncol=2)
  dev.off()
}

saveRDS(markers, file=file.path(outdir, "clusters", "markers.rds"))
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_preprocess.rds"))

###########################################
#### 1.b Cell type annotation - ImmGen ####
# if(!exists("seu")) seu <- readRDS(file.path(datadir, "seurat_obj", "seu_preprocess.rds"))
if(!exists("seu"))seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))

dir.create(file.path(outdir, "annotation"))
getImmgen <- TRUE
overwrite <- TRUE

# RNA preprocess and cell-cycle score
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu,
                     normalization.method = "LogNormalize")
seu <- CellCycleScoring(seu, s.features = str_to_title(cc.genes$s.genes), 
                        g2m.features = str_to_title(cc.genes$g2m.genes),
                        assay='RNA', set.ident = TRUE)
seu <- FindVariableFeatures(seu, 
                            selection.method = "vst",
                            nfeatures = 2000, 
                            verbose = FALSE)
seu <- ScaleData(seu)


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
DefaultAssay(seu) <- 'integrated'
colnames(seu@meta.data) <- gsub("X$", "", colnames(seu@meta.data))

# Makes annotations simpler (e.g. "B Cells (B.Fo) -> B.Fo)
id_map <- unique(seu$immgen.fine.cluster) %>%
  gsub("^.*\\(", "", .) %>%
  gsub("\\)", "", .) %>%
  gsub("(^.*?\\..*)\\..*", "\\1", .) %>%
  setNames(., unique(seu$immgen.fine.cluster))
seu$immgen.fine.cluster.simple <- id_map[seu$immgen.fine.cluster]
seu$immgen.fine.cell.simple <- id_map[seu$immgen.fine.cell]

# Relabel clusters based on manually curated markers
relabel_map <- unique(seu@meta.data[,c('immgen.fine.cluster.simple', 'seurat_clusters')]) %>%
  as.data.frame %>% 
  tibble::remove_rownames(.) %>% 
  tibble::column_to_rownames(., 'seurat_clusters') %>%
  rename_with(., ~"celltype")
rename_key <- c('8'='T.Tregs.Klf2+', '10'='T.Tregs.Klf2+', '19'='T.Tregs', 
                '20'='Unknown', '7'='T.4MEM', '17'='T.8MEM')
relabel_map[names(rename_key), 'celltype'] <- rename_key
seu$`immgen.fine.cluster.simple` <- relabel_map[as.character(seu$seurat_clusters),'celltype']
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
  plot_grid(dpcl, dp1, ncol=2)
  plot_grid(dpcl, dp1b, ncol=2)
  plot_grid(dpcl, dp1c, ncol=2)
  plot_grid(dpcl, dp1d, ncol=2)
  plot_grid(dp1c, dp2, ncol=2)
  dev.off()
}


saveRDS(seu, file=file.path(datadir, "seurat_obj", "4_seu_anno.rds"))

#################################################
#### 1.c Subsetting for just T- and NK-cells ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "4_seu_anno.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

# Relabel cell types according Sara's annotations
manual_anno <- c('0'='Treg', '1'='Macrophages (Mertk-)','2'='Treg','3'='Treg (Cxcl10+)',
                 '4'='Treg', '5'='Macrophages (Mertk+)', '6'='Macrophages (Mertk-)',
                 '7'='Naive CD4 Th2+','8'='Treg (8)','9'='Treg','10'='Treg (10)',
                 '11'='DC (Ccr7+)', '12'='Treg','13'='Monocyte Derived Macrophages',
                 '14'='CD8 Effector', '15'='Treg (15)','16'='Macrophages (Mertk+)',
                 '17'='CD8 Naive', '18'='Treg (18)', '19'='Treg (19)','20'='Macrophages',
                 '21'='Bcells', '22'='Eosinophils', '23'='NK')
seu$anno <- manual_anno[seu$seurat_clusters]
saveRDS(seu, file = file.path(datadir, "seurat_obj", "5_seu_manual_anno.rds"))

# pype_from_seurat <- function (seurat) {
#   seurat_status <- requireNamespace("Seurat", quietly = TRUE) && 
#     requireNamespace("SeuratObject", quietly = TRUE)
#   if (!seurat_status) 
#     stop("Install Seurat to use this function.")
#   stopifnot(inherits(seurat, "Seurat"))
#   all_graphs <- names(seurat@graphs)
#   if (is.null(all_graphs)) {
#     seurat <- Seurat::FindNeighbors(seurat)
#     all_graphs <- names(seurat@graphs)
#   }
#   snn_graphs <- all_graphs[grepl("_snn$", all_graphs)]
#   graph_choices <- paste(c("integrated", "WNN", "SCT", "RNA"), "snn", sep = "_")
#   graph_choice <- graph_choices[graph_choices %in% snn_graphs][1]
#   if ("umap" %in% names(seurat@reductions)) {
#     dimension_names <- c("UMAP_1", "UMAP_2")
#   } else if ("tsne" %in% names(seurat@reductions)) {
#     dimension_names <- c("tSNE_1", "tSNE_2")
#   } else if (any(grepl("umap|tsne", names(seurat@reductions)))) {
#     use <- names(seurat@reductions)[grepl("umap|tsne|UMAP|TSNE", 
#                                           names(seurat@reductions))]
#     use <- use[1]
#     dimension_names <- colnames(seurat@reductions[[use]]@cell.embeddings)
#   } else {
#     stop("Neither UMAP nor tSNE found.")
#   }
#   list(raw = SeuratObject::GetAssayData(seurat, "counts"), 
#        neighbors = methods::as(seurat@graphs[[graph_choice]], 
#                                "dgCMatrix") > 0.1, embed = Seurat::FetchData(seurat, 
#                                                                              dimension_names), totalUMI = seurat$nCount_RNA)
# }
# 
# pype <- seu %>%
#   pype_from_seurat %>%
#   rule("B",           "Ms4a1",   ">", 1)                    %>%
#   rule("CD14+ Mono",  "Cd14",    ">", 1)                    %>%
#   rule("CD14+ Mono",  "Lyz1",     ">", 20)                   %>%
#   rule("DC",          "Fcer1a",  ">", 1)                    %>%
#   rule("Platelet",    "Ppbp",    ">", 30)                   %>%
#   rule("T",           "Cd3e",    ">", 3.5)                  %>% 
#   rule("CD8+ T",      "Cd8a",    ">", .8,  parent="T")      %>%
#   rule("CD4+ T",      "Cd4",     ">", .05, parent="T")      %>%
#   # rule("Treg",        "Cxcr3",     ">", 2, parent="T")      %>%
#   rule("Treg",        "Foxp3",     ">", 2, parent="T")      %>%
#   rule("Treg_Klf2+",   "Klf2",     ">", 2, parent="Treg")      %>%
#   rule("Naive CD4+",  "Ccr7",    ">", 1.5, parent="CD4+ T") 
# #%>%
#   #rule("Memory CD4+",  "S100a4", ">", 13,  parent="CD4+ T")
# 
# pdf("~/xfer/test.pdf")
# plot_classes(pype)+ggtitle("PBMCs annotated with cellpypes")
# dev.off()

subset_celltypes <- c('Treg', 'Treg (Cxcl10+)', 'Treg (8)', 'Treg (10)', 
                      'Treg (15)', 'Treg (19)', 'Treg (18)', 'CD8 Effector',
                      'CD8 Naive','Naive CD4 Th2+')
Idents(seu) <- 'anno'
seu_sub <- subset(seu, idents=subset_celltypes)
DefaultAssay(seu_sub) <- 'integrated'

#Confirm #PC's determined explain > 95% of variance
seu_sub <- RunPCA(seu_sub, npcs=500, verbose = FALSE)
var <- seu_sub@reductions$pca@stdev^2
dvar <- abs(diff(var))
PCNum <- c(min(which(cumsum(dvar)/sum(dvar) >= 0.99)) + 1,
           min(which(cumsum(var)/sum(var) >= 0.8)))
PCNum <- 213 #80% of variance

# Visualizing and selecting PCs
DefaultAssay(seu_sub) <- 'integrated'
seu_sub <- RunUMAP(object = seu_sub, dims = 1:PCNum[2])

# pdf("~/xfer/test_3.2.pdf", width = 12, height = 10)
# fs <- c('Cd3e', 'Itgax', 'Itgam', 'Foxp3', 'Cd8a', 'Cd4', 'Cd19', 'Ikzf2',
#         'Adgre1', 'H2-Ab1', 'Siglecf', 'Siglech', 'Ly6g', 'Mpo', 'Adgre4',
#         'Apoe', 'Cd36', 'Csf1r', 'Cx3cr1', 'Ccr2', 'Mertk', 'Klf2', 'Irf4', 'Il1b',
#         'Ly6c', 'Marco', 'Pparg', 'Mrc1', 'Cd2', 'Siglec1', 'C1q', 'Mafb', 'Arg1', 'Mgl2', 'Trem2',
#         'Ccr2', 'Ly6c2', 'F13q1', 'Cx3cr1', 'Nr4a1', 'Pglyrp1', 'Eno3', 'Ace')
# fs <- c('Cd3e', 'Cd8a', 'Cd4', 'Cxcr3', 'Klf2', 'Cxcl10', 'Foxp3')
# DefaultAssay(seu_sub) <- 'RNA'
# FeaturePlot(seu_sub, features=fs, pt.size=0.5, order=T, ncol=3, slot='counts')
# DefaultAssay(seu_sub) <- 'integrated'
# DimPlot(seu_sub, label=TRUE, repel=TRUE, reduction='umap',
#         split.by='orig.ident', pt.size=0.5, group.by='seurat_clusters', ncol=3)
# DimPlot(seu_sub, label=TRUE, repel=TRUE, reduction='umap',
#         split.by='orig.ident', pt.size=0.5, group.by='immgen.fine.cluster.simple', ncol=3)
# dev.off()

seu <- seu_sub
saveRDS(seu, file = file.path(datadir, "seurat_obj", "5_seu_anno_tcells.rds"))



#########################################################
#### 1.c Finding all markers for each seurat-cluster ####
if(!exists("seu")) seu <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_anno_tcells.rds"))
if(!exists("seu")) seu <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_manual_anno.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

markers_f <- file.path(outdir, "annotation", "rds", "seurat_markers_allcells.rds")
if(!file.exists(markers_f)){
  ## Cluster-level markers
  Idents(seu) <- 'anno'
  DefaultAssay(seu) <- 'integrated'
  clusters <- levels(Idents(seu))
  sample_pairs <- combn(unique(seu$orig.ident), m=2)
  cl_marker <- lapply(setNames(clusters, clusters), function(cluster){
    print(paste0("Finding DEGs for : ", cluster, "..."))
    cluster_markers_f <- gsub(".rds$", paste0("_", cluster, ".rds"), markers_f)
    if(!file.exists(cluster_markers_f)){
      cl_markers <- apply(sample_pairs, 2, function(ids){
        idx <- apply(sample_pairs, 2, function(i) all(ids == i)) %>% which
        print(paste0("DEG for: ", cluster, " (", idx, "/", ncol(sample_pairs), ")"))
        n_ident_1 <- sum(seu$orig.ident==ids[1] & as.character(Idents(seu))==cluster, na.rm=T)
        n_ident_2 <- sum(seu$orig.ident==ids[2] & as.character(Idents(seu))==cluster, na.rm=T)
        
        markers <- tryCatch({
          FindMarkers(seu, slot = "data",test.use = "wilcox", 
                      ident.1=ids[1], ident.2=ids[2], 
                      group.by='orig.ident', subset.ident=cluster) %>%
            mutate(n.ident_1 = n_ident_1,
                   n.ident_2 = n_ident_2,
                   ident_1=ids[1],
                   ident_2=ids[2],
                   cluster=cluster) %>%
            tibble::rownames_to_column(., "gene") %>%
            mutate(ensemble=sym2ens_ids[gene]) %>% 
            mutate(biotype=ens2biotype_ids[ensemble])
        }, error=function(e){
          c(rep(NA, 6), n_ident_1, n_ident_2,  ids[1], 
            ids[2], cluster, rep(NA,2)) %>%
            t %>% 
            as.data.frame %>%
            setNames(., c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2',
                          'p_val_adj', 'n.ident_1', 'n.ident_2', 'ident_1',
                          'ident_2', 'cluster', 'ensemble', 'biotype'))
        })
        return(markers)
      })
      names(cl_markers) <- apply(sample_pairs,2,paste, collapse="_")
      saveRDS(cl_markers, file=cluster_markers_f)
    } else {
      cl_markers <- readRDS(cluster_markers_f)
    }
    return(cl_markers)
  })
  
  ## Celltype-level markers between samples
  Idents(seu) <- 'immgen.fine.cluster.simple'
  DefaultAssay(seu) <- 'integrated'
  celltypes <- levels(Idents(seu))
  ct_marker <- lapply(setNames(celltypes, celltypes), function(celltype){
    print(paste0("Finding DEGs for : ", celltype, "..."))
    celltype_markers_f <- gsub(".rds$", paste0("_", celltype, ".rds"), markers_f)
    if(!file.exists(celltype_markers_f)){
      ct_markers <- apply(sample_pairs, 2, function(ids){
        idx <- apply(sample_pairs, 2, function(i) all(ids == i)) %>% which
        print(paste0("DEG for: ", celltype, " (", idx, "/", ncol(sample_pairs), ")"))
        n_ident_1 <- sum(seu$orig.ident==ids[1] & as.character(Idents(seu))==celltype, na.rm=T)
        n_ident_2 <- sum(seu$orig.ident==ids[2] & as.character(Idents(seu))==celltype, na.rm=T)
        
        markers <- tryCatch({
          FindMarkers(seu, slot = "data",test.use = "wilcox", 
                      ident.1=ids[1], ident.2=ids[2], 
                      group.by='orig.ident', subset.ident=celltype) %>%
            mutate(n.ident_1 = n_ident_1,
                   n.ident_2 = n_ident_2,
                   ident_1=ids[1],
                   ident_2=ids[2],
                   celltype=celltype) %>%
            tibble::rownames_to_column(., "gene") %>%
            mutate(ensemble=sym2ens_ids[gene]) %>% 
            mutate(biotype=ens2biotype_ids[ensemble])
        }, error=function(e){
          c(rep(NA, 6), n_ident_1, n_ident_2,  ids[1], 
            ids[2], celltype, rep(NA,2)) %>%
            t %>% 
            as.data.frame %>%
            setNames(., c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2',
                          'p_val_adj', 'n.ident_1', 'n.ident_2', 'ident_1',
                          'ident_2', 'celltype', 'ensemble', 'biotype'))
        })
        return(markers)
      })
      names(ct_markers) <- apply(sample_pairs,2,paste, collapse="_")
      saveRDS(ct_markers, file=celltype_markers_f)
    } else {
      ct_markers <- readRDS(celltype_markers_f)
    }
    return(ct_markers)
  })
  
  markers <- list("cluster"=cl_marker,
                  "celltype"=ct_marker)
  saveRDS(markers, file=markers_f)
} else {
  markers <- readRDS(markers_f)
}

# Concatenate markers
markers_dts <- lapply(names(markers), function(grp_type){
  marker_j <- markers[[grp_type]]
  merge_cols <- c('gene', grp_type, 'ensemble', 'biotype')
  
  lapply(marker_j, function(marker_i){
    # Remove ident_1/2 column and make non-merge column names unique
    marker_i <- lapply(seq_along(marker_i), function(midx){
      mi <- marker_i[[midx]]
      mi %>% 
        mutate(ids=paste0(ident_1, "-", ident_2),
               avg_log2FC=round(as.numeric(as.character(avg_log2FC)), 3)) %>%
        select(-c(ident_1, ident_2)) %>%
        relocate(c(merge_cols, 'ids')) %>%
        rename_with(.data =., 
                    .cols = -merge_cols,
                    ~paste0(., ".", midx))
    })
    
    marker_i %>% 
      purrr::reduce(left_join, by=merge_cols)
  })
})
names(markers_dts) <- names(markers)
markers_dts_unl <- unlist(markers_dts, recursive=F)

# Write out the DEGs per celltype/cluster and merged each sample-sample comparison
dir.create(file.path(outdir, "annotation", "tsv"), showWarnings = F)
lapply(names(markers_dts_unl), function(m_id){
  m_id_san <- gsub("[^[:alnum:]. ]", "", m_id)
  write.table(markers_dts_unl[[m_id]], 
              file=file.path(outdir, "annotation", "tsv", paste0(m_id_san, ".tcell.tsv")),
              sep="\t", col.names = T, row.names = F, quote = F)
})


####################################################
#### 2.tmp Subset B-cells and identify identity ####
dir.create(file.path(outdir, "markers"), showWarnings = F)
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "4_seu_anno.rds"))
immgen_label <- setNames('immgen.label.main', "main")
marker_file <- paste0("markers_per_celltype.", names(immgen_label), ".rds")

data("pbmc3k")
pbmc3k <- preprocessAzimuth(pbmc3k)
# Annotate
seu <- runAzimuth(ref=pbmc3k, query=seu)

runSingleR <- function(seu, dataset='blueprint',
                       dataset_dir='/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets',
                       overwrite=FALSE){
  
  DefaultAssay(seu) <- 'RNA'
  if(all(seu@assays$RNA@data == seu@assays$RNA@counts)){
    print("logNormalizing data...")
    seu <- NormalizeData(seu,
                         normalization.method = "LogNormalize")
  }
  
  # Read in the reference datasets  
  bed.se <- switch(dataset,
                   blueprint='BlueprintEncodeData.rds',
                   immgen='immgen.rds',
                   mair='MairPBMCData.rds',
                   monaco='Monaco.rds')
  bed.se <- readRDS(file.path(dataset_dir, bed.se))
  if(dataset=='monaco') rownames(bed.se) <- stringr::str_to_title(rownames(bed.se))
  
  # Annotate cell/cluster level using the fine/main levels
  for(lbl in c('label.fine', 'label.main')){
    for(clus in c(TRUE, FALSE)){
      rds_file <- paste0("celldex_", 
                         dataset, ".",
                         gsub("label.", "", lbl), ".", 
                         if(clus) 'cluster' else 'cell', 
                         ".rds")
      id <- gsub("^.*?\\.(.*).rds", paste0(dataset, ".\\1"), rds_file)
      print(paste0("> ", id, "..."))
      
      singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
                              ref=bed.se, 
                              assay.type.test=1, labels=bed.se[[lbl]],
                              de.method="wilcox", genes='sd', de.n=50,
                              clusters=if(clus) seu$seurat_clusters else NULL)
  
      
      if(clus){
        cluster_ids <- setNames(singler_anno$pruned.labels, 
                                as.character(rownames(singler_anno)))
        seu@meta.data[,id] <- cluster_ids[as.character(seu$seurat_clusters)]
      } else {
        seu@meta.data[,id] <- singler_anno$pruned.labels
      }
    }
  }
  return(seu)
}

getNumOfPcs <- function(seu, min_var=0.9){
  sdev <- seu@reductions$pca@stdev
  percent_var <- round(sdev^2/sum(sdev^2),5)
  max_pcs <- which(cumsum(percent_var) <= min_var) %>% max
  return(max_pcs)
}


Idents(seu) <- 'seurat_clusters'
seu_subset_raw <- subset(seu, idents=c('1', '5', '6', '16', '13', '3'))
seu_subset <- seu_subset_raw
seu_subset$old_clusters <- seu$seurat_clusters
DefaultAssay(seu_subset) <- 'integrated'
seu_subset <- RunPCA(seu_subset, npcs = 200, verbose = FALSE)
pcs <- getNumOfPcs(seu_subset, 0.8)
seu_subset <- FindNeighbors(object = seu_subset, dims = 1:(pcs+1))
seu_subset <- FindClusters(object = seu_subset, resolution = 1.8)
seu_subset <- RunUMAP(object = seu_subset, dims = 1:(pcs+1))
seu_subset <- runSingleR(seu_subset, dataset = 'monaco')

noise <- with(seu_subset@meta.data, table(seurat_clusters, monaco.fine.cell)) %>% 
  colSums
noise_idx <- which(seu_subset$monaco.fine.cell %in% names(which(noise < 50)))
seu_subset@meta.data[noise_idx,'monaco.fine.cell'] <- NA
seu_subset_raw$monaco.fine.cell <- seu_subset$monaco.fine.cell
seu_subset_raw$monaco.fine.cluster <- seu_subset$monaco.fine.cluster
seu_subset_raw$new_cluster <- seu_subset$seurat_clusters

pdf("~/xfer/sara_bcells.pdf", width = 17)
dimplot_fun <- function(seu, group){
  DimPlot(seu, reduction = "umap", group.by = group, 
          label = TRUE, label.size = 3,
          repel = TRUE, combine=T)
}
dp1 <- dimplot_fun(seu_subset_raw, 'seurat_clusters')
dp2 <- dimplot_fun(seu_subset_raw, 'immgen.fine.cell.simple')
dp3 <- dimplot_fun(seu_subset_raw, 'immgen.fine.cluster.simple')
plot_grid(dp1, dp2, dp3, ncol=3)

dp1 <- dimplot_fun(seu_subset, 'old_clusters')
dp2 <- dimplot_fun(seu_subset, 'seurat_clusters')
dp3 <- dimplot_fun(seu_subset, 'immgen.fine.cell.simple')
plot_grid(dp1, dp2, dp3, ncol=3)


dp1 <- dimplot_fun(seu_subset, 'immgen.fine.cell.simple')
dp2 <- dimplot_fun(seu_subset, 'monaco.fine.cell')
dp3 <- dimplot_fun(seu_subset, 'monaco.fine.cluster')
plot_grid(dp1, dp2, dp3, ncol=3)

dp1 <- dimplot_fun(seu_subset_raw, 'new_cluster')
dp2 <- dimplot_fun(seu_subset_raw, 'monaco.fine.cell')
dp3 <- dimplot_fun(seu_subset_raw, 'monaco.fine.cluster')
plot_grid(dp1, dp2, dp3, ncol=3)
dev.off()


##############################################################
#### 1.d Relabelling clusters based on Sara's annotations ####
dir.create(file.path(outdir, "markers"), showWarnings = F)
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "4_seu_anno.rds"))
seu$old_clusters <- seu$seurat_clusters
map_key <- as.character(levels(seu$seurat_clusters))
map_key <- setNames(map_key, map_key)
map_key[as.character(c(2,9,0,12,4))] <- '2'
map_key[as.character(c(5,16,1,6))] <- '1'

new_clusters <- factor(as.integer(map_key[seu$seurat_clusters]))
levels(new_clusters) <- c(0:length(levels(new_clusters)))
seu$seurat_clusters <- new_clusters

pdf("~/xfer/sara_new_cluster.pdf")
DimPlot(seu, reduction = "umap", group.by = 'seurat_clusters', 
        label = TRUE, label.size = 3,
        repel = TRUE, combine=T)
dev.off()
saveRDS(seu, file=file.path(datadir, "seurat_obj", "5_seu_anno_recluster.rds"))

##################################################
#### 1.e Differential markers of the clusters ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_anno_recluster.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

markers_f <- file.path(outdir, "annotation", "rds", "seurat_markers_allcells.reanno.rds")
if(!file.exists(markers_f)){
  ## Cluster-level markers
  Idents(seu) <- 'seurat_clusters'
  DefaultAssay(seu) <- 'integrated'
  clusters <- levels(Idents(seu))
  sample_pairs <- combn(unique(seu$orig.ident), m=2)
  cl_marker <- lapply(setNames(clusters, clusters), function(cluster){
    print(paste0("Finding DEGs for : ", cluster, "..."))
    cluster_markers_f <- gsub(".rds$", paste0("_", cluster, ".rds"), markers_f)
    if(!file.exists(cluster_markers_f)){
      cl_markers <- apply(sample_pairs, 2, function(ids){
        idx <- apply(sample_pairs, 2, function(i) all(ids == i)) %>% which
        print(paste0("DEG for: ", cluster, " (", idx, "/", ncol(sample_pairs), ")"))
        n_ident_1 <- sum(seu$orig.ident==ids[1] & as.character(Idents(seu))==cluster, na.rm=T)
        n_ident_2 <- sum(seu$orig.ident==ids[2] & as.character(Idents(seu))==cluster, na.rm=T)
        
        markers <- tryCatch({
          FindMarkers(seu, slot = "data",test.use = "wilcox", 
                      ident.1=ids[1], ident.2=ids[2], 
                      group.by='orig.ident', subset.ident=cluster) %>%
            mutate(n.ident_1 = n_ident_1,
                   n.ident_2 = n_ident_2,
                   ident_1=ids[1],
                   ident_2=ids[2],
                   cluster=cluster) %>%
            tibble::rownames_to_column(., "gene") %>%
            mutate(ensemble=sym2ens_ids[gene]) %>% 
            mutate(biotype=ens2biotype_ids[ensemble])
        }, error=function(e){
          c(rep(NA, 6), n_ident_1, n_ident_2,  ids[1], 
            ids[2], cluster, rep(NA,2)) %>%
            t %>% 
            as.data.frame %>%
            setNames(., c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2',
                          'p_val_adj', 'n.ident_1', 'n.ident_2', 'ident_1',
                          'ident_2', 'cluster', 'ensemble', 'biotype'))
        })
        return(markers)
      })
      names(cl_markers) <- apply(sample_pairs,2,paste, collapse="_")
      saveRDS(cl_markers, file=cluster_markers_f)
    } else {
      cl_markers <- readRDS(cluster_markers_f)
    }
    return(cl_markers)
  })
  
  ## Celltype-level markers between samples
  Idents(seu) <- 'immgen.fine.cluster.simple'
  DefaultAssay(seu) <- 'integrated'
  celltypes <- levels(Idents(seu))
  ct_marker <- lapply(setNames(celltypes, celltypes), function(celltype){
    print(paste0("Finding DEGs for : ", celltype, "..."))
    celltype_markers_f <- gsub(".rds$", paste0("_", celltype, ".rds"), markers_f)
    if(!file.exists(celltype_markers_f)){
      ct_markers <- apply(sample_pairs, 2, function(ids){
        idx <- apply(sample_pairs, 2, function(i) all(ids == i)) %>% which
        print(paste0("DEG for: ", celltype, " (", idx, "/", ncol(sample_pairs), ")"))
        n_ident_1 <- sum(seu$orig.ident==ids[1] & as.character(Idents(seu))==celltype, na.rm=T)
        n_ident_2 <- sum(seu$orig.ident==ids[2] & as.character(Idents(seu))==celltype, na.rm=T)
        
        markers <- tryCatch({
          FindMarkers(seu, slot = "data",test.use = "wilcox", 
                      ident.1=ids[1], ident.2=ids[2], 
                      group.by='orig.ident', subset.ident=celltype) %>%
            mutate(n.ident_1 = n_ident_1,
                   n.ident_2 = n_ident_2,
                   ident_1=ids[1],
                   ident_2=ids[2],
                   celltype=celltype) %>%
            tibble::rownames_to_column(., "gene") %>%
            mutate(ensemble=sym2ens_ids[gene]) %>% 
            mutate(biotype=ens2biotype_ids[ensemble])
        }, error=function(e){
          c(rep(NA, 6), n_ident_1, n_ident_2,  ids[1], 
            ids[2], celltype, rep(NA,2)) %>%
            t %>% 
            as.data.frame %>%
            setNames(., c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2',
                          'p_val_adj', 'n.ident_1', 'n.ident_2', 'ident_1',
                          'ident_2', 'celltype', 'ensemble', 'biotype'))
        })
        return(markers)
      })
      names(ct_markers) <- apply(sample_pairs,2,paste, collapse="_")
      saveRDS(ct_markers, file=celltype_markers_f)
    } else {
      ct_markers <- readRDS(celltype_markers_f)
    }
    return(ct_markers)
  })
  
  markers <- list("cluster"=cl_marker,
                  "celltype"=ct_marker)
  saveRDS(markers, file=markers_f)
} else {
  markers <- readRDS(markers_f)
}

# Concatenate markers
markers_dts <- lapply(names(markers), function(grp_type){
  marker_j <- markers[[grp_type]]
  merge_cols <- c('gene', grp_type, 'ensemble', 'biotype')
  
  lapply(marker_j, function(marker_i){
    # Remove ident_1/2 column and make non-merge column names unique
    marker_i <- lapply(seq_along(marker_i), function(midx){
      mi <- marker_i[[midx]]
      mi %>% 
        mutate(ids=paste0(ident_1, "-", ident_2),
               avg_log2FC=round(as.numeric(as.character(avg_log2FC)), 3)) %>%
        select(-c(ident_1, ident_2)) %>%
        relocate(c(merge_cols, 'ids')) %>%
        rename_with(.data =., 
                    .cols = -merge_cols,
                    ~paste0(., ".", midx))
    })
    
    marker_i %>% 
      purrr::reduce(left_join, by=merge_cols)
  })
})
names(markers_dts) <- names(markers)
markers_dts_unl <- unlist(markers_dts, recursive=F)

# Write out the DEGs per celltype/cluster and merged each sample-sample comparison
dir.create(file.path(outdir, "annotation", "reanno_tsv"), showWarnings = F)
lapply(names(markers_dts_unl), function(m_id){
  m_id_san <- gsub("[^[:alnum:]. ]", "", m_id)
  write.table(markers_dts_unl[[m_id]], 
              file=file.path(outdir, "annotation", "reanno_tsv", paste0(m_id_san, ".tcell.tsv")),
              sep="\t", col.names = T, row.names = F, quote = F)
})

############################################################
#### 1.f Proportion of samples per cluster and celltype ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_anno_recluster.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

# Absolute sample proportions: proportions of clusters per sample
abs_prop_sample <- sapply(split(seu$seurat_clusters, seu$old.ident),
       function(i){
         table(i)/length(i)
       })

# Delta sample proportions: proportions of clusters per sample, compared to the baseline
base_prop_sample <- table(seu$seurat_clusters) / ncol(seu)
rel_prop_sample <- apply(abs_prop_sample, 2, function(i) i-base_prop_sample)


# Absolute proportions: proportions of samples per cluster
abs_prop_cluster <- sapply(split(seu$old.ident, seu$seurat_clusters),
                          function(i){
                            table(factor(i, levels=unique(seu$old.ident)))/length(i)
                          })

# Delta proportions: proportions of samples compared to the baseline per cluster
base_prop_cluster <- table(seu$old.ident) / ncol(seu)
rel_prop_cluster <- apply(abs_prop_cluster, 2, function(i) i-base_prop_cluster)
prop_list <- list("rel_clus"=rel_prop_cluster,
                  "abs_clus"=abs_prop_cluster,
                  "rel_samp"=rel_prop_sample,
                  "abs_samp"=abs_prop_sample)

dir.create(file.path(outdir, "annotation", "prop"), showWarnings = F)
lapply(names(prop_list), function(pid){
  write.table(prop_list[[pid]], file=file.path(outdir, "annotation", "prop", paste0(pid, ".tsv")), 
              sep="\t", col.names = T, row.names = T, quote = F)
})

#################################################################
#### 1.g Trajectory analysis of revised clustered - MONOCLE3 ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_anno_recluster.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

DefaultAssay(seu) <- 'integrated'
Idents(seu) <- 'immgen.main.cluster'
seu_sub <- subset(seu, idents='T cells')
seu_sub <- seu 


pdf("~/xfer/1cg.pdf", width = 9)
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
patchwork::wrap_plots(p1, p2)
dev.off()

data <- as(as.matrix(seu_sub@assays$RNA@data), 'sparseMatrix')
data <- as(as.matrix(seu_sub@assays$integrated@scale.data), 'sparseMatrix')
pd <- seu_sub@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(expression_data=data,
                          cell_metadata  = seu_sub@meta.data,
                          gene_metadata = fData)

## Step 1: Normalize and pre-process the data
# cds <- preprocess_cds(cds, num_dim = 50, norm_method='log')
cds <- preprocess_cds(cds, num_dim = 100, norm_method='size_only', pseudo_count=0)
pdf("~/xfer/cb.pdf"); plot_pc_variance_explained(cds); dev.off()
## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "orig.ident")
## Step 3: Reduce the dimensions using UMAP
DefaultAssay(seu_sub) <- 'integrated'
reducedDims(cds)$UMAP <- seu_sub@reductions$umap@cell.embeddings
reducedDims(cds)$PCA <- seu_sub@reductions$pca@cell.embeddings[,1:50]
# cds <- reduce_dimension(cds)
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)
## Step 6: Order cells
# cds <- order_cells(cds, reduction_method='UMAP', root_pr_nodes='1')
pdf("~/xfer/c1i.pdf"); 
plot_cells(cds, show_trajectory_graph=F)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition", show_trajectory_graph=F)
plot_cells(cds, color_cells_by="nCount_RNA", show_trajectory_graph=F)
plot_cells(cds, color_cells_by="nFeature_RNA", show_trajectory_graph=F)
plot_cells(cds, color_cells_by="G2M.Score", show_trajectory_graph=F)
plot_cells(cds, color_cells_by="Phase", show_trajectory_graph=F)
plot_cells(cds, color_cells_by='seurat_clusters')
plot_cells(cds, color_cells_by='orig.ident', show_trajectory_graph=F)
plot_cells(cds, color_cells_by='immgen.fine.cell', show_trajectory_graph=F)
plot_cells(cds, genes=c("Cd8a", "Marco", "Cd4", "Cd3", "Foxp3", "Ly6c", "Ly6g"),
           show_trajectory_graph=F)
dev.off()

marker_clus_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
# marker_test_res <- top_markers(cds, group_cells_by="partition", 
#                                reference_cells=1000, cores=8)

seu_sub$monocle3_clusters <- cds@clusters$UMAP$clusters
seu_sub$monocle3_partitions <- cds@clusters$UMAP$partitions
DefaultAssay(seu_sub) <- 'RNA'
Idents(seu_sub) <- 'monocle3_clusters'
marker_clus_res <- FindAllMarkers(seu_sub, test.use = "wilcox", slot='data',
                                  assay='RNA', logfc.threshold = 0.25, 
                                  max.cells.per.ident=2000)
top_markers <- marker_clus_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2) %>%
  arrange(., cell_group, marker_test_q_value)

pdf("~/xfer/c2h.pdf"); 
plot_genes_by_group(cds,
                    unique(top_markers %>% pull(gene_id)),
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()



# With regression:
gene_fits <- fit_models(cds, model_formula_str = "~orig.ident")
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")
emb_time_terms <- emb_time_terms %>% mutate(q_value = p.adjust(p_value))
sig_genes <- emb_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)

# With graph autocorrelation:
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# change Matrix::rBind() to rbind()
pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", 
                          cores=4, verbose=T)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
head(pr_test_res %>% arrange(q_value))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)


cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$old_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pdf("~/xfer/c3a.pdf")
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(cds, color_cells_by='seurat_clusters', label_cell_groups=TRUE)

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(46, 14, 35, 22)),
           label_cell_groups=TRUE,
           show_trajectory_graph=TRUE)
dev.off()

#x <- list("cds"=cds, "pr_test_res"=pr_test_res)
saveRDS(x, file="cds_pr.seuUMAP.rds")

##################################################################
#### 1.g Trajectory analysis of revised clustered - SLINGSHOT ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_anno_recluster.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

DefaultAssay(seu) <- 'integrated'
Idents(seu) <- 'immgen.main.cluster'
seu_sub <- subset(seu, idents='T cells')

seu_i <- SCTransform(seu_sub, assay = 'RNA', new.assay.name = 'SCT2',
                     vars.to.regress = c('percent.mt', 'CC.Difference'),
                     conserve.memory = T)
DefaultAssay(seu_i) <- 'SCT2'

pto.pca <- slingshot(Embeddings(seu_i, "pca")[,1:40], 
                     clusterLabels = seu_sub$seurat_clusters, 
                 omega = TRUE, stretch = 0)
# pto.umap <- slingshot(Embeddings(seu_i, "umap"), clusterLabels = seu_sub$seurat_clusters, 
#                      omega = TRUE, stretch = 0)
pto <- pto.pca
# pto <- pto.umap
curves <- slingCurves(pto, as.df = TRUE)

# Get top highly variable genes
top_hvg <- HVFInfo(seu_i) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(residual_variance)) %>% 
  top_n(300, residual_variance) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(seu_i, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(pto)[,2], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

dat_split <- rsample::initial_split(dat_use_df)
dat_train <- rsample::training(dat_split)
dat_val <- rsample::testing(dat_split)

model <- parsnip::rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)

val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)
rsample::metrics(data = val_results, truth, estimate)

var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_gene_name <- names(var_imp)[1:12]



hcl_ord <- NULL
gg_hmps <- lapply(c(1:ncol(assay(pto))), function(lineage_idx){
  pst.ord <- rownames(assay(pto))[order(assay(pto)[,lineage_idx], na.last = NA)]
  
  if(is.null(hcl_ord)){
    print("hcl_ord is NULL")
    dat <- GetAssayData(seu_i, slot='scale.data')
    hcl <- dat[names(var_imp),
               pst.ord[sample(c(1:length(pst.ord)), size=1000, replace = F)]] %>%
      dist(., method='euclidean') %>%
      hclust(.)
    hcl_ord <<- hcl$order
  }
  
  DoHeatmap(seu_i, features=names(var_imp)[hcl_ord], cells=pst.ord, raster=T)
})
pdf("~/xfer/x.pdf", height = 30)
gg_hmps
dev.off()




# library(tradeSeq)
# # fit negative binomial GAM to model relationshipships between gex and pseudotime
# seu.gam <- fitGAM(counts=seu_sub@assays$RNA@counts, sds=pto.pca)
# 
# # test for dynamic expression
# ATres <- associationTest(sce)
# 
# # For visualization, top genes by pseudotime
# topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
# pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
# heatdata <- assays(sce)$counts[topgenes, pst.ord]
# heatclus <- sce$GMM[pst.ord]
# 
# heatmap(log1p(heatdata), Colv = NA,
#         ColSideColors = brewer.pal(9,"Set1")[heatclus])


slingMST2 <- function(x){
  dfs <- lapply(seq_along(metadata(x)$lineages), function(l){
    lin <- metadata(x)$lineages[[l]]
    mst <- metadata(x)$mst
    centers <- do.call(rbind, igraph::V(mst)$coordinates)
    rownames(centers) <- igraph::V(mst)$name
    return(data.frame(centers[lin,,drop=F], Order = seq_along(lin), 
                      Lineage = l, Cluster = lin))
  })
  return(do.call(rbind, dfs))
}
mst <- slingMST2(pto)
curves_bkup <-curves; mst_bkup <- mst
mst$Lineage <- as.character(mst$Lineage)
curves$Lineage <- as.numeric(curves$Lineage)
curves$colors <- c('black', 'orange', 'red')[curves$Lineage]
curves <- curves %>% arrange(Order)

pdf("~/xfer/test.pdf")
DimPlot(seu_i, label=TRUE, repel=TRUE, reduction='pca',
               group.by='immgen.fine.cluster.simple', pt.size=0.5)

ggp <- DimPlot(seu_i, label=TRUE, repel=TRUE, reduction='pca',
        group.by='seurat_clusters', pt.size=0.5, combine=F)
ggf <- FeaturePlot(seu_i, features=top_gene_name, reduction='pca',
               pt.size=0.5, order=T, ncol=4, combine=F)

ggp[[1]] + 
  geom_point(data = mst, aes(x=PC_1, y=PC_2), size = 2) +
  geom_path(data = curves,
            aes(x=PC_1, y=PC_2, group = Lineage),
            color=curves$color, size = 1)
  
ggp[[1]] + 
  geom_path(data = curves,
            aes(x=PC_1, y=PC_2, group = Lineage),
            color=curves$color, size = 1)

ggfs <- lapply(ggf, function(gg_i){
  gg_i + 
    geom_path(data = curves %>% arrange(Order),
              aes(x=PC_1, y=PC_2, group = Lineage)) 
})
cowplot::plot_grid(plotlist=ggfs, ncol=3)
  
dev.off()







plot(rd, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl])
lines(SlingshotDataSet(pto), type = 'l', lwd=2, col='black')
pdf("~/xfer/test.pdf")
ggplot(SlingshotDataSet(pto)) + geom_lines()
# lines(SlingshotDataSet(pto), lwd=2, col='red')
dev.off()
lin1 <- getLineages(rd, cl, start.clus = '1')
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')

# getCurves function:
# construct smooth curves based on all the cells elimitating the problem of
# cells projecting onto vertices of piece-wise linear trajectories
#
# When only a single lineage, it is the principal curve  through the center
# of the data.When two or more lineages, the algorithm averages curves near
# shared cells.Both lineages should agree on cells that have yet to differentiate
crv1 <- getCurves(lin1)

############################################
#### 2.a Split Tumor and LN; preprocess ####
dir.create(file.path(outdir, "dimplot"), showWarnings = F)
dir.create(file.path(outdir, "markers"), showWarnings = F)
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "4_seu_anno.rds"))
immgen_label <- setNames('immgen.label.main', "main")
marker_file <- paste0("markers_per_celltype.", names(immgen_label), ".rds")

#split into tumor-ln, re-run umap on this smaller subset and re-cluster
seu$LN_Tumor <- gsub("_.*", "", seu$orig.ident)
ln_tumor <- unique(seu$LN_Tumor)

# Re-run PCA and UMAP on the subsetted SCT assays
seu_ln_tumor <- lapply(setNames(ln_tumor, ln_tumor), function(ln_or_tumor){
  rds_file <- paste0("seu_subset.", ln_or_tumor, ".rds")
  if(file.exists(file.path(outdir, "annotation", rds_file))){
    bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 
    
    Idents(seu) <- 'LN_Tumor'
    seu_i <- subset(seu, idents=ln_or_tumor)
    DefaultAssay(seu_i) <- 'integrated'
    
    #Confirm #PC's determined explain > 95% of variance
    seu_i <- RunPCA(seu_i, npcs=50, verbose = FALSE)
    stdev <- seu_i@reductions$pca@stdev
    var <- stdev^2
    PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))
    seu_i <- FindNeighbors(object = seu_i, dims = 1:PCNum)
    seu_i <- RunUMAP(object = seu_i, dims = 1:PCNum, n.epochs=200L)
    seu_i <- FindClusters(object = seu_i, resolution = 1.2)
    
    for(lbl in c('label.fine', 'label.main')){
      singler_anno <- SingleR(test=GetAssayData(seu_i), ref=bed.se, 
                              assay.type.test=1, labels=bed.se[[lbl]],
                              clusters=seu_i$seurat_clusters)
      cluster_ids <- setNames(singler_anno$labels, 
                              as.character(rownames(singler_anno)))
      id <- paste0("immgen.", lbl)
      seu_i@meta.data[,id] <- cluster_ids[as.character(seu_i$seurat_clusters)]
    }
    saveRDS(seu_i, file=file.path(outdir, "annotation", rds_file))
  } else {
    seu_i <- readRDS(file=file.path(outdir, "annotation", rds_file))
  }
  
  return(seu_i)
})

saveRDS(seu_ln_tumor, file = file.path(datadir, "seurat_obj", "5_seu_tumorLn_split.rds"))

################################################################################
#### 2.c Split Tumor and LN; Calculate proportion of cells per group/sample ####
if(!exists("seu_ln_tumor")) seu_ln_tumor <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_tumorLn_split.rds"))
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "4_seu_anno.rds"))

clusters <- list('seurat_clusters'=c("orig.ident", 'seurat_clusters'), 
                 'immgen.label.main'=c("orig.ident", 'immgen.fine.cluster.simple'),
                 'clusters_cell'=c("orig.ident", 'immgen.fine.cluster.simple'))

# Calculate proportion of cell-types per sample for the entire sample
celltype_props <- lapply(clusters, function(cl){
  celltype_prop <- table(seu@meta.data[,c(cl[1], cl[2])]) %>% 
      apply(., 1, function(i) i / sum(i)) %>%
      round(., 6) %>% 
      as.data.frame() 
  celltype_prop[is.na(celltype_prop)] <- 0
  return(celltype_prop)
})


# Calculate proportion of cell-types per sample for the Tumor_LN split
clusters <- list('seurat_clusters'=c("orig.ident", 'seurat_clusters'), 
                 'immgen.label.main'=c("orig.ident", 'immgen.label.fine'),
                 'clusters_cell'=c("seurat_clusters", 'immgen.label.fine'))
celltype_props <- lapply(clusters, function(cl){
  celltype_prop <- lapply(seu_ln_tumor, function(seu_i){
    table(seu_i@meta.data[,c(cl[1], cl[2])]) %>% 
      apply(., 1, function(i) i / sum(i)) %>%
      round(., 6) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(., 'cellType')
  }) %>% 
    Reduce(function(x,y) merge(x,y, by='cellType', all=T), .) %>%
    tibble::column_to_rownames(., 'cellType')
  celltype_prop[is.na(celltype_prop)] <- 0
  return(celltype_prop)
})
celltype_props <- celltype_props[-3]


.plotIt <- function(x) {
  meltx <- x %>%
    round(., 3) %>%
    melt() %>%
    rename_with(., ~ c('Sample', 'cellType', 'proportion')) %>%
    left_join(ref_tbl, 
              by=c('Sample'='Sample')) %>%
    mutate(Sample=factor(Sample, levels=ref_lvls),
           ct_join=paste0(cellType, "_", join))
  
  auc_delta_val <- meltx %>%   
    group_by(ct_join) %>%
    mutate(delta=-1*diff(proportion)) %>%
    filter(V3=='PBS')
  auc_delta_vals <- split(auc_delta_val, auc_delta_val$join)
  
  .setVal <- function(x, val=0, dir='pos'){
    if(dir=='pos'){
      x[x<=val] <- val
    } else{
      x[x>=val] <- val
    }
    return(x)
  }
  .gg <- function(x, xlim, rm.y=FALSE){
    rm <- theme(axis.text.y = element_blank(), 
                axis.ticks.y = element_blank())
    gg <- ggplot(x, aes(x=delta, y=cellType)) +
      # facet_grid(cols=vars(join), scales='free') +
      geom_bar(stat='identity', position = 'dodge') +
      xlim(xlim) +
      scale_fill_gradientn(colors = 'black') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90),
            axis.title.y = element_blank(),
            axis.title.x = element_blank())
    if(rm.y) gg + rm else gg
  }
  
  pgrids <- lapply(auc_delta_vals, function(auc_i){
    delta_neg_gg <- auc_i %>%
      mutate(delta=.setVal(delta, dir='neg')) %>%
      .gg(., xlim=c(-0.5,0)) + 
      geom_vline(xintercept = 0) +
      theme(axis.line.y = element_blank())
    
    hm_gg <- ggplot(auc_i, aes(x=join, y=cellType, fill=proportion)) +
      geom_tile() +
      scale_fill_gradientn(colors = colorPal(10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = 'bottom') 
    delta_pos_gg <- auc_i %>%
      mutate(delta=.setVal(delta, dir='pos')) %>%
      .gg(., xlim=c(0,0.5), rm.y=T)
    
    plot_grid(delta_neg_gg, hm_gg, delta_pos_gg, nrow=1, 
              align='h', axis='bt',rel_widths = c(4.5, 1, 3))
  })
  
  plot_grid(plotlist=pgrids, nrow=1)
}
winsorize <- function(x, probs=c(0.05, 0.95)) {
  q <- quantile(as.numeric(x), probs)
  x[x<q[1]] <- q[1]
  x[x>q[2]] <- q[2]
  return(x)
}

col <- c('#feebe2', '#dd3497', '#7a0177')
colorPal <- grDevices::colorRampPalette(col)
# celltype_prop %>% t() %>%
ref_tbl <- colnames(celltype_props[[1]]) %>% 
  strsplit(., "_", .) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>%
  mutate(Sample=colnames(celltype_props[[1]]))
ref_tbl$join <- with(ref_tbl, paste0(V1, "_", V2))
ref_tbl$V1 <- factor(ref_tbl$V1, c('LN', 'Tumor'))
ref_tbl$V2 <- factor(ref_tbl$V2, c('WT', 'KO'))
ref_tbl$V3 <- factor(ref_tbl$V3, c('PBS', '72h'))
ref_lvls <- ref_tbl[with(ref_tbl, order(V1, V2, V3)),]$Sample

gg_cps <- celltype_props[[2]]  %>% 
  t() %>% .plotIt() + ggtitle("Absolute proportions")

pdf(file.path(outdir, "dimplot", "cellType_proportions.pdf"), width = 15)
gg_cps
dev.off()

pdf(file.path(outdir, "dimplot", "dimplot_subsample.pdf"), width = 15)
# pdf(file.path("~/xfer", "dimplot_subsample.pdf"), width = 15)
DimPlot(seu_ln_tumor$Tumor, label=TRUE, repel=TRUE, reduction='umap',
        group.by='immgen.label.main', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$Tumor, label=TRUE, repel=TRUE, reduction='umap',
        group.by='immgen.label.fine', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$Tumor, label=TRUE, repel=TRUE, reduction='umap',
        group.by='seurat_clusters', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$LN, label=TRUE, repel=TRUE, reduction='umap',
        group.by='immgen.label.main', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$LN, label=TRUE, repel=TRUE, reduction='umap',
        group.by='immgen.label.fine', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$LN, label=TRUE, repel=TRUE, reduction='umap',
        group.by='seurat_clusters', pt.size=0.5, split.by='orig.ident')
dev.off()


######################################################################
#### 2.b Split Tumor and LN; differential expression per celltype ####
if(!exists("seu_ln_tumor")) seu_ln_tumor <- readRDS(file = file.path(datadir, "seurat_obj", "5_seu_tumorLn_split.rds"))
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "4_seu_anno.rds"))

# Find all markers that identify each individual `seurat_cluster`, using wilcoxon
seu_markers_ln_tumor <- lapply(seu_ln_tumor, function(seu_i){
  markers <- FindAllMarkers(object = seu_i, only.pos = TRUE, 
                            min.pct = 0.25, thresh.use = 0.25, 
                            test.use='wilcox', slot='data')
  markers
})

# Find markers that separate Treatment condition for each organ and cell-type
# immgen_label <- setNames('immgen.label.fine', "fine")
seu_grp_markers <- lapply(seu_ln_tumor, function(seu_i){
  if(file.exists(file.path(outdir, "clusters", marker_file))) stop("markers exists")
  # Assign metadata (e.g timepoint, treatment, and tissue-specific treatment)
  seu_i$timepoint <- sapply(strsplit(seu_i$orig.ident, split="_"), function(i) i[3])
  seu_i$treatment <- sapply(strsplit(seu_i$orig.ident, split="_"), function(i) i[2])
  seu_i$tissue_treatment <- with(seu_i@meta.data, paste0(LN_Tumor, "_", treatment))
  seu_i$tissue_timepoint <- with(seu_i@meta.data, paste0(LN_Tumor, "_", timepoint))
  
  # Split based on tissue-treatment (e.g. LN_KO and LN_WT to compare DEG betwee timepoints)
  DefaultAssay(seu_i) <- 'RNA'
  Idents(seu_i) <- immgen_label
  tissue_treatment <- c(unique(seu_i$tissue_treatment),
                        unique(seu_i$tissue_timepoint))
  ct_markers <- lapply(setNames(tissue_treatment,tissue_treatment), function(sub){
    col_id <- if(grepl("72h|PBS", sub)) 'tissue_timepoint' else 'tissue_treatment'
    
    ids <- table(seu_i@meta.data[,c(col_id, 'orig.ident')])[sub,]
    ids <- ids[which(ids!=0)]
    
    # Iterate through each celltype (e.g. DEG between 72hr and BL, for T.Tregs)
    celltypes <- levels(Idents(seu_i))
    ct_marker <- lapply(setNames(celltypes, celltypes), function(celltype){
      print(paste0(" > ", col_id, ": ", celltype, "..."))
      tryCatch({
        markers <- FindMarkers(seu_i, slot = "data",test.use = "wilcox", 
                               ident.1=names(ids)[1], ident.2=names(ids)[2], 
                               group.by='orig.ident', subset.ident=celltype) %>%
          mutate(!!paste0("n.", names(ids[1])) := ids[1],
                 !!paste0("n.", names(ids[2])) := ids[2])
        markers_sig <- markers %>%
          filter(p_val_adj <= 0.05)
        
        gsea_val <- iterateMsigdb(species='Mus musculus', 
                                  lfc_v=setNames(markers$avg_log2FC, sym2entrez_ids[rownames(markers)]), 
                                  fun=gseaFun,
                                  msig_lvls=list('H'=list(NULL),
                                                 'C2'=list('CP:REACTOME'),
                                                 'C5'=list('GO:BP', 'GO:CC', 'GO:MF')))
        gsea_val <- unlist(gsea_val, recursive=F)
        gsea_val <- lapply(gsea_val, as.data.frame) %>%
          do.call(rbind, .) %>% 
          select(c(4:11))
        gsea_val <- cbind(gsea_val, 
              do.call(rbind, strsplit(rownames(gsea_val), split="\\.")) %>%
                as.data.frame() %>%
                rename_with(., ~ c("Cmain", "Csub", "geneset")))
        return(list("markers"=markers, "marker_sig"=markers_sig, "gsea"=gsea_val))
      }, error=function(e){NULL})
    })
  })
  return(ct_markers)
})
if(!file.exists(file.path(outdir, "clusters", marker_file))){
  saveRDS(seu_grp_markers,  file=file.path(outdir, "clusters", marker_file))
} else {
  seu_grp_markers <- readRDS(file.path(outdir, "clusters", marker_file))
}

# Output a csv of DEGs, split by Tumor and LN
site_degs <- lapply(seu_grp_markers, function(deg){
  deg_conds <- lapply(deg, function(deg_cond){
    # Append celltypes
    deg_ct <- lapply(deg_cond, function(ct) {
      if(is.null(ct$markers)) return(NULL)
      ct$markers %>% 
        as.data.frame() %>%
        tibble::rownames_to_column(., "gene")
    }) %>%
      data.table::rbindlist(., idcol='celltype')
    return(deg_ct)
  })
  
  # Combine the KO and WT conditions
  deg_cond <- data.table::rbindlist(lapply(deg_conds,function(i){
    i %>% as.data.frame() %>%
      rename_with(.,  ~gsub("\\..*?_", "_", .), starts_with="n.") %>%
      mutate(avg_log2FC=round(avg_log2FC, 3))
  }), idcol='condition', fill=TRUE) %>%
    mutate(condition=gsub("^.*_", "", condition))
  return(deg_cond)
})
# Combine site type (LN and Tumor)
degs <- data.table::rbindlist(site_degs, idcol='site')
degs$ensemble <- sym2ens_ids[degs$gene]
degs$biotype <- ens2biotype_ids[degs$ensemble]
write.table(degs, file=file.path(outdir, "markers", "degs_per_celltype.csv"),
            sep=",", row.names = F, col.names = T, quote = F)

# Output a csv of GSEA values per celltype, split by Tumor and LN
site_gseas <- lapply(seu_grp_markers, function(gsea){
  gsea_conds <- lapply(gsea, function(gsea_cond){
    # Append celltypes
    gsea_ct <- lapply(gsea_cond, function(ct) {
      if(is.null(ct$gsea)) return(NULL)
      ct$gsea %>% 
        as.data.frame() %>%
        select(-core_enrichment)
    }) %>%
      data.table::rbindlist(., idcol='celltype')
    return(gsea_ct)
  })
  
  # Combine the KO and WT conditions
  gsea_cond <- data.table::rbindlist(lapply(gsea_conds,function(i){
    i %>% as.data.frame() %>%
      rename_with(.,  ~gsub("\\..*_?", "_", .), starts_with="n.") %>%
      mutate(NES=round(NES, 3),
             enrichmentScore=round(enrichmentScore, 3))
  }), idcol='condition') 
  return(gsea_cond)
})
# Combine site type (LN and Tumor)
gseas <- data.table::rbindlist(site_gseas, idcol='site') %>% 
  select(-qvalues,-rank)
write.table(gseas, file=file.path(outdir, "markers", "gseas_per_celltype.csv"),
            sep=",", row.names = F, col.names = T, quote = F)


####################################################
#### 2.b Split Tumor and LN: Visualize the DEGs #####
if(!exists("seu_grp_markers")) seu_grp_markers <- readRDS(file.path(outdir, "clusters", marker_file))
if(!exists("seu_ln_tumor")) {
  ln_tumor <- c("LN", "Tumor")
  seu_ln_tumor <- lapply(setNames(ln_tumor, ln_tumor), function(ln_or_tumor){
    rds_file <- paste0("seu_subset.", ln_or_tumor, ".rds")
    seu_i <- readRDS(file=file.path(outdir, "annotation", rds_file))
    return(seu_i)
  })
}
  
# Remove the PBS/72hr WT-KO comparisons (temporary)
seu_grp_markers <- lapply(seu_grp_markers, function(i) i[c('LN_KO', 'LN_WT')])

# Get all unique celltypes
celltypes <- sapply(seu_ln_tumor, function(seu_i){
  unique(seu_i@meta.data[,immgen_label])
}) %>% 
  unlist %>%  sort %>% unique

# Get top_n significnat gene sets per group
top_gs_all <- lapply(seu_grp_markers, function(lt_marker){
  gseas <- lapply(unlist(lt_marker, recursive=F), function(i) {
    if(is.null(i$gsea)) return(NULL)
    i$gsea %>% 
      as.data.frame() %>%
      mutate("C"=paste0(Cmain, ".", Csub),
             geneset=gsub("_", " ", geneset))
  }) 
  top_gs <- lapply(gseas, function(g){
    if(is.null(g)) return(NULL)
    sig_g <- g %>%
      group_by(., C) %>% 
      top_n(., 5, wt=-p.adjust) %>%
      filter(., p.adjust <= 0.05)
    return(split(sig_g$geneset, sig_g$C))
  }) 
  gs_groups <- sapply(gseas, function(i) unique(i$C)) %>% unlist() %>% unique()
  
  return(list("gs"=top_gs, "grps"=gs_groups))
})

gs_groups <- top_gs_all[[1]]$grps
top_gs_sel <- lapply(setNames(gs_groups, gs_groups), function(gsg){
  sapply(unlist(lapply(top_gs_all, function(i) i$gs), recursive=F), function(i) i[[gsg]]) %>% 
    unlist() %>%
    table() %>% 
    sort() %>% 
    tail(., 7) %>%
    names()
})
top_gs_sel <- unlist(top_gs_sel)

# Visualize the barplot of differential gene-sets per organ+celltype
gg_plots <- lapply(seu_grp_markers, function(lt_marker){
  gseas <- lapply(unlist(lt_marker, recursive=F), function(i) {
    if(is.null(i$gsea)) return(NULL)
    i$gsea %>% 
      as.data.frame() %>%
      mutate("C"=paste0(Cmain, ".", Csub),
             geneset=gsub("_", " ", geneset),)
  }) 
  
  # # Get top_n significnat gene sets per group
  # top_gs <- lapply(gseas, function(g){
  #   if(is.null(g)) return(NULL)
  #   sig_g <- g %>%
  #     group_by(., C) %>% 
  #     top_n(., 5, wt=-p.adjust) %>%
  #     filter(., p.adjust <= 0.05)
  #   return(split(sig_g$geneset, sig_g$C))
  # }) 
  # gs_groups <- sapply(gseas, function(i) unique(i$C)) %>% unlist() %>% unique()
  # top_gs_sel <- lapply(setNames(gs_groups, gs_groups), function(gsg){
  #   sapply(top_gs, function(i) i[[gsg]]) %>% 
  #     unlist() %>%
  #     table() %>% 
  #     sort() %>% 
  #     tail(., 5) %>%
  #     names()
  # })
  # top_gs_sel <- unlist(top_gs_sel)
  
  # Subset on the gsea values for only the top selected genesets
  gseas_sel <- gseas %>% 
    data.table::rbindlist(., idcol='celltype') %>%
    mutate('condition'=gsub("^(.*?)\\..*", "\\1", celltype),
           'celltype'=factor(gsub("^.*?\\.", "", celltype) %>%
                               gsub(" \\(.*", "", .), levels=celltypes),
           "sig_p"=as.integer(p.adjust <= 0.05)) %>%
    filter(geneset %in% top_gs_sel) %>%
    select(celltype, geneset, NES, C, condition, sig_p)
  
  # fill in missing celltypes
  keymap <- with(gseas_sel, table(C, gsub(" .*", "", geneset)))
  keymap <- which(t(keymap !=0), arr.ind=T) %>%
    as.data.frame() %>%
    mutate(Cid=rownames(keymap)[C]) %>%
    select(Cid)
  if(any(table(gseas_sel$celltype)==0)) {
    missing_ct <- which(table(gseas_sel$celltype)==0)
    missing_tbl <- with(gseas_sel, table(celltype, geneset))==0
    missing_tbl_idx <- which(missing_tbl, arr.ind=T)[c(1:length(missing_ct)),,drop=F]

    ct_placeholder <- gseas_sel[c(1:length(missing_ct)),,drop=F] %>%
      mutate(NES=NA,
             sig_p=0,
             celltype=rownames(missing_tbl)[missing_tbl_idx[,1]],
             geneset=colnames(missing_tbl)[missing_tbl_idx[,2]]) %>%
      mutate(C=keymap[gsub(" .*", "", geneset),'Cid'])
    gseas_sel <- rbind(gseas_sel, ct_placeholder)
  }
  
  # Fill in missing values and plot
  b <- c(-4, -2,  0, 2, 4)
  ggp <- gseas_sel  %>%
    mutate(NES=DescTools::Winsorize(NES, minval=min(b), maxval=max(b))) %>%
    tidyr::complete(., condition, celltype, tidyr::nesting(geneset, C)) %>%
    mutate(sig_p=replace_na(sig_p, 0)) %>%
    ggplot(., aes(x=celltype, y=geneset, fill=NES)) +
      facet_grid(C ~ condition, scales = "free_y", space='free', switch = "y") +
      geom_tile() +
      geom_point(data=. %>% as.data.frame(), aes(x=celltype, y=geneset, alpha=sig_p)) +
      scale_alpha(range = c(0, 1)) +
      theme_classic() +
      scale_fill_gradientn(limits = c(min(b),max(b)),
                           colours=c("#4575b4", "#d1e5f0", "white", "#fddbc7", "#b2182b"),
                           breaks=b, labels=format(b)) +
      theme(axis.text.y=element_text(size=8),
            axis.text.x=element_text(angle=90, size=9),
            strip.text.y.left = element_text(angle=90),
            panel.spacing.x = unit(2, "lines")) +
      scale_y_discrete(labels =  scales::wrap_format(26)) +
      scale_x_discrete(drop=FALSE)
  return(ggp)  
})
pdf(file.path(outdir, "markers", "gsea_per_celltype.pdf"), width = 10, height = 15)
gg_plots
dev.off()




pdf(file.path(outdir, "markers", "markers_deg.pdf"), width = 15, width = 10)
# pdf(file.path("~/xfer", "markers_deg.pdf"), height = 20, width = 10)
lapply(names(seu_ln_tumor), function(id){
  markers_i <- seu_markers_ln_tumor[[id]]
  seu_i <- seu_ln_tumor[[id]]
  
  top.markers <- do.call(rbind, lapply(split(markers_i, markers_i$cluster), head, n=12L))
  DoHeatmap(seu_i, features = top.markers$gene, raster=T, 
            group.by='seurat_clusters') + ggtitle(id)

})
dev.off()





########################################################################
#### 2.b Split by cell-type and calculate differential between sets ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "seu_anno.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

markers_f <- file.path(outdir, "markers", "ct_split_markers.rds")
if(!file.exists(markers_f)){
  Idents(seu) <- 'bp.fine.cell2'
  celltypes <- na.omit(as.character(unique(Idents(seu))))
  ct_markers <- lapply(celltypes, function(ct){
    seu_ct <- subset(seu, ident=ct)
    Idents(seu_ct) <- 'orig.ident'
    keep_idx <- (table(Idents(seu_ct)) > 50)
    if(any(keep_idx)){
      seu_ct2 <- subset(seu_ct, ident=names(keep_idx)[which(keep_idx)])
      DefaultAssay(seu_ct2) <- 'RNA'
      markers <- FindAllMarkers(object = seu_ct2, only.pos = FALSE, 
                                min.pct = 0.01, thresh.use = 0.01, 
                                test.use='wilcox', slot='data')
    } else {
      markers <- matrix(nrow=0, ncol=0)
    }
    return(markers)
  })
  names(ct_markers) <- celltypes
  
  ct_markers %>%
    do.call(rbind, .) %>%
    mutate(celltype=rep(names(ct_markers), unlist(sapply(ct_markers, nrow)))) %>%
    arrange(celltype, cluster, p_val_adj) %>%
    rename(., avg_logFC=avg_log2FC) %>% 
    write.table(., sep=",", file=gsub(".rds$", ".csv", markers_f),
                quote=F, row.names = F, col.names = T)
  saveRDS(ct_markers, file=markers_f)
} else {
  ct_markers <- readRDS(markers_f)
}


ct_markers_sig <- lapply(ct_markers, function(i){
  tryCatch({i[i$p_val_adj < 0.05,]}, error=function(e){NULL})
})
nullOrRow <- function(x) (is.null(x) || nrow(x)==0)
melt_ct_markers <- setNames(melt(ct_markers_sig[which(!sapply(ct_markers_sig, nullOrRow))]),
                            c("sample", "gene", "variable", "value", "celltype")) %>%
  filter(variable == 'p_val_adj')



if(visualize){
  ct_hms <- lapply(split(melt_ct_markers, melt_ct_markers$celltype), function(ct_m){
    Idents(seu) <- 'bp.fine.cell2'
    top_n_markers <- ct_m %>% 
      group_by(sample) %>%
      top_n(., -5, value) %>%
      as.data.frame() %>% 
      select(gene) %>% 
      unique()
    ct <- unique(ct_m$celltype)
    
    seu_ct <- subset(seu, ident=ct)
    Idents(seu_ct) <- 'orig.ident'
    hm_ct <- DoHeatmap(seu_ct, features = top_n_markers[,1], raster=F,
                       size=3, angle=15, hjust=0) +
      theme(legend.position='none') +
      ggtitle(ct)
  })
  
  table(seu@meta.data[,c('orig.ident', 'bp.fine.cell2')])
  
  pdf(gsub(".rds$", ".pdf", markers_f), height = 15, width = 15)
  # pdf("~/xfer/test.pdf", height = 16, width = 12)
  plot_grid(plotlist=ct_hms, ncol=1)
  dev.off()
}


###########################################
#### 2.c Slingshot trajectory analysis ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "seu_anno.rds"))
dir.create(file.path(outdir, "trajectory"), showWarnings = F)

DefaultAssay(seu) <- 'RNA'
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- RunPCA(seu, npcs=200, verbose = FALSE)

sce <- as.SingleCellExperiment(seu)
sce <- slingshot(sce, clusterLabels = 'orig.ident', reducedDim = "PCA",
                 allow.breaks = FALSE)








################################################################
#### OUTDATED - 7.a) DEG on Markers identifying INF for each cell type ####
library(org.Hs.eg.db)
library(msigdbr)

# Load in msigdb gene sets
m_df <- msigdbr(species = "Homo sapiens")
msig_ds <- list("immune"=list("C7", NULL),
                "pathway"=list("C2", "CP:REACTOME"),
                "gobp"=list("C5", "GO:BP"),
                "gomf"=list("C5", "GO:MF"))
msig_gsl <- lapply(msig_ds, function(mds){
  msigdbr(species = "Homo sapiens", category = mds[[1]], subcategory = mds[[2]]) %>% 
    dplyr::select(gs_name, entrez_gene)
})

# Create mapping from symbol to NCBI
dir.create(file.path(outdir, "differential"), showWarnings = F)
genome_gse <- org.Hs.eg.db
txby <- keys(genome_gse, 'SYMBOL')
gene_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                   keytype='SYMBOL', multiVals="first")
txby <- keys(genome_gse, 'ENTREZID')
entrez_ids <- mapIds(genome_gse, keys=txby, column='SYMBOL',
                     keytype='ENTREZID', multiVals="first")


# Massage the seurat object
# if(!exists("seu")) seu <- SeuratDisk::LoadH5Seurat(file.path(outdir, "seurat_obj", "scIntegrated_preproc.h5seurat"))
if(!exists("seu")) seu <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
anno_type <- 'bp.fine.cell'
Idents(seu) <- 'orig.ident'
seu_inf <- subset(x = seu, idents = c('INFB', 'unstimulated'))
Idents(seu_inf) <- anno_type
DefaultAssay(seu_inf) <- 'SCT'
anno_type <- 'bp.fine.cluster'
cell_types <- unique(sort(seu_inf@meta.data[[anno_type]]))
seu_inf$orig.ident <- seu_inf$old.ident

# Find DEGs between INFB and Unstimulated
degs <- lapply(cell_types, function(ct){
  print(paste0(ct, "..."))
  cnts <- split(seu_inf@meta.data[[anno_type]] == ct, unique(as.character(seu_inf$orig.ident)))
  
  if(!all(sapply(cnts, sum) > 50)){
    return(NULL)
  }
  
  markers <- FindMarkers(seu_inf, ident.1 = 'INFB', 
                         group.by = 'orig.ident', subset.ident=ct)
  markers0 <- markers[which(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > log2(2)),]
  markers_spl <- split(markers0, f=markers0$avg_log2FC > 0)
  dirnm <- setNames(c('UP', 'DN'), c('TRUE', 'FALSE'))
  names(markers_spl) <- dirnm[names(markers_spl)]
  markers_spl
})
names(degs) <- cell_types

# Form the IFN-based geneset per celltype for GSEA
deg_gsl <- lapply(degs, function(deg){
  deg1 <- lapply(names(deg), function(deg_dir){
    deg0 <- deg[[deg_dir]]
    tryCatch({
      data.frame("gs_name"=rep(paste0("INFB-vs-unstimulated_", deg_dir), nrow(deg0)),
                 "entrez_gene"=as.character(gene_ids[rownames(deg0)]))
    }, error=function(e){NULL})
  })
  do.call(rbind, deg1)
})
names(deg_gsl)
dir.create(file.path(outdir, "infb-unstimulated"), recursive = F, showWarnings = F)
saveRDS(deg_gsl, file=file.path(outdir, "infb-unstimulated", "degs.rds"))

################################################################
#### OUTDATED - 7.a) DEG on Markers identifying INF for each cell type ####
if(!exists("deg_gsl")) deg_gsl <- readRDS(file=file.path(outdir, "infb-unstimulated", "degs.rds"))
if(!exists("seu")) seu <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
anno_type <- 'bp.fine.cluster'
Idents(seu) <- 'old.ident'

seu_inf <- subset(x = seu, idents = c('CCP149', 'CCP228', 'CCP033', 'CCP246'))
Idents(seu_inf) <- anno_type
DefaultAssay(seu_inf) <- 'SCT'
cell_types <- names(deg_gsl)
seu_inf$orig.ident <- seu_inf$old.ident

## Calculate the DEGs and run GSEA on the pathways that define them
degs <- lapply(cell_types, function(ct){
  print(paste0(ct, "..."))
  cnts <- split(seu_inf@meta.data[[anno_type]] == ct, as.character(seu_inf$orig.ident))
  if(!all(sapply(cnts, sum) > 50)){
    return(NULL)
  }
  
  markers <- FindMarkers(seu_inf, ident.1 = 'lowISS', 
                         group.by = 'group', subset.ident=ct,
                         features=na.omit(entrez_ids[deg_gsl[[ct]]$entrez_gene]),
                         logfc.threshold=0, min.pct=0)
  markers_universe <- FindMarkers(seu_inf, ident.1 = 'lowISS', 
                                  group.by = 'group', subset.ident=ct,
                                  logfc.threshold=0)
  markers_merge <- rbind(markers, markers_universe[-which(rownames(markers_universe) %in% rownames(markers)),])
  
  # Make gene list
  gene_list <- markers_merge$avg_log2FC
  names(gene_list) <- gene_ids[rownames(markers_merge)]
  gene_list<-sort(na.omit(gene_list), decreasing = TRUE)
  msig <- tryCatch({
    GSEA(gene_list, TERM2GENE = deg_gsl[[ct]], 
         pvalueCutoff = 1, minGSSize=1, maxGSSize=2000)
  }, error=function(e){NULL})
  
  return(msig)
})
names(degs) <- cell_types
deg_df <- do.call(rbind, lapply(degs, function(i) tryCatch({i@result[,-c(2, 8, 11)]}, error=function(e){NULL})))
saveRDS(degs, file=file.path(outdir, "infb-unstimulated", "deg_infb-unstim.rds"))

degs <- readRDS(file.path(outdir, "infb-unstimulated", "deg_infb-unstim.rds"))
deg_df <- do.call(rbind, lapply(degs, function(i) tryCatch({i@result[,-c(2, 8, 11)]}, error=function(e){NULL})))


## RPS Genes


# Reduce GSEA result to a ggplot2 friendly format
godeg <- lapply(names(degs), function(id){
  i <- degs[[id]]
  tryCatch({
    x <- i$go@result
    x$cell_type <- id
    x
  }, error=function(e){NULL})
})
godeg <- do.call(rbind, godeg)
godeg_melt <- godeg[,c('Description', 'NES', 'qvalues', 'cell_type')]
godeg_melt$direction <- godeg$NES > 0
n <- 10
desc_oi <- lapply(split(godeg_melt, godeg_melt$direction), function(i){
  head(i$Description[order(sapply(split(i$NES, i$Description), mean))], n)
})
godeg_melt$Description <- factor(godeg_melt$Description, levels=unique(unlist(desc_oi)))
godeg_melt <- godeg_melt[-which(is.na(godeg_melt$Description)),]
godeg_melt$cell_type <- factor(godeg_melt$cell_type, levels=names(degs))

# Do the plotties
pdf(file.path(outdir, "differential", "deg_per-cell.pdf"), width = 25)
ggplot(data=godeg_melt, aes(x=Description, y=NES, fill=qvalues)) +
  geom_bar(stat='identity') + 
  facet_grid(cols=vars(cell_type),
             scales = "free_y", drop=F) +
  geom_hline(yintercept=0, size=1) +
  coord_flip() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

###############################################
#### 7.a) Add Motif info the seurat object #### 
dir.create(file.path(outdir, "differential"), showWarnings = F)
rerun_analysis <- T

# Massage the seurat object
if(!exists("seu")) seu <- readRDS(file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
anno_type <- 'bp.fine.cluster'
Idents(seu) <- anno_type
DefaultAssay(seu) <- 'peaks'
cell_types <- unique(sort(seu@meta.data[[anno_type]]))

# Add motif information from JASPAR2020
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", species = 'Homo sapiens'))
seu <- AddMotifs(seu, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)


###########################
#### 7.b) Run ChromVAR ####
## Run ChromVAR
chromvar_obj_f <- file.path(outdir, "seurat_obj", "scAtacMotif.rds")
if(!file.exists(chromvar_obj_f)){
  library(BiocParallel)
  BiocParallel::register(MulticoreParam(workers = 1))
  seu <- RunChromVAR(
    object = seu,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    new.assay.name = "chromvar"
  )
  saveRDS(seu, file = chromvar_obj_f)
} else {
  seu <- readRDS(chromvar_obj_f)
}

###########################################################
#### 7.c) Identifying Differential Peaks per cell type ####
chromvar_obj_f <- file.path(outdir, "seurat_obj", "scAtacMotif.rds")
if(!exists("seu")) seu <- readRDS(chromvar_obj_f)
rerun_analysis <- FALSE
anno_type <- 'bp.fine.cluster'
cell_types <- unique(sort(seu@meta.data[[anno_type]]))

## Calculate Differential Accessible regions
#    search for DNA motifs that are overrepresented in a set of peaks that are 
#    differentially accessible between cell types
group_by <- 'group'
grp_smps <- unique(seu@meta.data[,c('group', 'dataset')])
# list("ID" = [sampes-to-use], [column-to-compare], [group-to-compare])
# e.g. list("testA" = c("sampleA", "sampleB"), "treatment", "post") 
#  will compare sampleA (post) to sampleB (treat), based on the "treatment" category
comp_grps <- list("high-low"=list(idents=c(grp_smps[grp_smps$group=='highISS',]$dataset,
                                           grp_smps[grp_smps$group=='lowISS',]$dataset),
                                  'group', 'lowISS'),
                  "high-healthy"=list(idents=c(grp_smps[grp_smps$group=='highISS',]$dataset,
                                               grp_smps[grp_smps$group=='healthy',]$dataset),
                                      'group', 'healthy'),
                  "low-healthy"=list(idents=c(grp_smps[grp_smps$group=='lowISS',]$dataset,
                                              grp_smps[grp_smps$group=='healthy',]$dataset),
                                     'group', 'healthy'),
                  "highlow-healthy"=list(idents=c(grp_smps[grp_smps$group=='lowISS',]$dataset,
                                                  grp_smps[grp_smps$group=='highISS',]$dataset,
                                                  grp_smps[grp_smps$group=='healthy',]$dataset),
                                         'group', 'healthy'))

dpeaks_f <- file.path(outdir, "differential", "da_peaks.rds")
if((!file.exists(dpeaks_f)) | rerun_analysis){
  # Iterate through the group-comparisons given
  das_grps <- lapply(names(comp_grps), function(grp_id){
    print(paste(cg$idents, collapse=","))
    Idents(seu) <- 'dataset' # Extract the highISS and lowTSS
    seu_grp <- subset(seu, idents = cg$idents)
    Idents(seu_grp) <- anno_type
    
    # Call differential peaks per cell-type by the groups given
    das <- lapply(cell_types, function(ct){
      print(paste0(ct, "..."))
      
      cnts <- split(seu_grp@meta.data[[anno_type]] == ct, 
                    seu_grp@meta.data[[group_by]])
      if(!all(sapply(cnts, sum) > 50)){
        return(NULL)
      }
      
      # Differential accessibility test
      all_grps <- unique(seu_grp@meta.data[,cg[[2]]])
      grps_given <- cg[[3]]
      grps_inverse <- all_grps[!all_grps %in% cg[[3]]]
      da_peaks <- FindMarkers(object = seu_grp,
                              only.pos = F,
                              ident.1 = cg[[3]], 
                              group.by = cg[[2]], 
                              subset.ident=ct,
                              test.use = "LR", 
                              min.pct = 0.05,
                              latent.vars = 'nCount_ATAC')
      
      ## Split differential peaks into two lists, one for 
      # more-accessible in each group
      da_peaks_spl <- split(da_peaks, f=da_peaks$avg_log2FC >0)
      grps_map <- setNames(c(paste(grps_inverse, collapse="_"),
                             paste(grps_given, collapse="_")),
                           c('FALSE', 'TRUE'))
      names(da_peaks_spl) <- grps_map[names(da_peaks_spl)]
      return(da_peaks_spl)
    })
    names(das) <- cell_types
    return(das)
  })
  saveRDS(das_grps, file=dpeaks_f)
} else {
  das_grps <- readRDS(file=dpeaks_f)
}

#######################################################################
#### 7.d) Using the differential peaks to identify enriched motifs ####
dpeaks_f <- file.path(outdir, "differential", "da_peaks.rds")
if(!exists("seu")) seu <- readRDS(chromvar_obj_f)
if(!exists("das_grps")) das_grps <- readRDS(dpeaks_f)

dmotif_f <- file.path(outdir, "differential", "da_motifs.rds")
if((!file.exists(dmotif_f)) | rerun_analysis){
  # Iterate through groups
  dmotif_grps <- lapply(names(comp_grps), function(grp_id){
    cg <- comp_grps[[grp_id]]
    Idents(seu) <- 'dataset' # Extract the highISS and lowTSS
    seu_grp <- subset(seu, idents = cg$idents)
    Idents(seu_grp) <- anno_type
    
    # Iterate through cell types
    das <- lapply(cell_types, function(ct){
      print(paste0(ct, "..."))
      
      cnts <- split(seu_grp@meta.data[[anno_type]] == ct, 
                    seu_grp@meta.data[[group_by]])
      if(!all(sapply(cnts, sum) > 50)){
        return(NULL)
      }
      da_peaks_spl <- das_grps[[grp_id]][[ct]]
      
      ## Calculate the enriched motifs within the more-accessible peaks
      das_l <- lapply(names(da_peaks_spl), function(grp){
        da_peaks_i <- da_peaks_spl[[grp]]
        da_peaks_id <- rownames(da_peaks_i[which(da_peaks_i$p_val_adj < 0.05),])
        
        #isolating all accessible peaks in cell type: ct
        open.peaks <- AccessiblePeaks(seu_grp, idents = ct)
        
        # match the overall GC content in the peak set
        meta.feature <- GetAssayData(seu_grp, assay = "peaks", slot = "meta.features")
        enriched.motifs <- tryCatch({
          peaks.matched <- MatchRegionStats(meta.feature = meta.feature[open.peaks, ],
                                            query.feature = meta.feature[da_peaks_id, ],
                                            n = 50000)
          ## Matching GC.percent distribution to identify enrichment 
          FindMotifs(object = seu_grp,
                     features = da_peaks_id,
                     background=peaks.matched)
        }, error=function(e){NULL})
        
        return(list("posAccGrp"=grp, "da"=da_peaks_i, "motif"=enriched.motifs))
      })
      names(das_l) <- names(da_peaks_spl)
      return(das_l)
    })
    names(das) <- cell_types
    return(das)
  })
  names(dmotif_grps) <- names(comp_grps)
  saveRDS(dmotif_grps, file=dmotif_f)
} else {
  dmotif_grps <- readRDS(dmotif_f)
}

## Visualization
# Reduce differential motifs to a ggplot2 friendly format
# and return the enriched motifs within the differential peaks
n <- 20 #top n-features to extract
head_ct <-TRUE # extract the top n-features from each cell type rather than overall
da_motifs <- lapply(names(dmotif_grps), function(das_id){
  das <- dmotif_grps[[das_id]]
  grps <- unlist(sapply(das, function(i) sapply(i, function(j) j$posAccGrp)))
  grp_dir <- setNames(c(-1, 1), unique(sort(grps)))
  da <- lapply(names(das), function(id){
    i <- das[[id]]
    # Aggregate the bi-directional fold-enrichment and assign directions
    i_motifs <- lapply(i, function(i_dir){
      i_motif <- i_dir$motif
      i_motif$fold.enrichment <-  grp_dir[i_dir$posAccGrp] * i_motif$fold.enrichment
      i_motif$posAccGrp <- i_dir$posAccGrp
      return(i_motif)
    })
    i_motifs <- do.call(rbind, i_motifs)
    if(!is.null(i_motifs)) i_motifs$cell_type <- id
    i_motifs
  })
  
  da <- as.data.frame(do.call(rbind, da))
  # da <- da[which(da$observed >= 3),]
  da$q <- p.adjust(da$pvalue, method='bonferroni')
  da_filtq <- da[which(da$q < 0.05),]
  da_filt <- da_filtq[,c('motif.name', 'observed', 
                         'fold.enrichment', 'cell_type', "q")]
  if(head_ct){
    da_filt_l <- lapply(split(da_filt, f=da_filt$cell_type), head, n=n)
    da_filt <- as.data.frame(do.call(rbind, da_filt_l))
  } else {
    da_filt <- head(da_filt[order(da_filt$q),], n)
  }
  
  #sorted fold enrichment per motif
  desc_oi <- sort(sapply(split(da_filt[,'fold.enrichment'], 
                               f=da_filt[,'motif.name']), 
                         function(i) mean(abs(i))),
                  decreasing = T)
  
  # Remove motifs that are not to be plotted
  da_filt[,'motif.name'] <- factor(da_filt[,'motif.name'], 
                                   levels=unique(names(desc_oi)))
  # da_filt <- da_filt[-which(is.na(da_filt[,'motif.name'])),]
  da_filt$cell_type <- factor(da_filt$cell_type, levels=names(das))
  da_filt$grp <- names(grp_dir)[as.integer(da_filt$fold.enrichment > 0)+1]
  
  # Do the plotties
  pdf(file.path(outdir, "differential", paste0(das_id, "_da.pdf")), height=17, width = 25)
  # pdf("~/xfer/test3.pdf", height=14, width=25)
  b <- c(-4, -2, 0, 2, 4)
  gp <-ggplot(data=da_filt, aes(x=motif.name, y=fold.enrichment, 
                                fill=grp)) +
    geom_bar(stat='identity') + 
    geom_point(aes(x=motif.name, y=fold.enrichment, size=observed)) +
    facet_grid(cols=vars(cell_type), drop=F) +
    geom_hline(yintercept=0, size=1) +
    coord_flip() + 
    theme_minimal() + 
    scale_y_continuous(breaks=b, labels=abs(b), limits = c(min(b), max(b))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(gp)
  dev.off()
  
  return(da_filtq)
})

################################################################
#### 7.e) Annotate the genes within the differential peaks  ####
dpeaks_f <- file.path(outdir, "differential", "da_peaks.rds")
if(!exists("seu")) seu <- readRDS(chromvar_obj_f)
if(!exists("das_grps")) das_grps <- readRDS(dpeaks_f)

# Make TxDB object
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

# Creating mapping of Gene Ensembl to SYMBOLs
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

dpeak_anno_f <- file.path(outdir, "differential", "da_peak_anno.rds")
if(!file.exists(dpeak_anno_f)){
  das_anno_grps <- lapply(das_grps, function(da_grp){
    lapply(da_grp, function(da_ct){
      lapply(da_ct, function(da_i){
        da_i$chr <- gsub("-.*$", "", rownames(da_i))
        da_i$start <- gsub("^.*-(.*)-.*$", "\\1", rownames(da_i))
        da_i$end <- gsub("^.*-", "", rownames(da_i))
        da_gr <- makeGRangesFromDataFrame(da_i, keep.extra.columns = T)
        seqlevelsStyle(da_gr) <- 'NCBI'
        
        # Annotate the genomic features of differential accessible regions
        anno <- annotatePeak(da_gr, TxDb=txdb, tssRegion = c(-3000, 3000))
        anno_gr <- anno@anno
        anno_gr$symbol <- gene_ids[anno_gr$geneId]
        anno_gr$symbol[which(is.na(anno_gr$symbol))] <- anno_gr$geneId[which(is.na(anno_gr$symbol))]
        return(anno_gr)
      })
    })
  })
  saveRDS(das_anno_grps, file=dpeak_anno_f)
} else {
  das_anno_grps <- readRDS(dpeak_anno_f)
}

txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
gr <- GRanges(BSgenome.Hsapiens.UCSC.hg38@seqinfo)
gr <- keepStandardChromosomes(gr, pruning.mode = 'coarse')
seqlevels(gr) <- gsub("chr", "", seqlevels(gr))

TxDb <- ChIPseeker:::loadTxDb(txdb)
features <- ChIPseeker:::getGene(txdb, by = "transcript")
anno <- ChIPseeker:::getGenomicAnnotation(features, distance=0, tssRegion = c(-3000, 3000), TxDb=txdb,
                                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
anno <- annotatePeak(gr, TxDb=txdb, 
                     tssRegion = c(-3000, 3000))

################################
#### 8. CD4-subset analysis ####
cd4obj_f <- file.path(outdir, "seurat_obj", "sc_CD4_subset.rds")
if(!exists("seu")) seu <- readRDS(cd4obj_f)
seu$phenotype_detailed <- gsub(" ", "_", seu$phenotype_detailed)
seu$phenotype_detailed2 <- seu$phenotype_detailed
seu$phenotype_detailed2[which(seu$phenotype_detailed2 %in% c("CD4_N_mem", "CD4_TEMRA"))] <- 'CD4'

cd4outdir <- file.path(outdir, "CD4")
dir.create(cd4outdir, showWarnings = F)

# Make TxDB file from GTF
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
TxDb <- ChIPseeker:::loadTxDb(txdb)
features <- ChIPseeker:::getGene(txdb, by="gene")

GR <- transcripts(txdb)
PR <- promoters(txdb, upstream=2000, downstream=400)

# Add motif information from JASPAR2020
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", species = 'Homo sapiens'))
# Creating mapping of Gene Ensembl to SYMBOLs
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")
# Get mapping of JASPAR motif id to name
motif_names <- GetMotifData(object = seu, assay = 'peaks', 
                            slot = "motif.names")

# Subset based on the cells (idents) and test for differences between groups (grps)
comp_grps <- list("cd4s"=list(idents=c('CD4_TEM', 'CD4_N_mem', 'CD4_TEMRA'),
                              idents.id=c('phenotype_detailed'),
                              grps=c('highISS', 'lowISS', 'healthy'),
                              grps.id=c('group')),
                  "cd4n_temra"=list(idents=c('CD4'),
                                    idents.id=c('phenotype_detailed2'),
                                    grps=c('highISS', 'lowISS', 'healthy'),
                                    grps.id=c('group')))

### a) Differential accessible regions (DAR) between Hi-Lo ----
dir.create(file.path(cd4outdir, "dar"), showWarnings = F)
cd4_peaks_f <- file.path(cd4outdir, "dar", "da_peaks.rds")
rerun_analysis <- FALSE

# Subset based on the cells (idents) and test for differences between groups (grps)
comp_grps <- list("cd4s"=list(idents=c('CD4_TEM', 'CD4_N_mem', 'CD4_TEMRA'),
                              idents.id=c('phenotype_detailed'),
                              grps=c('highISS', 'lowISS', 'healthy'),
                              grps.id=c('group')),
                  "cd4n_temra"=list(idents=c('CD4'),
                                    idents.id=c('phenotype_detailed2'),
                                    grps=c('highISS', 'lowISS', 'healthy'),
                                    grps.id=c('group')))

if((!file.exists(cd4_peaks_f)) | rerun_analysis){
  # Iterate through the group-comparisons given
  das_grps <- lapply(names(comp_grps), function(grp_id){
    cg <- comp_grps[[grp_id]]
    print(paste(cg$idents, collapse=","))
    
    # Subset the seurat object based ona set of cells
    Idents(seu) <- cg$idents.id 
    cg_das <- lapply(cg$idents, function(cg_ident){
      print(paste0(cg_ident, "..."))
      # Subset for a specific cell-type
      seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
      
      # Differential accessibility test between every group
      grp_combn <- combn(cg$grps, m=2, simplify = T)
      colnames(grp_combn) <- apply(grp_combn, 2, paste, collapse="_")
      
      das <- apply(grp_combn, 2, function(grp){
        da_peaks <- FindMarkers(object = seu_grp,
                                only.pos = F,
                                ident.1 = grp[1], 
                                ident.2 = grp[2],
                                group.by = cg$grps.id, 
                                test.use = "LR", 
                                min.pct = 0.05,
                                logfc.threshold = 0.1,
                                latent.vars = 'nCount_peaks')
        return(da_peaks)
      })
      return(das)
    })
    names(cg_das) <- cg$idents
    return(cg_das)
  })
  das_grps <- unlist(das_grps, recursive = F)
  saveRDS(das_grps, file=cd4_peaks_f)
} else {
  das_grps <- readRDS(file=cd4_peaks_f)
}

### b) Genes & Motifs found within DARs ----
dir.create(file.path(cd4outdir, "dar"), showWarnings = F)
cd4_peaks_f <- file.path(cd4outdir, "dar", "da_peaks.rds")
das_grps <- readRDS(file=cd4_peaks_f)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
rerun_analysis <- FALSE

if((!file.exists(cd4_peaks_anno_f)) | rerun_analysis){
  # Annotate the genes and motifs within each DA peak
  das_anno <- lapply(das_grps, function(das){
    print(names(das))
    peak_motif <- lapply(das, function(da_i){
      print(".")
      if(is.null(da_i)) return(NULL)
      
      da_i$chr <- gsub("-.*$", "", rownames(da_i))
      da_i$start <- gsub("^.*-(.*)-.*$", "\\1", rownames(da_i))
      da_i$end <- gsub("^.*-", "", rownames(da_i))
      da_gr <- makeGRangesFromDataFrame(da_i, keep.extra.columns = T)
      seqlevelsStyle(da_gr) <- 'NCBI'
      
      # Annotate the genomic features of differential accessible regions
      anno <- annotatePeak(da_gr, TxDb=txdb, tssRegion = c(-3000, 3000),
                           level='gene', addFlankGeneInfo = TRUE, flankDistance = 500,
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", 
                                                         "Intron", "Downstream", "Intergenic"))
      anno_gr <- anno@anno
      anno_gr$symbol <- gene_ids[anno_gr$geneId]
      anno_gr$symbol[which(is.na(anno_gr$symbol))] <- anno_gr$geneId[which(is.na(anno_gr$symbol))]
      anno_gr$symbol_fl <- gene_ids[anno_gr$flank_geneIds]
      anno_gr$symbol_fl[which(is.na(anno_gr$symbol_fl))] <- anno_gr$flank_geneIds[which(is.na(anno_gr$symbol_fl))]
      
      # Annotate motifs for each fragment
      seqlevelsStyle(da_gr) <- 'UCSC'
      motif_ix <- matchMotifs(pfm, da_gr, genome = BSgenome.Hsapiens.UCSC.hg38,
                              out = "positions", p.cutoff = 5e-05) 
      motif_sc <- matchMotifs(pfm, da_gr, genome = BSgenome.Hsapiens.UCSC.hg38,
                              out = "scores", p.cutoff = 5e-05) 
      
      return(list("gene"=anno_gr, "motif_pos"=motif_ix, "motif_score"=motif_sc))
    })
  })
  
  # Add in the DA data to the gene and motif annotation
  das_grps <- lapply(names(das_anno), function(das_id){
    da_x <- lapply(names(das_anno[[das_id]]), function(da){
      x <- das_anno[[das_id]][[da]]
      x[['da']] <- das_grps[[das_id]][[da]]
      return(x)
    })
    names(da_x) <- names(das_anno[[das_id]])
    return(da_x)
  })
  names(das_grps) <- names(das_anno)
  
  saveRDS(das_grps, file=cd4_peaks_anno_f)
} else {
  das_grps <- readRDS(file=cd4_peaks_anno_f)
}

### c) TF motif enrichment across all Peaks (ChromVAR) ----
dir.create(file.path(cd4outdir, "tf"), showWarnings = F)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
das_grps <- readRDS(file=cd4_peaks_anno_f)
cd4_tf_f <- file.path(cd4outdir, "tf", "dtf.rds")
rerun_analysis <- FALSE

if((!file.exists(cd4_tf_f)) | rerun_analysis){
  # Iterate through the group-comparisons given
  dtfs_grps <- lapply(names(comp_grps), function(grp_id){
    cg <- comp_grps[[grp_id]]
    print(paste(cg$idents, collapse=","))
    
    # Subset the seurat object based ona set of cells
    Idents(seu) <- cg$idents.id 
    cg_dtfs <- lapply(setNames(cg$idents,cg$idents), function(cg_ident){
      print(paste0(cg_ident, "..."))
      # Subset for a specific cell-type
      seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
      DefaultAssay(seu_grp) <- 'chromvar'
      
      # Differential accessibility test between every group
      grp_combn <- combn(cg$grps, m=2, simplify = T)
      colnames(grp_combn) <- apply(grp_combn, 2, paste, collapse="_")
      
      dtfs <- apply(grp_combn, 2, function(grp){
        # motif_x <- names(motif_names[motif_names=='STAT1'])
        # mx <- GetAssayData(seu_grp)[motif_x,which(seu_grp@meta.data$group == grp[1])]
        # my <- GetAssayData(seu_grp)[motif_x,which(seu_grp@meta.data$group == grp[2])]
        # rbind(summary(mx), summary(my))
        
        dtf <- FindMarkers(object = seu_grp,
                           slot = "data",
                           only.pos = F,
                           ident.1 = grp[1], 
                           ident.2 = grp[2],
                           group.by = cg$grps.id, 
                           test.use = "wilcox", 
                           logfc.threshold=0.01,
                           min.pct = 0.05)
        dtf$MA <- rownames(dtf)
        dtf$motif <- motif_names[dtf$MA]
        return(dtf)
      })
      return(dtfs)
    })
    return(cg_dtfs)
  })
  dtfs_grps <- unlist(dtfs_grps, recursive = F)
  saveRDS(dtfs_grps, file=cd4_tf_f)
} else {
  dtfs_grps <- readRDS(file=cd4_tf_f)
}



# Idents(seu) <- 'phenotype_detailed' 
# cg_ident <- 'CD4_TEM'
# seu_grp <- subset(seu, idents = cg_ident) 
# DefaultAssay(seu_grp) <- 'chromvar'
# 
# # Differential accessibility test between every group
# grp <- c('highISS', 'lowISS')
# dtf <- FindMarkers(object = seu_grp,
#                    slot = "data",
#                    only.pos = F,
#                    ident.1 = grp[1], 
#                    ident.2 = grp[2],
#                    group.by = 'group', 
#                    test.use = "wilcox", 
#                    logfc.threshold = 0.01,
#                    min.pct = 0.05)
# dtf$MA <- rownames(dtf)
# dtf$motif <- motif_names[dtf$MA]
# dtf[which(dtf$motif %in% c('REL', 'FOS::JUNB')),]



### d) ChromVAR DARs: TF motif enrichment within DARs ----
dir.create(file.path(cd4outdir, "tf"), showWarnings = F)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
das_grps <- readRDS(file=cd4_peaks_anno_f)
cd4_tf_f <- file.path(cd4outdir, "tf", "dtf.rds")
chromvar_dtf <- readRDS(file=cd4_tf_f)
cd4_tf_dar_f <- file.path(cd4outdir, "tf", "dtf_dar.rds")
rerun_analysis <- FALSE

cg_ident <- 'CD4_TEM'
res_ident <- 'highISS_lowISS'

da <- das_grps[[cg_ident]][[res_ident]]$da
da_sig_spl <- split(da, f=da$p_val < 0.05)
names(da_sig_spl) <- setNames(c('nonDAR', 'DAR'), 
                              c('FALSE', 'TRUE'))[as.character(names(da_sig_spl))]

features <- list(chromvar_DAR=rownames(da_sig_spl$DAR),
                 chromvar_nonDAR=rownames(da_sig_spl$nonDAR))

## Run the customized chromVar function for calling peaks
for(assay_id in names(features)){
  print(assay_id)
  feature_i <- features[[assay_id]]
  
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 1))
  DefaultAssay(seu) <- 'peaks'
  seu <- RunChromVar_subset(
    object = seu,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    features=feature_i,
    new.assay.name=assay_id
  )
}

## Find differential motifs from chromVAR scores
dar_dtfs <- lapply(setNames(names(features),names(features)), 
                   function(assay_id){
                     seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
                     DefaultAssay(seu_grp) <- assay_id
                     
                     # Differential accessibility test between every group
                     grp_combn <- combn(cg$grps, m=2, simplify = T)
                     colnames(grp_combn) <- apply(grp_combn, 2, paste, collapse="_")
                     
                     grp <- strsplit(res_ident, split="_")[[1]]
                     dtf <- FindMarkers(object = seu_grp,
                                        slot = "data",
                                        only.pos = F,
                                        ident.1 = grp[1], 
                                        ident.2 = grp[2],
                                        group.by = cg$grps.id, 
                                        test.use = "wilcox", 
                                        logfc.threshold=0.1,
                                        min.pct = 0.05)
                     dtf$MA <- rownames(dtf)
                     dtf$motif <- motif_names[dtf$MA]
                     return(dtf)
                   })
saveRDS(dar_dtfs, file=cd4_tf_dar_f)

dar_dtfs <- readRDS(file=cd4_tf_dar_f)



cg_ident <- 'CD4_TEM'
res_ident <- 'highISS_lowISS'
res_ids <- strsplit(res_ident, split = "_")[[1]]
cells_split <- split(colnames(seu), f=list(seu$group, seu$phenotype_detailed))
cells_grp1 <- cells_split[[paste0(res_ids[1], ".", cg_ident)]]
cells_grp2 <- cells_split[[paste0(res_ids[2], ".", cg_ident)]]

motif_ids <- data.frame(id=rownames(seu_assay),
                        name=as.character(motif_names[rownames(seu_assay)]))
motif_i <- motif_ids$name


m_chromvars <- lapply(chromvar_ids, function(chromvar_id){
  DefaultAssay(seu) <- chromvar_id #'chromvar'
  seu_assay <- GetAssayData(seu)
  
  m_chromvars <- lapply(motif_i, function(mi){
    motif_id_i <- motif_ids[motif_ids$name == mi,]$id
    if(length(motif_id_i)==0) return(NULL)
    tval <- t.test(x=seu_assay[motif_id_i, cells_grp1],
                   y=seu_assay[motif_id_i, cells_grp2])
    m_chromvar <- data.frame("stat"=tval$statistic,
                             "p"=tval$p.value,
                             "grp1"=mean(seu_assay[motif_id_i, cells_grp1]),
                             "grp2"=mean(seu_assay[motif_id_i, cells_grp2]),
                             "motif"=mi)
    # print(paste0(chromvar_id, " - ", mi, ": ", w_pval))
    return(m_chromvar)
  })
  m_chromvars_seu <- as.data.frame(do.call(rbind, m_chromvars))
  m_chromvars_seu$chromvar_obj <- chromvar_id
  m_chromvars_seu$pval_adj <- p.adjust(m_chromvars_seu$p, method='fdr')
  m_chromvars_seu <- m_chromvars_seu[order(m_chromvars_seu$pval_adj),]
  return(m_chromvars_seu)
})

df <- as.data.frame(do.call(rbind, m_chromvars))
df$pval_adj <- -1*log10(as.numeric(df$pval_adj))
sig_subset <- df %>% 
  filter(pval_adj > (-1*log10(0.1))) %>%
  group_by(chromvar_obj) %>%
  slice_max(order_by = pval_adj, n = 50)
sig_subset2 <- df %>%
  filter(pval_adj > (-1*log10(0.1))) %>%
  mutate(direction=stat>0) %>%
  group_by(direction) %>%
  slice_max(order_by = pval_adj, n = 25)

sig_motifs <- unique(c(sig_subset$motif))
sig_motifs <- unique(c(sig_subset$motif, sig_subset2$motif))
sig_df <- df[which(df$motif %in% sig_motifs),]
sig_df$sig <- sig_df$pval_adj > (-1*log10(0.1))
sig_df$motif <- factor(sig_df$motif, levels=sig_motifs)
sig_df$chromvar_obj <- factor(sig_df$chromvar_obj,
                              levels=paste0("chromvar", c("", "_DAR", "_nonDAR")))

pdf("~/xfer/test2.pdf", height = 12)
ggplot(sig_df, aes(y=motif, x=stat, fill=chromvar_obj, group=chromvar_obj)) +
  geom_bar(stat='identity', position='dodge') +
  geom_bar(aes(x=stat, y=motif, group=chromvar_obj, alpha=sig), 
           stat='identity', position='dodge', fill='grey') +
  scale_alpha_discrete(range = c(1, 0)) +
  xlab("t-statistic") + 
  theme_classic()
dev.off()

### e) Correlation between chromVAR DMs ---------------
cd4obj_chromvar_f <- file.path(outdir, "seurat_obj", "sc_chromvar_CD4_subset.rds")
dir.create(file.path(cd4outdir, "tf"), showWarnings = F)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
das_grps <- readRDS(file=cd4_peaks_anno_f)
cd4_tf_f <- file.path(cd4outdir, "tf", "dtf.rds")
chromvar_dtf <- readRDS(file=cd4_tf_f)
cd4_tf_dar_f <- file.path(cd4outdir, "tf", "dtf_dar.rds")
rerun_analysis <- FALSE

DefaultAssay(seu) <- 'peaks'
peaks_mat <- GetAssayData(seu)
motifs_mat <- Motifs(seu)@data
DefaultAssay(seu) <- 'SCT'
expr_mat <- GetAssayData(seu)

seu <- SyncExprAndPeaks(
  object=seu, 
  peaks_mat=peaks_mat, 
  motifs_mat=motifs_mat,
  expr_mat=expr_mat,
  txdb=makeTxDbFromGFF(file = gtf_file, format = "gtf"),
  new.assay.id='pr_expr_peaks'
)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 1))
DefaultAssay(seu) <- 'pr_expr_peaks'
seu <- RunChromVAR_subset(
  object = seu,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix=f_motifs_mat,
  new.assay.name='chromvar_expr'
)

## ChromVAR Matrix [ motif x cells ] (chromvar_mat)
DefaultAssay(seu) <- 'chromvar'
peak_chromvar_mat <- GetAssayData(seu)

DefaultAssay(seu) <- 'chromvar_expr'
expr_chromvar_mat <- GetAssayData(seu)

cg_ident <- 'CD4_TEM'
res_ident <- 'highISS_lowISS'
res_ids <- strsplit(res_ident, split = "_")[[1]]
cells_split <- split(colnames(seu), f=list(seu$group, seu$phenotype_detailed))

motif_i <- rownames(peak_chromvar_mat)
m_chromvars <- lapply(res_ids, function(grp_id){
  cells_grp <- cells_split[[paste0(grp_id, ".", cg_ident)]]
  
  m_chromvars <- lapply(motif_i, function(mi){
    corval <- cor.test(x=peak_chromvar_mat[mi, cells_grp],
                       y=expr_chromvar_mat[mi, cells_grp], 
                       method='spearman')
    corval$estimate
    corval$p.value
    m_chromvar <- data.frame("r"=corval$estimate,
                             "p"=corval$p.value,
                             "peak_mean"=mean(peak_chromvar_mat[mi, cells_grp]),
                             "expr_mean"=mean(expr_chromvar_mat[mi, cells_grp]),
                             "motif"=as.character(motif_names[mi]))
    # print(paste0(chromvar_id, " - ", mi, ": ", w_pval))
    return(m_chromvar)
  })
  m_chromvars_seu <- as.data.frame(do.call(rbind, m_chromvars))
  m_chromvars_seu$group <- grp_id
  m_chromvars_seu$pval_adj <- p.adjust(m_chromvars_seu$p, method='fdr')
  m_chromvars_seu <- m_chromvars_seu[order(m_chromvars_seu$pval_adj),]
  
  return(m_chromvars_seu)
})
m_chromvar <- as.data.frame(do.call(rbind, m_chromvars))
m_chromvar$sig <- m_chromvar$pval_adj < 0.01
m_chromvar$log10q <- -1*log10(m_chromvar$pval_adj)

pdf("~/xfer/test3.pdf")
ggplot(m_chromvar, aes(x=r, y=log10q, color=group)) +
  geom_point() +
  geom_point(aes(x=r, y=log10q, alpha=sig), color='grey') +
  scale_alpha_discrete(range = c(1, 0)) +
  geom_text_repel(data=m_chromvar[which(m_chromvar$sig),], 
                  aes(x=r, y=log10q, label=motif), 
                  size=3, min.segment.length = 0, box.padding = 0.1) +
  theme_classic() +
  xlim(-0.25, 0.25) +
  xlab("Spearman correlation") + ylab("-log10(q)")
dev.off()
write.table(m_chromvar, file="~/xfer/test.csv",
            sep=",", quote=F, col.names = T, row.names = F)

saveRDS(seu, file=cd4obj_chromvar_f)



### f) Motif-specific TSS Plot between High and LowISS groups ----
motif_tags_f <- file.path(cd4outdir, "tf", "motifs_tssplot.rds")
rerun_analysis <- FALSE

if((!file.exists(motif_tags_f)) | rerun_analysis){
  DefaultAssay(seu) <- 'peaks'
  peaks_mat <- GetAssayData(seu)
  motifs_mat <- Motifs(seu)@data
  colnames(motifs_mat) <- motif_names[colnames(motifs_mat)]
  motif_tags <- getMotifTags(peaks_mat, motifs_mat, txdb, 
                             grps_col='group', grps=c('highISS', 'lowISS'),
                             return_gg=TRUE)
  saveRDS(motif_tags, file=motif_tags_f)
} else {
  motif_tags <- readRDS(motif_tags_f)
}


# pdf("~/xfer/test.pdf")
# lapply(motif_tags[c('STAT1', 'REL')], function(motif_gg){
#   plot_grid(plotlist = motif_gg, ncol=2)
# })
# dev.off()

### g) Investigate enrichment between DEGs and motif-specific peaks ----
# Given a set of DEG [g], we can create a 
# boolean matrix of motifs M [genes x motifs] focused only on whether there is 
# a peaks for motif_i in the promoter of g_x. We can then do something like 
# a propprtions test to compare the ratio of genes with motif_i to 
# genes not with motif_i
deg_outpath <- file.path(cd4outdir, "deg_motif")
dir.create(deg_outpath, showWarnings = F)

degs_f <- file.path(deg_outpath, "degs.rds")
peaks_f <- file.path(deg_outpath, "peaks.rds")
out_stats_f <- file.path(deg_outpath, "stats_p.csv")
rerun_analysis <- FALSE

## Differential expression of CD4 TEM cells high-lowISS group
if((!file.exists(degs_f)) | rerun_analysis){
  ## Get differentially expressed genes using DESeq2
  cg_ident <- 'CD4_TEM'  
  grp <- c('highISS', 'lowISS')
  group_id <- 'group'
  Idents(seu) <- 'phenotype_detailed'
  seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
  DefaultAssay(seu_grp) <- 'RNA'
  
  deg <- FindMarkers(object = seu_grp,
                     only.pos = F,
                     ident.1 = grp[1], 
                     ident.2 = grp[2],
                     group.by = group_id, 
                     test.use = "DESeq2", 
                     min.pct = 0.05)
  deg$gene <- rownames(deg)
  
  saveRDS(deg, file=degs_f)
} else {
  deg <- readRDS(degs_f)
}

## Annotation of peaks based on motifs, selected for genes with peaks in promoters
if((!file.exists(peaks_f)) | rerun_analysis){
  ## Annotate and grab all the peaks for motif_i that are in the promoters
  DefaultAssay(seu) <- 'peaks'
  peaks_mat <- GetAssayData(seu)
  motifs_mat <- Motifs(seu)@data
  colnames(motifs_mat) <- motif_names[colnames(motifs_mat)]
  
  # Transform peaks into GRanges object
  peaks_id <- do.call(rbind, strsplit(rownames(peaks_mat), split="-"))
  peaks_gr <- GRanges(seqnames=peaks_id[,1], 
                      ranges=IRanges(start = as.integer(peaks_id[,2]), 
                                     end = as.integer(peaks_id[,3])))
  # Get a peak-list matched index of motifs
  motif_idx <- apply((motifs_mat!=0), 2, which)
  
  # Go motif by motif to find the tag-matrix
  motif_annos <- lapply(names(motif_idx), function(motif_name){
    # separate out peaks that contain motif_i
    print(paste0(motif_name, " (", match(motif_name, names(motif_idx)), 
                 "/", length(motif_idx), ")"))
    motif_i <- motif_idx[[motif_name]]
    peaks_i <- peaks_gr[motif_i,]
    seqlevelsStyle(peaks_i) <- 'NCBI'
    peak_i_anno <- annotatePeak(peaks_i, tssRegion=c((-1*promoter_flank), promoter_flank),
                                TxDb=txdb, level="gene")
    promoter_i_anno <- peak_i_anno@anno %>%
      as.data.frame() %>%
      filter(grepl("Promoter", annotation))
    return(list("promoter"=promoter_i_anno, "peak"=peak_i_anno))
  })
  names(motif_annos) <- names(motif_idx)
  
  saveRDS(motif_annos, file=peaks_f)
} else {
  deg <- readRDS(peaks_f)
}

# Select only the promoter-subsetted peaks
prom_annos <- lapply(motif_annos, function(i) {
  p <- i$promoter
  p$gene <- gene_ids[p$geneId]
  return(p)
})

# Crappy hack together dataframe of genes by DEG and peaks in Motif_i
deg_genes <- deg[which(deg$p_val_adj < 0.05),]$gene
non_deg_genes <- deg[which(deg$p_val_adj > 0.05),]$gene
all_motif_genes <- unique(sort(unlist(sapply(prom_annos, function(p) p$gene))))

gene_df <- data.frame("gene"=unique(sort(c(deg_genes, non_deg_genes, all_motif_genes))))
gene_df$deg <- FALSE
gene_df$deg[which(gene_df$gene %in% deg_genes)] <- 'DEG'
gene_df$deg[which(gene_df$gene %in% non_deg_genes)] <- 'nonDEG'
gene_df$motif <- FALSE
gene_df$motif[which(gene_df$gene %in% all_motif_genes)] <- 'allMotif'

stats_l <- lapply(prom_annos, function(p){
  motif_genes <- p$gene
  gene_df$motif[which(gene_df$gene %in% motif_genes)] <- 'motif'
  contigency_mat <- table(gene_df[,c('deg', 'motif')])[-2, -2]
  p <- fisher.test(contigency_mat)$p.val
  residuals <- chisq.test(contigency_mat)$residuals
  stdres <- chisq.test(contigency_mat)$stdres
  
  return(list('p'=p, 'residuals'=residuals))
})

stats_df <- data.frame('p'=sapply(stats_l, function(i) i$p)) %>%
  mutate(padj = p.adjust(p, method='fdr')) %>% 
  tibble::rownames_to_column("Gene") %>%
  mutate(residual_deg_motif = sapply(stats_l, function(i) i$residuals['DEG', 'motif'])) %>%
  dplyr::arrange(padj)
write.table(stats_df, file=out_stats_f, sep=",",
            col.names = T, row.names = F, quote = F)



