renv::load("~/downloads/renvs")
## sara Tumor/LN KO and WT samples
library(AUCell)
library(monocle3)
library(parsnip)
# library(slingshot)
# library(tradeSeq)
library(tidyr)
library(msigdbr)
library(BiocParallel)
library(GenomicFeatures)
# library(harmony)
library(org.Mm.eg.db)
library(SingleR)
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
# library(cellpypes)
library(SeuratWrappers)
# library(velocyto.R)
library(ProjecTILs)
library(SCENIC)
library(SCopeLoomR)

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

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor'
setwd(PDIR)

gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
groups <- list.files(datadir)

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
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")
gm <- geneMap(species='Mus musculus')

# Last updated Dec-5-2022; Still works as of Dec-5-2022
recode_map <- function(x){
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

#### Other Functions ####

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

###########################################
#### 1.b Cell type annotation - ImmGen ####
if(!exists("seu"))seu <- readRDS(file=file.path(datadir, "seurat_obj", "2_seu_integratedX.rds"))
seu$ilc2 <- ifelse(seu$immgen.fine.cell == 'ILC (ILC2)', 'ILC2', 'nonILC2')


pdf("~/xfer/test.pdf", width = 9, height = 6)
DimPlot(seu, group.by=c('ilc2', 'manual_anno'))
FeaturePlot(seu, )
dev.off()

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
  plot_grid(dpcl, dp1, ncol=2)
  plot_grid(dpcl, dp1b, ncol=2)
  plot_grid(dpcl, dp1c, ncol=2)
  plot_grid(dpcl, dp1d, ncol=2)
  plot_grid(dp1c, dp2, ncol=2)
  dev.off()
}


saveRDS(seu, file=file.path(datadir, "seurat_obj", "2_seu_annoX.rds"))
##########################################
#### 1.c Cell type annotation - Other ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
# if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))

sctype_dir <- '/cluster/projects/mcgahalab/bin/sc-type/R'
sctype_db <- 'ScTypeDB_full.xlsx'
sctype_tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 


#--- i) ScType -------------------------------------------------------
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

#--- ii) scSorter ----------------------------------------------------
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

DefaultAssay(seu) <- 'RNA'
DefaultAssay(refobj_custom) <- 'RNA'
seu <- ProjecTILs.classifier(seu, refobj_custom, ncores = 1, 
                              split.by = "orig.ident", filter.cells=F)
# seu$functional.cluster.tit <- seu$functional.cluster
# seu$functional.cluster.cd4.2 <- seu$functional.cluster
# seu$functional.cluster.cd8.2 <- seu$functional.cluster
# seu$functional.cluster.dc.2 <- seu$functional.cluster
seu$functional.cluster.custom <- seu$functional.cluster
seu$functional.cluster.treg <- seu$functional.cluster
pdf("~/xfer/test9.pdf", width = 16, height = 12)
dp1 <- DimPlot(seu, group.by='seurat_clusters', raster=T, label = T)
dp2 <- DimPlot(seu, group.by='functional.cluster.custom', raster=T, label = T)
dp1 + dp2 

Idents(seu) <- 'functional.cluster.custom'
dp_clusters <- lapply(unique(seu$functional.cluster.custom), function(clid){
  scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = "red",
                         background_color = "lightgray", raster=T) + 
    NoLegend() +
    ggtitle(clid)
})
cowplot::plot_grid(plotlist = dp_clusters, ncol=7)

Idents(seu) <- 'functional.cluster.treg'
dp_clusters <- lapply(unique(seu$functional.cluster.treg), function(clid){
  scCustomize::Cluster_Highlight_Plot(seu, cluster_name = clid, highlight_color = "red",
                                      background_color = "lightgray", raster=T) + 
    NoLegend() +
    ggtitle(clid)
})
cowplot::plot_grid(plotlist = dp_clusters, ncol=7)

dev.off()
# table(seu$functional.cluster.custom2, seu$seurat_clusters)
# DimPlot(seu, group.by='functional.cluster.momac', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.cd4', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.cd8', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.cd4.2', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.dc.2', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.dc', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.tit', raster=T, label = T)
# DimPlot(seu, group.by='functional.cluster.custom', raster=T, label = T)
# dev.off()

saveRDS(seu, file=file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))

#################################################################
#### 1.d Trajectory analysis of integrated dataset - MONOCLE3 ####
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

#################################################################
#### 2. Recluster Tumor and TDLN independently and integrate ####
## Take the seurat object and split them into a tumor and TDLN subset,
# then reprocess them all indepedently, recluster, revisualize, and then
# integrate them together using MNN
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))

seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()
DefaultAssay(seu) <- 'RNA'

seu$splitid <- gsub("_.*", "", seu$orig.ident)
Idents(seu) <- 'splitid'
seul <- lapply(unique(seu$splitid), function(sid){
  seu_sub <- subset(seu, ident=sid)
  seux <- preprocessSeu(seu_sub, ncount_min=0, ncount_max=Inf,
                        nfeature_min=0, nfeature_max=Inf,
                        mt_max=100, org='mouse', numpcs=30, getPCs=FALSE,
                        variable_features=2000, res=0.8)
  seux[['integrated']] <- NULL
  seux@assays$RNA@scale.data <- as.matrix(0)
  seux@assays$SCT@scale.data <- as.matrix(0)
  return(seux)
})

seu.reint <- RunFastMNN(object.list = seul, 
                        features = 2000, assay='SCT')
seu.reint <- RunUMAP(seu.reint, reduction = "mnn", dims = 1:30, n.neighbors=30L,
                     min.dist = 0.1, return.model=TRUE)
seu.reint <- FindNeighbors(seu.reint, reduction = "mnn", dims = 1:30)
seu.reint <- FindClusters(seu.reint, resolution = 0.9, graph.name='SCT_snn')
DefaultAssay(seu.reint) <- 'RNA'
seu.reint <- NormalizeData(seu.reint,
                           normalization.method = "LogNormalize") %>%
  FindVariableFeatures(., selection.method = "vst",
                       nfeatures = 3000, verbose = FALSE) %>% 
  ScaleData(.)
seu.reint <- RunPCA(seu.reint, features=VariableFeatures(seu.reint), ndims.print=1:30,
                    nfeatures.print=5, npcs=30)
saveRDS(seu.reint, file=file.path(datadir, "seurat_obj", "3_seu_split.integrated.rds"))
saveRDS(seul, file=file.path(datadir, "seurat_obj", "3_seu_split.separate.rds"))


if(do.plot){
  dp1a <- DimPlot(seul[[1]], group.by='seurat_clusters', label=T, raster=T) + NoLegend() + ggtitle("LN")
  dp1b <- DimPlot(seul[[1]], group.by='manual_anno', label=T, raster=T) + NoLegend()
  dp1c <- DimPlot(seul[[1]], group.by='functional.cluster.custom', label=T, raster=T) + NoLegend()
  dp2a <- DimPlot(seul[[2]], group.by='seurat_clusters', label=T, raster=T) + NoLegend() + ggtitle("Tumor")
  dp2b <- DimPlot(seul[[2]], group.by='manual_anno', label=T, raster=T) + NoLegend()
  dp2c <- DimPlot(seul[[2]], group.by='functional.cluster.custom', label=T, raster=T) + NoLegend()
  dpinta <- DimPlot(seu.reint, group.by='seurat_clusters', reduction='umap', label=T,raster=T) + NoLegend() + ggtitle("LN_Tumor_Integrated")
  dpintb <- DimPlot(seu.reint, group.by='manual_anno', reduction='umap', label=T,raster=T) + NoLegend()
  dpintc <- DimPlot(seu.reint, group.by='functional.cluster.custom', reduction='umap', label=T,raster=T) + NoLegend()
  
  pdf("~/xfer/seu.split_and_reintegrate.pdf", width=16)
  dp1a + dp1b + dp1c
  dp2a + dp2b + dp2c
  dpinta + dpintb + dpintc
  dev.off()
}


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


files <- list.files()
xdf <- lapply(files, function(i){
  f <- list.files(file.path(i))
  df <- read.table(file.path(i, f))
  return(df)
})
xdf2 <- purrr::reduce(xdf, left_join, by='V1') %>% 
  tibble::column_to_rownames(., "V1") %>% 
  rename_with(., ~files)
write.table(round(xdf2, 2), file=file.path("~/xfer/", "robbie.paad.csv"),
            sep=",", quote = F, col.names = T, row.names = T)

##############################################
#### 4.a Differential expression analysis ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
if(!exists("cds"))  cds <- readRDS(file = file.path(datadir, "seurat_obj", "3_cds_degPrepX.rds"))
subset_of_seu <- TRUE
outdirdeg <- file.path(outdir, "degs")
dir.create(outdirdeg, recursive = F)
pval_sig <- 0.05
fc_sig <- 0.5
cols=c("KO.down_WT.down" = "#fdb462", 
       "KO.down_WT.up" = "#80b1d3", 
       "KO.up_WT.up" = "#fdb462", 
       "KO.up_WT.down" = "#fb8072")

seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()

seu$tissue_cond <- gsub("_72h|_PBS" , "", seu$orig.ident)
seu$kowt_comp <- gsub("WT_|KO_", "", seu$orig.ident)
seu.list <- SplitObject(seu, split.by = "tissue_cond") # Compares 72h timepoint to BL
seu.list2 <- SplitObject(seu, split.by = "kowt_comp")  # Compares KO to WT at each timepoint
seu.lists <- c(seu.list2, seu.list)

#--- a) Differential expression between groups across all clusters ----
## Differential expression per cluster (seurat_clusters) or celltype (manual_anno)
anno_ident <- 'seurat_clusters' #'manual_anno', 'seurat_clusters'
ct_markers_all <- lapply(names(seu.lists), function(seu_id){
  seu_i <- seu.lists[[seu_id]]
  # Idents(seu_i) <- 'monocle3_partitions'
  # seu_i <- subset(seu_i, idents='1')
  
  colid <- c('tissue_cond', 'kowt_comp')
  colid <- colid[grep(seu_id, seu_i@meta.data[,colid])]
  idents <- switch(colid,
                   kowt_comp={
                     id <- strsplit(seu_id, split="_")[[1]]
                     c(paste(c(id[1], "KO", id[2]), collapse="_"), 
                       paste(c(id[1], "WT", id[2]), collapse="_"))
                   },
                   tissue_cond={
                     c(paste0(seu_id, "_72h"),  paste0(seu_id, "_PBS"))
                   })
  
  # DEG across all cells (not cluster specific)
  Idents(seu_i) <- 'orig.ident'
  marker_partition <- FindMarkers(seu_i, assay = "RNA", test.use='wilcox',
                                  ident.1= idents[1], ident.2= idents[2],
                                  verbose = FALSE,
                                  logfc.threshold = 0,
                                  recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
    tibble::rownames_to_column(., "symbol")
  
  # DEG between groups within individual clusters  (e.g. LN_KO_72h vs LN_WT_72h in Cluster12)
  Idents(seu_i) <- anno_ident
  celltypes <- (table(Idents(seu_i)) > 50) %>% which %>% names
  celltypes <- c('all', celltypes)
  
  ct_markers <- lapply(setNames(celltypes, celltypes), function(ct){
    print(paste0(ct, ": ", paste(idents, collapse=", ")))
    seu_j <- if(ct=='all') seu_i else subset(seu_i, idents=ct)
    Idents(seu_j) <- 'orig.ident'
    ct_marker <- tryCatch({
      FindMarkers(seu_j, assay = "RNA", test.use='wilcox',
                  ident.1= idents[1], ident.2= idents[2],
                  verbose = FALSE,
                  logfc.threshold = 0,
                  recorrect_umi =  if(subset_of_seu) FALSE else TRUE) %>% 
        tibble::rownames_to_column(., "symbol")
    }, error=function(e){NULL})
    return(ct_marker)
  })
  return(ct_markers)
})
names(ct_markers_all) <- names(seu.lists)
saveRDS(ct_markers_all, file=file.path(outdirdeg, "celltypes_degs.rds"))
saveRDS(ct_markers_all2, file=file.path(outdirdeg, "clusters_degs.rds"))

if(!exists("ct_markers_all2")) ct_markers_all2 <- readRDS(file.path(outdirdeg, "clusters_degs.rds"))
if(!exists("ct_markers_all")) ct_markers_all <- readRDS(file.path(outdirdeg, "celltypes_degs.rds"))


# pdf("~/xfer/sara_volcano.cluster.pdf", height = 7, width = 15)
ct_markers_X <- ct_markers_all
all_cts <- sapply(ct_markers_X, function(i) names(i)) %>%
  unlist %>% as.character %>% unique
# marker_partition


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

volc_comps <- lapply(all_cts, function(ct_i){
  print(ct_i)
  volc_comp <- lapply(names(ct_markers_X), function(comp_id){
    i <- tryCatch({
      class_id <- gsub("^.*_", "", comp_id)
      idx <- which(seu.lists[[comp_id]]$manual_anno == ct_i)
      n  <- table(seu.lists[[comp_id]][,idx]$orig.ident) %>%
        as.data.frame 
      ids <- do.call(rbind, strsplit(as.character(n$Var1), split="_"))
      n <- setNames(n$Freq, c(setdiff(ids[1,], ids[2,]),
                              setdiff(ids[2,], ids[1,])))
      ct_markers_X[[comp_id]][[ct_i]] %>%
        mutate(sig=case_when((avg_log2FC >= fc_sig) & (p_val_adj <= qval_sig) ~ "up",
                             (avg_log2FC <= (-1*fc_sig)) & (p_val_adj <= qval_sig) ~ "down",
                             TRUE ~ "ns"),
               log10p=(-1*log10(p_val_adj)),
               group=comp_id,
               n=paste(paste0("n_", names(n), "=", n), collapse=";")) %>%
        select(symbol, avg_log2FC, sig, log10p, n) %>%
        rename_with(., ~paste0(., ".", class_id), .cols = -symbol)
    }, error=function(e){NULL})
    return(i)
  })
  names(volc_comp) <- names(ct_markers_X)
  if(any(sapply(volc_comp[c("LN_KO", "LN_WT", 
                            "Tumor_KO", "Tumor_WT")], is.null))){
    print("null")
    return(NULL)
  }
  
  ln_lfc <- volc_comp[c('LN_KO', 'LN_WT')] %>% 
    purrr::reduce(full_join, by=c('symbol')) %>% 
    mutate(sig=paste0("KO.", sig.KO, "_WT.", sig.WT),
           celltype="LN") %>%
    select(-c(sig.KO, sig.WT)) 
  
  tumor_lfc <- volc_comp[c('Tumor_KO', 'Tumor_WT')] %>% 
    purrr::reduce(full_join, by=c('symbol')) %>% 
    mutate(sig=paste0("KO.", sig.KO, "_WT.", sig.WT),
           celltype="Tumor") %>%
    select(-c(sig.KO, sig.WT)) 
  
  ggscatter <- function(df, cols=NULL){
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
      ggtitle(ct_i) + 
      NoLegend()
  }
  gg_scatter <- ggscatter(ln_lfc) + ggscatter(tumor_lfc)
  gg_scatter_ln <- ggscatter(ln_lfc)
  gg_scatter_tumor <- ggscatter(tumor_lfc)
  
  gg_volcano <- ggplot(rbind(ln_lfc, tumor_lfc) %>% 
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
})
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


.printSymbolOv(volc_comps[['Tregs; Intermediate']]$ln_lfc)
.printSymbolOv(volc_comps[['Tregs; Cxcl10_hi']]$ln_lfc)

#--- b) Differnetial expression of Cluster 18 and 12 ----
Idents(seu) <- 'seurat_clusters'

.getTopGenes <- function(deg, n=100, dir="up", verbose=T){
  xsel <- deg %>% 
    arrange(p_val_adj) %>%
    filter(if(dir=='up') avg_log2FC >=0 else avg_log2FC < 0) %>%
    filter(p_val_adj <= 0.05) %>%
    slice_head(., n=n) %>% 
    select(symbol)
  if(verbose) cat(paste(c(xsel[,1], ""), collapse="\n"))
  return(xsel)
}

# Gets the DEGs between cluster X and all other cells
clusters_test <- c('18', '12')
mps <- lapply(setNames(clusters_test,clusters_test), function(clid){
  FindMarkers(seu, assay = "RNA", test.use='wilcox',
              ident.1= clid, verbose = FALSE,
              recorrect_umi = FALSE) %>% 
    tibble::rownames_to_column(., "symbol")
})
mps_topn <- lapply(names(mps), function(mp_id){
  # Writes out the top 100 upregulated and downregulated genes for cluster X
  print(paste0(mp_id, " - up..."))
  up <- .getTopGenes(mps[[mp_id]], n=100, dir='up', verbose=F)
  print(paste0(mp_id, " - down..."))
  dn <- .getTopGenes(mps[[mp_id]], n=100, dir='down', verbose=F)
  return(list("up"=up, "down"=dn))
})
names(mps_topn) <- names(mps)
mps_topn_mat <- do.call(cbind, unlist(mps_topn, recursive=F))
colnames(mps_topn_mat) <- names(unlist(mps_topn, recursive=F))

write.table(mps_topn_mat, file=file.path(outdirdeg, "cluster.18_12.topn.csv"),
            sep=",", quote = F, col.names = T, row.names = F)


# Finds the DEGs between Cluster 12 and 18, 8, and 21+10
refcl <- '12'
compclusters <- list('18'='18', '8'='8', '21_10'=c('21', '10'))
comp_mps <- lapply(compclusters, function(compid){
  FindMarkers(seu, assay = "RNA", test.use='wilcox',
                ident.1= refcl, ident.2=compid, verbose = FALSE,
                recorrect_umi = FALSE) %>% 
    tibble::rownames_to_column(., "symbol")
})
names(comp_mps) <- names(compclusters)
comp_topn <- lapply(names(comp_mps), function(mp_id){
  # Writes out the top 100 upregulated and downregulated genes for cluster X
  print(paste0(mp_id, " - up..."))
  up <- .getTopGenes(comp_mps[[mp_id]], n=100, dir='up', verbose=T)
  print(paste0(mp_id, " - down..."))
  dn <- .getTopGenes(comp_mps[[mp_id]], n=100, dir='down', verbose=T)
  return(list("up"=up, "down"=dn))
})
names(comp_topn) <- names(comp_mps)
comp_topn_mat <- do.call(cbind, unlist(comp_topn, recursive = F))
colnames(comp_topn_mat) <- names(unlist(comp_topn, recursive = F))

write.table(comp_topn_mat, file=file.path(outdirdeg, "cluster.12_vs_18_8_2110.topn.csv"),
            sep=",", quote = F, col.names = T, row.names = F)

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
celltype='Tregs; Cycling' # names(gsea_markers$celltypes$LN_KO)
celltypes <- grep("Treg", names(gsea_markers$celltypes$LN_KO), value=T)
  
## Create GSEA data with the delta-NES between KO and WT, and calculate p-values
dNES_ct_dat <- lapply(setNames(celltypes,celltypes), function(celltype){
  lapply(setNames(databases,databases), function(db){
    # Extract the original NES data for KO and WT
    celltypeid='LN_KO' # 'Tumor_KO'  # just used to initialize, not important
    dat <- as.data.frame(gsea_markers$celltypes[[celltypeid]][[celltype]][[db]]) %>% 
      select(-c(NES, pvalue, p.adjust,qvalues))
    databaseid <- switch(db,
                         "H.base"='HALLMARK',
                         "C2.CP:REACTOME"='REACTOME',
                         "C5.GO:BP"='GOBP',
                         "C5.GO:CC"='GOCC',
                         "C5.GO:MF"='GOMF')
    
    # Extract the database-paired delta-NES between KO-and-WT
    dnes <- dgsea$celltypes$LN[[celltype]] %>%
      filter(Database == databaseid)
    
    # Monte Carlo randomization of NES values to generate a null deltaNES distribution
    set.seed(1234)
    ridx1 <- sample(1:nrow(dnes), size=100000, replace = T)
    ridx2 <- sample(1:nrow(dnes), size=100000, replace = T)
    null_dnes <- na.omit(dnes[ridx1,4]-dnes[ridx2,3])

    # Replace the NES of the original data with the deltaNES and calculate the
    # p-value of the deltaNES using the Monte Carlo null distribution
    dat_dnes <- left_join(dat, dnes %>% select(ID, dNES), by='ID') %>% 
      rename_with(., ~'NES', .cols='dNES')
    pval <- sapply(dat_dnes$NES, function(i){
      num <- ifelse(sign(i) == -1,
                    sum(null_dnes < i),
                    sum(null_dnes > i))
      denom <- length(null_dnes)
      (num/denom)*2
    })
    data.frame(dat_dnes$NES, pval)
    dat_dnes <- dat_dnes %>% 
      mutate(pvalue=pval,
             p.adjust=p.adjust(pval, method='BH'),
             qvalues=p.adjust)
    return(dat_dnes)
  })
})

# convert the GSEA data into a single cytoscape-format gsea df
for(celltype in names(dNES_ct_dat)){
  print(paste(celltype, "..."))
  dNES_dat <- dNES_ct_dat[[celltype]]
  dNES_df <- do.call(rbind, dNES_dat) %>% 
    as.data.frame %>% 
    mutate(p.adjust=p.adjust(pvalue, method='bonferroni')) %>%
    arrange(p.adjust, desc(abs(NES)))
  gsea_dat_i <-  gsea2CytoscapeFormat(dNES_df, gs_map)
  
  dir.create(file.path(outdir, "gse", "cytoscape", "ko_wt"), showWarnings = F)
  write.table(gsea_dat_i$up, 
              file=file.path(outdir, "gse", "cytoscape", "ko_wt",
                             paste0("gsea.", gsub(" ", "_", celltype), ".up.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(gsea_dat_i$down, 
              file=file.path(outdir, "gse", "cytoscape", "ko_wt",
                             paste0("gsea.", gsub(" ", "_", celltype), ".down.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  
}


#########################################################################
#### 6. SCENIC Regulon analysis between 72hr-vs-PBS for Tumor and LN ####
## Take the seurat object and split them into a tumor and TDLN subset,
# then reprocess them all indepedently, recluster, revisualize, and then
# integrate them together using MNN
outdirregulon <- file.path(outdir, "regulons")
dir.create(outdirregulon, showWarnings = F)
setwd(outdirregulon)
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))
seu$manual_anno <- seu$seurat_clusters %>%
  recode_map()
DefaultAssay(seu) <- 'RNA'

#--- a) Pre-pySCENIC ----
# Pre-filtering of the gene expression matrix
exprMat <- GetAssayData(seu, slot="data")
loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
exprMat_filt <- as.matrix(exprMat[loci1, ])
cellInfo <- seu@meta.data

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
loom <- build_loom(file.path(outdirregulon, "seu.loom"), dgem=exprMat_filt)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

#--- b) Post-pySCENIC ----
loom <- open_loom(file.path(outdirregulon, 'seu_scenic.loom'))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC), 
        file = file.path(outdirregulon, 'seu_regulon.rds'))
#--- c) Differential regulon per cluster ----
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

# Combine the delta LFC (KO-WT) into one matrix [regulon x cluster]
clusts <- names(delta_fcs[[1]])
clust_grp_lfcpadj <- lapply(setNames(clusts,clusts), function(clust_category){
  # Iterate through each clustering (manual_anno, seurat_clusters)
  
  tumorln <- names(delta_fcs)
  lapply(setNames(tumorln,tumorln), function(delta_grpid){
    delta_grp <- delta_fcs[[delta_grpid]]
    # Iterate through each sampletype (LN Tumor)
    nullidx <- sapply(delta_grp[[clust_category]], is.null)
    
    mergedeltafcs <- function(dfc, var='l2fc'){
      lapply(dfc, function(i){
        i$deltal2fc %>%
          as.data.frame %>% 
          tibble::rownames_to_column(., "gene") %>%
          select(gene, !!rlang::sym(var))
      })  %>%
        purrr::reduce(., .f=full_join, by='gene') %>% 
        tibble::column_to_rownames(., "gene") %>%
        rename_with(., ~names(dfc))
    }
    fc_mat <- mergedeltafcs(delta_grp[[clust_category]][-which(nullidx)], var='l2fc')
    padj_mat <- mergedeltafcs(delta_grp[[clust_category]][-which(nullidx)], var='padj')
    sigidx <- apply(padj_mat, 1, min, na.rm=T) < 0.1
    
    fc_padj_mat <- full_join(melt(t(fc_mat)), 
                             melt(t(padj_mat)), 
                             by=c('Var1', 'Var2')) %>% 
      rename_with(., ~c('Cluster', 'Regulon', 'L2FC', 'padj')) %>%
      mutate(SampleType=delta_grpid)
    return(list("mat"=fc_padj_mat, "sig"=names(which(sigidx))))
  })
})


pdf("~/xfer/regulon_KOvsWT_postNorm.pdf")
for(cluster in c('seurat_clusters', 'manual_anno')){
  b <- c(-0.15, 0.15)
  # cluster <- 'seurat_clusters' # 'seurat_clusters', 'manual_anno'
  lfcpadj <- lapply(clust_grp_lfcpadj[[cluster]], function(i) i$mat) %>%
    do.call(rbind,.) %>%  as.data.frame %>% 
    filter(Regulon %in% unlist(sapply(clust_grp_lfcpadj[[cluster]], function(i) i$sig))) %>%
    mutate(Cluster=as.character(Cluster))
  gp <- ggplot(lfcpadj, aes(x=Cluster, y=Regulon, color=L2FC, size=padj)) +
    facet_grid(.~SampleType, scales = 'free', space='free') +
    geom_point() + 
    theme_bw() +
    scale_color_gradientn(limits = c(min(b),max(b)),
                          colours=colorRampPalette(c('darkblue', "blue", "grey", "red", 'darkred'))(7),
                          breaks=seq(b[1], b[2], length.out=5), 
                          labels=seq(b[1], b[2], length.out=5)) +
    scale_size(trans='log10', range=c(10,0.5), limits = c(NA, 1)) +
    ylab("") + xlab("") + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
  plot(gp)
}
dev.off()

for(cluster in c('seurat_clusters', 'manual_anno')){
  lfcpadj <- lapply(clust_grp_lfcpadj[[cluster]], function(i) i$mat) %>%
    do.call(rbind,.) %>%  as.data.frame 
  outf <- file.path("~/xfer", paste0(cluster, ".lfc_padj_mat.csv"))
  write.table(lfcpadj, file=outf,
              sep=",", quote = F, col.names = T, row.names = F)
  
  cat(paste0("\nxfer ", basename(outf), "\n"))
}

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



##### XXX #####


{
  lfc_mat <- lapply(ct_markers, function(i){
    i %>% filter((p_val_adj < 0.01) & (abs(avg_log2FC) > 1)) %>%
      select(symbol, avg_log2FC)
  }) %>% 
    purrr::reduce(., full_join, by='symbol') %>%
    tibble::column_to_rownames(., "symbol") %>%
    rename_with(., ~celltypes)
  
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"))
  hm <- ComplexHeatmap::Heatmap(as.matrix(lfc_mat), cluster_rows=F,
                                cluster_columns=F, row_names_gp = grid::gpar(fontsize = 6),
                                col=col_fun, 
                                column_title=paste0(unique(seu_i$tissue_cond),
                                                    " 72hr vs PBS"))
  return(hm)
}
pdf("~/xfer/i.pdf", width = 5, height = 9); 
lfc_heatmaps
dev.off()


  
cds_groups <- lapply(seu.list, function(seu_i){
  Idents(seu_i) <- 'monocle3_partitions'
  seu_i <- subset(seu_i, idents='1')
  
  DefaultAssay(seu_i) <- 'integrated'
  seu_i <- RunPCA(seu_i, verbose = FALSE)
  seu_i <- RunUMAP(seu_i, reduction = "pca", dims = 1:30, n.neighbors=15L,
                 min.dist = 0.1)
  
  data <- as(as.matrix(seu_i@assays$integrated@scale.data), 'sparseMatrix')
  #data <- as(as.matrix(seu_i@assays$RNA@data), 'sparseMatrix')
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  cds <- new_cell_data_set(expression_data=data,
                           cell_metadata  = seu_i@meta.data,
                           gene_metadata = fData)
  
  cds <- preprocess_cds(cds, num_dim = 100, norm_method='size_only', pseudo_count=0)
  # cds <- preprocess_cds(cds, num_dim = 50, norm_method='log')
  # cds <- align_cds(cds, alignment_group = "orig.ident")
  # cds <- reduce_dimension(cds)
  DefaultAssay(seu_i) <- 'integrated'
  reducedDims(cds)$UMAP <- seu_i@reductions$umap@cell.embeddings
  reducedDims(cds)$PCA <- seu_i@reductions$pca@cell.embeddings[,1:30]
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  return(cds)
})


pdf("~/xfer/c1i.pdf", width = 14); 
lapply(cds_groups, function(cds){
  p1 <- plot_cells(cds, color_cells_by='manual_anno', group_cells_by='cluster', 
             show_trajectory_graph=T, cell_size=1, label_cell_groups=F)
  p2 <- plot_cells(cds, color_cells_by='orig.ident', group_cells_by='cluster', 
             show_trajectory_graph=T, cell_size=1, label_cell_groups=F)
  p1+p2
  # plot_cells(cds, color_cells_by='cluster',show_trajectory_graph=F, label_cell_groups=F)
  # plot_cells(cds, color_cells_by='partition',show_trajectory_graph=F, label_cell_groups=F)
  # plot_cells(cds, color_cells_by='seurat_clusters',show_trajectory_graph=F, label_cell_groups=F)
  # plot_cells(cds, color_cells_by='orig.ident',show_trajectory_graph=F, label_cell_groups=F)
  # plot_cells(cds, show_trajectory_graph=T)
})
dev.off()


