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
library(DESeq2)
library(ggrepel)
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
groups <- list.files(file.path(PDIR, "data", "datasets"))

# Make TxDB object
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

# doublet_barcodes <- file.path(outdir, "doubletfinder", "doublet_barcodes.rds")
# doublet_quantile_cutoff <- 0.95

###################
#### Functions ####
getRoughTpm <- function(cnt, deg, genes_txdb, sample_idx){
  tpm_idx <- which(rownames(deg) %in% genes_txdb$hugo)
  
  tpm <- rep(NA, nrow(deg))
  gene_length <- width(genes_txdb)[ which(genes_txdb$hugo %in% rownames(deg))]
  norm_cnts <- cnt[rownames(deg), sample_idx][tpm_idx,] /gene_length
  tpm[tpm_idx] <- (t(t(norm_cnts) / colSums(norm_cnts)) * 1e6) %>%
    rowMeans() %>%
    round(.,2)
  return(tpm)
}

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

######################################################
#### 0.a Create Data-type specific Seurat object  ####
seus_rds <- file.path(seurat_dir, "sc_atac_rna.rds")

library(EnsDb.Hsapiens.v86)
genome <- EnsDb.Hsapiens.v86
dir.create(file.path(outdir, "seurat_obj"), 
           recursive = TRUE, showWarnings = FALSE)


## Merge peaks to get reduced-combined peak 
# https://satijalab.org/signac/articles/merging.html 
gr_pks <- lapply(groups, function(grp){
  samples <- list.files(file.path(PDIR, "data", "datasets", grp))
  pks <- lapply(samples, function(s){
    pk_path <- file.path(PDIR, "data", "datasets", grp, s)
    pk <- read.table(file = file.path(pk_path, "atac_peaks.bed"),
                     stringsAsFactors = FALSE, check.names=F,
                     sep="\t", comment.char="#", header=F,
                     col.names = c("chr", "start", "end"))
    grpk <- makeGRangesFromDataFrame(pk)
    grpk <- keepStandardChromosomes(grpk, pruning.mode = 'coarse')
    return(grpk)
  })
  names(pks) <- samples
  return(pks)
}) %>%
  setNames(., groups)
combined_pks <- reduce(unlist(as(unlist(gr_pks, recursive = T), "GRangesList")))
pk_width <- width(combined_pks)
combined_pks <- combined_pks[pk_width  < 10000 & pk_width > 100]


dir.create(file.path(outdir, "seurat_obj", "samples"), showWarnings = F)
seus <- lapply(setNames(groups,groups), function(grp){
  samples <- list.files(file.path(PDIR, "data", "datasets", grp))
  seus_s <- lapply(setNames(samples, samples), function(s){
    rds_f <- file.path(PDIR, "data", "seurat_obj", "samples", paste0(s, ".rds"))
    if(file.exists(rds_f)){
      print(paste0("Reading in ", grp, "-", s, "..."))
      seu <- readRDS(rds_f)
      print("Read in!")
    } else {
      print(paste0("Could not find ", grp, "-", s, "..."))
      pk_path <- file.path(PDIR, "data", "datasets", grp, s)
      seu <- list()
      
      ## > Read in scRNA data
      {
        mtx <- Read10X(data.dir = file.path(pk_path, "filtered_feature_bc_matrix"), strip.suffix=TRUE)
        seu[['GEX']] <- CreateSeuratObject(counts = mtx[['Gene Expression']], 
                                           project = s, assay='RNA')
        seu[['GEX']]$group <- grp
      }
      
      ## > Read in scATAC data
      {
        # Read in metadata: 
        meta_atac <- read.csv(file = file.path(pk_path, "per_barcode_metrics.csv"), 
                              header = TRUE, row.names = 1)
        col_ids <- colnames(meta_atac)[grep("atac", colnames(meta_atac))]
        
        # create fragment objects
        frags_obj <- CreateFragmentObject(path = file.path(pk_path, "atac_fragments.tsv.gz"),
                                          cells = rownames(meta_atac), 
                                          max.lines = NULL)
        # Quantify peaks in each dataset using the combined reference
        frag_cnts <- FeatureMatrix(fragments = frags_obj,
                                   features = combined_pks,
                                   cells = rownames(meta_atac))
        chrom_assay <- CreateChromatinAssay(counts = frag_cnts, #mtx[['Peaks']],
                                            sep = c(":", "-"), 
                                            genome = seqinfo(genome), 
                                            fragments = frags_obj,
                                            min.cells = 10, 
                                            min.features = 200)
        seu_chrom <- CreateSeuratObject(chrom_assay, assay = "ATAC", 
                                        meta.data=meta_atac[,col_ids])
        seu[['Peaks']] <- RenameCells(seu_chrom, 
                                      new.names=gsub("-1$", "", Cells(seu_chrom)))
      }

      saveRDS(seu, file=rds_f)
    }
    print("I'm outski")
    return(seu)
  })
  return(seus_s)
})
if(!file.exists(seus_rds)) saveRDS(seus, file=seus_rds)

# Group together all the samples/groups into one seurat object
if(!file.exists(seus_rds)){
  sample_grp_map <- sapply(seus, names) %>% unlist
  sample_grp_map <- setNames(names(sample_grp_map), sample_grp_map) %>% 
    gsub("[0-9]$", "", .)
  seus <- unlist(seus, recursive = F)
  names(seus) <- as.character(grp_sample_map)
  saveRDS(seus, file=seus_rds) 
} else{
  seus <- readRDS(file=seus_rds)
}

## Merge the scRNA data and save
datatype <- 'GEX'
seus_rna <- lapply(names(seus), function(id) { appendDatasetId(seus, id, datatype) })

seu <- merge(x=seus_rna[[1]], y = seus_rna[-1], 
             add.cell.ids = names(seus), 
             project = dataset)
SeuratDisk::SaveH5Seurat(
  seu, filename = file.path(
    outdir, "seurat_obj",paste0(gsub("\\s", "_", datatype), ".h5seurat")
  ), 
  overwrite=TRUE
)
seu_gex <- seu


## Merge the scATAC data and save
datatype <- 'Peaks'
seus_peaks <- lapply(names(seus), function(id) { appendDatasetId(seus, id, datatype) })
seu <- merge(x=seus_peaks[[1]], y = seus_peaks[-1], 
             add.cell.ids = names(seus), 
             project = dataset)
saveRDS(seu, file=file.path(outdir, "seurat_obj",paste0(gsub("\\s", "_", datatype), ".rds")))
seu_peaks <- seu; rm(seu);
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
  seu_i <- SCTransform(seu_i, assay = 'RNA', new.assay.name = 'SCT',
                       vars.to.regress = c('percent.mt'), conserve.memory = T) #nCount_RNA regressed in vst
  seu_i <- CellCycleScoring(seu_i, s.features = str_to_title(cc.genes$s.genes), 
                            g2m.features = str_to_title(cc.genes$g2m.genes),
                            assay='SCT', set.ident = TRUE)
  seu_i$CC.Difference <- seu_i$S.Score - seu_i$G2M.Score
  seu_i <- SCTransform(seu_i, assay = 'RNA', new.assay.name = 'SCT',
                       vars.to.regress = c('percent.mt', 'CC.Difference'),
                       conserve.memory = T)
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
saveRDS(seu, file=file.path(datadir, "seurat_obj", "2_seu_integrated.rds"))


############################################
#### 1.b Preprocessing the counts - RNA ####
# Preprocessing of the RNA counts
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "2_seu_integrated.rds"))

## Find Neighbours and Cluster 
DefaultAssay(seu) <- 'integrated'
seu <- RunPCA(seu, npcs=500, verbose = FALSE)

#Confirm #PC's determined explain > 95% of variance
var <- seu@reductions$pca@stdev^2
dvar <- abs(diff(var))
PCNum <- c(min(which(cumsum(dvar)/sum(dvar) >= 0.99)) + 1,
           min(which(cumsum(var)/sum(var) >= 0.8)))
PCNum <- 213 #80% of variance

# Visualizing and selecting PCs
visualize_pcs <- FALSE
if(visualize_pcs){
  pdf(file.path(outdir, "clusters", "pca_dimplot.pdf"))
  .getVar <- function(x) paste0(x, " (", round(max(cumsum(var)[1:x], na.rm=T)/sum(var) * 100,1), "%)")
  for(dim in c(5, 10, 30, 60, 100, PCNum)){
    seu <- RunUMAP(object = seu, dims = 1:dim, n.epochs=200L)
    DimPlot(seu, label = FALSE, reduction = "umap", group.by="orig.ident") + ggtitle(.getVar(dim)) 
  }
  dev.off()
}
DefaultAssay(seu) <- 'integrated'
seu <- RunUMAP(object = seu, dims = 1:PCNum[1])
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:PCNum)
seu <- FindClusters(seu, resolution = 1.5)
seu <- BuildClusterTree(seu, reduction='pca')

saveRDS(seu, file=file.path(datadir, "seurat_obj", "3_seu_sct.rds"))

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


###############################################################
#### 7.a Separating out Monocyte population and preprocess ####
if(!exists("seu"))  seu <- readRDS(file = file.path(seurat_dir, "scIntegrated_anno.rds"))
Idents(seu) <- 'group'
seu <- subset(seu, idents=c('healthy', 'lowISS', 'highISS'))
Idents(seu) <- 'bp.fine.cluster'
seu_mono <- subset(seu, idents='Monocytes')


DefaultAssay(seu_mono) <- "RNA"
seu_mono <- SCTransform(seu_mono, assay = 'RNA', new.assay.name = 'SCT',
                        vars.to.regress = c('percent.mt', 'CC.Difference', 'orig.ident'),
                        conserve.memory = T)

## Find Neighbours and Cluster with HARMONY
dir.create(file.path(outdir, "harmony"), showWarnings = FALSE)
pdf(file.path(outdir, "harmony", "harmony_convergence.pdf"))
seu_mono <- RunPCA(seu_mono, verbose = FALSE)
seu_mono <- RunHarmony(seu_mono, c('orig.ident'), assay.use='SCT', plot_convergence = TRUE)
dev.off()

#Confirm #PC's determined explain > 95% of variance
stdev <- seu_mono@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.90))

DefaultAssay(seu_mono) <- "SCT"
seu_mono <- FindNeighbors(object = seu_mono, dims = 1:if(PCNum>40) 40 else PCNum)
seu_mono <- FindClusters(object = seu_mono, resolution = 0.6)
# Tweaked the UMAP parameters here
seu_mono <- RunUMAP(object = seu_mono, dims = 1:if(PCNum>40) 40 else PCNum,
                     n.neighbors = 30L, n.components = 2L, n.epochs=200L, min.dist=0.3)

# # build a joint neighbor graph using both assays
# DefaultAssay(seu_mono) <- 'SCT'
# seu_mono <- FindMultiModalNeighbors(object = seu_mono,
#                                      reduction.list = list("pca", "lsi"), 
#                                      dims.list = list(1:40, 2:40),
#                                      modality.weight.name = "RNA.weight",
#                                      verbose = TRUE)
# 
# # build a joint UMAP visualization
# seu_mono <- RunUMAP(object = seu_mono,
#                      nn.name = "weighted.nn",
#                      assay = "SCT", 
#                      verbose = TRUE)

saveRDS(seu_mono, file=file.path(seurat_dir, "scIntegrated_mono.rds"))

##########################################
#### 7.b Monocyte dim and featureplot ####
if(!exists("seu_mono"))  seu_mono <- readRDS(file=file.path(seurat_dir, "scIntegrated_mono.rds"))
monocyte_markers <- setNames(c('CD86', 'CD33', 'FCGR1A', 'CCR2', 'HLA-DRA', 'IL1RL1', 
                               'CD14', 'FCGR3A', 'HLA-DRA'),
                             c('CD86', 'CD33', 'CD64', 'CCR2', 'HLA-DR', 'ST2',  
                               'CD14', 'CD16', 'HLA-DRA'))

DefaultAssay(seu_mono) <- 'SCT'
# clusters
dp_clus <- DimPlot(seu_mono, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='seurat_clusters', pt.size=0.5, shuffle=TRUE)
dp_anno <- DimPlot(seu_mono, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='bp.fine.cell', pt.size=0.5, shuffle=TRUE)
dp_sample <- DimPlot(seu_mono, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='orig.ident', pt.size=0.5, shuffle=TRUE)
dp_group <- DimPlot(seu_mono, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='group', pt.size=0.5, shuffle=TRUE)

pdf(file.path(outdir, "dimplot", "mono_dimplots.pdf"), width=15)
dp_clus + dp_anno
dp_sample + dp_group
dev.off()

# Feature plots
DefaultAssay(seu_mono) <- 'SCT'
featp <- FeaturePlot(seu_mono, features = monocyte_markers)
pdf(file.path(outdir, "dimplot", "mono_featplots.pdf"), width=13, height = 16)
featp
dev.off()

########################################
#### 7.c Labelling Monocyte subsets ####
# Identify populations rich in monocyte subtypes
monocytes <- c('CD14', 'FCGR3A', 'HLA-DRA')
Idents(seu_mono) <- 'seurat_clusters'
DefaultAssay(seu_mono) <- 'RNA'
cts <- GetAssayData(seu_mono)

cl_cnt <- sapply(unique(as.character(seu_mono@meta.data[,c('seurat_clusters')])), function(cl){
  cl_idx <- which(Idents(seu_mono) == as.character(cl))
  round(rowSums(cts[monocytes,cl_idx]) / length(cl_idx),2)
}) %>% 
  t() %>%
  as.matrix() %>% as.data.frame() %>%
  tibble::rownames_to_column(., var="cluster")

hc <- hclust(dist(cl_cnt[,-1]))
cl_cnt[hc$order,]

monocyte_subtype_map <- list(
  '23'='Unk', '18'='Unk', '23'='Unk', '22'='Unk', '4'='Unk',
  '8'='ITM',
  '6'='NC', '20'='NC', '17'='NC', '19'='NC', '11'='NC',
  '10'='CL', '13'='CL', '12'='CL', '15'='CL', '0'='CL',
  '1'='CL', '9'='CL', '3'='CL', '14'='CL', '2'='CL',
  '5'='CL', '16'='CL', '21'='CL', '7'='CL'
)
mono_id <- monocyte_subtype_map[as.character(seu_mono@meta.data$seurat_clusters)]
seu_mono$monocyte_subtype <- as.character(mono_id)


DefaultAssay(seu_mono) <- 'SCT'
# clusters
dp_clus <- DimPlot(seu_mono, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='seurat_clusters', pt.size=0.5, shuffle=TRUE)
dp_mono <- DimPlot(seu_mono, label = TRUE, repel = TRUE, reduction = "umap",
                    group.by='monocyte_subtype', pt.size=0.5, shuffle=TRUE)
featp <- FeaturePlot(seu_mono, features = monocytes)

pdf(file.path(outdir, "dimplot", "mono_subtypes_dimplots.pdf"), width=15)
dp_clus + dp_mono
featp
dev.off()

saveRDS(seu_mono, file=file.path(seurat_dir, "scIntegrated_mono.rds"))
##################################################
#### 7.d PCA of Monocyte subtypes and samples ####
if(!exists("seu_mono"))  seu_mono <- readRDS(file=file.path(seurat_dir, "scIntegrated_mono.rds"))

seu_mono$group2 <- seu_mono$group %>% 
  gsub("highISS", "pretrt", .) %>%
  gsub("lowISS", "pretrt", .)
ids <- apply(seu_mono@meta.data[,c('orig.ident', 'monocyte_subtype', 'group2')], 
             1, paste, collapse="_") 
uids <- unique(ids)

DefaultAssay(seu_mono) <- 'RNA'
cnts <- GetAssayData(seu_mono)
cnts_grps <- sapply(setNames(uids,uids), function(id){
  rowMeans(cnts[,which(ids %in% id)]) %>%
    ceiling() %>%
    as.integer()
})

coldata <- colnames(cnts_grps) %>%
  as.data.frame() %>%
  rename_with(~"id") %>%
  mutate("sample"=gsub("_.*", "", id),
         "subtype"=gsub("^.*_(.*)_.*$", "\\1", id),
         "group"=gsub("^.*_", "", id))

# obtain normalized counts
keep_idx <- grep("[^Unk]", coldata$subtype)
dds <- DESeqDataSetFromMatrix(countData = cnts_grps[,keep_idx],
                              colData = coldata[keep_idx,],
                              design= ~ subtype)
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
summary(rowSums(counts(dds)))
vst_cnts <- varianceStabilizingTransformation(dds, blind=FALSE)
tcnts <- as.data.frame(t(assay(vst_cnts)))
pca <- prcomp(tcnts, scale=F)
percent_var <- round(pca$sdev^2/sum(pca$sdev^2),2)
names(percent_var) <- paste0("PC", c(1:length(percent_var)))

max_pc <- 6
pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc))],
                             "condition"=as.character(vst_cnts$subtype),
                             "sample"=as.character(vst_cnts$sample),
                             "group"=as.character(vst_cnts$group)))
for(id in paste0("PC", c(1:max_pc))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}
pca_x$Name <- rownames(pca_x)
pca_y <- pca_x 

pcs <- paste0("PC", c(1:max_pc))
pcs_l <- lapply(c(2:length(pcs)), function(i){c(pcs[[i-1]], pcs[[i]])})
ggps <- lapply(pcs_l, function(pc){
  ggplot(data=pca_y, aes_string(x=pc[1], y=pc[2], color="condition", 
                                fill='group', shape="group")) +
    geom_point() + 
    xlab(paste0(pc[1], " (", percent_var[pc[1]], ")")) +
    ylab(paste0(pc[2], " (", percent_var[pc[2]], ")")) +
    theme_classic()
})

dir.create(file.path(outdir, "pca"), showWarnings = F)
pdf(file.path(outdir, "pca", "pca_meta.pdf"), width = 9, height = 10)
plot_grid(plotlist = ggps, ncol = 2)
dev.off()

################################################
#### 7.e DEG between each monocyte subtypes ####
if(!exists("seu_mono"))  seu_mono <- readRDS(file=file.path(seurat_dir, "scIntegrated_mono.rds"))

# Mapping of ensemlble to hugo
genome <- org.Hs.eg.db
species <- 'Homo sapiens'
txby <- keys(genome, 'ENSEMBL')
ens2symb <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")
genes_txdb <- genes(txdb)
genes_txdb$hugo <- ens2symb[genes_txdb$gene_id]

# Subsetting monocytes for groups and removing samples
rm_idents <- c('CCP033', 'CCP246')
seu_mono$group2 <- seu_mono$group %>% 
  gsub("highISS", "pretrt", .) %>%
  gsub("lowISS", "pretrt", .)
monocyte_subtypes <- c("CL", "NC", "ITM")
names(monocyte_subtypes) <- monocyte_subtypes

DefaultAssay(seu_mono) <- 'RNA'
Idents(seu_mono) <- 'monocyte_subtype'
monocyte_degs <- lapply(monocyte_subtypes, function(mono_sub){
  print(mono_sub)
  seu_sub <- subset(seu_mono, idents = mono_sub)
  rm_idx <- which(seu_sub$orig.ident %in% rm_idents)
  seu_sub <- seu_sub[,-rm_idx]
  
  # Normalize/preprocess data
  seu_sub <- NormalizeData(seu_sub, normalization.method = "LogNormalize")
  seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst", nfeatures = 2000)
  seu_sub <- ScaleData(seu_sub)
  
  # Get matrices and counts
  dat <- GetAssayData(seu_sub, slot='data')
  cnt <- GetAssayData(seu_sub, slot='count')
  grp_idx1 <- which(seu_sub$group2 %in%  'pretrt')
  grp_idx2 <- which(seu_sub$group2 %in%  'healthy')
  deg <- FindMarkers(object = seu_sub,
                     only.pos = F,
                     ident.1 = 'pretrt', 
                     ident.2 = 'healthy',
                     group.by = 'group2', 
                     min.pct = 0.5, 
                     logfc.threshold = log2(2),
                     method='DESeq2')
  deg <- deg %>%
    mutate("grp1"="pretrt",
           "grp2"="healthy",
           "mean.1" = rowMeans(expm1(dat[rownames(deg), grp_idx1])) %>% round(.,3),
           "mean.2" = rowMeans(expm1(dat[rownames(deg), grp_idx2])) %>% round(.,3),
           "mean.tpm.1"=getRoughTpm(cnt, deg, genes_txdb, grp_idx1),
           "mean.tpm.2"=getRoughTpm(cnt, deg, genes_txdb, grp_idx2),
           "pseudocount"=1,
           "avg_FC" = 2^avg_log2FC %>% round(.,3))
  deg_sig <- deg[which(deg$p_val_adj <= 0.05),]

  deg_heatmap <- DoHeatmap(seu_sub, features=rownames(deg_sig), 
                           group.by='group2', raster=FALSE) +
    ggtitle(mono_sub)
  
  # log2(mean(expm1(dat[gene, idx2]))+1) - log2(mean(expm1(dat[gene, idx1]))+1)
  return(list("deg"=deg, "heatmap"=deg_heatmap))
})

dir.create(file.path(outdir, "monocyte_subset"), showWarnings = F)
pdf(file.path(outdir, "monocyte_subset", "monocyte_heatmaps.pdf"), width = 9, height = 10)
lapply(monocyte_degs, function(i) i$heatmap)
dev.off()

dir.create(file.path(outdir, "monocyte_subset"), showWarnings = F)
degs <- plyr::rbind.fill(lapply(monocyte_degs, function(i) i$deg)) %>%
  mutate("monocyte_subtype"=rep(names(monocyte_degs),
                                sapply(monocyte_degs, function(i) nrow(i$deg))))
write.table(degs, file=file.path(outdir, "monocyte_subset", "monocyte_heatmaps.csv"), 
            col.names = T, row.names = F, quote=F, sep=",")

###############################################
#### 7.f TESTING DEG between all monocyte subtypes ####
if(!exists("seu_mono"))  seu_mono <- readRDS(file=file.path(seurat_dir, "scIntegrated_mono.rds"))

# Mapping of ensemlble to hugo
genome <- org.Hs.eg.db
species <- 'Homo sapiens'
txby <- keys(genome, 'ENSEMBL')
ens2symb <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")
genes_txdb <- genes(txdb)
genes_txdb$hugo <- ens2symb[genes_txdb$gene_id]

# Subsetting monocytes for groups and removing samples
rm_idents <- c('CCP033', 'CCP246')
seu_mono$group2 <- seu_mono$group %>% 
  gsub("highISS", "pretrt", .) %>%
  gsub("lowISS", "pretrt", .)
monocyte_subtypes <- c("CL", "NC", "ITM")

DefaultAssay(seu_mono) <- 'RNA'
rm_idx <- which(seu_mono$orig.ident %in% rm_idents)
rm_idx <- c(rm_idx, which(!seu_mono$monocyte_subtype %in% monocyte_subtypes))
seu_sub <- seu_mono[,-rm_idx]

# Do the preprocessies
seu_sub <- NormalizeData(seu_sub, normalization.method = "LogNormalize")
seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst", nfeatures = 2000)
seu_sub <- ScaleData(seu_sub)

# Extracting the counts and normalized data
dat <- GetAssayData(seu_sub, slot='data')
cnt <- GetAssayData(seu_sub, slot='count')
grp_idx1 <- which(seu_sub$group2 %in%  'pretrt')
grp_idx2 <- which(seu_sub$group2 %in%  'healthy')
deg <- FindAllMarkers(object = seu_sub,
                      only.pos = F,
                      min.pct = 0.5, 
                      logfc.threshold = log2(2),
                      method='DESeq2')
deg <- deg %>%
  mutate("grp1"="pretrt",
         "grp2"="healthy",
         "mean.1" = rowMeans(expm1(dat[rownames(deg), grp_idx1])) %>% round(.,3),
         "mean.2" = rowMeans(expm1(dat[rownames(deg), grp_idx2])) %>% round(.,3),
         "mean.tpm.1"=getRoughTpm(cnt, deg, genes_txdb, grp_idx1),
         "mean.tpm.2"=getRoughTpm(cnt, deg, genes_txdb, grp_idx2),
         "pseudocount"=1,
         "avg_FC" = 2^avg_log2FC %>% round(.,3))
deg_sig <- deg[which(deg$p_val_adj <= 0.05),]

deg_heatmap <- DoHeatmap(seu_sub, features=rownames(deg_sig), 
                         group.by='group2', raster=FALSE) +
  ggtitle(mono_sub)






print(mono_sub)
seu_sub <- subset(seu_mono, idents = mono_sub)
rm_idx <- which(seu_sub$orig.ident %in% rm_idents)
seu_sub <- seu_sub[,-rm_idx]
dat <- GetAssayData(seu_sub, slot='data')
cnt <- GetAssayData(seu_sub, slot='count')
grp_idx1 <- which(seu_sub$group2 %in%  'pretrt')
grp_idx2 <- which(seu_sub$group2 %in%  'healthy')