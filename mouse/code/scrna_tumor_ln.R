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


doublet_quantile_cutoff <- 0.95

####################
#### Annotation ####
## Human
anno_geneset <- list("DC"=list("up"=c("CLE4C", "CLEC10A", "CST3", "FCER1A", "HLA-DRA"),
                               "down"=c("CD14")),
                     "DC1"=list("up"=c("CLEC9A", "HLA-DRA", "XCR1"),
                                "down"=c("CD14")),
                     "Endothelial"=list("up"=c("CDH5", "CLEC14A", "PECAM1", "VWF"),
                          "down"=c()),
                     "Eosinophil"=list("up"=c("RNASE2", "CLC", "ADGRE1", "CFD", "CCL4", "PNPLA6",
                                              "MARCKSL1", "SLC29A1", "PTGDR2"),
                                       "down"=c()),
                     "Erythrocyte"=list("up"=c("HBA2", "HBB"),
                                        "down"=c()),
                     "Exhaustion"=list("up"=c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
                                       "down"=c()),
                     "ILC1"=list("up"=c("IL7R", "KLRB1"),
                                 "down"=c("CD4E", "NKG7")),
                     "ILC2"=list("up"=c("GATA3", "KIT", "MS4A2", "TPSAB1", "TPSB2"),
                                 "down"=c("CD4E", "NKG7")),
                     "ILC3"=list("up"=c("KLRB1", "NCR1", "NCR2", "RORC"),
                                 "down"=c("CD4E", "NKG7")),
                     "Mast cell"=list("up"=c("GATA2", "KIT", "MS4A2", "TPSAB1", "TPSB2"),
                                      "down"=c()),
                     "Monocyte"=list("up"=c("CD14", "CD300E", "DNM2", "FCN1", "S100A8", "VCAN"),
                                     "down"=c("C1QC")),
                     "Neutrophil"=list("up"=c("CPD", "CSF3R", "FCGR3B", "S100A8", "S100A9"),
                                       "down"=c()),
                     "NK"=list("up"=c("GNLY", "KLRD1", "NCAM1", "NKG7"),
                               "down"=c("CD3E")),
                     "pDC"=list("up"=c("GZMB", "HLA-DRA", "IL3RA", "TSPAN13"),
                                "down"=c("CD14")),
                     "Plasma cell"=list("up"=c("IGHA1", "IGHG2", "IGHG1", "JCHAIN", "TNFRSF17"),
                                        "down"=c()),
                     "Platelet"=list("up"=c("PPBP"),
                                     "down"=c()),
                     "TAM"=list("up"=c("CD14", "CD68", "LYZ", "HLA-DRA"),
                                "down"=c()),
                     "TAM M2-like"=list("up"=c("CD16", "CD68", "CD163", "HLA-DRA", "LYZ"),
                                        "down"=c()))
## Mouse
anno_geneset <- list("DC"=list("up"=c("CLE4C", "CLEC10A", "CST3", "FCER1A", "H2-Eb1"),
                               "down"=c("CD14")),
                     "DC1"=list("up"=c("CLEC9A", "H2-Eb1", "XCR1"),
                                "down"=c("CD14")),
                     "Endothelial"=list("up"=c("CDH5", "CLEC14A", "PECAM1", "VWF"),
                                        "down"=c()),
                     "Eosinophil"=list("up"=c("RNASE2", "CLC", "ADGRE1", "CFD", "CCL4", "PNPLA6",
                                              "MARCKSL1", "SLC29A1", "PTGDR2"),
                                       "down"=c()),
                     "Erythrocyte"=list("up"=c("HBA2", "HBB"),
                                        "down"=c()),
                     "Exhaustion"=list("up"=c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
                                       "down"=c()),
                     "ILC1"=list("up"=c("IL7R", "KLRB1"),
                                 "down"=c("CD4E", "NKG7")),
                     "ILC2"=list("up"=c("GATA3", "KIT", "MS4A2", "TPSAB1", "TPSB2"),
                                 "down"=c("CD4E", "NKG7")),
                     "ILC3"=list("up"=c("KLRB1", "NCR1", "NCR2", "RORC"),
                                 "down"=c("CD4E", "NKG7")),
                     "Mast cell"=list("up"=c("GATA2", "KIT", "MS4A2", "TPSAB1", "TPSB2"),
                                      "down"=c()),
                     "Monocyte"=list("up"=c("CD14", "CD300E", "DNM2", "FCN1", "S100A8", "VCAN"),
                                     "down"=c("C1QC")),
                     "Neutrophil"=list("up"=c("CPD", "CSF3R", "FCGR3B", "S100A8", "S100A9"),
                                       "down"=c()),
                     "NK"=list("up"=c("GNLY", "KLRD1", "NCAM1", "NKG7"),
                               "down"=c("CD3E")),
                     "pDC"=list("up"=c("GZMB", "H2-Eb1", "IL3RA", "TSPAN13"),
                                "down"=c("CD14")),
                     "Plasma cell"=list("up"=c("IGHA1", "Ighg2b", "IGHG1", "JCHAIN", "TNFRSF17"),
                                        "down"=c()),
                     "Platelet"=list("up"=c("PPBP"),
                                     "down"=c()),
                     "TAM"=list("up"=c("CD14", "CD68", "LYZ1", "H2-Eb1"),
                                "down"=c()),
                     "TAM M2-like"=list("up"=c("FCGR3", "CD68", "CD163", "H2-Eb1", "LYZ1"),
                                        "down"=c()))

###################
#### Functions ####
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")

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
seu_doublet_list <- lapply(setNames(samples,samples), function(sample_i){
  seu_i <- subset(seu, idents=sample_i)
  doublet_rds_f <- file.path(doublet_dir, paste0(sample_i, ".rds"))
  
  if(file.exists(doublet_rds_f)){
    print(paste0("Loading existing analysis: ", sample_i))
    doublet <- readRDS(doublet_rds_f)
    return(doublet)
  } 
  
  # pre-process
  seu_i <- SCTransform(seu_i)
  seu_i <- RunPCA(seu_i, npcs=200)
  
  stdev <- seu_i@reductions$pca@stdev
  var <- stdev^2
  PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))
  
  seu_i <- FindNeighbors(object = seu_i, dims = 1:PCNum)
  seu_i <- FindClusters(object = seu_i, resolution = 0.6)
  seu_i <- RunUMAP(seu_i, dims = 1:PCNum)
  
  
  # cell-type identification
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/BlueprintEncodeData.rds") 
  lbl <- 'label.fine' #label.main
  mouse.bed.se <- bed.se
  rownames(mouse.bed.se) <- str_to_title(rownames(bed.se))
  blueprint_anno <- SingleR(test=GetAssayData(seu_i), ref=mouse.bed.se, 
                            assay.type.test=1, labels=mouse.bed.se[[lbl]],
                            clusters=seu_i$seurat_clusters)
  cluster_ids <- setNames(blueprint_anno$labels, 
                          as.character(rownames(blueprint_anno)))
  seu_i@meta.data[,'blueprint'] <- cluster_ids[as.character(seu_i$seurat_clusters)]
  
  # pK identification
  sweep_res_list <- paramSweep_v3(seu_i, PCs = 1:PCNum, sct = TRUE)
  sweep_seu <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_seu)
  pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # pk <- bcmvn$pK[which.max(bcmvn[as.numeric(as.character(bcmvn$pK))<0.2,]$BCmetric)]
  # pk <- as.numeric(as.character(pk))
  # pk <- 0.11
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(seu_i$blueprint)    
  nExp_poi <- round(0.016*nrow(seu_i@meta.data))  # assuming 3200 cells loaded https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seu_i <- doubletFinder_v3(seu_i, PCs = 1:PCNum, pN = 0.25, pK = pk,  # first bump under 0.1 from bcmvn_seu
                            nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  # seu_i@meta.data <- seu_i@meta.data[,-grep("28_28", colnames(seu_i@meta.data))]
  seul <- list("seu"=seu_i, "sweep"=sweep_res_list, 
               "pk"=pk, "nExp_poi"=nExp_poi, "bcmvn"=bcmvn)
  saveRDS(seul, file=doublet_rds_f)
  return(seu_i)
})

if(visualize){
  ggps <- lapply(seu_doublet_list, function(seu_i){
    seu_i <- seu_i$seu
    dp_dub <- DimPlot(seu_i, label = TRUE, repel = TRUE, reduction = "umap",
                      group.by=grep("DF", colnames(seu_i@meta.data), value = T), pt.size=0.5, shuffle=TRUE)
    dp_bp <- DimPlot(seu_i, label = TRUE, repel = TRUE, reduction = "umap",
                     group.by="blueprint", pt.size=0.5, shuffle=TRUE)
    plot_grid(plotlist=list(dp_dub+
                              theme(legend.position='none') + 
                              ggtitle(unique(seu_i$orig.ident)), 
                            dp_bp+theme(legend.position='none')), ncol=2)
  })
  
  pdf(file.path(outdir, "doublets", "doublet_umap.pdf"), height=18, width=7)
  plot_grid(plotlist=ggps, nrow=length(ggps))
  dev.off()
}

doublets <- sapply(seu_doublet_list, function(seu_i){
  seu_i <- seu_i$seu
  dub_idx <- seu_i@meta.data[,grep("DF", colnames(seu_i@meta.data), value = T)] == 'Doublet'
  colnames(seu_i)[which(dub_idx)]
}) %>% 
  unlist() %>% as.character()

seu$doublets <- 'Singlet'
seu@meta.data[which(colnames(seu) %in% doublets),'doublets'] <- 'Doublet'
seu_dubsumm <- table(seu@meta.data[,c('orig.ident', 'doublets')])

Idents(seu) <- 'doublets'
seu <- subset(seu, idents='Singlet')
saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_filt_dubrm.rds"))

############################################
#### 1.a Preprocessing the counts - RNA ####
# Preprocessing of the RNA counts
if(!exists("seu")) seu <- readRDS(file=file.path(datadir, "seurat_obj", "seu_filt_dubrm.rds"))

DefaultAssay(seu) <- "RNA"
# https://github.com/satijalab/seurat/issues/1739
seu <- SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT',
                   vars.to.regress = c('percent.mt', 'orig.ident'), conserve.memory = T) #nCount_RNA regressed in vst
seu <- CellCycleScoring(seu, s.features = str_to_title(cc.genes$s.genes), 
                        g2m.features = str_to_title(cc.genes$g2m.genes),
                        assay='SCT', set.ident = TRUE)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT',
                   vars.to.regress = c('percent.mt', 'CC.Difference', 'orig.ident'),
                   conserve.memory = T)

## Find Neighbours and Cluster with HARMONY
seu <- RunPCA(seu, npcs=200, verbose = FALSE)

#Confirm #PC's determined explain > 95% of variance
stdev <- seu@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))

seu <- FindNeighbors(object = seu, dims = 1:PCNum)
seu <- RunUMAP(object = seu, dims = 1:PCNum, n.epochs=200L)

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
if(!exists("seu")) seu <- readRDS(file.path(datadir, "seurat_obj", "seu_preprocess.rds"))

dir.create(file.path(outdir, "annotation"))
getBlueprint <- TRUE
if(getBlueprint){
  #ref <- BlueprintEncodeData()
  #saveRDS(ref, file="~/downloads/BlueprintEncodeData.rds")
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/immgen.rds") 

  for(lbl in c('label.fine', 'label.main')){
    for(clus in c(TRUE, FALSE)){
      rds_file <- paste0("celldex_immgen.", gsub("label.", "", lbl),
                         ".", if(clus) 'cluster' else 'cell', ".rds")
      id <- gsub("^.*immgen(.*).rds", "bp\\1", rds_file)
      if(!file.exists(file.path(outdir, "annotation", rds_file))){
        print(rds_file)
        singler_anno <- SingleR(test=GetAssayData(seu), ref=bed.se, 
                                  assay.type.test=1, labels=bed.se[[lbl]],
                                  clusters=if(clus) seu$seurat_clusters else NULL)
        saveRDS(singler_anno, file=file.path(outdir, "annotation", rds_file))
      } else {
        singler_anno <- readRDS(file.path(outdir, "annotation", rds_file))
      }
      
      if(clus){
        cluster_ids <- setNames(singler_anno$labels, 
                                as.character(rownames(singler_anno)))
        seu@meta.data[,id] <- cluster_ids[as.character(seu$seurat_clusters)]
      } else {
        seu@meta.data[,id] <- singler_anno$labels
      }
    }
  }
} 

id_map <- unique(seu$bp.fine.cluster) %>%
  gsub("^.*\\(", "", .) %>%
  gsub("\\)", "", .) %>%
  gsub("(^.*?\\..*)\\..*", "\\1", .) %>%
  setNames(., unique(seu$bp.fine.cluster))
seu$bp.fine.cluster2 <- id_map[seu$bp.fine.cluster]
seu$bp.fine.cell2 <- id_map[seu$bp.fine.cell]

if(visualize){
  dpcl <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="seurat_clusters2", pt.size=0.5, shuffle=TRUE)
  dp1 <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="bp.main.cell", pt.size=0.5, shuffle=TRUE)
  dp1b <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="bp.main.cluster", pt.size=0.5, shuffle=TRUE)
  dp1c <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="bp.fine.cluster2", pt.size=0.5, shuffle=TRUE)
  dp1d <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="bp.fine.cell2", pt.size=0.5, shuffle=TRUE)
  dp2 <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="orig.ident", pt.size=0.5, shuffle=TRUE)
  pdf(file.path(outdir, "annotation", "dimplot_immgen_anno.pdf"), width=15)
  plot_grid(dpcl, dp1, ncol=2)
  plot_grid(dpcl, dp1b, ncol=2)
  plot_grid(dpcl, dp1c, ncol=2)
  plot_grid(dpcl, dp1d, ncol=2)
  plot_grid(dpcl, dp2, ncol=2)
  dev.off()
}


saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_anno.rds"))

#########################################################
#### 1.c Finding all markers for each seurat-cluster ####
if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "seu_anno.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

markers_f <- file.path(outdir, "annotation", "seurat_markers.rds")
if(!file.exists(markers_f)){
  Idents(seu) <- 'seurat_clusters2'
  DefaultAssay(seu) <- 'RNA'
  markers <- FindAllMarkers(object = seu, only.pos = FALSE, 
                            min.pct = 0.01, thresh.use = 0.01, 
                            test.use='wilcox', slot='data')
  markers %>%
    rename(., avg_logFC=avg_log2FC) %>% 
    write.table(., sep=",", file=file.path(outdir, "markers", "markers.csv"),
                quote=F, row.names = F, col.names = T)
  saveRDS(markers, file=markers_f)
} else {
  markers <- readRDS(markers_f)
}

markers_cl <- split(markers, markers$cluster)
x <- lapply(anno_geneset, function(anno_gs){
  lapply(names(anno_gs), function(dir){
    anno_dir <- str_to_title(anno_gs[[dir]])
    lfc <- sapply(markers_cl, function(mcl){
      if(length(anno_dir)==0) anno_dir <- ''
      setNames(mcl[anno_dir,]$avg_log2FC,
               anno_dir)
    }) %>%
      round(., 2)
    if(!is.matrix(lfc)) lfc <-  matrix(lfc, nrow=1, dimnames=list(anno_dir, names(markers_cl)))
    lfc[is.na(lfc)] <- 0
    return(lfc[,order(colSums(lfc))])
  })
})
markers_cl <- split(markers, f=markers$cluster)

marker_id <- list("DC"=c("Cst3"),
                  "Mast cells"=c("Gata2", "Ms4a2", "Tpsab1", "Tpsb2"),
                  "Neutrophil"=c("Cpd", "Csf3r", "S100a8", "S100a9"),
                  "NK"=c("Klrd1", "Nkg7", "Cd3e"),
                  "Plasma cell"=c("Ighg2", "Ighg1", "Jchain"),
                  "TAM"=c("Cd68", "Fcgr3", "Lyz", "Cd14"))

marker_lbl <- c("16"="Mast cells",
                "17"="DC",
                "1"="Neutrophil",
                "15"="NK",
                "2"="Plasma cell",
                "13"="TAM",
                "10"="TAM")
seu$marker_labels <- marker_lbl[seu$seurat_clusters2]


if(visualize){
  dpcl <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by="seurat_clusters2", pt.size=0.5, shuffle=TRUE)
  dp1 <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="bp.fine.cell2", pt.size=0.5, shuffle=TRUE)
  dp2 <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                 group.by="marker_labels", pt.size=0.5, shuffle=TRUE)
  
  pdf(file.path(outdir, "annotation", "dimplot_manual_anno.pdf"), width=15)
  plot_grid(dpcl, dp1, ncol=2)
  plot_grid(dpcl, dp2, ncol=2)
  FeaturePlot(seu, features = unlist(marker_id))
  dev.off()
}

saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_anno.rds"))

###################################################
#### 2.a Split Tumor and LN; Visualize samples ####
dir.create(file.path(outdir, "dimplot"), showWarnings = F)

if(!exists("seu"))  seu <- readRDS(file = file.path(datadir, "seurat_obj", "seu_anno.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

#split into tumor-ln, re-run umap on this smaller subset and re-cluster
seu$LN_Tumor <- gsub("_.*", "", seu$orig.ident)
ln_tumor <- unique(seu$LN_Tumor)

# seu_ln <- subset(seu, idents='LN')
# seu_tumor <- subset(seu, idents='Tumor')

# Re-run PCA and UMAP on the subsetted SCT assays
seu_ln_tumor <- lapply(setNames(ln_tumor, ln_tumor), function(ln_or_tumor){
  Idents(seu) <- 'LN_Tumor'
  seu_i <- subset(seu, idents=ln_or_tumor)
  DefaultAssay(seu_i) <- 'SCT'
  
  #Confirm #PC's determined explain > 95% of variance
  seu_i <- RunPCA(seu_i, npcs=50, verbose = FALSE)
  stdev <- seu_i@reductions$pca@stdev
  var <- stdev^2
  PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))
  seu_i <- FindNeighbors(object = seu_i, dims = 1:PCNum)
  seu_i <- RunUMAP(object = seu_i, dims = 1:PCNum, n.epochs=200L)
  seu_i <- FindClusters(object = seu_i, resolution = 1.2)
  
  return(seu_i)
})

seu_markers_ln_tumor <- lapply(seu_ln_tumor, function(seu_i){
  markers <- FindAllMarkers(object = seu_i, only.pos = TRUE, 
                            min.pct = 0.25, thresh.use = 0.25, 
                            test.use='wilcox', slot='data')
  markers
})

seu_grp_markers <- lapply(seu_ln_tumor, function(seu_i){
  seu_i$timepoint <- sapply(strsplit(seu_i$orig.ident, split="_"), function(i) i[3])
  seu_i$treatment <- sapply(strsplit(seu_i$orig.ident, split="_"), function(i) i[2])
  seu_i$tissue_treatment <- with(seu_i@meta.data, paste0(LN_Tumor, "_", treatment))
  
  DefaultAssay(seu_i) <- 'RNA'
  Idents(seu_i) <- 'bp.fine.cell2'
  tissue_treatment <- unique(seu_i$tissue_treatment)
  ct_markers <- lapply(setNames(tissue_treatment,tissue_treatment), function(sub){
    
    ids <- table(seu_i@meta.data[,c('tissue_treatment', 'orig.ident')])[sub,]
    ids <- ids[which(ids!=0)]
    
    celltypes <- levels(Idents(seu_i))
    lapply(setNames(celltypes, celltypes), function(celltype){
      print(paste0(celltype, "..."))
      tryCatch({
        markers <- FindMarkers(seu_i, slot = "data",test.use = "wilcox", 
                               ident.1=names(ids)[1], ident.2=names(ids)[2], 
                               group.by='orig.ident', subset.ident=celltype) %>%
          filter(p_val_adj < 0.05)
        
        c2 <- iterateMsigdb(species='Mus musculus', 
                            lfc_v=setNames(markers$avg_log2FC, sym2entrez_ids[rownames(markers)]), 
                            fun=gseaFun,
                            msig_lvls=list('C2'=list('CP:REACTOME')))
        return(list("markers"=markers, "gsea"=as.data.frame(c2[[1]][[1]][,4:10])))
      }, error=function(e){NULL})
    })
  })
  return(ct_markers)
})

pdf(file.path(outdir, "clusters", "markers_deg.pdf"), width = 15, width = 10)
pdf(file.path("~/xfer", "markers_deg.pdf"), height = 20, width = 10)
lapply(names(seu_ln_tumor), function(id){
  markers_i <- seu_markers_ln_tumor[[id]]
  seu_i <- seu_ln_tumor[[id]]
  
  top.markers <- do.call(rbind, lapply(split(markers_i, markers_i$cluster), head, n=12L))
  DoHeatmap(seu_i, features = top.markers$gene, raster=T, 
            group.by='seurat_clusters') + ggtitle(id)

})
dev.off()

# Calculate proportion of cell-types per sample
clusters <- list('seurat_clusters'=c("orig.ident", 'seurat_clusters'), 
                 'bp.fine.cell2'=c("orig.ident", 'bp.fine.cell2'),
                 'clusters_cell'=c("seurat_clusters", 'bp.fine.cell2'))
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
  x %>%
    round(., 3) %>%
    melt() %>%
    rename_with(., ~ c('Sample', 'cellType', 'proportion')) %>%
    left_join(ref_tbl, 
              by=c('Sample'='Sample')) %>%
    mutate(Sample=factor(Sample, levels=ref_lvls)) %>%
    ggplot(., aes(x=Sample, y=cellType, fill=proportion)) +
    facet_grid(cols=vars(join), scales='free') +
    geom_tile() +
    scale_fill_gradientn(colors = colorPal(10)) +
    theme(axis.text.x = element_text(angle = 90))
}
winsorize <- function(x, probs=c(0.05, 0.95)) {
  q <- quantile(as.numeric(x), probs)
  x[x<q[1]] <- q[1]
  x[x>q[2]] <- q[2]
  return(x)
}

pdf(file.path(outdir, "dimplot", "cellType_proportions.pdf"), width = 15)
col <- c('#feebe2', '#dd3497', '#7a0177')
colorPal <- grDevices::colorRampPalette(col)
# celltype_prop %>% t() %>%
ref_tbl <- colnames(celltype_props) %>% 
  strsplit(., "_", .) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>%
  mutate(Sample=colnames(celltype_prop))
ref_tbl$join <- with(ref_tbl, paste0(V1, "_", V2))
ref_tbl$V1 <- factor(ref_tbl$V1, c('LN', 'Tumor'))
ref_tbl$V2 <- factor(ref_tbl$V2, c('WT', 'KO'))
ref_tbl$V3 <- factor(ref_tbl$V3, c('PBS', '72h'))
ref_lvls <- ref_tbl[with(ref_tbl, order(V1, V2, V3)),]$Sample
  
# apply(celltype_prop, 1, function(i) scales::rescale(i, to=c(0,1))) %>%
gg_cps <- lapply(celltype_props, function(celltype_prop){
  mad_gg <- apply(celltype_prop, 1, function(i) {
    mad_i <- (i-median(i))/mad(i)
    if(any(is.nan(mad_i) | is.infinite(mad_i))){
      mad_i[which(is.nan(mad_i) | is.infinite(mad_i))] <- 0
    }
    return(mad_i)
  }) %>% 
    winsorize(., probs=c(0, 0.99)) %>% 
    .plotIt() + ggtitle("MAD proportions")
 
 abs_gg <- celltype_prop  %>% 
    t() %>% .plotIt() + ggtitle("Absolute proportions") + theme()
 
 plot_grid(mad_gg, abs_gg, nrow=1)
})
gg_cps
dev.off()

pdf(file.path(outdir, "dimplot", "dimplot_subsample.pdf"), width = 15)
DimPlot(seu_ln_tumor$Tumor, label=TRUE, repel=TRUE, reduction='umap',
        group.by='bp.fine.cluster2', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$Tumor, label=TRUE, repel=TRUE, reduction='umap',
        group.by='seurat_clusters', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$LN, label=TRUE, repel=TRUE, reduction='umap',
        group.by='bp.fine.cluster2', pt.size=0.5, split.by='orig.ident')
DimPlot(seu_ln_tumor$LN, label=TRUE, repel=TRUE, reduction='umap',
        group.by='seurat_clusters', pt.size=0.5, split.by='orig.ident')
dev.off()





saveRDS(seu, file=file.path(datadir, "seurat_obj", "seu_anno.rds"))

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





FindMarkers(seu,
            )



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



