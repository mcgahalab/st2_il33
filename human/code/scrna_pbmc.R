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