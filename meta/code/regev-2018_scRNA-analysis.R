library(Seurat)
library(harmony)
library(dplyr)
library(Matrix)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(celldex)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/melanoma_meta'
outdir <- 'output/Regev_2018'
setwd(PDIR)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path("data/scRNA/Regev_2018", "GSE115978_cell.annotations.csv.gz"),
                 header=TRUE)
cnts <- read.csv(file.path("data/scRNA/Regev_2018", "GSE115978_counts.csv.gz"),
                 header=TRUE)

rownames(expr) <- expr[,1]
expr <- expr[,-1]
list.gene <- c('IL1RL1', 'IL33', 'PDCD1', 'CTLA4', 'TOX', 'KLRG1')

######################
#### 1. Functions ####
getMarkers <- function(database='panglaodb', 
                       dbpath='/cluster/projects/mcgahalab/ref/scrna/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv',
                       n=2, species='homo_sapiens', organ='Immune system'){
  db <- read.table(dbpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if(database=='panglaodb'){
    db_spl <- split(db, db$organ)
    if(organ=='all'){
      db_spl <- split(db, list(db$organ, db$cell.type))
      names(db_spl) <- gsub("\\.", "_", names(db_spl))
      db_spl <- db_spl[which(sapply(db_spl, nrow)>0)]
    } else {
      print(" > Organs: ")
      print(names(db_spl))
      db_spl <- split(db_spl[[organ]], db_spl[[organ]]$cell.type)
    }
    
    db_markers <- lapply(db_spl, function(db){
      # db <- db[which(db$canonical.marker ==1),]
      if(species=='homo_sapiens'){
        # with(db, ubiquitousness.index < 0.1 & sensitivity_human > 0.5 & specificity_human > 0.2)
        # (1-ubiquitousness.index) + 
        ord <- order(with(db, sensitivity_human + (1-specificity_human)), 
                     decreasing = TRUE)
        head(db[ord,]$official.gene.symbol,n)
      } else if(species=='mus_musculus'){
        ord <- order(with(db, (1-ubiquitousness.index) + sensitivity_mouse + specificity_mouse), 
                     decreasing = TRUE)
        firstup(tolower(head(db[ord,]$official.gene.symbol,n)))
      } else {
        stop("species must be either 'homo_sapiens' or 'mus_musculus'")
      }
      
    })
    names(db_markers) <- gsub(" ", "_", names(db_markers))
  }
  unlist(db_markers)
}

genDEGplots <- function(seu, stat_test='wilcox', celltypes='NA'){
  Idents(seu) <- seu$treatment.group
  all_degs <- FindMarkers(seu, ident.1 = "post.treatment", 
                          ident.2 = "treatment.naive", test.use = stat_test)
  deg_markers <- FindMarkers(seu, features=list.gene,
                             ident.1 = "post.treatment", ident.2 = "treatment.naive", 
                             test.use = stat_test, min.pct = 0.01, logfc.threshold = 0.01)
  all_degs$gene <- rownames(all_degs)
  deg_markers$gene <- rownames(deg_markers)
  ggdeg <- ggplot(all_degs, aes(x=avg_log2FC, y=p_val_adj)) +
    geom_point() +
    scale_y_continuous(log10_trans()) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_minimal()
  
  ggdeg2 <- ggdeg + 
    geom_point(data=deg_markers, aes(x=avg_log2FC, y=p_val_adj, colour=gene)) + 
    geom_text(data=deg_markers, aes(colour = factor(gene)), label=deg_markers$gene, 
              check_overlap = TRUE, hjust = 0)
  
  vln <- VlnPlot(seu, features = list.gene, ncol = 3) +
    labs(subtitle=celltypes)
  tbl <- tableGrob(round(deg_markers[,-6], 6), theme = mytheme, rows = rownames(deg_markers))
  plot_grid(vln, ggdeg2, tbl, nrow=3, ncol=1, rel_heights=c(2,1,1))
}

#################################
#### 2. Create Seurat Object ####
rownames(cnts) <- cnts[,1]
cnts <- cnts[,-1]
seu <- CreateSeuratObject(counts = cnts, min.cells = 3, min.genes = 200, assay = "RNA", project = "regev_2018")
seu@meta.data <- cbind(seu@meta.data, meta[match(rownames(seu@meta.data),meta$cells),])
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

## preprocess the data
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu_filt <- FindVariableFeatures(seu, mean.function = ExpMean, 
                                 dispersion.function = LogVMR, nfeatures = 2000)
seu_filt <- CellCycleScoring(seu_filt, s.features = cc.genes$s.genes, 
                             g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
seu_filt$CC.Difference <- seu_filt$S.Score - seu_filt$G2M.Score
seu_filt <- ScaleData(object = seu_filt, vars.to.regress = c("CC.Difference"))
seu_filt <- RunPCA(seu_filt, npcs = 30, features = VariableFeatures(object = seu_filt), verbose = FALSE)

#Confirm #PC's determined explain > 90% of variance
stdev <- seu_filt@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))

seu_filt <- FindNeighbors(seu_filt, reduction = "pca", dims = 1:PCNum)
seu_filt <- FindClusters(seu_filt, resolution = 1.2)
seu_filt <- RunUMAP(seu_filt, reduction = "pca", dims = 1:PCNum)

## Find Neighbours and Cluster with HARMONY
dir.create(file.path(outdir, "harmony"), showWarnings = FALSE)
pdf(file.path(outdir, "harmony", "harmony_convergence.pdf"))
seu_harmony <- seu_filt %>% 
  RunHarmony("samples", plot_convergence = TRUE)
dev.off()

seu_harmony <- FindNeighbors(object = seu_harmony, dims = 1:PCNum, reduction ="harmony")
seu_harmony <- FindClusters(object = seu_harmony, resolution = 1.2, reduction ="harmony")
# Tweaked the UMAP parameters here
seu_harmony <- RunUMAP(object = seu_harmony, dims = 1:PCNum, reduction = "harmony",
                       n.neighbors = 30L, n.components = 2L, n.epochs=400L, min.dist=0.2)
SeuratDisk::SaveH5Seurat(seu_harmony, filename = file.path(outdir, "seurat_obj", 
                                                           paste0("regev2018_harmony.h5seurat")), overwrite = TRUE)

#############################
#### 3. Marker Detection ####
seu_harmony <- SeuratDisk::LoadH5Seurat(file.path(outdir, "seurat_obj", "regev2018_harmony.h5seurat"))
datatype <- 'Gene_Expression'
if(!file.exists(file.path(outdir, "harmony", paste0(datatype,"_markers.rda")))){
  Idents(seu_harmony) <- seu_harmony$seurat_clusters
  seu_harmony_markers <- FindAllMarkers(object = seu_harmony, only.pos = TRUE, min.pct = 0.25, 
                                        thresh.use = 0.25)
  # seu.harmony.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
  col_idx <- which(colnames(seu_harmony_markers) == 'avg_log2FC')
  colnames(seu_harmony_markers)[col_idx] <- 'avg_logFC'
  write.table(seu_harmony_markers, file=file.path(outdir, "harmony", "seurat_markers.csv"),
              sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
  save(seu_harmony_markers, file=file.path(outdir, "harmony", paste0(datatype,"_markers.rda")))
} else {
  load(file.path(outdir, "harmony", paste0(datatype,"_markers.rda")))
}

###################################################
#### 4. a) Assigning cell type clusters - SCSA ####
# Cell type annotation was done using SCSA, run indepedently of this 
# script on the findAllMarkers results
Idents(seu_harmony) <- seu_harmony$seurat_clusters
cluster_keep <- length(table(seu_harmony_markers$cluster))
nL <- R.utils::countLines(file.path(outdir, "harmony", 'annotate.txt'))
df <- read.csv(file.path(outdir, "harmony", 'annotate.txt'), header=FALSE, skip=nL-cluster_keep,
               check.names = FALSE, stringsAsFactors = FALSE)

df$V1 <-  gsub("\\[", "", df$V1)
for (coli in seq_along(df)) {
  df[,coli] <- gsub("\\'", "", df[,coli])
}
df$sortedID <- sapply(strsplit(df$V3, split="\\|"), function(i) paste(sort(gsub("^ ", "", i)), collapse="/"))
cluster_ids <- setNames(df$sortedID, df$V1)

# manually resolve mixed annotations
cluster_ids['0'] <- strsplit(cluster_ids['0'], '/')[[1]][1] #T-cell (CD4)
cluster_ids['1'] <- strsplit(cluster_ids['1'], '/')[[1]][1]
cluster_ids['17'] <- strsplit(cluster_ids['17'], '/')[[1]][2]
cluster_ids['17'] <- strsplit(cluster_ids['17'], '/')[[1]][2]
cluster_ids['19'] <- strsplit(cluster_ids['19'], '/')[[1]][1]
cluster_ids['14'] <- strsplit(cluster_ids['14'], '/')[[1]][1]
cluster_ids['6'] <- strsplit(cluster_ids['6'], '/')[[1]][2]
cluster_ids['19'] <- strsplit(cluster_ids['19'], '/')[[1]][2]

seu_harmony$anno_clusters <- cluster_ids[as.character(seu_harmony$seurat_clusters)]

#######################################################
#### 4. b) Assigning cell type clusters - scRNAseq ####
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
      if(!file.exists(file.path(outdir, "annotation", rds_file))){
        print(rds_file)
        pred.seu <- SingleR(test=GetAssayData(seu_harmony), ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            clusters=if(clus) seu_harmony$seurat_clusters else NULL)
        saveRDS(pred.seu, file=file.path(outdir, "annotation", rds_file))
      } else {
        blueprint_anno <- readRDS(file.path(outdir, "annotation", rds_file))
        id <- gsub("^.*blueprint(.*).rds", "bp\\1", rds_file)
        if(clus){
          cluster_ids <- setNames(blueprint_anno$labels, 
                                  as.character(rownames(blueprint_anno)))
          seu_harmony@meta.data[,id] <- cluster_ids[as.character(seu_harmony$seurat_clusters)]
        } else {
          seu_harmony@meta.data[,id] <- blueprint_anno$labels
        }
      }
    }
  }
} 

################################
#### 5. Cell type fractions ####
meltAnno <- function(cpt){
  anno_breakdown <- as.matrix(sapply(cpt, function(i) {
    table(factor(i, levels=unique(sort(unlist(cpt)))))
  }))
  anno_fracdown  <- apply(anno_breakdown, 2, function(i) i/sum(i))
  anno_fracdown  <- melt(anno_fracdown)
  colnames(anno_fracdown) <- c('celltype', 'treatment', 'frac')
  
  ggplot(anno_fracdown, aes(fill=treatment, y=frac, x=celltype)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
dir.create(file.path(outdir, "cellfrac"))

st2_pos <- factor(as.numeric(GetAssayData(seu_harmony['IL1RL1',])) > 0)
levels(st2_pos) <- c('ST2neg', 'ST2pos')
Idents(seu_harmony) <- st2_pos
seu_st2pos <- subset(x = seu_harmony, idents = 'ST2pos')
seu_st2neg <- subset(x = seu_harmony, idents = 'ST2neg')

pdf(file.path(outdir, "cellfrac", "bp_fine_cluster.pdf"), height=15)
plot_grid(meltAnno(split(seu_harmony$bp.fine.cluster, f=seu_harmony$treatment.group)),
          meltAnno(split(seu_harmony$bp.fine.cell, f=seu_harmony$treatment.group)),
          meltAnno(split(seu_harmony$cell.types, f=seu_harmony$treatment.group)),
          nrow=3, ncol=1)

plot_grid(meltAnno(split(seu_st2pos$bp.fine.cluster, f=seu_st2pos$treatment.group)),
          meltAnno(split(seu_st2pos$bp.fine.cell, f=seu_st2pos$treatment.group)),
          meltAnno(split(seu_st2pos$cell.types, f=seu_st2pos$treatment.group)),
          nrow=3, ncol=1)

plot_grid(meltAnno(split(seu_st2neg$bp.fine.cluster, f=seu_st2neg$treatment.group)),
          meltAnno(split(seu_st2neg$bp.fine.cell, f=seu_st2neg$treatment.group)),
          meltAnno(split(seu_st2neg$cell.types, f=seu_st2neg$treatment.group)),
          nrow=3, ncol=1)
dev.off()


############################
#### 6. Dimension Plots ####
# Umap visualization, separated by clusters
dir.create(file.path(outdir, "dimplot"), showWarnings = FALSE)
pdf(file.path(outdir, "dimplot", "dimplot_clusters.global.pdf"), width=17)
p_cl <- DimPlot(object = seu_harmony, reduction = "umap", group.by='seurat_clusters', 
                label = TRUE, pt.size = 0.5)
p_lbl <- DimPlot(object = seu_harmony, reduction = "umap", group.by='cell.types', 
                 label = TRUE, pt.size = 0.5)
p_scsa <- DimPlot(object = seu_harmony, reduction = "umap", group.by='anno_clusters', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpml <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.main.cluster', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpfl <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.fine.cluster', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpmc <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.main.cell', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpfc <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.fine.cell', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_tx <- DimPlot(object = seu_harmony, reduction = "umap", group.by='treatment.group', 
                label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_cl + p_lbl
p_lbl + p_scsa
p_lbl + p_bpfl
p_bpfc + p_bpfl
p_bpmc + p_bpml
p_lbl + p_tx
p_bpml + p_tx
dev.off()

############################
#### 7. Feature Plots ####
# Umap visualization, separated by clusters
db_markers <- getMarkers(n=4)
dbl <- split(db_markers, gsub("[0-9]*$", "", names(db_markers)))
dbl$GOI <- list.gene
dbl$Monocytes <- c(dbl$Monocytes[-length(dbl$Monocytes)], 'CD14')
dbl$monocytes <- c('CD14', 'ITGAM', 'CCR2', 'FCGR3A')
dbl$extravascular_monocytes <- c('THBD', 'ITGAX', 'HLA-DRA', 'CCR7')


dir.create(file.path(outdir, "featplot"), showWarnings = FALSE)
fplots <- lapply(names(dbl), function(cell_type){
  print(cell_type)
  db_marker <- dbl[[cell_type]]
  tryCatch({
    FeaturePlot(seu_harmony, features = db_marker, pt.size = 0.5, ncol = 2,
                cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
                label=FALSE, raster=TRUE) +
      labs(subtitle=cell_type)
  }, error=function(e){NULL})
})
names(fplots) <- names(dbl)

pdf(file.path(outdir, "featplot", "featplot_clusters.pdf"), width=12)
lapply(fplots, print)
dev.off()
pdf(file.path(outdir, "featplot", "featplot_clusters.goi.pdf"), width=12, height = 11)
fplots[['GOI']]
dev.off()

############################
#### 8. Dotplots of GOI ####
# Dotplot visualization of expression and %expressed per celltype
Idents(seu_harmony) <- seu_harmony$bp.fine.cell
dp_bfc <- DotPlot(seu_harmony, features =  list.gene,
                  cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 12, y.text = 12)

Idents(seu_harmony) <- seu_harmony$bp.fine.cluster
dp_bfl <- DotPlot(seu_harmony, features =  list.gene,
                  cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 12, y.text = 12)

Idents(seu_harmony) <- seu_harmony$cell.types
dp_clu <- DotPlot(seu_harmony, features =  list.gene,
                  cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 12, y.text = 12)

dir.create(file.path(outdir, "dotplot"))
pdf(file.path(outdir, "dotplot", "cluster_dotplots.pdf"), height = 18, width = 12)
plot_grid(dp_bfc, dp_bfl, dp_clu, nrow=3, rel_heights=c(3,1,1))
dev.off()


#############################################
#### 9. DiffExp between Treatment and Naive ####
stat_test <- 'wilcox'
mytheme <- gridExtra::ttheme_default(
  core = list(padding = unit(c(2.5, 2.5), "mm")))
dir.create(file.path(outdir, "deg"))


# All cells treat-vs-naive
pdf(file.path(outdir, "deg", paste0(stat_test, "_all.deg.pdf")), height = 12)
genDEGplots(seu_harmony, stat_test=stat_test, celltypes = 'All')
dev.off()

clusters <- c('cell.types', 'bp.fine.cell', 'bp.fine.cluster')
for(cl in clusters){
  print(paste0(cl, "..."))
  
  # cell-type specific treat-vs-naive
  Idents(seu_harmony) <- seu_harmony@meta.data[,cl]
  all_celltypes <- "Monocytes"
  all_celltypes <- na.omit(unique( seu_harmony@meta.data[,cl]))
  vt_plts <- lapply(all_celltypes, function(celltypes){
    print(celltypes)
    seu_celltype <- subset(x = seu_harmony, idents = celltypes)
    tryCatch({
      genDEGplots(seu_celltype, stat_test=stat_test, celltypes = celltypes)
    }, error=function(e){NULL})
  })
  
  pdf(file.path(outdir, "deg", paste0(stat_test, "_", cl, ".deg.pdf")), height = 12)
  print(vt_plts)
  dev.off()
}

####################################
#### 10.a Subset Monocyte cells ####
Idents(seu_harmony) <- seu_harmony$bp.fine.cluster
seu_monocytes <- subset(x = seu_harmony, idents = 'Monocytes')
seu_monocytes <- FindVariableFeatures(seu_monocytes, mean.function = ExpMean, 
                                      dispersion.function = LogVMR, nfeatures = 2000)
seu_monocytes <- RunPCA(seu_monocytes, npcs = 30, features = VariableFeatures(object = seu_monocytes), verbose = FALSE)
stdev <- seu_monocytes@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))

## Find Neighbours and Cluster with HARMONY
dir.create(file.path(outdir, "harmony"), showWarnings = FALSE)
pdf(file.path(outdir, "harmony", "harmony_convergence.monocytes.pdf"))
seu_monocytes <- seu_monocytes %>% 
  RunHarmony("samples", plot_convergence = TRUE)
dev.off()

seu_monocytes <- FindNeighbors(object = seu_monocytes, dims = 1:PCNum, reduction ="harmony")
seu_monocytes <- FindClusters(object = seu_monocytes, resolution = 1.2, reduction ="harmony")
# Tweaked the UMAP parameters here
seu_monocytes <- RunUMAP(object = seu_monocytes, dims = 1:PCNum, reduction = "harmony",
                         n.neighbors = 20L, n.components = 2L, n.epochs=400L, min.dist=0.2)

bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/BlueprintEncodeData.rds") 
lbl <- 'label.fine'
for(clus in c(TRUE, FALSE)){
  print(as.character(clus))
  pred.seu <- SingleR(test=GetAssayData(seu_monocytes), ref=bed.se, 
                      assay.type.test=1, labels=bed.se[[lbl]],
                      clusters=if(clus) seu_monocytes$seurat_clusters else NULL)
  id <- paste0("bp.", gsub("label.", "", lbl), ".", if(clus) 'cluster' else 'cell')
  if(clus){
    cluster_ids <- setNames(pred.seu$labels, 
                            as.character(rownames(pred.seu)))
    seu_monocytes@meta.data[,id] <- cluster_ids[as.character(seu_monocytes$seurat_clusters)]
  } else {
    seu_monocytes@meta.data[,id] <- pred.seu$labels
  }
}

##############################
#### 10.a Do the plotties ####
# DEG & ST2 counting
trimmed_clusters <- seu_monocytes$seurat_clusters
levels(trimmed_clusters) <- c(rep("M1", 3), "M2", rep("M1", 2), 
                              "M3", rep("M1", 5))
Idents(seu_monocytes) <- trimmed_clusters
degs <- FindAllMarkers(seu_monocytes)

st2_pos <- factor(as.numeric(GetAssayData(seu_monocytes['IL1RL1',])) > 0)
table(st2_pos)

# Dimplot
dir.create(file.path(outdir, "dimplot"), showWarnings = FALSE)
pdf(file.path(outdir, "dimplot", "dimplot_clusters.monocytes.pdf"), width=17)
p_sc <- DimPlot(object = seu_monocytes, reduction = "umap", group.by='seurat_clusters', 
                label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpfl <- DimPlot(object = seu_monocytes, reduction = "umap", group.by='bp.fine.cluster', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpfc <- DimPlot(object = seu_monocytes, reduction = "umap", group.by='bp.fine.cell', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_sc + p_bpfc
p_bpfl + p_bpfc
dev.off()

# Featplot
pdf(file.path(outdir, "featplot", "featplot_clusters.monocytes.pdf"), width=12, height = 11)
fp <- FeaturePlot(seu_monocytes, features = list.gene, pt.size = 0.5, ncol = 2,
                  cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
                  label=FALSE, raster=TRUE)
fp
dev.off()
