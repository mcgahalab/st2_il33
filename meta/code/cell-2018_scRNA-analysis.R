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
library(scuttle)
library(celldex)
library(SingleR)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/melanoma_meta'
proj <- 'Cell_2018'
outdir <- file.path('output', proj)
setwd(PDIR)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
list.gene <- c('IL1RL1', 'IL33', 'PDCD1', 'CTLA4', 'TOX', 'KLRG1')

##############################
#### 0.a) Format Metadata ####
# Clean up the metadata df formatting
extra_meta <- read.csv(file.path("data/scRNA/Cell_2018", "metadata.csv"), header=TRUE)
annotations <- read.csv(file.path("data/scRNA/Cell_2018", "annotations.csv"), header=TRUE)

meta <- read.table(file.path("data/scRNA/Cell_2018", "GSE120575_patient_ID_single_cells.txt.gz"),
                   header=TRUE, stringsAsFactors = FALSE, check.names = FALSE, skip = 19, sep="\t")
meta <- meta[,1:7]
colnames(meta) <- gsub("characteristics: ", "", colnames(meta))
colnames(meta)[5] <- 'patient_ID'
meta <- meta[c(1:16291),]
meta$state <- gsub("_.*", "", meta$patient_ID)
meta$patient_ID <- gsub("^.*?_", "", meta$patient_ID) %>%
  gsub("_.*", "", .)


# combnine the metadata with the emetadata from the stables
meta2 <- merge(meta, annotations, by.x='title', by.y='Cell_Name', all=TRUE)
meta3 <- merge(meta2, extra_meta[,c('Patient_ID', 'Gender', 'Age', 'RECIST', 'Overall_survival', 'Status')], 
               by.x='patient_ID', by.y='Patient_ID', all=TRUE)
na_idx <- which(is.na(meta3$patient_ID) | is.na(meta3$Cluster_number))
meta <- meta3[-na_idx,]

# Format the metadata to match TCGA style
dead_idx <- which(meta$Status==1)
alive_idx <- which(meta$Status==0)
meta$vital_status <- 'Alive'
meta$vital_status[dead_idx] <- 'Dead'
meta$days_to_death <- meta$Overall_survival
meta$days_to_death[alive_idx] <- -Inf
meta$days_to_last_follow_up <- meta$Overall_survival
meta$days_to_last_follow_up[dead_idx] <- -Inf
meta$bcr_patient_barcode <- meta$title

#################################
#### 0.b) Read in TPM matrix ####
expr <- read.table(file.path("data/scRNA/Cell_2018", "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"),
                   header=TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep="\t", fill=TRUE)
colnames(expr) <- c(colnames(expr)[-1], ' ')
expr <- expr[,-ncol(expr)]
meta2 <- data.frame("samples"=colnames(expr),
                    "anno"=as.character(expr[1,]))
expr <- expr[-1,]

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
  Idents(seu) <- seu$treatment
  all_degs <- FindMarkers(seu, ident.1 = "post", 
                          ident.2 = "pre", test.use = stat_test)
  deg_markers <- FindMarkers(seu, features=list.gene,
                             ident.1 = "post", ident.2 = "pre", 
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
dir.create(file.path(outdir, "seurat_obj"))
seu <- CreateSeuratObject(counts = expr, min.cells = 3, min.genes = 200, assay = "RNA", project = "Cell_2018")
meta_merge <- meta[match(rownames(seu@meta.data), meta$bcr_patient_barcode),]
seu@meta.data <- cbind(seu@meta.data, meta_merge)
SeuratDisk::SaveH5Seurat(seu, filename = file.path(outdir, "seurat_obj", paste0(proj, ".h5seurat")), 
                         overwrite = TRUE)

## preprocess the data
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, 
                             g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- FindVariableFeatures(seu, mean.function = ExpMean, 
                            dispersion.function = LogVMR, nfeatures = 2000)
seu <- ScaleData(object = seu, do.scale=F, do.center=F,
                 vars.to.regress = c("CC.Difference", "percent.mt"))
seu <- RunPCA(seu, npcs = 30, features = VariableFeatures(object = seu), verbose = FALSE)

#Confirm #PC's determined explain > 90% of variance
stdev <- seu@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))

## Find Neighbours and Cluster with HARMONY
dir.create(file.path(outdir, "harmony"), showWarnings = FALSE)
pdf(file.path(outdir, "harmony", "harmony_convergence.pdf"))
seu_harmony <- seu %>% 
  RunHarmony("patient_ID", plot_convergence = TRUE)
dev.off()

seu_harmony <- FindNeighbors(object = seu_harmony, dims = 1:PCNum, reduction ="harmony")
seu_harmony <- FindClusters(object = seu_harmony, resolution = 1.2, reduction ="harmony")
# Tweaked the UMAP parameters here
seu_harmony <- RunUMAP(object = seu_harmony, dims = 1:PCNum, reduction = "harmony",
                       n.neighbors = 30L, n.components = 2L, n.epochs=400L, min.dist=0.2)
SeuratDisk::SaveH5Seurat(seu_harmony, filename = file.path(outdir, "seurat_obj", 
                                                           paste0(proj, "_harmony.h5seurat")), overwrite = TRUE)

###################################################
#### 3. Assigning cell type clusters: scRNAseq ####
seu_harmony <- SeuratDisk::LoadH5Seurat(file.path(outdir, "seurat_obj", 
                                                  paste0(proj, "_harmony.h5seurat")))

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
seu_st2pos <- subset(x = seu_harmony, idents = 'ST2pos') # 130/16291
seu_st2neg <- subset(x = seu_harmony, idents = 'ST2neg') # 16161/16291

pdf(file.path(outdir, "cellfrac", "bp_fine_cluster.pdf"), height=15)
plot_grid(meltAnno(split(seu_harmony$bp.fine.cluster, f=seu_harmony$response)),
          meltAnno(split(seu_harmony$bp.fine.cell, f=seu_harmony$response)),
          meltAnno(split(seu_harmony$Cluster_number2, f=seu_harmony$response)),
          nrow=3, ncol=1)

plot_grid(meltAnno(split(seu_st2pos$bp.fine.cluster, f=seu_st2pos$response)),
          meltAnno(split(seu_st2pos$bp.fine.cell, f=seu_st2pos$response)),
          meltAnno(split(seu_st2pos$Cluster_number2, f=seu_st2pos$response)),
          nrow=3, ncol=1)

plot_grid(meltAnno(split(seu_st2neg$bp.fine.cluster, f=seu_st2neg$response)),
          meltAnno(split(seu_st2neg$bp.fine.cell, f=seu_st2neg$response)),
          meltAnno(split(seu_st2neg$Cluster_number2, f=seu_st2neg$response)),
          nrow=3, ncol=1)
dev.off()

############################
#### 6. Dimension Plots ####
# Umap visualization, separated by clusters
# seu_harmony_bkup <- seu_harmony
# seu_harmony@reductions$umap@cell.embeddings[,1] <- seu_harmony_bkup$UMAP1
# seu_harmony@reductions$umap@cell.embeddings[,2] <- seu_harmony_bkup$UMAP2

dir.create(file.path(outdir, "dimplot"), showWarnings = FALSE)
pdf(file.path(outdir, "dimplot", "dimplot_clusters.global.pdf"), width=17)
p_cl <- DimPlot(object = seu_harmony, reduction = "umap", group.by='seurat_clusters', 
              label = TRUE, pt.size = 0.5)
p_lbl <- DimPlot(object = seu_harmony, reduction = "umap", group.by='Cluster_number2', 
              label = TRUE, pt.size = 0.5)
p_scsa <- DimPlot(object = seu_harmony, reduction = "umap", group.by='Cluster_number', 
              label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpml <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.main.cluster', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpfl <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.fine.cluster', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpmc <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.main.cell', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_bpfc <- DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.fine.cell', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_tx <- DimPlot(object = seu_harmony, reduction = "umap", group.by='response', 
              label = TRUE, pt.size = 0.5, shuffle = TRUE)
p_cl + p_lbl
p_lbl + p_scsa
p_lbl + p_bpfl
p_bpfc + p_bpfl
p_bpmc + p_bpml
p_lbl + p_tx
p_bpml + p_tx
dev.off()

# seu_harmony <- seu_harmony_bkup

##########################
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
DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.fine.cell', 
                  label = TRUE, pt.size = 0.5, shuffle = TRUE)
fplots
dev.off()

pdf(file.path(outdir, "featplot", "featplot_clusters.goi.pdf"), width=12, height = 11)
DimPlot(object = seu_harmony, reduction = "umap", group.by='bp.fine.cell', 
        label = TRUE, pt.size = 0.5, shuffle = TRUE)
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

Idents(seu_harmony) <- seu_harmony$Cluster_number2
dp_clu <- DotPlot(seu_harmony, features =  list.gene,
                  cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 12, y.text = 12)

dir.create(file.path(outdir, "dotplot"))
pdf(file.path(outdir, "dotplot", "cluster_dotplots.pdf"), height = 15, width = 12)
plot_grid(dp_bfc, dp_bfl, dp_clu, nrow=3, rel_heights=c(2,1,2))
dev.off()


####################################
#### 10.a Subset Monocyte cells ####
Idents(seu_harmony) <- seu_harmony$bp.fine.cell
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
  RunHarmony("patient_ID", plot_convergence = TRUE)
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
summary(as.numeric(GetAssayData(seu_monocytes['IL1RL1',]))[which(as.logical(st2_pos))])

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