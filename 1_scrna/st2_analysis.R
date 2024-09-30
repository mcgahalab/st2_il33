# renv::load("/cluster/projects/mcgahalab/envs/renvs/seuratv5_v2")
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
# Enrichment Analysis
library(msigdbr)
library(SCPA)
library(AUCell)
library(enrichplot)
# Regulon Analysis
library(SCopeLoomR)
library(SCENIC)
# QC 
library(scater)
library(DoubletFinder)
# Trajcectory/Velocity
library(parsnip)
library(tradeSeq)
library(condiments)
library(reticulate)
library("clusterProfiler")
library(ComplexHeatmap)
library(plyr)
library(grid)
library(GO.db)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e9)

visualize <- FALSE
seed <- 1234
set.seed(seed)
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d'
setwd(PDIR)

gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf' # Mus_musculus.GRCm38.102.gtf.gz
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
groups <- list.files(datadir, pattern="^(CD45|ST2|LN|Tumor)")

## Sets up geneset object
msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets
species <- 'Mus musculus'

min.geneset.size=10
geneset_l <- lapply(names(msig_lvls),function(mlvl){
  lapply(msig_lvls[[mlvl]], function(sublvl){
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
      dplyr::select(gs_name, entrez_gene, gene_symbol, gs_exact_source) %>%
      as.data.frame()
    keep_gs <- names(which(table(msig_ds$gs_name) >= min.geneset.size))
    msig_ds <- msig_ds %>%
      dplyr::filter(gs_name %in% keep_gs)
    return(msig_ds)
  }) %>% 
    setNames(., unlist(msig_lvls[[mlvl]]))
})
names(geneset_l) <- names(msig_lvls)
geneset_l <- unlist(geneset_l, recursive=F)
geneset_df <-  do.call(rbind, geneset_l) %>%
  mutate(entrez_gene = gene_symbol)

###################
#### Functions ####
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/list2df.R")
source("~/git/mini_projects/mini_functions/singlecell/writeSeuToFlat.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
source("~/git/mini_projects/mini_functions/singlecell/publicationDimPlot.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/bug_fixes/dietseurat.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")
gm <- geneMap(version='v2', species=species)

# Python functions for reticulate
# reticulate::use_python("/cluster/home/quever/miniconda3/bin/python3.9")
# python_pacmap <- import("pacmap")
# python_pandas <- import("pandas")
# python_numpy <- import("numpy")

##-- Recurring functions 1 ----
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
  seu_i$monocle3_partitions <- cds@clusters$UMAP$partitions
  return(seu_i)
}

## function used to create all unique 2-sample comparisons
.createGrps <- function(spl){
  grps <- lapply(spl, combn, m=2) %>%
    do.call(cbind, .)
  labels <- apply(grps, 2, function(i){
    table(c(strsplit(i[1], split = "\\.")[[1]],
            strsplit(i[2], split = "\\.")[[1]])) %>%
      sort %>% tail(., 1) %>% names %>%
      gsub("[^a-zA-Z0-9_ ]", "", .) %>%
      gsub(" ", "_", .)
  })
  grps %>%
    magrittr::set_colnames(make.unique(labels, sep = "."))
}

## function that uses regex to to help isolate which comparisons to be made
.selRows <- function(x, sample_regex, simple_regex, relevel_regex){
  x %>%
    mutate(s1=gsub(sample_regex, "\\1\\3", V1),
           s2=gsub(sample_regex, "\\1\\3", V2)) %>%
    dplyr::filter(s1 == s2) %>%
    dplyr::select(-c('s1', 's2')) %>%
    mutate(baselvl = paste0(c(grep(relevel_regex, V1, value=T, invert=T),
                              grep(relevel_regex, V2, value=T, invert=T)))) %>%
    mutate(id1=gsub(simple_regex, "\\2", V1),
           id2=gsub(simple_regex, "\\2", V2))
}

.relabelid <- function(x){
  if(any(grepl("Tumor", x))){
    x %>% 
      gsub("72h", "3d", .) %>%
      gsub("^Tumor", "B1_ST2_Tumor", .) %>%
      gsub("^LN", "B1_ST2_LN", .) %>%
      gsub("^ST2", "B2_ST2", .) %>%
      gsub("^CD45", "B2_CD45", .) %>%
      gsub("PBS$", "Un", .)
  } else {
    x %>% 
      gsub("72h", "3d", .) %>%
      gsub("^LN", "B1_ST2_LN", .) %>%
      gsub("^ST2", "B2_ST2", .) %>%
      gsub("^CD45", "B2_CD45", .) %>%
      gsub("PBS$", "Un", .)
  }
}

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
    AUCell::AUCell_run(expr_mat, msig_l)
  }, error=function(e){NULL})
  return(auc)
}

addLoomCellAnnotation <- function(loom, cellAnnotation){
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation))){
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation)){
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

sapply.allbyall <- function(x, fun, ...){
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(x), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  cnt <- 0
  sapply(x, function(i){
    setTxtProgressBar(pb, cnt)
    cnt <<- cnt + 1
    sapply(x, function(j){
      fun(i,j,...)
    })
  })
  close(pb) 
}

t_col <- function(color, percent = 0.50, name = NULL) {
  sapply(seq_along(color), function(colidx){
    if(length(percent)==1){
      percenti = percent
    } else if(length(percent) != length(color)){
      warning("length(percent) != length(color); using first entry")
      percenti = percent[1]
    } else {
      percenti <- percent[colidx]
    }
    coli <- color[colidx]
    
    rgb.val <- col2rgb(coli) * (((percenti/3)*2) + 0.33)
    rgb.val <- col2rgb(coli)
    rgb.val <- ((rgb.val/255) + percenti) * 255
    rgb.val[rgb.val > 255] <- 255
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 # alpha = (100 - percenti) * 255 / 100,
                 names = name)
    invisible(t.col)
  })
}

###########################################
#### 1. Preprocessing and integration  ####
## Description: 
# seus_preproc.rds = list of seurat object (CD45 and ST2) with doublets and singlets identified but no cells removed;
#   layers are split based on samples and each sample is individually normalized
# seus_integrated.rds = list of seurat object (CD45 and ST2) with doublets removed and samples integrated
#   using MNN. The layers were joined together and the integrated datasets were normalized

#--- 1.a) Create Seurat Objects ----
dir.create(file.path(datadir, "seurat_obj"), 
           recursive = TRUE, showWarnings = FALSE)

seus <- lapply(groups, function(grp){
  print(paste0(grp, "..."))
  mtx <- Read10X(data.dir = file.path(datadir, grp), strip.suffix=TRUE)
  seu <- CreateSeuratObject(counts = mtx, project = grp)
  return(seu)
})
names(seus) <- .relabelid(groups)


#--- 1.b QC - RNA data ----
mt_cutoff <- 10

seus_l <- split(seus, grepl("ST2", names(seus))) %>%
  setNames(., c('CD45', 'ST2'))
rm(seus); gc()

### Preprocess and outliers
seus_preproc_l <- lapply(setNames(names(seus_l),names(seus_l)), function(id){
  #--- Merge ----
  seus <- seus_l[[id]]
  seus <- lapply(setNames(names(seus),names(seus)), function(seuid){ 
    seus[[seuid]]$orig.ident <- seuid; 
    seus[[seuid]]
  })
  seu <- merge(seus[[1]], y = seus[-1], add.cell.ids = names(seus), project = id)
  seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mt")
  seu <- PercentageFeatureSet(seu, pattern = "^Rp[sl]", col.name= "percent.rps")
  seu@meta.data <- cbind(seu@meta.data, 
                         strsplit(seu$orig.ident, split="_") %>% 
                           do.call(rbind, .) %>%
                           as.data.frame %>% 
                           magrittr::set_colnames(c('Batch', 'Sort', 'Tissue', 'Condition', 'Timepoint')))
  
  #--- Split sample layers and preprocess ----
  seu <- JoinLayers(seu)
  DefaultAssay(seu) <- 'RNA'
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
  
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures(.)  %>%
    ScaleData(.)  %>%
    RunPCA(.) %>%
    FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
    FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
    RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  #--- Scan for outliers and doublets ----
  sces <- as.SingleCellExperiment.Seurat5(seu)
  mt_miqc <- FALSE
  
  #a) flexmix (miQC) method of mixed-modelling - Selected Method
  metrics <- seu@meta.data
  model_mt <- flexmix::flexmix(percent.mt ~ (nFeature_RNA), 
                               data = metrics,  k = 2)
  # sapply(split(posterior(model_mt)[,2]>0.997, f=seu$orig.ident), table) %>%t
  mt_outlier_miqc <- flexmix::posterior(model_mt)[,2]>0.997
  
  #b) 3-MAD/Mindiff from Median as outliers
  mt_outlier_scuttle <- isOutlier(seu$percent.mt, nmads=3, 
                                  min.diff=(mt_cutoff-median(seu$percent.mt)), 
                                  batch=seu$orig.ident, share.mads=T, 
                                  share.medians=T, type="higher")
  sapply(split(mt_outlier_scuttle, f=seu$orig.ident), table)
  
  seu$mt_outlier <- if(mt_miqc) mt_outlier_miqc else mt_outlier_scuttle
  # sapply(split(seu$percent.mt, f=seu$mt_outlier), summary)
  return(seu)
})


### Doublet calling
seus_preproc_l <- lapply(setNames(names(seus_preproc_l),names(seus_preproc_l)), function(seuid){
  seu <- seus_preproc_l[[seuid]]
  sce <- as.SingleCellExperiment.Seurat5(seu)
  
  print("Generating...")
  sce_cxds <- lapply(sce, runCxds, n=25)
  sce_bcds <- lapply(sce, runBcds, n=25)

  cxds_z <- lapply(sce_cxds, function(x) x$cxds %>% tibble::rownames_to_column('cell')) %>% do.call(rbind,.)
  bcds_score <- lapply(names(sce_bcds), function(id) {
    sce_bcds[[id]]$bcds %>%
      mutate(cell=colnames(sce[[id]])) %>%
      magrittr::set_colnames(., c('score', 'cell'))
  }) %>% do.call(rbind, .)
  
  seu$cxds_z <- setNames(cxds_z$z_cxds, cxds_z$cell)[Cells(seu)]
  seu$bcds_score <- setNames(bcds_score$score, bcds_score$cell)[Cells(seu)]
  seu$cxds_doublet <- seu$cxds_z > 3
  seu$bcds_doublet <- seu$bcds_score > 0.95
  
  
  seus_post <- SplitObject(seu, split.by='orig.ident')
  sing_cell_idx <- lapply(seus_post, function(seu_i){
    id <- unique(as.character(seu_i$orig.ident))
    print(id)
    print(ncol(seu_i))
    if(is.null(doublet_clusters)){
      idx <- c()
      dubidx <- c()
    } else {
      idx <- (!seu_i$seurat_clusters %in% doublet_clusters[[id]])
    }
    
    if(length(idx) == 0) idx <- rep(TRUE, ncol(seu_i))
    prop_doublets <- split(with(seu_i@meta.data, cxds_doublet + bcds_doublet) > 0,
                           seu_i$seurat_clusters) %>% 
      sapply(., function(i) round(sum(i)/length(i),2))
    if(!is.null(doublet_clusters)){
      dubidx <- which(names(prop_doublets) %in% doublet_clusters[[id]])
      if(any(prop_doublets[dubidx] < doublet.cutoff)){
        warning(paste0(id, "-CL", paste(doublet_clusters[[id]], collapse=","), ": is below doublet cutoff proportion"))
        warning(paste(prop_doublets[dubidx], collapse=","))
      }
    }
    dubidx <- if(length(dubidx) == 0) 10000 else  dubidx
    
    if(any(na.omit(prop_doublets[-dubidx] > doublet.cutoff))){
      exceed_doublet_clus <- names(which(prop_doublets[-dubidx] > doublet.cutoff))
      warning(paste0(id, "-CL", paste(exceed_doublet_clus, collapse=","), ": are above doublet cutoff proportion"))
      warning(paste(prop_doublets[exceed_doublet_clus], collapse=","))
    }
    return(setNames(idx, Cells(seu_i)))
  })
  sing_cell_idx <- stack(sing_cell_idx) %>% 
    tibble::rownames_to_column('barcode') %>%
    pull(values, name = barcode)
  seu$singlet <- sing_cell_idx[Cells(seu)]
  return(seu)
})

#--- 1.c) Integrating samples ----
percent.mt.cutoff <- 10
ncount.min <- 500
ncount.max <- Inf
nfeature.min <- 0
nfeature.max <- 6000

seus <- lapply(seus_preproc_l, function(seu){
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
  seu <- PercentageFeatureSet(seu, pattern = "^RP[SL]", col.name= "percent.rps")
  seu <- subset(seu, subset = 
                  nCount_RNA > ncount.min &
                  nCount_RNA < ncount.max &
                  nFeature_RNA > nfeature.min & 
                  nFeature_RNA < nfeature.max & 
                  percent.mt < percent.mt.cutoff & 
                  singlet)
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures(.)  %>%
    ScaleData(.)  %>%
    RunPCA(.) %>%
    FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
    FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
    RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  DefaultAssay(seu) <- 'RNA'
  seu_integ <- IntegrateLayers(
    object = seu, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = T,
    group.by='orig.ident'
  )
  seu_integ <- seu_integ %>% 
    FindNeighbors(., reduction = "integrated.mnn", dims = 1:30) %>%
    FindClusters(., resolution = 0.9, cluster.name = "mnn_clusters") %>% 
    RunUMAP(., reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
  seu_integ <- JoinLayers(seu_integ)
  return(seu_integ)
})
saveRDS(seus, file=file.path(PDIR, "data", "seurat_obj", "seus_integrated.rds"))

#######################
#### 2. Annotation ####
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_integrated.rds"))

#--- a) SingleR: Monaco ----
getImmgen <- TRUE
overwrite <- TRUE
references <- list("immgen"="immgen.rds",
                   "mouse"="MouseRNAseq.rds")

# RNA preprocess and cell-cycle score
seus_anno <- lapply(setNames(names(seus),names(seus)), function(seuid){
  seu <- seus[[seuid]]
  DefaultAssay(seu) <- 'RNA'
  
  
  for(database_id in names(references)){
    bed.se <- readRDS(file.path("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets", 
                                references[[database_id]])) 
    
    for(lbl in c('label.fine', 'label.main')){
      for(clus in c(TRUE, FALSE)){
        rds_file <- paste0("celldex_", database_id, ".", 
                           gsub("label.", "", lbl), ".", 
                           if(clus) 'cluster' else 'cell', ".",
                           seuid, ".rds")
        id <- gsub(paste0("^.*", database_id, "(.*).rds"), 
                   paste0(database_id, "\\1"), 
                   rds_file)
        if(!file.exists(file.path(outdir, "annotation2", rds_file)) | overwrite){
          print(paste0(" Annotating: ", rds_file))
          singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
                                  ref=bed.se, 
                                  assay.type.test=1, labels=bed.se[[lbl]],
                                  de.method="wilcox", genes='sd',
                                  clusters=if(clus) seu$seurat_clusters else NULL)
          saveRDS(singler_anno, file=file.path(outdir, "annotation2", rds_file))
        } else {
          print(paste0(" Reading in annotation: ", rds_file))
          singler_anno <- readRDS(file.path(outdir, "annotation2", rds_file))
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
  DefaultAssay(seu) <- 'RNA'
  
  
  pbmc_multimodal <- '/cluster/projects/mcgahalab/ref/scrna/pbmc_dataset/pbmc_multimodal_2023.rds'
  reference <- readRDS(pbmc_multimodal)
  anchor <- FindTransferAnchors(
    reference = reference,
    query = seu,
    reference.reduction = "spca",
    normalization.method = "SCT",
    dims = 1:50
  )
  seu <- MapQuery(
    anchorset = anchor,
    query = seu,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2"
    ),
    reduction.model = "wnn.umap"
  )
  return(seu)
})
saveRDS(seus_anno, file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))

###################################################
#### 3. Split ST2 Tumor and LN and remove doublets ####
if(!exists('seus_st2')) seus_st2 <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))

#--- a) Remove doublets and generate dimplots ----
if(!exists('seus')) seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
if(!exists("seus_st2")){
  seust2 <- seus$ST2
  seus_st2 <- split(seust2, f=seust2$orig.ident)
  Idents(seus_st2) <- 'orig.ident'
  uids <- as.character(unique(Idents(seus_st2)))
}

seus <- lapply(c('LN', 'Tumor'), function(sampleid){
  seu <- subset(seus_st2, ident=grep(sampleid, uids, value=T))
  seu <- DietSeurat(seu, layers=grep("count", Layers(seu), value=T),
                    assays='RNA')
  Idents(seu) <- 'manual_anno2'
  ids <- grep("^doublet", as.character(unique(Idents(seu))), 
              ignore.case = T, value=T, invert = T) %>%
    grep("^$", ., value=T, invert = T)
  seu <- subset(seu, ident=ids)
  
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures(.)  %>%
    ScaleData(.)  %>%
    RunPCA(.) %>%
    FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
    FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
    RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  DefaultAssay(seu) <- 'RNA'
  Idents(seu) <- 'orig.ident'
  seu_integ <- IntegrateLayers(
    object = seu, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = T
  )
  seu_integ <- seu_integ %>% 
    FindNeighbors(., reduction = "integrated.mnn", dims = 1:30) %>%
    FindClusters(., resolution = 0.9, cluster.name = "mnn_clusters") %>% 
    RunUMAP(., reduction = "integrated.mnn", dims = 1:30, 
            reduction.name = "umap.mnn",
            n.neighbors = 30L,
            min.dist = 0.3)
  seu_integ <- JoinLayers(seu_integ)
  for(nneigh in seq(30, 90, by=30)){
    for(mind in seq(0.1, 0.7, by=0.2)){
      print(paste0("umap.mnn_nn", nneigh, "_md", mind))
      seu_integ <- seu_integ %>% 
        RunUMAP(., reduction = "integrated.mnn", dims = 1:30, 
              reduction.name = paste0("umap.mnn_nn", nneigh, "_md", mind),
              n.neighbors = nneigh,
              min.dist = mind)
    }
  }
  return(seu_integ)
})
names(seus) <- c('LN', 'Tumor')
seus$LN@meta.data <- seus$LN@meta.data[Cells(seus$LN),]
seus$Tumor@meta.data <- seus$Tumor@meta.data[Cells(seus$Tumor),]
##############################
#### 4. Cell type subsets ####
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
#--- b) ST2 Treg cells ----
#--- b.i) Preprocess the reference TReg data ----
treg_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/gse_dat/treg_human_PBMID'
treg_rda <- '10X_all_data_preQC.rdata'
metaf <- 'all_data_cl_metadata.csv'

treg_meta <- read.csv(file.path(treg_refdir, metaf)) %>%
  tibble::column_to_rownames('X')
load(file.path(treg_refdir, treg_rda))
treg_seu <- UpdateSeuratObject(all_data)
rm(all_data); gc()
new.features <- gm$ENSEMBL$SYMBOL[Features(treg_seu)]
rmidx <- duplicated(new.features)
new.features <- as.character(make.unique(new.features))
new.features[which(is.na(new.features))] <- 'NArm'
new.features[which(new.features=='')] <- 'NArm'
treg_seu <- Seurat.utils::RenameGenesSeurat(treg_seu, new.features)
treg_seu <- subset(treg_seu, features = setdiff(Features(treg_seu), c('NArm', new.features[rmidx])))
treg_seu = Seurat::AddMetaData(treg_seu, metadata = treg_meta)


# Remove celltypes not wanted in analysis
Idents(treg_seu) <- 'tissue_sp'
treg_seu <- subset(treg_seu, ident=c('bLN', 'colon', 'mLN', 'spleen'))
Idents(treg_seu) <- 'cl_annot'
idents <- c(paste0("Treg_", c('lymphoid', 'effector', 'NLT', 'Stat1', 'suppressive',
                              'NLTlike', 'LTlike')), 
            'MitoStress')
treg_seu <- subset(treg_seu, ident=idents)



DefaultAssay(treg_seu) <- 'RNA'
treg_seul <- split(treg_seu, f=treg_seu$orig.ident)
treg_seu <- NormalizeData(treg_seul, normalization.method = "LogNormalize", 
                          scale.factor = 10000) %>%
  FindVariableFeatures(., num.bin = 100, binning.method = "equal_frequency")  %>%
  ScaleData(.,  model.use = "negbinom", do.center = T, do.scale = F, use.umi = T)  %>%
  RunPCA(., pcs.compute = 50, pcs.print = 1:3, genes.print = 10, pc.genes = rownames(treg_seul@data)) %>%
  FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
  FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DefaultAssay(treg_seu) <- 'RNA'
Idents(treg_seu) <- 'orig.ident'
seu_integ <- IntegrateLayers(
  object = treg_seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = T
)
seu_integ <- seu_integ %>% 
  FindNeighbors(., reduction = "integrated.mnn", dims = 1:30) %>%
  FindClusters(., resolution = 0.9, cluster.name = "mnn_clusters")
seu_integ <- JoinLayers(seu_integ)
for(nneigh in seq(30, 60, by=30)){
  for(mind in seq(0.1, 0.5, by=0.2)){
    print(paste0("umap.mnn_nn", nneigh, "_md", mind))
    seu_integ <- seu_integ %>% 
      RunUMAP(., reduction = "integrated.mnn", dims = 1:30, 
              reduction.name = paste0("umap.mnn_nn", nneigh, "_md", mind),
              n.neighbors = nneigh,
              min.dist = mind)
  }
}
seu_integ <- Seurat::AddMetaData(seu_integ, metadata = treg_meta)
treg_seu <- seu_integ
treg_seu <- treg_seu %>% 
  RunTSNE(., reduction = "integrated.mnn", dims.use = 1:30, seed.use = 1)
rm(seu_integ); gc()
saveRDS(treg_seu, file=file.path(treg_refdir, "seu_anno.normalized.mnn.umap.rds"))

#--- b.ii) Use transferAnchors to map our tregs to the reference dataset ----
treg_seu <- FindVariableFeatures(treg_seu,  num.bin = 100, binning.method = "equal_frequency")
treg_seu <- ScaleData(treg_seu)
seus_treg <- lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  Idents(seu) <- 'manual_anno2'
  ident <- unique(grep("treg", Idents(seu), ignore.case = T, value=T))
  seu_subset <- subset(seu, ident=ident)
  
  DefaultAssay(treg_seu) <- 'RNA'
  anchor <- FindTransferAnchors(
    reference = treg_seu,
    query = seu_subset,
    reference.reduction = "integrated.mnn",
    normalization.method = "LogNormalize",
    dims = 1:30
  )
  seu_subset <- MapQuery(
    anchorset = anchor,
    query = seu_subset,
    reference = treg_seu,
    refdata = list(
      cl_annot = "cl_annot"
    )
  )
  
  seu_subset <- seu_subset %>%
    RunUMAP(., reduction = "integrated.mnn", dims = 1:30, 
            reduction.name = 'umap_treg',
            n.neighbors = 20,
            min.dist = 0.4)
  
  seu_tmp <- seu_subset
  seu_tmp@meta.data <- seu_tmp@meta.data[,c('orig.ident', 'predicted.cl_annot',
                                            'manual_anno2')]
  seu_tmp[['RNA2']] <- CreateAssayObject(counts=seu_tmp@assays$RNA$counts)
  DefaultAssay(seu_tmp) <- 'RNA2'
  return(seu_subset)
}) %>% 
  setNames(., names(seus))

map <- c(setNames(seus$LN@meta.data[,'treg_anno'], Cells(seus$LN)),
         setNames(seus$Tumor@meta.data[,'treg_anno'], Cells(seus$Tumor)))
seus_treg$LN$treg_anno <- map[Cells(seus_treg$LN)]
seus_treg$Tumor$treg_anno <- map[Cells(seus_treg$Tumor)]
saveRDS(seus_treg, file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))

#--- b.iii) Re-integrate and dimensional reduction of just the TRegs cluster ----
if(!exists("seus_treg")) seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
dim<-15
seus_treg <- c(seus_treg, list("LN_Tumor"=JoinLayers(merge(seus_treg[[1]], seus_treg[-1]))))
seus_treg2 <- lapply(seus_treg, function(seu){
  Idents(seu) <- 'treg_anno'
  seu <- subset(seu, ident=setdiff(unique(Idents(seu)), 'NA'))
  # seu <- subset(seu, ident=c('eTreg_NLT-like', 'Treg-STAT1'))
  
  seu <- DietSeurat(seu, layers=grep("count", Layers(seu), value=T),
                    assays='RNA')
  seu <- split(seu, f=seu$orig.ident)
  
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures(.)  %>%
    ScaleData(.)  %>%
    RunPCA(.) %>%
    FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
    FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
    RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  DefaultAssay(seu) <- 'RNA'
  Idents(seu) <- 'orig.ident'
  seu_integ <- IntegrateLayers(
    object = seu, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = T
  )
  seu_integ <- seu_integ %>% 
    FindNeighbors(., reduction = "integrated.mnn", dims = 1:dim) %>%
    FindClusters(., resolution = 0.5, cluster.name = "mnn_clusters") %>% 
    RunUMAP(., reduction = "integrated.mnn", dims = 1:dim, 
            reduction.name = "umap.mnn",
            n.neighbors = 30L,
            min.dist = 0.3)
  seu_combined <- JoinLayers(seu_integ)
 
  return(seu_combined)

})
seus_treg <- seus_treg2

# Fix to treg_anno based on local Tumor or LN only population
ln_anno <- file.path(PDIR, "ref", "Sara.anno_final_CD45ST2Treg_LN.v2.csv")
tumor_anno <- file.path(PDIR, "ref", "Sara.anno_final_CD45ST2_TUMOR_Treg.v2.csv")
anno <- lapply(c(ln_anno, tumor_anno), read.csv, col.names=c('Barcode', 'Anno')) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  mutate(Anno = gsub(" $", "", Anno) %>% 
           gsub(" ", "_", .) %>%
           gsub("eTregNLT", "eTreg_NLT", .))
anno_map <- with(anno, setNames(Anno, Barcode))
seus_treg$LN$treg_anno <- as.character(anno_map[Cells(seus_treg$LN)])
seus_treg$Tumor@meta.data <- seus_treg$Tumor@meta.data[Cells(seus_treg$Tumor),]
seus_treg$Tumor$treg_anno <- as.character(anno_map[Cells(seus_treg$Tumor)])
for(id in c('LN', 'Tumor')){
  seus_treg[[id]] <- FindNeighbors(seus_treg[[id]], reduction = "integrated.mnn", 
                                   dims = 1:30, return.neighbor=T)
  idx <- which(anno$Anno == '')
  missing_anno <- sapply(grep(id, anno$Barcode[idx], value=T), function(cell_i){
    if(!cell_i %in% Cells(seus_treg[[id]])) return('NA')
    cells <- Seurat::TopNeighbors(seus_treg[[id]]@neighbors$RNA.nn, 
                                  cell=cell_i,
                                  n=20)
    majority_celltype <- anno %>% 
      dplyr::filter(Barcode %in% cells, Anno !='') %>% 
      pull(Anno) %>% table %>% sort %>% tail(., 1)
    if(length(majority_celltype)==0) majority_celltype <- table(cell_i)
    ifelse(majority_celltype < 10, 'NA', names(majority_celltype))
  })
  names(missing_anno) <- grep(id, anno$Barcode[idx], value=T)
  seus_treg[[id]]@meta.data[names(missing_anno), 'treg_anno'] <- as.character(missing_anno)
}

idx <- which(seus_treg$Tumor$treg_anno == 'Treg-STAT1' | 
               seus_treg$Tumor$treg_anno == 'NA')
seus_treg$Tumor$treg_anno[idx] <- 'eTreg_NLT-like'
anno_map2 <- c(seus_treg$LN$treg_anno, seus_treg$Tumor$treg_anno)
seus_treg$LN_Tumor$treg_anno <- anno_map2[Cells(seus_treg$LN_Tumor)]
seus_treg$LN_Tumor$treg_anno2 <- paste0(gsub(".*(LN|Tumor)_.*", "\\1", seus_treg$LN_Tumor$orig.ident),
                                         "_", seus_treg$LN_Tumor$treg_anno)
saveRDS(seus_treg, file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))


#--- b.iv) Visualize the tregs ----
seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))

grp <- 'treg_anno'
grps <- sapply(seus_treg, function(i) unique(i@meta.data[,grp])) %>% 
  unlist %>% unique
colors_use <- scCustomize::scCustomize_Palette(num_groups = length(grps)+1, 
                                               ggplot_default_colors = FALSE, 
                                               color_seed = 231)[-2] %>%
  setNames(., grps)
pdps <- lapply(names(seus_treg), function(seuid){
  print(seuid)
  seu <- seus_treg[[seuid]]
  colors_use_i <- as.character(colors_use[unique(seu@meta.data[,grp])])
  seu@reductions[['UMAP']] <- seu@reductions[['umap.mnn']]
  colnames(seu@reductions[['UMAP']]@cell.embeddings) <- paste0("UMAP_", c(1:2))
  seu@reductions[['UMAP']]@key <- 'UMAP_'
  seu@meta.data[,seuid] <- seu@meta.data[,grp]
  publicationDimPlot(seu, grp=seuid, simplify_labels=TRUE, 
                     colors_use=colors_use_i, pt.size=1.25, aspect.ratio=1,
                     reduction='UMAP',
                     legend.text = element_text(size=10),
                     labs(title=seuid))
})

pdf(file.path(PDIR, "results", "umap", "umap_CD45ST2_tumor_ln.TRegs_mnn.pdf"), width=7, height=5)
pdps
dev.off()

#########################
#### 5. DEG analysis ####
outdir <- file.path(PDIR, "results", "degs", "st2_analysis")
dir.create(outdir, recursive = T, showWarnings = F)
method <- 'wilcox'
sigdig <- 4

#--- a) ST2: All cells ----
grouping <- 'st2_all_cells'
dir.create(file.path(outdir, grouping), recursive = T, showWarnings = F)
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
overwrite <- F

## FindAllMarkers on everything and run DotPlot
allmarkers <- lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  Idents(seu) <- 'manual_anno2'
  seuid_markers <- FindAllMarkers(object = seu) %>%
    scCustomize::Add_Pct_Diff() 
  return(seuid_markers)
})
names(allmarkers) <-names(seus)

for(seuid in names(seus)){
  amarkers <- allmarkers[[seuid]] %>%
    dplyr::filter(pct_diff > 0.6)
  top_markers <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 20, 
                                     named_vector = FALSE, make_unique = TRUE)
  top_markers_tbl <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 50, 
                                     named_vector = FALSE, make_unique = TRUE)
  seu <- seus[[seuid]]
  Idents(seu) <- 'manual_anno2'
  genes <- read.csv(file.path(PDIR, "ref", "dotplots", paste0("ST2_", seuid, ".csv")),
           header = T, sep=",")
  
  write.table(amarkers %>% dplyr::filter(gene %in% top_markers_tbl),
              file=file.path("~/xfer", paste0("ST2_", seuid, ".csv")), 
              sep=",", col.names = T, row.names = F, quote = F)
  
  pdf(file.path("~/xfer", paste0("ST2_", seuid, ".pdf")), height = 20)
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
  
  cat(paste0("xfer ST2_", seuid, ".csv\n"))
  cat(paste0("xfer ST2_", seuid, ".pdf\n"))
  
}


# Output a rds containing all DEGs
# and a csv file containing the tables for the comparisons
method <- 'wilcox'
sigdig <- 3
gg_degumaps <- lapply(names(seus), function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$manual_anno2)
  uids <- unique(seu$group)
  x <- split(uids, f=gsub("^.*\\.", "", uids))
  xt <- .createGrps(x) %>% 
    t %>% as.data.frame 
  
  # list of 3d/7d vs Untreated samples per celltype
  xt1 <- xt %>%
    .selRows("(^.*)?(3d|7d|Un)\\..*", "(^.*)?((KO|WT)_(3d|7d|Un))\\..*", "Un")
  # list of KO vs WT
  xt2 <- xt %>%
    .selRows("(^.*)?(KO|WT)(.*)?\\..*", "(^.*)?((KO|WT)(.*)?)\\..*", "WT")
  # list of Batch UN comparisons
  xt3 <- xt %>%
    .selRows("^B[0-9]_(.*)?(KO_Un)\\..*", "^()(B1.*KO|B2.*KO)(.*?)\\..*", "B2")
  xt4 <- xt %>%
    .selRows("^B[0-9]_(.*)?(WT_Un)\\..*", "^()(B1.*WT|B2.*WT)(.*?)\\..*", "B2")
  xt_all <- do.call(rbind, list(xt1, xt2, xt3, xt4))
  xt_all <- do.call(rbind, list(xt2))
  if(!file.exists(file.path(outdir, grouping, paste0(seuid, ".wilcox_deg.rds"))) | overwrite){
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
    saveRDS(markers_all, file=file.path(outdir, grouping, paste0(seuid, ".wilcox_deg.rds")))
  } else {
    markers_all <- readRDS(file.path(outdir, grouping, paste0(seuid, ".wilcox_deg.rds")))
  }
  
  dir.create(file.path(outdir, grouping, seuid))
  ids <- sapply(names(markers_all), function(id){
    dat <- markers_all[[id]]
    if(is.null(dat)) return(c(NA, NA))
    grp1 <- paste0(gsub("_.*", "", unique(dat$fullid1)), "_", unique(dat$id1))
    grp2 <- paste0(gsub("_.*", "", unique(dat$fullid2)), "_", unique(dat$id2))
    grpids <- data.frame("grp1"=grp1, "grp2"=grp2)
    fileid <- paste0(gsub("\\.[0-9]*$", "", id), ".", grp1, "_vs_", grp2) %>%
      gsub("(B[12])_B[12]", "\\1", .)
    if(file.exists(file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv")))){
      print(paste0(id, " - ", fileid))
      return(grpids)
    }
    write.table(dat, 
                file=file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
    return(grpids)
  }) %>%
    t %>% as.data.frame
  ggs <- lapply(c('7d', '3d', '[3d|7d]'), function(grp){
    idx <- which(apply(ids, 1, function(i) all(grepl(paste0("[KO|WT]_", grp), i))))
    markers_sub <- lapply(markers_all[idx], function(i){
      i %>% 
        dplyr::filter(wilcox.p_val_adj <= 0.05,
                      abs(avg_log2FC) >= 0.5,
                      biotype == 'protein_coding',
                      !grepl("Rik$", gene),
                      !grepl("^Gm", gene))
    })
    names(markers_sub) <- gsub("\\.[0-9]*$", "", names(markers_sub))
    idmap <- list('Treg'='Treg',
                  'Treg1'='Treg1', #grep('treg', names(markers_sub), value=T, ignore.case = T),
                  'B'='B_cell', 'Cd4'='Cd4', 'Cd3_DN'='Cd3_DN',
                  'Cd8'=c('Cd8', 'Cd8_cycling'),
                  'ILC2'='ILC2',
                  'contamination'=grep('contam', names(markers_sub), value=T, ignore.case = T),
                  'NK'='Nk_cells',
                  'moMacrophages'=grep('monocyte|macrophage', names(markers_sub), value=T, ignore.case = T)) %>% 
      list2df
    idmap.v <- with(idmap, setNames(Var1, Var2))
    # ctcnts <- sapply(split(markers_sub, idmap.v[names(markers_sub)]), function(i) sum(sapply(i, nrow)))
    ctcnts <- sapply(split(markers_sub, names(markers_sub)), function(i) sum(sapply(i, nrow)))
    
    seu$newid <- as.character(idmap.v[gsub(" ", "_", seu$manual_anno2) %>% gsub("[-+]", "", .)])
    seu$newid <- seu$manual_anno2
    
    
    maxval <- if(seuid == 'LN') 2200 else 800
    maxval <- if(seuid == 'LN') 1000 else 550
    if(maxval <= max(ctcnts)) maxval <- max(ctcnts) * 1.05
    cols <- setNames(viridis_plasma_dark_high[ceiling(ctcnts * 1/(maxval/250))+1],
                     names(ctcnts))
    Idents(seu) <- 'newid'
    pdp <- publicationDimPlot(seu, grp='newid', colors_use=as.character(cols[gsub(" ", "_", levels(Idents(seu)))]), 
                              reduction='umap.mnn', return='list')
    
    
    ggleg <- ggplot(data=data.frame(id=seq_len(maxval)), 
                    aes(x=id, y=id, color=id)) +
      geom_point() +
      scale_color_gradientn(colors = scCustomize::viridis_plasma_dark_high) +
      labs(color='# of DEGs')
    leg <- ggpubr::get_legend(ggleg)
   
    pdp$plot <- pdp$plot + ggtitle(paste0(seuid, " - ", grp))
    plot_figure <- pdp$plot + pdp$axis_plot + ggpubr::as_ggplot(leg) +
      patchwork::plot_layout(design = c(
        patchwork::area(t = 1, l = 2, b = 11, r = 11),
        patchwork::area(t = 10, l = 1, b = 12, r = 2),
        patchwork::area(t = 1, l = 11, b = 4, r = 12))) & 
      theme(aspect.ratio=pdp$aspect.ratio)
    plot_figure
  })
  
  return(ggs)
})
pdf("~/xfer/ST2_LN_degUmaps.pdf")
publicationDimPlot(seus$LN, grp='manual_anno2', reduction='umap.mnn')
ggs
dev.off()

pdf("~/xfer/ST2_Tumor_degUmaps.pdf")
publicationDimPlot(seus$Tumor, grp='manual_anno2', reduction='umap.mnn')
ggs # gg_degumaps
dev.off()


#--- b) ST2: TRegs ----
grouping <- 'st2_tregs'
dir.create(file.path(outdir, grouping), recursive = T, showWarnings = F)
seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
overwrite <- F

## FindAllMarkers on everything and run DotPlot
grp <- 'treg_anno'
allmarkers <- lapply(names(seus_treg), function(seuid){
  seu <- seus_treg[[seuid]]
  Idents(seu) <- grp
  seuid_markers <- FindAllMarkers(object = seu) %>%
    scCustomize::Add_Pct_Diff() 
  return(seuid_markers)
})
names(allmarkers) <- names(seus_treg)

for(seuid in names(seus_treg)[1:2]){
  amarkers <- allmarkers[[seuid]] %>%
    dplyr::filter(pct_diff > 0.2)
  top_markers <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 15, 
                                     named_vector = FALSE, make_unique = TRUE)
  top_markers_tbl <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 50, 
                                     named_vector = FALSE, make_unique = TRUE)
  seu <- seus_treg[[seuid]]
  Idents(seu) <- grp
  seu <- subset(seu, idents=grep(".+", unique(Idents(seu)), value=T))
  
  genes <- read.csv(file.path(PDIR, "ref", "dotplots", paste0("ST2_", seuid, "_TReg.csv")),
                    header = T, sep=",")
  
  write.table(amarkers %>% dplyr::filter(gene %in% top_markers_tbl),
              file=file.path("~/xfer", paste0("ST2treg_", seuid, ".csv")), 
              sep=",", col.names = T, row.names = F, quote = F)
  
  pdf(file.path("~/xfer", paste0("ST2treg_", seuid, ".pdf")), height = 15)
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
  
  pdf(file.path("~/xfer", paste0("ST2treg_", seuid, "_umap.pdf")), width = 7, height = 7)
  publicationDimPlot(seu, grp=grp, reduction='umap.mnn', 
                     simplify_labels=if(grp=='treg_anno2') T else F) %>% plot()
  dev.off()
  
  cat(paste0("xfer ST2treg_", seuid, ".csv\n"))
  cat(paste0("xfer ST2treg_", seuid, "_umap.pdf\n"))
  cat(paste0("xfer ST2treg_", seuid, ".pdf\n"))
  
}


# Output a rds containing all DEGs
# and a csv file containing the tables for the comparisons
method <- 'wilcox'
sigdig <- 3
gg_degumaps <- lapply(names(seus_treg), function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus_treg[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno2)
  uids <- unique(seu$group)
  x <- split(uids, f=gsub("^.*\\.", "", uids))
  xt <- .createGrps(x) %>% 
    t %>% as.data.frame 
  
  # list of KO vs WT
  xt2 <- xt %>%
    .selRows("(^.*)?(KO|WT)(.*)?\\..*", "(^.*)?((KO|WT)(.*)?)\\..*", "WT")
  xt_all <- do.call(rbind, list(xt2)) %>%
    dplyr::filter(!grepl("Un$", id1))
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
  
  dir.create(file.path(outdir, grouping, seuid))
  ids <- sapply(names(markers_all), function(id){
    dat <- markers_all[[id]]
    if(is.null(dat)) return(c(NA, NA))
    grp1 <- paste0(gsub("_.*", "", unique(dat$fullid1)), "_", unique(dat$id1))
    grp2 <- paste0(gsub("_.*", "", unique(dat$fullid2)), "_", unique(dat$id2))
    grpids <- data.frame("grp1"=grp1, "grp2"=grp2)
    fileid <- paste0(gsub("\\.[0-9]*$", "", id), ".", grp1, "_vs_", grp2) %>%
      gsub("(B[12])_B[12]", "\\1", .)
    if(file.exists(file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv")))){
      print(paste0(id, " - ", fileid))
      return(grpids)
    }
    write.table(dat, 
                file=file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
    return(grpids)
  }) %>%
    t %>% as.data.frame
  ggs <- lapply(c('7d', '3d', '[3d|7d]'), function(grp){
    idx <- which(apply(ids, 1, function(i) all(grepl(paste0("[KO|WT]_", grp), i))))
    markers_sub <- lapply(markers_all[idx], function(i){
      i %>% 
        dplyr::filter(wilcox.p_val_adj <= 0.05,
                      abs(avg_log2FC) >= 0.5,
                      biotype == 'protein_coding',
                      !grepl("Rik$", gene),
                      !grepl("^Gm", gene))
    })
    names(markers_sub) <- gsub("\\.[0-9]*$", "", names(markers_sub))
    # ctcnts <- sapply(split(markers_sub, idmap.v[names(markers_sub)]), function(i) sum(sapply(i, nrow)))
    ctcnts <- sapply(split(markers_sub, names(markers_sub)), function(i) sum(sapply(i, nrow)))
    seu$newid <- seu$treg_anno2
    
    
    maxval <- if(seuid == 'LN') 500 else 500
    cols <- setNames(viridis_plasma_dark_high[ceiling(ctcnts * 1/(maxval/250))+1],
                     names(ctcnts))
    Idents(seu) <- 'newid'
    subids <- setdiff(gsub("-", "", levels(Idents(seu))), names(ctcnts))
    setdiff(letters[1:3], letters[1:3])
    if(length(subids) > 0){
      seu <- subset(seu, cells=Cells(seu)[gsub("-", "", seu$treg_anno2) %in% names(ctcnts)])
    }
    
    pdp <- publicationDimPlot(seu, grp='newid', colors_use=as.character(cols[gsub("-", "", levels(Idents(seu)))]), 
                              reduction='umap.mnn', return='list')
    
    
    ggleg <- ggplot(data=data.frame(id=seq_len(maxval)), 
                    aes(x=id, y=id, color=id)) +
      geom_point() +
      scale_color_gradientn(colors = scCustomize::viridis_plasma_dark_high) +
      labs(color='# of DEGs')
    leg <- ggpubr::get_legend(ggleg)
    
    pdp$plot <- pdp$plot + ggtitle(paste0(seuid, " - ", grp))
    plot_figure <- pdp$plot + pdp$axis_plot + ggpubr::as_ggplot(leg) +
      patchwork::plot_layout(design = c(
        patchwork::area(t = 1, l = 2, b = 11, r = 11),
        patchwork::area(t = 10, l = 1, b = 12, r = 2),
        patchwork::area(t = 1, l = 11, b = 4, r = 12))) & 
      theme(aspect.ratio=pdp$aspect.ratio)
    plot_figure
  })
  
  return(ggs)
})
names(gg_degumaps) <- names(seus_treg)


#############################
#### 6. Pathway analysis ####
outdir_pway <- file.path(PDIR, "results", "pathway")
dir.create(outdir_pway, recursive = F)
method <- 'aucell'
geneset_opt <- 'msigdb'
analyze.custom <- FALSE
mode <- 'treg'
runth17 <- FALSE

sigdig <- 4
pval_sig <- 0.05
fc_sig <- 0.5
overwritepathway <- F
min.geneset.size = 10

if(mode == 'treg'){
  seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
} else if(mode == 'all'){
  seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
}

th17dir <- file.path(PDIR, "ref", "th17_pathways")
th17_gs <- lapply(list.files(th17dir), function(id){
  read.csv(file.path(th17dir, id), col.names = 'entrez_gene') %>%
    mutate(gs_name = gsub(".csv", "", id)) %>%
    relocate(., gs_name)
}) %>% do.call(rbind, .)

geneset_l <- lapply(names(msig_lvls),function(mlvl){
  lapply(msig_lvls[[mlvl]], function(sublvl){
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
      dplyr::select(gs_name, entrez_gene) %>%
      as.data.frame()  %>%
      mutate(entrez_gene = gm$ENTREZ$SYMBOL[entrez_gene]) %>%
      dplyr::filter(!is.na(entrez_gene)) 
    keep_gs <- names(which(table(msig_ds$gs_name) >= min.geneset.size))
    msig_ds <- msig_ds %>%
      dplyr::filter(gs_name %in% keep_gs)
    return(msig_ds)
  }) %>% 
    setNames(., unlist(msig_lvls[[mlvl]]))
})
names(geneset_l) <- names(msig_lvls)
geneset_l <- unlist(geneset_l, recursive=F)
names(geneset_l)

if(runth17) geneset_l[['TH17']] <- th17_gs

#--- a) AUCell ----
if(mode=='treg'){
  seu <- seus$LN_Tumor
  DefaultAssay(seu) <- 'RNA'
} else if(mode=='all'){
  seu <- merge(seus[[1]], seus[-1])
  DefaultAssay(seu) <- 'RNA'
  seu <- JoinLayers(seu)
}
expr <- GetAssayData(seu, slot='counts')
expr <- expr[rowSums(expr)>=50, ]

# geneset_all <- do.call(rbind, geneset_l)
scores_l <- list()
for(id in names(geneset_l)){
  print(paste0(id, "..."))
  dir.create(file.path(outdir_pway, method), recursive = T, showWarnings = F)
  outf <- file.path(outdir_pway, method, paste0(method, ".", geneset_opt, id, ".", mode, ".rds"))
  # msig_ds <- data.frame(gs_name=gsub("[0-9]*$", "", names(unlist(geneset_l[id]))),
  #                       entrez_gene=as.character(unlist(geneset_l[id]))) %>%
  
  if(!file.exists(outf)){
    print("Running aucell..")
    if(method=='aucell'){
      scores <- aucellFun(msig_ds = geneset_l[[id]], 
                          expr_mat=expr, gm,
                          mapfrom='SYMBOL', mapto='SYMBOL')
    } else if(method=='pagoda2'){
      scores <- cal_pagoda2(expr, msig_l[[id]])
    }
    saveRDS(scores, file=outf)
  } else {
    print("Reading in aucell..")
    scores <- readRDS(outf)
  }
  scores_l[[id]] <- scores
}
# scores <- readRDS(file.path(outdir_pway, method, paste0("aucell.Custom.rds")))
# scores_l[['Custom']] <- scores
saveRDS(scores_l, file=file.path(outdir_pway, method, paste0("aucell.", geneset_opt, ".", mode, ".rds")))



#--- b) SCPA ----
data=as(cd8seu[['RNA']]$data, 'dgCMatrix')
cd8seu[['RNA2']] <- CreateAssayObject(data[which(rowSums(data) > 50),])

msig_l_scpa <- msig_ds %>% 
  rename_with(., ~'Pathway', .cols='gs_name') %>% 
  dplyr::mutate('Genes'=msig_ds_genename) %>% 
  filter(!is.na(Genes))
msig_l_scpa <- split(msig_l_scpa %>% 
                       dplyr::select(c('Pathway', 'Genes')),
                     f=msig_l_scpa$Pathway)
custom_l_scpa <- customgs %>%
  magrittr::set_colnames(c('Genes', 'Pathway')) %>%
  split(., f=.$Pathway)

for(compid in names(comp_lists)){
  print(compid)
  spl <- comp_lists[[compid]]
  grpspl <- .createGrps(spl)
  if(analyze.custom){
    outf <- file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".custom.rds"))
    pathways_scpa <- custom_l_scpa
  } else {
    outf <- file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".", geneset_opt, ".rds"))
    pathways_scpa <- msig_l_scpa
  }
  
  
  if(!file.exists(outf)){
    print(paste0("Analyzing in ", compid))
    markers <- lapply(colnames(grpspl), function(grpid){
      ids <- grpspl[,grpid]
      print(paste(ids, collapse="//"))
      
      scpa_out <- compare_seurat(cd8seu,
                                 assay='RNA2',
                                 group1 = 'group', 
                                 group1_population = ids,
                                 pathways = pathways_scpa,
                                 max_genes=500,
                                 parallel=F) %>%
        mutate(Database = gsub(":.*", "", Pathway),
               id1=ids[1],
               id2=ids[2],
               celltype=grpid) 
    }) %>%
      do.call(rbind, .) %>% 
      mutate(celltype=gsub("\\.[0-9]$", "", celltype))
    
    saveRDS(markers, file=outf)
  } else {
    print(paste0("Reading in ", compid))
    markers <- readRDS(file=outf)
  }
  # outf.gp <- file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".", 'gprofiler', ".rds"))
  # markers.gp <- readRDS(file=outf.gp)
  # outf.ms <- file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".", 'msigdb', ".rds"))
  # markers.ms <- readRDS(file=outf.ms)
  # markers <- rbind(markers.gp, markers.ms) %>%
  #   filter(Pathway %in% unique(msig_ds$gs_name)) %>%
  #   mutate(adjPval=p.adjust(Pval))
  # saveRDS(markers, file=outf)
}


for(compid in names(comp_lists)){
  markers <-  readRDS(file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".", geneset_opt, ".rds")))
  markers_custom <- readRDS(file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".custom.rds")))
  markers <- rbind(markers, markers_custom)
  dupidx <- with(markers, duplicated(paste0(Pathway, id1, id2)))
  if(any(dupidx)){
    stop(paste0(compid, ": Duplicated indices found at ", paste(which(dupidx), collapse=",")))
  } else {
    saveRDS(markers, file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".", geneset_opt, ".rds")))
  }
}



#--- c) Differential AUCell scores ----
dir.create(file.path(outdir_pway, method, "differential"), showWarnings = F)
scores_l <- readRDS(file=file.path(outdir_pway, method, paste0("aucell.", geneset_opt, ".", mode, ".rds")))
if(geneset_opt == 'gprofiler'){
  scores <- lapply(scores_l[-grep("HP|MIRNA", names(scores_l))], assay) %>% 
    do.call(rbind, .)
} else {
  scores <- lapply(scores_l, assay) %>% 
    do.call(rbind, .)
}



## Simple group by group comparison using wilcoxon test
lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu[['AUCell']] <- CreateAssayObject(data = scores[,Cells(seu)])
  DefaultAssay(seu) <- 'AUCell'

  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  x <- split(uids, f=gsub("^.*\\.", "", uids))
  xt <- .createGrps(x) %>% 
    t %>% as.data.frame 
  
  # list of 3d/7d vs Untreated samples per celltype
  xt1 <- xt %>%
    .selRows("(^.*)?(3d|7d|Un)\\..*", "(^.*)?((KO|WT)_(3d|7d|Un))\\..*", "Un")
  # list of KO vs WT
  xt2 <- xt %>%
    .selRows("(^.*)?(KO|WT)(.*)?\\..*", "(^.*)?((KO|WT)(.*)?)\\..*", "WT")
  # list of Batch UN comparisons
  xt3 <- xt %>%
    .selRows("^B[0-9]_(.*)?(KO_Un)\\..*", "^()(B1.*KO|B2.*KO)(.*?)\\..*", "B2")
  xt4 <- xt %>%
    .selRows("^B[0-9]_(.*)?(WT_Un)\\..*", "^()(B1.*WT|B2.*WT)(.*?)\\..*", "B2")
  xt_all <- do.call(rbind, list(xt1, xt2, xt3, xt4))
  if(!file.exists(file.path(outdir_pway, method, "differential", paste0(seuid, ".wilcox_degeneset.rds"))) | overwrite){
    markers_all <- apply(xt_all, 1, function(ids){
      # ids <- xt_all %>% dplyr::filter(id1 == 'KO_7d' & id2 == 'WT_7d' & grepl('eTreg_NLT-like', baselvl)) %>% unlist
      DefaultAssay(seu) <- 'RNA'
      print(paste0("Testing ", ids[4], " - ", ids[5]))
      ids_test <- factor(c(ids[1], ids[2]))
      ids_test <- relevel(ids_test, ids[3])
      ord <- order(ids_test)
      
      Idents(seu) <- 'group'
      print(paste0("Processing: ", paste(ids, collapse=",")))
      M <- tryCatch({
        FindMarkers(seu, assay='AUCell', 
                    group.by='group',
                    ident.1=ids[1], ident.2=ids[2], 
                    logfc.threshold=0, method='wilcox') %>%
          tibble::rownames_to_column("pathway") %>%
          mutate(id.1=ids[c(4,5)[ord[1]]],
                 id.2=ids[c(4,5)[ord[2]]],
                 fullid1=ids[c(1,2)[ord[1]]],
                 fullid2=ids[c(1,2)[ord[2]]],
                 p_val = round(p_val, sigdig),
                 avg_log2FC=round(avg_log2FC, sigdig),
                 p_val_adj = ifelse(p_val_adj < 10^(-1*sigdig), 
                                    p_val_adj, 
                                    round(p_val_adj, sigdig))) %>%
          relocate(., c('id.1', 'id.2'), .after='pathway')
      }, error=function(e){NULL})
      if(is.null(M)){
        MD <- NULL
      } else {
        cells <- Seurat:::IdentsToCells(object = seu, 
                                        ident.1 = ids[1], 
                                        ident.2 = ids[2], 
                                        cellnames.use = Cells(seu)) 
        D <- apply(GetAssayData(seu, assay='AUCell'), 1, function(i){
          x <- effsize::cohen.d(i[cells$cells.1], i[cells$cells.2])
          setNames(c(x$conf.int, x$estimate), c('D.lower', 'D.upper', 'D'))
        }) %>%
          t %>% as.data.frame %>%
          tibble::rownames_to_column("pathway") %>%
          mutate(D.lower = round(D.lower, sigdig),
                 D.upper =round(D.upper, sigdig),
                 D =round(D, sigdig))
        
        MD <- left_join(M, D, by='pathway')
      }
      
      return(MD)
      
    })
    saveRDS(markers_all, file=file.path(outdir_pway, method, "differential", paste0(seuid, ".wilcox_degeneset.rds")))
  } else {
    markers_all <- readRDS(file.path(outdir_pway, method, "differential", paste0(seuid, ".wilcox_degeneset.rds")))
  }
  
  dir.create(file.path(outdir_pway, method, "differential", "tables"), showWarnings = F)
  for(id in names(markers_all)){
    dat <- markers_all[[id]]
    if(is.null(dat)) next()
    grp1 <- paste0(gsub("_.*", "", unique(dat$fullid1)), "_", unique(dat$id.1))
    grp2 <- paste0(gsub("_.*", "", unique(dat$fullid2)), "_", unique(dat$id.2))
    fileid <- paste0(gsub("\\.[0-9]*$", "", id), ".", grp1, "_vs_", grp2) %>%
      gsub("(B[12])_B[12]", "\\1", .)
    # if(file.exists(file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv")))){
    #   print(paste0(id, " - ", fileid))
    #   stop()
    # }
    outf <- file.path(outdir_pway, method, "differential", "tables", 
                      paste0(seuid, ".DEgeneset_", fileid, ".csv"))
    if(overwrite | !file.exists(outf)){
      write.table(dat, file=outf,
                  sep=",", col.names = T, row.names = F, quote = F)
    }
  }
})


## Custom grouped analysis that uses the two-way RM ANOVA to test between groups
treg_markers_all <- lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu[['AUCell']] <- CreateAssayObject(data = scores[,Cells(seu)]) #th17 c(5073:5078)
  DefaultAssay(seu) <- 'AUCell'
  
  if(mode=='treg'){
    seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
    uids <- unique(seu$group)
  } else if(mode == 'all'){
    seu$group <- paste0(seu$orig.ident, ".", seu$manual_anno2)
    uids <- unique(seu$group)
    uids <- grep("Cd8|Cd4|B cell|DC", uids, value=T, ignore.case = T)
  }
  
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))
  
  outf <- file.path(outdir_pway, method, "differential", paste0(seuid, ".", mode, ".rmaov_degeneset.rds"))
  if(!file.exists(outf) | overwrite){
    # markers_all <- lapply(xt_all, 1, function(ids){
    markers_all <- lapply(meta_l, function(meta){
      message(paste0(unique(meta$celltype), ": ", unique(meta$batch)))
      # Make all comparisons treated vs Untreated
      meta$IDs <- factor(as.character(meta$IDs), 
                         levels=meta$IDs[order(grepl("Un", meta$IDs))])
      idents <- lapply(split(meta$IDs, f=meta$Condition), function(i) as.character(sort(i)))
      
      Idents(seu) <- 'group'
      tryCatch({
        FindMarkers.grouped(seu,meta=meta, idents=idents, assay='AUCell',
                            fc.method='cohensd') %>%
          tibble::rownames_to_column(., "Pathway") %>%
          relocate(., c(Pathway))
      }, error=function(e){NULL})
    })
    saveRDS(markers_all, file=outf)
  } else {
    markers_all <- readRDS(outf)
  }
  
  
  for(id in names(markers_all)){
    dat <- markers_all[[id]]
    if(is.null(dat)) next()
    grp1 <- 'KO'
    grp2 <- 'WT'
    fileid <- paste0(gsub("\\/", "-", id), ".", grp1, "_vs_", grp2)
    outf <- file.path(outdir_pway, method, "differential", "tables",
                      paste0(seuid, ".", mode, ".DEregulon_", fileid, ".csv"))
    if(overwrite | !file.exists(outf)){
      write.table(dat, file=outf,
                  sep=",", col.names = T, row.names = F, quote = F)
      # write.table(list2df(markers_all),
      #             file=file.path(outdir_pway, method, "differential", "tables", paste0(seuid, ".DEpathway_th17.csv")),
      #             col.names = T, row.names = F, quote = F)
      # file.copy(file.path(outdir_pway, method, "differential", "tables",paste0(seuid, ".DEpathway_th17.csv")),
      #           to="~/xfer", overwrite = T)
    }
  }
  
  return(markers_all)
})
names(treg_markers_all) <- names(seus_treg)[1:2]

#--- d) Cytoscape visualization of differential pathways ----
cytodir <- file.path(PDIR, "results", "cytoscape")
geneset_id_spl <- split(geneset_l, f=gsub(":.*", "", names(geneset_l)))
gs_map <- lapply(geneset_id_spl, function(i){
  data.frame("ID"=gsub("^(.*?)_.*$", "\\1", names(i)) %>%
               gsub(":[0-9]+$", "", .),
             "Desc"=gsub("^.*?_(.*)$", "\\1", names(i)),
             "setSize"=sapply(i, length)) %>%
    magrittr::set_rownames(gsub("_", "-", names(i))) %>%
    mutate(ID = gsub("REAC:.*", "REAC", ID))
}) %>% do.call(rbind, .)
dup <- duplicated(rownames(gs_map))
if(any(dup)) gs_map <- gs_map[-which(dup),]
rownames(gs_map) <- gsub("^.*?\\.", "", rownames(gs_map)) %>% gsub("_", "-", .)


for(compid in rev(names(comp_lists))){
  print(paste0(">>>    ", compid, "    <<<"))
  pathway_markers_f <- paste0("Pathway_", compid, ".", geneset_opt, ".rds")
  markers <- readRDS(file.path(outdir_pway, compid, pathway_markers_f))
  
  for(label in names(markers)){
    fileid <- .getFileid(grpid = label, df=markers[[label]])
    dat_raw=markers[[label]] %>%
      magrittr::set_rownames(., markers[[label]]$pathway)
    
    dat <- findmarkers2CytoscapeFormat(dat=dat_raw, gs_map=gs_map) 
    dat$up[is.na(dat$up)] <- 0
    dat$down[is.na(dat$down)] <- 0
    
    dir.create(file.path(cytodir,  compid, geneset_opt), showWarnings = F, recursive = T)
    write.table(dat$up, 
                file=file.path(cytodir,  compid, geneset_opt,
                               paste0("cytoscape.",label, "_", fileid, ".up.xls")),
                sep="\t", col.names = T, row.names = F, quote = F)
    write.table(dat$down, 
                file=file.path(cytodir,  compid, geneset_opt,
                               paste0("cytoscape.",label, "_", fileid, ".down.xls")),
                sep="\t", col.names = T, row.names = F, quote = F)
    write.table(dat_raw %>%
                  tibble::rownames_to_column("Pathway") %>%
                  dplyr::relocate(),
                file=file.path(cytodir,  compid, geneset_opt,
                               paste0("pathway_deg.",label, "_", fileid, ".csv")),
                sep=",", col.names = T, row.names = F)
  }
}


## Create the reference data for cytoscape plots
geneset_spl <- split((geneset_l), f=gsub(":.*", "",names(geneset_l)))
lapply(names(geneset_spl), function(msig_id){
  print(msig_id)
  dir.create(file.path(cytodir, "reference"), showWarnings=F)
  CON <- file(file.path(cytodir, "reference", 
                        paste0(geneset_opt, "2023.", msig_id, ".gmt")), "a")    #Open connection to append
  lapply(names(geneset_spl[[msig_id]]), function(geneset_id){
    line <- c(gsub("[_]", "-", toupper(geneset_id)),
              gsub("[_]", "-", toupper(geneset_id)) %>% gsub("^.*?:.*?[0-9]+-", "", .),
              geneset_spl[[msig_id]][[geneset_id]]) %>%
      as.character %>% 
      paste(., collapse="\t")
    writeLines(paste0(line), CON)
  })
  close(CON)       
})

#######################################
#### 7. SCENIC upstream regulators ####
outdir <- file.path(PDIR, "results", "regulons", "loom")
dir.create(outdir, recursive=T, showWarnings = F)
overwrite_regulon <- F

sigdig <- 4
pval_sig <- 0.05
fc_sig <- 0.5
overwritepathway <- F
min.geneset.size = 10
mode <- 'treg'

if(mode == 'treg'){
  seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
  seu <- seus$LN_Tumor
  DefaultAssay(seu) <- 'RNA'
} else if(mode == 'all'){
  seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
  seu <- merge(seus[[1]], seus[-1])
  DefaultAssay(seu) <- 'RNA'
  seu <- JoinLayers(seu)
} else if(mode == 'cd8'){
  seu <- readRDS(file=file.path(datadir, "seurat_obj", "CD45_Tumor_CD8.seuratobj.rds"))
  seu$treg_anno <- seu$manual_anno
  seu@reductions$umap.mnn <- seu@reductions$umap
}

#--- a) Generate loom file for pySCENIC ----
# code from https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html#run-scenic-vis-pyscenic
exprMat <- seu@assays$RNA$data
cellInfo <- seu@meta.data

loci1 <- which(rowSums(exprMat>0) > 0.01*ncol(exprMat))
exprMat_filter <- exprMat[loci1, ]

loomf <- file.path(outdir, paste0("st2.", mode, ".seu2loom.loom"))
if(!file.exists(loomf)){
  loom <- build_loom(loomf, dgem=exprMat_filter)
  loom <- addLoomCellAnnotation(loom, cellInfo)
  close_loom(loom)
}

#--- b) Post pySCENIC analysis ----
outdirregulon <- file.path(PDIR, "results", "regulons", "pyscenic")
dir.create(outdirregulon, showWarnings = F)
regulon_f <- file.path(PDIR, "results", "regulons", "regulons.csv")

looms <- lapply(c('LN', 'Tumor'), function(grp){
  open_loom(file.path(PDIR, "results", "regulons", "pyscenic", paste0(grp, '.seu_scenic.loom')))
}) %>% setNames(c('LN', 'Tumor'))

# Get the regulon genesets
regulons <- lapply(looms, function(loom){
  regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
  regulonsToGeneLists(regulons_incidMat)
})
regulon_ids <- unique(unlist(sapply(regulons, names)))
regulons <- lapply(setNames(regulon_ids, regulon_ids), function(id){
  unique(c(regulons$LN[[id]], regulons$Tumor[[id]]))
})
regulon_df <- data.frame('regulon'=rep(names(regulons), sapply(regulons, length)),
                         'gene'=as.character(unlist(regulons)))
write.table(regulon_df, file=regulon_f, sep=",", 
            col.names = T, row.names = F, quote = F)

# Get regulon AUC scores
regulonAUC <- lapply(looms, function(loom){
  assay(get_regulons_AUC(loom, column.attr.name='RegulonsAUC')) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., "regulons")
}) %>%
  purrr::reduce(., .f=full_join, by='regulons') %>%
  tibble::column_to_rownames(., "regulons") %>%
  as.matrix
regulonAUC <- new(Class = "aucellResults", 
      SummarizedExperiment::SummarizedExperiment(assays = list(AUC = regulonAUC)))
# seu[['regulons']] <- CreateAssayObject(data=assay(regulonAUC[,Cells(seu)]))

#--- c) ORA of regulons ----
regulon_f <- file.path(PDIR, "results", "regulons", "regulons.csv")
regulons <- read.csv(regulon_f, header=T)
regulons <- split(regulons$gene, f=regulons$regulon)

geneset_l <- lapply(names(msig_lvls),function(mlvl){
  lapply(msig_lvls[[mlvl]], function(sublvl){
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
      dplyr::select(gs_name, entrez_gene, gene_symbol) %>%
      as.data.frame()  %>%
      mutate(entrez2symbol_gene = gm$ENTREZ$SYMBOL[as.character(entrez_gene)],
             entrez_gene=gene_symbol) %>%
      dplyr::filter(!is.na(entrez_gene)) 
    keep_gs <- names(which(table(msig_ds$gs_name) >= min.geneset.size))
    msig_ds <- msig_ds %>%
      dplyr::filter(gs_name %in% keep_gs)
    return(msig_ds)
  }) %>% 
    setNames(., unlist(msig_lvls[[mlvl]]))
})
names(geneset_l) <- names(msig_lvls)
geneset_l <- unlist(geneset_l, recursive=F)


geneset_sel <- geneset_l[['C5.GO:CC']][grep("GOCC_CORTICAL_CYTOSKELETON", geneset_l[['C5.GO:CC']]$gs_name),]
geneset_sel3 <- geneset_l[['C5.GO:CC']][grep("GOCC_CORTICAL_CYTOSKELETON", geneset_l[['C5.GO:CC']]$gs_name),]
sublvl <- grep("cortical cytoskeleton$", names(gprof_ds), value=T)
geneset_sel2 <- data.frame("gs_name"=as.character(rep(sublvl, length(gprof_ds[[sublvl]]))),
                           "entrez_gene"=as.character(gm$ENSEMBL$SYMBOL[gprof_ds[[sublvl]]]))
table(geneset_sel3$gene_symbol %in% geneset_sel2$entrez_gene)



geneset_df <-  do.call(rbind, geneset_l)
# geneset_df <-   lapply(names(gprof_ds), function(sublvl){
#   data.frame("gs_name"=sublvl,
#              "entrez_gene"=gprof_ds[[sublvl]])
# }) %>% do.call(rbind,.) %>%
#   mutate(entrez_gene = gm$ENSEMBL$SYMBOL[entrez_gene]) %>%
#   dplyr::filter(!is.na(entrez_gene))

ora_res <- lapply(regulons, function(regulon_i){
  clusterProfiler::enricher(
    gene = regulon_i, # A vector of your genes of interest
    pAdjustMethod = "BH", # Method to be used for multiple testing correction
    universe = Features(seu), # A vector containing your background set genes
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.05,
    TERM2GENE = geneset_df
  )
})
saveRDS(ora_res, file=file.path(PDIR, "results", "regulons", "regulon_ora.rds"))

all_ora_res <- lapply(names(ora_res), function(id){
  ora_res[[id]]@result %>% 
    mutate(ID = id)
}) %>% do.call(rbind, .)
write.table(all_ora_res, file=file.path(PDIR, "results", "regulons", "regulon_ora.csv"),
            sep=",", col.names = T, row.names = F, quote = F)
file.copy(file.path(PDIR, "results", "regulons", "regulon_ora.csv"), 
          to="~/xfer", overwrite = T)

#--- d) Differential Regulon scores ----
dir.create(file.path(outdir, "..", "differential"), showWarnings = F)


looms <- lapply(c('LN', 'Tumor'), function(grp){
  open_loom(file.path(PDIR, "results", "regulons", "pyscenic", paste0(grp, '.seu_scenic.loom')))
}) %>% setNames(c('LN', 'Tumor'))
regulonAUC <- lapply(looms, function(loom){
  assay(get_regulons_AUC(loom, column.attr.name='RegulonsAUC')) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., "regulons")
}) %>%
  purrr::reduce(., .f=full_join, by='regulons') %>%
  tibble::column_to_rownames(., "regulons") %>%
  as.matrix
regulonAUC <- new(Class = "aucellResults", 
                  SummarizedExperiment::SummarizedExperiment(assays = list(AUC = regulonAUC)))

## Measure individual group comparisons (e.g. WT_7d vs WT_Un)
lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu[['regulons']] <- CreateAssayObject(data=assay(regulonAUC[,Cells(seu)]))
  DefaultAssay(seu) <- 'regulons'
  
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  
  x <- split(uids, f=gsub("^.*\\.", "", uids))
  xt <- .createGrps(x) %>% 
    t %>% as.data.frame 
  
  # list of 3d/7d vs Untreated samples per celltype
  xt1 <- xt %>%
    .selRows("(^.*)?(3d|7d|Un)\\..*", "(^.*)?((KO|WT)_(3d|7d|Un))\\..*", "Un")
  # list of KO vs WT
  xt2 <- xt %>%
    .selRows("(^.*)?(KO|WT)(.*)?\\..*", "(^.*)?((KO|WT)(.*)?)\\..*", "WT")
  # list of Batch UN comparisons
  xt3 <- xt %>%
    .selRows("^B[0-9]_(.*)?(KO_Un)\\..*", "^()(B1.*KO|B2.*KO)(.*?)\\..*", "B2")
  xt4 <- xt %>%
    .selRows("^B[0-9]_(.*)?(WT_Un)\\..*", "^()(B1.*WT|B2.*WT)(.*?)\\..*", "B2")
  xt_all <- do.call(rbind, list(xt1, xt2, xt3, xt4))
  if(!file.exists(file.path(outdir, "..", "differential", paste0(seuid, ".wilcox_deregulon.rds"))) | overwrite){
    markers_all <- apply(xt_all, 1, function(ids){
      # ids <- xt_all %>% dplyr::filter(id1 == 'KO_7d' & id2 == 'WT_7d' & grepl('eTreg_NLT-like', baselvl)) %>% unlist
      DefaultAssay(seu) <- 'regulons'
      print(paste0("Testing ", ids[4], " - ", ids[5]))
      ids_test <- factor(c(ids[1], ids[2]))
      ids_test <- relevel(ids_test, ids[3])
      ord <- order(ids_test)
      
      Idents(seu) <- 'group'
      print(paste0("Processing: ", paste(ids, collapse=",")))
      M <- tryCatch({
        FindMarkers(seu, assay='regulons', 
                    group.by='group',
                    ident.1=ids[1], ident.2=ids[2], 
                    logfc.threshold=0, method='wilcox') %>%
          tibble::rownames_to_column("regulons") %>%
          mutate(id.1=ids[c(4,5)[ord[1]]],
                 id.2=ids[c(4,5)[ord[2]]],
                 fullid1=ids[c(1,2)[ord[1]]],
                 fullid2=ids[c(1,2)[ord[2]]],
                 p_val = round(p_val, sigdig),
                 avg_log2FC=round(avg_log2FC, sigdig),
                 p_val_adj = ifelse(p_val_adj < 10^(-1*sigdig), 
                                    p_val_adj, 
                                    round(p_val_adj, sigdig))) %>%
          relocate(., c('id.1', 'id.2'), .after='regulons')
      }, error=function(e){NULL})
      if(is.null(M)){
        MD <- NULL
      } else {
        cells <- Seurat:::IdentsToCells(object = seu, 
                                        ident.1 = ids[1], 
                                        ident.2 = ids[2], 
                                        cellnames.use = Cells(seu)) 
        D <- apply(GetAssayData(seu, assay='regulons'), 1, function(i){
          x <- effsize::cohen.d(i[cells$cells.1], i[cells$cells.2])
          setNames(c(x$conf.int, x$estimate), c('D.lower', 'D.upper', 'D'))
        }) %>%
          t %>% as.data.frame %>%
          tibble::rownames_to_column("regulons") %>%
          mutate(D.lower = round(D.lower, sigdig),
                 D.upper =round(D.upper, sigdig),
                 D =round(D, sigdig))
        
        MD <- left_join(M, D, by='regulons')
      }
      
      return(MD)
      
    })
    saveRDS(markers_all, file=file.path(outdir, "..", "differential", paste0(seuid, ".wilcox_deregulon.rds")))
  } else {
    markers_all <- readRDS(file.path(outdir, "..", "differential", paste0(seuid, ".wilcox_deregulon.rds")))
  }
  
  dir.create(file.path(outdir, "..", "differential", "tables"), showWarnings = F)
  for(id in names(markers_all)){
    dat <- markers_all[[id]]
    if(is.null(dat)) next()
    grp1 <- paste0(gsub("_.*", "", unique(dat$fullid1)), "_", unique(dat$id.1))
    grp2 <- paste0(gsub("_.*", "", unique(dat$fullid2)), "_", unique(dat$id.2))
    fileid <- paste0(gsub("\\.[0-9]*$", "", id), ".", grp1, "_vs_", grp2) %>%
      gsub("(B[12])_B[12]", "\\1", .)
    # if(file.exists(file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv")))){
    #   print(paste0(id, " - ", fileid))
    #   stop()
    # }
    outf <- file.path(outdir, "..", "differential", "tables", 
                      paste0(seuid, ".DEregulon_", fileid, ".csv"))
    if(overwrite | !file.exists(outf)){
      write.table(dat, file=outf,
                  sep=",", col.names = T, row.names = F, quote = F)
    }
  }
})


## Measure two-way group comparisons (e.g. [KO_treated-vs-untreated] vs [WT_treated-vs-untreated])
treg_markers_all <- lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu[['regulons']] <- CreateAssayObject(data=assay(regulonAUC[,Cells(seu)]))
  DefaultAssay(seu) <- 'regulons'
  
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))
  
  outf <- file.path(outdir, "..", "differential", paste0(seuid, ".rmaov_deregulon.rds"))
  if(!file.exists(outf) | overwrite){
    # markers_all <- lapply(xt_all, 1, function(ids){
    markers_all <- lapply(meta_l, function(meta){
      message(paste0(unique(meta$celltype), ": ", unique(meta$batch)))
      # Make all comparisons treated vs Untreated
      meta$IDs <- factor(as.character(meta$IDs), 
                         levels=meta$IDs[order(grepl("Un", meta$IDs))])
      idents <- lapply(split(meta$IDs, f=meta$Condition), function(i) as.character(sort(i)))
      
      Idents(seu) <- 'group'
      tryCatch({
        FindMarkers.grouped(seu,meta=meta, idents=idents, assay='regulons',
                            fc.method='cohensd') %>%
          tibble::rownames_to_column(., "Pathway") %>%
          relocate(., c(Pathway))
      }, error=function(e){NULL})
    })
    saveRDS(markers_all, file=outf)
  } else {
    markers_all <- readRDS(outf)
  }
  
  for(id in names(markers_all)){
    dat <- markers_all[[id]]
    if(is.null(dat)) next()
    grp1 <- 'KO'
    grp2 <- 'WT'
    fileid <- paste0(gsub("\\/", "-", id), ".", grp1, "_vs_", grp2)
    outf <- file.path(outdir, "..", "differential", "st2_treg", "tables",
                      paste0(seuid, ".DEregulon_", fileid, ".csv"))
    if(overwrite | !file.exists(outf)){
      write.table(dat, file=outf,
                  sep=",", col.names = T, row.names = F, quote = F)
    }
  }
  
  return(markers_all)
})
names(treg_markers_all) <- names(seus_treg)[1:2]

#################################################################
#### 8. Trajectory analysis - MONOCLE3 ####
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

##############################################
#### 9. Velocity and Trajectory analysis ####
mode <- 'cd8' #'treg', 'cd8'
if(mode == 'cd8'){
  seu <- readRDS(file=file.path(datadir, "seurat_obj", "CD45_Tumor_CD8.seuratobj.rds"))
  seu$treg_anno <- seu$manual_anno
  seu@reductions$umap.mnn <- seu@reductions$umap
  seus <- list("Tumor"=seu)
} else if(mode == 'treg'){
  seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
  seus$LN_Tumor$treg_anno <- seus$LN_Tumor$treg_anno2
  seus <- seus[c('LN', 'Tumor')]
}
outdir <- file.path(PDIR, "results", "scVelo", mode)
dir.create(outdir, recursive=T, showWarnings = F)


outrds <- file.path(datadir, "velocyto", "seu_velocyto.rds")
dir.create(file.path(outdir, "scVelo"), showWarnings = F)
make.small <- FALSE


#--- a) scVelo ----
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
raw_velocyto_f <- gsub(".rds", "_raw.rds", outrds)
if(!file.exists(raw_velocyto_f)){
  sample_ids <- list.files(velocyto_dir, pattern = "[^merge]") %>% 
    grep("loom_files|.rds", ., invert = T, value=T)
  seu_velo_raw <- lapply(setNames(sample_ids,sample_ids), function(sid){
    print(paste0(sid, " => ", .relabelid(sid)))
    f <- list.files(file.path(velocyto_dir, sid))
    loom.data <- ReadVelocity(file = file.path(velocyto_dir, sid, f))
    loom.data <- lapply(loom.data, function(loom.mat){
      rownames(loom.mat) <- make.unique(rownames(loom.mat))
      loom.mat
    })
    loom.data <- as.Seurat(loom.data)
    loom.data$orig.ident <- rep(as.character(sid), nrow(loom.data))
    
    new.names <- gsub("^.*\\:", paste0(.relabelid(sid), "_"), Cells(loom.data)) %>% 
      gsub("x$", "", .)
    loom.data <- RenameCells(loom.data, 
                             new.names=new.names)
    return(loom.data)
  })
  saveRDS(seu_velo_raw, file=raw_velocyto_f)
} else {
  seu_velo_raw <- readRDS(file=raw_velocyto_f)
}

## Preprocess the seurat velocyto files
names(seu_velo_raw) <- .relabelid(names(seu_velo_raw))
seu_velo <- merge(seu_velo_raw[[1]], seu_velo_raw[-1], 
                  project = 'ST2_TDLN_LN')
if(mode=='treg'){
  seu_velo_l <- list("LN_Tumor"=seu_velo,
                     "Tumor"=seu_velo_raw[grep("ST2_Tumor", names(seu_velo_raw), value=T)],
                     "LN"=seu_velo_raw[grep("ST2_LN", names(seu_velo_raw), value=T)])
  cols <- c('orig.ident', 'treg_anno', 'predicted.cl_annot', 'mnn_clusters')
} else if(mode == 'cd8'){
  seu_velo_l <- list("LN_Tumor"=seu_velo,
                     "Tumor"=seu_velo_raw[grep("CD45_Tumor", names(seu_velo_raw), value=T)],
                     "LN"=seu_velo_raw[grep("CD45_LN", names(seu_velo_raw), value=T)])
  cols <- c('orig.ident', 'treg_anno', 'seurat_clusters')
}
rm(seu_velo, seu_velo_raw); gc()
# seu_velo <- seu_velo_raw[[1]]

## Generate the flat files used in scVelo
x <- lapply(names(seu_velo_l), function(seu_velo_id){
  print(paste0(seu_velo_id, " - ", mode))
  seu <- seus[[seu_velo_id]]
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

  ## Unify the cells between the velocyto and seurat obj
  newlabel <- .relabelid(Cells(seu))
  seu <- RenameCells(seu, new.names=newlabel)
  seu_small <- subset(seu,
                      cells=Cells(seu_velo))
  seu_velo_small <- subset(seu_velo, cells=Cells(seu_small))
  
  ### output for scVelo velocity analysis
  print("saving...")
  # flatfile method: https://smorabit.github.io/tutorials/8_velocyto/
  id <- paste0(seu_velo_id, ".", mode)
  DefaultAssay(seu_small) <- 'RNA'
  seu_small <- tryCatch({
    JoinLayers(seu_small)
  }, error=function(e){seu_small})
  seu_small@meta.data <- seu_small@meta.data[,cols]
  if(!"pca" %in% Reductions(seu_small)){
    seu_small <- RunPCA(seu_small)
  }
  writeSeuToFlat(seu_small, assay='RNA',reduction = 'umap.mnn',
                 out_metaf=file.path(outdir, paste0(id, ".metadata.csv")),
                 out_cntsf=file.path(outdir, paste0(id, ".counts.mtx")),
                 out_pcaf=file.path(outdir, paste0(id, ".pca.csv")),
                 out_featuref=file.path(outdir, paste0(id, ".genes.csv")))
})

#--- b) Condiments ----
## Condiments with monocle3 package
if(mode == 'cd8'){
  # trajs <- list('Tumor'=list('P1'=c('Y_2', NA)))
  trajs <- list('Tumor'=list('P1'=c('Y_14', NA)))
  comps <- list("CD45_KO"=c("CD45", "KO"),
                "CD45_WT"=c("CD45", "WT"))
} else if (mode =='treg'){
  trajs <- list('LN'=list('P1all'=c('Y_100',NA),
                          'P2all'=c('Y_252',NA)),
                'Tumor'=list('P1all'=c('Y_120', NA),
                             'P2all'=c('Y_260', NA)))
  comps <- list("B1_WT"=c("B1", "WT"),
                "B1_KO"=c("B1", "KO"),
                "B2_WT"=c("B2", "WT"),
                "B2_KO"=c("B2", "KO"))
}

seu_traj_imbalance <- lapply(names(seus), function(seuid){
  stopifnot(!file.exists(file.path(outdir, "seu_traj_imbalance.rds")))
  ### Generate trajectory analysis
  seu <- seus[[seuid]]
  Idents(seu) <- 'orig.ident'

  
  set.seed(1234)
  ncenter <- ifelse(seuid =='LN',  300, 275)
  cds <- runMonocle3(seu, reduction='umap.mnn', ncenter=ncenter)
  
  message("[Traj.Analysis] Adding monocle3 nodes")
  cds <- runMonocle3(cds, nodes_l=trajs[[seuid]])
  
  message("[Traj.Analysis] Calculating imbalance scores")
  scores_all <- lapply(names(comps), function(compid){
    message(compid)
    comp <- comps[[compid]]
    idx <- grep(paste0("^", comp[1], ".*", comp[2]), cds$conditions)
    
    message("Running imbalance")
    scores <- condiments::imbalance_score(Object = reducedDims(cds$cds[,idx])$UMAP,
                                          conditions =  colData(cds$cds[,idx])$orig.ident,
                                          k = 20, smooth = 40)
    message("Adding imbalance")
    df <- as.data.frame(reducedDims(cds$cds[,idx])$UMAP) %>%
      mutate("score"=as.numeric(round(scores$scores,3)),
             "scaled_scores"=as.numeric(round(scores$scaled_scores,3)),
             'scaled_scores2' = scales::rescale(scores$scaled_scores, to=c(0,1))) 
    return(df)
  }) %>% setNames(., names(comps))
  df <- do.call(rbind, scores_all)
  rownames(df) <- gsub("^.*\\.", "", rownames(df))
  compids <- setNames(rep(names(scores_all), sapply(scores_all, nrow)), 
                          as.character(unlist(sapply(scores_all, rownames))))

  cds$cds$imbalance_score <- df[colnames(cds$cds),]$scaled_scores
  cds$cds$imbalance_score_size <- abs(cds$cds$imbalance_score) %>%
    ifelse(. > 2, ., 0) %>%
    scales::rescale(., to=c(0.35, 2))
  cds$cds$group <- compids[colnames(cds$cds)]
  
  ### Create a dataframe of pseudotime and imbalance score, annotated per comparison
  message("[Traj.Analysis] Generating imbalance score per pseudotime bin")
  resdf <- lapply(names(comps), function(compid){
    message(compid)
    comp <- comps[[compid]]
    conditions <- unique(grep(paste0("^", comp[1], ".*", comp[2]), cds$conditions, value=T))
    all_pseudotime <- unlist(as.numeric(cds$pseudotime))
    pseudotime_breaks <- seq(0, max(all_pseudotime[!is.infinite(all_pseudotime)]), by=0.2)
    
    lapply(conditions, function(condition){
      lapply(colnames(cds$pseudotime), function(resid){
        idx <- which(cds$conditions == condition)
        resdf <- cds$pseudotime[idx,resid,drop=F]
        cells_x <- intersect(rownames(resdf),
                             rownames(cds$cellWeights)[which(cds$cellWeights[,resid,drop=F] == 1)])
        # cells_x <- rownames(resdf)[which(!is.infinite(resdf[,1]))]
        
        resdf[cells_x,,drop=F] %>%
          as.data.frame %>%
          magrittr::set_colnames(., "pseudotime") %>%
          mutate(conditions=condition, 
                 partition=resid,
                 pseudotime_breaks = cut(pseudotime, pseudotime_breaks),
                 lineages="lineage1",
                 sample=factor(colData(cds$cds[,cells_x])$orig.ident),
                 anno=factor(colData(cds$cds[,cells_x])$treg_anno),
                 imbalance_score=colData(cds$cds[,cells_x])$imbalance_score,
                 comp=compid) 
      })  %>% 
        do.call(rbind, .)
    }) %>%
      do.call(rbind, .)
  }) %>% do.call(rbind, .)
  resdf$conditions <- gsub("^.*(3d|7d|Un).*", "\\1", resdf$conditions)
  
  message("[Traj.Analysis] Getting proportion of celltypes per pseudotime bin")
  break_prop <- lapply(split(resdf, f=resdf$partition), .breakPseudotime,
                       torange=c(0,1), calculate='anno') %>%
    do.call(rbind, .)
  
  message("[Traj.Analysis] Getting proportion of samples per pseudotime bin")
  break_sample_prop <- lapply(split(resdf, f=resdf$partition), .breakPseudotime,
                              torange=NULL, calculate='sample', 
                              method='chisq', ref_proportions=(table(colData(cds$cds)$orig.ident) / ncol(cds$cds))) %>%
    do.call(rbind, .)
  
  message("[Traj.Analysis] Testing the differences of imbalance scores between KO and WT")
  batches <- if(mode=='cd8') c('CD45') else c('B1', 'B2')
  imbalance_breaks <- lapply(batches, function(batchid){
    koid <- paste0(batchid, "_KO")
    wtid <- paste0(batchid, "_WT")
    
    lapply(split(resdf, f=resdf$partition), function(partspl){
      partitionid <- unique(partspl$partition)
      compspl <- split(partspl, f=partspl$comp)
      
      ## Get null distribution for difference between imbalance scores
      nullrange <- sapply(c(1:1000), function(seed){
        set.seed(seed)
        
        sample(compspl[[koid]]$imbalance_score) - sample(compspl[[wtid]]$imbalance_score)
      })
      null_ecdf <- ecdf(as.numeric(nullrange))
      
      ## Calculate a pvalue and sumarrize statistics
      kospl <- split(compspl[[koid]], compspl[[koid]]$pseudotime_breaks)
      wtspl <- split(compspl[[wtid]], compspl[[wtid]]$pseudotime_breaks)
      uids <- unique(c(names(kospl), names(wtspl)))
      test_stats <- sapply(uids, function(id_i){
        x1 <- kospl[[id_i]]$imbalance_score
        x2 <- wtspl[[id_i]]$imbalance_score
        tryCatch({
          delta <- mean(x1)-mean(x2)
          res <- t.test(x1, x2)
          c('delta'=delta, 
            unlist(res[c('conf.int', 'p.value', 'estimate')]),
            'p.null'=null_ecdf(delta))
        }, error=function(e){
          rep(NA, 7) %>% 
            setNames(., c('delta', 'conf.int1',  'conf.int2',
                          'p.value', 'estimate.mean of x',  
                          'estimate.mean of y', 'p.null'))
        })
      }) %>% t %>% as.data.frame %>%
        mutate(partition=partitionid,
               batchid=batchid) %>%
        tibble::rownames_to_column(., "pseudotime_breaks")
      return(test_stats)
    }) %>%
      do.call(rbind, .)
  }) %>% do.call(rbind, .)  %>%
    mutate(pseudotime=.break2number(pseudotime_breaks),
           twotailed.p = (0.5 - abs(p.null - 0.5))*2,
           p=.to_simple_p(twotailed.p))
    
  imbalance_breaks <- imbalance_breaks %>%
    dplyr::filter(!is.na(delta))
  imbalance_breaks$pseudotime <- as.numeric(imbalance_breaks$pseudotime)

  return(list("cds"=cds,
              "imbalance_breaks"=imbalance_breaks,
              "break_prop"=break_prop,
              "break_sample_prop"=break_sample_prop,
              "resdf"=resdf))
})
seu_traj_imbalance[[seuid]] <- list("cds"=cds,
                                    "imbalance_breaks"=imbalance_breaks,
                                    "break_prop"=break_prop,
                                    "break_sample_prop"=break_sample_prop,
                                    "resdf"=resdf)
# saveRDS(seu_traj_imbalance, file=file.path(outdir, "seu_traj_imbalance.rds"))
seu_traj_imbalance <- readRDS(file=file.path(outdir, "seu_traj_imbalance.rds"))

imbalance_comp_gg <- lapply(seu_traj_imbalance, function(seuid_traj){
  sig.cols <- c('n.s.'='black', 'p<0.001'='#e31a1c', 'p<0.05'='#fd8d3c', 'p<0.1'='#fed976')
  imbalance_breaks <- seuid_traj$imbalance_breaks
  resdf <- seuid_traj$resdf
  break_prop <- seuid_traj$break_prop
  break_sample_prop <- seuid_traj$break_sample_prop
  break_sample_prop <- lapply(split(break_sample_prop, f=break_sample_prop$celltype), function(i){
    loessMod <- loess(proportion ~ pseudotime, data=i, span=0.25) # 25% smoothing span
    i$smoothed.proportion <- predict(loessMod) 
    return(i)
  }) %>% do.call(rbind, .)

  ## Visualize
  gg_imbalance <- ggplot(imbalance_breaks, aes(x=pseudotime)) +
    geom_ribbon(aes(ymin=conf.int1, ymax=conf.int2), linetype=2, alpha=0.2) +
    geom_line(aes(y=delta, group=1, color=p)) +
    facet_grid(batchid~partition, space='free', scale='free_x') +
    scale_color_manual(values=sig.cols) +
    cowplot::theme_cowplot() +
    ylim(-6,6) +
    geom_hline(yintercept=0, linetype=2) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2,"line"),
          legend.box="vertical", legend.margin=margin(), 
          axis.text.x=element_blank(),
          axis.title.x = element_blank()) +
    guides(fill=guide_legend(nrow=1))
  gg_sampleprop <- ggplot(break_sample_prop %>%
                            mutate(sig = abs(proportion) > 2), 
                          aes(x=pseudotime, y=smoothed.proportion, group=celltype, color=sig)) +
    geom_line() +
    facet_grid(celltype~partition, space='free', scale='free_x') +
    cowplot::theme_cowplot() +
    geom_hline(yintercept=0, linetype=2) +
    scale_color_manual(values=c('TRUE'='red', 'FALSE'='grey')) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2,"line"),
          legend.box="vertical", legend.margin=margin(), 
          axis.text.x=element_blank(),
          axis.title.x = element_blank()) +
    guides(fill=guide_legend(nrow=1))
  gg <- ggplot(resdf, aes(x = pseudotime)) +
    geom_density(data=resdf, 
                 aes(y=..scaled.., fill=conditions),
                 alpha = .70) + 
    # geom_label(data=prog_res, aes(label=label), x=Inf, y=Inf, hjust=1, vjust=1) +
    scale_fill_brewer(type = "qual") +
    facet_grid(comp~partition, space='free', scale='free_x') +
    cowplot::theme_cowplot() +
    ylim(0, 1) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2,"line"),
          legend.box="vertical", legend.margin=margin()) +
    guides(fill=guide_legend(nrow=1))
  
  ggbar <-  ggplot(break_prop, aes(x = pseudotime, y=proportion, fill=celltype)) +
    geom_bar(position='stack', stat='identity') +
    # ylim(0,1) +
    cowplot::theme_cowplot() +
    facet_grid(1~partition, space='free', scale='free_x') +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2,"line"),
          legend.box="vertical", legend.margin=margin(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    guides(fill=guide_legend(ncol=4, byrow=T))
  
  list("ggimbalance_group"=plot_grid(plotlist = list(gg_imbalance, ggbar), ncol=1, rel_heights = c(1.5,1), align='v'),
       "ggimbalance_individual"=plot_grid(plotlist = list(gg, ggbar), ncol=1, rel_heights = c(2,1), align='v') )
})
pdf("~/xfer/m3_imbalance_score.compare.pdf", height = 5, width = 6)
imbalance_comp_gg
dev.off()
  
imbalance_umap <- lapply(names(seu_traj_imbalance), function(cdsid){
  print(cdsid)
  cds <- seu_traj_imbalance[[cdsid]]$cds
  b <- c(-6, 0, 6)
  gg <- plot_cells(cds$cds, color_cells_by = "treg_anno", show_trajectory_graph=T,
                   label_cell_groups=T,
                   label_principal_points=F) +
    ggtitle(cdsid)
  gg2 <- plot_cells(cds$cds, color_cells_by = "imbalance_score", show_trajectory_graph=T,
                    label_cell_groups=T,
                    label_principal_points=F) +
    scale_color_gradientn(limits = c(min(b),max(b)),
                          colours=c("blue", "grey", "red"),
                          breaks=b, labels=format(b)) #+
    #geom_label(data=prog_all, aes(x=x, label=label), y=Inf, hjust=1, vjust=1)
  gg2$layers[[2]]$aes_params$size <- NULL
  gg2$layers[[2]]$aes_params$stroke <- 0
  gg2$layers[[2]]$aes_params$alpha <- 0.5
  gg2$layers[[2]]$mapping$size <- quo(abs(imbalance_score_size))
  
  colData(cds$cds) <- cbind(colData(cds$cds), (seu_traj_imbalance[[cdsid]]$cds$pseudotime))
  gg <- lapply(grep("all", colnames(seu_traj_imbalance[[cdsid]]$cds$pseudotime), value=T), function(pseudox){
    plot_cells(cds$cds, color_cells_by = pseudox, show_trajectory_graph=T,
               label_cell_groups=T,
               label_principal_points=F)
  }) 
  
  title <- ggdraw() + 
    draw_label(cdsid, fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  mainplot <- cowplot::plot_grid(plotlist=c(list(gg2), gg), 
                                 labels=c('Imbalance_score', 'Pseudotime', 'Pseudotime'),
                                 nrow = 1)
  return(plot_grid(title, mainplot,ncol = 1, rel_heights = c(0.1, 1)))
})
pdf("~/xfer/m3_imbalance_score.umap.pdf", width = 8)
imbalance_umap
dev.off()

pdf("~/xfer/m3_selected_paths.pdf")
lapply(names(seu_traj_imbalance), function(cdsid){
  cds <- seu_traj_imbalance[[cdsid]]$cds
  apply(cds$cellWeights, 2, function(i){
    plot_cells(cds$cds[,names(which(i==1))], 
               color_cells_by = "treg_anno", show_trajectory_graph=T,
               label_cell_groups=T,
               label_principal_points=F)
  })
})
lapply(names(seu_traj_imbalance), function(cdsid){
  cds <- seu_traj_imbalance[[cdsid]]$cds
  lapply(cds$cds_all, function(cdsi){
    plot_cells(cdsi, 
               color_cells_by = "pseudotime", show_trajectory_graph=T,
               label_cell_groups=T,
               label_principal_points=F)
  })
})
dev.off()



## Run GAM for Differential expression along trajectory
palette <- hcl.colors(30, palette = "blue-red")
select.genes <- T # Select what genes to plot

if(mode == 'treg'){
  ## Load in regulons
  outdirregulon <- file.path(PDIR, "results", "regulons", "pyscenic")
  looms <- lapply(c('LN', 'Tumor'), function(grp){
    open_loom(file.path(outdirregulon, paste0(grp, '.seu_scenic.loom')))
  }) %>% setNames(c('LN', 'Tumor'))
  regulonAUC <- lapply(looms, function(loom){
    assay(get_regulons_AUC(loom, column.attr.name='RegulonsAUC')) %>%
      as.matrix
  })
  
  ## Load in pathways
  outdir_pway <- file.path(PDIR, "results", "pathway")
  scores_l <- readRDS(file=file.path(PDIR, "results", "pathway", 'aucell', 
                                     paste0("aucell.", 'msigdb', ".", 'all', ".rds")))
  scores <- lapply(scores_l, assay) %>% 
    do.call(rbind, .)
} else {
  regulonAUC <- NULL
  scores <- NULL
}





method = 'rna' # 'rna', 'regulon', 'pathway'
treg_gg_gams <- lapply(names(seu_traj_imbalance), function(seuid){
  ### Generate trajectory analysis
  traj_i <- seu_traj_imbalance[[seuid]]
  if(mode=='treg') traj_i$cds$pseudotime <- traj_i$cds$pseudotime[,grep("all$", colnames(traj_i$cds$pseudotime))]
  break_prop <- traj_i$break_prop
  
  gams <- lapply((names(comps)), function(compid){
    message(compid)
    comp <- comps[[compid]]
    batchsamples <- unique(grep(paste0("^", comp[1]), traj_i$cds$conditions, value=T))
    conditions <- unique(grep(paste0("^", comp[1], ".*", comp[2]), traj_i$cds$conditions, value=T))
    
    ## Test each condition in each trajectory
    partition_gams <- lapply(colnames(traj_i$cds$pseudotime), function(pt_i){
      gammodelf <- file.path(outdir, "gamModels", paste0("gam.", seuid, "_", method, ".", compid, ".", pt_i, ".rds"))
      if(file.exists(gammodelf) & !overwrite){
        message(paste0("Reading in gam model for ", compid, " - ", pt_i))
        gamModels <- readRDS(file=gammodelf)
      } else {
        message("Generating GAM")
        non.inf.idx <- which(traj_i$cds$cellWeights[,pt_i]!=0)
        n_genes <- 3000
        qval_thresh <- 0.1
        
        message(paste0("Selecting for batch samples: ", paste(batchsamples, collapse=", ")))
        batch.idx <- unlist(sapply(batchsamples, grep, x=rownames(traj_i$cds$cellWeights)))
        comp.idx <- unlist(sapply(conditions, grep, x=rownames(traj_i$cds$cellWeights)))
        batch.noninf.idx <- intersect(non.inf.idx, batch.idx)
        comp.noninf.idx <- intersect(non.inf.idx, comp.idx)
        
        message(paste0("Running monocle3 for: ", seuid, " - ", pt_i, "..."))
        rowcnts <- rowSums(assay(traj_i$cds$cds))
        genes <- rownames(traj_i$cds$cds)[which(rowcnts > ncol(traj_i$cds$cds)*0.05)]
        
        monocle3_graphf <- file.path(outdir, "gamModels", paste0("m3.", seuid, "_", method, ".", comp[1], ".", pt_i, ".rds"))
        cds <- switch(method,
                      'rna'=traj_i$cds$cds[genes,],
                      'regulon'=.clone_cds_params(new.assay=regulonAUC[[seuid]],
                                                  ref.cds=traj_i$cds$cds[genes,]),
                      'pathway'=.clone_cds_params(new.assay=scores,
                                                  ref.cds=traj_i$cds$cds[genes,]))
        if(!file.exists(monocle3_graphf) | overwrite){
          cds_pr_test_res <- graph_test(cds[,batch.noninf.idx], 
                                        neighbor_graph="principal_graph", cores=4,
                                        expression_family=if(method=='rna') 'quasipoisson' else 'binomial')
          saveRDS(cds_pr_test_res, file=monocle3_graphf)
        } else {
          cds_pr_test_res <- readRDS(file=monocle3_graphf)
        }
        pr_deg_ids <- cds_pr_test_res %>% 
          arrange(q_value, desc(morans_I)) %>% head(., n_genes)
        if(!any(pr_deg_ids$q_value < qval_thresh)){
          pr_deg_ids <- pr_deg_ids %>%
            head(., 50) %>%
            pull(gene_short_name)
        } else {
          pr_deg_ids <- pr_deg_ids %>%
            dplyr::filter(q_value < qval_thresh) %>%
            pull(gene_short_name)
        }
        
        
        message(paste0("Running tradeSEQ using monocle3 trajectory genes for: ", seuid, " - ", pt_i, "..."))
        gamModels <- fitGAM(counts = as.matrix(assays(cds[pr_deg_ids,comp.noninf.idx])$counts),
                            pseudotime = traj_i$cds$pseudotime[comp.noninf.idx,pt_i,drop=F],
                            cellWeights = traj_i$cds$cellWeights[comp.noninf.idx,pt_i,drop=F],
                            conditions = factor(traj_i$cds$conditions[comp.noninf.idx]),
                            nknots = 3,
                            family=if(method=='rna') 'nb' else 'betar')  # 'betar', 'gaussian'
        # eps <- 1e-10
        saveRDS(gamModels, file=gammodelf)
      }
      return(gamModels)
    })
    names(partition_gams) <- colnames(traj_i$cds$pseudotime)
    
    return(partition_gams)
  })
  names(gams) <- names(comps)
  
  
  # Visualize the GAM differential expressions
  ngenes2report <- 100
  batches <- unique(gsub("_.*", "", names(comps)))
  gamdf <- lapply(batches, function(batchid){
    comp_ids <- grep(batchid, names(comps), value=T)
    
    lapply(names(trajs[[seuid]]), function(partitionid){
      message(paste0("Marker 1: ", seuid, ", ", batchid, " - ", partitionid))
      condresf <- file.path(outdir, "gamModels", "condres", paste0("gam.", seuid, "_", method, ".", batchid, ".", partitionid, ".rds"))
      dir.create(dirname(condresf), recursive = T, showWarnings = F)
      if(!file.exists(condresf) | overwrite){
        ## Assess differential expression patterns between batch-controlled treated and untreated (e.g. 3d vs Un)
        condres_genes <- lapply(comp_ids, function(compid){
          message(paste0("  >> ", compid))
          condRes_monocle <- conditionTest(gams[[compid]][[partitionid]], l2fc = log2(2))
          condRes <- condRes_monocle %>%
            mutate(padj = p.adjust(pvalue, "fdr"),
                   gene = rownames(.))
          
          waldstat_ordered_genes <- condRes %>% 
            dplyr::filter(!is.na(waldStat)) %>%
            arrange(padj, desc(waldStat)) %>%
            pull(gene)
          conditionGenes <- intersect(waldstat_ordered_genes, rownames(condRes)[condRes$pvalue <= 0.1])
          conditionGenes <- conditionGenes[!is.na(conditionGenes)]
          
          return(list("genes"=conditionGenes,
                      "condRes"=condRes))
        })
        names(condres_genes) <- comp_ids
        saveRDS(condres_genes, file=condresf)
      } else {
        condres_genes <- readRDS(file=condresf)
      }
      conditionGenes <- sapply(condres_genes, function(i) i$genes) %>%
        unlist %>% as.character %>% unique
      
      # Aggregate the KO-WT genes and scale across both differential conditions
      condres_res <- lapply(names(condres_genes), function(compid){
        if(length(conditionGenes)==0) return(NULL)
        condRes <- condres_genes[[compid]]$condRes
        # seu_traj_imbalance$LN$cds$cds
        pseudotime_breaks <- traj_i$imbalance_breaks %>% 
          dplyr::filter(partition == partitionid,
                        batchid == gsub("_.*", "", compid)) %>% 
          dplyr::select(pseudotime_breaks, pseudotime)
        yhatSmooth <- tradeseq.predictSmooth_conditions(models=gams[[compid]][[partitionid]], 
                                                        gene = conditionGenes, pseudotime_breaks=pseudotime_breaks, tidy=F)
        yhatSmooth <- log1p(yhatSmooth)
        yhatSmooth[is.infinite(yhatSmooth)] <- max(yhatSmooth[!is.infinite(yhatSmooth)])
        yhatSmoothScaled <- yhatSmooth
        # yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
        
        id <- comps[[compid]]
        ids <- grep(paste0(id[1], ".*?", id[2]), colnames(yhatSmoothScaled), value=T) 
        uids <- gsub("^.*condition", "", ids) %>% gsub("_point[0-9]*$", "", .) %>% unique
        
        # gene_ord <- intersect(rownames(yhatSmoothScaled),waldstat_ordered_genes[seq_len(ngenes2report)])
        # yhatSmoothScaled <- yhatSmoothScaled[gene_ord,]
        # alldat <- yhatSmoothScaled[,as.character(sapply(uids, grep, x=colnames(yhatSmoothScaled), value=T))]
        # hc <- hclust(dist(alldat))
        
        
        heatmaps <- lapply(uids, function(uid_i){
          dat <- yhatSmoothScaled[, grep(uid_i, colnames(yhatSmoothScaled), value=T)]
          dat_long <- as.data.frame(dat) %>% tibble::rownames_to_column("gene") %>% 
            tidyr::pivot_longer(., cols = !gene, names_to = 'timepoint', 
                                values_to = 'yHatSmoothScaled') %>% 
            mutate(timepoint = gsub("^.*point", "", timepoint),
                   pseudotime = pseudotime_breaks[as.integer(timepoint),'pseudotime'],
                   pseudotime.st = pseudotime - median(diff(pseudotime)),
                   ymin=as.numeric(as.factor(gene))-0.5,
                   ymax=as.numeric(as.factor(gene))+0.5,
                   uid=uid_i,
                   partition=partitionid,
                   comp=compid)
          return(dat_long)
          
          # ComplexHeatmap::Heatmap(dat,
          #                         cluster_rows=F, cluster_columns = FALSE, column_title = uid_i, 
          #                         show_row_names = T, show_column_names = F) %>% 
          #   draw %>% grid.grabExpr
        }) %>% do.call(rbind, .)
        return(heatmaps)
      })
      do.call(rbind, condres_res)
    }) %>%
      do.call(rbind, .)
  }) %>% do.call(rbind, .)
  gamdf$batch <- gsub("_.*", "", gamdf$comp)
  
  ## Selects for genes to plot and the order
  if(select.genes){
    genes <- read.table(file.path(PDIR, "ref", "trajectory_genesets", paste0(mode, ".", seuid, ".txt")), 
                        header=F, stringsAsFactors = F) %>%
      mutate(V1 = str_to_title(V1))
    genes <- genes[!duplicated(genes),]
    gamdf <- left_join(genes, gamdf, by=c('V1'='gene',
                                          'V2'='partition',
                                          'V3'='batch'), 
                       relationship = "many-to-many")
    colnames(gamdf)[1:3] <- c('gene', 'partition', 'batch')
  }
  
  
  ## Calculate difference between yhat for Treated and Untread, and between compids (KO vs WT)
  # And scale across KO and WT condition
  gamdf_spl <- split(gamdf, f=list(gamdf$partition, gamdf$comp))
  gams_delta_df <- lapply(names(trajs[[seuid]]), function(partition_j){
    gams_comps <- lapply(batches, function(batchid){
      message(paste0("Marker 2: ", batchid,  " - ", partition_j))
      comp_ids <- grep(batchid, names(comps), value=T)
      gamdf_wides <- lapply(comp_ids, function(comp_i){
        message(paste0("  >> ", comp_i))
        listidx <- paste(partition_j, comp_i, sep=".")
        gamdf_wide <- tidyr::pivot_wider(gamdf_spl[[listidx]], names_from = uid, values_from=yHatSmoothScaled)
        id = comps[[comp_i]]
        colids <- grep(paste0(id[1], ".*", id[2]), colnames(gamdf_wide), value=T)
        colids <- colids[order(grepl("Un", colids))]
        
        if(nrow(gamdf_wide) == 0){
          gamdf_wide <- NULL 
        }else {
          gamdf_wide <- gamdf_wide %>%
            mutate(delta_treated = as.numeric(unlist(!!rlang::sym(colids[1]) - !!rlang::sym(colids[2])))) %>%
            dplyr::select(-c(!!rlang::sym(colids[1]), !!rlang::sym(colids[2]), comp))
        }
        return(gamdf_wide)
      })
      names(gamdf_wides) <- comp_ids
      if(any(sapply(gamdf_wides, is.null)) | all(sapply(gamdf_wides, function(i) is.na(i$pseudotime)))) return(NULL)
      
      gamdf_wide <- purrr::reduce(gamdf_wides, full_join, 
                                  by=c('gene', 'timepoint', 'pseudotime', 
                                       'pseudotime.st', 'ymin', 'ymax', 'partition'),
                                  suffix=paste0(".", comp_ids))
      
      scaled_vals <- gamdf_wide %>% 
        dplyr::select(c(gene, timepoint, grep("delta_treated", colnames(.), value=T))) %>%
        tidyr::pivot_wider(., names_from = timepoint, values_from=grep("delta_treated", colnames(.), value=T)) %>%
        tibble::column_to_rownames('gene') %>%
        apply(., 1, scales::rescale) %>% t %>% as.data.frame %>%
        tibble::rownames_to_column(., "gene") %>% 
        tidyr::pivot_longer(., cols = !gene, names_to = 'tempids', 
                            values_to = 'delta_scaled')
      
      scaled_vals <- scaled_vals %>%
        mutate(timepoint=as.character(gsub("^.*?([0-9]*)$", "\\1", tempids)),
               comp=gsub("^.*\\.(.*?)_[0-9]*$", "\\1", tempids),
               batchid=batchid) %>%
        dplyr::select(-tempids)
      scaled_vals <- split(scaled_vals, f=scaled_vals$comp)
      
      lapply(comp_ids, function(i){
        left_join(gamdf_wides[[i]], scaled_vals[[i]], by=c('gene', 'timepoint'))
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
    return(gams_comps)
  }) %>% do.call(rbind, .)
  
  # gg_gams <- lapply(names(comps), function(comp_i){
  gg_gams <- lapply(batches, function(batch_i){
    comp_ids <- grep(batch_i, names(comps), value=T)
    
    
    gg_comp_gams <- lapply(names(trajs[[seuid]]), function(partition_j){
      gamdf_ij <- gams_delta_df %>%
        dplyr::filter(partition == partition_j & 
                        batchid == batch_i) %>%
        mutate(yHatSmoothScaled=delta_scaled)
      if(nrow(gamdf_ij)==0) return(NULL)
      dist_gam <- gamdf_ij %>% 
        dplyr::select(yHatSmoothScaled, timepoint, gene, comp) %>%
        group_split(., comp) %>% 
        lapply(., function(i) {
          pivot_wider(i %>% dplyr::select(-comp), 
                      names_from = timepoint, values_from = yHatSmoothScaled)   %>% 
            tibble::column_to_rownames('gene')
        }) %>%
        do.call(cbind, .) %>% 
        dist
      dist_gam[is.na(dist_gam)] <- max(dist_gam, na.rm=T)
      gene_ord <- rownames(as.matrix(dist_gam))[hclust(dist_gam)$order]
      gene_ord <- rev(unique(gamdf_ij$gene))
      
      coefs_comps <- lapply(comp_ids, function(comp_i){
        print(comp_i)
        gamdf_ij <- gams_delta_df %>%
          dplyr::filter(partition == partition_j &
                          comp == comp_i) %>% # batchid == batch_i) %>%
          mutate(yHatSmoothScaled=delta_scaled)
        coefs <- sapply(split(gamdf_ij, f=gamdf_ij$gene), function(i){
          tryCatch({
            lmx <- lm(data=i, yHatSmoothScaled~as.numeric(pseudotime))
            lmx$coefficients[2]
          }, error=function(e){ NA } )
        })
      }) %>%
        setNames(., comp_ids)
      # coefs <- sort(coefs_comps[[1]] - coefs_comps[[2]])
      # gene_ord <- gsub("\\..*", "", names(coefs))
      
      
      ggdat <- lapply(comp_ids, function(comp_i){
        gamdf_dat <- gams_delta_df %>%
          dplyr::filter(partition == partition_j & 
                          comp == comp_i) %>%
          mutate(yHatSmoothScaled=delta_scaled)
        
        gamdf_dat$geneord <- factor(as.character(gamdf_dat$gene), levels=gene_ord)
        return(gamdf_dat)
        
      }) %>%
        do.call(rbind, .)
      
      
      cowheatmap <- ggplot(ggdat, aes(fill=yHatSmoothScaled)) +
        facet_grid(.~comp, space='free', scales = 'free_y') +
        geom_rect(aes(xmin=pseudotime.st, xmax=pseudotime,
                      ymin=as.factor(geneord), ymax=as.numeric(as.factor(geneord))-1),
                  lineend = "butt") +
        cowplot::theme_cowplot() +
        scale_fill_gradientn(colours = palette, values = seq(0,1, length.out=30)) +
        labs(title=paste0(seuid, " - ", partition_j)) +
        theme(legend.position = "top",
              legend.key.size = unit(0.2,"line"),
              legend.box="vertical", legend.margin=margin(),
              axis.text.x.bottom = element_blank())

      break_prop_j <- break_prop %>% 
        dplyr::filter(partition == partition_j)
      break_prop_j <- rbind(break_prop_j,break_prop_j) %>%
        mutate(comp = sort(rep(comp_ids, nrow(break_prop_j))))
      ggbar <-  ggplot(break_prop_j, aes(x = pseudotime, y=proportion, fill=celltype)) +
        geom_bar(position='stack', stat='identity') +
        # ylim(0,1) +
        cowplot::theme_cowplot() +
        facet_grid(.~comp, space='free') +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.2,"line"),
              legend.box="vertical", legend.margin=margin(),
              strip.background = element_blank(),
              strip.text.x = element_blank()) +
        guides(fill=guide_legend(ncol=4, byrow=T))
      
      pgrid = plot_grid(plotlist = list(cowheatmap, ggbar), ncol=1, 
                        rel_heights = c(4,1), align='v') 
      list("grid"=pgrid,
           "gene_dir"=lapply(coefs_comps, function(coefs) split(gsub("\\..*", "", names(coefs)), coefs > 0)))
      
    })
    names(gg_comp_gams) <- names(trajs[[seuid]])
    gg_comp_gams
  })
  unlist(gg_gams, recursive=F)
})
names(treg_gg_gams) <- names(seu_traj_imbalance)
pdf("~/xfer/m3_differential_expression.pdf", height = 10, width=8)
unlist(gg_gams, recursive=F)
# treg_gg_gams
dev.off()


## Run ORA for Differential expressed genes along trajectory based on direction
ora_res <- lapply(names(treg_gg_gams), function(seuid) {
  lapply(names(treg_gg_gams[[seuid]]), function(part_comp_i) {
    message(paste0(seuid, " - ", part_comp_i))
    i <- treg_gg_gams[[seuid]][[part_comp_i]]
    names(i$gene_dir) <- c('TRUE'='up', 'FALSE'='down')[names(i$gene_dir)]
    lapply(names(i$gene_dir), function(direction){
      res <- clusterProfiler::enricher(
        gene = i$gene_dir[[direction]], # A vector of your genes of interest
        pAdjustMethod = "BH", # Method to be used for multiple testing correction
        universe = Features(seu), # A vector containing your background set genes
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.05,
        TERM2GENE = geneset_df
      )@result %>%
        mutate(seuid = seuid,
               comp_partition=part_comp_i,
               direction=direction) %>%
        dplyr::select(-c(Description)) %>%
        tibble::remove_rownames()
      return(res)
    })  %>% 
      do.call(rbind, .)
  }) %>% 
    do.call(rbind, .)
}) %>% 
  do.call(rbind, .)

write.table(ora_res, file=file.path("~/xfer", "m3_differential_expression.ora.csv"),
            sep=",", row.names = F, col.names = T, quote = F)





lapply(names(seu_velo_l), function(seu_velo_id){
  ### Generate trajectory analysis
  seu <- seus[[seu_velo_id]]
  Idents(seu) <- 'orig.ident'
  comps <- list("B1_WT"=c("B1", "WT"),
                "B1_KO"=c("B1", "KO"),
                "B2_WT"=c("B2", "WT"),
                "B2_KO"=c("B2", "KO"))
  
  # set.seed(1234)
  # cds <- runMonocle3(seu, reduction='umap.mnn')
  
  comp_cds_l <- lapply(names(comps), function(compid){
    message(compid)
    comp <- comps[[compid]]
    seusub <- subset(seu, ident=unique(grep(paste0("^", comp[1], ".*", comp[2]), Idents(seu), value=T)))
    # seusub <- .splitAndIntegrate(seusub)
    
    #>>> Monocle3
    set.seed(1234)
    cds <- runMonocle3(seusub)
    message("Adding nodes")
    cds <- runMonocle3(cds, nodes_l=trajs[[seu_velo_id]][[compid]])
    idx <- grep(paste0("^", comp[1], ".*", comp[2]), cds$conditions)
    
    message("Running imbalance")
    scores <- condiments::imbalance_score(Object = reducedDims(cds$cds[,idx])$UMAP,
                                          conditions =  colData(cds$cds[,idx])$orig.ident,
                                          k = 20, smooth = 40)
    message("Adding imbalance")
    df <- as.data.frame(reducedDims(cds$cds[,idx])$UMAP) %>%
      mutate("score"=as.numeric(round(scores$scores,3)),
             "scaled_scores"=as.numeric(round(scores$scaled_scores,3)),
             'scaled_scores2' = scales::rescale(scores$scaled_scores, to=c(0,1))) 
    cds$imbalance <- df
    
    cds$cds$imbalance_score <- 0
    cds$cds$group <- NA
    cds$cds[,idx]$imbalance_score <- df$scaled_scores
    cds$cds[,idx]$group <- compid
    return(cds)
  })
  saveRDS(comp_cds_l, file=file.path(PDIR, "results", "trajectory", "monocle3", paste0("ST2_TReg.condiments.", seu_velo_id, ".rds")))
 
   ## Plot imbalances 
  pdf("~/xfer/m3_imbalance_score.test.pdf", width = 10); 
  lapply(names(comp_cds_l), function(cdsid){
    cds <- comp_cds_l[[cdsid]]
    b <- c(-6, 0, 6)
    gg <- plot_cells(cds$cds, color_cells_by = "treg_anno", show_trajectory_graph=T,
                     label_cell_groups=T,
                     label_principal_points=F) +
      ggtitle(cdsid)
    gg2 <- plot_cells(cds$cds, color_cells_by = "imbalance_score", show_trajectory_graph=T,
                      label_cell_groups=T,
                      label_principal_points=F) +
      scale_color_gradientn(limits = c(min(b),max(b)),
                           colours=c("blue", "grey", "red"),
                           breaks=b, labels=format(b)) 
    gg2$layers[[2]]$aes_params$size <- NULL
    gg2$layers[[2]]$mapping$size <- quo(abs(imbalance_score))

    
    pdf("~/xfer/m3_imbalance_score.test.pdf", width = 10)
    gg + gg2
    dev.off()
  })
  dev.off()
  
  
  ## Plot progression and run tests
  all_prog_dat <- lapply(comp_cds_l, function(cds){
    prog_dat <- plotProgressionAndCelltype(cds)
    return(prog_dat)
  })
  pdf("~/xfer/m3_progression_score.tumor.pdf")
  lapply(names(all_prog_dat),function(id) all_prog_dat[[id]]$gg + ggtitle(id))
  dev.off()
  
})

## Run the Condiments workflow to compare between conditions
balance.groups <- F
lapply(names(seu_velo_l), function(seu_velo_id){
  ### Generate trajectory analysis
  seu <- seus[[seu_velo_id]]
  Idents(seu) <- 'orig.ident'
  comps <- list("B1_WT"=c("B1", "WT"),
                "B1_KO"=c("B1", "KO"),
                "B2_WT"=c("B2", "WT"),
                "B2_KO"=c("B2", "KO"))
  ps <- lapply(names(comps), function(compid){
    comp <- comps[[compid]]
    seusub <- subset(seu, ident=unique(grep(paste0("^", comp[1], ".*", comp[2]), Idents(seu), value=T)))
    seusub <- .splitAndIntegrate(seusub) %>%
      .monocle3Cluster(.)
    
    # Remove Cells where there are too few celltypes per group in a partition
    cnttbl <- table(seusub$treg_anno, seusub$orig.ident, seusub$monocle3_partitions)
    partition_keepids <- lapply(c(1:dim(cnttbl)[3]), function(tblidx){
      # rmidx <- which(cnttbl[,,tblidx] <= (sum(cnttbl[,,tblidx])*0.01), arr.ind=T)
      rmidx <- which(cnttbl[,,tblidx] <= 3, arr.ind=T)
      setdiff(rownames(cnttbl[,,tblidx]), unique(rownames(cnttbl[,,tblidx])[rmidx[,1]]))
    })
    names(partition_keepids) <- dimnames(cnttbl)[[3]]
    
    print(paste0(">>> ", compid))
    cells <- sapply(names(partition_keepids), function(partitionid){
      sel <- (seusub$monocle3_partitions == as.character(partitionid) & 
                seusub$treg_anno %in% partition_keepids[[as.character(partitionid)]])
      as.character(Cells(seusub)[sel])
    }) %>% unlist %>% as.character
    message(paste0(length(cells), "/", length(Cells(seusub)), " cells..."))
    seusub <- subset(seusub, cells=cells)
    print(paste0("<<< ", compid))
      
    ## Slingshot
    Idents(seusub) <- 'treg_anno'
    seusub@reductions <- seusub@reductions[c('umap')]
    partitions <- as.character(unique(seusub$monocle3_partitions))
    partition_dat <- lapply(setNames(partitions, partitions), function(partition){
      print(paste0(compid, " - ", partition))
      Idents(seusub) <- 'monocle3_partitions'
      sce <- as.SingleCellExperiment(subset(seusub, ident=as.character(partition)),
                                     assay='RNA')
      
      ids <- names(which(table(colData(sce)$treg_anno) > 5))
      sce <- sce[,which(sce$treg_anno %in% ids)]
      sce <- slingshot(sce, reducedDim = 'UMAP', 
                       clusterLabels = colData(sce)$treg_anno,
                       approx_points = ceiling(0.1* ncol(sce)))
      curves <- slingCurves(sce, as.df = TRUE)
      curve_gg <- geom_path(data = curves %>% arrange(Order),
                            aes(group = Lineage)) 
      
      ## Imbalance score
      scores <- condiments::imbalance_score(Object = reducedDims(sce)$UMAP,
                                            conditions =  colData(sce)$orig.ident,
                                            k = 20, smooth = 40)
      
      df <- as.data.frame(reducedDims(sce)$UMAP) %>%
        mutate("score"=as.numeric(round(scores$scores,3)),
               "scaled_scores"=as.numeric(round(scores$scaled_scores,3)))
      p <- ggplot(df, aes(x = umap_1, y = umap_2, col = scaled_scores)) +
        geom_point() +
        scale_color_viridis_c(option = "C") + 
        ylim(-12,12) + xlim(-12,12) +
        ggtitle(paste0("Treg ST2: ", seu_velo_id, 
                       ", compid: ", compid, 
                       ", Partition: ", partition)) +
        geom_path(data = curves %>% arrange(Order),
                  aes(group = Lineage), col='black')
        
      
      return(list("sce"=sce, "curves"=curves, "df"=df))
    })
    
    df <- lapply(partition_dat, function(dat){
      dat$df %>% mutate(scaled_scores = scales::rescale(scaled_scores, to=c(0,1)))
    }) %>% do.call(rbind, .)
    curves <- lapply(names(partition_dat), function(pid){
      partition_dat[[pid]]$curves %>% mutate(Lineage=pid)
    }) %>% do.call(rbind, .)
    
    pdf("~/xfer/test.pdf"); 
    gg <- ggplot(df, aes(x = umap_1, y = umap_2, col = scaled_scores)) +
      geom_point() +
      scale_color_viridis_c(option = "C") + 
      ylim(-12,12) + xlim(-12,12) +
      ggtitle(paste0("Treg ST2: ", seu_velo_id, 
                     ", compid: ", compid)) +
      geom_path(data = curves %>% arrange(Order),
                aes(group = Lineage), col='black')
    dev.off()
    
    return(list("sce"=lapply(partition_dat, function(i)i$sce), 
                "seu"=seusub,
                "gg"=gg, "curves"=curves, "df"=df))
  })
  names(ps) <- names(comps)
  
  pdf("~/xfer/x22.pdf")
  lapply(names(ps), function(id) DimPlot(ps[[id]]$seu, group.by='treg_anno', split.by='orig.ident') + ggtitle(id))
  lapply(names(ps), function(id) DimPlot(ps[[id]]$seu, group.by='treg_anno') + ggtitle(id))
  dev.off()
  
  pdf("~/xfer/y22.pdf")
  lapply(ps, function(i) i$gg)
  dev.off()
  
  toptests <- lapply(ps, function(grp_i){
    lapply(grp_i$sce, function(sce_j){
      topologyTest(sds = sce_j, 
                   conditions = colData(sce_j)$orig.ident,
                   parallel=F)
    })
  })
  
  
  #--- Differential Progression ----
  progtests  <- lapply(ps, function(grp_i){
    res_dat <- lapply(names(grp_i$sce), function(partitionid){
      sce_j <- grp_i$sce[[partitionid]]
      psts <- slingPseudotime(sce_j) %>%
        as.data.frame() %>%
        mutate(cells = rownames(.),
               conditions = colData(sce_j)$orig.ident) %>%
        pivot_longer(starts_with("Lineage"), 
                     values_to = "pseudotime", 
                     names_to = "lineages") %>%
        mutate("partition" = partitionid,
               'anno'=colData(sce_j[,cells])$treg_anno)
      prog_res <- progressionTest(sce_j, conditions = colData(sce_j)$orig.ident, 
                                  global = TRUE, lineages = TRUE) %>%
        mutate("partition"=partitionid)
      return(list("res"=prog_res, "df"=psts))
    })
    
    resdf <- lapply(res_dat, function(i) i$df) %>% 
      do.call(rbind, .) %>% as.data.frame %>%
      mutate(pseudotime_breaks = cut(pseudotime, breaks=100),
             anno=factor(anno))
    break_prop <- lapply(split(resdf, f=resdf$partition), function(partition_i){
      lapply(split(partition_i, f=partition_i$lineages), function(lineage_ij){
        lapply(split(lineage_ij, lineage_ij$pseudotime_breaks), function(pseudotime_x){
          as.data.frame(t(as.matrix(tail(sort(table(pseudotime_x$anno)/nrow(pseudotime_x)), n=3))))
        }) %>%
          rbind.fill(.) %>%
          mutate(lineages=unique(lineage_ij$lineages),
                 partition=unique(lineage_ij$partition),
                 pseudotime_breaks= names(split(lineage_ij, lineage_ij$pseudotime_breaks)),
                 pseudotime = .break2number(pseudotime_breaks, ret='median')) %>%
          as.data.frame
      }) %>% rbind.fill(.)
    })  %>% rbind.fill(.)
    break_prop <- pivot_longer(break_prop, 
                               cols=!c(lineages, partition, pseudotime_breaks, pseudotime),
                               names_to='celltype',
                               values_to = 'proportion') %>%
      mutate(proportion = scales::rescale(proportion, to = c(0, -0.25)))
    
    
    gg <- ggplot(resdf, aes(x = pseudotime, y=..scaled.., fill = conditions)) +
        geom_density(alpha = .70) + 
      geom_bar(data=break_prop, 
               aes(x=pseudotime, y=proportion, fill=celltype),
               position='stack', stat='identity') +
        scale_fill_brewer(type = "qual") +
        facet_grid(partition~lineages, space='free') +
        cowplot::theme_cowplot() +
        ylim(-0.25, 1) +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.2,"line"),
              legend.box="vertical", legend.margin=margin()) +
        guides(fill=guide_legend(nrow=4))
    res <- lapply(res_dat, function(i) i$res) %>% do.call(rbind, .)
    return(list("gg"=gg, "res"=res))
  })
   
  lapply(progtests, function(i) i$res)
  
  pdf("~/xfer/z2.pdf")
  lapply(progtests, function(i) i$gg)
  dev.off()
  
  
  #--- Differential expression ----
  library(tradeSeq)
  progtests  <- lapply(ps, function(grp_i){
    res_dat <- lapply(names(grp_i$sce), function(partitionid){
      sce_j <- grp_i$sce[[partitionid]]
      # get number of knots: 
      icMat <- evaluateK(counts = as.matrix(assays(sce_j)$counts),
                         pseudotime = pathStats(colData(sce_j)$slingshot)$pseudotime,
                         cellWeights = pathStats(colData(sce_j)$slingshot)$weights,
                         conditions = factor(colData(sce_j)$orig.idents),
                         nGenes = 300,
                         k = 3:7,
                         parallel = TRUE,
                         BPPARAM = BiocParallel::SerialParam())
      
      set.seed(3)
      sce_j <- fitGAM(counts = sce_j, nknots = 3,
                      conditions = factor(colData(sce_j)$orig.ident))
      mean(rowData(tgfb)$tradeSeq$converged)
      
      cond_genes <- conditionTest(sds)
      cond_genes$padj <- p.adjust(cond_genes$pvalue, "fdr")
      
      
    })
  })
  
  
  
})

## Run Monocle3 to create the trajectory datat
lapply(names(seu_velo_l), function(seu_velo_id){
  ### Generate trajectory analysis
  data <- as(as.matrix(seusub@assays$RNA$counts), 'sparseMatrix')
  pd <- seusub@meta.data
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  
  cds <- new_cell_data_set(expression_data=data,
                           cell_metadata  = seusub@meta.data,
                           gene_metadata = fData)
  
  ## Step 1: Normalize and pre-process the data
  cds <- preprocess_cds(cds, num_dim = 50, norm_method='log')
  ## Step 2: Remove batch effects with cell alignment
  cds <- align_cds(cds, alignment_group = "orig.ident")
  ## Step 3: Reduce the dimensions using UMAP
  reducedDims(cds)$UMAP <- Embeddings(seusub, 'umap')
  reducedDims(cds)$PCA <- Embeddings(seusub, 'pca')[,1:50]
  # cds <- reduce_dimension(cds)
  ## Step 4: Cluster the cells
  cds <- cluster_cells(cds)
  ## Step 5: Learn a graph
  cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE,
                     list(ncenter=275))
  pdf("~/xfer/x3.pdf")
  plot_cells(cds, color_cells_by = "treg_anno", show_trajectory_graph=T,
             label_cell_groups=F,
             label_principal_points=T)
  dev.off()
  
  ## Step 6: Order cells
  cds <- order_cells(cds, reduction_method='UMAP', root_pr_nodes='Y_93')
  cds1 <- choose_graph_segments(cds,
                                reduction_method = "UMAP",
                                starting_pr_node = 'Y_93',
                                ending_pr_nodes = 'Y_101',
                                clear_cds=FALSE)
  pseudotime(cds1, reduction_method='UMAP')
  
  
  
  
  seu_small$monocle3_clusters <- cds@clusters$UMAP$clusters
  seu_small$monocle3_partitions <- cds@clusters$UMAP$partitions
  cds_tumor <- choose_graph_segments(cds,
                                     reduction_method = "UMAP",
                                     starting_pr_node = 'Y_93',
                                     ending_pr_nodes = 'Y_101',
                                     clear_cds=FALSE)
  cds_ln <- choose_graph_segments(cds,
                                  reduction_method = "UMAP",
                                  starting_pr_node = 'Y_139',
                                  ending_pr_nodes = 'Y_113',
                                  clear_cds=FALSE)
  
  res <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
  subset_res <- graph_test(cds_tumor, neighbor_graph="principal_graph", cores=1)
  pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
  
  pdf(file.path(outdir, "trajectory", paste0("monocle3_traj.", seu_velo_id, ".pdf")), width=10);
  dp_mc <- DimPlot(seu_small, label = TRUE, repel = TRUE, reduction = "umap.mnn",
                   group.by="monocle3_partitions", pt.size=0.5, shuffle=TRUE, raster=T)
  dp_mp <- DimPlot(seu_small, label = TRUE, repel = TRUE, reduction = "umap.mnn",
                   group.by="monocle3_clusters", pt.size=0.5, shuffle=TRUE, raster=T)
  cowplot::plot_grid(dp_mc, dp_mp, ncol=2)
  
  publicationDimPlot(seu_small, grp='treg_anno', pt.size=1, reduction='umap.mnn')
  
  plot_cells(cds, color_cells_by='treg_anno',show_trajectory_graph=T,  label_principal_points = TRUE)
  plot_cells(cds, color_cells_by = "pseudotime",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5)
  plot_cells(cds_tumor, color_cells_by = "treg_anno", show_trajectory_graph=T)
  plot_cells(cds_ln, color_cells_by = "treg_anno",show_trajectory_graph=T)
  dev.off()
  file.copy(file.path(outdir, "trajectory", paste0("monocle3_traj.", seu_velo_id, ".pdf")),
            to="~/xfer", overwrite = T)
  cat(paste0("xfer monocle3_traj.", seu_velo_id, ".pdf\n"))
})


#--- a) TRegs or CD8 only ----
if(!exists("seul")) seul <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln2.rds"))
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
# if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
seul <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))

cols <- setNames(c('#1b9e77','#d95f02'), c('Treg_1','Treg_2'))
seul <- seul[c('LN', 'Tumor')]
seul$LN$treg_anno2[seul$LN$treg_anno2 == 'cTreg'] <- 'Treg-LT_Stat1'
seul <- lapply(seul, function(seu){ 
  seu$treg_anno <- paste0(seu$Tissue, ".", seu$treg_anno2)
  seu[['umap']] <- seu[['umap.mnn']]
  seu <- .monocle3Cluster(seu)
  Idents(seu) <- 'monocle3_partitions'
  seu <- subset(seu, ident='1')
  return(seu)
})
seul$LN.Tumor <- JoinLayers(merge(seul[[1]], seul[-1]))


prepSubData <- function(seurat_obj, genes = NULL){
  sub_data = seurat_obj
  DefaultAssay(sub_data) <- 'RNA'
  sub_data <- NormalizeData(object = sub_data, 
                            scale.factor = 10000, display.progress = F) %>% 
    Seurat::FindVariableFeatures(., display.progress = F, 
                                 num.bin = 100,
                                 binning.method = "equal_frequency")
  
  data_expr <- sub_data@assays$RNA@layers$data
  rownames(data_expr) <-  Features(sub_data@assays$RNA)
  colnames(data_expr) <- Cells(sub_data@assays$RNA)
  exp_mat = if(is.null(genes)){
    t(as.matrix(data_expr[VariableFeatures(sub_data@assays$RNA),]))
  } else{
    t(as.matrix(data_expr[genes,]))
  }
  save_mat = merge(sub_data@meta.data[,c("manual_anno", "seurat_clusters")],
                   as.matrix(exp_mat), by = 0)
  rownames(save_mat) = save_mat[,1]
  save_mat = save_mat[,-1]
  
  return(list("seurat" = sub_data, "mat" = save_mat))
}


sampleid <- 'LN.Tumor'
seu <- seul[[sampleid]]
seu$manual_anno <- seu$treg_anno
Idents(seu) <- 'orig.ident'
seu <- prepSubData(seu)
outf <- paste0("./results/bgplvm/", sampleid, ".varfeat.csv")
if(!file.exists(outf)){
  write.csv(seu$mat, file = outf,
          row.names = T, quote = F, col.names = T)
}
seul_celltype <- list(seu) %>% setNames(., sampleid)


#{{ Run GPy: runGpy.py }}


lv <- read.csv(paste0("./results/bgplvm/", sampleid, ".latent_var.csv"), header = T, row.names = 1)
ard <- read.csv(paste0("./results/bgplvm/", sampleid, ".ARD.csv"), header = T, row.names = 1)
lv_i <- 'LV1'

newidmap <- c('Tumor.eTreg_cycling'='Tumor.eTreg_LT', # 'Tumor.eTreg_LT_cycling',
              'Tumor.eTreg_NLT'='Tumor.eTreg_NLT',
              'Tumor.eTreg_NLT-like'='Tumor.eTreg_LT',
              'Tumor.Treg-STAT1'= 'Tumor.nTreg', # 'Tumor.NA', # 'Tumor.nTreg',
              'Tumor.eTreg1_cycling'='Tumor.NA',
              'LN.eTreg_NLT'='LN.eTreg',
              'LN.eTreg_NLT-like'='LN.eTreg_NLT-like',
              'LN.eTreg_suppressive_cycling'='LN.eTreg_NLT-like', # 'LN.eTreg_NLT-like_cycling',
              'LN.Treg_LT/NLT_like'='LN.',
              'LN.Treg-LT_Stat1'='LN.cTreg')
lvdf <- seu$seurat@meta.data[,'manual_anno', drop=F] %>% 
  tibble::rownames_to_column(., "cell") %>% 
  left_join(., lv, by='cell') %>%
  mutate(manual_anno = newidmap[manual_anno]) %>%
  dplyr::filter(!manual_anno %in% c('LN.', 'Tumor.NA')) %>%
  dplyr::filter(!is.na(manual_anno))
logexp <- GetAssayData(seu$seurat, slot='data')[,lv$cell]
treg_ids2 <- treg_ids <- setNames(RColorBrewer::brewer.pal(n = 8, name = 'Dark2'),
                                  unique(lvdf$manual_anno))
names(treg_ids2) <- gsub("^.*?\\.", "", names(treg_ids2))


sig <- read.csv(file.path(PDIR, "ref", "gpy_signature", "gpy_signature.csv"), header=F)
expr <- GetAssayData(seu$seurat, slot='counts')[,lvdf$cell]
auc <- AUCell::AUCell_run(expr, list("Signature"=stringr::str_to_title(sig$V1)))
lvdf <- cbind(lvdf, 'auc'=assay(auc)[1,])
.addSummAuc <- function(lvdf, lv_i, nbreaks=10){
  lvdf <- lvdf %>% dplyr::arrange(manual_anno, !!rlang::sym(lv_i))
  
  cuts <- cut(lvdf[,lv_i], nbreaks)
  mean.cuts <- sapply(split(lvdf$auc, cuts), mean)
  lvdf$auc_mean <- mean.cuts[cuts]
  return(lvdf)
}

visualize=T
if(visualize){
  seusub <- subset(seu$seurat, cells=lvdf$cell)
  pdf("~/xfer/gpy_signature_heatmap.pdf", height = 8)
  seusub$manual_anno2 <- NA
  seusub$manual_anno2[lvdf$cell] <- lvdf$manual_anno
  Idents(seusub) <- 'manual_anno2'
  
  ScaleData(seusub, features=stringr::str_to_title(sig$V1)) %>% 
    # DoHeatmap(., cells = lvdf$cell, features=stringr::str_to_title(sig$V1),
    #           group.by='manual_anno2', size=4)
    DotPlot_scCustom(seurat_object = ., features = stringr::str_to_title(sig$V1), 
                     colors_use = viridis_plasma_dark_high, group.by='manual_anno2',
                     flip_axes = T, x_lab_rotate = TRUE)
  dev.off()
  
  lvs <- c('LV0', 'LV1', 'LV2')
  pdf(paste0("~/xfer/", sampleid, ".allLVs.pdf"), height = 9)
  
  for(lv_i in lvs){
    for(lv_j in lvs){
      ggp <- ggplot(lvdf, aes_string(x = lv_i, y = lv_j, colour = 'manual_anno'))+
        scale_color_manual(values=treg_ids) +
        geom_point(shape = 19, size = 0.7, alpha=0.5)+
        theme_classic()+
        theme(legend.position = "none", plot.margin = unit(c(0,0.1,0,0), "cm")) + 
        ggtitle(sampleid)
      
      lvdf.l <- split(lvdf, f=gsub("\\..*", "", lvdf$manual_anno))
      ggl <- lvdf %>%
        mutate(manual_anno = gsub("^.*?\\.", "", manual_anno)) %>%
        ggplot(. , aes_string(x = lv_i, fill = 'manual_anno'))+
        scale_fill_manual(values=treg_ids2) +
        geom_density(alpha = 0.5, colour = "grey25")+
        theme_cowplot()+
        theme(legend.position = "none", plot.margin = unit(c(0,0.1,0,0), "cm")) + 
        ylim(0, 1.5)
      
      ggl.LN <- ggplot(lvdf.l$LN, aes_string(x = lv_i, fill = 'manual_anno'))+
        scale_fill_manual(values=treg_ids) +
        geom_density(alpha = 0.5, colour = "grey25")+
        theme_cowplot()+
        theme(legend.position = "none", plot.margin = unit(c(0,0.1,0,0), "cm")) + 
        ylim(0, 1.5)
      ggl.T <- ggplot(lvdf.l$Tumor, aes_string(x = lv_i, fill = 'manual_anno'))+
        scale_fill_manual(values=treg_ids) +
        geom_density(alpha = 0.5, colour = "grey25")+
        theme_cowplot()+
        theme(legend.position = "none", plot.margin = unit(c(0,0.1,0,0), "cm")) + 
        ylim(0, 1.5)
      
      lvdf <- .addSummAuc(lvdf, lv_i)
      gg.loess <- lvdf %>%
        mutate(manual_anno = gsub("^.*?\\.", "", manual_anno)) %>%
        ggplot(., aes_string(x=lv_i, y='auc', col='manual_anno', fill='manual_anno')) +
        geom_smooth(method = "loess", se = TRUE, level = 0.95) +
        scale_fill_manual(values=treg_ids) +
        # geom_boxplot(aes(group = cut_width(!!rlang::sym(lv_i), 0C.1)), outlier.alpha = 0.1)
        theme_cowplot() +
        theme(legend.position = "bottom", 
              # legend.key.width = unit(0.45, "lines"),
              # legend.key.height = unit(0.45, "lines"),
              # axis.text = element_text(size = 4.5, colour = "black"),
              # axis.title = element_text(size = 5.5, colour = "black"),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              # plot.margin = unit(c(0,0.1,0,0), "cm"), 
              legend.background = element_blank()) +
        ylim(0, 0.5)
      
      
      
      print(cowplot::plot_grid(ggp, ggl.LN, ggl.T, gg.loess, ncol=1, align = "v"))
      print(cowplot::plot_grid(ggp, ggl, gg.loess, ncol=1, align = "v"))
    }
  }
  dev.off()
  
  print(paste0( sampleid, ".allLVs.pdf"))
}

