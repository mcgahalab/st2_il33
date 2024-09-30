renv::load("/cluster/projects/mcgahalab/envs/renvs/seuratv5_v2")
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
library(tradeSeq)
library(condiments)
library(reticulate)
# library(harmony)
# # library("GOstats")
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
args <- commandArgs(trailingOnly = TRUE)
group <- args[1]

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d'
setwd(PDIR)

gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
groups <- list.files(datadir, pattern="^(CD45|ST2|LN|Tumor)")

# # Read in gtf file to map ENS->Biotype
# gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
# GTF <- rtracklayer::import(gtf)
# ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))
# sym2biotype_ids <- with(GTF, setNames(gene_biotype, gene_name))


doublet_quantile_cutoff <- 0.95

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

# GOTERMS <- (as.list(GOTERM))
# GOTERMS <- sapply(GOTERMS, Term)
# library(igraph)
# library(GO.db)
# BP <- toTable(GOBPPARENTS)
# CC <- toTable(GOCCPARENTS)
# MF <- toTable(GOMFPARENTS)
# GOg <- graph.data.frame( rbind(BP,CC,MF) )

# gg <- geneset_l[[4]]
# 
# goid <- sample(unique(gg$gs_exact_source), size=2)
# goid <- gg$gs_exact_source[sample(grep("INTERLEUKIN_6", gg$gs_name), size=2)]
# goid <- c('GO:0032675', 'GO:0001818')
# goid <- c('GO:0002730', 'GO:0002719')
# goid <- c('GO:0002704', 'GO:0002700')
# goid <- c('GO:0002704', 'GO:0002700', 'GO:0010629')
# 
# gg[match(goid, gg$gs_exact_source),]
# library(GO.db)
# library(igraph)
# BP <- toTable(GOBPPARENTS)
# CC <- toTable(GOCCPARENTS)
# MF <- toTable(GOMFPARENTS)
# g <- graph.data.frame( rbind(BP,CC,MF) )
# 
# 
# lca <- function(graph, ...) {
#   dots = c(...)
#   path = ego(graph, order=length(V(graph)), nodes=dots, mode="in")
#   V(graph)[[(Reduce(intersect, path))]]
# }
# V(g)[[lca(g, goid])]]
# 
# pn = 'GO:0002731'
# V(g)[[goid]]
# V(g)[[pn]]
# shortest_paths(g, from=pn, to=goid[1])
# shortest_paths(g, from=pn, to=goid[2])
# shortest_paths(g, from=pn, to=goid[3])
# shortest_paths(g, from='GO:0002719', to=goid[3])
#                
# BP[grep(goid[3], BP[,2]),]
# 
# 
# 
# 
# parent_node = attr(lca(g, goid), 'names')
# lca_id <- sapply(parent_node, function(pn){
#   c(length(shortest_paths(g, from=pn, to=goid[1])$vpath[[1]]),
#     length(shortest_paths(g, from=pn, to=goid[2])$vpath[[1]]),
#     length(shortest_paths(g, from=pn, to=goid[3])$vpath[[1]]))
# }) %>% apply(., 2, median) %>% which.min
# V(g)[[lca_id[1]]]
# 
# V(g)[[goid]]
# 
# 
# tatus_on <- filter(vertices, status == 1)$vertex  
# subcomponent(graph=g, v=goid, mode='in')
# vids <- ego(graph=g, order=1, nodes=goid, mode='in')
# subg <- subgraph(graph=g, path[[1]])
# 
# pdf("~/xfer/x.pdf")
# plot(subg)
# dev.off()
# 
# BP[match(goid, BP[,1]),]
# 
# x <- map(status_on, subcomponent, graph = g, mode = "in") %>%   
#   unlist() %>%   names()
# 
# gg[match(id, gg$gs_exact_source),]
# gg[match(path[[1]], gg$gs_exact_source),]


# # Read in gtf file to map ENS->Biotype
# gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
# GTF <- rtracklayer::import(gtf)
# ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

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

##-- Recurring functions 2 ----
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

.searchReference <- function(dat, searchterms){
  for(colx in names(searchterms)){
    # message(colx)
    dat <- dat %>%
      dplyr::filter(!!rlang::sym(colx) == searchterms[[colx]])
  }
  return(dat)
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
saveRDS(seus_preproc_l, file=file.path(PDIR, "data", "seurat_obj", "seus_preproc.rds"))


### Doublet calling
seus_preproc_l <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_preproc.rds"))
doublet.cutoff <- 0.30
doublet_clusters <- list("B1_ST2_LN_KO_3d"=as.character(c(25,28)),
                         "B1_ST2_LN_KO_Un"=as.character(c(7, 25,28)),
                         "B1_ST2_LN_WT_3d"=as.character(c(20,28,36)),
                         "B1_ST2_LN_WT_Un"=as.character(c(28)),
                         "B1_ST2_Tumor_WT_Un"=as.character(c(7,17)),
                         "B1_ST2_Tumor_WT_3d"=as.character(c()),
                         "B1_ST2_Tumor_KO_Un"=as.character(c()),
                         "B1_ST2_Tumor_KO_3d"=as.character(c()),
                         "B2_ST2_LN_KO_7d"=as.character(c(7, 28)),
                         "B2_ST2_LN_KO_Un"=as.character(c(7,17,28)),
                         "B2_ST2_LN_WT_7d"=as.character(c(28)),
                         "B2_ST2_LN_WT_Un"=as.character(c(16,28)),
                         "B2_ST2_Tumor_KO_7d"=as.character(c()),
                         "B2_ST2_Tumor_KO_Un"=as.character(c()),
                         "B2_ST2_Tumor_WT_7d"=as.character(c()),
                         "B2_ST2_Tumor_WT_Un"=as.character(c(28)),
                         "B2_CD45_LN_KO_7d"=as.character(c(13,31,33)),
                         "B2_CD45_LN_KO_Un"=as.character(c(13,29,31)),
                         "B2_CD45_LN_WT_7d"=as.character(c(13,33)),
                         "B2_CD45_LN_WT_Un"=as.character(c(13)),
                         "B2_CD45_Tumor_KO_7d"=as.character(c(31)),
                         "B2_CD45_Tumor_KO_Un"=as.character(c()),
                         "B2_CD45_Tumor_WT_7d"=as.character(c(31)),
                         "B2_CD45_Tumor_WT_Un"=as.character(c(31))
)


seus_preproc_l <- lapply(setNames(names(seus_preproc_l),names(seus_preproc_l)), function(seuid){
  seu <- seus_preproc_l[[seuid]]
  sce <- as.SingleCellExperiment.Seurat5(seu)
  cxds_f <- file.path(PDIR, "results", "doublets", paste0("cxds.", seuid, ".rds"))
  bcds_f <- file.path(PDIR, "results", "doublets", paste0("bcds.", seuid, ".rds"))
  if(!file.exists(cxds_f) & !file.exists(bcds_f)){
    print("Generating...")
    sce_cxds <- lapply(sce, runCxds, n=25)
    sce_bcds <- lapply(sce, runBcds, n=25)
    saveRDS(sce_cxds, file=cxds_f)
    saveRDS(sce_bcds, file=bcds_f)
  } else {
    print("Reading in...")
    sce_cxds <- readRDS(file=cxds_f)
    sce_bcds <- readRDS(file=bcds_f)
  }
  
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
saveRDS(seus_preproc_l, file=file.path(PDIR, "data", "seurat_obj", "seus_preproc.rds"))

#--- 1.c) Integrating samples ----
seus_preproc_l <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_preproc.rds"))

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

#--- b) Carrying over previous annotations ----
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.rds"))
map <- c(seul$Tumor$manual_anno,seul$LN$manual_anno)
map <- setNames(as.character(map), .relabelid(names(map)))
seus$ST2@meta.data[,'manual_anno'] <- map[Cells(seus$ST2)]
rm(seul);gc()
seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_CD45_7D.split.final.rds"))
map <- c(seul$Tumor$manual_anno,seul$LN$manual_anno)
map <- setNames(as.character(map), .relabelid(names(map)))
seus$CD45@meta.data[,'manual_anno'] <- map[Cells(seus$CD45)]
rm(seul);gc()

saveRDS(seus, file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
#--- c) TReg annotation ---
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))


#--- d) Visualizations ----
dps <- lapply(setNames(c('ST2', 'CD45'), c('ST2', 'CD45')), function(seuid){
  dp1 <- DimPlot(seus[[seuid]], label = TRUE, repel = TRUE, reduction = "umap.mnn",
                  group.by="seurat_clusters", pt.size=2, shuffle=TRUE,raster=T)
  dp2 <- DimPlot(seus[[seuid]], label = TRUE, repel = TRUE, reduction = "umap.mnn",
                  group.by="orig.ident", pt.size=2, shuffle=TRUE,raster=T)
  dp3 <- DimPlot(seus[[seuid]], label = TRUE, repel = TRUE, reduction = "umap.mnn",
                  group.by="manual_anno", pt.size=2, shuffle=TRUE,raster=T)
  list('dp1'=dp1, 'dp2'=dp2, 'dp3'=dp3)
})

pdf(file.path(outdir, "umap", "annotated.pdf"), width = 10, height = 7)
dps$ST2
dps$CD45
dev.off()
file.copy(file.path(outdir, "umap", "annotated.pdf"), to = "~/xfer", overwrite = T)



#--- e) Loupe ----
if(!exists('seus')) seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
doublets <- read.csv("~/xfer/sara_doublets.csv") %>%
  with(., setNames(Sara_anno, Barcode))

# seu_st2 <- .monocle3Cluster(seus$ST2)


dir.create(file.path(PDIR, "results", "cloupe"), showWarnings = F)
lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  seu$doublets <- NA
  seu@meta.data[, 'doublets'] <- doublets[Cells(seu)]
  cols <- c('orig.ident', 'Batch', 'Sort', 'Tissue', 'Condition', 'Timepoint',
            'seurat_clusters', 'mt_outlier', 'singlet', 'mnn_clusters', 
            'immgen.*cluster', 'mouse.*cluster', 'manual', 'doublets')
  cols <- sapply(cols, grep, x=colnames(seu@meta.data), ignore.case=T, value=T) %>%
    unlist %>% as.character
  meta <- seu@meta.data[,cols]
  for(i in seq_along(colnames(meta))){
    meta[,i] <- factor(as.character(meta[,i]))
  }
  seu@meta.data <- meta
  seu[['RNA2']] <- CreateAssayObject(counts=seu@assays$RNA$counts)
  DefaultAssay(seu) <- 'RNA2'
  # seu <- subset(seu, cells=sample(Cells(seu), size=1000))
  loupeR::create_loupe_from_seurat(seu, output_dir=file.path(PDIR, "results", "cloupe"),
                                   output_name=paste0("seu_", seuid), force=T)
  file.copy(file.path(PDIR, "results", "cloupe", paste0("seu_", seuid, ".cloupe")),
            to = "~/xfer", overwrite = T)
  cat(paste0("xfer seu_", seuid, ".cloupe\n"))
  
  # success <- loupeR:::create_loupe(
  #   count_mat=loupeR:::select_assay(seu)$RNA2@counts,
  #   clusters=loupeR:::select_clusters(seu, dedup=F),
  #   projections=loupeR:::select_projections(seu),
  #   output_dir=file.path(PDIR, "results", "cloupe"),
  #   output_name=paste0("seu_", seuid),
  #   force=T)
  # ok <- loupeR:::create_hdf5(count_mat=loupeR:::select_assay(seu)$RNA2@counts,
  #                            clusters=loupeR:::select_clusters(seu, dedup=F),
  #                            projections=loupeR:::select_projections(seu),
  #                            h5path=file.path(PDIR, "results", "cloupe", "test.h5"),
  #                            feature_ids=NULL,
  #                            seurat_obj_version='V5')
  # ok <- loupeR:::louper_create_cloupe(file.path(PDIR, "results", "cloupe", "test.h5"),
  #                            output_dir=file.path(PDIR, "results", "cloupe"),
  #                            output_name=paste0("seu_", seuid),
  #                            force = T)
  # input_flag <- sprintf("--input=%s", h5path)
  # output_flag <- sprintf("--output=%s", loupe_path)
  # args <- c("create", input_flag, output_flag)
  # '/cluster/home/quever/.local/share/R/loupeR/louper create --input test.h5 --output test.cloupe'
})


#--- f) Manual annotation ----
if(!exists('seus')) seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
manual_mode <- 'allCells'

# these annotations were created from sara's manual revision from the loupe file
if(manual_mode == 'allCells'){
  ln_anno <- file.path(PDIR, "ref", "CD45ST2_LN_annotations.csv")
  tumor_anno <- file.path(PDIR, "ref", "CD45ST2_Tumor_annotations.csv")
  colid <- 'manual_anno2'
} else if(manual_mode=='Treg'){
  ln_anno <- file.path(PDIR, "ref", "Sara.anno_final_CD45ST2Treg_LN.v2.csv")
  tumor_anno <- file.path(PDIR, "ref", "Sara.anno_final_CD45ST2_TUMOR_Treg.v2.csv")
  colid <- 'treg_anno'
}


anno <- lapply(c(ln_anno, tumor_anno), read.csv, col.names=c('Barcode', 'Anno')) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  mutate(Anno = gsub(" $", "", Anno) %>% 
           gsub(" ", "_", .) %>%
           gsub("eTregNLT", "eTreg_NLT", .))
  # mutate(Anno = gsub("doublet_in", "doublets_in", Anno) %>%
  #          gsub("monocyte-derived mac$", "monocyte-derived macrophages", .) %>%
  #          gsub("^CD3", "Cd3", .) %>%
  #          gsub("cycling $", "cycling", .) %>%
  #          gsub(" ", "_", .))
seus$ST2 <- FindNeighbors(seus$ST2, reduction = "integrated.mnn", 
                          dims = 1:30, return.neighbor=T)
idx <- which(anno$Anno == '')
missing_anno <- sapply(anno$Barcode[idx], function(cell_i){
  cells <- Seurat::TopNeighbors(seus$ST2@neighbors$RNA.nn, 
                       cell=cell_i,
                       n=20)
  majority_celltype <- anno %>% filter(Barcode %in% cells,
                                       Anno !='') %>% 
    pull(Anno) %>% table %>% sort %>% tail(., 1)
  if(length(majority_celltype)==0) majority_celltype <- table(cell_i)
  ifelse(majority_celltype < 10, 'NA', names(majority_celltype))
})
anno[idx,'Anno'] <- as.character(missing_anno)
seus$ST2@meta.data[,colid] <- 'NA'
seus$ST2@meta.data[anno$Barcode,colid] <- anno$Anno

saveRDS(seus, file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))

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
# seus$LN@meta.data[,'treg_anno'] <- map[rownames(seus$LN@meta.data)]
# seus$Tumor@meta.data[,'treg_anno'] <- map[rownames(seus$Tumor@meta.data)]

## Fix the APC doublets
seus$LN$manual_anno2 <- gsub("APC doublets_in ", "", seus$LN$manual_anno2)
seus$LN$manual_anno2 <- gsub("^DC$", "monocyte-derived DC", seus$LN$manual_anno2)
seus$Tumor$manual_anno2 <- gsub("APC doublets_in ", "", seus$Tumor$manual_anno2) %>%
  gsub("Macrophages", "macrophages", .) 
seus$Tumor$manual_anno2 <- gsub("^DC$", "monocyte-derived DC", seus$Tumor$manual_anno2)
saveRDS(seus, file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))

## Make cloupe files
lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  cols <- c('orig.ident', 'Batch', 'Sort', 'Tissue', 'Condition', 'Timepoint',
            'seurat_clusters', 'manual', 'treg_anno')
  cols <- sapply(cols, grep, x=colnames(seu@meta.data), ignore.case=T, value=T) %>%
    unlist %>% as.character
  meta <- seu@meta.data[,cols]
  meta[is.na(meta)] <- 'NA'
  for(i in seq_along(colnames(meta))){
    meta[,i] <- factor(as.character(meta[,i]))
  }
  seu@meta.data <- meta
  seu[['RNA2']] <- CreateAssayObject(counts=seu@assays$RNA$counts)
  DefaultAssay(seu) <- 'RNA2'
  # seu <- subset(seu, cells=sample(Cells(seu), size=1000))
  loupeR::create_loupe_from_seurat(seu, output_dir=file.path(PDIR, "results", "cloupe"),
                                   output_name=paste0("CD45ST2_", seuid), force=T)
  file.copy(file.path(PDIR, "results", "cloupe", paste0("CD45ST2_", seuid, ".cloupe")),
            to = "~/xfer", overwrite = T)
  cat(paste0("xfer CD45ST2_", seuid, ".cloupe\n"))
})


## Create the publication dimplots
grp <- 'manual_anno2'
grps <- sapply(seus, function(i) unique(i@meta.data[,grp])) %>% 
  unlist %>% unique
colors_use <- scCustomize::scCustomize_Palette(num_groups = length(grps)+1, 
                                               ggplot_default_colors = FALSE, 
                                               color_seed = 231)[-2] %>%
  setNames(., grps)
pdps <- lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  colors_use_i <- as.character(colors_use[unique(seu@meta.data[,grp])])
  seu@reductions[['UMAP']] <- seu@reductions[['umap.mnn_nn60_md0.3']]
  colnames(seu@reductions[['UMAP']]@cell.embeddings) <- paste0("UMAP_", c(1:2))
  seu@reductions[['UMAP']]@key <- 'UMAP_'
  seu@meta.data[,seuid] <- seu@meta.data[,grp]
  publicationDimPlot(seu, grp=seuid, simplify_labels=TRUE, 
                     colors_use=colors_use_i, pt.size=1, aspect.ratio=1,
                     reduction='UMAP',
                     legend.text = element_text(size=10),
                     labs(title=seuid))
})
pdf(file.path(PDIR, "results", "umap", "umap_CD45ST2_tumor_ln.allCells.pdf"), width=7, height=5)
pdps
dev.off()
file.copy(file.path(PDIR, "results", "umap", "umap_CD45ST2_tumor_ln.allCells.pdf"), to="~/xfer", overwrite = T)

#--- a) inferCNV tumor contaminant ----
library(infercnv)
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
infercnv_dir <- '/cluster/projects/mcgahalab/ref/scrna/infercnv'
infercnv_f <- "mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions"
seuid <- 'Tumor'
seu <- seus[[seuid]]
seu$highPga <- 'Low'
infercnv_obj = infercnv::CreateInfercnvObject(
  raw_counts_matrix=GetAssayData(seu, layer='counts', assay='RNA'),
  annotations_file=seu@meta.data[,'manual_anno2', drop=F],
  delim="\t",
  gene_order_file=file.path(infercnv_dir, infercnv_f),
  ref_group_names=NULL) 

dir.create(file.path(PDIR, "results", "infercnv", paste0("st2_", seuid)), showWarnings = F)
infercnv_obj = infercnv::run(
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=file.path(PDIR, "results", "infercnv", paste0("st2_", seuid)), 
  cluster_by_groups=TRUE, 
  denoise=TRUE,
  HMM=TRUE)
infercnv_obj <- readRDS(file.path(PDIR, "results", "infercnv", paste0("st2_", seuid), 'run.final.infercnv_obj'))
thresholds <- read.table(file.path(PDIR, "results", "infercnv", paste0("st2_", seuid), 'infercnv.heatmap_thresholds.txt'))
infercnv::plot_cnv(infercnv_obj, out_dir=file.path(PDIR, "results", "infercnv", paste0("st2_", seuid)),
                   output_filename = "st2_infercnv", output_format="pdf")

# pga <- colSums(abs(infercnv_obj@expr.data - 1))
pga <- colSums((infercnv_obj@expr.data <= min(thresholds$V1)) | 
  (infercnv_obj@expr.data >= max(thresholds$V1)))
seu@meta.data[names(which(pga >= quantile(pga, 0.9))),]$highPga <- 'High'

cols <- c('orig.ident', 'Batch', 'Sort', 'Tissue', 'Condition', 'Timepoint',
          'seurat_clusters', 'manual', 'highPga')
cols <- sapply(cols, grep, x=colnames(seu@meta.data), ignore.case=T, value=T) %>%
  unlist %>% as.character
seu@meta.data <- seu@meta.data[,cols]
seu[['RNA2']] <- CreateAssayObject(counts=seu@assays$RNA$counts)
DefaultAssay(seu) <- 'RNA2'
# seu <- subset(seu, cells=sample(Cells(seu), size=1000))
loupeR::create_loupe_from_seurat(seu, output_dir=file.path(PDIR, "results", "cloupe"),
                                 output_name=paste0("CD45ST2_", seuid, ".infercnv"), force=T)
file.copy(file.path(PDIR, "results", "cloupe", paste0("CD45ST2_", seuid, ".infercnv.cloupe")),
          to = "~/xfer", overwrite = T)
cat(paste0("xfer CD45ST2_", seuid, ".infercnv.cloupe\n"))

#--- b) Integrate with B16 ----
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
b16seu <- readRDS(file=file.path(dirname(PDIR), "b16_scrna", "data", "tmp", "2_rnamnn_seu.rds"))

seu <- DietSeurat(merge(seus$LN, b16seu), layers='counts', assay='RNA')
idx <- which(seu$orig.ident == 'SeuratProject')
seu$orig.ident[idx] <- seu$cell_type[idx]
seu <- JoinLayers(seu)
seu <- DietSeurat(seu, layers='counts')
seu <- split(seu, f=seu$orig.ident)
Idents(seu) <- 'orig.ident'
DefaultAssay(seu) <-'RNA'
seu <- NormalizeData(seu) %>%
  FindVariableFeatures(.)  %>%
  ScaleData(.)  %>%
  RunPCA(.) %>%
  FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
  FindClusters(., resolution = 0.5, cluster.name = "unintegrated_clusters_std") %>%
  RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated_std")

Idents(seu) <- 'orig.ident'
seu <- IntegrateLayers(
  object = seu, method = SeuratWrappers::FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = T
)
seu <- seu %>% 
  FindNeighbors(., reduction = "integrated.mnn", dims = 1:30) %>%
  FindClusters(., resolution = 0.9, cluster.name = "mnn_clusters") %>% 
  RunUMAP(., reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn",
          n.neighbors = 25L, min.dist = 0.3)
DefaultAssay(seu) <- 'RNA'
seu <- JoinLayers(seu)

mat <- GetAssayData(seu, assay='RNA', layer='counts')
clusters <- seu@meta.data[,c('manual_anno2', 'Cluster1', 'Cluster2', 'Cell_cycle', 'Kdm5b', 'cell_type')]
clusters[is.na(clusters)] <- 'NA'
projections <- loupeR:::select_projections(seu)
projections <- c(projections, list('pca'=seu@reductions$pca@cell.embeddings[,c(1:6)]))
makeLoupe(mat=mat, meta=clusters, projections=projections,
        output_dir=file.path(PDIR, "results", "cloupe"), 
        output_name='b16_ln')
file.copy(file.path(PDIR, "results", "cloupe", "b16_ln.cloupe"), to="~/xfer", overwrite = T)


#####################################
#### 4. Pathway activity scoring ####
# seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))

outdir_pway <- file.path(PDIR, "results", "pathway")
dir.create(outdir_pway, recursive = F)

# seul <- lapply(seul, recode_list, newcol='manual_clusters', grp='orig.ident', anno=F)
# seul <- lapply(seul, recode_list, newcol='manual_anno', grp='orig.ident', anno=T)
method <- 'aucell'
anno_ident <- 'manual_anno2' #'manual_anno', 'manual_clusters'
pval_sig <- 0.05
fc_sig <- 0.5
sigdig <- 4


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

#--- a) AUCell ----
# custom_l <- list('Custom'=split(customgs$Gene, f=customgs$Geneset))

msig_l <- lapply(split(gm$ENSEMBL$SYMBOL[msig_ds$entrez_gene], msig_ds$gs_name), function(i){
  i[!is.na(i)]
})
msig_l <- split(msig_l, f=gsub(":.*", "", names(msig_l)))

# msig_l <- custom_l

lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  DefaultAssay(seu) <- 'RNA'
  expr <- GetAssayData(seu, slot='counts')
  expr <- expr[rowSums(expr)>=50, ]
  
  scores_l <- list()
  for(id in names(msig_l)){
    print(paste0(id, "..."))
    dir.create(file.path(outdir_pway, method, seuid), recursive = T, showWarnings = F)
    outf <- file.path(outdir_pway, method, seuid, paste0("gprofiler.", id, ".rds"))
    msig_ds <- data.frame(gs_name=gsub("[0-9]*$", "", names(unlist(msig_l[id]))),
                          entrez_gene=as.character(unlist(msig_l[id]))) %>%
      filter(entrez_gene != '') %>%
      mutate(gs_name = gsub("^.*?\\.", "", gs_name) %>%
               gsub("\\.ENSMUSG", "", .))
    
    if(!file.exists(outf)){
      print("Running aucell..")
      if(method=='aucell'){
        scores <- aucellFun(msig_ds = msig_ds, 
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
  # scores <- readRDS(file.path(outdir_pway, method, "aucell.Custom.rds"))
  # scores_l[['Custom']] <- scores
  saveRDS(scores_l, file=file.path(outdir_pway, method, seuid, paste0("aucell.rds")))
})


#--- b) SCPA ----
msig_l_scpa <- msig_ds %>% 
  rename_with(., ~'Pathway', .cols='gs_name') %>% 
  mutate(Genes=gm$ENSEMBL$SYMBOL[entrez_gene]) %>% 
  filter(!is.na(Genes))
msig_l_scpa <- split(msig_l_scpa %>% 
                       dplyr::select(c('Pathway', 'Genes')),
                     f=msig_l_scpa$Pathway)
# custom_l_scpa <- customgs %>%
#   magrittr::set_colnames(c('Genes', 'Pathway')) %>%
#   split(., f=.$Pathway)

for(seuid in names(seus)){
  print(paste0(" >>> ", seuid, " <<< "))
  outf <- file.path(outdir_pway, 'scpa', paste0(seuid, ".scpa.rds"))
  
  seu <- seus[[seuid]]
  data=as(seu[['RNA']]$data, 'dgCMatrix')
  seu[['RNA2']] <- CreateAssayObject(data[which(rowSums(data) > 50),])
  
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
  
  if(!file.exists(outf)){
    markers_all <- apply(xt_all[1:50,], 1, function(ids){
      print(paste0("Testing ", ids[4], " - ", ids[5]))
      ids_test <- factor(c(ids[1], ids[2]))
      ids_test <- relevel(ids_test, ids[3])
      ord <- order(ids_test)
      
      Idents(seu) <- 'group'
      scpa_out <- tryCatch({
        compare_seurat(seu,
                       assay='RNA2',
                       group1 = 'group', 
                       group1_population = ids[1:2],
                       pathways = msig_l_scpa,  # custom_l_scpa,
                       max_genes=1000,
                       parallel=F) %>%
          mutate(Database = gsub(":.*", "", Pathway),
                 id1=ids[4],
                 id2=ids[5],
                 fullid1=ids[1],
                 fullid2=ids[2],
                 celltype=gsub("^.*\\.", "", ids[1])) 
      }, error=function(e){print("Error"); NULL})
      return(scpa_out)
    })
    saveRDS(markers_all, file=outf)
  } else {
    markers_all <- readRDS(outf)
  }
  
  dir.create(file.path(outdir_pway, 'scpa', seuid))
  for(id in names(markers_all)){
    write.table(markers_all[[id]], 
                file=file.path(outdir_pway, 'scpa', seuid, paste0("SCPA_", id, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  }
}


for(compid in names(comp_lists)){
  print(compid)
  spl <- comp_lists[[compid]]
  grpspl <- .createGrps(spl)
  # outf <- file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".custom.rds"))
  outf <- file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".rds"))
  
  if(!file.exists(outf)){
    markers <- lapply(colnames(grpspl), function(grpid){
      ids <- grpspl[,grpid]
      print(paste(ids, collapse="//"))
      
      scpa_out <- compare_seurat(cd8seu,
                                 assay='RNA2',
                                 group1 = 'group', 
                                 group1_population = ids,
                                 pathways = msig_l_scpa,  # custom_l_scpa,
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
    markers <- readRDS(file=outf)
  }
}


for(compid in names(comp_lists)){
  markers <-  readRDS(file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".rds")))
  markers_custom <- readRDS(file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".custom.rds")))
  markers <- rbind(markers, markers_custom)
  dupidx <- with(markers, duplicated(paste0(Pathway, id1, id2)))
  if(any(dupidx)){
    stop(paste0(compid, ": Duplicated indices found at ", paste(which(dupidx), collapse=",")))
  } else {
    saveRDS(markers, file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".rds")))
  }
}



#--- c) Visualization and Differential AUCell scores ----
list.files(file.path(outdir_pway, method))
scores_l <- readRDS(file=file.path(outdir_pway, method, "aucell.rds"))
scores <- lapply(scores_l, assay) %>% do.call(rbind, .)
cd8seu[['AUCell']] <- CreateAssayObject(data = scores)
DefaultAssay(cd8seu) <- 'AUCell'
Idents(cd8seu) <- 'group'

for(compid in rev(names(comp_lists))){
  print(paste0(">>>    ", compid, "    <<<"))
  spl <- comp_lists[[compid]]
  # Read in SCPA data
  scpa_scores <-  readRDS(file=file.path(outdir_pway, "scpa", paste0("SCPA_", compid, ".rds"))) %>%
    dplyr::select(c(Pathway, id1, id2, Pval, adjPval)) %>%
    magrittr::set_colnames(., c('pathway', 'id.1', 'id.2', 'scpa.pval', 'scpa.adjPval')) %>%
    mutate(pathway=gsub("_", "-", pathway))
  
  outdir <- file.path(PDIR, "results", "pathway", compid)
  
  if(!file.exists(file.path(outdir, paste0("Pathway_", compid, ".rds")))){
    markers <- apply(.createGrps(spl), 2, function(ids){
      print(paste0("Processing: ", paste(ids, collapse=",")))
      M <- FindMarkers(cd8seu, assay='AUCell', 
                       group.by='group',
                       ident.1=ids[1], ident.2=ids[2], 
                       logfc.threshold=0) %>%
        tibble::rownames_to_column("pathway") %>%
        mutate(id.1=ids[1],
               id.2=ids[2],
               p_val = round(p_val, sigdig),
               avg_log2FC=round(avg_log2FC, sigdig),
               p_val_adj = ifelse(p_val_adj < 10^(-1*sigdig), 
                                  p_val_adj, 
                                  round(p_val_adj, sigdig))) %>%
        relocate(., c('id.1', 'id.2'), .after='pathway')
      cells <- Seurat:::IdentsToCells(object = cd8seu, 
                                      ident.1 = ids[1], 
                                      ident.2 = ids[2], 
                                      cellnames.use = Cells(cd8seu)) 
      D <- apply(GetAssayData(cd8seu, assay='AUCell'), 1, function(i){
        x <- effsize::cohen.d(i[cells$cells.1], i[cells$cells.2])
        setNames(c(x$conf.int, x$estimate), c('D.lower', 'D.upper', 'D'))
      }) %>%
        t %>% as.data.frame %>%
        tibble::rownames_to_column("pathway") %>%
        mutate(D.lower = round(D.lower, sigdig),
               D.upper =round(D.upper, sigdig),
               D =round(D, sigdig))
      
      MD <- left_join(M, D, by='pathway')
      return(MD)
    })
    
    saveRDS(markers, file=file.path(outdir, paste0("Pathway_", compid, ".rds")))
  } else {
    markers <- readRDS(file=file.path(outdir, paste0("Pathway_", compid, ".rds")))
  }
  
  for(grpid in names(markers)){
    m <- left_join(markers[[grpid]], scpa_scores, by=c('pathway', 'id.1', 'id.2'))
    dir.create(outdir, recursive = T, showWarnings = F)
    write.table(m, file=file.path(outdir, paste0("Pathway_", grpid, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  }
}



ensure.same.wilcox.results <- T
if(ensure.same.wilcox.results){
  ## Tests to make sure the FindMarkers p-value gives the same value as a simple wilcox.test
  ass <- GetAssayData(cd8seu_sub, assay='AUCell')
  rowidx <- which(rownames(ass) == rownames(M)[1])
  cells1 <- Cells(cd8seu_sub)[cd8seu_sub$group == 'P14-ISO_Tpex']
  cells2 <- Cells(cd8seu_sub)[cd8seu_sub$group == 'P14-DB_Tpex']
  c(wilcox.test(ass[rowidx,cells1], ass[rowidx,cells2])$p.value,
    M[1,]$p_val)
}


mat <- GetAssayData(cd8seu, assay='AUCell')
makeLoupe(mat=mat, meta=clusters, projections=projections,
          output_dir=file.path(PDIR, "results", "loupe"), 
          output_name='p14_pathway_activity')
file.copy(file.path(PDIR, "results", "loupe", "p14_pathway_activity.cloupe"), to="~/xfer", overwrite = T)

effsize::cohen.d(rnorm(100, mean=0, sd=1), rnorm(100, mean=0, sd=1))


#--- a) CD45: Pagoda2 activity score ####
cellgrp <- 'CD45' # 'ST2', 'CD45'
DefaultAssay(seus[[toupper(cellgrp)]]) <- 'RNA'
expr <- GetAssayData(seus[[toupper(cellgrp)]], slot='counts')
expr <- expr[rowSums(expr)>=50, ]
method <- 'aucell'

for(id in names(msig_l)[4]){
  print(paste0(id, "..."))
  dir.create(file.path(outdir_pway, "pathways", method), recursive = T, showWarnings = F)
  outf <- file.path(outdir_pway, "pathways", method, paste0(method, ".", cellgrp, ".", id, ".rds"))
  msig_ds <- data.frame(gs_name=gsub("[0-9]*$", "", names(unlist(msig_l[id]))),
                        entrez_gene=as.character(unlist(msig_l[id]))) %>%
    filter(entrez_gene != '') %>%
    mutate(gs_name = gsub("^.*?\\.", "", gs_name) %>%
             gsub("\\.ENSMUSG", "", .))
  
  if(!file.exists(outf)){
    print("Running aucell..")
    if(method=='aucell'){
      scores <- aucellFun(msig_ds = msig_ds, 
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
}


## 
scores_l <- lapply(c('ST2', 'CD45'), function(cellgrp){
  lapply(names(msig_l), function(id){
    tryCatch({
      readRDS(file=file.path(outdir_pway, "pathways", method, 
                             paste0(method, ".", cellgrp, ".", id, ".rds")))
    }, error=function(e){NULL})
  })
})
names(scores_l) <- c('ST2', 'CD45')

############################
#### 5. SCENIC analysis ####
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
outdir <- file.path(PDIR, "results", "regulons", "loom")
dir.create(outdir, recursive=T, showWarnings = F)
overwrite_regulon <- F

#--- a) Generate loom file for pySCENIC ----
for(seuid in names(seus)){
  # code from https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html#run-scenic-vis-pyscenic
  exprMat <- seus[[seuid]]@assays$RNA$data
  cellInfo <- seus[[seuid]]@meta.data
  
  loci1 <- which(rowSums(exprMat>0) > 0.01*ncol(exprMat))
  exprMat_filter <- exprMat[loci1, ]
  
  add_cell_annotation <- function(loom, cellAnnotation){
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
  
  loomf <- file.path(outdir, paste0(seuid, ".seurat2loom.loom"))
  if(!file.exists(loomf)){
    loom <- build_loom(loomf, dgem=exprMat_filter)
    loom <- add_cell_annotation(loom, cellInfo)
    close_loom(loom)
  }
}

#--- b) Post pySCENIC analysis ----
loom <- open_loom(file.path(PDIR, "results", "regulons", "pyscenic", 'cd8seu_scenic.loom'))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC), 
        file = file.path(outdirregulon, 'seu_regulon.rds'))
if((!'regulons' %in% names(cd8seu@assays)) | overwrite_regulon){
  cd8seu[['regulons']] <- CreateAssayObject(data=assay(regulonAUC))
  cd8seu <- saveRDS(cd8seu, file.path(PDIR, "data", "P14_seurs_with_geneset_enrichment_scores.rds"))
}


##############################
#### 6. Cell type subsets ####
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "seus_annotated.rds"))
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))

outdir_pway <- file.path(PDIR, "results", "pathway")
method <- 'aucell'
ids <- if(exists("msig_l")) names(msig_l) else c('GO', 'HP', 'MIRNA', 'REAC', 'WP')
scores_l <- lapply(c('ST2', 'CD45'), function(cellgrp){
  lapply(ids, function(id){
    tryCatch({
      readRDS(file=file.path(outdir_pway, "pathways", method, 
                             paste0(method, ".", cellgrp, ".", id, ".rds")))
    }, error=function(e){NULL})
  }) %>%
    setNames(.,ids)
})
names(scores_l) <- c('ST2', 'CD45')

#--- a) CD45 CD8 cells ----
cytodir <- file.path(PDIR, "results", "pathway", "pathways", "cytoscape")
cd45_cd8_seu <- readRDS(file=file.path(PDIR, "cellsubset", "CD45_Tumor_CD8.seuratobj.rds"))
newlabel <- .relabelid(Cells(cd45_cd8_seu))
cd45_cd8_seu <- RenameCells(cd45_cd8_seu, new.names=newlabel)


## Subsetting for a specific celltype
seu_cd45 <- subset(seus$CD45, cells=Cells(cd45_cd8_seu))
cd45_cd8_seu <- subset(cd45_cd8_seu, cells=Cells(seu_cd45))
cd45_cd8_reac <- scores_l$CD45$REAC[,Cells(seu_cd45)]
seu_cd45[['umap.cd8']] <- CreateDimReducObject(embeddings = Embeddings(cd45_cd8_seu, 'umap'), 
                                               key = 'cd8umap_', assay = 'RNA')

reds <- c('umap.unintegrated', 'umap.mnn', 'umap.cd8') %>%
  setNames(., .)
cd45_cd8_reac <- (cd45_cd8_reac + abs(min(cd45_cd8_reac)))
f <- makeLoupe(mat=cd45_cd8_reac, 
               meta=seu_cd45@meta.data %>%
                 select(c(orig.ident, Condition, Timepoint, immgen.fine.cluster.CD45, 
                          mouse.fine.cluster.CD45, manual_anno)),
               projections=lapply(reds, function(i) Embeddings(seu_cd45, i)),
               output_dir=file.path(PDIR, "results", "cloupe"), 
               output_name="cd45_cd8_pathways")
file.copy(f, to="~/xfer", overwrite=T)

## FindMarkers on the pathways
seu_cd45[[method]] <- if(method=='pagoda2'){
  CreateAssayObject(data=as(cd45_cd8_reac, 'dgCMatrix'))
} else if(method=='aucell'){
  CreateAssayObject(data=as(assay(cd45_cd8_reac), 'dgCMatrix'))
}
DefaultAssay(seu_cd45) <- method
seu_cd45$ident <- gsub("^.*Tumor_", "", seu_cd45$orig.ident)
Idents(seu_cd45) <- 'ident'
combos <- combn(as.character(unique(Idents(seu_cd45))), m=2)
markers <- apply(combos, 2, function(i){
  FindMarkers(seu_cd45, 
              ident.1=i[1],
              ident.2=i[2],
              logfc.threshold = 0)
}) %>%
  setNames(., apply(combos, 2, paste, collapse="."))
names(markers) <- paste0("All.", names(markers))

Idents(seu_cd45) <- 'manual_anno'
markers_by_celltype <- lapply(levels(Idents(seu_cd45)), function(celltype){
  seu_cd45_celltype <- subset(seu_cd45, ident=celltype)
  apply(combos, 2, function(i){
    Idents(seu_cd45_celltype) <- 'ident'
    FindMarkers(seu_cd45_celltype, 
                ident.1=i[1],
                ident.2=i[2],
                logfc.threshold = 0)
  }) %>%
    setNames(., apply(combos, 2, paste, collapse="."))
})
names(markers_by_celltype) <- levels(Idents(seu_cd45))
markers_all <- c(markers, unlist(markers_by_celltype, recursive=F)) 


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
  )#,
  #   reduction.model = "wnn.umap"
  # )
  
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
  loupeR::create_loupe_from_seurat(seu_tmp, output_dir=file.path(PDIR, "results", "cloupe"),
                                   output_name=paste0("CD45ST2_TReg_", seuid), force=T)
  file.copy(file.path(PDIR, "results", "cloupe", paste0("CD45ST2_TReg_", seuid, ".cloupe")),
            to = "~/xfer", overwrite = T)
  cat(paste0("xfer CD45ST2_TReg_", seuid, ".cloupe\n"))
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
  pdf("~/xfer/stat1_etreg.pdf")
  DimPlot(seu_combined, group.by='treg_anno', reduction='umap.mnn')
  dev.off()
  
  
  
  
  # 
  # #### 
  # seu <- subset(seu, ident=setdiff(unique(Idents(seu)), 'NA'))
  # seu <- SplitObject(seu, split.by = "orig.ident")
  # seu <- lapply(seu, SCTransform, vars.to.regress = "percent.mt")
  # features  <- SelectIntegrationFeatures(seu, nfeatures = 3000)
  # seu <- PrepSCTIntegration(seu, anchor.features = features)
  # anchors <- FindIntegrationAnchors(seu, normalization.method = "SCT", anchor.features = features)
  # seu_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
  #                                  k.weight = 46)
  # seu_combined <- seu_combined %>%
  #   RunPCA %>%
  #   RunUMAP(reduction = "pca", dims = 1:dim, reduction_name=paste0('umap_treg.', dim),
  #           n.neighbors=20L, min.dist=0.5)
  
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
# seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln2.rds"))
# embeds <- lapply(seus_treg, Embeddings, reduction='umap')
seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
# seus_treg$LN[['umap.anchors']] <- CreateDimReducObject(embeddings = embeds$LN, key = 'umap_', assay = 'RNA')
# seus_treg$Tumor[['umap.anchors']] <- CreateDimReducObject(embeddings = embeds$Tumor, key = 'umap_', assay = 'RNA')
# saveRDS(seus_treg, file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))





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
file.copy(file.path(PDIR, "results", "umap", "umap_CD45ST2_tumor_ln.TRegs_mnn.pdf"),
          to="~/xfer", overwrite = T)


#--- b.v) Output the cloupe file ----
if(!exists("treg_seu")) treg_seu <- readRDS(file=file.path(treg_refdir, "seu_anno.normalized.mnn.umap.rds"))
if(!exists("treg_seu")) seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
treg_refdir <- file.path(PDIR, "results", 'cloupe')
lapply(names(seus_treg), function(seuid){
  print(seuid)
  treg_seu <- seus_treg[[seuid]]
  treg_seu@meta.data <- treg_seu@meta.data[Cells(treg_seu),]
  Idents(treg_seu) <- 'orig.ident'
  treg_seu@meta.data <- treg_seu@meta.data[,c('orig.ident', 'treg_anno', 'treg_anno2')]
        #'cell_type_sort', 'cl_allCells', 'cl_annot', 'tissue_gen', 'tissue_sp', 'cl_tiss','cl_ct', )]
  treg_seu[['RNA2']] <- CreateAssayObject(counts=treg_seu@assays$RNA$counts)
  DefaultAssay(treg_seu) <- 'RNA2'
  treg_seu@reductions$umap.unintegrated <- NULL
  treg_seu@reductions$umap.anchors <- NULL
  # treg_seu <- subset(treg_seu, cells=sample(Cells(treg_seu), size=1000))
  loupeR::create_loupe_from_seurat(treg_seu, output_dir=treg_refdir,
                                   output_name=paste0("st2_treg.", seuid, ".mnn"), force=T)
  file.copy(file.path(treg_refdir, paste0("st2_treg.", seuid, ".mnn.cloupe")),
            to = "~/xfer", overwrite = T)
  cat(paste0("xfer st2_treg.", seuid, ".mnn.cloupe\n"))
})




#########################
#### 7. DEG analysis ####
outdir <- file.path(PDIR, "results", "degs", "st2_analysis")
dir.create(outdir, recursive = T, showWarnings = F)
method <- 'wilcox'
sigdig <- 4

#--- a) ST2: All cells ----
grouping <- 'st2_all_cells'
dir.create(file.path(outdir, grouping), recursive = T, showWarnings = F)
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
overwrite <- F
lapply(seus, function(i) write.table(table(i$manual_anno2, i$orig.ident), sep=",", quote=F))


# Remove b cells and relabel DC cells
seus$LN$manual_anno2[seus$LN$manual_anno2 == 'monocyte-derived DC'] <- 'DC'
seus$LN$manual_anno2[seus$LN$manual_anno == 'CD4_CCR7hi' & seus$LN$manual_anno2 == 'DC'] <- 'DCc'
seus$LN$manual_anno2[seus$LN$manual_anno2 == 'monocyte-derived macrophages'] <- 'moMacrophages'
seus$LN$manual_anno2[seus$LN$manual_anno == 'cDC1'] <- 'moDC'
seus$LN$manual_anno2[seus$LN$manual_anno2 == 'Treg cycling'] <- 'Treg'
seus$LN$manual_anno2[seus$LN$manual_anno2 == 'Treg1 cycling'] <- 'Treg1'

seus$LN <- subset(seus$LN, cells=Cells(seus$LN)[-which(seus$LN$manual_anno2 == 'B cell')])
# Remove b cells and relabel DC cells
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Treg1'] <- 'Treg'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Treg1_cycling'] <- 'Treg'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'macrophages'] <- 'Monocytes'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'monocyte-derived DC'] <- 'DC'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Cd8 cycling'] <- 'Cd8'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Treg2'] <- 'Treg1'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Treg2 cycling'] <- 'Treg1'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Cd19'] <- 'moMacrophages'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'monocyte-derived DC cycling'] <- 'moMacrophages'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'DC_Cd4Cd8pos'] <- 'moMacrophages'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'monocyte-derived macrophages'] <- 'moMacrophages'
seus$Tumor$manual_anno2[seus$Tumor$manual_anno2 == 'Treg1 cycling'] <- 'Treg'
seus$Tumor <- subset(seus$Tumor, cells=Cells(seus$Tumor)[grep("contaminat", seus$Tumor$manual_anno2, invert = T)])
pdf("~/xfer/x.pdf")
publicationDimPlot(seus$LN, grp='manual_anno2', reduction='umap.mnn')
publicationDimPlot(seus$Tumor, grp='manual_anno2', reduction='umap.mnn')
dev.off()
# 
# id <- 'LN'
# seux <- seus[[id]]@meta.data %>% 
#   dplyr::select(c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, 
#                   Batch, Sort, Tissue, Condition, Timepoint,
#                   manual_anno2)) %>%
#   rename_with(., ~gsub("manual_anno2", "anno_allCells", .))
# seux <- seus[[id]]@reductions$umap.mnn@cell.embeddings  %>%
#   as.data.frame %>% 
#   rename_with(., ~c('umap_allCells_1', 'umap_allCells_2')) %>%
#   cbind(seux, .) %>%
#   tibble::rownames_to_column(., "barcode")
# 
# seuy <- seus_treg[[id]]@meta.data %>% 
#   dplyr::select(treg_anno2) %>%
#   cbind(.,  seus_treg[[id]]@reductions$umap.mnn@cell.embeddings %>% 
#           as.data.frame %>% 
#           rename_with(., ~c('umap_Tregs_1', 'umap_Tregs_2')))  %>%
#   tibble::rownames_to_column(., "barcode")
# 
# seuxy <- left_join(seux, seuy, by='barcode')
# dir='/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/geo_submission/scrna/processed'
# write.table(seuxy, file=file.path(dir, paste0("ST2_", id, "_meta.csv")),
#             sep=",", col.names = T, row.names = F, quote = F)
# cnts <- GetAssayData(seus[[id]], assay='RNA', slot='counts')
# as.data.frame(cnts) %>% tibble::rownames_to_column(., "feature") %>%
#   write.table(., file=file.path(dir, paste0("ST2_", id, "_cnts.csv")),
#             sep=",", col.names = T, row.names = F, quote = F)



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





seus_treg.lt <- seus_treg[[3]]
seus_treg <- seus_treg[-3]
seus_treg$LN$treg_anno2[seus_treg$LN$treg_anno2 == 'cTreg'] <- 'Treg-LT_Stat1'

tregmap <- sapply(seus_treg[-3], function(treg){
  setNames(treg$treg_anno2, Cells(treg))
}) %>% unlist
names(tregmap) <- gsub("(LN|Tumor)\\.", "", names(tregmap))
seus_treg.lt@meta.data$treg_anno3 <- NA
seus_treg.lt@meta.data[names(tregmap),'treg_anno3'] <- as.character(tregmap)

seus_treg$LN <- subset(seus_treg$LN, cells=Cells(seus_treg$LN)[!seus_treg$LN$treg_anno2 %in% c('Treg_LT/NLT_like', 'eTreg_cycling', '')])
seus_treg$Tumor <- subset(seus_treg$Tumor, cells=Cells(seus_treg$Tumor)[!seus_treg$Tumor$treg_anno2 %in% c('eTreg_NLT/NLT_like', 'eTreg1_cycling', 'NA')])

seus_treg.lt[['umap']] <- seus_treg.lt[['umap.mnn']]
seus_treg.lt <- .monocle3Cluster(seus_treg.lt)

pdf("~/xfer/y.pdf")
publicationDimPlot(seus_treg$LN, grp='treg_anno2', reduction='umap.mnn')
publicationDimPlot(seus_treg$Tumor, grp='treg_anno2', reduction='umap.mnn')
publicationDimPlot(seus_treg.lt, grp='treg_anno3', reduction='umap.mnn')
publicationDimPlot(seus_treg.lt, grp='monocle3_partitions', reduction='umap.mnn')
dev.off()

## Violin plot between the two Treg clusters for LN and Tumor
seus_treg.lt$monocle3_partitions <- as.character(paste0("Treg_", seus_treg.lt$monocle3_partitions))


cols <- setNames(c('#1b9e77','#d95f02'), c('Treg_1','Treg_2'))
features <- c("S100a4", "S100a6", "Itgb1", "Ccr2", "Icos")
expr <- as.data.frame(t(GetAssayData(seus_treg.lt, assay='RNA', slot='data')[features,])) %>% 
  mutate(monocle3_partitions=seus_treg.lt$monocle3_partitions,
         Tissue = seus_treg.lt$Tissue) %>%
  pivot_longer(., cols=!c(monocle3_partitions, Tissue))
exprcnt <- as.data.frame(t(GetAssayData(seus_treg.lt, assay='RNA', slot='counts')[features,])) %>% 
  mutate(monocle3_partitions=seus_treg.lt$monocle3_partitions,
         Tissue = seus_treg.lt$Tissue) %>%
  pivot_longer(., cols=!c(monocle3_partitions, Tissue))
exprcnt$value <- log2(exprcnt$value)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



pdf("~/xfer/treg_1vs2.pdf", width = 5, height = 7)
Stacked_VlnPlot(seurat_object = seus_treg.lt, features = features, 
                colors_use = cols, x_lab_rotate = TRUE, plot_spacing = 0.3, 
                group.by='Tissue', split.by = "monocle3_partitions", plot_legend=T)
ggplot(expr, aes(x=Tissue, y=value, fill = monocle3_partitions)) + 
  geom_split_violin() +
  facet_grid(name ~., space='free', switch = "y") + 
  ylim(0,5) +
  ylab("") +
  cowplot::theme_cowplot()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")
ggplot(exprcnt, aes(x=Tissue, y=value, fill = monocle3_partitions)) + 
  geom_split_violin() +
  facet_grid(name ~., space='free', switch = "y") + 
  ylim(0,10) +
  ylab("") +
  cowplot::theme_cowplot()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")

dev.off()


## FindAllMarkers on everything and run DotPlot
grp <- 'treg_anno2'
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

pdf("~/xfer/ST2treg_LN_degUmaps.pdf")
publicationDimPlot(seus_treg$LN, grp='treg_anno2', reduction='umap.mnn')
gg_degumaps$LN # gg
dev.off()

pdf("~/xfer/ST2treg_Tumor_degUmaps.pdf")
publicationDimPlot(seus_treg$Tumor, grp='treg_anno2', reduction='umap.mnn')
gg_degumaps$Tumor # gg_degumaps
dev.off()




## Do a custom grouped analysis that uses the two-way RM ANOVA to test between groups
treg_markers_all <- lapply(names(seus_treg)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus_treg[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))
  
  outf <- file.path(outdir, grouping, paste0(seuid, ".rmaov_deg.rds"))
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
        FindMarkers.grouped(seu,meta=meta, idents=idents, assay='RNA') %>%
          tibble::rownames_to_column(., "gene") %>%
          mutate(ens=gm$SYMBOL$ENSEMBL[gene],
                 # biotype=ens2biotype_ids[ens],
                 biotype=gm$ENSEMBL$gene_biotype[ens]) %>%
          relocate(., c(gene, ens, biotype))
      }, error=function(e){NULL})
      
    })
    saveRDS(markers_all, file=outf)
  } else {
    markers_all <- readRDS(outf)
  }

  dir.create(file.path(outdir, grouping, seuid), showWarnings = F)
  for(id in names(markers_all)){
    dat <- markers_all[[id]]
    if(is.null(dat)) next()
    grp1 <- 'KO'
    grp2 <- 'WT'
    fileid <- paste0(gsub("\\/", "-", id), ".", grp1, "_vs_", grp2)
    outf <- file.path(outdir, grouping, seuid, paste0("DEG_", fileid, ".csv"))
    if(overwrite | !file.exists(outf)){
      write.table(dat, file=outf,
                  sep=",", col.names = T, row.names = F, quote = F)
    }
  }
  
  return(markers_all)
})
names(treg_markers_all) <- names(seus_treg)[1:2]

## Run ORA for Differential expressed genes along trajectory based on direction
lapply(names(treg_markers_all), function(seuid){
  treg_seuid <- treg_markers_all[[seuid]]
  
  outrds <- file.path(outdir, grouping, seuid, paste0("DEGora_all_", seuid, ".rds"))
  outf <- file.path(outdir, grouping, seuid, paste0("DEGora_all_", seuid, ".csv"))
  
  if(file.exists(outrds)){
    oras <- readRDS(outrds)
  } else {
    oras <- lapply(names(treg_seuid), function(treg_subtype_id){
      message(treg_subtype_id)
      treg_subtype <- treg_seuid[[treg_subtype_id]]
      if(is.null(treg_subtype)) return(NULL)
      sig <- treg_subtype %>% 
        dplyr::filter(padj_int < 0.1) %>%
        dplyr::filter(biotype == 'protein_coding')
      sigl <- split(sig, f=sig$KO.avg_log2FC > sig$WT.avg_log2FC)
      sigmap <- c('TRUE'='KO>WT', 'FALSE'='KO<WT')
      names(sigl) <- sigmap[names(sigl)]
      
      lapply(names(sigl), function(dir){
        res <- clusterProfiler::enricher(
          gene = pull(sigl[[dir]], gene), # A vector of your genes of interest
          pAdjustMethod = "BH", # Method to be used for multiple testing correction
          universe = Features(seu), # A vector containing your background set genes
          pvalueCutoff = 0.05,
          minGSSize = 10,
          maxGSSize = 500,
          qvalueCutoff = 0.05,
          TERM2GENE = geneset_df
        )@result %>%
          mutate(seuid = seuid,
                 treg_subtype=treg_subtype_id,
                 direction=dir) %>%
          dplyr::select(-c(Description)) %>%
          tibble::remove_rownames()
        return(res)
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
    # oras %>% dplyr::filter(p.adjust < 0.05) %>% 
    #   group_split(treg_subtype)
    saveRDS(oras, file=outf)
  }
  
  
  if(overwrite | !file.exists(outf)){
    write.table(oras, file=outf,
                sep=",", col.names = T, row.names = F, quote = F)
  }
  file.copy(outf, to="~/xfer", overwrite = T)
  message(paste0("xfer DEGora_all_", seuid, ".csv"))
})


##### Temp: GSEA to cytoscape ####
gs_map <- lapply(msig_l, function(i){
  data.frame("ID"=gsub("^(.*?)_.*$", "\\1", names(i)),
             "Desc"=gsub("^.*?_(.*)$", "\\1", names(i)),
             "setSize"=sapply(i, length)) %>%
    magrittr::set_rownames(gsub("_", "-", names(i)))
}) %>% do.call(rbind, .)
dup <- duplicated(rownames(gs_map))
if(any(dup)){
  gs_map <- gs_map[-which(dup),]
}
rownames(gs_map) <- gsub("^.*?\\.", "", rownames(gs_map)) %>% gsub("_", "-", .)

cytoscape_dat <- lapply(names(markers_all), function(comparison){
  print(comparison)
  dat <- markers_all[[comparison]]
  if(is.null(dat)) return(NULL)
  dat <- findmarkers2CytoscapeFormat(dat=dat, gs_map=gs_map)
  dat$up[is.na(dat$up)] <- 0
  dat$down[is.na(dat$down)] <- 0
  dir.create(file.path(outdir, "gse", "cytoscape",  "st2_tumor"), 
             showWarnings = F, recursive = T)
  write.table(dat$up, 
              file=file.path(outdir, "gse", "cytoscape",  "st2_tumor",
                             paste0("cytoscape_", comparison, ".up.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(dat$down, 
              file=file.path(outdir, "gse", "cytoscape", "st2_tumor", 
                             paste0("cytoscape_", comparison, ".down.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(markers_all[[comparison]] %>%
                tibble::rownames_to_column("Pathway") %>%
                dplyr::relocate(),
              file=file.path(outdir, "gse", "cytoscape", "st2_tumor",
                             paste0("pathway_deg_", comparison, ".csv")),
              sep=",", col.names = T, row.names = F)
})

lapply(names(msig_l), function(msig_id){
  CON <- file(file.path(outdir, "gse", "cytoscape", paste0("goprofiler.", msig_id, ".gmt")), "a")    #Open connection to append
  lapply(names(msig_l[[msig_id]]), function(geneset_id){
    line <- c(gsub("[_]", "-", toupper(geneset_id)),
              gsub("[_]", "-", toupper(geneset_id)) %>% gsub("^.*?:.*?[0-9]+-", "", .),
              msig_l[[msig_id]][[geneset_id]]) %>%
      as.character %>% 
      paste(., collapse="\t")
    writeLines(paste0(line), CON)
  })
  close(CON)       
})

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
    dplyr::filter(Database == databaseid)
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





#############################
#### 8. Pathway analysis ####
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



#--- b) SCPA - untested ----
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

#--- d) Loupe visualization ----
mat <- GetAssayData(cd8seu, assay='AUCell')
clusters <- cd8seu@meta.data[,c('orig.ident', 'sample', 'sample', 'cell_type', 'group')]
projections <- loupeR:::select_projections(cd8seu)
makeLoupe(mat=mat, meta=clusters, projections=projections,
          output_dir=file.path(PDIR, "results", "loupe"), 
          output_name='p14_pathway_activity')
file.copy(file.path(PDIR, "results", "loupe", "p14_pathway_activity.cloupe"), to="~/xfer", overwrite = T)

#--- e) Cytoscape visualization of differential pathways ----
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

#--- f) ST2: TRegs - suppressive Signature ----
grouping <- 'st2_tregs'
seus_treg <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.rds"))
overwrite <- F


seus_treg <- seus_treg[-3]
suppressive_sig <- read.csv(file.path(PDIR, "ref", "treg_suppressive_signature.csv"), header = F) %>% 
  magrittr::set_colnames(., "entrez_gene") %>%
  mutate(gs_name='suppressive_signature') %>%
  relocate(gs_name)

suppdf <- lapply(names(seus_treg), function(seuid){
  seu <- seus_treg[[seuid]]
  expr <- GetAssayData(seu, slot='counts')
  expr <- expr[rowSums(expr>0)>=(ncol(seu) > 0.05), ]
  
  scores <- aucellFun(msig_ds = suppressive_sig, 
                      expr_mat=expr, gm,
                      mapfrom='SYMBOL', mapto='SYMBOL')
  melt(split(assay(scores), seu$treg_anno2)) %>%
    mutate(seutype=seuid)
}) %>% do.call(rbind, .) %>% 
  dplyr::filter((L1 != '' & !is.na(L1) & L1 != 'NA'))


saveRDS(suppdf, file=file.path("~/xfer", "suppresive_scores.rds"))

suppdf <- readRDS("suppresive_scores.rds")
order <- sapply(split(suppdf, suppdf$L1), function(i){
  min(sapply(split(i$value, i$seutype), median) )
}) %>% sort
suppdf$L1 <- factor(as.character(suppdf$L1), levels=names(order))



  
pdf("~/xfer/suppresive_scores.pdf")
suppdf %>% 
  ggplot(., aes(x=L1, y=value, fill=seutype)) +
  facet_grid(seutype ~ ., space = 'free') +
  geom_violin() + 
  ggbeeswarm::geom_quasirandom(size=0.5, alpha=0.2) +
  cowplot::theme_cowplot() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust=1, vjust = 1))
dev.off()
#######################################
#### 9. SCENIC upstream regulators ####
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



#--- c) OLD Differential regulon scores ----

for(compid in (names(comp_lists))){
  print(paste0(">>>    ", compid, "    <<<"))
  spl <- comp_lists[[compid]]
  
  outdir <- file.path(PDIR, "results", "regulons", compid)
  dir.create(outdir, showWarnings = F, recursive = T)
  
  if(!file.exists(file.path(outdir, paste0("scenic_", compid, ".rds")))){
    markers <- apply(.createGrps(spl), 2, function(ids){
      print(paste0("Processing: ", paste(ids, collapse=",")))
      M <- FindMarkers(cd8seu, assay='regulons', 
                       group.by='group',
                       ident.1=ids[1], ident.2=ids[2], 
                       logfc.threshold=0) %>%
        tibble::rownames_to_column("regulon") %>%
        mutate(id.1=ids[1],
               id.2=ids[2],
               p_val = round(p_val, sigdig),
               avg_log2FC=round(avg_log2FC, sigdig),
               p_val_adj = ifelse(p_val_adj < 10^(-1*sigdig), 
                                  p_val_adj, 
                                  round(p_val_adj, sigdig))) %>%
        relocate(., c('id.1', 'id.2'), .after='regulon')
      cells <- Seurat:::IdentsToCells(object = cd8seu, 
                                      ident.1 = ids[1], 
                                      ident.2 = ids[2], 
                                      cellnames.use = Cells(cd8seu)) 
      D <- apply(GetAssayData(cd8seu, assay='regulons'), 1, function(i){
        x <- effsize::cohen.d(i[cells$cells.1], i[cells$cells.2])
        setNames(c(x$conf.int, x$estimate), c('D.lower', 'D.upper', 'D'))
      }) %>%
        t %>% as.data.frame %>%
        tibble::rownames_to_column("regulon") %>%
        mutate(D.lower = round(D.lower, sigdig),
               D.upper =round(D.upper, sigdig),
               D =round(D, sigdig))
      
      MD <- left_join(M, D, by='regulon')
      return(MD)
    })
    
    saveRDS(markers, file=file.path(outdir, paste0("scenic_", compid, ".rds")))
  } else {
    markers <- readRDS(file=file.path(outdir, paste0("scenic_", compid, ".rds")))
  }
  
  for(grpid in names(markers)){
    m <- markers[[grpid]] %>%
      mutate(regulon = gsub(",", "-", regulon))
    write.table(m, file=file.path(outdir, paste0("Pathway_", grpid, ".csv")),
                sep=",", col.names = T, row.names = F, quote = F)
  }
}

#--- d) loupe visualization ----
mat <- GetAssayData(cd8seu, assay='regulons')
clusters <- cd8seu@meta.data[,c('orig.ident', 'sample', 'sample', 'cell_type', 'group')]
projections <- loupeR:::select_projections(cd8seu)
makeLoupe(mat=mat, meta=clusters, projections=projections,
          output_dir=file.path(PDIR, "results", "loupe"), 
          output_name='p14_regulons')
file.copy(file.path(PDIR, "results", "loupe", "p14_regulons.cloupe"), to="~/xfer", overwrite = T)


#############################################################
#### 10. Integrating Differential Gene, Regulon, Pathway ####
outdir <- file.path(PDIR, "results", "gpr_integration")
dir.create(outdir, recursive=T, showWarnings = F)
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
}

outdir_pway <- file.path(PDIR, "results", "pathway")
method='aucell'
geneset_opt='msigdb'
scores_l <- readRDS(file=file.path(outdir_pway, method, paste0("aucell.", geneset_opt, ".", mode, ".rds")))
if(geneset_opt == 'gprofiler'){
  scores <- lapply(scores_l[-grep("HP|MIRNA", names(scores_l))], assay) %>% 
    do.call(rbind, .)
} else {
  scores <- lapply(scores_l, assay) %>% 
    do.call(rbind, .)
}


refined_networks <- read.csv(file.path(outdir, "saras_refined_networks.v2.csv"), 
                             header = T, sep=",")
jaccf <- file.path(PDIR, "ref", "jacc_gsmat.rds")
gosimf <- file.path(PDIR, "ref", "semantic_gsmat.rds")

if(!file.exists(gosimf)){
  goids <- c('GO:BP', 'GO:MF', 'GO:CC')
  semantic_mats <- lapply(goids, function(goid){
    message(goid)
    geneset_dfi <- geneset_l[[grep(goid, names(geneset_l))]]
    goterms <- unique(geneset_dfi$gs_exact_source)
    X = sapply(goterms, .semanticsim, y=goterms, g=GOg)
  })
  names(semantic_mats) <- goids
  
  semantic_mats2 <- lapply(semantic_mats, function(sem.mat){
    max.dist <- max(c(max(sem.mat, na.rm=T), 15))
    message(max.dist)
    sem.mat[sem.mat==0 | is.na(sem.mat)] <- max.dist
    return(sem.mat)
  })

  saveRDS(semantic_mats2, file=gosimf)
} else {
  semantic_mats <- readRDS(file=gosimf)
}

if(!file.exists(jaccf)){
  jacc_gsmat <- sapply.allbyall(with(geneset_df, split(gene_symbol, f=gs_name)), .jaccard)
  saveRDS(jacc_gsmat, file=jaccf)
} else {
  jacc_gsmat <- readRDS(file=jaccf)
}
geneset_size <- sapply(with(geneset_df, split(gene_symbol, f=gs_name)), length)
geneset_size <- geneset_size/300
# jacc_gsmat <- 1-jacc_gsmat

moi_list <- list("LN"=c('eTreg_NLT-like.B1', 'eTreg_suppressive_cycling.B1', "Treg_LT/NLT_like.B1",
                        'eTreg_NLT.B2', 'eTreg_NLT-like.B2', 'eTreg_suppressive_cycling.B2'),
                 "Tumor"=c('eTreg_NLT-like.B1', 'eTreg_cycling.B1',
                           'eTreg_NLT-like.B2', 'eTreg_NLT.B2', 'eTreg_cycling.B2'))
dirmap <- c('TRUE'='KO>WT', 'FALSE'='KO<WT')


#### LCA Test section ####
treg_markers_all <- lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))
  
  
  regulon_outf <- file.path(PDIR, "results", "regulons", "differential", paste0(seuid, ".rmaov_deregulon.rds"))
  regulon_markers_all <- readRDS(regulon_outf)
  
  regulonora_outf <- file.path(PDIR, "results", "regulons", "regulon_ora.rds")
  regulonora_markers_all <- readRDS(regulonora_outf)
  
  pathway_outf <- file.path(PDIR, "results", "pathway", 'aucell', "differential", paste0(seuid, ".rmaov_degeneset.rds"))
  pathway_markers_all <- readRDS(pathway_outf)
  
  deg_outf <- file.path(PDIR, "results", "degs", "st2_analysis", 'st2_tregs', paste0(seuid, ".rmaov_deg.rds"))
  deg_markers_all <- readRDS(deg_outf)
  
  degora_outf <- file.path(PDIR, "results", "degs", "st2_analysis", 'st2_tregs', seuid, paste0("DEGora_all_", seuid, ".rds"))
  degora_markers_all <- readRDS(degora_outf)
  
  min.max <- 2:10
  padj_threshold <- 0.05
  set.seed(1234)
  goterm <- 'GO:BP'
  verbose=F
  
  gomap <- with(geneset_df, setNames(gs_exact_source, gs_name))
  gomap.rev <- setNames(names(gomap), gomap)
  markers.of.interest <- moi_list[[seuid]]
  
  dir_mod3.integrate <- lapply(dirmap, function(dir){
    outf <- file.path(outdir, paste0("mod3gosemantic.", seuid, ".", dir, ".rds"))
    if(!file.exists(outf) | overwrite){
      message("Generating in...")
      mod3.integrate  = lapply(markers.of.interest, function(subgroup){
        print(subgroup)
        
        # if(is.null(regulon_markers_all[[subgroup]])) return(NULL)
        deg_pathways <- degora_markers_all %>%
          dplyr::filter(., direction==dir) %>% 
          dplyr::filter(., treg_subtype == subgroup) %>%
          dplyr::filter(., p.adjust < padj_threshold) %>% 
          pull(ID)
        
        regulon_sel <- regulon_markers_all[[subgroup]] %>%  
          mutate(direction = dirmap[as.character(KO.D > WT.D)]) %>%
          dplyr::filter(., padj_int < padj_threshold) %>%  
          dplyr::filter(., direction < dir) %>%  
          pull(Pathway)
        regulon_pathways <- sapply(regulon_sel, function(regulon_i){
          dplyr::filter(regulonora_markers_all[[regulon_i]]@result, p.adjust < 0.01) %>%
            pull(ID)
        }) %>% unlist %>% as.character
        
        pathways <- pathway_markers_all[[subgroup]] %>% 
          mutate(direction = dirmap[as.character(KO.D > WT.D)]) %>%
          dplyr::filter(., padj_int < padj_threshold)  %>% 
          dplyr::filter(., direction==dir) %>% 
          pull(Pathway)
        
        mod.pathways <- list("gene"=deg_pathways,
                             "regulon"=regulon_pathways,
                             "pathway"=pathways)
        
        goterm.simple <- gsub(":", "", goterm)
        mod.pathways_go <- lapply(mod.pathways, grep, pattern=goterm.simple, value=T)
        mod.pathways_goid <- lapply(mod.pathways_go, function(i) {
          go.i <- gomap[gsub("-", "_", i)]
          go.i[go.i %in% V(GOg)$name]
        })
        
        
        
        sem.mat <- semantic_mats[[goterm]]
        sem.mat <- abs(1-round((sem.mat-1)/(max(sem.mat)-1), 2)) * 100
        set.seed(1234)
        mod.min.parents <- lapply(mod.pathways_goid, function(goid.mod){
          message(paste0("goid.mod length: ", length(goid.mod)))
          if(length(goid.mod) == 0) return(NULL)
          pamk.split <- .fconn_networks(sem.mat[goid.mod,goid.mod])
          # pamk.split <- pamk.split[sapply(pamk.split, length)>1]
          
          idx <- 1
          min.parents <- lapply(pamk.split, function(i){
            # print(idx)
            idx <<- idx + 1
            if(length(i) == 1){
              lca_id <- data.frame("parent"=i, "sub"=i)
            } else {
              lca_id <- .getLca(GOg, i, M=sem.mat, verbose=verbose, get.sublvls=T, level=3)
            }
            return(lca_id)
          })
          min.parents <- min.parents[!sapply(min.parents, is.null)]
          names(min.parents) <- sapply(min.parents, function(i) unique(i$parent))
          # GOTERMS[names(min.parents)]
          # cbind(GOTERMS[min.parents], GOTERMS[as.character(min.parents.bkup)])
          return(min.parents)
        })
        # x <- list2df(lapply(mod.min.parents, function(i) do.call(rbind, i)))
        # xdf <- x %>% 
        #   mutate(dir = dir, 
        #          parentGO=GOTERMS[parent],
        #          subGO=GOTERMS[sub]) 
        # write.table(xdf, file=file.path("~/xfer", "xdf.csv"), sep=",", col.names = T,
        #             row.names = F, quote = F)
        
        
        parentids <- lapply(mod.min.parents, function(i) na.omit(names(i)))
        
        aggregate_go <- semantic.combine(golist=mod.pathways_goid, g=GOg, M=sem.mat,
                                         min.cnt=1, get.sublvls=T)
        # aggregate_go <- semantic.combine(golist=parentids, g=GOg, M=sem.mat, 
        #                                  min.cnt=1, get.sublvls=T)
        # subids <- lapply(mod.min.parents, function(i) unique(na.omit(do.call(rbind, i)$sub)))
        # aggregate_go.sub <- semantic.combine(golist=subids, g=GOg, min.cnt=1, get.sublvls=T)
        if(class(aggregate_go) == 'list') aggregate_go <- do.call(rbind, aggregate_go)
        # GOTERMS[unlist(aggregate_go)]
        # GOTERMS[aggregate_go[,1]]
        return(aggregate_go)
      })
      names(mod3.integrate) <- markers.of.interest
      saveRDS(mod3.integrate, file=outf)
    } else {
      message("Reading in...")
      mod3.integrate <- readRDS(file=outf)
    }
    return(mod3.integrate)
  })
  df <- list2df(dir_mod3.integrate) %>% 
    mutate(Var.1 = dirmap[Var.1], 
           parentGO=GOTERMS[parent],
           subGO=GOTERMS[sub]) 
  write.table(df, file=file.path("~/xfer", paste0(seuid, ".GO_semantic_int.csv")),
              sep=",", col.names = T, row.names = F)
  paste0("xfer ", seuid, ".GO_semantic_int.csv")
  
  GOTERMS[semantic.combine(dir_mod3.integrate[[1]], g=GOg)]
  GOTERMS[semantic.combine(dir_mod3.integrate[[2]], g=GOg)]
  GOTERMS[semantic.combine(unlist(dir_mod3.integrate), g=GOg)]
  
  #### Test make barplots ####
  ##  Check for interferon gamma TYPE II INTERFERON
  ## INTERLEUKIN 6, 12, 17
  
  cbind(sapply(godat[[1]], function(i) gomap.rev[i]),
        sapply(godat[[1]], function(i) GOTERMS[i]))
  
  godat <- dir_mod3.integrate[[1]]
  goterms <- na.omit(as.character(sapply(godat[[1]], function(i) gomap.rev[i])))
  
  seu <- seus[[seuid]]
  expr <- GetAssayData(seu, slot='counts')
  expr <- expr[rowSums(expr)>=50, ]
  geneset_i <- geneset_df %>% 
    dplyr::filter(gs_name %in% goterms)
  scores <- aucellFun(msig_ds = geneset_i, 
                      expr_mat=expr, gm,
                      mapfrom='SYMBOL', mapto='SYMBOL')
        
  
  seu[['AUCell']] <- CreateAssayObject(data = assay(scores))
  DefaultAssay(seu) <- 'AUCell'
  seu$group <- paste0(seu$treg_anno2, ".", seu$Batch)
  Idents(seu) <- 'group'
  seusub <- subset(seu, ident=markers.of.interest)
  
  
  
    
  Idents(seu) <- 'treg_anno'
  seu$group <- paste0(seu$Condition, ".", gsub("[37]d", "Trx", seu$Timepoint))
  seu$group <- factor(as.character(seu$group), levels=names(colors_list))
  pt.size <- 0
  ggdat <- Stacked_VlnPlot2(
    seurat_object = seu,
    features = Features(seu), 
    split.by='group',
    x_lab_rotate = TRUE,
    colors_use = colors_list,
    pt.size=pt.size
  )
  seuidx <- 1
  ggdat <- lapply(ggdat, function(ggobj){
    ggobj <- ggobj + scale_fill_manual(values=colors_list) +
      theme(plot.margin = margin(t = 1, b = 1)) + 
      ylim(c(0,0.2))
    if(seuidx != 1) ggobj <- ggobj  + theme(axis.text.y=element_blank(),
                                            axis.title.y=element_blank())
    return(ggobj)
  })
  ggdat[[1]] <- ggdat[[1]] + ggtitle(seuid) + theme(plot.title=element_text())
  
  add.legend=T
  pdf("~/xfer/x.pdf", height = 10, width=20)
  wrp <- patchwork::wrap_plots(plotlist = ggdat, ncol = 1)
  if(add.legend){
    wrp <- wrp & theme(legend.position = "right")
    wrp <- wrp + patchwork::plot_layout(guides = 'collect')
  }
  wrp
  dev.off()
  #### END - make barplots ####
})


#### End test ####

treg_markers_all <- lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))
  
  
  regulon_outf <- file.path(PDIR, "results", "regulons", "differential", paste0(seuid, ".rmaov_deregulon.rds"))
  regulon_markers_all <- readRDS(regulon_outf)
  
  regulonora_outf <- file.path(PDIR, "results", "regulons", "regulon_ora.rds")
  regulonora_markers_all <- readRDS(regulonora_outf)
  
  pathway_outf <- file.path(PDIR, "results", "pathway", 'aucell', "differential", paste0(seuid, ".rmaov_degeneset.rds"))
  pathway_markers_all <- readRDS(pathway_outf)
  
  deg_outf <- file.path(PDIR, "results", "degs", "st2_analysis", 'st2_tregs', paste0(seuid, ".rmaov_deg.rds"))
  deg_markers_all <- readRDS(deg_outf)
  
  degora_outf <- file.path(PDIR, "results", "degs", "st2_analysis", 'st2_tregs', seuid, paste0("DEGora_all_", seuid, ".rds"))
  degora_markers_all <- readRDS(degora_outf)
  
  # subgroup <- 'Treg_LT/NLT_like.B2'
  min.match <- 2
  min.max <- 2:10
  padj_threshold <- 0.05
  pattern.match <- '^GOBP|^REACTOME|^HALLMARK'
  set.seed(1234)
  
  
  dir_mod3.integrate <- lapply(dirmap, function(dir){
    markers.of.interest <- moi_list[[seuid]]
    outf <- file.path(outdir, paste0("mod3integrate.", seuid, ".", dir, ".rds"))
    if(!file.exists(outf)){
      message("Generating in...")
      mod3.integrate  = lapply(markers.of.interest, function(subgroup){
        print(subgroup)
        
        if(is.null(regulon_markers_all[[subgroup]])) return(NULL)
        deg_pathways <- degora_markers_all %>%
          dplyr::filter(., direction==dir) %>% 
          dplyr::filter(., treg_subtype == subgroup) %>%
          dplyr::filter(., p.adjust < padj_threshold) %>% 
          pull(ID)
        
        regulon_sel <- regulon_markers_all[[subgroup]] %>%  
          mutate(direction = dirmap[as.character(KO.D > WT.D)]) %>%
          dplyr::filter(., padj_int < padj_threshold) %>%  
          dplyr::filter(., direction < dir) %>%  
          pull(Pathway)
        regulon_pathways <- sapply(regulon_sel, function(regulon_i){
          dplyr::filter(regulonora_markers_all[[regulon_i]]@result, p.adjust < 0.01) %>%
            pull(ID)
        }) %>% unlist %>% as.character
        
        pathways <- pathway_markers_all[[subgroup]] %>% 
          mutate(direction = dirmap[as.character(KO.D > WT.D)]) %>%
          dplyr::filter(., padj_int < padj_threshold)  %>% 
          dplyr::filter(., direction==dir) %>% 
          pull(Pathway)
        
        mod.pathways <- list("gene"=deg_pathways,
                             "regulon"=regulon_pathways,
                             "pathway"=pathways)
        
        mod.networks <- lapply(mod.pathways, function(uids){
          uids <- gsub("-", "_", uids)
          uids <- grep(pattern.match, uids, value=T)
          message(paste0(length(uids), "; ", min(min.max), " - ", max(min.max)))
          tryCatch({
            pamk.res <- fpc::pamk(jacc_gsmat[uids,uids], krange=min.max)
            pamk.spl <- .reducePamkClusters(pamk.res, mat=jacc_gsmat, min.max=min.max)
            return(pamk.spl)
          }, error=function(e){NULL})
        })
        if(sum(sapply(mod.networks, is.null))>=2) return(NULL)
        
        pathway_universe <- unique(as.character(unlist(mod.networks)))
        mod.networks.unl <- unlist(mod.networks, recursive=F)
        lapply(mod.networks.unl, function(i) grep("t.?helper", i, value=T, ignore.case = T))
        
        phyper_mat <- sapply.allbyall(mod.networks.unl, .avgjaccard.genes, 
                                      reflist=with(geneset_df, split(gene_symbol, f=gs_name)))
        pamk.res <- fpc::pamk(phyper_mat, krange=min.max)
        pamk.split <- .reducePamkClusters(pamk.res, max.clus.size=10, 
                                          mat=phyper_mat, min.max=min.max)
        
        pamk.split.idx <- sapply(pamk.split, function(i) {
          sum(c(as.integer(any(grepl("gene", i))),
                as.integer(any(grepl("regulon", i))),
                as.integer(any(grepl("pathway", i))))) >= 2
        })
        pamk.genesets <- lapply(pamk.split[pamk.split.idx], function(i){
          (mod.networks.unl[i])
        })
        names(pamk.genesets) <- make.unique(rep("NoN", length(pamk.genesets)))
        # mod3.integrate[[subgroup]] <- pamk.genesets
        return(pamk.genesets)
      }) %>% 
        setNames(., markers.of.interest)
      saveRDS(mod3.integrate, file=outf)
    } else {
      message("Reading in...")
      mod3.integrate <- readRDS(outf)
    }
    return(mod3.integrate)
  })
  return(dir_mod3.integrate)
})
names(treg_markers_all) <- names(seus)[1:2]
mod3.df <- list2df(treg_markers_all) %>% 
  magrittr::set_colnames(c('gpr_network', 'pathway', 'non', 
                           'celltype.batch', 'direction', 'tissue')) %>% 
  mutate(direction=dirmap[as.character(direction)])

write.table(mod3.df, file=file.path(outdir, "integrated.network_of_networks.csv"),
            sep=",", col.names = T, row.names = F, quote = F)
file.copy(file.path(outdir, "integrated.network_of_networks.csv"), to="~/xfer",
          overwrite = T)
    

ids <- list(list('KO>WT', 'eTreg_NLT-like.B1', 'LN', paste0("NoN.", c(11, 16, 37, 53))),
            list('KO>WT', 'Treg_LT/NLT_like.B1', 'LN', paste0("NoN.", c(8,9))),
            list('KO>WT', 'eTreg_NLT-like.B2', 'LN', paste0("NoN.", c(8,17,37,67))),
            list('KO>WT', 'eTreg_NLT.B2', 'LN', c('NoN', paste0("NoN.", c(2,4,13,26,52,68,74,78,107)))),
            list('KO>WT', 'eTreg_suppressive cycling.B2', 'LN', paste0("NoN.", c(8,29,28,72))),
            list('KO>WT', 'eTreg_NLT-like.B1', 'Tumor', paste0("NoN.", c(4,20,23,27,55,100))),
            list('KO>WT', 'eTreg_NLT-like.B2', 'Tumor', paste0("NoN.", c(9,20,43,55))),
            list('KO>WT', 'eTreg_NLT.B2', 'Tumor', c('NoN', paste0("NoN.", c(28,46,27,29,60,13,34)))),
            list('KO>WT', 'eTreg_cycling.B2', 'Tumor', paste0("NoN.", c(13,17))))
x <- lapply(ids, function(i){
  mod3.df %>% 
    dplyr::filter(direction == i[[1]]) %>% 
    dplyr::filter(celltype.batch == i[[2]]) %>% 
    dplyr::filter(tissue == i[[3]]) %>% 
    dplyr::filter(non %in% i[[4]]) %>%
    mutate(uid = paste0(non,tissue,celltype.batch))
}) %>% do.call(rbind, .)
term_doc <- xtabs(~ uid + pathway, data=x, sparse = TRUE)
co_occur <- crossprod(term_doc, term_doc)
diag(co_occur) <- 0
pamk.res <- fpc::pamk(co_occur, krange=c(2:20))

     
     


#--- x. Make Pathway networks ----
networks.to.plot <- refined_networks %>% 
  dplyr::select(-c(gpr_network, pathway)) %>%
  dplyr::filter(!duplicated(.)) %>% 
  collapse::rsplit(.,
                   ~ tissue + celltype.batch + direction)
# networks.to.label <- list('Tumor'=list(
#   'eTreg_NLT-like.B1'=list(
#     'TRUE'=list('NoN.4'=c('REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE')),
#     'FALSE'=list('NoN'=c('REACTOME_TRANSLATION'))
#   )
# ))
dir.colors <- c('TRUE'='red', 'FALSE'='blue')
   
network_plots <- lapply(names(networks.to.plot), function(seuid){
  lapply(names(networks.to.plot[[seuid]]), function(group){
    lapply(names(networks.to.plot[[seuid]][[group]]), function(direction){
      dir <- setNames(names(dirmap), dirmap)[direction]
      networks <- networks.to.plot[[seuid]][[group]][[direction]]
      lapply(networks$NETWORK, function(network_id){
        print(paste(c(seuid, group, dir, network_id), collapse="_"))
        network_i <- treg_markers_all[[seuid]][[dir]][[group]][[network_id]]
        searchterms <- list('NETWORK'=network_id,
                            'celltype.batch'=group,
                            'tissue'=seuid,
                            'direction'=as.character(direction))
        refineddat <- .searchReference(refined_networks, searchterms)
        network_i <- lapply(split(network_i, gsub("\\..*", "", names(network_i))), function(i){
          i <- unique(as.character(unlist(i)))
          j <- i[which(i %in% unique(as.character(refineddat$pathway)))]
          return(j)
        })
        # label <- networks.to.label[[seuid]][[group]][[dir]][[network_id]]
        
        # Make igraph network
        genes.in.network <- as.character(unlist(network_i))
        uids <- unique(genes.in.network)
        matraw <- mat <- jacc_gsmat[uids,uids]
        pnormval <- 2
        connect.nodes <- TRUE
        while(connect.nodes){
          pnormval <- pnormval - 0.01
          mat <- matraw
          mat[mat<quantile(mat, pnorm(pnormval))] <- 0
          g <- graph_from_adjacency_matrix(mat, weighted=T, mode="undirected", diag=F)
          
          # check if 90% of vertices are connected
          connected.vertices <- (sum(degree(g) == 0) <= min(c(3, ceiling(length(uids)*0.1))))
          # check how dense the connections are
          median.degree <- quantile(degree(g), 0.25) >= 2
          # check the components
          if(components(g)$no == 1) connect.nodes <- FALSE
          # if(connected.vertices & median.degree & evratio) connect.nodes <- FALSE
          if(pnormval <= 0) connect.nodes <- FALSE
        }

        # Set colors
        gene.cnts <- table(genes.in.network)[uids]
        gene.alpha <- (1-round(gene.cnts * (1/3),2)) * 1.15
        gene.alpha <- 0
        my_color <- t_col(rep(dir.colors[dir], length(uids)), gene.alpha)
        # # set labels
        # idx <- grepl(label, uids)
        # uids.label <- rep("", length(uids))
        # uids.label[which(idx)] <- uids[idx]
        
        g$size <-  0.5+(geneset_size[uids]*3)
        ggp <- GGally::ggnet2(g,
                              # label=uids.label,
                              label=uids,
                              node.size = g$size, 
                              max_size=4.5,
                              node.color = my_color, 
                              edge.size = E(g)$weight*2) +
          labs(title=paste0(seuid, ": ", group),
               subtitle=network_id) +
          theme_classic() + theme(legend.position = "none")
        return(ggp)
      })
    })
  })
})

pdf("~/xfer/networks.refined.labelled.pdf")
network_plots
dev.off()




#--- x. Make DEG heatmap ----
# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation

deg_ggs <- lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))

  deg_outf <- file.path(PDIR, "results", "degs", "st2_analysis", 'st2_tregs', paste0(seuid, ".rmaov_deg.rds"))
  deg_markers_all <- readRDS(deg_outf)
  
  gene_sel <- sample(deg_markers_all[[1]]$gene, size=20)
  gene_hms <- lapply(names(deg_markers_all), function(treg_subtype_id){
    idx <- match(treg_subtype_id, names(deg_markers_all))
    is.last <- ifelse(idx == length(deg_markers_all), TRUE, FALSE) 
    message(idx)
    
    # treg_subtype_id <- names(deg_markers_all)[10]  #### TESTING DATA
    treg_subtype <- deg_markers_all[[treg_subtype_id]]
    if(is.null(treg_subtype)) return(NULL)
    sigmap <- c('TRUE'='KO>WT', 'FALSE'='KO<WT')
    sig <- treg_subtype %>% 
      dplyr::filter(padj_int < 0.1) %>%
      dplyr::filter(biotype == 'protein_coding')
    sig$direction <- sigmap[as.character(with(sig, KO.avg_log2FC > WT.avg_log2FC))]
    sig <- sig[order(sig$KO.avg_log2FC, decreasing = T),]
    sigmat <- sig %>%
      dplyr::select(c(gene, KO.avg_log2FC, WT.avg_log2FC)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames('gene') %>%
      rename_with(., ~gsub("\\..*", "", .))

    
    ha = rowAnnotation(foo = anno_mark(at = which(rownames(sigmat) %in% gene_sel), 
                                       labels = gene_sel[which(gene_sel %in% rownames(sigmat))],
                                       side='left'))
    gg <- Heatmap(as.matrix(sigmat), name = "mat", cluster_rows = FALSE,
            left_annotation = ha, row_names_gp = gpar(fontsize = 4),
            show_row_names = FALSE, #row_names_side = "left",
            show_column_dend=FALSE,
            show_column_names = is.last,
            show_heatmap_legend = is.last,
            row_split=factor(as.character(sig$direction), levels=c('KO>WT', 'KO<WT')),
            width=unit(10, "cm"),
            row_title_rot = 0,
            height=ncol(sigmat)*unit(10, "mm"),
            row_title = treg_subtype_id) %>% 
      draw() %>% 
      grid.grabExpr()
    return(gg)
  })
  gene_hms <- gene_hms[which(!sapply(gene_hms, is.null))]
  
  patchgg <- patchwork::wrap_plots(gene_hms, ncol = 1)
  return(patchgg)
})
  
pdf("~/xfer/degs_all.pdf", height=20)
deg_ggs
dev.off()


#--- x. Make Regulon heatmap ----
regulon_ggs <- lapply(names(seus)[1:2], function(seuid){
  print(paste0(">> ", seuid))
  seu <- seus[[seuid]]
  seu$group <- paste0(seu$orig.ident, ".", seu$treg_anno)
  uids <- unique(seu$group)
  
  df <- data.frame(IDs=uids,
                   celltype=gsub("^.*\\.", "", uids),
                   batch=gsub("_.*", "", uids)) %>% 
    mutate(Treatment = gsub("^.*(7d|3d|Un).*", "\\1", uids),
           Condition = gsub("^.*(KO|WT).*", "\\1", uids))
  meta_l <- split(df, f=list(df$celltype, df$batch))
  
  regulon_outf <- file.path(PDIR, "results", "regulons", "differential", paste0(seuid, ".rmaov_deregulon.rds"))
  regulon_markers_all <- readRDS(regulon_outf)
  regulon_sel <- sample(regulon_markers_all[[1]]$Pathway, size=5)
  
  regulon_hms <- lapply(names(regulon_markers_all), function(treg_subtype_id){
    idx <- match(treg_subtype_id, names(regulon_markers_all))
    is.last <- ifelse(idx == length(regulon_markers_all), TRUE, FALSE) 
    message(idx)
    
    # treg_subtype_id <- names(regulon_markers_all)[10]  #### TESTING DATA
    treg_subtype <- regulon_markers_all[[treg_subtype_id]]
    if(is.null(treg_subtype)) return(NULL)
    if(min(treg_subtype$padj_int) >= 0.1) return(NULL)
    sigmap <- c('TRUE'='KO>WT', 'FALSE'='KO<WT')
    sig <- treg_subtype %>% 
      dplyr::filter(padj_int < 0.1)
    sig <- sig[order(sig$KO.D, decreasing = T),]
    sig$direction <- sigmap[as.character(with(sig, KO.D > WT.D))]
    sigmat <- sig %>%
      dplyr::select(c(Pathway, KO.D, WT.D)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames('Pathway')
    
    ha = rowAnnotation(foo = anno_mark(at = which(rownames(sigmat) %in% regulon_sel), 
                                       labels = regulon_sel[which(regulon_sel %in% rownames(sigmat))],
                                       side='right'))
    gg <- Heatmap(as.matrix(sigmat), name = "mat", cluster_rows = FALSE,
                  right_annotation = ha, row_names_gp = gpar(fontsize = 4),
                  show_row_names = FALSE, #row_names_side = "left",
                  show_column_dend=FALSE,
                  show_column_names = is.last,
                  show_heatmap_legend = is.last,
                  row_split=factor(as.character(sig$direction), levels=c('KO>WT', 'KO<WT')),
                  width=unit(10, "cm"),
                  row_title_rot = 0,
                  height=ncol(sigmat)*unit(10, "mm"),
                  row_title = treg_subtype_id) %>% 
      draw() %>% 
      grid.grabExpr()
    return(gg)
  })
  regulon_hms <- regulon_hms[which(!sapply(regulon_hms, is.null))]
  
  patchgg <- patchwork::wrap_plots(regulon_hms, ncol = 1)
  return(patchgg)
})

pdf("~/xfer/regulons_all.pdf", height=20)
regulon_ggs
dev.off()


  
#############################################
#### 11. Velocity and Trajectory analysis - inprogres ####
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
#   
#   # Normalize and cluster
#   # workflow from https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.md
#   set.seed(seed)
#   seu_velo_small[["RNA"]] <- seu_velo_small[["spliced"]]
#   seu_velo_small <- seu_velo_small %>%
#     SCTransform(.) %>%
#     RunPCA %>%
#     RunUMAP(., dims=1:30) %>%
#     FindNeighbors(., dims = 1:30) %>%
#     FindClusters(.)
#   DefaultAssay(seu_velo_small) <- "RNA"
#   seu_velo_small <- JoinLayers(seu_velo_small)
#   seu_velo_small[['umap_orig']] <- CreateDimReducObject(embeddings = Embeddings(seu_small, reduction='umap'),
#                                                         key = 'UMAP_', assay = 'RNA')
#   seu_velo_small$functional_cluster <- as.character(seu_small$treg_anno)
#   
#   print("saving...")
#   # flatfile method: https://smorabit.github.io/tutorials/8_velocyto/
#   id <- paste0(seu_velo_id, ".", mode)
#   writeSeuToFlat(seu_velo_small, assay='RNA',reduction = 'umap_orig',
#                  out_metaf=file.path(outdir, "scVelo", paste0(id, ".metadata.csv")),
#                  out_cntsf=file.path(outdir, "scVelo", paste0(id, ".counts.mtx")),
#                  out_pcaf=file.path(outdir, "scVelo", paste0(id, ".pca.csv")),
#                  out_featuref=file.path(outdir, "scVelo", paste0(id, ".genes.csv")))
#   
#   # zellkonverter method: https://github.com/satijalab/seurat/discussions/7402#discussioncomment-7918919
#   sce_obj <- as.SingleCellExperiment(seu_velo_small, assay = c("RNA"))
#   library(zellkonverter)
#   zellkonverter::writeH5AD(sce_obj, file.path(outdir, "scVelo", paste0(id, ".sce_obj.h5ad")), 
#                            X_name = 'counts')
#   
#   # Old method doesnt work
#   # SeuratDisk::SaveH5Seurat(seu_velo_small, 
#   #                          filename = file.path(outdir, "scVelo", 
#   #                                               paste0("seu_velo.", seu_velo_id, ".", path_j, ".h5Seurat")),
#   #                          overwrite=T)
#   # cwd <- getwd()
#   # setwd(file.path(outdir, "scVelo"))
#   # SeuratDisk::Convert(source = paste0("seu_velo.", seu_velo_id, ".", path_j, ".h5Seurat"), 
#   #                     dest = "h5ad", overwrite=T)
#   # setwd(cwd)
#   return(NULL)
# })



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


#--- b) TRegs BGPLVM analysis ----
if(!exists("seul"))  seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))

celltype <- c('cd8'="^T_CD8")
celltype <- c('tregs'="^TReg")

prepSubData <- function(seurat_obj, genes = NULL){
  sub_data = seurat_obj
  DefaultAssay(sub_data) <- 'RNA'
  sub_data <- NormalizeData(object = sub_data, 
                            scale.factor = 10000, display.progress = F) #%>% 
    # Seurat::FindVariableFeatures(., display.progress = F, 
    #                              num.bin = 100,
    #                              binning.method = "equal_frequency")
  
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
    seu$manual_anno <- seu$treg_anno
    Idents(seu) <- 'manual_anno'
    seu <- subset(seu, ident=setdiff(unique(grep(celltype, Idents(seu), value=T, ignore.case=T)),
                                     unique(grep("remove|cycling", Idents(seu), value=T, ignore.case = T))))
    
    # Subset for a specific sample
    Idents(seu) <- 'orig.ident'
    seu <- prepSubData(subset(seu, ident=sampleid))
    write.csv(seu$mat, file = paste0("./results/bgplvm/", sampleid, ".varfeat.csv"),
              row.names = T, quote = F, col.names = T)
    return(seu)
  }) %>% setNames(., sampleids)
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
      dplyr::select(manual_anno, lv_lvl)  %>%
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

seul_celltype2 <- lapply(seul_celltype, function(seu){
  merge(seu[[1]]$seurat, lapply(seu[-1], function(i) i$seurat))
})

lvdf <- seul_celltype2$LN@meta.data[,'manual_anno', drop=F] %>% 
  tibble::rownames_to_column(., "cell") %>% 
  left_join(., lv, by='cell')
logexp <- as.matrix(seul_celltype2$LN@assays$RNA@data[,lv$cell])

avg_grp_lv0 <- sapply(split(lvdf$LV0, lvdf$manual_anno), mean, na.rm=T)
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
  dplyr::select(-c(geneID, Description)) %>%
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
    dplyr::filter(p.adjust < 0.001) %>%
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

#################################
#### 12. Violin Treg Subsets ####
require(Seurat)
require(tidyverse)
require(patchwork)
require(scCustomize)
source("~/git/mini_projects/mini_functions/singlecell/Stacked_VlnPlot2.R")

mode <- 'treg' #'treg', 'cd8'
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.diet.rds"))
# seus2 <- lapply(seus, function(seu){
#   seu <- DietSeurat(seu, layers=c('data', 'counts'))
#   seu@meta.data <- seu@meta.data[,c('orig.ident', 'Condition', 'Timepoint', 'treg_anno')]
#   return(seu)
# })
# saveRDS(seus2, file=file.path("~/xfer", "CD45ST2_TReg_tumor_ln.diet.rds"))
outdir <- file.path(PDIR, "results", "scVelo", mode)
dir.create(outdir, recursive=T, showWarnings = F)

gene_list_plot <- c('Ccr8', 'Cxcr4')
pt.size <- 0.15
add.legend <- TRUE
colors_list <- c('KO.Trx'='#b2182b', 'KO.Un'='#ef8a62', 'WT.Trx'='#4d4d4d', 'WT.Un'='#999999')
seugg <- lapply(names(seus), function(seuid){
  seuidx <- match(seuid, names(seus))
  seu <- seus[[seuid]]
  Idents(seu) <- 'treg_anno'
  seu$group <- paste0(seu$Condition, ".", gsub("[37]d", "Trx", seu$Timepoint))
  seu$group <- factor(as.character(seu$group), levels=names(colors_list))
  
  ggdat <- Stacked_VlnPlot2(
    seurat_object = seu,
    features = gene_list_plot, 
    split.by='group',
    x_lab_rotate = TRUE,
    colors_use = colors_list,
    pt.size=pt.size
    )
  ggdat <- lapply(ggdat, function(ggobj){
    ggobj <- ggobj + scale_fill_manual(values=colors_list) +
      theme(plot.margin = margin(t = 1, b = 1)) + 
      ylim(c(0,5))
    if(seuidx != 1) ggobj <- ggobj  + theme(axis.text.y=element_blank(),
                                            axis.title.y=element_blank())
    return(ggobj)
  })
  ggdat[[1]] <- ggdat[[1]] + ggtitle(seuid) + theme(plot.title=element_text())
  
    
  return(ggdat)
})


wrp.ln <- wrap_plots(plotlist = seugg[[1]], ncol = 1)
wrp.tumor <- wrap_plots(plotlist = seugg[[2]], ncol = 1)
if(add.legend){
  wrp.tumor <- wrp.tumor & theme(legend.position = "right")
  wrp.tumor <- wrp.tumor + plot_layout(guides = 'collect')
}
pdf("~/xfer/x.pdf")
wrp.ln | wrp.tumor 
dev.off()









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
