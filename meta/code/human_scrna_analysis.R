## Human Analysis
## sara Tumor/LN KO and WT samples
renv::load("/cluster/home/quever/downloads/renvs/")
# Core
library(monocle3)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
# Annotation
library(pagoda2)

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))
sym2biotype_ids <- with(GTF, setNames(gene_biotype, gene_name))
split.datasets <- T


###################
#### Functions ####
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
# source("~/git/mini_projects/mini_functions/singlecell/publicationDimPlot.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
source("~/git/mini_projects/mini_functions/singlecell/aucell_helper.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")
gm <- geneMap(species='Homo sapiens')

#### Other Functions ####
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



#### Params ####
visualize <- FALSE
humanpdir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/meta_analysis/human_scrna/'
pathway_f <- file.path("ref", "pathways.csv")
ortholog_map <- '/cluster/projects/mcgahalab/ref/maps/jax_human_mouse_genemap.csv'
setwd(humanpdir)


#### Main ####
#--- Preprocess and integrate datasets ----
dir.create(file.path("datasets", "processed"), showWarnings = F)
datasets <- list(
  "GSE131907_luad_ln"=list(
    "data"=file.path("datasets", "GSE131907_luad_ln", "GSE131907_Lung_Cancer_raw_UMI_matrix.rds"),
    "meta"=file.path("datasets", "GSE131907_luad_ln", "GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
    "samples"=c(paste0("EBUS_", c("10","12", "13", "15", "19", "51")), "BRONCHO_11", 
                paste0("LN_0", c(1:8)), "LN_11", "LN_12"),
    "out"=file.path("datasets", "processed", "GSE131907_luad_ln.raw.rds")
  ),
  "GSE212797_hnscc_ln"=list(
    "data"=file.path("datasets", "GSE212797_hnscc_ln", "GSE212797_adata.h5ad"),
    "meta"=file.path("datasets", "GSE212797_hnscc_ln", "GSE212797_IPIHNSC118_comb_TCR_filtered_contig_annotations.csv.gz"),
    "samples"=c('117_lymph node', '118_lymph node', '126_lymph node', 'multi_lymph node'),
    "out"=file.path("datasets", "processed", "GSE212797_hnscc_ln.raw.rds")
  )
)

sampleid <- 'GSE212797_hnscc_ln'
.processGSE212797 <- function(datasets){
  cat("Processing GSE212797_hnscc_ln dataset\n")
  ## Read in the h5ad data 
  ad <- reticulate::import("anndata", convert = FALSE)
  anndata <- ad$read_h5ad(datasets[[sampleid]]$data)
  
  ## Parse the counts, metadata and reduction from the data
  cnts <- reticulate::py_to_r(anndata$raw$X) %>%
    as.data.frame %>% 
    magrittr::set_rownames(as.vector(anndata$obs_names$values)) %>%
    magrittr::set_colnames(as.vector(anndata$var_names$values)) %>%
    t %>% as.data.frame
  meta <- reticulate::py_to_r(anndata$obs) 
  
  # The raw data appears to be library-size normalized and log2(x+1) scaled, reversing it 
  # to get raw counts
  log2ScaledToCounts <- function(cnt, meta){
    ncount <- meta$n_counts
    rawcnt <- exp(cnt)-1
    scalefactor <- ncount / colSums(rawcnt)
    rawcnt2 <- rawcnt
    for(idx in seq_along(colnames(rawcnt))){
      rawcnt2[,idx] <- round(rawcnt[,idx] * scalefactor[idx], 0)
    }
    return(rawcnt2)
  }
  rawcnts <- log2ScaledToCounts(cnts, meta)
  rawcnts <- rawcnts %>% 
    magrittr::set_colnames(as.vector(anndata$obs_names$values)) %>%
    magrittr::set_rownames(as.vector(anndata$raw$var_names$values)) %>%
    as.data.frame
  
  meta <- meta %>%
    filter(sample_name %in% datasets[[sampleid]]$samples) %>%
    mutate(sample_name = as.character(sample_name))
  umap <- reticulate::py_to_r(anndata$obsm[['X_umap']]) %>%
    as.matrix %>% 
    magrittr::set_rownames(as.vector(anndata$obs_names$values)) %>%
    magrittr::set_colnames(c('UMAP_1', 'UMAP_2')) 
  
  ## Assemble and save the seurat object
  seu <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(rawcnts[,colnames(rawcnts) %in% rownames(meta)]),sparse = T),
                            meta.data=meta)
  seu[["umap_orig"]] <- CreateDimReducObject(embeddings = umap[rownames(umap) %in% rownames(meta),], 
                                             key = "UMAP_", 
                                             assay = DefaultAssay(seu))
  saveRDS(seu, file=datasets[[sampleid]]$out)
  return(seu)
}
.processGSE131907 <- function(datasets){
  cat("Processing GSE131907_luad_ln dataset\n")
  ## Read in the count matrix and metadtaa
  dat <- readRDS(datasets[[sampleid]]$data)
  meta <- read.table(datasets[[sampleid]]$meta, sep="\t", header = T, stringsAsFactors = F) %>%
    tibble::column_to_rownames("Index") %>%
    filter(Sample %in% datasets[[sampleid]]$samples)
  dat <- dat[,colnames(dat) %in% rownames(meta)]
  
  ## Assemble the seurat object
  seu <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(dat),sparse = T),
                            meta.data=meta)
  saveRDS(seu, file=datasets[[sampleid]]$out)
  return(seu)
}

sampleids <- c('GSE212797_hnscc_ln', 'GSE131907_luad_ln')
seul <- lapply(setNames(sampleids,sampleids), function(sampleid){
  if(!file.exists(datasets[[sampleid]]$out)){
    switch(sampleid,
           GSE131907_luad_ln=.processGSE131907(datasets),
           GSE212797_hnscc_ln=.processGSE212797(datasets))
  } else {
    readRDS(file=datasets[[sampleid]]$out)
  }
})


seul <- lapply(setNames(names(seul), names(seul)), function(id){
  seul[[id]]$dataset <- id
  if(id=='GSE212797_hnscc_ln'){
    samplecol <- 'sample_name'
    annocol <- 'cell_type'
  } else if(id=='GSE131907_luad_ln'){
    samplecol <- 'Sample'
    annocol <- 'Cell_subtype'
  }
  seul[[id]]$orig.ident <- gsub(" ", ".", as.character(seul[[id]]@meta.data[,samplecol]))
  seul[[id]]$functional_clusters <- gsub(" ", ".", as.character(seul[[id]]@meta.data[,annocol]))
  
  return(seul[[id]])
})

if(!split.datasets){
  seul <- list(merge(seul[[1]], seul[-1], add.cell.ids=sampleids))
}


seuls_integ <- lapply(seul, function(seu){
  seu <- tryCatch({
    JoinLayers(seu)
  }, error=function(e){seu})
  
  DefaultAssay(seu) <- 'RNA'
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
  
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures(.)  %>%
    ScaleData(.)  %>%
    RunPCA(.) %>% 
    FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
    FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
    RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
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
    group.by='orig.ident'
  )
  seu_integ <- seu_integ %>% 
    FindNeighbors(., reduction = "integrated.mnn", dims = 1:30) %>%
    FindClusters(., resolution = 0.9, cluster.name = "mnn_clusters") %>% 
    RunUMAP(., reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
  
  seu_integ <- JoinLayers(seu_integ)
  
  return(seu_integ)
})



saveRDS(seu, file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.rds"))
saveRDS(seuls_integ, file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.splitdataset.rds"))
# seu <- readRDS(file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.rds"))

if(visualize){
  pdf("~/xfer/x.pdf", width = 9, height = 9)
  print(publicationDimPlot(seu_integ, grp='cell_type', reduction='umap.mnn', pt.size=2))
  print(publicationDimPlot(seu_integ, grp='Cell_type.refined', reduction='umap.mnn', pt.size=2))
  print(publicationDimPlot(seu_integ, grp='Cell_subtype', reduction='umap.mnn', pt.size=2))
  print(publicationDimPlot(seu_integ, grp='dataset', reduction='umap.mnn', pt.size=2))
  print(publicationDimPlot(seu_integ, grp='mnn_clusters', reduction='umap.mnn', pt.size=2))
  dev.off()
}



# seu_integ <- IntegrateLayers(
#   object = seu, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = T
# )


#--- a) Pathway annotation ----
dir.create(file.path(humanpdir, "results", "pathways"), recursive = F, showWarnings = F)
# seu <- readRDS(file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.rds"))
seul <- readRDS(file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.splitdataset.rds"))
ortho_map <- read.csv(ortholog_map, header=T, sep=",")
ortho_map <- with(ortho_map, setNames(Symbol.Human, Symbol.Mouse))
pagodards <- ifelse(split.datasets, 
                    "pagoda2.custom_activitiy.splitdataset.rds", 
                    "pagoda2.custom_activitiy.rds")
aucellrds <- ifelse(split.datasets, 
                    "aucell.custom_activitiy.splitdataset.rds", 
                    "aucell.custom_activitiy.rds")

pathways <- read.csv(pathway_f, header=T)
msig_l <- as.list(pathways) %>%
  lapply(., function(i) {
    ortho_map[i[which(i != '')]]
  })
msig_ds <- data.frame(gs_name=gsub("[0-9]*$", "", names(unlist(msig_l))) %>%
                        gsub("\\..*$", "", .),
                      entrez_gene=as.character(unlist(msig_l))) %>%
  filter(entrez_gene != '')


aucell_scores_l <- lapply(seul, function(seu){
  print(unique(seu$orig.ident))
  DefaultAssay(seu) <- 'RNA'
  expr <- seu[['RNA']]$counts
  expr <- expr[rowSums(expr)>=50, ]
  
  aucell_scores_l <- aucellFun(msig_ds = msig_ds, 
                               expr_mat=expr, gm,
                               mapfrom='SYMBOL', mapto='SYMBOL')
})
saveRDS(aucell_scores_l, file=file.path(humanpdir, "results", "pathways", aucellrds))

pagoda2_scores_l <- lapply(seul, function(seu){
  DefaultAssay(seu) <- 'RNA'
  expr <- seu[['RNA']]$counts
  expr <- expr[rowSums(expr)>=50, ]
  pagoda2_scores <- cal_pagoda2(expr, msig_l, max_gset_size=1500)
  return(pagoda2_scores)
})
saveRDS(pagoda2_scores_l, file=file.path(humanpdir, "results", "pathways", pagodards))
                                         

aucell_scores_l <- readRDS(file=file.path(humanpdir, "results", "pathways", aucellrds))
pagoda2_scores_l <- readRDS(file=file.path(humanpdir, "results", "pathways", pagodards))

lapply(names(seul), function(seuid){
  seu <- seul[[seuid]]  
  pagoda2_scores <- scores_l[[seuid]]
  rownames(pagoda2_scores) <- gsub("_", ".", rownames(pagoda2_scores))
  
  # Add assay and scale features
  pagoda2_scores[is.na(pagoda2_scores)] <- 0
  seu[['pagoda2']] <- CreateAssayObject(data=pagoda2_scores)
  DefaultAssay(seu) <- 'pagoda2'
  seu <- ScaleData(seu, features=gsub("_", ".", names(msig_l)))
  
  
  # Re-annotate the cells based on their clusters/cell-annotation overlap
  tbl <- table(seu$seurat_clusters, seu$functional_clusters)
  clus_map <- apply(tbl, 1, function(i) round(i/sum(i), 2)) %>%
    apply(., 2, function(i){
      if(max(i) < 0.5){
        idx <- head(order(i, decreasing = T), 2)
      } else {
        idx <- which.max(i)
      }
      paste(colnames(tbl)[idx], collapse="/")
    })
  seu@meta.data$aggregate_clusters <- clus_map[as.character(seu$seurat_clusters)]
  seu@meta.data$aggregate_clusters <- seu$functional_clusters
  
  Idents(seu) <- 'aggregate_clusters'
  subset_ids <- grep("Mac", unique(seu$aggregate_clusters), ignore.case = T, value = T)
  # seumac <- subset(seu, ident=c(subset_ids))
  
  if(seuid == 'GSE131907_luad_ln'){
    stages <- c("EBUS_10"="IV","BRONCHO_11"="IV","EBUS_12"="IV","EBUS_13"="IV","EBUS_15"="IIIA",
                "EBUS_19"="IV","EBUS_51"="IV","LN_01"="IIB","LN_02"="IB","LN_03"="IIB","LN_04"="IA2",
                "LN_05"="IA3","LN_06"="IA3","LN_07"="IA","LN_08"="IB","LN_11"="IB","LN_12"="IA")
    seu@meta.data$stages <- stages[as.character(seu$orig.ident)]
    
    comparisons <- list('nLN.mLN'=c('Sample_Origin', 'nLN', 'mLN'),
                        'stages'=c('stages', 'kruskal-wallis'))
  } else if(seuid == 'GSE212797_hnscc_ln'){
    comparisons <- list('KO.D7.trx_vs_un'=c('B2_LN_KO_7d', 'B2_LN_KO_Un'),
                        'WT.D7.trx_vs_un'=c('B2_LN_WT_7d', 'B2_LN_WT_Un'),
                        'KO.D3.trx_vs_un'=c('B1_LN_KO_3d', 'B1_LN_KO_Un'),
                        'WT.D3.trx_vs_un'=c('B1_LN_WT_3d', 'B1_LN_WT_Un'),
                        'D7.trx.wt_vs_ko'=c('B2_LN_WT_7d', 'B2_LN_KO_7d'),
                        'D3.trx.wt_vs_ko'=c('B1_LN_WT_3d', 'B1_LN_KO_3d'))
  }
  
  roughcounts <- seu@assays$pagoda2$data + abs(min(seu@assays$pagoda2$data))
  roughcounts <- round(roughcounts*10, 0)
  seu@assays$pagoda2$counts <- as(roughcounts, 'dgCMatrix')
  seu2loupe <- seu
  seu2loupe@meta.data <- if(seuid=='GSE131907_luad_ln'){
    seu2loupe@meta.data[,c('orig.ident', 'Sample_Origin', 'Cell_subtype', 'Cell_type.refined', 'stages')]
  } else if (seuid == 'GSE212797_hnscc_ln'){
    seu2loupe@meta.data[,c('orig.ident', 'patient_tissue', 'cell_type', 'functional_clusters')]
  }
  dir.create(file.path(humanpdir, "results", "cloupe"), showWarnings = F)
  create_loupe_from_seurat(seu2loupe, output_dir=file.path(humanpdir, "results", "cloupe"),
                           output_name=paste0(seuid, ".pagoda2"), force=T)
  file.copy(file.path(humanpdir, "results", "cloupe", paste0(seuid, ".pagoda2.cloupe")), to = "~/xfer", overwrite = T)
  cat(paste0("xfer ", seuid, ".pagoda2.cloupe\n"))
             
  seu2loupe[['RNA2']] <- CreateAssayObject(counts=seu2loupe@assays$RNA$counts)
  DefaultAssay(seu2loupe) <- 'RNA2'
  create_loupe_from_seurat(seu2loupe, output_dir=file.path(humanpdir, "results", "cloupe"),
                           output_name=paste0(seuid, ".rna"), force=T)
  file.copy(file.path(humanpdir, "results", "cloupe", paste0(seuid, ".rna.cloupe")), to = "~/xfer", overwrite = T)
  cat(paste0("xfer ", seuid, ".rna.cloupe\n")) 
  
  
  Seurat:::FindMarkers.default
  degs <- lapply(subset_ids, function(sub_i){
    seu_mac_i <- subset(seu, ident=sub_i)
    lapply(comparisons, function(comp_j){
      Idents(seu_mac_i) <- comp_j[1]
      deg <- tryCatch({
        FindMarkers(seu_mac_i, 
                    ident.1=comp_j[2],
                    ident.2=comp_j[3],
                    logfc.threshold=0)
      }, error=function(e){NULL})
      return(deg)
    })
  })
  
  
})

#--- b) Identify genes used in AUCell scoring ----
dir.create(file.path(humanpdir, "results", "pathways"), recursive = F, showWarnings = F)
# seu <- readRDS(file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.rds"))
seul <- readRDS(file=file.path(humanpdir, "results", "seuratobj", "seu_integrated.splitdataset.rds"))
ortho_map <- read.csv(ortholog_map, header=T, sep=",")
ortho_map <- with(ortho_map, setNames(Symbol.Human, Symbol.Mouse))
pagodards <- ifelse(split.datasets, 
                    "pagoda2.custom_activitiy.splitdataset.rds", 
                    "pagoda2.custom_activitiy.rds")
aucellrds <- ifelse(split.datasets, 
                    "aucell.custom_activitiy.splitdataset.rds", 
                    "aucell.custom_activitiy.rds")

pathways <- read.csv(pathway_f, header=T)
msig_l <- as.list(pathways) %>%
  lapply(., function(i) {
    ortho_map[i[which(i != '')]]
  })
msig_ds <- data.frame(gs_name=gsub("[0-9]*$", "", names(unlist(msig_l))) %>%
                        gsub("\\..*$", "", .),
                      entrez_gene=as.character(unlist(msig_l))) %>%
  filter(entrez_gene != '')

aucell_genes_l <- lapply(seul, function(seu){
  print(unique(seu$orig.ident))
  DefaultAssay(seu) <- 'RNA'
  expr <- seu[['RNA']]$counts
  expr <- expr[rowSums(expr)>=50, ]
  ####  Test Area
  mapped_id <- msig_ds$entrez_gene
  exprMat <- expr
  geneSets =  lapply(split(mapped_id, msig_ds$gs_name), function(i){
    i[!is.na(i)]
  })
  cnt_mat <- getGenesUsedInAucell(exprMat, geneSets)
  return(cnt_mat)
})
for(id in names(aucell_genes_l)){
  write.table(aucell_genes_l[[id]], 
              file=file.path(humanpdir, "results", "pathways", paste0("aucell_gene_count.", id, ".csv")),
              sep=",", quote = F, col.names = T, row.names = F)
  file.copy(file.path(humanpdir, "results", "pathways", paste0("aucell_gene_count.", id, ".csv")),
            to="~/xfer", overwrite = T)
  cat(paste0("xfer aucell_gene_count.", id, ".csv\n"))
}

#--- Comparison 2 ----
pathway_caller <- 'aucell'
aucell_scores_l <- readRDS(file=file.path(humanpdir, "results", "pathways", aucellrds))
scores_l <- lapply(aucell_scores_l, assay)
make_loupe <- FALSE
min.pct <- 0.1

group_pcts_l <- lapply(names(seul), function(seuid){
  seu <- seul[[seuid]]  
  scores <- scores_l[[seuid]]
  rownames(scores) <- gsub("_", ".", rownames(scores))
  
  # Add assay and scale features
  scores[is.na(scores)] <- 0
  seu[[pathway_caller]] <- CreateAssayObject(data=scores)
  DefaultAssay(seu) <- pathway_caller
  seu <- ScaleData(seu, features=gsub("_", ".", names(msig_l)))
  
  if(seuid == 'GSE131907_luad_ln'){
    Idents(seu) <- 'Cell_subtype'
    subset_ids <- grep("(mo.Mac|mono)", (unique(seu$Cell_subtype)), ignore.case = T, value = T)
    # stages <- c("EBUS_10"="IV","BRONCHO_11"="IV","EBUS_12"="IV","EBUS_13"="IV","EBUS_15"="IIIA",
    #             "EBUS_19"="IV","EBUS_51"="IV","LN_01"="IIB","LN_02"="IB","LN_03"="IIB","LN_04"="IA2",
    #             "LN_05"="IA3","LN_06"="IA3","LN_07"="IA","LN_08"="IB","LN_11"="IB","LN_12"="IA")
    # seu@meta.data$stages <- stages[as.character(seu$orig.ident)]
    
    comparisons <- list('nLN.mLN'=c('Sample_Origin', 'nLN', 'mLN'))
  } else if(seuid == 'GSE212797_hnscc_ln'){
    Idents(seu) <- 'cell_type'
    subset_ids <- grep("Mac", unique(seu$cell_type), ignore.case = T, value = T)
    comparisons <- list('117_ln'=c('sample_name', '117_lymph node', NULL),
                        '118_ln'=c('sample_name', '118_lymph node', NULL),
                        '126_ln'=c('sample_name', '126_lymph node', NULL),
                        'multi_ln'=c('sample_name', 'multi_lymph node', NULL))
  }
  
  ## Check how many genes in the geneset are expressed in greater than 10% of each group
  group_pcts <- lapply(subset_ids, function(sub_i){
    seu_mac_i <- subset(seu, ident=sub_i)
    lapply(comparisons, function(comp_j){
      Idents(seu_mac_i) <- comp_j[1]
      group_pct <- lapply(msig_l, function(i) {
        DefaultAssay(seu_mac_i) <- 'RNA'
        fc <- FoldChange(seu_mac_i, 
                         features=intersect(i, Features(seu_mac_i)), 
                         ident.1=comp_j[2],
                         ident.2=if(is.na(comp_j[3])) NULL else  comp_j[3])
        fc_all <- FoldChange(seu_mac_i, 
                             features=Features(seu_mac_i), 
                             ident.1=comp_j[2],
                             ident.2=if(is.na(comp_j[3])) NULL else  comp_j[3])
        fc1 <- (fc$`pct.1` > min.pct)
        fc2 <- (fc$`pct.2` > min.pct)
        fcall <- (fc_all$`pct.1` > min.pct) | (fc_all$`pct.2` > min.pct) 
        list("group1_genes"=rownames(fc)[fc1],
             "group2_genes"=rownames(fc)[fc2],
             "expr"=c('group1'=table(fc1)['TRUE']/length(fc1),
                      'group2'=table(fc2)['TRUE']/length(fc2),
                      "background"=table(fcall)['TRUE']/nrow(seu_mac_i)))
      })
      genes <- sapply(group_pct, function(i) c(i$group1_genes, i$group2_genes)) %>%
        unlist 
      
      group_pct_df <- lapply(group_pct, function(i) i$expr) %>% do.call(rbind, .)
      
      # return(list("group"=group_pct_df, "genes"=genes))
      return(group_pct_df)
    })
  }) %>%
    setNames(., subset_ids)
  return(group_pcts)
}) %>%
  setNames(., names(seul))

degs_all <- lapply(names(seul), function(seuid){
  seu <- seul[[seuid]]  
  scores <- scores_l[[seuid]]
  rownames(scores) <- gsub("_", ".", rownames(scores))
  
  # Add assay and scale features
  scores[is.na(scores)] <- 0
  seu[[pathway_caller]] <- CreateAssayObject(data=scores)
  DefaultAssay(seu) <- pathway_caller
  seu <- ScaleData(seu, features=gsub("_", ".", names(msig_l)))
  
  ## Make loupe file
  if(make_loupe){
    roughcounts <- seu@assays[[pathway_caller]]$data + abs(min(seu@assays[[pathway_caller]]$data))
    roughcounts <- round(roughcounts*(100/max(roughcounts)), 0)
    seu@assays[[pathway_caller]]$counts <- as(roughcounts, 'dgCMatrix')
    seu2loupe <- seu
    seu2loupe@meta.data <- if(seuid=='GSE131907_luad_ln'){
      seu2loupe@meta.data[,c('orig.ident', 'Sample_Origin', 'Cell_subtype', 'Cell_type.refined')]
    } else if (seuid == 'GSE212797_hnscc_ln'){
      seu2loupe@meta.data[,c('orig.ident', 'patient_tissue', 'cell_type', 'functional_clusters')]
    }
    dir.create(file.path(humanpdir, "results", "cloupe"), showWarnings = F)
    loupeR::create_loupe_from_seurat(seu2loupe, output_dir=file.path(humanpdir, "results", "cloupe"),
                                     output_name=paste0(seuid, ".", pathway_caller), force=T)
    file.copy(file.path(humanpdir, "results", "cloupe", paste0(seuid, ".", pathway_caller, ".cloupe")), to = "~/xfer", overwrite = T)
    cat(paste0("xfer ", seuid, ".", pathway_caller, ".cloupe\n"))
  }
  # seumac <- subset(seu, ident=c(subset_ids))
  
  if(seuid == 'GSE131907_luad_ln'){
    Idents(seu) <- 'Cell_subtype'
    subset_ids <- grep("(mo.Mac|mono)", (unique(seu$Cell_subtype)), ignore.case = T, value = T)
    # stages <- c("EBUS_10"="IV","BRONCHO_11"="IV","EBUS_12"="IV","EBUS_13"="IV","EBUS_15"="IIIA",
    #             "EBUS_19"="IV","EBUS_51"="IV","LN_01"="IIB","LN_02"="IB","LN_03"="IIB","LN_04"="IA2",
    #             "LN_05"="IA3","LN_06"="IA3","LN_07"="IA","LN_08"="IB","LN_11"="IB","LN_12"="IA")
    # seu@meta.data$stages <- stages[as.character(seu$orig.ident)]
    
    comparisons <- list('nLN.mLN'=c('Sample_Origin', 'nLN', 'mLN'))
  } else if(seuid == 'GSE212797_hnscc_ln'){
    Idents(seu) <- 'cell_type'
    subset_ids <- grep("Mac", unique(seu$cell_type), ignore.case = T, value = T)
    comparisons <- list('117_ln'=c('sample_name', '117_lymph node', NULL),
                        '118_ln'=c('sample_name', '118_lymph node', NULL),
                        '126_ln'=c('sample_name', '126_lymph node', NULL),
                        'multi_ln'=c('sample_name', 'multi_lymph node', NULL))
  }
  
  ## Differential between the two groups for the pathway
  degs <- lapply(subset_ids, function(sub_i){
    seu_mac_i <- subset(seu, ident=sub_i)
    lapply(comparisons, function(comp_j){
      Idents(seu_mac_i) <- comp_j[1]
      
      DefaultAssay(seu_mac_i) <- pathway_caller
      deg <- tryCatch({
        cells <- Seurat:::IdentsToCells(object = seu_mac_i,  ident.1=comp_j[2],
                               ident.2=if(is.na(comp_j[3])) NULL else  comp_j[3], 
                               cellnames.use = Cells(seu_mac_i))
        ass <- GetAssayData(object = seu_mac_i, slot = "data")
        data.1 <- log(rowMeans(ass[Features(seu_mac_i), cells$`cells.1`]), base=2)
        data.2 <- log(rowMeans(ass[Features(seu_mac_i), cells$`cells.2`]), base=2)
        
        FindMarkers(seu_mac_i, 
                    ident.1=comp_j[2],
                    ident.2=if(is.na(comp_j[3])) NULL else  comp_j[3],
                    logfc.threshold=0) %>%
          mutate(avg_log2FC=(data.1 - data.2),
                 dataset=seuid,
                 celltype=sub_i,
                 ident.1=comp_j[2],
                 ident.2=if(is.na(comp_j[3])) NA else comp_j[3]) %>%
          tibble::rownames_to_column("pathway") %>%
          relocate(c("pathway", "dataset", "celltype", "ident.1", "ident.2"))
      }, error=function(e){NULL})
      return(deg)
    })
  })
  names(degs) <- subset_ids
  return(degs)
})
degs_all <- do.call(rbind, unlist(unlist(degs_all, recursive = F), recursive=F)) %>%
  tibble::remove_rownames(.)

outf <- file.path(humanpdir, "results", "pathways", "diffexp_macMono.csv")
write.table(degs_all, file=outf,
            sep=",", col.names = T, row.names = F, quote = F)
file.copy(outf, to="~/xfer", overwrite = T)
cat(paste0("xfer ", basename(outf), "\n"))