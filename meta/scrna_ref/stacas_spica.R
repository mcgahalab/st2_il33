library(ProjecTILs)
library(Seurat)
library(dplyr)
library(STACAS)
library(ggplot2)

source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/geneMap.R")


###########################
#### Global Parameters ####
doublet_quantile_cutoff <- 0.95
visualize_qc <- FALSE

mt_miqc <- FALSE # whether to use miQC for percent.mt filter, else use isOutlier
call_doublets <- FALSE # whether to call doublets

# Annotation
atlas_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils'
scgate_refdir <- file.path(atlas_refdir, 'scgate_model')
projectils_refdir <- file.path(atlas_refdir, "ProjecTILs")
projectils_dir <- file.path(atlas_refdir, "spica")

gm <- geneMap(species='Homo sapiens')
##############
#### Main ####
# ---- Read and integrate SPICA datasets ----
## Read in pre-made annotated reference datasets from SPICA
dir.create(file.path(projectils_dir, 'custom'), showWarnings = F)
refobj_dc <- orthogene.scgate(refobj_dc, from='human', 'mouse')
refobj_dc <- readRDS(file.path(projectils_dir, 'DC_human_ref_v1.rds'))
refobj_cd4 <- readRDS(file.path(projectils_dir, 'CD4T_human_ref_v1.rds'))
refobj_cd8 <- readRDS(file.path(projectils_dir, 'CD8T_human_ref_v1.rds'))

refobj_lcmv <- readRDS(file.path(projectils_dir, 'ref_LCMV_Atlas_mouse_v1.rds'))
refobj_momac <- readRDS(file.path(projectils_dir, "APC_atlas_v1_SPICA.rds"))

## Integrate all the SPICA reference datasets together
refobj_int <- Run.STACAS(list("dc"=refobj_dc, "cd4"=refobj_cd4, "cd8"=refobj_cd8),
                         anchor.features=3000) # , "momac"=refobj_momac)
refobj_int <- Run.STACAS(list("lcmv"=refobj_lcmv, 'momac'=refobj_momac),
                         anchor.features=3000) # , "momac"=refobj_momac)

## Run PCA, UMAP and make ProjecTILs reference
ndim = 30
seed = 1234
refobj_int <- RunPCA(refobj_int, features = VariableFeatures(refobj_int), ndims.print = 1:5,
                     nfeatures.print = 5, npcs = ndim)
refobj_int <- RunUMAP(refobj_int, reduction = "pca", dims = 1:ndim, seed.use = seed)
ref.refobj_int <- make.reference(ref = refobj_int, ndim = ndim, seed = seed, recalculate.umap = TRUE,
                                 annotation.column = "functional.cluster")

saveRDS(ref.refobj_int, file=file.path(projectils_dir, 'custom', 'dc.cd4.cd8.rds'))

# ---- Read and integrate custom datasets ----
# Reads in a metadata file containing the GSE, GSM, and dataset label
# to be used when annotating the datasets downstream
ref_seu_dir <- file.path(atlas_refdir, "gse_dat", "teresa_MDSC")
ref_celltypes <- read.csv(file.path(ref_seu_dir, "celltypes.csv"), 
                          check.names = F, stringsAsFactors = F) %>%
  tibble::column_to_rownames("GSM") %>% 
  filter(grepl("(Mouse_ILC2|Mouse_Bcell)", ID))
id_map <- c('GSE168158_1'='GSM5130034',
            'GSE168158_2'='GSM51300342')

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/teresa_ly6/scrna'
dir='/cluster/projects/mcgahalab/ref/scrna/projectils/gse_dat/Treg_and_Tmem_scRNA'
outrds_int2 <- file.path(PDIR, "data", "seurat_obj", 'ref1b_filtered.rds')
outrds_int <- file.path(PDIR, "data", "seurat_obj", "ref2_dubfiltered.intermediate.rds")
if(!exists('seus_post')) seus_post <- readRDS(outrds_int)
if(!exists('seus_preprocessed')) seus_preprocessed <- readRDS(outrds_int2)
if(!exists("seu_treg")) seu_treg <- readRDS(file.path(dir, "10X_all_data_preQC.reprocessed.rdata"))

names(seus_preprocessed) <- id_map[names(seus_preprocessed)]
seus_post <- seus_post[intersect(rownames(ref_celltypes), names(seus_post))]
seus_preprocessed <- seus_preprocessed[intersect(rownames(ref_celltypes), names(seus_preprocessed))]
seu_treg <- subsetSeu(seu_treg, colid='cl_annot')

custom_int <- lapply(c(seus_post, seus_preprocessed, seu_treg), function(seu_i){
  DefaultAssay(seu_i) <- 'RNA'
  return(seu_i)
}) %>%
  Run.STACAS(., anchor.features=3000)
custom_int$functional.cluster <- ref_celltypes[custom_int$orig.ident,]$ID
ident_idx <- which(!is.na(custom_int$ident))
custom_int$functional.cluster[ident_idx] <- custom_int$ident[ident_idx]
ident_idx <- which(is.na(custom_int$functional.cluster))
custom_int$functional.cluster[ident_idx] <- custom_int$cl_annot[ident_idx]

## Read in SPICA mouse data and integrate with custom
refobj_tit <- readRDS(file.path(projectils_dir, "ref_TILAtlas_mouse_v1.rds")) # tumor-infiltrating T cells
# refobj_lcmv <- readRDS(file.path(projectils_dir, 'ref_LCMV_Atlas_mouse_v1.rds')) # virus-specific CD8 T cells
refobj_momac <- readRDS(file.path(projectils_dir, "APC_atlas_v1_SPICA.rds")) # mononuclear phagocytes - MoMacDC
custom_all_int <- Run.STACAS(list(refobj_tit, refobj_momac, custom_int),  #refobj_lcmv, 
                             anchor.features=3000)

## Run PCA, UMAP and make ProjecTILs reference
## Run PCA, UMAP and make ProjecTILs reference
ndim = 30
seed = 1234
custom_int <- RunPCA(custom_int, features = VariableFeatures(custom_int), ndims.print = 1:5,
                     nfeatures.print = 5, npcs = ndim)
custom_int <- RunUMAP(custom_int, reduction = "pca", dims = 1:ndim, seed.use = seed)
ref.custom_int <- make.reference(ref = custom_int, ndim = ndim, seed = seed, recalculate.umap = TRUE,
                                 annotation.column = "functional.cluster")

saveRDS(ref.custom_int, file=file.path(projectils_dir, 'custom', 'custom_ilc2.treg.rds'))
saveRDS(ref.custom_all_int, file=file.path(projectils_dir, 'custom', 'custom_ilc2.mouse_atlasB.rds'))

pdf("~/xfer/test7.pdf")
DimPlot(custom_int, group.by = "functional.cluster", label = T, repel = T, label.size = 4, raster=T) +
  theme(aspect.ratio = 1) + NoLegend()
DimPlot(refobj_tit, group.by = "functional.cluster", label = T, repel = T, label.size = 4, raster=T) +
  theme(aspect.ratio = 1) + NoLegend()
DimPlot(refobj_momac, group.by = "functional.cluster", label = T, repel = T, label.size = 4, raster=T) +
  theme(aspect.ratio = 1) + NoLegend()
DimPlot(ref.custom_all_int, group.by = "functional.cluster", label = T, repel = T, label.size = 4, raster=T) +
  theme(aspect.ratio = 1) + NoLegend()
dev.off()

# ---- Tregs ----
dir='/cluster/projects/mcgahalab/ref/scrna/projectils/gse_dat/Treg_and_Tmem_scRNA'
treg_meta <- read.csv(file.path(dir, "all_data_cl_metadata.csv"), header = T) %>%
  tibble::column_to_rownames(., "X")

load(file.path(dir, "10X_all_data_preQC.rdata")) #all_data
rawdata <- attributes(all_data)$raw.data
symbol_ids <- gm$ENSEMBL$SYMBOL[rownames(rawdata)]
rm_idx <- which(duplicated(symbol_ids) | is.na(symbol_ids))
rawdata <- rawdata[-rm_idx,]
rownames(rawdata) <- symbol_ids[-rm_idx]
idx <- which(colnames(rawdata) %in% rownames(treg_meta))

seu_treg <- CreateSeuratObject(rawdata[,idx], project = "Treg", assay = "RNA",
                               min.cells = 0, min.features = 0, names.field = 1,
                               names.delim = "_", meta.data = treg_meta)
Idents(seu_treg) <- 'tissue_sp'
seu_treg <- subset(seu_treg, idents=c('bLN', 'mLN', 'skin'))
seu_treg <- preprocessSeu(seu_treg, ncount_min=0, ncount_max=Inf,
                          nfeature_min=0, nfeature_max=Inf,
                          mt_max=100, org='mouse', numpcs=30, getPCs=FALSE,
                          variable_features=2000, res=1.2)
saveRDS(seu_treg, file=file.path(dir, "10X_all_data_preQC.reprocessed.rdata"))

pdf("~/xfer/treg_dimplot.pdf", width = 12)
DimPlot(seu_treg, group.by='seurat_clusters', raster=T, label=T)
DimPlot(seu_treg, group.by='cl_annot', raster=T, label=T)
DimPlot(seu_treg, group.by='cl_allCells', raster=T, label=T)
dev.off()


seu_treg <- all_data



# ---- x ----
corexpr <- function(obj, genes, clusterid='seurat_clusters', treatmentid='orig.ident'){
  DefaultAssay(obj) <- 'RNA'
  dat <- GetAssayData(obj, slot='data')
  
  # Iterate through all possible gene pairs
  cordat <- apply(combn(genes, 2), 2, function(gene_pair){
    genes_dat <- t(as.matrix(dat[gene_pair,]))
    
    
    
    # Iterate through all possible clusters
    clusters <- unique(obj@meta.data[,clusterid])
    cl_cor <- lapply(clusters, function(cl){
      clus_idx <- which(obj@meta.data[[clusterid]] == cl)

      # Iterate through all treatment conditions
      ids <- unique(obj@meta.data[,treatmentid])
      sapply(ids, function(id_i){
        id_idx <- which(obj@meta.data[[treatmentid]] == id_i)
        cells_idx <- Cells(obj)[intersect(clus_idx, id_idx)]
        genes_dat_cl_id <- genes_dat[cells_idx,]
        
        r <- tryCatch({
          cor.test(genes_dat_cl_id[,1], genes_dat_cl_id[,2])
        }, error=function(e){
          list("estimate"=NA, "p.value"=NA)
        })
        c('cluster'=cl, 'id'=id_i, 'g1'=gene_pair[1], 'g2'=gene_pair[2], 
          'r'=r$estimate, 'p'=r$p.value,  'n'=length(cells_idx))
      }) %>% t %>% as.data.frame
    })
    names(cl_cor) <- clusters
    return(cl_cor)
  })
  names(cordat) <- apply(combn(genes, 2), 2, paste, collapse="_")
  return(cordat)
}
corexpr(seu, genes=c("Mtor","Pten","Ahr" ))


# ---- New ----
seu$lntumor <- gsub("_.*", "", seu$orig.ident)
Idents(seu) <- 'lntumor'
seus <- lapply(c('t'='Tumor', 'ln'='LN'), function(id) subset(seu, ident=id))
seus <- lapply(seus, function(seu_i){
  preprocessSeu(seu_i, ncount_min=0, nfeature_min=0, nfeature_max=Inf,
                mt_max=100, numpcs=30, res=0.8)
})

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
library(SeuratWrappers)
library(STACAS)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
# Annotation
library(SingleR)
library(ProjecTILs)
# QC 
library(scater)
library(DoubletFinder)

###############
#### Setup ####

visualize <- FALSE
seed <- 1234
set.seed(seed)


outdir <- file.path(PDIR, "results")
groups <- list.files(file.path(datadir, 'filtered_feature_bc_matrices'))

# Create ENZ -> SYMBOL mapping key
# genome_gse <- org.Mm.eg.db
# txby <- keys(genome_gse, 'SYMBOL')

# Read in gtf file to map ENS->Biotype
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf_file)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

###########################
#### Global Parameters ####
doublet_quantile_cutoff <- 0.95
visualize_qc <- FALSE

mt_miqc <- FALSE # whether to use miQC for percent.mt filter, else use isOutlier
call_doublets <- FALSE # whether to call doublets

# Preprocessed datasets:
#     Must contain a '*counts.csv.gz' file and '*metadata.csv.gz' file
#     in the /cluster/projects/mcgahalab/ref/scrna/projectils/gse_dat/GSEXXXX/
#     directory
pp_datasets <- c('GSE216166', 'GSE117131') #, 'GSE184423')

# Annotation
atlas_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils'
scgate_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/scgate_model'
projectils_refdir <- file.path(atlas_refdir, "ProjecTILs")

###################
#### Functions ####
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/STACAS_fix.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")



#######################################
#### 1.a Make ProjecTILs reference ####
# This section of the code reads in all the 10x cellranger
# files and loads them into seurat objects. It saves an
# RDS object which contains a list of seurat objects, one
# for each dataset
ref_groups <- list.files(ref_seu_dir, pattern = "^GSM")
outrds <- file.path(datadir, "seurat_obj", "ref1_unfiltered.rds")
relabel_id <- gsub("^GSM[0-9]*_", "", ref_groups) %>%
  setNames(., ref_groups)

if(file.exists(outrds)){
  cat("Reading in existing 'seus' object...\n")
  seus <- readRDS(file=outrds)
} else {
  seus <- NULL
}

seus <- lapply(ref_groups, function(grp){
  if(relabel_id[grp] %in% names(seus)){
    cat(paste0("Existing ID: ", grp, "\n"))
    seu <- seus[[relabel_id[grp]]]
  } else {
    gsm_id <- gsub("_.*", "", grp)
    print(paste0(grp, "..."))
    seu <- tryCatch({
      mtx <- Read10X(data.dir = file.path(ref_seu_dir, grp), 
                     strip.suffix=TRUE)
      if(class(mtx)=='list') mtx <- mtx[['Gene Expression']]
      CreateSeuratObject(counts = mtx, project = gsm_id)
    }, error=function(e){cat("Failed\n"); NULL})
  }
  return(seu)
})
names(seus) <- relabel_id
seus <- seus[which(!sapply(seus, is.null))]

# For each sampel, drop all cells where there are at least 1 gene detected
sapply(seus,dim) %>% as.data.frame.matrix %>% t
seus <- lapply(seus, function(seu_i){
  subset(seu_i, subset = nFeature_RNA > 50)
})
sapply(seus,dim) %>% as.data.frame.matrix %>% t

saveRDS(seus, file=outrds)

##########################################################
#### 1.b Load in preprocessed/packaged seurat objects ####
# This section of the code loads in all the data that was
# released on GEO in a non-standard format (i.e. not features,
# barcodes, matrix.mtx.gz files). Since each dataset 
# had specifics on how they released their data, each 
# dataset was manually hardcoded to address and standardize
# the formatting.
outrds <- file.path(datadir, "seurat_obj", "ref1b_filtered.rds")

if(file.exists(outrds)){
  cat("Reading in existing 'seus' object...\n")
  seus_preprocessed <- readRDS(file=outrds)
} else {
  seus_preprocessed <- NULL
}

.relabel_origident <- function(seu, map, colid='orig.ident'){
  idx <- which(seu@meta.data[,colid] %in% names(map))
  seu$tmp <- NA
  seu$tmp[idx] <- as.character(map[seu@meta.data[idx,colid]])
  seu$orig.ident <- seu$tmp
  seu$tmp <- NULL
  return(seu)
}
.filterGSE216166 <- function(seu_ppi){
  cat("Fixing GSE216166...\n")
  seu_gse216166_map <- c('S1163175'='GSM6659730',
                         'S1163177'='GSM6659732',
                         'S1163179'='GSM6659734')
  
  Idents(seu_ppi) <- 'tissue_detail'
  seu_ppi <- subset(seu_ppi, idents='naive')
  seu_ppi$percent.mt <- seu_ppi$percent.mito
  seu_ppi <- .relabel_origident(seu_ppi, seu_gse216166_map)
  seu_ppi
}
.filterGSE117131 <- function(seu_ppi, ident_keep='Healthy'){
  cat("Fixing GSE117131...\n")
  seu_gse117131_map <- c('Mouse_Healthy1_001'='GSM3271853',
                         'Mouse_Healthy1_002'='GSM3271859',
                         'Mouse_Healthy2_001'='GSM3271855',
                         'Mouse_Healthy2_002'='GSM3271861',
                         'Mouse_Healthy3_001'='GSM3271857',
                         'Mouse_Healthy3_002'='GSM3271863',
                         'Mouse_Tumor1_001'='GSM3271852', 
                         'Mouse_Tumor1_002'='GSM3271858',
                         'Mouse_Tumor2_001'='GSM3271854',
                         'Mouse_Tumor2_002'='GSM3271860',
                         'Mouse_Tumor3_001'='GSM3271856',
                         'Mouse_Tumor3_002'='GSM3271862')
  
  seu_ppi$tissue_detail <- gsub("[0-9]_[0-9]*$", "", seu_ppi$Sample.Name) %>%
    gsub("Mouse_", "", .)
  Idents(seu_ppi) <- 'tissue_detail'
  # seu_ppi <- subset(seu_ppi, idents=ident_keep)
  seu_ppi <- PercentageFeatureSet(seu_ppi, pattern = "^mt-", col.name = "percent.mt")
  seu_ppi <- .relabel_origident(seu_ppi, seu_gse117131_map, colid='Sample.Name')
  seu_ppi
}

seus_preprocessed <- lapply(pp_datasets, function(pp_i){
  if(pp_i %in% names(seus_preprocessed)){
    cat(paste0("Using existing preprocessed data for ", pp_i, "\n"))
    seu_ppi <- seus_preprocessed[[pp_i]]
  } else {
    cat(paste0("Loading in count data and preprocessing ", pp_i, "\n"))
    dirx <- file.path(atlas_refdir, "gse_dat", pp_i)
    cts_f <- grep("counts.(txt|tsv|csv).gz$", list.files(dirx), value=T, perl=T)
    meta_f <- grep("(metadata|annotation).(txt|tsv|csv).gz$", list.files(dirx), value=T, perl=T)
    
    if(grepl("txt|tsv", cts_f)){
      cts <- read.table(file.path(dirx, cts_f), sep="\t", header = T)
      meta <- read.table(file.path(dirx, meta_f), sep="\t", header = T)
    } else {
      cts <- read.csv(file.path(dirx, cts_f))
      meta <- read.csv(file.path(dirx, meta_f))
    }
    if('Barcode' %in% colnames(meta)){
      meta <- meta  %>% 
        tibble::column_to_rownames('Barcode')
    }
    if('X' %in% colnames(cts)){
      cts <- cts  %>% 
        tibble::column_to_rownames('X')
    }
    
    seu_ppi <-   CreateSeuratObject(
      counts = cts, project = pp_i, 
      assay = "RNA", names.field = 1,
      names.delim = "_", meta.data = meta,
    )
    seu_ppi <- switch(pp_i,
                      GSE216166=.filterGSE216166(seu_ppi),
                      GSE117131=.filterGSE117131(seu_ppi),
                      seu_ppi)
    
    seu_ppi <- preprocessSeu(
      seu_ppi, ncount_min=1000, ncount_max=Inf, 
      nfeature_min=0, nfeature_max=Inf, 
      mt_max=Inf, org='mouse', numpcs=30, getPCs=FALSE
    )
  }
  return(seu_ppi)
})
names(seus_preprocessed) <- pp_datasets

saveRDS(seus_preprocessed, file=outrds)

###########################################
#### 2. Preprocess to remove doublets  ####
dir.create(file.path(outdir, "gex", "qc"), showWarnings = F, recursive = T)
dir.create(file.path(outdir, "seurat_obj", "per_sample"), showWarnings = F,)
if(!exists("seus")) seus <- readRDS(file.path(datadir, "seurat_obj", "ref1_unfiltered.rds"))
if(!exists("seus_preprocessed")) seus_preprocessed <- readRDS(file.path(datadir, "seurat_obj", "ref1b_filtered.rds"))

outrds_mnn <- file.path(datadir, "seurat_obj", "ref2_dubfiltered.mnn.rds")
outrds_stacas <- file.path(datadir, "seurat_obj", "ref2_dubfiltered.stacas.rds")

### Preprocess and Doublet Calling
names(seus) <- as.character(sapply(seus, function(i) unique(i$orig.ident)))
seus <- lapply(seus, function(seu_i){
  PercentageFeatureSet(seu_i, pattern = "^mt-", col.name = "percent.mt")
})
# seu <- merge(seus[[1]], y = seus[-1], add.cell.ids = names(seus), project = "myeloid_reference")

#.... 1) Identify mitochondrial/nCount/nFeature outliers: ----
outliers <- lapply(seus, function(seu_i){
  # 3-MAD/Mindiff from Median as outliers
  # mitochondrial outlier: scater::isOutlier
  mt_outlier_scuttle <- isOutlier(seu_i$percent.mt, nmads=3, 
                                  min.diff=(10-median(seu_i$percent.mt)), 
                                  batch=seu_i$orig.ident, share.mads=T, 
                                  share.medians=T, type="higher")
  # nCount outlier: scater::isOutlier
  ncnt_outlier_scuttle <- isOutlier(seu_i$nCount_RNA, nmads=4, 
                                    # min.diff=median(seu_i$nCount_RNA), 
                                    batch=seu_i$orig.ident, share.mads=F, 
                                    share.medians=F, type="higher")
  # nFeature outlier: scater::isOutlier
  nft_outlier_scuttle <- isOutlier(seu_i$nFeature_RNA, nmads=3, 
                                   # min.diff=median(seu$nFeature_RNA), 
                                   batch=seu_i$orig.ident, share.mads=F, 
                                   share.medians=F, type="higher")
  
  nmt.max <- min(seu_i$percent.mt[mt_outlier_scuttle]) # 6530
  nfeature.max <- min(seu_i$nFeature_RNA[nft_outlier_scuttle]) # 6530
  ncnt.max <- min(seu_i$nCount_RNA[ncnt_outlier_scuttle]) # 37252
  c('mt'=max(c(nmt.max, 10)), 
    'nft'=max(c(nfeature.max, 5000)), 
    'ncnt'=max(c(ncnt.max, 17000)))
})

pdf(file.path(outdir, "gex", "qc", "ref_qc_featureScatters.pdf"))
# pdf(file.path("~/xfer", "ref_qc_featureScatters.raw.pdf"))
fsl <- lapply(seus, function(seu_i){
  fs1 <- FeatureScatter(seu_i, 'nCount_RNA', 'nFeature_RNA', shuffle=TRUE, raster=TRUE)
  fs2 <- FeatureScatter(seu_i, 'nCount_RNA', 'percent.mt', shuffle=TRUE, raster=TRUE)
  cowplot::plot_grid(fs1, fs2, ncol=2)
})
print(fsl)
dev.off()

#.... 2) Flag Potential Doublets ----
if(file.exists(outrds_int)) {
  seus_post <- readRDS(outrds_int)
  # names(seus_post) <- names(seus)
} else{
  seus_post <- NULL
}

#a) Preprocess for Doublet Calling
s.features <- stringr::str_to_title(cc.genes$s.genes)
g2m.features <- stringr::str_to_title(cc.genes$g2m.genes)

seu_ids <- setdiff(names(seus), names(seus_post))
seus_pre <- lapply(seu_ids, function(seu_id){
  if(seu_id %in% names(seus_post)){
    cat(paste0("Loading in existing preprocessed data: ", seu_id, "\n"))
    seu_i <- seus_post[[seu_id]]
  } else {
    cat(paste0("Preprocessing and adding on ", seu_id, "\n"))
    seu_i <- seus[[seu_id]]
    
    gsm_id <- as.character(unique(seu_i$orig.ident))
    outlier_i <- outliers[[gsm_id]]
    seu_i <- preprocessSeu(seu_i, ncount_min=1000, ncount_max=Inf, 
                           nfeature_min=50, nfeature_max=outlier_i['nft'], 
                           mt_max=outlier_i['mt'], org='mouse', numpcs=30, getPCs=FALSE)
  }
  return(seu_i)
})

#b) Doublet Calling
if(call_doublets){
  cat("Calling doublets...\n")
  seus_pre2 <- lapply(seus_pre, function(seu_i){
    ## DoubletDetection
    sce_cxds <- runCxds(seu_i, n=25)  # Looks for genes that are mutually exclusive, scores their expression
    seu_i$cxds_z <- sce_cxds$cxds$z_cxds
    seu_i$cxds_doublet <- (seu_i$cxds_z>3)
    
    sce_bcds <- runBcds(seu_i, n=25) # Creates artificial doublets, trains a classifier to predicts them
    seu_i$bcds_score <- sce_bcds$bcds[,1]
    seu_i$bcds_doublet <- (seu_i$bcds_score>0.95)
    
    seu_i
  })
  seus_post <- c(seus_post, seus_pre2)
} else {
  cat("Skipping doublet calling...\n")
  new_ids <- c(names(seus_post), seu_ids)
  seus_post <- c(seus_post, seus_pre)
}
names(seus_post) <- new_ids
saveRDS(seus_post, file=outrds_int)
if(!exists('seus_post')) seus_post <- readRDS(outrds_int)

#.... 3) Integrate datasets ----

# Integrate the data: FastMNN
seus_post <- seus_post[which(!duplicated(names(seus_post)))]
seus_post2 <- subsetSeu(c(seus_post, seus_preprocessed), n=1000)
seu.mnn <- RunFastMNN(object.list = seus_post2, 
                      features = 2000, assay='SCT')
seu.mnn <- RunUMAP(seu.mnn, reduction = "mnn", dims = 1:50, n.neighbors=30L,
                   min.dist = 0.1, return.model=TRUE)
seu.mnn <- FindNeighbors(seu.mnn, reduction = "mnn", dims = 1:50)
seu.mnn <- FindClusters(seu.mnn, resolution = 0.9, graph.name='SCT_snn')
DefaultAssay(seu.mnn) <- 'RNA'
seu.mnn <- NormalizeData(seu.mnn,
                         normalization.method = "LogNormalize") %>%
  FindVariableFeatures(., selection.method = "vst",
                       nfeatures = 3000, verbose = FALSE) %>% 
  ScaleData(.)
seu.mnn <- RunPCA(seu.mnn, features=VariableFeatures(seu.mnn), ndims.print=1:50,
                  nfeatures.print=5, npcs=50)

seu.mnn$dataset_celltype <- ref_celltypes[seu.mnn$orig.ident,]$ID
pdf("~/xfer/x.pdf", width = 13, height = 12)
DimPlot(seu.mnn, group.by='dataset_celltype', raster=T, reduction='umap')
lapply(unique(seu.mnn$dataset_celltype), function(celltype){
  seu.mnn$celltype <- ifelse(seu.mnn$dataset_celltype == celltype, celltype, 'NA')
  DimPlot(seu.mnn, group.by='celltype', raster=T, order=celltype, 
          cols=setNames(c('grey', '#de2d26'), c('NA', celltype)),
          reduction='umap') + 
    ggtitle('') + ylab('') + xlab(celltype) + 
    theme(legend.position='none', 
          axis.text.x = element_blank(),
          axis.text.y = element_blank())
}) %>%
  cowplot::plot_grid(plotlist=., ncol = 5)
dev.off()
saveRDS(seu.mnn, file=outrds_mnn)
# saveRDS(seu.mnn, file=gsub(".rds$", ".sara.rds", outrds_mnn))


# Integrate the data: STACAS
seu.stacas <- lapply(seus_post2, function(seu_i){
  DefaultAssay(seu_i) <- 'RNA'
  return(seu_i)
}) %>%
  Run.STACAS(., future.maxSize=25)

seu.stacas <- RunUMAP(seu.stacas, dims = 1:30, n.neighbors=30L,
                      min.dist = 0.1)
seu.stacas <- FindNeighbors(seu.stacas, dims = 1:30)
seu.stacas <- FindClusters(seu.stacas, resolution = 1.2)

seu.stacas$dataset_celltype <- ref_celltypes[seu.stacas$orig.ident,]$ID
# saveRDS(seu.stacas, file=gsub(".rds$", ".sara.rds", outrds_stacas))


#.... 4) Save and Print Dimplots ----

rm(seus, seus_pre, seus_post)
# rm(, seu.mnn, seu.stacas)

# DimPlot(seu.stacas, group.by = c("Method","CellType")) 
pdf(file.path(outdir, "gex", "qc", "ref_integration_dimplot.pdf"), width = 10, height = 8)
# pdf(file.path("~/xfer", "ref_integration_dimplot.pdf"), width = 12, height = 8)
DimPlot(seu.stacas, group.by = c("orig.ident", "dataset_celltype"), raster=TRUE) +
  ggtitle("STACAS")
DimPlot(seu.stacas, split.by = "dataset_celltype", raster=TRUE, ncol=4) +
  ggtitle("STACAS")
DimPlot(seu.mnn, group.by = c("orig.ident", "dataset_celltype"), raster=TRUE)  +
  ggtitle("MNN")
DimPlot(seu.mnn,  split.by = "dataset_celltype", raster=TRUE, ncol=4) +
  ggtitle("MNN")
dev.off()

###########################################
#### 3. Create a ProjecTILs reference  ####
dir.create(projectils_refdir, showWarnings = F, recursive = F)
if(!exists("seu.stacas")) seu.stacas <- readRDS(file.path(datadir, "seurat_obj", "ref2_dubfiltered.stacas.rds"))
if(!exists("seu.mnn")) seu.mnn <- readRDS(file.path(datadir, "seurat_obj", "ref2_dubfiltered.mnn.rds"))

DefaultAssay(seu.stacas) <- 'integrated'
ref.stacas <- make.reference(ref = seu.stacas, ndim = 30, seed = seed, recalculate.umap = FALSE,
                             annotation.column = "dataset_celltype")

seu.stacas.small <- subsetSeu(seu.stacas, 'dataset_celltype', 1000, seed=seed)
ref.stacas.small <- make.reference(ref = seu.stacas.small, assay='integrated',
                                   ndim = 25, seed = seed, recalculate.umap = FALSE,
                                   annotation.column = "dataset_celltype")

seu.mnn <- RenameAssays(seu.mnn, mnn.reconstructed='integrated')
DefaultAssay(seu.mnn) <- 'integrated'
ref.mnn <- make.reference(ref = seu.mnn, assay='integrated',
                          ndim = 25, seed = seed, recalculate.umap = FALSE,
                          annotation.column = "dataset_celltype")

seu.mnn.small <- subsetSeu(seu.mnn, 'dataset_celltype', 1000, seed=seed)
ref.mnn.small <- make.reference(ref = seu.mnn.small, assay='integrated',
                                ndim = 25, seed = seed, recalculate.umap = FALSE,
                                annotation.column = "dataset_celltype")
# query <- readRDS("/cluster/projects/mcgahalab/data/mcgahalab/teresa_ly6/scrna/data/seurat_obj/2_dubfiltered.mnn.rds")
# x <- ref.mnn.small@misc$pca_object
# x$rotation %>% head

# x <- subset(ref.mnn, cells=sample(Cells(ref.mnn), size=1000))  # Makes sure the renaming didn't botch the object

saveRDS(ref.stacas, file=file.path(projectils_refdir, "seu.stacas.rds"))
saveRDS(ref.mnn, file=file.path(projectils_refdir, "seu.mnn.rds"))
# saveRDS(ref.mnn, file=file.path(projectils_refdir, "seu.mnn.small.sara.rds"))
saveRDS(ref.mnn.small, file=file.path(projectils_refdir, "seu.mnn.small.rds"))
# saveRDS(ref.stacas, file=file.path(projectils_refdir, "seu.stacas.small.sara.rds"))
saveRDS(ref.stacas.small, file=file.path(projectils_refdir, "seu.stacas.small.rds"))
# ref.mnn <- readRDS(file=file.path(projectils_refdir, "seu.mnn.rds"))

pdf("~/xfer/dimplot_ref.pdf", width = 9)
DimPlot(ref.stacas, label = T, repel = T, label.size = 4, reduction='umap')
DimPlot(ref.mnn, label = T, repel = T, label.size = 4, reduction='umap')
pdf("~/xfer/dimplot_ref3.pdf", width = 9)
DimPlot(ref.mnn.small, label = T, repel = T, label.size = 4, reduction='umap')
dev.off()
