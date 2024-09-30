# renv::load("/cluster/home/quever/downloads/renvs/")

library(tidyverse)
library(ComplexHeatmap)
library(ggrastr)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(umap)
library(cowplot)
library(ggrepel)
library(GSVA)
library(org.Hs.eg.db)
library(icellnet)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(gridExtra)
library(jetset)
library(RColorBrewer)
library(clusterProfiler)


colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
species <- 'Mus musculus'
pdir <- file.path('/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_MSM_SSM_123/results')

gprofiler_dir <- '/cluster/projects/mcgahalab/ref/gprofiler'
gprofiler_f <- file.path(gprofiler_dir, 'gprofiler_full_mmusculus.ENSG.gmt')
# 
barcodes_f="~/git/mini_projects/ref/barcodes.tsv"

#### Functions ####
# code cloned from https://github.com/quevedor2/mini_projects
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/makeLoupe.R")
gm <- geneMap(species=species)

###################
#### Functions ####
ssGseaFun <- function(msig_ds, lfc_v, ss_method='ssgsea'){
  require(GSVA)
  ssgsea <- tryCatch({
    sig_ens_gs <- split(setNames(msig_ds$entrez_gene, msig_ds$entrez_gene), 
                        f=msig_ds$gs_name)
    gsva(lfc_v, sig_ens_gs, verbose=FALSE, method=ss_method)
  }, error=function(e){NULL})
  return(ssgsea)
}

getDEGedgeRwt <- function(cts, meta, group){
  idx <- match(rownames(meta), colnames(cts))
  if(any(is.na(idx))){
    na_idx <- which(is.na(idx))
    # print(paste0("NA idx: ", paste(na_idx, collapse=",")))
    idx <- idx[-na_idx]
    meta <- meta[-na_idx,]
  }
  
  ## Get EdgeR Results
  se <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(cts[,idx])), 
                             colData = meta)
  
  
  edgey <- tryCatch({
    edgey <- DGEList(counts=cts[,idx],
                     samples=rownames(meta),
                     group=meta[,group])
    # design <- with(meta, model.matrix(as.formula(formula)))
    design <- model.matrix(~ meta[,group])
    keep <- filterByExpr(edgey, design)
    edgey <- edgey[keep, , keep.lib.sizes=FALSE]
    edgey <- calcNormFactors(edgey, method = "TMM")
    edgey <- estimateDisp(edgey)
    edgey
  }, error=function(e){NULL})
  
  # Calculate CPM/TMM values
  edge_tmm <- cpm(edgey)
  edge_tmm_spl <- edge_tmm %>% t %>% as.data.frame %>%
    split(., meta[,group]) 
  
  # Differential testing
  er_dat <- exactTest(edgey)
  et_res <- er_dat$table  %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("ensemble") %>% 
    mutate(padj=p.adjust(PValue, method='BH')) %>%
    mutate(sig=ifelse(padj < 0.05, "er_sig", "er_ns"))
  
  glvl <- levels(factor(meta[,group]))
  wilcox_res <- sapply(seq_along(edge_tmm_spl[[1]]), function(idx){  
    wt <- wilcox.test(edge_tmm_spl[[glvl[1]]][,idx], edge_tmm_spl[[glvl[2]]][,idx])
    fc <- mean(edge_tmm_spl[[glvl[1]]][,idx],na.rm=T) /  mean(edge_tmm_spl[[glvl[2]]][,idx], na.rm=T)
    c('W'=wt$statistic, 'pval'=wt$p.value, "FC"=fc, 'Log2FC'=log2(fc+1),
      'ensemble'=colnames(edge_tmm_spl[[1]][idx]))
  }) %>% 
    t %>% as.data.frame %>%
    mutate(padj=p.adjust(pval, method='BH'),
           FC=as.numeric(FC),
           Log2FC=as.numeric(Log2FC)) %>%
    mutate(sig=ifelse(padj < 0.05, "wt_sig", "wt_ns"))
  
  # pdf("~/xfer/dds2.pdf") 
  # # ggplot(dds_et, aes(x=log2FoldChange, y=logFC)) +
  # ggplot(dds_et, aes(x=log2FoldChange, y=logFC, color=sig, group=sig)) +
  #   geom_point()
  # ggplot(wt_et, aes(x=Log2FC, y=logFC, color=sig, group=sig)) +
  #   geom_point()
  # dev.off()
  return(list("edger"=et_res,  "wilcox"=wilcox_res))
}

rmVstBatchEffects <- function(vsd, dds, samples_meta, 
                              condition='condition', batchcol='batch'){
  mat <- assay(vsd)
  mm <- model.matrix(as.formula(paste0("~", condition)), colData(dds))
  mat <- limma::removeBatchEffect(mat, design=mm, 
                                  batch=samples_meta[colnames(dds),batchcol])
  assay(vsd) <- mat
  return(vsd)
}


##############
#### Main ####
dir.create(file.path(pdir, "manual", "objs"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "pca"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "loupe"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "deg"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "gsea"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "ssgsea"), showWarnings = F, recursive = T)

rm.samples <- c('XS_A', 'XS_9', 'XS_11')

deganaldir <- file.path(pdir, "manual", "differential_expression")
file <- 'all.tsv'
min_expr <- 3
expr_frac_samples <- 0.2

# Load bulk RNAseq data
data=read.table(file.path(pdir, "counts", file), header = T, check.names=FALSE, 
                stringsAsFactors = FALSE, na.strings = "")
genes <- gm$ENSEMBL$SYMBOL[data$gene]
genes[is.na(genes)] <- 'NA'
if(any(genes=='')){
  genes[genes==''] <- names(genes[genes==''])
}


data <- data %>% 
  dplyr::select(-c(gene)) %>% 
  split(., f=genes) %>%
  sapply(., colMeans)  %>% 
  t  %>% as.data.frame %>%
  rename_with(., ~gsub("genesresults$", "", .)) 

if(length(rm.samples)>0){
  rm.samples <- sapply(rm.samples, grep, x=colnames(data), value=T)
  keep.samples <- setdiff(colnames(data), rm.samples)
  data <- data[,keep.samples]
}



# Remove genes where [expr_frac_sample] of the dataset has fewer than [min_expr] 
# reads linked to that gene
# e.g. Remove: LOC101929116 is expressed at higher than 3 counts in only 2/42 samples
low_expr_idx <- which(rowSums(data >= min_expr) >= 2 ) #(ncol(data) * expr_frac_samples))
data <- data[low_expr_idx,]



#### 1. Create metadata ####
# Create metadata data frame
metadata <- colnames(data)
metad <- strsplit(gsub("_S[0-9]*$", "", metadata) %>%
                    gsub("(.*)([0-9])$", "\\1_\\2", .), split="_") %>% 
  do.call(rbind, .) %>% 
  magrittr::set_colnames(., c('Macrophage', 'Location', 'Treatment', 'Replicate')) %>%
  as.data.frame %>% 
  dplyr::mutate('sample'=as.character(metadata),
                'condition'=paste(Macrophage, Location, Treatment, sep="_"))
saveRDS(metad, file.path(pdir, "manual", "objs", "metad.rds"))

data <- as.matrix(data)
storage.mode(data) <- 'integer'
dds_all <- DESeqDataSetFromMatrix(countData=data,
                                  colData=metad,
                                  design=as.formula("~condition"))
saveRDS(dds_all, file=file.path(pdir, "manual", "objs", "deseq_all.rds"))

