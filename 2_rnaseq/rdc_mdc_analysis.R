## Sara's mouse_DC project
# Comparing rDC and mDC between Cis and PBS treated samples in the WT and Tumor
library(cowplot)
library(igraph)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(org.Mm.eg.db)
library(dplyr)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)

# SCENIC params
scenic_org_code <- 'mgi'
rcis_db <- '/cluster/projects/mcgahalab/ref/scenic/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

# Params
pdir <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_DCs"
dedir <- file.path(pdir, "results", "diffexp")
outdir <- file.path(pdir, "results", "manual")
dir.create(outdir, recursive = F, showWarnings = F)
dir.create(file.path(outdir, "rds"), showWarnings = F)
setwd(pdir)

dds_main <- readRDS(file.path("results", "deseq2", "all.rds"))
rm_idx <- which(colnames(dds_main) %in% rm_samples)
dds_main <- dds_main[,-rm_idx]

#######################
#### Functions Git ####
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/st2_il33/functions/msm_ssm_functions.R")

gm <- geneMap("Mus musculus")


#############################################
#### 0. Create DESeq object and Metadata ####
cts <- read.table(file.path(pdir, "results", "counts", "mdc_rdc.counts.tsv"),
                  header = T, check.names = F, stringsAsFactors = F) %>%
  tibble::column_to_rownames(., "gene") %>%
  rename_with(., ~gsub("genesresults", "", .))
meta <- read.table(file.path(pdir, "config", "mdc_rdc.samples.tsv"),
                  header=T, check.names=F, stringsAsFactors = F)  %>%
  mutate("group"=gsub("^.*_([ABCDE]).*", "\\1", sample_name),
         "modifier"=gsub("^.*_[ABCDE]_?(.*)", "\\1", sample_name),
         "dc"=gsub("_.*", "", sample_name),
         "treatment"=gsub("^.*_", "", condition),
         "wt.t"=gsub("^.*_(.*)_.*", "\\1", condition)) %>%
  tibble::column_to_rownames(., "sample_name")
meta <- meta[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=meta,
                              design=as.formula('~condition'))

# remove uninformative columns
min.cnt <- 5
min.n <- 0.1
dds <- dds[ which(rowSums(counts(dds)>min.cnt) > (ncol(dds)*min.n)),]
# normalization and preprocessing
dds <- DESeq(dds)

dir.create(file.path(pdir, "results", "deseq2"))
dir.create(file.path(pdir, "results", "metadata"))
saveRDS(dds, file=file.path(pdir, "results", "deseq2", "all.rds"))
saveRDS(dds, file=file.path(pdir, "results", "metadata", "meta.rds"))

############################
#### 1.a) PCA analysis  ####
# Run the main PCA analysis on vst-counts
grp_ggps <- lapply(seq_along(sample_l), function(sids_idx){
  sids <- sample_l[[sids_idx]]
  sample_id <- names(sample_l)[sids_idx]
  print(sample_id)
  sample_idx <- which(rownames(colData(dds_main)) %in% sids)
  dds <- dds_main[,sample_idx]
  
  # obtain normalized counts
  counts <- vst(dds, blind=T)
  tcnts <- as.data.frame(t(assay(counts)))
  pca <- prcomp(tcnts, scale=F)
  max_pc_i <- min(c(max_pc, ncol(pca$x)))
  pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc_i))],
                               "condition"=as.character(counts$condition)))
  for(id in paste0("PC", c(1:max_pc_i))){
    pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
  }
  
  # Format in any factors/metadata you want to plot by
  pca_y <- pca_x %>% 
    tibble::rownames_to_column(., "Name") %>%
    left_join(samples_meta %>%
                tibble::rownames_to_column(., "samples"), 
              by=c('Name'='samples')) %>%
    mutate(tissue_treatment=paste0(dc, "_", tissue, "_", treatment))


  pcs <- paste0("PC", c(1:max_pc_i))
  pcs_l <- lapply(c(2:length(pcs)), function(i){c(pcs[[i-1]], pcs[[i]])})
  
  # grp <- 'sample_ct'  # Treatment specific triangulation
  # if(sids_idx == length(sample_l)) grp <- 'Sample' # Subject specific triangulation
  ggps <- lapply(pcs_l, function(pc){
    ggplot(data=pca_y, aes_string(x=pc[1], y=pc[2], shape='treatment',
                                  color='tissue_treatment', fill='tissue_treatment',
                                  label='id')) +
      geom_point() + 
      geom_text_repel() +
      scale_shape_manual(values=c(21,24)) +
      geom_line(aes(group = tissue_treatment), alpha=0.3)  +
      ggtitle(paste0(sample_id, ": ", paste(pc, collapse="-")))
  })
  
  return(list("ggp"=ggps, "pca"=pca, "pca_y"=pca_y))
})
names(grp_ggps) <- names(sample_l)

gg_pcas <- lapply(grp_ggps, function(grp) grp$ggp)
pcas <- lapply(grp_ggps, function(grp) grp$pca)
pca_ct <- lapply(grp_ggps, function(grp) grp$pca_y) # PCA vals + meta

pdf(file.path(pca_dir, "pca.pdf"))
gg_pcas
dev.off()

# Investigating PC components that best predict metadata components
pc_mat <- pca_ct$All %>%
  tibble::column_to_rownames(., var="Name") %>%
  select(grep("^PC", colnames(.)))
meta_mat <- pca_ct$All %>%
  tibble::column_to_rownames(., var="Name") %>%
  select(grep("^PC", colnames(.), invert = T)) %>% 
  rename_with(., ~ gsub("[[:punct:]]", "", .x))
if(any(is.na(meta_mat))) meta_mat[is.na(meta_mat)] <- 'NA'


saveRDS(pc_mat, file=file.path(pca_dir, "pcmat.rds"))
saveRDS(meta_mat, file=file.path(pca_dir, "metamat.rds"))
saveRDS(pcas, file=file.path(pca_dir, "pcas.rds"))

##############################################################################
#### 2.a) Differential expression controlling for interaction terms - Raw ####
cntsdir <- file.path(pdir, "results", "counts")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(file.path(cntsdir, "all.tsv"), header=TRUE, 
                  row.names="gene", check.names=FALSE)
colnames(cts) <- gsub("B_1", "B1", colnames(cts)) %>% gsub("B_2", "B2",. )
if(any(!colnames(cts) %in% rownames(samples_meta))){
  # ensure coldata and cts are in the same sample order
  cts <- cts[,match(rownames(samples_meta), colnames(cts))]
}

samples_meta$tissue <-factor(samples_meta$tissue, c("WT", "Tumor"))
samples_meta$treatment <- factor(samples_meta$treatment, c('PBS', 'CIS'))
samples_meta$dc <- factor(samples_meta$dc, c('rDC', 'mDC'))
samples_meta_l <- c(split(samples_meta, samples_meta$dc), 
               list("All"=samples_meta))
samples_meta_treat_ln_l <- split(samples_meta, list(samples_meta$treatment, samples_meta$tissue))

celltypes <- setNames(names(samples_meta_l), names(samples_meta_l))
resl <- lapply(celltypes, function(celltype){
  ## INTERACTION TERMS
  # if the log2 fold change attributable to a given condition is different 
  # based on another factor, for example if the treatment effect differs 
  # across lnstatus
  samples_meta_l[[celltype]]$treatment <- relevel(samples_meta_l[[celltype]]$treatment, "PBS")
  samples_meta_l[[celltype]]$tissue <- relevel(samples_meta_l[[celltype]]$tissue, "WT")
  
  dds_raw <- DESeqDataSetFromMatrix(countData=cts[,rownames(samples_meta_l[[celltype]])],
                                    colData=samples_meta_l[[celltype]],
                                    design=as.formula('~tissue+treatment+tissue:treatment'))
  
  # remove uninformative columns
  dds <- dds_raw[ rowSums(counts(dds_raw)) > 1, ]
  # normalization and preprocessing
  dds <- DESeq(dds)
  
  ## MAIN TREATMENT EFFECT 1
  # treatment_CIS_vs_PBS
  # the main treatment effect only represents the effect of treatment 
  # for the reference level of lymph-node status.
  #
  # Relates to the treatment effect of CIS vs PBS for LN
  # - Need to relevel to evaluate CIS vs PBS for TDLN
  resultsNames(dds)
  coef <- 'treatment_CIS_vs_PBS'
  overall_reslfc_ln <- lfcShrink(dds, coef=coef, type="apeglm")
  overall_ma_ln <- plotMA(overall_reslfc_ln, ylim=c(-3,3), cex=.8)
  abline(h=c(-1,1), col="dodgerblue", lwd=2)
  overall_reslfc_ln <- overall_reslfc_ln %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".ln_cis-pbs")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## INTERPRETING INTERACTION TERMS
  # lnstatusTDLN.treatmentCIS
  # interaction terms lnstatusTDLN.treatmentCIS gives the DIFFERENCE between
  # the treatment effect for a given lymph-node status, and the treatment 
  # effect for the reference lymph-node status
  #
  # Relates to the treatment effect for TDLN vs LN
  # - whether the effect of treatment is different between TDLN and LN
  resultsNames(dds)
  coef <- 'tissueTumor.treatmentCIS'
  interact_reslfc <- lfcShrink(dds, coef=coef, type="apeglm")
  interact_reslfc <- interact_reslfc %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".interaction")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## MAIN TREATMENT EFFECT 2
  # Relates to the treatment effect of CIS vs PBS for TDLN
  # - By setting the reference level to TDLN, we can evaluate CIS vs PBS
  # differences in terms of TDLN status
  samples_meta_l[[celltype]]$tissue <- relevel(samples_meta_l[[celltype]]$tissue, "Tumor")
  dds2 <- DESeqDataSetFromMatrix(countData=cts[,rownames(samples_meta_l[[celltype]])],
                                 colData=samples_meta_l[[celltype]],
                                 design=as.formula('~tissue+treatment+tissue:treatment'))
  
  # remove uninformative columns
  dds2 <- dds2[ rowSums(counts(dds2)) > 1, ]
  dds2 <- DESeq(dds2)
  coef <- 'treatment_CIS_vs_PBS'
  overall_reslfc_tdln <- lfcShrink(dds2, coef=coef, type="apeglm")
  overall_reslfc_tdln <-  overall_reslfc_tdln %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".tumor_cis-pbs")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## MAIN TREATMENT EFFECT 3
  # Relates to the treatment effect of TDLN-CIS vs LN-CIS
  # - By setting the reference level to LN and subsetting for only the CIS samples, 
  # we can evaluate TDLN-CIS vs LN-CIS differences
  samples_meta_l[[celltype]]$tissue <- relevel(samples_meta_l[[celltype]]$tissue, "WT")
  cis_samples_meta <- samples_meta_l[[celltype]] %>% 
    filter(treatment=='CIS')
  dds3 <- DESeqDataSetFromMatrix(countData=cts[,rownames(cis_samples_meta)],
                                 colData=cis_samples_meta,
                                 design=as.formula('~tissue'))
  
  # remove uninformative columns
  dds3 <- dds3[ rowSums(counts(dds3)) > 1, ]
  dds3 <- DESeq(dds3)
  coef <- resultsNames(dds3)[2]
  overall_reslfc_cis <- lfcShrink(dds3, coef=coef, type="apeglm")
  overall_reslfc_cis <- overall_reslfc_cis %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".cis_tumor-wt")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## MAIN TREATMENT EFFECT 4
  # Relates to the treatment effect of PBS-TDLN vs PBS-LN
  # - By setting the reference level to LN and subsetting for only the PBS samples, 
  # we can evaluate PBS-TDLN vs PBS-LN differences
  samples_meta_l[[celltype]]$tissue <- relevel(samples_meta_l[[celltype]]$tissue, "WT")
  pbs_samples_meta <- samples_meta_l[[celltype]] %>% 
    filter(treatment=='PBS')
  dds4 <- DESeqDataSetFromMatrix(countData=cts[,rownames(pbs_samples_meta)],
                                 colData=pbs_samples_meta,
                                 design=as.formula('~tissue'))
  
  # remove uninformative columns
  dds4 <- dds4[ rowSums(counts(dds4)) > 1, ]
  dds4 <- DESeq(dds4)
  coef <- resultsNames(dds4)[2]
  overall_reslfc_pbs <- lfcShrink(dds4, coef=coef, type="apeglm")
  overall_reslfc_pbs <- overall_reslfc_pbs %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".pbs_tumor-wt")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## Identify the differential genes
  res <- Reduce(function(x,y) merge(x,y,by='ens'), list(overall_reslfc_tdln, 
                                                        overall_reslfc_ln, 
                                                        overall_reslfc_cis,
                                                        overall_reslfc_pbs,
                                                        interact_reslfc))
  res$gene <- gene_ids[res$ens]
  res$gene[is.na(res$gene)] <- res$ens[is.na(res$gene)]
  
  res_sdigits <- res %>%
    mutate(biotype=ens2biotype_ids[res$ens]) %>%
    relocate(c('gene', 'biotype'), .after=ens)
  for(col_i in colnames(res)){
    if(is.numeric(res[,col_i])){
      res_i <- res[,col_i]
      sdig <- if(mean(abs(res_i), na.rm=T) > 10) 2 else 6
      res_sdigits[,col_i] <- round(res[,col_i], sdig)
    }
  }
  
  write.table(res_sdigits, file=file.path(outdir, "degs", paste0(celltype, "_degs.csv")), 
              col.names = T,row.names = F, quote=F, sep=",")
  saveRDS(res, file=file.path(outdir, "rds", paste0(celltype, "_degs.rds")))
  return(list("res"=res, "dds_lnBase"=dds, "dds_tdlnBase"=dds2))
})


## Compares mDC vs rDC for each grouping (e.g. WT_PBS or WT_CIS)
resl_msm_v_ssm <- lapply(samples_meta_treat_ln_l, function(col_grp){
  col_grp$dc <- factor(col_grp$dc)
  col_grp$dc <- relevel(col_grp$dc, "rDC")
  
  dds <- DESeqDataSetFromMatrix(countData=cts[,rownames(col_grp)],
                                colData=col_grp,
                                design=as.formula('~dc'))
  
  # remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  
  ## MAIN TREATMENT EFFECT 1
  # celltype_MSM_vs_SSM
  # Compares the differences between MSM and SSM celltypes under every condition
  lbl <- with(col_grp, paste(tissue, treatment, sep="-")) %>% unique
  reslfc_ct <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") %>%
    as.data.frame %>%
    rename_with(., ~ paste0(., ".dc_", lbl)) %>%
    tibble::rownames_to_column(var='ens') 
  return(list("dds"=dds, "res"=reslfc_ct))
})
res_msm_ssm <- lapply(resl_msm_v_ssm, function(i) i$res) %>% 
  purrr::reduce(., full_join, by='ens') %>% 
  mutate("symbol"=gm$ENSEMBL$SYMBOL[ens]) %>%
  relocate(ens, symbol)
resl[['Celltype']] <- c(lapply(resl_msm_v_ssm, function(i) i$res),
                        list("res"=res_msm_ssm))

saveRDS(resl, file=file.path(outdir, "rds", "resl.rds"))
saveRDS(samples_meta, file=file.path(outdir, "rds", "samples_meta.rds"))


### DEMO: Explaining how to read the results from the DEG ----
do_demo <- FALSE
if(do_demo){
  gene <- rownames(resl$MSM$res)[1]
  dds <- resl$MSM$dds_lnBase
  dl <- lapply(colnames(dds@colData), function(i){
    d <- plotCounts(dds, gene=gene, intgroup=i, 
                    returnData=TRUE)
    d$ID <- rownames(d)
    d
  })
  d <- Reduce(function(x,y) merge(x,y,by=c('ID', 'count')), dl)
  d$count <- log2(d$count+1)
  d$treatment <- relevel(d$treatment, "PBS")
  pdf("~/xfer/test.pdf")
  ggplot(d, aes(x=treatment, y=count)) +
    facet_grid(cols=vars(lnstatus)) +
    geom_point()
  dev.off()
  
  x <- sapply(split(d$count, f=d$condition), mean)
  (x[4]-x[3]) - (x[2]-x[1])
}


######################################
#### 3. GSEA analysis of the DEGs ####
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
min_baseMean <- 5

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) # ontology gene sets
                  

gseas_all <- lapply(setNames(c('rDC', 'mDC'), c('rDC', 'mDC')), function(res_id){
  # get LFC table
  print(res_id)
  res <- lfc_df <- res_dds[[res_id]]$res 
  # if(filter_genes) res <- res %>% filter(ens %in% genes_to_include$V1)
  
  lfc_df <- res %>% 
    filter(baseMean.interaction > min_baseMean) %>%
    dplyr::select(grep("ens|log2FoldChange", colnames(.), value=T)) %>%
    tibble::column_to_rownames(., "ens") #%>% 
    #select(grep("interaction", colnames(.), value=T))
  
  # get GSEA table
  gseas <- apply(lfc_df, 2, function(lfc_v){
    print("...")
    lfc_v <- setNames(lfc_v,
                      gm$ENSEMBL$ENTREZ[rownames(lfc_df)])
    
    iterateMsigdb(species='Mus musculus', msig_lvls=msig_lvls, 
                  fun=gseaFun, lfc_v=lfc_v)
  })
  return(gseas)
})
saveRDS(gseas_all, file=file.path(outdir, 'rds', "gseas_all.rds"))
gseas_all <- readRDS(file=file.path(outdir, 'rds', "gseas_all.rds"))

# Write out the GSEA to cytoscape format
gs_map <- .mapGsToExactSource(species="Mus musculus")
for(celltype in names(gseas_all)){
  dir.create(file.path(outdir, "gse", "cytoscape"), recursive = T, showWarnings = F)
  gsea_l <- unlist(gseas_all[[celltype]]$log2FoldChange.interaction[-4], recursive=F)
  gsea_df <- do.call(rbind, lapply(gsea_l, as.data.frame))
  gsea_dat_i <- gsea2CytoscapeFormat(gsea_df, gs_map)
  write.table(gsea_dat_i$up, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("gsea.", celltype, ".up.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(gsea_dat_i$down, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("gsea.", celltype, ".down.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
}

## continue with merging gsea into one dataframe
gsea_dcs <- lapply(gseas_all, function(gseas){  
  # Aggregate all the GSEA data into one dataframe
  gsea <- lapply(gseas, unlist, recursive=F)
  gsea_df <- lapply(names(gsea), function(gsea_id) {
    gsea_i <- gsea[[gsea_id]]
    gsea_id <- gsub("log2.*\\.", "", gsea_id)
    
    gsea_df_i <- do.call(rbind, lapply(gsea_i, as.data.frame))
    gsea_df_i %>% 
      select(ID, enrichmentScore, NES, pvalue, p.adjust) %>% 
      rename_with(., ~paste0(gsea_id, ".", .), 
                  .cols = c(enrichmentScore, NES, pvalue, p.adjust))
  }) %>% 
    purrr::reduce(full_join, by='ID') %>% 
    mutate("Dataset"=gsub("_.*", "", ID),
           "Celltype"=res_id) %>%
    relocate(Dataset, Celltype)
  
  return(gsea_df)
})
saveRDS(gsea_dcs, file=file.path(outdir, 'rds', "gsea_dcs.rds"))

for(id in names(gsea_dcs)){
  write.table(gsea_dcs[[id]], 
              file=file.path("~/xfer/", paste0("gsea.", id, ".csv")),
              sep=",", quote = F, col.names = T, row.names = F)
}




#### 4) Final Figures: ####
#--- SFig2.b: Volcano Plot of rDC and mDC - DEGs for Tumor-CIS vs Tumor-PBS ----
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))

lfc_lim <- 2
qval_lim <- 0.01
maxq <- 100
cols <- c("sigup"="#fc8d59", "sigdown"="#91bfdb")

# comparisons <- names(res_dds[[1]])
# lapply(setNames(comparisons, comparisons), function(comp_i){
#   comp_name <- switch(comp_i,
#                       res='Interaction',
#                       dds_lnBase='Ref:LN',
#                       dds_tdlnBase='Ref:TDLN')
# })
gg_volcanos <- lapply(c('rDC', 'mDC'), function(id){
  res_dds_i <- res_dds[[id]]
  res_df <- res_dds_i$res %>% 
    dplyr::select(ens, gene, `log2FoldChange.interaction`, `padj.interaction`) %>%
    rename_with(., ~c('ens', 'gene', 'logFC', 'FDR')) %>%
    mutate(id=id,
           log10FDR=(-1*log10(FDR)),
           dir=ifelse(FDR<=qval_lim, 
                      ifelse(logFC>=lfc_lim, 
                             "sigup", 
                             ifelse(logFC<=(-1*lfc_lim),
                                    "sigdown", 
                                    "nonsig")),
                      "nonsig"))
  res_df$log10FDR[res_df$log10FDR> maxq] <- maxq
  return(res_df)
}) %>% do.call(rbind, .)

pdf(file.path(outdir, "degs", "rdc_mdc.volcano.pdf"), height = 4)
ggplot(as.data.frame(gg_volcanos), aes(x=logFC, y=log10FDR, col=dir)) +
  facet_wrap(vars(id), ncol=2) +
  ggrastr::rasterise(geom_point(alpha=0.5)) +
  scale_color_manual(values=cols) +
  theme_classic() +
  xlim(-20,20) +
  geom_hline(yintercept = (-1*log10(qval_lim)), linetype = "dashed") +
  geom_vline(xintercept = c((-1*lfc_lim), lfc_lim), linetype = "dashed") +
  xlab("log2(Fold Change)") + ylab("-log10(q)") + ggtitle("Interaction") +
  theme(legend.position='none')
dev.off()
