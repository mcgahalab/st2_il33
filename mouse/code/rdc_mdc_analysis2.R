# renv::load("~/downloads/renvs")
## Sara's mouse_DC project
# Comparing rDC and mDC between Cis and PBS treated samples in the WT and Tumor
# library(WGCNA)
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
library(SCENIC)

# SCENIC params
scenic_org_code <- 'mgi'
rcis_db <- '/cluster/projects/mcgahalab/ref/scenic/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

# Params
pdir <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_DCs"
old_mdc_data <- file.path("data_old", "Sara_DC_R15_Raw.csv")
old_rdc_data <- file.path("data_old", "Sara_DC_R9_Raw.csv")
old_cpm <- file.path("data_old", "180221_Sara_RNAseq_CPM_formatted_Andreas_MSM_SSM.csv")
data_dedup <- file.path(file.path("results/star/se/mDC_Tumor_CIS_C-merged/rsem",
                                  "mDC_Tumor_CIS_C.genes.results"))
dedir <- file.path(pdir, "results", "diffexp")
outdir <- file.path(pdir, "results", "manual")
dir.create(outdir, recursive = F, showWarnings = F)
dir.create(file.path(outdir, "rds"), showWarnings = F)
setwd(pdir)

cts_main <- readRDS(file=file.path(pdir, "results", "counts", "cts.rds"))
dds_main <- readRDS(file.path("results", "deseq2", "all.main.rds"))
meta_main <- readRDS(file=file.path(pdir, "results", "metadata", "meta.rds"))
meta <- meta_main[colnames(dds_main),]

#######################
#### Functions Git ####
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/st2_il33/functions/msm_ssm_functions.R")

gm <- geneMap("Mus musculus")

.getDds <- function(cts, meta, formula='~wt.t+treatment+wt.t:treatment',
                    min.cnt=5, min.n=0.1){
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData=meta,
                                design=as.formula(formula))
  dds <- dds[ which(rowSums(counts(dds)>min.cnt) > (ncol(dds)*min.n)),]
  DESeq(dds)
}

#############################################
#### 0. Create DESeq object and Metadata ####
# remove uninformative genes
min.cnt <- 5
min.n <- 0.1

cts <- read.table(file.path(pdir, "results", "counts", "all.tsv"),
                  header = T, check.names = F, stringsAsFactors = F) %>%
  tibble::column_to_rownames(., "gene") %>%
  rename_with(., ~gsub("genesresults", "", .))
meta <- read.table(file.path(pdir, "config", "samples.tsv"),
                  header=T, check.names=F, stringsAsFactors = F)  %>%
  mutate("group"=gsub("^.*_([ABCDE]).*", "\\1", sample_name),
         "modifier"=gsub("^.*_[ABCDE]_?(.*)", "\\1", sample_name),
         "dc"=gsub("_.*", "", sample_name),
         "treatment"=gsub("^.*_", "", condition),
         "wt.t"=gsub("^.*_(.*)_.*", "\\1", condition)) %>%
  tibble::column_to_rownames(., "sample_name")
meta <- meta[colnames(cts),]
meta$batch <- ifelse(meta$modifier == '', 'batch1', 'batch2')
keep.ids <- meta %>% filter(modifier != 'merge')  %>%
  rownames %>% 
  grep("rDC_WT_PBS_C", ., value=T, invert = T)


dds_main <- .getDds(cts[,keep.ids], meta[keep.ids,], formula="~condition")
dds_interaction <- .getDds(cts[,keep.ids], meta[keep.ids,], formula='~wt.t+treatment+wt.t:treatment')

## Get the DeSeq2 object for each celltype and  batch/merged-batch independently
# metax <- meta[grep("rDC_WT_PBS_C", rownames(meta), value=T, invert = T),]
# meta_l <- split(metax, f=list(metax$dc, metax$wt.t))
meta_l <- split(meta, f=meta$dc) # default method
meta_l <- list("all"=meta)
dds_l <- lapply(names(meta_l), function(dc.id){
  ids_l <- split(meta_l[[dc.id]], meta_l[[dc.id]]$modifier)
  if(any(names(ids_l) == 'merge')){
    samples <- gsub("^(.*)_.*?$", "\\1", rownames(ids_l$merge))
    refidx <- which(sapply(ids_l, function(i) all(i$modifier=='')))
    ref_ids <- ids_l[[refidx]][-match(samples, rownames(ids_l[[refidx]])),]
    ids_l[[refidx]] <- ids_l[[refidx]][match(samples, rownames(ids_l[[refidx]])),]
    ids_l <- lapply(ids_l, function(i) c(rownames(ref_ids), rownames(i)))
    names(ids_l) <- gsub("^$", "batch1", names(ids_l))
  } else {
    ids_l <- list("batch1"=rownames(ids_l[[1]]),
                  "merge"=NA,
                  "reseq"=NA)
  }
  meta$treatment <- factor(meta$treatment)
  meta$treatment <- relevel(meta$treatment, ref='PBS')
  lapply(ids_l, function(ids){
    tryCatch({
      if(length(unique(meta[ids,]$batch)) == 1){
        .getDds(cts[,ids], meta[ids,], formula="~condition")
      } else {
        meta$condition <- factor(meta$condition)
        meta$condition <- relevel(meta$condition, ref='mDC_T_PBS')
        # .getDds(cts[,ids], meta[ids,], formula="~batch + condition")
        .getDds(cts[,ids], meta[ids,], formula="~condition")
      }
    }, error=function(e){NULL})
  })
}) %>% setNames(., names(meta_l))



dir.create(file.path(pdir, "results", "deseq2"))
dir.create(file.path(pdir, "results", "metadata"))
saveRDS(cts, file=file.path(pdir, "results", "counts", "cts.rds"))
saveRDS(dds_main, file=file.path(pdir, "results", "deseq2", "all.main.rds"))
saveRDS(dds_interaction, file=file.path(pdir, "results", "deseq2", "all.interaction.rds"))
saveRDS(meta, file=file.path(pdir, "results", "metadata", "meta.rds"))
saveRDS(dds_l, file=file.path(pdir, "results", "deseq2", "all.batches.rds"))

############################
#### 1.a) PCA analysis  ####
dds <- dds_main
max_pc <- 10
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
  left_join(meta %>%
              tibble::rownames_to_column(., "samples"), 
            by=c('Name'='samples')) %>%
  mutate(wtt_treatment=paste0(wt.t, "_", treatment),
         group_mod=paste0(group, ".", modifier))


pcs <- paste0("PC", c(1:max_pc_i))
pcs_l <- lapply(c(2:length(pcs)), function(i){c(pcs[[i-1]], pcs[[i]])})

ggps <- lapply(pcs_l, function(pc){
  ggplot(data=pca_y, aes_string(x=pc[1], y=pc[2], shape='wtt_treatment',
                                color='wtt_treatment', fill='wt.t',
                                label='group_mod')) +
    facet_grid(.~dc, scales = 'free') +
    geom_point() + 
    ggrepel::geom_text_repel() +
    geom_line(aes(group = wtt_treatment), alpha=0.3)  +
    ggtitle(paste(pc, collapse="-"))
})
pdf(file.path(pdir, "results", "pca", "pca_batch12.pdf"),#"pca_filtered.pdf"),
    width = 12, height = 7)
ggps
dev.off()
file.copy(file.path(pdir, "results", "pca", "pca_batch12.pdf"), to="~/xfer", overwrite = T)


#############################################
#### 2. Differential Expression analysis ####
lfc_lim <- 1
qval_lim <- 0.01
maxq <- 30
cols <- c("sigup"="#fc8d59", "sigdown"="#91bfdb")

#--- a) Comparison for each batch/merge ----
dds_l <- readRDS(file=file.path(pdir, "results", "deseq2", "all.batches.rds"))

all_volanco_plots <- lapply(names(dds_l[[1]]), function(group.id){
  deg_res <- lapply(names(dds_l), function(dc.id){
    print(paste0("> ", group.id, " - ", dc.id))
    dds <- dds_l[[dc.id]][[group.id]]
    if(is.null(dds)) return(NULL)
    coef <- resultsNames(dds)[length(resultsNames(dds))]# 'wt.tWT.treatmentPBS'
    coef <- grep("mDC_T_CIS", resultsNames(dds), value=T)
    overall_reslfc_ln <- lfcShrink(dds, coef=coef, type="apeglm")
    overall_reslfc_ln <- overall_reslfc_ln %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='ens') %>%
      mutate(symbol=gm$ENSEMBL$SYMBOL[ens])
    
    
    res_df <- overall_reslfc_ln %>% 
      dplyr::select(ens, symbol, `log2FoldChange`, `padj`) %>%
      rename_with(., ~c('ens', 'gene', 'logFC', 'FDR')) %>%
      mutate(id=paste0(dc.id, ".", group.id),
             log10FDR=(-1*log10(FDR)),
             dir=ifelse(FDR<=qval_lim, 
                        ifelse(logFC>=lfc_lim, 
                               "sigup", 
                               ifelse(logFC<=(-1*lfc_lim),
                                      "sigdown", 
                                      "nonsig")),
                        "nonsig"),
             test=coef)
    res_df$log10FDR[res_df$log10FDR> maxq] <- maxq
    return(res_df)
  }) %>% do.call(rbind, .)
  
  if(is.null(deg_res)){
    ggp <- NULL
  } else {
    ggp <- ggplot(as.data.frame(deg_res), aes(x=logFC, y=log10FDR, col=dir)) +
      facet_wrap(vars(id), ncol=2) +
      ggrastr::rasterise(geom_point(alpha=0.5)) +
      scale_color_manual(values=cols) +
      theme_classic() +
      xlim(-20,20) +
      geom_hline(yintercept = (-1*log10(qval_lim)), linetype = "dashed") +
      geom_vline(xintercept = c((-1*lfc_lim), lfc_lim), linetype = "dashed") +
      xlab("log2(Fold Change)") + ylab("-log10(q)") + 
      ggtitle(paste0(group.id, ":", unique(deg_res$test))) +
      theme(legend.position='none')
  }
  return(list("gg"=ggp, "res"=deg_res))
})
names(all_volanco_plots) <- names(dds_l[[1]])

pdf(file.path(pdir, "results", "diffexp", "volcano_interaction_batchControlCondition.pdf"),
    width = 10, height = 5)
all_volanco_plots
dev.off()
file.copy(file.path(pdir, "results", "diffexp", "volcano_interaction_batchControlCondition.pdf"),
          to="~/xfer", overwrite = T)

#--- b) Comparison batch1 and batch2 data ---- 
keep.ids <- colnames(dds_main)
cts <- cts_main
meta <- meta_main
meta$condition_modifier <- with(meta, paste0(condition, modifier))

### Compare the samples from Batch2 vs Batch1 
meta$condition_modifier <- factor(meta$condition_modifier)
reflvls <- c('rDC_WT_PBSreseq', 'mDC_T_PBSreseq')
res_all <- lapply(reflvls, function(reflvl_i){
  print(reflvl_i)
  altlvl_i <- gsub("reseq", "", reflvl_i)
  meta$condition_modifier <- relevel(meta$condition_modifier, ref=reflvl_i)
  dds <- .getDds(cts[,keep.ids], meta[keep.ids,], formula='~condition_modifier')
  
  
  coef <- grep(paste0(altlvl_i, "_"), resultsNames(dds), value=T)
  overall_reslfc_ln <- lfcShrink(dds, coef=coef, type="apeglm") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var='ens') %>%
    mutate(symbol=gm$ENSEMBL$SYMBOL[ens])
  
  
  res_df <- overall_reslfc_ln %>% 
    dplyr::select(ens, symbol, `log2FoldChange`, `padj`) %>%
    rename_with(., ~c('ens', 'gene', 'logFC', 'FDR')) %>%
    mutate(id=paste0(dc.id, ".", group.id),
           log10FDR=(-1*log10(FDR)),
           dir=ifelse(FDR<=qval_lim, 
                      ifelse(logFC>=lfc_lim, 
                             "sigup", 
                             ifelse(logFC<=(-1*lfc_lim),
                                    "sigdown", 
                                    "nonsig")),
                      "nonsig"),
           test=coef)
  res_df$log10FDR[res_df$log10FDR> maxq] <- maxq
  return(list("res"=res_df, "dds"=dds))
}) %>% setNames(., reflvls)


### volcano plots of all the comparisons made
gg_volcanos <- lapply(names(res_all), function(resid){
  deg_res <- res_all[[resid]]
  ggp <- ggplot(as.data.frame(deg_res$res), aes(x=logFC, y=log10FDR, col=dir)) +
    # facet_wrap(vars(id), ncol=2) +
    ggrastr::rasterise(geom_point(alpha=0.9), dpi = 300) +
    scale_color_manual(values=cols) +
    theme_classic() +
    xlim(-10,10) +
    geom_hline(yintercept = (-1*log10(qval_lim)), linetype = "dashed") +
    geom_vline(xintercept = c((-1*lfc_lim), lfc_lim), linetype = "dashed") +
    xlab("log2(Fold Change)") + ylab("-log10(q)") + 
    labs(subtitle = gsub("^.*?_", "", unique(deg_res$res$test)) %>%
           gsub("_vs_", " vs ", .)) +
    theme(legend.position='none')
  ggp
})

pdf(file.path(pdir, "results", "diffexp", "volcano_batch1vs2.pdf"),
    width = 6, height = 5)
gg_volcanos
dev.off()
file.copy(file.path(pdir, "results", "diffexp", "volcano_batch1vs2.pdf"),
          to="~/xfer", overwrite = T)



#--- c) Comparison for filtered data ---- 
keep.ids <- colnames(dds_main)
cts <- cts_main
meta <- meta_main

### Compare all [CIS vs PBS] for each DC and treatment group
meta$condition <- factor(meta$condition)
reflvls <- grep("PBS", levels(meta$condition), value=T)
res_cispbs <- lapply(reflvls, function(reflvl_i){
  print(reflvl_i)
  altlvl_i <- gsub("PBS", "CIS", reflvl_i)
  meta$condition <- relevel(meta$condition, ref=reflvl_i)
  dds <- .getDds(cts[,keep.ids], meta[keep.ids,], formula='~condition')
  
  
  coef <- grep(altlvl_i, resultsNames(dds), value=T)
  overall_reslfc_ln <- lfcShrink(dds, coef=coef, type="apeglm") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var='ens') %>%
    mutate(symbol=gm$ENSEMBL$SYMBOL[ens])
  
  
  res_df <- overall_reslfc_ln %>% 
    dplyr::select(ens, symbol, `log2FoldChange`, `padj`, baseMean, pvalue, lfcSE) %>%
    rename_with(., ~c('ens', 'gene', 'logFC', 'FDR',
                      'baseMean', 'pvalue', 'lfcSE')) %>%
    mutate(id=paste0("ref_", reflvl_i),
           log10FDR=(-1*log10(FDR)),
           dir=ifelse(FDR<=qval_lim, 
                      ifelse(logFC>=lfc_lim, 
                             "sigup", 
                             ifelse(logFC<=(-1*lfc_lim),
                                    "sigdown", 
                                    "nonsig")),
                      "nonsig"),
           test=coef)
  res_df$log10FDR[res_df$log10FDR> maxq] <- maxq
  return(list("res"=res_df, "dds"=dds))
}) %>% setNames(., reflvls)

### Compare the change of [CIS vs PBS] between [TUMOR vs WT] for each DC
meta$dc_wt.t <- factor(with(meta, paste0(dc, "_", wt.t)))
meta$treatment <- factor(meta$treatment)
meta$treatment <- relevel(meta$treatment, ref='PBS')
reflvls <- grep("WT", levels(meta$dc_wt.t), value=T)
res_int <- lapply(reflvls, function(reflvl_i){
  print(reflvl_i)
  altlvl_i <- gsub("WT", "T", reflvl_i)
  meta$dc_wt.t <- relevel(meta$dc_wt.t, ref=reflvl_i)
  dds <- .getDds(cts[,keep.ids], meta[keep.ids,], formula='~dc_wt.t+treatment+dc_wt.t:treatment')
  
  
  coef <- grep(paste0(altlvl_i, "\\."), resultsNames(dds), value=T)
  overall_reslfc_ln <- lfcShrink(dds, coef=coef, type="apeglm") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var='ens') %>%
    mutate(symbol=gm$ENSEMBL$SYMBOL[ens])
  
  
  res_df <- overall_reslfc_ln %>% 
    dplyr::select(ens, symbol, `log2FoldChange`, `padj`, baseMean, pvalue, lfcSE) %>%
    rename_with(., ~c('ens', 'gene', 'logFC', 'FDR',
                      'baseMean', 'pvalue', 'lfcSE')) %>%
    mutate(id=paste0("ref_", reflvl_i),
           log10FDR=(-1*log10(FDR)),
           dir=ifelse(FDR<=qval_lim, 
                      ifelse(logFC>=lfc_lim, 
                             "sigup", 
                             ifelse(logFC<=(-1*lfc_lim),
                                    "sigdown", 
                                    "nonsig")),
                      "nonsig"),
           test=coef)
  res_df$log10FDR[res_df$log10FDR> maxq] <- maxq
  return(list("res"=res_df, "dds"=dds))
}) %>% setNames(., reflvls)
names(res_int) <- paste0("Interaction_", names(res_int))

### Write out DEGs into a table
lapply(names(res_all), function(resid){
  res_all[[resid]]$res %>%
    mutate(biotype=ens2biotype_ids[ens]) %>%
    write.table(., file=file.path(pdir, "results", "diffexp", paste0("deg.", resid, ".csv")),
                sep=",", col.names = T, row.names = F,quote = F)
})

### volcano plots of all the comparisons made
res_all <- c(res_cispbs, res_int)
gg_volcanos <- lapply(names(res_all), function(resid){
  deg_res <- res_all[[resid]]
  uid <- unique(deg_res$res$test)
  if(grepl("condition", uid)){
    main_title <- gsub("condition_", "", uid) %>%
      gsub("_vs.*", "", .) %>%
      gsub("_[a-zA-Z0-9]*$", "", .)
    subtitle <- gsub("^.*?_", "", uid) %>%
      gsub("_vs_", " vs ", .) %>%
      gsub(paste0(main_title, "_"), "", .)
  } else {
    main_title <- 'Interaction'
    subtitle <- gsub("dc_wt.t|treatment", "", uid) %>%
      paste0(., " vs ", gsub("_T", "_WT", .) %>%
               gsub("CIS", "PBS", .))
  }
  n <- reshape2::melt(t(data.frame("sigup"=paste0("Sig.up = ", table(as.data.frame(deg_res$res)$dir)['sigup']),
                                 "sigdown"=paste0("Sig.down = ", table(as.data.frame(deg_res$res)$dir)['sigdown'])))) %>%
    mutate(x=c(10, -10),
           y=c(30, 30),
           cols=c(cols['sigup'], cols['sigdown']),
           hjust=c(1, 0))
  
  ggp <- ggplot(as.data.frame(deg_res$res), aes(x=logFC, y=log10FDR, col=dir)) +
    # facet_wrap(vars(id), ncol=2) +
    ggrastr::rasterise(geom_point(alpha=0.9), dpi = 300) +
    scale_color_manual(values=cols) +
    theme_classic() +
    xlim(-10,10) + ylim(0, 30) +
    geom_hline(yintercept = (-1*log10(qval_lim)), linetype = "dashed") +
    geom_vline(xintercept = c((-1*lfc_lim), lfc_lim), linetype = "dashed") +
    xlab("log2(Fold Change)") + ylab("-log10(q)") + 
    labs(title=main_title,
         subtitle = subtitle) +
    theme(legend.position='none')
  
  ggp + geom_text(data=n, aes(x=x, y=y, col=Var1, label=value, hjust=hjust))
})

pdf(file.path(pdir, "results", "diffexp", "volcano_plots_all.filtered.pdf"),
    width = 12, height = 11)
cowplot::plot_grid(plotlist=gg_volcanos[1:4], ncol=2,  nrow = 2)
cowplot::plot_grid(plotlist=gg_volcanos[5:6], ncol=2, nrow = 2)
dev.off()
file.copy(file.path(pdir, "results", "diffexp", "volcano_plots_all.filtered.pdf"),
          to="~/xfer", overwrite = T)

saveRDS(res_all, file=file.path(pdir, "results", "diffexp", "res_all.rds"))
# res_all <- readRDS(file=file.path(pdir, "results", "diffexp", "res_all.rds"))

#--- d) Overlap DEGs with MSM/SSM data ----
res_all <- readRDS(file=file.path(pdir, "results", "diffexp", "res_all.rds"))

msmdir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_MSM_SSM/results/manual'
res_dds <- readRDS(file.path(msmdir, "rds", "resl.rds"))
sm_groups <- grep("\\.", colnames(res_dds$MSM$res), value=T) %>%
  gsub("^.*\\.", "", .) %>% gsub("tdln_|ln_", "", .) %>%
  grep("tdln-ln", ., invert = T, value=T) %>%
  unique

pcutoff <- 0.01


lapply(sm_groups, function(compgroup){
  .getmdat <- function(res, pattern){
    res %>% 
      select(ens, gene, grep(pattern, colnames(res), value=T)) %>%
      rename_with(., ~gsub(paste0(".", pattern), "", .)) %>%
      rename_with(., ~gsub("log2FoldChange", "logFC", .) %>% gsub("padj", "FDR", .)) %>%
      mutate(dir=ifelse(FDR<=qval_lim, 
                        ifelse(logFC>=lfc_lim, 
                               "sigup", 
                               ifelse(logFC<=(-1*lfc_lim),
                                      "sigdown", 
                                      "nonsig")),
                        "nonsig")) %>%
      filter(dir != 'nonsig')
  }
  
  if(compgroup == 'cis-pbs'){
    res_dat <- list("rDC"=filter(res_all$rDC_T_PBS$res, dir != 'nonsig'),
                    "mDC"=filter(res_all$mDC_T_PBS$res, dir != 'nonsig'),
                    'MSM'=.getmdat(res_dds$MSM$res, 'tdln_cis-pbs'),
                    'SSM'=.getmdat(res_dds$SSM$res, 'tdln_cis-pbs'))
  } else if(compgroup=='interaction') {
    res_dat <- list("rDC"=filter(res_all$Interaction_rDC_WT$res, dir != 'nonsig'),
                    'mDC'=filter(res_all$Interaction_mDC_WT$res, dir != 'nonsig'),
                    'MSM'=.getmdat(res_dds$MSM$res, 'interaction'),
                    'SSM'=.getmdat(res_dds$SSM$res, 'interaction'))
  }
  
  sig_df <- lapply(res_dat, function(i) select(i, c(ens, dir))) %>% 
    purrr::reduce(full_join, by='ens') %>%
    tibble::column_to_rownames(., "ens") %>%
    magrittr::set_colnames(., names(res_dat))
  
  sigdirs <- c('sigup', 'sigdown')
  ms <- lapply(sigdirs, function(sig_i){
    sig_df_sel <- (sig_df==sig_i)
    sig_df_sel[is.na(sig_df_sel)] <- 0
    storage.mode(sig_df_sel) <- 'integer'
    sig_df_sel <- sig_df_sel[which(rowSums(sig_df_sel) > 0),]

    # Make sure all unique combinations are represented
    combinations <- function(size, choose) {
      d <- do.call("expand.grid", rep(list(0:1), size))
      d[rowSums(d) == choose,]
    }
    sig_df_sel <- lapply(c(0:4), combinations, size=4) %>% 
      do.call(rbind,.) %>%
      magrittr::set_colnames(., colnames(sig_df_sel)) %>%
      rbind(sig_df_sel, .)
    m = ComplexHeatmap::make_comb_mat(sig_df_sel, 
                                      mode = "intersect")
    
    newvals <- log10(attr(m, 'comb_size')-0.99)
    newvals[newvals<0] <- 0
    attr(m, 'comb_size') <- newvals
    return(m)
  }) %>% setNames(., sigdirs)
  
  
  combcol_map <- c('macrophage'="#1b9e77", 'intersect'="#636363", 'dc'="#7570b3")
  combcols <- rep("intersect", length(comb_name(ms$sigdown)))
  combcols[grepl("^00[01]*$", comb_name(ms$sigdown))] <- 'macrophage'
  combcols[grepl("^[01][01]00$", comb_name(ms$sigdown))] <- 'dc'
  
  pdf(file.path(pdir, "results", "diffexp", paste0("macro_dc_compSig.", compgroup, ".pdf")))
  p <- ComplexHeatmap::UpSet(
    t(ms$sigdown), 
    set_order = names(res_dat), 
    comb_order = order(combcols, comb_size(ms$sigdown)),
    comb_col = combcol_map[combcols],
    left_annotation = upset_left_annotation(t(ms$sigdown),
                                            gp = gpar(fill = cols['sigdown']),
                                            width=unit(6,'cm')),
    right_annotation = upset_right_annotation(t(ms$sigup),
                                              gp = gpar(fill = cols['sigup']),
                                              width=unit(6,'cm')))
  print(p)
  dev.off()
  file.copy(file.path(pdir, "results", "diffexp", paste0("macro_dc_compSig.", compgroup, ".pdf")), 
            to = "~/xfer", overwrite = T)
  cat(paste0("xfer macro_dc_compSig.", compgroup, ".pdf\n"))
  
})

######################################
#### 4. GSEA analysis of the DEGs ####
res_all <- readRDS(file=file.path(pdir, "results", "diffexp", "res_all.rds"))
min_baseMean <- 5

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets
                  #'C7'=list('IMMUNESIGDB'))             # immunologic signature gene sets

gsea_all <- lapply(names(res_all), function(compid){
  res <- res_all[[compid]]$res
  
  lfc_df <- res %>% 
    filter(baseMean > min_baseMean) %>%
    dplyr::select(grep("ens|logFC", colnames(.), value=T)) %>%
    tibble::column_to_rownames(., "ens")
  
  # get GSEA table
  lfc_v <- setNames(lfc_df$logFC,
                    gm$ENSEMBL$ENTREZ[rownames(lfc_df)])
  names(lfc_v) <- make.unique(names(lfc_v))
  
  ggseas <- lapply(names(msig_lvls), function(mlvl){
    sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
      print(paste0(">", compid, "-", mlvl, ":", sublvl, "..."))
      msig_ds <- msigdbr(species = 'Mus musculus', category = mlvl, subcategory = sublvl) %>%
        dplyr::select(gs_name, entrez_gene) %>%
        as.data.frame()
      
      # overrepresentation analysis
      for(attempt in c(1:5)){
        print(paste0("start - ", attempt))
        msig_gsea <- tryCatch({
          GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = msig_ds, pvalueCutoff = 1)
        }, error=function(e){print("fail"); NULL})
        if(!is.null(msig_gsea)) break()
      }
      return(msig_gsea)
    }) %>% setNames(., unlist(msig_lvls[[mlvl]]))
  }) %>%
    setNames(., names(msig_lvls))
  
  ggseas <- unlist(ggseas, recursive = F)
  gsea_df <- lapply(ggseas, function(i) i@result) %>%
    do.call(rbind, .) %>%
    mutate(qvalue_total=p.adjust(pvalue, method = 'BH'))
  
  gsea_df %>% arrange(qvalue_total) %>%
    select(-core_enrichment)
  return(list("obj"=ggseas, "df"=gsea_df))
})
names(gsea_all) <-  names(res_all)
saveRDS(gsea_all, file=file.path(outdir, 'rds', "gseas_all.rds"))
gsea_all <- readRDS(file=file.path(outdir, 'rds', "gseas_all.rds"))

# Write out GSEA
dir.create(file.path(outdir, "gse", "tables"), showWarnings = F)
lapply(names(gsea_all), function(id){
  df <- gsea_all[[id]]$df
  df <- df %>% 
    select(-Description) %>%
    mutate(leading_edge=gsub(",", "_", leading_edge),
           ID=gsub(",", "_", ID))
  write.table(df, file=file.path(outdir, "gse", "tables", paste0(id, ".csv")),
              sep=",", col.names = T, row.names = F, quote = F)
})
          
# Write out the GSEA to cytoscape format
gs_map <- .mapGsToExactSource(species="Mus musculus")
for(celltype in names(gsea_all)){
  dir.create(file.path(outdir, "gse", "cytoscape"), recursive = T, showWarnings = F)
  # gsea_l <- unlist(gsea_all[[celltype]]$log2FoldChange.interaction[-4], recursive=F)
  # gsea_df <- do.call(rbind, lapply(gsea_l, as.data.frame))
  gsea_df <- gsea_all[[celltype]]$df %>%
    mutate(qvalues=qvalue_total)
  gsea_dat_i <- gsea2CytoscapeFormat(gsea_df, gs_map) 
  write.table(gsea_dat_i$up, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("gsea.", celltype, ".up.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(gsea_dat_i$down, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("gsea.ref", celltype, ".down.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
}


