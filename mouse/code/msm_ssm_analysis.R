## Sara's MSM_SSM project
library(WGCNA)
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

# Create ENZ -> SYMBOL mapping key
genome_gse <- org.Mm.eg.db
txby <- keys(genome_gse, 'ENSEMBL')
gene_ids <- mapIds(genome_gse, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")
ens2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                         keytype='ENSEMBL', multiVals="first")

txby <- keys(genome_gse, 'SYMBOL')
ens_ids <- mapIds(genome_gse, keys=txby, column='ENSEMBL',
                   keytype='SYMBOL', multiVals="first")
sym2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                  keytype='SYMBOL', multiVals="first")

# Params
pdir <- "/cluster/projects/mcgahalab/data/mcgahalab/sara_MSM_SSM/results"
dedir <- file.path(pdir, "diffexp")
outdir <- file.path(pdir, "manual")
dir.create(outdir, recursive = F, showWarnings = F)
dir.create(file.path(outdir, "rds"), showWarnings = F)
setwd(pdir)

celltypes <- c('MSM', 'SSM')
celltypes <- setNames(celltypes,celltypes)

###################
#### Functions ####
# Takes a DESeq results data frame, orders the data based on padjusted values,
# splits it based on Log2FC greater-than or less-than 0, then selects the top N
# genes in both directions.
getBidirSigGenes <- function(res, fc_col='log2FoldChange', padj_col='padj',
                             qthresh=0.05, topn=25){
  res <- res[order(res[,padj_col[1]]),]
  resl <- lapply(split(res, f=res[,fc_col]>0), function(i) {
    # Select all genes below the q-cutoff threshold
    iq <- if(length(padj_col)>1){
      rowSums(sapply(padj_col, function(p_id) i[,p_id] < qthresh), na.rm=T) == length(padj_col)
    } else {
      i[,padj_col] < qthresh
    }
    if(sum(iq,na.rm=T) > 0){
      i <- i[which(iq),]
    } else {
      i <- i[0,]
    }
    # Select the top X genes under the q-cutoff
    if(nrow(i) > topn) i <- i[1:topn,]
    return(i)
  })
}

# Takes a DESeq expression matrix (EnsGenes x Samples), as well as a set of genes.
# It will subset the data for the genes that overlap, and create a pheatmap using
# the metafile provided for annotation
degHeatmap <- function(exprmat, genes, meta_df=NULL, genes_style='ENS', 
                       samples=NULL,
                       sample_order=NULL, cluster_rows=FALSE, 
                       cluster_cols=FALSE, title='', scale_score=TRUE,
                       scale_time='pre', max_rownames=30, ...){
  if(any(duplicated(genes))) genes <- genes[-which(duplicated(genes))]
  ens_genes <- if(genes_style=='SYMBOL') ens_ids[genes] else genes
  sym_genes <- if(genes_style=='SYMBOL') genes else gene_ids[genes]
  
  # Keep only genes found in Gene-set and expression-matrix
  gene_m_idx <- which(ens_genes %in% rownames(exprmat))
  ens_genes <- ens_genes[gene_m_idx]
  sym_genes <- sym_genes[gene_m_idx]
  
  non_na_idx <- !is.na(ens_genes)
  ens_genes <- ens_genes[which(non_na_idx)]
  sym_genes <- sym_genes[which(non_na_idx)]
  if(any(is.na(sym_genes))){
    idx <- which(is.na(sym_genes))
    sym_genes[idx] <- ens_genes[idx]
  }
  
  # Create gene expression matrix subset
  if(scale_time == 'pre'){
    print('scaling prior to subsetting')
    gene_heatmap0 <- if(scale_score) as.data.frame(t(apply(exprmat[ens_genes,,drop=F], 1, scale))) else exprmat[ens_genes,,drop=F]
    colnames(gene_heatmap0) <- colnames(exprmat)
    gene_heatmap0 <- if(is.null(samples)) gene_heatmap0 else gene_heatmap0[,samples,drop=F]
    
    if(is.null(sample_order)) sample_order <- c(1:ncol(exprmat))
    gene_heatmap2 <- gene_heatmap0[,sample_order,drop=F]
    colnames(gene_heatmap2) <- colnames(gene_heatmap0)[sample_order]
    
  } else {
    print('scaling after subsetting')
    if(is.null(sample_order)) sample_order <- c(1:ncol(exprmat))
    gene_heatmap <- exprmat[ens_genes,sample_order]
    gene_heatmap2 <- if(scale_score) as.data.frame(t(apply(gene_heatmap, 1, scale))) else gene_heatmap
    colnames(gene_heatmap2) <- colnames(gene_heatmap)
  }
  rownames(gene_heatmap2) <- sym_genes
  
  if(any(is.na(gene_heatmap2))) gene_heatmap2[is.na(gene_heatmap2)] <- 0
  if(!is.null(meta_df)){
    meta_df <- meta_df[sample_order,,drop=F]
  }
  phm <- pheatmap(gene_heatmap2, 
                  cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                  show_rownames=if(nrow(gene_heatmap2) > max_rownames) F else T, 
                  show_colnames=F,
                  annotation_col=if(!is.null(meta_df)) meta_df else NULL,
                  main=title,
                  ...)
  return(phm)
}

# Functions to run the GSEA function and do basic visualization of:
#  - Barplot of top N significant genesets in either direction where X-axis == NES
#  - Heatmap of LFC for gene-sets by genes
getGSEA <- function(geneset_df, lfc_v, return.df=FALSE){
  gsea_out <- tryCatch({
    GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = geneset_df, pvalueCutoff = 1)
  })
  gsea_df <- if(return.df) as.data.frame(gsea_out)[,1:10] else gsea_out
  return(gsea_df)
}
vizGSEA.bp <- function(gsea_df, topn=15, title=NULL){
  # Barplot of NES for the GSEA terms, ordered by NES
  b <- c(0, 0.05, 0.1)
  if(!any(gsea_df$p.adjust < 0.1)){
    return(NULL)
  }
  gsea_df %>%
    filter(p.adjust < 0.1) %>%
    mutate(direction=NES>0) %>%
    group_by(direction) %>%
    slice(1:topn) %>%
    enrichplot:::barplot.enrichResult(x='NES', showCategory=(topn*2)) +
    scale_fill_gradientn(limits = c(min(b),max(b)),
                         colours=colorRampPalette(c("red", "orange", "yellow"))(7),
                         breaks=b, labels=format(b)) +
    ggtitle(title) +
    theme(text = element_text(size=8),
          axis.text.y = element_text(size=8))  +
    xlim(-6,6) +
    xlab("NES") + theme_bw()
}
vizGSEA.hm <- function(gsea, lfc_v, topn=15){
  # Heatmap of LFC for the gene-set
  gsea_x <- setReadable(gsea, org.Mm.eg.db, 'ENTREZID')
  colors <- c('navyblue', 'lightgrey', 'darkred')
  maxval <- max(abs(lfc_v))
  b <- c(-ceiling(maxval), 0, ceiling(maxval))
  
  gsea_x %>%
    heatplot(showCategory=topn, foldChange=lfc_v) +
    theme_bw() +
    scale_fill_gradientn(limits = c(min(b),max(b)), colors = colors,
                         breaks = b, labels = format(b)) +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
}

#####################################################################
#### 1) Comparison of LFC between LN (cis-PBS) to TDLN (cis-PBS) ####
# Plot LFC of [LN Cis-to-PBS] compared to LFC of [TDLN Cis-to-PBS]
ggs <- lapply(c("msm", "ssm"), function(ttype){
  tdln <- read.table(file.path(dedir, paste0(ttype, "-tdln-cis-vs-pbs.diffexp.tsv")),
                     header=T, check.names = F, stringsAsFactors = F)
  ln <- read.table(file.path(dedir, paste0(ttype, "-ln-cis-vs-pbs.diffexp.tsv")),
                   header=T, check.names = F, stringsAsFactors = F)
  cols <- c('gene', 'baseMean', 'symbol', 'log2FoldChange', 'lfcSE')
  m_ln <- merge(tdln[,cols], ln[,cols], by="gene", all=T)
  colnames(m_ln) <- c('gene', 'tdln_base', 'tdln_symbol', 'tdln_lfc', 'tdln_lfcSE',
                      'ln_base', 'ln_symbol', 'ln_lfc', 'ln_lfcSE')
  m_ln$avg_base <- log10(rowMeans(m_ln[,c('tdln_base', 'ln_base')])+1)
  m_ln0 <- m_ln[which(m_ln$ln_base > 50),]
  m_ln0 <- m_ln0[which(m_ln0$ln_base < quantile(m_ln0$ln_base, 0.99)),]
  with(m_ln0, cor(tdln_lfc, ln_lfc), method='spearman')

  gg<- ggplot(m_ln0, aes(x=ln_lfc, y=tdln_lfc, color=avg_base)) +
    geom_point() +
    ggtitle(ttype) +
    theme_minimal()
  return(gg)
})
pdf(file.path(outdir, "LFC_ln-vs-tdln.pdf"))
lapply(ggs, function(i) i)
dev.off()

######################################################################
#### 2) Differential expression controlling for interaction terms ####
cntsdir <- file.path(pdir, "counts")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(file.path(cntsdir, "all.tsv"), header=TRUE, 
                  row.names="gene", check.names=FALSE)
coldata <- read.table(file.path(pdir, "..", 'config', 'samples.tsv'), 
                      header=TRUE, row.names="sample_name", check.names=FALSE)
if(any(colnames(cts) != rownames(coldata))){
  # ensure coldata and cts are in the same sample order
  coldata <- coldata[match(colnames(cts), rownames(coldata)),,drop=FALSE]
}

coldata$lnstatus <-factor(gsub("-.*", "", coldata$condition), c("LN", "TDLN"))
coldata$treatment <- factor(gsub(".*-(.*)_.*", "\\1", coldata$condition), c('PBS', 'CIS'))
coldata$celltype <- gsub(".*_", "", coldata$condition)
coldata_l <- c(split(coldata, coldata$celltype), 
               list("All"=coldata))

celltypes <- setNames(names(coldata_l), names(coldata_l))
resl <- lapply(celltypes, function(celltype){
  ## INTERACTION TERMS
  # if the log2 fold change attributable to a given condition is different 
  # based on another factor, for example if the treatment effect differs 
  # across lnstatus
  coldata_l[[celltype]]$treatment <- relevel(coldata_l[[celltype]]$treatment, "PBS")
  coldata_l[[celltype]]$lnstatus <- relevel(coldata_l[[celltype]]$lnstatus, "LN")
  
  dds <- DESeqDataSetFromMatrix(countData=cts[,rownames(coldata_l[[celltype]])],
                                colData=coldata_l[[celltype]],
                                design=as.formula('~lnstatus+treatment+lnstatus:treatment'))
  
  # remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  # normalization and preprocessing
  dds <- DESeq(dds)
  
  ## MAIN TREATMENT EFFECT
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
  coef <- 'lnstatusTDLN.treatmentCIS'
  interact_reslfc <- lfcShrink(dds, coef=coef, type="apeglm")
  interact_reslfc <- interact_reslfc %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".interaction")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## MAIN TREATMENT EFFECT 1
  # Relates to the treatment effect of CIS vs PBS for TDLN
  # - By setting the reference level to TDLN, we can evaluate CIS vs PBS
  # differences in terms of TDLN status
  coldata_l[[celltype]]$lnstatus <- relevel(coldata_l[[celltype]]$lnstatus, "TDLN")
  dds2 <- DESeqDataSetFromMatrix(countData=cts[,rownames(coldata_l[[celltype]])],
                                 colData=coldata_l[[celltype]],
                                 design=as.formula('~lnstatus+treatment+lnstatus:treatment'))
  
  # remove uninformative columns
  dds2 <- dds2[ rowSums(counts(dds2)) > 1, ]
  dds2 <- DESeq(dds2)
  coef <- 'treatment_CIS_vs_PBS'
  overall_reslfc_tdln <- lfcShrink(dds2, coef=coef, type="apeglm")
  overall_reslfc_tdln <-  overall_reslfc_tdln %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".tdln_cis-pbs")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## MAIN TREATMENT EFFECT 2
  # Relates to the treatment effect of TDLN-CIS vs LN-CIS
  # - By setting the reference level to LN and subsetting for only the CIS samples, 
  # we can evaluate TDLN-CIS vs LN-CIS differences
  coldata_l[[celltype]]$lnstatus <- relevel(coldata_l[[celltype]]$lnstatus, "LN")
  cis_coldata <- coldata_l[[celltype]] %>% 
    filter(treatment=='CIS')
  dds3 <- DESeqDataSetFromMatrix(countData=cts[,rownames(cis_coldata)],
                                 colData=cis_coldata,
                                 design=as.formula('~lnstatus'))
  
  # remove uninformative columns
  dds3 <- dds3[ rowSums(counts(dds3)) > 1, ]
  dds3 <- DESeq(dds3)
  coef <- resultsNames(dds3)[2]
  overall_reslfc_cis <- lfcShrink(dds3, coef=coef, type="apeglm")
  overall_reslfc_cis <- overall_reslfc_cis %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".cis_tdln-ln")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## Identify the differential genes
  res <- Reduce(function(x,y) merge(x,y,by='ens'), list(overall_reslfc_tdln, 
                                                        overall_reslfc_ln, 
                                                        overall_reslfc_cis,
                                                        interact_reslfc))
  res$gene <- gene_ids[res$ens]
  res$gene[is.na(res$gene)] <- res$ens[is.na(res$gene)]
  
  res_sdigits <- res
  for(col_i in colnames(res)){
    if(is.numeric(res[,col_i])){
      res_i <- res[,col_i]
      sdig <- if(mean(abs(res_i), na.rm=T) > 10) 2 else 6
      res_sdigits[,col_i] <- round(res[,col_i], sdig)
    }
  }
  
  write.table(res_sdigits, file=file.path(outdir, paste0(celltype, "_degs.csv")), 
              col.names = T,row.names = F, quote=F, sep=",")
  saveRDS(res, file=file.path(outdir, "rds", paste0(celltype, "_degs.rds")))
  return(list("res"=res, "dds_lnBase"=dds, "dds_tdlnBase"=dds2))
})
saveRDS(resl, file=file.path(outdir, "rds", "resl.rds"))
saveRDS(coldata, file=file.path(outdir, "rds", "coldata.rds"))


### DEMO: Explaining how to read the results from the DEG
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


#################################
#### 3) Gene set of interest ####
resl <- readRDS(file.path(outdir, "rds", "resl.rds"))

st4_inflammatory_goi <- read.csv("/cluster/projects/mcgahalab/data/mcgahalab/sara_MSM_SSM/ref/st4_inflammatory.csv",
                                 sep=",", header = F) %>%
  unlist() %>% as.character() %>% str_to_title()
st4_antiinflammatory_goi <- read.csv("/cluster/projects/mcgahalab/data/mcgahalab/sara_MSM_SSM/ref/st4_antiinflammatory.csv",
                                 sep=",", header = F) %>%
  unlist() %>% as.character() %>% str_to_title()

list2_goi <- read.csv("/cluster/projects/mcgahalab/data/mcgahalab/sara_MSM_SSM/ref/list2.csv",
                      sep=",", header = F) %>%
  unlist() %>% as.character() %>% str_to_title()

inflammatory_goi <- c("CCR2", "IL18", "IL15", "IRF5", "CXCL11", "CXCL10", "IL6", 
                      "IL1B", "IL1A", "PTGES", "IL15RA", "RELA", "TNF", "CD86", "CD80", 
                      "SOCS3", "NLRP3", "SELL", "STAT1", "CD40", "Slc7a2", 
                      "Serpinb2", "Slc7a11", "STAT3", "HIF1A") 
inf2 <- c("Il1rap", "Jun", "Plpp1", "Nos2", "Ptgs2", 
          "Ifng", "Il23a", "Il12a", "Il12b")
inflammatory_goi <- str_to_title(c(inflammatory_goi, inf2))

anti_inflammatory_goi <- c("IL10RA", "MMP13", "MSR1", "CCL17", "CCL22", "IL4RA", 
                           "CLEC7A", "CLEC2i", "clec4a2", "pdgfa", "cd36", "cxcl2", "cxcl3", 
                           "pdgfc", "mmp9", "mrc1", "fpr1", "trem2", "retnla", 
                           "pdgfb", "stab1", "f13a1", "irf4", "stat3", "ccl12", "ccl7", "clec10a", 
                           "mgl2", "stat6", "igf1", "ptgs2", "IL10",  
                           "SOCS1", "GATA3", "IL4", "IL13", "VEGFA", "Ahrr", 
                           "ATF6b", "atf4", "atf1", "atf2", "atf3", "atf4", "atf7", 
                           "arg1", "arg2", "Ahr", "Cyp1a1", "Ido1", "Ido2", 
                           "Ddit3", "il33")
ainf2 <- c('Tgfb1', 'Chil3', 'Tnfsf10', 'Socs2', 'Nr1h3', 'Pparg', 'Retnla', 'cxcl15', 
           'Atf5', 'Atf6', 'cyp1b1', 'Gcn1')
anti_inflammatory_goi <- str_to_title(c(anti_inflammatory_goi, ainf2))

gois <- list("anti-inflammatory"=anti_inflammatory_goi, 
             "inflammatory"=inflammatory_goi,
             "st4_antiinflammatory_goi"=st4_antiinflammatory_goi,
             "st4_inflammatory_goi"=st4_inflammatory_goi,
             "list2"=list2_goi)

inflammatory_goi[which(!inflammatory_goi %in% gene_ids)]
anti_inflammatory_goi[which(!anti_inflammatory_goi %in% gene_ids)]
list2_goi[which(!list2_goi %in% gene_ids)]
lapply(names(resl), function(ct){
  lapply(names(gois), function(goi_id){
    res <- resl[[ct]]$res
    res_goi <- res[res$gene %in% gois[[goi_id]],]
    write.table(res_goi, 
                file=file.path(outdir, paste0(ct, "_degs_", goi_id, ".csv")), 
                col.names = T, row.names = F, quote=F, sep=",")
  })
})

saveRDS(gois,
        file=file.path(outdir, "goi.rds"))

##################################
#### 4 vizDEG: Visualize DEGs ####
celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
gois <- readRDS(file.path(outdir, "goi.rds"))

## Params
top_genes <- 50
print_n <- 1000
q_threshold <- 0.05
get_interaction=TRUE
get_delta_tdln=TRUE
subset_sig_goi=FALSE # Subset the gois for just the genes that are significant in our deg analysis
scale_time='post' # z-scale the gene expression 'pre' or 'post' subsetting for samples to plot

## Create heatmaps of DEGs, based on the interaction between:
# 1) DEG of CIS: TDLN-LN
# 2) Genes with significant differences based on the treatment effect
annotation_colors <- list("treatment"=c("PBS"="black", "CIS"="red"))
phms <- lapply(names(res_dds), function(celltype, get_delta_tdln, get_interaction){
  print(celltype)
  res <- res_dds[[celltype]]$res        # DESeq results data frame
  dds <- res_dds[[celltype]]$dds_lnBase # DESeq object
  
  # Get cnts table 
  cnts0 <- counts(dds)
  vsd <- vst(dds)
  cnts <- counts(dds,normalized=TRUE)
  meta_df <- as.data.frame(colData(dds)[,c("lnstatus","treatment", "celltype")])
  order_idx <- with(meta_df, order(celltype, lnstatus, treatment))
  meta_df_tdln <- split(meta_df, meta_df$lnstatus)$TDLN
  order_idx_tdln <- with(meta_df_tdln, order(celltype, lnstatus, treatment))

  # ((t(cnts0[ens_ids[c('Il10ra', 'Rela', 'Il12b')],]) / colSums(cnts0)) * 1e6)
  # t(cnts[ens_ids[c('Il10ra', 'Rela', 'Il12b')],])
  
  # Subset gois if interested in only genes that are significant according to 
  # our DEG analysis
  if(subset_sig_goi){
    res_sig <- res[which(with(res, (`padj.interaction` < q_threshold) & 
                                (`padj.cis_tdln-ln` < q_threshold))),]
    gois_sub <- lapply(names(gois), function(goi){
      gois[[goi]][which(gois[[goi]] %in% res_sig$gene)]
    })
    names(gois_sub) <- names(gois)
  }else {
    gois_sub <- gois
  }
  
  
  phm_tops <- list()
  # Get top X significant genes from DESeq results; split based on direction
  if(get_interaction){
    ## Top X genes based on the interaction term with LN-PBS as the base
    resl <- getBidirSigGenes(res, fc_col=c('log2FoldChange.interaction'),
                             padj_col=c('padj.interaction', 
                                        'padj.cis_tdln-ln'),
                             qthresh=q_threshold, topn=top_genes)
    genes <- unlist(lapply(resl, function(i) i$ens))
    symbols <- unlist(lapply(resl, function(i) i$gene))
    if(any(is.na(symbols))) symbols[is.na(symbols)] <- genes[is.na(symbols)]
    
    ## DEG Heatmap for top 40 genes in either direction
    exprmat_fc <- res %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var='ens') %>%
      select(grep("log2Fold.*\\.", colnames(res), value=T) %>% 
               grep("interaction", ., value=T,invert=T)) %>%
      rename_with(~ gsub("log2.*\\.", "", .))
    
    meta_df_fc <- data.frame("celltype"=rep(celltype, ncol(exprmat_fc)),
                             "group"=colnames(exprmat_fc), 
                             row.names = colnames(exprmat_fc))
    phm_tops[['interaction_sample']] <- degHeatmap(
      exprmat=assay(vsd), samples=rownames(meta_df_tdln), meta_df=meta_df_tdln,  
      genes=genes, genes_style='ENS', sample_order=order_idx_tdln, 
      title=paste0("CIS(TDLN-LN) & Treatment: ", celltype, "-  top_", top_genes),
      annotation_colors=annotation_colors, max_rownames = print_n,
      scale_time=scale_time
    )
    phm_tops[['interaction_fc']] <- degHeatmap(
      exprmat=exprmat_fc, meta_df=meta_df_fc,  
      genes=genes, genes_style='ENS', sample_order=c(1:ncol(exprmat_fc)), 
      title=paste0("CIS(TDLN-LN) & Treatment: ", celltype, "-  top_", top_genes),
      scale_score = FALSE, annotation_colors=annotation_colors, 
      max_rownames = print_n, scale_time=scale_time
    )
  }
  
  if(get_delta_tdln){
    ## Top X genes based on differnetially expressed between TDLN-CIS to TDLN-PBS
    res2 <- res[which(res$padj.int < 0.05),]
    resl <- getBidirSigGenes(res2, fc_col='log2FoldChange.cis_tdln-ln', 
                             padj_col='padj.cis_tdln-ln', qthresh=q_threshold, 
                             topn=top_genes)
    genes <- unlist(lapply(resl, function(i) i$ens))
    symbols <- unlist(lapply(resl, function(i) i$gene))
    if(any(is.na(symbols))) symbols[is.na(symbols)] <- genes[is.na(symbols)]
    
    ## DEG Heatmap for top 40 genes in either direction
    phm_tops[['cis_tdln-ln']] <- degHeatmap(
      exprmat=assay(vsd), samples=rownames(meta_df_tdln), 
      meta_df=meta_df_tdln,  genes=genes, 
      genes_style='ENS', sample_order=order_idx_tdln, 
      title=paste0("CIS_TDLN-LN: ", celltype, "- top_", top_genes),
      annotation_colors=annotation_colors, max_rownames = print_n, 
      scale_time=scale_time)
  }
  
  ## DEG Heatmap for geneset of interest
  phm_gois <- lapply(names(gois_sub), function(goi){
    if(length(gois_sub[[goi]]) <= 1) return(NULL)
    degHeatmap(exprmat=assay(vsd), samples=rownames(meta_df_tdln), 
               meta_df=meta_df_tdln,  sample_order=order_idx_tdln, 
               genes=gois_sub[[goi]], genes_style='SYMBOL', 
               cluster_rows=F, title=paste0(celltype, ": ", goi),
               annotation_colors=annotation_colors, max_rownames = print_n, 
               scale_time=scale_time)
  })
  ## DEG Heatmap for top 40 genes in either direction
  if(get_interaction){
    
    phm_gois_interaction <- lapply(names(gois_sub), function(goi){
      degHeatmap(
        exprmat=exprmat_fc, meta_df=meta_df_fc,  
        genes=ens_ids[gois_sub[[goi]]], genes_style='ENS', sample_order=c(1:ncol(exprmat_fc)), 
        cluster_rows=F, scale_score = FALSE,
        title=paste0(celltype, ": interaction - ", goi),
        annotation_colors=annotation_colors, max_rownames = print_n, 
        scale_time=scale_time
      )
    })
    phm_gois <- c(phm_gois, phm_gois_interaction)
  }
  
  return(list("res"=do.call(rbind, resl), "phm_tops"=phm_tops, "phms_goi"=phm_gois))
}, get_interaction=get_interaction, get_delta_tdln=get_delta_tdln)
names(phms) <- names(res_dds)

pdf(file.path(outdir, paste0("heatmap_top", top_genes, ".subset_", subset_sig_goi, ".pdf")), 
    width = 12, height = 20)
phm_tops <- unlist(lapply(phms, function(i) i$phm_tops), recursive = F)
lapply(phm_tops, function(i){
  grid::grid.newpage()
  grid::grid.draw(i)
})

phms_goi <- unlist(lapply(phms, function(i) i$phms_goi), recursive = F)
lapply(phms_goi, function(i){
  grid::grid.newpage()
  grid::grid.draw(i)
})
dev.off()

##############################################
#### 5) GSEA: enrichment analysis on DEGs ####
celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
gois <- readRDS(file.path(outdir, "goi.rds"))

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                  # 'C7'=list('IMMUNESIGDB'),             # immunologic signature gene sets
                  'C8'=list(NULL),                      # cell type signature gene sets
                  'custom'=list('custom'=gois))           # custom genesets

ct_res_gras <- lapply(names(res_dds), function(celltype){
  print(celltype)
  res <- res_dds[[celltype]]$res        # DESeq results data frame

  resl <- list()
  # LFC for interaction terms
  resl[['interaction']] <- list(res=res, 
                                lfc_col="log2FoldChange.int")
  
  res_gras <- lapply(resl, function(res_i){
    res <- res_i$res
    lfc_col <- res_i$lfc_col
    lfc_v <- setNames(res[,lfc_col],
                      ens2entrez_ids[res$ens])
    
    gra <- lapply(names(msig_lvls), function(mlvl){
      sub_gra <- lapply(msig_lvls[[mlvl]], function(sublvl){
        # mlvl <- 'H'; sublvl <- NULL
        # get Msigdb gene sets
        print(paste0(">", mlvl, ":", sublvl, "..."))
        msig_ds <- if(mlvl == 'custom'){
          data.frame(gs_name=gsub("[0-9]*$", "", names(unlist(gois))),
                     entrez_gene=sym2entrez_ids[unlist(gois)]) %>%
            filter(!is.na(entrez_gene))
        } else {
          msigdbr(species = 'Mus musculus', category = mlvl, 
                  subcategory = sublvl) %>%
            dplyr::select(gs_name, entrez_gene) %>%
            as.data.frame()
        }
        
        # run GSEA on the gene set given
        msig_gsea <- getGSEA(msig_ds, lfc_v, return.df=FALSE)
        msig_gsea_df <- as.data.frame(msig_gsea)[,1:10]
        
        # Barplot of NES for significant GSEA terms, ordered by NES
        title <- paste0(mlvl, " [", sublvl, "] : ", celltype)
        gsea_bp <- vizGSEA.bp(msig_gsea_df, topn=25, title = title)
        
        # Heatmap visualizing the LFC of Genesets by Genes
        # gsea_heat <- vizGSEA.hm(msig_gsea, lfc_v=lfc_v, topn=15)
        
        return(list("gsea-viz"=list(gsea_bp), #, gsea_heat), ### dp, gg_gsea, ridge),
                    "gsea-tbl"=msig_gsea_df))
      })
      
      if(mlvl!='custom'){
        sublvls <- as.character(unlist(msig_lvls[[mlvl]]))
        if(length(sublvls)==0) sublvls <- 'NA'
        names(sub_gra) <- sublvls
      } else {
        names(sub_gra) <- names(msig_lvls[[mlvl]])
      }
      return(sub_gra)
    })
    names(gra) <- names(msig_lvls)
    return(gra)
  })
  
  return(res_gras)
})
names(ct_res_gras) <- names(res_dds)
saveRDS(ct_res_gras, file=file.path(outdir, "GSEA-response.rds"))

## Write out GSEA results
gsea_tbls <- lapply(unlist(unlist(unlist(ct_res_gras, recursive = F), recursive = F), recursive=F), 
                    function(i){
  i$`gsea-tbl`
})

drop_cols <- c('Description', 'leading_edge')
gsea_tbl <- do.call(rbind, gsea_tbls) %>%
  mutate("LE.tags"=gsub(".*tags=([0-9]*%).*", "\\1", leading_edge)) %>%
  mutate("LE.list"=gsub(".*list=([0-9]*%).*", "\\1", leading_edge)) %>%
  mutate("LE.signal"=gsub(".*signal=([0-9]*%).*", "\\1", leading_edge)) %>%
  select(-one_of(drop_cols)) %>%
  mutate(grp_id=rep(names(gsea_tbls), sapply(gsea_tbls, nrow)))
write.table(gsea_tbl, file=file.path(outdir, "gsea_tbl.csv"), sep=",", 
            col.names = T, row.names = F, quote = F)

## Visualization of GSEA barplots
pdf(file.path(outdir, "gsea_response.pdf"), width = 12)
gsea_vizs <- lapply(unlist(unlist(unlist(ct_res_gras, recursive = F), recursive = F), recursive=F), 
                    function(i){
                      i$`gsea-viz`
                    })
gsea_vizs
dev.off()

#########################################
#### 6) WGCNA: Coexpression analysis ####
celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))

## A) Setup the filtered TPM matrix and metadata -------------------------------
## Create the DESeq object based on TPM matrix
cts <- read.table(file.path("counts", "all_tpm.tsv"), header=TRUE,
                  row.names="gene", check.names=FALSE, stringsAsFactors = F)
cts_tmp <- apply(cts, 2, as.numeric)
rownames(cts_tmp) <- rownames(cts)
coldata <- readRDS(file=file.path(outdir, "rds", "coldata.rds"))

## Remove low expressed genes
cts_filt <- cts_tmp %>%
  as.data.frame() %>%
  dplyr::filter(rowSums(.) >= 50)


## Create DESeq2 matrix
dds <- DESeqDataSetFromMatrix(countData = cts_filt,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)

## B) Identifying power for network modules ------------------------------------
## Pick a soft power threshold  for linking the Similarity Matrix to the Adjacency Matrix
normalized_array <- t(assay(dds))
storage.mode(normalized_array) <- 'numeric'
sft <- pickSoftThreshold(normalized_array,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

## Computer the network modules and eigengenes from the adjacency matrix using TOM signed
pwr <- 16
cor <- WGCNA::cor
bwnet <- blockwiseModules(normalized_array,
                          maxBlockSize = 5000, 
                          TOMType = "signed",
                          power=pwr,
                          numericLabels = TRUE,
                          randomSeed = 1234)
cor<-stats::cor
module_eigengenes <- bwnet$MEs

# # plot adjacency network (https://www.biostars.org/p/402720/)
# adjacency <- adjacency(normalized_array, power=pwr, type="signed")
# tom <- TOMsimilarity(adjacency)
# tom[tom > 0.1] = 1
# tom[tom != 1] = 0
# network <- graph.adjacency(tom)
# network <- simplify(network)  # removes self-loops
# V(network)$color <- bwnet$colors
# network <- delete.vertices(network, degree(network)==0)
# plot(network, layout=layout.fruchterman.reingold(network), 
#      edge.arrow.size = 0.2)

## C) Associating eigengene modules with significant interaction DEGs ------------------------------------
gene_modules <- lapply(names(res_dds), function(res_id){
  res <- res_dds[[res_id]]$res
  
  # Merge eigengene modules with p-values per MSM interaction
  gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
    dplyr::mutate(module = paste0("ME", module)) %>% 
    mutate("hugo"=gene_ids[gene]) %>%
    mutate("entrez"=ens2entrez_ids[gene]) %>%
    inner_join(res %>%
                 select(c("ens", "log2FoldChange.interaction",
                          "log2FoldChange.cis_tdln-ln",
                          "padj.interaction", "padj.cis_tdln-ln")),
               by=c("gene" = "ens")) 
  
  # Summarize the p-adjusted values per eigengene module
  module_deg <- gene_module_key %>% 
    group_by(module) %>%
    dplyr::summarise(mean_lfc_cis=mean(`log2FoldChange.cis_tdln-ln`, na.rm=T),
                     sd_lfc_cis=sd(`log2FoldChange.cis_tdln-ln`, na.rm=T),
                     mean_lfc_int=mean(`log2FoldChange.interaction`, na.rm=T),
                     sd_lfc_int=sd(`log2FoldChange.interaction`, na.rm=T),
                     q95_cis=quantile(`padj.cis_tdln-ln`, 0.95, na.rm=T),
                     q95_int=quantile(padj.interaction, 0.95, na.rm=T),
                     n=n()) %>%
    arrange(q95_int+q95_cis)
  
  # Subset for the gene_modules that have a Q95 padj value < 0.05
  sig_modules <- gene_module_key %>%
    group_by(module) %>%
    filter(module %in% (module_deg %>% 
                          filter(q95_int < 0.2 & q95_cis < 0.2) %>%
                          select(module) %>% 
                          unlist())) %>%
    group_split()
  names(sig_modules) <- sapply(sig_modules, function(i) unique(i$module))
  
  
  ## D) Overrepresentation analysis of modules -----------------------------------
  oraFun <- function(msig_ds, entrez_genes){
    # overrepresentation analysis
    sig_ora <- tryCatch({
      enricher(gene = na.omit(entrez_genes), TERM2GENE = msig_ds)@result
    }, error=function(e){NULL})
    return(sig_ora)
  }
  
  ## Select the genes of module X
  module_ids <- sapply(sig_modules, function(i) unique(i$module))
  module_oras <- lapply(setNames(module_ids,module_ids), function(me){
    print(paste0(res_id, ": ", me, "..."))
    module_genes <- gene_module_key %>%
      dplyr::filter(module == me)
    
    oras <- iterateMsigdb(species='Mus musculus', fun=oraFun, 
                          entrez_genes=module_genes$entrez)
    oras <- lapply(oras, summarizeOra, keep_genes=TRUE)
    
    # Get the LFC for all the genes composing that geneset
    oras <- lapply(oras, function(ora){
      if(nrow(ora)==0) return(ora)
      
      lfc_df <- sapply(ora$geneID, function(genes){
        lfcs <- genes %>% 
          strsplit(., split="/") %>%
          as.data.frame() %>%
          rename_with(~ "entrez") %>%
          inner_join(gene_module_key,
                     by=c("entrez" = "entrez")) %>%
          summarise(mean_lfc_cis=mean(`log2FoldChange.cis_tdln-ln`, na.rm=T),
                    sd_lfc_cis=sd(`log2FoldChange.cis_tdln-ln`, na.rm=T),
                    mean_lfc_int=mean(`log2FoldChange.interaction`, na.rm=T),
                    sd_lfc_int=sd(`log2FoldChange.interaction`, na.rm=T)) %>%
          as.numeric()
      }) %>% 
        t() %>% 
        as.data.frame() %>%
        rename_with(~ c('mean_lfc_cis', 'sd_lfc_cis', 'mean_lfc_int', 'sd_lfc_int'))
      
      return(as.data.frame(cbind(ora, lfc_df)))
    })
    oras <- plyr::rbind.fill(oras) %>%
      select(-c(geneID))
    return(oras)
  })
  saveRDS(module_oras, file=file.path(outdir, "rds", paste0(res_id, "_coexpression_modules_intAndCis.rds")))
  
  
  ## Get frequent word combinations of Gene sets for each module
  n <- 6
  freqWords <- function(textv, split="_", word_pairs=F, n=5){
    if(length(textv)==0) return('')
    words <- unlist(strsplit(textv, split=split))
    if(word_pairs) words <- paste(words, lead(words), sep=":")
    paste(names(
      head(sort(table(words), decreasing = T), n)), 
      collapse=",")
  }
  freq_word_pairs <- sapply(module_oras, function(i){
    freqWords(textv=i$ID, split="_", word_pairs=TRUE, n=n)
  })
  freq_words <- sapply(module_oras, function(i){
    freqWords(textv=i$ID, split="_", word_pairs=FALSE, n=n)
  })
  
  module_deg_anno <- module_deg %>% 
    filter(module %in% names(sig_modules)) %>%
    mutate("freq_w"=freq_words[names(sig_modules)]) %>%
    mutate("freq_w_pairs"=freq_word_pairs[names(sig_modules)])
  
  module_ora <- do.call(rbind, module_oras) %>%
    as.data.frame() %>%
    mutate("module"=rep(names(module_oras), sapply(module_oras, nrow)))
  
  write.table(module_ora, file=file.path(outdir, paste0(res_id, "coexpression_module_intAndCis.csv")),
              sep=",", quote=F, row.names = F, col.names = T)
  return(list("ora"=module_oras, "modules"=sig_modules,
              "anno"=module_deg_anno))
})
names(gene_modules) <- names(res_dds)


gene_modules_df <- unlist(gene_modules, recursive=F)
gene_modules_df <- gene_modules_df[grep("anno$", names(gene_modules_df))] 
gene_modules_df <- gene_modules_df %>%
  rbind.fill() %>% 
  as.data.frame() %>%
  mutate("celltype"=rep(gsub(".anno", "", names(gene_modules_df)), 
                        sapply(gene_modules_df, nrow)))

pdf(file.path(outdir, "coexpression_ORA_genesets.pdf"), height = 8)
lapply(c("MSM", "SSM"), function(ct){
  gene_modules_df %>% 
    select(c('module', 'mean_lfc_cis', 'mean_lfc_int', 'freq_w_pairs', 'celltype')) %>%
    melt() %>%
    filter(celltype==ct) %>%
    mutate(variable=gsub("mean_lfc_", "", variable) %>%
             gsub("cis", "CIS_TDLN-LN", .) %>%
             gsub("int", "CIS-Treatment", .)) %>%
    mutate(freq_w_pairs=gsub(",", "\n", freq_w_pairs)) %>%
    ggplot(., aes(y=freq_w_pairs, x=value, fill=variable)) +
      geom_bar(position="dodge", stat="identity") + 
      facet_grid(rows=vars(module), scales = "free") +
      xlim(-5,5) +
      theme(axis.text.x = element_text(size=5)) +
      ggtitle(ct) + 
      xlab("mean-LFC") + ylab("") +
      theme_classic()
})
dev.off()

##################################################################################
#### 6. single-sample GSEA on each cell-type using the top Response gene-sets ####
gras <- readRDS(file.path(PDIR, "results", "manual", "GSEA-response.rds"))
deg_res <- readRDS(file=file.path(PDIR, "results", "manual", "deg.rds"))

## Get all significant gene-sets to assess using ssGSEA
q_cutoff <- 0.2
category='H'; subcategory='NA'

gras_unlist <- unlist(lapply(gras, function(i) i[[category]][[subcategory]]), recursive=F)
geneset_tbl <- lapply(gras_unlist, function(gra_i) {
  tbl <- gra_i$`gsea-tbl`
  tbl[which(tbl$p.adjust < q_cutoff),]
})
geneset_tbl <- as.data.frame(do.call(rbind, geneset_tbl))[,c('ID', 'p.adjust', 'NES')]

## Preprocessing counts data for variance stabilization
counts <- assay(vst(dds_main, blind=T))

## Isolating significant genesets and assessing using GSVA
sig_gs <- msigdbr(species = species, category = category) %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame()
sig_gs <- sig_gs[which(sig_gs$gs_name %in% unique(geneset_tbl$ID)),]
sig_gs$ens <- entrez2ens_ids[as.character(sig_gs$entrez_gene)]
sig_ens_gs <- split(setNames(sig_gs$ens, sig_gs$entrez_gene), f=sig_gs$gs_name)


## calculate single-sample enrichment score
ss_method <- 'ssgsea'
gsva_es <- round(gsva(counts, sig_ens_gs, verbose=FALSE, method=ss_method),2)
gsva_es <- as.data.frame(t(gsva_es))
gsva_es <- merge(samples_meta, gsva_es, 
                 by.x='samples', by.y=0, all=T)

## Call the ssGSVA score for the BL and Treated samples and assembled
# a melted data frame for the GSVA scores by samples by treatments
gsva_trx <- split(gsva_es, f=gsva_es$treatment)
gs_score <- lapply(gsva_trx, function(gsva_xi){
  gs_melts <- lapply(colnames(gsva_xi)[-c(1:4)], function(geneset){
    gs_mat <- dcast(gsva_xi[,c('patient', 'celltype', geneset)],
                    celltype~patient)
    rownames(gs_mat) <- gs_mat$celltype
    gs_melt <- melt(t(gs_mat[,-1]))
    gs_melt$geneset <- geneset
    gs_melt
  })
  gs_melt <- as.data.frame(do.call(rbind, gs_melts))
  return(gs_melt)
})
gs_score <- merge(gs_score[[1]], gs_score[[2]], by=c('Var1', 'Var2', 'geneset'))
colnames(gs_score) <- c('patient', 'celltype', 'geneset', names(gsva_trx))
gs_score <- melt(gs_score)
gs_score_anno <- merge(gs_score, meta, by.x='patient', 'TH_ID', all=T)
gs_score_anno_spl <- split(gs_score_anno, f=gs_score_anno$geneset)

## Calculate a summary statistic of ssGSEA change per cell type
# split on geneset
gs_stats <- lapply(names(gs_score_anno_spl), function(gs_id){  
  gsx <- gs_score_anno_spl[[gs_id]]
  # Split the genesets based on cell types
  gsx_resp <- lapply(split(gsx, f=gsx$celltype), function(ct_i) {
    # Split the geneset-celltype based on response variables
    resp_df <- lapply(split(ct_i, f=ct_i$Response), function(gsx_i){
      # Aggregate the BL-Treat for each cell-type > response
      ss_stat <- dcast(data=gsx_i[,c('patient', 'variable', 'value')], 
                       formula=as.formula("patient ~ variable"))
      delta <- ss_stat$treat - ss_stat$BL
      return(round(data.frame("mean"=mean(delta),
                              "sd"=sd(delta),
                              "cov"=cov(ss_stat[,-1])[1,2],
                              "p_corr"=cor(ss_stat[,-1], method='pearson')[1,2],
                              "s_corr"=cor(ss_stat[,-1], method='spearman')[1,2],
                              "mae"=mae(ss_stat[,-1]),
                              "mbe"=mbe(ss_stat[,-1]),
                              "rmse"=rmse(ss_stat[,-1])), 3))
    })
    resp_df <- as.data.frame(do.call(rbind, resp_df))
    resp_df$celltype <- unique(as.character(ct_i$celltype))
    resp_df$response <- rownames(resp_df)
    return(resp_df)
  })
  gsx_resp <- as.data.frame(do.call(rbind, gsx_resp))
  gsx_resp$geneset <- gs_id
  return(gsx_resp)
})
gs_stats_df <- as.data.frame(do.call(rbind,gs_stats))
gs_stats_l <- split(gs_stats_df, f=gs_stats_df$response)
gs_stats_df <- Reduce(function(x,y) merge(x,y,by=c('geneset', 'celltype')), gs_stats_l)
write.table(gs_stats_df, file=file.path(PDIR, "results", "manual", paste0("GSEA_ssgsva.", ss_method, ".tsv")),
            sep=",", row.names = F, col.names = T, quote = F)

gs_stats_df$delta <- apply(gs_stats_df[,grep("mbe", colnames(gs_stats_df))], 1, function(i){
  max(i) - min(i)
})
gs_stats_df[order(gs_stats_df$delta),]  

## Do the plotties
pdf(file=file.path(PDIR, "results", "manual", paste0("GSEA_ssgsva.", ss_method, ".pdf")))
lapply(names(gs_score_anno_spl), function(geneset){
  gs_score_anno <- gs_score_anno_spl[[geneset]]
  lims <- c(min(gs_score_anno$value, na.rm=T), 
            max(gs_score_anno$value, na.rm=T))
  
  ggplot(gs_score_anno, aes(x=variable, y=value, color=patient)) +
    geom_point() +
    facet_grid(rows=vars(Response), cols=vars(celltype), scales = "free") +
    geom_line(aes(group = patient), alpha=0.3) +
    ggtitle(geneset) +
    ylim(lims) +
    theme_classic()
})
dev.off()

# 
# gs_id <- 'HALLMARK_TNFA_SIGNALING_VIA_NFKB'
# ct <- 'CM'
# x <- gs_score_anno_spl[[gs_id]][gs_score_anno_spl[[gs_id]]$Response=='Partial Response',]
# xl <- split(x, f=list(x$variable, x$celltype))
# ids <- c('patient', 'variable', 'value')
# xli <- merge(xl[[paste0("BL.", ct)]][,ids], xl[[paste0("treat.", ct)]][,ids], by='patient')
# xli






#### X #####
