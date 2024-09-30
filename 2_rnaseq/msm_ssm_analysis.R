## Sara's MSM_SSM project
library(DOSE)
library(extrafont)
library(WGCNA)
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

###############
#### Setup ####
# SCENIC params
scenic_org_code <- 'mgi'
rcis_db <- '/cluster/projects/mcgahalab/ref/scenic/cistarget_db/version1/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'

# Create ENZ -> SYMBOL mapping key
genome_gse <- org.Mm.eg.db
txby <- keys(genome_gse, 'ENSEMBL')
gene_ids <- mapIds(genome_gse, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")
ens2sym_ids <- gene_ids
ens2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                         keytype='ENSEMBL', multiVals="first")
entrez2ens_ids <- setNames(names(ens2entrez_ids), ens2entrez_ids)

txby <- keys(genome_gse, 'SYMBOL')
ens_ids <- mapIds(genome_gse, keys=txby, column='ENSEMBL',
                  keytype='SYMBOL', multiVals="first")
sym2ens_ids <- ens_ids
sym2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                         keytype='SYMBOL', multiVals="first")
entrez2sym_ids <- setNames(names(sym2entrez_ids), sym2entrez_ids)

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))
# GTF_gene <- GTF[which(GTF$type=='gene' & GTF$gene_biotype=='protein_coding'),]
# gene_length <- setNames(width(GTF_gene), GTF_gene$gene_name)


# Params
pdir <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_MSM_SSM/results"
dedir <- file.path(pdir, "diffexp")
outdir <- file.path(pdir, "manual")
dir.create(outdir, recursive = F, showWarnings = F)
dir.create(file.path(outdir, "rds"), showWarnings = F)
setwd(pdir)

delta_col <- c('delta_int'='#f03b20',
               'delta_cis'='#636363')
celltypes <- c('MSM', 'SSM')
celltypes <- setNames(celltypes,celltypes)

#######################
#### Functions Git ####
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/st2_il33/functions/msm_ssm_functions.R")
gm <- geneMap(species='Mus musculus') 

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

##########################################
#### 0) Read in Counts and Other Data ####
  cntsdir <- file.path(pdir, "counts")
  
  # colData and countData must have the same sample order, but this is ensured
  # by the way we create the count matrix
  cts <- read.table(file.path(cntsdir, "msm_ssm.counts.tsv"), header=TRUE, 
                    row.names="gene", check.names=FALSE)
  tpm <- read.table(file.path(cntsdir, "msm_ssm.tpm.tsv"), header=TRUE, 
                    row.names="gene", check.names=FALSE)
  coldata <- read.table(file.path(pdir, "..", 'config', 'msm_ssm.samples.tsv'), 
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
  
#####################################################################
#### 1) Differential expression controlling for interaction terms ####
cntsdir <- file.path(pdir, "counts")

coldata_treat_ln_l <- split(coldata, list(coldata$treatment, coldata$lnstatus))

celltypes <- setNames(names(coldata_l), names(coldata_l))
resl <- lapply(celltypes, function(celltype){
  ## INTERACTION TERMS
  # if the log2 fold change attributable to a given condition is different 
  # based on another factor, for example if the treatment effect differs 
  # across lnstatus
  coldata_l[[celltype]]$treatment <- relevel(coldata_l[[celltype]]$treatment, "PBS")
  coldata_l[[celltype]]$lnstatus <- relevel(coldata_l[[celltype]]$lnstatus, "LN")
  
  dds_raw <- DESeqDataSetFromMatrix(countData=cts[,rownames(coldata_l[[celltype]])],
                                colData=coldata_l[[celltype]],
                                design=as.formula('~lnstatus+treatment+lnstatus:treatment'))
  
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
  coef <- 'lnstatusTDLN.treatmentCIS'
  interact_reslfc <- lfcShrink(dds, coef=coef, type="apeglm")
  interact_reslfc <- interact_reslfc %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".interaction")) %>%
    tibble::rownames_to_column(var='ens') 
  
  ## MAIN TREATMENT EFFECT 2
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
  
  ## MAIN TREATMENT EFFECT 3
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
  
  ## MAIN TREATMENT EFFECT 4
  # Relates to the treatment effect of PBS-TDLN vs PBS-LN
  # - By setting the reference level to LN and subsetting for only the PBS samples, 
  # we can evaluate PBS-TDLN vs PBS-LN differences
  coldata_l[[celltype]]$lnstatus <- relevel(coldata_l[[celltype]]$lnstatus, "LN")
  pbs_coldata <- coldata_l[[celltype]] %>% 
    filter(treatment=='PBS')
  dds4 <- DESeqDataSetFromMatrix(countData=cts[,rownames(pbs_coldata)],
                                 colData=pbs_coldata,
                                 design=as.formula('~lnstatus'))
  
  # remove uninformative columns
  dds4 <- dds4[ rowSums(counts(dds4)) > 1, ]
  dds4 <- DESeq(dds4)
  coef <- resultsNames(dds4)[2]
  overall_reslfc_pbs <- lfcShrink(dds4, coef=coef, type="apeglm")
  overall_reslfc_pbs <- overall_reslfc_pbs %>%
    as.data.frame() %>%
    rename_with(~ paste0(., ".pbs_tdln-ln")) %>%
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
  
  write.table(res_sdigits, file=file.path(outdir, paste0(celltype, "_degs.csv")), 
              col.names = T,row.names = F, quote=F, sep=",")
  saveRDS(res, file=file.path(outdir, "rds", paste0(celltype, "_degs.rds")))
  return(list("res"=res, "dds_lnBase"=dds_raw, "dds_tdlnBase"=dds2))
})
resl_msm_v_ssm <- lapply(coldata_treat_ln_l, function(col_grp){
  col_grp$celltype <- factor(col_grp$celltype)
  col_grp$celltype <- relevel(col_grp$celltype, "SSM")
  
  dds <- DESeqDataSetFromMatrix(countData=cts[,rownames(col_grp)],
                                colData=col_grp,
                                design=as.formula('~celltype'))
  
  # remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  
  ## MAIN TREATMENT EFFECT 1
  # celltype_MSM_vs_SSM
  # Compares the differences between MSM and SSM celltypes under every condition
  lbl <- with(col_grp, paste(lnstatus, treatment, sep="-")) %>% unique
  reslfc_ct <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") %>%
    as.data.frame %>%
    rename_with(., ~ paste0(., ".celltype_", lbl)) %>%
    tibble::rownames_to_column(var='ens') 
  return(list("dds"=dds, "res"=reslfc_ct))
})
res_msm_ssm <- lapply(resl_msm_v_ssm, function(i) i$res) %>% 
  purrr::reduce(., full_join, by='ens') %>% 
  mutate("symbol"=ens2sym_ids[ens]) %>%
  relocate(ens, symbol)
resl[['Celltype']] <- c(lapply(resl_msm_v_ssm, function(i) i$res),
                        list("res"=res_msm_ssm))

saveRDS(resl, file=file.path(outdir, "rds", "resl.rds"))
saveRDS(coldata, file=file.path(outdir, "rds", "coldata.rds"))


### DEMO: Explaining how to read the results from the DEG
resl <- readRDS(file=file.path(outdir, "rds", "resl.rds"))
coldata <- readRDS(file=file.path(outdir, "rds", "coldata.rds"))

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


######################################################
#### 2) Overlap of sigGenes between MSM and SSM ####
celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
qval <- 0.001  # 0.05, 0.01, 0.001
rm_gene_regex <- "^Gm[0-9]+$|^Fam[0-9]|^Olfr|^BC0|Rik$"


# Find the MSM/SSM-only genes and those that intersect based on a sig.qvalue
msm_res <- res_dds$MSM$res %>%
  filter(!grepl(rm_gene_regex, gene))
ssm_res <- res_dds$SSM$res %>%
  filter(!grepl(rm_gene_regex, gene))

.getSigGenes <- function(x){
  x %>% 
    filter(padj.interaction < qval,
           abs(log2FoldChange.interaction) > 2) %>%
    pull(ens)
}
msm_genes <- .getSigGenes(msm_res)
ssm_genes <- .getSigGenes(ssm_res)
all_genes <- unique(sort(c(msm_genes, ssm_genes)))

msm_only_genes <- setdiff(msm_genes, ssm_genes)
ssm_only_genes <- setdiff(ssm_genes, msm_genes)
intersect_genes <- intersect(msm_genes, ssm_genes)

# Create a new dataframe of the interaction MSM/SSM with their labels of overlap
msm_ssm_overlap_res <- msm_res %>% 
  filter(ens %in% all_genes) %>%
  select(ens, baseMean.interaction, log2FoldChange.interaction, 
         pvalue.interaction, padj.interaction) %>%
  left_join(., ssm_res %>% 
              filter(ens %in% all_genes) %>%
              select(ens, baseMean.interaction, log2FoldChange.interaction, 
                     pvalue.interaction, padj.interaction),
            by='ens', suffix=c('.MSM', '.SSM')) %>%
  mutate(gene=gm$ENSEMBL$SYMBOL[ens], 
         biotype=ens2biotype_ids[ens],
         overlap=ifelse(ens %in% msm_only_genes, 'MSM_only',
                        ifelse(ens %in% ssm_only_genes, 'SSM_only', 'intersect'))) %>%
  relocate(c(gene, biotype), .after=ens)

  
write.table(msm_ssm_overlap_res, file=file.path(outdir, paste0("overlap_interaction_degs.q", qval, ".csv")),
            col.names = T,row.names = F, quote=F, sep=",")


##############################################
#### 3) GSEA: enrichment analysis on DEGs ####
celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
gois <- readRDS(file.path(outdir, "goi.rds"))
genes_to_include <- read.table(file.path(pdir, "..", "ref", "msm_ssm_gsea_genes", "msm_ssm.txt"),
                               header = F, check.names = F) %>%
  dplyr::filter(!duplicated(V1))
filter_genes <- FALSE

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF'))#, # ontology gene sets
                  # 'C7'=list('IMMUNESIGDB'))             # immunologic signature gene sets
                  # 'C8'=list(NULL),                      # cell type signature gene sets
                  # 'custom'=list('custom'=gois))           # custom genesets
q_thresh <- 0.2  # Threshold to select top genesets from

lfcs <- list("grp1"=c("interaction"='log2FoldChange.interaction',
                      "cis_tdln_ln"="log2FoldChange.cis_tdln-ln"),
             "grp2"=c("tdln_cis_pbs"="log2FoldChange.tdln_cis-pbs",
                      "pbs_tdln_ln"="log2FoldChange.pbs_tdln-ln"))
lfcs <- list("grp1"=c("interaction"='log2FoldChange.interaction'))
ct_res_gsea <- lapply(lfcs, function(lfc_cols){
  # lfc_cols <- c("interaction"='log2FoldChange.interaction',
  #               "cis_tdln_ln"="log2FoldChange.cis_tdln-ln")
  ids <- gsub("_.*", "", names(lfc_cols))
  
  ct_res_gsea <- lapply(celltypes, function(celltype){
    print(celltype)
    res <- res_dds[[celltype]]$res #%>%
    # mutate(biotype=ens2biotype_ids[ens]) %>%
    # filter(biotype == 'protein_coding') # DESeq results data frame
    if(filter_genes) res <- res %>% filter(ens %in% genes_to_include$V1)
    
    # Run GSEA using the LFC for the interaction and the CIS-treatment terms
    gsea <- lapply(lfc_cols, function(lfc_col){
      lfc_v <- setNames(res[,lfc_col],
                        ens2entrez_ids[res$ens])
      
      gsea_v <- iterateMsigdb(species='Mus musculus', msig_lvls=msig_lvls, 
                              fun=gseaFun, lfc_v=lfc_v)
      # return(list("gsea"=gsea_v, "lfc"=lfc_v))  # for GSEA2cytoscape
      return(gsea_v)  # for standard GSEA
    })
    
    # Merge the GSEA values between the two terms, select values that are 
    # significant below a certain threshold for each term
    gsea <- lapply(gsea, unlist, recursive=F)
    gsea_merge <- lapply(names(gsea[[1]]), function(geneset_lvl){
      int_gsea <- gsea[[names(lfc_cols)[1]]][[geneset_lvl]] # Interaction GSEA
      cis_gsea <- gsea[[names(lfc_cols)[2]]][[geneset_lvl]] # CIS_TDLN-LN GSEA
      
      .isolateGsea <- function(x, id=NULL){
        x %>%
          select(ID, setSize, NES, p.adjust, rank, leading_edge) %>%
          rename_with(., ~ paste0(., ".", id), .cols=matches("[^ID]", perl=T))
      }
      gsea_m <- tryCatch({
        full_join(int_gsea %>% 
                    as.data.frame() %>%
                    .isolateGsea(., ids[1]),
                  cis_gsea %>% 
                    as.data.frame() %>%
                    .isolateGsea(., ids[2]),
                  by='ID') %>% 
          mutate(geneset=geneset_lvl,
                 celltype=celltype)
      }, error=function(e){NULL})
      return(gsea_m)
    })
    names(gsea_merge) <- names(gsea[[1]])
    
    gsea_merge_sig <- lapply(gsea_merge, function(gsea_m){
      tryCatch({
        gsea_m %>%
          filter(!!as.symbol(paste0("p.adjust.", ids[1])) < q_thresh & 
                   !!as.symbol(paste0("p.adjust.", ids[2])) < q_thresh)
      }, error=function(e){NULL})
    })
    
    return(list("gsea"=gsea_merge, "gsea_sig"=gsea_merge_sig))
  })
})
saveRDS(ct_res_gsea, file=file.path(outdir, "GSEA-response.allGenes.2cyto.rds"))
ct_res_gsea_bkup <- readRDS(file.path(outdir, "GSEA-response.allGenes.2cyto.rds"))

gs_map <- .mapGsToExactSource(species="Mus musculus")
for(celltype in names(ct_res_gsea_bkup)){
  print(celltype)
  dir.create(file.path(outdir, "gse", "cytoscape"), recursive = T, showWarnings = F)
  
  ## Format the GSEA Ranks
  ranks <- ct_res_gsea_bkup[[celltype]]$interaction$lfc
  ranks <- data.frame(ranks) %>%
    dplyr::rename_with(., ~"rank") %>%
    mutate(GeneName=gm$ENTREZID$SYMBOL[names(ranks)]) %>% 
    relocate(GeneName) %>% 
    filter(!duplicated(GeneName),
           !is.na(GeneName)) %>%
    arrange(desc(rank))
  
  ## Format the TPM for the included samples
  coldata_i <- coldata %>% filter(celltype==!!celltype)
  tpm_i <- tpm %>% 
    select(rownames(coldata_i)) %>%
    tibble::rownames_to_column("ens") %>% 
    mutate(GeneName=gm$ENSEMBL$SYMBOL[ens],
           Description=ens2biotype_ids[ens]) %>% 
    filter(Description=='protein_coding') %>% 
    select(-c(ens)) %>% 
    relocate(c('GeneName', 'Description')) %>% 
    filter(!duplicated(GeneName),
           !is.na(GeneName)) %>%
    left_join(x = ranks[,1,drop=F], 
              y=., by='GeneName') %>%
    rename_with(., .fn = ~c("Name", 'Description', make.unique(coldata_i$condition)))
  tpm_i[is.na(tpm_i)] <- 0
  
  ## Format the GSEA Enrichment Scores
  gsea_l <- unlist(ct_res_gsea_bkup[[celltype]]$interaction$gsea, recursive=F)
  gsea_df <- do.call(rbind, lapply(gsea_l[grep("REACTOME", names(gsea_l))], as.data.frame))
  gsea_df <- do.call(rbind, lapply(gsea_l, as.data.frame))
  gsea_dat_i <- gsea2CytoscapeFormat(gsea_df, gs_map)
  fileConn<-file(file.path(outdir, "gse", "cytoscape", 
                           paste0("meta.", celltype, ".cls")))
  # writeLines(paste(c(nrow(coldata_i), 
  #                    length(unique(coldata_i$condition)),
  #                    "1"), collapse="\t"), fileConn)
  # writeLines(paste0("#", paste(unique(coldata_i$condition), collapse="\t")), fileConn)
  writeLines(paste(coldata_i$condition, collapse="\t"), fileConn)
  close(fileConn)
  write.table(tpm_i, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("tpm.", celltype, ".expression.txt")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(ranks, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("l2fc.", celltype, ".rnk")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(gsea_dat_i$up, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("gsea.", celltype, ".up.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
  write.table(gsea_dat_i$down, 
              file=file.path(outdir, "gse", "cytoscape", 
                             paste0("gsea.", celltype, ".down.xls")),
              sep="\t", col.names = T, row.names = F, quote = F)
}



## Collapse GSEA results to a single table per celltype
grp_gsea_tbls_merge <- lapply(ct_res_gsea, function(grps){
  gsea_tbls <- lapply(grps, function(gsea_i) {
    lapply(gsea_i[which(!sapply(gsea_i, is.null))], function(i) {
      padj_cols <- grep("p.adjust", colnames(i), value=T)
      i %>% 
        top_n(wt=(!!as.symbol(padj_cols[1])) + (!!as.symbol(padj_cols[2])), 5)
    }) %>%
      do.call(rbind, .)
  })
  gsea_tbls_merge <- do.call(rbind, gsea_tbls) %>% 
    as.data.frame()
  return(gsea_tbls_merge)
})


## Visualize the GSEA results
gg_gseas <- lapply(grp_gsea_tbls_merge, function(gsea_tbls_merge){
  
  melt_gsea <- gsea_tbls_merge %>% 
    select(ID, grep("^NES", colnames(.), value=T), celltype, geneset) %>%
    reshape2::melt() %>%
    mutate(geneset=gsub("^.*?\\.", "", geneset) %>% 
             gsub("base", "Hallmark", .),
           ID=gsub("_", " ", ID) %>%
             gsub("HALLMARK |GOBP |REACTOME |GOCC |GOM F", "", .),
           variable=gsub("NES.cis", "NES (CIS)", variable) %>%
             gsub("NES.int", "NES (TDLN-LN)", .) %>%
             gsub("NES.tdln", "NES (TDLN_CIS-PBS)", .) %>%
             gsub("NES.pbs", "NES (PBS_TDLN-LN)", .)) %>%
    filter(geneset!='IMMUNESIGDB')
  gg_gsea <- ggplot(melt_gsea, aes(x=value, y=ID, fill=variable, group=variable)) +
    facet_grid(geneset ~ celltype, scales = "free_y", space='free',switch = "y") +
    geom_bar(stat='identity', position='dodge') +
    # scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
    scale_fill_manual(values=as.character(rev(delta_col))) +
    theme_classic() + ylab("") +
    geom_hline(aes(yintercept=ID), linetype='dashed', color='grey') + 
    theme(axis.text.y = element_text(size = 7,
                                     hjust=0), 
          strip.text.y.left = element_text(angle=90)) +
    scale_y_discrete(labels =  scales::wrap_format(26))
  gg_gsea
})

pdf(file.path("~/xfer", "gsea_response.pdf"), width = 10, height = 12)
gg_gseas
dev.off()




## Write out GSEA results
gsea_tbls <- lapply(ct_res_gsea, function(x) do.call(rbind, x))
gsea_tbls$MSM$celltype <- 'MSM'
gsea_tbls$SSM$celltype <- 'SSM'


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

##############################
#### 4) SCENIC - Regulons ####
## Initialize
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
coldata <- readRDS(file=file.path(outdir, "rds", "coldata.rds"))
dir.create(file.path(outdir, "regulon"), showWarnings = F)

# dir.create("SCENIC_msm_tdln", recursive = F, showWarnings = F)
setwd(file.path(pdir, "SCENIC_msm_tdln"))

dds <- res_dds$All$dds_lnBase
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org=scenic_org_code, 
                                  dbDir=dirname(rcis_db), dbs=basename(rcis_db), 
                                  datasetTitle='msm_ssm', nCores=10) 
scenicOptions@settings$defaultTsne <- list("perpl"=3,
                                           "dims"=30,
                                           "aucType"="AUC")

# Pre-filtering of the gene expression matrix
samples <- as.data.frame(colData(dds)) %>%
  dplyr::filter(celltype == 'MSM',
                lnstatus == 'TDLN') %>%
  rownames

exprMat <- assay(dds)[which(!is.na(ens2sym_ids[rownames(dds)])),samples]
rownames(exprMat) <- ens2sym_ids[rownames(exprMat)] 
exprMat <- exprMat[which(!duplicated(rownames(exprMat))),]
genesKept <- geneFiltering(exprMat, scenicOptions,
                           minCountsPerGene=10,
                           minSamples=3)
exprMat_filtered <- exprMat[genesKept, ]

# GENIE3 and SCENIC analysis wrapper
if(!file.exists(file.path("int", "3.4_regulonAUC.Rds"))){
  runCorrelation(tpmx_filtered, scenicOptions)
  exprMat_filtered_log <- log2(tpmx_filtered+1) 
  runGenie3(exprMat_filtered_log, scenicOptions)
  
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) #, coexMethod=c("top5perTarget")) #** Only for toy run!!
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, tpmx_filtered)
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
}


# Aug 15 test
# dir.create("SCENIC_msm_ln", recursive = F, showWarnings = F)
setwd(file.path(pdir, "SCENIC_msm_tdln"))
ilfamily <- list(
  IL1 = c("Il1a", "Il1b", "Il1rn"),
  IL2 = c("Il2", "Il4", "Il5", "Il9", "Il13", "Il15"),
  IL6 = c("Il6", "Il11", "Il27", "Il31", "Il33", "Il35", "Il36a", "Il36b", "Il36g"),
  IL10 = c("Il10", "Il19", "Il20", "Il22", "Il24", "Il26"),
  IL17 = c("Il17a", "Il17b", "Il17c", "Il17d", "Il17f"),
  Other = c("Il3", "Il7", "Cxcl8", "Il12a", "Il12b", "Il14", "Il16", "Il18", "Il34")
)  

data_groups <- c('msm_tdln'='SCENIC_msm_tdln', 'msm_ln'='SCENIC_msm_ln')
tf_grps <- lapply(data_groups, function(grp){
  setwd(file.path(pdir, grp))
  
  ## 0. Set up the SCENIC env
  data(list="motifAnnotations_mgi_v9", package="RcisTarget")
  motifAnnotations_mgi <- motifAnnotations_mgi_v9
  scenicOptions <- initializeScenic(org=scenic_org_code, 
                                    dbDir=dirname(rcis_db), dbs=basename(rcis_db), 
                                    datasetTitle='msm_ssm', nCores=10) 
  scenicOptions@settings$defaultTsne <- list("perpl"=3,
                                             "dims"=30,
                                             "aucType"="AUC")
  
  # Pre-filtering of the gene expression matrix
  samples <- as.data.frame(colData(dds)) %>%
    dplyr::filter(celltype == 'MSM',
                  lnstatus == 'TDLN') %>%
    rownames
  
  motifAnnot <- getDbAnnotations(scenicOptions) # Motif annotations
  allTFs=getDbTfs(scenicOptions) # all TFs
  rnkName = getDatabases(scenicOptions)
  il33_ranks <- importRankings(rnkName, columns='Il33')
  motifs_AUC <- loadInt(scenicOptions, 'motifs_AUC')
  
  ## 1) Load in the correlation matrix and extra cytokine links
  # corrMat <- loadInt(scenicOptions, 'corrMat')
  
  ## 2.a) Get the GENIE3 TF-modules that have IL33 within it:
  linkList <- loadInt(scenicOptions, "genie3ll")
  il33 <- linkList[grep("Il33", linkList$Target,ignore.case=T),] %>%  arrange(weight)
  ils_ll <- lapply(unlist(ilfamily), function(i){
    linkList[grep(paste0("^", i, "$"), linkList$Target,ignore.case=T),] %>% 
      dplyr::select(-Target) %>%
      magrittr::set_colnames(., c('TF', i))
  }) %>%
    purrr::reduce(full_join, by='TF') %>%
    tibble::column_to_rownames(., "TF")
  il33_vs_ilf <- apply(ils_ll, 1, function(i){
    ecdf(i)(i['Il33'])
  }) %>% sort
  
  ## 2.b) Load in the list of TF-modules and their genes that meet the cutoff for high correlation
  ## However, not all of these TFs have a matching motif that is found in Il33 promoter region
  tfModules_forEnrichment <- loadInt(scenicOptions, 'tfModules_forEnrichment')
  motifEnrichment_full <-  loadInt(scenicOptions, 'motifEnrichment_full')
  il33_all_ids <- names(which(sapply(tfModules_forEnrichment, function(i) any(grepl("Il33", i)))))
  il33_ids <- list("top"=grep(pattern="top.*Target", il33_all_ids, value = T),
                   "corr"=grep(pattern="w0", il33_all_ids, value = T))
  
  ## 3. For each of those modules, check if any of them are defined by a motif
  motifEnrichment<- loadInt(scenicOptions, 'motifEnrichment_full')
  TFs <- motifEnrichment %>% 
    dplyr::filter(geneSet %in% c(il33_ids$top)) %>% 
    dplyr::filter(TFinDB == '**') 
  TF_wMotifs <- unique(TFs$geneSet) # unique(gsub("_.*", "", TFs$geneSet))
  TF_wMotifs_genes <- unique(gsub("_.*", "", TFs$geneSet))
  
  ## 4. Cross check the significant TF_wMotifs and see if they contain Il33 or if they regulate
  ## TFs that do conatin Il33 (from section 2.b)
  onedeg_il33_ids <- sapply(TF_wMotifs_genes, function(gene_j){
    sapply(tfModules_forEnrichment[il33_ids$top], function(i) {
      gene_j %in% i
    })
  })
  onedeg_il33_tfs <- apply(onedeg_il33_ids, 2, function(i){
    setNames(gsub("_.*", "", names(which(i))),
             names(which(i))) %>% 
      sort %>% list
  }) %>% unlist(., recursive=F)
  onedeg_il33_genes <- lapply(onedeg_il33_tfs, function(i) {
    tfModules_forEnrichment[names(i)]
  })
  
  ## 5. Extract AUC rankings to be used for downstream calcAUC
  aucellRankings <- loadInt(scenicOptions, "aucell_rankings")
  
  return(list("TF"=onedeg_il33_tfs,
              "Genes"=onedeg_il33_genes,
              'Il33_vAll'=il33_vs_ilf,
              "aucellRankings"=aucellRankings))
})


## ResA] Find Il33-linked TFs that are matched between TDLN and LN,
## and compare Il33 weights to all other IL-family weights
tf_join <- full_join(tibble::rownames_to_column(as.data.frame(tf_grps$msm_tdln$Il33_vAll), "TF"),
          tibble::rownames_to_column(as.data.frame(tf_grps$msm_ln$Il33_vAll), "TF"),
          by='TF') %>%
  magrittr::set_colnames(., c('Tf', 'TDLN', 'LN')) %>%
  tibble::column_to_rownames(., 'Tf')
tf_join0 <- tf_join[which(rowSums(is.na(tf_join)) == 0),] %>%
  round(., 2) %>% 
  arrange(TDLN, desc(LN))


## ResB] Isolate for TFs that are supported by motifs and top targets for
## the TDLN and LN, within 1 degree of separation
## Then trim down the joined Il33-linked TFs based on this list
sig_tf_motifs <- rbind(reshape2::melt(tf_grps$msm_tdln$TF) %>%
                         mutate("group"="TDLN"),
                       reshape2::melt(tf_grps$msm_ln$TF) %>% 
                         mutate("group"="LN")) %>%
  magrittr::set_colnames(c("tf_corr", "tf_motif", "group"))
sig_tf_motifs_df <- tf_join[sig_tf_motifs[,1],] %>%
  cbind(sig_tf_motifs[,c(3:1)], .) %>%
  mutate("link"=ifelse(tf_motif == tf_corr, "primary", "1_deg"),
         "id"=ifelse(link=='primary', tf_motif, paste0(tf_motif, " > ", tf_corr))) %>%
  unique %>% 
  arrange(desc(group), tf_motif, desc(link), desc(TDLN))
.names2ypos <- function(namesdf){
  # namesdf <- sig_tf_motifs_df[,c('tf_motif', 'tf_corr')]
  namesdf[,1] <- factor(namesdf[,1], levels=unique(namesdf[,1]))
  ypos1 <- as.integer(namesdf[,1])
  namesl <- split(namesdf, namesdf[,1]) 
  ypos <- 0
  for(ni in seq_along(namesl)){
    i <- namesl[[ni]]
    ypos <- c(ypos,
              c(max(ypos)+2, (seq(nrow(i)-1) + max(ypos) + 1)))
  }
  return(cbind(namesdf, "ypos"=(ypos[-1] - 1)))
}
sig_tf_motifs_df <- cbind(sig_tf_motifs_df,
                          'ypos'=.names2ypos(sig_tf_motifs_df[,c('tf_motif', 'tf_corr')])$ypos)
sig_tf_motifs_df$xpos <- as.integer(factor(sig_tf_motifs_df$link,
                                           levels=c('primary', '1_deg')))
sig_tf_motifs_df$xpos <- sig_tf_motifs_df$xpos / (2*max(sig_tf_motifs_df$xpos))
sig_tf_motifs_df$ypos <- sig_tf_motifs_df$ypos / max(sig_tf_motifs_df$ypos)

sig_tf_motifs_df <- split(sig_tf_motifs_df, sig_tf_motifs_df$tf_motif) %>%
  lapply(., function(i) {
    i$xpos0 <- i$xpos[1]
    i$ypos0 <-  i$ypos[1]
    return(i)
    }) %>% do.call(rbind, .)

X <- melt(sig_tf_motifs_df, measure.vars=c('TDLN', 'LN'))  %>%
  mutate(id = factor(as.character(id), levels=unique(id)),
         blank=1-value)



library(scatterpie)
p <- ggplot(X, aes(x=xpos, y=ypos, label=tf_corr)) +
  geom_segment(aes(x=xpos0, y=ypos0, xend=xpos, yend=ypos)) +
  scatterpie::geom_scatterpie(data=X, aes(x=xpos, y=ypos), 
                              cols=c('value', 'blank'), pie_scale=2.5) +
  geom_text(data=dplyr::filter(X, link=='1_deg'),
            nudge_x=0.03, nudge_y=0, hjust=0, vjust=0) +
  geom_text(data=dplyr::filter(X, link=='primary'),
            nudge_x=-0.03, nudge_y=0, hjust=1, vjust=0) +
  coord_equal() +
  scale_fill_manual(values=c('value'='black', 'blank'='white')) +
  xlim(0, 1.0) +
  facet_grid(variable ~ ., switch = "y") +
  cowplot::theme_cowplot() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        legend.position = 'none',
        axis.line = element_blank())
pdf("~/xfer/y.pdf", height = 15)
p
p1
dev.off()

## ResC] Run AUC using all the genesets from both datasets (TDLN and LN)
## 5. Use the extended list of TFs from section 4 to check how Cis scores
## compare to PBS score
library(AUCell)
gene_list <- lapply(tf_grps, function(i) i$Genes) %>%
  unlist(., recursive=F)
for(id in names(tf_grps)){
  aucellRankings <- tf_grps[[id]]$aucellRankings
  set.seed(1234)
  regulonAUC <- AUCell_calcAUC(unlist(gene_list, recursive=F), 
                               aucellRankings, 
                               aucMaxRank=aucellRankings@nGenesDetected["1%"], 
                               nCores=1) 
  tf_grps[[id]]$regulon_auc <- assay(regulonAUC)
}

## ResD] Finally, For each of the tf_motifs in ResB, generate the boxplot
## comparisons between Cis and PBS treated samples
regulon_auc_df <- lapply(c("tdln", "ln") %>% setNames(.,.), function(grp){
  df <- tf_grps[[paste0("msm_", grp)]]$regulon_auc %>%
    as.data.frame
  cpidx <- grepl("Cis", colnames(df), ignore.case = T)
  
  apply(df, 1, function(row_i){
    X <- effsize::cohen.d(row_i[which(cpidx)], row_i[which(!cpidx)])
    c('estimate'=X$estimate, X$conf.int)
  }) %>% t %>% as.data.frame %>%
    mutate(group=toupper(grp),
           geneset_group=gsub("^(.*?)\\..*", "\\1", rownames(.)),
           # tf_motif=gsub("\\..*", "", rownames(.)),
           tf_corr=gsub("(.*)_.*", "\\1", rownames(.)) %>% 
             gsub("^.*\\.", "", .)) %>%
    tibble::remove_rownames() %>%
    dplyr::filter(!duplicated(tf_corr))
}) %>% do.call(rbind, .) %>%
  arrange(desc(geneset_group)) %>%
  mutate(tf_corr = factor(as.character(tf_corr), levels=unique(tf_corr)),
         group = factor(as.character(group), levels=c('TDLN', 'LN')))

p <- ggplot(regulon_auc_df, aes(y=tf_corr, x=estimate, fill=geneset_group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,
                position=position_dodge(.9)) +
  facet_grid(group ~ ., space='free', scales='free') +
  cowplot::theme_cowplot() +
  xlim(-140, 140) +
  geom_vline(xintercept = 0, linetype='dashed') +
  ylab("TF_module") + xlab("Cohen's D (AUC)") +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust=1, vjust=1))

pdf("~/xfer/x.pdf", height=12, width = 4)
p
dev.off()

###################
#### 5) ssGSEA ####
## Initialize
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
coldata <- readRDS(file=file.path(outdir, "rds", "coldata.rds"))
tpm <- read.table(file.path(outdir, "..", "counts", "all_tpm.tsv"), header=TRUE,
                  row.names="gene", check.names=FALSE, stringsAsFactors = F)
dir.create(file.path(outdir, "ssgsea"), showWarnings = F)
setwd(file.path(outdir, "ssgsea"))
# dir.create("SCENIC_msm_ln", recursive = F, showWarnings = F)


dds <- res_dds$All$dds_lnBase
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org=scenic_org_code, 
                                  dbDir=dirname(rcis_db), dbs=basename(rcis_db), 
                                  datasetTitle='msm_ssm', nCores=10) 
scenicOptions@settings$defaultTsne <- list("perpl"=3,
                                           "dims"=30,
                                           "aucType"="AUC")

# Pre-filtering of the gene expression matrix
samples_df <- as.data.frame(colData(dds)) %>%
  dplyr::filter(celltype == 'MSM') 
samples <- rownames(samples_df)
samplel <- split(samples, f=samples_df$condition)
  #               lnstatus == 'TDLN') %>%
  # rownames

exprMat <- assay(dds)[which(!is.na(gm$ENSEMBL$SYMBOL[rownames(dds)])),samples]
rownames(exprMat) <- gm$ENSEMBL$SYMBOL[rownames(exprMat)] 
exprMat <- exprMat[which(!duplicated(rownames(exprMat))),]
genesKept <- geneFiltering(exprMat, scenicOptions,
                           minCountsPerGene=10,
                           minSamples=3)
exprMat_filtered <- exprMat[genesKept, ]

tpmx <- tpm[which(!is.na(gm$ENSEMBL$SYMBOL[rownames(tpm)])),samples,]
tpmx <- tpmx[which(!duplicated(gm$ENSEMBL$SYMBOL[rownames(tpmx)])),]
rownames(tpmx) <- gm$ENSEMBL$SYMBOL[rownames(tpmx)] 
tpmx_filtered <- as.matrix(tpmx[genesKept,])


msig_lvls <- list('C2'=list('CP:REACTOME'))#              # curated gene sets
mlvl <- 'C2'
sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
  msig_ds <- msigdbr(species = 'Mus musculus', category = mlvl, subcategory = sublvl) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    as.data.frame()
  
  library(AUCell)
  cells_rankings <- AUCell::AUCell_buildRankings(exprMat_filtered, plotStats=FALSE)
  gene_list <- split(msig_ds$gene_symbol, msig_ds$gs_name)
  # gene_list <- gene_list[grep("(TLR9|STAT[0-9]|NFKB)", gsub("_", "", names(gene_list)))]
  set.seed(1234)
  genesetAUC <- AUCell::AUCell_calcAUC(gene_list, 
                                       cells_rankings, 
                                       # aucMaxRank=aucellRankings@nGenesDetected["1%"], 
                                       nCores=1) 
  ssgseaparam <- GSVA::ssgseaParam(tpmx_filtered, gene_list)
  ssgsea_dat_tpm<- GSVA::gsva(ssgseaparam)
  # X <- ssgsea_dat_tpm[order(apply(ssgsea_dat_tpm,1,median), decreasing = T),]
  
  df <- as.data.frame(assay(genesetAUC))
  Y <- df[order(apply(df, 1, median), decreasing = T),]
  top5perc <- ceiling(nrow(Y) * 0.05)
  
  pattern <- "(TLR|STAT[0-9]|NFKB|IRF|Jun|AP1)"
  ids <- setNames(grep(pattern, gsub("_", "", rownames(Y))),
                  rownames(Y)[grep(pattern, gsub("_", "", rownames(Y)))])
  ids <- ids[which(ids < top5perc)]
  
  .allbyall <- function(X, fun, ...){
    lapply(X, function(i) {
      lapply(X, function(j){
        fun(i,j, ...)
      }) %>% do.call(rbind, .) %>%
        as.data.frame %>%
        tibble::rownames_to_column(., "Alt") 
    }) %>%
      do.call(rbind, .) %>%
      mutate(Ref=as.character(sapply(names(X), rep, times=length(X)))) %>%
      dplyr::relocate(., "Ref") %>% 
      arrange(p)
  }
  .diff <- function(i,j,row_i){
    X <- effsize::cohen.d(row_i[i], row_i[j])
    p <- tryCatch({t.test(row_i[i], row_i[j])$p.value}, error=function(e){NA})
    c('estimate'=X$estimate, X$conf.int, 'p'=p)
  }
  
  # df <- as.data.frame(ssgsea_dat_tpm)
  geneset_auc_df <- apply(df, 1, function(row_i){
    .allbyall(samplel, .diff, row_i=row_i)
  }) 
  gsselect_auc_df <- do.call(rbind, geneset_auc_df[names(ids)]) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., 'geneset') %>%
    dplyr::mutate(geneset=gsub("\\..*", "", geneset)) %>%
    dplyr::filter(Ref == 'TDLN-CIS_MSM' & 
                    (Alt == 'LN-CIS_MSM' | Alt == 'TDLN-PBS_MSM'))
  
  pdf("~/xfer/x.pdf", width = 10)
  gsselect_auc_df %>% 
    dplyr::mutate(geneset = gsub("REACTOME_", "", geneset) %>%
                    gsub("_", " ", .)) %>%
    ggplot(aes(y= geneset, x = estimate, fill = Alt)) + 
    geom_bar(stat="identity", alpha=0.5, 
             position=position_dodge()) +
    geom_errorbar(aes(xmin=lower, xmax=upper), width=.2, colour="orange", 
                  position=position_dodge(.9)) +
    xlab("Cohen's D") + ylab("") + 
    cowplot::theme_cowplot() +
    ggtitle("TDLN-CIS_MSM vs ...")
  dev.off()
  # %>% t %>% as.data.frame %>%
  #   mutate(geneset_group=gsub("^REACTOME_", "", rownames(.))) %>%
  #   tibble::remove_rownames() %>%
  #   arrange(p)
  return(geneset_auc_df)
})


###########################
#### 6) Final Figures: ####
#--- June 3 2024 - Fig1.F: Heatmap for MSM and SSM ----
order_type <- 'alphabetical'
merge_replicates <- FALSE

dir.create(file.path(outdir, "final_figures"), showWarnings = F)
goi <- c('Mertk', 'Axl', 'Cd36', 'Il10', 'Tgfb3', 'Ido1', 'Ido2',
         'Cyp1a1', 'Cyp1b1', 'Ddit3', 'Pdcd1lg2', 'Timd4', 'Irf4',
         'Lepr', 'Tslp', 'Mrc1', 'Ccl17', 'Cxcl2', 'Ccl4', 'Il6',
         'Il1b', 'Gzmb', 'Vegfc', 'Il17a', 'Il1rn', 'Cd80', 'Cd40lg')
celltype_order <- c('SSM', 'MSM')
lnstatus_order <- c('LN', 'TDLN')
treatment_order <- c('PBS', 'CIS')

celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
cts <- read.table(file.path("counts", "all_tpm.tsv"), header=TRUE,
                  row.names="gene", check.names=FALSE, stringsAsFactors = F)

annotation_colors=list(
  LNstatus = c('TDLN'='#7fbc41', 'LN'='#de77ae'), # pink, green
  Treatment = c('Cisplatin'='#006d2c', 'Untreated'='#a1d99b'), # light-green, dark-green
  Celltype=c('SSM'='white', 'MSM'='black')  
)
res_celltype <- res_dds[[c('All')]]

cat(paste0("Params: \n\t order_type = ", order_type, 
           "\n\t merge_replicates = ", merge_replicates, "\n"))
# row_ord <- NULL # Set to NULL if you want to cluster the genes
row_ord <- (goi) # Specify an order if you want to keep that gene-order
heatmaps <- lapply(res_dds[c('All')], function(res_celltype){
  dds <- res_celltype$dds_lnBase 
  celltype_id <- unique(colData(dds)$celltype) %>%
    paste(., collapse="_")
  
  cnt_mat <- assay(dds)
  cnt_mat_goi <- cnt_mat %>%
    as.data.frame %>% 
    tibble::rownames_to_column(., "id") %>%
    mutate(id=ens2sym_ids[id]) %>%
    # filter(id %in% goi) %>% 
    dplyr::filter(!duplicated(id),
           !is.na(id)) %>%
    tibble::column_to_rownames(., "id")
  cat(paste0("Genes not in your GOI: ", 
             paste(goi[!goi %in% rownames(cnt_mat_goi)], collapse="\n")), "\n")
  
  vst_assay <- vst(cnt_mat, blind=T)
  vst_assay_goi <- vst_assay %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "id") %>%
    mutate(id=ens2sym_ids[id]) %>%
    dplyr::filter(id %in% goi) %>%
    dplyr::filter(!duplicated(id),
                  !is.na(id)) %>%
    tibble::column_to_rownames(., "id")
  vst_assay_mat <- as.matrix(cnt_mat_goi)
  vst_assay_mat <- as.matrix(vst_assay_goi)
  storage.mode(vst_assay_mat) <- 'numeric'
  # cor_vst <- cor(t(vst_assay_mat), method='spearman')
  # cor_vst['Il33',] %>% round(., 3) %>% sort %>% tail(., 50) %>% as.matrix
  
  vst_grp_goi_list <- split(as.data.frame(t(vst_assay_goi)),#[goi,])), 
                       colData(dds)$condition)
  if(merge_replicates){
    # Merge biological replicates into one group average
    vst_grp_goi <- vst_grp_goi_list %>%
      lapply(., colMeans) %>%
      do.call(rbind, .) %>% t
  } else {
    # Keep biological replicates separate and annotate accordingly
    vst_grp_goi <- vst_grp_goi_list %>% 
      do.call(rbind, .) %>% 
      t %>% as.data.frame %>%
      rename_with(., ~make.unique(gsub("\\..*", "", .)))
  }
  
  # setdiff(goi, rownames(vst_grp_goi))  # Check if any goi are not found in the matrix
  annotation_col <- colnames(vst_grp_goi) %>% 
    gsub("\\..*", "", .) %>% 
    strsplit(., split="-|_") %>% 
    do.call(rbind, .) %>% as.data.frame  %>%
    rename_with(., ~c('LNstatus', 'Treatment', 'Celltype')) %>%
    mutate(Celltype = factor(Celltype, levels=celltype_order),
           LNstatus = factor(LNstatus, levels=lnstatus_order),
           Treatment = factor(Treatment, levels=treatment_order))
  treatment_anno_order <- c('Untreated', 'Cisplatin')
  if(treatment_order[1] == 'CIS') treatment_anno_order <-  rev(treatment_anno_order)
  levels(annotation_col$Treatment) <- treatment_anno_order
  order_idx <- with(annotation_col, order(Celltype, LNstatus, Treatment))

  
  
  print(row_ord)
  if(is.null(row_ord)) {
    row_ord <<- switch(order_type,
                       cis_pbs_diff=(order(vst_grp_goi[,'TDLN-CIS_MSM'] - vst_grp_goi[,'TDLN-PBS_MSM'])),
                       alphabetical=(order(rownames(vst_grp_goi))))
  }
  
  return(list("matrix"=vst_grp_goi[row_ord, order_idx], "annotation_col"=annotation_col[order_idx,]))
}) 

# Annotation Meta
annotation_col <- do.call(rbind, lapply(heatmaps, function(i) i$annotation_col))
annotation_col$split <- with(annotation_col, paste0(Celltype, LNstatus))
annotation_col$split <- factor(annotation_col$split, levels=unique(annotation_col$split))

# Normalize for a condition (i.e. Cis-PBS)
lnidx <- grep("^LN", annotation_col$LNstatus)
tdlnidx <- grep("^TDLN", annotation_col$LNstatus)
raw_vst <- heatmaps$All$matrix
colnames(raw_vst) <- with(annotation_col, paste(LNstatus, Celltype, Treatment, sep="_"))
scaled_vst <- t(apply(heatmaps$All$matrix, 1, scale))
colnames(scaled_vst) <- with(annotation_col, paste(LNstatus, Celltype, Treatment, sep="_"))
scaled_grp_goi <- heatmaps$All$matrix[,tdlnidx] - heatmaps$All$matrix[,lnidx]
scaled_grp_goi <- t(apply(scaled_grp_goi, 1, scale))
# scaled_grp_goi <- scaled_grp_goi[which(goi %in% rownames(scaled_grp_goi)),]
colnames(scaled_grp_goi) <- with(annotation_col[lnidx,], 
                                 paste(Celltype, Treatment, sep="_")) %>% 
  make.unique

###### Important segment: Outputs the res datastructure for Sara to make heatmaps
# res <- list("raw_vst"=raw_vst,
#             "scaled_vst"=scaled_vst,
#             "cnt_mat"=scaled_grp_goi,
#             "ens2sym_ids"=ens2sym_ids,
#             "raw_condition"=annotation_col,
#             "condition"=annotation_col[lnidx,])
# saveRDS(res, file=file.path("~/xfer/res.rds"))

gg_phm <- ComplexHeatmap::pheatmap(as.matrix(scaled_grp_goi), scale = 'none',
                         color = colorRampPalette(c("#053061", "#4393c3", "#f7f7f7", 
                                                    "#d6604d", "#67001f"))(20),
                         annotation_col = annotation_col[lnidx,],
                         column_split= annotation_col[lnidx,]$split,
                         annotation_colors = annotation_colors,
                         fontsize_row = 14,show_colnames = FALSE,
                         use_raster=FALSE, cluster_cols = FALSE, cluster_rows = FALSE
) %>% 
  ggplotify::as.ggplot(.) #+
  #theme_minimal(base_family = "Arial")
  
# pdf(file.path(outdir, "final_figures", "msm_ssm_heatmap_gois.pdf"), width =9)
pdf(file.path("~/xfer", "msm_ssm_heatmap_gois.y.pdf"), height = 12, width =4.5)
gg_phm
dev.off()

#--- Fig1.b: GSEA for MSM and SSM ----
plot_classical_gsea <- FALSE

## Read in gene sets
dirto <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_MSM_SSM/ref/genesets'
bus.inf <- read.table(file.path(dirto, "buscher_inflamm.txt"), header = F) %>%
  mutate(V1=str_to_title(V1))
rahul_dn <- read.table(file.path(dirto, "rahul_up.csv"), header = F) %>%
  mutate(V1=str_to_title(V1)) # Rahul's test conditions are reversed from ours
rahul_up <- read.table(file.path(dirto, "rahul_down.csv"), header = F) %>%
  mutate(V1=str_to_title(V1)) # Rahul's test conditions are reversed from ours

genesets <- list("Buscher.Inflammatory"=bus.inf$V1, 
                    "Rahul_down"=rahul_dn$V1, "Rahul_up"=rahul_up$V1)

geneset_df <- lapply(names(genesets), function(id){
  data.frame("geneset"=id, "genes"=genesets[[id]])
}) %>% 
  do.call(rbind, .)
msig_lvls <- list('custom'=genesets)

res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
gsea_res <- lapply(c('SSM', 'MSM'), function(celltype){
  res <- res_dds[[celltype]]$res %>%
    mutate(biotype=ens2biotype_ids[ens]) %>%
    filter(biotype=='protein_coding')
  gene_id <- res$gene
  # abslfc <- abs(res$log2FoldChange.interaction) * (-1*log10(res$padj.interaction + 0.0001))
  lfc <- res$log2FoldChange.interaction#  * (-1*log10(res$padj.interaction + 0.0001))
  lfc_genes <- setNames(lfc, gene_id)
  
  gseaFun <- function(msig_ds, lfc_v){
    gsea <- tryCatch({
      GSEA(sort(na.omit(lfc_v), decreasing = T), 
           TERM2GENE = msig_ds, pvalueCutoff = 1, maxGSSize=1000)
    }, error=function(e){NULL})
    return(gsea)
  }
  oras <- iterateMsigdb(species='Mus musculus', fun=gseaFun, 
                        lfc_v=lfc_genes, msig_lvls=msig_lvls)
  
  .mkgseval <- function(x){
    do.call(rbind, x) %>% 
      mutate(celltype=celltype,
             group=ifelse(grepl("Inflammatory", Description), "Inflammatory", "Tolerogenic"))
  }
  gseval <- lapply(oras$custom, function(i) i@result) %>% .mkgseval
  
  # merge up and down tolerogenic signatures
  {
    gseval_l <- split(gseval, f=gseval$group)
    tolero <- gseval_l$Tolerogenic
    upidx <- sapply(c('up', 'down'), grep, x=tolero$ID)
    tolero$enrichmentScore <- as.numeric(tolero$enrichmentScore)
    tolero$NES <- as.numeric(tolero$NES)
    tolero_es <- rbind(tolero[upidx['up'], c('enrichmentScore', 'NES')],
                       (-1*tolero[upidx['down'], c('enrichmentScore', 'NES')])) %>%
      colMeans
    pvals <- c(with(tolero, if(NES[upidx['up']] > 0) pvalue[upidx['up']] else 0.99), 
               with(tolero, if(NES[upidx['down']] < 0) pvalue[upidx['down']] else 0.99))
    tolero_p <- metap::sumz(p=as.numeric(pvals))$p %>%
      as.numeric
                
    
    tolero[3,] <- with(tolero,
                       c(rep('Rahul.Tolerogenic', 2), 
                         sum(as.numeric(setSize)), tolero_es, 
                         rep(tolero_p, 2),
                         NA, NA, NA, paste(core_enrichment, collapse="/"),
                         unique(celltype), unique(group)))
    gseval <- rbind(gseval_l$Inflammatory, tolero[3,,drop=F])
  }                
  
  if(plot_classical_gsea){
    pdf(file.path("~/xfer/", paste0("sara_gsea_plots.", celltype, ".pdf"))); 
    for(tit in names(oras$custom)){
      x <- gseaplot2(oras$custom[[tit]], geneSetID = 1, 
                     title=tit, pvalue_table=F)
      print(x)
    }
    gseaplot2(oras_tolerogenic$custom[['Rahul.Tolerogenic']], geneSetID = 1, 
              title='Tolerogenic', pvalue_table=F)
    gseaplot2(oras_inflammatory$custom[['Buscher.Inflammatory']], geneSetID = 1, 
              title='Inflammatory', pvalue_table=F)
    dev.off()
  }
  
  return(gseval)
})

gsea_resm <- do.call(rbind, gsea_res) %>%
  mutate(qvalues=round(p.adjust(pvalue, method='BH'),2),
         celltype=factor(celltype, levels=c('SSM', 'MSM')),
         vjust=ifelse(NES>0, 1.5,0),
         NES=as.numeric(NES)) %>%
  mutate(sig=ifelse((qvalues > 0.05), "", ifelse(qvalues > 0.01, " * ", "***")))

pdf("~/xfer/sara_gsea_nes.pdf", height = 4)
ggplot(gsea_resm, aes(x=NES, y=ID, color=celltype, fill=celltype, group=celltype, label=sig)) +
  facet_grid(group~., scales = 'free', space = 'free', switch='y') +
  geom_bar(stat='identity', position='dodge', alpha=0.8) +
  geom_text(aes(vjust=vjust), position = position_dodge(width = 1), 
            hjust=0.5, angle = 90, size = 10 / .pt) +
  xlim(-2,2) + ylab("") +
  theme_classic() +
  scale_fill_manual(name="",  values=c('MSM'='gray25', 'SSM'='gray75')) +
  scale_color_manual(name="",  values=c('MSM'='gray25', 'SSM'='gray75')) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()
write.table(gsea_resm, file=file.path("~/xfer", "sara_gsea_nes.csv"),
            quote = F, row.names = F, col.names = T, sep=",")


#--- Fig2.b: Waterfall of select genes ----
# Goal: Waterfall plot, showing the Log2FC MSM vs SSM for the LN-PBS for select genes
dir.create(file.path(outdir, "final_figures"), showWarnings = F)
goi <- c('Usp4', 'Sigirr', 'Pdlim1', 'Hes1', 'Arr1', 'Socs3', 
         'Arrb2', 'Zfp36', 'Klf4', 'Dusp1')
goi <- factor(goi, levels=goi)

celltypes <- setNames(celltypes,celltypes)
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
lfc_goi <- res_dds$Celltype$res %>%
  dplyr::filter(symbol %in% goi) %>%
  dplyr::select(symbol, `log2FoldChange.celltype_LN-PBS`, `padj.celltype_LN-PBS`) %>%
  rename_with(., ~c('symbol', 'logFC', 'padj')) %>%
  mutate(symbol=factor(symbol, levels=goi),
         sig=(padj < 0.05), 
         dir=(logFC>0))
b <- c(0, 5, 10, 20)

pdf(file.path(outdir, "final_figures", "waterfall_celltypes_ln_pbs.pdf"), height = 3.4)
pdf(file.path("~/xfer", "waterfall_celltypes_ln_pbs.pdf"), height = 3.4)
ggplot(lfc_goi, aes(x=symbol, y=logFC, fill=logFC)) +
  geom_bar(stat='identity', position='stack') +
  geom_text(aes(label = ifelse(sig, "*", "")), 
            position = position_dodge(width = .9), 
            vjust = ifelse(lfc_goi$dir, 0.4, 1.2), size = 20 / .pt) +
  scale_fill_gradientn(limits = c(-2,2),
                       colours=c("#08519c", "grey", "#a50f15"),
                       breaks=c(-2,0,2)) +
  xlab("") +
  theme_bw()
dev.off()
  


#--- SFig2.b: Volcano Plot of SSM and MSM - DEGs for TDLN-CIS vs TDLN-PBS ----
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))

lfc_lim <- 2
qval_lim <- 0.01
maxq <- 100
cols <- c("sigup"="#fc8d59", "sigdown"="#91bfdb")

gg_volcanos <- lapply(c('MSM', 'SSM'), function(id){
  res_dds_i <- res_dds[[id]]
  res_dds_i$res %>%
    dplyr::select(grep("padj", colnames(.), value=T)) %>%
    apply(., 2, summary) %>% round(., 4)
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
  
pdf(file.path(outdir, "degs", "msm_ssm.volcano.pdf"), height = 4)
# pdf("~/xfer/msm_ssm.volcano.pdf", height = 4)
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

lfc_goi <- res_dds$Celltype$res %>%
  filter(symbol %in% goi) %>%
  dplyr::select(symbol, `log2FoldChange.celltype_LN-PBS`, `padj.celltype_LN-PBS`) %>%
  rename_with(., ~c('symbol', 'logFC', 'padj')) %>%
  mutate(symbol=factor(symbol, levels=goi),
         sig=(padj < 0.05), 
         dir=(logFC>0))
b <- c(0, 5, 10, 20)