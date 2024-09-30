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
  cts <- read.table(file.path(cntsdir, "all.tsv"), header=TRUE, 
                    row.names="gene", check.names=FALSE)
  tpm <- read.table(file.path(cntsdir, "all_tpm.tsv"), header=TRUE, 
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

#####################################################################
#### 2) Differential expression controlling for interaction terms ####
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
#### 2.b) Overlap of sigGenes between MSM and SSM ####
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

# inters <- make.unique(rep("intersect", 1508))
# venn_x <- list("MSM"=c(make.unique(rep("MSM", 2555)), inters),
#                "SSM"=c(make.unique(rep("SSM", 2937)), inters))
# library(ggvenn)
# pdf("~/Projects/mcgaha/st2_il33/results/mouse/sara_MSM-SSM/results/degs/overlap_MSM_SSM/ggvenn.pdf",
#     height = 3)
# ggvenn(
#   venn_x, fill_color = c('#7b3294', '#008837'),
#   stroke_size = 0.5, set_name_size = 4
# )
# dev.off()

#################################
#### 3) Gene set of interest ####
resl <- readRDS(file.path(outdir, "rds", "resl.rds"))
ref_dir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/mouse_MSM_SSM/ref'

st4_inflammatory_goi <- read.csv(file.path(pdir, "..", "ref", "st4_inflammatory.csv"),
                                 sep=",", header = F) %>%
  unlist() %>% as.character() %>% str_to_title()
st4_antiinflammatory_goi <- read.csv(file.path(pdir, "..", "ref", "st4_antiinflammatory.csv"),
                                 sep=",", header = F) %>%
  unlist() %>% as.character() %>% str_to_title()

list2_goi <- read.csv(file.path(pdir, "..", "ref", "list2.csv"),
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
#   })  # for GSEA2cytoscape
#   return(ct_res_gsea)  # for GSEA2cytoscape
# })  # for GSEA2cytoscape
    
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
# saveRDS(ct_res_gsea_bkup, file=file.path(outdir, "GSEA-response2.rds"))
# saveRDS(ct_res_gsea, file=file.path(outdir, "GSEA-response.rds"))
saveRDS(ct_res_gsea, file=file.path(outdir, "GSEA-response.selectedGenes.2cyto.rds"))
saveRDS(ct_res_gsea, file=file.path(outdir, "GSEA-response.allGenes.2cyto.rds"))
# ct_res_gsea <- readRDS(file.path(outdir, "GSEA-response.rds"))
ct_res_gsea <- readRDS(file.path(outdir, "GSEA-response.selectedGenes.2cyto.rds"))
# ct_res_gsea_bkup <- readRDS(file.path(outdir, "GSEA-response2.rds"))
ct_res_gsea_bkup <- readRDS(file.path(outdir, "GSEA-response.selectedGenes.2cyto.rds"))
ct_res_gsea_bkup <- readRDS(file.path(outdir, "GSEA-response.allGenes.2cyto.rds"))
ct_res_gsea_bkup <- ct_res_gsea_bkup$grp1

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


# pdf(file.path(outdir, "gsea_response2.pdf"), width = 10, height = 12)
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
#### 6) SCENIC - Regulons ####
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
{
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
  
  
  
  
}
# Read in Regulon information
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulon_genesets <- readRDS("int/2.6_regulons_asGeneSet.Rds")
regulon_genesets_df <- lapply(names(regulon_genesets), function(i){
  data.frame("Regulon"=i, "Genes"=regulon_genesets[[i]])
}) %>% do.call(rbind, .)
write.table(regulon_genesets_df, file=file.path("~/xfer/regulons_df.csv"),
            sep=",", row.names = F, col.names = T, quote = F)

# Set the sample-order plotting
sample_ord <- coldata[with(coldata, order(celltype, lnstatus, treatment)),]

# Generate the Jenson-Shannon Divergence value (0 = identical, 1 = divergent) between a cell-type
# and the entire population, where RSS = 1- sqrt(JSD)
cols <- c('#969696', '#ef3b2c', '#67000d')
rssPlots <- lapply(colnames(coldata), function(col_i){
  rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=coldata[colnames(regulonAUC),col_i])
  list("plot"=plotRSS(rss[,unique(sample_ord[,col_i])],
                      col.low=cols[1], col.mid=cols[2], col.high=cols[3]), 
       "rss"=reshape2::melt(rss))
})
names(rssPlots) <- colnames(coldata)

# Identify regulons that are descriptive of the celltype,lnstatus,group and plot the 
# normalized AUC values +/- sd for each group
auc_l <- calcMeanAUC(getAUC(regulonAUC), 
                     cellAnnotation=coldata[colnames(regulonAUC),]$condition)
thresh_auc <- auc_l$auc
# [LN|Tumor] _ [Cis|PBS] _ [MSM|SSM]
#fg1: = paste(V3_V1): paste(MSM|SSM_LN|Tumor)
#fg1: = paste(V2_V3): paste(CIS|PBS_MSM|SSM)
site_auc_grp <- thresh_auc %>% .splitFg(., 'fg3') # Constant CIS_MSM: Compare between TDLN/LN
trt_auc_grp <- thresh_auc %>% .splitFg(., 'fg1') # Constant MSM_LN: compare between PBS/Cis

#--- A1) Calculate delta-avgAUC (difference) between TDLN and LN for CIS_[MS]SM:
# Returns: A list of MSM and SSM (both CIS)
#          each dataframe (e.g. MSM) contains the LN and TDLN AUC, and the delta (TDLN-LN)
cis_ids <- c("MSM"='CIS_MSM', "SSM"='CIS_SSM')
cis_sm_delta <- lapply(cis_ids, function(grp_id){
  # Create a regulon by LN,Tumor matrix; calculate Delta
  reshape2::dcast(site_auc_grp[[grp_id]],
                                Regulon ~ V1, value.var='mean') %>%
    mutate(delta=TDLN-LN)  %>% 
    rename_with(., ~ paste0(., ".auc"), matches("[^Regulon]"))
})

#--- A2) Calculate delta-avgAUC (difference) between TDLN and LN  for PBS_[MS]SM:
# Returns: A list of MSM and SSM (both PBS)
#          each dataframe (e.g. MSM) containts the LN and TDLN AUC, and the delta (TDLN-LN)
pbs_ids <- c("MSM"='PBS_MSM', "SSM"='PBS_SSM')
pbs_sm_delta <- lapply(pbs_ids, function(grp_id){
  reshape2::dcast(site_auc_grp[[grp_id]],
                  Regulon ~ V1, value.var='mean') %>%
    mutate(delta=TDLN-LN)  %>% 
    rename_with(., ~ paste0(., ".auc"), matches("[^Regulon]"))
})

#--- B1) Calculate delta-avgAUC (difference) between Cis and PBS for [MS]SM_TDLN:
# Returns: A list of MSM and SSM (both TDLN)
#          each dataframe (e.g. MSM) contains the Cis and PBS AUC, and the delta (Cis-PBS)
tdln_ids <- c("MSM"='MSM_TDLN', "SSM"='SSM_TDLN')
tdln_sm_delta <- lapply(tdln_ids, function(grp_id){
  reshape2::dcast(trt_auc_grp[[grp_id]],
                  Regulon ~ V2, value.var='mean') %>%
    mutate(delta=CIS-PBS)  %>% 
    rename_with(., ~ paste0(., ".auc"), matches("[^Regulon]"))
})

#--- B2) Calculate delta-avgAUC (difference) between Cis and PBS for [MS]SM_LN:
# Returns: A list of MSM and SSM (both LN)
#          each dataframe (e.g. MSM) contains the Cis and PBS AUC, and the delta (Cis-PBS)
ln_ids <- c("MSM"='MSM_LN', "SSM"='SSM_LN')
ln_sm_delta <- lapply(ln_ids, function(grp_id){
  reshape2::dcast(trt_auc_grp[[grp_id]],
                  Regulon ~ V2, value.var='mean') %>%
    mutate(delta=CIS-PBS)  %>% 
    rename_with(., ~ paste0(., ".auc"), matches("[^Regulon]"))
})


#--- C) Calculate delta-avgAUC (difference) between Cis and PBS for MSM_[TD]LN,
# With that, calculate the distance between MSM_TDLN and MSM_LN corrected for Cis:
# Returns: A list of MSM and SSM (contains both TLDN and LN)
#   each dataframe contains:
#     [a] TDLN: Cis, PBS, and delta(Cis-PBS)
#     [b]   LN: Cis, PBS, and delta(Cis-PBS)
#        Delta: [a] - [b]; TDLN_Delta(Cis-PBS) - LN_Delta(Cis-PBS) 
int_ids <- list("MSM"=c('MSM_TDLN', 'MSM_LN'),
                "SSM"=c('SSM_TDLN', 'SSM_LN'))
int_sm_delta <- lapply(int_ids, function(grp_id){
  auc_tdln <- reshape2::dcast(trt_auc_grp[[grp_id[1]]],
                         Regulon ~ V2, value.var='mean') %>%
    mutate(delta=CIS-PBS) %>% 
    rename_with(., ~ paste0(., ".tdln"), matches("[^Regulon]"))
  auc_ln <- reshape2::dcast(trt_auc_grp[[grp_id[2]]],
                              Regulon ~ V2, value.var='mean') %>%
    mutate(delta=CIS-PBS) %>% 
    rename_with(., ~ paste0(., ".ln"), matches("[^Regulon]"))
  
  full_join(auc_tdln, auc_ln, by='Regulon') %>%
    mutate("delta.tdln_ln"=delta.tdln - delta.ln)
})

regulon_aucs <- list("interaction"=int_sm_delta,   # C [B1 - B2]
                     "tdln_cis_pbs"=tdln_sm_delta, # B1 TDLN Cis-x-PBS
                     "ln_cis_pbs"=ln_sm_delta,     # B2 LN Cis-x-PBS
                     "pbs_tdln_ln"=pbs_sm_delta,   # A2 PBS TDLN-x-LN
                     "cis_tdln_ln"=cis_sm_delta)   # A1 Cis TDLN-x-LN
saveRDS(regulon_aucs, file=file.path(outdir, "regulon", "regulon_auc.rds"))

regulon_aucs <- readRDS(file.path(outdir, "regulon", "regulon_auc.rds"))
# [Cis TDLN-x-LN] + [Interaction] regulon AUC plots   [A1, C]
{
grp1_regulons <- .getRegulons(intx=reshape2::melt(regulon_aucs$interaction$MSM),
                              cisx=reshape2::melt(regulon_aucs$cis_tdln_ln$MSM)) %>%
  tail(., n=10)
gg_grp1_auc_plot <- lapply(setNames(celltypes,celltypes), function(ct){
  aucPlot(int_i=regulon_aucs$interaction[[ct]], 
          grp1_i=regulon_aucs$cis_tdln_ln[[ct]], grp2_i=regulon_aucs$interaction[[ct]], 
          new_ids=c('(Cis_TDLN-LN)', '(Interaction)'), 
          celltype=ct, regulons = grp1_regulons,
          delta_grp1="delta.auc", delta_grp2='delta.tdln_ln',
          delta_col=c('paleturquoise3','gray30'))
})
}

# [PBS TDLN-x-LN] + [TDLN Cis-x-PBS] regulon AUC plots  [A2, B1]
{
grp2_regulons <- .getRegulons(cisx=reshape2::melt(regulon_aucs$pbs_tdln_ln$MSM),
                              intx=reshape2::melt(regulon_aucs$tdln_cis_pbs$MSM),
                              delta_int='delta.auc') %>%
  tail(., n=10)
gg_grp2_auc_plot <- lapply(setNames(celltypes,celltypes), function(ct){
  aucPlot(int_i=regulon_aucs$interaction[[ct]],
          grp1_i=regulon_aucs$pbs_tdln_ln[[ct]], grp2_i=regulon_aucs$tdln_cis_pbs[[ct]], 
          new_ids=c('(PBS TDLN-LN)', '(TDLN Cis-PBS)'), 
          celltype=ct, regulons = grp2_regulons,
          delta_grp1="delta.auc", delta_grp2='delta.auc',
          delta_col=c('royalblue2', 'indianred1'))
})
}

# [TDLN Cis-x-PBS] + [LN Cis-x-PBS] regulon AUC plots  [B1, B2]
{
grp3_regulons <- .getRegulons(cisx=reshape2::melt(regulon_aucs$tdln_cis_pbs$MSM),
                              intx=reshape2::melt(regulon_aucs$ln_cis_pbs$MSM),
                              delta_int='delta.auc') %>%
  tail(., n=10)
gg_grp3_auc_plot <- lapply(setNames(celltypes,celltypes), function(ct){
  aucPlot(int_i=regulon_aucs$interaction[[ct]],
          grp1_i=regulon_aucs$tdln_cis_pbs[[ct]], grp2_i=regulon_aucs$ln_cis_pbs[[ct]], 
          new_ids=c('(TDLN Cis-PBS)', '(LN Cis-PBS)'), 
          celltype=ct, regulons = grp3_regulons,
          delta_grp1="delta.auc", delta_grp2='delta.auc',
          delta_col=c('indianred1', 'orangered2'))
})
}

# [Cis TDLN-x-LN] + [PBS TDLN-x-LN] regulon AUC plots  [A1, A2]
{
grp4_regulons <- .getRegulons(cisx=reshape2::melt(regulon_aucs$cis_tdln_ln$MSM),
                              intx=reshape2::melt(regulon_aucs$pbs_tdln_ln$MSM),
                              delta_int='delta.auc') %>%
  tail(., n=10)
gg_grp4_auc_plot <- lapply(setNames(celltypes,celltypes), function(ct){
  aucPlot(int_i=regulon_aucs$interaction[[ct]],
          grp1_i=regulon_aucs$cis_tdln_ln[[ct]], grp2_i=regulon_aucs$pbs_tdln_ln[[ct]], 
          new_ids=c('(Cis TDLN-LN)', '(PBS TDLN-LN)'), 
          celltype=ct, regulons = grp4_regulons,
          delta_grp1="delta.auc", delta_grp2='delta.auc',
          delta_col=c('paleturquoise3','royalblue2'))
})
}


# Do the plotties
## Regulon-AUC values and their corresponding delta's between CIS and Interaction
pdf(file.path(outdir, "regulon", "scenic_regulons.pdf"), width = 13, height = 6)
# pdf(file.path("~/xfer", "scenic_regulons2.pdf"), width = 13, height = 6)
gg_grp1_aucs <- plot_grid(plotlist=gg_grp1_auc_plot, nrow=1, align='h', axis='bt')
gg_grp2_aucs <- plot_grid(plotlist=gg_grp2_auc_plot, nrow=1, align='h', axis='bt')
gg_grp3_aucs <- plot_grid(plotlist=gg_grp3_auc_plot, nrow=1, align='h', axis='bt')
gg_grp4_aucs <- plot_grid(plotlist=gg_grp4_auc_plot, nrow=1, align='h', axis='bt')
plot_grid(plot_grid(demoAucPlot('A1', 'C'), ncol=2), 
          gg_grp1_aucs, nrow=2, align='v', axis='lr', rel_heights = c(1,3.5))
plot_grid(plot_grid(demoAucPlot('A2', 'B1'), ncol=2), 
          gg_grp2_aucs, nrow=2, align='v', axis='lr', rel_heights = c(1,3.5))
plot_grid(plot_grid(demoAucPlot('B1', 'B2'), ncol=2), 
          gg_grp3_aucs, nrow=2, align='v', axis='lr', rel_heights = c(1,3.5))
plot_grid(plot_grid(demoAucPlot('A1', 'A2'), ncol=2), 
          gg_grp4_aucs, nrow=2, align='v', axis='lr', rel_heights = c(1,3.5))
dev.off()

write.table(int_sm_delta$MSM, 
            file=file.path("~/xfer", "MSM.regulon_delta.csv"),
            sep=",", quote=F, row.names = F, col.names = T)
write.table(int_sm_delta$SSM, 
            file=file.path("~/xfer", "SSM.regulon_delta.csv"),
            sep=",", quote=F, row.names = F, col.names = T)
## Regulon-RSS values marking their group-specificity
pdf(file.path(outdir, "regulon", "scenic_regulons_RSS.pdf"))
lapply(rssPlots, function(i) i$plot$plot)
dev.off()


###################
#### 7) ssGSEA ####
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


msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'))#              # curated gene sets
# 'C5'=list('GO:BP', 'GO:CC', 'GO:MF'))#, # ontology gene sets
# 'C7'=list('IMMUNESIGDB'))             # immunologic signature gene sets
# 'C8'=list(NULL),                      # cell type signature gene sets
# 'custom'=list('custom'=gois))           # custom genesetsmsig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
mlvl <- 'C2'
sublvl <- "CP:REACTOME"
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


######################################
#### 7) Pairing Regulons with DEG ####
## Initialize
res_dds <- readRDS(file.path(outdir, "rds", "resl.rds"))
coldata <- readRDS(file=file.path(outdir, "rds", "coldata.rds"))
regulon_aucs <- readRDS(file=file.path(outdir, "regulon", "regulon_auc.rds"))
dir.create(file.path(outdir, "regulon_deg"), showWarnings = F)

# Read in Regulon information
regulon_targets <- readRDS(file.path("int", "2.5_regulonTargetsInfo.Rds"))

# For each celltype, merge the gene-LFC/pvals and the regulon-AUC/Deltas
regulon_targets_cts <- lapply(celltypes, function(celltype){
  # celltype <- 'MSM'
  # Gene-level LFC for DEGs
  lfcs <- res_dds[[celltype]]$res %>%
    select(gene, grep("^padj", colnames(.), value=T), 
           grep("log2FoldChange", colnames(.), value=T)) %>%
    rename_with(., ~ gsub("-", "_", .))
  
  # Regulon-level delta-AUCs
  regulons <- lapply(names(regulon_aucs), function(reg_id){
    reg_i <- regulon_aucs[[reg_id]]
    reg_i[[celltype]] %>%
      mutate(TF=gsub("(_extended)? \\(.*", "", Regulon)) %>%
      filter(!duplicated(TF)) %>% 
      select(TF, grep("delta", colnames(.))) %>%
      rename_with(., ~ gsub("tdln_ln$|auc$", reg_id, .), .cols=matches("delta", perl=T))
  }) %>%
    Reduce(function(x,y) full_join(x,y,by='TF'), .)
  
  # Merge the regulon-target with the LFCs and regulon-deltas
  regulon_targets_full <- regulon_targets %>% 
    left_join(lfcs, by='gene') %>%
    left_join(regulons, by='TF')
  
  return(regulon_targets_full)
})

# Set the list of regulons interested in
regulons <- .getRegulons(intx=reshape2::melt(regulon_aucs$interaction$MSM),
                         cisx=reshape2::melt(regulon_aucs$cis_tdln_ln$MSM)) %>%
  tail(., n=10) %>% 
  gsub("(_extended)? \\(.*", "", .) %>%
  unique()

grps <- list("grp1"=c('interaction', 'cis_tdln_ln'))
regulon_network <- lapply(names(regulon_targets_cts), function(celltype){
  ct_reg <- regulon_targets_cts[[celltype]]
  
  # ct_reg <- regulon_targets_cts$MSM
  ct_reg_subset <- ct_reg %>%
    filter(TF %in% regulons,
           highConfAnnot==TRUE) %>%
    select(TF, gene, NES, spearCor, grep(paste(grps[[1]], collapse="|"), colnames(.),value=T)) %>% 
    rename_with(., .fn=~paste0(., ".", celltype), .cols=matches("padj|log2|delta"))
  
  gpr <- ct_reg_subset %>% 
    group_by(TF) %>% 
    summarise(genes=paste(gene, collapse=","))
  nodes <- ct_reg_subset %>%
    select(gene, grep(paste(grps[[1]], collapse="|"), colnames(.),value=T)) %>%
    unique()
  edges <- ct_reg_subset %>% 
    select(TF, gene, spearCor, NES)
  
  return(list("edge"=edges, "node"=nodes, "genes_per_regulon"=gpr))
})
names(regulon_network) <- names(regulon_targets_cts)

for(celltype in celltypes){
  
  
  write.table(regulon_network$MSM$edge, 
              file=file.path(outdir, "regulon_deg", "regulon_deg.edge.tsv"),
              sep="\t", col.names=T, row.names = F, quote = F)
  full_join(regulon_network$MSM$node,
            regulon_network$SSM$node, by=c("gene")) %>%
    write.table(., 
              file=file.path(outdir, "regulon_deg", "regulon_deg.node.tsv"),
              sep="\t", col.names=T, row.names = F, quote = F)
}

###########################
#### 8) Final Figures: ####
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
# inf <- read.table(file.path(dirto, "inflamm_gene_2.txt"), header = F) %>%
#   mutate(V1=str_to_title(V1))
# antiinf <- read.table(file.path(dirto, "antiinflamm_gene_2.txt"), header = F) %>%
#   mutate(V1=str_to_title(V1))
# bus.antiinf <- read.table(file.path(dirto, "buscher_antiinflamm.txt"), header = F) %>%
#   mutate(V1=str_to_title(V1))
# mariejo_up <- read.table(file.path(dirto, "mariejo_up.csv"), header = F) %>%
#   mutate(V1=str_to_title(V1))
# mariejo_dn <- read.table(file.path(dirto, "mariejo_down.csv"), header = F) %>%
#   mutate(V1=str_to_title(V1))
# barton <- read.table(file.path(dirto, "barton_2017_fig3d.csv"), header = F) %>%
#   mutate(V1=str_to_title(V1))
bus.inf <- read.table(file.path(dirto, "buscher_inflamm.txt"), header = F) %>%
  mutate(V1=str_to_title(V1))
rahul_dn <- read.table(file.path(dirto, "rahul_up.csv"), header = F) %>%
  mutate(V1=str_to_title(V1)) # Rahul's test conditions are reversed from ours
rahul_up <- read.table(file.path(dirto, "rahul_down.csv"), header = F) %>%
  mutate(V1=str_to_title(V1)) # Rahul's test conditions are reversed from ours

genesets <- list("Buscher.Inflammatory"=bus.inf$V1, 
                    "Rahul_down"=rahul_dn$V1, "Rahul_up"=rahul_up$V1)
# genesets <- list("Buscher.Inflammatory"=bus.inf$V1, 
#                  "Rahul.Tolerogenic"=unique(c(rahul_dn$V1, rahul_up$V1)))
                 #"Buscher.Anti-inflammatory"=bus.antiinf$V1,
                 #"MarieJo_M2_down"=mariejo_dn$V1, "MarieJo_M2_up"=mariejo_up$V1,
                 #"Barton_2017"=barton$V1

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


#########################################
#### 7) WGCNA: Coexpression analysis ####
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
# normalized_array <- t(assay(dds))
normalized_array <- t(vst(assay(dds)))
storage.mode(normalized_array) <- 'numeric'
sft <- pickSoftThreshold(normalized_array,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

## Computer the network modules and eigengenes from the adjacency matrix using TOM signed
pwr <- 18
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


## C) Correlating eigengenes to metadata ---------------------------------------
## Create the design matrix from the `time_point` variable
group <- 'treatment' # lnstatus, treatment, celltype, condition
des_mat <- model.matrix(as.formula(paste0("~ coldata$", group)))
stats_df <- linearModelEigen(t(module_eigengenes), des_mat)

## Combine eigenegenes to the metadata
module_df <- module_eigengenes %>%
  tibble::rownames_to_column("Sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(coldata %>%
                      tibble::rownames_to_column("Sample") %>%
                      dplyr::select(Sample, lnstatus, treatment, celltype, condition),
                    by = c("Sample" = "Sample"))


## D1) Associating eigengene modules with significant interaction DEGs ---------
res_module_oras <- lapply(names(res_dds), function(res_id){
  res <- res_dds[[res_id]]$res
  # Merge eigengene modules with p-values per MSM interaction
  gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
    dplyr::mutate(module = paste0("ME", module)) %>% 
    mutate("hugo"=ens2sym_ids[gene]) %>%
    mutate("entrez"=ens2entrez_ids[gene]) %>%
    inner_join(res %>%
                 as.data.frame() %>%
                 select(c("ens", grep("(padj)|(log2Fold)", colnames(.), value=T))),
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
  
  
  ## Select the genes of module X
  module_ids <- sapply(sig_modules, function(i) unique(i$module))
  module_oras <- lapply(setNames(module_ids,module_ids), function(me){
    print(paste0(res_id, ": ", me, "..."))
    module_genes <- gene_module_key %>%
      dplyr::filter(module == me)
    
    oras <- iterateMsigdb(species='Mus musculus', fun=oraFun, 
                          entrez_genes=module_genes$entrez)
    oras <- lapply(oras, summarizeOra, keep_genes=TRUE, qcutoff=0.15)
    
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
  
  do.call(rbind, module_oras) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., "module") %>%
    mutate(module = gsub("\\.[0-9]*$", "", module))
})
names(res_module_oras) <- names(res_dds)
saveRDS(res_module_oras, file=file.path(outdir, "network", "lfc_module_oras.rds"))

## D2) Get frequent word combinations of Gene sets for each module ---------
if(!exists("res_module_oras")) res_module_oras <- readRDS(file.path(outdir, "network", "lfc_module_oras.rds"))
x <- split(res_module_oras$MSM, res_module_oras$MSM$module)
names(x)

n <- 6
freqWords <- function(textv, split="_", word_pairs=F, n=5){
  if(length(textv)==0) return('')
  words <- unlist(strsplit(textv, split=split))
  if(word_pairs) words <- paste(words, lead(words), sep=":")
  paste(names(
    head(sort(table(words), decreasing = T), n)), 
    collapse=",")
}

pairwise_words <- lapply(res_module_oras, function(module_i){
  if(nrow(module_i)==0) return(NULL)
  spl_module_i <- split(module_i, module_i$module)
  freq_word_pairs <- sapply(spl_module_i, function(i){
    freqWords(textv=i$ID, split="_", word_pairs=TRUE, n=n)
  })
  freq_words <- sapply(spl_module_i, function(i){
    freqWords(textv=i$ID, split="_", word_pairs=FALSE, n=n)
  })
  lfc_delta <- sapply(spl_module_i, function(i){ c("mean"=mean(i$mean_lfc_int), 
                                                   "sd"=mean(i$mean_lfc_int))})
  t(lfc_delta) %>%
    as.data.frame() %>%
    mutate("freq_w"=freq_words,
           "freq_w_pairs"=freq_word_pairs)
})


## E) Cytoscape export from network ------------
lfc <- lapply(names(res_dds), function(res_id) {
  .addID <- function(x) paste0(x, ".", res_id)
  res_dds[[res_id]]$res %>%
    dplyr::select(., grep("(log2Fold)|(ens)|(gene)", colnames(.), value=T)) %>%
    rename_with(., .addID, starts_with("log2"))
}) %>% 
  Reduce(function(x,y) merge(x,y,by=c('ens', 'gene'), all=T), .)

quantile_edge <- 0.9
TOM = TOMsimilarityFromExpr(normalized_array, power = pwr)
modules <- colnames(bwnet$MEs) %>%
  gsub("ME", "", .) %>%
  setNames(., colnames(bwnet$MEs))

cyts <- lapply(modules, function(module){
  print(module)
  in_module <- which(bwnet$colors == as.integer(module))

  datexpr = normalized_array[,in_module]
  TOM_module = TOMsimilarityFromExpr(datexpr, power = pwr, networkType = "signed", TOMType="signed");
  probes = colnames(datexpr)
  dimnames(TOM_module) = list(probes, probes)
  
  sym_ids <- setNames(ens2sym_ids[probes], probes)
  sym_ids[which(is.na(sym_ids))] <- probes[which(is.na(sym_ids))]
  
  cyt = exportNetworkToCytoscape(TOM_module,
                                 weighted = TRUE,
                                 threshold = 0.10)
  cyt = exportNetworkToCytoscape(TOM_module,
                                 edgeFile = file.path("manual", "network", 
                                                      paste("cytoscape-edges-", 
                                                            paste(module, collapse="-"), 
                                                            ".txt", sep="")),
                                 nodeFile = file.path("manual", "network", 
                                                      paste("cytoscape-nodes-", 
                                                            paste(module, collapse="-"), 
                                                            ".txt", sep="")),
                                 weighted = TRUE,
                                 threshold = quantile(cyt$edgeData$weight, quantile_edge),
                                 nodeNames = probes,
                                 altNodeNames = sym_ids,
                                 nodeAttr = lfc[match(probes, lfc$ens), 
                                                    grep("log2FoldChange", colnames(lfc))] %>%
                                   round(., 3));
  return(cyt)
})

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






