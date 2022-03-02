library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/melanoma_meta'
outdir <- 'output/Regev_2018'
setwd(PDIR)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path("data/scRNA/Regev_2018", "GSE115978_cell.annotations.csv.gz"),
                 header=TRUE)
expr <- read.csv(file.path("data/scRNA/Regev_2018", "GSE115978_tpm.csv.gz"),
                 header=TRUE)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
list.gene <- c('IL1RL1', 'IL33', 'PDCD1', 'CTLA4', 'TOX', 'KLRG1')

######################################################################
#### Plot Gene Expression by treatment for malignant/nonmalignant ####
## Split expression by treatment status and malignant
meta_celltype <- split(meta, meta$cell.types == 'Mal')
celltype_legend <- c('FALSE'='NonMal', 'TRUE'='Mal')
names(meta_celltype) <- celltype_legend[as.character(names(meta_celltype))]
meta_treat <- lapply(meta_celltype, function(m) split(m, m$treatment.group))

## Create a melted data structure of expression by gene and treatment group
expr_treat <- lapply(meta_treat, function(meta_ct){
  melt_expr <- lapply(names(meta_ct), function(meta_id){
    m <- meta_ct[[meta_id]]
    melt_e <- melt(t(expr[list.gene,m$cells]))
    colnames(melt_e) <- c('cell', 'gene', 'expr')
    melt_e$treatment <- rep(meta_id, nrow(melt_e))
    return(melt_e)
  })
  return(as.data.frame(do.call(rbind, melt_expr)))
})

gps <- lapply(expr_treat, function(melt_e){
  melt_e$treatment <- factor(melt_e$treatment, c('treatment.naive', 'post.treatment'))
  melt_e$expr <- log2(as.numeric(melt_e$expr)+1)
  expr_dp <- ggplot(melt_e, aes(x=gene, y=expr, fill=treatment)) + 
    geom_violin(trim=FALSE, scale = "width") +
    labs(title="Gene expression per cell type and treatment",
         x="Genes", y = "Log2(TPM + 1)") +
    #geom_boxplot(width=0.2, fill="white") +
    ylim(0,5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(expr_dp)
})

pdf(file.path(outdir, "scRNA-expr.pdf"))
gps[[1]]
gps[[2]]
dev.off()

##############################################################
#### Plot ST2/IL33 expression by treatment and cell types ####
## Split expression by treatment status and malignant
meta_celltype <- split(meta, meta$cell.types)
meta_treat <- lapply(meta_celltype, function(m) split(m, m$treatment.group))

## Create a melted data structure of expression by gene and treatment group
expr_treat <- lapply(names(meta_treat), function(celltype_id){
  meta_ct <- meta_treat[[celltype_id]]
  melt_expr <- lapply(names(meta_ct), function(meta_id){
    m <- meta_ct[[meta_id]]
    melt_e <- melt(t(expr[list.gene,m$cells]))
    colnames(melt_e) <- c('cell', 'gene', 'expr')
    melt_e$treatment <- rep(meta_id, nrow(melt_e))
    return(melt_e)
  })
  melt_expr <- as.data.frame(do.call(rbind, melt_expr))
  melt_expr$celltype <- rep(celltype_id, nrow(melt_expr))
  
  return(melt_expr)
})
expr_treat <- as.data.frame(do.call(rbind, expr_treat))
expr_gene <- split(expr_treat, expr_treat$gene)

gps <- lapply(expr_gene, function(melt_e){
  gene <- as.character(unique(melt_e$gene))
  
  ## Preformat the data
  melt_e$treatment <- factor(melt_e$treatment, c('treatment.naive', 'post.treatment'))
  melt_e$celltype <- factor(melt_e$celltype, levels=names(meta_celltype))
  melt_e$expr <- log2(as.numeric(melt_e$expr)+1)
  
  ## Plot the expression of the celltypes
  expr_dp <- ggplot(melt_e, aes(x=celltype, y=expr, fill=treatment)) + 
    geom_violin(trim=FALSE, scale = "width") +
    labs(title=gene,
         x="Cell type", y = "Log2(TPM + 1)") +
    #geom_boxplot(width=0.2, fill="white") +
    ylim(0,5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ## Plot the statistics of the celltype groupings
  df <- lapply(split(melt_e, melt_e$celltype), function(ct){
    df <- lapply(split(ct, ct$treatment), function(treat){
      data.frame('metric'=c('mean', 'median', 'pct'),
                 'value'=c(mean(treat$expr), 
                           median(treat$expr), 
                           sum(treat$expr > 0)/nrow(treat)))
    })
    df <- as.data.frame(do.call(rbind, df))
    df$treatment <- gsub("\\.[0-9]+$", "", rownames(df))
    rownames(df) <- NULL
    return(df)
  })
  df <- as.data.frame(do.call(rbind, df))
  df$celltype <- gsub("\\.[0-9]+$", "", rownames(df))
  stat_dp <- ggplot(df, aes(x=treatment, y=value, fill=metric)) + 
    facet_wrap(~celltype, scale="free") +
    geom_bar(stat='identity', position=position_dodge()) +
    labs(title=gene, x="treatment", y = "value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  return(list("violin"=expr_dp, "barplot"=stat_dp))
})

pdf(file.path(outdir, "scRNA-celltype.pdf"), width = 12)
plot_grid(gps[[1]]$violin, gps[[1]]$barplot, nrow=1, ncol=2)
plot_grid(gps[[2]]$violin, gps[[2]]$barplot, nrow=1, ncol=2)
plot_grid(gps[[3]]$violin, gps[[3]]$barplot, nrow=1, ncol=2)
plot_grid(gps[[4]]$violin, gps[[4]]$barplot, nrow=1, ncol=2)
plot_grid(gps[[5]]$violin, gps[[5]]$barplot, nrow=1, ncol=2)
plot_grid(gps[[6]]$violin, gps[[6]]$barplot, nrow=1, ncol=2)
dev.off()

######################################################
#### Plot GenesOfInterest expression per celltype ####
df <- sapply(split(expr_treat,list(expr_treat$gene, 
                             expr_treat$celltype,
                             expr_treat$treatment)), function(gene_expr){
  data.frame("gene"=unique(as.character(gene_expr$gene)),
             "treatment"=unique(gene_expr$treatment),
             "celltype"=unique(gene_expr$celltype),
             "expr"=mean(gene_expr$expr))
})
df <- as.data.frame(t(df))
df$gene <- as.character(unlist(df$gene))
df$treatment <- (unlist(df$treatment))
df$expr <- unlist(df$expr)
df$celltype <- unlist(df$celltype)
rownames(df) <- NULL

pdf(file.path(outdir, "scRNA-gene-expr.pdf"), width = 12)
ggplot(df, aes(gene, celltype, fill= expr)) + 
  facet_wrap(~treatment, scale="free") +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu")
dev.off()