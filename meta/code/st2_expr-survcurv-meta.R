###########################
#### Read in TCGA Data ####
library(reshape2)
library(ggplot2)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(gridExtra)
library(cowplot)
library(scales)
library(MASS)
library(taRifx)
source("~/git/mini_projects/ST2_IL33/code/TCGAanalyze_SurvivalKM2.R")

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcga'
outdir <- file.path(PDIR, 'results')
icdo3_codes <- file.path(PDIR, "ref", "sitetype.icdo3.d20210607.csv")
icd10_codes <- file.path(PDIR, "ref", "icd10_codes.csv")
tcgadir <- file.path(PDIR, 'input/data/fpkm/obj')
clindir <- file.path(PDIR, 'input/data/clinical')
visualize <- FALSE

icdo3 <- read.csv(icdo3_codes, sep=",", header=TRUE, fill=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
icd10 <- read.csv(icd10_codes, sep=",", header=FALSE, fill=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
icd10$V3 <- paste0(icd10$V1, ".", icd10$V2)
projects <- c("TCGA-BRCA", "TCGA-GBM", "TCGA-OV", "TCGA-LUAD", "TCGA-UCEC",
              "TCGA-KIRC", "TCGA-HNSC", "TCGA-LGG", "TCGA-THCA", "TCGA-LUSC",
              "TCGA-PRAD", "TCGA-SKCM", "TCGA-COAD", "TCGA-STAD", "TCGA-BLCA",
              "TCGA-LIHC", "TCGA-CESC", "TCGA-KIRP", "TCGA-SARC", "TCGA-LAML", 
              "TCGA-PAAD", "TCGA-ESCA", "TCGA-PCPG", "TCGA-READ", "TCGA-TGCT",
              "TCGA-THYM", "TCGA-KICH", "TCGA-ACC", "TCGA-MESO", "TCGA-UVM",
              "TCGA-DLBC", "TCGA-UCS", "TCGA-CHOL")
workflow_type <- list("counts"="HTSeq - Counts", "fpkm"="HTSeq - FPKM-UQ")

list.gene <- c('IL1RL1' = 'ENSG00000115602', # IL1RL1
               'PDCD1_1' = "ENSG00000188389", # PDCD1
               #'PDCD1_2' = "ENSG00000276977", # PDCD1
               'CTLA4' = "ENSG00000163599", # CTLA4
               'TOX' = "ENSG00000198846", # TOX
               'KLRG1' = "ENSG00000139187") # KLRG1
assay_type <- "HTSeq - FPKM-UQ"

## Read in TCGA expression and metadata
setwd(tcgadir)
tcga_files <- list.files(path = tcgadir, pattern="obj.rds$")
objs <- lapply(tcga_files, readRDS)
names(objs) <- gsub("\\..*", "", tcga_files)
exprs <- lapply(objs, function(i) assay(i, assay_type))

## Read in TCGA clinical data
setwd(clindir)
clin_files <- list.files(path = clindir, pattern="clinical.rds$")
clins <- lapply(clin_files, readRDS)
names(clins) <- gsub("\\..*", "", clin_files)

###################
#### Functions ####
annotateIcdCode <- function(p_df, icd_df, icd_code_col, icd_anno_col){
  m_idx <- match(p_df$Var2, icd_df[,icd_code_col])
  p_df$Var2 <- paste0("[", p_df$Var2, "] ",
                      icd_df[m_idx, icd_anno_col])
  
  m_idx <- match(p_df$Var1, icd_df[,icd_code_col])
  p_df$Var1 <- paste0("[", p_df$Var1, "] ",
                      icd_df[m_idx, icd_anno_col])
  p_df
}

###########################
#### Survival Analysis ####
## Obtain the survival curves for each cancer type in TCGA
survivals <- lapply(projects, function(proj){
  # get Log2 FPKM + 1
  print(paste0(proj, " (", grep(proj, projects), "/", length(projects), ")"))
  log2exprs <- log2(exprs[[proj]] + 1)
  
  # Run survival analysis and get log-rank p-values
  
  dir.create(file.path(outdir, "survival_curves"), showWarnings = FALSE)
  if(visualize) pdf(file.path(outdir, "survival_curves", paste0(proj, "_skm.pdf")))
  skm <- TCGAanalyze_SurvivalKM2(clins[[proj]], log2exprs, list.gene, Survresult=TRUE,
                         threshcuts = c(0.25, 0.50, 0.75), caption=proj, p.cut = 1)
  if(visualize) dev.off()
  skm <- as.data.frame(do.call(rbind, skm))
  rownames(skm) <- list.gene
  skm$HUGO <- names(list.gene[match(rownames(skm), list.gene)])
  return(skm)
})
# Format the data structures
names(survivals) <- projects
st2_pvals <- as.data.frame(t(sapply(survivals, function(s) s[list.gene['IL1RL1'],])))
save(st2_pvals, file=file.path(outdir, "survival_curves", "survival_pvals.rda"))

################################
#### Co-expression analysis ####
# Plot IL1RL1 expression against all other marker genes and return correlation + p.value
coexprs <- lapply(projects, function(proj){
  # get Log2 FPKM + 1
  print(paste0(proj, " (", grep(proj, projects), "/", length(projects), ")"))
  log2exprs <- log2(exprs[[proj]] + 1)
  
  dir.create(file.path(outdir, "coexpression"), showWarnings = FALSE)
  if(visualize) pdf(file.path(outdir, "coexpression", paste0(proj, "_coexpr.pdf")))
  if(visualize) par(mfrow=c(2,2))
  
  k=6
  my.cols <-  colorRampPalette(c("#f03b20", "#feb24c", "#ffeda0"))(k)
  
  st2_expr <- as.numeric(log2exprs[list.gene['IL1RL1'],])  
  # Compute correlation between IL1RL1 and all other marker genes
  cors <- sapply(names(list.gene[-grep("IL1RL1", names(list.gene))]), function(gene_y){
    gene_y_expr <- as.numeric(log2exprs[list.gene[gene_y],])
    z <- kde2d(st2_expr, gene_y_expr, n=50) # Computer 2D KDE
    cor <- cor.test(st2_expr, gene_y_expr)
      
    # Scatterplot of IL1RL1 to marker gene, overlaid with a contour plot
    plot(x=st2_expr, y=gene_y_expr, pch=16,
         col=alpha("black", 0.4), xlab='IL1RL1', ylab=gene_y,
         main=paste0("cor=", round(cor$estimate, 3), 
                     "; p=",round(cor$p.value, 3)))
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    return(c("r"=cor$estimate, "p"=cor$p.value))
  })
  if(visualize) dev.off()
  return(as.data.frame(t(cors)))
})
names(coexprs) <- projects

# Separate correlations by marker genes
cor_genes <- lapply(names(list.gene)[-1], function(geneid){
  sapply(coexprs, function(i){ i[geneid,] })
})
names(cor_genes) <-  names(list.gene)[-1]

# Generate a barplot of all correlations for each cancer type, separated by marker genes
if(visualize) {
  pdf(file.path(outdir, "coexpression", "cor_ctypes.pdf"))
  par(mfrow=c(5,1))
  sapply(names(cor_genes), function(geneid){
    par(mar=c(1, 4.1, 2, 2.1))
    my.cols <-  c("<0.001"="#d73027", "0.001-0.05"="#fc8d59", 
                  "0.05-0.1"="#e0f3f8", "nosnsig"="grey")
    p_idx <- as.integer(cut(as.numeric(cor_genes[[geneid]]['p',]), 
                            c(0, 0.001, 0.05, 0.1, 1)))
    barplot(unlist(cor_genes[[geneid]]['r.cor',]),  col=my.cols[p_idx], xaxt='n',
            ylim=c(-1,1), las=2, cex.names = 0.6, ylab="cor", main=geneid)
  })
  par(mar=c(10.1, 4.1, 0, 2.1))
  barplot(setNames(rep(0, length(projects)), projects), ylim=c(0,0.1), las=2)
  dev.off()
}  

#################################################################
#### Isolate Metadata and GOI Expression across cancer types ####
load(file.path(outdir, "survival_curves", "survival_pvals.rda"))
groups <- c('ajcc_pathologic_stage.x', 'tumor_stage.x', 'morphology.x',
            'ajcc_pathologic_t.x', 'ajcc_pathologic_n.x', 'ajcc_clinical_m.x',
            'classification_of_tumor.x', 'icd_10_code.x', 'tumor_grade.x',
            'progression_or_recurrence.x', 'race.x', 'gender.x')
            # 'state', 'tissue_type', 'ajcc_pathologic_stage', 'tumor_stage', 
            # 'ajcc_pathologic_t', 'ajcc_pathologic_n', 'ajcc_pathologic_m',
            # 'morphology', 'tumor_grade', 'progression_or_recurrence', 
            # 'paper_Tumor.Type', 'paper_pathologic_stage', 'paper_Tumor_Grade')
goi <-'IL1RL1'
st2_dat <- lapply(setNames(names(objs), names(objs)), function(id) {
  i <- objs[[id]]
  clin_i        <- clins[[id]]
  metadata      <- colData(i)
  clin_metadata <- merge(clin_i, metadata, by.x='submitter_id', 
                         by.y='patient', all=TRUE)
  metadata      <- clin_metadata
    
  expr <- assay(i, assay_type)
  expr <- expr[match(list.gene[goi], rownames(expr)),]
  expr <- data.frame("id"=names(expr),
                     "expr"=log2(expr + 1))
  expr_metadata <- merge(expr, metadata[,which(colnames(metadata) %in% c('barcode', groups))], 
                         by.x='id', by.y='barcode')
  expr_metadata$index <- c(1:nrow(expr_metadata))
  expr_metadata$ctype <- id
  expr_metadata
})
# Reduce expression and metadat a by cancer types
st2_exprs <- Reduce(function(x,y) merge(x,y,by='index', all=TRUE), lapply(st2_dat, function(i) i[,c('index', 'expr')]))
st2_exprs <- st2_exprs[,-1]
colnames(st2_exprs) <- names(st2_dat)
st2_exprs <- st2_exprs[,rownames(st2_pvals)]
st2_dat <- lapply(st2_dat, function(i) i[,-grep("index", colnames(i))])

###########################################
#### Visualize Expression and Survival ####
# Reorder cancer types based on survival analysis pvals
melt_st2_exprs <- melt(as.data.frame(st2_exprs))
melt_st2_exprs <- melt_st2_exprs[-which(is.na(melt_st2_exprs$value)),]
melt_st2_exprs$value <- as.numeric(as.character(melt_st2_exprs$value))

# Visualize gene expression in FPKM+1 ordered by survival pvals
scaleFUN <- function(x) sprintf("%.1f", x)
expr_dp <- ggplot(melt_st2_exprs, aes(x=variable, y=value)) + 
  geom_violin(trim=FALSE, fill="#74a9cf", scale = "width") +
  labs(title="IL1RL1 gene expression",x="TCGA Projects", y = "Log2(FPKM + 1)") +
  geom_boxplot(width=0.2, fill="white") +
  theme_minimal() +
  scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

if(visualize){
  pdf(file.path("~/xfer", "expr.pdf"), width=12, height=5)
  expr_dp
  dev.off()
}


# Visualize survival pval ordered by survival pvals
st2_pvals$groups <- factor(rownames(st2_pvals), levels=rownames(st2_pvals))
st2_pvals$sig <- st2_pvals$pvalue < 0.05
surv_dp <- ggplot(st2_pvals, aes(x=groups, y=pvalue, fill=sig)) + 
  geom_bar(stat="identity", color="black") + 
  scale_fill_manual(values=c("#ef8a62", "#67a9cf")) +
  labs(title="IL1RL1 Survival p-values",x="TCGA Projects", y = "Log-rank p-value") +
  theme_minimal() +
  coord_cartesian(ylim=c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

if(visualize){
  pdf(file.path(outdir, "expr_skm.pdf"), width=12, height=9)
  grid.arrange(expr_dp + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank()), 
               surv_dp, 
               ncol = 1, nrow = 2)
  dev.off()
}  


##########################################
#### Expression and Metadata analysis ####
# Report all metadata categories
# x <- lapply(st2_dat, function(st2){
#   grpids <- colnames(st2)[-c(1:2)]
#   categories <- lapply(setNames(grpids, grpids), function(g){
#     table(st2[,g])
#   })
#   return(categories)
# })

# Run an ANOVA analysis of IL1RL1 categorized by each metadata column
aov_pvals <- lapply(st2_dat, function(st2){
  grpids <- colnames(st2)[-c(1:2)]
  aovs <- lapply(setNames(grpids, grpids), function(g){
    res.aov <- tryCatch({
      aov(as.formula(paste0("expr ~ ", g)), data = st2)
    }, error=function(e){NA})
    res.aov
  })
  pval <- sapply(aovs, function(i)  tryCatch(summary(i)[[1]]$`Pr(>F)`[1], error=function(e){NA}))
  data.frame("category"=grpids, "p"=pval)
})
pvals <- Reduce(function(x,y) as.data.frame(merge(x,y,by='category', all=TRUE)), aov_pvals)
pvals <- matrix(unlist(pvals), ncol=ncol(pvals))
rownames(pvals) <- pvals[,1]
pvals <- pvals[,-1]
colnames(pvals) <-names(st2_dat)
storage.mode(pvals) <- 'numeric'
pvals <- round(pvals, 4)
qvals <- round(matrix(p.adjust(pvals, method='fdr'), ncol=ncol(pvals)), 4)
colnames(qvals) <- colnames(pvals)
rownames(qvals) <- rownames(pvals)

## Heatmap visualization of ANOVA analysis
meta_pvals <- qvals[,rownames(st2_pvals)]
melt_meta_pvals <- melt(meta_pvals)
meta_dp <- ggplot(melt_meta_pvals, aes(Var2, Var1, fill= value)) + 
  geom_tile() +
  labs(title="ANOVA p-value for IL1RL1 expression",x="TCGA Projects", y = "Metadata") +
  scale_fill_gradientn(colours = c("#cb181d","#2171b5", "#6baed6", "#c6dbef", "#f7fbff"), 
                       values = rescale(c(0, 0.05, 0.1, 0.15, 1)),
                       guide = "colorbar", limits=c(0,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

if(visualize){
  pdf(file.path("~/xfer", "expr_skm_meta.pdf"), width=12, height=14)
  plot_grid(expr_dp + theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()), 
            surv_dp + theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()),
            meta_dp, 
            ncol=1, align="v")
  dev.off()
}


###############################################
#### Pairwise tests of sig. ANOVA metadata ####
sig_meta <- melt_meta_pvals[which(melt_meta_pvals$value < 0.05),]
sig_meta$Var1 <- as.character(sig_meta$Var1)
sig_meta$Var2 <- as.character(sig_meta$Var2)
rownames(sig_meta) <- paste0(sig_meta$Var2, "_", sig_meta$Var1)

sig_pairwise <- apply(sig_meta, 1, function(sig_i){
  ctype <- sig_i['Var2']
  metad <- sig_i['Var1']
  
  res_aov <- aov(as.formula(paste0("expr ~ ", metad)), data = st2_dat[[ctype]])
  pairwise_tukey  <- TukeyHSD(res_aov)[[1]]
  pairwise_t      <- pairwise.t.test(st2_dat[[ctype]]$expr, st2_dat[[ctype]][,metad], 
                                     p.adjust.method = "fdr")
  
  
  # Combine tukey, pairwise-t, and group sizes
  pairwise_tukey <- as.data.frame(pairwise_tukey)
  tukey_ids <- do.call(rbind, strsplit(rownames(pairwise_tukey), split="-"))
  colnames(tukey_ids) <- c('Var1', 'Var2')
  pairwise_tukey <- cbind(tukey_ids, pairwise_tukey)
  
  n_table <- melt(table(st2_dat[[ctype]][,metad]))
  
  pairwise_p <- merge(pairwise_tukey, n_table, by.x='Var2', by.y='Var1')
  pairwise_p <- merge(pairwise_p, n_table, by.x='Var1', by.y='Var1')
  nc <- ncol(pairwise_p)
  colnames(pairwise_p)[c(nc-1,nc)] <- c('n_v2', 'n_v1')
  pairwise_p <- pairwise_p[order(pairwise_p[,'p adj']),]
  
  if(grepl("icd_10", metad)){
    pairwise_p <- annotateIcdCode(pairwise_p, icd10, 'V3', 'V4')
  } else if(grepl("morphology", metad)){
    pairwise_p <- annotateIcdCode(pairwise_p, icdo3, 'Histology/Behavior', 
                                  'Histology/Behavior Description')
  }

  return(pairwise_p)
})


dir.create(file.path(outdir, "meta"), showWarnings = FALSE)
sig_q <- lapply(setNames(names(sig_pairwise),names(sig_pairwise)), function(id){
  raw_p <- sig_pairwise[[id]]
  sig_p <-raw_p[raw_p[,'p adj'] < 0.1,]
  
  write.table(raw_p, file=file.path(outdir, "meta", paste0("RAW_", gsub(".x$", "", id), ".csv")),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  write.table(sig_p, file=file.path(outdir, "meta", paste0("SIG_", gsub(".x$", "", id), ".csv")),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  sig_p
})


# 
# kws <- lapply(setNames(groups, groups), function(g){
#   res.kw <- tryCatch({
#     kruskal.test(as.formula(paste0("expr ~ ", g)), data = df)
#   }, error=function(e){NA})
#   res.kw
# })
# aovs <- lapply(setNames(groups, groups), function(g){
#   res.aov <- tryCatch({
#     aov(as.formula(paste0("expr ~ ", g)), data = df)
#   }, error=function(e){NA})
#   res.aov
# })
# pvals <- sapply(aovs, function(g){
#   tryCatch({
#     summary(g)[[1]]$`Pr(>F)`[1]
#   }, error=function(e){NA})
# })
# paired <- lapply(aovs, function(g){
#   tryCatch({
#     TukeyHSD(res.aov)
#   }, error=function(e){NA})
# })
# paired_t <- lapply(setNames(groups, groups), function(g){
#   pairwise.t.test(df$expr, df[,g], p.adjust.method = "fdr")
# })
# paired[names(which(pvals < 0.05))]
# paired_t[names(which(pvals < 0.05))]
# 
# summary(res.aov)
# 
# # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
# BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)