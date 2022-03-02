library(reshape2)
library(ggplot2)
library(TCGAbiolinks)
library(gridExtra)
library(cowplot)
source("~/git/mini_projects/ST2_IL33/code/TCGAanalyze_SurvivalKM2.R")

PDIR='/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/melanoma_meta/'
outdir <- file.path("output", "Cell_2017")
setwd(PDIR)

expr <- read.csv(file.path("data/RNA/Cell_2017", "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"),
                 header=TRUE)
meta <- read.csv(file.path("data/RNA/Cell_2017", "stable2.csv"), header=TRUE)
list.gene <- c('IL1RL1'='9173',  'IL33'='90865',
               'PDCD1_1' = "5133", 'CTLA4' = "1493", 'TOX' = "9760", 'KLRG1' = "10219")

## Preformat of expression samples
rownames(expr) <- expr$X
expr <- expr[,-1]
ids <- data.frame("ids"=colnames(expr))
ids$Patient <- gsub("_.*", "", ids$ids)

## Format the metadata to match TCGA style
dead_idx <- which(meta$Dead_Alive)
alive_idx <- which(!meta$Dead_Alive)
meta$vital_status <- 'Alive'
meta$vital_status[dead_idx] <- 'Dead'
meta$days_to_death <- meta$Time_to_death
meta$days_to_death[alive_idx] <- -Inf
meta$days_to_last_follow_up <- meta$Time_to_death
meta$days_to_last_follow_up[dead_idx] <- -Inf

## Filter metadata and split based on pre or on-treatment samples
meta_sub <- merge(meta, ids, by='Patient', all=TRUE)
meta_sub$bcr_patient_barcode <- meta_sub$ids
meta_sub <- meta_sub[-which(is.na(meta_sub$ids)),]
meta_sub$treatment <- gsub("^.*_(On|Pre)_.*$", "\\1", meta_sub$ids)
meta_sub_spl <- split(meta_sub, is.na(meta_sub$Cohort))
meta_treatment <- split(meta_sub_spl$`FALSE`, meta_sub_spl$`FALSE`$treatment)
meta_patients <- split(meta_sub_spl$`FALSE`, meta_sub_spl$`FALSE`$Patient)


## Generate survival curve analysis
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
skms <- lapply(meta_treatment, function(meta_t){
  expr_t <- log2(expr[,match(meta_t$ids, colnames(expr))]+1)
  
  pdf(file.path(outdir, paste0("SC_", unique(meta_t$treatment), "-treatment.pdf")))
  skm <- TCGAanalyze_SurvivalKM2(meta_t, expr_t, list.gene, Survresult=TRUE,
                                 threshcuts = c(0.33, 0.66), caption=unique(meta_t$treatment),
                                 p.cut = 1, dataset='Cell_2017')
  dev.off()
  names(skm) <- names(list.gene)[match(names(skm),list.gene)]
  return(skm)
})
skms <- lapply(skms, function(i) do.call(rbind, i))
save(skms, file=file.path(outdir, "skms.rda"))


expr_melt <- lapply(meta_treatment, function(meta_t){
  # isolate the genes and response of interest
  meta_resp <- split(meta_t, meta_t$Response)
  expr_resp <- lapply(meta_resp, function(i){
    expr[as.character(list.gene),i$ids,drop=FALSE]
  })
  
  # Create a melted data structure of gene expression by response by gene
  expr_gene_resp <- lapply(list.gene, function(gene){
    gene_expr <- sapply(expr_resp, function(expr_t){
      as.numeric(expr_t[as.character(gene),])
    })
    data.frame("class"=rep(names(gene_expr), sapply(gene_expr, length)),
               "expr"=log2(unlist(gene_expr)+1),
               "gene"=rep(gene, length(unlist(gene_expr))))
  })
  expr_gene_resp <- as.data.frame(do.call(rbind, expr_gene_resp))
  
  expr_gene_resp
})


## Compute pairwise statistics using Tukeys HSD
pairwise_stats <- lapply(expr_melt, function(expr_i){
  pairwise_stat <- lapply(split(expr_i, expr_i$gene), function(ei){
    res_aov <- aov(as.formula(paste0("expr ~ class")), data = ei)
    pairwise_tukey  <- TukeyHSD(res_aov)[[1]]
    pairwise_tukey
  })
  names(pairwise_stat) <- names(list.gene)[match(names(pairwise_stat), list.gene)]
  pairwise_stat
})

## Plot the Expression per response and the pairwise tables
#/INTERACTIVE-SESSION/#
expr_id <- names(expr_melt)[1]  # 'Pre' or 'On'
expr_melt_i <- expr_melt[[expr_id]]
scaleFUN <- function(x) sprintf("%.1f", x)

## Generate violin+boxplots of expression per RECIST criteria
dps <- lapply(split(expr_melt_i, expr_melt_i$gene), function(expr_i){
  gene <- unique(expr_i$gene)
  gene <- names(list.gene)[match(gene, list.gene)]
  
  expr_i$class <- factor(expr_i$class, levels=c('NE', 'CR', 'PR', 'SD', 'PD'))
  expr_dp <- ggplot(expr_i, aes(x=class, y=expr)) + 
    geom_violin(trim=FALSE, fill="#74a9cf", scale = "width") +
    labs(title=paste0(gene, " gene expression"),
         x="RECIST", y = "Log2(FPKM + 1)") +
    geom_boxplot(width=0.2, fill="white") +
    theme_minimal() +
    scale_y_continuous(labels=scaleFUN) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(expr_dp)
})

## Generate tables of tukeys HSD pariwise stats for RECIST criteria
ord_pairwise <- pairwise_stats[[expr_id]][match(names(dps), 
                                                list.gene[names(pairwise_stats[[expr_id]])])]
tps <- lapply(ord_pairwise, function(stats_i){
  mytheme <- gridExtra::ttheme_default(
    core = list(padding = unit(c(2.5, 2.5), "mm")))
  stats_i <- stats_i[order(stats_i[,4]),]
  tbl <- tableGrob(round(stats_i, 3), theme = mytheme, rows = rownames(stats_i))
  return(tbl)
})

pdf(file.path(outdir, paste0(expr_id, "_recist-gene-expr.pdf")))
plot_grid(dps[[1]],tps[[1]], nrow=2, ncol=1)
plot_grid(dps[[2]],tps[[2]], nrow=2, ncol=1)
plot_grid(dps[[3]],tps[[3]], nrow=2, ncol=1)
plot_grid(dps[[4]],tps[[4]], nrow=2, ncol=1)
plot_grid(dps[[5]],tps[[5]], nrow=2, ncol=1)
plot_grid(dps[[6]],tps[[6]], nrow=2, ncol=1)
dev.off()