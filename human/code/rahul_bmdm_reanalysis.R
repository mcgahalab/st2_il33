## Sara's re-analysis of Rahul's BMDM vs apBMDM 2016 data
library(cowplot)
library(igraph)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(ggrepel)

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

# Params
pdir <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/rahul_bmdm"
data_dedup <- file.path("results", "rsem")
dedir <- file.path(pdir, "results", "diffexp")
outdir <- file.path(pdir, "results", "manual")
dir.create(file.path(outdir, "rds"), showWarnings = F, recursive = F)
setwd(pdir)

# Read in deseq2 data
dds_main <- readRDS(file.path("results", "deseq2", "all.rds"))
dds <- dds_main

#######################
#### Functions Git ####
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/st2_il33/functions/msm_ssm_functions.R")
source("~/git/mini_projects/mini_functions/geneMap.R")

###########################
#### 0. Metadata setup ####
gm <- geneMap('Mus musculus')

pca_dir <- file.path(outdir, "pca")
dir.create(pca_dir, showWarnings = F, recursive = T)
max_pc <- 10  # Plot PC component 1 to X
visualize <- TRUE

# Set up the metadata
sample_ids <- colnames(dds_main) %>%
  gsub("^(.*)_([123]_[123])$", "\\1.\\2", .)
samples_meta <- strsplit(sample_ids, split="_") %>% do.call(rbind,.) %>%
  as.data.frame %>%
  rename_with(., ~c('group', 'id')) %>%
  mutate(sample=gsub("\\.", "_", sample_ids),
         group=gsub("\\.[123]", "", group)) %>%
  tibble::column_to_rownames(., 'sample')

############################
#### 1.a) PCA analysis  ####
dir.create(file.path(outdir, "qc"), showWarnings = F)

# obtain normalized counts
counts <- vst(dds, blind=T)
tcnts <- as.data.frame(t(assay(counts)))
pca <- prcomp(tcnts[grep("_2_", rownames(tcnts)),], scale=F)
max_pc_i <- min(c(max_pc, ncol(pca$x)))

# Format in any factors/metadata you want to plot by
pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc_i))],
                             "condition"=as.character(counts$condition)))
for(id in paste0("PC", c(1:max_pc_i))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}
pca_y <- pca_x %>% 
  tibble::rownames_to_column(., "Name")

# pair PCs against each other
pcs <- paste0("PC", c(1:max_pc_i))
pcs_l <- lapply(c(2:length(pcs)), function(i){c(pcs[[i-1]], pcs[[i]])})

# plot the pc"eh"
ggps <- lapply(pcs_l, function(pc){
  ggplot(data=pca_y, aes_string(x=pc[1], y=pc[2], color='condition', 
                                fill='condition', label='Name')) +
    geom_point() + 
    geom_text_repel() +
    theme_minimal() +
    scale_shape_manual(values=c(21,24)) +
    geom_line(aes(group = condition), alpha=0.3)
})
# pdf(file.path(pca_dir, "pca.pdf"), width = 8, height = 12)
pdf(file.path("~/xfer", "pca.pdf"), width = 8, height = 12)
cowplot::plot_grid(plotlist=ggps, ncol=2)
dev.off()
saveRDS(pca_x, file=file.path(pca_dir, "pca_x.rds"))




##############################################################################
#### 2.a) Differential expression controlling for interaction terms - Raw ####
dir.create(file.path(outdir, "degs"), showWarnings = F)
cntsdir <- file.path(pdir, "results", "counts")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(file.path(cntsdir, "all.tsv"), header=TRUE, 
                  row.names="gene", check.names=FALSE)

## MAIN TREATMENT EFFECT: BMDM vs apBMDM
# Relates to the treatment effect of CIS vs PBS for TDLN
# - By setting the reference level to TDLN, we can evaluate CIS vs PBS
# differences in terms of TDLN status
dds2 <- DESeqDataSetFromMatrix(countData=cts[,rownames(samples_meta)],
                               colData=samples_meta,
                               design=as.formula('~group'))

# remove uninformative columns
dds2 <- DESeq(dds2[ rowSums(counts(dds2)) > 5, ])
coef <- resultsNames(dds2)[2] # 'group_BMDM_vs_apBMDM'
res <- lfcShrink(dds2, coef=coef, type="apeglm")
res_df <- res %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var='ens') %>%
  mutate(gene=gm$ENSEMBL$SYMBOL[ens],
         biotype=ens2biotype_ids[ens],
         score=log2FoldChange * (-1*log10(padj))) %>%
  relocate(c(gene, biotype), .after=ens) %>%
  arrange(desc(score))
res_df %>% filter(padj < 0.05) %>% head


write.table(res_df, file=file.path(outdir, "degs", "degs.csv"), 
            col.names = T,row.names = F, quote=F, sep=",")
resl <- list("res"=res, "res_df"=res_df)
saveRDS(resl, file=file.path(outdir, "degs", "degs.rds"))


######################################
#### 3. GSEA analysis of the DEGs ####
dir.create(file.path(outdir, "gsea"), showWarnings = F)

res_dds <- readRDS(file.path(outdir, "degs", "degs.rds"))
min_baseMean <- 5

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF'))#, # ontology gene sets
# 'C7'=list('IMMUNESIGDB'))             # immunologic signature gene sets

lfc_df <- res_dds$res_df %>% 
  filter(baseMean > min_baseMean,
         biotype=='protein_coding') %>% 
  dplyr::select(grep("ens|log2FoldChange|score", colnames(.), value=T)) %>%
  tibble::column_to_rownames(., "ens")

# get GSEA table
gseas <- apply(lfc_df, 2, function(lfc_v){
  print("...")
  lfc_v <- setNames(lfc_v,
                    gm$ENSEMBL$ENTREZID[rownames(lfc_df)])
  
  iterateMsigdb(species='Mus musculus', msig_lvls=msig_lvls, 
                fun=gseaFun, lfc_v=lfc_v)
})

saveRDS(gseas, file=file.path(outdir, "gsea", "gseas.rds"))
gseas <- readRDS(file.path(outdir, "gsea", "gseas.rds"))
gseas <- lapply(gseas, function(i) unlist(i, recursive = F))
gseas_df <- lapply(gseas, function(gsea_i){
  do.call(rbind, lapply(gsea_i, as.data.frame)) %>% 
    select(ID, setSize, enrichmentScore, NES, pvalue, qvalues, leading_edge) %>%
    tibble::rownames_to_column('group') %>%
    mutate(leading_edge=gsub("[,%]", "", leading_edge) %>%
             gsub(" ", "_", .),
           group=strsplit(group, split="\\.") %>% 
             sapply(., function(i) paste0(i[1], ".", i[2])))# %>%
  # filter(qvalues < 0.2) %>% 
  # arrange(qvalues)
})

for(id in names(gseas_df)){
  write.table(gseas_df[[id]], file=file.path(outdir, "gsea", paste0("gsea.", id, ".csv")),
              sep=",", col.names = T, row.names = F, quote = F)
}

