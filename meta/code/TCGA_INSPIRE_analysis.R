renv::load("/cluster/projects/mcgahalab/envs/renvs/seuratv5_v2")
### Analyze Xena-processed TCGA/GTEx/TARGET data using a custom treg signature
library(tidyverse)
library(ggplot2)
library(GSVA)
library(survival)
library(ggmap)
library(immunedeconv)


source("~/git/mini_projects/mini_functions/TCGAanalyze_SurvivalKM2.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
gm <- geneMap()

# INS.PDIR='/cluster/projects/mcgahalab/data/mcgahalab/INSPIRE'
INS.PDIR='/cluster/projects/mcgahalab/ext_data/inspire/Source Data'
INS.DATADIR=file.path(INS.PDIR, "SourceData_Fig4")
GTEX.PDIR='/cluster/projects/mcgahalab/ext_data/gtex'
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcga_treg_sig'
DATADIR <- file.path(PDIR, "data")
OUTDIR <- file.path(PDIR, "results")
gtex_f <- file.path(GTEX.PDIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
gtex_rds_f <- file.path(GTEX.PDIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.rds")
inscounts_f <- file.path(INS.DATADIR, "gene-expression-matrix-TPM-final.tsv")
inscounts_rds_f <- file.path(DATADIR, "rna", "inspire.rds")
normcounts_f <- file.path(DATADIR, "rna", "TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz")
tpm_f <- file.path(DATADIR, "rna", "TcgaTargetGtex_rsem_gene_tpm.gz")
probemap_f <- file.path(DATADIR, "rna", "probeMap%2Fgencode.v23.annotation.gene.probemap")
survival_f <- file.path(DATADIR, "phenotype", "TCGA_survival_data.txt")
tissue_f <- file.path(DATADIR, "phenotype", "TCGA_GTEX_category.txt")
phenotype_f <- file.path(DATADIR, "phenotype", "TcgaTargetGTEX_phenotype.txt.gz")
signature_f <- file.path(PDIR, "ref", 'treg_signature/tcga_treg_signature.csv')

if(file.exists(gsub(".gz$", ".rds", normcounts_f))){
  message("Reading in rds files...")
  expr <- readRDS(file=gsub(".gz$", ".rds", normcounts_f))
  tpm <- readRDS(file=gsub(".gz$", ".rds", tpm_f))
  ins.expr <- readRDS(file=inscounts_rds_f)
  gtex.tpm <- readRDS(file=gtex_rds_f)
} else {
  message("Reading in txt file...")
  expr <- read.table(normcounts_f, sep="\t", header = T)
  tpm <- read.table(tpm_f, sep="\t", header = T)
  ins.expr <- read.table(inscounts_f, sep="\t", header=T)
  gtex.tpm <- read.table(gtex_f, sep="\t", header=T, skip = 2)
  biotypes <- gm$ENSEMBL$gene_biotype[gsub("\\.[0-9]*$", "", gtex.tpm$Name)]
  idx <- which(biotypes == 'protein_coding' &
                 !duplicated(gm$ENSEMBL$SYMBOL[gsub("\\.[0-9]*$", "", gtex.tpm$Name)]))
  gtex.tpm <- gtex.tpm[idx,] %>%
    tibble::remove_rownames() %>%
    mutate(symbol=gm$ENSEMBL$SYMBOL[gsub("\\.[0-9]*$", "", gtex.tpm$Name)]) %>%
    tibble::column_to_rownames(., "symbol") %>%
    dplyr::select(-c(Name, Description))
  # gtex.expr <- read.table(gsub("tpm", "reads", gtex_f), sep="\t", header=T, skip = 2)
  # gtex.expr <- gtex.expr[idx,] %>%
  #   tibble::remove_rownames() %>%
  #   mutate(symbol=gm$ENSEMBL$SYMBOL[gsub("\\.[0-9]*$", "", gtex.expr$Name)]) %>%
  #   tibble::column_to_rownames(., "symbol") %>%
  #   dplyr::select(-c(Name, Description))
  
  saveRDS(round(gtex.tpm, 3), file=gtex_rds_f)
  saveRDS(ins.expr, file=inscounts_rds_f)
  saveRDS(expr, file=gsub(".gz$", ".rds", normcounts_f))
  saveRDS(tpm, file=gsub(".gz$", ".rds", tpm_f))
}
probemap <- read.table(probemap_f, sep="\t", header = T)
probemap_filt <- probemap %>% 
  mutate(biotype =  gm$ENSEMBL$gene_biotype[gsub("\\.[0-9]$", "", id)]) %>%
  dplyr::filter(biotype == 'protein_coding' | is.na(biotype))  %>%
  dplyr::filter(!duplicated(gene))
survival <- read.table(survival_f, sep="\t", header = T)
phenotype <- read.table(phenotype_f, sep="\t", header = T)
signature <- read.table(signature_f, header=F, col.names = c('direction', 'gene'))
signature <- rbind(signature, 
                   read.table(gsub("signature.csv", "benchmark_signature.csv", signature_f), 
                              header=F, col.names = c('direction', 'gene'))) %>%
  mutate(gene=toupper(gene))
# signature[which(!signature$gene %in% rownames(expr_mat)),]

.cleanExpr <- function(expr, probemap, samplecol='X',
                       id1='gene', id2='id'){
  probemap_v <- setNames(probemap[,id1], probemap[,id2])
  if(samplecol=='gene') names(probemap_v) <- gsub("\\.[0-9]*$", "", names(probemap_v))
  
  expr_mat <- as.data.frame(expr) %>%
    tibble::remove_rownames() %>% 
    dplyr::filter(!!rlang::sym(samplecol) %in% names(probemap_v)) 
  expr_mat[,samplecol] <- probemap_v[expr_mat[,samplecol]]
  expr_mat <- expr_mat %>%
    tibble::column_to_rownames(samplecol) %>% 
    dplyr::rename_with(., ~gsub("\\.", "-", .)) %>% 
    as.matrix
  return(expr_mat)
}
expr_mat <- .cleanExpr(expr, probemap_filt)
ins.expr_mat <- .cleanExpr(ins.expr %>% tibble::rownames_to_column("gene"), 
                           probemap_filt, samplecol='gene', id2='gene')
# gtex.expr_mat <- .cleanExpr(gtex.expr %>% tibble::rownames_to_column("gene"), 
#                            probemap_filt, samplecol='gene', id2='gene')
gtex.tpm_mat <- .cleanExpr(gtex.tpm %>% tibble::rownames_to_column("gene"), 
                            probemap_filt, samplecol='gene', id2='gene')
tpm_mat <- .cleanExpr(tpm, probemap_filt, samplecol='sample')
signature <- with(signature, split(gene, f=direction))
ids <- split(names(signature), f=gsub("_(up|down)$", "", names(signature),ignore.case = T))
signature2 <- lapply(ids[sapply(ids, length)>1], function(id_i){
  unique(as.character(unlist(signature[id_i])))
})
signature <- c(signature, signature2)

##################################
#### Read in GTEx metadata ####
gtex.meta <- read.table(file.path(GTEX.PDIR, "gtex_meta.txt"),
                        sep="\t", header=T, check.names = F, stringsAsFactors = F, comment.char = "")

gtex_meta <- gtex.meta %>% 
  dplyr::filter(SAMPID %in% colnames(gtex.tpm_mat)) %>%
  dplyr::filter(SMTS %in% 'Skin')
gtex_tpm <- gtex.tpm_mat[,gtex_meta$SAMPID]

##################################
#### Read in INSPIRE metadata ####
ins.meta <- lapply(list.files(file.path(INS.PDIR, "csv")), function(f){
  read.csv(file.path(INS.PDIR, "csv", f))
}) %>% 
  purrr::reduce(., full_join, by='PATIENT_ID') %>% 
  dplyr::select(-grep("\\.[xy]$", colnames(.), value=T))

inspire_expr <- as.data.frame(ins.expr_mat) %>% 
  dplyr::select(., grep("-ST", colnames(.), value=T)) %>%
  magrittr::set_colnames(., gsub("-ST$","", colnames(.)))

inspire_meta <- ins.meta %>% 
  dplyr::filter(PATIENT_ID %in% colnames(inspire_expr)) %>%
  dplyr::rename_with(., ~gsub("event$", "", .)) %>%
  dplyr::rename_with(., ~gsub("TIME_Months_since_C3$", ".time", .)) %>%
  dplyr::rename_with(., ~gsub("PFS", "PFI", .)) %>%
  dplyr::rename_with(., ~gsub("PATIENT_ID", "sample", .))
  
inspire_expr <- as.matrix(inspire_expr[,inspire_meta$sample])
  
################################
#### TCGA specific analysis ####
tcga_expr <- expr_mat[,grep("TCGA", colnames(expr_mat))]
tcga_tpm <- tpm_mat[,grep("TCGA", colnames(tpm_mat))]
tcga_phenotype <- phenotype %>% dplyr::filter(X_study == 'TCGA')

## Load extended TCGA clinical ##
tcgadir <- '/cluster/projects/mcgahalab/ref/TCGA/data/clinical'
tcgaclin <- lapply(list.files(tcgadir, pattern=".rds$"), function(f){
  readRDS(file.path(tcgadir, f))
}) %>% plyr::rbind.fill(.)

tcga_phenotype_filt  <- tcga_phenotype %>% 
  dplyr::filter(!X_sample_type %in% c('Control analyte', 'Additional Metastatic', 
                                      'Solid Tissue Normal', 'Primary Blood Derived Cancer - Peripheral Blood')) %>%
  dplyr::filter(!X_primary_site %in% c('White blood cell'))
tcga_phenotype_filt  <- tcga_phenotype %>% 
  dplyr::filter(!X_sample_type %in% c('Solid Tissue Normal'))

tcga_phenotype_trx <- tcga_phenotype_filt %>%
  mutate(submitter_id=gsub("-[0-9]*$", "", sample)) %>%
  left_join(., tcgaclin[,c('submitter_id', 'treatments_radiation_treatment_or_therapy', 'treatments_pharmaceutical_treatment_or_therapy')])
tcga_phenotype_trx$treatment <- (tcga_phenotype_trx$treatments_radiation_treatment_or_therapy == 'yes') + 
  (tcga_phenotype_trx$treatments_pharmaceutical_treatment_or_therapy == 'yes')
id1 <- tcga_phenotype_trx$treatments_radiation_treatment_or_therapy %>%
  ifelse(.=='yes', 2, .) %>% 
  ifelse(.=='no', 1, .) %>% 
  ifelse(.=='not reported', 0, .)
id2 <- tcga_phenotype_trx$treatments_pharmaceutical_treatment_or_therapy %>%
  ifelse(.=='yes', 2, .) %>% 
  ifelse(.=='no', 1, .) %>% 
  ifelse(.=='not reported', 0, .)
tcga_phenotype_trx$treatment_code <- paste0(id1, id2)

tcga_phenotype_negtrx  <- tcga_phenotype_trx %>% 
  dplyr::filter(treatment == 0)
tcga_phenotype_postrx  <- tcga_phenotype_trx %>% 
  dplyr::filter(treatment > 0)

#########################################
#### Extract TCGA Tumor-normal pairs ####
normal_tcga <- tcga_phenotype %>% 
  dplyr::filter(X_sample_type %in% c('Solid Tissue Normal'))
tn_idx <- sapply(gsub("-[0-9]*$", "", normal_tcga$sample), function(i) grep(i, tcga_phenotype$sample))
tn_idx <- tn_idx[which(sapply(tn_idx, length) > 1)]
tcga_phenotype_tn <- lapply(tn_idx, function(i) tcga_phenotype[i,])


######################
#### ssGSEA main  ####
ssgseaparam <- GSVA::ssgseaParam(tcga_tpm, signature)
ssgsea_dat <- GSVA::gsva(ssgseaparam)

ssgseaparam <- GSVA::ssgseaParam(inspire_expr, signature)
ssgsea_dat_inspire <- GSVA::gsva(ssgseaparam)

ssgseaparam <- GSVA::ssgseaParam(gtex_tpm, signature)
ssgsea_dat_gtex <- GSVA::gsva(ssgseaparam)

# gsvaparam <- GSVA::gsvaParam(tcga_expr, signature)
# gsva_dat <- gsva(gsvaparam)

# x <- split(as.data.frame(t(ssgsea_dat)), 
#            f=with(tcga_phenotype, setNames(primary.disease.or.tissue, sample))[colnames(ssgsea_dat)])
# xb <- sapply(x, function(i){colMeans(i)})
# lapply(setNames(rownames(xb),rownames(xb)), function(i) t(t(xb[i, order(xb[i,])])))
# 




pdf("~/xfer/tcga_inspire_all.pdf", width = 11)
metric <- 'PFI'
for(geneseti in rownames(ssgsea_dat)){
  quant_spl <- split(gsub("\\.", "-", colnames(ssgsea_dat)), 
                     f=cut(ssgsea_dat[geneseti,], breaks=quantile(ssgsea_dat[geneseti,], c(0, 0.25, 0.5, 0.75, 1))))
  x=XenaTCGAanalyze_SurvivalKM(clinical_patient=survival, data_grp=quant_spl,
                               Survresult = T, metric=metric, caption=paste0("All_TCGA [", metric, "]: ", geneseti))
  print(x)
}

for(geneseti in rownames(ssgsea_dat)){
  quant_spl <- split(colnames(ssgsea_dat_inspire), 
                     f=cut(ssgsea_dat_inspire[geneseti,], 
                           breaks=quantile(ssgsea_dat_inspire[geneseti,], c(0, 0.25, 0.5, 0.75, 1))))
  y=XenaTCGAanalyze_SurvivalKM(clinical_patient=inspire_meta, data_grp=quant_spl,
                               Survresult = T, metric=metric, 
                               caption=paste0("All_INSPIRE [", metric, "]: ", geneseti),
                               skipvalidation=T)
  print(y)
}
dev.off()

########################################
#### INSPIRE: ssGSEA per cancerType ####
phenospl <- with(inspire_meta, split(sample, f=COHORT))

Survresult <- T
breakid <- 'median' # 'median', 'Q14'
breaks <- if(breakid == 'median') c(0, 0.5, 1) else c(0, 0.25, 0.75, 1)
nbins <- 2
metric <- 'PFI'
ggs <- lapply(names(phenospl)[], function(phenoid){
  phenospl_i <- phenospl[[phenoid]]
  phenospl_i <- phenospl_i[(phenospl_i %in% colnames(inspire_expr))]
  gsvapar_i <- GSVA::ssgseaParam(inspire_expr[,phenospl_i], signature)
  ssgsea_dat_i <- gsva(gsvapar_i)
  
  
  gg_scs <- lapply(rownames(ssgsea_dat_i), function(geneseti){
    bins <- cut(ssgsea_dat_i[geneseti,], breaks=round(quantile(ssgsea_dat_i[geneseti,], breaks),2), 
                include.lowest=T)
    span <- ceiling(nbins/2)
    bin.n <- length(levels(bins))
    low.idx <- 1:span
    high.idx <- (bin.n-span+1):bin.n
    levels(bins)[setdiff(c(1:bin.n), c(low.idx, high.idx))] <- NA
    
    quant_spl <- split(gsub("\\.", "-", colnames(ssgsea_dat_i)), 
                       f=bins) # f=cut(ssgsea_dat_i[geneseti,], breaks=quantile(ssgsea_dat_i[geneseti,], c(0, 0.25, 0.5, 0.75, 1))))
    ggsc=tryCatch({
      XenaTCGAanalyze_SurvivalKM(clinical_patient=inspire_meta, data_grp=quant_spl,
                                 Survresult = T, metric=metric, add.pvaltbl=F, ret.pvaltbl=T,
                                 caption=paste0("INSPIRE ", phenoid, " [", metric, "]: ", geneseti),
                                 skipvalidation=T)
    }, error=function(e){NULL})
    return(ggsc)
  })
  
  
  return(gg_scs)
})
pdf(paste0("~/xfer/inspire_cancer_specific.", metric, ".", breakid, ".pdf"), width=12)
ggs
dev.off()
cat(paste0("xfer inspire_cancer_specific.", metric, ".", breakid, ".pdf\n"))

#####################################
#### TCGA: ssGSEA per cancerType ####
phenospl <- with(tcga_phenotype_filt, split(sample, f=primary.disease.or.tissue))
phenospl <- phenospl[grep("(Breast|Melanoma)", names(phenospl), value=T)]

Survresult <- T
breakid <- 'median' # 'median', 'Q14'
breaks <- if(breakid == 'median') c(0, 0.5, 1) else c(0, 0.25, 0.75, 1)
nbins <- 2
metric <- 'OS'
ggs <- lapply(names(phenospl)[], function(phenoid){
  message(phenoid)
  phenospl_i <- phenospl[[phenoid]]
  phenospl_i <- phenospl_i[(phenospl_i %in% colnames(tcga_expr))]
  gsvapar_i <- GSVA::ssgseaParam(tcga_expr[,phenospl_i], signature)
  ssgsea_dat_i <- gsva(gsvapar_i)
  
  
  
  
  gg_scs <- lapply(rownames(ssgsea_dat_i), function(geneseti){
    bins <- cut(ssgsea_dat_i[geneseti,], breaks=round(quantile(ssgsea_dat_i[geneseti,], breaks),2), 
                include.lowest=T)
    span <- ceiling(nbins/2)
    bin.n <- length(levels(bins))
    low.idx <- 1:span
    high.idx <- (bin.n-span+1):bin.n
    levels(bins)[setdiff(c(1:bin.n), c(low.idx, high.idx))] <- NA
    
    quant_spl <- split(gsub("\\.", "-", colnames(ssgsea_dat_i)), 
                       f=bins) # f=cut(ssgsea_dat_i[geneseti,], breaks=quantile(ssgsea_dat_i[geneseti,], c(0, 0.25, 0.5, 0.75, 1))))
    ggsc=XenaTCGAanalyze_SurvivalKM(clinical_patient=survival, data_grp=quant_spl,
                                 Survresult = T, metric=metric, add.pvaltbl=F, ret.pvaltbl=T,
                                 caption=paste0("TCGA ", phenoid, " [", metric, "]: ", geneseti))
    return(ggsc)
  })
  
  
  return(gg_scs)
})
pdf(paste0("~/xfer/tcga_cancer_specific.", metric, ".", breakid, ".pdf"), width=12)
ggs
dev.off()
cat(paste0("xfer tcga_cancer_specific.", metric, ".", breakid, ".pdf\n"))


#############################################################
#### TCGA: ssGSEA per cancerType - with treatment status ####
phenospl <- split(tcga_phenotype_trx,tcga_phenotype_trx$primary.disease.or.tissue)
phenospl <- phenospl[grep("(Breast|Melanoma)", names(phenospl), value=T)]


metric <- 'OS'
breakid <- 'Q14' # 'median', 'Q14'
Survresult <- T
breaks <- if(breakid == 'median') c(0, 0.5, 1) else c(0, 0.25, 0.75, 1)
nbins <- 2
ggs <- lapply(names(phenospl)[], function(phenoid){
  # ggs <- lapply(names(phenospl)['Breast Invasive Carcinoma'], function(phenoid){
  message(phenoid)
  phenospl_i <- phenospl[[phenoid]] %>%
    dplyr::filter(sample %in% colnames(tcga_expr))
  gsvapar_i <- GSVA::ssgseaParam(tcga_expr[,phenospl_i$sample], signature)
  ssgsea_dat_i <- gsva(gsvapar_i)
  
  gg_scs <- lapply(rownames(ssgsea_dat_i), function(geneseti){
    bins <- cut(ssgsea_dat_i[geneseti,], breaks=round(quantile(ssgsea_dat_i[geneseti,], breaks),2), 
                include.lowest=T)
    span <- ceiling(nbins/2)
    bin.n <- length(levels(bins))
    low.idx <- 1:span
    high.idx <- (bin.n-span+1):bin.n
    levels(bins)[setdiff(c(1:bin.n), c(low.idx, high.idx))] <- NA
    
    quant_spl <- split(gsub("\\.", "-", phenospl_i$sample), 
                       list(phenospl_i$treatment_code, bins))
    ggsc=tryCatch({
      XenaTCGAanalyze_SurvivalKM(clinical_patient=survival, data_grp=quant_spl,
                                   Survresult = Survresult, metric=metric, add.pvaltbl=F,
                                 caption=paste0("TCGA ", phenoid, " [", metric, "]: ", geneseti))
    }, error=function(e){NULL})
    return(ggsc)
  })
  
  
  return(gg_scs)
})

pdf(paste0("~/xfer/tcga_cancer_specific.treatment.", metric, ".", breakid, ".pdf"), height = 10, width=15)
ggs
dev.off()
cat(paste0("xfer tcga_cancer_specific.treatment.", metric, ".", breakid, ".pdf\n"))



ggs2 <- lapply(ggs2, function(i){
  names(i) <- names(signature)
  return(i)
})
names(ggs2) <- names(phenospl)[-1]
pval_df <- lapply(names(ggs2), function(ctid) lapply(names(ggs2[[ctid]]), function(genesetid) {
  reshape2::melt(ggs2[[ctid]][[genesetid]]$pval) %>% 
    magrittr::set_colnames(., c("Start", "End", "p.value")) %>%
    mutate("Celltype"=ctid,
           "Geneset"=genesetid)
}) %>% do.call(rbind, .)) %>% do.call(rbind, .)
write.table(pval_df, file="~/xfer/tcga_os.csv", 
            quote = F, sep="\t", col.names = T, row.names = F)

###########################################################
#### TCGA: Comparing ssGSEA between Tumor-Normal pairs ####
tn.ssgsea <- sapply(tcga_phenotype_tn, function(phenotype_i){
  print(phenotype_i$sample[1])
  if(!all(phenotype_i$sample %in% colnames(ssgsea_dat))) return(NULL)
  df <- as.data.frame(t(ssgsea_dat[,phenotype_i$sample])) %>%
    mutate(sample_type=phenotype_i$X_sample_type)
  
  
  cols <- c(1:nrow(ssgsea_dat))
  tid <- if(any(grepl("Tumor", df$sample_type))){'Tumor'} else {'Metastatic'}
  dfsumm <- (df[grep(tid, df$sample_type),cols] - df[grep("Normal", df$sample_type),cols]) %>%
    mutate(disease=phenotype_i$primary.disease.or.tissue[1],
           sample=gsub("-[0-9]*$", "", phenotype_i$sample[1]))  %>%
    relocate(., c(sample, disease))
  return(dfsumm)
})
tn.ssgsea.df <- tn.ssgsea %>% 
  do.call(rbind, .) %>%
  tidyr::pivot_longer(., cols=!c('sample', 'disease'))
pdf("~/xfer/tcga_tn.pdf", width = 15, height = 10)
ggplot(tn.ssgsea.df, aes(x=disease, y=value, fill=disease)) +
  facet_grid(name~., space = 'free') +
  geom_violin() +
  ylim(-0.2, 0.2) +
  # ggbeeswarm::geom_beeswarm() + 
  cowplot::theme_cowplot() + 
  ylab("Delta [ Tumor - Normal ]") +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
  geom_hline(yintercept = 0, linetype='dotted', col='grey')
dev.off()
########################################################
#### TCGA: Comparing ssGSEA-scores to tcga metadata ####
categories <- c('ajcc_pathologic_stage', 'tumor_stage', 'ajcc_clinical_m', 
                'ajcc_pathologic_n', 
                'treatments_pharmaceutical_treatment_or_therapy',
                'treatments_radiation_treatment_or_therapy') %>%
  setNames(., .)
genesetids <- rownames(ssgsea_dat) %>%
  setNames(., .)
ssgsea_df <- as.data.frame(t(ssgsea_dat)) %>% 
  tibble::rownames_to_column('sample') %>%
  mutate(submitter_id = gsub("^(.*)-([0-9]*)$", "\\1", sample))
ssgsea_clin_df <- ssgsea_df %>% 
  left_join(., tcgaclin[,c('submitter_id', 'disease', categories)], by='submitter_id')
diseases <- unique(ssgsea_clin_df$disease)

aov_res <- lapply(diseases, function(d.k){
  ssgsea_clin_df.k <- split(ssgsea_clin_df, f=ssgsea_clin_df$disease)[[d.k]]
  lapply(genesetids, function(gs.i){
    lapply(categories, function(cat.j){
      aovres <- tryCatch({
        aov(as.formula(paste0(gs.i,  " ~ ", cat.j)), data = ssgsea_clin_df.k)
      }, error=function(e){NULL})
      if(is.null(aovres)) return(aovres)
      
      raw.dat <- split(ssgsea_clin_df.k[,gs.i], ssgsea_clin_df.k[,cat.j]) 
      tukeydf <- as.data.frame(TukeyHSD(aovres)[[1]]) %>% 
        rename_with(., ~gsub("p adj", "p.adj", .)) %>% 
        tibble::rownames_to_column(., "comparisons") %>%
        arrange(`p.adj`) %>%
        mutate(grp1=gsub("-.*?$", "", comparisons),
               grp2=gsub("^.*?-", "", comparisons)) %>%
        mutate(grp1.mean = sapply(raw.dat, mean)[grp1],
               grp2.mean = sapply(raw.dat, mean)[grp2],
               grp1.n = sapply(raw.dat, length)[grp1],
               grp2.n = sapply(raw.dat, length)[grp2])

      aov.pval <- data.frame("comparisons"="ANOVA", "p.adj"=summary(aovres)[[1]]$`Pr(>F)`[1])
      plyr::rbind.fill(aov.pval, tukeydf) %>% 
        mutate(Geneset=gs.i,
               Category=cat.j,
               Disease=d.k)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

aov_res <- aov_res %>%
  tibble::remove_rownames()
apply(aov_res, 2, function(i) grepl(",", i)) %>% colSums
write.table(aov_res, file=file.path("~/xfer/tcga_aov_res.csv"),
            sep=",", col.names = T, row.names = F, quote = F)

#### INSPIRE: Comparin
###########################################################
#### INSPIRE: Comparing ssGSEA-scores to clinical metadata ####
# ANOVA and Tukeys HSD for categorical data
categories <- c('BEST_OVERALL_RECIST') %>%
  setNames(., .)
genesetids <- rownames(ssgsea_dat_inspire) %>%
  setNames(., .)
ssgsea_df <- as.data.frame(t(ssgsea_dat_inspire)) %>% 
  tibble::rownames_to_column('sample') 
ssgsea_clin_df <- ssgsea_df %>% 
  left_join(., inspire_meta, by='sample')
diseases <- unique(ssgsea_clin_df$COHORT)

aov_res <- lapply(diseases, function(d.k){
  ssgsea_clin_df.k <- split(ssgsea_clin_df, f=ssgsea_clin_df$COHORT)[[d.k]]
  lapply(genesetids, function(gs.i){
    lapply(categories, function(cat.j){
      aovres <- tryCatch({
        aov(as.formula(paste0(gs.i,  " ~ ", cat.j)), data = ssgsea_clin_df.k)
      }, error=function(e){NULL})
      if(is.null(aovres)) return(aovres)
      
      raw.dat <- split(ssgsea_clin_df.k[,gs.i], ssgsea_clin_df.k[,cat.j]) 
      tukeydf <- as.data.frame(TukeyHSD(aovres)[[1]]) %>% 
        rename_with(., ~gsub("p adj", "p.adj", .)) %>% 
        tibble::rownames_to_column(., "comparisons") %>%
        arrange(`p.adj`) %>%
        mutate(grp1=gsub("-.*?$", "", comparisons),
               grp2=gsub("^.*?-", "", comparisons)) %>%
        mutate(grp1.mean = sapply(raw.dat, mean)[grp1],
               grp2.mean = sapply(raw.dat, mean)[grp2],
               grp1.n = sapply(raw.dat, length)[grp1],
               grp2.n = sapply(raw.dat, length)[grp2])
      low.n.idx <- with(tukeydf, (grp1.n < 3)  | (grp2.n < 3))
      if(any(low.n.idx)) tukeydf$p.adj[which(low.n.idx)] <- NA
      
      aov.pval <- data.frame("comparisons"="ANOVA", "p.adj"=summary(aovres)[[1]]$`Pr(>F)`[1])
      plyr::rbind.fill(aov.pval, tukeydf) %>% 
        mutate(Geneset=gs.i,
               Category=cat.j,
               Disease=d.k)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

aov_res <- aov_res %>%
  tibble::remove_rownames()
apply(aov_res, 2, function(i) grepl(",", i)) %>% colSums
write.table(aov_res, file=file.path("~/xfer/inspire_aov_res.csv"),
            sep=",", col.names = T, row.names = F, quote = F)




categories2 <- c('BEST_RECIST_TM_TARGET', 'TARGETED_LESION_TM_EARLY', 
                'change_ctDNA', 'TMB_n72') %>%
  setNames(., .)

cor_res <- lapply(diseases, function(d.k){
  ssgsea_clin_df.k <- split(ssgsea_clin_df, f=ssgsea_clin_df$COHORT)[[d.k]]
  lapply(genesetids, function(gs.i){
    lapply(categories2, function(cat.j){
      x <- ssgsea_clin_df.k[,gs.i]
      y <- as.numeric(ssgsea_clin_df.k[,cat.j])
      cor.res <- cor.test(x, y, use='complete.obs', method='spearman')
      data.frame("r"=cor.res$estimate, "p"=cor.res$p.value,
                 "category"=cat.j, "disease"=d.k, "geneset"=gs.i,
                 "n"=min(c(sum(!is.na(x)), sum(!is.na(y)))))
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  as.data.frame %>% 
  mutate(padj=p.adjust(p))

write.table(cor_res, file=file.path("~/xfer/inspire_cor_res.csv"),
            sep=",", col.names = T, row.names = F, quote = F)

###################
#### CIBERSORT ####
ggplotit <- function(res){
  cols <- rainbow(n=nrow(res), s = 0.5, v=0.6)
  cols <- sample(cols, length(cols), replace = F)
  
  res %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_manual(values=cols) +
    scale_x_discrete(limits = rev(levels(res)))
}

# CIBERSORT
cibersort_path <- '/cluster/projects/mcgahalab/ref/immunedeconv/cibersort'
set_cibersort_binary(file.path(cibersort_path, "CIBERSORT.R"))
set_cibersort_mat(file.path(cibersort_path, "LM22.txt"))

#--- a) TCGA ----
# CIBERSORT per cancertype
phenospl <- with(tcga_phenotype, split(sample, f=primary.disease.or.tissue))

cibersort_dat <- lapply(names(phenospl)[-1], function(phenoid){
  outid <- gsub(" ", "_", phenoid) %>%
    gsub("&", "and", .) %>% gsub("-", "_", .)
  message(outid)
  
  dir.create(file.path(OUTDIR, "cibersort"), showWarnings = F)
  outf <- file.path(OUTDIR, "cibersort", paste0(outid, ".rds"))
  if(!file.exists(outf)){
    file.create(outf)
    phenospl_i <- phenospl[[phenoid]]
    phenospl_i <- phenospl_i[(phenospl_i %in% colnames(tcga_tpm))]
    phenospl_i <- setdiff(phenospl_i, 'TCGA-25-1870-01')
    cibersort <- deconvolute(tcga_tpm[,phenospl_i], "cibersort")
    saveRDS(cibersort, file=outf)
  } else {
    cibersort <- tryCatch({readRDS(outf)}, error=function(e){NULL})
  }
  return(cibersort)
})
names(cibersort_dat) <- names(phenospl)[-1]
cibersort_dat  <- cibersort_dat[!sapply(cibersort_dat, is.null)]
cibersort_mat <- as.data.frame(do.call(cbind, cibersort_dat)) %>% 
  dplyr::rename_with(., ~make.unique(gsub("^.*?\\.", "", .))) %>%
  tibble::column_to_rownames(., "cell_type") 
  

cibersort_cor <- lapply(names(phenospl)[-1], function(phenoid){
  phenospl_i <- phenospl[[phenoid]]
  phenospl_i <- phenospl_i[(phenospl_i %in% colnames(tcga_tpm))]
  
  apply(ssgsea_dat, 1, function(i){
    int_ids <- purrr::reduce(list(names(i), colnames(cibersort_mat), phenospl_i), .f=intersect)
    if(length(int_ids)==0) return(NULL)
    apply(cibersort_mat[,int_ids], 1, function(j){
      cor.res <- cor.test(i[int_ids], as.numeric(j))
      # ifelse(cor.res$p.value <= 0.05, cor.res$estimate, NA)
      return(cor.res)
    })
  })
})
names(cibersort_cor) <- names(phenospl)[-1]
cibersort_cor <- cibersort_cor[!sapply(cibersort_cor, is.null)]

genesets <- c('Sign_Tumor_V002_DOWN', 'Sign_Tumor_V002')
cancertypes <- c('Skin Cutaneous Melanoma', 'Breast Invasive Carcinoma')
cor_df <- lapply(cancertypes, function(cancertype_i) {
  lapply(genesets, function(geneset_j){
    print(paste0(cancertype_i, " - ", geneset_j))
    i <- cibersort_cor[[cancertype_i]]
    data.frame("geneset"=geneset_j, 
               "cancertype"=cancertype_i,
               "r"=sapply(i[[geneset_j]], function(x) x$estimate),
               "p"=sapply(i[[geneset_j]], function(x) x$p.value)) %>%
      mutate("log10padj"=-1*log10(p.adjust(p, method='BH')),
             "sig" = case_when(r >= 0.2 & log10padj >= (-1*log10(0.01)) ~ "up",
                               r <= -0.2 & log10padj >= (-1*log10(0.01)) ~ "down",
                               TRUE ~ "ns")) 
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
cor_df <- cor_df %>%
  tibble::rownames_to_column(., "celltype") %>%
  mutate(celltype = gsub("\\.cor.*", "", celltype))
sig_up <- cor_df %>%
  dplyr::filter(sig %in% 'up')
sig_down <- cor_df %>%
  dplyr::filter(sig %in% 'down')

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
final_plot <- ggplot(cor_df, aes(x=r, y=log10padj)) +
  geom_point(aes(colour = sig), 
             alpha = 0.7, 
             shape = 16,
             size = 3) + 
  geom_point(data = sig_up,
             shape = 21,
             size = 4, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = sig_down,
             shape = 21,
             size = 4, 
             fill = "steelblue", 
             colour = "black") + 
  facet_grid(geneset ~ cancertype, space='free') +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.2, 0.2),
             linetype = "dashed") +
  ggrepel::geom_label_repel(data = rbind(sig_up, sig_down),   
                   aes(label = celltype),
                   size=3,
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  xlim(-0.5, 0.5) +
  ylim(0, 40) +
  # labs(title = "Gene expression changes in diseased versus healthy samples",
  #      x = "log2(fold change)",
  #      y = "-log10(adjusted P-value)",
  #      colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

pdf("~/xfer/tcga_volcano_plots.pdf", width = 10, height = 10); final_plot; dev.off()



cor_df <- reshape2::melt(cibersort_cor) %>%
  magrittr::set_colnames(., c('Celltype', 'Geneset', 'corr', 'Disease'))
write.table(cor_df, file="~/xfer/cibersort_corr.csv",
            sep=",", quote = F, col.names = T, row.names = F)

#--- b) INSPIRE ----
phenospl <- with(inspire_meta, split(sample, f=COHORT))


cibersort_dat <- lapply(names(phenospl), function(phenoid){
  outid <- gsub(" ", "_", phenoid) %>%
    gsub("&", "and", .) %>% gsub("-", "_", .) %>%
    gsub("^.*?\\:_", "INSPIRE.", .)
  message(outid)
  
  dir.create(file.path(OUTDIR, "cibersort"), showWarnings = F)
  outf <- file.path(OUTDIR, "cibersort", paste0(outid, ".rds"))
  if(!file.exists(outf)){
    file.create(outf)
    phenospl_i <- phenospl[[phenoid]]
    phenospl_i <- phenospl_i[(phenospl_i %in% colnames(inspire_expr))]
    cibersort <- deconvolute(inspire_expr[,phenospl_i], "cibersort")
    saveRDS(cibersort, file=outf)
  } else {
    cibersort <- tryCatch({readRDS(outf)}, error=function(e){NULL})
  }
  return(cibersort)
})
cibersort_dat  <- cibersort_dat[!sapply(cibersort_dat, is.null)]
cibersort_mat <- as.data.frame(do.call(cbind, cibersort_dat)) %>% 
  tibble::column_to_rownames(., "cell_type") 


cibersort_cor <- lapply(names(phenospl), function(phenoid){
  phenospl_i <- phenospl[[phenoid]]
  phenospl_i <- phenospl_i[(phenospl_i %in% colnames(inspire_expr))]
  
  apply(ssgsea_dat_inspire, 1, function(i){
    int_ids <- purrr::reduce(list(names(i), colnames(cibersort_mat), phenospl_i), .f=intersect)
    if(length(int_ids)==0) return(NULL)
    apply(cibersort_mat[,int_ids], 1, function(j){
      cor.res <- cor.test(i[int_ids],j)
      ifelse(cor.res$p.value <= 0.05, cor.res$estimate, NA)
    })
  })
})
names(cibersort_cor) <- names(phenospl)
cibersort_cor <- cibersort_cor[!sapply(cibersort_cor, is.null)]

cor_df <- reshape2::melt(cibersort_cor) %>%
  magrittr::set_colnames(., c('Celltype', 'Geneset', 'corr', 'Disease'))
write.table(cor_df, file="~/xfer/cibersort_corr.inspire.csv",
            sep=",", quote = F, col.names = T, row.names = F)

#--- c) GTEx ----
phenospl <- with(gtex_meta, split(SAMPID, f=SMTSD))#SMTS))
x <- split(gtex_meta, f=gtex_meta$SMTSD)
sapply(x, nrow)
apply(x[[2]][,-1], 2, table)

cibersort_dat <- lapply(names(phenospl)[-1], function(phenoid){
  outid <- gsub(" ", "_", phenoid) %>%
    gsub("&", "and", .) %>% gsub("-", "_", .) %>%
    gsub("_+", "_", .)
  message(outid)
  
  dir.create(file.path(OUTDIR, "cibersort"), showWarnings = F)
  outf <- file.path(OUTDIR, "cibersort", paste0(outid, ".rds"))
  if(!file.exists(outf)){
    file.create(outf)
    phenospl_i <- phenospl[[phenoid]]
    phenospl_i <- phenospl_i[(phenospl_i %in% colnames(gtex_tpm))]
    cibersort <- deconvolute(gtex_tpm[,phenospl_i], "cibersort")
    saveRDS(cibersort, file=outf)
  } else {
    cibersort <- tryCatch({readRDS(outf)}, error=function(e){NULL})
  }
  return(cibersort)
})
cibersort_dat  <- cibersort_dat[!sapply(cibersort_dat, is.null)]
cibersort_mat <- as.data.frame(do.call(cbind, cibersort_dat)) %>% 
  tibble::column_to_rownames(., "cell_type") 


cibersort_cor <- lapply(names(phenospl), function(phenoid){
  phenospl_i <- phenospl[[phenoid]]
  phenospl_i <- phenospl_i[(phenospl_i %in% colnames(gtex_tpm))]
  
  apply(ssgsea_dat_gtex, 1, function(i){
    int_ids <- purrr::reduce(list(names(i), colnames(cibersort_mat), phenospl_i), .f=intersect)
    if(length(int_ids)==0) return(NULL)
    apply(cibersort_mat[,int_ids], 1, function(j){
      cor.res <- cor.test(i[int_ids], as.numeric(j))
      # ifelse(cor.res$p.value <= 0.05, cor.res$estimate, NA)
      return(cor.res)
    })
  })
})
names(cibersort_cor) <- names(phenospl)
cibersort_cor <- cibersort_cor[!sapply(cibersort_cor, is.null)]

genesets <- c('Sign_Tumor_V002_DOWN', 'Sign_Tumor_V002')
cancertypes <- names(cibersort_cor)
cor_df <- lapply(cancertypes, function(cancertype_i) {
  lapply(genesets, function(geneset_j){
    i <- cibersort_cor[[cancertype_i]]
    data.frame("geneset"=geneset_j, 
               "cancertype"=cancertype_i,
               "r"=sapply(i[[geneset_j]], function(x) x$estimate),
               "p"=sapply(i[[geneset_j]], function(x) x$p.value)) %>%
      mutate("log10padj"=-1*log10(p.adjust(p, method='BH')),
             "sig" = case_when(r >= 0.2 & log10padj >= (-1*log10(0.01)) ~ "up",
                               r <= -0.2 & log10padj >= (-1*log10(0.01)) ~ "down",
                               TRUE ~ "ns")) 
    
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
cor_df <- cor_df %>%
  tibble::rownames_to_column(., "celltype") %>%
  mutate(celltype = gsub("\\.cor.*", "", celltype))
sig_up <- cor_df %>%
  dplyr::filter(sig %in% 'up')
sig_down <- cor_df %>%
  dplyr::filter(sig %in% 'down')

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
final_plot <- ggplot(cor_df, aes(x=r, y=log10padj)) +
  geom_point(aes(colour = sig), 
             alpha = 0.7, 
             shape = 16,
             size = 3) + 
  geom_point(data = sig_up,
             shape = 21,
             size = 4, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = sig_down,
             shape = 21,
             size = 4, 
             fill = "steelblue", 
             colour = "black") + 
  facet_grid(geneset ~ cancertype, space='free') +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.2, 0.2),
             linetype = "dashed") +
  ggrepel::geom_label_repel(data = rbind(sig_up, sig_down),   
                            aes(label = celltype),
                            size=3,
                            force = 2,
                            nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  xlim(-0.5, 0.5) +
  ylim(0,40) +
  # labs(title = "Gene expression changes in diseased versus healthy samples",
  #      x = "log2(fold change)",
  #      y = "-log10(adjusted P-value)",
  #      colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

pdf("~/xfer/gtex_volcano_plots.pdf", width = 10, height = 10); final_plot; dev.off()



cor_df <- reshape2::melt(cibersort_cor) %>%
  magrittr::set_colnames(., c('Celltype', 'Geneset', 'corr', 'Disease'))
write.table(cor_df, file="~/xfer/cibersort_corr.csv",
            sep=",", quote = F, col.names = T, row.names = F)

