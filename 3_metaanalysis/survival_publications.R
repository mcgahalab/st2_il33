# renv::load("/cluster/home/quever/downloads/renvs/")
# renv::load("/cluster/projects/mcgahalab/envs/renvs/seuratv5_v2")
library(Seurat)
library(ggplot2)
library(cowplot)
library(harmony)
library(scuttle)
library(tidyverse)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e9)


library(GSVA)
library(survival)

PDIR <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/meta_analysis/treg_signature"
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
proj <- 'ipmn_pdac'
setwd(PDIR)

response_score <- list()


###################
#### Functions ####
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/as.SingleCellExperiment.Seurat5.R")
source("~/git/mini_projects/mini_functions/makeLoupe.R")
source("~/git/mini_projects/mini_functions/TCGAanalyze_SurvivalKM2.R")
gm <- geneMap(species='Homo sapiens')

##-- GSE functions ----
aucellFun <- function(msig_ds, expr_mat, gm, mapfrom='SYMBOL', mapto='ENTREZID'){
  mapped_id <- if(mapfrom==mapto){
    msig_ds$entrez_gene
  } else {
    gm[[mapfrom]][[mapto]][as.character(msig_ds$entrez_gene)]
  }
  
  msig_l <- lapply(split(mapped_id, msig_ds$gs_name), function(i){
    i[!is.na(i)]
  })
  auc <- tryCatch({
    AUCell::AUCell_run(expr_mat, msig_l)#, ...)
  }, error=function(e){NULL})
  return(auc)
}

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
###################
#### 0. Prelim ####
treg_sig <- read.csv(file.path(PDIR, "ref", "treg_signature_paper.human.csv"), header = F) %>% 
  magrittr::set_colnames(., c("gs_name", "entrez_gene", "direction"))
treg_sig.dir <- treg_sig
treg_sig.dir$gs_name <- with(treg_sig.dir, paste0(gs_name, "_", direction))
treg_sig <- rbind(treg_sig, treg_sig.dir)

#########################################################
#### 1. PMID 30388456: Sade-Feldman et al. Cell 2018 ####
##  Defining T Cell States Associated with Response to ##
## Checkpoint Immunotherapy in Melanoma                ##
proj = 'pmid_30388456'
setwd(file.path(PDIR, proj))

mtx <- read.table("All_cells_counts_GEO.txt", header = T, sep="\t",
                  check.names=F, stringsAsFactors = F, fill=T)
colids <- colnames(mtx)[-1]
metadata <- as.character(mtx[1,-ncol(mtx)])
mtx <- mtx[-1,-ncol(mtx)]
colnames(mtx) = colids
mtx <- as.matrix(mtx)
storage.mode(mtx) <- 'numeric'

metadata3 <- read.table("GSE120575_patient_ID_single_cells.txt", header=T,
                        sep="\t", check.names = F, stringsAsFactors = F)
metadata2 <- read.csv("cluster_cell_name_mapping.csv", header=T,
                        check.names = F, stringsAsFactors = F) %>%
  mutate("treatment" = metadata)
metadata2 <- cbind(metadata2, metadata3)
cluster_label <- c('1'='Bcells',
                   '2'='Plasma_cells',
                   '3'='Monocyte.Macrophages',
                   '4'='Dendritic_cells',
                   '5'='Lymphocytes',
                   '6'='Exhausted_CD8_Tcells',
                   '7'='Regulatory_Tregs',
                   '8'='Cytotoxicity_Lymphocytes',
                   '9'='Exhausted.HS_CD8_Tcells',
                   '10'='Memory_Tcells',
                   '11'='Lymphocytes_exhausted.Cellcycle')
metadata2$Clusters <- cluster_label[as.character(metadata2$`Cluster Number`)]

seu <- CreateSeuratObject(counts = mtx, project = proj)
seu <- Seurat::AddMetaData(seu, metadata2)

saveRDS(seu, file=file.path(".", "seuobj.rds"))



seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
seu <- SCTransform(seu,
                   assay = 'RNA',
                   new.assay.name = 'SCT',
                   vst.flavor='v2',
                   vars.to.regress = c('percent.mt'))

seu <- CellCycleScoring(seu, 
                        s.features = (cc.genes$s.genes), 
                        g2m.features = (cc.genes$g2m.genes), 
                        set.ident = TRUE,
                        assay='SCT')

# seu <- split(seu, seu$sample)
seu <- SCTransform(seu, vst.flavor = "v2", verbose = TRUE,
                   assay='RNA', 
                   new.assay.name='SCT',
                   # variable.features.n = 1000,
                   vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 
                                       "S.Score", "G2M.Score"))
seu <- seu %>% 
  ScaleData(.)  %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE, 
          reduction.name='umap.unintegrated')

saveRDS(seu, file=file.path(".", "seuobj.rds"))



### TReg suppressive signature
seu <- readRDS(file=file.path(".", "seuobj.rds"))
expr <- GetAssayData(seu, slot='counts')
expr <- expr[rowSums(expr>0)>=(ncol(seu) > 0.05), ]


scores <- aucellFun(msig_ds = treg_sig, 
                    expr_mat=expr, gm,
                    mapfrom='SYMBOL', mapto='SYMBOL')
seu_aucell <- CreateSeuratObject(counts = assay(scores), project = proj)
seu_aucell <- Seurat::AddMetaData(seu_aucell, seu@meta.data)
seu_aucell@assays$RNA$data <- seu_aucell@assays$RNA$counts

id='Regulatory_Tregs'
seu_aucell_id <- subset(seu_aucell, ident=id)
# features <- gsub("_", "-", rownames(scores)) 
# df <- FindMarkers(seu_aucell_id, features=features, 
#                   group.by='response',
#                   ident.1='Responder', ident.2='Non-responder',
#                   min.pct=0, logfc.threshold=0)
pdf("~/xfer/aucell_treg.pdf")
df <- as.data.frame(t(GetAssayData(seu_aucell_id))) %>%
  dplyr::mutate("response"=seu_aucell_id$response)
my_comparisons <- list(unique(df$response))
ggpubr::ggviolin(reshape2::melt(df), x = "response", y = "value", fill = "response",
                 add = "boxplot", add.params = list(fill = "white"),
                 facet.by='variable',
                 title=id, ylab='AUCell_score', xlab='')+
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  ggpubr::stat_compare_means(label.y = 0.5) + 
  theme(axis.text.x = element_blank())
dev.off()

seu <- PrepSCTFindMarkers(seu, assay = "SCT", verbose = TRUE)
Idents(seu) <- 'Clusters'
seu_i <- subset(seu, ident='Regulatory_Tregs')
markers <- FindMarkers(object = seu_i, features='IL1RL1',
                       ident.1 = "Responder", ident.2="Non-responder",
                       group.by='response', assay = "SCT")
pdf("~/xfer/il1rl1.pdf")
scCustomize::VlnPlot_scCustom(seurat_object = seu_i, features = "IL1RL1", 
                              group.by='response', )
GetAssayData(seu_i)['IL1RL1',]
dev.off()


seu <- readRDS(file=file.path(".", "seuobj.rds"))

seutmp <- seu
cols <- c('orig.ident', 'Cell Name', 'Cluster Number', 'treatment', 'Phase',
          'Clusters', 'response')
meta <- seu@meta.data[,cols]
seutmp@meta.data <- meta
seutmp[['RNA']] <- CreateAssayObject(counts=seutmp@assays$SCT$counts)
DefaultAssay(seutmp) <- 'RNA'
obj <- loupeR:::select_assay(seutmp)$RNA
counts = loupeR:::counts_matrix_from_assay(obj)
counts <- rbind(counts, rbind(assay(scores)))
projections = loupeR:::select_projections(seutmp)
makeLoupe(counts, meta, projections,
          output_dir = ".", output_name = "pmid_30388456", 
          tmpdir="./")
file.copy(file.path("pmid_30388456.cloupe"),
          to = "~/xfer", overwrite = T)
cat("xfer pmid_30388456.cloupe\n")



#########################################################
#### 2. PMID31792460: Liu et al. Nature Med 2019     ####
## Integrative molecular and clinical modeling of      ##
## clinical outcomes to PD1 blockade in patients with  ##
## metastatic melanoma                                 ##
proj = 'pmid_31792460'
setwd(file.path(PDIR, proj))

tpm <- read.table(file.path(".", "41591_2019_654_MOESM3_ESM.txt"),
                  header=T, sep="\t") %>%
  tibble::column_to_rownames(., "X") %>% 
  t %>% as.data.frame
metadata <- read.table(file.path(".", "41591_2019_654_MOESM4_ESM.csv"),
                       header=T, sep=",") %>%
  tibble::column_to_rownames(., "Patient")
survival <- metadata %>%
  tibble::rownames_to_column(., "sample") %>%
  dplyr::select(sample, progressed, PFS, dead, OS) %>%
  dplyr::mutate(dss=0, dss2=0, x=0, x1=0) %>%
  magrittr::set_colnames(., c('sample', 'PFI', 'PFI.time', 'OS', 'OS.time',  
                              'DSS', 'DSS.time','DFI', 'DFI.time'))
  
ssgseaparam <- GSVA::ssgseaParam(as.matrix(tpm), 
                                 with(treg_sig, split(entrez_gene, gs_name)))
ssgsea_dat <- GSVA::gsva(ssgseaparam)


metrics <- c('PFI', 'OS')
for(measure in c('med', 'q')){
  msplit <- if(measure=='q') c(0, 0.25, 0.5, 0.75, 1) else c(0, 0.5, 1.0)
  surv <- lapply(metrics, function(metric){
    pdf(paste0("./", proj, ".", tolower(metric), "_scurve.", measure, ".pdf"), width = 4.5)
    pval <- list()
    for(geneseti in rownames(ssgsea_dat)){
      quant_spl <- split(colnames(ssgsea_dat), 
                         f=cut(ssgsea_dat[geneseti,], breaks=quantile(ssgsea_dat[geneseti,], msplit)))
      quant_spl <- quant_spl[c(1, length(quant_spl))]
      x=XenaTCGAanalyze_SurvivalKM(clinical_patient=survival, data_grp=quant_spl,
                                   Survresult = T, metric=metric, 
                                   caption=paste0(proj, " [", metric, "]: ", geneseti),
                                   xticks=7)
      pval[[geneseti]] <- x$pval
      print(x$tab)
    }
    dev.off()
    return(pval)
  })
  names(surv) <- metrics
  for(metric in metrics){
    file.copy(paste0("./", proj, ".", tolower(metric), "_scurve.", measure, ".pdf"), 
              to="~/xfer", overwrite = T)  
    message(paste0("xfer ", proj, ".", tolower(metric), "_scurve.", measure, ".pdf\n"))
  }
}
saveRDS(surv, file="surv.rds")

response_score[[proj]] <- data.frame(
  'group'=proj,
  "score"=ssgsea_dat[1,],
  'response'=metadata[colnames(ssgsea_dat),]$BR)

#########################################################
#### 3. PMID30753825: Gide et al. Cancer Cell 2019   ####
## Distinct Immune Cell Populations Define Response    ##
## to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4  ##
## Combined Therapy                                    ##
proj <- 'pmid_30753825'
setwd(file.path(PDIR, proj))

# Load in preiminary metadata data
metadata <- read.table(file.path(".", "1-s2.0-S1535610819300376-mmc2.csv"),
                       header=T, sep=",") %>%
  dplyr::mutate(Patient=paste0("Patient", `Patient.no.`),
                pfi=ifelse(Progressed=='Yes', 1, 0), 
                os=ifelse(`Last.Followup.Status`=='Dead',1 , 0)) %>%
  tibble::column_to_rownames(., "Patient") 

## Conver count matrix to a rough estimate of TPM based on median gene length
# genelength file created using GTFtools v0.90 on GRCh38 GTF file
genelength_f <- '/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/GTF/genome.gtf.genelength'
genelength <- read.table(genelength_f, sep="\t", header=T)
ensgene_map <- read.table(file.path("Gide_Quek_CancerCell2019", 
                              "cancercell_normalized_counts_genenames.txt"),
                  header=T, sep="\t")
ensgene_map <- with(ensgene_map, setNames(Gene, ID))
cnts <- read.table(file.path("Gide_Quek_CancerCell2019",
                             "PD1-IPIPD1_counts.txt"),
                   header=T, sep="\t") %>%
  tibble::column_to_rownames(., "ID")
genelength <- genelength[match(rownames(cnts), genelength$gene),]
na.idx <- which(is.na(genelength$gene))
genelength <- genelength[-na.idx,]
cnts <- cnts[-na.idx,]

tpm=DGEobj.utils::convertCounts(as.matrix(cnts),'TPM',genelength$median)
mapped_ids <- ensgene_map[rownames(tpm)]
non.na.idx <- which(!is.na(mapped_ids))
tpm <- tpm[non.na.idx,] %>%
  magrittr::set_rownames(., ensgene_map[rownames(tpm)][non.na.idx])


## Refine the metadata based on the actual sample names
mtmp <- strsplit(colnames(tpm), split='_') %>% do.call(rbind, .) %>%
  as.data.frame %>% 
  magrittr::set_colnames(., c('pd1_ipi', 'Patient.no.', 'pre_edt'))  %>%
  dplyr::mutate(ID=colnames(tpm))
metadata<- left_join(mtmp, 
                     dplyr::mutate(metadata, `Patient.no.`=as.character(`Patient.no.`)), 
                     by='Patient.no.')

survival <- metadata %>%
  dplyr::select(ID, pfi, Progression.Free.Survival..Days., 
                os, Overall.Survival..Days.) %>%
  dplyr::mutate(dss=0, dss2=0, x=0, x1=0) %>%
  magrittr::set_colnames(., c('sample', 'PFI', 'PFI.time', 'OS', 'OS.time',  
                              'DSS', 'DSS.time','DFI', 'DFI.time'))



ssgseaparam <- GSVA::ssgseaParam(as.matrix(tpm), 
                                 with(treg_sig, split(entrez_gene, gs_name)))
ssgsea_dat <- GSVA::gsva(ssgseaparam)
ssgsea_l <- split(as.data.frame(t(ssgsea_dat)), 
                  f=grepl("EDT$", colnames(ssgsea_dat)))
ssgsea_l2 <- split(ssgsea_l[['FALSE']], 
                   f=grepl("^ipiPD1", rownames(ssgsea_l[['FALSE']])))
ssgsea_l <- list("PRE"=t(ssgsea_l[['FALSE']]), 
                 "PRE_ipiPD1"=t(ssgsea_l2[['TRUE']]), 
                 "PRE_PD1"=t(ssgsea_l2[['FALSE']]), 
                 "EDT"=t(ssgsea_l[['TRUE']]))


measure='med'
metric <- 'PFI'
for(metric in c('OS', 'PFI')){
  for(measure in c('med', 'q')){
    msplit <- if(measure=='q') c(0, 0.25, 0.5, 0.75, 1) else c(0, 0.5, 1.0)
    pdf(paste0("./", proj, ".", tolower(metric), "_scurve.", measure, ".pdf"), width = 3.5)
    tukeys <- lapply(names(ssgsea_l), function(id){
      ssgsea_dat <- ssgsea_l[[id]]
      
      survival_i <- survival %>%
        dplyr::filter(sample %in% colnames(ssgsea_dat))
      
      pval <- list()
      tukey <- list()
      for(geneseti in rownames(ssgsea_dat)){
        # model=lm( ssgsea_dat[geneseti,] ~ metadata_i$Response )
        # tukey[[geneseti]] <- TukeyHSD(x=aov(model), 'metadata_i$Response', conf.level=0.95)
        quant_spl <- split(colnames(ssgsea_dat), 
                           f=cut(ssgsea_dat[geneseti,], breaks=quantile(ssgsea_dat[geneseti,], msplit)))
        quant_spl <- quant_spl[c(1, length(quant_spl))]
        x=XenaTCGAanalyze_SurvivalKM(clinical_patient=survival_i, data_grp=quant_spl,
                                     Survresult = T, metric=metric, 
                                     caption=paste0(id, " [", metric, "]: ", geneseti),
                                     xticks=7)
        print(x$tab)
        pval[[geneseti]] <- x$pval
      }
      return(list("surv"=pval))
    }) %>%
      setNames(., names(ssgsea_l))
    dev.off()
    file.copy(paste0("./", proj, ".", tolower(metric), "_scurve.", measure, ".pdf"), 
              to="~/xfer", overwrite = T)
    cat(paste0("xfer ", proj, ".", tolower(metric), "_scurve.", measure, ".pdf\n"))
    
    surv <- lapply(tukeys, function(i) i$surv)
  }
}
saveRDS(surv, file="surv.q.rds")


response_score[[paste0(proj, "_PrePD1")]] <- data.frame(
  'group'=paste0(proj, "_PrePD1"),
  "score"=ssgsea_dat[1,grep("^PD1.*PRE", colnames(ssgsea_dat))],
  'response'=metadata[match(grep("^PD1.*PRE", colnames(ssgsea_dat), value=T),metadata$ID),]$Progressed)
response_score[[paste0(proj, "_PreipiPD1")]] <- data.frame(
  'group'=paste0(proj, "_PreipiPD1"),
  "score"=ssgsea_dat[1,grep("^ipiPD1.*PRE", colnames(ssgsea_dat))],
  'response'=metadata[match(grep("^ipiPD1.*PRE", colnames(ssgsea_dat), value=T),metadata$ID),]$Progressed)

#########################################################
#### 4. PMID29033130: Riaz et al. Cell 2017          ####
## Tumor and Microenvironment Evolution during         ##
## Immunotherapy with Nivolumab                        ##
proj <- 'pmid_29033130'
setwd(file.path(PDIR, proj))

fpkm <- read.csv(file.path(".", "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"),
                    header=T, sep=",")
fpkm.tmp <- fpkm %>% 
  dplyr::mutate(gene=gm$ENTREZID$SYMBOL[as.character(X)],
                ens=gm$ENTREZID$ENSEMBL[as.character(X)],
                biotype=gm$ENSEMBL$gene_biotype[as.character(ens)]) %>%
  dplyr::filter(biotype == 'protein_coding',
                !is.na(gene))
fpkm.spl <- split(dplyr::select(fpkm.tmp, -c(X, gene, ens, biotype)),
                  f=fpkm.tmp$gene)
fpkm <- lapply(fpkm.spl, colMeans, na.rm=T) %>%
  do.call(rbind, .)

fpkm_on <- as.data.frame(fpkm) %>% dplyr::select(., contains("_On_")) %>%
  dplyr::rename_with(., ~make.unique(gsub("_.*", "", .))) %>%
  dplyr::select(-ends_with(".1"))
fpkm_pre <- as.data.frame(fpkm) %>% dplyr::select(., contains("_Pre_")) %>%
  dplyr::rename_with(., ~make.unique(gsub("_.*", "", .))) %>%
  dplyr::select(-ends_with(".1"))

fpkm_l <- list("pre"=fpkm_pre, "on"=fpkm_on)

metadata <- read.table(file.path(".", "NIHMS907788-supplement-9.csv"),
                       header=T, sep=",") %>%
  dplyr::mutate(os=ifelse(`Dead_or_Alive`==TRUE,1 , 0),
                OS.time=as.integer(Time.to.Death.weeks * 7)) %>%
  tibble::column_to_rownames(., "Patient") 




metric <- 'OS'
for(measure in c('med', 'q')){
msplit <- if(measure=='q') c(0, 0.25, 0.5, 0.75, 1) else c(0, 0.5, 1.0)
pdf(paste0("./", proj, ".", tolower(metric), "_scurve.", measure, ".pdf"), width = 3.5)
tukeys <- lapply(names(fpkm_l), function(fpkm_id){
  fpkm <- fpkm_l[[fpkm_id]]
  ssgseaparam <- GSVA::ssgseaParam(as.matrix(fpkm), 
                                   with(treg_sig, split(entrez_gene, gs_name)))
  ssgsea_dat <- GSVA::gsva(ssgseaparam)
  
  metadata_i <- metadata[match(colnames(fpkm), rownames(metadata)),]
  survival_i <- metadata_i %>%
    tibble::rownames_to_column(., "sample") %>%
    dplyr::select(sample, os, OS.time) %>%
    dplyr::mutate(pfi=0, pfit=0, dss=0, dsst=0, dfi=0, dfit=0) %>%
    magrittr::set_colnames(., c('sample', 'OS', 'OS.time',  'PFI', 'PFI.time', 
                                'DSS', 'DSS.time','DFI', 'DFI.time'))
  
  pval <- list()
  tukey <- list()
  for(geneseti in rownames(ssgsea_dat)){
    
    model=lm( ssgsea_dat[geneseti,] ~ metadata_i$Response )
    tukey[[geneseti]] <- TukeyHSD(x=aov(model), 'metadata_i$Response', conf.level=0.95)
    
    quant_spl <- split(colnames(ssgsea_dat), 
                       f=cut(ssgsea_dat[geneseti,], breaks=quantile(ssgsea_dat[geneseti,], msplit)))
    quant_spl <- quant_spl[c(1, length(quant_spl))]
    x=XenaTCGAanalyze_SurvivalKM(clinical_patient=survival_i, data_grp=quant_spl,
                                 Survresult = T, metric=metric, 
                                 caption=paste0(fpkm_id, " [", metric, "]: ", geneseti),
                                 xticks=4.5)
    print(x$tab)
    pval[[geneseti]] <- x$pval
  }
  return(list("tukey"=tukey, "surv"=pval))
}) %>%
  setNames(., names(fpkm_l))
dev.off()
file.copy(paste0("./", proj, ".", tolower(metric), "_scurve.", measure, ".pdf"), 
          to="~/xfer", overwrite = T)
cat(paste0("xfer ", proj, ".", tolower(metric), "_scurve.", measure, ".pdf\n"))

surv <- lapply(tukeys, function(i) i$surv) 
}
saveRDS(surv, file="surv.q.rds")



fpkm <- fpkm_l[['on']]
ssgseaparam <- GSVA::ssgseaParam(as.matrix(fpkm), 
                                 with(treg_sig, split(entrez_gene, gs_name)))
ssgsea_dat <- GSVA::gsva(ssgseaparam)
response_score[[paste0(proj)]] <- data.frame(
  'group'=paste0(proj),
  "score"=ssgsea_dat[1,],
  'response'=metadata[colnames(ssgsea_dat),]$Response)

################################
#### 5. Manual calculations ####

# Meta-P Fishers method
# Treg_Sig, median
pvals <- c("riaz_on"=0.31,
           "liu"=0.059,
           "gide_pre_pd1"=0.01,
           "gide_pre_ipipd1"=0.78)
metap::sumlog(pvals[-4]) ## w/o ipiPD1: 0.008531024; w/ ipiPD1: 0.02350837 

# Treg_Sig, top-bottom quantile
pvals <- c("riaz_on"=0.054,
           "liu"=0.028,
           "gide_pre_pd1"=0.023,
           "gide_pre_ipipd1"=0.23)
metap::sumlog(pvals) ## w/o ipiPD1: 0.002224549, w/ipiPD1 0.002807715 



## Visualizing response score in context of recist progression
X <- as.data.frame(do.call(rbind, response_score)) %>%
  dplyr::filter(!is.na(response))
pdf("~/xfer/response_score.pdf")
lapply(split(X, f=grepl("pmid_30753825", X$group)), function(Xi){
  print("A")
  my_comparisons <- combn(names(which(table(Xi$response) > 2)), 2) %>% 
    apply(., 2, list) %>% unlist(., recursive=F)
  # my_comparisons <- list(combn(unique(Xi$response), 2))
  ggpubr::ggviolin(Xi, x = "response", y = "score", fill = "response",
                   add = "boxplot", add.params = list(fill = "white"),
                   facet.by='group',
                   ylab='AUCell_score', xlab='')+
    ylim(0,1.75) +
    ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
    ggpubr::stat_compare_means(label.y = 1.6) + 
    theme(axis.text.x = element_blank())
}) %>%
  cowplot::plot_grid(plotlist=., nrow=2)
dev.off()

#############################################
#### 6. TCGA: Melanoma and Breast cancer ####
DATADIR <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcga_treg_sig/data"
OUTDIR <- file.path(DATADIR, "..", "results")
survival_f <- file.path(DATADIR, "phenotype", "TCGA_survival_data.txt")
survival <- read.table(survival_f, sep="\t", header = T)

## Read in TPM
tpm_f <- file.path(DATADIR, "rna", "TcgaTargetGtex_rsem_gene_tpm.gz")
probemap_f <- file.path(DATADIR, "rna", "probeMap%2Fgencode.v23.annotation.gene.probemap")

probemap <- read.table(probemap_f, sep="\t", header = T)
probemap_filt <- probemap %>% 
  mutate(biotype =  gm$ENSEMBL$gene_biotype[gsub("\\.[0-9]$", "", id)]) %>%
  dplyr::filter(biotype == 'protein_coding' | is.na(biotype))  %>%
  dplyr::filter(!duplicated(gene))

tpm <- read.table(tpm_f, sep="\t", header = T)
tpm_mat <- .cleanExpr(tpm, probemap_filt, samplecol='sample')
tcga_tpm <- tpm_mat[,grep("TCGA", colnames(tpm_mat))]

## Read in TCGA phenotypes
phenotype_f <- file.path(DATADIR, "phenotype", "TcgaTargetGTEX_phenotype.txt.gz")
phenotype <- read.table(phenotype_f, sep="\t", header = T)
tcga_phenotype <- phenotype %>% dplyr::filter(X_study == 'TCGA')

tcga_phenotype_filt  <- tcga_phenotype %>% 
  dplyr::filter(!X_sample_type %in% c('Control analyte', 'Additional Metastatic', 
                                      'Solid Tissue Normal', 'Primary Blood Derived Cancer - Peripheral Blood')) %>%
  dplyr::filter(!X_primary_site %in% c('White blood cell'))
tcga_phenotype_filt  <- tcga_phenotype %>% 
  dplyr::filter(!X_sample_type %in% c('Solid Tissue Normal'))
phenospl <- with(tcga_phenotype_filt, split(sample, f=primary.disease.or.tissue))
phenospl <- phenospl[grep("(Breast|Melanoma)", names(phenospl), value=T)]


#--- A) Survival Curves ----
## TCGA survival curves
Survresult <- T
breakid <- 'median' # 'median', 'Q14'
breaks <- if(breakid == 'median') c(0, 0.5, 1) else c(0, 0.25, 0.75, 1)
nbins <- 2
metric <- 'OS'
ggs <- lapply(names(phenospl)[-3], function(phenoid){
  message(phenoid)
  phenospl_i <- phenospl[[phenoid]]
  phenospl_i <- phenospl_i[(phenospl_i %in% colnames(tcga_tpm))]
  gsvapar_i <- GSVA::ssgseaParam(tcga_tpm[,phenospl_i], 
                                 with(treg_sig, split(entrez_gene, gs_name)))
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
                                    caption=paste0("TCGA ", phenoid, " [", metric, "]: ", geneseti),
                                    xticks=7)
    return(ggsc$tab)
  })
  
  
  return(gg_scs)
})
pdf(paste0("~/xfer/tcga_cancer_specific.", metric, ".", breakid, ".pdf"), width=4.5)
ggs
dev.off()
cat(paste0("xfer tcga_cancer_specific.", metric, ".", breakid, ".pdf\n"))

#--- B) Volcano Plots ----
# CIBERSORT
cibersort_path <- '/cluster/projects/mcgahalab/ref/immunedeconv/cibersort'
immunedeconv::set_cibersort_binary(file.path(cibersort_path, "CIBERSORT.R"))
immunedeconv::set_cibersort_mat(file.path(cibersort_path, "LM22.txt"))

ssgseaparam <- GSVA::ssgseaParam(tcga_tpm, with(treg_sig, split(entrez_gene, gs_name)))
ssgsea_dat <- GSVA::gsva(ssgseaparam)

# CIBERSORT per cancertype
tcga_phenotype_filt  <- tcga_phenotype %>% 
  dplyr::filter(!X_sample_type %in% c('Solid Tissue Normal'))
phenospl <- with(tcga_phenotype_filt, split(sample, f=primary.disease.or.tissue))
phenospl <- phenospl[grep("(Breast|Melanoma)", names(phenospl), value=T)]

cibersort_dat <- lapply(names(phenospl)[-3], function(phenoid){
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
names(cibersort_dat) <- names(phenospl)[-3]
cibersort_dat  <- cibersort_dat[!sapply(cibersort_dat, is.null)]
cibersort_mat <- as.data.frame(do.call(cbind, cibersort_dat)) %>% 
  dplyr::rename_with(., ~make.unique(gsub("^.*?\\.", "", .))) %>%
  tibble::column_to_rownames(., "cell_type") 


cibersort_cor <- lapply(names(phenospl)[-3], function(phenoid){
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
names(cibersort_cor) <- names(phenospl)[-3]
cibersort_cor <- cibersort_cor[!sapply(cibersort_cor, is.null)]

genesets <- c('Treg_Sig')
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
  xlim(-0.51, 0.51) +
  ylim(0, 70) +
  # labs(title = "Gene expression changes in diseased versus healthy samples",
  #      x = "log2(fold change)",
  #      y = "-log10(adjusted P-value)",
  #      colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position='none') 

pdf("~/xfer/tcga_volcano_plots.pdf", width = 10, height = 5); final_plot; dev.off()



cor_df <- reshape2::melt(cibersort_cor) %>%
  magrittr::set_colnames(., c('Celltype', 'Geneset', 'corr', 'Disease'))
write.table(cor_df, file="~/xfer/cibersort_corr.csv",
            sep=",", quote = F, col.names = T, row.names = F)

#####################################
#### 7. GTEx: non-Malignant Skin #### 
GTEX.PDIR='/cluster/projects/mcgahalab/ext_data/gtex'
DATADIR <- "/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcga_treg_sig/data"
OUTDIR <- file.path(DATADIR, "..", "results")
probemap_f <- file.path(DATADIR, "rna", "probeMap%2Fgencode.v23.annotation.gene.probemap")

probemap <- read.table(probemap_f, sep="\t", header = T)
probemap_filt <- probemap %>% 
  mutate(biotype =  gm$ENSEMBL$gene_biotype[gsub("\\.[0-9]$", "", id)]) %>%
  dplyr::filter(biotype == 'protein_coding' | is.na(biotype))  %>%
  dplyr::filter(!duplicated(gene))


gtex_f <- file.path(GTEX.PDIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
gtex_rds_f <- file.path(GTEX.PDIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.rds")
gtex.tpm <- read.table(gtex_f, sep="\t", header=T, skip = 2)
biotypes <- gm$ENSEMBL$gene_biotype[gsub("\\.[0-9]*$", "", gtex.tpm$Name)]
idx <- which(biotypes == 'protein_coding' &
               !duplicated(gm$ENSEMBL$SYMBOL[gsub("\\.[0-9]*$", "", gtex.tpm$Name)]))
gtex.tpm <- gtex.tpm[idx,] %>%
  tibble::remove_rownames() %>%
  mutate(symbol=gm$ENSEMBL$SYMBOL[gsub("\\.[0-9]*$", "", gtex.tpm[idx,]$Name)]) %>%
  tibble::column_to_rownames(., "symbol") %>%
  dplyr::select(-c(Name, Description))

gtex.tpm_mat <- .cleanExpr(gtex.tpm %>% tibble::rownames_to_column("gene"), 
                           probemap_filt, samplecol='gene', id2='gene')

gtex.meta <- read.table(file.path(GTEX.PDIR, "gtex_meta.txt"),
                        sep="\t", header=T, check.names = F, stringsAsFactors = F, comment.char = "")

gtex_meta <- gtex.meta %>% 
  dplyr::filter(SAMPID %in% colnames(gtex.tpm_mat)) %>%
  dplyr::filter(SMTS %in% 'Skin')
gtex_tpm <- gtex.tpm_mat[,gtex_meta$SAMPID]

# CIBERSORT
cibersort_path <- '/cluster/projects/mcgahalab/ref/immunedeconv/cibersort'
immunedeconv::set_cibersort_binary(file.path(cibersort_path, "CIBERSORT.R"))
immunedeconv::set_cibersort_mat(file.path(cibersort_path, "LM22.txt"))

ssgseaparam <- GSVA::ssgseaParam(gtex_tpm, with(treg_sig, split(entrez_gene, gs_name)))
ssgsea_dat_gtex <- GSVA::gsva(ssgseaparam)


#--- A) Volcano plot ----
phenospl <- with(gtex_meta, split(SAMPID, f=SMTSD))#SMTS))

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

genesets <- c('Treg_Sig')
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
  xlim(-0.51, 0.51) +
  ylim(0,70) +
  # labs(title = "Gene expression changes in diseased versus healthy samples",
  #      x = "log2(fold change)",
  #      y = "-log10(adjusted P-value)",
  #      colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position='none') 

pdf("~/xfer/gtex_volcano_plots.pdf", width = 10, height = 5); final_plot; dev.off()

# 
# 
# cor_df <- reshape2::melt(cibersort_cor) %>%
#   magrittr::set_colnames(., c('Celltype', 'Geneset', 'corr', 'Disease'))
# write.table(cor_df, file="~/xfer/cibersort_corr.csv",
#             sep=",", quote = F, col.names = T, row.names = F)
# 
