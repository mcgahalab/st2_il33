library(ggplot2)
library(reshape2)
library(oligo)

# dir='/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/meta_analysis/st2_il33_paper'
# rma_file <- file.path(dir, "results", "st2_il33.rma_norm.rds")

## Set the path to the rds file
rma_file <- "/path/to/st2_il33.rma_norm.rds"
## Set the path for the output PDF
outpdf <- "~/Dekstop/rnaarray.pdf"
## Set what genes should be plotted
genes <- c('Kit', 'Ccr1', 'Gata3', 'Adam8', 'Gpr15', 'Itga9', 'Myo1g',
           'Ccl8', 'Ednrb', 'Sell', 'Cxcr5', 'Ccl20', 'Selp', 'Ffar2',
           'Cnr2', 'B4galt1', 'Prkca', 'Vav1', 'Nup85', 'F2rl1', 'Itga6', 'S1pr1')
  
###################
#### Functions ####
getLog2FC <- function(dat, genes){
  ## Idnetifying indices of matcheed mouse data for ST2- and ST2+
  pdata_ord <- pData(dat)[with(pData(dat), order(mouseID)),]
  wt_idx <- pdata_ord[which(pdata_ord$STwt),]$index
  ko_idx <- pdata_ord[which(!pdata_ord$STwt),]$index
  
  ## Visualize the data
  
  genes <- rev(genes)
  probeid <- fData(dat)[match(genes, fData(dat)$SYMBOL),]$PROBEID
  log2expr_sel <- (exprs(dat)[probeid,wt_idx] - exprs(dat)[probeid,ko_idx]) %>%
    as.data.frame %>%
    mutate(ids=genes) %>% 
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames('ids')
  return(log2expr_sel)
}

visHeatmap <- function(dat){
  ggplot(reshape2::melt(t(dat)), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(limits = c(-3,3),
                         colours=c("cyan", "blue", "black", "yellow", "red"))
}

####################################################
#### Preprocess: Don't need to modify, just run ####
if(process_data){
  library(affycoretools)
  library(clariomsmousetranscriptcluster.db)
  library(limma)
  
  print("Processing data")
  dir='/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/meta_analysis/st2_il33_paper'
  setwd(file.path(dir, "data"))
  
  celfiles <- list.files(file.path(dir, "data"), pattern="CEL")
  ab <- read.celfiles(celfiles)
  eset <- rma(ab, background=TRUE, normalize=TRUE, subset=NULL)
  eset <- annotateEset(eset, clariomsmousetranscriptcluster.db, 
                       columns = c("ACCNUM","PROBEID", "ENTREZID", "SYMBOL", "GENENAME","ENSEMBL"))
  sampleNames(eset) <- gsub(".CEL$", "", sampleNames(eset))
  
  
  eset2 <- eset[!is.na(fData(eset)$SYMBOL),]
  avemat <- avereps(eset2, fData(eset2)$SYMBOL)
  ## and put back into an ExpressionSet
  eset3 <- eset2[match(row.names(avemat), fData(eset2)$SYMBOL),]
  row.names(avemat) <- row.names(eset3) # relabel symbols with probeids
  exprs(eset3) <- avemat
  
  ## Add in metadata
  meta <- read.table(file.path(dir, "data", 'E-MTAB-6842.sdrf.txt'), sep="\t",
                     header = T, check.names = F, stringsAsFactors = F,
                     comment.char = '', fill=T)
  pData(eset3) <- cbind(pData(eset3), meta[match(meta$`Source Name`, sampleNames(eset3)),]) %>% 
    as.data.frame %>% 
    dplyr::select(c('index', 'Characteristics[phenotype]', 'Characteristics[individual]')) %>%
    rename_with(., ~c('index', 'phenotype', 'mouseID')) %>%
    mutate("STwt"=grepl("ST2\\+", phenotype))
  
  ## Limma DEG
  design <- cbind(WT=1, ST2PosvsNeg=pData(eset3)$STwt)
  fit <- lmFit(eset3, design)
  fit <- eBayes(fit)
  # topTable(fit, coef="ST2PosvsNeg")
  saveRDS(eset3, file=rma_file)
} else {
  print("Reading in existing RMA data...")
  eset3 <- readRDS(rma_file)
}

#######################
#### Visualization ####
log2expr_sel <- getLog2FC(eset3, genes)

pdf(outpdf)
visHeatmap(log2expr_sel)
dev.off()
