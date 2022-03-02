library(BiocParallel)
library(GenomicFeatures)
library(harmony)
library(org.Hs.eg.db)
library(chromVAR)
library(ChIPseeker)
library(SingleR)
library(reshape2)
library(enrichplot)
library(Signac)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(plyr)
library(cowplot)
library(ggrastr)
# library("GOstats")
library(ggplot2)
library(GenomicRanges)
library("clusterProfiler")
library("enrichplot")
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
# library(ArchR)

seed <- 1234
set.seed(seed)
args = commandArgs(trailingOnly=TRUE)

sidx <- args[1] # 1  # Where to start preprocessing individual samples from
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/giselle_scATAC'
dataset <- 'pgmc'
macs2_path <- '/cluster/home/quever/miniconda3/envs/r4/bin/macs2'
pbmc_annot <- '/cluster/projects/mcgahalab/ref/scrna/pbmc_dataset/pbmc_multimodal.h5seurat'
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/GTF/genome.gtf'
panglao_dbf <- '/cluster/projects/mcgahalab/ref/scrna/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv'

count_blacklist_frags <- FALSE # counting balcklist from 10x fragments is extremely lengthy process

out_id <- 'hilo_iss'
outdir <- file.path(PDIR, "results", out_id)
# outdir <- file.path(outdir, "scatac")
dir.create(outdir, recursive = T, showWarnings = F)
groups <- list.files(file.path(PDIR, "datasets", dataset))
setwd(PDIR)

doublet_barcodes <- file.path(outdir, "doubletfinder", "doublet_barcodes.rds")
doublet_quantile_cutoff <- 0.95

###################
#### Functions ####
checkObj <- function(seu){
  seu$test <- rep(FALSE, ncol(seu))
  seu$test[1:5] <- TRUE
  subset(seu, subset=test)
}

appendDatasetId <- function(seus, id, datatype){
  i <- seus[[id]]
  seu_p <- i[[datatype]]
  seu_p$dataset <- id
  return(seu_p)
}
 
######################################################
#### 0.a Create Data-type specific Seurat object  ####
library(EnsDb.Hsapiens.v86)
genome <- EnsDb.Hsapiens.v86
dir.create(file.path("results", dataset, "seurat_obj"), 
           recursive = TRUE, showWarnings = FALSE)

## Merge peaks to get reduced-combined peak 
# https://satijalab.org/signac/articles/merging.html 
gr_pks <- lapply(groups, function(grp){
  samples <- list.files(file.path(PDIR, "datasets", dataset, grp))
  pks <- lapply(samples, function(s){
    pk_path <- file.path("datasets", dataset, grp, s)
    pk <- read.table(file = file.path(pk_path, "atac_peaks.bed"),
                     stringsAsFactors = FALSE, check.names=F,
                     sep="\t", comment.char="#", header=F,
                     col.names = c("chr", "start", "end"))
    grpk <- makeGRangesFromDataFrame(pk)
    return(grpk)
  })
  names(pks) <- samples
  return(pks)
})
names(gr_pks) <- groups
combined_pks <- reduce(unlist(as(unlist(gr_pks, recursive = T), "GRangesList")))
pk_width <- width(combined_pks)
combined_pks <- combined_pks[pk_width  < 10000 & pk_width > 20]

## Read in the 10x count data for each group and sample
seus_rds <- file.path(outdir, "seurat_obj", "samples", "seus.rds")
dir.create(file.path(outdir, "seurat_obj", "samples"), 
           recursive = T, showWarnings = F)

  
seus <- if(!file.exists(seus_rds)) lapply(groups[sidx:length(groups)], function(grp){
  samples <- list.files(file.path(PDIR, "datasets", dataset, grp))
  seus <- lapply(samples, function(s){
    rds_f <- file.path(outdir, "seurat_obj", "samples", paste0(s, ".rds"))
    if(file.exists(rds_f)){
      seu <- readRDS(rds_f)
    } else {
      pk_path <- file.path("datasets", dataset, grp, s)
      mtx <- Read10X(data.dir = file.path(pk_path, "filtered_feature_bc_matrix"), strip.suffix=TRUE)
      
      seu <- list()
      ## Gene expression seurat object
      {
        seu[['GEX']] <- CreateSeuratObject(counts = mtx[['Gene Expression']], project = s)
        seu[['GEX']]$group <- grp
      }
      
      ## ATAC seurat object
      {
        # Fix mismatch between fragment and counts
        # colnames(mtx[['Peaks']]) <- paste0(colnames(mtx[['Peaks']]), "-1") 
        
        # Read in metadata: 
        meta_atac <- read.csv(file = file.path(pk_path, "per_barcode_metrics.csv"), 
                              header = TRUE, row.names = 1)
        col_ids <- colnames(meta_atac)[grep("atac", colnames(meta_atac))]
        
        # create fragment objects
        frags_obj <- CreateFragmentObject(path = file.path(pk_path, "atac_fragments.tsv.gz"),
                                          cells = rownames(meta_atac), 
                                          max.lines = NULL)
        # Quantify peaks in each dataset using the combined reference
        frag_cnts <- FeatureMatrix(fragments = frags_obj,
                                   features = combined_pks,
                                   cells = rownames(meta_atac))
        chrom_assay <- CreateChromatinAssay(counts = frag_cnts, #mtx[['Peaks']],
                                            sep = c(":", "-"), 
                                            genome = seqinfo(genome), 
                                            fragments = frags_obj,
                                            min.cells = 10, 
                                            min.features = 200)
        seu_chrom <- CreateSeuratObject(chrom_assay, assay = "ATAC", 
                                        meta.data=meta_atac[,col_ids])
        seu[['Peaks']] <- RenameCells(seu_chrom, 
                                      new.names=gsub("-1$", "", Cells(seu_chrom)))
        
        # #Create chromatin assay from fragments
        # chrom_assay <- CreateChromatinAssay(counts = mtx[['Peaks']], 
        #                                     sep = c(":", "-"), 
        #                                     genome = seqinfo(genome), 
        #                                     fragments = file.path(pk_path, "atac_fragments.tsv.gz"), 
        #                                     min.cells = 10, 
        #                                     min.features = 200)
        # # Append cellRanger QC data to the seurat object
        # chrom_assay <- CreateSeuratObject(counts = chrom_assay,
        #                                   assay = "ATAC",
        #                                   meta.data = meta_atac[,col_ids])
        
        # Rename cells here so the Fragments barcodes also rename their mapping
        # 123456-1 => [group]_123456
        # seu[['Peaks']] <- RenameCells(chrom_assay, 
        #                               new.names=gsub("-1$", "", Cells(chrom_assay)))
      }
      
      saveRDS(seu, file=rds_f)
    }
    return(seu)
  })
})
if(!file.exists(seus_rds)){
  names(seus) <- groups
  seus <- unlist(seus, recursive = F)
  
  # Create sample<->group mapping key
  grp_sample_map <- unlist(sapply(groups, function(grp) list.files(file.path(PDIR, "datasets", dataset, grp))))
  sample_grp_map <- setNames(names(grp_sample_map), grp_sample_map)
  sample_grp_map <- gsub("[1-9]*$", "", sample_grp_map)
  names(seus) <- as.character(grp_sample_map)
  saveRDS(seus, file=seus_rds) 
} else{
  seus <- readRDS(file=seus_rds)
}


## Merge the scRNA data and save
dir.create(file.path(outdir, "seurat_obj"), showWarnings = F)
datatype <- 'GEX'
seus_rna <- lapply(names(seus), function(id) { appendDatasetId(seus, id, datatype) })

seu <- merge(x=seus_rna[[1]], y = seus_rna[-1], 
             add.cell.ids = names(seus), 
             project = dataset)
SeuratDisk::SaveH5Seurat(
  seu, filename = file.path(
    outdir, "seurat_obj",paste0(gsub("\\s", "_", datatype), ".h5seurat")
  ), 
  overwrite=TRUE
)
seu_gex <- seu


## Merge the scATAC data and save
datatype <- 'Peaks'
seus_peaks <- lapply(names(seus), function(id) { appendDatasetId(seus, id, datatype) })
seu <- merge(x=seus_peaks[[1]], y = seus_peaks[-1], 
             add.cell.ids = names(seus), 
             project = dataset)
saveRDS(seu, file=file.path(outdir, "seurat_obj",paste0(gsub("\\s", "_", datatype), ".rds")))
seu_peaks <- seu; rm(seu);

#############################################################
#### 0.b Doublet Finding on a per-sample basis: RUN ONCE ####
if(file.exists(doublet_barcodes)) warning("Doublet files already identified")
seus_rds <- file.path(outdir, "seurat_obj", "samples", "seus.rds")
if(!exists("seus")) seus <- readRDS(file=seus_rds)

## Apply DoubletFinder to every scRNA sample independently
datatype <- 'GEX'
seus_rna <- lapply(names(seus), function(id) { appendDatasetId(seus, id, datatype) })
dir.create(file.path(outdir, "doubletfinder", "samples"), recursive = T, showWarnings = F)
seus_rna_doublet <- lapply(rev(seus_rna), function(seu){
  id <- as.character(unique(seu$orig.ident))
  doublet_rds_f <- file.path(outdir, "doubletfinder", "samples", paste0(id, ".rds"))
  if(file.exists(doublet_rds_f)){
    print("Loading existing analysis")
    doublet <- readRDS(doublet_rds_f)
    return(doublet)
  } 
  
  # pre-process
  print(paste0("> ", id, "..."))
  seu <- SCTransform(seu)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(object = seu, dims = 1:30)
  seu <- FindClusters(object = seu, resolution = 1.2)
  seu <- RunUMAP(seu, dims = 1:10)
  
  # cell-type identification
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/BlueprintEncodeData.rds") 
  lbl <- 'label.fine' #label.main
  blueprint_anno <- SingleR(test=GetAssayData(seu), ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            clusters=seu$seurat_clusters)
  cluster_ids <- setNames(blueprint_anno$labels, 
                          as.character(rownames(blueprint_anno)))
  seu@meta.data[,'blueprint'] <- cluster_ids[as.character(seu$seurat_clusters)]
  
  # pK identification
  # sweep_res_list <- paramSweep_v3(seu, PCs = 1:10, sct = TRUE)
  # sweep_seu <- summarizeSweep(sweep_res_list, GT = FALSE)
  # bcmvn <- find.pK(sweep_seu)
  # pk <- bcmvn$pK[which.max(bcmvn[as.numeric(as.character(bcmvn$pK))<0.15,]$BCmetric)]
  # pk <- as.numeric(as.character(pk))
  pk <- 0.11
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(seu$blueprint)    
  nExp_poi <- round(0.075*nrow(seu@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
   
  seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pk,  # first bump under 0.1 from bcmvn_seu
                          nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  ## Assemble proportions of doublets
  meta <- seu@meta.data
  df_col <- grep("^DF\\.", colnames(meta), value=T)
  doublet_prop <- table(meta[,c('seurat_clusters', df_col)])
  doublet_df <- data.frame("prop"=round(doublet_prop[,1] / rowSums(doublet_prop),2),
                           "ncells"=rowSums(doublet_prop))
  cluster_df <- unique(meta[,c('seurat_clusters', 'blueprint')])
  rownames(cluster_df) <- cluster_df$seurat_clusters
  doublet_df$quantile <- round(ecdf(sort(doublet_df$prop))(doublet_df$prop),2)
  doublet_df <- cbind(doublet_df, cluster_df[rownames(doublet_df),-1])
  doublet_df <- doublet_df[order(doublet_df$prop),]
  
  ## Plot the dimension plot of doublets by cluster by cell type
  dp_doublet <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                     group.by=df_col, pt.size=0.5, shuffle=TRUE)
  dp_anno <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                       group.by='blueprint', pt.size=0.5, shuffle=TRUE)
  dp_clus <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                     group.by='seurat_clusters', pt.size=0.5, shuffle=TRUE)
  gg_dp <- dp_doublet + dp_anno + dp_clus
  
  doublet_list <- list("seu"=seu, "doublet"=doublet_df, "gg"=gg_dp)
  saveRDS(doublet_list, file=doublet_rds_f)
  # pdf("~/xfer/id.pdf", width = 17)
  # gg_dp
  # dev.off()
  return(doublet_list)
})

## Plot out the DimPlot of each sample - seuratCluster,Annotation,Doublets
ids <- sapply(seus_rna_doublet, function(i) unique(i$seu$orig.ident))
names(seus_rna_doublet) <- ids
pdf(file.path(outdir, "doubletfinder", "doublets_samples.pdf"), width=20)
lapply(names(seus_rna_doublet), function(id) {
  seus_rna_doublet[[id]]$gg[[1]] <- seus_rna_doublet[[id]]$gg[[1]] + ggtitle(id)
  seus_rna_doublet[[id]]$gg
})
dev.off()

## Calculate a total ECDF of proportion of doublets, identify
## the total quantile for every sample-specific proportion
## and draw a threshold to identify clusters of high doublet proportions
doublets_l <- lapply(seus_rna_doublet, function(i) i$doublet)
doublet_ecdf <- ecdf(do.call(rbind, doublets_l)[,1])
doublets_l <- lapply(names(doublets_l), function(id){
  i <- doublets_l[[id]]
  i$cluster_num <- rownames(i)
  i$quantile_total <- doublet_ecdf(i$prop)
  i$sample <- id
  i[order(i$quantile_total),]
})
doublets_df <- as.data.frame(do.call(rbind, doublets_l))
write.table(doublets_df, file=file.path(outdir, "doubletfinder", "doublets_prop.csv"),
            sep=",", quote=F, col.names=T, row.names=F)

pdf(file.path(outdir, "doubletfinder", "doublets_prop.pdf"))
ggplot(doublets_df, aes(x=prop, y=quantile_total, size=ncells, color=sample)) +
  geom_point() + 
  geom_hline(yintercept=doublet_quantile_cutoff, color="red") +
  theme_classic() 
dev.off()

## Identify the cell-identities of clusters with high doublets
dubz_wins <- doublets_df[doublets_df$quantile_total > doublet_quantile_cutoff,]
samples <- split(dubz_wins, f=dubz_wins$sample)
barcodes_rm <- lapply(names(samples), function(s){
  dub_clusters <- samples[[s]]$cluster_num
  seu_s <- seus_rna_doublet[[s]]$seu
  
  # subset seurat for clusters to remove and return cell barcodes
  Idents(seu_s) <- 'seurat_clusters'
  seu_rm <- subset(seu_s, idents=dub_clusters)
  return(Cells(seu_rm))
})
names(barcodes_rm) <- names(samples)
saveRDS(barcodes_rm, file=doublet_barcodes)

######################################
#### 1.a Merge scATAC with scRNA  ####
library(EnsDb.Hsapiens.v86)
genome <- EnsDb.Hsapiens.v86
seu_obj_path <- file.path(outdir, "seurat_obj")
if(!exists("seu_gex")) seu_gex <- SeuratDisk::LoadH5Seurat(file.path(seu_obj_path, paste0('GEX', ".h5seurat")))
if(!exists("seu_peaks")) seu_peaks <- readRDS(file.path(seu_obj_path, paste0('Peaks', ".rds")))

# Combine the scRNA data with the scATAC data
seu_gex$hasATAC <- Cells(seu_gex) %in% Cells(seu_peaks)
seu_gex_tmp <- subset(seu_gex, subset=hasATAC)
seu_peaks$hasATAC <- Cells(seu_peaks) %in% Cells(seu_gex)
seu_peaks_tmp <- subset(seu_peaks, subset=hasATAC)

seu_integrated <- seu_gex_tmp
seu_integrated[['ATAC']] <- seu_peaks_tmp@assays$ATAC
seu_integrated <- AddMetaData(seu_integrated, metadata=seu_peaks_tmp@meta.data)
seu_integrated$orig.ident <- seu_gex_tmp$orig.ident
rm(seu_peaks_tmp, seu_gex_tmp)

# Annotate the scATAC data
DefaultAssay(seu_integrated) <- "ATAC"
annotations <- GetGRangesFromEnsDb(ensdb = genome)
annotations <- renameSeqlevels(annotations, c(paste0("chr", seqlevels(annotations))))
Annotation(seu_integrated) <- annotations

# Rename seqlevels [NCBI] to keep consistent with fragments [UCSC] (seqlevelsStyle needs NCBIlinker)
seq_idx <- match(c(1:22, "X", "Y", "MT"), seqnames(seu_integrated))
seq_ids <- seqnames(seu_integrated)
seq_ids[seq_idx] <- paste0("chr", seq_ids)
seu_integrated <- renameSeqlevels(seu_integrated, seq_ids)

saveRDS(seu_integrated, file = file.path(outdir, "seurat_obj", "scIntegrated.rds"))
rm(seu_gex, seu_peaks)

#############################
#### 1.b Remove Doublets ####
if(!file.exists(doublet_barcodes)) {
  stop("Run DoubletFinder in Step 0.b before merging and pre-processing data")
} else {
  barcodes_rm <- readRDS(doublet_barcodes)
}
if(!exists("seu_integrated")) seu_integrated <- readRDS(file.path(outdir, "seurat_obj", "scIntegrated.rds"))

# Remove doublets
cells_to_remove <- sapply(names(barcodes_rm), function(id){
  paste0(id, "_", barcodes_rm[[id]])
})
cells_to_remove <- as.character(unlist(cells_to_remove))

length(colnames(seu_integrated))
seu_integrated$rm <- 'singlet'
seu_integrated$rm[match(cells_to_remove, rownames(seu_integrated@meta.data))] <- 'doublet'
Idents(seu_integrated) <- 'rm'
seu_integrated <- subset(seu_integrated, idents='singlet')

# saveRDS(seu_integrated, file = file.path(outdir, "seurat_obj", "scIntegrated.rds"))

############################
#### 2.a QC - ATAC data ####
dir.create(file.path(outdir, "qc"), showWarnings = F)
# seu_test<- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
if(!exists("seu_integrated")) seu_integrated<- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated.rds"))
checkObj(seu_integrated)

# Mononucleosome signal and TSS enrichment
DefaultAssay(seu_integrated) <- "ATAC"
seu_integrated <- NucleosomeSignal(object=seu_integrated, verbose=T) # strength of the nucleosome signal per cell [Ratio of (mononucleosome)/(nucleosome-free) (147-294bp fragments/ <147 bp)]
seu_integrated <- TSSEnrichment(seu_integrated, fast=T)

# add blacklist ratio and fraction of reads in peaks
seu_integrated$pct_reads_in_peaks <- seu_integrated$atac_peak_region_fragments / seu_integrated$atac_fragments * 100
if(count_blacklist_frags) seu_integrated$blacklist_ratio <- seu_integrated$blacklist_region_fragments / seu_integrated$peak_region_fragments
seu_integrated$high.tss <- ifelse(seu_integrated$TSS.enrichment > 2, 'High', 'Low')


## ATAC Viz
Idents(seu_integrated) <- seu_integrated$orig.ident
vp <- VlnPlot(object = seu_integrated,
              features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal",
                           'atac_peak_region_fragments'),
              ncol = 3, pt.size = 0, combine = F) 
pdf(file.path(outdir, "qc", "vlnplot.pdf"))
plot_grid(plotlist=lapply(vp, function(i) i + theme(text = element_text(size = 7),
                                                    axis.title.x=element_blank())), 
          ncol=2)
dev.off()

# tssp <- TSSPlot(seu_integrated, group.by = 'high.tss') + NoLegend()
# pdf(file.path(outdir, "qc", "tssenrich.pdf"))
# tssp
# dev.off()

# saveRDS(seu_integrated, file = file.path(outdir, "seurat_obj", "scIntegrated.rds"))

###########################
#### 2.b QC - RNA data ####
DefaultAssay(seu_integrated) <- "RNA"
seu_integrated <- PercentageFeatureSet(seu_integrated, pattern = "^MT-", col.name = "percent.mt")
# Plot percent.mt per sample
pdf(file.path(outdir, "qc", "percent_mt.pdf"))
VlnPlot(object = seu_integrated, features = 'percent.mt', split.by = 'orig.ident')
dev.off()

pdf(file.path(outdir, "qc", "featureScatter.pdf"), width = 14, height = 7)
plot1 <- FeatureScatter(seu_integrated, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt", shuffle = TRUE)
plot2 <- FeatureScatter(seu_integrated, feature1 = "nCount_ATAC", 
                        feature2 = "nCount_RNA", shuffle = TRUE)
plot_grid(plot1, plot2, nrow=1)
dev.off()

## https://www.nature.com/articles/s41591-021-01323-8#Sec15
# All cells expressing <200 or >6,000 genes were removed
# cells that contained <400 unique molecular identifiers (UMIs) 
# >15% mitochondrial counts


# Filter out high and low nCount ATAC data, or low TSS signal cells
## I think this may be causing a downstream issue, unsure why, but removing for now...
unfilt_tbl <- table(seu_integrated$orig.ident)
seu_integ <- subset(x = seu_integrated,
                    subset = nCount_ATAC < quantile(seu_integrated$nCount_ATAC, 0.99) &
                      nCount_ATAC > 1000 &
                      nucleosome_signal < 2 &
                      TSS.enrichment > 1)
atac_filt_tbl <- table(seu_integ$orig.ident)

seu_integ <- subset(seu_integ,
                    subset = nFeature_RNA > 200 & 
                      nCount_RNA > 400 &
                      nFeature_RNA < 6000 & 
                      percent.mt < 20)
rna_filt_tbl <- table(seu_integ$orig.ident)

# Create a count of the filtered cells
filtered_tbl <- t(rbind(unfilt_tbl, atac_filt_tbl, rna_filt_tbl))
write.table(filtered_tbl, file.path(outdir, "qc", "filtered_cell_cnt.tsv"), 
            sep="\t", quote = F, row.names = T, col.names = T)

rm(seu_integrated)
saveRDS(seu_integ, file = file.path(outdir, "seurat_obj", "scIntegrated_filtered.rds"))

########################
#### 3.Peak Calling ####
dir.create(file.path(outdir, "macs2"), recursive = T, showWarnings = F)
if(!exists("seu_integ")) seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_filtered.rds"))
DefaultAssay(seu_integ) <- 'ATAC'

# call peaks using MACS2
if(!file.exists(file.path(outdir, "macs2", "peaks_cell.rds"))){
  peaks_cells <- CallPeaks(seu_integ, macs2.path = macs2_path,
                           outdir=file.path(outdir, "macs2"), 
                           name='peaks_cells', cleanup=F)
  saveRDS(peaks_cells, file=file.path(outdir, "macs2", "peaks_cell.rds"))
} else {
  peaks_cells <- readRDS(file.path(outdir, "macs2", "peaks_cell.rds"))
}

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks_cells, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE) #ENCODE blacklist

# quantify counts in each peak
if(!file.exists(file.path(outdir, "macs2", "feat_mat.rds"))){
  macs2_counts <- FeatureMatrix(fragments = Fragments(seu_integ),
                                features = peaks,
                                cells = colnames(seu_integ))
  saveRDS(macs2_counts, file=file.path(outdir, "macs2", "feat_mat.rds"))
} else {
  macs2_counts <- readRDS(file.path(outdir, "macs2", "feat_mat.rds"))
}

# create a new assay using the MACS2 peak set and add it to the Seurat object
seu_integ[['peaks']] <- CreateChromatinAssay(counts = macs2_counts,
                                             fragments = Fragments(seu_integ),
                                             annotation = Annotation(seu_integ))
saveRDS(seu_integ, file=file.path(outdir, "seurat_obj", "scIntegrated_peaks.rds"))
#############################################
#### 4.a Preprocessing the counts - ATAC ####
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(outdir, "seurat_obj", "scIntegrated_peaks.rds"))
checkObj(seu_integ)

# Preprocessing of the ATAC peaks
DefaultAssay(seu_integ) <- "peaks"
seu_integ <- FindTopFeatures(seu_integ, min.cutoff = 5)
seu_integ <- RunTFIDF(seu_integ)
seu_integ <- RunSVD(seu_integ)

DefaultAssay(seu_integ) <- "ATAC"
seu_integ <- FindTopFeatures(seu_integ, min.cutoff = 5)
seu_integ <- RunTFIDF(seu_integ)
seu_integ <- RunSVD(seu_integ)

saveRDS(seu_integ, file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))


############################################
#### 4.b Preprocessing the counts - RNA ####
# Preprocessing of the RNA counts
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))

DefaultAssay(seu_integ) <- "RNA"
# https://github.com/satijalab/seurat/issues/1739
seu_integ <- SCTransform(seu_integ, assay = 'RNA', new.assay.name = 'SCT',
                         vars.to.regress = c('percent.mt'), conserve.memory = T) #nCount_RNA regressed in vst
seu_integ <- CellCycleScoring(seu_integ, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,
                        assay='SCT', set.ident = TRUE)
seu_integ$CC.Difference <- seu_integ$S.Score - seu_integ$G2M.Score
seu_integ <- SCTransform(seu_integ, assay = 'RNA', new.assay.name = 'SCT',
                         vars.to.regress = c('percent.mt', 'CC.Difference'),
                         conserve.memory = T)

## Find Neighbours and Cluster with HARMONY
dir.create(file.path(outdir, "harmony"), showWarnings = FALSE)
pdf(file.path(outdir, "harmony", "harmony_convergence.pdf"))
seu_integ <- RunPCA(seu_integ, verbose = FALSE)
seu_integ <- RunHarmony(seu_integ, c('orig.ident'), assay.use='SCT', plot_convergence = TRUE)
dev.off()

#Confirm #PC's determined explain > 95% of variance
stdev <- seu_integ@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.95))

seu_integ <- FindNeighbors(object = seu_integ, dims = 1:if(PCNum>30) 30 else PCNum, reduction ="harmony")
seu_integ <- FindClusters(object = seu_integ, resolution = 1.2, reduction ="harmony")
# Tweaked the UMAP parameters here
seu_integ <- RunUMAP(object = seu_integ, dims = 1:if(PCNum>30) 30 else PCNum, reduction = "harmony",
               n.neighbors = 30L, n.components = 2L, n.epochs=200L, min.dist=0.3)

# build a joint neighbor graph using both assays
DefaultAssay(seu_integ) <- 'SCT'
seu_integ <- FindMultiModalNeighbors(object = seu_integ,
                                     reduction.list = list("pca", "lsi"), 
                                     dims.list = list(1:30, 2:40),
                                     modality.weight.name = "RNA.weight",
                                     verbose = TRUE)

# build a joint UMAP visualization
seu_integ <- RunUMAP(object = seu_integ,
                     nn.name = "weighted.nn",
                     assay = "SCT", reduction = "harmony",
                     verbose = TRUE)

saveRDS(seu_integ, file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))

#####################################################
#### 5.a Cell type annotation - ENCODE/Blueprint ####
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))

dir.create(file.path(outdir, "annotation"))
getBlueprint <- TRUE
if(getBlueprint){
  #ref <- BlueprintEncodeData()
  #saveRDS(ref, file="~/downloads/BlueprintEncodeData.rds")
  bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/BlueprintEncodeData.rds") 
  for(lbl in c('label.fine', 'label.main')){
    for(clus in c(TRUE, FALSE)){
      rds_file <- paste0("celldex_blueprint.", gsub("label.", "", lbl),
                         ".", if(clus) 'cluster' else 'cell', ".rds")
      id <- gsub("^.*blueprint(.*).rds", "bp\\1", rds_file)
      if(!file.exists(file.path(outdir, "annotation", rds_file))){
        print(rds_file)
        blueprint_anno <- SingleR(test=GetAssayData(seu_integ), ref=bed.se, 
                            assay.type.test=1, labels=bed.se[[lbl]],
                            clusters=if(clus) seu_integ$seurat_clusters else NULL)
        saveRDS(blueprint_anno, file=file.path(outdir, "annotation", rds_file))
      } else {
        blueprint_anno <- readRDS(file.path(outdir, "annotation", rds_file))
      }
      
      if(clus){
        cluster_ids <- setNames(blueprint_anno$labels, 
                                as.character(rownames(blueprint_anno)))
        seu_integ@meta.data[,id] <- cluster_ids[as.character(seu_integ$seurat_clusters)]
      } else {
        seu_integ@meta.data[,id] <- blueprint_anno$labels
      }
    }
  }
} 


saveRDS(seu_integ, file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))

###########################################
#### 5.b Cell type annotation - PBMC3k ####
# load PBMC reference
if(!exists("seu_integ")) seu_integ <- readRDS(file=file.path(outdir, "seurat_obj", "scIntegrated_preproc.rds"))
checkObj(seu_integ)

reference <- SeuratDisk::LoadH5Seurat(pbmc_annot)
DefaultAssay(seu_integ) <- "SCT"

# transfer cell type labels from reference to query
tr_anchors <- FindTransferAnchors(reference = reference,
                                  query = seu_integ,
                                  normalization.method = "SCT",
                                  reference.reduction = "spca",
                                  recompute.residuals = FALSE,
                                  dims = 1:50)

preds <- TransferData(anchorset = tr_anchors, 
                      refdata = reference$celltype.l2,
                      weight.reduction = seu_integ[['pca']],
                      dims = 1:50)

seu_integ <- AddMetaData(object = seu_integ,
                         metadata = preds)
saveRDS(seu_integ, file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))

#########################################################
#### 5.c Finding all markers for each seurat-cluster ####
if(!exists("seu_integ"))  seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
dir.create(file.path(outdir, "markers"), showWarnings = F)

markers_f <- file.path(outdir, "markers", "seurat_markers.rds")
if(!file.exists(markers_f)){
  Idents(seu_integ) <- 'seurat_clusters'
  markers <- FindAllMarkers(seu_integ)
  saveRDS(markers, file=markers_f)
} else {
  markers <- readRDS(markers_f)
}
markers_cl <- split(markers, f=markers$cluster)

#############################################################
#### 5.d Building/Plotting cluster tree Dimplot [MANUAL] ####
if(!exists("seu_integ"))  seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))

seu_integ <- BuildClusterTree(object = seu_integ)
# Breaking apart and customizing the PlotClusterTree function
getMajorityAnnoFromCells <- function(meta, cluster='seurat_clusters', anno=NULL){
  meta_cl <- split(meta, f=meta[,cluster])
  cluster_prop <- lapply(meta_cl, function(mcl_x) {
    sort(round(table(mcl_x[,anno]) / nrow(mcl_x),2), decreasing = T)
  })
  setNames(
    paste0(sapply(cluster_prop, function(i) names(i)[1]), " [", 
           sapply(cluster_prop, function(i) i[1]), "]"),
    names(meta_cl)
  )
}
anno_cls <- list('pbmc'='predicted.id', 
                 'blueprint_fine'='bp.fine.cell',
                 'blueprint_main'='bp.main.cell')
pdf(file.path(outdir, "dimplot", "clustertree.pdf"))
cl_map <- lapply(anno_cls, function(acl){
  clkey <- getMajorityAnnoFromCells(meta=seu_integ@meta.data, 
                           cluster='seurat_clusters',
                           anno=acl)
  data.tree <- Tool(object = seu_integ, slot = "BuildClusterTree")
  data.tree$tip.label <- paste0(data.tree$tip.label, " - ", clkey[data.tree$tip.label])
  ape::plot.phylo(x = data.tree, direction = 'rightwards', cex=0.75)
  return(clkey)
})
dev.off()

## Reduce the annotated clusters to a dataframe of anno + proportion for printing 
cl_map_l <- lapply(cl_map, function(i){
 x <- strsplit(i, split=" \\[")
 x <- as.data.frame(do.call(rbind, x))
 x[,2] <- as.numeric(gsub("]$", "", x[,2]))
 x
})
cl_map_df <- as.data.frame((do.call(cbind, cl_map_l)))
write.table(cl_map_df, file=file.path(outdir, "dimplot", "clust_map.tsv"),
            sep="\t", col.names=T, row.names = F, quote=F)
#########################################################################
#### 5.e Aggregating clusters based on manual revision of annotation ####
if(!exists("seu_integ"))  seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
old_clusters <- seu_integ$seurat_clusters

res_dimplot_f <- file.path(outdir, "dimplot", "resolution_dimplot.pdf")
if(!file.exists(res_dimplot_f)){
  res_dps <- lapply(seq(0.1, 1.5, by=0.1), function(resolution){
    seu_integ <- FindClusters(object = seu_integ, resolution = resolution, reduction ="harmony")
    dp_clus <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
                       group.by='seurat_clusters', pt.size=0.5, shuffle=TRUE)
    return(dp_clus)
  })
  pdf(res_dimplot_f)
  lapply(res_dps, print)
  dev.off()
} 

res <- 0.3
seu_integ <- FindClusters(object = seu_integ, resolution = res, reduction ="harmony")
bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/BlueprintEncodeData.rds") 
for(lbl in c('label.fine', 'label.main')){
  rds_file <- paste0("res", res,
                     ".celldex_blueprint.", gsub("label.", "", lbl),
                     ".", 'cluster', ".rds")
  id <- gsub("^.*blueprint(.*).rds", "bp\\1", rds_file)
  if(!file.exists(file.path(outdir, "annotation", rds_file))){
    print(rds_file)
    blueprint_anno <- SingleR(test=GetAssayData(seu_integ), ref=bed.se, 
                              assay.type.test=1, labels=bed.se[[lbl]],
                              clusters=seu_integ$seurat_clusters)
    saveRDS(blueprint_anno, file=file.path(outdir, "annotation", rds_file))
  } else {
    blueprint_anno <- readRDS(file.path(outdir, "annotation", rds_file))
  }
  
  cluster_ids <- setNames(blueprint_anno$labels, 
                          as.character(rownames(blueprint_anno)))
  seu_integ@meta.data[,id] <- cluster_ids[as.character(seu_integ$seurat_clusters)]
}

seu_integ$old_clusters <- old_clusters
saveRDS(seu_integ, file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))

################################
#### 6.a Multimodal Dimplot ####
if(!exists("seu_integ"))  seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))

dir.create(file.path(outdir, "dimplot"), showWarnings = F)

DefaultAssay(seu_integ) <- 'SCT'
called_doublets <- grepl("^DF.classifications", colnames(seu_integ@meta.data))
if(called_doublets){
  doublet_col <- grep("^DF.classifications", colnames(seu_integ@meta.data))
  dp_doublet <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
                        group.by=doublet_col, pt.size=0.5, shuffle=TRUE)
}
dp_anno <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
               group.by='bp.fine.cell', pt.size=0.5, shuffle=TRUE)
dp_annocl <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='bp.fine.cluster', pt.size=0.5, shuffle=TRUE)
dp_pred <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
                   group.by='predicted.id', pt.size=0.5, shuffle=TRUE)
dp_clus <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
               group.by='seurat_clusters', pt.size=0.5, shuffle=TRUE)
dp_lbl <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
               group.by='orig.ident', pt.size=0.5, shuffle=TRUE)
dp_anno_lbl <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by='bp.fine.cluster', split.by='orig.ident', pt.size=0.5, shuffle=TRUE)
dp_anno_grp <- DimPlot(seu_integ, label = TRUE, repel = TRUE, reduction = "umap",
                  group.by='bp.fine.cluster', split.by='group', pt.size=0.5, shuffle=TRUE)



pdf(file.path(outdir, "dimplot", "dimplot.pdf"), width=17)
dp_lbl + dp_clus
dp_anno + dp_annocl
dp_anno + dp_pred
plot_grid(dp_anno_lbl, ncol=1)
plot_grid(dp_anno_grp, ncol=1)
if(called_doublets) dp_clus + dp_doublet
dev.off()

################################################################
#### OUTDATED - 7.a) DEG on Markers identifying INF for each cell type ####
library(org.Hs.eg.db)
library(msigdbr)

# Load in msigdb gene sets
m_df <- msigdbr(species = "Homo sapiens")
msig_ds <- list("immune"=list("C7", NULL),
                "pathway"=list("C2", "CP:REACTOME"),
                "gobp"=list("C5", "GO:BP"),
                "gomf"=list("C5", "GO:MF"))
msig_gsl <- lapply(msig_ds, function(mds){
  msigdbr(species = "Homo sapiens", category = mds[[1]], subcategory = mds[[2]]) %>% 
    dplyr::select(gs_name, entrez_gene)
})

# Create mapping from symbol to NCBI
dir.create(file.path(outdir, "differential"), showWarnings = F)
genome_gse <- org.Hs.eg.db
txby <- keys(genome_gse, 'SYMBOL')
gene_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                   keytype='SYMBOL', multiVals="first")
txby <- keys(genome_gse, 'ENTREZID')
entrez_ids <- mapIds(genome_gse, keys=txby, column='SYMBOL',
                   keytype='ENTREZID', multiVals="first")


# Massage the seurat object
# if(!exists("seu_integ")) seu_integ <- SeuratDisk::LoadH5Seurat(file.path(outdir, "seurat_obj", "scIntegrated_preproc.h5seurat"))
if(!exists("seu_integ")) seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
anno_type <- 'bp.fine.cell'
Idents(seu_integ) <- 'orig.ident'
seu_inf <- subset(x = seu_integ, idents = c('INFB', 'unstimulated'))
Idents(seu_inf) <- anno_type
DefaultAssay(seu_inf) <- 'SCT'
anno_type <- 'bp.fine.cluster'
cell_types <- unique(sort(seu_inf@meta.data[[anno_type]]))
seu_inf$orig.ident <- seu_inf$old.ident

# Find DEGs between INFB and Unstimulated
degs <- lapply(cell_types, function(ct){
  print(paste0(ct, "..."))
  cnts <- split(seu_inf@meta.data[[anno_type]] == ct, unique(as.character(seu_inf$orig.ident)))
  
  if(!all(sapply(cnts, sum) > 50)){
    return(NULL)
  }
  
  markers <- FindMarkers(seu_inf, ident.1 = 'INFB', 
                         group.by = 'orig.ident', subset.ident=ct)
  markers0 <- markers[which(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > log2(2)),]
  markers_spl <- split(markers0, f=markers0$avg_log2FC > 0)
  dirnm <- setNames(c('UP', 'DN'), c('TRUE', 'FALSE'))
  names(markers_spl) <- dirnm[names(markers_spl)]
  markers_spl
})
names(degs) <- cell_types

# Form the IFN-based geneset per celltype for GSEA
deg_gsl <- lapply(degs, function(deg){
  deg1 <- lapply(names(deg), function(deg_dir){
    deg0 <- deg[[deg_dir]]
    tryCatch({
      data.frame("gs_name"=rep(paste0("INFB-vs-unstimulated_", deg_dir), nrow(deg0)),
                 "entrez_gene"=as.character(gene_ids[rownames(deg0)]))
    }, error=function(e){NULL})
  })
  do.call(rbind, deg1)
})
names(deg_gsl)
dir.create(file.path(outdir, "infb-unstimulated"), recursive = F, showWarnings = F)
saveRDS(deg_gsl, file=file.path(outdir, "infb-unstimulated", "degs.rds"))

################################################################
#### OUTDATED - 7.a) DEG on Markers identifying INF for each cell type ####
if(!exists("deg_gsl")) deg_gsl <- readRDS(file=file.path(outdir, "infb-unstimulated", "degs.rds"))
if(!exists("seu_integ")) seu_integ <- readRDS(file = file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
anno_type <- 'bp.fine.cluster'
Idents(seu_integ) <- 'old.ident'

seu_inf <- subset(x = seu_integ, idents = c('CCP149', 'CCP228', 'CCP033', 'CCP246'))
Idents(seu_inf) <- anno_type
DefaultAssay(seu_inf) <- 'SCT'
cell_types <- names(deg_gsl)
seu_inf$orig.ident <- seu_inf$old.ident

## Calculate the DEGs and run GSEA on the pathways that define them
degs <- lapply(cell_types, function(ct){
  print(paste0(ct, "..."))
  cnts <- split(seu_inf@meta.data[[anno_type]] == ct, as.character(seu_inf$orig.ident))
  if(!all(sapply(cnts, sum) > 50)){
    return(NULL)
  }
  
  markers <- FindMarkers(seu_inf, ident.1 = 'lowISS', 
                         group.by = 'group', subset.ident=ct,
                         features=na.omit(entrez_ids[deg_gsl[[ct]]$entrez_gene]),
                         logfc.threshold=0, min.pct=0)
  markers_universe <- FindMarkers(seu_inf, ident.1 = 'lowISS', 
                                  group.by = 'group', subset.ident=ct,
                                  logfc.threshold=0)
  markers_merge <- rbind(markers, markers_universe[-which(rownames(markers_universe) %in% rownames(markers)),])
  
  # Make gene list
  gene_list <- markers_merge$avg_log2FC
  names(gene_list) <- gene_ids[rownames(markers_merge)]
  gene_list<-sort(na.omit(gene_list), decreasing = TRUE)
  msig <- tryCatch({
      GSEA(gene_list, TERM2GENE = deg_gsl[[ct]], 
           pvalueCutoff = 1, minGSSize=1, maxGSSize=2000)
    }, error=function(e){NULL})

  return(msig)
})
names(degs) <- cell_types
deg_df <- do.call(rbind, lapply(degs, function(i) tryCatch({i@result[,-c(2, 8, 11)]}, error=function(e){NULL})))
saveRDS(degs, file=file.path(outdir, "infb-unstimulated", "deg_infb-unstim.rds"))

degs <- readRDS(file.path(outdir, "infb-unstimulated", "deg_infb-unstim.rds"))
deg_df <- do.call(rbind, lapply(degs, function(i) tryCatch({i@result[,-c(2, 8, 11)]}, error=function(e){NULL})))


## RPS Genes


# Reduce GSEA result to a ggplot2 friendly format
godeg <- lapply(names(degs), function(id){
  i <- degs[[id]]
  tryCatch({
    x <- i$go@result
    x$cell_type <- id
    x
  }, error=function(e){NULL})
})
godeg <- do.call(rbind, godeg)
godeg_melt <- godeg[,c('Description', 'NES', 'qvalues', 'cell_type')]
godeg_melt$direction <- godeg$NES > 0
n <- 10
desc_oi <- lapply(split(godeg_melt, godeg_melt$direction), function(i){
  head(i$Description[order(sapply(split(i$NES, i$Description), mean))], n)
})
godeg_melt$Description <- factor(godeg_melt$Description, levels=unique(unlist(desc_oi)))
godeg_melt <- godeg_melt[-which(is.na(godeg_melt$Description)),]
godeg_melt$cell_type <- factor(godeg_melt$cell_type, levels=names(degs))

# Do the plotties
pdf(file.path(outdir, "differential", "deg_per-cell.pdf"), width = 25)
ggplot(data=godeg_melt, aes(x=Description, y=NES, fill=qvalues)) +
  geom_bar(stat='identity') + 
  facet_grid(cols=vars(cell_type),
             scales = "free_y", drop=F) +
  geom_hline(yintercept=0, size=1) +
  coord_flip() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

###############################################
#### 7.a) Add Motif info the seurat object #### 
dir.create(file.path(outdir, "differential"), showWarnings = F)
rerun_analysis <- T

# Massage the seurat object
if(!exists("seu_integ")) seu_integ <- readRDS(file.path(outdir, "seurat_obj", "scIntegrated_anno.rds"))
anno_type <- 'bp.fine.cluster'
Idents(seu_integ) <- anno_type
DefaultAssay(seu_integ) <- 'peaks'
cell_types <- unique(sort(seu_integ@meta.data[[anno_type]]))

# Add motif information from JASPAR2020
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", species = 'Homo sapiens'))
seu_integ <- AddMotifs(seu_integ, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)


###########################
#### 7.b) Run ChromVAR ####
## Run ChromVAR
chromvar_obj_f <- file.path(outdir, "seurat_obj", "scAtacMotif.rds")
if(!file.exists(chromvar_obj_f)){
  library(BiocParallel)
  BiocParallel::register(MulticoreParam(workers = 1))
  seu_integ <- RunChromVAR(
    object = seu_integ,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    new.assay.name = "chromvar"
  )
  saveRDS(seu_integ, file = chromvar_obj_f)
} else {
  seu_integ <- readRDS(chromvar_obj_f)
}

###########################################################
#### 7.c) Identifying Differential Peaks per cell type ####
chromvar_obj_f <- file.path(outdir, "seurat_obj", "scAtacMotif.rds")
if(!exists("seu_integ")) seu_integ <- readRDS(chromvar_obj_f)
rerun_analysis <- FALSE
anno_type <- 'bp.fine.cluster'
cell_types <- unique(sort(seu_integ@meta.data[[anno_type]]))

## Calculate Differential Accessible regions
#    search for DNA motifs that are overrepresented in a set of peaks that are 
#    differentially accessible between cell types
group_by <- 'group'
grp_smps <- unique(seu_integ@meta.data[,c('group', 'dataset')])
# list("ID" = [sampes-to-use], [column-to-compare], [group-to-compare])
# e.g. list("testA" = c("sampleA", "sampleB"), "treatment", "post") 
#  will compare sampleA (post) to sampleB (treat), based on the "treatment" category
comp_grps <- list("high-low"=list(idents=c(grp_smps[grp_smps$group=='highISS',]$dataset,
                                           grp_smps[grp_smps$group=='lowISS',]$dataset),
                                  'group', 'lowISS'),
                  "high-healthy"=list(idents=c(grp_smps[grp_smps$group=='highISS',]$dataset,
                                               grp_smps[grp_smps$group=='healthy',]$dataset),
                                      'group', 'healthy'),
                  "low-healthy"=list(idents=c(grp_smps[grp_smps$group=='lowISS',]$dataset,
                                              grp_smps[grp_smps$group=='healthy',]$dataset),
                                     'group', 'healthy'),
                  "highlow-healthy"=list(idents=c(grp_smps[grp_smps$group=='lowISS',]$dataset,
                                                  grp_smps[grp_smps$group=='highISS',]$dataset,
                                                  grp_smps[grp_smps$group=='healthy',]$dataset),
                                         'group', 'healthy'))

dpeaks_f <- file.path(outdir, "differential", "da_peaks.rds")
if((!file.exists(dpeaks_f)) | rerun_analysis){
  # Iterate through the group-comparisons given
  das_grps <- lapply(names(comp_grps), function(grp_id){
    print(paste(cg$idents, collapse=","))
    Idents(seu_integ) <- 'dataset' # Extract the highISS and lowTSS
    seu_grp <- subset(seu_integ, idents = cg$idents)
    Idents(seu_grp) <- anno_type

    # Call differential peaks per cell-type by the groups given
    das <- lapply(cell_types, function(ct){
      print(paste0(ct, "..."))
      
      cnts <- split(seu_grp@meta.data[[anno_type]] == ct, 
                    seu_grp@meta.data[[group_by]])
      if(!all(sapply(cnts, sum) > 50)){
        return(NULL)
      }
      
      # Differential accessibility test
      all_grps <- unique(seu_grp@meta.data[,cg[[2]]])
      grps_given <- cg[[3]]
      grps_inverse <- all_grps[!all_grps %in% cg[[3]]]
      da_peaks <- FindMarkers(object = seu_grp,
                              only.pos = F,
                              ident.1 = cg[[3]], 
                              group.by = cg[[2]], 
                              subset.ident=ct,
                              test.use = "LR", 
                              min.pct = 0.05,
                              latent.vars = 'nCount_ATAC')

      ## Split differential peaks into two lists, one for 
      # more-accessible in each group
      da_peaks_spl <- split(da_peaks, f=da_peaks$avg_log2FC >0)
      grps_map <- setNames(c(paste(grps_inverse, collapse="_"),
                             paste(grps_given, collapse="_")),
                           c('FALSE', 'TRUE'))
      names(da_peaks_spl) <- grps_map[names(da_peaks_spl)]
      return(da_peaks_spl)
    })
    names(das) <- cell_types
    return(das)
  })
  saveRDS(das_grps, file=dpeaks_f)
} else {
  das_grps <- readRDS(file=dpeaks_f)
}

#######################################################################
#### 7.d) Using the differential peaks to identify enriched motifs ####
dpeaks_f <- file.path(outdir, "differential", "da_peaks.rds")
if(!exists("seu_integ")) seu_integ <- readRDS(chromvar_obj_f)
if(!exists("das_grps")) das_grps <- readRDS(dpeaks_f)

dmotif_f <- file.path(outdir, "differential", "da_motifs.rds")
if((!file.exists(dmotif_f)) | rerun_analysis){
  # Iterate through groups
  dmotif_grps <- lapply(names(comp_grps), function(grp_id){
    cg <- comp_grps[[grp_id]]
    Idents(seu_integ) <- 'dataset' # Extract the highISS and lowTSS
    seu_grp <- subset(seu_integ, idents = cg$idents)
    Idents(seu_grp) <- anno_type
    
    # Iterate through cell types
    das <- lapply(cell_types, function(ct){
      print(paste0(ct, "..."))
      
      cnts <- split(seu_grp@meta.data[[anno_type]] == ct, 
                    seu_grp@meta.data[[group_by]])
      if(!all(sapply(cnts, sum) > 50)){
        return(NULL)
      }
      da_peaks_spl <- das_grps[[grp_id]][[ct]]
      
      ## Calculate the enriched motifs within the more-accessible peaks
      das_l <- lapply(names(da_peaks_spl), function(grp){
        da_peaks_i <- da_peaks_spl[[grp]]
        da_peaks_id <- rownames(da_peaks_i[which(da_peaks_i$p_val_adj < 0.05),])
        
        #isolating all accessible peaks in cell type: ct
        open.peaks <- AccessiblePeaks(seu_grp, idents = ct)
        
        # match the overall GC content in the peak set
        meta.feature <- GetAssayData(seu_grp, assay = "peaks", slot = "meta.features")
        enriched.motifs <- tryCatch({
          peaks.matched <- MatchRegionStats(meta.feature = meta.feature[open.peaks, ],
                                            query.feature = meta.feature[da_peaks_id, ],
                                            n = 50000)
          ## Matching GC.percent distribution to identify enrichment 
          FindMotifs(object = seu_grp,
                     features = da_peaks_id,
                     background=peaks.matched)
        }, error=function(e){NULL})
        
        return(list("posAccGrp"=grp, "da"=da_peaks_i, "motif"=enriched.motifs))
      })
      names(das_l) <- names(da_peaks_spl)
      return(das_l)
    })
    names(das) <- cell_types
    return(das)
  })
  names(dmotif_grps) <- names(comp_grps)
  saveRDS(dmotif_grps, file=dmotif_f)
} else {
  dmotif_grps <- readRDS(dmotif_f)
}

## Visualization
# Reduce differential motifs to a ggplot2 friendly format
# and return the enriched motifs within the differential peaks
n <- 20 #top n-features to extract
head_ct <-TRUE # extract the top n-features from each cell type rather than overall
da_motifs <- lapply(names(dmotif_grps), function(das_id){
  das <- dmotif_grps[[das_id]]
  grps <- unlist(sapply(das, function(i) sapply(i, function(j) j$posAccGrp)))
  grp_dir <- setNames(c(-1, 1), unique(sort(grps)))
  da <- lapply(names(das), function(id){
    i <- das[[id]]
    # Aggregate the bi-directional fold-enrichment and assign directions
    i_motifs <- lapply(i, function(i_dir){
      i_motif <- i_dir$motif
      i_motif$fold.enrichment <-  grp_dir[i_dir$posAccGrp] * i_motif$fold.enrichment
      i_motif$posAccGrp <- i_dir$posAccGrp
      return(i_motif)
    })
    i_motifs <- do.call(rbind, i_motifs)
    if(!is.null(i_motifs)) i_motifs$cell_type <- id
    i_motifs
  })
  
  da <- as.data.frame(do.call(rbind, da))
  # da <- da[which(da$observed >= 3),]
  da$q <- p.adjust(da$pvalue, method='bonferroni')
  da_filtq <- da[which(da$q < 0.05),]
  da_filt <- da_filtq[,c('motif.name', 'observed', 
                        'fold.enrichment', 'cell_type', "q")]
  if(head_ct){
    da_filt_l <- lapply(split(da_filt, f=da_filt$cell_type), head, n=n)
    da_filt <- as.data.frame(do.call(rbind, da_filt_l))
  } else {
    da_filt <- head(da_filt[order(da_filt$q),], n)
  }
  
  #sorted fold enrichment per motif
  desc_oi <- sort(sapply(split(da_filt[,'fold.enrichment'], 
                               f=da_filt[,'motif.name']), 
                         function(i) mean(abs(i))),
                  decreasing = T)
  
  # Remove motifs that are not to be plotted
  da_filt[,'motif.name'] <- factor(da_filt[,'motif.name'], 
                                   levels=unique(names(desc_oi)))
  # da_filt <- da_filt[-which(is.na(da_filt[,'motif.name'])),]
  da_filt$cell_type <- factor(da_filt$cell_type, levels=names(das))
  da_filt$grp <- names(grp_dir)[as.integer(da_filt$fold.enrichment > 0)+1]
  
  # Do the plotties
  pdf(file.path(outdir, "differential", paste0(das_id, "_da.pdf")), height=17, width = 25)
  # pdf("~/xfer/test3.pdf", height=14, width=25)
  b <- c(-4, -2, 0, 2, 4)
  gp <-ggplot(data=da_filt, aes(x=motif.name, y=fold.enrichment, 
                                fill=grp)) +
    geom_bar(stat='identity') + 
    geom_point(aes(x=motif.name, y=fold.enrichment, size=observed)) +
    facet_grid(cols=vars(cell_type), drop=F) +
    geom_hline(yintercept=0, size=1) +
    coord_flip() + 
    theme_minimal() + 
    scale_y_continuous(breaks=b, labels=abs(b), limits = c(min(b), max(b))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(gp)
  dev.off()
  
  return(da_filtq)
})

################################################################
#### 7.e) Annotate the genes within the differential peaks  ####
dpeaks_f <- file.path(outdir, "differential", "da_peaks.rds")
if(!exists("seu_integ")) seu_integ <- readRDS(chromvar_obj_f)
if(!exists("das_grps")) das_grps <- readRDS(dpeaks_f)

# Make TxDB object
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

# Creating mapping of Gene Ensembl to SYMBOLs
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

dpeak_anno_f <- file.path(outdir, "differential", "da_peak_anno.rds")
if(!file.exists(dpeak_anno_f)){
  das_anno_grps <- lapply(das_grps, function(da_grp){
    lapply(da_grp, function(da_ct){
      lapply(da_ct, function(da_i){
        da_i$chr <- gsub("-.*$", "", rownames(da_i))
        da_i$start <- gsub("^.*-(.*)-.*$", "\\1", rownames(da_i))
        da_i$end <- gsub("^.*-", "", rownames(da_i))
        da_gr <- makeGRangesFromDataFrame(da_i, keep.extra.columns = T)
        seqlevelsStyle(da_gr) <- 'NCBI'
        
        # Annotate the genomic features of differential accessible regions
        anno <- annotatePeak(da_gr, TxDb=txdb, tssRegion = c(-3000, 3000))
        anno_gr <- anno@anno
        anno_gr$symbol <- gene_ids[anno_gr$geneId]
        anno_gr$symbol[which(is.na(anno_gr$symbol))] <- anno_gr$geneId[which(is.na(anno_gr$symbol))]
        return(anno_gr)
      })
    })
  })
  saveRDS(das_anno_grps, file=dpeak_anno_f)
} else {
  das_anno_grps <- readRDS(dpeak_anno_f)
}

txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
gr <- GRanges(BSgenome.Hsapiens.UCSC.hg38@seqinfo)
gr <- keepStandardChromosomes(gr, pruning.mode = 'coarse')
seqlevels(gr) <- gsub("chr", "", seqlevels(gr))

TxDb <- ChIPseeker:::loadTxDb(txdb)
features <- ChIPseeker:::getGene(txdb, by = "transcript")
anno <- ChIPseeker:::getGenomicAnnotation(features, distance=0, tssRegion = c(-3000, 3000), TxDb=txdb,
                                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
anno <- annotatePeak(gr, TxDb=txdb, 
                     tssRegion = c(-3000, 3000))

################################
#### 8. CD4-subset analysis ####
cd4obj_f <- file.path(outdir, "seurat_obj", "sc_CD4_subset.rds")
if(!exists("seu")) seu <- readRDS(cd4obj_f)
seu$phenotype_detailed <- gsub(" ", "_", seu$phenotype_detailed)
seu$phenotype_detailed2 <- seu$phenotype_detailed
seu$phenotype_detailed2[which(seu$phenotype_detailed2 %in% c("CD4_N_mem", "CD4_TEMRA"))] <- 'CD4'

cd4outdir <- file.path(outdir, "CD4")
dir.create(cd4outdir, showWarnings = F)

# Make TxDB file from GTF
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
TxDb <- ChIPseeker:::loadTxDb(txdb)
features <- ChIPseeker:::getGene(txdb, by="gene")

GR <- transcripts(txdb)
PR <- promoters(txdb, upstream=2000, downstream=400)

# Add motif information from JASPAR2020
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", species = 'Homo sapiens'))
# Creating mapping of Gene Ensembl to SYMBOLs
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")
# Get mapping of JASPAR motif id to name
motif_names <- GetMotifData(object = seu, assay = 'peaks', 
                            slot = "motif.names")

# Subset based on the cells (idents) and test for differences between groups (grps)
comp_grps <- list("cd4s"=list(idents=c('CD4_TEM', 'CD4_N_mem', 'CD4_TEMRA'),
                              idents.id=c('phenotype_detailed'),
                              grps=c('highISS', 'lowISS', 'healthy'),
                              grps.id=c('group')),
                  "cd4n_temra"=list(idents=c('CD4'),
                                    idents.id=c('phenotype_detailed2'),
                                    grps=c('highISS', 'lowISS', 'healthy'),
                                    grps.id=c('group')))

### a) Differential accessible regions (DAR) between Hi-Lo ----
dir.create(file.path(cd4outdir, "dar"), showWarnings = F)
cd4_peaks_f <- file.path(cd4outdir, "dar", "da_peaks.rds")
rerun_analysis <- FALSE

# Subset based on the cells (idents) and test for differences between groups (grps)
comp_grps <- list("cd4s"=list(idents=c('CD4_TEM', 'CD4_N_mem', 'CD4_TEMRA'),
                              idents.id=c('phenotype_detailed'),
                              grps=c('highISS', 'lowISS', 'healthy'),
                              grps.id=c('group')),
                  "cd4n_temra"=list(idents=c('CD4'),
                                    idents.id=c('phenotype_detailed2'),
                                    grps=c('highISS', 'lowISS', 'healthy'),
                                    grps.id=c('group')))

if((!file.exists(cd4_peaks_f)) | rerun_analysis){
  # Iterate through the group-comparisons given
  das_grps <- lapply(names(comp_grps), function(grp_id){
    cg <- comp_grps[[grp_id]]
    print(paste(cg$idents, collapse=","))
    
    # Subset the seurat object based ona set of cells
    Idents(seu) <- cg$idents.id 
    cg_das <- lapply(cg$idents, function(cg_ident){
      print(paste0(cg_ident, "..."))
      # Subset for a specific cell-type
      seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
      
      # Differential accessibility test between every group
      grp_combn <- combn(cg$grps, m=2, simplify = T)
      colnames(grp_combn) <- apply(grp_combn, 2, paste, collapse="_")
      
      das <- apply(grp_combn, 2, function(grp){
        da_peaks <- FindMarkers(object = seu_grp,
                                only.pos = F,
                                ident.1 = grp[1], 
                                ident.2 = grp[2],
                                group.by = cg$grps.id, 
                                test.use = "LR", 
                                min.pct = 0.05,
                                logfc.threshold = 0.1,
                                latent.vars = 'nCount_peaks')
        return(da_peaks)
      })
      return(das)
    })
    names(cg_das) <- cg$idents
    return(cg_das)
  })
  das_grps <- unlist(das_grps, recursive = F)
  saveRDS(das_grps, file=cd4_peaks_f)
} else {
  das_grps <- readRDS(file=cd4_peaks_f)
}

### b) Genes & Motifs found within DARs ----
dir.create(file.path(cd4outdir, "dar"), showWarnings = F)
cd4_peaks_f <- file.path(cd4outdir, "dar", "da_peaks.rds")
das_grps <- readRDS(file=cd4_peaks_f)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
rerun_analysis <- FALSE

if((!file.exists(cd4_peaks_anno_f)) | rerun_analysis){
  # Annotate the genes and motifs within each DA peak
  das_anno <- lapply(das_grps, function(das){
    print(names(das))
    peak_motif <- lapply(das, function(da_i){
      print(".")
      if(is.null(da_i)) return(NULL)
      
      da_i$chr <- gsub("-.*$", "", rownames(da_i))
      da_i$start <- gsub("^.*-(.*)-.*$", "\\1", rownames(da_i))
      da_i$end <- gsub("^.*-", "", rownames(da_i))
      da_gr <- makeGRangesFromDataFrame(da_i, keep.extra.columns = T)
      seqlevelsStyle(da_gr) <- 'NCBI'
      
      # Annotate the genomic features of differential accessible regions
      anno <- annotatePeak(da_gr, TxDb=txdb, tssRegion = c(-3000, 3000),
                           level='gene', addFlankGeneInfo = TRUE, flankDistance = 500,
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", 
                                                         "Intron", "Downstream", "Intergenic"))
      anno_gr <- anno@anno
      anno_gr$symbol <- gene_ids[anno_gr$geneId]
      anno_gr$symbol[which(is.na(anno_gr$symbol))] <- anno_gr$geneId[which(is.na(anno_gr$symbol))]
      anno_gr$symbol_fl <- gene_ids[anno_gr$flank_geneIds]
      anno_gr$symbol_fl[which(is.na(anno_gr$symbol_fl))] <- anno_gr$flank_geneIds[which(is.na(anno_gr$symbol_fl))]
      
      # Annotate motifs for each fragment
      seqlevelsStyle(da_gr) <- 'UCSC'
      motif_ix <- matchMotifs(pfm, da_gr, genome = BSgenome.Hsapiens.UCSC.hg38,
                              out = "positions", p.cutoff = 5e-05) 
      motif_sc <- matchMotifs(pfm, da_gr, genome = BSgenome.Hsapiens.UCSC.hg38,
                              out = "scores", p.cutoff = 5e-05) 
      
      return(list("gene"=anno_gr, "motif_pos"=motif_ix, "motif_score"=motif_sc))
    })
  })
  
  # Add in the DA data to the gene and motif annotation
  das_grps <- lapply(names(das_anno), function(das_id){
    da_x <- lapply(names(das_anno[[das_id]]), function(da){
      x <- das_anno[[das_id]][[da]]
      x[['da']] <- das_grps[[das_id]][[da]]
      return(x)
    })
    names(da_x) <- names(das_anno[[das_id]])
    return(da_x)
  })
  names(das_grps) <- names(das_anno)
  
  saveRDS(das_grps, file=cd4_peaks_anno_f)
} else {
  das_grps <- readRDS(file=cd4_peaks_anno_f)
}

### c) TF motif enrichment across all Peaks (ChromVAR) ----
dir.create(file.path(cd4outdir, "tf"), showWarnings = F)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
das_grps <- readRDS(file=cd4_peaks_anno_f)
cd4_tf_f <- file.path(cd4outdir, "tf", "dtf.rds")
rerun_analysis <- FALSE

if((!file.exists(cd4_tf_f)) | rerun_analysis){
  # Iterate through the group-comparisons given
  dtfs_grps <- lapply(names(comp_grps), function(grp_id){
    cg <- comp_grps[[grp_id]]
    print(paste(cg$idents, collapse=","))
    
    # Subset the seurat object based ona set of cells
    Idents(seu) <- cg$idents.id 
    cg_dtfs <- lapply(setNames(cg$idents,cg$idents), function(cg_ident){
      print(paste0(cg_ident, "..."))
      # Subset for a specific cell-type
      seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
      DefaultAssay(seu_grp) <- 'chromvar'
      
      # Differential accessibility test between every group
      grp_combn <- combn(cg$grps, m=2, simplify = T)
      colnames(grp_combn) <- apply(grp_combn, 2, paste, collapse="_")
      
      dtfs <- apply(grp_combn, 2, function(grp){
        # motif_x <- names(motif_names[motif_names=='STAT1'])
        # mx <- GetAssayData(seu_grp)[motif_x,which(seu_grp@meta.data$group == grp[1])]
        # my <- GetAssayData(seu_grp)[motif_x,which(seu_grp@meta.data$group == grp[2])]
        # rbind(summary(mx), summary(my))
        
        dtf <- FindMarkers(object = seu_grp,
                           slot = "data",
                           only.pos = F,
                           ident.1 = grp[1], 
                           ident.2 = grp[2],
                           group.by = cg$grps.id, 
                           test.use = "wilcox", 
                           logfc.threshold=0.01,
                           min.pct = 0.05)
        dtf$MA <- rownames(dtf)
        dtf$motif <- motif_names[dtf$MA]
        return(dtf)
      })
      return(dtfs)
    })
    return(cg_dtfs)
  })
  dtfs_grps <- unlist(dtfs_grps, recursive = F)
  saveRDS(dtfs_grps, file=cd4_tf_f)
} else {
  dtfs_grps <- readRDS(file=cd4_tf_f)
}



# Idents(seu) <- 'phenotype_detailed' 
# cg_ident <- 'CD4_TEM'
# seu_grp <- subset(seu, idents = cg_ident) 
# DefaultAssay(seu_grp) <- 'chromvar'
# 
# # Differential accessibility test between every group
# grp <- c('highISS', 'lowISS')
# dtf <- FindMarkers(object = seu_grp,
#                    slot = "data",
#                    only.pos = F,
#                    ident.1 = grp[1], 
#                    ident.2 = grp[2],
#                    group.by = 'group', 
#                    test.use = "wilcox", 
#                    logfc.threshold = 0.01,
#                    min.pct = 0.05)
# dtf$MA <- rownames(dtf)
# dtf$motif <- motif_names[dtf$MA]
# dtf[which(dtf$motif %in% c('REL', 'FOS::JUNB')),]



### d) ChromVAR DARs: TF motif enrichment within DARs ----
dir.create(file.path(cd4outdir, "tf"), showWarnings = F)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
das_grps <- readRDS(file=cd4_peaks_anno_f)
cd4_tf_f <- file.path(cd4outdir, "tf", "dtf.rds")
chromvar_dtf <- readRDS(file=cd4_tf_f)
cd4_tf_dar_f <- file.path(cd4outdir, "tf", "dtf_dar.rds")
rerun_analysis <- FALSE

cg_ident <- 'CD4_TEM'
res_ident <- 'highISS_lowISS'

da <- das_grps[[cg_ident]][[res_ident]]$da
da_sig_spl <- split(da, f=da$p_val < 0.05)
names(da_sig_spl) <- setNames(c('nonDAR', 'DAR'), 
                              c('FALSE', 'TRUE'))[as.character(names(da_sig_spl))]

features <- list(chromvar_DAR=rownames(da_sig_spl$DAR),
                 chromvar_nonDAR=rownames(da_sig_spl$nonDAR))

## Run the customized chromVar function for calling peaks
for(assay_id in names(features)){
  print(assay_id)
  feature_i <- features[[assay_id]]
  
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 1))
  DefaultAssay(seu) <- 'peaks'
  seu <- RunChromVar_subset(
    object = seu,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    features=feature_i,
    new.assay.name=assay_id
  )
}

## Find differential motifs from chromVAR scores
dar_dtfs <- lapply(setNames(names(features),names(features)), 
                   function(assay_id){
  seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
  DefaultAssay(seu_grp) <- assay_id
  
  # Differential accessibility test between every group
  grp_combn <- combn(cg$grps, m=2, simplify = T)
  colnames(grp_combn) <- apply(grp_combn, 2, paste, collapse="_")
  
  grp <- strsplit(res_ident, split="_")[[1]]
  dtf <- FindMarkers(object = seu_grp,
                     slot = "data",
                     only.pos = F,
                     ident.1 = grp[1], 
                     ident.2 = grp[2],
                     group.by = cg$grps.id, 
                     test.use = "wilcox", 
                     logfc.threshold=0.1,
                     min.pct = 0.05)
  dtf$MA <- rownames(dtf)
  dtf$motif <- motif_names[dtf$MA]
  return(dtf)
})
saveRDS(dar_dtfs, file=cd4_tf_dar_f)

dar_dtfs <- readRDS(file=cd4_tf_dar_f)



cg_ident <- 'CD4_TEM'
res_ident <- 'highISS_lowISS'
res_ids <- strsplit(res_ident, split = "_")[[1]]
cells_split <- split(colnames(seu), f=list(seu$group, seu$phenotype_detailed))
cells_grp1 <- cells_split[[paste0(res_ids[1], ".", cg_ident)]]
cells_grp2 <- cells_split[[paste0(res_ids[2], ".", cg_ident)]]

motif_ids <- data.frame(id=rownames(seu_assay),
                        name=as.character(motif_names[rownames(seu_assay)]))
motif_i <- motif_ids$name


m_chromvars <- lapply(chromvar_ids, function(chromvar_id){
  DefaultAssay(seu) <- chromvar_id #'chromvar'
  seu_assay <- GetAssayData(seu)
  
  m_chromvars <- lapply(motif_i, function(mi){
    motif_id_i <- motif_ids[motif_ids$name == mi,]$id
    if(length(motif_id_i)==0) return(NULL)
    tval <- t.test(x=seu_assay[motif_id_i, cells_grp1],
                   y=seu_assay[motif_id_i, cells_grp2])
    m_chromvar <- data.frame("stat"=tval$statistic,
                             "p"=tval$p.value,
                             "grp1"=mean(seu_assay[motif_id_i, cells_grp1]),
                             "grp2"=mean(seu_assay[motif_id_i, cells_grp2]),
                             "motif"=mi)
    # print(paste0(chromvar_id, " - ", mi, ": ", w_pval))
    return(m_chromvar)
  })
  m_chromvars_seu <- as.data.frame(do.call(rbind, m_chromvars))
  m_chromvars_seu$chromvar_obj <- chromvar_id
  m_chromvars_seu$pval_adj <- p.adjust(m_chromvars_seu$p, method='fdr')
  m_chromvars_seu <- m_chromvars_seu[order(m_chromvars_seu$pval_adj),]
  return(m_chromvars_seu)
})

df <- as.data.frame(do.call(rbind, m_chromvars))
df$pval_adj <- -1*log10(as.numeric(df$pval_adj))
sig_subset <- df %>% 
  filter(pval_adj > (-1*log10(0.1))) %>%
  group_by(chromvar_obj) %>%
  slice_max(order_by = pval_adj, n = 50)
sig_subset2 <- df %>%
  filter(pval_adj > (-1*log10(0.1))) %>%
  mutate(direction=stat>0) %>%
  group_by(direction) %>%
  slice_max(order_by = pval_adj, n = 25)

sig_motifs <- unique(c(sig_subset$motif))
sig_motifs <- unique(c(sig_subset$motif, sig_subset2$motif))
sig_df <- df[which(df$motif %in% sig_motifs),]
sig_df$sig <- sig_df$pval_adj > (-1*log10(0.1))
sig_df$motif <- factor(sig_df$motif, levels=sig_motifs)
sig_df$chromvar_obj <- factor(sig_df$chromvar_obj,
                              levels=paste0("chromvar", c("", "_DAR", "_nonDAR")))

pdf("~/xfer/test2.pdf", height = 12)
ggplot(sig_df, aes(y=motif, x=stat, fill=chromvar_obj, group=chromvar_obj)) +
  geom_bar(stat='identity', position='dodge') +
  geom_bar(aes(x=stat, y=motif, group=chromvar_obj, alpha=sig), 
           stat='identity', position='dodge', fill='grey') +
  scale_alpha_discrete(range = c(1, 0)) +
  xlab("t-statistic") + 
  theme_classic()
dev.off()

### e) Correlation between chromVAR DMs ---------------
cd4obj_chromvar_f <- file.path(outdir, "seurat_obj", "sc_chromvar_CD4_subset.rds")
dir.create(file.path(cd4outdir, "tf"), showWarnings = F)
cd4_peaks_anno_f <- file.path(cd4outdir, "dar", "da_peaks_anno.rds")
das_grps <- readRDS(file=cd4_peaks_anno_f)
cd4_tf_f <- file.path(cd4outdir, "tf", "dtf.rds")
chromvar_dtf <- readRDS(file=cd4_tf_f)
cd4_tf_dar_f <- file.path(cd4outdir, "tf", "dtf_dar.rds")
rerun_analysis <- FALSE

DefaultAssay(seu) <- 'peaks'
peaks_mat <- GetAssayData(seu)
motifs_mat <- Motifs(seu)@data
DefaultAssay(seu) <- 'SCT'
expr_mat <- GetAssayData(seu)

seu <- SyncExprAndPeaks(
  object=seu, 
  peaks_mat=peaks_mat, 
  motifs_mat=motifs_mat,
  expr_mat=expr_mat,
  txdb=makeTxDbFromGFF(file = gtf_file, format = "gtf"),
  new.assay.id='pr_expr_peaks'
)
  
BiocParallel::register(BiocParallel::MulticoreParam(workers = 1))
DefaultAssay(seu) <- 'pr_expr_peaks'
seu <- RunChromVAR_subset(
  object = seu,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix=f_motifs_mat,
  new.assay.name='chromvar_expr'
)

## ChromVAR Matrix [ motif x cells ] (chromvar_mat)
DefaultAssay(seu) <- 'chromvar'
peak_chromvar_mat <- GetAssayData(seu)

DefaultAssay(seu) <- 'chromvar_expr'
expr_chromvar_mat <- GetAssayData(seu)

cg_ident <- 'CD4_TEM'
res_ident <- 'highISS_lowISS'
res_ids <- strsplit(res_ident, split = "_")[[1]]
cells_split <- split(colnames(seu), f=list(seu$group, seu$phenotype_detailed))

motif_i <- rownames(peak_chromvar_mat)
m_chromvars <- lapply(res_ids, function(grp_id){
  cells_grp <- cells_split[[paste0(grp_id, ".", cg_ident)]]
  
  m_chromvars <- lapply(motif_i, function(mi){
    corval <- cor.test(x=peak_chromvar_mat[mi, cells_grp],
                       y=expr_chromvar_mat[mi, cells_grp], 
                       method='spearman')
    corval$estimate
    corval$p.value
    m_chromvar <- data.frame("r"=corval$estimate,
                             "p"=corval$p.value,
                             "peak_mean"=mean(peak_chromvar_mat[mi, cells_grp]),
                             "expr_mean"=mean(expr_chromvar_mat[mi, cells_grp]),
                             "motif"=as.character(motif_names[mi]))
    # print(paste0(chromvar_id, " - ", mi, ": ", w_pval))
    return(m_chromvar)
  })
  m_chromvars_seu <- as.data.frame(do.call(rbind, m_chromvars))
  m_chromvars_seu$group <- grp_id
  m_chromvars_seu$pval_adj <- p.adjust(m_chromvars_seu$p, method='fdr')
  m_chromvars_seu <- m_chromvars_seu[order(m_chromvars_seu$pval_adj),]
  
  return(m_chromvars_seu)
})
m_chromvar <- as.data.frame(do.call(rbind, m_chromvars))
m_chromvar$sig <- m_chromvar$pval_adj < 0.01
m_chromvar$log10q <- -1*log10(m_chromvar$pval_adj)

pdf("~/xfer/test3.pdf")
ggplot(m_chromvar, aes(x=r, y=log10q, color=group)) +
  geom_point() +
  geom_point(aes(x=r, y=log10q, alpha=sig), color='grey') +
  scale_alpha_discrete(range = c(1, 0)) +
  geom_text_repel(data=m_chromvar[which(m_chromvar$sig),], 
                  aes(x=r, y=log10q, label=motif), 
                  size=3, min.segment.length = 0, box.padding = 0.1) +
  theme_classic() +
  xlim(-0.25, 0.25) +
  xlab("Spearman correlation") + ylab("-log10(q)")
dev.off()
write.table(m_chromvar, file="~/xfer/test.csv",
            sep=",", quote=F, col.names = T, row.names = F)

saveRDS(seu, file=cd4obj_chromvar_f)



### f) Motif-specific TSS Plot between High and LowISS groups ----
motif_tags_f <- file.path(cd4outdir, "tf", "motifs_tssplot.rds")
rerun_analysis <- FALSE

if((!file.exists(motif_tags_f)) | rerun_analysis){
  DefaultAssay(seu) <- 'peaks'
  peaks_mat <- GetAssayData(seu)
  motifs_mat <- Motifs(seu)@data
  colnames(motifs_mat) <- motif_names[colnames(motifs_mat)]
  motif_tags <- getMotifTags(peaks_mat, motifs_mat, txdb, 
                             grps_col='group', grps=c('highISS', 'lowISS'),
                             return_gg=TRUE)
  saveRDS(motif_tags, file=motif_tags_f)
} else {
  motif_tags <- readRDS(motif_tags_f)
}


# pdf("~/xfer/test.pdf")
# lapply(motif_tags[c('STAT1', 'REL')], function(motif_gg){
#   plot_grid(plotlist = motif_gg, ncol=2)
# })
# dev.off()

### g) Motif-specific TSS Plot between High and LowISS groups ----
# Given a set of DEG [g], we can create a 
# boolean matrix of motifs M [genes x motifs] focused only on whether there is 
# a peaks for motif_i in the promoter of g_x. We can then do something like 
# a propprtions test to compare the ratio of genes with motif_i to 
# genes not with motif_i
deg_outpath <- file.path(cd4outdir, "deg_motif")
dir.create(deg_outpath, showWarnings = F)

degs_f <- file.path(deg_outpath, "degs.rds")
peaks_f <- file.path(deg_outpath, "peaks.rds")
out_stats_f <- file.path(deg_outpath, "stats_p.csv")
rerun_analysis <- FALSE

## Differential expression of CD4 TEM cells high-lowISS group
if((!file.exists(degs_f)) | rerun_analysis){
  ## Get differentially expressed genes using DESeq2
  cg_ident <- 'CD4_TEM'  
  grp <- c('highISS', 'lowISS')
  group_id <- 'group'
  Idents(seu) <- 'phenotype_detailed'
  seu_grp <- subset(seu, idents = cg_ident) # cg_ident <- 'CD4_TEM'
  DefaultAssay(seu_grp) <- 'RNA'
  
  deg <- FindMarkers(object = seu_grp,
                     only.pos = F,
                     ident.1 = grp[1], 
                     ident.2 = grp[2],
                     group.by = group_id, 
                     test.use = "DESeq2", 
                     min.pct = 0.05)
  deg$gene <- rownames(deg)
  
  saveRDS(deg, file=degs_f)
} else {
  deg <- readRDS(degs_f)
}

## Annotation of peaks based on motifs, selected for genes with peaks in promoters
if((!file.exists(peaks_f)) | rerun_analysis){
  ## Annotate and grab all the peaks for motif_i that are in the promoters
  DefaultAssay(seu) <- 'peaks'
  peaks_mat <- GetAssayData(seu)
  motifs_mat <- Motifs(seu)@data
  colnames(motifs_mat) <- motif_names[colnames(motifs_mat)]
  
  # Transform peaks into GRanges object
  peaks_id <- do.call(rbind, strsplit(rownames(peaks_mat), split="-"))
  peaks_gr <- GRanges(seqnames=peaks_id[,1], 
                      ranges=IRanges(start = as.integer(peaks_id[,2]), 
                                     end = as.integer(peaks_id[,3])))
  # Get a peak-list matched index of motifs
  motif_idx <- apply((motifs_mat!=0), 2, which)
  
  # Go motif by motif to find the tag-matrix
  motif_annos <- lapply(names(motif_idx), function(motif_name){
    # separate out peaks that contain motif_i
    print(paste0(motif_name, " (", match(motif_name, names(motif_idx)), 
                 "/", length(motif_idx), ")"))
    motif_i <- motif_idx[[motif_name]]
    peaks_i <- peaks_gr[motif_i,]
    seqlevelsStyle(peaks_i) <- 'NCBI'
    peak_i_anno <- annotatePeak(peaks_i, tssRegion=c((-1*promoter_flank), promoter_flank),
                                TxDb=txdb, level="gene")
    promoter_i_anno <- peak_i_anno@anno %>%
      as.data.frame() %>%
      filter(grepl("Promoter", annotation))
    return(list("promoter"=promoter_i_anno, "peak"=peak_i_anno))
  })
  names(motif_annos) <- names(motif_idx)
  
  saveRDS(motif_annos, file=peaks_f)
} else {
  deg <- readRDS(peaks_f)
}

# Select only the promoter-subsetted peaks
prom_annos <- lapply(motif_annos, function(i) {
  p <- i$promoter
  p$gene <- gene_ids[p$geneId]
  return(p)
})

# Crappy hack together dataframe of genes by DEG and peaks in Motif_i
deg_genes <- deg[which(deg$p_val_adj < 0.05),]$gene
non_deg_genes <- deg[which(deg$p_val_adj > 0.05),]$gene
all_motif_genes <- unique(sort(unlist(sapply(prom_annos, function(p) p$gene))))

gene_df <- data.frame("gene"=unique(sort(c(deg_genes, non_deg_genes, all_motif_genes))))
gene_df$deg <- FALSE
gene_df$deg[which(gene_df$gene %in% deg_genes)] <- 'DEG'
gene_df$deg[which(gene_df$gene %in% non_deg_genes)] <- 'nonDEG'
gene_df$motif <- FALSE
gene_df$motif[which(gene_df$gene %in% all_motif_genes)] <- 'allMotif'

stats_l <- lapply(prom_annos, function(p){
  motif_genes <- p$gene
  gene_df$motif[which(gene_df$gene %in% motif_genes)] <- 'motif'
  contigency_mat <- table(gene_df[,c('deg', 'motif')])[-2, -2]
  p <- fisher.test(contigency_mat)$p.val
  residuals <- chisq.test(contigency_mat)$residuals
  stdres <- chisq.test(contigency_mat)$stdres
  
  return(list('p'=p, 'residuals'=residuals))
})

stats_df <- data.frame('p'=sapply(stats_l, function(i) i$p)) %>%
  mutate(padj = p.adjust(p, method='fdr')) %>% 
  tibble::rownames_to_column("Gene") %>%
  mutate(residual_deg_motif = sapply(stats_l, function(i) i$residuals['DEG', 'motif'])) %>%
  dplyr::arrange(padj)
write.table(stats_df, file=out_stats_f, sep=",",
            col.names = T, row.names = F, quote = F)
  


