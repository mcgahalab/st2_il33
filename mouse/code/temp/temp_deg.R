## sara Tumor/LN KO and WT samples
# Core
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(monocle3)
library(BiocParallel)
library(GenomicFeatures)
library(tidyr)
library(Seurat)
library(dplyr)
library(plyr)
library(reshape2)
library(stringr)
# Preprocess
library(STACAS)
library(SeuratWrappers)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
library(scCustomize)
# Annotation
library(SingleR)
library(ProjecTILs)
# library(cellpypes)
# Enrichment Analysis
library(msigdbr)
library(SCPA)
library(AUCell)
library(GSEABase)
library(enrichplot)
library(pagoda2)
# Regulon Analysis
library(SCopeLoomR)
library(SCENIC)
# QC 
library(scater)
library(DoubletFinder)
# Trajcectory/Velocity
library(parsnip)
# library(velocyto.R)
# library(slingshot)
library(tradeSeq)
library(condiments)
library(reticulate)
# library(harmony)
# # library("GOstats")
library("clusterProfiler")
library(ComplexHeatmap)
library(plyr)
library(grid)
library(GO.db)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e9)

visualize <- FALSE
seed <- 1234
set.seed(seed)
args <- commandArgs(trailingOnly = TRUE)
group <- args[1]

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d'
setwd(PDIR)

gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
groups <- list.files(datadir, pattern="^(CD45|ST2|LN|Tumor)")

# # Read in gtf file to map ENS->Biotype
# gtf <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
# GTF <- rtracklayer::import(gtf)
# ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))
# sym2biotype_ids <- with(GTF, setNames(gene_biotype, gene_name))


doublet_quantile_cutoff <- 0.95

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets
species <- 'Mus musculus'

min.geneset.size=10
geneset_l <- lapply(names(msig_lvls),function(mlvl){
  lapply(msig_lvls[[mlvl]], function(sublvl){
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
      dplyr::select(gs_name, entrez_gene, gene_symbol, gs_exact_source) %>%
      as.data.frame()
    keep_gs <- names(which(table(msig_ds$gs_name) >= min.geneset.size))
    msig_ds <- msig_ds %>%
      dplyr::filter(gs_name %in% keep_gs)
    return(msig_ds)
  }) %>% 
    setNames(., unlist(msig_lvls[[mlvl]]))
})
names(geneset_l) <- names(msig_lvls)
geneset_l <- unlist(geneset_l, recursive=F)
geneset_df <-  do.call(rbind, geneset_l) %>%
  mutate(entrez_gene = gene_symbol)

GOTERMS <- (as.list(GOTERM))
GOTERMS <- sapply(GOTERMS, Term)
library(igraph)
library(GO.db)
BP <- toTable(GOBPPARENTS)
CC <- toTable(GOCCPARENTS)
MF <- toTable(GOMFPARENTS)
GOg <- graph.data.frame( rbind(BP,CC,MF) )

###################
#### Functions ####
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/list2df.R")
source("~/git/mini_projects/mini_functions/singlecell/writeSeuToFlat.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
source("~/git/mini_projects/mini_functions/singlecell/publicationDimPlot.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/bug_fixes/dietseurat.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")
gm <- geneMap(version='v2', species=species)


##-- Recurring functions 1 ----
# PacMAP clustering
runPacmap <- function(seu, reduction='pca'){
  dat_pd <- reticulate::r_to_py(as.data.frame(seu@reductions[[reduction]]@cell.embeddings))
  nparray <- dat_pd$values
  nparray <- nparray$astype(python_numpy$float)
  tfmr <- python_pacmap$PaCMAP(n_neighbors=10L, MN_ratio=0.5, FP_ratio=2.0) 
  X_transformed <- tfmr$fit_transform(nparray)#, init="pca")
  seu[["pacmap"]] <- CreateDimReducObject(embeddings = data.frame(X_transformed) %>%
                                            magrittr::set_rownames(., Cells(seu)) %>%
                                            as.matrix, 
                                          key = "PacMAP_", 
                                          assay = DefaultAssay(seu))
  return(seu)
}

.monocle3Cluster <- function(seu_i, ...){
  set.seed(1234)
  data <- as(as.matrix(seu_i@assays$mnn.reconstructed@data), 'sparseMatrix')
  pd <- seu_i@meta.data
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  cds <- new_cell_data_set(expression_data=data,
                           cell_metadata  = seu_i@meta.data,
                           gene_metadata = fData)
  cds <- preprocess_cds(cds, num_dim = 30, norm_method='size_only', pseudo_count=0)
  reducedDims(cds)$UMAP <- seu_i@reductions$umap@cell.embeddings
  # reducedDims(cds)$PCA <- seu_i@reductions$pca@cell.embeddings[,1:30]
  cds <- cluster_cells(cds, ...)
  seu_i$monocle3_clusters <- cds@clusters$UMAP$clusters
  seu_i$monocle3_partitions <- cds@clusters$UMAP$partitions
  return(seu_i)
}

## function used to create all unique 2-sample comparisons
.createGrps <- function(spl){
  grps <- lapply(spl, combn, m=2) %>%
    do.call(cbind, .)
  labels <- apply(grps, 2, function(i){
    table(c(strsplit(i[1], split = "\\.")[[1]],
            strsplit(i[2], split = "\\.")[[1]])) %>%
      sort %>% tail(., 1) %>% names %>%
      gsub("[^a-zA-Z0-9_ ]", "", .) %>%
      gsub(" ", "_", .)
  })
  grps %>%
    magrittr::set_colnames(make.unique(labels, sep = "."))
}

## function that uses regex to to help isolate which comparisons to be made
.selRows <- function(x, sample_regex, simple_regex, relevel_regex){
  x %>%
    mutate(s1=gsub(sample_regex, "\\1\\3", V1),
           s2=gsub(sample_regex, "\\1\\3", V2)) %>%
    dplyr::filter(s1 == s2) %>%
    dplyr::select(-c('s1', 's2')) %>%
    mutate(baselvl = paste0(c(grep(relevel_regex, V1, value=T, invert=T),
                              grep(relevel_regex, V2, value=T, invert=T)))) %>%
    mutate(id1=gsub(simple_regex, "\\2", V1),
           id2=gsub(simple_regex, "\\2", V2))
}

##-- Recurring functions 2 ----
.relabelid <- function(x){
  if(any(grepl("Tumor", x))){
    x %>% 
      gsub("72h", "3d", .) %>%
      gsub("^Tumor", "B1_ST2_Tumor", .) %>%
      gsub("^LN", "B1_ST2_LN", .) %>%
      gsub("^ST2", "B2_ST2", .) %>%
      gsub("^CD45", "B2_CD45", .) %>%
      gsub("PBS$", "Un", .)
  } else {
    x %>% 
      gsub("72h", "3d", .) %>%
      gsub("^LN", "B1_ST2_LN", .) %>%
      gsub("^ST2", "B2_ST2", .) %>%
      gsub("^CD45", "B2_CD45", .) %>%
      gsub("PBS$", "Un", .)
  }
}

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
    AUCell::AUCell_run(expr_mat, msig_l)
  }, error=function(e){NULL})
  return(auc)
}

addLoomCellAnnotation <- function(loom, cellAnnotation){
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation))){
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation)){
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

sapply.allbyall <- function(x, fun, ...){
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(x), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  cnt <- 0
  sapply(x, function(i){
    setTxtProgressBar(pb, cnt)
    cnt <<- cnt + 1
    sapply(x, function(j){
      fun(i,j,...)
    })
  })
  close(pb) 
}

t_col <- function(color, percent = 0.50, name = NULL) {
  sapply(seq_along(color), function(colidx){
    if(length(percent)==1){
      percenti = percent
    } else if(length(percent) != length(color)){
      warning("length(percent) != length(color); using first entry")
      percenti = percent[1]
    } else {
      percenti <- percent[colidx]
    }
    coli <- color[colidx]
    
    rgb.val <- col2rgb(coli) * (((percenti/3)*2) + 0.33)
    rgb.val <- col2rgb(coli)
    rgb.val <- ((rgb.val/255) + percenti) * 255
    rgb.val[rgb.val > 255] <- 255
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 # alpha = (100 - percenti) * 255 / 100,
                 names = name)
    invisible(t.col)
  })
}

.searchReference <- function(dat, searchterms){
  for(colx in names(searchterms)){
    # message(colx)
    dat <- dat %>%
      dplyr::filter(!!rlang::sym(colx) == searchterms[[colx]])
  }
  return(dat)
}



#########################
#### 7. DEG analysis ####
outdir <- file.path(PDIR, "results", "degs", "st2_analysis")
dir.create(outdir, recursive = T, showWarnings = F)
method <- 'wilcox'
sigdig <- 4

#--- a) ST2: All cells ----
grouping <- 'st2_all_cells'
dir.create(file.path(outdir, grouping), recursive = T, showWarnings = F)
seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_tumor_ln.rds"))
overwrite <- F

## FindAllMarkers on everything and run DotPlot
allmarkers <- lapply(names(seus), function(seuid){
  seu <- seus[[seuid]]
  Idents(seu) <- 'manual_anno2'
  seuid_markers <- FindAllMarkers(object = seu) %>%
    scCustomize::Add_Pct_Diff() 
  return(seuid_markers)
})
names(allmarkers) <-names(seus)

for(seuid in names(seus)){
  amarkers <- allmarkers[[seuid]] %>%
    dplyr::filter(pct_diff > 0.6)
  top_markers <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 20, 
                                     named_vector = FALSE, make_unique = TRUE)
  top_markers_tbl <- amarkers %>%
    scCustomize::Extract_Top_Markers(marker_dataframe = ., num_genes = 50, 
                                     named_vector = FALSE, make_unique = TRUE)
  seu <- seus[[seuid]]
  Idents(seu) <- 'manual_anno2'
  genes <- read.csv(file.path(PDIR, "ref", "dotplots", paste0("ST2_", seuid, ".csv")),
                    header = T, sep=",")
  
  write.table(amarkers %>% dplyr::filter(gene %in% top_markers_tbl),
              file=file.path("~/xfer", paste0("ST2_", seuid, ".csv")), 
              sep=",", col.names = T, row.names = F, quote = F)
  
  pdf(file.path("~/xfer", paste0("ST2_", seuid, ".pdf")), height = 20)
  scCustomize::DotPlot_scCustom(seurat_object = seu, 
                                features = unique(genes$gene),
                                flip_axes = T,
                                x_lab_rotate = TRUE)
  scCustomize::Clustered_DotPlot(seurat_object = seu, 
                                 features = unique(genes$gene),
                                 k=length(unique(Idents(seu))),
                                 x_lab_rotate = TRUE, print_exp_quantiles = T,
                                 exp_color_min = -1, exp_color_max = 3)
  
  scCustomize::DotPlot_scCustom(seurat_object = seu, 
                                features = top_markers,
                                flip_axes = T,
                                x_lab_rotate = TRUE)
  scCustomize::Clustered_DotPlot(seurat_object = seu, 
                                 features = top_markers,
                                 k=length(unique(Idents(seu))),
                                 x_lab_rotate = TRUE, print_exp_quantiles = T,
                                 exp_color_min = -1, exp_color_max = 3)
  dev.off()
  
  cat(paste0("xfer ST2_", seuid, ".csv\n"))
  cat(paste0("xfer ST2_", seuid, ".pdf\n"))
  
}

for(id in names(allmarkers)[1:2]){
  allmarkers[[id]]$cluster <- gsub(",", ".", allmarkers[[id]]$cluster)
  allmarkers[[id]]$gene <- gsub(",", ".", allmarkers[[id]]$gene)
  write.table(allmarkers[[id]],
              file=file.path(outdir, grouping, paste0("findallmarkers.", grouping, ".", id, ".csv")),
              sep=",", col.names = T, row.names = F, quote = F)
  file.copy(file.path(outdir, grouping, paste0("findallmarkers.", grouping, ".", id, ".csv")), 
            to="~/xfer", overwrite = T)
  message(paste0("xfer findallmarkers.", grouping, ".", id, ".csv"))
}






