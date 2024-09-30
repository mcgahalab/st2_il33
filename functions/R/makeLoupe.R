## This is a helper function designed to output .cloupe files using the 
## loupeR package on non-seurat or single-cell data
# 
# An example of how to use this function is illustrated in the 
# .demo() function
#
# NOTE: You require a barcodes.tsv file containing column list of all acceptable 
# 10x cellranger barcodes
#

.demo <- function(){
  message("This is a demonstation piece of code on how to use the makeLoupe() function")
  message("You require a barcodes.tsv file containing column list of all acceptable
          10x cellranger barcodes")
  # You require 3 pieces of information, a count matrix, metadata in a dataframe,
  # and a projection (e.g. pca, umap, tsne, ...) that you want to encode
  mat <- as.data.frame(matrix(1:200, ncol=10)) %>%
    magrittr::set_colnames(paste0("Sample", c(1:10))) %>%
    magrittr::set_rownames(paste0("Feature", c(1:20)))
  
  meta <- data.frame("cluster"=c(rep("A",5), rep("B", 5)),
                     "group"=sample(c("X", "Y", "Z"), size=10, replace=T)) %>%
    magrittr::set_rownames(colnames(mat))
  
  # Either of the following methods works. If you want to include multiple reductions,
  # use a list of projections. If you just want to use a single projection matrix,
  # you can input just a single matrix.
  projections <- data.frame("PCA_1"=rnorm(n=10, mean=5, sd=2),
                            "PCA_2"=rnorm(n=10, mean=5, sd=2),
                            "PCA_3"=rnorm(n=10, mean=5, sd=2)) %>%
    magrittr::set_rownames(colnames(mat))
  projections = list("UMAP"=data.frame("umap_1"=rnorm(n=10, mean=5, sd=2),
                                       "umap_2"=rnorm(n=10, mean=5, sd=2)) %>%
                       magrittr::set_rownames(colnames(mat)),
                     "tSNE"=data.frame("tsne_1"=rnorm(n=10, mean=5, sd=2),
                                       "tsne_2"=rnorm(n=10, mean=5, sd=2)) %>%
                       magrittr::set_rownames(colnames(mat)))
  
  makeLoupe(mat=mat, meta=meta, projections=projections,
            output_dir="~/xfer", output_name="demo")
  # the file will be outputed as: /path/to/output_directory/demo.cloupe
}

# mat=assay(counts)
# meta=metad2
# projections=list("gene"=pca_x,"ssgsea"=pca_pathway)
# output_dir=file.path(pdir, "manual", "loupe")
# output_name="pca_all"
# barcodes_f=barcodes_f

makeLoupe <- function(mat, meta, projections,
                      output_dir, output_name,
                      barcodes_f="~/mcgahalab/ref/scrna/barcodes.tsv",
                      tmpdir=NULL){
  mat_l <- .formatMat(mat, barcodes_f)
  mat <- mat_l$mat
  meta <- .formatMeta(meta, mat, barcodemapping=mat_l$barcodemapping)
  projections <- .formatProjections(projections, mat, barcodemapping=mat_l$barcodemapping)
  
  if(any(rownames(mat) == '')){
    warning("One of your genes had no name: '' ")
    mat <- mat[-which(rownames(mat) == ''),]
  }
  loupeR:::create_loupe(count_mat=mat, clusters = meta, projections = projections, 
                        output_dir = output_dir, output_name = output_name, 
                        executable_path = NULL, force = T, seurat_obj_version = NULL,
                        tmpdir=tmpdir)
  print(paste0("Output to: ", 
               file.path(output_dir, paste0(output_name, ".cloupe"))))
  return(file.path(output_dir, paste0(output_name, ".cloupe")))
}

.formatMat <- function(mat, barcodes_f){
  ## Data validation
  mat_prelim <- as.matrix(mat)
  if(is.character(mat_prelim)){
    stop("Please insert a numeric matrix with column and row-names")
  }
  storage.mode(mat_prelim) <- 'numeric'
  
  if(any(is.na(mat_prelim))){
    warning("NA values will be converted to 0")
    mat_prelim[is.na(mat_prelim)] <- 0
  }
  mat_range <- c(min(mat_prelim), max(mat_prelim))
  if(max(mat_range) < 1){
    warning("Max value in matrix is too low, must be >1")
  } else if(min(mat_range) < 0){
    stop("Min value in matrix is too low, must be >0")
  }

  ## Mapping matrix and barcodes
  mat <- as(round(mat_prelim* 100, 1), 'dgCMatrix')
  barcodes <- read.table(barcodes_f, header=T)
  random_barcodes <- sample(barcodes$x, ncol(mat))
  barcodemapping <- setNames(
    paste(random_barcodes, colnames(mat), sep="-"),
    colnames(mat)
  )
  colnames(mat) <- as.character(barcodemapping[colnames(mat)])
  return(list(mat=mat, barcodemapping=barcodemapping))
}
.formatProjections <- function(projections, mat, barcodemapping){
  if(class(projections) %in% c('matrix', 'data.frame')){
    warning("Converting projection data into a list")
    projections=list("proj"=projections)
  }
  
  projections <- lapply(projections, function(proj_i){
    mapped_ids <- barcodemapping[rownames(proj_i)]
    if(any(is.na(mapped_ids))){
      stop("Your projection data contains samples that were not found in your count matrix")
    }
    if(length(mapped_ids) != ncol(mat)){
      stop("Your projection data is missing a few samples that were found in your count matrix")
    }
    rownames(proj_i) <- mapped_ids
    
    return(proj_i)
  })
  
  non_num_check <- any(sapply(projections, function(proj_i){
    is.character(proj_i)
  }))
  if(any(non_num_check)){
    stop("Your projection data should be a list of dataframes, where each column
         corresponds to the spatial location in the projected space")
  }
  
  projections <- lapply(projections, function(proj_i){
    if(ncol(proj_i) > 2){
      warning("Splitting projection into two-dimensional space")
      npcs <- ncol(proj_i)
      proj_mats <- lapply(2:npcs, function(pc_j){
        pc_i <- pc_j-1
        proj_i[,c(pc_i, pc_j),drop=F]
      })
      names(proj_mats) <- paste0(c(1:(npcs-1)), ".", c(2:npcs))
      proj_i <- lapply(proj_mats, function(i) round(as.matrix(i), 2))
    } else {
      proj_i <- list(proj_i)
    }
    return(proj_i)
  }) 
  
  if(any(sapply(projections, class) == 'list')) projections <- unlist(projections, recursive = F)
    
  
  if(length(names(projections))==0){
    warning("List of projections is unnamed, setting arbitrary names for projections")
    names(projections) <- paste0("P", c(1:length(projections)))
  }
  projections <- lapply(projections, as.matrix)
  return(projections)
}
.formatMeta <- function(meta, mat, barcodemapping){
  ## Data validation
  if(!any(class(meta) %in% c('data.frame', 'DFrame'))){
    stop("Your meta file must be a dataframe or related data structure with the columns
         you wish to include in your export, and rownames matching your count
         matrix columns")
  }
  meta <- as.data.frame(meta)
  
  mapped_ids <- barcodemapping[rownames(meta)]
  if(any(is.na(mapped_ids))){
    stop("Your metadata contains samples not found in your count matrix")
  }
  if(length(mapped_ids) != ncol(mat)){
    stop("Your metadata is missing a few samples that were found in your count matrix")
  }
  rownames(meta) <- mapped_ids
  
  ## Factorize a named vector and return a list of clusters
  meta <- apply(meta, 2, simplify = F, function(i){
    factor(setNames(i, rownames(meta)))
  })
  return(meta)
}

# 
# 
# h5path <- sprintf("%s.h5", tempfile())
# h5path <- '~/xfer/tmp.h5'
# feature_ids <- NULL
# ok <- loupeR:::create_hdf5(count_mat, clusters, projections, h5path, 
#                   feature_ids, seurat_obj_version)
# f <- hdf5r::H5File$new(h5path, mode = "w")
# loupeR:::write_mat(f, count_mat, feature_ids)
# loupeR:::write_clusters(f, clusters)
# loupeR:::write_projections(f, projections)
# loupeR:::write_metadata(f, seurat_obj_version)
# f$close()
# SUCCESS
# 
# ok <- loupeR:::louper_create_cloupe(h5path, 
#                            output_dir = output_dir, 
#                            output_name = output_name, 
#                            executable_path = NULL, 
#                            force = force)
# loupe_path <- sprintf("%s.cloupe", file.path(output_dir, 
#                                              output_name))
# h5path <- normalizePath(path.expand(h5path))
# loupe_path <- suppressWarnings(normalizePath(path.expand(loupe_path)))
# input_flag <- sprintf("--input=%s", h5path)
# output_flag <- sprintf("--output=%s", loupe_path)
# args <- c("create", input_flag, output_flag, "--force")
# executable_path <- loupeR:::find_executable()
# args=c('create', '--input=/cluster/home/quever/xfer/tmp.h5', '--output=/cluster/home/quever/xfer/demo.cloupe', '--force')
# executable_path='/cluster/home/quever/.local/share/R/loupeR/louper'
# 
# system2(command = executable_path, args = args)
# system(command = paste(c(executable_path, args), collapse=" "))
# # /cluster/home/quever/.local/share/R/loupeR/louper \
# # create \
# # --input=/cluster/home/quever/xfer/tmp.h5 \
# # --output=/cluster/home/quever/xfer/demo.cloupe \
# # --force
# if (status == 0) {
#   return(SUCCESS)
# }
# else {
#   return(err(sprintf("Louper executable failed: status code %d", 
#                      status)))
# }
