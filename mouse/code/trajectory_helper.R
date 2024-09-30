.break2number <- function(i, ret='median'){
  strl <- strsplit(gsub('[()]|\\[|\\]', "", as.character(i)), split=",")
  switch(ret,
         median=sapply(strl, function(j) median(as.numeric(j))),
         min=sapply(strl, function(j) min(as.numeric(j))),
         max=sapply(strl, function(j) max(as.numeric(j))))
}

.to_simple_p <- function(p){
  breaks <- c(0, 0.001, 0.05, 0.1, 1)
  pcut <- cut(c(breaks, p), breaks=breaks,
              include.lowest = T)
  levels(pcut) <- c('p<0.001', 'p<0.05', 'p<0.1', 'n.s.')
  return(as.character(pcut)[-seq_along(breaks)])
}

# .clone_cds_params(new.assay=regulonAUC[[seuid]],
#                   ref.cds=traj_i$cds$cds[genes,batch.noninf.idx])
.clone_cds_params <- function(new.assay, ref.cds, reduction_method=NULL){
  data <- as(new.assay[,colnames(ref.cds)], 'sparseMatrix')
  pd <- colData(ref.cds)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  cds <- new_cell_data_set(expression_data=data,
                           cell_metadata  = pd,
                           gene_metadata = fData)
  reducedDims(cds) <- reducedDims(ref.cds)
  if(is.null(reduction_method)) reduction_method <- names(ref.cds@principal_graph)
  cds@principal_graph[[reduction_method]] <- ref.cds@principal_graph[[reduction_method]]
  cds@principal_graph_aux[[reduction_method]] <- ref.cds@principal_graph_aux[[reduction_method]]
  return(cds)
}


runMonocle3 <- function(dat, nodes_l=NULL, reduction='umap', ncenter=275, condition='orig.ident',
                        celltype_group='treg_anno', comp_group='orig.ident'){
  ## Preprocess the seurat to monocle cds
  if(class(dat)=='Seurat'){
    seu <- dat
    data <- as(as.matrix(seu@assays$RNA$counts), 'sparseMatrix')
    pd <- seu@meta.data
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    
    cds <- new_cell_data_set(expression_data=data,
                             cell_metadata  = seu@meta.data,
                             gene_metadata = fData) %>%
      estimate_size_factors(.) %>%
      preprocess_cds(., num_dim = 50, norm_method='log') %>%
      align_cds(., alignment_group = comp_group)
    reducedDims(cds)$UMAP <- Embeddings(seu, reduction)
    cds <- cds %>% 
      cluster_cells(.) %>%
      learn_graph(., use_partition = TRUE, close_loop = TRUE,
                  list(ncenter=ncenter))
  } else {
    cds <- dat
  }
  
  ## order cells and add in trajectory information
  if(is.null(nodes_l)){
    pdf("~/xfer/principal_nodes.pdf")
    plot_cells(cds, color_cells_by = celltype_group, show_trajectory_graph=T,
               label_cell_groups=F,
               label_principal_points=T) %>%
      print
    dev.off()
    message("Principal nodes umap plotted to ~/xfer/principal_nodes.pdf")
    res <- cds
  } else {
    trajdat <- lapply(nodes_l, function(nodes_i){
      print(nodes_i)
      start_node <- nodes_i[1]
      end_node <- nodes_i[-1]
      cds <- order_cells(cds, reduction_method='UMAP', root_pr_nodes=start_node)
      pseudotime <- matrix(pseudotime(cds))
      if(all(is.na(end_node))){
        df <- data.frame(weights = as.numeric(colnames(cds) %in% Cells(cds)[!is.infinite(pseudotime)]))
      } else {
        path <- choose_graph_segments(cds,
                                      reduction_method = "UMAP",
                                      starting_pr_node = start_node,
                                      ending_pr_nodes = end_node,
                                      return_list=T)
        df <- data.frame(weights = as.numeric(colnames(cds) %in% path$cells))
      }
      
      return(list('weights'=df, 'pseudotime'=pseudotime, 'cds'=cds))
    }) 
    cellWeights <- lapply(trajdat, function(i) i$weights) %>% 
      do.call(cbind, .) %>%
      magrittr::set_colnames(., names(nodes_l)) %>%
      magrittr::set_rownames(., colnames(cds))
    pseudotime <- lapply(trajdat, function(i) i$pseudotime) %>% 
      do.call(cbind, .) %>%
      magrittr::set_colnames(., names(nodes_l)) %>%
      magrittr::set_rownames(., colnames(cds))
    cds_all <- lapply(trajdat, function(i) i$cds)
    
    res <- list("pseudotime" = pseudotime,
                "cellWeights" = cellWeights,
                "conditions" = colData(cds)[,condition],
                "cds"=cds,
                'cds_all'=cds_all)
  }
  return(res)
}

.splitAndIntegrate <- function(seu, splitid='orig.ident',
                               min.dist=0.2, n.neighbors=30, ndim=30, res=0.9){
  seusub <- DietSeurat(seu, assays='RNA', layers='counts')
  seusub <- split(seusub, f=seusub@meta.data[,splitid])
  seusub <- NormalizeData(seusub, normalization.method = "LogNormalize", 
                          scale.factor = 10000) %>%
    FindVariableFeatures(., num.bin = 100, binning.method = "equal_frequency")  %>%
    ScaleData(.,  model.use = "negbinom", do.center = T, do.scale = F, use.umi = T)  %>%
    RunPCA(., pcs.compute = 50, pcs.print = 1:3, genes.print = 10, pc.genes = rownames(seusub@data)) %>%
    FindNeighbors(., dims = 1:30, reduction = "pca")
  
  DefaultAssay(seusub) <- 'RNA'
  Idents(seusub) <- splitid
  seu_integ <- IntegrateLayers(
    object = seusub, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = T
  )
  seusub <- JoinLayers(seu_integ)
  seusub <- seusub %>% 
    FindNeighbors(., reduction = "integrated.mnn", dims = 1:ndim) %>%
    FindClusters(., resolution = res, cluster.name = "mnn_clusters") %>%
    RunUMAP(., reduction = "integrated.mnn", dims = 1:ndim, 
            n.neighbors = n.neighbors, min.dist = min.dist)
  
  return(seusub)
}

plotProgressionAndCelltype <- function(cds, pseudoval=0.01){
  ## Aggregate the pseudotime and weights per partition
  res_dat <- lapply(seq_len(ncol(cds$pseudotime)), function(idx){
    df <- data.frame(pseudotime = cds$pseudotime[,idx],
                     conditions = cds$conditions) %>%
      mutate(partition=colnames(cds$pseudotime)[idx],
             anno=colData(cds$cds)$treg_anno) %>%
      dplyr::filter(!is.infinite(pseudotime))
    return(df)
  }) %>%
    setNames(., colnames(cds$pseudotime))
  
  ## Run condiments::progressionTest
  prog_res <- lapply(names(res_dat), function(resid){
    print(resid)
    resdf <- res_dat[[resid]]
    nontraj_idx <- which(cds$cellWeights[,resid,drop=F] != 0)
    cellids <- intersect(rownames(resdf), Cells(cds$cds)[nontraj_idx])
    if(resid == 'Tpex_all'){
      cellids <- rownames(resdf)
      cds$cellWeights[cellids,resid] <- 1
    }
    res <- progressionTest(cds$pseudotime[cellids,resid,drop=F], 
                           cds$cellWeights[cellids,resid,drop=F], 
                           conditions = cds$conditions[cellids])
    data.frame("statistic"=res$statistic, 
               "p"=res$p.value, 
               "partition"=resid) %>%
      mutate(label=paste0("statistic=", round(statistic, 3), 
                          "; ", .to_simple_p(p)))
  }) %>% do.call(rbind, .)
  
  
  ## Create the proportion of celltype per pseudotime per partition
 resdf <- do.call(rbind, res_dat) %>% 
    as.data.frame  %>%
    mutate(pseudotime_breaks = cut(pseudotime, breaks=100),
           anno=factor(anno))
  break_prop <- .breakPseudotime(resdf)
  
  
  ## Visualize
  gg <- ggplot(resdf, aes(x = pseudotime)) +
    geom_density(data=resdf, 
                 aes(y=..scaled.., fill=conditions),
                 alpha = .70) + 
    geom_label(data=prog_res, aes(label=label), x=Inf, y=Inf, hjust=1, vjust=1) +
    geom_bar(data=break_prop, 
             aes(x=pseudotime, y=proportion, fill=celltype),
             position='stack', stat='identity') +
    scale_fill_brewer(type = "qual") +
    facet_grid(partition~., space='free') +
    cowplot::theme_cowplot() +
    ylim(-0.25, 1) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.2,"line"),
          legend.box="vertical", legend.margin=margin()) +
    guides(fill=guide_legend(nrow=4))
  
  return(list("gg"=gg, "res"=resdf, "test"=prog_res))
}


.breakPseudotime <- function(resdf, colids=c('partition', 'pseudotime_breaks', 'pseudotime', 'lineages'),
                             torange=c(0, -0.25), calculate='anno', method='montecarlo',
                             ref_proportions=NULL, simulation_n=1000){
  ## Create the proportion of celltype per pseudotime per partition
  stopifnot(all(c('partition', 'pseudotime_breaks', 'lineages', 'pseudotime') %in% colnames(resdf)))
  .tbl2df <- function(x){
    as.data.frame(t(as.matrix(x)))
  }
  
  break_prop <- lapply(split(resdf, f=resdf$partition), function(lineage_ij){
    if(calculate == 'anno'){
      # calculate <- 'celltype_proportion'
      datlist <- lapply(split(lineage_ij, lineage_ij$pseudotime_breaks), function(pseudotime_x){
        .tbl2df(tail(sort(table(pseudotime_x[,calculate])/nrow(pseudotime_x)), n=3))
      })
    } else if(calculate == 'sample'){
      stopifnot(!is.null(ref_proportions))
      datlist <- lapply(split(lineage_ij, lineage_ij$pseudotime_breaks), function(pseudotime_x){
        sample_cnt <- table(pseudotime_x[,calculate])
        if(all(sample_cnt==0)){
          stdres <- .tbl2df(sample_cnt)
        } else {
          if(method=='montecarlo'){
            mc_cnt <- sapply(seq_len(simulation_n), function(i){
              table(factor(
                sample(names(ref_proportions), 
                       size=sum(sample_cnt), 
                       prob = ref_proportions, replace=T),
                levels=names(sample_cnt)))
            }) 
            sd <- rowSds(mc_cnt)
            mu <- rowMeans(mc_cnt)
            stdres <- .tbl2df((sample_cnt - mu) / sd)
          } else if(method=='chisq'){
            stdres <- .tbl2df(chisq.test(sample_cnt, p=ref_proportions[names(sample_cnt)])$stdres)
          } else {
            stop("method must be either chisq or montecarlo")
          }
        }
        return(stdres)
      })
    } else if(calculate == 'imbalance_score') {
      datlist <- lapply(split(lineage_ij, lineage_ij$pseudotime_breaks), function(pseudotime_x){
        # data.frame("mean"=mean(pseudotime_x$imbalance_score),
        #            "sd"=sd(pseudotime_x$imbalance_score),
        #            "mean"=median(pseudotime_x$imbalance_score))
        t(pseudotime_x$imbalance_score) %>%
          as.data.frame %>% 
          magrittr::set_colnames(make.unique(rep("imbalance", nrow(pseudotime_x))))
          
      })
    }
    
    df <- datlist %>%
      rbind.fill(.) %>%
      mutate(pseudotime_breaks= names(split(lineage_ij, lineage_ij$pseudotime_breaks)),
             pseudotime = .break2number(pseudotime_breaks, ret='median')) %>%
      as.data.frame
    
    addedcols <- setdiff(colids, c(colnames(df), calculate))
    for(eachcol in addedcols){
      if(length(unique(resdf[,eachcol])) > 1) stop(paste0(eachcol, " has too many unique values"))
      df[,eachcol] <- rep(unique(resdf[,eachcol]), nrow(df))
    }
    return(df)
  })  %>% rbind.fill(.)
  break_prop <- pivot_longer(break_prop, 
                             cols=!all_of(colids),
                             names_to='celltype',
                             values_to = 'proportion') 
  if(!is.null(torange)) break_prop <- dplyr::mutate(break_prop, proportion = scales::rescale(proportion, to = torange))
  return(break_prop)
}

tradeseq.predictSmooth_conditions <- function(models, gene, pseudotime_breaks, tidy){
  ## Class-level preprocessing
  # models = gams[[compid]][[partitionid]]
  id = gene
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  beta <- as.matrix(rowData(models)$tradeSeq$beta[[1]][id,])
  pseudotime <- as.data.frame(colData(models)$crv) %>% 
    dplyr::select(grep("pseudotime", colnames(.), value=T)) %>%
    as.matrix
  conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
  ## End of class level preprocessing 
  
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  nConditions <- nlevels(conditions)
  nPoints <- nrow(pseudotime_breaks)
  
  # get predictor matrix
  if (tidy) out <- list()
  for (jj in seq_len(nCurves)) {
    if (tidy) out_cond <- list()
    for(kk in seq_len(nConditions)){
      df <- tradeSeq:::.getPredictRangeDf(dm, lineageId = jj, conditionId = kk,
                                          nPoints = nPoints)
      df$t1 <- pseudotime_breaks$pseudotime
      
      
      Xdf <- tradeSeq:::predictGAM(lpmatrix = X,
                                   df = df,
                                   pseudotime = pseudotime,
                                   conditions = conditions)
      if(kk == 1) XallCond <- Xdf
      if(kk > 1) XallCond <- rbind(XallCond, Xdf)
      if (tidy) {
        out_cond[[kk]] <- data.frame(lineage = jj, time = df[, paste0("t",jj)],
                                     condition = levels(conditions)[kk])
      }
    }
    if (jj == 1) Xall <- XallCond
    if (jj > 1) Xall <- rbind(Xall, XallCond)
    if (tidy) out[[jj]] <- do.call(rbind, out_cond)
  }
  if (tidy) outAll <- do.call(rbind, out)
  
  # loop over all genes
  yhatMat <- matrix(NA, nrow = length(gene), ncol = nCurves * nConditions * nPoints)
  rownames(yhatMat) <- gene
  pointNames <- expand.grid(1:nCurves, 1:nConditions)
  baseNames <- paste0("lineage", rep(1:nCurves, each=nConditions), "_condition",
                      levels(conditions)[rep(1:nConditions, nCurves)])
  colnames(yhatMat) <- c(sapply(baseNames, paste0, "_point",1:nPoints))
  for (jj in 1:length(gene)) {
    yhat <- c(exp(t(Xall %*% t(beta[as.character(gene[jj]), ,
                                    drop = FALSE])) +
                    df$offset[1]))
    yhatMat[jj, ] <- yhat
  }
  ## return output
  if (!tidy) {
    return(yhatMat)
  } else {
    outList <- list()
    for (gg in seq_len(length(gene))){
      curOut <- outAll
      curOut$gene <- gene[gg]
      curOut$yhat <- yhatMat[gg,]
      outList[[gg]] <- curOut
    }
    return(do.call(rbind, outList))
  }
}