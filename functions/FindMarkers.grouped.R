# 
# FindMarkers.grouped(seu,meta=meta, idents=idents, assay='RNA')
# 
# 

FindMarkers.grouped <- function(object, meta, idents, 
                                reduction = NULL, 
                                fc.method='default',
                                assay = NULL, slot = "data", 
                                pseudocount.use = 1, base = 2,
                                features = NULL, min.pct = 0.01, min.diff.pct = -Inf){
  .getcells <- function(object, ident_i, cellnames.use){
    cells <- Seurat:::IdentsToCells(object = object, ident.1 = ident_i[1], 
                                    ident.2 = ident_i[2], cellnames.use = cellnames.use)
    sapply(X = cells, FUN = intersect, y = cellnames.use, 
           simplify = FALSE, USE.NAMES = TRUE)
  }
  
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    assayobj <- object[[assay]]
    assaydat <- GetAssayData(object, assay=assay)
    cellnames.use <- colnames(x = assayobj)
  } else {
    assayobj <- object[[reduction]]
    assaydat <- GetAssayData(object, assay=assay)
    cellnames.use <- rownames(x = assayobj)
  }
  features <- features %||% rownames(x = object)
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && 
                     is.null(x = reduction)) {
    Command(object = object, command = norm.command, value = "normalization.method")
  } else {
    NULL
  }
  
  data.slot <- 'data'
  data.use <- GetAssayData(object = assayobj, slot = data.slot)
  counts <- switch(EXPR = data.slot, scale.data = GetAssayData(object = data.use, 
                                                               slot = "counts"), numeric())
  
  ## Identify features to analyze based on minimum alphas
  features.alpha <- lapply(idents, function(ident_i){
    thresh.min <- 0
    cells <- .getcells(object, ident_i, cellnames.use)
    pct.1 <- round(x = rowSums(
      x = data.use[features, cells$cells.1, drop = FALSE] > thresh.min)/length(x = cells$cells.1), digits = 3
    )
    pct.2 <- round(x = rowSums(
      x = data.use[features, cells$cells.2, drop = FALSE] > thresh.min)/length(x = cells$cells.2), digits = 3
    )
    alpha.min <- setNames(pmax(pct.1, pct.2), features)
    alpha.diff <- setNames(alpha.min - pmin(pct.1, pct.2), features)
    features.i <- names(x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct))
    return(features.i)
  })
  features <- purrr::reduce(features.alpha, intersect)
  
  # Find the cell barcodes for each identity
  cells <- sapply(meta$IDs, function(id){
    Cells(object)[which(Idents(object) %in% id)]
  })
  
  # Assemble the groupings metadata
  dat <- data.frame("cell"=unlist(cells), row.names = NULL) %>%
    mutate(ID=gsub("_(KO|WT).*", "", cell))
  for(colid in colnames(meta)[-1]){
    dat[,colid]=rep(meta[,colid], sapply(cells, length))
  }
  
  ## Two-way repeated measures ANOVA for condition:treatment interaction
  message("Running two-way repeated measures ANOVA...")
  aov.p <- sapply(features, function(f){
    dat$expr <- data.use[f,dat$cell]
    aov.res <- rstatix::anova_test(data=dat, expr ~ Condition * Treatment)
    with(aov.res, setNames(p, Effect))
  }) %>% 
    t %>% as.data.frame
  aov.p$padj_int <- p.adjust(aov.p$`Condition:Treatment`, method='BH')
    
  fc.res <- lapply(names(idents), function(id){
    ident_i <- idents[[id]]
    cells <- .getcells(object, ident_i, cellnames.use)
    if(fc.method == 'default'){
      message("Calculating effect size: Log2FC...")
      res <- FoldChange(object = assayobj, slot = data.slot, 
                 cells.1 = cells$cells.1, cells.2 = cells$cells.2, 
                 features = features, pseudocount.use = pseudocount.use, 
                 base = base, norm.method = norm.method) %>%
        rename_with(., ~paste0(id, ".", .))
    } else if (fc.method=='cohensd'){
      message("Calculating effect size: Cohens D...")
      res <- apply(assaydat, 1, function(i){
        x <- effsize::cohen.d(i[cells$cells.1], i[cells$cells.2])
        setNames(c(x$conf.int, x$estimate), c('D.lower', 'D.upper', 'D'))
      }) %>%
        t %>% as.data.frame %>%
        mutate(D.lower = round(D.lower, sigdig),
               D.upper =round(D.upper, sigdig),
               D =round(D, sigdig)) %>%
        rename_with(., ~paste0(id, ".", .))
    }
    return(res)
  })  %>%
    do.call(cbind, .)
  
  
  aov.fc.res <- left_join(tibble::rownames_to_column(aov.p, "dataid"),
            tibble::rownames_to_column(fc.res, "dataid"),
            by='dataid') %>% 
    tibble::column_to_rownames(., "dataid") %>% 
    arrange(padj_int)
  return(aov.fc.res)
}

