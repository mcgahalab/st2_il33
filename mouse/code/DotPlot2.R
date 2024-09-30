# DotPlot2: Customized functions for Sara's ST2 TDLN/Tumor project
# to split the LN and Tumor DotPlot, and plot the log2ratio
# of the 72hr vs PBS avg expression and prcent expressed
.compare72toPbs <- function(DPdat, operations='log2ratio', batch='b1'){
  x <- DPdat$data
  # add on extra columns separating KO-WT and 72h-PBS
  if(batch=='b1'){
    x$LnTumor <- gsub("_.*", "", x$id)
    x$id <- gsub("^.*?_", "", x$id)
    metdat <- do.call(rbind, strsplit(x$id, split="_")) %>%
      as.data.frame %>% 
      rename_with(~c('KoWt', 'Pbs72'))
    x <- cbind(x, metdat) %>%
      mutate(Pbs72=factor(Pbs72, levels=c('72h', 'PBS')))
  } else {
    x$LnTumor <- gsub("^.*(LN|Tumor)_.*", "\\1", x$id, ignore.case = T)
    x$id <- gsub("_(LN|Tumor)", "", x$id)
    metdat <- do.call(rbind, strsplit(x$id, split="_")) %>%
      as.data.frame %>% 
      rename_with(~c('Batch', 'KoWt', 'Pbs72'))
    x <- cbind(x, metdat) %>%
      mutate(Pbs72=factor(Pbs72, levels=c('3d', '7d', 'Un')))
  }
  
 
  
  # create new metric summarizing the avg.exp and pct.exp
  xnew <- sapply(split(x, f=list(x$LnTumor, x$features.plot, x$KoWt)), function(i, operations){
    i <- arrange(i, Pbs72)
    colid <- c('avg.exp', 'pct.exp', 'avg.exp.scaled')
    delta <- if(operations=='diff'){
      i[1,colid] - i[2,colid] 
    } else if (operations=='ratio'){
      i[1,colid] / i[2,colid] 
    } else if (operations=='log2ratio'){
      log2(i[1,colid] / i[2,colid])
    }
    # if(operations != 'diff') delta['avg.exp.scaled'] <- delta['avg.exp']
    inew <- i[1,,drop=F]
    inew[,colid] <- delta
    return(as.matrix(inew))
  }, operations=operations) %>% 
    t %>% as.data.frame %>% 
    rename_with(., ~colnames(x)) %>%
    mutate(avg.exp=as.numeric(avg.exp),
           pct.exp=abs(as.numeric(pct.exp)),
           avg.exp.scaled=as.numeric(avg.exp.scaled))
  
  DPdat$data <- xnew
  return(DPdat)
}
.groupCompare72toPbs <- function(DPdat, operations='ratio',
                                 group.by, split.by, seu){
  # Initialize data and vars
  x <- DPdat$data
  groups <- unique(as.character(seu@meta.data[,group.by]))
  splits <- unique(as.character(seu@meta.data[,split.by]))
  mapdf <- as.data.frame(expand.grid(groups, splits)) %>%
    rename_with(., ~c('Group', 'Split')) %>%
    mutate(id=paste0(Group, "_", Split))
  
  # Label the group and split ids per uniqueid
  x <- left_join(x, mapdf, by='id')
  
  # add on extra columns separating KO-WT and 72h-PBS
  x$LnTumor <- gsub("_.*", "", x$Group)
  x$Group <- gsub("^.*?_", "", x$Group)
  metdat <- as.data.frame(stringr::str_split_fixed(x$Group, pattern="_", 2)) %>%
    rename_with(~c('KoWt', 'Pbs72'))
  x <- cbind(x, metdat) %>%
    mutate(Pbs72=factor(Pbs72, levels=c('72h', 'PBS')),
           id=Group)
  
  ## Creates new metric summarizing the avg.exp and pct.exp
  # First Pass: Subtracts 72hr from PBS
  xnew <- lapply(split(x, f=list(x$LnTumor, x$features.plot, x$Split)), function(i, operations){
   
    .preprocess <- function(ix, arrangeby, splitby, operations, log=FALSE){
      ix <- dplyr::arrange(ix, !!rlang::sym(arrangeby))  # puts 72h in 1st row, PBS in 2nd
      ixspl <- if(is.null(splitby)){
        list(ix)
      } else {
        split(ix, f=ix[,splitby])
      }

      all_ij <- lapply(ixspl, function(ij){
        if(nrow(ij)!=2) return(NULL)
        colid <- c('avg.exp', 'pct.exp', 'avg.exp.scaled')
        delta <- if(operations=='diff'){
          1 + as.numeric(ij[1,colid]) - as.numeric(ij[2,colid])
        } else if (operations=='ratio'){
          as.numeric(ij[1,colid]) / as.numeric(ij[2,colid])
        } 
        if(log) delta <- log2(delta)
        
        # if(operations != 'diff') delta['avg.exp.scaled'] <- delta['avg.exp']
        inew <- ij[1,,drop=F]
        inew[,colid] <- delta
        return(as.matrix(inew))
      }) 
      null_idx <- sapply(all_ij, is.null)
      if(any(null_idx)) {
        all_ij <- NULL
      } else {
        all_ij <- all_ij %>% 
          do.call(rbind, .) %>%
          as.data.frame 
      }
      return(all_ij)
    }
    
    i_kowt <- .preprocess(i, arrangeby='Pbs72', splitby='KoWt', 
                          operations=operations, log=FALSE)
    ifin <- if(is.null(i_kowt)){
      NULL
    } else {
      .preprocess(i_kowt, arrangeby='KoWt', splitby=NULL, 
                  operations='diff', log=FALSE)
    }
    return(ifin)
  }, operations=operations) %>% 
    do.call(rbind, .) %>% as.data.frame %>% 
    rename_with(., ~colnames(x)) %>%
    mutate(avg.exp=as.numeric(avg.exp),
           pct.exp=abs(as.numeric(pct.exp)),
           avg.exp.scaled=as.numeric(avg.exp.scaled))
  
  DPdat$data <- xnew
  return(DPdat)
}

DotPlot2 <- function(dat, cols=c("cyan", "blue", "grey", "red", "yellow"),
                     scl_range=c(2, 6), cap=5, comparison=NULL){
  
  b <- seq((-1*cap), cap, by=1)
  cap_idx <- (abs(dat$data$avg.exp)) > cap
  if(any(cap_idx)){
    dat$data$avg.exp[which(dat$data$avg.exp > cap)] <- cap
    dat$data$avg.exp[which(dat$data$avg.exp < (-1*cap))] <- (-1*cap)
  }
  
  yvar <- if(comparison == 'Pbs72'){
    'id'
  } else if (comparison == 'KoWt'){
    'Split'
  } else {
    stop("Comparison must be either 'KoWt' or 'Pbs72'")
  }
  plot <- ggplot(data = dat$data, mapping = aes_string(x = 'features.plot', y = yvar)) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(x = 'Features', y = 'Identity') +
    theme_cowplot() + 
    scale_radius(range = scl_range, 
                 limits = c(NA, NA)) +
    scale_color_gradientn(limits = c(min(b),max(b)), colours=cols,
                          breaks=b, labels=format(b))
  plot + facet_grid(.~LnTumor, scales='free', space = 'free') + 
    coord_flip() +
    theme(axis.text.x = element_text(angle=90))
}

if(FALSE){
  # Seurat object has all the samples included:
  # LN_KO_72h LN_KO_PBS LN_WT_72h LN_WT_PBS Tumor_KO_72h Tumor_KO_PBS Tumor_WT_72h Tumor_WT_PBS
  genes <- c('Mtor', 'Pten', 'Ahr')
  DefaultAssay(seu) <- 'RNA'
  Idents(seu) <- 'orig.ident'
  
  dp <- DotPlot(seu,features = genes, group.by = 'orig.ident',  scale=F)
  dp2 <- .compare72toPbs(dp, operations='log2ratio')
  
  dp <- DotPlot(seu,features = genes, group.by = 'orig.ident', split.by='manual_anno', 
                scale=F, cols=rep("black", 100))
  dp3 <- .groupCompare72toPbs(dp, operations='ratio', seu=seu,
                              group.by = 'orig.ident', split.by='manual_anno')
  
  pdf("~/xfer/b.pdf")
  DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
  DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
  dev.off()
}

