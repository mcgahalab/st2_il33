.reducePamkClusters <- function(pamk.res, mat, max.clus.size=50,
                                min.max=2:10){
  pamk.spl <- split(names(pamk.res$pamobject$clustering), 
                    pamk.res$pamobject$clustering)
  large.idx <- (sapply(pamk.spl, length)>max.clus.size)
  
  idx <- 0
  while(any(large.idx)){
    idx <<- idx + 1
    large.clus <- names(which(large.idx))[1]
    large.uids <- pamk.spl[[large.clus]]#names(which(pamk.res$pamobject$clustering == large.clus))
    if(length(large.uids) > max(min.max)){
      min.max.i <- min.max
    } else {
      min.max.i <- c(min(min.max):(length(large.uids)-1))
    }
    message(paste0(large.clus, " - ", length(large.uids)))
    pamk.res <- fpc::pamk(mat[large.uids,large.uids], krange=min.max.i)
    if(!is.null(pamk.res$pamobject)){
      pamk.sub.spl <- split(large.uids, pamk.res$pamobject$clustering)
      pamk.spl <- c(pamk.spl[-match(large.clus, names(pamk.spl))], 
                    pamk.sub.spl)
      names(pamk.spl) <- make.unique(names(pamk.spl))
      large.idx <- (sapply(pamk.spl, length)>max.clus.size)
    } else {
      print(idx)
      large.idx[which(large.idx)[1]] <- FALSE
    }
  }
  return(pamk.spl)
}

.hypertest <- function(x,y, universe){
  universe <- length(universe)
  enriched.in.y <- length(y)
  enriched.in.x <- length(x)
  intersect.xy <- length(intersect(x,y))
  phyper((intersect.xy-1), enriched.in.y, (universe-enriched.in.y), 
         enriched.in.x, lower.tail=FALSE)
}

.avgjaccard.genes <- function(x,y, reflist){
  # mean(refmat[x,y])
  # reflist = with(geneset_df, split(gene_symbol, f=gs_name))
  xgenes = as.character(unique(unlist(reflist[x])))
  ygenes = as.character(unique(unlist(reflist[y])))
  .jaccard(xgenes, ygenes)
}

.jaccard <- function(x,y){
  length(intersect(x, y))/length(unique(c(x,y)))
}

.semanticsim <- function(g, x, y){
  all.len <- setNames(rep(NA, length(y)), y)
  if(x %in% V(g)$name) {
    keep.idx <- which(y %in% V(g)$name)
    semlen <- sapply(shortest_paths(g, from=x, to=y[keep.idx])$vpath, length)
    all.len[keep.idx] <- semlen
  }
  return(all.len)
}

.getSubLca <- function(g, M, vertices){
  tbl <- table(M[vertices, vertices])
  message("Pruning network to get sub-fully connected networks")
  if(length(tbl) > 3){
    sub.networks <- lapply(seq_len(length(tbl)-3), function(idx){
      minval <- as.integer(names(tbl)[idx + 1])
      M.lvl <- ceiling(M[vertices,vertices]-minval)
      M.lvl[M.lvl<0] <- 0
      unlist(sapply(.fconn_networks(M.lvl), .getLca, g=g) %>% sapply(., names))
    })
  } else {
    sub.networks = NULL
  } 
  return(sub.networks)
}
  
.getLca <- function(g, vertices, get.sublvls=T, 
                    M=NULL, force.subnetwork=T, verbose=F){
  if(length(vertices)==1) return(NULL)
  parent_nodes = attr(lca(g, vertices, M=M, 
                          get.sublvls=get.sublvls, 
                          force.subnetwork=force.subnetwork), 'names')
  
  if(verbose) message(paste(GOTERMS[parent_nodes], collapse="\n"))
  if(verbose) message(paste(as.character(GOTERMS[i]), collapse="\n"))
  if(length(parent_nodes) >0){
    sapply(parent_nodes, function(pn){
      sapply(shortest_paths(g, from=pn, to=vertices)$vpath, length)
    }) %>% colSums %>% which.min
  } else {
    NULL
  }
}

lca <- function(g, vertices, get.sublvls=F, M=NULL, force.subnetwork=T) {
  pathi = ego(g, order=length(V(g)), nodes=vertices, mode='in') # original 'in'
  res <- V(g)[[(Reduce(intersect, pathi))]]
  
  if((length(res) == 0)){
    if(force.subnetwork & !is.null(M)){
    subgo <- .getSubLca(g, M, vertices)
    res <- as.character(subgo[[which.max(!sapply(subgo, is.null))]])
    # cnts <- sapply(pathi, attr, which='names') %>% unlist %>% table
    # res <- names(cnts[cnts==max(cnts)])
    attr(res, 'names') = res
    }  else {
    res <- NA
    }
  }
  return(res)
}


## Get all fully connected subnetworks from a binary matrix
.fconn_networks <- function(X){
  g <- graph_from_adjacency_matrix(X)
  cg <- components(g)
  split(names(cg$membership), cg$membership)
}


# Given a list of GOIDs from different modalities, tries to combine
# them using semantic similarity and then report the LCA for each
# group
semantic.combine <- function(golist, g, min.cnt = 1){
  # Integrate the 3 modalities
  goterms <- golist[!sapply(golist, is.null)] %>%
    lapply(., function(i){
    names(i) <- make.unique(rep("XX", length(i)))
    return(i)
  }) %>% unlist(., recursive=T)
  X <- sapply(goterms, function(i){
    sapply(shortest_paths(g, from=i, to=goterms)$vpath, length) %>%
      setNames(., names(goterms))
  })
  X[X==0] <- max(X)+1
  X <- abs(1-round((X-1)/(max(X)-1), 2)) * 100
  
  pamk.spl <- .fconn_networks(X)
  # pamk.res <- fpc::pamk(X,
  #                       krange=min.max)
  # pamk.spl <- .reducePamkClusters(pamk.res, max.clus.size=10,
  #                                 mat=X,
  #                                 min.max=min.max)
  
  multi.idx <- sapply(pamk.spl, function(i) {
    length(unique(gsub("\\.XX.*$", "", i))) > min.cnt
  })
  pamk.spl <- pamk.spl[multi.idx]
  
  aggregate_go <- sapply(seq_along(pamk.spl), function(i){
    Xi <- X[pamk.spl[[i]],pamk.spl[[i]],drop=F]
    names(.getLca(g, goterms[pamk.spl[[i]]]))
    
    # # Remove clusters of no-association
    # diag(Xi) <- 1
    # bool <- any(Xi < 1)
    # if(bool) {
    #   # print(GOTERMS[goterms[pamk.spl[[i]]]])
    #   # print(GOTERMS[names(print(.getLca(g, goterms[pamk.spl[[i]]])))])
    #   # print(Xi)
    #   names(.getLca(g, goterms[pamk.spl[[i]]]))
    # } else {
    #   NULL
    # }
  })
  return(aggregate_go)
}