.imbalance_score <- function(rd, conditions, k = 10, smooth = k) {
  # Code inspired from the monocle3 package
  # https://github.com/cole-trapnell-lab/monocle3/blob/9becd94f60930c2a9b51770e3818c194dd8201eb/R/cluster_cells.R#L194
  if (length(conditions) != nrow(rd)) {
    stop("The conditions and reduced dimensions do not contain the same cells")
  }
  
  props <- as.vector(table(conditions) / length(conditions))
  groups <- unique(conditions)
  n_groups <- length(groups)
  if (n_groups == 1) stop("conditions should have at least 2 classes")
  
  # Get the graph
  # We need to add 1 because by default, nn2 counts each cell as its own
  # neighbour
  tmp <- RANN::nn2(rd, rd, k + 1, searchtype = "standard")
  neighborMatrix <- tmp[[1]]
  cdMatrix <- matrix(factor(conditions)[neighborMatrix], ncol = k + 1)
  
  # Get the smoothed scores
  scores <- .multinomial.test(cdMatrix, groups, props)
  scores <- unlist(scores)
  names(scores) <- rownames(rd)
  formula <- paste0("scores ~ s(",
                    paste0("rd[, ", seq_len(ncol(rd)), "], ", collapse = ""),
                    "k = smooth)")
  mm <- mgcv::gam(stats::as.formula(formula))
  scaled_scores <- mgcv::predict.gam(mm, type = "response")
  
  return(list("scores" = scores, "scaled_scores" = scaled_scores))
}



.multinomial.test <- function (cdMatrix, groups, props) {
  size <- ncol(cdMatrix)
  eventMat <- .findVectors(length(groups), size)
  
  if(length(groups)==2){
    # Added 2 group comparisons
    # binomial distribution for pvalue of 2 categories, focus on preserving direction
    pseudoval <- 5
    pvalues <- apply(cdMatrix, 1, function(conds) {
      real <- as.vector(table(factor(conds, levels = groups)))
      p.value <- stats::pbinom(real, size, props, lower.tail = T)
      p.value[p.value==1] <- pnorm(pseudoval)
      p.value[p.value==0] <- pnorm(-1*pseudoval)
      return(p.value)
    }) %>% t
    res <- stats::qnorm(pvalues[,1])
  } else {
    eventProb <- apply(eventMat, 1, function(x) {
      dmultinom(x, size = size, prob = props)
    })
    
    # standard condiments:::.dmultinom test
    pvalues <- apply(cdMatrix, 1, function(conds) {
      real <- as.vector(table(factor(conds, levels = groups)))
      pObs <- stats::dmultinom(real, size, props)
      p.value <- sum(eventProb[eventProb <= pObs])
      return(p.value)
    })
    res <- -stats::qnorm(unlist(pvalues)/2)
  }
  return(res)
}