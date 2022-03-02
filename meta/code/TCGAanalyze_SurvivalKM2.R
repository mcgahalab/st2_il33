# clinical_patient=meta_t
# dataGE=expr_t
# Genelist=list.gene
# Survresult=TRUE
# threshcuts = c(0.25, 0.50, 0.75)
# caption=unique(meta_t$treatment)
# p.cut = 1

TCGAanalyze_SurvivalKM2 <- function (clinical_patient, dataGE, Genelist, Survresult = FALSE, 
          threshcuts = c(0.25, 0.5, 0.75), p.cut = 0.05, caption='NULL', dataset='TCGA',
          group1, group2) {
  require(plotrix)
  TCGAbiolinks:::check_package("survival")
  Genelist <- intersect(rownames(dataGE), as.character(Genelist))
  dataCancer <- dataGE[as.character(Genelist), group2, drop = FALSE]
  dataNormal <- dataGE[as.character(Genelist), group1, drop = FALSE]
  #dataCancer <- dataGE[as.character(Genelist), , drop = FALSE]
  #dataNormal <- dataGE[as.character(Genelist), , drop = FALSE]
  if(dataset=='TCGA'){
    colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)
    cfu <- clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% 
                              substr(colnames(dataCancer), 1, 12), ]
  } else {
    cfu <- clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% colnames(dataCancer), ]
  }
  if ("days_to_last_followup" %in% colnames(cfu)) 
    colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <- "days_to_last_follow_up"
  cfu <- as.data.frame(subset(cfu, select = c("bcr_patient_barcode", 
                                              "days_to_death", "days_to_last_follow_up", "vital_status")))
  if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 
      0) 
    cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), 
        "days_to_death"] <- "-Inf"
  if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 
      0) 
    cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), 
        "days_to_last_follow_up"] <- "-Inf"
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  followUpLevel <- FALSE
  tabSurv_Matrix <- matrix(0, nrow(as.matrix(rownames(dataNormal))), 
                           8)
  colnames(tabSurv_Matrix) <- c("mRNA", "pvalue", "Cancer Deaths", 
                                "Cancer Deaths with Top", "Cancer Deaths with Down", 
                                "Mean Tumor Top", "Mean Tumor Down", "Mean Normal")
  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <- as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"]
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  cfu_complete <- cfu
  ngenes <- nrow(as.matrix(rownames(dataNormal)))
  tabSurvKM <- lapply(1:nrow(as.matrix(rownames(dataNormal))), function(i) {
    cat(paste0((ngenes - i), "."))
    gene <- rownames(dataNormal)[i]
    mRNAselected <- as.matrix(rownames(dataNormal))[i]
    mRNAselected_values <- dataCancer[rownames(dataCancer) == 
                                        mRNAselected, ]
    mRNAselected_values_normal <- dataNormal[rownames(dataNormal) == 
                                               mRNAselected, ]
    if (all(mRNAselected_values == 0)) 
      next
    tabSurv_Matrix[i, "mRNA"] <- mRNAselected
    
    # get quantile values of cutoffs
    mRNAselected_values_ordered <- sort(mRNAselected_values, 
                                        decreasing = TRUE)
    mRNAselected_cuts <- as.numeric(quantile(as.numeric(mRNAselected_values_ordered), 
                                             c(0, threshcuts, 1)))
    mRNAselected_values_newvector <- mRNAselected_values
    if (!is.na(mRNAselected_values_ordered)) {
      # Subset the values above and below the quantile thresholds
      numberOfSamples <- length(mRNAselected_values_ordered)
      mRNA_cuts <- cut(as.numeric(mRNAselected_values_ordered), unique(mRNAselected_cuts), 
                       include.lowest = TRUE)
      samples_split_mRNA_selected <- split(names(mRNAselected_values_ordered), mRNA_cuts)
      
      # subset the metadata to match the groups
      cfu_split <- lapply(samples_split_mRNA_selected, function(ids){
        cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% ids, ]
      })
      
      # count the number of deaths
      ttime <- as.numeric(cfu[, "days_to_death"])
      deads_complete <- sum(status <- ttime > 0)    # number of dead [ALL]
      deads_spl <- sapply(cfu_split, function(cfu_x){
        if (dim(cfu_x)[1] >= 1) {
          ttime_x <- cfu_x[, "days_to_death"]  
          deads_x <- sum(ttime_x > 0)  # number of dead [BOTTOM]
        } else {
          deads_x <- 0
        }
        deads_x
      })
      tabSurv_c <- setNames(c(deads_complete, deads_spl), 
                            c('Cancer_Deaths', paste0('Cancer_Deaths', names(deads_spl))))
      
      # Subset the cancer samples for top and bottom expressed gene
      dataCancer_spl <- lapply(samples_split_mRNA_selected, function(ids){
        dataCancer[rownames(dataCancer) == mRNAselected, ids, drop=FALSE]
      })
      
      # Report the mean expressino for top and bottom expression for gene
      mean_expr <- sapply(dataCancer_spl, function(i) mean(as.numeric(i)))
      tabSurv_c <- c(tabSurv_c, setNames(round(mean_expr, 2), 
                                         paste0('Mean_Tumor', names(mean_expr))))
      
      # estimate survival p-value and chisquare value
      sc_stats <- lapply(names(cfu_split), function(grp_i){
        sapply(names(cfu_split), function(grp_j){
          if(grp_i != grp_j){
            cfu <- do.call(rbind, cfu_split[c(grp_i, grp_j)])
            ttime <- as.numeric(cfu[, "days_to_death"])
            status <- ttime > 0
            
            ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
            ttime[which(ttime == -Inf)] <- 0
            ttime <- survival::Surv(time=ttime, event=status)
            rownames(ttime) <- rownames(cfu)
            
            tabSurv_pvalue <- tryCatch({
              
              tabSurv <- survival::survdiff(ttime ~ rep(c(grp_i, grp_j), sapply(cfu_split[c(grp_i, grp_j)], nrow))  )
              tabSurv_chis <- unlist(tabSurv)$chisq
              tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), 
                                                      df = 1))
              list('pval'=tabSurv_pvalue, 'chis'=tabSurv_chis, 'ttime'=ttime)
            }, error = function(e) {
              return(list('pval'=NA, 'chis'=NA, 'ttime'=NA))
            })
          } else {
            return(c('pval'=NA, 'chis'=NA, 'ttime'=NA))
          }
        })
      })
      names(sc_stats) <- names(cfu_split)
      
      pval_mat <- round(sapply(sc_stats, function(i) unlist(i['pval',])),3)
      pval_mat[is.na(pval_mat)] <- 1
      id_mat <- sapply(colnames(pval_mat), function(i){
        sapply(rownames(pval_mat), function(j) paste0(i, "_", j))
      })
      tabSurv_c <- c(tabSurv_c,setNames(pval_mat[lower.tri(pval_mat)],
                                        paste0("p_", id_mat[lower.tri(id_mat)])))
      
      # do the plotty plotties
      if (Survresult == TRUE) {
        cols <- RColorBrewer::brewer.pal(length(cfu_split)+1, "YlOrRd")[-1]
        par(mfrow=c(2,1))
        par(mar=c(3, 4.1, 4.1, 2.1))
        ## Plot the survival curve
        titlePlot <- paste0("Kaplan-Meier (", gene, "): ", caption)
        grp_i <- names(cfu_split)[1]
        for(grp_j in names(cfu_split)[-1]){
          ttime <- sc_stats[[grp_i]]['ttime',grp_j][[1]]
          fit <- survival::survfit(ttime ~ rep(c(grp_i, grp_j), sapply(cfu_split[c(grp_i, grp_j)], nrow)))
          if(grp_j == names(cfu_split)[[2]]){
            plot(fit, main = titlePlot, las=1,
                 col = cols[c(2,1)], xlab = "Days", ylab = "Survival")
          } else {
            lines(fit, col=c(cols[match(grp_j, names(cfu_split))], cols[1]))
          }
        }
        legend("topright", bty='n', legend = paste0("expr ", names(cfu_split)), 
               cex=0.6, col = cols[c(1:length(cfu_split))], pch = 15)
        
        ## Plot the survival-curve p-values
        par(mar=c(9.1, 4.1, 0, 2.1))
        plot.new()
        plotrix::addtable2plot(0,0,pval_mat, bty='o', cex=0.6,ypad=2, xpad=1,
                      display.colnames=TRUE, display.rownames = TRUE,
                      hlines = TRUE, vlines = TRUE, bg='lightgrey')
      }
      
      tabSurv_c
    } else {
      tabSurv_c <- NULL
    }
  
    return(tabSurv_c)
  })
  
  names(tabSurvKM) <- rownames(dataNormal)
  return(tabSurvKM)
}