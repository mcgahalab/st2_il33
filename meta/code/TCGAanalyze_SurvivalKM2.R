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

# clinical_patient <- meta_surv
# data_grp <- patient_grp
TCGAanalyze_SurvivalKM3 <- function (clinical_patient, data_grp, 
                                     Survresult = FALSE, ...) {
  require(plotrix)
  TCGAbiolinks:::check_package("survival")
  all_samples <- as.character(unlist(data_grp))

  cfu <- clinical_patient[clinical_patient$bcr_patient_barcode %in% all_samples, ] %>%
    select(bcr_patient_barcode, days_to_death,
           days_to_last_follow_up, vital_status) %>%
    as.data.frame()
  cfu$time <- apply(cfu[,c('days_to_death', 'days_to_last_follow_up')], 1, function(i){
    i <- as.numeric(i)
    i[which(!is.infinite(i))]
  })
  cfu$status <- setNames(c(1,0), c("Dead", "Alive"))[cfu$vital_status]
  
  if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0) 
    cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <- "-Inf"
  if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0) 
    cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <- "-Inf"
  
  cfu_complete <- cfu %>%
    filter(!is.na(days_to_last_follow_up),
           !is.na(days_to_death)) %>%
    mutate(days_to_death = as.numeric(days_to_death),
           days_to_last_follow_up = as.numeric(as.character(days_to_last_follow_up))) %>%
    tibble::column_to_rownames(., var='bcr_patient_barcode')

  # subset the metadata to match the groups
  data <- matrix(unlist(data_grp)) %>%
    as.data.frame() %>%
    rename_with(., ~ "samples") %>%
    mutate(group=rep(names(data_grp), sapply(data_grp, length))) %>%
    filter(!duplicated(samples),
           samples %in% rownames(cfu_complete))
  
  cfu_complete <- cfu_complete[match(data$samples,rownames(cfu_complete)),]
  cfu_complete$grp <- data$group
  cfu_split <- split(cfu_complete, f=data$group)
  
  
  # count the number of deaths
  ttime <- as.numeric(cfu_complete[, "days_to_death"])
  deads_complete <- sum(status <- ttime > 0)    # number of dead [ALL]
  deads_spl <- sapply(cfu_split, cntDead)
  
  tabSurv_c <- setNames(c(deads_complete, deads_spl), 
                        c('Cancer_Deaths', paste0('Cancer_Deaths', names(deads_spl))))
  
  # estimate survival p-value and chisquare value
  sc_stats <- lapply(names(cfu_split), function(grp_i){
    sapply(names(cfu_split), function(grp_j){
      if(grp_i != grp_j){
        cfu <- do.call(rbind, cfu_split[c(grp_i, grp_j)]) %>%
          as.data.frame() %>% 
          mutate(grp=rep(c(grp_i, grp_j), sapply(cfu_split[c(grp_i, grp_j)], nrow)))
        get_surv_pval(cfu)
      } else {
        c('pval'=NA, 'chis'=NA, 'ttime'=NA)
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
    # plotSurvival(cfu_split, pval_mat, sc_stats, ...)
    print(ggPlotSurvival(cfu_complete, pval_mat, ...))
  } else {
    tabSurv_c <- NULL
  }
  
  return(tabSurv_c)
}

get_surv_pval <- function(cfu){
  ttime <- as.numeric(cfu[, "days_to_death"])
  status <- ttime > 0
  ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
  ttime[which(ttime == -Inf)] <- 0
  ttime <- survival::Surv(time=ttime, event=status)
  rownames(ttime) <- rownames(cfu)
  cfu$ttime <- ttime
  
  tabSurv_pvalue <- tryCatch({
    tabSurv <- survival::survdiff(ttime ~ cfu$grp)
    tabSurv_chis <- unlist(tabSurv)$chisq
    tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
    list('pval'=tabSurv_pvalue, 'chis'=tabSurv_chis, 'ttime'=ttime)
  }, error = function(e) {
    list('pval'=NA, 'chis'=NA, 'ttime'=NA)
  })
  return(tabSurv_pvalue)
}

cntDead <- function(cfu_x){
  if (dim(cfu_x)[1] >= 1) {
    ttime_x <- cfu_x[, "days_to_death"]  
    deads_x <- sum(ttime_x > 0)  # number of dead [BOTTOM]
  } else {
    deads_x <- 0
  }
  deads_x
}

ggPlotSurvival <- function(cfu, mytable, caption=''){
  require(gridExtra)
  require(survminer)
  fit <- survfit(Surv(time, status) ~ grp, data = cfu)
  
  # grps <- c('2>1', '1>2')
  # grps <- c("FALSE", "TRUE")
  # fit.cox <- coxph(Surv(time, status) ~ grp, data = cfu[which(cfu$grp %in% grps),])
  # exp(fit.cox$coefficients)
  # 
  # rr <- seq(-3, 3, by=0.1)
  # hazard_rr <- 2^(rr)
  # res <-  sapply(setNames(hazard_rr,rr), function(i){
  #   powerCT(formula = Surv(time, status) ~ grp,
  #                dat = cfu[which(cfu$grp %in% grps),] %>%
  #                  mutate(grp=gsub(grps[1], "E", grp)) %>%
  #                  mutate(grp=gsub(grps[2], "C", grp)),
  #                nE = fit$strata[paste0("grp=", grps[2])],
  #                nC = fit$strata[paste0("grp=", grps[1])],
  #                RR = i,
  #                # RR = exp(fit.cox$coefficients),
  #                alpha = 0.05)$power %>%
  #     as.numeric()
  # })
  # 
  # pdf("~/xfer/test.pdf", height = 2, width = 5)
  # res %>% 
  #   as.data.frame() %>%
  #   rename_with(., ~ "power") %>%
  #   tibble::rownames_to_column(., var="hazard_ratio") %>%
  #   mutate(hazard_ratio=as.numeric(hazard_ratio)) %>%
  #   ggplot(aes(x=hazard_ratio, y=power, grp=1)) +
  #     geom_point() + geom_line() +
  #     theme_classic() + xlab("log2(hazard ratio)") +
  #     geom_vline(xintercept = log2(exp(fit.cox$coefficients))) +
  #     ylim(0,1)
  # dev.off()
  
  # pdf("~/xfer/test.pdf")
  # ggsurvplot(survfit(fit.cox, newdata = cfu), data=cfu,
  #            ggtheme = theme_minimal())
  # 
  # dev.off()
  
  time_br <- ceiling(max(cfu$time)/10)
  ggsurv <- ggsurvplot(
    fit, data = cfu,
    conf.int = FALSE, pval = TRUE,
    risk.table = TRUE, risk.table.col = "strata",
    risk.table.height = 0.25, 
    ncensor.plot = TRUE, ncensor.plot.height = 0.15,
    break.time.by = if(time_br > 50) 100 else time_br, 
    ggtheme = theme_bw()
  )
  ggtbl <- ggplot(data.frame(0)) + 
    xlim(0,10) + ylim(0,10) +
    annotation_custom(tableGrob(mytable, theme = ttheme_default(base_size = 5)),
                      xmin=1, xmax=8, ymin=0, ymax=10) +
    theme_nothing()
  ggsurv_grid <- plot_grid(ggsurv$plot + ggtitle(caption), 
                           ggsurv$table + theme(legend.position='none'), 
                           ggsurv$ncensor.plot + theme(legend.position='none'),
                           nrow=3, rel_heights = c(3,2,1))
  plot_grid(ggsurv_grid, ggtbl, ncol = 2, rel_widths = c(3,1))
}

plotSurvival <- function(cfu_split, pval_mat, sc_stats, 
                         cols=NULL, caption=''){
  if(is.null(cols)){
    cols <- RColorBrewer::brewer.pal(length(cfu_split)+1, "Set1")[-1]
  }
  layout(matrix(c(1,1,3,
                  1,1,2,
                  1,1,3), nrow = 3, ncol = 3, byrow = TRUE))  
  par(mar=c(5.1, 4.1, 4.1, 1))
  ## Plot the survival curve
  titlePlot <- paste0("Kaplan-Meier: ", caption)
  grp_i <- names(cfu_split)[1]
  for(grp_j in names(cfu_split)[-1]){
    ttime <- sc_stats[[grp_i]]['ttime',grp_j][[1]]
    fit <- survival::survfit(ttime ~ rep(c(grp_i, grp_j), 
                                         sapply(cfu_split[c(grp_i, grp_j)], nrow)))
    if(grp_j == names(cfu_split)[[2]]){
      plot(fit, main = titlePlot, las=1,
           col = c(cols[match(grp_i, names(cfu_split))],
                   cols[match(grp_j, names(cfu_split))]), 
           xlab = "Days", ylab = "Survival",
           mark.time=TRUE, censor=TRUE)
    } else {
      print(match(grp_j, names(cfu_split)))
      lines(fit, col=scales::alpha(c(1,cols[match(grp_j, names(cfu_split))]),
                                   c(0,1)),
            mark.time=TRUE, censor=TRUE)
    }
  }
  
  legend("topright", bty='n', legend = paste0("expr ", names(cfu_split)), 
         cex=0.6, col = cols[c(1:length(cfu_split))], pch = 15)
  
  ## Plot the survival-curve p-values
  par(mar=c(0, 0, 0, 1))
  plot.new()
  plotrix::addtable2plot(0,0,pval_mat, bty='o', cex=0.6,ypad=2, xpad=1,
                         display.colnames=TRUE, display.rownames = TRUE,
                         hlines = TRUE, vlines = TRUE, bg='lightgrey')
}
