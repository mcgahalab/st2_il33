## Functions for the MSM_SSM RNAseq analysis, split by analysis-type

##########################
#### Regulon activity ####
#' Code used plot the activity of regulons inferred from SCENIC

# Calculate mean+/-sd AUC and normalized AUC (z-scale) based on cell annotation
calcMeanAUC <- function(AUC, cellAnnotation){
  normAUC <- apply(AUC, 1, scale) %>% t() %>% as.data.frame() %>%
    rename_with(., ~colnames(AUC))
  
  .getMeanSd <- function(auc, cellAnnotation){
    auc_l <- split(as.data.frame(t(auc)), cellAnnotation)
    mean_auc <- sapply(auc_l, colMeans) %>% 
      reshape2::melt() %>% 
      rename_with(., ~c('Regulon', 'cellType', 'mean'))
    sd_auc <- sapply(auc_l, function(i) apply(i, 2, sd)) %>% 
      reshape2::melt() %>% 
      rename_with(., ~c('Regulon', 'cellType', 'sd')) %>%
      mutate(sd = max(sd)-sd)
    mean_auc %>% 
      full_join(., sd_auc, by=c('Regulon', 'cellType')) %>%
      left_join(., 
                mean_auc$cellType %>% 
                  as.character() %>%
                  strsplit(., split="-|_") %>%
                  do.call(rbind, .) %>% 
                  as.data.frame() %>%
                  mutate(cellType=mean_auc$cellType) %>%
                  unique(), by='cellType')
  }
  
  list("auc"=.getMeanSd(AUC, cellAnnotation), "norm_auc"=.getMeanSd(normAUC, cellAnnotation))
}

# A function that uses ggplot to create a barplot of AUC-scores per regulon 
# and the delta-values for both the Interaction-term and the Cis_TDLN-LN 
# condition, given a matrix of both.
aucPlot <- function(int_i, grp2_i, grp1_i=NULL, regulons=NULL,  new_ids=NULL,
                    celltype=NULL, auc_pattern="CIS|PBS|LN|TDLN|DAB|DMSO",  
                    delta_grp1="delta.tdln_ln", delta_grp2='delta.auc', 
                    delta_col=NULL){
  if(is.null(delta_col)){
    delta_col <- c('delta_int'='#f03b20',
                   'delta_cis'='#636363')
  }
  
  
  intx <- reshape2::melt(int_i) %>%
    mutate(Regulon=as.character(Regulon))
  
  if(is.null(grp1_i)) grp1_i <- int_i
  grp1x <- reshape2::melt(grp1_i) %>%
    mutate(Regulon=as.character(Regulon))
  grp2x <- reshape2::melt(grp2_i) %>%
    mutate(Regulon=as.character(Regulon))
  if(is.null(regulons)){
    regulons <- .getRegulons(grp1x, grp2x)
  }
  
  gg_auc <- intx %>%
    filter(grepl(auc_pattern, variable)) %>%
    filter(Regulon %in% regulons) %>%
    mutate(Regulon=factor(as.character(Regulon), levels=regulons)) %>%
    ggplot(., aes(x=value, y=Regulon, fill=variable)) +
      geom_bar(position='dodge', stat='identity') +
      theme_classic() +
      xlim(0, 1) + xlab("AUC") +
      scale_fill_manual(values = c("#2c7fb8", "#41b6c4", 
                                   "#c51b8a", '#f768a1'), 
                        name = "") +
      ggtitle(celltype) +
      theme(legend.position='bottom',
            legend.justification = "right") 
  
  .pltDelta <- function(x, delta_var, new_id=' (Cis)', xlim=0.5, col='darkgrey'){
    ylab <- delta_var %>%
      gsub("delta", "Delta", .) %>%
      gsub("\\.auc", paste0(" ", new_id), .) %>%
      gsub("\\.tdln_ln", " (TDLN-LN)", .)
    x %>%
      filter(grepl(delta_var, variable)) %>%
      filter(Regulon %in% regulons) %>%
      mutate(Regulon=factor(as.character(Regulon), levels=regulons)) %>%
      ggplot(., aes(x=value, y=Regulon)) +
        geom_bar(position='dodge', stat='identity', fill=col) +
        theme_classic() + 
        xlim((-1 * xlim), xlim) + xlab(ylab) + 
        geom_vline(xintercept=0) +
        theme(axis.text.y=element_blank(),
              axis.line.y=element_blank(),
              axis.title.y=element_blank(),
              axis.title.x = element_text(size=8),
              axis.text.x=element_text(angle=90))
  }
  
  gg_grp1_delta <- .pltDelta(grp1x, delta_var=delta_grp1, new_id=new_ids[1], col=delta_col[1], xlim=0.5)
  gg_grp2_delta <- .pltDelta(grp2x, delta_var=delta_grp2, new_id=new_ids[2], col=delta_col[2], xlim=0.5)
  
  
  plot_grid(gg_auc, gg_grp1_delta, gg_grp2_delta, nrow=1,
            align = "h", axis='bt', rel_widths = c(3,1,1))
}

aucPlotSingle <- function(int_i, regulons=NULL,  new_ids=NULL,
                          celltype=NULL, auc_pattern="CIS|PBS|LN|TDLN|DAB|DMSO",  
                          delta_var="delta.auc"){
  delta_col <- c('delta_int'='#f03b20',
                 'delta_cis'='#636363')
  
  intx <- reshape2::melt(int_i) %>%
    mutate(Regulon=as.character(Regulon))
  
  if(is.null(regulons)){
    # regulons <- .getRegulons(grp1x, grp2x)
  }
  
  gg_auc <- intx %>%
    filter(grepl(auc_pattern, variable)) %>%
    filter(Regulon %in% regulons) %>%
    mutate(Regulon=factor(as.character(Regulon), levels=regulons)) %>%
    ggplot(., aes(x=value, y=Regulon, fill=variable)) +
      geom_bar(position='dodge', stat='identity') +
      theme_classic() +
      xlim(0, 1) + xlab("AUC") +
      scale_fill_manual(values = c("#2c7fb8", "#41b6c4"), 
                        name = "") +
      ggtitle(celltype) +
      theme(legend.position='bottom',
            legend.justification = "right") 
  
  .pltDelta <- function(x, delta_var, new_id=' (Cis)', xlim=0.5, col='darkgrey'){
    ylab <- delta_var %>%
      gsub("delta", "Delta", .) %>%
      gsub("\\.auc", paste0(" ", new_id), .) %>%
      gsub("\\.tdln_ln", " (TDLN-LN)", .)
    x %>%
      filter(grepl(delta_var, variable)) %>%
      filter(Regulon %in% regulons) %>%
      mutate(Regulon=factor(as.character(Regulon), levels=regulons)) %>%
      ggplot(., aes(x=value, y=Regulon)) +
      geom_bar(position='dodge', stat='identity', fill=col) +
      theme_classic() + 
      xlim((-1 * xlim), xlim) + xlab(ylab) + 
      geom_vline(xintercept=0) +
      theme(axis.text.y=element_blank(),
            axis.line.y=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x = element_text(size=8),
            axis.text.x=element_text(angle=90))
  }
  
  gg_intx_delta <- .pltDelta(intx, delta_var=delta_var, new_id=new_ids[1], col=delta_col[1], xlim=0.5)

  
  plot_grid(gg_auc, gg_intx_delta, nrow=1,
            align = "h", axis='bt', rel_widths = c(3,1))
}

# Stand alone function to create a demo plot fo the barplots and how to 
# interpret the results
demoAucPlot <- function(grp1='A1', grp2='B1'){
  variables <- c("CIS.tdln", "PBS.tdln", "CIS.ln", "PBS.ln")
  demo_int <- data.frame("Regulon"="TF-X (10g)",
                         "variable"=factor(variables, levels=variables),
                         "value"=c(0.6, 0.4, 0.7, 0.2))
  
  grps <- setNames(c(grp1, grp2), c('grp1', 'grp2'))
  
  # -- Set up all the comparison lines and comparison barplots
  grps_barplots <- list()
  grps_cols <- c()
  grps_lbl <- c()
  delta_lines <- as.data.frame(matrix(nrow=0, ncol=4)) %>%
    magrittr::set_colnames(., c('Regulon', 'variable', 'value', 'group'))
  for(grp_i in grps){
    grp_idx <- match(grp_i, grps)
    if(grp_i=='B1'){
      # B1: tdln_sm_delta: difference between Cis and PBS for [MS]SM_TDLN
      grps_barplots[[paste0("grp", grp_idx)]] <- data.frame("Regulon"="TF-X (10g)",
                                                            "value"=0.2)
      delta_line_tmp <- demo_int[demo_int$variable %in% c("CIS.tdln", "PBS.tdln"),]
      grps_cols <- c(grps_cols, 'indianred1')
      grps_lbl <- c(grps_lbl, "TDLN [Cis-PBS]")
    } else if(grp_i=='B2'){
      # B2: ln_sm_delta: difference between Cis and PBS for [MS]SM_LN
      grps_barplots[[paste0("grp", grp_idx)]] <- data.frame("Regulon"="TF-X (10g)",
                                                            "value"=0.5)
      delta_line_tmp <- demo_int[demo_int$variable %in% c("CIS.ln", "PBS.ln"),]
      grps_cols <- c(grps_cols, 'orangered2')
      grps_lbl <- c(grps_lbl, "LN [Cis-PBS]")
    } else if(grp_i=='A1'){
      # A1: cis_sm_delta: difference between TDLN and LN for CIS_[MS]SM
      grps_barplots[[paste0("grp", grp_idx)]] <- data.frame("Regulon"="TF-X (10g)",
                                                            "value"=-0.1)
      delta_line_tmp <- demo_int[demo_int$variable %in% c("CIS.tdln", "CIS.ln"),]
      grps_cols <- c(grps_cols, 'paleturquoise3')
      grps_lbl <- c(grps_lbl, "Cis [TDLN-LN]")
    } else if(grp_i=='A2'){
      # A2: pbs_sm_delta: difference between TDLN and LN for PBS_[MS]SM
      grps_barplots[[paste0("grp", grp_idx)]] <- data.frame("Regulon"="TF-X (10g)",
                                                            "value"=0.2)
      delta_line_tmp <- demo_int[demo_int$variable %in% c("PBS.tdln", "PBS.ln"),]
      grps_cols <- c(grps_cols, 'royalblue2')
      grps_lbl <- c(grps_lbl, "PBS [TDLN-LN]")
    } else if(grp_i=='C') {
      # C: int_sm_delta: difference between Cis and PBS for MSM_[TD]LN,
      #     [a] TDLN: delta(Cis-PBS)
      #     [b]   LN: delta(Cis-PBS)
      #        Delta: [a] - [b]; TDLN_Delta(Cis-PBS) - LN_Delta(Cis-PBS) 
      grps_barplots[[paste0("grp", grp_idx)]] <- data.frame("Regulon"="TF-X (10g)",
                                                            "value"=(0.2 - 0.5))
      delta_line_tmp <- demo_int[demo_int$variable %in% c('CIS.tdln', 'CIS.ln', "PBS.tdln", "PBS.ln"),]
      grps_cols <- c(grps_cols, c('gray30', 'gray50'))
      grps_lbl <- c(grps_lbl, "TDLN[Cis-PBS] - LN[Cis-PBS]")
      grp_idx <- c(rep(grp_idx,2), rep(grp_idx+1,2))
    }
    delta_line_tmp$group <- grp_idx
    delta_lines <- rbind(delta_lines, delta_line_tmp)
  }
  delta_lines$group <- factor(delta_lines$group)
  
  # -- Plot all the barplots and comparisons
  gg_auc <- ggplot(demo_int, aes(x=value, y=variable, fill=variable)) +
    geom_bar(position='dodge', stat='identity') +
    theme_classic() +
    xlim(0, 1) + xlab("AUC") + ylab("Regulon") +
    scale_y_discrete("Regulon",
                     labels = c("PBS.ln" = "", "CIS.ln" = "", 
                                "PBS.tdln" = "TF-X (10g)","CIS.tdln" = "")) +
    scale_fill_manual(values = c("#2c7fb8", "#41b6c4", 
                                 "#c51b8a", '#f768a1'), 
                      name = "") +
    theme(legend.position='none',
          axis.title.y=element_text(margin = margin(t = 0, r = 50, b = 0, l = 0)))
  gg_auc <- gg_auc + 
    geom_point(data = delta_lines, aes(x=value, y=variable)) +
    geom_line(data = delta_lines, aes(x=value, y=variable, group=group, col=group), size=2) +
    scale_color_manual(values=grps_cols)
  
  .pltDelta <- function(x, ylab, xlim, col){
    ggplot(x, aes(x=value, y=Regulon)) +
      geom_bar(position='dodge', stat='identity', fill=col) +
      theme_classic() + 
      xlim((-1*xlim), xlim) + xlab(ylab) + 
      geom_vline(xintercept=0) +
      theme(axis.text.y=element_blank(),
            axis.line.y=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=8),
            axis.text.x=element_text(angle=90))
  }
  gg_cis_delta <- .pltDelta(grps_barplots[[1]], grps_lbl[1], 0.5, grps_cols[1]) # '#636363')
  gg_int_delta <- .pltDelta(grps_barplots[[2]], grps_lbl[2], 0.5, grps_cols[2]) # '#f03b20')
  
  plot_grid(gg_auc, gg_cis_delta, gg_int_delta, nrow=1,
            align = "h", axis='bt', rel_widths = c(3,1,1))
  
}


# Small function to take the Interaction-term delta and the Cis_TDLN-LN
# delta, combine the values into their mean-absolute-value per regulon, and 
# then order in ascending order
.getRegulons <- function(intx, cisx,  
                         delta_int="delta.tdln_ln", delta_cis='delta.auc',
                         return_df=FALSE){
  delta_comb <- intx %>%
    filter(grepl(delta_int, variable)) %>%  # select only delta-values between condition 1-2
    select(Regulon, value) %>%              # select only the Regulon + delta-values
    full_join(cisx %>%                      # full_join with intx with cisx based on Regulon
                filter(grepl(delta_cis, variable)) %>%
                select(Regulon, value),
              by=c('Regulon')) %>%
    tibble::column_to_rownames(., var = 'Regulon') %>% 
    mutate(absDelta=rowMeans(abs(.))) %>%
    arrange(absDelta)
  
  if(return_df){
    delta_comb
  } else {
    rownames(delta_comb)
  }
}

# Highly specialized function meant to create pairings of columns by the name 
# of V1, V2, and V3 (melted)
.splitFg <- function(x, col_x='fg2'){
  x %>% mutate(fg1=paste0(V3, "_", V1),
               fg2=paste0(V2, "_", V1),
               fg3=paste0(V2, "_", V3)) %>%
    split(., .[,col_x])
}
