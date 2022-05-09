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
                    celltype=NULL, auc_pattern="CIS|PBS|LN|TDLN",  
                    delta_grp1="delta.tdln_ln", delta_grp2='delta.auc'){
  delta_col <- c('delta_int'='#f03b20',
                 'delta_cis'='#636363')
  
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
  
  
  plot_grid(gg_auc, gg_grp2_delta, gg_grp1_delta, nrow=1,
            align = "h", axis='bt', rel_widths = c(3,1,1))
}

# Stand alone function to create a demo plot fo the barplots and how to 
# interpret the results
demoAucPlot <- function(){
  variables <- c("CIS.tdln", "PBS.tdln", "CIS.ln", "PBS.ln")
  demo_int <- data.frame("Regulon"="TF-X (10g)",
                         "variable"=factor(variables, levels=variables),
                         "value"=c(0.6, 0.4, 0.7, 0.2))
  df <- rbind(demo_int[c(1,3),], demo_int)
  df$group <- c('1', '1', '2', '2', '3', '3')
  demo_delta_cis <- data.frame("Regulon"="TF-X (10g)",
                               "value"=0.4)
  demo_delta_int <- data.frame("Regulon"="TF-X (10g)",
                               "value"=-0.3)
  
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
    geom_point(data = df, aes(x=value, y=variable)) +
    geom_line(data = df, aes(x=value, y=variable, group=group, col=group), size=2) +
    scale_color_manual(values=c('#636363', '#e31a1c', '#fd8d3c'))
  
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
  gg_cis_delta <- .pltDelta(demo_delta_cis, "Delta (Cis_TDLN-LN)", 0.5, '#636363')
  gg_int_delta <- .pltDelta(demo_delta_int, "Delta (TDLN-LN)", 0.5, '#f03b20')
  
  plot_grid(gg_auc, gg_cis_delta, gg_int_delta, nrow=1,
            align = "h", axis='bt', rel_widths = c(3,1,1))
  
}

# Small function to take the Interaction-term delta and the Cis_TDLN-LN
# delta, combine the values into their mean-absolute-value per regulon, and 
# then order in ascending order
.getRegulons <- function(intx, cisx,  
                         delta_int="delta.tdln_ln", delta_cis='delta.auc'){
  delta_comb <- intx %>%
    filter(grepl(delta_int, variable)) %>%
    select(Regulon, value) %>%
    full_join(cisx %>%
                filter(grepl(delta_cis, variable)) %>%
                select(Regulon, value),
              by=c('Regulon')) %>%
    tibble::column_to_rownames(., var = 'Regulon')
  rownames(delta_comb[order(rowMeans(abs(delta_comb))),])
}

# Highly specialized function meant to create pairings of columns by the name 
# of V1, V2, and V3 (melted)
.splitFg <- function(x, col_x='fg2'){
  x %>% mutate(fg1=paste0(V3, "_", V1),
               fg2=paste0(V2, "_", V1),
               fg3=paste0(V2, "_", V3)) %>%
    split(., .[,col_x])
}
