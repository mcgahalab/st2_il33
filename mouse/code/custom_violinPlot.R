#################################
#### 12. Violin Treg Subsets ####
require(Seurat)
require(tidyverse)
require(patchwork)
require(scCustomize)

source("/path/to/Stacked_VlnPlot2.R")
seus <- readRDS(file="/path/to/CD45ST2_TReg_tumor_ln.diet.rds")

# seus <- readRDS(file=file.path(PDIR, "data", "seurat_obj", "CD45ST2_TReg_tumor_ln.diet.rds"))
# seus2 <- lapply(seus, function(seu){
#   seu <- DietSeurat(seu, layers=c('data', 'counts'))
#   seu@meta.data <- seu@meta.data[,c('orig.ident', 'Condition', 'Timepoint', 'treg_anno')]
#   return(seu)
# })
# saveRDS(seus2, file=file.path("~/xfer", "CD45ST2_TReg_tumor_ln.diet.rds"))


gene_list_plot <- c('Ccr8', 'Cxcr4')
pt.size <- 0.15
add.legend <- TRUE
colors_list <- c('KO.Trx'='#b2182b', 'KO.Un'='#ef8a62', 'WT.Trx'='#4d4d4d', 'WT.Un'='#999999')
outfile <- "~/Desktop/x.pdf"


### DON'T MODIFY BELOW, JUST RUN IT ###
#######################################
seugg <- lapply(names(seus), function(seuid){
  seuidx <- match(seuid, names(seus))
  seu <- seus[[seuid]]
  Idents(seu) <- 'treg_anno'
  seu$group <- paste0(seu$Condition, ".", gsub("[37]d", "Trx", seu$Timepoint))
  seu$group <- factor(as.character(seu$group), levels=names(colors_list))
  
  ggdat <- Stacked_VlnPlot2(
    seurat_object = seu,
    features = gene_list_plot, 
    split.by='group',
    x_lab_rotate = TRUE,
    colors_use = colors_list,
    pt.size=pt.size
  )
  ggdat <- lapply(ggdat, function(ggobj){
    ggobj <- ggobj + scale_fill_manual(values=colors_list) +
      theme(plot.margin = margin(t = 1, b = 1)) + 
      ylim(c(0,5))
    if(seuidx != 1) ggobj <- ggobj  + theme(axis.text.y=element_blank(),
                                            axis.title.y=element_blank())
    return(ggobj)
  })
  ggdat[[1]] <- ggdat[[1]] + ggtitle(seuid) + theme(plot.title=element_text())
  
  
  return(ggdat)
})


wrp.ln <- wrap_plots(plotlist = seugg[[1]], ncol = 1)
wrp.tumor <- wrap_plots(plotlist = seugg[[2]], ncol = 1)
if(add.legend){
  wrp.tumor <- wrp.tumor & theme(legend.position = "right")
  wrp.tumor <- wrp.tumor + plot_layout(guides = 'collect')
}

pdf(outfile)
wrp.ln | wrp.tumor 
dev.off()
