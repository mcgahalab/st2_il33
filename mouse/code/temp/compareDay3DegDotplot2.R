#https://www.notion.so/ST2-Investigate-Day3-discrepancies-e50ae2b785074559b0d9d52b36bf0f20?pvs=4
day3deg <- c("Tox", "Tnfrsf18", "Mki67", "Maf", "Lgals1", "Klrg1", "Itgae", "Il4i1", "Il2rb", "Il2ra", "Ikzf4", "Ikzf2", "Icos", "Gzmb", "Gata3", "Foxp3", "Cxcr3", "Ctla4", "Cd69", "Ccr7", "Ccr2")

datadir <- file.path('/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor', "data")
seu_d3 <- readRDS(file = file.path(datadir, "seurat_obj", "3_seu_degPrepX.rds"))

recode_map <- function(x){
  x %>% 
    recode("0"='Tregs',
           "1"="B",
           "2"='Tregs',
           "3"="B",
           "4"='Tregs',
           "5"='Tregs; Cycling',
           "6"='Tregs; Cxcl10_hi',
           "7"='Tregs; Cycling',
           "8"="DC",
           "9"='Tregs; Intermediate',
           "10"='CD4 Th2',
           "11"='CD8 Effector',
           "12"='Monocytes',
           "13"='CD8 Naive',
           "14"="B",
           "15"='Tregs; highMt',
           "16"='Tregs; Intermediate',
           "17"='CD4 Naive',
           "18"='Tregs; Unknown',
           "19"='B.GC',
           "20"='Eosinophils',
           "21"='CD4 Th2',
           "22"='NKcells')
}


datadir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor_7d/data'
seul <- readRDS(file = file.path(datadir, "seurat_obj", "seu_integ_B3D7D.split.final.rds"))
seu_d37 <- seul$LN

############
#### D3 ####
seu <- seu_d3
seu$manual_anno <- seu$seurat_clusters %>% 
  recode_map()

idents <- unique(grep("treg", seu$manual_anno, ignore.case = T, value=T))
Idents(seu) <- 'manual_anno'
seu_treg <- subset(seu, ident=idents)

Idents(seu_treg) <- 'orig.ident'
idents <- unique(grep("LN.*72", seu_treg$orig.ident, ignore.case = T, value=T))
seu_ln_treg <- subset(seu_treg, ident=idents)


pdf("~/xfer/vlnplots_day3.pdf", width = 15, height = 25)
Idents(seu_ln_treg) <- 'manual_anno'
scCustomize::Stacked_VlnPlot(seu_ln_treg, features=day3deg, 
                             split.by='orig.ident', raster=T, x_lab_rotate=T,
                             plot_legend=T)
seu_ln_treg$all <- 'Treg'
Idents(seu_ln_treg) <- 'all'
scCustomize::Stacked_VlnPlot(seu_ln_treg, features=day3deg, 
                             split.by='orig.ident', raster=T, x_lab_rotate=T,
                             plot_legend=T)
dev.off()

dp <- DotPlot(seu_treg,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio')
dp <- DotPlot(seu_treg, features = day3deg, group.by = 'orig.ident', split.by='manual_anno', 
              scale=F, cols=rep("black", 100))
dp3 <- .groupCompare72toPbs(dp, operations='ratio', seu=seu_treg,
                            group.by = 'orig.ident', split.by='manual_anno')

pdf("~/xfer/day3_dotplot2.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()


d3_mapping <- setNames(seu_treg$manual_anno, Cells(seu_treg))
d3_cells <- Cells(seu_treg)
saveRDS(d3_cells, file="~/xfer/d3_cells.rds")

#a) --- Overlap with Day3and7-overlap cells -----
d3_d37_ov_cells <- readRDS(file="~/xfer/d3_d37_ov_cells.rds")

table(Cells(seu_treg) %in% d3_d37_ov_cells)
seu_tregx <- subset(seu_treg, cells=Cells(seu_treg)[which(Cells(seu_treg) %in% d3_d37_ov_cells)])
table(seu_tregx$orig.ident)

dp <- DotPlot(seu_tregx,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')


pdf("~/xfer/day3_dotplot2.x.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()

#################
#### Day 3+7 ####
seu <- seu_d37
idents <- unique(grep("treg", seu$manual_anno, ignore.case = T, value=T))
Idents(seu) <- 'manual_anno'
seu_treg <- subset(seu, ident=idents)

Idents(seu_treg) <- 'orig.ident'
idents <- unique(grep("B1.*3d", seu_treg$orig.ident, ignore.case = T, value=T))
seu_ln_treg <- subset(seu_treg, ident=idents)





pdf("~/xfer/vlnplots_day3and7.pdf", width = 15, height = 25)
Idents(seu_ln_treg) <- 'manual_anno'
scCustomize::Stacked_VlnPlot(seu_ln_treg, features=day3deg, 
                             split.by='orig.ident', raster=T, x_lab_rotate=T,
                             plot_legend=T)
seu_ln_treg$all <- 'Treg'
Idents(seu_ln_treg) <- 'all'
scCustomize::Stacked_VlnPlot(seu_ln_treg, features=day3deg, 
                             split.by='orig.ident', raster=T, x_lab_rotate=T,
                             plot_legend=T)
dev.off()

dp <- DotPlot(seu_treg,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')


pdf("~/xfer/day3and7_dotplot2.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()

#a) --- Overlap with Day3 cells -----
d3_cells <- readRDS(file="~/xfer/d3_cells.rds")
table(Cells(seu_treg) %in% d3_cells)
seu_tregx <- subset(seu_treg, cells=Cells(seu_treg)[which(Cells(seu_treg) %in% d3_cells)])
seu_x <- subset(seu, cells=Cells(seu)[which(Cells(seu) %in% d3_cells)])
table(seu_tregx$orig.ident)

dp <- DotPlot(seu_tregx,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')


pdf("~/xfer/day3and7_dotplot2.x.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()
saveRDS(Cells(seu_tregx), file="~/xfer/d3_d37_ov_cells.rds")


#b) --- Include all original unremoved cells -----
seu <- seul$LN

intersect_cells <- intersect(Cells(seu), d3_cells)
# length(intersect_cells); length(d3_cells)
seuln_intersect <- subset(seu, cells=intersect_cells)
seuln_intersect$d3_manual_anno <- d3_mapping[Cells(seuln_intersect)]
seul_intersect <- SplitObject(seuln_intersect, split.by='orig.ident')
table(seuln_intersect$d3_manual_anno, seuln_intersect$manual_anno)
lapply(seul_intersect, function(seu_i){
  table(seu_i$d3_manual_anno, seu_i$manual_anno)
})


dp <- DotPlot(seuln_intersect,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')

pdf("~/xfer/day3and7_dotplot2.intersect.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()

#c) --- Remove only the Bcells -----
bcells <- Cells(seuln_intersect)[seuln_intersect$manual_anno == 'B']
seuln_intersect_bcellrm <- subset(seuln_intersect, cells=setdiff(Cells(seuln_intersect), bcells))


dp <- DotPlot(seuln_intersect_bcellrm,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')

pdf("~/xfer/day3and7_dotplot2.intersect.bremove.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()

#d) --- Remove non TReg cells -----
treg_cells <- Cells(seuln_intersect_bcellrm)[grep("Treg", seuln_intersect_bcellrm$manual_anno, ignore.case=T)]
seuln_intersect_bcellrm_tregonly <- subset(seuln_intersect_bcellrm, cells=treg_cells)

dp <- DotPlot(seuln_intersect_bcellrm_tregonly,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')

pdf("~/xfer/day3and7_dotplot2.intersect.bremove.treg_only.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()


#e) --- Remove TREG_REMOVE cells -----
remove_cells <- Cells(seuln_intersect_bcellrm_tregonly)[grep("remove", seuln_intersect_bcellrm_tregonly$manual_anno, ignore.case=T)]
seuln_intersect_bcellrm_tregonly_tregrm <- subset(seuln_intersect_bcellrm_tregonly, cells=setdiff(Cells(seuln_intersect_bcellrm_tregonly), remove_cells))

dp <- DotPlot(seuln_intersect_bcellrm_tregonly_tregrm,features = day3deg, group.by = 'orig.ident',  scale=F)
dp2 <- .compare72toPbs(dp, operations='log2ratio', batch='b2')

pdf("~/xfer/day3and7_dotplot2.intersect.bremove.treg_only.tregrem.pdf")
DotPlot2(dp2, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='Pbs72')
# DotPlot2(dp3, cols=rev(RColorBrewer::brewer.pal(n=7, name='RdBu')), comparison='KoWt')
dev.off()
