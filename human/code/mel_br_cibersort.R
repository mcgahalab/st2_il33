# Custom plotting of CIBERSORT results for MEL-BR human samples
# Warm (red) colours are used for myeloid population
# Cool (blue) colours are used for lymphoid population

library(dplyr)
library(ggplot2)

pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/human_MEL_BR/results/immune'
setwd(pdir)

immune <- readRDS("immunedeconv.rds")

sapply(immune, function(i) levels(i$Var2))

min_frac <- 0.05
immune_i <- immune$cibersort


reduced_immune <- lapply(immune, function(immune_i){
  immune_spl <- immune_i %>%
    group_by(Var1) %>%
    filter(value > min_frac) %>% 
    as.data.frame %>% 
    group_split(Var1)
  
  # Fill in the sum of all celltypes under the min_frac as "other"
  lapply(immune_spl, function(i){
    rbind(i, 
          data.frame("Var1"=unique(i$Var1), "Var2"="Other", "value"=1 - sum(i$value)))
  }) %>% 
    plyr::rbind.fill(.) %>%
    as.data.frame
})


# lymphoid = 'blue'
# myeloid = 'red'
colors <- list("cibersort"=c('NK cell activated' = '#66a9e3', # light blue
                             'NK cell resting' = '#4ca0e1', 
                             'B cell memory' = '#346888', #dark blue-ish
                             'B cell naive' = '#7aa6c2',
                             'B cell plasma' = '#5886a5',
                             'T cell CD4+ memory activated' = '#b394c8', #purpleish
                             'T cell CD8+' = '#7a509a',
                             'T cell follicular helper' ='#935498', #purple/piunkish
                             'T cell regulatory (Tregs)' = '#c295c4',
                             'Eosinophil' = '#dec464', # yellow
                             'Mast cell activated' = '#e36748', #orange
                             'Mast cell resting' = '#e36748',
                             'Macrophage M1' = '#de425b', #reddish
                             'Macrophage M2' = '#e35c6c', 
                             'Monocyte' = '#ec9ea1', # light red
                             'Myeloid dendritic cell resting' ='#ac825a', #brownish
                             'Neutrophil' = '#e2985a', #yellow-orangeish
                             'Other'='grey'))
immune_i <- reduced_immune$cibersort
immune_i$celltype <- gsub("[0-9]+.*$", "", immune_i$Var1)
immune_i$Var2 <- factor(as.character(immune_i$Var2), 
                        levels=names(colors$cibersort))
mat_i <- dcast(immune_i,  Var2 ~ Var1) %>% 
  as.data.frame %>% 
  tibble::column_to_rownames(., 'Var2')  %>% t 
mat_i[is.na(mat_i)] <- 0
dist_i <- dist(mat_i)
hcl_i <-hclust(dist_i)
sample_ord <- colnames(as.matrix(dist_i))[hcl_i$order]
immune_i$Var1 <- factor(as.character(immune_i$Var1),
                        levels=sample_ord)

immune_j <- as.data.frame(melt(mat_i))
immune_j$celltype <- gsub("[0-9]+.*$", "", immune_j$Var1)
median_celltype <- immune_j %>% 
  # filter(grepl("^BR", Var1)) %>%
  group_by(Var2) %>%
  # summarise(mean=(1/mean(1/value))) %>%
  # summarise(mean=prod(value)^(1/length(value))) %>%
  summarise(mean=mean(value)) %>%
  arrange(desc(mean))
immune_j$Var2 <- factor(as.character(immune_j$Var2),
                        levels=median_celltype$Var2)

pdf(file.path(pdir, "cibersort_custom.pdf"))
# pdf(file.path("~/xfer", "test.pdf"))
ggplot(immune_i, aes(y=Var1, x=value, fill=Var2)) +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_grid(celltype ~ ., scales = 'free', space = 'free') +
  theme_minimal() +
  scale_fill_manual(values=colors$cibersort) +
  xlab("Estimated fraction") + ylab("") + 
  guides(fill=guide_legend(title="Cell type"))
dev.off()

pdf(file.path(pdir, "cibersort_grouped.pdf"), width = 10, height = 6)
ggplot(immune_j, aes(x=Var2, y=value, group=Var1, fill=Var2)) +
  facet_grid(. ~ celltype, scales='free', space='free') +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=colors$cibersort) +
  theme_classic() +
  xlab("") + ylab("Estimated fraction") +
  guides(fill=guide_legend(title="Cell type")) +
  theme(axis.text.x=element_blank())
dev.off()