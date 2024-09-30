# Takes the geneSet data frame outputted from clusterProfiler::GSEA()@result
# and converts it into the UP and DOWN files with the corresponding GO/REACTOME
# identifiers in the GS.DETAILS section for plotting in enrichmentMap of 
# cytoscape
gsea2CytoscapeFormat <- function(dat, gs_map, ...){
  if(missing(gs_map)){
    gs_map <- .mapGsToExactSource(...)
  }
  
  dat <- dat %>% 
    mutate(Details=gs_map[ID]) %>%
    dplyr::select(ID, Description, Details, setSize, enrichmentScore, NES, pvalue,
           p.adjust, qvalues, rank, leading_edge) %>%
    rename_with(., ~c('NAME', 'GS.br..follow.link.to.MSigDB', 'GS.DETAILS', 
                      'SIZE', 'ES', 'NES', 'NOM.p.val', 'FDR.q.val', 
                      'FWER.p.val', 'RANK.AT.MAX', 'LEADING.EDGE'))
  datl <- split(dat, (dat$NES > 0))
  return(list("up"=datl[['TRUE']],
              "down"=datl[['FALSE']]))
}

# Takes the geneSet data frame outputted from clusterProfiler::GSEA()@result
# and converts it into the UP and DOWN files with the corresponding GO/REACTOME
# identifiers in the GS.DETAILS section for plotting in enrichmentMap of 
# cytoscape
findmarkers2CytoscapeFormat <- function(dat, gs_map, ...){
  dat <- dat %>% 
    tibble::rownames_to_column('ID') %>%
    mutate(Details=toupper(gs_map[ID,]$ID),
           Description=toupper(gs_map[ID,]$Desc),
           Details=NA,
           setSize=gs_map[ID,]$setSize,
           rank=NA,
           leading_edge=NA,
           FC=2^avg_log2FC,
           qval=p_val_adj) %>%
    mutate(ID=toupper(ID)) %>%
    dplyr::select(ID, Description, Details, setSize, FC, avg_log2FC, p_val,
                  p_val_adj, qval, rank, leading_edge) %>%
    rename_with(., ~c('NAME', 'GS.br..follow.link.to.MSigDB', 'GS.DETAILS', 
                      'SIZE', 'ES', 'NES', 'NOM.p.val', 'FDR.q.val', 
                      'FWER.p.val', 'RANK.AT.MAX', 'LEADING.EDGE'))
  datl <- split(dat, (dat$NES > 0))
  return(list("up"=datl[['TRUE']],
              "down"=datl[['FALSE']]))
}


# Creates a named vector mapping of msigdbr genesets to unique exact code
# e.g. c('GO:BP_Pathway_X'='GO_123456')
.mapGsToExactSource <- function(msig_lvls, species){
  if(missing(msig_lvls)){
    warning("'msig_lvls' undefined: Defaulting to H, C2, C5, and C8")
    msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                      'C2'=list('CP:REACTOME'),             # curated gene sets
                      'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                      'C8'=list(NULL)) 
  }
  if(missing(species)){
    warning("'species' undefined: Defaulting to Homo sapiens")
    species <- 'Homo sapiens'
  }
  
  main_obj <- lapply(names(msig_lvls), function(mlvl){
    sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
      print(paste0(">", mlvl, ":", sublvl, "..."))
      msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>% 
        as.data.frame %>%
        dplyr::select(gs_name, gs_exact_source) %>% 
        unique
      with(msig_ds, setNames(gs_exact_source, gs_name))
    })
  }) %>% unlist %>% 
    gsub("[ :-?]", "_", .)
}

