# Returns a mapping from key to value (e.g. SYMBOL -> ENSEMBL == map$SYMBOL$ENSEMBL)
geneMap <- function(version='v2', ...){
  if(version=='v1'){
    geneMap.v1(...)
  } else if(version=='v2'){
    geneMap.v2(...)
  }
}

geneMap.v1 <- function(species='Homo sapiens'){
  genome_gse <- switch(species,
                       "Homo sapiens"={
                         library(org.Hs.eg.db)
                         org.Hs.eg.db
                       },
                       "Mus musculus"={
                         library(org.Mm.eg.db)
                         org.Mm.eg.db
                       })
  db <- c('SYMBOL', 'ENSEMBL', 'ENTREZID')
  db <- setNames(db,db)
  
  gene_map <- lapply(db, function(db_i){
    lapply(db[-match(db_i, db)], function(db_j){
      txby <- keys(genome_gse, db_i)
      i2j_ids <- mapIds(genome_gse, keys=txby, column=db_j,
                        keytype=db_i, multiVals="first")
      i2j_ids
    })
  })
  
  return(gene_map)
}

geneMap.v2 <- function(species='Homo sapiens'){
  gene_info <- switch(species,
                       "Homo sapiens"={
                         library(EnsDb.Hsapiens.v86)
                         ensembldb::genes(EnsDb.Hsapiens.v86, 
                                                       return.type='data.frame')
                       },
                       "Mus musculus"={
                         library(EnsDb.Mmusculus.v79)
                         ensembldb::genes(EnsDb.Mmusculus.v79 ,
                                          return.type='data.frame')
                       })
  entrezid <- unlist(gene_info$entrezid)
  reps <- sapply(gene_info$entrezid, length)
  gene_info2 <- data.frame("SYMBOL"=rep(gene_info$gene_name, reps),
                           "ENSEMBL"=rep(gene_info$gene_id, reps),
                           "gene_biotype"=rep(gene_info$gene_biotype, reps),
                           'ENTREZID'=entrezid)
  
  db <- c('SYMBOL', 'ENSEMBL', 'ENTREZID')
  db <- setNames(db,db)
  dbextra <- c('gene_biotype'='gene_biotype')
  
  
  gene_map <- lapply(db, function(db_i){
    lapply(c(db[-match(db_i, db)], dbextra), function(db_j){
      gmap <- setNames(gene_info2[,db_j], gene_info2[,db_i])
      idx <- which(duplicated(names(gmap)) | is.na(names(gmap)))
      gmap[-idx]
    })
  })
  
  return(gene_map)
}