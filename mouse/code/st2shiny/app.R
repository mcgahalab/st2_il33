library(dplyr)
library(shiny)
library(ggplot2)
library(reshape2)
library(DESeq2)
y <- data.frame(colnames(dat$cnt_mat), dat$condition, colnames(vst_assay_goi))
View(x)
  
# setwd("~/git/st2_il33/mouse/code/st2shiny")
# write.table(dat$cnt_mat, file=file.path("data", "rsem_mat.tsv"),
#             sep="\t", col.names = T, row.names = T)
# write.table(meta, file=file.path("data", "meta.tsv"),
#             sep="\t", col.names = T, row.names = F)
# meta <- dat$condition %>%
#   strsplit(., split=delims) %>%
#   do.call(rbind, .) %>% as.data.frame %>%
#   rename_with(., ~c('LNstatus', 'Treatment', 'Celltype')) %>%
#   mutate(Treatment=recode(Treatment,
#                           "PBS"="Untreated",
#                           "CIS"="Cisplatin")) %>%
#   mutate(ID=colids,
#          Sample=colnames(vst_assay_goi))


dat <- readRDS(file.path("data/res.rds"))
# meta <- readRDS(file.path("data/meta.rds"))

#### UI ####
ui <- fluidPage(
  titlePanel("MSM/SSM Heatmap Plotter"),
  sidebarLayout(
    sidebarPanel(
      uiOutput('extraParams'),
      
      fileInput("cntmat", "Count matrix (ENS-genes by Samples)",
                multiple = FALSE, accept = ".tsv", placeholder = "count_matrix.tsv"),
      fileInput("metamat", "Metadata data-frame",
                multiple=FALSE, accept = ".tsv", placeholder = "meta.tsv"),
      
      radioButtons("genome", "Genome", choices = c('Mouse', 'Human'), selected = 'Mouse'),
      
      checkboxInput('vst', 'VST', value=T),
      checkboxInput('inferanno', 'Infer metadata', value=F),
      
      textInput('meta_sample_col', 'Sample column ID for metadata',
                value=c('Sample')),
      textInput('groupid', 'Metadata column specifying groups to merge samples across',
                value=c('ID')),
      textInput('goi', 'Genes of Interest',
                value=c('Il10, Cd274, Pdcd1lg2, Tgfb3, Ccl17, Mrc1')),
      
      checkboxInput(inputId = 'merge_replicates',
                    label = 'Whether to merge biological replicates',
                    value=FALSE),
      checkboxGroupInput(inputId = "keepgrp",
                         label="Identify which samples to keep in the plot",
                         choices = list("SSM"="SSM",
                                        "MSM"="MSM",
                                        "LN"="LN",
                                        "TDLN"="TDLN",
                                        "PBS"="Untreated",
                                        "CIS"="Cisplatin"),
                         select=c('SSM', 'MSM', 'LN', 'TDLN', 
                                  'Untreated', 'Cisplatin'))
    ),
    mainPanel(
      plotOutput(outputId = "msmSsmHeatmap")
    )
  )
)

#### Server #####
server <- function(input, output){
  # Hardcoded parameters [ Can Change ] 
  annotation_colors=list(
    LNstatus = c('TDLN'='#7fbc41', 'LN'='#de77ae'), # pink, green
    Treatment = c('Cisplatin'='#006d2c', 'Untreated'='#a1d99b'), # light-green, dark-green
    Celltype=c('SSM'='white', 'MSM'='black')  
  )
  
  # Hardcoded parameters [ Don't Change ] 
  celltype_order <- c('SSM', 'MSM')
  lnstatus_order <- c('LN', 'TDLN')
  treatment_order <- c('Untreated', 'Cisplatin')
  
  #... Load the Ens2Symbol Mapping ----
  mapInput <- reactive({
    genome_gse <- switch(input$genome,
                         Mouse={library(org.Mm.eg.db); org.Mm.eg.db},
                         Human={library(org.Hs.eg.db); org.Hs.eg.db})
    txby <- keys(genome_gse, 'ENSEMBL')
    gene_ids <- mapIds(genome_gse, keys=txby, column='SYMBOL',
                       keytype='ENSEMBL', multiVals="first")
    gene_ids
  })
  
  #... Load in (and transform) the count matrix ----
  dataInput <- reactive({
    ens2sym_ids <- mapInput()
    inFile <- input$cntmat
    if(is.null(inFile)) return(NULL)
    cnt_mat <- read.table(file=inFile$datapath, header=T, stringsAsFactors = F,
                          check.names = F, row.names=1)
    
    if(input$vst){
      vst_assay <- vst(as.matrix(cnt_mat), blind=T)
    } else {
      vst_assay <- as.matrix(cnt_mat)
    }
    vst_assay_goi <- vst_assay %>% 
      as.data.frame %>%
      tibble::rownames_to_column(., "id") %>%
      mutate(id=ens2sym_ids[id]) %>%
      filter(!duplicated(id),
             !is.na(id)) %>%
      tibble::column_to_rownames(., "id")
    
    vst_assay_goi
  })
  
  #... Infer/Load the metadata annotations ----
  metaInput <- reactive({
    if(input$inferanno){
      return(NULL)
      # delims <- gsub("(.)", "\\1|", input$annosep) %>%
      #   gsub(".$", "", .) %>% gsub("\\.", "\\\\.", .)
      # 
      # meta_headers <- gsub(" " , "", strsplit(input$metaheaders, ",")[[1]])
      # 
      # colids <- colnames(dataInput())
      # annotation_col <- colids %>% 
      #   strsplit(., split=delims) %>% 
      #   do.call(rbind, .) %>% as.data.frame
      # colnames(annotation_col) <- if(length(meta_headers) == ncol(annotation_col)){
      #   meta_headers 
      # } else {
      #   letters[1:ncol(annotation_col)]
      # }
      #   # mutate(Celltype = factor(Celltype, levels=celltype_order),
      #   #        LNstatus = factor(LNstatus, levels=lnstatus_order),
      #   #        Treatment = factor(Treatment, levels=treatment_order))
    } else {
      inData <- input$metamat
      if(is.null(inData)) return(NULL)
      annotation_col <- read.table(inData$datapath, header=T, stringsAsFactors = F,
                 check.names = F) %>% 
        tibble::column_to_rownames(., input$meta_sample_col)
    }
    
    annotation_col
  })
  
  #... Take the Metadata and infer the parameters to tune ----
  output$extraParams <- renderUI({
    meta <- metaInput()
    if(is.null(meta)) return(NULL)

    tagList(
      selectInput("metadata_col", "Metadata Column ID", choices=colnames(meta)),
      selectInput("metadata_ord", "Metadata Order", choices=NULL)
    )
  })
  
  output$msmSsmHeatmap <- renderPlot({
    # Read in the data matrix
    vst_assay_goi <- dataInput() 
    if(is.null(vst_assay_goi)) return(NULL)
    
    # Read in the metadata
    annotation_col <- metaInput()
    if(is.null(annotation_col)) return(NULL)
    
    # Order the metadata and expression matrix
    annotation_col <- annotation_col[match(colnames(vst_assay_goi), 
                                           rownames(annotation_col)),]
    ordidx <- order(annotation_col[,input$groupid])
    annotation_col <- annotation_col[ordidx,]
    vst_assay_goi <- vst_assay_goi[,ordidx]

    # Reactive values
    goi <- gsub(" " , "", strsplit(input$goi, ",")[[1]])
    if(any(!goi %in% rownames(vst_assay_goi))){
      warning(paste0("Genes not in your GOI: ", 
                     paste(goi[!goi %in% rownames(vst_assay_goi)], collapse="\n")), "\n")
    }
    
    
    # Merge matrices and decided whether to merge or separate samples
    vst_grp_goi_list <- split(as.data.frame(t(vst_assay_goi[goi,])), 
                              as.character(annotation_col[,input$groupid]))
    row_ord <- seq_along(goi) # Specify an order if you want to keep that gene-order
    
    if(input$merge_replicates){
      # Merge biological replicates into one group average
      vst_grp_goi <- vst_grp_goi_list %>%
        lapply(., colMeans) %>%
        do.call(rbind, .) %>% t
    } else {
      # Keep biological replicates separate and annotate accordingly
      vst_grp_goi <- vst_grp_goi_list %>% 
        do.call(rbind, .) %>% 
        t %>% as.data.frame %>%
        rename_with(., ~make.unique(gsub("\\..*", "", .)))
    }
    
    
    # Generate metadata order
    updateMetaCol <- reactive({
      req(input$metadata_col)
    })
    observeEvent(updateMetaCol(), {
      choices <- unique(annotation_col[,input$metadata_col])
      updateSelectInput(inputId = "metadata_ord", choices = choices) 
    })
    
    # Update the annotation plotting order based on the selection
    updateMetaOrder <- reactive({
      req(input$metadata_ord)
      annotation_col[,input$metadata_col] <- factor(annotation_col[,input$metadata_col])
      annotation_col[,input$metadata_col] <- relevel(annotation_col[,input$metadata_col], 
                                                     input$metadata_ord)
      order_idx <- do.call("order", annotation_col[colnames(annotation_col)])
      return(order_idx)
    })
    order_idx <- updateMetaOrder()
    
    # If no row order is given, do some custom ordering based on difference
    heatmaps <- list("matrix"=vst_grp_goi[row_ord, order_idx], 
                     "annotation_col"=annotation_col[order_idx,])
    
    # Filter out celltypes
    annotation_col <- heatmaps$annotation_col
    keep_idx <- apply(annotation_col[,1:3], 2, function(i){
      which(i %in% input$keepgrp)
    }) %>%
      Reduce(function(x,y) intersect(x,y), .)
    if(length(keep_idx) == 0 | is.null(keep_idx)) keep_idx <- c(1:nrow(annotation_col))
    
    # Z-scale the count matrix
    scaled_grp_goi <- t(apply(heatmaps$matrix[,keep_idx], 1, scale))
    
    # Prepare final ordering for heatmap
    annotation_col <- annotation_col[keep_idx,]
    annotation_col$split <- with(annotation_col, paste0(Celltype, LNstatus))
    annotation_col$split <- factor(annotation_col$split, levels=unique(annotation_col$split))
    
    # Plot heatmap
    gg_phm <- ComplexHeatmap::pheatmap(as.matrix(scaled_grp_goi), scale = 'none',
                                       color = colorRampPalette(c("#053061", "#4393c3", "#f7f7f7", 
                                                                  "#d6604d", "#67001f"))(20),
                                       annotation_col = annotation_col,
                                       column_split= annotation_col$split,
                                       annotation_colors = annotation_colors,
                                       fontsize_row = 14,show_colnames = FALSE,
                                       use_raster=FALSE, cluster_cols = FALSE, cluster_rows = FALSE
    ) %>% 
      ggplotify::as.ggplot(.) #+
      #theme_minimal(base_family = "Arial")
    
    if(FALSE){
      keep_idx <<- keep_idx
      heatmaps <<- heatmaps
      annotation_col <<- annotation_col 
      scaled_grp_goi <<- scaled_grp_goi
      vst_assay_goi <<- vst_assay_goi
      order_idx <<- order_idx
    }
    
    gg_phm
  })
}

shinyApp(ui = ui, server = server)

# rm(annotation_col, scaled_grp_goi, vst_assay_goi, order_idx, heatmaps, keep_idx)
# Il10, Tgfb3, Mrc1, Ido1, Ido2, Cyp1a1, Cyp1b1, Asns, Ddit3, Il18bp, Pdcd1lg2, Cd274, Ccl17, Il1rn, Vsir, Il4, Tnfaip3, Havcr2, Lilrb4a, Ahr, Il1b, Ccr2, Il17f, Tnfsf4, Ccl3, Ccl4, Ccl2, Nlrp3, Il6, Cd80, Il12ra, Il12rb1, Il2, Trem1, Cxcl10, Il17ra, Il12b
