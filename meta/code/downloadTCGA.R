########################
#### Download Files ####
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(reshape2)
library(ggplot2)

# Defining a list of Projects to download
projects <- c("TCGA-BRCA", "TCGA-GBM", "TCGA-OV", "TCGA-LUAD", "TCGA-UCEC",
              "TCGA-KIRC", "TCGA-HNSC", "TCGA-LGG", "TCGA-THCA", "TCGA-LUSC",
              "TCGA-PRAD", "TCGA-SKCM", "TCGA-COAD", "TCGA-STAD", "TCGA-BLCA",
              "TCGA-LIHC", "TCGA-CESC", "TCGA-KIRP", "TCGA-SARC", "TCGA-LAML", 
              "TCGA-PAAD", "TCGA-ESCA", "TCGA-PCPG", "TCGA-READ", "TCGA-TGCT",
              "TCGA-THYM", "TCGA-KICH", "TCGA-ACC", "TCGA-MESO", "TCGA-UVM",
              "TCGA-DLBC", "TCGA-UCS", "TCGA-CHOL")
workflow_type <- list("counts"="HTSeq - Counts", "fpkm"="HTSeq - FPKM-UQ")
# Query platform Illumina HiSeq with a list of barcode 
proj_queries <- lapply(projects, function(proj){
  print(proj)
  query <-GDCquery(project = proj, 
                   data.category = 'Transcriptome Profiling', 
                   data.type = 'Gene Expression Quantification', 
                   experimental.strategy = "RNA-Seq",
                   workflow.type = workflow_type$fpkm,
                   legacy = FALSE)
  
  GDCdownload(query, directory = "GDCdata")
  return(query)
})

proj_queries <- lapply(projects, function(proj){
  print(paste0(proj, " (", match(proj, projects), "/", length(projects), ")"))
  query <-GDCquery(project = proj,
                   data.category = 'Transcriptome Profiling',
                   data.type = 'Gene Expression Quantification',
                   experimental.strategy = "RNA-Seq",
                   workflow.type = workflow_type$counts,
                   legacy = FALSE)
  
  gdc_obj <- GDCprepare(query, directory = "GDCdata")
  gdc_matrix <- assay(gdc_obj,"HTSeq - Counts") # "HTSeq - FPKM-UQ"
  clinical_patient <- GDCquery_clinic(proj,"clinical")
  
  saveRDS(clinical_patient, file=file.path("GDCdata", "dat", paste0(proj, ".clinical.rds")))
  saveRDS(gdc_obj, file=file.path("GDCdata", "dat", paste0(proj, ".", gsub(" ", "", workflow_type$counts), ".obj.rds")))
  saveRDS(gdc_matrix, file=file.path("GDCdata", "dat", paste0(proj, ".", gsub(" ", "", workflow_type$counts), ".assay.rds")))
  return(NULL)
})
