# To work with the split LN/Tumor list of seurat objects
recode_list <- function(seu, newcol, grp=NULL, ...){
  if(!is.null(grp)){
    grp <- table(strsplit(unique(seu@meta.data[,grp]), split="_") %>% unlist)
    grp <- names(which.max(grp))
    print(paste0("grp = ", grp))
  }
  seu@meta.data[,newcol] <- seu@meta.data$seurat_clusters %>% 
    recode_map(., grp=grp, ...)
  return(seu)
}

# Relabel the orig.ident column
.relabelid <- function(x){
  if(any(grepl("Tumor", x))){
    x %>% 
      gsub("72h", "3d", .) %>%
      gsub("^Tumor", "B1_Tumor", .) %>%
      gsub("^LN", "B1_LN", .) %>%
      gsub("^ST2", "B2", .) %>%
      gsub("^CD45", "B2_CD45", .) %>%
      gsub("PBS$", "Un", .)
  } else {
    x %>% 
      gsub("72h", "3d", .) %>%
      gsub("^LN", "B1_LN", .) %>%
      gsub("^ST2", "B2", .) %>%
      gsub("^CD45", "B2_CD45", .) %>%
      gsub("PBS$", "Un", .)
  }
}

# Updated [Apr 24-2023] to annotate clusters and annotations
recode_map <- function(x, grp=NULL, anno=TRUE){
  if(grepl('^ln$', grp, ignore.case = T)){
    x %>% 
      recode("0"=ifelse(anno, 'Treg', "0"),
             "1"=ifelse(anno, "B", "1"),
             "2"=ifelse(anno, 'B', "1"),
             "3"=ifelse(anno, "DC_CCR7hi", "3"),
             "4"=ifelse(anno, 'Treg', "4"),
             "5"=ifelse(anno, 'Treg', "5"),
             "6"=ifelse(anno, 'Treg', "6"),
             "7"=ifelse(anno, 'Treg', "7"),
             "8"=ifelse(anno, "CD8_naive", "8"),
             "9"=ifelse(anno, 'Treg', "9"),
             "10"=ifelse(anno, 'B', "1"),
             "11"=ifelse(anno, 'CD4_CCR7hi', "11"),
             "12"=ifelse(anno, 'B', "12"),
             "13"=ifelse(anno, 'Treg', "13"),
             "14"=ifelse(anno, "CD8_proliferating", "14"),
             "15"=ifelse(anno, 'CD4_Tmem', "15"),
             "16"=ifelse(anno, 'MoMac', "16"),
             "17"=ifelse(anno, 'B', "17"),
             "18"=ifelse(anno, 'B', "1"),
             "19"=ifelse(anno, 'CD4_Tmem', "19"),
             "20"=ifelse(anno, 'CD4_Tmem', "20"),
             "21"=ifelse(anno, 'cDC1', "21"),
             "22"=ifelse(anno, 'B', "12"),
             '23'=ifelse(anno, 'DC_CCR7hi', "3"),
             "24"=ifelse(anno, 'B', "12"),
             "25"=ifelse(anno, "Treg", "4"),
             '26'=ifelse(anno, 'DC_CCR7hi', "3"))
  } else if(grepl('^tumor$', grp, ignore.case = T)){
    x %>% 
      recode("0"=ifelse(anno, "Treg", "0"),
             "1"=ifelse(anno, "Treg", "1"),
             "2"=ifelse(anno, "Treg", "1"),
             "3"=ifelse(anno, "Treg", "3"),
             "4"=ifelse(anno, "Treg", "4"),
             "5"=ifelse(anno, "Treg", "5"),
             "6"=ifelse(anno, "CD4_Tcell", "6"),
             "7"=ifelse(anno, "Treg", "7"),
             "8"=ifelse(anno, "MoMac", "8"),
             "9"=ifelse(anno, "MoMac", "8"),
             "10"=ifelse(anno, "CD8_Tcell", "10"),
             "11"=ifelse(anno, "Treg", "11"),
             "12"=ifelse(anno, "Patrolling_monocyte", "12"),
             "13"=ifelse(anno, "cDC2", "13"),
             "14"=ifelse(anno, "Treg", "1"),
             "15"=ifelse(anno, "CD8_Tcell", "15"),
             "16"=ifelse(anno, "ILC2", "16"),
             "17"=ifelse(anno, "DC", "17"),
             "18"=ifelse(anno, "CD4", "18"),
             "19"=ifelse(anno, "Eosinophil", "19"),
             "20"=ifelse(anno, "Macrophage", "20"),
             "21"=ifelse(anno, "Treg", "0"),
             "22"=ifelse(anno, "NK", "22"),
             "23"=ifelse(anno, "Treg", "3"))
  } else {
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
}

cd45_recode_map <- function(x, grp=NULL, anno=TRUE){
  if(grepl('^ln$', grp, ignore.case = T)){
    x %>% 
      recode("24"=ifelse(anno, 'DC', "0"),
             "20"=ifelse(anno, 'Monocyte_derived_macrophages', "1"),
             "27"=ifelse(anno, 'Monocyte_derived_macrophages', "1"),
             "11"=ifelse(anno, 'B', "2"),
             "31"=ifelse(anno, 'B', "2"),
             "1"=ifelse(anno, 'B', "3"),
             "23"=ifelse(anno, 'B', "3"),
             "13"=ifelse(anno, 'B', "3"),
             "0"=ifelse(anno, 'B', "3"),
             "7"=ifelse(anno, 'B', "3"),
             "28"=ifelse(anno, 'B', "3"),
             "10"=ifelse(anno, 'B', "4"),
             "9"=ifelse(anno, 'B', "4"),
             "6"=ifelse(anno, 'Remove', "-1"))
  } else if(grepl('^tumor$', grp, ignore.case = T)){
    x %>% 
      recode("0"=ifelse(anno, 'NK', "0"),
             "1"=ifelse(anno, 'CD8', "1"),
             "2"=ifelse(anno, 'CD8', "1"),
             "3"=ifelse(anno, 'CD8', "1"),
             "4"=ifelse(anno, 'Treg', "4"),
             "5"=ifelse(anno, 'Monocyte_derived_macrophage', "5"),
             "6"=ifelse(anno, 'Macrophage', "6"),
             "7"=ifelse(anno, 'CD4', "7"),
             "8"=ifelse(anno, 'Monocyte_derived_macrophage', "5"),
             '9'=ifelse(anno, 'CD8', "9"),
             '10'=ifelse(anno, 'DC', "10"),
             '11'=ifelse(anno, 'CD8', "1"),
             '12'=ifelse(anno, 'Tgd', "12"),
             '13'=ifelse(anno, 'CD8_proliferating', "13"),
             '14'=ifelse(anno, 'DC_Cd4pos_Cd8pos', "14"),
             '15'=ifelse(anno, 'ILC3', "15"),
             '16'=ifelse(anno, 'NK', "0"),
             '17'=ifelse(anno, 'Treg', "17"),
             '18'=ifelse(anno, 'CD8', "1"),
             '19'=ifelse(anno, 'Treg', "4"),
             '20'=ifelse(anno, 'Treg_proliferating', "20"),
             '21'=ifelse(anno, 'DC', "21"),
             '22'=ifelse(anno, 'NK_proliferating', "22"),
             '23'=ifelse(anno, 'B', "23"),
             '24'=ifelse(anno, 'CD3_DN_Naive', "24"),
             '25'=ifelse(anno, 'Macrophage', "25"),
             '26'=ifelse(anno, 'Monocyte', "26"))
  } 
}

# ---- Tcells ----
tcell_recode_map <- function(x, grp=NULL, day=NULL, anno=TRUE){
  if(grepl('^ln$', grp, ignore.case = T)){
    if(grepl("Day3", day, ignore.case=T)){
      x %>% 
        recode("1"=ifelse(anno, 'CD4_naive', "0"),
               "10"=ifelse(anno, 'CD4_naive', "0"),
               "2"=ifelse(anno, 'CD8_naive', "1"),
               "3"=ifelse(anno, 'CD8_naive', "1"),
               "4"=ifelse(anno, 'CD8_cycling', "2"),
               "5"=ifelse(anno, 'CD8_cycling', "2"),
               "6"=ifelse(anno, 'CD8_cycling', "2"),
               "8"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "7"=ifelse(anno, 'CD4_Isc.Tfh', "4"),
               "0"=ifelse(anno, 'CD4_Isc.Tfh', "4"),
               "9"=ifelse(anno, 'CD4_Isc.Tfh', "4"))
    } else if (grepl("Day7", day, ignore.case=T)){
      x %>% 
        recode("5"=ifelse(anno, 'CD4_naive', "0"),
               "3"=ifelse(anno, 'CD8_naive', "1"),
               "4"=ifelse(anno, 'CD8_naive', "1"),
               "10"=ifelse(anno, 'CD8_cycling', "2"),
               "6"=ifelse(anno, 'CD8_cycling', "2"),
               "7"=ifelse(anno, 'CD8_cycling', "2"),
               "9"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "0"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "11"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "1"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "2"=ifelse(anno, 'CD4_Isc.Tfh', "4"),
               "8"=ifelse(anno, 'CD4_Isc.Tfh', "4"))
    }
  } else if(grepl('^tumor$', grp, ignore.case = T)){
    if(grepl("Day3", day, ignore.case=T)){
      x %>%
        recode("2"=ifelse(anno, 'CD8_cycling', "2"),
               "4"=ifelse(anno, 'CD8_cycling', "2"),
               "7"=ifelse(anno, 'CD8_cycling', "2"),
               "8"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "6"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "3"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "5"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "9"=ifelse(anno, 'CD8_effector', "4"),
               "1"=ifelse(anno, 'CD8_effector', "4"),
               "0"=ifelse(anno, 'CD8_effector', "4"))
    } else if(grepl("Day7", day, ignore.case=T)){
      x %>%
        recode("9"=ifelse(anno, 'CD8_cycling', "2"),
               "11"=ifelse(anno, 'CD8_cycling', "2"),
               "5"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "1"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "4"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "7"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "8"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "3"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "10"=ifelse(anno, 'CD4_CCR7pos_TCF7pos', "3"),
               "0"=ifelse(anno, 'CD8_effector', "4"),
               "2"=ifelse(anno, 'CD8_effector', "4"),
               "12"=ifelse(anno, 'CD8_effector', "4"),
               "13"=ifelse(anno, 'CD8_effector', "4"),
               "6"=ifelse(anno, 'Remove', "-1"),
               "14"=ifelse(anno, 'Remove', "-1"))
    }
  }
}

cd45_tcell_recode_map <- function(x, grp=NULL, day=NULL, anno=TRUE){
  if(grepl('^ln$', grp, ignore.case = T)){
    x %>% 
      recode("9"=ifelse(anno, 'Remove', "-1"),
             "7"=ifelse(anno, 'Remove', "-1"),
             "20"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', "2"),
             "3"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', "2"),
             "4"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', "2"),
             "2"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', "2"),
             "19"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', "2"),
             "23"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_002', "2"),
             "17"=ifelse(anno, 'CD8_memory', "3"),
             "14"=ifelse(anno, 'CD8_memory', "3"),
             "13"=ifelse(anno, 'CD8_memory', "3"),
             "24"=ifelse(anno, 'CD8_memory', "3"),
             '0'=ifelse(anno, 'CD8_naive_memory', "4"),
             '6'=ifelse(anno, 'CD8_naive_memory', "4"),
             '1'=ifelse(anno, 'CD8_naive_memory', "4"),
             '5'=ifelse(anno, 'CD8_naive_memory', "4"),
             '12'=ifelse(anno, 'CD4_naive', "5"),
             '16'=ifelse(anno, 'CD4_naive', "5"),
             '10'=ifelse(anno, 'CD8_naive', "6"),
             '21'=ifelse(anno, 'Remove', "-1"),
             '22'=ifelse(anno, 'Remove', "-1"),
             '25'=ifelse(anno, 'Remove', "-1"),
             '8'=ifelse(anno, 'TReg', "-2"),
             '11'=ifelse(anno, 'TReg', "-2"),
             '15'=ifelse(anno, 'TReg', "-2"),
             '18'=ifelse(anno, 'TReg', "-2"))
  } else if(grepl('^tumor$', grp, ignore.case = T)){
    x %>% 
      recode("4"=ifelse(anno, 'CD3_DN_activated', "0"),
             "16"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_001', "1"),
             "5"=ifelse(anno, 'CD4_intermediate_TCF7pos_CCR7pos_TCF7pos_001', "1"),
             "11"=ifelse(anno, 'CD8_effector_cycling', "2"),
             "17"=ifelse(anno, 'CD8_effector_cycling', "2"),
             "20"=ifelse(anno, 'CD8_effector_cycling', "2"),
             "8"=ifelse(anno, 'CD8_effector', "3"),
             "3"=ifelse(anno, 'CD8_effector', "3"),
             "2"=ifelse(anno, 'CD8_effector', "3"),
             "10"=ifelse(anno, 'CD8_memory_precursor_effector', "4"),
             "7"=ifelse(anno, 'CD8_memory_precursor_effector', "4"),
             "14"=ifelse(anno, 'CD8_memory_precursor_effector', "4"),
             "13"=ifelse(anno, 'CD8_memory', "5"),
             "6"=ifelse(anno, 'CD8_memory', "5"),
             "1"=ifelse(anno, 'CD8_memory', "5"),
             "22"=ifelse(anno, 'CD8_memory', "5"),
             "15"=ifelse(anno, 'CD8_memory2', "6"),
             "9"=ifelse(anno, 'CD8_memory2', "6"),
             "21"=ifelse(anno, 'TReg', "-2"),
             "12"=ifelse(anno, 'TReg', "-2"),
             "0"=ifelse(anno, 'TReg', "-2"),
             "19"=ifelse(anno, 'TReg', "-2"),
             "18"=ifelse(anno, 'TReg', "-2"),
             )
  }
}


# ---- TRegs ----
# Finalized [June 14-2023] to annotate clusters and annotations
treg_recode_map <- function(x, grp=NULL, day=NULL, anno=TRUE){
  if(grepl('^ln$', grp, ignore.case = T)){
    if(grepl("Day3", day, ignore.case=T)){
      x %>% 
        recode("0"=ifelse(anno, 'NLTlike.Effector.Central', "0"),
               "1"=ifelse(anno, 'NLTlike.Effector.Central', "0"),
               "13"=ifelse(anno, 'NLTlike.Effector.Central', "0"),
               "14"=ifelse(anno, 'NLTlike.Effector.Central', "0"),
               "2"=ifelse(anno, 'NLTlike.Effector.Central', "1"),
               "15"=ifelse(anno, 'NLTlike.Effector.Central', "1"),
               "6"=ifelse(anno, 'Central', "2"),
               "9"=ifelse(anno, 'NLTlike.STAT1', "3"),
               "7"=ifelse(anno, 'Effector', "4"),
               "5"=ifelse(anno, 'Effector.NLTlike_cycling', "5"),
               "10"=ifelse(anno, 'Effector.NLTlike_cycling', "5"),
               "3"=ifelse(anno, 'NLTlike.Effector_cycling', "6"),
               "4"=ifelse(anno, 'NLTlike.Effector_cycling', "6"),
               "8"=ifelse(anno, 'NLTlike_cycling', "7"),
               "12"=ifelse(anno, 'NLTlike_cycling', "7"),
               "16"=ifelse(anno, 'NLTlike_cycling', "7"),
               "11"=ifelse(anno, 'Remove', "-1"))
    } else if (grepl("Day7", day, ignore.case=T)){
      x %>% 
        recode("2"=ifelse(anno, 'NLTlike.Effector_cycling', "0"),
               "6"=ifelse(anno, 'NLTlike.Effector_cycling', "0"),
               "11"=ifelse(anno, 'NLTlike_cycling', "1"),
               "17"=ifelse(anno, 'NLTlike_cycling', "1"),
               "8"=ifelse(anno, 'Effector.NLTlike_cycling', "2"),
               "4"=ifelse(anno, 'Effector', "3"),
               "0"=ifelse(anno, 'NLT.NLTlike', "4"),
               "3"=ifelse(anno, 'NLT.NLTlike', "4"),
               "20"=ifelse(anno, 'NLT.NLTlike', "4"),
               "1"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "7"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "9"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "12"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "13"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "15"=ifelse(anno, 'Central', "6"),
               "14"=ifelse(anno, 'NLT', "7"),
               "16"=ifelse(anno, 'NLT', "7"),
               "10"=ifelse(anno, 'NLTlike.Effector.Central', "10"),
               "5"=ifelse(anno, 'Remove', "-1"),
               "18"=ifelse(anno, 'Remove', "-1"),
               "19"=ifelse(anno, 'Remove', "-1"))
    }
  } else if(grepl('^tumor$', grp, ignore.case = T)){
    if(grepl("Day3", day, ignore.case=T)){
      x %>%
        recode("10"=ifelse(anno, 'NLTlike_cycling', "1"),
               "14"=ifelse(anno, 'NLTlike_cycling', "2"),
               "13"=ifelse(anno, 'NLTlike_cycling', "3"),
               "7"=ifelse(anno, 'Effector', "4"),
               "9"=ifelse(anno, 'Effector', "4"),
               "8"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "6"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "1"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "2"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "4"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "5"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "11"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "15"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "16"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "0"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "3"=ifelse(anno, 'NLTlike.STAT1', "6"),
               "12"=ifelse(anno, 'NLTlike.STAT1', "6"))
    } else if(grepl("Day7", day, ignore.case=T)){
      x %>%
        recode("12"=ifelse(anno, 'NLTlike_cycling', "1"),
               "20"=ifelse(anno, 'NLTlike_cycling', "1"),
               "14"=ifelse(anno, 'NLTlike_cycling', "2"),
               "18"=ifelse(anno, 'NLTlike_cycling', "3"),
               "19"=ifelse(anno, 'Effector', "4"),
               "10"=ifelse(anno, 'Effector', "4"),
               "3"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "15"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "5"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "7"=ifelse(anno, 'NLTlike.Effector.Central', "5"),
               "9"=ifelse(anno, 'NLTlike.Effector.Central', "6"),
               "17"=ifelse(anno, 'NLTlike.Effector.Central', "6"),
               "2"=ifelse(anno, 'NLTlike.STAT1', "7"),
               "6"=ifelse(anno, 'NLTlike.STAT1', "7"),
               "8"=ifelse(anno, 'NLT', "8"),
               "13"=ifelse(anno, 'NLT.NLTlike', "9"),
               "16"=ifelse(anno, 'NLT.NLTlike', "9"),
               "4"=ifelse(anno, 'NLT.NLTlike', "9"),
               "0"=ifelse(anno, 'NLT.NLTlike', "9"),
               "22"=ifelse(anno, 'NLT.NLTlike', "9"),
               "1"=ifelse(anno, 'NLT.NLTlike', "9"),
               "11"=ifelse(anno, 'NLT.NLTlike', "9"),
               "21"=ifelse(anno, 'Remove', "-1"))
    }
  }
}

# Finalized [June 14-2023] to annotate clusters and annotations
cd45_treg_recode_map <- function(x, grp=NULL, day=NULL, anno=TRUE){
  if(grepl('^ln$', grp, ignore.case = T)){
    x %>% 
      recode("3"=ifelse(anno, 'NLTlike.Central', "0"),
             "12"=ifelse(anno, 'NLTlike.Central', "0"),
             "14"=ifelse(anno, 'NLTlike.Central', "0"),
             "1"=ifelse(anno, 'NLTlike', "1"),
             "8"=ifelse(anno, 'NLTlike_cycling', "2"),
             "0"=ifelse(anno, 'Central', "3"),
             "2"=ifelse(anno, 'Central', "3"),
             "5"=ifelse(anno, 'Central', "3"),
             "10"=ifelse(anno, 'Effector', "4"),
             "13"=ifelse(anno, 'Effector', "4"),
             "4"=ifelse(anno, 'Remove', "-1"),
             "6"=ifelse(anno, 'Remove', "-1"),
             "7"=ifelse(anno, 'Remove', "-1"),
             "9"=ifelse(anno, 'Remove', "-1"),
             "11"=ifelse(anno, 'Remove', "-1"),
             "15"=ifelse(anno, 'Remove', "-1"),
             "16"=ifelse(anno, 'Remove', "-1"))
  } else if(grepl('^tumor$', grp, ignore.case = T)){
    x %>%
      recode("5"=ifelse(anno, 'NLTlike', "0"),
             "12"=ifelse(anno, 'NLTlike', "0"),
             "10"=ifelse(anno, 'NLTlike_cycling', "1"),
             "8"=ifelse(anno, 'NLTlike_cycling', "1"),
             "2"=ifelse(anno, 'NLTlike.NLT', "2"),
             "11"=ifelse(anno, 'NLTlike.NLT', "2"),
             "13"=ifelse(anno, 'NLTlike.NLT', "2"),
             "0"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "1"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "3"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "4"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "6"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "7"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "9"=ifelse(anno, 'NLTlike.NLT.Effector', "3"),
             "14"=ifelse(anno, 'NLTlike.NLT.Effector', "3"))
  }
}