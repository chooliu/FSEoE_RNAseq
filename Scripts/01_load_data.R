# ==============================================================================
# 01_load_data.R
# import metadata, counts
# ==============================================================================




# sample mapping ===============================================================
# (deidentified ID --> clinical/demographic information)

metadata <-
  read_excel("./Data/SubjectInfo_RNAseq-042020.xlsx") %>%
  mutate(Phenotype = str_sub(RNASeqID, 0, 1) %>% as_factor()) %>%
  set_names(., names(.) %>% gsub("/| |-", "_", .))

# RNA extraction metadata (for RIN)
aliquots <-
  read_excel("./Data/SeqAliquotPlan-Menard-Katcher.xlsx")
aliquot_summary <-
  aliquots %>%
  group_by(`Subject ID`) %>%
  summarize(RIN = mean(c(RIN...7, RIN...12, RIN...17), na.rm = T))


metadata$RNASeqID %in% aliquot_summary$`Subject ID`
aliquot_summary$`Subject ID` %in% metadata$RNASeqID




# RNA-seq counts (processed by M.D.G.) =========================================

counts <-
  read_tsv("./Data/EoE2_raw_counts.txt")

# check same order
identical(metadata$RNASeqID %>% toupper,
          counts %>% names %>% .[-1] %>% toupper)

counts_matrix <- 
  counts %>%
  select(-1) %>%
  as.matrix %>%
  set_rownames(counts$Geneid)




# combine for final metadata ===================================================
metadata <-
  metadata %>%
  mutate(
    
    # distensibility labeled as "mean"
    distensibility = as.numeric(Mean),
    
    # alternative grouping: use of STCs
    PhenotypeSTC = case_when(
      Phenotype == "C" ~ "Cntl",
      Phenotype != "C" & Medication_STCs == "No" ~ "No",
      Phenotype != "C" & Medication_STCs == "Yes" ~ "Yes"),
    
    # name of groups in final manuscript
    Group = fct_recode(
      Phenotype, Control = "C", `non-fs EoE` = "E", `FS EoE` = "F")
  
  ) %>%
  left_join(aliquot_summary, c("RNASeqID" = "Subject ID"))



