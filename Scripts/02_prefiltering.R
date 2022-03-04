# ==============================================================================
# 02_prefiltering.R
# QC participants, genes
# ==============================================================================




# generate filter criteria for participants & genes ============================

# remove one participant with acute esophagitis* 
# in later EDA appears as clear outlier
# (* abnormal esophageal biopsy, but not EoE)
filter_subjects <-
  metadata$RNASeqID != "C106"

# genes detected in (i.e., non-zero) at least 2/3rds of dataset
filter_genes <-
  counts_matrix %>%
  apply(., 1, function(x) { sum(x == 0) < 9 })




# apply filtering for final datasets ===========================================

# 17190 genes x 27 participants
metadata_final <-
  metadata[filter_subjects, ]
counts_final <-
  counts_matrix[filter_genes, filter_subjects]
  
# for plotting/sensitivity analysis, optional filter to examine cases only
filter_subjects_nocntl <-
  metadata_final$Phenotype != "C"

# 17190 x 20
metadata_final_nocntl <- 
  metadata_final[filter_subjects_nocntl, ]
counts_final_nocntl <- 
  counts_final[ , filter_subjects_nocntl ]


