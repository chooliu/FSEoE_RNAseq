# ==============================================================================
# 04_deseq_categorical.R
# compare categorical phenotype (cntl, non-fs EoE, FS EoE)
# ==============================================================================




# run naive versions of DESeq model ============================================
# includes no RUV correction model for sensitivity analysis / normalized counts

preruv_formula <-
  "~ Phenotype" %>%
  as.formula

preruv_deseq_obj <-
  DESeqDataSetFromMatrix(
    countData = counts_final,
    design = preruv_formula,
    colData = metadata_final
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions()

plotDispEsts(preruv_deseq_obj)

preruv_deseq_obj <- nbinomWaldTest(preruv_deseq_obj)

resultstable_noRUV_EvC_all <-
  results(preruv_deseq_obj, contrast = c("Phenotype", "E", "C"),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)
resultstable_noRUV_FvC_all <-
  results(preruv_deseq_obj, contrast = c("Phenotype", "F", "C"),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)
resultstable_noRUV_FvE_all <-
  results(preruv_deseq_obj, contrast = c("Phenotype", "F", "E"),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)

resultstable_noRUV_EvC_signif <- filter(resultstable_noRUV_EvC_all, padj < 0.05)
resultstable_noRUV_FvC_signif <- filter(resultstable_noRUV_FvC_all, padj < 0.05)
resultstable_noRUV_FvE_signif <- filter(resultstable_noRUV_FvE_all, padj < 0.05)




# run RUVr =====================================================================

# use rlog counts from above
preruv_rlog_counts <- rlog(preruv_deseq_obj)

# 27 x 17190
preruv_rlog_residuals <-
  preruv_rlog_counts %>%
  assay %>%
  apply(., 1,
        function(y) {
          tmp <-
            metadata_final %>%
            bind_cols(y = y)
          lm(update(preruv_formula, y ~ .),
             tmp) %>% 
            resid() }
  )

# elbow plot to check number of RUV components
check_ruv_kvals <-
  function(k) {
    ruv_tmp <-
      RUVr(x = assay(preruv_rlog_counts),
           k = k,
           residuals = preruv_rlog_residuals %>% t,
           isLog = T)$W
    adonis(preruv_rlog_residuals ~ ruv_tmp,
           method = "manhattan"
    ) %>%
      .$aov.tab %>%
      .$R2 %>%
      .[1]
  }

r2_by_k <- sapply(1:20, check_ruv_kvals)
ggplot(data = NULL, aes(x = 1:20, y = r2_by_k)) +
  geom_point() +
  geom_line() +
  xlab("# of RUVr Components (k)") +
  ylab("% of Residual Variance Explained\n(Multivariate PERMANOVA)") +
  theme_classic()

# run RUVr with k = 3 components
ruv_soln <-
  RUVr(x = assay(preruv_rlog_counts),
       k = 3,
       residuals = preruv_rlog_residuals %>% t,
       isLog = T)$W

# add to metadata
metadata_final$RUV1 <- ruv_soln[ , 1]
metadata_final$RUV2 <- ruv_soln[ , 2]
metadata_final$RUV3 <- ruv_soln[ , 3]

# check for correlations with clinical features
check_continuous_with_RUV <-
  function(y) {
    sapply(paste0("RUV", 1:3),
           function(ruv_factor) {
             cor(metadata_final[ , ruv_factor, drop = T],
                 y,
                 method = "spearman", use = "complete")
           })
  }

metadata_final %>%
  select_if(is.numeric) %>%
  summarize_all(check_continuous_with_RUV)




# run DESeq2 ===================================================================
# (final analysis presented in manuscript, adjusting for RUVr components)

# final DESeq2 object
test_formula <-
  update(preruv_formula, ~  . + RUV1 + RUV2 + RUV3) %>%
  as.formula

deseq_obj <-
  DESeqDataSetFromMatrix(
    countData = counts_final,
    design = test_formula,
    colData = metadata_final
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions()

plotDispEsts(preruv_deseq_obj)
rlog_counts <- rlog(deseq_obj)

# Wald testing / contrasts
deseq_obj <- nbinomWaldTest(deseq_obj)
resultstable_EvC_all <-
  results(deseq_obj, contrast = c("Phenotype", "E", "C"),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)
resultstable_FvC_all <-
  results(deseq_obj, contrast = c("Phenotype", "F", "C"),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)
resultstable_FvE_all <-
  results(deseq_obj, contrast = c("Phenotype", "F", "E"),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)

resultstable_EvC_signif <- filter(resultstable_EvC_all, padj < 0.05)
resultstable_FvC_signif <- filter(resultstable_FvC_all, padj < 0.05)
resultstable_FvE_signif <- filter(resultstable_FvE_all, padj < 0.05)




# custom fxn ===================================================================

# to clean environment in between statistical tests:
# used in next two scripts to make sure no leftovers from testing prev phenotype
clean_environment <- function() {
  rm(preruv_formula, preruv_deseq_obj, preruv_rlog_counts,
     ruv_soln, test_formula, deseq_obj, rlog_residuals,
     envir = .GlobalEnv
  )
}