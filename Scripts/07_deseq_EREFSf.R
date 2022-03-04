# ==============================================================================
# 07_deseq_EREFSf.R
# DE with continuous phenotype #2; for comments, see script 04 (similar code)
# ==============================================================================




# run naive versions of DESeq model ============================================

clean_environment()

preruv_formula <-
  "~ EREFf" %>%
  as.formula

preruv_deseq_obj <-
  DESeqDataSetFromMatrix(
    countData = counts_final,
    design = ~1,
    colData = metadata_final
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions()

plotDispEsts(preruv_deseq_obj)

preruv_rlog_counts <- rlog(preruv_deseq_obj)

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


# run RUVr =====================================================================

ruv_soln <-
  RUVr(x = assay(preruv_rlog_counts),
       k = 3,
       residuals = preruv_rlog_residuals %>% t,
       isLog = T)$W

metadata_final$RUV1 <- ruv_soln[ , 1]
metadata_final$RUV2 <- ruv_soln[ , 2]
metadata_final$RUV3 <- ruv_soln[ , 3]




# run DESeq2 ===================================================================

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

deseq_obj <- nbinomWaldTest(deseq_obj)

resultstable_EREFf_all <-
  results(deseq_obj, name = "EREFf",
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue)

resultstable_EREFf_signif <-
  filter(resultstable_EREFf_all, padj < 0.05)

rlog_residuals <-
  rlog_counts %>%
  assay %>%
  apply(., 1,
        function(y) {
          tmp <-
            metadata_final %>%
            bind_cols(y = y)
          lm(update(test_formula, y ~ . - EREFf),
             tmp) %>% 
            resid() }
  )




# plotting DE results (Fig 1) ==================================================

rlog_residuals %>%
  as_tibble() %>%
  bind_cols(ID = metadata_final$RNASeqID, .) %>%
  write_tsv("20211228_resids_EREFf.tsv")

plot_by_EREFf <-
  function(genename, residuals = T, padding = NULL) {
    
    if (residuals) {
      Expression <- rlog_residuals[ , genename ]
      yaxis_label <- "Residualized Expression" 
    } else {
      Expression <- assay(rlog_counts)[genename, ]
      yaxis_label <- "Expression (rlog transform)"
    }
    
    yaxisrange <- 
      Expression %>% range(na.rm = T)
    if( is.null(padding) ) {
      padding <- Expression %>% sd(na.rm = T) }
    
    
    metadata_final %>%
      bind_cols(Expression) %>%
      ggplot(data = .,
             aes(x = EREFf, y = Expression,
                 color = Group, shape = Group)) +
      geom_point(alpha = 0.8, size = 2) +
      theme_few() +
      geom_smooth(aes(group = 1),
                  method = "lm", se = F, color = "black", alpha = 0.8) +
      scale_x_continuous("EREFf Severity") +
      scale_y_continuous(paste0(genename, " Expression"),
                         breaks = pretty_breaks(n = 6), 
                         limits = c(yaxisrange[1] - padding/3,
                                    yaxisrange[2] + padding/3)) + 
      theme(legend.position = "right") +
      scale_color_manual(values = palette_color_phenotype) +
      scale_shape_manual(values = palette_shape_phenotype)
  }

plot_by_EREFf("CXCL1", residuals = T)
ggsave("Fig4_CXCL1_EREFSf.pdf", width = 5, height = 3.5)
plot_by_EREFf("NEFL", residuals = T)
ggsave("Fig4_NEFL_EREFSf.pdf", width = 5, height = 3.5)





# check euler for all three phenotypes =========================================

set.seed(1234)
euler_fit <-
  euler(
    list(FvC = resultstable_FvC_signif$row,
         Dysphagia = resultstable_dysphagia_signif$row,
         EREFf = resultstable_EREFf_signif$row),
    shape = "ellipse", 
  )

plot(euler_fit, labels = T, quantities = T, cex = 0.6)





