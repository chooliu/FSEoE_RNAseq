# ==============================================================================
# 05_figures_Venn_PCA_heatmap_DE.R
# some major figures / processing of primary analysis (DEGs by phenotype)
# ==============================================================================




# Venn diagram of DE results (Fig 2) ===========================================

# use eulerr package to generate Euler diagram ---------------------------------
# (Venn diagram but with scaled circle areas)
set.seed(1234)
euler_fit <-
  euler(
    list(FvC = resultstable_FvC_signif$row,
         EvC = resultstable_EvC_signif$row,
         FvE = resultstable_FvE_signif$row),
    shape = "circle"
  )

pdf(file = "Fig2_VennDiagram.pdf", width = 5, height = 5)
plot(euler_fit, labels = F, quantities = T, cex = 0.6)
dev.off()

# check counts manually
intersect(resultstable_FvC_signif$row, resultstable_EvC_signif$row) %>% length
setdiff(resultstable_EvC_signif$row, resultstable_FvC_signif$row) %>% length
setdiff(resultstable_FvC_signif$row, resultstable_EvC_signif$row) %>% length
setdiff(resultstable_EvC_signif$row, resultstable_FvC_signif$row) %>% length




# check directionality ---------------------------------------------------------
# TRUE if non-fs EoE expression levels are intermediate btwn cntls and FS EoE

left_join(resultstable_FvC_signif, resultstable_EvC_all, by = "row") %>%
  mutate(interm = if_else(
    log2FoldChange.x < 0,
    log2FoldChange.x < log2FoldChange.y & log2FoldChange.y < 0,
    log2FoldChange.x > log2FoldChange.y & log2FoldChange.y > 0)) %>%
  .$interm %>%
  table(useNA = "ifany")

7320/(7320+777) # 0.904




# check # of signif genes FC > 2 -----------------------------------------------
# for consistency with previous work, note the # of genes DE with |FC| > 2

resultstable_FvC_signif %>% filter(abs(log2FoldChange) > sqrt(2)) %>% nrow
resultstable_EvC_signif %>% filter(abs(log2FoldChange) > sqrt(2)) %>% nrow
resultstable_FvE_signif %>% filter(abs(log2FoldChange) > sqrt(2)) %>% nrow




# PCA (Fig 2) ==================================================================

PCAsoln <-
  prcomp(assay(rlog_counts) %>% t,
         center = T, scale = F)

PCAsoln_percentvar <-
  PCAsoln$sdev^2 %>% `/`(., sum(.)) %>% `*`(100) %>% .[1:2] %>% format(digits = 3)

Fig2_PCA <-
  metadata_final %>%
  bind_cols(PC1 = PCAsoln$x[ , 1], PC2 = PCAsoln$x[ , 2]) %>%
  ggplot(data = ., aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 2) +
  theme_few() +
  scale_shape_manual(values = palette_shape_phenotype) +
  scale_color_manual(values = palette_color_phenotype) +
  xlab(paste0("PC1 (", PCAsoln_percentvar[1], "%)")) +
  ylab(paste0("PC2 (", PCAsoln_percentvar[2], "%)")) 

ggsave("Fig2_PCA.pdf", width = 5, height = 3.5)
dev.off()



Fig2_heatmap <-
  Heatmap(use_raster = F, show_row_dend = F,
  name = "Expression\n(scaled rlog)", col = viridis(200),
  assay(rlog_counts)[
    union(resultstable_FvC_signif %>% filter(abs(log2FoldChange) > 2) %>% .$row,
          resultstable_EvC_signif %>% filter(abs(log2FoldChange) > 2) %>% .$row) %>%
      union(resultstable_FvE_signif %>% filter(abs(log2FoldChange) > 2) %>% .$row), ] %>%
    apply(., 1, scale) %>% t %>%
    set_colnames(NULL) %>% set_rownames(NULL),
  bottom_annotation =
    columnAnnotation(Phenotype = metadata_final$Group,
                  col = list(Phenotype = palette_color_phenotype))
 )

pdf("Fig2_heatmap.pdf", width = 4, height = 7)
draw(Fig2_heatmap,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()








# plotting DE results (Fig 3) ==================================================

# RUV-adjusted resid expression ------------------------------------------------

# 27 x 17190
# residuals (remove RUV effect) for plotting
rlog_residuals <-
  rlog_counts %>%
  assay %>%
  apply(., 1,
        function(y) {
          tmp <-
            metadata_final %>%
            bind_cols(y = y)
          lm(update(test_formula, y ~ . - Phenotype),
             tmp) %>% 
            resid() }
  )

# export rlog counts for GEO / collaborators
rlog_counts %>%
  assay %>%
  t %>%
  bind_cols(ID = metadata_final$RNASeqID, .) %>%
  write_tsv("rlogvals.tsv")

rlog_residuals %>%
  as_tibble() %>%
  bind_cols(ID = metadata_final$RNASeqID, .) %>%
  write_tsv("resids_categorical.tsv")

# plotting functions -----------------------------------------------------------
plot_by_categorical <-
  function(genename, residuals = T, unfiltered_counts = F, padding = NULL) {
    
    if (unfiltered_counts) {
      Expression <- counts %>% filter(Geneid %in% genename) %>%
        select(-1) %>% t %>% .[filter_subjects, ]
      residuals <- F
      yaxis_label <- "Raw Counts"
    } else {
      if (residuals) {
        Expression <- rlog_residuals[ , genename ]
        yaxis_label <- "Residualized Expression" 
      } else {
        Expression <- assay(rlog_counts)[genename, ]
        yaxis_label <- "Expression (rlog transform)"
      }
    }
    
    yaxisrange <- 
      Expression %>% range(na.rm = T)
    if( is.null(padding) ) {
      padding <- Expression %>% sd(na.rm = T) }
    
    metadata_final %>%
      bind_cols(Expression) %>%
      ggplot(data = .,
             aes(x = Group, y = Expression, color = Group, shape = Group)) +
      geom_quasirandom(alpha = 0.8, size = 2, width = 0.1, bandwidth = 0.2) +
      theme_few() +
      stat_summary(fun.data = function(x) {
        c(quantile(x, 0.25, na.rm = T), median(x, na.rm = T),
          quantile(x, 0.75, na.rm = T)) %>%
          set_names(c("ymin", "y", "ymax")) },
        geom = "errorbar", color = "black",
        width = 0.4, size = 0.8, alpha = 0.8) +
      stat_summary(fun.data = function(x) {
        rep(median(x, na.rm = T), 3) %>%
          set_names(c("ymin", "y", "ymax")) },
        geom = "errorbar", color = "black",
        width = 0.2, size = 1.5, alpha = 1) +
      scale_x_discrete("Phenotype") + 
      scale_y_continuous(paste0(genename, " Expression"),
                         breaks = pretty_breaks(n = 6), 
                         limits = c(yaxisrange[1] - padding/3,
                                    yaxisrange[2] + padding/3)) + 
      theme(legend.position = "none") +
      scale_shape_manual(values = palette_shape_phenotype) +
      scale_color_manual(values = palette_color_phenotype)
  }


# Figure 3: DE genes fig -------------------------------------------------------
plot_grid(
  plot_by_categorical("TSPAN12"),
  plot_by_categorical("COL8A2"),
  plot_by_categorical("RABGAP1L"),
  plot_by_categorical("NPM1"))
ggsave("Fig3_DEGs.pdf", width = 6, height = 6)



