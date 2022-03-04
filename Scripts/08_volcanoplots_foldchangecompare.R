# ==============================================================================
# 08_volcanoplots_foldchangecompare.R
# generate volcano plot for all phenotypes; check FvC and EvC effect size
# ==============================================================================




# generic volcano plot function ================================================

palette_color_signif_volcano <-
  c(pos = "#5ab4ac", neg = "#8c510a", NS = "grey")

palette_shape_signif_volcano <-
  c(pos = 19, neg = 19, NS= 1)

make_volcano_plot <-
  function(volcano_tmp, nbreaksx = 6, nbreaksy = 6) {
    
     # which genes to label
  volcano_tmp_symbols_to_include <-
    c(
      volcano_tmp %>% filter(log2FoldChange > 0) %>% .$row %>% .[1:5],
      volcano_tmp %>% filter(log2FoldChange < 0) %>% .$row %>% .[1:5],
      volcano_tmp %>% filter(padj < 0.05) %>% arrange(log2FoldChange) %>%
        .$row %>% .[1:5],
      volcano_tmp %>% filter(padj < 0.05) %>% arrange(-log2FoldChange) %>%
        .$row %>% .[1:5]
    )
  
  volcano_tmp$row[!volcano_tmp$row %in% volcano_tmp_symbols_to_include] <- ""
  
  # add coloring, transformation on log10
  volcano_tmp <-
    volcano_tmp %>%
    mutate(
      Signif = case_when(
        padj < 0.05 & log2FoldChange > 0 ~ "pos",
        padj < 0.05 & log2FoldChange < 0 ~ "neg",
        T ~ "NS"
      ),
      neglog = -log10(pvalue)
    )
  
  # tidy breaks
  xaxisrange <- volcano_tmp$log2FoldChange %>% range(na.rm = T)
    xpadding <- volcano_tmp$log2FoldChange %>% sd(na.rm = T)
    yaxisrange <- volcano_tmp$neglog %>% range(na.rm = T)
    ypadding <- volcano_tmp$neglog %>% sd(na.rm = T)
  
    # plot
    ggplot(
      data = volcano_tmp,
      aes(
        x = log2FoldChange, y = neglog,
        color = Signif, shape = Signif
      )
    ) +
    geom_point(alpha = 0.6, size = 0.8) +
    geom_text_repel(aes(label = row),
      color = "black",
      box.padding = unit(0.01, "cm"), max.overlaps = 20,
      size = 2.5
    ) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_color_manual(values = palette_color_signif_volcano) +
    scale_shape_manual(values = palette_shape_signif_volcano) +
    theme_few() +
    theme(legend.position = "none") +
    scale_x_continuous(expression(log[2] ~ "(Fold Change)"),
                       breaks = pretty_breaks(n = nbreaksx), 
                       limits = c(xaxisrange[1] - xpadding,
                                  xaxisrange[2] + xpadding)) +
      scale_y_continuous(expression(-log[10] ~ "(p-value)"),
         breaks = pretty_breaks(n = nbreaksy), 
         limits = c(yaxisrange[1] - ypadding/3,
                    yaxisrange[2] + ypadding/3))

}




# make and save volcano plots ==================================================

make_volcano_plot(resultstable_EvC_all)
ggsave("Fig2_Volcano_EvC.pdf", width = 5, height = 4)

make_volcano_plot(resultstable_FvC_all, 8)
ggsave("Fig2_Volcano_FvC.pdf", width = 5, height = 4)

make_volcano_plot(resultstable_FvE_all, 7)
ggsave("Fig2_Volcano_FvE.pdf", width = 5, height = 4)

make_volcano_plot(resultstable_dysphagia_all)
ggsave("Fig4_Volcano_Dysphagia.pdf", width = 5, height = 4)

make_volcano_plot(resultstable_EREFf_all)
ggsave("Fig4_Volcano_EREFSf.pdf", width = 5, height = 4)




# check relative effect size of FvC vs EvC =====================================
# FvC magnitudes tend to be more extreme, same direction as EvC

left_join(resultstable_FvC_signif,
          resultstable_EvC_all,
          by = "row", suffix = c(".FvC", ".EvC")) %>%
  ggplot(data = ., aes(log2FoldChange.EvC, log2FoldChange.FvC)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_hex(bins = 100, alpha = 0.75) +
  geom_abline(lty = 2, intercept = 0, slope = 1) +
  scale_fill_viridis_c(name = "# Genes",
                       breaks = c(1, 10, 50, 100, 250, 500),
                       trans = pseudo_log_trans(sigma = 0.001)) +
  theme_few() +
  scale_x_continuous(expression(log[2] ~ "(FC): FS EoE vs. Cntl"),
                     breaks = seq(-8, 12, 4), limits = c(-8, 12)) +
  scale_y_continuous(expression(log[2] ~ "(FC): non-fs EoE vs. Cntl"),
                     breaks = seq(-8, 12, 4), limits = c(-8, 12))

ggsave("Fig3_logFC_Comparison.pdf", width = 5, height = 4)



