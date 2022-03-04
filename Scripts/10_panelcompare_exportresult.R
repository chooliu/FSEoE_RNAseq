# ==============================================================================
# 10_panelcompare_exportresult.R
# compare to 96-gene panel referred to throughout manuscript
# & export final results tables in tidier format (compile by gene)
# ==============================================================================




# import .xlsx from C.M.K. =====================================================
# targeted 96-gene panel known in (all/total) EoE versus controls

# load genes
panelset <-
  read_excel("./Data/Rothnberg Gene Panel.xlsx") %>%
  group_by(Gene) %>% 
  transmute(Gene,
            FC = as.numeric(FC),
            TaqManID = `TaqMan ID`,
            About = c(`...4`, `...5`, `...6`) %>%
              .[!is.na(.)] %>% paste(., collapse = " ")) %>%
  bind_rows(.,
            tribble(~Gene, ~FC, ~TaqManID, ~About, # [*]
                    "CXCL8", 27, "", "=IL8",
                    "FCGR3A", 8.6, "", "=FCGR3A/B",
                    "FCGR3B", 8.6, "", "=FCGR3A/B",
                    "TPSB2", 9.9, "", "=TPSB2/AB1",
                    "TPSAB1", 9.9, "", "=TPSB2/AB1",
                    "C7orf68", -16.8, "", "=HILPDA"
            )
            ) %>%
  mutate(FC = if_else(FC < 0, 1/abs(FC), FC))

# [*] some manual additions from checking
# which panel genes were not tested in our RNA-seq dataset 
panelset %>% filter(!Gene %in% resultstable_EvC_all$row)

# IL8 --> CXCL8
# IL5    [quantified, not tested]
# CRISP2    [quantified, not tested]
# FCGR3A/B --> FCGR3A & FCGR3B
# TPSB2/AB1 --> TPSB2 & TPSAB1
# C7orf68 --> HILPDA
# IL4    [quantified, not tested]
# CCL8   [quantified, not tested]
# EPX   [quantified, not tested]
# 18S  [not quantified]

# tested, different names
resultstable_EvC_all$row %>% .[grepl("CXCL8|FCGR3|TPSB2|TPSAB1|HILPDA", .)]

# quantified, not tested: very low counts
counts$Geneid %>% .[grepl("^IL5|^CRISP2|^IL4|^CCL8|^EPX", .)] %>% sort
counts$Geneid %>% .[grepl("RNA18SN5", .)] %>% sort

plot_by_categorical("CRISP2", unfiltered_counts = T)
plot_by_categorical("IL4", unfiltered_counts = T)
plot_by_categorical("IL5", unfiltered_counts = T)
plot_by_categorical("CCL8", unfiltered_counts = T)
plot_by_categorical("EPX", unfiltered_counts = T)

  


# summarize results by gene ====================================================
# compiling statistical tests across multiple phenotypes into one table
prep_results_for_gene_merge <- 
  function(x, prefix) {
    
    x %>%
      select(row, log2FoldChange, padj) %>%
      transmute(row, FC = 2^log2FoldChange,
                padj, if_else(padj < 0.05, "*", "")) %>%
      set_names(., c("Gene", paste0(prefix, c(": FC", ": FDR", ": Signif"))))
    
  }

summary_of_DE_results_by_gene <-
  prep_results_for_gene_merge(resultstable_EvC_all, "EvC") %>%
  inner_join(prep_results_for_gene_merge(resultstable_FvC_all, "FvC")) %>%
  inner_join(prep_results_for_gene_merge(resultstable_FvE_all, "FvE")) %>%
  inner_join(prep_results_for_gene_merge(resultstable_dysphagia_all, "Dysphagia")) %>%
  inner_join(prep_results_for_gene_merge(resultstable_EREFf_all, "EREFf")) %>%
  inner_join(prep_results_for_gene_merge(resultstable_noRUV_EvC_all, "EvC_noRUV")) %>%
  inner_join(prep_results_for_gene_merge(resultstable_noRUV_FvC_all, "FvC_noRUV")) %>%
  inner_join(prep_results_for_gene_merge(resultstable_noRUV_FvE_all, "FvE_noRUV")) %>%
  mutate(`In Rothnberg Panel` = if_else(Gene %in% panelset$Gene, "*", "")) %>%
  left_join(., panelset %>% select(1:2) %>% set_names(c("Gene", "Rothnberg Panel FC"))) %>%
  mutate(`Rothnberg Panel FC` = if_else(is.na(`Rothnberg Panel FC`), "", as.character(`Rothnberg Panel FC`)))





# questions about targeted panel genes =========================================

panelset_comparison <-
  panelset %>%
  select(1:2) %>%
  set_names(c("Gene", "FC (Rothnberg Panel)")) %>%
  arrange(Gene) %>%
  left_join(., summary_of_DE_results_by_gene, by = "Gene") %>%
  filter(!is.na(`EvC: FC`))

# were genes in the panel DEGs in our EvC and FvC comparisons?
panelset_comparison$`EvC: Signif` %>% table
panelset_comparison$`FvC: Signif` %>% table

# did the genes in the panel have consistent effect directions
# between cases and controls? (regardless of statistical significance)
panelset_comparison $`FC (Rothnberg Panel)` %>% is.na %>% sum # 15

panelset_comparison %>%
  filter((log2(`FC (Rothnberg Panel)`) * log2(`EvC: FC`)) < 0)  %>% .$Gene
panelset_comparison %>%
  filter((log2(`FC (Rothnberg Panel)`) * log2(`FvC: FC`)) < 0)  %>% .$Gene
plot_by_categorical("FCGR3A", residuals = T)
plot_by_categorical("FCGR3B", residuals = T)

# are any of the Rothnberg panel genes FS-EoE vs. non-fs EoE DEGs (FvE)?
panelset_comparison %>%
  filter(`FvE: Signif` == "*") %>% .$Gene %>% paste0(collapse = ", ")

panelset_comparison %>%
  filter(`FvE: Signif` == "*") %>%
  filter(abs(log2(`FvC: FC`)) < abs(log2(`EvC: FC`)))






# tidy results formatting ======================================================

# format DESeq2 results,
# adding columns that show stat signif in all other phenotype comparisons
export_deseq2_tables <-
  function(x) {
    
    x %>%
      transmute(
        Gene = row,
        Log2FC = log2FoldChange,
        SE = lfcSE,
        P = pvalue %>% format(digits = 3, format = "g"),
        FDR = padj %>% format(digits = 3, format = "g")
      ) %>%
      mutate(
        Log2FC = Log2FC %>% formatC(digits = 3),
        SE = SE %>% formatC(digits = 3)
      ) %>%
      mutate(
        `Signif (EvC)` = if_else(Gene %in% resultstable_EvC_signif$row, "*", ""),
        `Signif (FvC)` = if_else(Gene %in% resultstable_FvC_signif$row, "*", ""),
        `Signif (FvE)` = if_else(Gene %in% resultstable_FvE_signif$row, "*", ""),
        `Signif (Dysphagia)` = if_else(Gene %in% resultstable_dysphagia_signif$row, "*", ""),
        `Signif (EREFf)` = if_else(Gene %in% resultstable_EREFf_signif$row, "*", ""),
        `Signif (EvC, no RUV)` = if_else(Gene %in% resultstable_noRUV_EvC_signif$row, "*", ""),
        `Signif (FvC, no RUV)` = if_else(Gene %in% resultstable_noRUV_FvC_signif$row, "*", ""),
        `Signif (FvE, no RUV)` = if_else(Gene %in% resultstable_noRUV_FvE_signif$row, "*", ""),
        `Signif (Rothnberg Panel Set)` = if_else(Gene %in% panelset$Gene, "*", "")
      )
    
  }

# export final results to .xlsx
# (tab 1-10 are organized mainly by comparison/test,
# tab 11-12 compile all tests to a by-gene basis)
write_xlsx(
  list(
    `1_DE_EoE_v_Cntl` = export_deseq2_tables(resultstable_EvC_all),
    `2_DE_FSEoE_v_Cntl` = export_deseq2_tables(resultstable_FvC_all),
    `3_DE_FSEoE_v_EoE` = export_deseq2_tables(resultstable_FvE_all),
    `4_fGSEA_EoE_v_Cntl` = format_fgsea_for_export(fgsea_out$EvC, resultstable_EvC_all),
    `5_fGSEA_FSEoE_v_Cntl` = format_fgsea_for_export(fgsea_out$FvC, resultstable_FvC_all),
    `6_fGSEA_FSEoE_v_EoE` = format_fgsea_for_export(fgsea_out$FvE, resultstable_FvE_all),
    `7_DE_Dysphagia` = export_deseq2_tables(resultstable_dysphagia_all),
    `8_DE_EREFf` = export_deseq2_tables(resultstable_EREFf_all),
    `9_fGSEA_Dysphagia` = format_fgsea_for_export(fgsea_out$Dysph, resultstable_dysphagia_all),
    `10_fGSEA_EREFf` = format_fgsea_for_export(fgsea_out$EREFf, resultstable_EREFf_all),
    
    `11_AllDEResultsByGene` = summary_of_DE_results_by_gene,
    `12_ComparisonToRothnberg` = panelset_comparison
  ),
  path = paste0("./Output/FSEoE_RNAseq-Results.xlsx"),
  format_headers = F)




