# ==============================================================================
# 09_fGSEA.R
# process GO data --> fGSEA testing
# ==============================================================================




# GO ontology (retrieved June 15th 2020) =======================================

# processing via ontologyIndex package
go_obo <-
  get_ontology("./Data/Annotations/go-basic.obo")
go_ont <-
  enframe(go_obo$name) %>%
  set_names(c("pathway", "Description"))

go_annotations <-
  fread("Data/Annotations/goa_human.gaf", header = F) %>%
  .[ V9 == "P", ] %>% # biological process
  .[ , c("V3", "V5")] %>%
  set_names(c("Symbol", "GOTerm")) %>%
  unique %>%
  filter(Symbol %in% rownames(counts_final))

# formatting for fGSEA, include categories with at least 20 genes tested
go_annotations_as_list <-
  go_annotations %>%
  split(.$GOTerm) %>%
  map(., function(x) {x$Symbol %>% unique}) %>%
  .[map_dbl(., length) %>% `>=`(., 20)]




# run fGSEA ====================================================================

fgsea_out <-
  list(
    FvE = resultstable_FvE_all$log2FoldChange %>%
      set_names(resultstable_FvE_all$row),
    FvC = resultstable_FvC_all$log2FoldChange %>%
      set_names(resultstable_FvE_all$row),
    EvC = resultstable_EvC_all$log2FoldChange %>%
      set_names(resultstable_FvE_all$row),
    Dysph = resultstable_dysphagia_all$log2FoldChange %>%
      set_names(resultstable_dysphagia_all$row),
    EREFf = resultstable_EREFf_all$log2FoldChange %>%
      set_names(resultstable_EREFf_all$row)
    ) %>%
  map(.,
      function(inputlist) {
        set.seed(1234)
        fgsea(
          pathways = go_annotations_as_list,
          stats = inputlist,
          eps = 1e-30)
        })




# functions for export =========================================================

extract_top_ten_symbol <- function(listgenes, resultstable) {
  test <-
    unlist(listgenes)
  resultstable %>%
    filter(row %in% test) %>%
    mutate(row = if_else(padj < 0.05, paste0(row, "*"), row)) %>%
    arrange(-abs(log2FoldChange)) %>%
    .$row %>%
    .[1:min(length(.), 8)] %>%
    paste(collapse = ", ")
}

format_fgsea_for_export <- function(x, resultstable) {

  x %>%
    left_join(go_ont) %>%
    arrange(padj, -abs(NES)) %>%
    rowwise() %>%
    transmute(`GO Term` = pathway,
              `GO Description` = Description,
              `# Genes` = size,
              NES = formatC(NES, digits = 2, format = "f"),
              `P-Value` = formatC(pval, digits = 2, format = "g"),
              FDR = formatC(padj, digits = 2, format = "g"),
              `Signif (FDR < 0.05)` = if_else(padj < 0.05, "*", ""),
              `Leading Edge` = extract_top_ten_symbol(leadingEdge, resultstable))
}



