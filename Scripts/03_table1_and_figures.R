# ==============================================================================
# 03_table1.R
# demographic / clinical features by group
# ==============================================================================




# Table 1 fxns ================================================================
# existing fxns i frequently use to summarize demographics

summaryCountPercent <-
  function(x, values, count_NAs = T,
           digits = 1,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = F,
           NA_count_in_brackets_alt = T,
           fuzzy = F, inverse = F) {
    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100
    
    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }
    
    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }
    
    # exact match (default)
    n <- sum(x %in% values, na.rm = T)
    tot <- ifelse(count_NAs, length(x), length(x) - n_NA)
    
    # fuzzy matching
    if (fuzzy) {
      search_term <-
        paste(values, collapse = "|")
      n <- grepl(search_term, x, ignore.case = T) %>% sum()
    }
    
    # inverse
    if (inverse) {
      n <- tot - n
    }
    
    output <- paste0(
      n, " (", formatC(n / tot * 100, format = "f", digits = digits),
      "%)"
    ) 
    
    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }
    
    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }
    
    return(output)
  }

summaryMedianIQR <-
  function(x, digits = 1, na.rm = F,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = F,
           NA_count_in_brackets_alt = T,
           symbol = NULL) {
    if (!(typeof(x) %in% c("double", "integer"))) {
      warning("input vector not numeric")
      x <- as.numeric(x)
    }
    
    tot <- length(x)
    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100
    
    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }
    
    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }
    
    output <- quantile(
      as.numeric(x),
      c(0.25, 0.5, 0.75),
      na.rm = na.rm
    )
    
    output <- formatC(output, format = "f", digits = digits)
    
    output <- ifelse(is.null(symbol),
                     paste0(output[2], " (", output[1], ", ", output[3], ")"),
                     paste0(output[2], symbol, " (", output[1],
                            symbol, ", ", output[3], symbol, ")")
    )
    
    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }
    
    if (NA_count_in_brackets_alt & !NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [n = ", tot - n_NA, "]")
    }
    
    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }
    
    return(output)
  }




# calc summary stats ===========================================================

Table1 <-
  metadata_final %>%
  bind_rows(., metadata_final_nocntl %>% mutate(Group = "Total EoE")) %>%
  bind_rows(., metadata_final %>% mutate(Group = "All")) %>%
  group_by(Group) %>%
  summarize(N = n() %>% as.character,
            `Gender (Male)` = summaryCountPercent(Gender, "Male"),
            Age = summaryMedianIQR(Age),
            Dilation = summaryCountPercent(Dilation, "Yes"),
            Distensibility = summaryMedianIQR(distensibility, na.rm = T),
            `Dysphagia Severity` = summaryMedianIQR(Dysphagia_Severity),
            EREFi = summaryMedianIQR(EREFi),
            EREFf = summaryMedianIQR(EREFf),
            `Eos/hpf (Min)` = summaryMedianIQR(min_eos_hpf),
            `Eos/hpf (Average)` = summaryMedianIQR(average_eos_hpf),
            `Eos/hpf (Peak)` = summaryMedianIQR(peak_eos_hpf),
            `Fibrotic Score` = summaryMedianIQR(Fibrotic_Score),
            `Food Impaction` = summaryCountPercent(Food_impaction, "Yes"),
            `Medication-PPI` = summaryCountPercent(Medication_PPI, "Yes"),
            `Medication-STCs` = summaryCountPercent(Medication_STCs, "Yes")) %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  pivot_wider(id_cols = name, values_from = value, names_from = Group) %>%
  rename(name = "Feature")




# 3-group and 2-group testing ==================================================

categorical_test_3group <-
  function(x) {
    fisher.test(metadata_final$Group,
                metadata_final[ , x, drop = T])$p.value
  }

continuous_test_3group <-
  function(x) {
    kruskal.test(metadata_final[ , x, drop = T],
                 metadata_final$Group)$p.value
  }

categorical_test_2group <-
  function(x) {
    fisher.test(metadata_final$Group[filter_subjects_nocntl],
                metadata_final[filter_subjects_nocntl, x, drop = T])$p.value
  }
continuous_test_2group <-
  function(x) {
    kruskal.test(metadata_final[filter_subjects_nocntl, x, drop = T],
                 metadata_final$Group[filter_subjects_nocntl])$p.value
  }

fmt_scientific <-
  function(x) {
    if_else(x < 0.001, formatC(x, format = "e", digits = 2),
            formatC(x, digits = 2)) %>%
      paste0(if_else(x < 0.05, paste0(" *"), "")) 
  }




# run statistical testing ======================================================

pval_three_group <-
  c(NA_real_,
    categorical_test_3group("Gender"),
    continuous_test_3group("Age"),
    categorical_test_3group("Dilation"),
    continuous_test_3group("distensibility"),
    continuous_test_3group("Dysphagia_Severity"),
    continuous_test_3group("EREFi"),
    continuous_test_3group("EREFf"),
    continuous_test_3group("min_eos_hpf"),
    continuous_test_3group("average_eos_hpf"),
    continuous_test_3group("peak_eos_hpf"),
    continuous_test_3group("Fibrotic_Score"),
    categorical_test_3group("Food_impaction"),
    categorical_test_3group("Medication_PPI"),
    categorical_test_3group("Medication_STCs")
  )

pval_two_group <-
  c(NA_real_,
    categorical_test_2group("Gender"),
    continuous_test_2group("Age"),
    categorical_test_2group("Dilation"),
    continuous_test_2group("distensibility"),
    continuous_test_2group("Dysphagia_Severity"),
    continuous_test_2group("EREFi"),
    continuous_test_2group("EREFf"),
    continuous_test_2group("min_eos_hpf"),
    continuous_test_2group("average_eos_hpf"),
    continuous_test_2group("peak_eos_hpf"),
    continuous_test_2group("Fibrotic_Score"),
    categorical_test_2group("Food_impaction"),
    categorical_test_2group("Medication_PPI"),
    categorical_test_2group("Medication_STCs")
  )




# export Table 1  ==============================================================

Table1 <-
  Table1 %>%
  bind_cols(`P (3-group)` = pval_three_group %>% fmt_scientific,
            `P (EoE vs. FS-EoE)` = pval_two_group %>% fmt_scientific) %>%
  mutate_all(function(x) { gsub("NANA", "", x) } ) # remove NAs

write_xlsx(Table1,
           path = "./Table1.xlsx")




# palettes =====================================================================

palette_shape_phenotype <- c(
  Control = 1, `non-fs EoE` = 9, `FS EoE` = 15
)

palette_color_phenotype <- c(
  Control = "#e41a1c", `non-fs EoE` = "#33a02c", `FS EoE` = "#1f78b4"
)



# plots of clinical indices by group for Figure 1A =============================

plot_clinical_by_group <-
  function(yvariable, ylabel = yvariable, padding = NULL) {
    
    yaxisrange <- 
      metadata_final[ , yvariable, drop = T] %>% range(na.rm = T)
    if( is.null(padding) ) {
        padding <- metadata_final[ , yvariable, drop = T] %>% sd(na.rm = T) }
    
  ggplot(metadata_final,
         aes_string("Group", yvariable, color = "Group", shape = "Group")) +
    geom_quasirandom(size = 2.5, alpha = 0.6, width = 0.25, bandwidth = 0.2) +
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
    theme_few() +
      scale_y_continuous(ylabel, breaks = pretty_breaks(n = 5), 
                         limits = c(yaxisrange[1] - padding/3,
                                    yaxisrange[2] + padding/3)) +
      xlab("Phenotype") +
    scale_color_manual(values = palette_color_phenotype) +
    scale_shape_manual(values = palette_shape_phenotype) +
    theme(legend.position = "none")
  }

plot_grid(
  plot_clinical_by_group("distensibility", "Distensibility (mm)"),
  plot_clinical_by_group("EREFf", "EREFSf"),
  plot_clinical_by_group("Fibrotic_Score", "Fibrotic Score"),
  plot_clinical_by_group("Dysphagia_Severity", "Dysphagia Severity"),
  nrow = 1
  )

ggsave("Fig1A.pdf", width = 12, height = 3.5)





# spearman cor of clinical indices for Figure 1B ===============================

plot_correlations <-
  function(xvariable, yvariable,
           xlabel = gsub("_", " ", xvariable),
           ylabel = gsub("_", " ", yvariable),
           xpadding = NULL, ypadding = NULL) {
    
    cortest <-
      cor.test(metadata_final_nocntl[ , xvariable, drop = T], 
               metadata_final_nocntl[ , yvariable, drop = T],
               method = "spearman")
    
    corlabel <-
      paste0("Spearman r = ", formatC(cortest$estimate, digits = 2),
             "; p = ", cortest$p.value %>% fmt_scientific())
    
    xaxisrange <- 
      metadata_final[ , xvariable, drop = T] %>% range(na.rm = T)
    if( is.null(xpadding) ) {
      xpadding <- metadata_final[ , xvariable, drop = T] %>% sd(na.rm = T) }    
    yaxisrange <- 
      metadata_final[ , yvariable, drop = T] %>% range(na.rm = T)
    if( is.null(ypadding) ) {
      ypadding <- metadata_final[ , yvariable, drop = T] %>% sd(na.rm = T) }
    
    ggplot(metadata_final_nocntl,
           aes_string(xvariable, yvariable,
                      color = "Group", shape = "Group")) +
      geom_quasirandom(alpha = 0.6, size = 2) +
      theme_few() +
      scale_shape_manual(values = palette_shape_phenotype) +
      scale_color_manual(values = palette_color_phenotype) +
      scale_x_continuous(xlabel, breaks = pretty_breaks(n = 5), 
                         limits = c(xaxisrange[1] - xpadding/3,
                                    xaxisrange[2] + xpadding/3)) +
      scale_y_continuous(ylabel, breaks = pretty_breaks(n = 5), 
                         limits = c(yaxisrange[1] - ypadding/3,
                                    yaxisrange[2] + ypadding/3)) +
      ggtitle("", corlabel) +
      theme(legend.position = "none")
}


plot_grid(
  plot_correlations("Dysphagia_Severity", "Fibrotic_Score"),
  plot_correlations("Dysphagia_Severity", "distensibility",
                    ylabel = "Distensibility (mm)"),
  plot_correlations("Dysphagia_Severity", "EREFf", ylabel = "EREFSf"),
  plot_correlations("Dysphagia_Severity", "EREFi", ylabel = "EREFSi"),
  nrow = 1)
ggsave("Fig1B.pdf", width = 12, height = 3.5)


