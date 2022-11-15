Analysis code accompanying "Fibroproliferative gene signatures define patients with fibrostenotic eosinophilic esophagitis", now published in _Allergy_.

We profile esophageal mucosal transcriptomes using RNA-sequencing, and detect gene expression differences between pediatric/adolescent patients diagnosed with eosinophilic esophagitis (EoE), a fibrostenotic EoE phenotype (FS-EoE), and controls.

```
File Structure
├─ README.md                                # this readme file
├─ Scripts/                                 # R code; run in sequential order
│   └─ 00_setup.R                           # - library load (also see sessionInfo)
│   └─ 01_load_data.R                       # - counts (from M.D.G.) & metadata
│   └─ 02_prefiltering.R                    # - QC on counts, partiicpants
│   └─ 03_table1_and_figures.R              # - check metadata by phenotype (categorical group)
│   └─ 04_deseq_categorical.R               # - diff expression by group
│   └─ 05_figures_DE_PCA_heatmap.R          # - visualize diffs (PCA, venn, scatter)
│   └─ 06_deseq_dysphagia.R                 # - diff exp, continuous FS-EoE phenotype definitions
│   └─ 07_deseq_EREFSf.R                    # - " "
│   └─ 08_volcanoplots_foldchangecompare.R  # - check logFCs
│   └─ 09_fGSEA.R                           # - pathway analyses
│   └─ 10_panelcompare_exportresult.R       # - summarize all, export
└─ sessionInfo.txt                          # library version control
```

Deidentified participant level metadata and raw/processed RNA-seq data are publicly available at GEO, accession [GSE197702](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197702).

**Citation:** Menard-Katcher C, Liu C, Galbraith MD, Benson T, Burger C, Dobias D, Larsen L, O'Brien C, Spencer LA, Furuta GT, Masterson JC. Fibrostenotic Eosinophilic Esophagitis phenotype is defined by a proliferative gene signature. Allergy. 2022 Oct 22.

**Manuscript:** [Wiley Allergy](https://onlinelibrary.wiley.com/doi/10.1111/all.15557) \
**PMID:** [36273270](https://pubmed.ncbi.nlm.nih.gov/36273270/) \
**doi:** [doi:10.1111/all.15557](https://doi.org/10.1111/all.15557)
