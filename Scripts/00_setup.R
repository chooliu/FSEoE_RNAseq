# ==============================================================================
# 00_setup.R
# load R libraries
# ==============================================================================



# libraries 
library(magrittr)
library(tidyverse)

# modeling/stats
library(DESeq2)
library(readxl)
library(RUVSeq)
library(vegan)
library(eulerr)

# GO/pathway analysis
library(ontologyIndex)
library(fgsea)
library(data.table)

# plotting
library(ggthemes)
library(ggrepel)
library(cowplot)
library(viridis)
library(ggbeeswarm)
library(scales)
library(writexl)
library(ComplexHeatmap)


# version control
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



