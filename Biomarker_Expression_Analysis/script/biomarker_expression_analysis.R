# ==========================================================
#
# RNA-Seq Analysis using Linux and R
#
# Module:
# Biomarker Expression Analysis
#
# Description:
# Explore stage-specific biomarkers of Toxoplasma gondii
# using TPM values, heatmaps, Z-score normalization,
# paired statistical tests and biological interpretation.
#
# Dataset:
# RNA-Seq from infected and non-infected samples
#
# Author:
# Jean Resende
#
# ==========================================================

# ==========================================================
# 1. Load required packages
# ==========================================================

library(SummarizedExperiment)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

# ==========================================================
# 2. Load the SummarizedExperiment object
# ==========================================================

load("/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/se_toxo_gene.RData")

# ==========================================================
# 3. Explore the dataset
# ==========================================================
# What assays are available?
assayNames(se_gene)

# How many genes and samples?
dim(se_gene)

# Sample information
colData(se_gene)

# Gene annotation
rowData(se_gene)

# Available metadata
metadata(se_gene)

# ==========================================================
# 4. Extract TPM values
# ==========================================================
tpm <- assay(se_gene,"abundance")
head(tpm)

# ==========================================================
# 5. Select infected samples
# ==========================================================
# We are interested in parasite gene expression.
# Therefore, only infected samples are used.

infected_samples <- c(
  "quant_INF1",
  "quant_INF2",
  "quant_INF3")

tpm_inf <- tpm[, infected_samples]
dim(tpm_inf)

# ==========================================================
# 6. Import biomarker list
# ==========================================================
markers <- read.csv(
  "/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/genes_markers.csv",
  stringsAsFactors = FALSE
)

markers

# ==========================================================
# 7. Select biomarker genes
# ==========================================================
marker_matrix <- tpm_inf[markers$Gene_ID, ]

rownames(marker_matrix) <- markers$Gene

marker_matrix

# ==========================================================
# Question 1
#
# Are Bradyzoite markers more expressed than
# Tachyzoite markers?
#
# ==========================================================
bag1_sag1 <- data.frame(
    Sample = infected_samples,
    BAG1 = marker_matrix["BAG1",],
    SAG1 = marker_matrix["SAG1",]
  )

bag1_sag1

############################################################

# Before calculating any statistics,
# always inspect the raw data.
#
# Ask yourself:
#
# • Is one gene consistently more expressed?
# • Is there high variability?
# • Are there possible outliers?

############################################################

# ==========================================================
# Descriptive statistics
# ==========================================================

summary(bag1_sag1)

colMeans(bag1_sag1[ ,c("BAG1","SAG1")])

apply(bag1_sag1[, c("BAG1", "SAG1")], 2, sd)

# ==========================================================
# Fold Change
# ==========================================================

mean_bag1 <- mean(bag1_sag1$BAG1)
mean_sag1 <- mean(bag1_sag1$SAG1)

fold_change <- mean_bag1 / mean_sag1
fold_change

log2_fc <- log2(fold_change)
log2_fc

############################################################

# Fold Change > 1
#
# BAG1 shows higher average expression.
#
# Fold Change < 1
#
# SAG1 shows higher average expression.

############################################################

# ==========================================================
# Boxplot
# ==========================================================

library(tidyr)

bag1_sag1_long <- pivot_longer(
  bag1_sag1,
  cols = c(BAG1, SAG1),
  names_to = "Gene",
  values_to = "TPM")

ggplot(bag1_sag1_long, aes(x = Gene, y = TPM, fill = Gene)) +
  geom_boxplot(width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 3) +
  theme_classic(base_size = 14) +
  labs(
    title = "Expression of BAG1 and SAG1",
    y = "TPM")

############################################################

# IMPORTANT
#
# With only three biological replicates,
# the boxplot provides limited information.
#
# The individual observations are often more informative
# than the box itself.
#
############################################################

# ==========================================================
# Paired Plot
# ==========================================================
paired_data <- pivot_longer(
  bag1_sag1,
  cols = c(BAG1, SAG1),
  names_to = "Gene",
  values_to = "TPM")

ggplot(paired_data, aes(x = Gene, y = TPM, group = Sample, color = Sample)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_classic(base_size = 14) +
  labs(
    title = "Paired comparison between BAG1 and SAG1",
    y = "TPM")

############################################################

# Why is this plot important?
#
# Each line represents one biological replicate.
#
# If all lines point in the same direction,
# the difference is consistent among replicates.
#
# Consistency across biological replicates
# is often more informative than simply comparing means.
#
############################################################

