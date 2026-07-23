# ==========================================================
#
# RNA-Seq Analysis using Linux and R
#
# Module:
# Gene Expression Statistics
#
# Description:
# Statistical inference for paired gene expression analysis
# using biological replicates.
#
# Dataset:
# Toxoplasma gondii RNA-Seq
#
# ==========================================================

# ==========================================================
# 1. Load required packages
# ==========================================================

library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ==========================================================
# 2. Load the SummarizedExperiment object
# ==========================================================

load("/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/se_toxo_gene.RData")

# ==========================================================
# 3. Extract TPM matrix
# ==========================================================

tpm <- assay(se_gene, "abundance")

# ==========================================================
# 4. Select infected samples
# ==========================================================

infected_samples <- c(
  "quant_INF1",
  "quant_INF2",
  "quant_INF3")

tpm_inf <- tpm[, infected_samples]

# ==========================================================
# 5. Import biomarker list
# ==========================================================

markers <- read.csv(
  "/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/genes_markers.csv",
  stringsAsFactors = FALSE)

# ==========================================================
# 6. Select biomarker genes
# ==========================================================

marker_matrix <- tpm_inf[markers$Gene_ID, ]

rownames(marker_matrix) <- markers$Gene

# ==========================================================
# 7. Create paired dataset
# ==========================================================

bag1_sag1 <- data.frame(
  Sample = infected_samples,
  BAG1 = as.numeric(marker_matrix["BAG1", ]),
  SAG1 = as.numeric(marker_matrix["SAG1", ]))

bag1_sag1

# ==========================================================
# Experimental Design
# ==========================================================

# Before performing any statistical analysis,
# we must understand the experimental design.

# INF1, INF2 and INF3 represent biological replicates.

# Biological replicates are independent biological
# samples subjected to the same experimental condition.

# Biological replication allows estimation of
# biological variability.

############################################################

# IMPORTANT
#
# BAG1 and SAG1 are NOT biological replicates.
#
# They are two different genes measured in the
# same biological samples.
#
# Therefore, observations are naturally paired.
#
############################################################

# ==========================================================
# Statistical Hypotheses
# ==========================================================

# H0:
#
# The average paired difference between BAG1
# and SAG1 is equal to zero.
#
# H1:
#
# The average paired difference is different
# from zero.

# ==========================================================
# Calculate paired differences
# ==========================================================

bag1_sag1$Difference <- bag1_sag1$BAG1 - bag1_sag1$SAG1

bag1_sag1

############################################################

# A paired statistical test does not compare
# two independent groups.

# Instead, it evaluates whether the paired
# differences are consistently different from zero.

############################################################

# ==========================================================
# Descriptive statistics
# ==========================================================

summary(bag1_sag1)
mean(bag1_sag1$Difference)
sd(bag1_sag1$Difference)
median(bag1_sag1$Difference)

# ==========================================================
# Plot paired differences
# ==========================================================

ggplot(bag1_sag1, aes(x = Sample, y = Difference)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red") +
  geom_point(size = 4) +
  geom_line(group = 1) +
  theme_classic(base_size = 14) +
  labs(
    title = "Paired differences (BAG1 - SAG1)",
    y = "Difference (TPM)",
    x = "Biological replicate")

############################################################

# Positive values indicate higher BAG1 expression.

# Negative values indicate higher SAG1 expression.

# If all biological replicates remain on the
# same side of zero, the expression difference
# is highly consistent.

############################################################

# ==========================================================
# Test normality of paired differences
# ==========================================================

shapiro.test(bag1_sag1$Difference)

# ==========================================================
# QQ Plot
# ==========================================================

qqnorm(bag1_sag1$Difference, pch = 19)

qqline(bag1_sag1$Difference, col = "red", lwd = 2)

############################################################

# Shapiro-Wilk evaluates whether paired
# differences follow a normal distribution.

# It does NOT test BAG1 separately.

# It does NOT test SAG1 separately.

############################################################

############################################################

# IMPORTANT

# With only three biological replicates,
# the Shapiro-Wilk test has very low power.

# Therefore, graphical inspection is equally
# important when evaluating normality.

############################################################

# ==========================================================
# Interpretation of the normality assessment
# ==========================================================

############################################################

# The Shapiro-Wilk test returned p > 0.05.

# Therefore, we do not reject the assumption
# of normality for the paired differences.

# However, because we only have three biological
# replicates, this result should be interpreted
# with caution.

# The Shapiro-Wilk test has very low statistical
# power for such small sample sizes.

# Consequently, graphical inspection (QQ plot)
# and biological interpretation remain essential.

############################################################

# ==========================================================
# Paired Student's t-test
# ==========================================================

paired_ttest <- t.test(bag1_sag1$BAG1, 
                       bag1_sag1$SAG1,
                       paired = TRUE,
                       conf.level = 0.95)

paired_ttest

# ==========================================================
# Interpretation
# ==========================================================

if(paired_ttest$p.value < 0.05){
  message("The paired t-test detected a statistically significant difference.")
}else{
  message("The paired t-test did not detect a statistically significant difference.")
}

# ==========================================================
# Wilcoxon signed-rank test
# ==========================================================

paired_wilcox <-  wilcox.test(bag1_sag1$BAG1,
                              bag1_sag1$SAG1,
                              paired = TRUE,
                              exact = FALSE)

paired_wilcox

############################################################

# The Wilcoxon signed-rank test does not assume
# normality.

# It evaluates whether the median paired
# difference differs from zero.

############################################################

# ==========================================================
# Compare statistical tests
# ==========================================================

results <- data.frame(
    Test = c(
      "Paired t-test",
      "Wilcoxon signed-rank"),
    
    Statistic = c(
      paired_ttest$statistic,
      paired_wilcox$statistic),
    
    P_value = c(
      paired_ttest$p.value,
      paired_wilcox$p.value))

results

write.csv(results, "Biomarker_Expression_Analysis/results/statistical_tests.csv",
          row.names = FALSE)

# ==========================================================
# Effect size (Cohen's d)
# ==========================================================
mean_difference <- mean(bag1_sag1$Difference)
sd_difference <- sd(bag1_sag1$Difference)
cohens_d <- mean_difference / sd_difference
cohens_d

############################################################

# Cohen's d for paired samples is calculated as

# mean(differences) / sd(differences)

############################################################

if(abs(cohens_d) < 0.2){
  interpretation <- "Negligible"
}else if(abs(cohens_d) < 0.5){
  interpretation <- "Small"
}else if(abs(cohens_d) < 0.8){
  interpretation <- "Medium"
}else{
  interpretation <- "Large"
}

interpretation


############################################################

# IMPORTANT
#
# Statistical significance and biological relevance
# are not the same thing.
#
# The p-value indicates how compatible the data are
# with the null hypothesis.
#
# Cohen's d quantifies the magnitude of the observed
# biological effect.
#
# Therefore, both metrics should always be interpreted
# together.

############################################################