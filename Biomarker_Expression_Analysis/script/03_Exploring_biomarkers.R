# ==========================================================
# Module 3
# Multi-marker Expression Analysis
#
# Objective:
# Explore the expression of multiple biomarkers and
# construct a molecular signature representing the
# biological stage of Toxoplasma gondii.
#
# Dataset:
# T. gondii RNA-seq
#
# Author:
# Jean Resende
# ==========================================================



############################################################
##
## 1. Load packages
##
############################################################

library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(rstatix)



############################################################
##
## 2. Load the RNA-seq dataset
##
############################################################

load("/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/se_toxo_gene.RData")



############################################################
##
## 3. Extract TPM matrix
##
############################################################

tpm <- assay(
  se_gene,
  "abundance")



############################################################
##
## 4. Select infected samples
##
############################################################

infected_samples <- colnames(se_gene)[colData(se_gene)$sample_type == "infected"]

tpm_inf <- tpm[, infected_samples]


############################################################
##
## 5. Load biomarker annotation
##
############################################################

markers <- read.csv(
  "/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/genes_markers.csv",
  stringsAsFactors = FALSE)



############################################################
##
## 6. Build biomarker expression matrix
##
############################################################

marker_matrix <- tpm_inf[markers$Gene_ID,]

rownames(marker_matrix) <-  markers$Gene



############################################################
##
## 7. Inspect the biomarker matrix
##
############################################################

marker_matrix

dim(marker_matrix)

head(marker_matrix)

############################################################

# Biological question

# Can a single biomarker accurately describe
# parasite differentiation?

# Or should we analyze all biomarkers together?

############################################################

############################################################
##
## 8. Descriptive statistics
##
############################################################

marker_stats <- data.frame(
    Gene = rownames(marker_matrix),
    Mean = rowMeans(marker_matrix),
    Median =
      apply(marker_matrix, 1, median),
    SD = apply(marker_matrix, 1, sd),
    Minimum = apply(marker_matrix, 1, min),
    Maximum = apply(marker_matrix, 1, max))



marker_stats$CV <- (marker_stats$SD / marker_stats$Mean)*100

marker_stats

############################################################
##
## 9. Rank biomarkers
##
############################################################

marker_stats <- marker_stats |>
  arrange(desc(Mean))


marker_stats

write.csv(marker_stats, "Biomarker_Expression_Analysis/results/marker_statistics.csv", 
          row.names=FALSE)

############################################################
##
## 10. Convert to long format
##
############################################################

marker_long <-
  marker_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("Gene") |>
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "TPM"
  )

############################################################
##
## 11. Biomarker expression
##
############################################################

ggplot(marker_long, aes(Gene, TPM, fill=Sample))+
  geom_col(position="dodge")+
  theme_classic(base_size=14)+
  theme(axis.text.x= element_text(angle=45, hjust=1))+
  labs(
    title="Expression of stage-specific biomarkers",
    y="TPM")



# ==========================================================
# Module 3
# Part 2
#
# Building a Molecular Signature
#
# Objective:
# Construct stage-specific molecular signatures by
# summarizing the coordinated expression of multiple
# biomarkers.
# ==========================================================

############################################################
##
## Biological Question
##
## Until now we analyzed each biomarker individually.
##
## However, parasite differentiation is regulated by
## coordinated changes in multiple genes.
##
## Can we summarize the expression of all biomarkers
## into a single molecular signature?
##
############################################################

############################################################
##
## 1. Separate biomarkers by parasite stage
##
############################################################

brady_genes <- markers |>
  filter(Stage == "Bradyzoite") |>
  pull(Gene)

tachy_genes <-
  markers |>
  filter(Stage == "Tachyzoite") |>
  pull(Gene)

length(brady_genes)
length(tachy_genes)
brady_genes
tachy_genes


############################################################
##
## 2. Calculate signature scores
##
############################################################

signature <- data.frame(
  Sample = colnames(marker_matrix),
  Bradyzoite = colMeans(marker_matrix[brady_genes, ]),
  Tachyzoite = colMeans(marker_matrix[tachy_genes, ]))

signature

############################################################
##
## 3. Calculate paired differences
##
############################################################

signature$Difference <- signature$Bradyzoite - signature$Tachyzoite

signature

############################################################
##
## 4. Descriptive statistics
##
############################################################

summary(signature)

apply(signature[,2:4], 2, mean)

apply(signature[,2:4], 2, sd)

############################################################
##
## 5. Convert to long format
##
############################################################

signature_long <-
  signature |>
  pivot_longer(
    cols = c(Bradyzoite, Tachyzoite),
    names_to = "Stage",
    values_to = "Expression")

############################################################
##
## 6. Molecular signature
##
############################################################

ggplot(signature_long, aes(Stage, Expression, group = Sample, color = Sample))+
  geom_line(size=1)+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(
    title = "Stage-specific molecular signatures",
    y = "Mean TPM")

############################################################
##
## 7. Statistical inference
##
############################################################

signature

############################################################
##
## 8. Descriptive statistics
##
############################################################

summary(signature)

apply(signature[,2:4], 2, mean)

apply(signature[,2:4], 2, median)

apply(signature[,2:4], 2, sd)

############################################################
##
## 9. Paired differences
##
############################################################

ggplot(signature, aes(Sample, Difference, group = 1))+
  geom_line(color="steelblue")+
  geom_point(size=3)+
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "red"
  )+
  theme_classic(base_size = 14)+
  labs(
    title = "Difference between molecular signatures",
    y = "Bradyzoite - Tachyzoite"
  )

############################################################
##
## 10. Shapiro-Wilk test
##
############################################################

shapiro_result <- shapiro.test(signature$Difference)

shapiro_result

############################################################
##
## 11. QQ Plot
##
############################################################

qqnorm(signature$Difference)

qqline(
  signature$Difference,
  col="red",
  lwd=2
)

############################################################
##
## 12. Paired Student's t-test
##
############################################################

t_result <- t.test(
  signature$Bradyzoite,
  signature$Tachyzoite,
  paired = TRUE)

t_result

############################################################
##
## 13. Wilcoxon signed-rank test
##
############################################################

wilcox_result <-
  wilcox.test(
    signature$Bradyzoite,
    signature$Tachyzoite,
    paired = TRUE
  )

wilcox_result

############################################################
##
## 14. Cohen's d
##
############################################################

cohens_d <-
  mean(signature$Difference) /
  sd(signature$Difference)

cohens_d

############################################################
##
## 15. Summary table
##
############################################################

results <- data.frame(
  
  Test = c(
    "Shapiro-Wilk",
    "Paired t-test",
    "Wilcoxon"
  ),
  
  Statistic = c(
    shapiro_result$statistic,
    t_result$statistic,
    wilcox_result$statistic
  ),
  
  P_value = c(
    shapiro_result$p.value,
    t_result$p.value,
    wilcox_result$p.value
  )
  
)

results

