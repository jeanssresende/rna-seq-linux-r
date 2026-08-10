# ==========================================================
#
# RNA-Seq Analysis using Linux and R
#
# Module 5
# Transcriptomic Signature Analysis
#
# Description:
# Quantify and compare Bradyzoite and Tachyzoite
# transcriptomic signatures in infected samples.
#
# Objectives:
# • Calculate transcriptomic signature scores
# • Compare stage-specific signatures
# • Apply paired statistical tests
# • Estimate effect size
# • Export reproducible results
#
# Dataset:
# Toxoplasma gondii RNA-Seq
#
# Author:
# Jean Resende
#
# ==========================================================


############################################################
##
## Biological Question
##
## Do infected samples exhibit a stronger
## Bradyzoite transcriptomic signature than
## Tachyzoite transcriptomic signature?
##
############################################################



############################################################
##
## 1. Load packages
##
############################################################

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)



############################################################
##
## 2. Load transcriptomic signature matrices
##
############################################################

load("Biomarker_Expression_Analysis/results/transcriptomic_signature_matrices.RData")

# Objects available:
# brady_matrix
# tachy_matrix



############################################################
##
## 3. Explore matrix dimensions
##
############################################################

dim(brady_matrix)

dim(tachy_matrix)

head(brady_matrix)

head(tachy_matrix)



############################################################
##
## 4. Explore expression distributions
##
############################################################

summary(as.vector(brady_matrix))

summary(as.vector(tachy_matrix))

range(brady_matrix)

range(tachy_matrix)



############################################################
##
## 5. Convert matrices to long format
##
############################################################

brady_long <-
  brady_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("Gene") |>
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "TPM"
  ) |>
  mutate(Stage = "Bradyzoite")

tachy_long <-
  tachy_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("Gene") |>
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "TPM"
  ) |>
  mutate(Stage = "Tachyzoite")

expression_long <- bind_rows(brady_long, tachy_long)



############################################################
##
## 6. Visualize global expression distributions
##
############################################################

ggplot(expression_long, aes(Stage, TPM, fill = Stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of transcriptomic signature expression",
    y = "TPM (log10 scale)"
  )



############################################################
##
## 7. Violin plot
##
############################################################

ggplot(expression_long, aes(Stage, TPM, fill = Stage)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(
    title = "Violin plot of transcriptomic signatures",
    y = "TPM (log10 scale)"
  )



############################################################
##
## 8. Calculate transcriptomic signature scores
##
############################################################

brady_score <- colMeans(brady_matrix)

tachy_score <- colMeans(tachy_matrix)



############################################################
##
## 9. Build signature table
##
############################################################

signature <- data.frame(
  Sample = names(brady_score),
  Bradyzoite = brady_score,
  Tachyzoite = tachy_score
)

signature



############################################################
##
## 10. Calculate paired differences
##
############################################################

signature$Difference <-
  signature$Bradyzoite - signature$Tachyzoite

signature



############################################################
##
## 11. Descriptive statistics
##
############################################################

summary(signature)

apply(signature[, 2:4], 2, mean)

apply(signature[, 2:4], 2, median)

apply(signature[, 2:4], 2, sd)



############################################################
##
## 12. Convert signature table to long format
##
############################################################

signature_long <-
  signature |>
  pivot_longer(
    cols = c(Bradyzoite, Tachyzoite),
    names_to = "Stage",
    values_to = "Expression"
  )



############################################################
##
## 13. Paired signature plot
##
############################################################

ggplot(signature_long,
       aes(Stage, Expression,
           group = Sample,
           color = Sample)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_classic(base_size = 14) +
  labs(
    title = "Paired transcriptomic signature comparison",
    y = "Mean TPM"
  )



############################################################
##
## 14. Difference plot
##
############################################################

ggplot(signature,
       aes(Sample, Difference, group = 1)) +
  geom_line(color = "steelblue") +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "red"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Bradyzoite - Tachyzoite signature difference",
    y = "Difference"
  )



############################################################
##
## 15. Shapiro-Wilk normality test
##
############################################################

shapiro_result <-
  shapiro.test(signature$Difference)

shapiro_result



############################################################
##
## 16. QQ Plot
##
############################################################

qqnorm(signature$Difference, pch = 19)

qqline(signature$Difference, col = "red", lwd = 2)



############################################################
##
## 17. Paired Student's t-test
##
############################################################

t_result <-
  t.test(
    signature$Bradyzoite,
    signature$Tachyzoite,
    paired = TRUE,
    conf.level = 0.95
  )

t_result



############################################################
##
## 18. Wilcoxon signed-rank test
##
############################################################

wilcox_result <-
  wilcox.test(
    signature$Bradyzoite,
    signature$Tachyzoite,
    paired = TRUE,
    exact = FALSE
  )

wilcox_result



############################################################
##
## 19. Effect size (Cohen's d)
##
############################################################

mean_difference <- mean(signature$Difference)

sd_difference <- sd(signature$Difference)

cohens_d <- mean_difference / sd_difference

cohens_d



############################################################
##
## 20. Interpret effect size
##
############################################################

if(abs(cohens_d) < 0.2){
  effect_interpretation <- "Negligible"
}else if(abs(cohens_d) < 0.5){
  effect_interpretation <- "Small"
}else if(abs(cohens_d) < 0.8){
  effect_interpretation <- "Medium"
}else{
  effect_interpretation <- "Large"
}

effect_interpretation



############################################################
##
## 21. Confidence interval
##
############################################################

t_result$conf.int



############################################################
##
## 22. Summary table
##
############################################################

results <- data.frame(
  Test = c(
    "Shapiro-Wilk",
    "Paired t-test",
    "Wilcoxon",
    "Cohen_d"
  ),
  Statistic = c(
    unname(shapiro_result$statistic),
    unname(t_result$statistic),
    unname(wilcox_result$statistic),
    cohens_d
  ),
  P_value = c(
    shapiro_result$p.value,
    t_result$p.value,
    wilcox_result$p.value,
    NA
  )
)

results



############################################################
##
## 23. Export results
##
############################################################

write.csv(
  signature,
  "../results/signature_scores.csv",
  row.names = FALSE
)

write.csv(
  results,
  "../results/signature_statistics.csv",
  row.names = FALSE
)



############################################################
##
## 24. Biological interpretation
##
############################################################

if(mean(signature$Difference) > 0){
  message(
    "Bradyzoite signature predominates in infected samples."
  )
}else{
  message(
    "Tachyzoite signature predominates in infected samples."
  )
}



############################################################
##
## End of Module 5
##
############################################################