# ============================================================
# AULA 06: Z-score, Heatmap e Assinatura Transcriptômica Normalizada
# ============================================================

rm(list = ls())

setwd("/home/jean/Documentos/Projetos/rna-seq-linux-r/Biomarker_Expression_Analysis/script")

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

path_results <- "../results"
path_data <- "/home/jean/Documentos/Projetos/rna-seq-linux-r/Biomarker_Expression_Analysis/Data"

if (!dir.exists(path_results)) {
  dir.create(path_results, recursive = TRUE)
}

# ============================================================
# 2. CARREGAR DADOS
# ============================================================

load(file.path(path_data, "se_toxo_gene.RData"))
load(file.path(path_results, "transcriptomic_signature_matrices.RData"))

tpm_matrix <- assay(se_gene, "abundance")
infected_samples <- c("quant_INF1", "quant_INF2", "quant_INF3")

tpm_infected <- tpm_matrix[, infected_samples]

genes_brady <- rownames(brady_matrix)
genes_tachy <- rownames(tachy_matrix)
genes_all <- unique(c(genes_brady, genes_tachy))

genes_all <- intersect(genes_all, rownames(tpm_infected))
tpm_subset <- tpm_infected[genes_all, ]

# ============================================================
# 3. EXPLORAR DISTRIBUIÇÃO DOS TPMs BRUTOS
# ============================================================

pdf(file.path(path_results, "01_tpm_distribution_raw.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))

boxplot(tpm_subset,
        main = "Distribuição de TPM por Amostra (Bruto)",
        ylab = "TPM",
        las = 2,
        cex.axis = 1.2)

hist(tpm_subset, breaks = 50,
     main = "Histograma dos TPMs (Bruto)",
     xlab = "TPM",
     ylab = "Frequência",
     col = "skyblue")

par(mfrow = c(1, 1))
dev.off()

# ============================================================
# 4. TRANSFORMAÇÃO LOGARÍTMICA
# ============================================================

log_tpm <- log2(tpm_subset + 1)

pdf(file.path(path_results, "02_tpm_distribution_log.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))

boxplot(log_tpm,
        main = "Distribuição de log2(TPM+1) por Amostra",
        ylab = "log2(TPM+1)",
        las = 2,
        cex.axis = 1.2)

hist(log_tpm, breaks = 50,
     main = "Histograma dos log2(TPM+1)",
     xlab = "log2(TPM+1)",
     ylab = "Frequência",
     col = "lightcoral")

par(mfrow = c(1, 1))
dev.off()

# ============================================================
# 5. IDENTIFICAR GENES SEM VARIABILIDADE
# ============================================================

sd_log_genes <- apply(log_tpm, 1, sd)
threshold_sd <- 1e-10
genes_no_var <- which(sd_log_genes < threshold_sd)

log_tpm_filtered <- log_tpm[-genes_no_var, ]

# ============================================================
# 6. CRIAR FUNÇÃO Z-SCORE CUSTOMIZADA
# ============================================================

z_score <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  
  if (is.na(sd_x) || sd_x == 0) {
    return(rep(NaN, length(x)))
  }
  
  z <- (x - mean_x) / sd_x
  return(z)
}

# ============================================================
# 7. APLICAR Z-SCORE GENE POR GENE
# ============================================================

# Aplicar z_score sem transpor
zscore_list <- apply(log_tpm_filtered, 1, z_score)

# zscore_list é uma matriz 3 x 720 (amostras x genes)
# Transpor para genes x amostras
zscore_matrix <- t(zscore_list)

# Restaurar nomes
rownames(zscore_matrix) <- rownames(log_tpm_filtered)
colnames(zscore_matrix) <- colnames(log_tpm_filtered)

# ============================================================
# 8. VALIDAR A TRANSFORMAÇÃO Z-SCORE
# ============================================================

pdf(file.path(path_results, "03_zscore_validation.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))

hist(zscore_matrix, breaks = 50,
     main = "Distribuição dos Z-scores",
     xlab = "Z-score",
     ylab = "Frequência",
     col = "lightgreen")
abline(v = 0, col = "red", lty = 2, lwd = 2)

boxplot(zscore_matrix,
        main = "Z-scores por Amostra",
        ylab = "Z-score",
        las = 2,
        cex.axis = 1.2)
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1))
dev.off()

# Verificação rápida
cat("Validação do Z-score:\n")
cat("Média:", round(mean(zscore_matrix), 16), "\n")
cat("SD:", round(mean(apply(zscore_matrix, 1, sd)), 4), "\n")
cat("Z-score calculado corretamente!\n")

# ============================================================
# 9. CONSTRUIR HEATMAP (COM TOP GENES MAIS VARIADOS)
# ============================================================

# Calcular variabilidade (SD) dos Z-scores por gene
sd_zscore <- apply(zscore_matrix, 1, sd, na.rm = TRUE)

# Selecionar top 25 genes mais variados
top_n <- 25
top_genes <- names(sort(sd_zscore, decreasing = TRUE)[1:top_n])
zscore_matrix_top <- zscore_matrix[top_genes, ]

color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

pdf(file.path(path_results, "04_heatmap_zscore.pdf"), width = 12, height = 10)
pheatmap(zscore_matrix_top,
         main = paste("Heatmap de Z-scores: Top", top_n, "Genes Mais Variados"),
         color = color_palette,
         breaks = seq(-3, 3, length.out = 51),
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize = 9,
         fontsize_row = 10,
         fontsize_col = 14,
         cellwidth = 50,
         cellheight = 8)
dev.off()

# ============================================================
# 10. CALCULAR SCORES DE ASSINATURA
# ============================================================

genes_brady_filtered <- intersect(genes_brady, rownames(zscore_matrix))
genes_tachy_filtered <- intersect(genes_tachy, rownames(zscore_matrix))

bradyzoite_score <- colMeans(zscore_matrix[genes_brady_filtered, ], na.rm = TRUE)
tachyzoite_score <- colMeans(zscore_matrix[genes_tachy_filtered, ], na.rm = TRUE)

# ============================================================
# 11. CALCULAR STAGE INDEX
# ============================================================

stage_index <- bradyzoite_score - tachyzoite_score

# ============================================================
# 12. VISUALIZAR SCORES
# ============================================================

scores_df <- data.frame(
  Sample = names(bradyzoite_score),
  Bradyzoite = bradyzoite_score,
  Tachyzoite = tachyzoite_score,
  Stage_Index = stage_index,
  row.names = NULL
)

pdf(file.path(path_results, "05_scores_comparison.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

barplot(rbind(scores_df$Bradyzoite, scores_df$Tachyzoite),
        beside = TRUE,
        names.arg = scores_df$Sample,
        legend.text = c("Bradyzoite", "Tachyzoite"),
        main = "Scores de Assinatura por Amostra",
        ylab = "Score (média dos Z-scores)",
        col = c("steelblue", "coral"),
        cex.axis = 1.2,
        cex.names = 1.2)
abline(h = 0, col = "black", lty = 2)

barplot(scores_df$Stage_Index,
        names.arg = scores_df$Sample,
        main = "Stage Index por Amostra",
        ylab = "Stage Index",
        col = ifelse(scores_df$Stage_Index > 0, "steelblue", "coral"),
        cex.axis = 1.2,
        cex.names = 1.2)
abline(h = 0, col = "black", lty = 2)

par(mfrow = c(1, 1))
dev.off()

# ============================================================
# 13. EXPORTAR RESULTADOS
# ============================================================

save(zscore_matrix, file = file.path(path_results, "zscore_matrix.RData"))

write.csv(scores_df, file = file.path(path_results, "scores_summary.csv"), row.names = FALSE)

zscore_stats <- data.frame(
  Gene = rownames(zscore_matrix),
  Mean_Zscore = rowMeans(zscore_matrix, na.rm = TRUE),
  SD_Zscore = apply(zscore_matrix, 1, sd, na.rm = TRUE),
  Min_Zscore = apply(zscore_matrix, 1, min, na.rm = TRUE),
  Max_Zscore = apply(zscore_matrix, 1, max, na.rm = TRUE),
  Assinatura = ifelse(rownames(zscore_matrix) %in% genes_brady_filtered, "Bradyzoite",
                      ifelse(rownames(zscore_matrix) %in% genes_tachy_filtered, "Tachyzoite", "Ambas"))
)

write.csv(zscore_stats, file = file.path(path_results, "zscore_statistics.csv"), row.names = FALSE)

# ============================================================
# 14. ANÁLISE ESTATÍSTICA (EXPLORATÓRIA)
# ============================================================

cat("\n=== ANÁLISE DESCRITIVA DOS SCORES ===\n\n")

# Resumo descritivo
cat("Bradyzoite Score:\n")
print(summary(scores_df$Bradyzoite))

cat("\nTachyzoite Score:\n")
print(summary(scores_df$Tachyzoite))

cat("\nStage Index:\n")
print(summary(scores_df$Stage_Index))

# Diferença pareada
diferenca <- scores_df$Bradyzoite - scores_df$Tachyzoite
cat("\nDiferença (Brady - Tachy):\n")
print(diferenca)
cat("Média da diferença:", mean(diferenca), "\n")
cat("SD da diferença:", sd(diferenca), "\n")

# Teste de Wilcoxon (não-paramétrico)
cat("\n=== TESTE NÃO-PARAMÉTRICO ===\n")
wilcox_result <- wilcox.test(scores_df$Bradyzoite, 
                             scores_df$Tachyzoite, 
                             paired = TRUE)
print(wilcox_result)

cat("\n LIMITAÇÕES (n=3):\n")
cat("- Poder estatístico muito baixo\n")
cat("- Resultados exploratórios apenas\n")
cat("- Validação experimental necessária\n")

# Análise por gene
gene_stats <- data.frame(
  Gene = rownames(zscore_matrix),
  Mean = rowMeans(zscore_matrix),
  SD = apply(zscore_matrix, 1, sd),
  Range = apply(zscore_matrix, 1, max) - apply(zscore_matrix, 1, min)
)

gene_stats_sorted <- gene_stats[order(gene_stats$Range, decreasing = TRUE), ]

cat("\nTOP 15 GENES MAIS INFORMATIVOS:\n")
print(head(gene_stats_sorted, 15))

write.csv(gene_stats_sorted, 
          file = file.path(path_results, "gene_variability_ranking.csv"), 
          row.names = FALSE)