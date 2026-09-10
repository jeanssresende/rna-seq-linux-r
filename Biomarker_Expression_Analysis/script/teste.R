# ============================================================
# AULA 08: Comparação de Marcadores Não-Correlacionáveis
#          com Dados Pareados
# ============================================================
# Objetivo: Comparar dois marcadores (DBA e SAG1) que medem
#           aspectos diferentes do parasita (quantidade vs
#           aglomerado), entre dois grupos biológicos
#           (diferenciação vs proliferação)
#
# Pergunta biológica:
#   "Como DBA e SAG1 variam entre os grupos?"
#   "Qual marcador melhor discrimina os grupos?"
#   "Os marcadores se correlacionam dentro de cada grupo?"
#
# Conceitos abordados:
#   - Dados pareados
#   - Marcadores não-correlacionáveis
#   - Padronização (Z-score) para comparação
#   - Testes estatísticos pareados
#   - Visualizações comparativas
# ============================================================

rm(list = ls())

setwd("Biomarker_Expression_Analysis/script")

library(tidyverse)
library(ggplot2)
library(gridExtra)

path_results <- "../results"
path_data <- "../data"

if (!dir.exists(path_results)) {
  dir.create(path_results, recursive = TRUE)
}

# ============================================================
# 1. CARREGAR E EXPLORAR DADOS
# ============================================================

# Ler arquivo (cuidado com separadores)
data_raw <- read.csv(file.path(path_data, "Integrated_density_normalized_Class_new.csv"),
                     header = TRUE, stringsAsFactors = FALSE)

head(data_raw)
dim(data_raw)
colnames(data_raw)

# Limpar nomes das colunas
colnames(data_raw) <- c("Sample_diff", "DBA_diff", "SAG1_diff", 
                        "Sample_prol", "DBA_prol", "SAG1_prol")

# Converter para numérico (remover vírgulas se necessário)
data_raw$DBA_diff <- as.numeric(gsub(",", ".", data_raw$DBA_diff))
data_raw$SAG1_diff <- as.numeric(gsub(",", ".", data_raw$SAG1_diff))
data_raw$DBA_prol <- as.numeric(gsub(",", ".", data_raw$DBA_prol))
data_raw$SAG1_prol <- as.numeric(gsub(",", ".", data_raw$SAG1_prol))

# Remover NAs
data_clean <- data_raw %>%
  filter(!is.na(DBA_diff), !is.na(SAG1_diff), 
         !is.na(DBA_prol), !is.na(SAG1_prol))

cat("Dados carregados:\n")
cat("  Dimensões:", dim(data_clean), "\n")
cat("  Amostras pareadas:", nrow(data_clean), "\n")

# ============================================================
# 2. EXPLORAR DISTRIBUIÇÃO DOS MARCADORES
# ============================================================

cat("\n--- Estatísticas descritivas ---\n")

# Grupo Diferenciação
cat("\nGrupo DIFERENCIAÇÃO:\n")
cat("  DBA:\n")
cat("    Média:", mean(data_clean$DBA_diff), "\n")
cat("    Mediana:", median(data_clean$DBA_diff), "\n")
cat("    SD:", sd(data_clean$DBA_diff), "\n")
cat("    Min:", min(data_clean$DBA_diff), "\n")
cat("    Max:", max(data_clean$DBA_diff), "\n")

cat("  SAG1:\n")
cat("    Média:", mean(data_clean$SAG1_diff), "\n")
cat("    Mediana:", median(data_clean$SAG1_diff), "\n")
cat("    SD:", sd(data_clean$SAG1_diff), "\n")
cat("    Min:", min(data_clean$SAG1_diff), "\n")
cat("    Max:", max(data_clean$SAG1_diff), "\n")

# Grupo Proliferação
cat("\nGrupo PROLIFERAÇÃO:\n")
cat("  DBA:\n")
cat("    Média:", mean(data_clean$DBA_prol), "\n")
cat("    Mediana:", median(data_clean$DBA_prol), "\n")
cat("    SD:", sd(data_clean$DBA_prol), "\n")
cat("    Min:", min(data_clean$DBA_prol), "\n")
cat("    Max:", max(data_clean$DBA_prol), "\n")

cat("  SAG1:\n")
cat("    Média:", mean(data_clean$SAG1_prol), "\n")
cat("    Mediana:", median(data_clean$SAG1_prol), "\n")
cat("    SD:", sd(data_clean$SAG1_prol), "\n")
cat("    Min:", min(data_clean$SAG1_prol), "\n")
cat("    Max:", max(data_clean$SAG1_prol), "\n")

# ============================================================
# 3. PROBLEMA: ESCALAS DIFERENTES
# ============================================================

cat("\n--- Problema: Escalas diferentes ---\n")

cat("\nDBA vs SAG1 (Diferenciação):\n")
cat("  Razão DBA/SAG1 (média):", mean(data_clean$DBA_diff) / mean(data_clean$SAG1_diff), "\n")
cat("  → DBA tem escala ~10x maior que SAG1\n")
cat("  → Não podem ser comparados diretamente\n")

# ============================================================
# 4. SOLUÇÃO: PADRONIZAR COM Z-SCORE
# ============================================================

cat("\n--- Padronizando com Z-score ---\n")

# Função Z-score
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Aplicar Z-score para cada marcador
data_zscore <- data_clean %>%
  mutate(
    # Z-score por marcador (não por grupo)
    DBA_diff_z = z_score(DBA_diff),
    SAG1_diff_z = z_score(SAG1_diff),
    DBA_prol_z = z_score(DBA_prol),
    SAG1_prol_z = z_score(SAG1_prol)
  )

cat("Dados padronizados com Z-score\n")
cat("  Média DBA_diff_z:", mean(data_zscore$DBA_diff_z), "\n")
cat("  SD DBA_diff_z:", sd(data_zscore$DBA_diff_z), "\n")

# ============================================================
# 5. COMPARAÇÃO DENTRO DE CADA GRUPO
# ============================================================

cat("\n--- Comparação dentro de cada grupo ---\n")

# Diferenciação: DBA vs SAG1
cat("\nGrupo DIFERENCIAÇÃO:\n")
cat("  DBA (z-score):\n")
cat("    Média:", mean(data_zscore$DBA_diff_z), "\n")
cat("    SD:", sd(data_zscore$DBA_diff_z), "\n")

cat("  SAG1 (z-score):\n")
cat("    Média:", mean(data_zscore$SAG1_diff_z), "\n")
cat("    SD:", sd(data_zscore$SAG1_diff_z), "\n")

# Teste pareado: DBA vs SAG1 em Diferenciação
test_diff <- t.test(data_zscore$DBA_diff_z, data_zscore$SAG1_diff_z, paired = TRUE)
cat("\n  Teste t pareado (DBA vs SAG1):\n")
cat("    t =", test_diff$statistic, "\n")
cat("    p-value =", test_diff$p.value, "\n")
if (test_diff$p.value < 0.05) {
  cat("    ✓ Significativo: DBA e SAG1 diferem em Diferenciação\n")
} else {
  cat("    ✗ Não significativo\n")
}

# Proliferação: DBA vs SAG1
cat("\nGrupo PROLIFERAÇÃO:\n")
cat("  DBA (z-score):\n")
cat("    Média:", mean(data_zscore$DBA_prol_z), "\n")
cat("    SD:", sd(data_zscore$DBA_prol_z), "\n")

cat("  SAG1 (z-score):\n")
cat("    Média:", mean(data_zscore$SAG1_prol_z), "\n")
cat("    SD:", sd(data_zscore$SAG1_prol_z), "\n")

# Teste pareado: DBA vs SAG1 em Proliferação
test_prol <- t.test(data_zscore$DBA_prol_z, data_zscore$SAG1_prol_z, paired = TRUE)
cat("\n  Teste t pareado (DBA vs SAG1):\n")
cat("    t =", test_prol$statistic, "\n")
cat("    p-value =", test_prol$p.value, "\n")
if (test_prol$p.value < 0.05) {
  cat("    ✓ Significativo: DBA e SAG1 diferem em Proliferação\n")
} else {
  cat("    ✗ Não significativo\n")
}

# ============================================================
# 6. COMPARAÇÃO ENTRE GRUPOS (MESMO MARCADOR)
# ============================================================

cat("\n--- Comparação entre grupos (mesmo marcador) ---\n")

# DBA: Diferenciação vs Proliferação
cat("\nMarcador DBA:\n")
cat("  Diferenciação (z-score):\n")
cat("    Média:", mean(data_zscore$DBA_diff_z), "\n")
cat("    SD:", sd(data_zscore$DBA_diff_z), "\n")
cat("    n =", sum(!is.na(data_zscore$DBA_diff_z)), "\n")

cat("  Proliferação (z-score):\n")
cat("    Média:", mean(data_zscore$DBA_prol_z), "\n")
cat("    SD:", sd(data_zscore$DBA_prol_z), "\n")
cat("    n =", sum(!is.na(data_zscore$DBA_prol_z)), "\n")

# Teste NÃO-PAREADO: DBA em Diferenciação vs Proliferação
test_dba <- t.test(data_zscore$DBA_diff_z, data_zscore$DBA_prol_z, paired = FALSE)
cat("\n  Teste t NÃO-PAREADO (Diferenciação vs Proliferação):\n")
cat("    t =", test_dba$statistic, "\n")
cat("    p-value =", test_dba$p.value, "\n")
cat("    Diferença de médias:", test_dba$estimate[1] - test_dba$estimate[2], "\n")
if (test_dba$p.value < 0.05) {
  cat("    ✓ Significativo: DBA varia entre grupos\n")
} else {
  cat("    ✗ Não significativo\n")
}

# SAG1: Diferenciação vs Proliferação
cat("\nMarcador SAG1:\n")
cat("  Diferenciação (z-score):\n")
cat("    Média:", mean(data_zscore$SAG1_diff_z), "\n")
cat("    SD:", sd(data_zscore$SAG1_diff_z), "\n")
cat("    n =", sum(!is.na(data_zscore$SAG1_diff_z)), "\n")

cat("  Proliferação (z-score):\n")
cat("    Média:", mean(data_zscore$SAG1_prol_z), "\n")
cat("    SD:", sd(data_zscore$SAG1_prol_z), "\n")
cat("    n =", sum(!is.na(data_zscore$SAG1_prol_z)), "\n")

# Teste NÃO-PAREADO: SAG1 em Diferenciação vs Proliferação
test_sag1 <- t.test(data_zscore$SAG1_diff_z, data_zscore$SAG1_prol_z, paired = FALSE)
cat("\n  Teste t NÃO-PAREADO (Diferenciação vs Proliferação):\n")
cat("    t =", test_sag1$statistic, "\n")
cat("    p-value =", test_sag1$p.value, "\n")
cat("    Diferença de médias:", test_sag1$estimate[1] - test_sag1$estimate[2], "\n")
if (test_sag1$p.value < 0.05) {
  cat("    ✓ Significativo: SAG1 varia entre grupos\n")
} else {
  cat("    ✗ Não significativo\n")
}
# ============================================================
# 7. CORRELAÇÃO ENTRE MARCADORES
# ============================================================

cat("\n--- Correlação entre marcadores ---\n")

# Correlação em Diferenciação
cor_diff <- cor.test(data_zscore$DBA_diff_z, data_zscore$SAG1_diff_z)
cat("\nGrupo DIFERENCIAÇÃO:\n")
cat("  Correlação (DBA vs SAG1):", cor_diff$estimate, "\n")
cat("  p-value:", cor_diff$p.value, "\n")
if (cor_diff$p.value < 0.05) {
  cat("  ✓ Correlação significativa\n")
} else {
  cat("  ✗ Sem correlação significativa\n")
}

# Correlação em Proliferação
cor_prol <- cor.test(data_zscore$DBA_prol_z, data_zscore$SAG1_prol_z)
cat("\nGrupo PROLIFERAÇÃO:\n")
cat("  Correlação (DBA vs SAG1):", cor_prol$estimate, "\n")
cat("  p-value:", cor_prol$p.value, "\n")
if (cor_prol$p.value < 0.05) {
  cat("  ✓ Correlação significativa\n")
} else {
  cat("  ✗ Sem correlação significativa\n")
}

# ============================================================
# 8. VISUALIZAR DISTRIBUIÇÕES
# ============================================================

cat("\n--- Gerando gráficos ---\n")

# Gráfico 1: Boxplot - Comparação entre marcadores (Diferenciação)
pdf(file.path(path_results, "01_markers_diff_boxplot.pdf"), width = 10, height = 6)

data_boxplot_diff <- data_zscore %>%
  select(DBA_diff_z, SAG1_diff_z) %>%
  pivot_longer(everything(), names_to = "Marcador", values_to = "Z_score") %>%
  mutate(Marcador = gsub("_diff_z", "", Marcador))

ggplot(data_boxplot_diff, aes(x = Marcador, y = Z_score, fill = Marcador)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribuição de Marcadores - Grupo DIFERENCIAÇÃO",
       x = "Marcador",
       y = "Z-score",
       subtitle = "Dados pareados") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

dev.off()
cat("Gráfico salvo: 01_markers_diff_boxplot.pdf\n")

# Gráfico 2: Boxplot - Comparação entre marcadores (Proliferação)
pdf(file.path(path_results, "02_markers_prol_boxplot.pdf"), width = 10, height = 6)

data_boxplot_prol <- data_zscore %>%
  select(DBA_prol_z, SAG1_prol_z) %>%
  pivot_longer(everything(), names_to = "Marcador", values_to = "Z_score") %>%
  mutate(Marcador = gsub("_prol_z", "", Marcador))

ggplot(data_boxplot_prol, aes(x = Marcador, y = Z_score, fill = Marcador)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribuição de Marcadores - Grupo PROLIFERAÇÃO",
       x = "Marcador",
       y = "Z-score",
       subtitle = "Dados pareados") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

dev.off()
cat("Gráfico salvo: 02_markers_prol_boxplot.pdf\n")

# Gráfico 3: Comparação entre grupos (DBA)
pdf(file.path(path_results, "03_dba_between_groups.pdf"), width = 10, height = 6)

data_dba_groups <- data.frame(
  Z_score = c(data_zscore$DBA_diff_z, data_zscore$DBA_prol_z),
  Grupo = c(rep("Diferenciação", nrow(data_zscore)), 
            rep("Proliferação", nrow(data_zscore)))
)

ggplot(data_dba_groups, aes(x = Grupo, y = Z_score, fill = Grupo)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Marcador DBA - Comparação entre Grupos",
       x = "Grupo",
       y = "Z-score",
       subtitle = "Dados pareados") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

dev.off()
cat("Gráfico salvo: 03_dba_between_groups.pdf\n")

# Gráfico 4: Comparação entre grupos (SAG1)
pdf(file.path(path_results, "04_sag1_between_groups.pdf"), width = 10, height = 6)

data_sag1_groups <- data.frame(
  Z_score = c(data_zscore$SAG1_diff_z, data_zscore$SAG1_prol_z),
  Grupo = c(rep("Diferenciação", nrow(data_zscore)), 
            rep("Proliferação", nrow(data_zscore)))
)

ggplot(data_sag1_groups, aes(x = Grupo, y = Z_score, fill = Grupo)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Marcador SAG1 - Comparação entre Grupos",
       x = "Grupo",
       y = "Z-score",
       subtitle = "Dados pareados") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

dev.off()
cat("Gráfico salvo: 04_sag1_between_groups.pdf\n")

# Gráfico 5: Scatter plot - Correlação em Diferenciação
pdf(file.path(path_results, "05_correlation_diff.pdf"), width = 8, height = 8)

ggplot(data_zscore, aes(x = DBA_diff_z, y = SAG1_diff_z)) +
  geom_point(size = 3, alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.2) +
  theme_minimal() +
  labs(title = "Correlação entre DBA e SAG1 - DIFERENCIAÇÃO",
       x = "DBA (Z-score)",
       y = "SAG1 (Z-score)",
       subtitle = paste("r =", round(cor_diff$estimate, 3), 
                        ", p-value =", round(cor_diff$p.value, 4))) +
  theme(plot.subtitle = element_text(size = 12, face = "italic"))

dev.off()
cat("Gráfico salvo: 05_correlation_diff.pdf\n")

# Gráfico 6: Scatter plot - Correlação em Proliferação
pdf(file.path(path_results, "06_correlation_prol.pdf"), width = 8, height = 8)

ggplot(data_zscore, aes(x = DBA_prol_z, y = SAG1_prol_z)) +
  geom_point(size = 3, alpha = 0.6, color = "coral") +
  geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.2) +
  theme_minimal() +
  labs(title = "Correlação entre DBA e SAG1 - PROLIFERAÇÃO",
       x = "DBA (Z-score)",
       y = "SAG1 (Z-score)",
       subtitle = paste("r =", round(cor_prol$estimate, 3), 
                        ", p-value =", round(cor_prol$p.value, 4))) +
  theme(plot.subtitle = element_text(size = 12, face = "italic"))

dev.off()
cat("Gráfico salvo: 06_correlation_prol.pdf\n")

# Gráfico 7: Paired plot - Mudança de cada amostra
pdf(file.path(path_results, "07_paired_changes.pdf"), width = 12, height = 8)

par(mfrow = c(2, 2))

# DBA
plot(1:nrow(data_zscore), data_zscore$DBA_diff_z, 
     main = "DBA: Mudança de Diferenciação para Proliferação",
     xlab = "Amostra pareada",
     ylab = "Z-score",
     pch = 16, col = "steelblue", ylim = c(min(data_zscore$DBA_diff_z, data_zscore$DBA_prol_z),
                                           max(data_zscore$DBA_diff_z, data_zscore$DBA_prol_z)))
points(1:nrow(data_zscore), data_zscore$DBA_prol_z, pch = 16, col = "coral")
for (i in 1:nrow(data_zscore)) {
  lines(c(i, i), c(data_zscore$DBA_diff_z[i], data_zscore$DBA_prol_z[i]), 
        col = "gray", lty = 2)
}
legend("topright", legend = c("Diferenciação", "Proliferação"), 
       col = c("steelblue", "coral"), pch = 16)

# SAG1
plot(1:nrow(data_zscore), data_zscore$SAG1_diff_z, 
     main = "SAG1: Mudança de Diferenciação para Proliferação",
     xlab = "Amostra pareada",
     ylab = "Z-score",
     pch = 16, col = "steelblue", ylim = c(min(data_zscore$SAG1_diff_z, data_zscore$SAG1_prol_z),
                                           max(data_zscore$SAG1_diff_z, data_zscore$SAG1_prol_z)))
points(1:nrow(data_zscore), data_zscore$SAG1_prol_z, pch = 16, col = "coral")
for (i in 1:nrow(data_zscore)) {
  lines(c(i, i), c(data_zscore$SAG1_diff_z[i], data_zscore$SAG1_prol_z[i]), 
        col = "gray", lty = 2)
}
legend("topright", legend = c("Diferenciação", "Proliferação"), 
       col = c("steelblue", "coral"), pch = 16)

par(mfrow = c(1, 1))

dev.off()
cat("Gráfico salvo: 07_paired_changes.pdf\n")

# ============================================================
# 9. EXPORTAR RESULTADOS
# ============================================================

cat("\n--- Exportando resultados ---\n")

# Tabela com dados originais e padronizados
results_table <- data_zscore %>%
  select(Sample_diff, DBA_diff, SAG1_diff, DBA_diff_z, SAG1_diff_z,
         Sample_prol, DBA_prol, SAG1_prol, DBA_prol_z, SAG1_prol_z)

write.csv(results_table, 
          file = file.path(path_results, "markers_analysis_results.csv"), 
          row.names = FALSE)
cat("Tabela salva: markers_analysis_results.csv\n")

# Tabela de testes estatísticos
stats_tests <- data.frame(
  Comparacao = c("DBA vs SAG1 (Diferenciação)",
                 "DBA vs SAG1 (Proliferação)",
                 "DBA (Diferenciação vs Proliferação)",
                 "SAG1 (Diferenciação vs Proliferação)",
                 "Correlação DBA-SAG1 (Diferenciação)",
                 "Correlação DBA-SAG1 (Proliferação)"),
  Estatistica = c(test_diff$statistic, test_prol$statistic, 
                  test_dba$statistic, test_sag1$statistic,
                  cor_diff$estimate, cor_prol$estimate),
  P_value = c(test_diff$p.value, test_prol$p.value,
              test_dba$p.value, test_sag1$p.value,
              cor_diff$p.value, cor_prol$p.value),
  Significativo = c(ifelse(test_diff$p.value < 0.05, "Sim", "Não"),
                    ifelse(test_prol$p.value < 0.05, "Sim", "Não"),
                    ifelse(test_dba$p.value < 0.05, "Sim", "Não"),
                    ifelse(test_sag1$p.value < 0.05, "Sim", "Não"),
                    ifelse(cor_diff$p.value < 0.05, "Sim", "Não"),
                    ifelse(cor_prol$p.value < 0.05, "Sim", "Não"))
)

write.csv(stats_tests, 
          file = file.path(path_results, "statistical_tests.csv"), 
          row.names = FALSE)
cat("Tabela salva: statistical_tests.csv\n")

# ============================================================
# 10. RESUMO FINAL
# ============================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("RESUMO DA AULA 08\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Pergunta biológica:\n")
cat("  'Como DBA e SAG1 variam entre diferenciação e proliferação?'\n\n")

cat("Dados:\n")
cat("  Amostras pareadas:", nrow(data_zscore), "\n")
cat("  Marcadores: DBA (densidade), SAG1 (aglomerado)\n")
cat("  Grupos: Diferenciação vs Proliferação\n\n")

cat("Conceitos-chave:\n")
cat("  1. Dados pareados: mesma amostra em dois estados\n")
cat("  2. Marcadores não-correlacionáveis: escalas diferentes\n")
cat("  3. Solução: Z-score para padronização\n")
cat("  4. Análise: testes t pareados\n\n")

cat("Resultados:\n")
if (test_dba$p.value < 0.05) {
  cat("  ✓ DBA varia significativamente entre grupos (p =", 
      round(test_dba$p.value, 4), ")\n")
} else {
  cat("  ✗ DBA não varia significativamente entre grupos (p =", 
      round(test_dba$p.value, 4), ")\n")
}

if (test_sag1$p.value < 0.05) {
  cat("  ✓ SAG1 varia significativamente entre grupos (p =", 
      round(test_sag1$p.value, 4), ")\n")
} else {
  cat("  ✗ SAG1 não varia significativamente entre grupos (p =", 
      round(test_sag1$p.value, 4), ")\n")
}

if (cor_diff$p.value < 0.05) {
  cat("  ✓ DBA e SAG1 correlacionam em Diferenciação (r =", 
      round(cor_diff$estimate, 3), ")\n")
} else {
  cat("  ✗ Sem correlação em Diferenciação\n")
}

if (cor_prol$p.value < 0.05) {
  cat("  ✓ DBA e SAG1 correlacionam em Proliferação (r =", 
      round(cor_prol$estimate, 3), ")\n")
} else {
  cat("  ✗ Sem correlação em Proliferação\n")
}

cat("\nArquivos gerados:\n")
cat("  - markers_analysis_results.csv\n")
cat("  - statistical_tests.csv\n")
cat("  - 01_markers_diff_boxplot.pdf\n")
cat("  - 02_markers_prol_boxplot.pdf\n")
cat("  - 03_dba_between_groups.pdf\n")
cat("  - 04_sag1_between_groups.pdf\n")
cat("  - 05_correlation_diff.pdf\n")
cat("  - 06_correlation_prol.pdf\n")
cat("  - 07_paired_changes.pdf\n\n")

cat("Fim da aula 08.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")