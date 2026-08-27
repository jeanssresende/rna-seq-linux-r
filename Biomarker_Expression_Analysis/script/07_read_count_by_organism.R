# ============================================================
# AULA 07: Contagem de Reads por Gene e Organismo
# ============================================================
# Objetivo: Quantificar reads mapeados para genes de cada
#           organismo (Homo sapiens vs Toxoplasma gondii)
#
# Pergunta biológica:
#   "Qual é a proporção de reads do hospedeiro vs parasita?"
#   "Qual amostra tem maior carga parasitária?"
#
# Conceitos abordados:
#   - Contagem de reads (assay)
#   - Identificação de organismo por gene
#   - Agregação de reads por organismo
#   - Visualização de proporções
# ============================================================

rm(list = ls())

setwd("/home/jean/Documentos/Projetos/rna-seq-linux-r/Biomarker_Expression_Analysis/script")

library(tidyverse)
library(SummarizedExperiment)

path_results <- "../results"
path_data <- "/home/jean/Documentos/Projetos/rna-seq-linux-r/Biomarker_Expression_Analysis/Data"

if (!dir.exists(path_results)) {
  dir.create(path_results, recursive = TRUE)
}

# ============================================================
# 1. CARREGAR DADOS
# ============================================================

# Carrega SummarizedExperiment com dados de Homo sapiens
load(file.path(path_data, "se_gene.RData"))
se_human <- se_gene
cat("SummarizedExperiment Homo sapiens carregado\n")
cat("  Dimensões:", dim(se_human)[1], "genes x", dim(se_human)[2], "amostras\n")

# Carrega SummarizedExperiment com dados de Toxoplasma gondii
load(file.path(path_data, "se_toxo_gene.RData"))
se_toxo <- se_gene
cat("\nSummarizedExperiment Toxoplasma gondii carregado\n")
cat("  Dimensões:", dim(se_toxo)[1], "genes x", dim(se_toxo)[2], "amostras\n")

# Extrair matrizes de contagem
counts_human <- assay(se_human, "counts")
counts_toxo <- assay(se_toxo, "counts")

cat("\nMatrizes de contagem extraídas\n")
cat("  Homo sapiens:", dim(counts_human), "\n")
cat("  Toxoplasma gondii:", dim(counts_toxo), "\n")

# ============================================================
# 2. EXPLORAR ESTRUTURA DOS DADOS
# ============================================================

cat("\n--- Explorando estrutura dos dados ---\n")

# Homo sapiens
cat("\nHomo sapiens:\n")
cat("  Nomes dos genes (primeiros 10):\n")
print(head(rownames(counts_human), 10))

cat("\n  Amostras:\n")
print(colnames(counts_human))

cat("\n  Estatísticas de contagem:\n")
cat("    Contagem total:", sum(counts_human), "reads\n")
cat("    Média por gene:", mean(rowSums(counts_human)), "\n")
cat("    Mediana por gene:", median(rowSums(counts_human)), "\n")

# Toxoplasma gondii
cat("\nToxoplasma gondii:\n")
cat("  Nomes dos genes (primeiros 10):\n")
print(head(rownames(counts_toxo), 10))

cat("\n  Amostras:\n")
print(colnames(counts_toxo))

cat("\n  Estatísticas de contagem:\n")
cat("    Contagem total:", sum(counts_toxo), "reads\n")
cat("    Média por gene:", mean(rowSums(counts_toxo)), "\n")
cat("    Mediana por gene:", median(rowSums(counts_toxo)), "\n")

# ============================================================
# 3. CONTAGEM TOTAL DE READS POR ORGANISMO E AMOSTRA
# ============================================================

cat("\n--- Contagem total de reads por organismo e amostra ---\n")

# Por amostra
reads_human_per_sample <- colSums(counts_human)
reads_toxo_per_sample <- colSums(counts_toxo)

# Criar tabela
read_count_summary <- data.frame(
  Sample = colnames(counts_human),
  Homo_sapiens = reads_human_per_sample,
  Toxoplasma_gondii = reads_toxo_per_sample,
  Total = reads_human_per_sample + reads_toxo_per_sample,
  Percent_Human = (reads_human_per_sample / (reads_human_per_sample + reads_toxo_per_sample)) * 100,
  Percent_Toxo = (reads_toxo_per_sample / (reads_human_per_sample + reads_toxo_per_sample)) * 100
)

row.names(read_count_summary) <- NULL

cat("\nTabela de contagem de reads:\n")
print(read_count_summary)

write.csv(read_count_summary, 
          file = file.path(path_results, "read_count_summary.csv"), 
          row.names = FALSE)
cat("\nTabela salva: read_count_summary.csv\n")

# ============================================================
# 4. CONTAGEM TOTAL DE READS POR ORGANISMO (AGREGADO)
# ============================================================

cat("\n--- Contagem total agregada por organismo ---\n")

total_human <- sum(counts_human)
total_toxo <- sum(counts_toxo)
total_all <- total_human + total_toxo

cat("Homo sapiens:\n")
cat("  Total de reads:", total_human, "\n")
cat("  Percentual:", (total_human / total_all) * 100, "%\n")

cat("\nToxoplasma gondii:\n")
cat("  Total de reads:", total_toxo, "\n")
cat("  Percentual:", (total_toxo / total_all) * 100, "%\n")

cat("\nTotal geral:\n")
cat("  Reads:", total_all, "\n")

# ============================================================
# 5. CONTAGEM POR GENE (TOP 20 GENES MAIS EXPRESSOS)
# ============================================================

cat("\n--- Top 20 genes mais expressos (Homo sapiens) ---\n")

# Soma de reads por gene
reads_per_gene_human <- rowSums(counts_human)
top_genes_human <- sort(reads_per_gene_human, decreasing = TRUE)[1:20]

top_genes_human_df <- data.frame(
  Gene = names(top_genes_human),
  Total_Reads = as.numeric(top_genes_human),
  Percent_of_Total = (as.numeric(top_genes_human) / total_human) * 100
)

cat("\n")
print(top_genes_human_df)

write.csv(top_genes_human_df, 
          file = file.path(path_results, "top_genes_human.csv"), 
          row.names = FALSE)

cat("\n--- Top 20 genes mais expressos (Toxoplasma gondii) ---\n")

reads_per_gene_toxo <- rowSums(counts_toxo)
top_genes_toxo <- sort(reads_per_gene_toxo, decreasing = TRUE)[1:20]

top_genes_toxo_df <- data.frame(
  Gene = names(top_genes_toxo),
  Total_Reads = as.numeric(top_genes_toxo),
  Percent_of_Total = (as.numeric(top_genes_toxo) / total_toxo) * 100
)

cat("\n")
print(top_genes_toxo_df)

write.csv(top_genes_toxo_df, 
          file = file.path(path_results, "top_genes_toxo.csv"), 
          row.names = FALSE)

# ============================================================
# 6. VISUALIZAR PROPORÇÃO DE READS POR AMOSTRA
# ============================================================

cat("\n--- Gerando gráficos ---\n")

# Gráfico 1: Barras empilhadas
pdf(file.path(path_results, "01_read_count_stacked_bar.pdf"), width = 10, height = 6)

data_plot <- read_count_summary[, c("Sample", "Homo_sapiens", "Toxoplasma_gondii")]
rownames(data_plot) <- data_plot$Sample
data_plot$Sample <- NULL

barplot(t(data_plot),
        main = "Contagem Total de Reads por Amostra e Organismo",
        xlab = "Amostra",
        ylab = "Número de Reads",
        legend.text = c("Homo sapiens", "Toxoplasma gondii"),
        col = c("steelblue", "coral"),
        beside = FALSE)

dev.off()
cat("Gráfico salvo: 01_read_count_stacked_bar.pdf\n")

# Gráfico 2: Proporção (%)
pdf(file.path(path_results, "02_read_count_proportion.pdf"), width = 10, height = 6)

data_prop <- read_count_summary[, c("Sample", "Percent_Human", "Percent_Toxo")]
rownames(data_prop) <- data_prop$Sample
data_prop$Sample <- NULL

barplot(t(data_prop),
        main = "Proporção de Reads por Organismo (%)",
        xlab = "Amostra",
        ylab = "Percentual (%)",
        legend.text = c("Homo sapiens", "Toxoplasma gondii"),
        col = c("steelblue", "coral"),
        beside = FALSE,
        ylim = c(0, 100))

dev.off()
cat("Gráfico salvo: 02_read_count_proportion.pdf\n")

# Gráfico 3: Top genes Homo sapiens
pdf(file.path(path_results, "03_top_genes_human.pdf"), width = 10, height = 8)

barplot(top_genes_human_df$Total_Reads,
        names.arg = top_genes_human_df$Gene,
        main = "Top 20 Genes Mais Expressos - Homo sapiens",
        xlab = "Gene",
        ylab = "Contagem de Reads",
        col = "steelblue",
        las = 2)

dev.off()
cat("Gráfico salvo: 03_top_genes_human.pdf\n")

# Gráfico 4: Top genes Toxoplasma gondii
pdf(file.path(path_results, "04_top_genes_toxo.pdf"), width = 10, height = 8)

barplot(top_genes_toxo_df$Total_Reads,
        names.arg = top_genes_toxo_df$Gene,
        main = "Top 20 Genes Mais Expressos - Toxoplasma gondii",
        xlab = "Gene",
        ylab = "Contagem de Reads",
        col = "coral",
        las = 2)

dev.off()
cat("Gráfico salvo: 04_top_genes_toxo.pdf\n")

# ============================================================
# 7. CRIAR TABELA COMPARATIVA COMPLETA
# ============================================================

cat("\n--- Criando tabela comparativa por gene e amostra ---\n")

# Combinar informações de ambos os organismos
comparison_table <- data.frame(
  Gene_Human = c(rownames(counts_human), rep(NA, nrow(counts_toxo))),
  Gene_Toxo = c(rep(NA, nrow(counts_human)), rownames(counts_toxo)),
  Organism = c(rep("Homo sapiens", nrow(counts_human)), 
               rep("Toxoplasma gondii", nrow(counts_toxo))),
  Total_Reads = c(rowSums(counts_human), rowSums(counts_toxo))
)

# Adicionar contagens por amostra
for (sample in colnames(counts_human)) {
  comparison_table[[paste0("Reads_", sample)]] <- c(
    counts_human[, sample],
    counts_toxo[, sample]
  )
}

# Ordenar por total de reads
comparison_table <- comparison_table[order(comparison_table$Total_Reads, decreasing = TRUE), ]

cat("\nTabela de comparação (primeiras 20 linhas):\n")
print(head(comparison_table, 20))

write.csv(comparison_table, 
          file = file.path(path_results, "gene_read_count_comparison.csv"), 
          row.names = FALSE)
cat("\nTabela salva: gene_read_count_comparison.csv\n")

# ============================================================
# 8. ESTATÍSTICAS DETALHADAS
# ============================================================

cat("\n--- Estatísticas detalhadas por organismo ---\n")

stats_summary <- data.frame(
  Organismo = c("Homo sapiens", "Toxoplasma gondii"),
  Total_Genes = c(nrow(counts_human), nrow(counts_toxo)),
  Total_Reads = c(total_human, total_toxo),
  Media_Reads_por_Gene = c(mean(rowSums(counts_human)), mean(rowSums(counts_toxo))),
  Mediana_Reads_por_Gene = c(median(rowSums(counts_human)), median(rowSums(counts_toxo))),
  Min_Reads_por_Gene = c(min(rowSums(counts_human)), min(rowSums(counts_toxo))),
  Max_Reads_por_Gene = c(max(rowSums(counts_human)), max(rowSums(counts_toxo))),
  Percentual_Total = c((total_human / total_all) * 100, (total_toxo / total_all) * 100)
)

cat("\n")
print(stats_summary)

write.csv(stats_summary, 
          file = file.path(path_results, "organism_statistics.csv"), 
          row.names = FALSE)
cat("\nTabela salva: organism_statistics.csv\n")

# ============================================================
# 9. RESUMO FINAL
# ============================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("RESUMO DA AULA 07\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Pergunta biológica:\n")
cat("  'Qual é a proporção de reads do hospedeiro vs parasita?'\n\n")

cat("Resultados principais:\n")
cat("  Homo sapiens:\n")
cat("    - Total de reads:", total_human, "\n")
cat("    - Percentual:", round((total_human / total_all) * 100, 2), "%\n")
cat("    - Genes:", nrow(counts_human), "\n\n")

cat("  Toxoplasma gondii:\n")
cat("    - Total de reads:", total_toxo, "\n")
cat("    - Percentual:", round((total_toxo / total_all) * 100, 2), "%\n")
cat("    - Genes:", nrow(counts_toxo), "\n\n")

cat("Interpretação:\n")
if ((total_toxo / total_all) * 100 < 5) {
  cat("Baixa carga parasitária (< 5% dos reads)\n")
} else if ((total_toxo / total_all) * 100 < 20) {
  cat("  ✓ Carga parasitária moderada (5-20% dos reads)\n")
} else {
  cat("  ✓✓ Alta carga parasitária (> 20% dos reads)\n")
}

cat("\nArquivos gerados:\n")
cat("  - read_count_summary.csv\n")
cat("  - top_genes_human.csv\n")
cat("  - top_genes_toxo.csv\n")
cat("  - gene_read_count_comparison.csv\n")
cat("  - organism_statistics.csv\n")
cat("  - 01_read_count_stacked_bar.pdf\n")
cat("  - 02_read_count_proportion.pdf\n")
cat("  - 03_top_genes_human.pdf\n")
cat("  - 04_top_genes_toxo.pdf\n\n")

cat("Fim da aula 07.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")