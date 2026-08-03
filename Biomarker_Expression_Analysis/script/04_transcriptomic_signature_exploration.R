# ==========================================================
#
# RNA-Seq Analysis using Linux and R
#
# Module 4
# Transcriptomic Signature Analysis
#
# Description:
# Construct and explore transcriptomic signatures
# representing the developmental stages of
# Toxoplasma gondii using RNA-Seq data.
#
# ==========================================================



############################################################
##
## Biological Question
##
## Biological processes are rarely controlled
## by a single gene.
##
## Instead, they result from coordinated changes
## in the expression of multiple genes.
##
## Therefore, we will analyze transcriptomic
## signatures rather than isolated biomarkers.
##
############################################################



############################################################
##
## 1. Load packages
##
############################################################

library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(pheatmap)
library(RColorBrewer)



############################################################
##
## 2. Load RNA-Seq object
##
############################################################

load("/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/se_toxo_gene.RData")



############################################################
##
## 3. Explore the SummarizedExperiment object
##
############################################################

# Available assays
assayNames(se_gene)

# Number of genes and samples
dim(se_gene)

# Sample metadata
colData(se_gene)

# Gene annotation
rowData(se_gene)

# Experiment metadata
metadata(se_gene)



############################################################
##
## 4. Extract TPM matrix
##
############################################################

tpm <- assay(
  se_gene,
  "abundance"
)

head(tpm)



############################################################
##
## 5. Select infected samples
##
############################################################

infected_samples <- colData(se_gene)$sample_type == "infected"

tpm_inf <- tpm[ , infected_samples]

############################################################

# Why only infected samples?
#
# Parasite genes are expected to be expressed
# only in infected samples.
#
# Therefore, only infected samples are
# biologically informative for this analysis.
#
############################################################



############################################################
##
## 6. Explore the TPM matrix
##
############################################################

dim(tpm_inf)

head(tpm_inf)

summary(as.vector(tpm_inf))

range(tpm_inf)



############################################################
##
## 7. Import transcriptomic signatures
##
############################################################

brady_signature <- read.csv("/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/UP_Brady_FC10.csv",
                            stringsAsFactors = FALSE)



tachy_signature <- read.csv("/media/jean/OneDrive/Bioinformatica_na_pratica/Aulas_particulares/Bioinformatica_vanessa/UP_Tachy_FC5.csv",
                            stringsAsFactors = FALSE)

############################################################
##
## 8. Explore transcriptomic signatures
##
############################################################

head(brady_signature)

head(tachy_signature)

dim(brady_signature)

dim(tachy_signature)

colnames(brady_signature)

colnames(tachy_signature)

############################################################

# Before using these signatures,
# always inspect their structure.
#
# Questions:
#
# • How many genes compose each signature?
#
# • Which columns contain the Gene IDs?
#
# • Which columns contain Gene Symbols?
#
############################################################



############################################################
##
## 9. Extract Gene IDs
##
############################################################

brady_genes <- unique(brady_signature$Gene.ID)

tachy_genes <- unique(tachy_signature$Gene.ID)

############################################################
##
## 10. Explore transcriptomic signatures
##
############################################################

length(brady_genes)

length(tachy_genes)

head(brady_genes)

head(tachy_genes)



############################################################
##
## 11. Verify overlap with RNA-Seq dataset
##
############################################################

sum(brady_genes %in%
      rownames(tpm_inf))

sum(tachy_genes %in%
      rownames(tpm_inf))

############################################################

# Not every gene from the published signature
# is necessarily detected in our RNA-Seq dataset.
#
# Therefore, checking the overlap between
# signatures and the expression matrix is
# an important quality control step.
#
############################################################



############################################################
##
## 12. Build expression matrices
##
############################################################

brady_matrix <- tpm_inf[rownames(tpm_inf) %in% brady_genes,]

tachy_matrix <- tpm_inf[rownames(tpm_inf) %in% tachy_genes,]

############################################################
##
## 13. Explore expression matrices
##
############################################################

dim(brady_matrix)

dim(tachy_matrix)

head(brady_matrix)

head(tachy_matrix)



############################################################

# At this point we have two expression matrices.
#
# Bradyzoite signature
#
#       Genes × Samples
#
#
# Tachyzoite signature
#
#       Genes × Samples
#
# The next step will be to summarize the
# coordinated expression of each signature.
#
############################################################

############################################################
##
## 14. Save expression matrices
##
############################################################

save(brady_matrix, tachy_matrix,
     file = "Biomarker_Expression_Analysis/results/transcriptomic_signature_matrices.RData")
