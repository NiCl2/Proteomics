#!/usr/bin/Rscript
#mass-spec_out.R
#process statistical output from mass spectrometry protein quantification
#Nicholas Clifton 
#5th Jan 2021

library(dplyr)

setwd("/Users/nclifton/Documents/Proteomics/Setd1a_Dev/")
ms_out <- read.csv("TimePointsStatsOut.csv", header = T, fill = T)

logP_cutoff <- -log10(0.05/nrow(ms_out))

getSigGenes <- function(timepoint1, timepoint2, genotype) {
  return(filter(ms_out, ms_out[, paste0("X.Log.Student.s.T.test.p.value.", timepoint1, "_", genotype, "_", timepoint2, "_", genotype)] > logP_cutoff)$Gene.names)
}

p35_p7 <- getSigGenes("P35", "P7", "WT")
p70_p35 <- getSigGenes("P70", "P35", "WT")


getHumanEnsembl <- function(mouseSymbols) {
  
  ensembl <- useMart("ensembl")
  mm <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
  hu <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  conv <- getLDS(attributes = c("mgi_symbol", "ensembl_gene_id"), 
                 filters="mgi_symbol", values = mouseSymbols, mart = mm, 
                 attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = hu,
                 uniqueRows = T)

}

p35_p7_hu <- getHumanEnsembl(p35_p7)



