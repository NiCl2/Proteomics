#!/usr/bin/Rscript
#mass-spec_stats.R
#perform statistical analysis of mass spectrometry protein quantification
#Nicholas Clifton 
#6th Jan 2021

library(biomaRt)
library(org.Mm.eg.db)
library(UniProt.ws)
library(dplyr)
library(SummarizedExperiment)
library(modelr)
library(limma)
# library(qvalue)
library(ggplot2)
library(ggrepel)
library(splines)
library(reshape2)
library(tidyr)
library(cowplot)

setwd("~/Documents/Proteomics/Setd1a_Dev/")

# read data and select relevant columns
ms_raw <- read.table("Canonical db unique pep quant/proteinGroups.txt", sep = "\t", header = T, fill = T) %>% 
  # rename certain colnames to allow grep selection of assay groups
  dplyr::rename(IBAQpeptides = iBAQ.peptides, SequenceCoveragePercent = Sequence.coverage....)

# read in meta data
ms_meta <- read.csv("sample_meta.csv", header = T) %>% arrange(sampleID)

# filter out contaminants and reverse database hits
ms_raw <- filter(ms_raw, !grepl("CON__", Protein.IDs) & !grepl("REV__", Protein.IDs))

# filter out proteins with 0 Intensity
ms_raw <- filter(ms_raw, Intensity != "0")

# get first protein from Majority.protein.IDs (first explains most of the data)
ms_raw$First.protein.ID <- sapply(ms_raw$Majority.protein.IDs, function(x) {
  unlist(strsplit(x, split = ";"))[1]
}, USE.NAMES = F)

# read in pre/post synaptic bias 
pre_post <- read.csv("Bayes2017_synaptic-enrichment.csv", header = T) 
pre_post$Enrichment <- apply(pre_post, 1, function(x) {
  if(x["PSD.enriched"] == "+") {"PSD"}
  else if (x["PSD.Depleted"] == "+") {"Syn"}
  else "NS"
})
pre_post <- dplyr::select(pre_post, ENSEBLE.ID, Enrichment)

#####################
### ID CONVERSION ###
#####################

# uniprot id mapping ftp file
Gene.IDs <- read.table("MOUSE_10090_idmapping_selected.tab", sep = "\t", na.strings = "") %>%
  dplyr::select(1:3, 19) %>% dplyr::rename(UNIPROTKB = V1, GENEID = V2, ENTREZID = V3, ENSEMBL = V19)
Gene.IDs$GENEID <- gsub("_MOUSE", "", Gene.IDs$GENEID)

# select first ENTREZ / ENSEMBL ID if many
Gene.IDs$ENSEMBL <- gsub(";.*", "", Gene.IDs$ENSEMBL)
Gene.IDs$ENTREZID <- gsub(";.*", "", Gene.IDs$ENTREZID)

# filter to match data
Gene.IDs <- filter(Gene.IDs, UNIPROTKB %in% ms_raw$First.protein.ID)

# get gene symbols from AnnotationDbi
Gene.symbols <- unlist(AnnotationDbi::mapIds(org.Mm.eg.db, 
                                      keys = Gene.IDs$ENSEMBL, 
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL", 
                                      multiVals = "first"))
Gene.symbols <- data.frame(ENSEMBL = names(Gene.symbols), SYMBOL = Gene.symbols)

# remove any 1:many conversions
Gene.symbols <- Gene.symbols %>% group_by(ENSEMBL) %>% filter(n() == 1) %>% ungroup

# join
Gene.IDs <- left_join(Gene.IDs, Gene.symbols, by = "ENSEMBL")

# replace missing symbols with Uniprot entry names
Gene.IDs$SYMBOL[is.na(Gene.IDs$SYMBOL)] <- Gene.IDs$GENEID[is.na(Gene.IDs$SYMBOL)]

# get human IDs from biomaRt
ensembl = useMart("ensembl")
mm = useDataset("mmusculus_gene_ensembl", mart=ensembl)
hu = useDataset("hsapiens_gene_ensembl", mart=ensembl)
conv = getLDS(attributes = c("ensembl_gene_id"), 
              filters="ensembl_gene_id", values = Gene.IDs$ENSEMBL, mart = mm, 
              attributesL = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), martL = hu)

# if duplicated, prioritise rows with HGNC symbols
conv <- conv %>% group_by(Gene.stable.ID) %>% filter(nchar(HGNC.symbol) == max(nchar(HGNC.symbol))) %>% ungroup

# join
Gene.IDs <- left_join(Gene.IDs, conv, by = c("ENSEMBL" = "Gene.stable.ID"))

# add pre/post synaptic bias
Gene.IDs <- left_join(Gene.IDs, pre_post, by = c("ENSEMBL" = "ENSEBLE.ID"))

# get unique protein entries
Gene.IDs_unique <- Gene.IDs[!duplicated(Gene.IDs$UNIPROTKB), ]
Gene.IDs_unique <- Gene.IDs_unique[!duplicated(Gene.IDs$SYMBOL), ]

# join to dataset
ms_raw <- left_join(ms_raw, Gene.IDs_unique, by = c("First.protein.ID" = "UNIPROTKB"))

# remove Gene ID (uniprot) duplicates, prioritising high intensity
ms_raw <- group_by(ms_raw, First.protein.ID) %>% filter(Intensity == max(Intensity)) %>% ungroup
ms_raw <- filter(ms_raw, !is.na(SYMBOL))
ms_raw <- filter(ms_raw, !duplicated(SYMBOL))

#############################
### SUMMARIZED EXPERIMENT ###
#############################

# check sample order matching
stopifnot(gsub(pattern = "LFQ.intensity.", replacement = "", colnames(ms_raw)[grepl("LFQ.intensity.", colnames(ms_raw))]) == ms_meta$sampleID)

# create ranged summarized experiment object
ms_rse <- SummarizedExperiment(
  assays = list(
    # peptides = ms_raw[, grepl("Peptides.", colnames(ms_raw))],
    # unique_peptides = ms_raw[, grepl("Unique.peptides.", colnames(ms_raw))],
    # unique_peptides = ms_raw[, grepl("Unique.peptides.", colnames(ms_raw))],
    # identification_type = ms_raw[, grepl("Identification.type.", colnames(ms_raw))],
    # sequence_coverage = ms_raw[, grepl("Sequence.coverage.", colnames(ms_raw))],
    # intensity = ms_raw[, grepl("Intensity.", colnames(ms_raw))],
    # iBAQ = ms_raw[, grepl("iBAQ.", colnames(ms_raw))],
    LFQ_intensity = ms_raw[, grepl("LFQ.intensity.", colnames(ms_raw))]#,
    # MSMS_count = ms_raw[, grepl("MS.MS.count.", colnames(ms_raw))]
  ),
  # row data is selected by those colnames without "_". Not future-proof
  rowData = ms_raw[, !grepl("_", colnames(ms_raw))],
  colData = ms_meta
)

# reset columns names
colnames(ms_rse) <- ms_rse$sampleID

#########################################
### DATA FILTERING AND TRANSFORMATION ###
#########################################

# log transformation
assays(ms_rse, withDimnames = F)$LFQ_intensity <- log2(assays(ms_rse, withDimnames = F)$LFQ_intensity)

# filter rows based on valid values in at least 4 samples from at least 1 group
filtering_index <- cbind(ms_meta, t(assays(ms_rse)$LFQ_intensity) > 0) %>% 
  group_by(age, genotype) %>% 
  summarise(across(`1`:tail(names(.), 1), sum)) %>%
  ungroup %>%
  dplyr::select(!c(age, genotype)) %>%
  t()

ms_rse_filtered <- ms_rse[rowSums(filtering_index >= 4) >= 1, ]

# set row names
rownames(ms_rse_filtered) <- rowData(ms_rse_filtered)$SYMBOL

# replace -Inf with NaN
assays(ms_rse_filtered, withDimnames = F)$LFQ_intensity[assays(ms_rse_filtered, withDimnames = F)$LFQ_intensity == -Inf] <- NaN

# subtract the median per column
assays(ms_rse_filtered)$LFQ_intensity <- apply(assays(ms_rse_filtered)$LFQ_intensity, 2, function(x) { x - median(x, na.rm = T)} )

# keep record of missing values before imputation
assays(ms_rse_filtered)$LFQ_missing <- sapply(colnames(ms_rse_filtered), function(x) !is.finite(assays(ms_rse_filtered)$LFQ_intensity[, x]))

# Impute missing values using normal distribution
impute_data = function(df, width = 0.3, downshift = 1.8) {
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  
  df.names = colnames(df)
  
  # Create new columns indicating whether the values are missing
  # df[paste(colnames(df), "missing", sep = "_")] = sapply(colnames(df), function(x) !is.finite(df[, x]))
  
  # Imputation
  df[df.names] = sapply(df.names,
                          function(x) {
                            temp = df[[x]]
                            temp[!is.finite(temp)] = NA
                            
                            temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
                            temp.mean = mean(temp, na.rm = TRUE) - 
                              downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
                            
                            n.missing = sum(is.na(temp))
                            temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                            return(temp)
                          })
  return(df)
}

set.seed(1)
assays(ms_rse_filtered)$LFQ_intensity <- impute_data(as.data.frame(assays(ms_rse_filtered)$LFQ_intensity))

# boxplot: intensities of all channels after data preprocessing and normalization
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=0.8)
boxplot(assays(ms_rse_filtered)$LFQ_intensity, main="Boxplot normalized Intensities")
plotMDS(assays(ms_rse_filtered)$LFQ_intensity)

# define the design
sample_groups <- factor(paste(colData(ms_rse_filtered)$age, colData(ms_rse_filtered)$genotype, sep = "."))
MODEL <- model.matrix(~ 0 + sample_groups)
colnames(MODEL) <- levels(sample_groups)
contrast_spec <- makeContrasts(
  # E18.WT-E14.WT,
  # P7.WT-E18.WT,
  # P35.WT-P7.WT,
  # P70.WT-P35.WT,
  E14.Het-E14.WT,
  E18.Het-E18.WT,
  P7.Het-P7.WT,
  P35.Het-P35.WT,
  P70.Het-P70.WT,
  # Diff.E18_E14 = (E18.Het-E14.Het)-(E18.WT-E14.WT),
  # Diff.P7_E18 = (P7.Het-E18.Het)-(P7.WT-E18.WT),
  # Diff.P35_P7 = (P35.Het-P7.Het)-(P35.WT-P7.WT),
  # Diff.P70_E35 = (P70.Het-P35.Het)-(P70.WT-P35.WT),
  # embryonic.Het-embryonic.WT,
  # postnatal.Het-postnatal.WT,
  levels = MODEL
)

# limma-type analysis
limma_analysis <- function(x, contrast_spec) {
  fit = lmFit(x, MODEL)
  fit_spec = contrasts.fit(fit, contrast_spec)
  # EITHER
  fit.eb = eBayes(fit_spec) # works with multiple simultaneous contrasts
  res.eb = topTable(fit.eb, number = Inf)
  res.eb$SYMBOL = rownames(res.eb)
  
  #OR
  # fit.eb = eb(fit_spec, lfc = log2(1.1)) # ensure minimum logfc. Only works with single contrasts
  # res.eb = topTreat(fit.eb, number = Inf, adjust.method = "BH")
  # res.eb$SYMBOL = rownames(res.eb)
  return(res.eb)
}

res.eb <- limma_analysis(x = assays(ms_rse_filtered)$LFQ_intensity, contrast_spec = contrast_spec)

# write to file
res.eb <- left_join(res.eb, Gene.IDs, by = "SYMBOL")
write.table(res.eb, "DEP_Setd1aHet-WT_Embryonic.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# write gene property file
res.gp <- res.eb %>% dplyr::select(Gene.stable.ID.1, "F") %>% filter(complete.cases(.)) %>% filter(!duplicated(Gene.stable.ID.1))
write.table(res.gp, "DEP_Setd1aHet-WT_AllTimepoints_gp.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# classic models
# colData(ms_rse_filtered)$genotype <- factor(colData(ms_rse_filtered)$genotype, levels = c("WT", "Het"))
# MODEL <- model.matrix(~synaptosome.conc + age + genotype, 
#                       data = colData(ms_rse_filtered))
# fit <- lmFit(assays(ms_rse_filtered)$LFQ_intensity, MODEL)
# fit.eb <- eBayes(fit)
# res.eb <- topTable(fit.eb, coef="genotypeHet", number = Inf)

# extract results (one comparison only)
# n <- dim(assays(ms_rse_filtered)$LFQ_intensity)[1]
# get_output <- function(eb) {
#   logFC = eb$coefficients[, 1]
#   df.r = eb$df.residual
#   df.0 = rep(eb$df.prior, n)
#   s2.0 = rep(eb$s2.prior, n)
#   s2 = (eb$sigma)^2
#   s2.post = eb$s2.post
#   t.ord = eb$coefficients[, 1]/eb$sigma/eb$stdev.unscaled[, 1]
#   t.mod = eb$t[, 1]
#   p.ord = 2*pt(-abs(t.ord), eb$df.residual)
#   p.mod = eb$p.value[, 1]
#   q.ord = qvalue(p.ord)$q
#   q.mod = qvalue(p.mod)$q
#   results.eb = data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
#   results.eb = results.eb[order(results.eb$p.mod), ]
#   return(results.eb)
# }

# res.eb <- get_output(fit.eb)
# res.eb$symbol <- rownames(res.eb)

# volcano plot
volcano_plot <- function(res.eb) {
  lava_cutoff <- 0.05
  ggplot(res.eb, aes(x = logFC, y = -log10(P.Value), colour = adj.P.Val < lava_cutoff)) + 
    scale_color_manual(values = c("black", "red")) +
    geom_text_repel(data = res.eb[res.eb$adj.P.Val < lava_cutoff,], aes(label=SYMBOL), show.legend = FALSE) +
    theme_bw() + 
    geom_point(show.legend = FALSE)
}

volcano_plot(res.eb)

# alternative volcano plot method
# volcanoplot(fit.eb, highlight = sum(res.eb$q.mod < lava_cutoff), names = rownames(fit.eb), hl.col = "red")

ggsave("volcano_P70WT-P35WT.png", width = 25, height = 25, units = "cm")

# Loop over contrasts

contrast_list <- c(
  # "E18.WT-E14.WT",
  # "P7.WT-E18.WT",
  # "P35.WT-P7.WT",
  # "P70.WT-P35.WT"
  "E14.Het-E14.WT",
  "E18.Het-E18.WT",
  "P7.Het-P7.WT",
  "P35.Het-P35.WT",
  "P70.Het-P70.WT"
  # "(E18.Het-E14.Het)-(E18.WT-E14.WT)",
  # "(P7.Het-E18.Het)-(P7.WT-E18.WT)",
  # "(P35.Het-P7.Het)-(P35.WT-P7.WT)",
  # "(P70.Het-P35.Het)-(P70.WT-P35.WT)"
)

# prime gene set output
genesets <- data.frame(GeneSetID = "All_synaptosome", EnsemblID = rowData(ms_rse_filtered)$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))

for (i in 1:length(contrast_list)) {
  # get each contrast specification
  each_contrast = contrast_list[i]
  print(each_contrast)
  contrast_spec = makeContrasts(
    each_contrast,
    levels = MODEL
  )
  
  # perform limma differential expression analysis
  each.res = limma_analysis(x = assays(ms_rse_filtered)$LFQ_intensity, contrast_spec = contrast_spec)
  
  # append meta data (Ensembl ID and pre/post synapse bias)
  each.res = left_join(each.res, Gene.IDs_unique, by = "SYMBOL")
  each.res <- filter(each.res, !duplicated(SYMBOL))
  
  write.table(each.res, paste0("DE_", each_contrast, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  volcano_plot(each.res)
  ggsave(paste0("volcano_", each_contrast, ".png"), units = "cm", width = 25, height = 25)
  
  # collate gene sets
  # if(any(each.res$adj.P.Val < 0.05, na.rm = T)){
  #   gs.all = data.frame(GeneSetID = paste0("DE_", each_contrast, "_All"), EnsemblID = filter(each.res, adj.P.Val < 0.05)$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.all)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC > 0, na.rm = T)){
  #   gs.up = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Up"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC > 0)$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.up)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC < 0, na.rm = T)){
  #   gs.down = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Down"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC < 0)$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.down)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC > 0 & each.res$Enrichment == "PSD", na.rm = T)){
  #   gs.up_psd = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Up_PSD-enriched"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC > 0 & Enrichment == "PSD")$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.up_psd)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC > 0 & each.res$Enrichment == "Syn", na.rm = T)){
  #   gs.up_syn = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Up_Syn-enriched"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC > 0 & Enrichment == "Syn")$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.up_syn)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC > 0 & each.res$Enrichment == "NS", na.rm = T)){
  #   gs.up_ns = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Up_not-enriched"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC > 0 & Enrichment == "NS")$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.up_ns)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC < 0 & each.res$Enrichment == "PSD", na.rm = T)){
  #   gs.down_psd = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Down_PSD-enriched"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC < 0 & Enrichment == "PSD")$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.down_psd)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC < 0 & each.res$Enrichment == "Syn", na.rm = T)){
  #   gs.down_syn = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Down_Syn-enriched"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC < 0 & Enrichment == "Syn")$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.down_syn)
  # }
  # if(any(each.res$adj.P.Val < 0.05 & each.res$logFC < 0 & each.res$Enrichment == "NS", na.rm = T)){
  #   gs.down_ns = data.frame(GeneSetID = paste0("DE_", each_contrast, "_Down_not-enriched"), EnsemblID = filter(each.res, adj.P.Val < 0.05 & logFC < 0 & Enrichment == "NS")$Gene.stable.ID.1) %>% filter(!is.na(EnsemblID))
  #   genesets <- rbind(genesets, gs.down_ns)
  # }
  
}

# add geneset for constituent synaptic proteins (for age contrasts)
genesets <- rbind(genesets, data.frame(GeneSetID = "synaptosome_no-changes", EnsemblID = setdiff(genesets$EnsemblID[genesets$GeneSetID == "All_synaptosome"], genesets$EnsemblID[genesets$GeneSetID != "All_synaptosome"])))

# write gene sets for downstream analyses
write.table(genesets, "DEP_Het-WT_EachTimepoint_up-down_enrichment_gs.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# write ENTREZ ID version
genesets_entrez <- genesets %>% left_join(Gene.IDs_unique, by = c("EnsemblID" = "Gene.stable.ID.1")) %>% dplyr::select(GeneSetID, NCBI.gene..formerly.Entrezgene..ID) %>%
  group_by(GeneSetID) %>% filter(!duplicated(NCBI.gene..formerly.Entrezgene..ID))
write.table(genesets_entrez, "DEP_Dev-steps_up-down_enrichment_genesets_entrez.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# Spline analysis

# setup a basis for a natural regression spline
spline.nat_reg <- ns(colData(ms_rse_filtered)$age_days, df = 4) # knots ~ df
genotype <- factor(colData(ms_rse_filtered)$genotype, levels = c("WT", "Het"))
spline.MODEL <- model.matrix(~genotype*spline.nat_reg)
# get genes with different time trends for het vs wt
spline.fit <- lmFit(assays(ms_rse_filtered)$LFQ_intensity, spline.MODEL)
spline.eb <- eBayes(spline.fit)
res.spline <- topTable(spline.eb, coef = 7:10, number = Inf)

# Developmental plots

# plot wt vs het developmental protein expression of selected genes
plot_select <- rownames(res.spline)[1:10]
plot_select <- rownames(res.eb)[1:10]
plot_select <- c("SYT2")
lfq.struct <- as.data.frame(t(assays(ms_rse_filtered)$LFQ_intensity[plot_select, ]))
lfq.struct$stage <- colData(ms_rse_filtered)$age
lfq.struct$genotype <- colData(ms_rse_filtered)$genotype
lfq.struct <- melt(lfq.struct, id = c("stage", "genotype"))
lfq.struct$stage <- factor(lfq.struct$stage, levels = c("E14", "E18", "P7", "P35", "P70"))
spline_plots <- list()
for (i in 1:length(plot_select)) {
  spline_plots[[i]] <- ggplot(filter(lfq.struct, variable == plot_select[i]), aes(x = stage, y = value, group = genotype, colour = genotype)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "loess", alpha = 0.15, aes(colour = genotype, fill = genotype), level = 0.95) + 
    scale_color_manual(values=c("purple3", "steelblue3")) +
    ggtitle(plot_select[i]) +
    theme_cowplot() +
    ylab("Normalized LFQ intensity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), 
          plot.title = element_text(size = 10),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 10))
}
plot_grid(plotlist = spline_plots, nrow = 2)
# ggsave("DevPlots_DE_Het-WT_AllTimepoints_top10.png", units = "cm", width = 30, height = 10)


## Cluster by developmental change

# filter for wildtype samples
LFQ.mean <- as.data.frame(assays(ms_rse_filtered[, colData(ms_rse_filtered)$genotype == "WT"])$LFQ_intensity)

# scale data
LFQ.mean <- as.data.frame(t(scale(t(LFQ.mean))))
LFQ.mean <- as.data.frame(apply(LFQ.mean, 1, function(x){ (x-mean(x))/sd(x-mean(x))}))

# average across stages
LFQ.mean$stage <- factor(colData(ms_rse_filtered[, colData(ms_rse_filtered)$genotype == "WT"])$age, 
                             levels = c("E18", "E14", "P7", "P35", "P70"))
LFQ.mean <- melt(LFQ.mean, id = "stage")
LFQ.mean <- dcast(LFQ.mean, variable ~ stage, mean)
rownames(LFQ.mean) <- LFQ.mean$variable
LFQ.mean$variable <- NULL

# scree plot
set.seed(1)
wss <- sapply(1:40, function(i) kmeans(LFQ.mean, centers = i, nstart = 20, iter.max = 20)$tot.withinss)
plot(1:40, wss, type = "b", xlab = "Number of Clusters", ylab = "Total within sum of squares")
text(1:40, wss, pos = 3, cex = 0.5)

# kmeans
set.seed(1)
LFQ.km <- kmeans(LFQ.mean, centers = 4, nstart = 20, iter.max = 20)
table(LFQ.km$cluster)

lineplot_LFQ <- function(d, c) {
  ggplot(data = d, aes(x = stage, y = LFQ, group = GeneSymbol)) +
    stat_smooth(geom = "line",
                size = 1, 
                method = "loess", 
                se = F, 
                span = 2, 
                color = "coral",
                alpha = 0.3) + 
    # scale_color_manual(values=c("#0F7FFE", "#FB0106")) +
    theme_cowplot() +
    ggtitle(paste("Cluster", c)) +
    theme(axis.text.x = element_text(size = 20), 
          plot.title = element_text(size = 20),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 20)) +
    xlab(NULL) +
    ylab(expression("Normalised LFQ intensity")) +
    ylim(-2.0, 2.0)
}

cluster_plots <- list()
for (c in 1: max(LFQ.km$cluster)) {
  LFQ.mean_cluster = LFQ.mean %>% mutate(GeneSymbol = rownames(LFQ.mean)) %>%
    filter(LFQ.km$cluster == c) %>%
    gather(stage, LFQ, 1:5, factor_key = T) 
  cluster_plots[[c]] <- lineplot_LFQ(LFQ.mean_cluster, c)
}

plot_grid(plotlist = cluster_plots, nrow = 2)
ggsave("clusters_wt_1-4_norm.png", units = "cm", width = 60, height = 30)

# create cluster genesets
gs_clusters <- data.frame(GeneSetID = paste0("cluster_", LFQ.km$cluster), 
                          GeneSymbol = rownames(LFQ.mean)) %>%
  left_join(Gene.IDs_unique, by = c("GeneSymbol" = "SYMBOL"))

gs_clusters.all <- gs_clusters %>%
  dplyr::select(GeneSetID, Gene.stable.ID.1) %>%
  filter(!is.na(Gene.stable.ID.1)) %>%
  arrange(GeneSetID)

# add pre/post synapse division
gs_clusters.pre_post <- gs_clusters
gs_clusters.pre_post$GeneSetID <- paste0(gs_clusters.pre_post$GeneSetID, "_", gs_clusters.pre_post$Enrichment)
gs_clusters.pre_post <- filter(gs_clusters.pre_post, !is.na(Enrichment)) %>% dplyr::select(GeneSetID, Gene.stable.ID.1) %>%
  filter(!is.na(Gene.stable.ID.1)) %>% 
  arrange(GeneSetID)

# add background
all_synaptosome <- data.frame(GeneSetID = "All_synaptosome", Gene.stable.ID.1 = rowData(ms_rse_filtered)$Gene.stable.ID.1) %>% filter(!is.na(Gene.stable.ID.1))
gs_clusters.magma <- rbind(gs_clusters.all, gs_clusters.pre_post, all_synaptosome)

# write ensembl and entrez versions
write.table(gs_clusters.magma, "synaptosome-protein_dev-clusters_wt_1-6.txt", sep = "\t", row.names = F, col.names = F, quote = F)

gs_clusters.magma_entrez <- gs_clusters.magma %>% left_join(Gene.IDs_unique, by = "Gene.stable.ID.1") %>% dplyr::select(GeneSetID, NCBI.gene..formerly.Entrezgene..ID) %>%
  group_by(GeneSetID) %>% filter(!duplicated(NCBI.gene..formerly.Entrezgene..ID))
write.table(gs_clusters.magma_entrez, "synaptosome-protein_dev-clusters_wt_1-6_entrez.txt", sep = "\t", row.names = F, col.names = F, quote = F)

