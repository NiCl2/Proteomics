#!/usr/bin/Rscript
#mass-spec_stats.R
#perform statistical analysis of mass spectrometry protein quantification
#Nicholas Clifton 
#6th Jan 2021

library(biomaRt)
library(dplyr)
library(SummarizedExperiment)
library(modelr)
library(limma)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(splines)

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

#####################
### ID CONVERSION ###
#####################

# Uniprot gene name conversion function
uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri)
  dat
}

# perform conversions (500 at a time, due to server limits)
setsOf500 <- floor(length(ms_raw$First.protein.ID) / 500)
if(setsOf500 > 0) {
  print(paste("Converting group 1 of", setsOf500 + 1, "protein groups"))
  Gene.names.uniprot <- uniprot_mapping(ms_raw$First.protein.ID[1:500])
  for (i in 2:setsOf500) {
    print(paste("Converting group", i, "of", setsOf500 + 1, "protein groups"))
    Gene.names.uniprot <- rbind(Gene.names.uniprot, uniprot_mapping(ms_raw$First.protein.ID[((i-1)*500 + 1):(i*500)]))
  }
  print(paste("Converting group", setsOf500 + 1, "of", setsOf500 + 1, "protein groups"))
  Gene.names.uniprot <- rbind(Gene.names.uniprot, uniprot_mapping(ms_raw$First.protein.ID[(setsOf500*500 + 1):length(ms_raw$First.protein.ID)]))
} else {
  Gene.names.uniprot <- uniprot_mapping(ms_raw$First.protein.ID)
}

# Create new column selecting the first gene name from 1:many conversions
Gene.names.uniprot$First.gene.name <- sapply(Gene.names.uniprot$Gene.names, function(x) {
  unlist(strsplit(x, split = " "))[1]
}, USE.NAMES = F)

# filter out unrequested entries and remove duplicates
Gene.names.uniprot <- filter(Gene.names.uniprot, Entry %in% ms_raw$First.protein.ID & !duplicated(Entry))

# Select columns of interest and rename
Gene.names.uniprot <- dplyr::select(Gene.names.uniprot, Entry, Protein.names, Gene.names, First.gene.name)
colnames(Gene.names.uniprot) <- paste("uniprot", colnames(Gene.names.uniprot), sep = ".")

# get human gene IDs from biomaRt
getHumanEnsembl <- function(mouseSymbols) {
  
  ensembl = useMart("ensembl")
  mm = useDataset("mmusculus_gene_ensembl", mart=ensembl)
  hu = useDataset("hsapiens_gene_ensembl", mart=ensembl)
  conv = getLDS(attributes = c("mgi_symbol", "ensembl_gene_id"), 
                filters="mgi_symbol", values = mouseSymbols, mart = mm, 
                attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = hu,
                uniqueRows = T)
  conv
}
Human.ensembl <- getHumanEnsembl(Gene.names.uniprot$uniprot.First.gene.name)
Human.ensembl <- group_by(Human.ensembl, MGI.symbol) %>% filter(nchar(HGNC.symbol) == max(nchar(HGNC.symbol))) %>% ungroup
Human.ensembl <- filter(Human.ensembl, !duplicated(MGI.symbol))
Gene.names.uniprot <- left_join(Gene.names.uniprot, Human.ensembl, by = c("uniprot.First.gene.name" = "MGI.symbol"))

# join
ms_raw <- left_join(ms_raw, Gene.names.uniprot, by = c("First.protein.ID" = "uniprot.Entry"))

# remove Gene name (uniprot) duplicates, prioritising high intensity
ms_raw <- group_by(ms_raw, uniprot.First.gene.name) %>% filter(Intensity == max(Intensity)) %>% ungroup
ms_raw <- filter(ms_raw, !is.na(uniprot.First.gene.name))

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
rownames(ms_rse_filtered) <- rowData(ms_rse_filtered)$uniprot.First.gene.name

# replace -Inf with NaN
assays(ms_rse_filtered, withDimnames = F)$LFQ_intensity[assays(ms_rse_filtered, withDimnames = F)$LFQ_intensity == -Inf] <- NaN

# subtract the median per columns
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
  set.seed(1)
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

assays(ms_rse_filtered)$LFQ_intensity <- impute_data(as.data.frame(assays(ms_rse_filtered)$LFQ_intensity))

# boxplot: intensities of all channels after data preprocessing and normalization
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=0.8)
boxplot(assays(ms_rse_filtered)$LFQ_intensity, main="Boxplot normalized Intensities")

# define the design
sample_groups <- factor(paste(colData(ms_rse_filtered)$age, colData(ms_rse_filtered)$genotype, sep = "."))
MODEL <- model.matrix(~ 0 + sample_groups)
colnames(MODEL) <- levels(sample_groups)
contrast_spec <- makeContrasts(
  # E18.WT-E14.WT,
  # P7.WT-E18.WT,
  # P35.WT-P7.WT,
  P70.WT-P35.WT,
  # E14.Het-E14.WT,
  # E18.Het-E18.WT,
  # P7.Het-P7.WT,
  # P35.Het-P35.WT,
  # P70.Het-P70.WT,
  # (P70.WT|P70.Het)-!(P70.WT|P70.Het),
  levels = MODEL
)

# limma-type analysis
n <- dim(assays(ms_rse_filtered)$LFQ_intensity)[1]
fit <- lmFit(assays(ms_rse_filtered)$LFQ_intensity, MODEL)
fit_spec <- contrasts.fit(fit, contrast_spec)
# fit.eb <- eBayes(fit_spec)
fit.treat <- treat(fit_spec, lfc = log2(1.1))

# extract results (one comparison only)
get_output <- function(eb) {
  logFC = eb$coefficients[, 1]
  df.r = eb$df.residual
  df.0 = rep(eb$df.prior, n)
  s2.0 = rep(eb$s2.prior, n)
  s2 = (eb$sigma)^2
  s2.post = eb$s2.post
  t.ord = eb$coefficients[, 1]/eb$sigma/eb$stdev.unscaled[, 1]
  t.mod = eb$t[, 1]
  p.ord = 2*pt(-abs(t.ord), eb$df.residual)
  p.mod = eb$p.value[, 1]
  q.ord = qvalue(p.ord)$q
  q.mod = qvalue(p.mod)$q
  results.eb = data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb = results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

# res.eb <- get_output(fit.eb)
# res.eb$symbol <- rownames(res.eb)
res.treat <- get_output(fit.treat)
res.treat$symbol <- rownames(res.treat)

# volcano plot
lava_cutoff <- 0.05
ggplot(res.treat, aes(x = logFC, y = -log10(p.mod), colour = q.mod < lava_cutoff)) + 
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(data = res.treat[res.treat$q.mod < lava_cutoff,], aes(label=symbol), show.legend = FALSE) +
  theme_bw() + 
  geom_point(show.legend = FALSE)

# volcanoplot(fit.treat, highlight = sum(res.treat$q.mod < lava_cutoff), names = rownames(fit.treat), hl.col = "red")

ggsave("volcano_P70WT-P35WT.png", width = 25, height = 25, units = "cm")

# Spline analysis



# two-sample t-tests
# age1 <- "P7"
# age2 <- "P7"
# geno1 <- "WT"
# geno2 <- "Het"
# t_out <- apply(assays(ms_rse_filtered)$LFQ_intensity, 1, function(x) {
#   # each_t = tryCatch(
#     each_t = t.test(x[colData(ms_rse_filtered)$age == age1 & colData(ms_rse_filtered)$genotype == geno1], x[colData(ms_rse_filtered)$age == age2 & colData(ms_rse_filtered)$genotype == geno2])
#     # error = function(e) { NULL })
#   # if(is.null(each_t)) {
#   #   return(c(NaN, NaN))
#   # } else {
#     return(c(each_t$statistic, each_t$p.value))
#   # }
# })
# t_out <- data.frame(gene = rowData(ms_rse_filtered)$Gene.names, protein = rowData(ms_rse_filtered)$Protein.names, t = t_out[1,], p = t_out[2,])
# t_out$fdr <- p.adjust(t_out$p, method = "fdr")
# t_out$log10p <- -log10(t_out$p)
# View(arrange(t_out, p))
# plot(t_out$t, t_out$log10p)
