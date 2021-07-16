# SimpliFi protein analysis

library(biomaRt)

setwd("~/Documents/Proteomics/Setd1a_Dev")

Q_ALPHA <- 0.05

read_data <- function(csv_file) {
  read.csv(csv_file) %>% 
    filter(p.value != "NCO", q.value < Q_ALPHA) %>% 
    dplyr::select(Accession.number, UniProt.description, p.value, q.value, fold.change) %>%
    mutate(First.protein.ID = gsub(";.*", "", Accession.number))
}

E14 <- read_data("SimpliFi E14_WT vs. E14_Het.csv")
E18 <- read_data("SimpliFi E18_WT vs. E18_Het.csv")
P7 <- read_data("SimpliFi P7_WT vs. P7_Het.csv")
P35 <- read_data("SimpliFi P35_WT vs. P35_Het.csv")
P70 <- read_data("SimpliFi P70_WT vs. P70_Het.csv")

# all_syn <- read.csv("SimpliFi E14_WT vs. E14_Het.csv") %>% mutate(GeneSetID = "all_synaptosome") %>% dplyr::select(GeneSetID, UniProt.description)
all_syn <- read.table("all_synaptosome_mm_ensembl.txt", sep = "\t", col.names = c("GeneSetID", "ENSEMBL"))

E14_up <- E14 %>% filter(fold.change > 1) %>% mutate(GeneSetID = "Het-WT_E14_up") %>% dplyr::select(GeneSetID, First.protein.ID)
E14_down <- E14 %>% filter(fold.change < 1) %>% mutate(GeneSetID = "Het-WT_E14_down") %>% dplyr::select(GeneSetID, First.protein.ID)
E18_up <- E18 %>% filter(fold.change > 1) %>% mutate(GeneSetID = "Het-WT_E18_up") %>% dplyr::select(GeneSetID, First.protein.ID)
E18_down <- E18 %>% filter(fold.change < 1) %>% mutate(GeneSetID = "Het-WT_E18_down") %>% dplyr::select(GeneSetID, First.protein.ID)
P7_up <- P7 %>% filter(fold.change > 1) %>% mutate(GeneSetID = "Het-WT_P7_up") %>% dplyr::select(GeneSetID, First.protein.ID)
P7_down <- P7 %>% filter(fold.change < 1) %>% mutate(GeneSetID = "Het-WT_P7_down") %>% dplyr::select(GeneSetID, First.protein.ID)
P35_up <- P35 %>% filter(fold.change > 1) %>% mutate(GeneSetID = "Het-WT_P35_up") %>% dplyr::select(GeneSetID, First.protein.ID)
P35_down <- P35 %>% filter(fold.change < 1) %>% mutate(GeneSetID = "Het-WT_P35_down") %>% dplyr::select(GeneSetID, First.protein.ID)
P70_up <- P70 %>% filter(fold.change > 1) %>% mutate(GeneSetID = "Het-WT_P70_up") %>% dplyr::select(GeneSetID, First.protein.ID)
P70_down <- P70 %>% filter(fold.change < 1) %>% mutate(GeneSetID = "Het-WT_P70_down") %>% dplyr::select(GeneSetID, First.protein.ID)

genesets <- rbind(E14_up, E14_down, E18_up, E18_down, P7_up, P7_down, P35_up, P35_down, P70_up, P70_down)

# uniprot id mapping ftp file
Gene.IDs <- read.table("MOUSE_10090_idmapping_selected.tab", sep = "\t", na.strings = "") %>%
  dplyr::select(1:3, 19) %>% dplyr::rename(UNIPROTKB = V1, GENEID = V2, ENTREZID = V3, ENSEMBL = V19)
Gene.IDs$GENEID <- gsub("_MOUSE", "", Gene.IDs$GENEID)
# select first ENTREZ / ENSEMBL ID if many
Gene.IDs$ENSEMBL <- gsub(";.*", "", Gene.IDs$ENSEMBL)
Gene.IDs$ENTREZID <- gsub(";.*", "", Gene.IDs$ENTREZID)
# filter to match data
Gene.IDs <- filter(Gene.IDs, UNIPROTKB %in% genesets$First.protein.ID)

genesets <- left_join(genesets, Gene.IDs, by = c("First.protein.ID" = "UNIPROTKB")) %>%
  dplyr::select(GeneSetID, ENSEMBL) %>%
  filter(!is.na(ENSEMBL)) %>%
  group_by(GeneSetID) %>%
  filter(!duplicated(ENSEMBL)) %>% ungroup

genesets <- rbind(genesets, all_syn)

# write mouse IDs to file
write.table(genesets, "SimpliFi_DEGs_EachTimepoint_up-down_mm.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# convert to human IDs
ensembl = useMart("ensembl")
mm = useDataset("mmusculus_gene_ensembl", mart=ensembl)
hu = useDataset("hsapiens_gene_ensembl", mart=ensembl)
conv = getLDS(attributes = c("ensembl_gene_id"), 
              filters="ensembl_gene_id", values = genesets$ENSEMBL, mart = mm, 
              attributesL = c("ensembl_gene_id"), martL = hu)

genesets_hu <- left_join(genesets, conv, by = c("ENSEMBL" = "Gene.stable.ID"))
genesets_hu <- filter(genesets_hu, !is.na(Gene.stable.ID.1)) %>% dplyr::select(GeneSetID, Gene.stable.ID.1)

write.table(genesets_hu, "SimpliFi_DEGs_EachTimepoint_up-down_hu.txt", sep = "\t", row.names = F, col.names = F, quote = F)
