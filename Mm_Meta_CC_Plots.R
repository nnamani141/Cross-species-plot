library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Set file path
fp <- "C:/Users/dnnamani/OneDrive - University of Tennessee/UTHSC/Personal_Ideas/"

# get to the data folder
dfp <- paste0(fp, "Diabetic_Liver_Data/Meta_analysis/") 


# 1. Collect the reproducible-core gene vectors you already saved per tier.
Obesity        <- read.csv(paste0(dfp, "Tier1_recurrence_core.csv"))[[1]]
Female_Obesity <- read.csv(paste0(dfp, "Tier1_Female_recurrence_core.csv"))[[1]]
STZ             <- read.csv(paste0(dfp, "Tier2_recurrence_core.csv"))[[1]]
DiabetesOnObs  <- read.csv(paste0(dfp, "Tier3a_recurrence_core.csv"))[[1]]
DietOnDiabetes <- read.csv(paste0(dfp, "Tier3b_recurrence_core.csv"))[[1]]

#    Use ENTREZ IDs; convert each with bitr
to_entrez <- function(sym) {
  bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
}

tier_gene_lists <- list(
  "M_Tier1"    = to_entrez(Obesity),
  "F_Tier1"    = to_entrez(Female_Obesity),
  "M_Tier2"   = to_entrez(STZ),
  "M_Tier3a" = to_entrez(DiabetesOnObs),
  "M_Tier3b"  = to_entrez(DietOnDiabetes)
)

# Get the universes rom the full tables.
t1_all  <- read.xlsx(file.path(dfp, "Obesity_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t1_F_all  <- read.xlsx(file.path(dfp, "Obesity_Females_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t2_all  <- read.xlsx(file.path(dfp, "Diabetes_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t3a_all <- read.xlsx(file.path(dfp, "DiabOnObesity_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t3b_all <- read.xlsx(file.path(dfp, "ObsOnDiabetes_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")

shared_universe_mm <- Reduce(intersect, list(t1_all$gene, t1_F_all$gene, t2_all$gene,
                                             t3a_all$gene, t3b_all$gene))

# 2. One call runs enrichment on every tier and keeps them aligned.
cc_go <- compareCluster(
  geneCluster   = tier_gene_lists,
  fun           = "enrichGO",
  universe      = to_entrez(shared_universe_mm),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cc_go1 <- compareCluster(
  geneCluster   = tier_gene_lists,
  universe      = to_entrez(shared_universe_mm),
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cc_go2 <- compareCluster(
  geneCluster   = tier_gene_lists,
  universe      = to_entrez(shared_universe_mm),
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)



# 3. The comparative dotplot: rows = pathways, columns = tiers.
dp <- dotplot(cc_go, showCategory = 8) 
# Adjust font sizes across different text elements
dp + theme(
  # Shrink pathway terms on the Y-axis (e.g., set to size 8 or 9)
  axis.text.y = element_text(size = 10),
)

dotplot(cc_go1, showCategory = 8) 
dpc <- dotplot(cc_go2, showCategory = 8) 
# Adjust font sizes across different text elements
dpc + theme(
  # Shrink pathway terms on the Y-axis (e.g., set to size 8 or 9)
  axis.text.y = element_text(size = 10),
)
# get df for the results
CC_BP_results <- as.data.frame(cc_go)
CC_CC_results <- as.data.frame(cc_go1)
CC_MF_results <- as.data.frame(cc_go2)


# AUTOMATIC EXPORT TO EXCEL (SUPPLEMENTARY TABLE S1)
wb <- createWorkbook()
# CC enrichment
addWorksheet(wb, "CC_BP")
addWorksheet(wb, "CC_CC")
addWorksheet(wb, "CC_MF")

# writedata for GSEA enrichment
writeData(wb, "CC_BP",  CC_BP_results)
writeData(wb, "CC_CC",    CC_CC_results)
writeData(wb, "CC_MF",  CC_MF_results)

# Save the workbook to your directory
output_file <- file.path(dfp, "CompareCluster_Mouse_Meta_DEGs.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

cat(paste0("\nFile successfully exported to: ", output_file, "\n"))


# cnetplot
cnetplot(cc_go, showCategory = 8)
cnetplot(cc_go1, showCategory = 8)
cnetplot(cc_go2, showCategory = 8)


#dotplot(cc_go, showCategory = 8) +
 # scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    #    axis.text.y = element_text(size = 8))



