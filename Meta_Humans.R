# ==============================================================================
# HUMAN SIDE: four stages
#   Steatosis        = MASL (GSE281797) + NAFL (GSE126848)   -> 2-cohort meta
#   Steatohepatitis  = MASH (GSE281797) + NASH (GSE126848)   -> 2-cohort meta
#   Obese            = GSE126848 only                        -> single dataset
#   Fibrosis         = GSE246221 F1/F2/F3-4                  -> single dataset,
#                                                               severity gradient
#   Progression      = MASH_vs_MASL + NASH_vs_NAFL           -> 2-cohort meta
# With k = 2, I2/Cochran's Q are uninformative (df = 1). Reproducibility is
# defined by cross-cohort directional concordance instead.
# ==============================================================================

library(dplyr)
library(openxlsx)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggridges)

fp  <- "C:/Users/dnnamani/OneDrive - University of Tennessee/UTHSC/Personal_Ideas/"
dfp <- paste0(fp, "Diabetic_Liver_Data/Meta_analysis/")

# Stage cores for compareCluster (human side)
# 1. Collect the reproducible-core gene vectors you already saved per tier
# Male
MASL        <- read.csv(paste0(dfp, "Tier_MASL_Human_recurrence_core.csv"))[[1]]
MASH <- read.csv(paste0(dfp, "Tier_MASH_Human_recurrence_core.csv"))[[1]]
Progression             <- read.csv(paste0(dfp, "Tier_Progression_Human_recurrence_core.csv"))[[1]]
Fibrosis  <- read.csv(paste0(dfp, "Tier_Fibrosis_Human_recurrence_core.csv"))[[1]]


# Female
Female_MASL <- read.csv(paste0(dfp, "Tier_Female_MASL_Human_recurrence_core.csv"))[[1]]
Female_MASH <- read.csv(paste0(dfp, "Tier_Female_MASH_Human_recurrence_core.csv"))[[1]]
Female_Progression <- read.csv(paste0(dfp, "Tier_Female_Progression_Human_recurrence_core.csv"))[[1]]
Female_Fibrosis <- read.csv(paste0(dfp, "Tier_Female_Fibrosis_Human_recurrence_core.csv"))[[1]]

#    Use ENTREZ IDs; convert each with bitr
to_entrez_hs <- function(sym) {
  bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}

# Combined
human_cores <- list(
  "M_MASL"           = to_entrez_hs(MASL),
  "M_MASH"           = to_entrez_hs(MASH),
  "M_Progresssion"   = to_entrez_hs(Progression),
  "M_Fibrosis"       = to_entrez_hs(Fibrosis),
  "F_MASL"         = to_entrez_hs(Female_MASL),
  "F_MASH"         = to_entrez_hs(Female_MASH),
  "F_Progresssion" = to_entrez_hs(Female_Progression),
  "F_Fibrosis"      = to_entrez_hs(Female_Fibrosis)
)

# Males
Male_human_cores <- list(
  "Male_MASL"           = to_entrez_hs(MASL),
  "Male_MASH"           = to_entrez_hs(MASH),
  "Male_Progresssion"   = to_entrez_hs(Progression),
  "Male_Fibrosis"       = to_entrez_hs(Fibrosis)
)

# Females
Female_human_cores <- list(
  "Female_MASL"         = to_entrez_hs(Female_MASL),
  "Female_MASH"         = to_entrez_hs(Female_MASH),
  "Female_Progresssion" = to_entrez_hs(Female_Progression),
  "FemaleFibrosis"      = to_entrez_hs(Female_Fibrosis)
)

# Get the universes rom the full tables.
t1_all  <- read.xlsx(file.path(dfp, "Human_Steatosis_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t2_all  <- read.xlsx(file.path(dfp, "Human_Steatohepatitis_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t3_all <- read.xlsx(file.path(dfp, "Human_Progression_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t4_all <- read.xlsx(file.path(dfp, "Human_FibrosisvF0_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t5_all  <- read.xlsx(file.path(dfp, "Human_Female_Steatosis_Highly_Reproducible_Meta_DEGs.xlsx"),
                       sheet = "All_Meta_Analysis_Results")
t6_all  <- read.xlsx(file.path(dfp, "Human_Female_Steatohepatitis_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t7_all  <- read.xlsx(file.path(dfp, "Human_Female_Progression_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")
t8_all  <- read.xlsx(file.path(dfp, "Human_Female_FibrosisvF0_Highly_Reproducible_Meta_DEGs.xlsx"),
                     sheet = "All_Meta_Analysis_Results")

# Get the universes rom the full tables.
shared_universe <- Reduce(intersect, list(t1_all$gene,  t2_all$gene, t3_all$gene, 
                                            t4_all$gene, t5_all$gene, t6_all$gene,
                                            t7_all$gene, t8_all$gene))

# Get the universes rom the full tables.
# Males
shared_universe_male <- Reduce(intersect, list(t1_all$gene,  t2_all$gene, t3_all$gene, 
                                          t4_all$gene))
# Get the universes rom the full tables.
# Females
shared_universe_female <- Reduce(intersect, list(t5_all$gene, t6_all$gene,
                                          t7_all$gene, t8_all$gene))



cc_go_hs <- compareCluster(
  geneCluster   = human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cc_go_hs1 <- compareCluster(
  geneCluster   = human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cc_go_hs2 <- compareCluster(
  geneCluster   = human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)


# Males
M_cc_go_hs <- compareCluster(
  geneCluster   = Male_human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe_male),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

M_cc_go_hs1 <- compareCluster(
  geneCluster   = Male_human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe_male),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

M_cc_go_hs2 <- compareCluster(
  geneCluster   = Male_human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe_male),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Females
F_cc_go_hs <- compareCluster(
  geneCluster   = Female_human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe_female),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

F_cc_go_hs1 <- compareCluster(
  geneCluster   = Female_human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe_female),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

F_cc_go_hs2 <- compareCluster(
  geneCluster   = Female_human_cores,
  fun           = "enrichGO",
  universe      = to_entrez_hs(shared_universe_female),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

dp <- dotplot(cc_go_hs, showCategory = 8) 
# Adjust font sizes across different text elements
dp + theme(
  # Shrink pathway terms on the Y-axis (e.g., set to size 8 or 9)
  axis.text.y = element_text(size = 10),
)

dotplot(cc_go_hs1, showCategory = 8)
dotplot(cc_go_hs2, showCategory = 8)
# Males
dotplot(M_cc_go_hs, showCategory = 8)
dotplot(M_cc_go_hs1, showCategory = 8)
dotplot(M_cc_go_hs2, showCategory = 8)
# Females
dotplot(F_cc_go_hs, showCategory = 8)
dotplot(F_cc_go_hs1, showCategory = 8)
dotplot(F_cc_go_hs2, showCategory = 8)




# get df for the results
CC_BP_results <- as.data.frame(cc_go_hs)
CC_CC_results <- as.data.frame(cc_go_hs1)
CC_MF_results <- as.data.frame(cc_go_hs2)
# Males
Male_CC_BP_results <- as.data.frame(M_cc_go_hs)
Male_CC_CC_results <- as.data.frame(M_cc_go_hs1)
Male_CC_MF_results <- as.data.frame(M_cc_go_hs2)
# Females
Female_CC_BP_results <- as.data.frame(F_cc_go_hs)
Female_CC_CC_results <- as.data.frame(F_cc_go_hs1)
Female_CC_MF_results <- as.data.frame(F_cc_go_hs2)

# AUTOMATIC EXPORT TO EXCEL (SUPPLEMENTARY TABLE S1)
wb <- createWorkbook()
# CC enrichment
addWorksheet(wb, "CC_BP")
addWorksheet(wb, "CC_CC")
addWorksheet(wb, "CC_MF")
# Male CC enrichment
addWorksheet(wb, "Male_CC_BP")
addWorksheet(wb, "Male_CC_CC")
addWorksheet(wb, "Male_CC_MF")
# Male CC enrichment
addWorksheet(wb, "Female_CC_BP")
addWorksheet(wb, "Female_CC_CC")
addWorksheet(wb, "Female_CC_MF")

# writedata for GSEA enrichment
writeData(wb, "CC_BP",  CC_BP_results)
writeData(wb, "CC_CC",    CC_CC_results)
writeData(wb, "CC_MF",  CC_MF_results)
writeData(wb, "CC_BP",  Male_CC_BP_results)
writeData(wb, "CC_CC",    Male_CC_CC_results)
writeData(wb, "CC_MF",  Male_CC_MF_results)
writeData(wb, "CC_BP",  Female_CC_BP_results)
writeData(wb, "CC_CC",    Female_CC_CC_results)
writeData(wb, "CC_MF",  Female_CC_MF_results)

# Save the workbook to your directory
output_file <- file.path(dfp, "CompareCluster_Human_Meta_DEGs.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

cat(paste0("\nFile successfully exported to: ", output_file, "\n"))

# cnetplot
cnetplot(cc_go_hs, showCategory = 8)
cnetplot(cc_go_hs1, showCategory = 8)
cnetplot(cc_go_hs2, showCategory = 8)

cnetplot(M_cc_go_hs, showCategory = 8)
cnetplot(M_cc_go_hs1, showCategory = 8)
cnetplot(M_cc_go_hs2, showCategory = 8)

cnetplot(F_cc_go_hs, showCategory = 8)
cnetplot(F_cc_go_hs1, showCategory = 8)
cnetplot(F_cc_go_hs2, showCategory = 8)

