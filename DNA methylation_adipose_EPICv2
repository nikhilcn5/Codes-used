setwd("E:/EPIC dataset/april 2025/AT_epic_april_2025/")

# Load following libraries
library(ChAMP)
library(minfi)
library(ChAMPdata)
library(ggplot2)
library(qqman)
library(ggrepel)
library (dplyr)
library(EnhancedVolcano)

# Use the following codes:
# Loading data
idat.dir_AT = "E:/EPIC dataset/DNL project/O580691_07_12_2023/"
myLoad_AT_overall_sex_outlier <- champ.load(idat.dir_AT,arraytype="EPICv2", filterXY = F)

# Results
  #Your data set contains 937055 Meth probes
    #Filtering probes with a detection p-value above 0.01.
      #Removing 5784 probes.
    #Filtering probes with a beadcount <3 in at least 5% of samples.
      #Removing 1656 probes
    #Filtering NoCG Start
      #Only Keep CpGs, removing 3629 probes from the analysis.
    #Filtering SNPs Start, Using general mask options
      #Removing 31266 probes from the analysis.
    #Filtering MultiHit Start, Filtering probes that align to multiple locations as identified in Nordlund et al
      #Removing 0 probes from the analysis.

# Normalization
myNorm_AT_overall_sex_outlier <- champ.norm(beta=myLoad_AT_overall_sex_outlier$beta,arraytype="EPICv2",cores=8)

# Batch effect correction
myCombat1_AT_overall_sex_outlier <- champ.runCombat(beta=myNorm_AT_overall_sex_outlier,pd=myLoad_AT_overall_sex_outlier$pd,batchname=c("Slide"))
myCombat2_AT_overall_sex_outlier <- champ.runCombat(beta=myCombat1_AT_overall_sex_outlier,pd=myLoad_AT_overall_sex_outlier$pd,batchname=c("Array"))
myCombat3_AT_overall_sex_outlier <- champ.runCombat(beta=myCombat2_AT_overall_sex_outlier,pd=myLoad_AT_overall_sex_outlier$pd,batchname=c("Sex"))

# Plot the SVD plot to see if the corrections have worked
#champ.SVD(beta = myCombat3_AT_sex, pd = myLoad_AT_sex$pd, resultsDir = "./CHAMP_SVDimages/SVDsummary_combat3_AT_sex.pdf")

# MDS plots
mdsPlot(myCombat2_AT_overall_sex_outlier, sampGroups = myLoad_AT_overall_sex_outlier$pd$Sample_Group, sampNames = myLoad_AT_overall_sex_outlier$pd$Sample_Name)
mdsPlot(myCombat2_AT_overall_sex_outlier, sampGroups = myLoad_AT_overall_sex_outlier$pd$Sex, sampNames = myLoad_AT_overall_sex_outlier$pd$Sample_Name)

# DMP results
myDMP_combat3_AT_overall_sex_outlier = champ.DMP(beta = myCombat3_AT_overall_sex_outlier, pheno = myLoad_AT_overall_sex_outlier$pd$Sample_Group, compare.group = c("Normal", "T2D"), arraytype = "EPICv2", adjPVal = 0.05, adjust.method = "fdr")
write.csv(myDMP_combat3_AT_overall_sex_outlier, "./AT_Overall_DMP_fdr_overall_sex_outlier_analysis.csv")

# To include all the probes without a cutoff of p Value - Combat4
myDMP_combat4_AT_overall_sex_outlier = champ.DMP(beta = myCombat3_AT_overall_sex_outlier, pheno = myLoad_AT_overall_sex_outlier$pd$Sample_Group, compare.group = c("Normal", "T2D"), arraytype = "EPICv2", adjPVal = 1, adjust.method = "none" )
write.csv(myDMP_combat4_AT_overall_sex_outlier, "./AT_Overall_combat4_DMP_overall_sex_outlier_correction-none_analysis.csv")

# To include all the probes without a cutoff of p Value AND no correction for Sex - Combat5
myDMP_combat5_AT_overall_sex_outlier = champ.DMP(beta = myCombat2_AT_overall_sex_outlier, pheno = myLoad_AT_overall_sex_outlier$pd$Sample_Group, compare.group = c("Normal", "T2D"), arraytype = "EPICv2", adjPVal = 1, adjust.method = "none" )
write.csv(myDMP_combat5_AT_overall_sex_outlier, "./AT_Overall_combat5_DMP_correction-none_overall_sex_outlier_analysis.csv")


# Manhattan plots

# 1. Subset to autosomes + X + Y
# for combat4 data
keep.chrs <- c(as.character(1:22), "X", "Y")
dmp_df <- subset(myDMP_combat4_all_AT_df, CHR %in% keep.chrs)

# 2. Recode X→23, Y→24, then coerce to numeric
dmp_df$CHR <- with(dmp_df,
                   ifelse(CHR == "X", 23,
                          ifelse(CHR == "Y", 24,
                                 as.numeric(as.character(CHR))
                          ))
)

# 3. Build the qqman data.frame
dmp_all_AT_4 <- data.frame(
  SNP = dmp_df$Name,
  CHR = dmp_df$CHR,
  BP  = as.integer(dmp_df$MAPINFO),
  P   = as.numeric(dmp_df$P.Value)
)
# 1) Open the PNG device
png(
  filename = "Manhattan_chr1-22_XY.png",
  width    = 6,    # width in inches
  height   = 4,    # height in inches
  units    = "in",
  res      = 300   # resolution in DPI
)

# 2) Re-run your manhattan plot

manhattan(
  dmp_all_AT_4,
  chr            = "CHR",
  bp             = "BP",
  p              = "P",
  snp            = "SNP",
  genomewideline = -log10(5e-8),
  suggestiveline = FALSE,
  col            = c("gray30","gray60"),
  main           = "Manhattan: chr1–22, X & Y",
  annotatePval   = TRUE
)

# 3) Close the device
dev.off()


#Volcano plot


df1 <- myDMP_combat4_all_AT_df %>%
  mutate(
    log10p      = -log10(P.Value),
    significance = case_when(
      P.Value < 5e-8 & abs(deltaBeta) >= 0 ~ "Significant",
      TRUE                                     ~ "Not Significant"
    )
  )

ggplot(df1, aes(x = deltaBeta, y = log10p, color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Significant" = "firebrick", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.0, 0.0), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Normal vs T2D (Adipose Tissue)",
    x = expression(Delta~"Beta (T2D - Normal)"),
    y = expression(-log[10](italic(p))),
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") 

#####################
# for combat 5 data
myDMP_combat5_overall_sex_outlier_df = myDMP_combat5_AT_overall_sex_outlier$Normal_to_T2D

myDMP_combat5_overall_sex_outlier_df$CHR <- sub("^chr", "", myDMP_combat5_overall_sex_outlier_df$CHR)

dmp_df_5 = myDMP_combat5_overall_sex_outlier_df
table(dmp_df_5$CHR)
# Recode X→23, Y→24, then coerce to numeric
dmp_df_5$CHR[dmp_df_5$CHR == "X"] <- "23"
dmp_df_5$CHR[dmp_df_5$CHR == "Y"] <- "24"
dmp_df_5$CHR <- as.numeric(dmp_df_5$CHR)
dmp_df_5 = subset(dmp_df_5, CHR != "0")

# 4. Check
table(dmp_df_5$CHR)


# 3. Build the qqman data.frame
dmp_df_5 <- data.frame(
  SNP = dmp_df_5$Name,
  CHR = as.numeric(dmp_df_5$CHR),
  BP  = as.integer(dmp_df_5$MAPINFO),
  P   = as.numeric(dmp_df_5$P.Value)
)
# 1) Open the PNG device
png(
  filename = "Manhattan_chr1-22_XY_combat5_sex-outliers-removed.png",
  width    = 6,    # width in inches
  height   = 4,    # height in inches
  units    = "in",
  res      = 300   # resolution in DPI
)

# 2) Re-run your manhattan plot

manhattan(
  dmp_df_5,
  chr            = "CHR",
  bp             = "BP",
  p              = "P",
  snp            = "SNP",
  genomewideline = -log10(5.9e-8),
  suggestiveline = FALSE,
  col            = c("gray30","gray60"),
  main           = "Manhattan: chr1–22, X & Y",
  annotatePval   = TRUE
)

# 3) Close the device
dev.off()


#Volcano plot
library(ggrepel)
library (dplyr)
library(EnhancedVolcano)

df1_5 <- myDMP_combat5_overall_sex_outlier_df %>%
  mutate(
    log10p      = -log10(P.Value),
    significance = case_when(
      P.Value < 5.9e-8 & abs(deltaBeta) >= 0 ~ "Significant",
      TRUE                                     ~ "Not Significant"
    )
  )

ggplot(df1_5, aes(x = deltaBeta, y = log10p, color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Significant" = "firebrick", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.0, 0.0), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(5.9e-8), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Normal vs T2D (Adipose Tissue)",
    x = expression(Delta~"Beta (T2D - Normal)"),
    y = expression(-log[10](italic(p))),
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")


