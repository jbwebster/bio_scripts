
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(plyr)
library(dplyr)
setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")






master <- readRDS(file="data/TCGA/master.rds")
df <- master
scaled <- readRDS(file="data/TCGA/pancancer.pam50.scaled.rds")
robust <- readRDS(file="data/TCGA/pancancer.pam50.robust.rds")
pam50.none <- readRDS(file="data/TCGA/pancancer.pam50.none.rds")
x <- subset(pam50.none, Cohort=="BRCA")
x$Patient <- gsub("\\.", "-", rownames(x))
x$Patient <- gsub("-[0-9A-Z]*$", "", x$Patient)
x <- subset(x, Patient %in% df$Patient)
table(df$PAM50.cbioportal)
table(subset(scaled, scaled$Cohort=="BRCA")$SubtypeCall)
table(subset(robust, robust$Cohort=="BRCA")$SubtypeCall)
table(subset(pam50.none, pam50.none$Cohort=="BRCA")$SubtypeCall)

plotOne <- function(df, subtype) {
  df <- subset(df, df$PAM50.cbioportal != "Normal" & !is.na(df$PAM50.cbioportal))
  df$BCAR4_log <- log2(df$BCAR4.fpkm + 1)
  med <- median(df$BCAR4_log)
  df$BCAR4_centered <- df$BCAR4_log - med
  x <- subset(df, df$PAM50.cbioportal==subtype & !is.na(df$IHC.HER2))
  x$IHC.HER2 <- factor(x$IHC.HER2, levels=c("Positive", "Equivocal", "Negative", "Indeterminate"))
  p <- ggplot(x, aes(x=IHC.HER2, y=BCAR4_centered)) +
    geom_boxplot() +
    ylim(-0.25, 3) +
    labs(x="IHC-HER2 Status",
         y="Normalized BCAR4",
         title=subtype) 
  p + stat_compare_means()
}
plotTwo <- function(df) {
  #df <- subset(df, df$PAM50.cbioportal)
  df <- subset(df, df$PAM50.cbioportal != "Normal" & !is.na(df$PAM50.cbioportal))
  df$BCAR4_log <- log2(df$BCAR4.fpkm + 1)
  med <- median(df$BCAR4_log)
  df$BCAR4_centered <- df$BCAR4_log - med
  df$HER2Group <- ifelse(df$IHC.HER2=="Positive", "HER2+",
                         ifelse(!is.na(df$IHC.HER2), "NOT-POSITIVE", "Other"))
  df$Group <- paste0(df$PAM50.cbioportal, "_", df$HER2Group)
  df <- subset(df, df$HER2Group != "Other")
  df$Group <- factor(df$Group, 
                     levels = c("Basal_HER2+",
                                "Her2_HER2+",
                                "LumA_HER2+",
                                "LumB_HER2+",
                                "Basal_NOT-POSITIVE",
                                "Her2_NOT-POSITIVE",
                                "LumA_NOT-POSITIVE",
                                "LumB_NOT-POSITIVE"))
  p_meds <- ddply(df, .(Group), summarise, med = round(median(BCAR4_centered), digits=3))
  p_counts <- ddply(df, .(Group), summarise, counts = length(BCAR4_centered))

  p <- ggplot(df, aes(x=Group,y=BCAR4_centered)) +
    geom_jitter(alpha=0.1) +
    geom_boxplot(outlier.shape=NA, alpha=0.9, aes(fill=PAM50.cbioportal)) +
    stat_compare_means() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(data = p_meds, aes(x = Group, y = 1, label = med), 
              size = 3, vjust = -1.5) +
    geom_text(data = p_counts, aes(x = Group, y = 0.9, label = paste0("n=",counts)), 
              size = 3, vjust = -1.5) + 
    labs(x="PAM50",
         y="Normalized BCAR4",
         title="BCAR4 Expression by PAM50 and IHC-HER2 Status",
         caption="FPKM values were log2-transformed and median-centered\n
         Excludes N/A's and samples labeled as 'Normal'") +
    scale_y_continuous(trans='log2')

}

plotThree <- function(df) {
  #df <- subset(df, df$PAM50.cbioportal)
  df <- subset(df, df$PAM50.cbioportal != "Normal" & !is.na(df$PAM50.cbioportal))
  df$BCAR4_log <- log2(df$BCAR4.fpkm + 1)
  med <- median(df$BCAR4_log)
  df$BCAR4_centered <- df$BCAR4_log - med
  df$HER2Group <- ifelse(df$IHC.HER2=="Positive", "HER2+",
                         ifelse(!is.na(df$IHC.HER2), "NOT-POSITIVE", "Other"))
  #df$Group <- paste0(df$PAM50.cbioportal, "_", df$HER2Group)
  #df <- subset(df, df$HER2Group != "Other")
  #df$Group <- factor(df$Group, 
  #                   levels = c("Basal_HER2+",
  #                              "Her2_HER2+",
  #                              "LumA_HER2+",
  #                              "LumB_HER2+",
  #                              "Basal_NOT-POSITIVE",
  #                              "Her2_NOT-POSITIVE",
  #                              "LumA_NOT-POSITIVE",
  #                              "LumB_NOT-POSITIVE"))
  p_meds <- ddply(df, .(HER2Group), summarise, med = round(median(BCAR4_centered), digits=3))
  p_counts <- ddply(df, .(HER2Group), summarise, counts = length(BCAR4_centered))
  
  p <- ggplot(df, aes(x=HER2Group,y=BCAR4_centered)) +
    geom_jitter(alpha=0.1) +
    geom_boxplot(outlier.shape=NA, alpha=0.9) +
    stat_compare_means() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(data = p_meds, aes(x = HER2Group, y = 1, label = med), 
              size = 3, vjust = -1.5) +
    geom_text(data = p_counts, aes(x = HER2Group, y = 0.9, label = paste0("n=",counts)), 
              size = 3, vjust = -1.5) + 
    labs(x="PAM50",
         y="Normalized BCAR4",
         title="BCAR4 Expression by PAM50 and IHC-HER2 Status",
         caption="FPKM values were log2-transformed and median-centered\n
         Excludes N/A's and samples labeled as 'Normal'") +
    scale_y_continuous(trans='log10')
  
}

basal <- plotOne(df, "Basal")
her2 <- plotOne(df, "Her2")
luma <- plotOne(df, "LumA")
lumb <- plotOne(df, "LumB")
grid.arrange(basal, her2, luma, lumb, 
             top="BCAR4 Expression, Grouped by PAM50 and IHC-HER2",
             bottom="Normalized, log2-transformed, and median-centered")
p <- plotTwo(master)
p
p <- plotThree(master)
p

df <- master
df <- subset(df, df$PAM50.cbioportal != "Normal" & !is.na(df$PAM50.cbioportal))
df$BCAR4_log <- log2(df$BCAR4.fpkm + 1)
med <- median(df$BCAR4_log)
df$BCAR4_centered <- df$BCAR4_log - med
df$BCAR4_centered <- df$BCAR4_log - med
df$HER2Group <- ifelse(df$IHC.HER2=="Positive", "HER2+",
                       ifelse(!is.na(df$IHC.HER2), "NOT-POSITIVE", "Other"))
df$Group <- paste0(df$PAM50.cbioportal, "_", df$HER2Group)
df <- subset(df, df$HER2Group != "Other")


x <- df %>%
  group_by(Group) %>%
  summarise(median = median(BCAR4_centered))

ggplot(df, aes(x=PAM50.cbioportal, y=BCAR4_centered)) +
  geom_boxplot()



