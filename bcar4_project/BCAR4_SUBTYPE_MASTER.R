

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("edgeR")

library(DESeq2)
library(edgeR)
library(data.table)
library(ggplot2)
library(forcats)
library(survminer)
library(survival)
library(dplyr)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")



################################################################################
########################### CREATION OF MASTER DF ##############################
# ###### IMPORT AND NORMALIZED READ COUNT DATA FOR ALL COHORTS ######
# # Using htseq count data donloaded from Xena Browser, in the GDC Data Hub
# files <- list.files(path="data/TCGA/pancancer_expression", 
#                     pattern = "counts.tsv.gz.tsv", full.names=T)
# 
# f = files[1]
# raw <- read.csv(f, sep="\t", header=T)
# df <- data.frame("Ensembl_ID"=raw$Ensembl_ID)
# for(f in files){
#   message(f)
#   raw <- read.csv(f, sep="\t", header=T)
#   #"ENSG00000262117.4" == BCAR4
#   # "ENSG00000141736.12" == ERBB2
#   df <- merge(df, raw, by="Ensembl_ID")
# }
# #saveRDS(df, file="data/TCGA/pancancer_expression/ALLCOHORTCOUNTS.rds")
# #rm(df)
# df <- readRDS(file="data/TCGA/pancancer_expression/ALLCOHORTCOUNTS.rds")
# # Randomly drop some genes to save memory and avoid errors
# # Keep a sufficiently large amount so as to maintain the assumption
# # that the data includes at least some non-differentially expressed genes
# set.seed(42)
# x <- df[sample(nrow(df), 9998), ]
# tmp <- subset(df, df$Ensembl_ID=="ENSG00000262117.4" | df$Ensembl_ID=="ENSG00000141736.12")
# x <-rbind(x,tmp)
# # Confirm that BCAR4 and ERBB2 didn't get dropped
# #"ENSG00000262117.4" %in% x$Ensembl_ID
# #"ENSG00000141736.12" %in% x$Ensembl_ID
# genes <- x$Ensembl_ID
# 
# x <- subset(x, select=-c(Ensembl_ID))
# meta <- data.frame("Samples" = colnames(x))
# meta$ExampleGroup <- c(rep("A", 5000), rep("B", 6057))
# #Undo the log2(counts + 1) normalization done by Xena to count data
# #https://www.biostars.org/p/383669/
# y <- round(((2^as.matrix(x)) - 1), 0)
# dds <- DESeqDataSetFromMatrix(countData=y,
#                              colData = meta,
#                              design = ~ ExampleGroup)
# dds <- estimateSizeFactors(dds)
# normalized_counts <- counts(dds, normalized=T)
# n <- data.frame(t(normalized_counts))
# colnames(n) <- genes
# n <- subset(n, select=c(ENSG00000262117.4, ENSG00000141736.12))
# colnames(n) <- c("BCAR4.deseq2.norm", "ERBB2.deseq2.norm")
# n$Patient <- rownames(n)
# saveRDS(n, file="data/TCGA/pancancer_expression/TCGA-deseq2-normalized.rds")
# 
# #Alternatively, try edgeR normalization, as it allows comparisons within a sample
# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# df <- readRDS(file="data/TCGA/pancancer_expression/ALLCOHORTCOUNTS.rds")
# genes <- df$Ensembl_ID
# df <- subset(df, select=-c(Ensembl_ID))
# d <- DGEList(counts=as.matrix(df))
# d <- calcNormFactors(d)
# normalized <- d$counts
# ndf <- data.frame(normalized)
# rownames(ndf) <- genes
# ndf <- ndf[c("ENSG00000262117.4", "ENSG00000141736.12"),]
# x <- data.frame(t(ndf))
# x$Patient <- rownames(x)
# colnames(x) <- c("BCAR4.exp.edgeR", "ERBB2.exp.edgeR", "Patient")
# v1 <- readRDS(file="data/TCGA/pancancer_expression/TCGA-deseq2-normalized.rds")
# m <- merge(x, v1)
# saveRDS(m, file="data/TCGA/pancancer_expression/TCGA-deseq2-edge-normalized.rds")
# 
# #### Gather clinical info
# files <- list.files(path="data/TCGA/pancancer_phenotype",
#                     pattern = "*.txt", full.names=T)
# clinical <- data.frame()
# for(f in files){
#   cohort <- gsub("_survival.txt","",gsub(".*phenotype/","",f))
#   message(cohort)
#   raw <- read.csv(f, sep="\t", header=T)
#   raw$cohort <- cohort
#   colnames(raw) <- c("Sample", "Patient", "OS", "OS.time",
#                      "DSS", "DSS.time", "DFI", "DFI.time",
#                      "PFI", "PFI.time", "Redaction", "Cohort")
#   raw <- raw[!duplicated(raw$Patient),]
#   clinical <- rbind(clinical, raw)
# }
# 
# ###### Combine all #####
# # DESeq2 custom normalized data
# exp <- readRDS("data/TCGA/pancancer_expression/TCGA-deseq2-edge-normalized.rds")
# exp$Patient <- gsub("-[A-Z0-9a-z]{3}$", "", gsub("\\.", "-", exp$Patient))
# colnames(exp) <- c("Patient", "BCAR4.exp.edgeR", "ERBB2.exp.edgeR", "BCAR4.norm.deseq2", "ERBB2.norm.deseq2")
# 
# # TCGA's BRCA PAM50
# library(TCGAbiolinks)
# dataSubt <- TCGAquery_subtype(tumor = "BRCA")
# brca_pam50 <- dataSubt[,c("patient", "BRCA_Subtype_PAM50")]
# colnames(brca_pam50) <- c("Patient", "BRCA.PAM50")
# 
# # Gistic CNV from Xena
# cnv <-  readRDS("data/TCGA/pancancer_gistic2/BCAR4.ERBB2.allCohorts.rds")
# colnames(cnv) <- c("Patient", "BCAR4.cnv", "ERBB2.cnv", "Cohort")
# cnv$Patient <- gsub("-[A-Z0-9a-z]{2}$", "", gsub("\\.", "-", cnv$Patient))
# 
# # Merge
# master <- merge(exp, clinical, by="Patient")
# master <- merge(master, cnv, by="Patient")
# master <- merge(brca_pam50, master, by="Patient", all=T)
# keep <- c("Patient", "Sample", "Cohort.x", "BCAR4.exp.edgeR",
#           "ERBB2.exp.edgeR", "BCAR4.norm.deseq2", "ERBB2.norm.deseq2",
#           "BCAR4.cnv", "ERBB2.cnv",
#           "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time",
#           "PFI", "PFI.time", "Redaction", "BRCA.PAM50")
# master <- master[,keep]
# colnames(master)[3] <- "Cohort"
# master$Cohort <- factor(master$Cohort)
# master$BRCA.PAM50 <- factor(master$BRCA.PAM50)
# fusions <- read.csv("data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=F)
# fusions <- subset(fusions, fusions$V7 %like% "BCAR4")
# fusions$Patient <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$", "", fusions$V19)
# master$hasBCAR4Fusion <- ifelse(master$Patient %in% fusions$Patient, 1, 0)
# master$hasBCAR4Fusion <- factor(master$hasBCAR4Fusion)
# 
# saveRDS(master, file="data/TCGA/MASTER/master.April11.rds")
# rm(list=ls())
################################################################################


################################################################################
########################### EXPRESSION PROFILES ################################

master <- readRDS("data/TCGA/MASTER/master.April11.rds")
figurepath <- "figures/ProperExpressionUnits/"

brca <- subset(master, master$Cohort=="BRCA")
set.seed(42)
top25v <- summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["3rd Qu."]
brca$IsTop25P <- ifelse(brca$BCAR4.exp.edgeR>=top25v, "Yes", "No")

#PanCanExpression.png
ggplot() +
  geom_jitter(data=subset(master, !is.na(master$Cohort)), 
              aes(x=Cohort, y=BCAR4.exp.edgeR, color=hasBCAR4Fusion),
              alpha=0.5) +
  geom_violin(data=subset(master, !is.na(master$Cohort)),
              aes(x=Cohort, y=BCAR4.exp.edgeR+1),
              alpha=0.5) +
  geom_hline(yintercept=(min(subset(brca, brca$IsTop25P=="Yes")$BCAR4.exp.edgeR))) +
  labs(x="TCGA Cohort",
       y="BCAR4 TMM Normalized Expression",
       title="BCAR4 Expression by Cohort",
       color="BCAR4 Fusion",
       caption="Line at top 25% of BCAR4 expressing BRCA samples") +
  scale_color_hue(labels=c("No", "Yes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#BCAR4ERBB2DistributionMedians.png
#Do 25% of samples fall in each quadrant?
brca <- subset(master, master$Cohort=="BRCA")
ggplot(brca) +
  geom_point(aes(x=ERBB2.exp.edgeR, y=ERBB2.cnv))


bmedian <- median(brca$BCAR4.exp.edgeR)
emedian <- median(brca$ERBB2.exp.edgeR)
middle <- nrow(brca) / 2
set.seed(42)
brca$BCAR4Rank <- rank(brca$BCAR4.exp.edgeR, ties.method="random")
brca$ERBB2Rank <- rank(brca$ERBB2.exp.edgeR, ties.method="random")
brca$Comparison <- ifelse(brca$BCAR4Rank > middle & brca$ERBB2Rank > middle, "TR",
                          ifelse(brca$BCAR4Rank > middle & brca$ERBB2Rank < middle, "LR",
                                 ifelse(brca$BCAR4Rank < middle & brca$ERBB2Rank > middle, "TL",
                                        ifelse(brca$BCAR4Rank < middle & brca$ERBB2Rank < middle, "LL", "Other"))))
ggplot() +
  geom_jitter(data=subset(brca, brca$Comparison!="Other"), aes(x=BCAR4.exp.edgeR, y=ERBB2.exp.edgeR, color=Comparison)) +
  geom_vline(xintercept=bmedian) +
  geom_hline(yintercept=emedian) +
  labs(title="BCAR4 vs ERBB2",
       caption="Lines represent medians")
brca$BComp <- ifelse(brca$BCAR4Rank>middle, "T", "L")
brca$EComp <- ifelse(brca$ERBB2Rank>middle, "T", "L")
brcasub <- subset(brca, brca$Comparison!="Other")
brcasub$Comparison <- factor(brcasub$Comparison)
table(brcasub$BComp, brcasub$EComp)
chisq.test(brcasub$BComp, brcasub$EComp)

tq <- nrow(brca) * 0.90
middle <- tq # For top quartile analysis
set.seed(42)
# Following is modified
brca$Comparison <- ifelse(brca$BCAR4Rank > middle & brca$ERBB2Rank > middle, "TR",
                          ifelse(brca$BCAR4Rank > middle & brca$ERBB2Rank < middle, "LR",
                                 ifelse(brca$BCAR4Rank < middle & brca$ERBB2Rank > middle, "TL",
                                        ifelse(brca$BCAR4Rank < middle & brca$ERBB2Rank < middle, "LL", "Other"))))
brca$BComp <- ifelse(brca$BCAR4Rank>middle, "T", "L")
brca$EComp <- ifelse(brca$ERBB2Rank>middle, "T", "L")
brcasub <- subset(brca, brca$Comparison!="Other")
brcasub$Comparison <- factor(brcasub$Comparison)
table(brcasub$BComp, brcasub$EComp)
chisq.test(brcasub$BComp, brcasub$EComp)

##
tq <- nrow(brca) * 0.90
bmiddle <- tq # For top quartile analysis
tq <- nrow(brca) * 0.75
emiddle <- tq
set.seed(42)
# Following is modified
brca$Comparison <- ifelse(brca$BCAR4Rank > bmiddle & brca$ERBB2Rank > emiddle, "TR",
                          ifelse(brca$BCAR4Rank > bmiddle & brca$ERBB2Rank < emiddle, "LR",
                                 ifelse(brca$BCAR4Rank < bmiddle & brca$ERBB2Rank > emiddle, "TL",
                                        ifelse(brca$BCAR4Rank < bmiddle & brca$ERBB2Rank < emiddle, "LL", "Other"))))
brca$BComp <- ifelse(brca$BCAR4Rank>bmiddle, "T", "L")
brca$EComp <- ifelse(brca$ERBB2Rank>emiddle, "T", "L")
brca$ERBB2.Amp <- factor(ifelse(brca$ERBB2.cnv>=1, "Yes", "No"))
brca$EComp <- brca$ERBB2.Amp
brcasub <- subset(brca, brca$Comparison!="Other")
brcasub$Comparison <- factor(brcasub$Comparison)
table(brcasub$BComp, brcasub$EComp)
chisq.test(brcasub$BComp, brcasub$EComp)
ggplot() +
  geom_jitter(data=subset(brca, brca$Comparison!="Other"), 
              aes(x=BCAR4.exp.edgeR, y=ERBB2.exp.edgeR, color=Comparison)) +
  geom_vline(xintercept=subset(brcasub, brcasub$BCAR4Rank==round(bmiddle))[,4]) +
  geom_hline(yintercept=subset(brcasub, brcasub$ERBB2Rank==round(emiddle))[,5]) +
  labs(title="BCAR4 vs ERBB2",
       caption="Lines represent medians")



#Instead create quadrants based on BCAR4 expression vs HER2 Amplification
tq <- nrow(brca) * 0.75
middle <- tq # For top quartile analysis
brca$ERBB2.Amp <- factor(ifelse(brca$ERBB2.cnv>=1, "Yes", "No"))
brca$Comparison <- ifelse(brca$BCAR4Rank>middle & brca$ERBB2.Amp=="Yes", "TR",
                          ifelse(brca$BCAR4Rank>middle & brca$ERBB2.Amp=="No", "TL",
                                 ifelse(brca$BCAR4Rank<middle & brca$ERBB2.Amp=="Yes", "LR",
                                        ifelse(brca$BCAR4Rank<middle & brca$ERBB2.Amp=="No", "LL", "Other"))))
brca$BComp <- ifelse(brca$BCAR4Rank>middle, "T", "L")
brcasub <- subset(brca, brca$Comparison!="Other")
table(brcasub$BComp, brcasub$ERBB2.Amp)
chisq.test(brcasub$BComp, brcasub$EComp)

ggplot() +
  geom_point(data=brca,
             aes(x=BCAR4.exp.edgeR, y=ERBB2.Amp),
             alpha=0.2) +
  geom_vline(xintercept=bmedian+1) +
  labs(title="BCAR4 Exp vs ERBB2 Amp",
       caption="Line is median BCAR4 exp")

# What if it was done using top quartile? Instead of median?
topquartile <- summary(brca$BCAR4.exp.edgeR)["3rd Qu."]
brca$Comparison <- ifelse(brca$BCAR4.exp.edgeR>topquartile & brca$ERBB2.Amp=="Yes", "TR",
                          ifelse(brca$BCAR4.exp.edgeR>topquartile & brca$ERBB2.Amp=="No", "TL",
                                 ifelse(brca$BCAR4.exp.edgeR<topquartile & brca$ERBB2.Amp=="Yes", "LR",
                                        ifelse(brca$BCAR4.exp.edgeR<topquartile & brca$ERBB2.Amp=="No", "LL", "Other"))))
obs <- summary(factor(subset(brca, brca$Comparison!="Other")$Comparison))
summary(subset(brca, brca$BCAR4.exp.edgeR!=2)$ERBB2.Amp) # To figure out expected values
exp <- c(476.5, 73.5, 476.4, 73.5)
chisq.test(obs,)


df <- subset(master, !is.na(master$Cohort) & master$BCAR4.exp.deseq2==0)
df2 <- subset(master, !is.na(master$Cohort) & master$BCAR4.exp.deseq2>0)
df$Group <- "None"
df2$BCAR4Rank <- rank(-df2$BCAR4.exp.deseq2)
df2$Group <- ifelse(df2$BCAR4.exp.deseq2 > (nrow(df2) / 2), "High", "Low")

#########
brca <- subset(master, master$Cohort=="BRCA")
set.seed(42)
top25v <- summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["3rd Qu."]
brca$IsTop25P <- ifelse(brca$BCAR4.exp.edgeR>=top25v, "Yes", "No")
m <- brca %>%
  group_by(BRCA.PAM50, IsTop25P) %>%
  summarise(count=n())
m <- na.omit(m)
n <- brca %>%
  group_by(BRCA.PAM50) %>%
  summarise(subtypeCount = n())
m <- merge(m, n)
o <- brca %>%
  group_by(IsTop25P) %>%
  summarise(GroupCount = n())
m <- merge(m, o)
m$Prop <- m$count / m$subtypeCount
m$GroupProp <- m$count / m$GroupCount
m$BRCA.PAM50 <- factor(m$BRCA.PAM50, levels=c("LumB", "Basal", "Her2", "LumA", "Normal"))
n <- subset(m, m$IsTop25P=="Yes")
ggplot(data=n) +
  geom_col(aes(x=BRCA.PAM50, y=Prop)) +
  geom_text(label=paste0(n$count, "/", n$subtypeCount), x=n$BRCA.PAM50, y=n$Prop+0.01) +
  scale_y_continuous(limits=c(0, 0.3)) +
  labs(title="BCAR4+ as a proportion of PAM50 subtypes",
       x="PAM50 Subtype",
       y="% of PAM50 Subtype Samples") 
ggplot(data=n) +
  geom_col(aes(x=BRCA.PAM50, y=GroupProp)) +
  labs(title="Proportion of BCAR4+ in each PAM50 subtype",
       x="PAM50 Subtype",
       y="% of BCAR4+ samples")
m$subtypeProp <- m$subtypeCount / 1196
m$X <- m$count / m$GroupCount
ggplot(data=subset(m, m$IsTop25P=="Yes")) +
  geom_point(aes(x=subtypeProp, y=GroupProp, color=BRCA.PAM50)) +
  geom_abline(slope=1, intercept=0) +
  scale_x_continuous(limits=c(0, 0.6)) +
  scale_y_continuous(limits=c(0, 0.6)) +
  labs(x="% of Cohort",
       y="% of BCAR4+ samples",
       caption="BCAR4=Top 25% of expressing samples",
       title="Subtype Enrichment")


ihc <- read.csv("data/TCGA/brca_tcga_clinical_data.tsv", sep="\t", header=T)
keep <- c("Patient.ID", "IHC.HER2",
          "ER.Status.By.IHC", "PR.status.by.ihc")
ihc <- ihc[,keep]
colnames(ihc) <- c("Patient",  "IHC.HER2",
                   "IHC.ER", "IHC.PR")
m <- merge(brca, ihc, all=T)

m$Group <- factor(paste0(m$BRCA.PAM50, "-", m$IHC.HER2))
n <- subset(m, m$IHC.HER2 != "Indeterminate" & !is.na(m$IHC.HER2) & !is.na(m$BRCA.PAM50))
n$Group <- factor(n$Group)
n$IHC.HER2 <- factor(n$IHC.HER2, levels=c("Positive", "Equivocal", "Negative"))
n$Mod.IHC.HER2 <- ifelse(n$IHC.HER2 == "Equivocal", "HER2-", "HER2+") # Like FGFR4 paper
#n$Mod.IHC.HER2 <- factor(n$Mod.IHC.HER2, levels=c("Positive", "Negative"))
ggplot() +
  geom_jitter(data=n, aes(x=Mod.IHC.HER2, y=BCAR4.exp.edgeR, color=BRCA.PAM50), height=0.1, width=0.1) +
  geom_boxplot(data=n, aes(x=Mod.IHC.HER2, y=BCAR4.exp.edgeR, color=BRCA.PAM50), alpha=0.4) +
  facet_wrap(~BRCA.PAM50) +
  labs(x="IHC-HER2",
       y="BCAR4 Exp")

n$ModGroup <-paste0(n$BRCA.PAM50, "_", n$Mod.IHC.HER2)
n$ModGroup <- factor(n$ModGroup, levels=c("Basal_HER2-", "Basal_HER2+",
                                          "Her2_HER2-", "Her2_HER2+",
                                          "LumA_HER2-", "LumA_HER2+",
                                          "LumB_HER2-", "LumB_HER2+",
                                          "Normal_HER2-", "Normal_HER2+"))
ggplot(data=subset(n, n$BRCA.PAM50!="Normal"), aes(x=ModGroup, y=BCAR4.exp.edgeR, color=BRCA.PAM50)) +
  geom_jitter(height=0.01,width=0.1) +
  geom_boxplot(data=subset(n,n$BRCA.PAM50!="Normal"),
               alpha=0.4, outlier.shape=NA, 
               aes(fill=BRCA.PAM50)) +
  scale_color_manual(values=c("red", "magenta", "blue", "lightblue")) +
  scale_fill_manual(values=c("red", "magenta", "blue", "lightblue")) +
  labs(x="",
       y="mRNA Expression",
       title="BCAR4") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "None")

ggplot(data=n, aes(x=Mod.IHC.HER2, y=BCAR4.exp.edgeR)) +
  geom_jitter(height=0,width=0.1) +
  geom_violin(alpha=0.4) +
  labs(x="IHC HER2 Status", 
       y="BCAR4 Expression",
       title="BCAR4 by clinical HER2 status")

ggplot(data=brca, aes(x=BCAR4.exp.edgeR, y=BCAR4.cnv)) +
  geom_jitter() +
  geom_smooth() +
  labs(x="BCAR4 Expression",
       y="BCAR4 Gistic Score",
       title="Expression vs Amplification of BCAR4")


################################################################################

################################################################################
########################### SURVIVAL CURVES ####################################

master <- readRDS("data/TCGA/MASTER/master.April11.rds")

pancan <- master
set.seed(42)
pancan$BRank <- rank(pancan$BCAR4.exp.edgeR, ties.method="random")
minrank <- min(subset(pancan, BCAR4.exp.edgeR>0)$BRank)
middle <- ((nrow(pancan)-minrank) / 2) + minrank
pancan$Group <- ifelse(pancan$BCAR4.exp.edgeR==0, "None",
                       ifelse(pancan$BRank<middle, "Low", "High"))

f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = pancan)
ggsurvplot(f, data = pancan, pval=T,
           risk.table=T) +
  labs(title="PanCancer BCAR4 Expression Tiers")

#
brca <- subset(master, master$Cohort=="BRCA" & !is.na(master$BCAR4.exp.edgeR))
brca$BRank <- rank(brca$BCAR4.exp.edgeR, ties.method="random")
minrank <- min(subset(brca, BCAR4.exp.edgeR>0)$BRank)
middle <- ((nrow(brca)-minrank) / 2) + minrank
brca$Group <- ifelse(brca$BCAR4.exp.edgeR==0, "None",
                     ifelse(brca$BRank<middle, "Low", "High"))
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = brca)
ggsurvplot(f, data = brca, pval=T,
           risk.table=T) +
  labs(title="BRCA BCAR4 Expression Tiers (50-50)")
top25v = summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["3rd Qu."]
top10r <- (((nrow(brca)-minrank) / 10)*9) + minrank
top10v <- subset(brca, brca$BRank==round(top10r))[,4]
ggplot(subset(brca, brca$BCAR4.exp.edgeR>0)) +
  geom_density(aes(BCAR4.exp.edgeR)) +
  geom_vline(xintercept=median(brca$BCAR4.exp.edgeR)) +
  geom_vline(xintercept=top25v, color="red") +
  geom_vline(xintercept=top10v, color="blue") +
  labs(title="BCAR4 Expression in non-0 samples", 
       caption="black=median, red=top 25%, blue=top 10%")
brca$Group <- ifelse(brca$BCAR4.exp.edgeR>top10v, "Top 10%",
                     ifelse(brca$BCAR4.exp.edgeR>top25v, "Top 25%", 
                     ifelse(brca$BCAR4.exp.edgeR>0, "Expressed", "None")))
brca$Group <- factor(brca$Group, levels=c("Top 10%", "Top 25%", "Expressed", "None"))
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = brca)
ggsurvplot(f, data = brca, pval=T,
           risk.table=T) +
  labs(title="BRCA BCAR4 Expression Tiers")

#
# Comparing BCAR4+ vs no expression in different subtypes
brca <- subset(master, master$Cohort=="BRCA")
set.seed(42)
top25v <- summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["3rd Qu."]
brca$IsTop25P <- ifelse(brca$BCAR4.exp.edgeR>=top25v, "Yes", "No")
brca$Group <- paste0(brca$BRCA.PAM50, "-", brca$IsTop25P)
s <- subset(brca, !is.na(brca$BRCA.PAM50))
s <- subset(s, s$IsTop25P=="Yes" | s$BCAR4.exp.edgeR==0)
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = s)
ggsurvplot(f, data = s, pval=T,
           risk.table=T) +
  labs(title="BRCA BCAR4 Expression Tiers")

s <- subset(brca, brca$BRCA.PAM50=="Her2" & (brca$IsTop25P=="Yes" | brca$BCAR4.exp.edgeR==0))
s$Group <- ifelse(s$IsTop25P=="Yes", "BCAR4+", "No Expression")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = s)
ggsurvplot(f, data = s, pval=T,
           risk.table=T) +
  labs(title="HER2")

s <- subset(brca, brca$BRCA.PAM50=="LumA" & (brca$IsTop25P=="Yes" | brca$BCAR4.exp.edgeR==0))
s$Group <- ifelse(s$IsTop25P=="Yes", "BCAR4+", "No Expression")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = s)
ggsurvplot(f, data = s, pval=T,
           risk.table=T) +
  labs(title="LumA")

s <- subset(brca, brca$BRCA.PAM50=="LumB" & (brca$IsTop25P=="Yes" | brca$BCAR4.exp.edgeR==0))
s$Group <- ifelse(s$IsTop25P=="Yes", "BCAR4+", "No Expression")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = s)
ggsurvplot(f, data = s, pval=T,
           risk.table=T) +
  labs(title="LumB")

s <- subset(brca, brca$BRCA.PAM50=="Basal" & (brca$IsTop25P=="Yes" | brca$BCAR4.exp.edgeR==0))
s$Group <- ifelse(s$IsTop25P=="Yes", "BCAR4+", "No Expression")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = s)
ggsurvplot(f, data = s, pval=T,
           risk.table=T) +
  labs(title="Basal")

s <- subset(brca, !is.na(brca$BRCA.PAM50))
s <- subset(s, s$IsTop25P=="Yes" | s$BCAR4.exp.edgeR==0)
s$Group <- ifelse(s$IsTop25P=="Yes", "BCAR4+", "No Expression")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = s)
ggsurvplot(f, data = s, pval=T,
           risk.table=T) +
  labs(title="BRCA BCAR4+ vs no BCAR4 Expression")

## IHC
brca <- subset(master, master$Cohort=="BRCA")
set.seed(42)
top25v <- summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["3rd Qu."]
brca$IsTop25P <- ifelse(brca$BCAR4.exp.edgeR>=top25v, "Yes", "No")

ihc <- read.csv("data/TCGA/brca_tcga_clinical_data.tsv", sep="\t", header=T)
keep <- c("Patient.ID", "IHC.HER2",
          "ER.Status.By.IHC", "PR.status.by.ihc")
ihc <- ihc[,keep]
colnames(ihc) <- c("Patient",  "IHC.HER2",
                   "IHC.ER", "IHC.PR")
m <- merge(brca, ihc, all=T)
m <- subset(m, !is.na(m$IsTop25P))
m$TripleNegative <- ifelse(m$IHC.ER=="Negative" & m$IHC.PR=="Negative" & m$IHC.HER2=="Negative",
                           "Yes", "No")
m$Group <- ifelse(m$IsTop25P=="Yes", "BCAR4+",
                  ifelse(m$TripleNegative=="Yes", "TripleNeg",
                         as.character(m$BRCA.PAM50)))
#m <- subset(m, m$Group != "Other")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = m)
ggsurvplot(f, data = m, pval=T,
           risk.table=T) +
  labs(title="BCAR4 Subtypes")

brca <- subset(master, master$Cohort=="BRCA")
set.seed(42)
top25v <- summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["3rd Qu."]
bottom25v <- summary(subset(brca, brca$BCAR4.exp.edgeR>0)$BCAR4.exp.edgeR)["1st Qu."]
brca$IsTop25P <- ifelse(brca$BCAR4.exp.edgeR>=top25v, "Yes", "No")
brca$IsBottom25P <- ifelse(brca$BCAR4.exp.edgeR==1, "Yes", "No")

ihc <- read.csv("data/TCGA/brca_tcga_clinical_data.tsv", sep="\t", header=T)
keep <- c("Patient.ID", "IHC.HER2",
          "ER.Status.By.IHC", "PR.status.by.ihc")
ihc <- ihc[,keep]
colnames(ihc) <- c("Patient",  "IHC.HER2",
                   "IHC.ER", "IHC.PR")
m <- merge(brca, ihc, all=T)
n <- subset(m, m$IsTop25P=="Yes" | m$IHC.HER2=="Positive")
n$Group <- ifelse(n$IsTop25P=="Yes", "BCAR4+", "IHC-HER2 Pos")
m$Group <- ifelse(m$IsTop25P=="Yes" & m$IHC.HER2=="Positive", "BCAR4+/IHC-HER2+",
                  ifelse(m$IsTop25P=="Yes" & m$IHC.HER2=="Negative", "BCAR4+ Only",
                         ifelse(m$BCAR4.exp.edgeR==0 & m$IHC.HER2=="Negative", "Neither",
                                ifelse(m$BCAR4.exp.edgeR==0 & m$IHC.HER2=="Positive", "IHC-HER2+ Only",
                                "Other"))))
m <- subset(m, m$Group !="Other")
n <- m
n$Group <- factor(n$Group, levels=c("BCAR4+/IHC-HER2+", "BCAR4+ Only", "IHC-HER2+ Only", "Neither"))
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = n)
ggsurvplot(f, data = n, pval=T,
           risk.table=T) +
  labs(title="BCAR4+ and IHC-HER2 Permutations")
nrow(subset(m, m$IHC.HER2=="Positive"))
nrow(subset(m, m$IHC.HER2=="Positive" & m$BRCA.PAM50=="Her2"))
nrow(subset(m, m$IsTop25P=="Yes"))
nrow(subset(m, m$IsTop25P=="Yes" & m$BRCA.PAM50=="Her2" & m$IHC.HER2=="Positive"))
nrow(subset(m, m$IsTop25P=="Yes" &  m$IHC.HER2=="Positive"))


n <- m
n$Group <- ifelse(n$BRCA.PAM50=="LumA" & n$IsTop25P=="Yes", "LumA/BCAR4+",
                  ifelse(n$BRCA.PAM50=="LumA" & n$BCAR4.exp.edgeR==1, "LumA/BCAR4-",
                         ifelse(n$BRCA.PAM50=="Her2" & n$BCAR4.exp.edgeR==1, "Her2/BCAR4-", "Other")))
n <- subset(n, n$Group!="Other")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = n)
ggsurvplot(f, data = n, pval=T,
           risk.table=T) +
  labs(title="LumA/BCAR4 vs HER2")



f <- survfit(Surv(OS.time / 365, OS) ~ BRCA.PAM50, data = brca)
ggsurvplot(f, data = brca, pval = T,
           risk.table = T) +
  labs(title="TCGA BRCA PAM50 Subtypes")

m <- merge(brca, ihc, all=T)
m$Group <- ifelse(m$IHC.HER2=="Positive", "HER2+",
                  ifelse(m$IHC.HER2=="Negative" | m$IHC.HER2=="Equivocal", "HER2-", "Other"))
m <- subset(m, m$Group != "Other")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = m)
ggsurvplot(f, data = m, pval = T,
           risk.table = T) +
  labs(title="TCGA IHC-HER2 Status")
