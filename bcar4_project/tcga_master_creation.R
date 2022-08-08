


library(dplyr)
library(TCGAbiolinks)
library(data.table)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")

# From xenabrowser.net
# Survival
survival <- read.csv("data/TCGA/TCGA-BRCA.survival.tsv", header=T, sep="\t")
colnames(survival) <- c("Sample", "OS", "Patient", "OS.time")
# Expression (ERBB2 and BCAR4 only - others available)
expression <- read.csv("data/TCGA/TCGA-BRCA.htseq_fpkm.roi.tsv", header=T, sep="\t")
expression_long <- data.frame(t(expression))
expression_long$Sample <- gsub("\\.", "-", rownames(expression_long))
expression_long$Sample <- gsub("-01[A-B]$", "", expression_long$Sample)
expression_long <- expression_long[-1,]
colnames(expression_long) <- c("ERBB2.fpkm", "BCAR4.fpkm", "Patient")

# Alterations
alterations <- read.csv("data/TCGA/cbioportal_alterations.tsv", sep="\t", header=T)
keep <- c("Patient.ID", "BCAR4..AMP", "ERBB2..MUT","ERBB2..AMP")
alterations <- alterations[,keep]
colnames(alterations) <- c("Patient", "BCAR4.amp.cbioportal", "ERBB2.mut.cbioportal", "ERBB2.amp.cbioportal")
# get subtype information
dataSubt <- TCGAquery_subtype(tumor = "BRCA")
keep <- c("patient", "BRCA_Subtype_PAM50")
dataSubt <- dataSubt[,keep]
colnames(dataSubt) <- c("Patient", "PAM50.cbioportal")
#IHC
ihc <- read.csv("data/TCGA/brca_tcga_clinical_data.tsv", sep="\t", header=T)
keep <- c("Patient.ID", "Disease.Free..Months.",
          "Disease.Free.Status", "IHC.HER2",
          "ER.Status.By.IHC", "PR.status.by.ihc",
          "Overall.Survival..Months.", "Overall.Survival.Status")
ihc <- ihc[,keep]
colnames(ihc) <- c("Patient", "Disease.Free.Months",
                   "Disease.Free.Status", "IHC.HER2",
                   "IHC.ER", "IHC.PR",
                   "Overall.Survival.Months", "Overall.Survival.Status")

a <- merge(expression_long, ihc)
a <- merge(a, alterations, all.x=T)
a <- merge(a, dataSubt, all.x=T)
master <- a
master$BCAR4.amp.cbioportal <- as.factor(master$BCAR4.amp.cbioportal)
master$ERBB2.amp.cbioportal <- as.factor(master$ERBB2.amp.cbioportal)
master$ERBB2.mut.cbioportal <- as.factor(master$ERBB2.mut.cbioportal)
master$IHC.HER2 <- as.factor(master$IHC.HER2)
master$IHC.ER <- as.factor(master$IHC.ER)
master$IHC.PR <- as.factor(master$IHC.PR)
master$PAM50.cbioportal <- as.factor(master$PAM50.cbioportal)
fusions <- c("TCGA-B6-A40B", "TCGA-BH-A0DG",
             "TCGA-C8-A1HL", "TCGA-LL-A740",
             "TCGA-AO-A0JL", "TCGA-A1-A0SO",
             "TCGA-D8-A27G", "TCGA-OL-A66N")
master$BCAR4_fusion <- ifelse(master$Patient %in% fusions, "Fusion", "NoFusion")
master$ERBB2.fpkm <- as.numeric(master$ERBB2.fpkm)
master$BCAR4.fpkm <- as.numeric(master$BCAR4.fpkm)
saveRDS(master, file="data/TCGA/master.rds")

