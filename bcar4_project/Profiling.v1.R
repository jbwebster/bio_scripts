

library(ggplot2)
library(ggrepel)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")

## Creating a BCAR4+ Profile


## Cell line data from lab
raw <- read.csv("data/Profiling/Microarray_results_All_genes.FromNicole.csv", 
                header=T, sep=",")
keep <- c("ProbeName", "SystematicName", "GeneName",
          "HME_CONTROL_1", "HME_CONTROL_2",
          "TUH_CONTROL_1", "TUH_CONTROL_2")

df <- raw[,keep]

pam50.centroids <- read.csv("scripts/PAM50/PAM50_R/pam50_centroids.txt",
                            header=T, sep="\t")
pam50.genes <- c(pam50.centroids$X, "BCAR4")
#NUF2 = CDCA1
#NDC80 = KNTC2
#ORC6 = ORC6L
# ??? = PRG
# ??? = TMEM45B
pam50.genes[which(pam50.genes=="CDCA1")] = "NUF2"
pam50.genes[which(pam50.genes=="KNTC2")] = "NDC80"
pam50.genes[which(pam50.genes=="ORC6L")] = "ORC6"

df.minimal <- subset(df, df$GeneName %in% pam50.genes)[3:7]
rownames(df.minimal) = df.minimal$GeneName
df.minimal = df.minimal[,keep[4:7]]
df.minimal$HME_Avg = (df.minimal$HME_CONTROL_1 + df.minimal$HME_CONTROL_2) / 2
df.minimal$TUH_Avg = (df.minimal$TUH_CONTROL_1 + df.minimal$TUH_CONTROL_2) / 2
df.minimal$HME_Avg = scale(df.minimal$HME_Avg - median(df.minimal$HME_Avg))
df.minimal$TUH_Avg = scale(df.minimal$TUH_Avg - median(df.minimal$TUH_Avg))
df.minimal$Gene = rownames(df.minimal)

ggplot(df.minimal) +
  geom_point(aes(x=HME_Avg, y=TUH_Avg)) +
  geom_text_repel(data=subset(df.minimal, abs(df.minimal$HME_Avg-df.minimal$TUH_Avg)>0.75), 
            aes(x=HME_Avg, y=TUH_Avg, label=Gene)) +
  theme(panel.grid.minor = element_blank())


exp <- read.csv("data/TCGA/TCGA-BRCA.htseq_fpkm.roi.tsv",sep="\t")
rownames(exp) <- exp$Ensembl_ID
texp <- data.frame(t(exp[,2:1218]))
ntexp <- subset(texp, texp$ENSG00000262117.4==0)


m = read.csv("data/TCGA/microarray/AgilentG4502A_07_3_BRCA.tsv", sep="\t")
rownames(m) <- m$sample
mt <- data.frame(t(m))



#############################
# Predictions from Python machine learning profiler
df = read.csv("data/TCGA/microarray/pred.tsv", sep=",", header=T)
colnames(df) <- c("Sample", "Prediction")

exp <- read.csv("data/TCGA/TCGA-BRCA.htseq_fpkm.roi.tsv",sep="\t")
rownames(exp) <- exp$Ensembl_ID
texp <- data.frame(t(exp[,2:1218]))
texp$Sample = gsub("[A-Z]$", "", gsub("\\.", "-", rownames(texp)))
ntexp <- subset(texp, texp$Sample %in% df$Sample)

df <- subset(df, df$Sample %in% ntexp$Sample)
data = merge(df, ntexp)


