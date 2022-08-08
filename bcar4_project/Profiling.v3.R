

library(caret)
library(data.table)
library(dplyr)
library(stats)
library(RColorBrewer)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")

set.seed(42)

################
# Data import
microarray = read.csv("data/TCGA/microarray/AgilentG4502A_07_3_BRCA.tsv", sep="\t")
#fusions = read.csv("data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=F)
#fpkm = read.csv("data/TCGA/TCGA-BRCA.htseq_fpkm.roi.tsv", sep="\t")
experiments = read.csv("data/Profiling/Microarray_results_All_genes.FromNicole.csv")


# Get samples of interest
original_cols = colnames(experiments)
keep_cols = c("GeneName",
              original_cols[original_cols %like% "status" & (original_cols %like% "HME" | original_cols %like% "TUH")],
              original_cols[(original_cols %like% "_1" | original_cols %like% "_2") &
                              (original_cols %like% "HME" | original_cols %like% "TUH")])
experiment_subset = experiments[,keep_cols]


# Count DE
status_cols = c(keep_cols[keep_cols %like% "status"])
de_exp = experiment_subset[,status_cols]
de_exp$FreqChange = ifelse(de_exp$TUH_B4siR268.status == "unchanged", 0, 1)
de_exp$FreqChange = ifelse(de_exp$TUH_B4siR271.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_ZCBOE.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_LBOE.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_LBMUT585.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_B4OE.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$GeneName = as.factor(experiment_subset$GeneName)

# Label genes that are in TCGA microarray, since we only want to use those
de_exp$InMicroarray <- ifelse(de_exp$GeneName %in% microarray$sample, 1, 0)

# Select gene candidates
experiment_gene_candidates = subset(de_exp, 
                                    de_exp$FreqChange==4 & de_exp$InMicroarray==1)

# Isolate expression data
exp_cols = c(keep_cols[keep_cols %like% "_1" | keep_cols %like% "_2"])
exp_data = data.frame(experiment_subset[,exp_cols])
exp_data$GeneName = as.factor(experiment_subset$GeneName)

# Select data for model
gene_cand_exp = subset(exp_data, exp_data$GeneName %in% experiment_gene_candidates$GeneName)


# Reformat
genes = gene_cand_exp$GeneName
data = data.frame(t(gene_cand_exp))
colnames(data) <- genes
data <-data[-c(3,4,5,6,17),] # Remove GeneName and LBMUT rows. Also B4OE without fusion

# B4 Status
status = c(0,0,1,1,1,1, # HME. 2 is overexpression of B4 without fusion
           1,1,0,0,0,0) # TUH
data$B4Status = status

Xdata = data[,1:145]
data$B4Status = ifelse(data$B4Status == 0, "No", "Yes")
Ydata = as.factor(data[,146])

# Optional
#Xdata = subset(Xdata, select = c(rownames(imp)[order(imp$Overall)[1:50]]))


fitControl <- trainControl(
  method = "LOOCV",
  number = 1,
  savePredictions = "final",
  classProbs = T,
  seed = as.list(rep(1,22)),
  summaryFunction = twoClassSummary
)

m = train(x = Xdata, 
          y = Ydata,
          method = "rf",
          tuneGrid=data.frame(mtry=1),
          trControl = fitControl)
m

m$finalModel
imp = data.frame(m$finalModel$importance)
colnames(imp) <- "Overall"



########
plot_df = data
plot_df$Sample = gsub("_.*", "", rownames(plot_df))
plot_df$Treatment = gsub("_[0-9]$", "", gsub("^[A-Z]*_", "", rownames(plot_df)))
plot_df$IsControl = ifelse(plot_df$Treatment == "CONTROL", "Control",
                           ifelse(plot_df$Sample == "HME", "OverExpressedB4",
                                  "siRNAB4"))

ggplot(data = plot_df) +
  geom_jitter(aes(x=B4Status, y=as.numeric(ADAMTS5), color=IsControl), width = 0.05) +
  geom_boxplot(aes(x=B4Status, y=as.numeric(ADAMTS5)), alpha=0.2) +
  facet_wrap(~Sample) +
  labs(x="B4 Fusion Status",
       y = "ADAMTS5 Expression",
       title="Most Important Feature: ADAMTS5")
  
ggplot(data = plot_df) +
  geom_jitter(aes(x=B4Status, y=as.numeric(RUNDC3B), color=IsControl), width = 0.05) +
  geom_boxplot(aes(x=B4Status, y=as.numeric(RUNDC3B)), alpha=0.2) +
  facet_wrap(~Sample) +
  labs(x="B4 Fusion Status",
       y = "RUNDC3B Expression",
       title="Least Important Feature: RUNDC3B")

tmp <- subset(data, select = c(RUNDC3B, B4Status))



####################
# Build model using RNA-Seq data from TCGA
#Look at all BCAR4 fusions
fusions = read.csv("data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=F)

# Find samples with fusions
bcar_fusions = subset(fusions, fusions$V7 %like% "BCAR4")
bcar_fusions$V19 = gsub("[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}$", "", bcar_fusions$V19)
fusion_samples = unique(bcar_fusions$V19)
fusion_samples.mod = gsub("-", ".", fusion_samples)
cohorts = unique(bcar_fusions$V20)

# Load RNA-Seq data for fusion samples and a random set of non-BCAR4 samples
# Units are log2(x+1) transformed RSEM normalized count
files = c()
for (cohort in cohorts) {
  if(cohort=="CRC"){
    cohort="COAD"
  }
  files = c(files, paste0("data/TCGA/pancancer_expression/TCGA-", cohort, ".HiSeqV2.txt"))
}

raw = read.csv(files[1], sep="\t", header=T)
fusion_rnaseq = data.frame(raw$sample)

set.seed(42)
for(f in files){
  message(f)
  raw = read.csv(f, sep="\t", header=T)
  rownames(raw) = raw$sample
  keep = c("sample", fusion_samples.mod)
  cols = colnames(raw) %in% keep
  n1 = sum(cols)
  message(paste0("Original keep ", sum(cols)))
  others_n = sum(cols) - 1
  kept = 0
  while (kept < others_n) {
    s = sample(colnames(raw), 1)
    if (raw["BCAR4",s] == 0) { # Only keep non-BCAR4 expressing samples
      keep = c(keep, s)
      kept = kept + 1
    }
  }
  cols2 = colnames(raw) %in% keep
  raw.sub = subset(raw, select = cols2)
  fusion_rnaseq = cbind(fusion_rnaseq, raw.sub)
}

repeated_cols = colnames(fusion_rnaseq) == "sample"
fusion_rnaseq = fusion_rnaseq[,!repeated_cols]         

# Subset down to the genes that were differentially expressed
# in our experimental data
de_exp = read.csv("data/Profiling/Microarray_results_All_genes.FromNicole.csv")
de_exp$FreqChange = ifelse(de_exp$TUH_B4siR268.status == "unchanged", 0, 1)
de_exp$FreqChange = ifelse(de_exp$TUH_B4siR271.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_ZCBOE.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_LBOE.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_LBMUT585.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
de_exp$FreqChange = ifelse(de_exp$HME_B4OE.status == "unchanged",
                           de_exp$FreqChange, de_exp$FreqChange + 1)
#de_exp$GeneName = as.factor(experiment_subset$GeneName)
de_exp = subset(de_exp, de_exp$FreqChange >= 5) # >= 5 is top 39

# Cuts down to top 22, as not all of them are present in RNASeq data
de_fusion = subset(fusion_rnaseq, fusion_rnaseq$raw.sample %in% de_exp$GeneName &
                     fusion_rnaseq$raw.sample != "BCAR4")

# Reformat for use with ML
rownames(de_fusion) <- de_fusion$raw.sample
t_fusion = data.frame(t(de_fusion[,-c(1)]))

t_fusion$Status = ifelse(rownames(t_fusion) %in% fusion_samples.mod, 1, 0)


# Build model
df = t_fusion

Xdata = df[,1:21]
df$Status = ifelse(df$Status == 0, "No", "Yes")
Ydata = as.factor(df[,22])

# Optional
#Xdata = subset(Xdata, select = c(rownames(imp)[order(imp$Overall)[1:50]]))


fitControl <- trainControl(
  method = "LOOCV",
  number = 1,
  savePredictions = "final",
  classProbs = T,
  seed = as.list(rep(1,87)),
  summaryFunction = twoClassSummary
)

m = train(x = Xdata, 
          y = Ydata,
          method = "rf",
          tuneGrid=data.frame(mtry=1),
          trControl = fitControl)

m$finalModel
m$finalModel$obsLevels


##Vis
mat = as.matrix(df[,1:21])
status_groups = ifelse(df$Status == "Yes", 2, 1)
colSide = brewer.pal(3,"Set1")[status_groups]
heatmap(mat, Colv = NA, RowSideColors=colSide)         
         