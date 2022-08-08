
library(caret)
library(data.table)
library(dplyr)
library(stats)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")

set.seed(42)

################
# Data import
microarray = read.csv("data/TCGA/microarray/AgilentG4502A_07_3_BRCA.tsv", sep="\t")
fusions = read.csv("data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=F)
fpkm = read.csv("data/TCGA/TCGA-BRCA.htseq_fpkm.roi.tsv", sep="\t")
experiments = read.csv("data/Profiling/Microarray_results_All_genes.FromNicole.csv")


########
# TCGA Data selection

# Find samples that have fusions 
bcar_fusion_samples = subset(fusions, fusions$V20=="BRCA" & fusions$V7 %like% "BCAR4")
bcar_fusion_samples$V19 = gsub("[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}$", "", bcar_fusion_samples$V19)

# Find samples with no BCAR4 expression from RNA-Seq
t_fpkm = t(fpkm)
colnames(t_fpkm) = c("HER2", "BCAR4") #c("ENSG00000141736.12", "ENSG00000262117.4")
t_fpkm = data.frame(t_fpkm[-c(1),])
t_fpkm$BCAR4 = as.numeric(t_fpkm$BCAR4)
t_fpkm$Samples = gsub("[A-Z]$", "", gsub("\\.", "-", rownames(t_fpkm)))
nobcar4 = subset(t_fpkm, t_fpkm$BCAR4 == 0)

# Get list of samples that have microarray data
micro_samples = gsub("\\.", "-", colnames(microarray))

bcar_fusion_samples = subset(bcar_fusion_samples, bcar_fusion_samples$V19 %in% micro_samples)
nobcar4 = subset(nobcar4, nobcar4$Samples %in% micro_samples)

samples_of_interest = c(bcar_fusion_samples$V19, nobcar4$Samples)
classification = c(rep("Fusion", 9), rep("No Fusion", 214))
tcga_samples = data.frame(Samples = samples_of_interest,
                          Classification = classification)
tcga_samples$SamplesMod = gsub("-", ".", tcga_samples$Samples)


#########
# Experimental data selection
# Keep only relevant columns
original_cols = colnames(experiments)
keep_cols = c("GeneName",
              original_cols[original_cols %like% "status" & (original_cols %like% "HME" | original_cols %like% "TUH")],
              original_cols[(original_cols %like% "_1" | original_cols %like% "_2") &
                              (original_cols %like% "HME" | original_cols %like% "TUH")])
experiment_subset = experiments[,keep_cols]

# Scale expression data (I think the same way that TCGA microarray data was scaled)
#exp_cols = c(keep_cols[keep_cols %like% "_1" | keep_cols %like% "_2"])
#exp_data = as.matrix(experiment_subset[,exp_cols])
#sc_ed = scale(exp_data, center = T, scale = T)
#sc_ed_d = data.frame(sc_ed)
#sc_ed_d$GeneName = as.factor(experiment_subset$GeneName)

# Count how often genes were differentially expressed
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

# Genes in microarray
micro_genes = microarray$sample
de_exp$InMicroarray = ifelse(de_exp$GeneName %in% micro_genes, 1, 0)

# Merge and 
# Select most frequently altered genes from experiments
# as possible candidate genes for ML model
de_scaled = merge(de_exp, sc_ed_d, by.x="GeneName", by.y="GeneName")
experiment_gene_candidates = subset(de_scaled, 
                                    de_scaled$FreqChange==4 & de_scaled$InMicroarray==1)

###########
# Gather data that will be used for model
# From experiments
cols = colnames(experiment_gene_candidates)
keep = c("GeneName",
         cols[cols %like% "_1" | cols %like% "_2"])
exper_exp = experiment_gene_candidates[,keep]
# From TCGA
tcga_micro_sub = subset(microarray, microarray$sample %in% exper_exp$GeneName)
set.seed(42)
keep = c("sample", 
         unique(subset(tcga_samples, tcga_samples$Classification == "Fusion")$SamplesMod),
         sample(unique(subset(tcga_samples, tcga_samples$Classification == "No Fusion")$SamplesMod), size=5))
tcga_exp = tcga_micro_sub[,keep]

tcga_exper = merge(x=exper_exp, y=tcga_exp,
                   by.x = "GeneName", by.y = "sample")


#########
# Reformat data
genes = tcga_exper$GeneName
t_t_e = data.frame(t(tcga_exper))
colnames(t_t_e) <- genes
t_t_e <-t_t_e[-c(1),]

#########
# Assign fusion status

# 0 = No BCAR4
# 1 = BCAR4 fusion
# 2 = BCAR4 overexpressed, but not fusion
# 3 = LBMUT585
status = c(0,0,2,2,3,3,1,1,1,1, # HME Controls, B4OE, LBMUT, LBOE, ZCBOE
           1,1,0,0,0,0, # TUH controls, si
           1,1,1,1, # TCGA samples with fusions
           0,0,0,0,0 # TCGA sapmles with no BCAR4 expression
           )
t_t_e$B4Status <- status

data <- subset(t_t_e, t_t_e$B4Status<=1)



########
# Check if additional scaling is needed

presc <- data %>%
  select_if(~ !any(is.na(.)))
presc$B4Status <- ifelse(presc$B4Status == 0, "No", "Yes")
presc$B4Status <- as.factor(presc$B4Status)

presc_exp <- data.frame(sapply(presc[,c(1:143)], as.numeric))
sc <- data.frame(scale(presc_exp))

sc_data <- sc
sc_data$B4Status <- presc$B4Status
rownames(sc_data) <- rownames(presc)
# None performed at this time. May be needed




#######
# Model creation
df <- sc_data %>%
  select_if(~ !any(is.na(.)))
#df <- data
#df$B4Status <- ifelse(df$B4Status == 0, "No", "Yes")
#df$B4Status <- as.factor(df$B4Status)

df_exp <- data.frame(sapply(df[,c(1:143)], as.numeric))
df_exp <- subset(df_exp, select = c(rownames(imp)[order(imp$Overall)[1:100]]))
df_exp <- df_exp[1:12,]
df_status <- df[,c("B4Status")]
df_status <- df[1:12,c("B4Status")]
fitControl <- trainControl(
  method = "LOOCV",
  number = 1,
  savePredictions = "final",
  classProbs = T,
  seed = as.list(rep(1,22)),
  summaryFunction = twoClassSummary
)

m = train(x = df_exp, 
      y = df_status,
      method = "rf",
      tuneGrid=data.frame(mtry=1),
      trControl = fitControl)

m
varImp(m)
imp <- varImp(m)
p <- m$pred

m$finalModel
fm <- m$finalModel$predicted
v = data.frame(m[["finalModel"]][["votes"]])
imp <- varImp(m$finalModel)
imp
o <- m$finalModel
o

fm
m
######
# General plotting
forplot
ggplot(data=df_exp[], aes(x=B4Status, y=as.numeric(TPRKB))) +
  geom_jitter(width=0.05) +
  geom_boxplot(alpha=0.2) +
  labs(x="BCAR4 Status",
       y="RRAGC Expression",
       title="Most important feature : RRAGC")

x = data.frame(df[,"RRAGC"])
rownames(x) <- rownames(df)


mede = median(scale(experiments$HME_B4OE_1))
medt = median(microarray$TCGA.C8.A1HL.01, na.rm = T)
ggplot() +
  geom_density(data=experiments, aes(x=scale(HME_B4OE_1) - mede), color="red") +
  geom_density(data=microarray, aes(x=TCGA.C8.A1HL.01), color="blue") +
  labs(x="Scaled value (mean-centered, divided by sd)",
       title="TCGA-C8-A1HL-01 (blue) vs HME_B4OE_1 (red)",
       caption="Both are B4+")



######
# Apply model to all microarray samples
m <- microarray



