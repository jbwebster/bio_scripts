
library(genefu)
library(breastCancerTRANSBIG) # For annotation info
library(ggplot2)
library(dplyr)




setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


## CREATE PAM50 SUBTYPE PREDICTIONS
# Load PAM50 centroid scaling data
#data(pam50.scale)
#data(pam50.robust)
data(pam50)
#pam50.robust <- pam50.scale # Trick genefu into using pam50.scale info, rather than robust
pam50.robust <- pam50
#data(pam50.robust)
#rm(pam50.scale)
rm(pam50)
# Load annotation data
dd <- get(data(list="transbig"))
dannot <- featureData(dd)@data

files <- list.files(path="data/TCGA/pancancer_expression",
                    pattern = "*.txt", full.names=T)
df <- data.frame()
for(f in files){
  cohort <- gsub(".HiSeqV2.txt","",gsub(".*TCGA-","",f))
  message(cohort)
  raw <- read.csv(f, sep="\t", header=T)
  # Two genes use alternate aliases in our data, update them accordingly
  raw$sample <- ifelse(raw$sample=="NUF2", "CDCA1",
                       ifelse(raw$sample=="NDC80", "KNTC2", raw$sample))
  rownames(raw) <- raw$sample
  t.raw <- t(raw[,-1])
  pam50_genes <- pam50.robust$centroids.map$probe
  pam50_genes_only <- t.raw[,colnames(t.raw) %in% pam50_genes]
  #pam50_genes_only <- apply(pam50_genes_only, 2, function(y)  (y - mean(y)))
  #pam50_genes_only <- apply(pam50_genes_only, 2, function(y) ((y - mean(y)) / sd(y)))
  
  PAM50Preds<-molecular.subtyping(sbt.model="pam50",data=pam50_genes_only,
                                  annot=dannot)
  cohort.predictions <- data.frame(PAM50Preds[["subtype.proba"]])
  cohort.predictions$SubtypeCall <- PAM50Preds[["subtype"]]
  cohort.predictions$Cohort <- cohort
  df <- rbind(df, cohort.predictions)
}

#saveRDS(df, "data/TCGA/pancancer.pam50.scaled.rds")
#saveRDS(df, "data/TCGA/pancancer.pam50.robust.rds")
saveRDS(df, "data/TCGA/pancancer.pam50.none.rds") # Using this one for now


##########################

##########################
n <- readRDS("data/TCGA/pancancer.pam50.none.rds")
df <- n
df$Cohort <- factor(df$Cohort)
df$SubtypeCall <- factor(df$SubtypeCall)

# Visualize PAM50 results
df$Count <- 1
x <- df %>%
  group_by(Cohort, SubtypeCall) %>%
  summarise(Count=n())
xx <- df %>%
  group_by(Cohort) %>%
  summarise(CohortSize=sum(Count))
yy <- merge(x,xx)
yy$SubtypeProportion <- yy$Count / yy$CohortSize

ggplot() +
  geom_col(data=yy, aes(x=Cohort, y=SubtypeProportion, fill=SubtypeCall)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  labs(x="TCGA Cohort",
       y="Frequency of Subtype",
       title="PanCancer PAM50 Subtyping",
       caption="genefu results. Scaling='none'. log2(RMSE+1)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))








## COMPARISON

# library(TCGAbiolinks)
# dataSubt <- TCGAquery_subtype(tumor = "BRCA")
# n <- cohort.predictions
# #n <- subset(readRDS("data/TCGA/pancancer.pam50.scaled.rds"), Cohort=="BRCA")
# n$Patient <- gsub("-[0-9A-Z]+$", "", gsub("\\.", "-", rownames(n)))
# n = n[!duplicated(n$Patient),]
# n <- subset(n, n$Patient %in% dataSubt$patient)
# dataSubt <- subset(dataSubt, dataSubt$patient %in% n$Patient)
# table(dataSubt$BRCA_Subtype_PAM50)
# table(n$SubtypeCall)
# 
# 
# x <- read.csv("data/TCGA/x.txt", header=T, sep="\t")
# x <- subset(x, x$Tumor.Type=="BRCA")
# x$Sample <- gsub("-[0-9A-Z]+$", "", x$Sample)
# x <- x[!duplicated(x$Sample),]
# x <- subset(x, x$Sample %in% dataSubt$patient)
# table(x$Subtype)
# 
# 
# x <- read.csv("data/TCGA/PAM50_nongenefu/BRCA/BRCA_pam50scores.txt", sep="\t")
# x$X <- gsub("-[0-9A-Z]+$", "", gsub("\\.", "-", x$X))
# x <- subset(x, x$X %in% dataSubt$patient)
# x <- x[!duplicated(x$X),]
# table(x$Call)
# 
# 
# 
# 
# x <- read.csv('scripts/PAM50/PAM50_R/mediansPerDataset_v2.txt', header=T,sep="\t")
# cen <-read.csv('scripts/PAM50/PAM50_R/pam50_centroids.txt', header=T, sep="\t")


xxx <- cohort.predictions
xxx$Patient <- gsub("-[0-9A-Z]{2}$", "", gsub("\\.", "-", rownames(xxx)))
xxx <- xxx[!duplicated(xxx$Patient),]
table(xxx$SubtypeCall)

yyy <- data.frame(pam50_genes_only)
mean(yyy$ERBB2)





### OTHER ATTEMPT... ###
# https://www.cell.com/cell/pdfExtended/S0092-8674(15)01195-2
# 1) downsampling to a similar ER distribution as the original training set, 
# 2) median center all samples based on the medians found in the downsampled cohort, 
# 3) apply Packer code as originally described
# Original training set from Packer:
#   ER Positive       ER Negative
#     114               77
# Ratio: 114/77 = 1.480519
brca_pheno <- read.csv("data/TCGA/TCGA.BRCA.sampleMap-BRCA_clinicalMatrix", sep="\t")
brca_pheno <- subset(brca_pheno)
table(brca_pheno$ER_Status_nature2012)
# ER Positive: 601, ER Negative: 179
# Keeping ER Negative at 179, we want only
# 179 * 1.480519 = 265
set.seed(42)
neg <- subset(brca_pheno, ER_Status_nature2012=="Negative") # 179 samples
all.pos <- subset(brca_pheno, ER_Status_nature2012=="Positive") # 601 samples
all.pos <- all.pos[!duplicated(all.pos$X_PATIENT),] # Still 601
pos <- sample(all.pos$X_PATIENT, size=265)
sub.pos <- subset(all.pos, all.pos$X_PATIENT %in% pos)
training.set <- rbind(neg, all.pos)


cohort <- "BRCA"
f=files[3]
data(pam50.robust)
raw <- read.csv(f, sep="\t", header=T)
# Two genes use alternate aliases in our data, update them accordingly
raw$sample <- ifelse(raw$sample=="NUF2", "CDCA1",
                     ifelse(raw$sample=="NDC80", "KNTC2", raw$sample))
rownames(raw) <- raw$sample
t.raw <- t(raw[,-1])
pam50_genes <- pam50.robust$centroids.map$probe
pam50_genes_only <- t.raw[,colnames(t.raw) %in% pam50_genes]
pam50_genes_only <- data.frame(pam50_genes_only)
pam50_genes_only$PATIENT <- gsub("-[0-9A-Z]{2}$", "", gsub("\\.", "-", rownames(pam50_genes_only)))
training.samples <- subset(pam50_genes_only, pam50_genes_only$PATIENT %in% training.set$X_PATIENT)
training.samples <- training.samples[!duplicated(training.samples$PATIENT),]
training.samples <- subset(training.samples, select=-c(PATIENT))
medians <- apply(training.samples, 2, function(y) (median(y))) # USE THESE TO CENTER



raw <- read.csv(f, sep="\t", header=T)
# Two genes use alternate aliases in our data, update them accordingly
raw$sample <- ifelse(raw$sample=="NUF2", "CDCA1",
                     ifelse(raw$sample=="NDC80", "KNTC2", raw$sample))
rownames(raw) <- raw$sample
t.raw <- t(raw[,-1])
pam50_genes <- pam50.robust$centroids.map$probe
pam50_genes_only <- t.raw[,colnames(t.raw) %in% pam50_genes]
pam50_genes_only <- data.frame(pam50_genes_only)
pam50_genes_only$PATIENT <- gsub("-[0-9A-Z]{2}$", "", gsub("\\.", "-", rownames(pam50_genes_only)))
pam50_genes_only <- pam50_genes_only[!duplicated(pam50_genes_only$PATIENT),]
pam50_genes_only <- subset(pam50_genes_only, select=-c(PATIENT))
if(sum(names(medians) == colnames(pam50_genes_only))==50){
  message("Names match")
}
xx <- pam50_genes_only
for(i in 1:ncol(xx)){
  coln <- colnames(xx)[i]
  xx[,coln] <- ((xx[,coln] - median(xx[,coln])) + medians[coln])
}
data(pam50)
#data(pam50.robust)
pam50.robust <- pam50

PAM50Preds<-molecular.subtyping(sbt.model="pam50",data=xx,
                                annot=dannot)
cohort.predictions <- data.frame(PAM50Preds[["subtype.proba"]])
cohort.predictions$SubtypeCall <- PAM50Preds[["subtype"]]
cohort.predictions$Cohort <- cohort

df <- data.frame(t(xx))
#dir.create(paste0("data/TCGA/PAM50_nongenefu/",cohort))
write.table(df, file=paste0("data/TCGA/PAM50_nongenefu/",cohort,"/",cohort,".txt"),
            quote=F, sep="\t", col.names=NA)

library(TCGAbiolinks)
dataSubt <- TCGAquery_subtype(tumor = "BRCA")
n <- res
n <- cohort.predictions
#n <- subset(readRDS("data/TCGA/pancancer.pam50.scaled.rds"), Cohort=="BRCA")
n$Patient <- gsub("-[0-9A-Z]+$", "", gsub("\\.", "-", rownames(n)))
n$Patient <- gsub("-[0-9A-Z]+$", "", gsub("\\.", "-", n$X))
n = n[!duplicated(n$Patient),]
n <- subset(n, n$Patient %in% dataSubt$patient)
dataSubt <- subset(dataSubt, dataSubt$patient %in% n$Patient)

table(n$Call)
table(cohort.predictions$Cohort)

p <- n[,c("Patient", "SubtypeCall")]
colnames(p) <- c("patient", "SubtypeCall")
tt <- dataSubt[,c("patient", "BRCA_Subtype_PAM50")]
m <- merge(p, tt)
ggplot(m) +
  geom_jitter(aes(x=SubtypeCall, y=BRCA_Subtype_PAM50), width=0.2, height=0.2)
  
  
  