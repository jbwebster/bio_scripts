
library(ggplot2) 
library(data.table)
library(dplyr)
library(forcats)
library(tidyr)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


## INITIAL DATA IMPORT
files <- list.files(path="data/TCGA/pancancer_expression", 
                    pattern = "*.txt", full.names=T)

# Expression values have units:  log2(x+1) transformed RSEM normalized count
f = files[1]
df <- data.frame()
for(f in files){
  raw <- read.csv(f, sep="\t", header=T)
  bcar4 <- subset(raw, raw$sample=="BCAR4")
  erbb2 <- subset(raw, raw$sample=="ERBB2")
  roi <- rbind(bcar4, erbb2)
  m <- reshape2::melt(roi)
  colnames(m) <- c("Gene", "Sample", "Expression")
  w <- m %>%
    spread(Gene, Expression)
  w$NumWithBCAR4 <- sum(w$BCAR4!=0)
  w$TotalCohort <- dim(w)[1]
  w$PropWithBCAR4 <- w$NumWithBCAR4 / w$TotalCohort
  cohort <- gsub(".HiSeqV2.txt","",gsub(".*TCGA-","",f))
  w$cohort <- cohort
  df <- rbind(df, w)
}
df$cohort <- as.factor(df$cohort)
saveRDS(df, "data/TCGA/pancancer_expression/BCAR4.ERBB2.allCohorts.rds")

cnv_files <- list.files(path="data/TCGA/pancancer_gistic2",
                        pattern = "*.txt", full.names=T)
f <- cnv_files[1]
cnv <- data.frame()
for(f in cnv_files) {
  raw <- read.csv(f, sep="\t", header=T)
  bcar4 <- subset(raw, raw$Gene.Symbol=="BCAR4")
  erbb2 <- subset(raw, raw$Gene.Symbol=="ERBB2")
  roi <- rbind(bcar4, erbb2)
  m <- reshape2::melt(roi)
  colnames(m) <- c("Gene", "Sample", "Gistic")
  w <- m %>%
    spread(Gene, Gistic)
  w$cohort <- gsub(".Gistic2.txt","",gsub(".*TCGA-","",f))
  cnv <- rbind(cnv, w)
}
cnv$cohort <- as.factor(cnv$cohort)
saveRDS(cnv, "data/TCGA/pancancer_gistic2/BCAR4.ERBB2.allCohorts.rds")



#####
## INITIAL DATA VIS
exp <- readRDS("data/TCGA/pancancer_expression/BCAR4.ERBB2.allCohorts.rds")
cnv <- readRDS("data/TCGA/pancancer_gistic2/BCAR4.ERBB2.allCohorts.rds")
colnames(exp) <- c("Sample", "BCAR4.exp", "ERBB2.exp", "BCAR4.Exp.Num", 
                   "CohortSize", "BCAR4.Exp.Prop", "Cohort")
colnames(cnv) <- c("Sample", "BCAR4.cnv", "ERBB2.cnv", "Cohort")
df <- merge(exp, cnv)
df$Sample <- gsub("\\.", "-", df$Sample)
df$Sample <- gsub("-[0-9A-Z]+$", "", df$Sample)
fusions <- read.csv("data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=F)
b4_fusions <- subset(fusions, V7 %like% "BCAR4")
b4_fusions$Patient <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*$", "", b4_fusions$V19)
df$hasBCAR4Fusion <- ifelse(df$Sample %in% b4_fusions$Patient, "hasFusion", "noFusion")

x <- exp %>%
  group_by(Cohort) %>%
  summarize(prop = mean(BCAR4.Exp.Prop), count=n())
xfac <- with(x, reorder(x$Cohort, prop, median, order=T))
x$Cohort <- factor(x$Cohort, levels=levels(xfac))
x$isBRCA <- ifelse(x$Cohort == "BRCA", "1", "0")
ggplot() +
  geom_col(data=x, aes(x=Cohort, y=prop, fill=isBRCA)) +
  geom_text(data=x, aes(x=Cohort, y=0.1, label = count)) +
  scale_x_discrete(limits = rev) +
  labs(y="% of Samples w/ BCAR4 Expression > 0",
       x="TCGA Cohort",
       title="Frequency of BCAR4 Expression") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot() +
  geom_point(data=df, aes(x=BCAR4.exp, y=BCAR4.cnv), alpha=0.2) +
  geom_smooth(data=df, aes(x=BCAR4.exp, y=BCAR4.cnv)) +
  labs(x="Expression (log2(normalized_count+1)",
       y="CNV (Gistic values)",
       title="Pan-Can BCAR4")

ggplot() +
  geom_point(data=df, aes(x=ERBB2.exp, y=ERBB2.cnv), alpha=0.2) +
  geom_smooth(data=df, aes(x=ERBB2.exp, y=ERBB2.cnv)) +
  labs(x="Expression (log2(normalized_count+1)",
       y="CNV (Gistic values)",
       title="Pan-Can ERBB2")

ggplot() +
  geom_jitter(data=df, aes(x=Cohort, y=BCAR4.exp), alpha=0.1) +
  geom_violin(data=df, aes(x=Cohort, y=BCAR4.exp), alpha=0.5) +
  labs(x="TCGA Cohort",
       y="BCAR4 Exp (log(normalized_count+1))",
       title="BCAR4 Expression by Cohort") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

x <- subset(df, BCAR4.exp>0)
x$Cohort = with(x, reorder(Cohort, BCAR4.exp, median))
y <- x %>%
  group_by(Cohort) %>%
  summarise(med=median(BCAR4.exp))
x$hasBCAR4Fusion <- ifelse(x$hasBCAR4Fusion=="hasFusion", "Yes", "No")
x$hasBCAR4Fusion <- factor(x$hasBCAR4Fusion, levels=c("Yes", "No"))
ggplot() +
  geom_jitter(data=x,
              aes(x=fct_rev(Cohort), y=BCAR4.exp, color=hasBCAR4Fusion), alpha=0.5) +
  geom_violin(data=x, aes(x=fct_rev(Cohort), y=BCAR4.exp), alpha=0.6) +
  labs(x="TCGA Cohort",
       y="BCAR4 Exp (log2(normalized_count+1))",
       title="BCAR4 Expression in TCGA",
       caption="Excluding samples with no BCAR4 expression. Ordered by median.",
       color="BCAR4 Fusion") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  
x <- df
x$ERBB2Amp <- ifelse(x$ERBB2.cnv > 1, "ERBB2_Amp", "No_ERBB2_Amp")
ggplot() +
  geom_jitter(data=x, aes(x=Cohort, y=ERBB2.exp, shape=ERBB2Amp), alpha=0.6) +
  geom_violin(data=x, aes(x=Cohort, y=ERBB2.exp), alpha=0.6) +
  scale_shape_manual(values = c(1,19)) +
  labs(x="TCGA Cohort",
       y="ERBB2 Exp (log(normalized_count+1))",
       title="HER2 Expression by Cohort") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  

ggplot() +
  geom_point(data=x, aes(x=BCAR4.exp, y=ERBB2.exp, 
                          color=Cohort, 
                          shape=ERBB2Amp)) +
  scale_shape_manual(values = c(1,19)) +
  labs(x="BCAR4 Expression",
       y="ERBB2 Expression",
       title="BCAR4 vs ERBB2 Pan-Can")

ggplot() +
  geom_point(data=x, aes(x=BCAR4.exp, y=ERBB2.exp,  
                         shape=ERBB2Amp), alpha=0.2) +
  scale_shape_manual(values = c(6,19)) +
  facet_wrap(~Cohort) +
  labs(x="BCAR4 Expression",
       y="ERBB2 Expression",
       title="BCAR4 vs ERBB2 Pan-Can")


# Criteria for elevated HER2 expression
ggplot() +
  geom_density(data=x, aes(ERBB2.exp)) +
  geom_vline(xintercept=14)

ggplot() + 
  geom_density(data=x,aes(BCAR4.exp))


### Mutual exclusivity attempt
exp <- readRDS("data/TCGA/pancancer_expression/BCAR4.ERBB2.allCohorts.rds")
cnv <- readRDS("data/TCGA/pancancer_gistic2/BCAR4.ERBB2.allCohorts.rds")
colnames(exp) <- c("Sample", "BCAR4.exp", "ERBB2.exp", "BCAR4.Exp.Num", 
                   "CohortSize", "BCAR4.Exp.Prop", "Cohort")
colnames(cnv) <- c("Sample", "BCAR4.cnv", "ERBB2.cnv", "Cohort")
df <- merge(exp, cnv)
df$Sample <- gsub("\\.", "-", df$Sample)
df$Sample <- gsub("-[0-9A-Z]+$", "", df$Sample)
fusions <- read.csv("data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=F)
b4_fusions <- subset(fusions, V7 %like% "BCAR4")
b4_fusions$Patient <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*$", "", b4_fusions$V19)
df$hasBCAR4Fusion <- ifelse(df$Sample %in% b4_fusions$Patient, "hasFusion", "noFusion")

df$HER2Amp <- ifelse(df$ERBB2.cnv > 1, 1, 0)
df$BCAR4Pos <- ifelse(df$BCAR4.exp > 2.5, 1, 0)
df$Both <- ifelse(df$HER2Amp ==1 & df$BCAR4Pos==1, 1, 0)
df$Neither <- ifelse(df$HER2Amp==0 & df$BCAR4Pos==0, 1, 0)
x <- df[!duplicated(df$Sample),]
neither <- sum(x$Neither)
anotb <- sum(x$BCAR4Pos) - sum(x$Both)
bnota <- sum(x$HER2Amp)-sum(x$Both)
ab <- sum(x$Both)

#neither <- 1000
#anotb <- 53
#bnota <- 128
#ab <- 10

counts <- matrix(
  c(ab,
    anotb,
    bnota,
    neither), nrow=2
)
result <- fisher.test(counts, alternative="less")
result
result.df <- data.frame(Cohort=c("PanCan"), fisher.p=c(round(result$p.value, 5)))
counts

for(cohort in unique(df$Cohort)){
  message(cohort)
  x <- subset(df, Cohort==cohort)
  neither <- sum(x$Neither)
  anotb <- sum(x$BCAR4Pos) - sum(x$Both)
  bnota <- sum(x$HER2Amp)-sum(x$Both)
  ab <- sum(x$Both)
  counts <- matrix(
    c(neither,
      anotb,
      bnota,
      ab), nrow=2
  )
  result <- fisher.test(counts, alternative="less")
  message(result$p.value)
  tmp <- data.frame(Cohort=c(cohort), fisher.p=c(round(result$p.value, 5)))
  result.df <- rbind(result.df, tmp)
}

xxx <- df %>%
  group_by(Cohort) %>%
  summarise(size=n(), neither=sum(Neither), 
            anotb=sum(BCAR4Pos)-sum(Both),
            bnota=sum(HER2Amp)-sum(Both),
            ab=sum(Both))
xx <- data.frame(Cohort="PanCan", size=nrow(df),
                 neither=sum(df$Neither),
                 anotb=sum(df$BCAR4Pos)-sum(df$Both),
                 bnota=sum(df$HER2Amp)-sum(df$Both),
                 ab=sum(df$Both))
xxx <- rbind(xxx, xx)
xxy <- merge(result.df, xxx)

### ONCOPRINT INPUT FILE CREATION
# PANCANCER
tcgaGenomicOncoprint <- function(df, cohort) {
  genomic <- data.frame(Patient=character(),
                        Gene=character(),
                        Alteration=character(),
                        Type=character(),
                        Track=character())
  if(cohort != "PANCAN"){
    df <- subset(df, Cohort==cohort)
  }
  for(i in 1:nrow(df)) {
    addedRow=FALSE
    patient <- gsub("\\.","-",df[i,"Sample"])
    if (df[i,"ERBB2.cnv"]>1){
      genomic[nrow(genomic)+1,]=c(patient, "ERBB2", "AMP", "CNA", "ERBB2_Amp")
      addedRow=TRUE
    }
    if (df[i,"ERBB2.exp"]>14){
      genomic[nrow(genomic)+1,]=c(patient, "ERBB2", "HIGH", "EXP", "ERBB2_High(>14)")
      addedRow=TRUE
    }
    if (df[i,"BCAR4.exp"]>0){
      genomic[nrow(genomic)+1,]=c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Expressed")
      addedRow=TRUE
    }
    if (df[i,"BCAR4.exp"]>3){
      genomic[nrow(genomic)+1,]=c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_High(>3)")
      addedRow=TRUE
    }
    if (!addedRow) {
      genomic[nrow(genomic) + 1,] = c(patient, "", "", "", "")
    }
  }
  write.table(genomic, file=paste0("data/TCGA/pancancer_oncoprints/TCGA-",cohort,".tsv"),
              row.names=FALSE, col.names=FALSE, sep="\t",
              quote=FALSE)
}

tcgaGenomicOncoprint(df, "PANCAN")
cohorts <- unique(as.character((df$Cohort)))
for(cohort in cohorts){
  print(cohort)
  tcgaGenomicOncoprint(df, cohort)
}
##

## Vis of oncoprint results
rm(list=ls())
onc <- read.csv("data/TCGA/pancancer_oncoprints/SummaryStats.csv", header=T)
s <- onc[,c("A","B","Cohort","NeitherPercent","ANotBPercent","BNotAPercent","BothPercent")]


onc_long <- melt(setDT(s), id.vars=c("A","B","Cohort"), variable.name="Value")


f <- subset(onc_long, ((A=="BCAR4_High(>3)" & B=="ERBB2_Amp") |
              (A=="ERBB2_Amp" & B=="BCAR4_High(>3)")) &
              (Value=="NeitherPercent"))
lev <- unique(c("PanCan",f$Cohort[unique(order(f$value))]))
onc_long$Cohort <- factor(onc_long$Cohort, 
                          levels=lev,
                          ordered=TRUE)
x <- subset(onc_long, (A=="BCAR4_High(>3)" & B=="ERBB2_Amp") |
              (A=="ERBB2_Amp" & B=="BCAR4_High(>3)"))
## Define A to be BCAR4 high and B to be ERBB2 amp

x$UpdatedValue <- ifelse(x$Value=="NeitherPercent", "NeitherPercent",
                         ifelse(x$Value=="BothPercent", "BothPercent",
                          ifelse(x$A=="ERBB2_Amp" & x$Value=="ANotBPercent", "BNotAPercent",
                            ifelse(x$A=="ERBB2_Amp" & x$Value=="BNotAPercent", "ANotBPercent",
                              ifelse(x$Value=="ANotBPercent", "ANotBPercent","BNotAPercent")))))
x$UpdatedValue <- factor(x$UpdatedValue, 
                         levels=c("NeitherPercent", "BothPercent", 
                                  "ANotBPercent", "BNotAPercent"))

### Need to modify this to confirm that A/B split is correct, since
# some samples flip which is which
ggplot() +
  geom_col(data = subset(x, UpdatedValue!="NeitherPercent"),
           aes(x=Cohort, y=100*value, fill=UpdatedValue)) +
  scale_fill_hue(labels = c("Both", "BCAR4High", "HER2Amp")) +
  scale_y_continuous(limits=c(0,25), expand=c(0,0)) +
  labs(x="TCGA Cohort",
       y="Percent of Cohort",
       title="HER2Amp and BCAR4High TCGA Comparison",
       caption="BCARHigh: log2(normalized_count+1) > 3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


y <- subset(onc_long, (A=="BCAR4_Expressed" & B=="ERBB2_Amp") |
              (A=="ERBB2_Amp" & B=="BCAR4_Expressed"))
## Define A to be BCAR4 high and B to be ERBB2 amp

y$UpdatedValue <- ifelse(y$Value=="NeitherPercent", "NeitherPercent",
                         ifelse(y$Value=="BothPercent", "BothPercent",
                                ifelse(y$A=="ERBB2_Amp" & x$Value=="ANotBPercent", "BNotAPercent",
                                       ifelse(y$A=="ERBB2_Amp" & x$Value=="BNotAPercent", "ANotBPercent",
                                              ifelse(y$Value=="ANotBPercent", "ANotBPercent","BNotAPercent")))))
y$UpdatedValue <- factor(y$UpdatedValue, 
                         levels=c("NeitherPercent", "BothPercent", 
                                  "ANotBPercent", "BNotAPercent"))

f <- subset(y, UpdatedValue=="NeitherPercent")
lev <- unique(c("PanCan",as.character(f$Cohort[unique(order(f$value))])))
y$Cohort <- factor(y$Cohort, 
                          levels=lev,
                          ordered=TRUE)

ggplot() +
  geom_col(data = subset(y, UpdatedValue!="NeitherPercent"),
           aes(x=Cohort, y=100*value, fill=UpdatedValue)) +
  scale_fill_hue(labels = c("Both", "BCAR4Expressed", "HER2Amp")) +
  scale_y_continuous(limits=c(0,60), expand=c(0,0)) +
  labs(x="TCGA Cohort",
       y="Percent of Cohort",
       title="HER2Amp and BCAR4Expressed TCGA Comparison") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## STAT TEST
xx <- data.frame(
  A=c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
  B=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)
)
yy <- data.frame(
  A=c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
  B=c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)
)
zz <- data.frame(
  A=c(rep(0,49), rep(1,11), rep(0,10)),
  B=c(rep(0,49), rep(0,11), rep(1,10))
)
f1 = fisher.test(xx$A, xx$B, alternative="g") 
f2 = fisher.test(yy$A, yy$B)
f3 = fisher.test(zz$A, zz$B, alternative="less")
f1$p.value
f2$p.value
f3$p.value

