

library(data.table)
library(survminer)
library(survival)
library(gridExtra)
library(ggrepel)
library(dplyr)

library(genefu)



setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


pam50.calls <- readRDS("data/TCGA/pancancer.pam50.none.rds")
pam50.calls$SubtypeCall <- as.character(pam50.calls$SubtypeCall)
pam50.calls$SubtypeCall <- ifelse(pam50.calls$Cohort=="BRCA", pam50.calls$SubtypeCall, NA)
exp <- readRDS("data/TCGA/pancancer_expression/BCAR4.ERBB2.allCohorts.rds")
cnv <-  readRDS("data/TCGA/pancancer_gistic2/BCAR4.ERBB2.allCohorts.rds")
colnames(exp) <- c("Sample", "BCAR4.exp", "ERBB2.exp", "BCAR4.Exp.Num", 
                   "CohortSize", "BCAR4.Exp.Prop", "Cohort")
colnames(cnv) <- c("Sample", "BCAR4.cnv", "ERBB2.cnv", "Cohort")
mol <- merge(exp, cnv)
mol$Sample <- gsub("\\.", "-", mol$Sample)

files <- list.files(path="data/TCGA/pancancer_phenotype",
                    pattern = "*.txt", full.names=T)
clinical <- data.frame()
for(f in files){
  cohort <- gsub("_survival.txt","",gsub(".*phenotype/","",f))
  message(cohort)
  raw <- read.csv(f, sep="\t", header=T)
  raw$cohort <- cohort
  colnames(raw) <- c("Sample", "Patient", "OS", "OS.time",
                     "DSS", "DSS.time", "DFI", "DFI.time",
                     "PFI", "PFI.time", "Redaction", "Cohort")
  raw <- raw[!duplicated(raw$Patient),]
  clinical <- rbind(clinical, raw)
}
pam50.calls$Sample <- gsub("\\.", "-", rownames(pam50.calls))
pam50dt <- data.table(pam50.calls)
moldt <- data.table(mol)

dt <- merge(pam50dt, moldt, by = 'Sample', all.x = TRUE)
clinicaldt <- data.table(clinical)
dt <- merge(dt, clinicaldt, by='Sample', all.x = TRUE)
dt <- dt[!duplicated(dt$Patient),]
saveRDS(dt, "data/TCGA/pancan.exp.cnv.pheno.pam50.rds")

###########
dt <- readRDS("data/TCGA/pancan.exp.cnv.pheno.pam50.rds")

dt$BCAR4.high <- ifelse(dt$BCAR4.exp>2, 1, 0)
dt$ERBB2.amp <- ifelse(dt$ERBB2.cnv>1, 1, 0)

performSurvivalComparisons = function(cohort, df, group1, group2) {
  b = subset(df, df$Cohort==cohort)
  b$Group = ifelse(b$ERBB2.amp==1 & b$BCAR4.high==1, "Both",
                    ifelse(b$ERBB2.amp==1 & b$BCAR4.high==0, "ERBB2.Amp",
                           ifelse(b$ERBB2.amp==0 & b$BCAR4.high==1, "BCAR4.High", "Neither")))
  bsub = subset(b, b$Group==group1 | b$Group==group2)
  if(group1 %in% bsub$Group & group2 %in% bsub$Group) {
    f = survfit(Surv(OS.time / 365, OS) ~ Group, data = bsub)
    p = surv_pvalue(f, data=bsub)
    return(as.numeric(p$pval))
  }
  return(NA)
}

performSurvivalComparisons2 = function(cohort, df, group1, group2) {
  b = subset(df, df$Cohort==cohort)
  b$Group = ifelse(b$BCAR4.high==1, 3,
                   ifelse(b$BCAR4.exp>0, 2, 1))
  bsub = subset(b, b$Group>=group1 | b$Group==group2)
  bsub$Group = ifelse(bsub$Group>=group1, "A", "B")
  if("A" %in% bsub$Group & "B" %in% bsub$Group) {
    f = survfit(Surv(OS.time / 365, OS) ~ Group, data = bsub)
    p = surv_pvalue(f, data=bsub)
    return(as.numeric(p$pval))
  }
  return(NA)
}

df <- dt
cohort="BRCA"
group1="BCAR4.High"
group2="ERBB2.Amp"

survival.pvals <- data.frame()

for(cohort in unique(dt$Cohort)) {
  message(cohort)
  a=performSurvivalComparisons(cohort=cohort, df=dt, group1="BCAR4.High", group2="ERBB2.Amp")
  b=performSurvivalComparisons(cohort=cohort, df=dt, group1="BCAR4.High", group2="Neither")
  cc=performSurvivalComparisons(cohort, dt, "Neither", "Both")
  d=performSurvivalComparisons(cohort, dt, "ERBB2.Amp", "Both")
  e=performSurvivalComparisons(cohort, dt, "ERBB2.Amp", "Neither")
  f=performSurvivalComparisons2(cohort, dt, 3, 1)
  g=performSurvivalComparisons2(cohort, dt, 2, 1)
  n1 <- nrow(subset(dt, dt$Cohort==cohort & dt$BCAR4.high==1))
  n2 <- nrow(subset(dt, dt$Cohort==cohort & dt$ERBB2.amp==1))
  n3 <- nrow(subset(dt, dt$Cohort==cohort))
  x <- data.frame(cohort=c(cohort),
                  bcar4high.erbb2amp=c(a),
                  bcar4high.neither=c(b),
                  neither.both=c(cc),
                  erbb2amp.both=c(d),
                  erbb2amp.neither=c(e),
                  bcar4high.bcar4none=c(f),
                  bcar4exp.bcar4none=c(g),
                  cohort.size=c(n3),
                  bcar4high.n=c(n1),
                  erbb2amp.n=c(n2))
  survival.pvals <- rbind(survival.pvals, x)
}

survival.pvals <- subset(survival.pvals, survival.pvals$cohort!="BRCA")
write.csv(survival.pvals, "data/TCGA/PancanSurvivalPvals_Bcar4HighGE2.csv")

### Plotting the significance of survival curves and some standout examples

OSSignificancePlot <- function(df, group, groupa, groupb) {
  df <- df[,c("cohort", group, "bcar4high.n", "erbb2amp.n", "cohort.size")]
  colnames(df) <- c("cohort", "pvals", "bcar4high.n", "erbb2amp.n", "cohort.size")
  p <- ggplot(data=df, aes(x=bcar4high.n/cohort.size, 
                                  y=erbb2amp.n/cohort.size)) + 
    geom_point(aes(size=(1-pvals)**4, color=cohort.size)) +
    geom_label_repel(data=subset(df,
                                 (bcar4high.n/cohort.size) > 0.02 &
                                   (erbb2amp.n/cohort.size) > 0.02),
                     label.size=NA,
                     size=7,
                     aes(label=paste0(cohort, ", p=", round(pvals, 3)))) +
    # geom_label_repel(data=subset(df,
    #                             (bcar4high.n/cohort.size) > 0.02 &
    #                               (erbb2amp.n/cohort.size) > 0.02 &
    #                               pvals<0.05),
    #                 label.size=NA,
    #                 size=7,
    #                 color="red",
    #                 aes(label=paste0(cohort, ", p=", round(pvals, 3)))) +
    theme_classic() +
    theme(plot.title = element_text(size=22)) +
    guides(size="none") +
    labs(x="% of cohort with BCAR4+",
         y="% of cohort with HER2 Amp",
         title=paste0("OS significance per cohort when comparing ", groupa, " vs ", groupb),
         caption="BCAR4+ : log2(rsme+1) >= 2") 
    
  return(p)
}

## BCAR4+ vs ERBB2-Amp
p1 <- OSSignificancePlot(survival.pvals, "bcar4high.erbb2amp", "BCAR4+", "HER2-Amp")
p2 <- OSSignificancePlot(survival.pvals, "bcar4high.neither", "BCAR4+", "Neither")
p3 <- OSSignificancePlot(survival.pvals, "neither.both", "Neither", "Both")
p4 <- OSSignificancePlot(survival.pvals, "erbb2amp.both", "HER2-Amp", "Both")
p5 <- OSSignificancePlot(survival.pvals, "erbb2amp.neither", "HER2-Amp", "Neither")

#grid.arrange(p1, p2, p3, p4, p5)
p1
p2
p3
p4
p5


ggplot(data=survival.pvals, aes(x=bcar4high.n/cohort.size, 
                                y=erbb2amp.n/cohort.size)) + 
  geom_point(aes(size=(1-bcar4high.erbb2amp)**4, color=cohort.size)) +
  geom_label_repel(data=subset(survival.pvals,
                               (bcar4high.n/cohort.size) > 0.02 &
                                (erbb2amp.n/cohort.size) > 0.02),
                   label.size=NA,
                   aes(label=paste0(cohort, ", p=", round(bcar4high.erbb2amp, 3)))) +
  theme_classic() +
  guides(size="none") +
  labs(x="% of cohort with BCAR4+",
       y="% of cohort with HER2 Amp",
       title="OS significance per cohort when comparing BCAR4+ vs ERBB2-Amp",
       caption="BCAR4+ : log2(rsme+1) >= 2")

#

### TESTING
# Plot BRCA just to make sure it looks as expected
cohort="BRCA"
b <- dt
b <- subset(dt, Cohort==cohort)
b$Group <- ifelse(b$BCAR4.exp > 0  & b$ERBB2.amp==1, "Both",
                  ifelse(b$BCAR4.exp> 0 & b$ERBB2.amp==0, "BCAR4+ Only",
                         ifelse(b$BCAR4.exp==0 & b$ERBB2.amp==1, "HER2-Amp Only", "Neither")))
b$Group <- ifelse(b$BCAR4.exp > 0.5, "BCAR4+",  
                  ifelse(b$BCAR4.exp == 0, "No B4 Expression", "Intermediate Expression"))

group1="BCAR4+"
group2="No B4 Expression"
bsub = subset(b, b$Group==group1 | b$Group==group2)
bsub = subset(bsub, bsub$SubtypeCall=="LumA")
bsub = b
bsub$Group <- factor(bsub$Group, levels=c("Both", "BCAR4+ Only", "HER2-Amp Only", "Neither"))
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = bsub)
ggsurvplot(f, data = bsub, pval=T,
           risk.table=T) +
  labs(title="PanCan")
p <- surv_pvalue(f)
p$pval.txt

###

old_master <- readRDS(file="data/TCGA/master.rds")
library(TCGAbiolinks)
dataSubt <- TCGAquery_subtype(tumor = "BRCA")
mm <- merge(old_master, dt, by="Patient")
dataSubt <- dataSubt[,c("patient", "BRCA_Subtype_PAM50", "days_to_death", "days_to_last_followup", "vital_status")]
dataSubt$days_to_death <- ifelse(dataSubt$days_to_death!="NA", dataSubt$days_to_death, dataSubt$days_to_last_followup)
colnames(dataSubt) <- c("Patient", "BRCA_Subtype_PAM50", "days_to_death", "days_to_followup", "vital_status")
mm <- merge(mm, dataSubt)
mm$TripleNegative <- ifelse(mm$IHC.ER=="Negative" &
                              mm$IHC.HER2=="Negative" &
                              mm$IHC.PR=="Negative", 1, 0)
mm$Group <- ifelse(mm$BCAR4.exp>0 & mm$IHC.HER2=="Positive", "Both",
                             ifelse(mm$BCAR4.exp>0 & mm$IHC.HER2!="Positive", "BCAR4+ Only",
                                    ifelse(mm$BCAR4.exp==0 & mm$IHC.HER2=="Positive", "HER2+ Only", "Neither")))
mm$Group <- ifelse(mm$BCAR4.exp>0 & mm$PAM50.cbioportal=="Her2", "Both",
                   ifelse(mm$BCAR4.exp>0 & mm$PAM50.cbioportal!="Her2", "BCAR4+ Only",
                          ifelse(mm$BCAR4.exp==0 & mm$PAM50.cbioportal=="Her2", "HER2+ Only", "Neither")))
mm$Group <- ifelse(mm$BCAR4.exp==0, "No Expression", 
                   ifelse(mm$BCAR4.exp>2, "High",
                          ifelse(mm$BCAR4.exp>0.75, "Intermediate", "Low")))
mm$Group <- ifelse(mm$BCAR4.exp==0, "No Expression", #  "Expressed")
                   ifelse(mm$BCAR4.exp>2, "High", "Low"))
mm$Group <- ifelse(mm$BCAR4.exp == 0, "No Expression", "Expressed")
mm$Group <- ifelse(mm$BCAR4.fpkm==0, "No Expression", "Expressed")
mm$Group <- ifelse(mm$BCAR4.exp == 0, "No Expresion", 
                   ifelse(mm$BCAR4.exp>0.5767, "Top 25%", "Other"))
nrow(subset(mm, mm$PAM50.cbioportal=="LumA"))
mn <- subset(mm, mm$PAM50.cbioportal=="LumA")
mn <- mn[order(mn$BCAR4.exp)[1:90],] # Top 20%
mo <- subset(mm, mm$PAM50.cbioportal=="LumA" & mm$BCAR4.exp==0)
mn$Group = "Top 20%"
mo$Group = "No Expression"
mm <- rbind(mn, mo)
nrow(subset(mm, mm$SubtypeCall=="LumA"))
mn <- subset(mm, mm$SubtypeCall=="LumA")
mn <- mn[order(mn$BCAR4.exp)[1:121],] # Top 20%
mo <- subset(mm, mm$SubtypeCall=="LumA" & mm$BCAR4.exp==0)
mn$Group = "Top 20%"
mo$Group = "No Expression"
mm <- rbind(mn, mo)

mm$Group <- factor(mm$Group, levels=c("Both", "BCAR4+ Only", "HER2+ Only", "Neither"))
mm <- subset(mm, mm$PAM50.cbioportal=="LumA")
mm <- subset(mm, mm$Group!="Other")
mm <- subset(mm, mm$BRCA_Subtype_PAM50=="LumB" | mm$BRCA_Subtype_PAM50=="LumA")
f <- survfit(Surv(OS.time / 365, OS) ~ Group, data = mm)
ggsurvplot(f, data = mm, pval=T,
           risk.table=T) +
  labs(title="BCAR4 Expression in BRCA LumA/B",
       x="Years")
mm$days_to_death <- as.numeric(mm$days_to_death)
mm$vital_status <- ifelse(mm$vital_status=="Alive", 0, 1)
f <- survfit(Surv(days_to_death, vital_status) ~ Group, data = mm)
ggsurvplot(f, data = mm, pval=T,
           risk.table=T) +
  labs(title="BCAR4 Expression in BRCA LumA/B",
       x="Years")

x <- subset(mm, mm$PAM50.cbioportal=="Her2")
nrow(subset(x, x$BCAR4.exp>2))

BR##
b <- dt
b$Group <- ifelse(b$BCAR4.exp > 2, "High", 
                  ifelse(b$BCAR4.exp == 0, "None", "Low"))
ggplot() +
  geom_density(data=b, aes(x=ERBB2.exp, fill=factor(ERBB2.amp)),
               alpha=0.2)
ggplot() +
  geom_density(data=b, aes(x=ERBB2.exp)) +
  geom_vline(xintercept=14) +
  geom_vline(xintercept=9) +
  labs(x="ERRB2 Expression",
       title="ERBB2 Grouping")

ggplot() +
  geom_density(data=subset(mm, Cohort=="BRCA"), aes(x=BCAR4.exp, group=PAM50.cbioportal)) +
  geom_vline(xintercept=2) +
  geom_vline(xintercept=0.75) +
  labs(x="BCAR4 Expression in BRCA",
       title="BCAR4 Grouping")

x <- b %>%
  group_by(ERBB2.amp, Group) %>%
  summarise(count=n())
x = x[c(1:6),]
x$n <- ifelse(x$ERBB2.amp==0, 8938, 371)
x$prop <- x$count / x$n

b$ERBB2.rank <- ifelse(b$ERBB2.exp>14, "High",
                       ifelse(b$ERBB2.exp>9, "Intermediate",
                              "Low"))
f <- survfit(Surv(OS.time / 365, OS) ~ ERBB2.rank, data = b)
ggsurvplot(f, data = b, pval=T,
           risk.table=T) +
  labs(title="PanCan")

ggplot() +
  geom_col(data=x, aes(x=factor(ERBB2.amp), y=prop, fill=Group))





f <- survfit(Surv(OS.time / 365, OS) ~ Group, data=b)
ggsurvplot(f, data=b, pval=T) +
  labs(title="BCAR4 expression tiers across 33 TCGA cohorts",
       x="Years")

mm$Group <- ifelse(mm$BCAR4.exp==0, "None", 
                   ifelse(mm$BCAR4.exp>2, "High", "Low"))
mm <- subset(mm, !is.na(mm$Group))
g1 <- mm %>%
  group_by(PAM50.cbioportal) %>%
  summarise(PAMCount=n())
g <- mm %>%
  group_by(PAM50.cbioportal, Group) %>%
  summarise(Count=n())
g3 <- merge(g1, g)
g3$Prop <- (g3$Count / g3$PAMCount) * 100
g3 <- subset(g3, !is.na(g3$PAM50.cbioportal))

g3$PAM50.cbioportal <- factor(g3$PAM50.cbioportal, levels=c("LumB", "Basal", "Her2", "LumA", "Normal"))
ggplot(data=g3, aes(x=PAM50.cbioportal, y=Prop, group=Group, fill=Group)) +
  geom_col() +
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(x="PAM50 Subtype",
       y="Proportion of Samples",
       title="BCAR4+ distribution among BRCA subtypes",
       fill="BCAR4 Expression")

xxx <- mm[,c("Patient","OS.time", "OS", "Overall.Survival.Months", "Overall.Survival.Status")]
ggplot() +
  geom_jitter(data=mm, aes(x=BRCA_Subtype_PAM50, y=PAM50.cbioportal)) +
  labs(x="TCGAbiolinks PAM50",
       y="cbioportal PAM50")
ggplot() +
  geom_point(data=mm, aes(x=OS.time, y=as.numeric(days_to_death))) +
  labs(x="TCGAbiolinks survival days",
       y="cbioportal survival days")
ggplot() +
  geom_point(data=mm, aes(x=OS.time, y=Overall.Survival.Months)) +
  labs(x="TCGAbiolinks survival days",
       y="Xena survival months")
ggplot() +
  geom_point(data=mm, aes(x=BCAR4.exp, y=BCAR4.fpkm)) +
  labs(x="Xena IlluminaHiSeq log2(RSME+1)",
       y="Xena HTSeq FPKM",
       title="BCAR4 Normalized Expression")
mn <- mm[,c("BCAR4.exp", "BCAR4.fpkm")]
mn$BCAR4.exp.rank <- rank(-mn$BCAR4.exp, ties.method="average")
mn$BCAR4.fpkm.rank <- rank(-mn$BCAR4.fpkm, ties.method="average")
ggplot(data=na.omit(mn), aes(x=BCAR4.exp.rank, BCAR4.fpkm.rank)) +
  geom_point(alpha=0.2) +
  geom_smooth() +
  geom_abline(slope=1, intercept=0, color="red") +
  scale_x_continuous(breaks=c(0,200,400,600,800,1000), limits=c(0,1000)) +
  scale_y_continuous(breaks=c(0,200,400,600,800,1000), limits=c(0,1000)) +
  labs(x="Xena IlluminaHiSeq Rank",
       y="Xena HTSeq FPKM Rank",
       title="BCAR4 Expression Rank Ordered",
       caption="Ties are given the same rank. Rank 1 is highest expression")
