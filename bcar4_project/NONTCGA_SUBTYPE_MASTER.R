
library(ggplot2)
library(survminer)
library(survival)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


########################################################
###### GSE6532

#Made with GSE6532_master_creation.R
master <- readRDS(file="data/GSE6532/master.rds")


###########################
#Basic exploration
m <- master
ggplot() +
  geom_density(data=m, aes(x=BCAR4_exp, fill=Treatment), alpha=0.5) +
  labs(x="BCAR4 Expression",
       title="GSE6532 B4 Expression Grouped by Treatment")
ggplot() +
  geom_density(data=m, aes(x=HER2_long_exp, fill=Treatment), alpha=0.5) +
  labs(x="HER2 Expression",
       title="GSE6532 HER2 Expression by Treatment")
#ggplot() +
#  geom_density(data=m, aes(x=HER2_short_exp, fill=Treatment), alpha=0.5) +
#  labs(x="HER2 Expression",
#       title="GSE6532 HER2 Expression by Treatment")

ggplot() +
  geom_point(data=m, aes(x=t.rfs, y=t.dmfs)) +
  labs(x="Relapse Free Survival (days)",
       y="Distant Metastasis Free Survival (days)",
       title="GSE6532 Comparison of Outcome Measurements")

ggplot() +
  geom_density(data=m, aes(x=HER2_long_exp, fill=PAM50.local.genefu), alpha=0.2) +
  labs(x="HER2 Expression",
       fill="PAM50 (In-House)",
       title="HER2 Expression By Subtype")

ggplot() +
  geom_density(data=m, aes(x=BCAR4_exp, fill=PAM50.local.genefu), alpha=0.2) +
  labs(x="BCAR4 Expression",
       fill="PAM50 (In-House)",
       title="BCAR4 Expression By Subtype")

m <- master
m$BCAR4Rank <- rank(m$BCAR4_exp)
q3 <- round((nrow(m) / 4 ) * 3, 0)
q3_threshold = subset(m,m$BCAR4Rank==q3)$BCAR4_exp
top10p <- round((nrow(m) / 10) * 9, 0)
top10p_threshold = subset(m, m$BCAR4Rank==top10p)$BCAR4_exp
ggplot() +
  geom_jitter(data=m, aes(x=PAM50.local.genefu, y=BCAR4_exp),
              height=0.05, width=0.05) +
  geom_hline(yintercept=q3_threshold) +
  geom_hline(yintercept=top10p_threshold) +
  labs(x="PAM50",
       y="BCAR4 Expression",
       title="BCAR4 Expression By Subtype",
       caption="Lines represent top 25% and top 10% thresholds")

ggplot() +
  geom_point(data=m, aes(x=HER2_long_exp, y=BCAR4_exp)) +
  geom_hline(yintercept=q3_threshold) +
  geom_hline(yintercept=top10p_threshold) +
  labs(x="HER2 Expression",
       y="BCAR4 Expression",
       title="HER2 vs BCAR4",
       caption="Lines represent top 25% and top 10% thresholds for BCAR4")

ggplot() +
  geom_point(data=master, aes(y=BCAR4_exp, x=t.rfs, color=factor(e.rfs))) +
  scale_color_manual(labels=c("LastUpdate", "Progressed"), values=c("blue", "red")) +
  labs(y="BCAR4 Expression",
       x="Relapse Free Survival",
       color="Event",
       title="Relapse Free Survival and BCAR4 Expression")

#########################
###Survival Curves

m <- master
q1 <- summary(m$BCAR4_exp)['1st Qu.']
q3 <- summary(m$BCAR4_exp)['3rd Qu.']
m$Group <- ifelse(m$BCAR4_exp>=q3, "TopQuartile",
                  ifelse(m$BCAR4_exp<=q1, "BottomQuartile", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 (Treated AND Untreated) BCAR4 Status")

m <- master
m$Rank <- rank(m$BCAR4_exp)
bottom_10p = round(nrow(m) / 10, 0)
top_10p = round((nrow(m)/10) * 9, 0)
m$Group <- ifelse(m$Rank <= bottom_10p, "Bottom10Percent",
                  ifelse(m$Rank >= top_10p, "Top10Percent", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 (Treated AND Untreated) BCAR4 Status")

m <- master
m$Rank <- rank(m$BCAR4_exp)
bottom_10p = round(nrow(m) / 10, 0)
top_10p = round((nrow(m)/10) * 9, 0)
m$Group <- ifelse(m$Rank <= bottom_10p, "Bottom10Percent",
                  ifelse(m$Rank >= top_10p, "Top10Percent", "Other"))
sm <- subset(m, m$Group!="Other" & m$Treatment=="Untreated")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 Untreated BCAR4 Status",
       caption="Percentiles calculated prior to splitting based on treatment")

m <- master
m$Rank <- rank(m$BCAR4_exp)
bottom_10p = round(nrow(m) / 10, 0)
top_10p = round((nrow(m)/10) * 9, 0)
m$Group <- ifelse(m$Rank <= bottom_10p, "Bottom10Percent",
                  ifelse(m$Rank >= top_10p, "Top10Percent", "Other"))
sm <- subset(m, m$Group!="Other" & m$Treatment=="Treated")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 Treated BCAR4 Status",
       caption="Percentiles calculated prior to splitting based on treatment")

m <- master
m$Rank <- rank(m$BCAR4_exp)
bottom_1q = round(nrow(m) / 4, 0)
top_3q = round((nrow(m)/4) * 3, 0)
m$Group <- ifelse(m$Rank <= bottom_1q, "BCAR4Bottom25%",
                  ifelse(m$Rank >= top_3q, "BCAR4Top25%", "Other"))
sm <- subset(m, m$Group!="Other" & m$Treatment=="Treated")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 Treated BCAR4 Status",
       caption="Percentiles calculated prior to splitting based on treatment")


m <- master
m$Rank <- rank(m$BCAR4_exp)
bottom_10p = round(nrow(m) / 10, 0)
top_10p = round((nrow(m)/10) * 9, 0)
m$BCAR4Group <- ifelse(m$Rank <= bottom_10p, "Bottom10Percent",
                  ifelse(m$Rank >= top_10p, "Top10Percent", "Other"))
m$Group <- ifelse(m$BCAR4Group=="Top10Percent" & m$PAM50.local.genefu!="Her2", "Top10PercentB4",
                  ifelse(m$BCAR4Group=="Bottom10Percent" & m$PAM50.local.genefu!="Her2", "Bottom10PercentB4", "Other"))
sm <- subset(m, m$Group!="Other" & m$Treatment=="Treated")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 Non-HER2 and Tamoxifen Treated")
# Splitting by 10th percentile is more statistically significant here (above)
# than splitting by quartiles

m <- master
m$BCAR4Rank <- rank(m$BCAR4_exp)
m$HER2Rank <- rank(m$HER2_long_exp)
bottom_1q <- round(nrow(m) / 4, 0)
top_3q <- round((nrow(m)/4)*3,0)
m$BCAR4Group <- ifelse(m$BCAR4Rank>=top_3q, "BCAR4TopQ",
                       ifelse(m$BCAR4Rank<=bottom_1q, "BCAR4BottomQ", "Other"))
m$HER2Group <- ifelse(m$HER2Rank>=top_3q, "HER2TopQ",
                     ifelse(m$HER2Rank<=bottom_1q, "HER2BottomQ", "Other"))
m$Group <- ifelse(m$BCAR4Group=="BCAR4TopQ" & m$HER2Group!="HER2TopQ", "BCAR4+hi/HER2-nothigh",
                  ifelse(m$BCAR4Group=="BCAR4BottomQ" & m$HER2Group=="HER2BottomQ", "BCAR4-low/HER2-low", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 BCAR4 and HER2",
       caption="BCAR4-hi=Top25%, HER2-nothigh=NotTop25%,\nBCAR4-low=Bottom25%, HER2-low=Bottom25%")

## This is how I did it when I originally looked at this data
m <- master
m$BCAR4Rank <- rank(m$BCAR4_exp)
top5p <- round((nrow(m)/20)*19, 0)
m$Group <- ifelse(m$BCAR4_exp>=1, "BCAR4>=1",
                  ifelse(m$BCAR4_exp<0 & m$HER2_long_exp<0, "BCAR4-low/HER2-low", "Other"))
sm <- subset(m, m$Group!="Other" & m$Treatment=="Treated")
f <- survfit(Surv(t.dmfs / 365, e.dmfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE6532 (Treated) BCAR4 and HER2",
       caption="Low = Lower than median. Distant Metastasis Free Survival")


m <- master
m$Rank <- rank(m$BCAR4_exp)
bottom_10p = round(nrow(m) / 10, 0)
top_10p = round((nrow(m)/10) * 9, 0)
m$BCAR4Group <- ifelse(m$Rank <= bottom_10p, "Bottom10Percent",
                       ifelse(m$Rank >= top_10p, "Top10Percent", "Other"))





########################################################
#### GSE4922
master <- readRDS("data/GSE4922/master.rds")




##########################
###### Expression Profile
m <- master

ggplot() +
  geom_density(data=m, aes(x=BCAR4_97.exp, fill=Cohort), alpha=0.5) +
  labs(x="BCAR4 Expression",
       title="BCAR4 Expression Both Cohorts")
ggplot() +
  geom_density(data=m, aes(x=HER2_long_96.exp, fill=Cohort), alpha=0.5) +
  labs(x="HER2 Expression",
       title="HER2 Expression Both Cohorts")

ggplot() +
  geom_point(data=m, aes(x=BCAR4_97.exp, y=HER2_long_96.exp, color=Cohort)) +
  geom_rug(data=m, aes(x=BCAR4_97.exp, y=HER2_long_96.exp),
           col=rgb(.5,0,0,alpha=.2)) +
  labs(x="BCAR4 Expression",
       y="HER2 Expression",
       title="BCAR4 vs HER2 by Cohort")


###########################
####### Survival Curves
colnames(master)[9] <- "DFS.Event"
m <- master

m$BCAR4Rank <- rank(m$BCAR4_97.exp)
top10p <- (nrow(m)/10)*9
bottom10p <- (nrow(m)/10)
m$Group <- ifelse(m$BCAR4Rank>=top10p, "BCAR4Top10%",
                  ifelse(m$BCAR4Rank<=bottom10p, "BCAR4Bottom10%", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(DFS.TIME..yrs. / 365, DFS.Event) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE4922 (All) BCAR4 Status")

m <- master
q1 <- summary(m$BCAR4_97.exp)['1st Qu.']
q3 <- summary(m$BCAR4_97.exp)['3rd Qu.']
m$Group <- ifelse(m$BCAR4_97.exp>=q3, "BCAR4Top25%",
                  ifelse(m$BCAR4_97.exp<=q3, "BCAR4Bottom25%", "Other"))
sm <- subset(m, m$Group!="Other" & m$ER...endocrine.therapy.only..1.included.in.survival.analysis.==1)
f <- survfit(Surv(DFS.TIME..yrs. / 365, DFS.Event) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE4922 (Endocrine Therapy Only) BCAR4 Status")

m <- master
m$BCAR4Rank <- rank(m$BCAR4_97.exp)
top10p <- (nrow(m)/10)*9
bottom10p <- (nrow(m)/10)
m$Group <- ifelse(m$BCAR4Rank>=top10p, "BCAR4Top10%",
                  ifelse(m$BCAR4Rank<=bottom10p, "BCAR4Bottom10%", "Other"))
sm <- subset(m, m$Group!="Other" & m$ER...endocrine.therapy.only..1.included.in.survival.analysis.==1)
f <- survfit(Surv(DFS.TIME..yrs. / 365, DFS.Event) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="GSE4922 (Endocrine Therapy Only) BCAR4 Status")





##################################################
### MERGE
gse6532 <- readRDS(file="data/GSE6532/master.rds")
gse4922 <- readRDS(file="data/GSE4922/master.rds")

# Get minimally comparable columns
a <- gse6532[,c("t.rfs", "e.rfs", "BCAR4_exp", "HER2_long_exp", "PAM50.local.genefu")]
a$Cohort <- "GSE6532"
colnames(a) <- c("t.rfs", "e.rfs", "BCAR4_exp", "HER2_exp", "PAM50", "Cohort")
gse4922$PAM50.local.genefu <- NA
b <- gse4922[,c("DFS.TIME..yrs.", "DFS.EVENT..0.censored..1.event.defined.as.any.type.of.recurrence..local..regional.or.distant..or.death.from.breast.cancer",
                "BCAR4_97.exp", "HER2_long_96.exp", "PAM50.local.genefu")]
b$Cohort <- "GSE4922"
colnames(b) <- c("t.rfs", "e.rfs", "BCAR4_exp", "HER2_exp", "PAM50", "Cohort")
bm <- median(b$BCAR4_exp)
b$BCAR4_exp <- b$BCAR4_exp - bm
hm <- median(b$HER2_exp)
b$HER2_exp <- b$HER2_exp - hm


merged <- rbind(a, b)

#############
## Check merge

ggplot() +
  geom_density(data=merged, aes(x=BCAR4_exp, fill=Cohort), alpha=0.5) +
  labs(x="BCAR4 Expression",
       title="BCAR4 Expression Post-Merge")
ggplot() +
  geom_density(data=merged, aes(x=HER2_exp, fill=Cohort), alpha=0.5) +
  labs(x="HER2 Expression",
       title="HER2 Expression Post-Merge")
## This plot might suggest that GSE4922 doesn't really have any PAM50 samples
## I already know GSE4922 is enriched for ER+ patients

m <- merged
m$Group <- ifelse(is.na(m$PAM50),  "GSE4922", as.character(m$PAM50))
m$Group <- factor(m$Group, levels=c("GSE4922", "Her2", "LumB", "LumA", "Basal", "Normal"))
ggplot() +
  geom_jitter(data=m, aes(x=Group, y=HER2_exp),
              height=0.05, width=0.05) +
  geom_violin(data=m, aes(x=Group, y=HER2_exp), 
              alpha=0.4) +
  labs(x="Group",
       y="HER2 Expression",
       caption="Showing cohort (GSE4922) or PAM50 (GSE6532)",
       title="HER2 Expression by Group")

############
## Merged survival curves

m <- merged

bottom10p <- nrow(m) / 10
top10p <- (nrow(m)/10)*9
m$Rank <- rank(m$BCAR4_exp)
m$Group <- ifelse(m$Rank>=top10p, "BCAR4Top10%", 
                  ifelse(m$Rank<=bottom10p, "BCAR4Bottom10%", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="Merged Cohorts BCAR4 Status")

m <- merged
q1 <- summary(m$BCAR4_exp)['1st Qu.']
q3 <- summary(m$BCAR4_exp)['3rd Qu.']
m$Group <- ifelse(m$BCAR4_exp>=q3, "TopQuartile",
                  ifelse(m$BCAR4_exp<=q1, "BottomQuartile", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="Merged Cohorts BCAR4 Status")


m <- merged
bottom10p <- nrow(m) / 20
top10p <- (nrow(m)/20)*19
m$Rank <- rank(m$BCAR4_exp)
m$Group <- ifelse(m$Rank>=top10p, "BCAR4Top10%", 
                  ifelse(m$Rank<=bottom10p, "BCAR4Bottom10%", "Other"))
sm <- subset(m, m$Group!="Other")
f <- survfit(Surv(t.rfs / 365, e.rfs) ~ Group, data = sm)
ggsurvplot(f, data = sm, pval=T,
           risk.table=T) +
  labs(title="Merged Cohorts BCAR4 Status")
