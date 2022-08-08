

library(ggplot2)
library(ggpmisc)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline/Data/ProstateTMB")

tmb <- read.csv("TCGA-TMB.tsv", sep="\t", comment.char="#")


median(tmb$TMB) 
#=0.63, consistent with published analyses
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6882947/
my.formula = y ~ x
tmb2 = subset(tmb, tmb$Sample!="TCGA-XK-AAIW-01")#Remove outlier
tmb3 = subset(tmb, tmb$TMB<5)
ggplot(tmb, aes(x=TotalMutations, y=MutationsInPanel)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  #annotate("text", x=4000, y=10, label="R^2 = 0.84") +
  #annotate("text", x=4000, y=5, label="Exome~=117xPanel") +
  labs(x="TCGA-PRAD Whole Exome\n(#mutations)",
       y="TCGA-PRAD Targeted Panel\n(#mutations)",
       title="Mutation load comparison",
       caption="Excluding 1 outlier (#mutations=5923)")
ggplot(tmb2, aes(x=TMB, y=PanelTMB)) +
  geom_point() +
  geom_smooth(method="lm")


su2c_tmb <- read.csv("su2c_TMB.tsv", sep="\t", header=T)
su2c_tmb$Ratio = su2c_tmb$TotalInPanel / su2c_tmb$TotalAllRegions
dim(subset(su2c_tmb, Ratio==0))

my.formula = y ~ x
ggplot(su2c_tmb, aes(x=TotalAllRegions, y=TotalInPanel)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(x="SU2C Whole Exome\n(# mutations)",
       y="SU2C Mutations in Panel Regions\n(# mutations)",
       title="Missense/Nonsense Mutations in SU2C Prostate Cohort")

ggplot(su2c_tmb) +
  geom_point(aes(x=TotalAllRegions, y=TotalInPanel))
