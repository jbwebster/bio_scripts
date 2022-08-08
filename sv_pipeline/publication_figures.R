
library(ggplot2)
library(dplyr)
library(gridExtra)
library(gt)
library(ggpubr)
library(ggsci)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline")

###
# May need to rename these depending on how the paper goes, but starting plan is:
# Fig1 will likely be a workflow schematic/diagram
# Fig2 will be PCa related results from JCOPO, with:
## a) recreation of AR events from cfDNA, b) simulation diagram and c) simulation results
# Fig3/Table1 results from Horizon Discovery
# Fig4 CRC simulation
# Fig5 HCC1395 results


################ Recreation of JCOPO Findings ################################

### Only LiquidSCAN and Aperture detect the AR/Enhancer Dups from the JCOPO cohort
# Output SVs were manually compared to breakpoints reported in original analysis
# Aperture reported many other SVs in the region that did not correspond to published SVs
res <- c(100, 100, 0, 0) # These numbers were collected manually from tool outputs
Tool <- c("Pipeline", "Aperture", "SViCT", "Factera")
Tool <- factor(Tool, levels=Tool)
df1 <- data.frame(res, Tool)
f2a <- 
  ggplot() +
  geom_col(data=df1, width=0.8,
           aes(x=Tool, y=res, fill=Tool)) +
  scale_fill_jco() +
  labs(x="Tool",
       y="Events Detected (%)",
       #tag="A",
       title="AR/Enhancer Duplications"
       ) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  theme_bw() +
  #theme(aspect.ratio = 1/1.5) +
  scale_x_discrete(expand = c(0.25,0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),  
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        plot.title = element_text(size=60),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=50),
        legend.text = element_text(size=40))
f2a

# Assuming 4 T2-ERG fusions
res <- c(100, 100, 75, 25) # These numbers were collected manually from tool outputs
Tool <- c("Pipeline", "Aperture", "SViCT", "Factera")
Tool <- factor(Tool, levels=Tool)
df2 <- data.frame(res, Tool)
f2b <- 
  ggplot() +
  geom_col(data=df2, width=0.8,
           aes(x=Tool, y=res, fill=Tool)) +
  scale_fill_jco() +
  labs(x="Tool",
       y="Events Detected (%)",
       #tag="B",
       title="TMPRSS2-ERG Fusions"
       ) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  theme_bw() +
  scale_x_discrete(expand = c(0.25,0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),  
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=30),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=25),
        legend.text = element_text(size=23))

f2a
f2b
grid.arrange(f2a, f2b, ncol=1)
#png with height=500 and width=670 used for thesis proposal
#ggarrange(f2a, f2b, ncol=2, labels=c("A", "B"))
################ JCOPO Prostate Simulation of T2-ERG #########################
data <- read.csv("data/simulation/T2-ERG.100Iterations.low-dilutions.tsv", 
                 sep="\t", header=T)
wide <- read.csv("data/simulation/T2-ERG.100Iterations.low-dilutions.wide.tsv",
                 sep="\t", header=T)

df2 <- data %>%
  group_by(Dilution, Tool) %>%
  summarise(Calls=sum(CallMade))
iterations <- 100 # CHANGE AS NEEDED
#PB078 goes up to 7.5% 
#PB079 goes up to 12.5%
#PB202 goes up to 30%
#PB242 goes up to 25%
df2$Total <- ifelse(df2$Dilution <= 0.075, iterations * 4,
                   ifelse(df2$Dilution > 0.075 & df2$Dilution<=0.125, iterations * 3,
                          ifelse(df2$Dilution>0.125 & df2$Dilution<=0.25, iterations * 2,
                                 ifelse(df2$Dilution>0.25, iterations * 1, 0))))
# Aggregated calls
df2$Percent <- df2$Calls / df2$Total
df2$TumorContent <- df2$Dilution * 100
#df2$Dilution <- as.factor(df2$Dilution)
#df2$TumorContent <- as.factor(df2$TumorContent)
x <- data.frame(c(0.3, 0.25), c("SViCT", "SViCT"), c(0, 0), c(100, 200), c(0, 0), c(30, 25))
names(x) <- c("Dilution", "Tool", "Calls", "Total", "Percent", "TumorContent")
#x$Dilution <- as.factor(x$Dilution)
#x$TumorContent <- as.factor(x$TumorContent)
df2 <- rbind(df2, x)
df2$Tool <- ifelse(df2$Tool == "LiquidSCAN", "Pipeline", df2$Tool)
df2 <- subset(df2, df2$Tool != "SViCT")
#df2$Tool <- factor(df2$Tool, levels=c("Pipeline", "Aperture", "SViCT", "Factera"))
df2$Tool <- factor(df2$Tool, levels=c("Pipeline", "Aperture", "Factera"))

jco_colors <- c("#0073C2FF", "#EFC000FF", "#CD534CFF")

f3b <- 
ggplot() +
  geom_line(data=subset(df2, as.numeric(as.character(Dilution))>=0.01),
            size=3,
            aes(x=TumorContent, y=Percent*100, 
                group=Tool, color=Tool)) +
  #scale_color_jco() +
  scale_color_manual(values=jco_colors) +
  labs(x="Simulated Tumor Content (%)",
       y="Sensitivity (%)",
       title="Prostate Cancer Cohort",
       tag="D") +
  scale_x_continuous(expand = c(0, 0), breaks = c(5,10,15,20,25,30)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.5)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()#,
        #plot.title = element_text(size=22)
        ) +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),  
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=30),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=25),
        legend.text = element_text(size=23))
  
f3b

################ Horizon dataset results #########################
Tool <- c("PACT", "Aperture", "Factera", "SViCT")
SLC34A2_ROS1 <- c("Yes", "Yes", "Yes", NA)
CCDC6_RET <- c("Yes", "Yes", "No", NA)
OtherCalls <- c(4, 1636, 13, NA)
horizon <- data.frame(Tool, SLC34A2_ROS1, CCDC6_RET, OtherCalls)
colnames(horizon) <- c("Tool", "SLC34A2-ROS1", "CCDC6-RET", "Other Calls")

Tool <- c("PACT", "Aperture", "Factera")
SLC34A2_ROS1 <- c("Yes", "Yes", "Yes")
CCDC6_RET <- c("Yes", "Yes", "No")
OtherCalls <- c(4, 1636, 13)
horizon <- data.frame(Tool, SLC34A2_ROS1, CCDC6_RET, OtherCalls)
colnames(horizon) <- c("Tool", "SLC34A2-ROS1", "CCDC6-RET", "Other Calls")

Tool <- c("PACT", "Aperture", "Factera")
SLC34A2_ROS1 <- c("Yes", "Yes", "Yes")
CCDC6_RET <- c("Yes", "Yes", "No")
horizon <- data.frame(Tool, SLC34A2_ROS1, CCDC6_RET)
colnames(horizon) <- c("Tool", "SLC34A2-ROS1", "CCDC6-RET")


table1 <- horizon %>%
  gt(rowname_col="Tool") %>% 
  tab_stubhead(label="Tool") %>%
  tab_header(
    title = "Horizon cfDNA Reference Data"
  ) %>%
  tab_spanner(label = "Validated Calls", columns = matches("SLC34A2-ROS1|CCDC6-RET")) %>%
  cols_align(
    align = "center",
    columns = c("SLC34A2-ROS1", "CCDC6-RET")
  ) 

table1

jco_colors <- c("#0073C2FF", "#EFC000FF", "#CD534CFF")
df = data.frame(Tool, OtherCalls)
df$Tool <- factor(df$Tool, levels=c("PACT", "Aperture", "Factera"))
fp = ggplot(df) +
  geom_col(aes(x=Tool,y=OtherCalls, fill=Tool)) +
  scale_fill_manual(values=jco_colors) +
  labs(y="# False Positives",
       title="False Positives in Reference Data",
       #tag="E"
       ) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()#,
        #plot.title = element_text(size=22)
  ) +
  theme(axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),  
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        plot.title = element_text(size=60),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=50),
        legend.text = element_text(size=40))
  

fp
# Accuracy = True positives / Total calls
# Would make other tools look terrible, but even ours would only have 33%, which sounds awful

################ CRC Simulation of 7 SVs #########################

df3 <- read.csv("data/crc/simulation_results/100iterations.tsv", sep="\t", header=T)

agg.grouped <- df3 %>%
  group_by(Dilution, Tool) %>%
  summarize(CallsMade = sum(CallMade))

agg.grouped$Percent <- ifelse(agg.grouped$Dilution==.3, (agg.grouped$CallsMade / 400) * 100,
                              (agg.grouped$CallsMade / 700) * 100)
agg.grouped$TumorContent <- agg.grouped$Dilution * 100
#agg.grouped$TumorContent <- as.factor(agg.grouped$TumorContent)
Dilution <- unique(agg.grouped$Dilution)
Tool <- rep("SViCT", 15)
CallsMade <- rep(0, 15)
Percent <- rep(0, 15)
tmp <- data.frame(Dilution, Tool, CallsMade, Percent)
tmp$TumorContent <- tmp$Dilution * 100
#tmp$TumorContent <- as.factor(tmp$Dilution * 100)
agg.grouped <- rbind(agg.grouped, tmp)
agg.grouped$Tool <- ifelse(agg.grouped$Tool == "LiquidSCAN", "Pipeline", agg.grouped$Tool)
agg.grouped <- subset(agg.grouped, agg.grouped$Tool != "SViCT")
#agg.grouped$Tool <- factor(agg.grouped$Tool, levels=c("Pipeline", "Aperture", "SViCT", "Factera"))
agg.grouped$Tool <- factor(agg.grouped$Tool, levels=c("Pipeline", "Aperture", "Factera"))

jco_colors <- c("#0073C2FF", "#EFC000FF", "#CD534CFF")

f3a <- 
ggplot() +
  geom_line(data=agg.grouped,
            size=3,
            aes(x=TumorContent, y=Percent, 
                group=Tool, color=Tool)) +
  #scale_color_jco() +
  scale_color_manual(values=jco_colors) +
  labs(x="Simulated Tumor Content (%)",
       y="Sensitivity (%)",
       title="Colorectal Cancer Cohort",
       tag="C") +
  scale_x_continuous(expand = c(0, 0), breaks=c(5,10,15,20,25,30)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.5)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()#,
        #plot.title = element_text(size=22
  ) +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),  
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=30),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=25),
        legend.text = element_text(size=23))
f3a



### Review figures
f2a
f2c
f2b
f3a
f3b
table1
grid.arrange(f3a,f3b, nrow=1)
grid.arrange(f2a,f2b,f3a,f3b,fp)
fp
# For thesis proposal
f3a <- f3a + labs(tag="C")
f3b <- f3b + labs(tag="D")
grid.arrange(f3a, f3b, ncol=1)
# Used height=450 width=450

grid.arrange(f2a, f2b, f3a, f3b, ncol=2)


#####Merging CRC and PCa simulations
# CRC
df3 <- read.csv("data/crc/simulation_results/100iterations.tsv", sep="\t", header=T)

agg.grouped <- df3 %>%
  group_by(Dilution, Tool) %>%
  summarize(CallsMade = sum(CallMade))
agg.grouped$Total <- ifelse(agg.grouped$Dilution == 0.3, 400, 700)
agg.grouped$Percent <- ifelse(agg.grouped$Dilution==.3, (agg.grouped$CallsMade / 400) * 100,
                              (agg.grouped$CallsMade / 700) * 100)
agg.grouped$TumorContent <- agg.grouped$Dilution * 100
#agg.grouped$TumorContent <- as.factor(agg.grouped$TumorContent)
Dilution <- unique(agg.grouped$Dilution)
Tool <- rep("SViCT", 15)
CallsMade <- rep(0, 15)
Percent <- rep(0, 15)
tmp <- data.frame(Dilution, Tool, CallsMade, Percent)
tmp$TumorContent <- tmp$Dilution * 100
#tmp$TumorContent <- as.factor(tmp$Dilution * 100)
agg.grouped <- rbind(agg.grouped, tmp)
agg.grouped$Tool <- ifelse(agg.grouped$Tool == "LiquidSCAN", "Pipeline", agg.grouped$Tool)
agg.grouped <- subset(agg.grouped, agg.grouped$Tool != "SViCT")
#agg.grouped$Tool <- factor(agg.grouped$Tool, levels=c("Pipeline", "Aperture", "SViCT", "Factera"))
agg.grouped$Tool <- factor(agg.grouped$Tool, levels=c("PACT", "Aperture", "Factera"))

# PCa
data <- read.csv("data/simulation/T2-ERG.100Iterations.low-dilutions.tsv", 
                 sep="\t", header=T)
wide <- read.csv("data/simulation/T2-ERG.100Iterations.low-dilutions.wide.tsv",
                 sep="\t", header=T)

df2 <- data %>%
  group_by(Dilution, Tool) %>%
  summarise(Calls=sum(CallMade))
iterations <- 100 # CHANGE AS NEEDED
#PB078 goes up to 7.5% 
#PB079 goes up to 12.5%
#PB202 goes up to 30%
#PB242 goes up to 25%
df2$Total <- ifelse(df2$Dilution <= 0.075, iterations * 4,
                    ifelse(df2$Dilution > 0.075 & df2$Dilution<=0.125, iterations * 3,
                           ifelse(df2$Dilution>0.125 & df2$Dilution<=0.25, iterations * 2,
                                  ifelse(df2$Dilution>0.25, iterations * 1, 0))))
# Aggregated calls
df2$Percent <- df2$Calls / df2$Total
df2$TumorContent <- df2$Dilution * 100
#df2$Dilution <- as.factor(df2$Dilution)
#df2$TumorContent <- as.factor(df2$TumorContent)
x <- data.frame(c(0.3, 0.25), c("SViCT", "SViCT"), c(0, 0), c(100, 200), c(0, 0), c(30, 25))
names(x) <- c("Dilution", "Tool", "Calls", "Total", "Percent", "TumorContent")
#x$Dilution <- as.factor(x$Dilution)
#x$TumorContent <- as.factor(x$TumorContent)
df2 <- rbind(df2, x)
df2$Tool <- ifelse(df2$Tool == "LiquidSCAN", "PACT", df2$Tool)
df2 <- subset(df2, df2$Tool != "SViCT")
#df2$Tool <- factor(df2$Tool, levels=c("Pipeline", "Aperture", "SViCT", "Factera"))
df2$Tool <- factor(df2$Tool, levels=c("PACT", "Aperture", "Factera"))

m <- cbind(df2, agg.grouped)
m$TotalPossible = m$Total...10 + m$Total...4
m$TotalCalls = m$Calls + m$CallsMade
m$TotalPercent = (m$TotalCalls / m$TotalPossible)

colnames(m) <- c("Dilution", "Tool", "Calls", "TotalA",
                 "PercentA", "TumorContentA", "DilutionB",
                 "ToolB", "CallsB", "PercentB", "TumorContentB",
                 "TotalB", "TotalCalls", "TotalPossible", "TotalPercent")

ggplot() +
  geom_line(data=m,
            size=3,
            aes(x=TumorContentA, y=TotalPercent*100, 
                group=Tool, color=Tool)) +
  #scale_color_jco() +
  scale_color_manual(values=jco_colors) +
  labs(x="Simulated Tumor Content (%)",
       y="Sensitivity (%)",
       title="In silico simulation",
       #tag="C"
       ) +
  scale_x_continuous(expand = c(0, 0), breaks=c(5,10,15,20,25,30)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.5)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()#,
        #plot.title = element_text(size=22
  ) +
  theme(axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),  
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        plot.title = element_text(size=60),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=50),
        legend.text = element_text(size=40))

###
# Merging JCOPO re-analysis findings
# AR (out of 5)
AR <- c(5, 5, 0, 0) # These numbers were collected manually from tool outputs
ARpossible <- c(5,5,5,5)
Tool <- c("PACT", "Aperture", "SViCT", "Factera")
Tool <- factor(Tool, levels=Tool)
df1 <- data.frame(AR, Tool, ARpossible)

# T2
TE <- c(4, 4, 2, 1) # These numbers were collected manually from tool outputs
TEpossible <- c(4,4,4,4)
Tool <- c("PACT", "Aperture", "SViCT", "Factera")
Tool <- factor(Tool, levels=Tool)
df2 <- data.frame(TE, Tool, TEpossible)

df <- merge(df1,df2)
df$TotalPossible <- df$ARpossible + df$TEpossible
df$Calls <- df$AR + df$TE 
df$Percent <- (df$Calls / df$TotalPossible) * 100

ggplot() +
  geom_col(data=df, width=0.8,
           aes(x=Tool, y=Percent, fill=Tool)) +
  scale_fill_jco() +
  labs(x="Tool",
       y="Events Detected (%)",
       #tag="A",
       title="Re-Analysis of Published cfDNA Results"
  ) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  theme_bw() +
  #theme(aspect.ratio = 1/1.5) +
  scale_x_discrete(expand = c(0.25,0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),  
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        plot.title = element_text(size=60),
        legend.key.height= unit(4, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.title = element_text(size=50),
        legend.text = element_text(size=40))
