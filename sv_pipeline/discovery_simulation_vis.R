

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(ggsci)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline")


## T2-ERG Calls
data <- read.csv("data/simulation/T2-ERG.100Iterations.low-dilutions.tsv", 
                   sep="\t", header=T)
wide <- read.csv("data/simulation/T2-ERG.100Iterations.low-dilutions.wide.tsv",
                 sep="\t", header=T)

df <- data %>%
  group_by(Dilution, Tool) %>%
  summarise(Calls=sum(CallMade))
iterations <- 100 # CHANGE AS NEEDED
#PB078 goes up to 7.5% 
#PB079 goes up to 12.5%
#PB202 goes up to 30%
#PB242 goes up to 25%
df$Total <- ifelse(df$Dilution <= 0.075, iterations * 4,
                    ifelse(df$Dilution > 0.075 & df$Dilution<=0.125, iterations * 3,
                    ifelse(df$Dilution>0.125 & df$Dilution<=0.25, iterations * 2,
                    ifelse(df$Dilution>0.25, iterations * 1, 0))))
# Aggregated calls
df$Percent <- df$Calls / df$Total
df$TumorContent <- df$Dilution * 100
df$Dilution <- as.factor(df$Dilution)
df$TumorContent <- as.factor(df$TumorContent)
df$Tool <- factor(df$Tool, levels=c("LiquidSCAN", "Aperture", "Factera", "SViCT"))
ggplot(subset(df, as.numeric(as.character(Dilution))>=0.01),
       aes(x=TumorContent, y=Percent*100, group=Tool, color=Tool)) +
  geom_line() +
  labs(x="Approx. Tumor Content (%)",
       y="Percent of Simulated Samples (%)",
       title="Detection of T2-ERG in Simulated Samples") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.5)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()#,
        #plot.title = element_text(size=22)
  )

individualPlot <- function(d, patient) {
  x <-subset(data, data$Patient == patient) %>%
    group_by(Dilution, Tool) %>%
    summarise(Calls=sum(CallMade))
  x$Percent = x$Calls / 100
  x$TumorContent <- x$Dilution * 100
  #x$TumorContent <- as.factor(x$TumorContent)
  x$Tool <- factor(x$Tool, levels=c("LiquidSCAN", "Aperture", "SViCT", "Factera"))
  p <- ggplot(x, aes(x=TumorContent, y=Calls, group=Tool, color=Tool)) +
    geom_line() +
    labs(x="Approx. Tumor Content (%)",
         y="Sensitivity (%)",
         title=paste0(patient, ": T2-ERG Detection")) +
    scale_x_continuous(expand = c(0, 0), breaks=c(5,10,15,20,25,30), limits=c(0,30)) +
    scale_y_continuous(limits = c(0,100), expand = c(0, 0.5)) +
    scale_color_jco() +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )
}
patient <- "PB078"
p1 <- individualPlot(data, "PB078")
p1
p2 <- individualPlot(data, "PB079")
p3 <- individualPlot(data, "PB202")
p4 <- individualPlot(data, "PB242")
grid.arrange(p1, p2, p3, p4)


# Read Support
df <- data
df$TumorContent <- df$Dilution * 100
df$Dilution <- as.factor(df$Dilution)
df$TumorContent <- as.factor(df$TumorContent)
df$Tool <- factor(df$Tool, levels=c("LiquidSCAN", "Aperture", "Factera", "SViCT"))
ggplot(df, aes(x=TumorContent, y=SR, fill=Tool)) +
  geom_boxplot() +
  ylim(0, 100) +
  labs(x="Approx. Tumor Content (%)",
       y="Split Read Support Reported",
       caption="Read support is reported as 0 if no T2-ERG call was made",
       title="Split Read Support for T2-ERG")
ggplot(df, aes(x=TumorContent, y=PE, fill=Tool)) +
  geom_boxplot() +
  ylim(0, 100) +
  labs(x="Approx. Tumor Content (%)",
       y="Paired-End Read Support Reported",
       caption="Read support is reported as 0 if no T2-ERG call was made",
       title="Paired-End Read Support for T2-ERG")
ggplot(subset(df, df$Tool %in% c("LiquidSCAN", "Aperture")), aes(x=PE, y=SR)) +
  geom_point() +
  facet_wrap(~Tool) +
  labs(x="Paired-End Support of T2-ERG",
       y="Split Read Support of T2-ERG",
       title="SR vs PE Support")
ggplot(df, aes(x=TumorContent, y=SR+PE, fill=Tool)) +
  geom_boxplot() +
  labs(title="Total Read Support")

# Per Sample Read Support 
wide <- subset(wide, wide$Dilution <= 0.3)
wide$L_Total = wide$L_PE + wide$L_SR
wide$A_Total = wide$A_PE + wide$A_SR
wide$TumorContent = wide$Dilution * 100
wide$TumorContent = as.factor(wide$TumorContent)
ggplot(wide, aes(x=L_Total, y=A_Total)) +
  geom_point(alpha=0.2, aes(color=TumorContent)) +
  geom_abline(intercept=0, slope=1, color="red") +
  #geom_smooth(method="lm") +
  labs(x="LiquidSCAN Total Read Support",
       y="Aperture Total Read Support",
       title="LiquidSCAN vs Aperture Read Support for T2-ERG",
       caption="Red line: y=x")


ggplot(wide, aes(x=L_SR, y=L_PE)) +
  geom_point()
ggplot(wide, aes(x=A_SR, y=A_PE)) +
  geom_point()

y <- subset(wide, wide$L_Call==0 & wide$A_Call==1)
yy <- subset(wide, wide$L_Call==1 & wide$A_Call==0)
x <- subset(wide, wide$L_Call == 1 & wide$A_Call==1)

wide$Difference <- wide$L_Total - wide$A_Total
ggplot(x, aes(x=Difference, fill=Patient)) + 
  geom_density(alpha=0.25) +
  facet_wrap(~TumorContent) +
  scale_x_continuous(breaks=c(-10,-5,0,5,10)) +
  geom_vline(xintercept=0) +
  #annotate("text", x = 5, y = 0.135, label = "More reads from LiquidSCAN", size=3) +
  #geom_segment(aes(x = 1, y = 0.13, xend = 10, yend = 0.13),
  #             arrow = arrow(length = unit(0.5, "cm"))) +
  labs(x="Difference in Total Read Support",
       title="Difference in T2-ERG Read Support",
       caption="Score = LiquidSCAN Reads - Aperture Reads\n
       Only plotting samples with a call made by both samples")



ggplot() +
  geom_jitter(data=subset(data, CallMade== 1 & (Tool=="Aperture")), #| Tool=="LiquidSCAN")), 
              alpha=0.05,
              aes(x=Dilution * 100, y=SR+PE, color=Patient)) +
  geom_smooth(data=subset(data, CallMade==1 & (Tool=="Aperture")), #| Tool=="LiquidSCAN")),
              alpha=-.25,
              method="gam", formula = y ~ s(x, bs = "cs", k=5),
              aes(x=Dilution * 100, y=SR+PE, color=Patient)) +
  facet_wrap(~Tool) +
  labs(x="Approx Tumor Content (%)",
       y="Total Read Support",
       title="PCa Simulation",
       caption="Showing only calls made")



