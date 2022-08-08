

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggvenn)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline")


df <- read.csv("data/crc/simulation_results/100iterations.tsv", sep="\t", header=T)

grouped <- df %>%
  group_by(Sample, Event, Dilution, Tool) %>%
  summarize(CallsMade = sum(CallMade))
grouped$Dilution <- grouped$Dilution * 100
#grouped$Dilution <- as.factor(grouped$Dilution * 100)
x <- data.frame(c("SRR11100269", "SRR11100320", "SRR11100350", "SRR11100358", "SRR11100358", "SRR11100416", "SRR11100416"),
                c(7, 3, 6, 5, 4, 1, 2),
                c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
                c("SViCT", "SViCT", "SViCT", "SViCT", "SViCT", "SViCT", "SViCT"),
                c(0, 0, 0, 0, 0, 0, 0),
                c(0, 0, 0, 0, 0, 0, 0))
colnames(x) <- colnames(grouped)
grouped <- rbind(grouped, x)
grouped$Tool <- factor(grouped$Tool, levels=c("LiquidSCAN", "Aperture", "SViCT", "Factera"))
grouped$Percent <- (grouped$CallsMade / 100 ) * 100

agg.grouped <- df %>%
  group_by(Dilution, Tool) %>%
  summarize(CallsMade = sum(CallMade))
agg.grouped$Dilution <- as.factor(agg.grouped$Dilution * 100)
agg.grouped$Tool <- factor(agg.grouped$Tool, levels=c("LiquidSCAN", "Aperture", "Factera"))
agg.grouped$Percent <- ifelse(agg.grouped$Dilution==30, (agg.grouped$CallsMade / 400) * 100,
                              (agg.grouped$CallsMade / 700) * 100)

ggplot() +
  geom_line(data=agg.grouped, aes(x=Dilution, y=Percent, , group=Tool, color=Tool)) +
  labs(x="Approx Tumor Content (%)",
       y="Percent of Simulated SV's Detected",
       caption="100 rounds of simulation",
       title="Aggregated") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot() +
  geom_col(data=agg.grouped, aes(x=Dilution, y=Percent, fill=Tool), 
           position="dodge") +
  labs(x="Approx Tumor Content (%)",
       y="Percent of Simulated Samples",
       caption="100 rounds of simulation",
       title="Aggregated") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


eventPlot <- function(grouped, event) {
  patient <- as.character(unique(subset(grouped, grouped$Event==event)$Sample))
  p <- ggplot() +
    #geom_col(data=subset(grouped, Event==event), aes(x=Dilution, y=Percent, fill=Tool), 
    #         position="dodge") +
    geom_line(data=subset(grouped, Event==event), aes(x=Dilution, y=Percent, 
                                                      group=Tool, color=Tool)) +
    labs(x="Approx Tumor Content (%)",
         y="Sensitivity (%)",
         title=paste0("Sample:", patient, " SV: ", event)) +
    scale_x_continuous(expand = c(0, 0), breaks=c(5,10,15,20,25,30), limits=c(0,30)) +
    scale_y_continuous(limits = c(0,100), expand = c(0, 0.5)) +
    scale_color_jco() +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )
  return(p)
}
p1 <- eventPlot(grouped, 1)
p2 <- eventPlot(grouped, 2)
p3 <- eventPlot(grouped, 3)
p4 <- eventPlot(grouped, 4)
p5 <- eventPlot(grouped, 5)
p6 <- eventPlot(grouped, 6)
p7 <- eventPlot(grouped, 7)
grid.arrange(p1,p2,p3,p4,p5,p6,p7)


# READ SUPPORT
l <- subset(df, Tool == "LiquidSCAN")
a <- subset(df, Tool == "Aperture")
f <- subset(df, Tool == "Factera")
colnames(l) <- c("Iteration", "Sample", "Event", "Dilution", "LCall", "LTool", "L_SR", "L_PE")
colnames(a) <- c("Iteration", "Sample", "Event", "Dilution", "ACall", "ATool", "A_SR", "A_PE")
colnames(f) <- c("Iteration", "Sample", "Event", "Dilution", "FCall", "FTool", "F_SR", "F_PE")
m <- merge(l, a)
m$Event <- as.factor(m$Event)
ggplot(m) +
  geom_point(aes(x=L_SR+L_PE, y=A_SR+A_PE, color=Event)) +
  geom_abline(slope=1, intercept=0) +
  labs(x="Total LiquidSCAN Read Support",
       y="Total Aperture Read Support",
       title="Read support for selected SVs")

ggplot(subset(df, df$Tool %in% c("LiquidSCAN", "Aperture")), aes(x=PE, y=SR,color=factor(Event))) +
  geom_point() +
  facet_wrap(~Tool) +
  labs(x="Paired-End Support of T2-ERG",
       y="Split Read Support of T2-ERG",
       title="SR vs PE Support")
m <- merge(l, a)
x <- list(LiquidSCAN = which(m$LCall==1), Aperture = which(m$ACall==1))
ggvenn(x)
m2 <- merge(l, f)
x <- list(LiquidSCAN = which(m2$LCall==1), Factera = which(m2$FCall==1))
ggvenn(x)
m3 <- merge(m, f)
x <- list(LiquidSCAN = which(m3$LCall==1), Aperture=which(m3$ACall==1), Factera=which(m3$FCall==1))
ggvenn(x) +
  labs(title="Overlap of Called Events",
       caption="100 iterations")

ggplot() +
  geom_point(data=subset(df, df$Event==2), aes(x=SR, y=PE, color=Tool)) +
  geom_abline(slope=1, intercept=0)


ggplot() +
  geom_point(data=subset(df, df$Tool=="Aperture"), aes(x=SR,y=PE), alpha=0.2) +
  facet_wrap(~Event) +
  geom_abline(slope=1, intercept=0) +
  labs(title="Aperture Read Support")

df$Event <- as.factor(df$Event)
ggplot() +
  geom_jitter(data=subset(df, CallMade== 1 & (Tool=="Aperture" | Tool=="LiquidSCAN")), 
              alpha=0.1,
              aes(x=Dilution * 100, y=SR+PE, color=Event)) +
  geom_smooth(data=subset(df, CallMade==1 & (Tool=="Aperture" | Tool=="LiquidSCAN")),
              alpha=-.25,
              aes(x=Dilution * 100, y=SR+PE, color=Event)) +
  facet_wrap(~Tool) +
  labs(x="Approx Tumor Content (%)",
       y="Total Read Support",
       title="CRC Simulation",
       caption="Showing only calls made")


