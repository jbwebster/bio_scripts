
library(ggplot2)
library(gridExtra)
library(dplyr)

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")



# Load master files
fpkm <- readRDS(file="data/5hmC/master_files/master.fpkm.rds")
ht <- readRDS(file="data/5hmC/master_files/ht_annot.rds")
liver2 <- readRDS(file="data/5hmC/master_files/Liver2_annot.rds")
livermet2 <- readRDS(file="data/5hmC/master_files/LiverMet2_annot.rds")
livermet7 <- readRDS(file="data/5hmC/master_files/LiverMet7_annot.rds")
lovo <- readRDS(file="data/5hmC/master_files/lovo_annot.rds")
mkn7 <- readRDS(file="data/5hmC/master_files/mkn7_annot.rds")
snu <- readRDS(file="data/5hmC/master_files/snu_annot.rds")
t47d <- readRDS(file="data/5hmC/master_files/t47d_annot.rds")
tuhr <- readRDS(file="data/5hmC/master_files/tuhr_annot.rds")
uninv <- readRDS(file="data/5hmC/master_files/Uninv_annot.rds")


getPromoters <- function(annot, promoters, fpkm_sub, cell) {
  colnames(fpkm_sub) <- c("geneSymbol", cell)
  fpkm_bcar4 <- subset(fpkm_sub, fpkm_sub$geneSymbol=="BCAR4")
  p <- subset(annot, annot$geneSymbol=="BCAR4" & annot$annotation=="Promoter")
  if(nrow(p)>0) {
    for(i in 1:nrow(p)) {
      promoters[nrow(promoters) + 1,] = c(cell, p$V5, fpkm_bcar4[1,2])
    }
  }
  else {
    promoters[nrow(promoters) + 1,] = c(cell, 0, fpkm_bcar4[1,2])
  }
  return(promoters)
}
promoters <- data.frame(Cell=character(),
                        PromoterSignal=numeric(),
                        fpkm=numeric())
promoters <- getPromoters(lovo, promoters, fpkm[,c("gene.name", "G26204.LoVo.2")], "LoVo")
promoters <- getPromoters(snu, promoters, fpkm[,c("gene.name", "G27481.SNU-308.2")], "SNU")
promoters <- getPromoters(mkn7, promoters, fpkm[,c("gene.name", "G25214.MKN7.1")], "MKN7")
promoters <- getPromoters(tuhr, promoters, fpkm[,c("gene.name", "G30558.TUHR14TKB.1")], "TUHR14")
promoters <- getPromoters(t47d, promoters, fpkm[,c("gene.name", "G41738.T-47D.5")], "T47D")
promoters$PromoterSignal <- as.numeric(promoters$PromoterSignal)
promoters$fpkm <- as.numeric(promoters$fpkm)
ggplot(promoters, aes(x=PromoterSignal, y=fpkm)) +
  geom_point() +
  geom_text(aes(label=Cell),hjust=0, vjust=0) +
  labs(x="PromoterSignalScore",
       y="CCLE FPKM",
       title="5hmC Promoter Signal vs FPKM")

createMultiPromoterPlot <- function(cell, fpkm.subset, cell.name) {
  colnames(fpkm.subset) <- c("gene.name", "FPKM")
  promoters <- subset(cell, cell$annotation=="Promoter")
  promSubset <- promoters[,c("V5", "geneSymbol")]
  colnames(promSubset) <- c("PromoterSignalScore", "gene.name")  
  promSubset <- promSubset %>%
    group_by(gene.name) %>%
    summarise(MeanPromoterSignalScore = mean(PromoterSignalScore))
  promoter.fpkm <- merge(promSubset, fpkm.subset)
  promoter.fpkm <- promoter.fpkm %>%
    group_by(gene.name) %>%
    summarise(MeanPromoterSignalScore=mean(MeanPromoterSignalScore), MeanFPKM=mean(FPKM))
  p <- ggplot(promoter.fpkm, aes(x=MeanPromoterSignalScore, y=log(MeanFPKM+1))) +
    geom_point() +
    geom_smooth() +
    labs(title=cell.name)
}

p.Lovo <- createMultiPromoterPlot(lovo, fpkm[,c("gene.name", "G26204.LoVo.2")], "LoVo")
p.snu <- createMultiPromoterPlot(snu,  fpkm[,c("gene.name", "G27481.SNU-308.2")], "SNU")
p.mkn7 <- createMultiPromoterPlot(mkn7,  fpkm[,c("gene.name", "G25214.MKN7.1")], "MKN7")
p.tuhr <- createMultiPromoterPlot(tuhr, fpkm[,c("gene.name", "G30558.TUHR14TKB.1")], "TUHR14")
p.t47d <- createMultiPromoterPlot(t47d, fpkm[,c("gene.name", "G41738.T-47D.5")], "T47D")
grid.arrange(p.Lovo, p.snu, p.mkn7, p.tuhr, p.t47d, top="5hmC Promoter Signal vs CCLE FPKM")


createMultiGenePlot <- function(cell, fpkm.subset, cell.name) {
  colnames(fpkm.subset) <- c("gene.name", "FPKM")
  peaks <- subset(cell, cell$annotation=="Promoter" | cell$annotation=="Downstream (<=300bp)")
  promSubset <- promoters[,c("V5", "geneSymbol")]
  colnames(promSubset) <- c("SignalScore", "gene.name")  
  promSubset <- promSubset %>%
    group_by(gene.name) %>%
    summarise(MeanSignalScore = mean(SignalScore))
  promoter.fpkm <- merge(promSubset, fpkm.subset)
  promoter.fpkm <- promoter.fpkm %>%
    group_by(gene.name) %>%
    summarise(MeanSignalScore=mean(MeanSignalScore), MeanFPKM=mean(FPKM))
  p <- ggplot(promoter.fpkm, aes(x=MeanSignalScore, y=log(MeanFPKM+1))) +
    geom_point() +
    geom_smooth() +
    labs(title=cell.name)
}

p.Lovo <- createMultiGenePlot(lovo, fpkm[,c("gene.name", "G26204.LoVo.2")], "LoVo")
p.snu <- createMultiGenePlot(snu,  fpkm[,c("gene.name", "G27481.SNU-308.2")], "SNU")
p.mkn7 <- createMultiGenePlot(mkn7,  fpkm[,c("gene.name", "G25214.MKN7.1")], "MKN7")
p.tuhr <- createMultiGenePlot(tuhr, fpkm[,c("gene.name", "G30558.TUHR14TKB.1")], "TUHR14")
p.t47d <- createMultiGenePlot(t47d, fpkm[,c("gene.name", "G41738.T-47D.5")], "T47D")
grid.arrange(p.Lovo, p.snu, p.mkn7, p.tuhr, p.t47d, 
             top="5hmC TSS/TES Signal vs CCLE FPKM",
             bottom="Using peaks labeled as 'Promoter' or 'Downstream (<=300bp)")

ht$Sample <- "ht"
liver2$Sample <- "liver2"
livermet2$Sample <- "livermet2"
livermet7$Sample <- "livermet7"
lovo$Sample <- "lovo"
mkn7$Sample <- "mkn7"
snu$Sample <- "snu"
t47d$Sample <- "t47d"
uninv$Sample <- "uninv"
tuhr$Sample <- "tuhr"

qc <- rbind(ht, liver2)
qc <- rbind(qc, livermet2)
qc <- rbind(qc, livermet7)
qc <- rbind(qc, lovo)
qc <- rbind(qc, mkn7)
qc <- rbind(qc, snu)
qc <- rbind(qc, t47d)
qc <- rbind(qc, uninv)
qc <- rbind(qc, tuhr)


ggplot(qc, aes(x=Sample)) +
  geom_bar() +
  labs(x="Sample",
       y="# of 5hmC Peaks",
       title="Peaks per Sample")
