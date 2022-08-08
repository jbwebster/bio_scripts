
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install('clusterProfiler')
BiocManager::install('annotables')
BiocManager::install('org.Hs.eg.db')
install.packages('igraph')
install.packages('clusterProfiler')
install.packages('annotables')
install.packages("devtools")
devtools::install_github("stephenturner/annotables")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(annotables)
library(org.Hs.eg.db)
library(dplyr)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


samplefiles1 <- list.files("data/5hmC/batch1/macs2-output", pattern= "chrms.bed", full.names=T)
samplefiles2 <- list.files("data/5hmC/batch2/macs2-output", pattern="chrms.bed", full.names=T)
samplefiles <- c(samplefiles1, samplefiles2)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Lovo", "MKN7",
                        "T47D", "1787_2_Uninv",
                        "1787_2_Liver", "1782_2_Liver_Met",
                        "2042_7_Liver_Met", "HT342CI",
                        "SNU-308", "TUHR14TKB")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE) # MODIFY THIS VALUE

plotAnnoPie(peakAnnoList[["TUHR14TKB"]])
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList, 
              title="Distribution of transcription factor-binding loci \n relative to TSS")


lovo <- as.data.frame(peakAnnoList[["Lovo"]]@anno)
mkn7 <- as.data.frame(peakAnnoList[["MKN7"]]@anno)
t47d <- as.data.frame(peakAnnoList[["T47D"]]@anno)
Uninv <- as.data.frame(peakAnnoList[["1787_2_Uninv"]]@anno)
Liver2 <- as.data.frame(peakAnnoList[["1787_2_Liver"]]@anno)
LiverMet2 <- as.data.frame(peakAnnoList[["1782_2_Liver_Met"]]@anno)
LiverMet7 <- as.data.frame(peakAnnoList[["2042_7_Liver_Met"]]@anno)
ht <- as.data.frame(peakAnnoList[["HT342CI"]]@anno)
snu <- as.data.frame(peakAnnoList[["SNU-308"]]@anno)
tuhr <- as.data.frame(peakAnnoList[["TUHR14TKB"]]@anno)

getIDs <- function(df) {
  entrezids <- unique(df$geneId)
  entrez2gene <- grch37 %>%
    filter(entrez %in% entrezids) %>%
    dplyr::select(entrez, symbol)
  m <- match(df$geneId, entrez2gene$entrez)
  df_annot <- cbind(df[,1:14], geneSymbol=entrez2gene$symbol[m], df[,15:16])
}

lovo_annot <- getIDs(lovo)
mkn7_annot <- getIDs(mkn7)
t47d_annot <- getIDs(t47d)
Uninv_annot <- getIDs(Uninv)
Liver2_annot <- getIDs(Liver2)
LiverMet2_annot <- getIDs(LiverMet2)
LiverMet7_annot <- getIDs(LiverMet7)
ht_annot <- getIDs(ht)
snu_annot <- getIDs(snu)
tuhr_annot <- getIDs(tuhr)

saveRDS(lovo_annot, "data/5hmC/master_files/lovo_annot.rds")
saveRDS(mkn7_annot, "data/5hmC/master_files/mkn7_annot.rds")
saveRDS(t47d_annot, "data/5hmC/master_files/t47d_annot.rds")
saveRDS(Uninv_annot, "data/5hmC/master_files/Uninv_annot.rds")
saveRDS(Liver2_annot, "data/5hmC/master_files/Liver2_annot.rds")
saveRDS(LiverMet2_annot, "data/5hmC/master_files/LiverMet2_annot.rds")
saveRDS(LiverMet7_annot, "data/5hmC/master_files/LiverMet7_annot.rds")
saveRDS(ht_annot, "data/5hmC/master_files/ht_annot.rds")
saveRDS(snu_annot, "data/5hmC/master_files/snu_annot.rds")
saveRDS(tuhr_annot, "data/5hmC/master_files/tuhr_annot.rds")
