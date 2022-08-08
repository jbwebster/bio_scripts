


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


#crc <- readRDS("data/5hmC/ccle-expression/CRC.fpkm.rds")
#stad <- readRDS("data/5hmC/ccle-expression/STAD.fpkm.rds")
#brca <- readRDS("data/5hmC/ccle-expression/BRCA.fpkm.rds")
#lihc <- readRDS("data/5hmC/ccle-expression/LIHC.fpkm.rds")
#kirc <- readRDS("data/5hmC/ccle-expression/KIRC.fpkm.rds")

# LoVo.2 in crc - batch1
# SNU-308 in lihc - batch2
#MKN7 in stad (5hmC data is labeled as MNK7, rather than MKN7) - batch1
#T-47D in brca - batch1
#TUHR14TKB in kirc - batch2

#grabSubset <- function(ccle, cell_line) {
#  keep <- c("gene", "gene.name",  cell_line)
#  df <- ccle[,keep]
#}
#lovo <- grabSubset(crc, "G26204.LoVo.2")
#snu308 <- grabSubset(lihc, "G27481.SNU-308.2")
#mkn7 <- grabSubset(stad, "G25214.MKN7.1")
#t47d <- grabSubset(brca, "G41738.T-47D.5")
#tuhr.1 <- grabSubset(kirc, "G30558.TUHR14TKB.1")
#tuhr.2 <- grabSubset(kirc, "C836.TUHR14TKB.2")

#ccle_master <- merge(lovo, snu308)
#ccle_master <- merge(ccle_master, mkn7)
#ccle_master <- merge(ccle_master, tuhr.1)
#ccle_master <- merge(ccle_master, tuhr.2)
#ccle_master <- merge(ccle_master, t47d)

saveRDS(ccle_master, file="data/5hmC/master_files/master.fpkm.rds")



