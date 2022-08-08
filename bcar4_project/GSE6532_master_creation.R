

library(genefu)



setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")





load("data/GSE6532/LUMINAL.RData")


# 230854_at = BCAR4
# 3 that might corresponds to HER2
# 216836_s_at = Long ERBB-2
# 210930_s_at = Short HER-2
# 234354_x_at = short thing..?
# http://www.genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A35109780%2D35138441&hgsid=1188966687_dQlOfSdCVqv5ZauSXSaNLT2KOloj


bcar_id = "230854_at"
her2_long_id = "216836_s_at"
her2_short_id = "210930_s_at"

bcar.tam <- data.frame(data.tam[,c(bcar_id)])
bcar.tam$Sample <- rownames(bcar.tam)
bcar.tam$Treatment <- "Treated"
colnames(bcar.tam) <- c("BCAR4_exp", "Sample", "Treatment")

bcar.unt <- data.frame(data.untreated[,c(bcar_id)])
bcar.unt$Sample <- rownames(bcar.unt)
bcar.unt$Treatment <- "Untreated"
colnames(bcar.unt) <- c("BCAR4_exp", "Sample", "Treatment")

bcar <- rbind(bcar.tam, bcar.unt)

her2.tam <- data.frame(data.tam[,c(her2_long_id, her2_short_id)])
her2.tam$Sample <- rownames(her2.tam)
her2.tam$Treatment <- "Treated"
colnames(her2.tam) <- c("HER2_long_exp", "HER2_short_exp", "Sample", "Treatment")

her2.unt <- data.frame(data.untreated[,c(her2_long_id, her2_short_id)])
her2.unt$Sample <- rownames(her2.unt)
her2.unt$Treatment <- "Untreated"
colnames(her2.unt) <- c("HER2_long_exp", "HER2_short_exp", "Sample", "Treatment")

her2 <- rbind(her2.tam, her2.unt)

expression <- merge(bcar, her2)

tam.demo <- demo.tam
tam.demo$Treatment <- "Treated"
unt.demo <- demo.untreated
unt.demo$Treatment <- "Untreated"
demo <- rbind(tam.demo, unt.demo)
colnames(demo) <- c("Sample", "id", "series", "age", "grade", "size", "er",
                   "pgr", "node", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs", "Treatment")

master <- merge(demo, expression)




# PAM50
data(pam50)
data(pam50.robust)
data(pam50.scale)
pam50.genes <- pam50$centroids.map$EntrezGene.ID
annot.test <- as.data.frame(annot.tam)
colnames(annot.test) <- c("id", "EntrezGene.ID", "probe", "Alignment.score",
                          "Length.of.probe", "NCBI.gene.symbol", "HUGO.gene.symbol",
                          "Cytoband", "Alternative.symbols", "Description")
annot.test$EntrezGene.ID <- gsub(" ", "", annot.test$EntrezGene.ID)
annot.test <- subset(annot.test, annot.test$EntrezGene.ID %in% pam50.genes)
annot.test <- as.matrix(annot.test)
tam.pam50 <- molecular.subtyping(sbt.model="pam50", data=data.tam, annot=annot.test,
                          do.mapping=T, verbose=T)
subtypes <- as.data.frame(tam.pam50)
subtypes$Sample <- rownames(subtypes)
subtypes <- subtypes[,c("subtype", "Sample")]
colnames(subtypes) <- c("PAM50.local.genefu", "Sample")

annot.test <- as.data.frame(annot.untreated)
colnames(annot.test) <- c("id", "EntrezGene.ID", "probe", "Alignment.score",
                          "Length.of.probe", "NCBI.gene.symbol", "HUGO.gene.symbol",
                          "Cytoband", "Alternative.symbols", "Description")
annot.test$EntrezGene.ID <- gsub(" ", "", annot.test$EntrezGene.ID)
annot.test <- subset(annot.test, annot.test$EntrezGene.ID %in% pam50.genes)
#annot.test <- as.matrix(annot.test)

unt.pam50 <- molecular.subtyping(sbt.model="pam50", data=data.untreated, annot=annot.test,
                                  do.mapping=T, verbose=T)
unt.subtypes <- as.data.frame(unt.pam50)
unt.subtypes$Sample <- rownames(unt.subtypes)
unt.subtypes <- unt.subtypes[,c("subtype", "Sample")]
colnames(unt.subtypes) <- c("PAM50.local.genefu", "Sample")

all.subtypes <- rbind(subtypes, unt.subtypes)

master <- merge(master, all.subtypes)
saveRDS(master, file="data/GSE6532/master.rds")
