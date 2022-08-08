
library(data.table)


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


series96.exp <- read.csv("data/GSE4922/GSE4922-GPL96_series_matrix.droppedMetadata.txt",
                     header=T, sep="\t")
series97.exp <- read.csv("data/GSE4922/GSE4922-GPL97_series_matrix.droppedMetadata.txt",
                         header=T, sep="\t")

## PROBES OF INTEREST
# GPL97 array: 230854_at == BCAR4
# GPL97 array: 234354_x_at == HER2 3 region alternatively spliced, 175nt partial
# GPL96 array: 210930_s_at == HER2 alternatively spliced, complete cds

# Extract probe data
bcar_97 = "230854_at"
her2_short_97 = "234354_x_at"
her2_long_96 = "210930_s_at"
probes_97 <- c(bcar_97, her2_short_97)
series97_poi <- subset(series97.exp, series97.exp$ID_REF %in% probes_97)
series96_poi <- subset(series96.exp, series96.exp$ID_REF == her2_long_96)
# Relabel
series97 <- data.frame(t(series97_poi))
colnames(series97) <- c("BCAR4_97", "HER2_short_97")
series97$Sample <- rownames(series97)
series97 <- series97[-1,]
series96 <- data.frame(t(series96_poi))
colnames(series96) <- c("HER2_long_96")
series96$Sample <- rownames(series96)
series96 <- series96[-1,]

clinical <- read.csv("data/GSE4922/GSE4922_Clinical_file_for_both_Uppsala_Singapore_Samples.txt",
                     header=T, sep="\t")
clinical <- separate(data=clinical,
                          col=GSM.ID..A.B.chip.,
                          into=c("ChipA", "ChipB")) 

m <- merge(series96, clinical, by.x="Sample", by.y="ChipA")
m <- merge(series97, m, by.x="Sample", by.y="ChipB")
cnames <- colnames(m)
cnames[1] <- "SampleName_97"  
cnames[4] <- "SampleName_96"
cnames[2] <- "BCAR4_97.exp"
cnames[3] <- "HER2_short_97.exp"
cnames[5] <- "HER2_long_96.exp"
colnames(m) <- cnames
newnames <- c(cnames[1], cnames[4], cnames[6], cnames[8:21], cnames[2], cnames[3], cnames[5])
m <- m[,newnames]
m$BCAR4_97.exp <- as.numeric(m$BCAR4_97.exp)
m$HER2_long_96.exp <- as.numeric(m$HER2_long_96.exp)
m$HER2_short_97.exp <- as.numeric(m$HER2_short_97.exp)

saveRDS(m, file="data/GSE4922/master.rds")

# Correlation between the two HER2 probes isn't super great
#ggplot(m, aes(x=rank(HER2_long_96.exp), y=rank(HER2_short_97.exp))) +
#  geom_point() +
#  geom_smooth()



