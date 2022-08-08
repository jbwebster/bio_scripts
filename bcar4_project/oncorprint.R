

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")


## CREATE ONCOPRINTS
tcgaGenomicOncoprint <- function(df, filename) {
  genomic <- data.frame(Patient=character(),
                        Gene=character(),
                        Alteration=character(),
                        Type=character(),
                        Track=character())
  df$BCAR4_ntile <- ntile(df$BCAR4.fpkm, 20)
  df$ERBB2_ntile <- ntile(df$ERBB2.fpkm, 10)
  for(i in 1:nrow(df)) {
    addedRow=FALSE
    patient <- df[i,"Patient"]
    if (!is.na(df[i,"BCAR4.amp.cbioportal"]) & df[i,"BCAR4.amp.cbioportal"]=="AMP") {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "AMP", "CNA", "BCAR4_Amp")
      addedRow=TRUE
    }
    if (!is.na(df[i,"ERBB2.amp.cbioportal"]) & df[i,"ERBB2.amp.cbioportal"]=="AMP (driver)") {
      genomic[nrow(genomic) + 1,] = c(patient, "ERBB2", "AMP", "CNA", "ERBB2_Amp")
      addedRow=TRUE
    }
    if (!is.na(df[i,"ERBB2.mut.cbioportal"]) & 
        grepl("driver", df[i,"ERBB2.mut.cbioportal"])) {
      genomic[nrow(genomic) + 1,] = c(patient, "ERBB2", df[i,"ERBB2.mut.cbioportal"], "MISSENSE_DRIVER", "ERBB2_Mutations")
      addedRow=TRUE
        }
    if (df[i,"BCAR4_fusion"]=="Fusion") {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "FUSION", "FUSION", "BCAR4_Fusion")
      addedRow=TRUE
    }
    if (df[i,"BCAR4_ntile"] >= 20) {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Top5%")
      addedRow=TRUE
    }
    if (df[i,"BCAR4_ntile"] >= 19) {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Top10%")
      addedRow=TRUE
    }
    if (df[i,"BCAR4_ntile"] >= 17) {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Top20%")
      addedRow=TRUE
    }
    if (df[i,"ERBB2_ntile"] >=10) {
      genomic[nrow(genomic) + 1,] = c(patient, "ERBB2", "HIGH", "EXP", "ERBB2_Top10%")
      addedRow=TRUE
    }
    if (!addedRow) {
      genomic[nrow(genomic) + 1,] = c(patient, "", "", "", "")
    }
  }
  write.table(genomic, file=filename,
              row.names=FALSE, col.names=FALSE, sep="\t",
              quote=FALSE)
}

tcgaClinicalOncoprint <- function(df, filename) {
  clinical <- data.frame(Sample=character(),
                         PAM50=character(),
                         IHC_HER2=character())
  clinical[nrow(clinical) + 1,] = c("sample", "PAM50(string)", "IHC-HER2(string)")
  for(i in 1:nrow(df)) {
    patient = df[i,"Patient"]
    ihc= "N/A"
    if (!is.na(df[i,"IHC.HER2"])) {
      ihc = as.character(df[i,"IHC.HER2"])
      if (ihc != "Positive" & ihc != "Negative") {
        #ihc="Negative" # Added this change to v2
        ihc = "N/A" # v3
      }
    }
    pam50 = "N/A"
    if (!is.na(df[i,"PAM50.cbioportal"])) {
      pam50 = as.character(df[i,"PAM50.cbioportal"])
    }
    clinical[nrow(clinical) + 1,] = c(patient, pam50, ihc)
  }
  write.table(clinical, file=filename,
              row.names=FALSE, col.names=FALSE, sep="\t",
              quote=FALSE)
}

gse6532GenomicOncoprint <- function(df, filename) {
  genomic <- data.frame(Patient=character(),
                        Gene=character(),
                        Alteration=character(),
                        Type=character(),
                        Track=character())
  df$BCAR4_ntile <- ntile(df$BCAR4_exp, 20)
  df$ERBB2_ntile <- ntile(df$HER2_long_exp, 10)
  for(i in 1:nrow(df)) {
    addedRow=FALSE
    patient <- df[i,"Sample"]
    if (df[i,"BCAR4_ntile"] >= 20) {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Top5%")
      addedRow=TRUE
    }
    if (df[i,"BCAR4_ntile"] >= 19) {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Top10%")
      addedRow=TRUE
    }
    if (df[i,"BCAR4_ntile"] >= 17) {
      genomic[nrow(genomic) + 1,] = c(patient, "BCAR4", "HIGH", "EXP", "BCAR4_Top20%")
      addedRow=TRUE
    }
    if (df[i,"ERBB2_ntile"] >=10) {
      genomic[nrow(genomic) + 1,] = c(patient, "ERBB2", "HIGH", "EXP", "ERBB2_Top10%")
      addedRow=TRUE
    }
    if (!addedRow) {
      genomic[nrow(genomic) + 1,] = c(patient, "", "", "", "")
    }
  }
  write.table(genomic, file=filename,
              row.names=FALSE, col.names=FALSE, sep="\t",
              quote=FALSE)
}

gse6532ClinicalOncoprint <- function(df, filename) {
  clinical <- data.frame(Sample=character(),
                         PAM50=character(),
                         er=character(),
                         pgr=character())
  clinical[nrow(clinical) + 1,] = c("sample", "PAM50(string)", "ER(string)", "PGR(string)")
  for(i in 1:nrow(df)) {
    patient = df[i,"Sample"]
    pam50 = "N/A"
    er = "N/A"
    pgr = "N/A"
    if (!is.na(df[i,"PAM50.local"])) {
      pam50 = as.character(df[i,"PAM50.local"])
    }
    if (!is.na(df[i,"er"])) {
      er = ifelse(df[i,"er"]==1,"Positive", "Negative")
    }
    if (!is.na(df[i,"pgr"])) {
      pgr = ifelse(df[i,"pgr"]==1,"Positive","Negative")
    }
  
    clinical[nrow(clinical) + 1,] = c(patient, pam50, er, pgr)
  }
  write.table(clinical, file=filename,
              row.names=FALSE, col.names=FALSE, sep="\t",
              quote=FALSE)
}

tcga <- readRDS("data/TCGA/master.rds")
tcgaGenomicOncoprint(tcga, "data/TCGA/oncoprints/genomic.v1_multip.txt")
tcgaClinicalOncoprint(tcga, "data/TCGA/oncoprints/clinical.v3.txt")

gse6532 <- readRDS("data/GSE6532/master.rds")
gse6532GenomicOncoprint(gse6532, "data/GSE6532/oncoprints/genomic.v1_multip.txt")
gse6532ClinicalOncoprint(gse6532, "data/GSE6532/oncoprints/clinical.v1.txt")
