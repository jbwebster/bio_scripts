


setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline")

# This version is what I originally shared with the group and were hypermutated samples
# were originally detected
originals = read.csv("data/tulane/snv/cleaned-tulane_validation.withHypermutation.tsv",
                     header=T, sep="\t")
hmrerun = read.csv("data/tulane/snv/cleaned-tulane_hmrerun.tsv",
                   header=T, sep="\t")

source <- originals %>%
  group_by(Source) %>%
  summarise(Count = n())

ggplot() +
  geom_col(data=source, aes(x=reorder(Source, -Count), y=Count)) +
  labs(x="Sample",
       y="# of SNVs",
       caption="Counting SNVs in genes of interest, passing all filters 
       and labeled as probably_damaging",
       title="Original detection of 'hypermutated' samples") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

hmrerun_source <- hmrerun %>%
  group_by(Source) %>%
  summarise(Count = n())
source$isHypermutated <- ifelse(source$Count>=50, "Yes", "No")  
source$wasHypermutated <- source$isHypermutated
was_hypermutated <- c(subset(source, isHypermutated=="Yes")[,1])
hmrerun_source$isHypermutated <- ifelse(hmrerun_source$Source, "Yes", "No")
hmrerun_source$wasHypermutated <- ifelse(hmrerun_source$Source %in% was_hypermutated$Source, "Yes", "No")
m <- rbind(source, hmrerun_source)


ggplot() +
  geom_col(data=subset(m, isHypermutated=="No"), aes(x=reorder(Source, -Count), y=Count, fill=wasHypermutated)) +
  labs(x="Sample",
       y="# of SNVs",
       caption="Counting SNVs in genes of interest, passing all filters 
       and labeled as probably_damaging",
       title="Reanalysis of 'hypermutated' samples") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
