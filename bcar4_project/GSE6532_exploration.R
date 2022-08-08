
library(ggplot2)
library(scales)



setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/BCAR4")





master <- readRDS("data/GSE6532/master.rds")



# PAM50 Subtyping Results
ggplot(master, aes(x=Treatment, fill=PAM50.local)) +
  geom_bar() +
  labs(x="Treatment Group",
       y="Sample Count",
       title="Generated PAM50 Status",
       caption="GSE6532")


# Which HER2 probe to use?
ggplot(master, aes(x=HER2_long_exp, y=HER2_short_exp)) +
  geom_point() +
  labs(x="Shorter Probe Expression",
       y="Longer Probe Expression",
       title="Short vs Long Probe Expression Scores - HER2")

df <- master
df <- subset(master, master$PAM50.local!="Normal")
df <- subset(df, df$Treatment=="Treated")
p_meds <- ddply(df, .(PAM50.local), summarise, med = round(median(BCAR4_exp), digits=3))
p_counts <- ddply(df, .(PAM50.local), summarise, counts = n())
ggplot(df, aes(x=PAM50.local,y=BCAR4_exp)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(outlier.shape=NA, alpha=0.9, aes(fill=PAM50.local)) +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(data = p_meds, aes(x = PAM50.local, y = 4, label = med), 
            size = 3, vjust = -1.5) +
  geom_text(data = p_counts, aes(x = PAM50.local, y = 3.5, label = paste0("n=",counts)), 
            size = 3, vjust = -1.5) + 
  labs(x="PAM50",
       y="Normalized BCAR4",
       title="BCAR4 Expression by PAM50",
       caption="GSE6532 - Treated Samples Only")

