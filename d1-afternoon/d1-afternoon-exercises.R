#Q1
library(tidyverse)

df <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")




#Q2
glimpse(df)




#Q3
RNA_seq <- df %>% 
  filter(SMGEBTCHT == "TruSeq.v1")




#Q4

#this task can be achieved in a single coding step...
#ggplot(data = RNA_seq, 
#       mapping = aes(x = SMTSD)) +
#  geom_bar()

#below is the same task, but split into two separate coding steps

RNA_seq_tissue <- RNA_seq %>%
  group_by(SMTSD) %>%
  summarize(sample_number = n())

RNA_seq_plot <- ggplot(data = RNA_seq_tissue, 
                       mapping = aes(x = SMTSD, y= sample_number)) +
  geom_bar(stat = "identity")


RNA_seq_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(
    title = "Number of Samples per Tissue Type",
    x = "Tissue",
    y = "Count")



#Q5

glimpse(RNA_seq)

ggplot(data = RNA_seq, 
       mapping = aes(x = SMRIN)) +
  geom_histogram(binwidth = 0.1) +
  labs(
    title = "RNA Integrity Number Distribution",
    x = "RNA Integrity Number",
    y = "Count")


#A histogram is best for visualizing the distribution of these data
#The shape is bimodal, with a major/broader peak in the center of the plot and a second peak right at "10"




#Q6

ggplot(data = RNA_seq, 
       mapping = aes(x = SMRIN, y = SMTSD)) +
  geom_boxplot() +
  labs(
    title = "RNA Integrity Number Distribution per Tissue Type",
    x = "RNA Integrity Number",
    y = "Tissue")

#A boxplot is best for visualizing the distribution of the RIN while separating by tissue type
#Samples taken from cultured cell lines tend to have greater RNA quality (greater RIN)
#Cells - Leukemia cell line (CML)
#Cells - EBV-transformed lymphocytes
#Cells - Cells - Cultured fibroblasts
#The high qualtiy in cultured cells likely stems from the fact that cultured cells are grown in stable and controlled conditions, allowing for collection of large quantities of pure RNA




#Q7

glimpse(RNA_seq)

ggplot(data = RNA_seq,
       mapping = aes(x = SMGNSDTC, y = SMTSD)) +
  geom_boxplot() +
  labs(
    title = "Gene Read Distribution per Tissue Type",
    x = "Number of Gene Reads",
    y = "Tissue")

#As in Q6, a boxplot is best for visualizing distribution of gene reads while separating by tissue type
#In general, there aren't any major differences between tissues
#However, the testis is an outlier with a much greater number of gene reads than there are in the other tissue types
#As per Xia et al., the greater-than-normal gene read number in the testis may be due to the fact that the testis expresses a greater number of genes than any other organ




#Q8

isch_RIN_plot <- ggplot(data = RNA_seq,
                        mapping = aes(x = SMTSISCH, y = SMRIN)) +
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(. ~ SMTSD) +
  geom_smooth(method = "lm") +
  labs(
    title = "Ischemic Time vs RNA Integrity Number per Tissue Type",
    x = "Ischemic Time",
    y = "RNA Integrity Number")

isch_RIN_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


#In many tissues, there doesn't appear to be a correlation between Ischemic Time and RIN
#In certain tissues (Vagina, Uterus, Liver, Lung, etc.), there is a negative correlation between the variables, where RIN decreases as Ischemic Time increases




#Q9

isch_RIN_plot_2 <- ggplot(data = RNA_seq,
                          mapping = aes(x = SMTSISCH, y = SMRIN)) +
  geom_point(aes(color = SMATSSCR), 
             size = 0.5, 
             alpha = 0.5) +
  facet_wrap(. ~ SMTSD) +
  geom_smooth(method = "lm") +
  labs(
    title = "Ischemic Time vs RNA Integrity Number per Tissue Type",
    x = "Ischemic Time",
    y = "RNA Integrity Number",
    color = "Autolysis Score"
  )

isch_RIN_plot_2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#In many tissue types, a greater autolysis score correlates to both a higher Ischemic Time and a lower RIN
#This relationship is tissue-dependent. The trend is observed in tissues like the thyroid and transverse colon, but not in tissues like the skeletal muscle

