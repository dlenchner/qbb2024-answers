library(tidyverse)

df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") #assign the text file data to the variable "df"

df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 ) #generate a SUBJECT column

#Arrange the SUBJECTs in ascending and descending order by "sample_number"
df %>% #search within the dataset "df"
  group_by(SUBJECT) %>% #group the data in "df" by SUBJECT
  summarize( sample_number=n() ) %>% #summarize the number of data points associated with each SUBJECT - variable name changed from "n()" to "sample_number"
  arrange(sample_number) %>% #arrange "sample_number in ascending order
  arrange(-sample_number) #arrange "sample_number" in descending order ("-" before the variable name reverses order)

#Two SUBJECTS with most samples: K-562 (217 samples) and GTEX-NPJ8 (72 samples)
#Two SUBJECTS with fewest samples: GTEX-1JMI6 (1 sample) and GTEX-1PAR6 (1 sample)


#Arrange the SMTSDs (tissue types) in ascending and descending order by "sample_number"
df %>%
  group_by(SMTSD) %>%
  summarize( sample_number=n() ) %>%
  arrange(sample_number) %>%
  arrange(-sample_number)

#Two SMTSDs with most samples: Whole Blood (3288 samples) and Muscle-Skeletal (1132 samples)
  #WHY: Easy to collect
#Two SMTSDs with fewest samples: Kidney-Medulla (4 samples) and Cervix-Ectocervix (9 samples)
  #WHY: Difficult to collect


df_npj8 <- filter(df, SUBJECT == "GTEX-NPJ8") #filter/retain data from SUBJECT "GTEX-NPJ8" within the dataset "df" and store the new dataset as "df_npj8"


#Arrange the SMTSDs (tissue types) from SUBJECT "GTEX-NPJ8" in descending order by "sample_number_npj8"
df_npj8 %>% #search within the dataset "df_npj8"
  group_by(SMTSD) %>% #group the data in "df_npj8" by SMTSD (tissue type)
  summarize( sample_number_npj8=n()) %>% #summarize the number of data points associated with each SMTSD - variable name changed from "n()" to "sample_number_npj8"
  arrange(-sample_number_npj8) #arrange "sample_number_npj8" in descending order

#Tissue with the most samples: Whole Blood (9 samples)

view(df_npj8)  

#Within the "Whole Blood" samples, what is different between the samples: Nucleic acid collection method and sequencing method

#Filter out any SMATSSCR value of "NA" from the "df" dataset
df %>%
  filter( !is.na(SMATSSCR) )















