library(tidyverse)
library(ggthemes)

# read in the AF.txt file and save as a new dataframe
df <- read_tsv("~/qbb2024-answers/week3/AF.txt")

# create a histogram demonstrating the distribution of alternate allele frequencies in the yeast samples that were sequenced
ggplot(data = df, 
       mapping = aes(x = allele_frequency)) +
  geom_histogram(bins=11, color = "black", fill = "lightgray") + 
  labs(
    title = "Distribution of SNP frequencies in a set of sample yeast sequences", 
    subtitle = "Each sample sequence was aligned to the sacCer3 reference genome",
    x = "Alternate Allele/SNP Frequency",
    y = "Count\n(Number of genomic positions with given allele frequency)"
    )

# save the histogram as a png file
ggsave(filename = "~/qbb2024-answers/week3/week3_allele_frequency_spectrum.png")


### Question 3.1 ###

# Allele frequency appears to follow a normal or Gaussian distribution, 
# centered around a frequency of around 0.5 (more like 0.4). This makes 
# sense given the coin-flipping analogy Rajiv provided in class, whereby 
# the frequency of tails in a series of coin-flipping experiments 
# with a given number of coin flips per experiment, 
# follows a similar distribution (since tails should occur 50% of the 
# time most frequently). In the coin flip 
# analogy, the number of coin flipping experiments is analogous 
# to the number of genome positions examined, the number of coin tosses 
# per experiment is analogous to the number of sequencing reads per 
# genome position, and the occurrence of tails is analogous to the 
# frequency of an alternate allele or SNP at the given position.





