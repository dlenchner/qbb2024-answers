library(tidyverse)


# read in the AF.txt file and save as a new dataframe
df <- read_tsv("~/qbb2024-answers/week3/AF.txt")

# create a histogram demonstrating the distribution of alternate allele frequencies in the yeast samples that were sequenced
ggplot(data = df, 
       mapping = aes(x = allele_frequency)) +
  geom_histogram(bins = 11, color = "black", fill = "lightgray") + 
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


# read in the DP.txt file and save as a new dataframe
df <- read_tsv("~/qbb2024-answers/week3/DP.txt")

# create a histogram demonstrating the distribution of read depths
ggplot(data = df, 
       mapping = aes(x = read_depth)) +
  geom_histogram(bins = 21, color = "black", fill = "lightgray") +
  xlim(0, 20) +
  labs(
    title = "Distribution of coverage across the yeast reference genome",
    x = "Coverage",
    y = "Count\n(Number of genomic positions with given coverage)"
  )

ggsave(filename = "~/qbb2024-answers/week3/week3_read_depth_distribution.png")


### Question 3.2 ###

# Coverage distribution appears to follow a Poisson distribution,
# centered around a coverage of 4x. The modal peak of this histogram 
# at 4x is as expected given the 4x estimated average coverage I 
# calculated with the reads/coverage in the A01_09 sample (see 
# Question 1.3). The shape of the distribution also makes sense when
# you consider that the number of samples (ten) and the number of 
# reads per sample (close to 700,000 for sample A01_09) are small 
# enough to make lower coverages across the genome (like 1-7x) more likely, 
# while making higher coverages (like 15x) possible but not likely.



