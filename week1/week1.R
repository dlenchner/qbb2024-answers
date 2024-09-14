library(tidyverse)


# read in the genome_coverage_10x file (Exercise 1.5)
genome_coverage <- read.delim("~/qbb2024-answers/week1/genome_coverage_10x.txt")

# save the genome size (1 Mbp) as a variable
genome_size = 1000000

# determine the maximum coverage of all the positions in the genome (the poisiton in the genome with the most overlapping reads)
max_coverage <- max(genome_coverage)

# calculate the poisson distribution of coverage values ranging from 0 to the maximum coverage value. Change lambda to 10 to account for the 10x coverage
poisson_estimates <- genome_size * dpois(0:max_coverage, 10)

# calculate the normal distribution of coverage values ranging from 0 to the maximum coverage value. Change the mean to 10 and the standard deviation to 3.16 (the sqrt of 10)
normal_estimates <- genome_size * dnorm(0:max_coverage, 10, 3.16)


# plot the histogram and the two distributions using ggplot
# geom_histogram - plotting data from the genome_coverage data frame created when reading in the file
# geom_smooth - each of these is used to plot the indicated distribution (normal or poisson)
ggplot() +
  geom_histogram(data = genome_coverage, mapping = aes(x = number_reads, color = "Genomic Positions per Read Number"), binwidth = 0.5) + 
  geom_smooth(mapping = aes(x = 0:max_coverage, y = poisson_estimates, color = "Poisson Distribution")) +
  geom_smooth(mapping = aes(x = 0:max_coverage, y = normal_estimates, color = "Normal Distribution")) +
  scale_color_manual(name = "Legend", values = c("Genomic Positions per Read Number" = "black", "Poisson Distribution" = "red", "Normal Distribution" = "blue")) + 
  labs(
    title = "Distribution of Coverage Across the Genome",
    subtitle = "Simulated 10x Coverage",
    x = "Number of Reads",
    y = "Number of Positions in the Genome"
  )

# save the plot as a png
ggsave(filename = "~/qbb2024-answers/week1/ex1_10x_cov.png")




