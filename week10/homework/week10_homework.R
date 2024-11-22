library(tidyverse)

# read in the file generated from my python script
signals <- read.delim("~/qbb2024-answers/week10/homework/week10_signals.tsv")

# fitler out any NA values
signals <- signals %>%
  filter(!is.na(log2ratio))

# plot the nascent RNA signal per gene as a violin plot
ggplot(data = signals, mapping = aes(x = Gene_ID, y = nascentRNA)) +
  geom_violin() +
  labs(
    title = "Nascent RNA signal per gene",
    x = "Gene",
    y = "Nascent RNA signal"
  )

# save the plot
ggsave(filename = "~/qbb2024-answers/week10/homework/nascentRNA_plot.png")


# plot the PCNA signal per gene as a violin plot
ggplot(data = signals, mapping = aes(x = Gene_ID, y = PCNA)) +
  geom_violin() +
  labs(
    title = "PCNA signal per gene",
    x = "Gene",
    y = "PCNA signal"
  )

# save the plot
ggsave(filename = "~/qbb2024-answers/week10/homework/PCNA_plot.png")


# plot the log2 ratio between nascent RNA signal and PCNA signal per gene as a violin plot
ggplot(data = signals, mapping = aes(x = Gene_ID, y = log2ratio)) +
  geom_violin() +
  labs(
    title = "Log2 ratio of nascent RNA signal to PCNA signal per gene",
    x = "Gene",
    y = "log2 ratio (Nascent RNA/PCNA)"
  )

# save the plot
ggsave(filename = "~/qbb2024-answers/week10/homework/log2ratio_plot.png")





