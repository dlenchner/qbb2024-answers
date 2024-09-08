library(tidyverse)
library(broom)

# Exercise 1.1 and 1.3

# read in the data
dnm <- read_csv(file = "~/qbb2024-answers/d4-afternoon/aau1043_dnm.csv")

ages <- read_csv(file = "~/qbb2024-answers/d4-afternoon/aau1043_parental_age.csv")


view(dnm)
view(ages)

# Exercise 1.2

# calculate and organize the number of maternal and paternal dnms, grouping by Proband_id
dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))


# Exercise 1.4

# combine the two datasets into one, using the Proband_id to attach each parent age to the correct proband
dnm_by_parental_age <- left_join(dnm_summary, ages, by = "Proband_id")



# Exercise 2.1

# Generate the maternal dnm vs. maternal age plot

ggplot(data = dnm_by_parental_age, 
       mapping = aes(x = Mother_age, 
                     y = n_maternal_dnm)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  labs(
    title = "Maternal de novo Mutations vs. Age of Mother",
    x = "Age of Mother",
    y = "Maternal de novo Mutations"
  )


# Generate the paternal dnm vs. paternal age plot

ggplot(data = dnm_by_parental_age, 
       mapping = aes(x = Father_age, 
                     y = n_paternal_dnm)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  labs(
    title = "Paternal de novo Mutations vs. Age of Father",
    x = "Age of Father",
    y = "Paternal de novo Mutations"
  )



# Exercise 2.2

# generate a linear regression model for the maternal age vs. maternal dnm data

maternal_model <- lm(data = dnm_by_parental_age,
                     formula = n_maternal_dnm ~ 1 + Mother_age)

summary(maternal_model)

# Questions:
  # The relationship between maternal age and the number of 
      # maternal dnms is relatively small (slope = 0.37757).
      # This means that age of the mother does positively 
      # correlate to maternal dnm number, but as age of the 
      # mother increases, the occurrence of maternal dnms 
      # doesn't increase very much. Roughly, for every three 
      # years added to the mother's age, the offspring is 
      # estimated to have 1 additional dnm. Overall, this 
      # matches the visual expectation I had after viewing 
      # the plot.
  # Regardless of the "size" of the relationship between 
      # maternal age and maternal dnm number, the relationship 
      # is significant, as indicated by the small p-value of 
      # < 2e-16. This means that we can have high confidence 
      # in the association between maternal age and maternal
      # dnm number


# Exercise 2.3

# generate a linear regression model for the paternal age vs. paternal dnm data

paternal_model <- lm(data = dnm_by_parental_age,
                     formula = n_paternal_dnm ~ 1 + Father_age)

summary(paternal_model)

# Questions:
  # The relationship between paternal age and the number 
      # of paternal dnms is much larger than was observed 
      # when comparing maternal age to maternal dnm counts. 
      # The larger slope of 1.35384 here suggests that for 
      # every 3 years added to the age of the father, the 
      # offspring is estimated to have an additional 4 dnms. 
      # Visually, this model matches the the trends observed
      # in the plot 
  # The relationship between paternal age and paternal dnm 
      # count is significant, indicated by the low p-value of 
      # < 2e-16. This means that we can have high confidence 
      # in the association between paternal age and paternal
      # dnm number



# Exercise 2.4


# Assign a new variable (new_father) for the 50.5 year-old father in a tibble
new_father <- tibble(
  Father_age = 50.5) 

# Run the predict function on the new_father using the paternal_model to 
    # calculate the estimated response variable (paternal dnm count)
predict(paternal_model, newdata = new_father)

# The proband of the new father is estimated to have 78 paternal dnms




# Exercise 2.5

ggplot() +
  geom_histogram(data = dnm_by_parental_age, mapping = aes(x = n_maternal_dnm, fill = "Maternal DNMs"), alpha = 0.5) + 
  geom_histogram(data = dnm_by_parental_age, mapping = aes(x = n_paternal_dnm, fill = "Paternal DNMs"), alpha = 0.5) +
  labs(
    title = "Distribution of Parental DNMs",
    x = "Number of DNMs",
    y = "Count",
    fill = NULL)




# Exercise 2.6


# Assign the maternal dnm and paternal dnm columns from the combined dataset to their own variables
n_maternal_dnm <- dnm_by_parental_age$n_maternal_dnm
n_paternal_dnm <- dnm_by_parental_age$n_paternal_dnm

# Run a Welch's t-test to compare the paternal and maternal dnm count, assuming unequal variance
t.test(n_paternal_dnm, n_maternal_dnm, alternative = "two.sided", var.equal = FALSE)


# I chose a Welch's t-test because the variances of maternal and paternal dnms are not equal 
    # (paternal dnms have greater variance)

# Yes, the result was significant, indicated by the small p-value of < 2.2e-16. 
    # Given the distribution of the data observed in the histogram in Exercise 
    # 2.5, this suggests that children are likely to inherit significantly more 
    # DNMs from their father than their mother.















