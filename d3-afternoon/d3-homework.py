#!/usr/bin/env python3


import sys

import numpy

# open file
fs = open(sys.argv[1], mode = 'r')

# skip 2 lines
fs.readline()
fs.readline()

# split column header by tabs, remove new line character, and skip first two entries
line = fs.readline()
fields = line.strip("\n").split("\t")
tissues = fields[2:]

# create way to hold gene names
# create way to hold gene IDs
# create way to hold expression values
gene_names = []
gene_IDs = []
expression = []

# for each line
for line in fs:
    fields = line.strip("\n").split("\t")   # split line and remove new line character
    gene_IDs.append(fields[0])              # save field 0 (column 1) into gene IDs
    gene_names.append(fields[1])            # save field 1 (column 2) into gene names
    expression.append(fields[2:])           # save field 2+ (columns 3 and on) into expression values

# close file
fs.close()

# translate each of the lists into arrays
tissues = numpy.array(tissues)
gene_IDs = numpy.array(gene_IDs)
gene_names = numpy.array(gene_names)
expression = numpy.array(expression, dtype = float)

# Why do you need to tell numpy what type of data the expression data is but not for the other data lists?
    # Each of the elements in each of the lists is a string. Expression data, visually, are numbers, but numpy does not recognize the numbers as dtype = float. 
    # So, we need to tell numpy that we want it to read each of the strings in the expression data list as a float.

# Question 4

# Calculate the mean (using numpy.mean) for the first 10 rows, but every column (expressed as [:10,:]), along the horizontal axis (axis = 1). 
# Save this calculation as variable expression_10_means.
expression_10_means = numpy.mean( expression[:10,:], axis = 1 )

#print(expression_10_means)


# Question 5

expression_means = numpy.mean( expression )
expression_medians = numpy.median( expression )

# mean expression = 16.557814350910945
# median expression = 0.0271075

# The low median suggests that a large chunk of the genes are very lowly expressed in a large subset of tissues. 
# The larger mean, however, suggests that a few of the genes have incredibly high gene expressions in certain tissue types.


# Question 6

# Add a pseudocount of 1 to each value in the expression array
expression_pseudocount = expression + 1

# Take the log base 2 of each value in the array
expression_log2 = numpy.log2( expression_pseudocount )

# Take the mean and median of the normalized values in the array expression_log2
expression_log2_means = numpy.mean( expression_log2 )
expression_log2_medians = numpy.median( expression_log2 ) 

# normalized mean expression = 1.1150342022364093
# normalized median expression = 0.03858718613570538


# Now, the mean and median are much closer to equal. However, the mean is still over 30-fold greater than the median, which reflects the same trend observed in the non-transformed values.
# Because of the log normalization, the median wasn't changed very much (0.039 instead of 0.027), but the mean was changed rather dramatically (1.12 instead of 16.56).


# Question 7

# Sort the log-adjusted expression array
expression_log2_sorted = numpy.sort(expression_log2, axis = 1)

# Generate a new 1D array (diff_array) by subtracting the second-to-last value in each row of the expression_log2_sorted array from the last value in each row of that array.
diff_array = expression_log2_sorted[ :, -1 ] - expression_log2_sorted[ :, -2]



# Question 8

# Threshold the array to filter out values that are less than 10. 
# Generate a new array whose values are the positions of values above the threshold in the original array.
diff_array_thresholded = numpy.where( diff_array >= 10 )

# Calculate the size of the thresholded array
high_expressing_genes = numpy.size(diff_array_thresholded)


# Print the number of high-expressing genes, calculated with the previous size method 
print(high_expressing_genes)

# There are 33 high-expressing genes