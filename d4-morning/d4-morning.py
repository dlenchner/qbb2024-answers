#!/usr/bin/env python3

import sys

import numpy

# Exercise 1 - maintained the original dictionary that uses the (gene, tissue) tuple as keys

# get target file, name the file, and open the file
filename = sys.argv[1]
fs = open(filename, mode = 'r')

# create an empty dictionary to hold samples for gene-tissue pairs
relevant_samples = {}

# step through file
for line in fs:
    fields = line.rstrip("\n").split("\t")      # split line into fields
    key = (fields[0], fields[2])                # create key from gene and tissue
    relevant_samples[key] = []                  # initialize dictionary from key with list to hold samples

fs.close()


# print(relevant_samples)                       # test the above code



# Exercise 2

# get metadata file, name the file, and open the file
filename = sys.argv[2]
fs = open(filename, mode = 'r')

# skip the first line
fs.readline()

# create an empty dictionary to hold samples for sample-tissue pairs
tissue_samples = {}

# step through file
for line in fs:
    fields = line.rstrip("\n").split("\t")      # split line into fields
    key = fields[6]                             # create key from tissue
    value = fields[0]                           # create value from sampleID
    tissue_samples.setdefault( key, [] )        # initialize dictionary from key with list to hold sampleIDs
    tissue_samples[key].append(value)           # using the key created, add the sampleID to the value list in the dictionary
             
fs.close()

# print(tissue_samples)                         # test the above code                


# Exercise 3-5 - these steps were covered by Mike in lecture and compiled into this streamlined method


# get raw data file, name the file, and open the file
filename = sys.argv[3]
fs = open(filename, mode = 'r')

# skip the first two lines
fs.readline()
fs.readline()

# take the third line in the file (the column names), strip the new line character and split by tabs. Assign to a new variable header
header = fs.readline().rstrip("\n").split("\t") 

header = header[2:]

# print(header)                                 # test the above code

# create an empty dictionary to hold indices associated with particular tissues in the raw data file
tissue_columns = {}

for tissue, samples in tissue_samples.items():      # use the items function to return both the keys and values from the tissue_samples dictionary. Store the keys as the variable tissue and the values as the variable samples
    tissue_columns.setdefault(tissue, [])           # in the empty tissue_columns dictionary, any time a new tissue is encountered, give it a value of a blank list
    for sample in samples:                          # step through each sample in the list of samples associated with the given tissue
        if sample in header:                        # check if each sample is in the header created in line 63. When that condition is satisfied -->
            position = header.index(sample)         # save the index as the variable position
            tissue_columns[tissue].append(position) # append the position (index) to the list created (line 71) for that given tissue

# fs.close()

# # print(tissue_columns)                             # test the above code



# # get raw data file, name the file, and open the file
# filename = sys.argv[3]
# fs = open(filename, mode = 'r')

# # skip the first three lines
# fs.readline()
# fs.readline()
# fs.readline()                                       # NOTE: lines 77-90 are unnecessary. The code runs fine (and in fact, a bit faster) if the program doesn't need to close and reopen the raw data file. I'm choosing to leave this extra code in place to show my initial logic.


# Exercise 6-7 - I combined these steps given the nested for loops shown below

# create the header for the resulting tsv file
print("Expression", "Gene_ID", "Tissue")

# step through the file
i = 0
for line in fs:
    fields = line.rstrip("\n").split("\t")                                          # split line into fields
    gene_ID = fields[0]                                                             # assign the value in the first column to the variable gene_ID
    for gene, tissue in relevant_samples.keys():                                    # use the keys function to return the keys from the relevant_samples dictionary. Store the gene in each key as the variable gene and the tissue in each key as the variable tissue
        if gene_ID == gene:                                                         # check if the gene_ID from the raw data file matches the gene from the relevant samples dictionary. When that condition is satisfied -->
            expression_values = []                                                  # create an empty list for the expression_values associated with the gene
            expression_values.append(fields[2:])                                    # append the expression values associated with the gene in the raw data file to the expression_values list
            expression_values_array = numpy.array(expression_values, dtype = float) # turn the expression_values list into a numpy array, specifying the data type as float
            sample_indices = tissue_columns.get(tissue)                             # use the get function on the tissue_columns dictionary to retrieve the indices (values) of interest for the tissue (key) associated with the gene being examined
            # print(tissue)                                                         # test the above code
            # print(sample_indices)                                                 # test the above code
            # print(expression_values_array)                                        # test the above code
            expression_val_interest = list(expression_values_array[0][sample_indices])  # from the expression_values_array, extract/list the expression values corresponding to the sample indices retrieved from the tissue_columns dictionary
            # print(expression_val_interest)                                        # test the above code
            for index in expression_val_interest:                                   # step through each index in the newly generated expression_val_interest list
                print(expression_val_interest[i], gene, tissue)                     # beneath the header printed in line 94, print the expression value at index i, the gene being examined, and the paired tissue
                i += 1                                                              # add 1 to i to continue stepping through the indices of expression_val_interest
            i = 0                                                                   # reset i to 0
                 
fs.close()




# Using Unix code, run this program and save the results as a tsv file, titled "genes_tissues_expression.tsv"






