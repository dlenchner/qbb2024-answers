#!/usr/bin/env python3

import sys

import numpy

# get target file name and open the file
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


# print(relevant_samples)





# get metadata file name and open the file
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
    tissue_samples.setdefault( key, [] )        # initialize dictionary from key with list to hold samples
    tissue_samples[key].append(value)        
             
fs.close()

# print(tissue_samples)


# get raw data file name and open the file
filename = sys.argv[3]
fs = open(filename, mode = 'r')

# skip the first two lines
fs.readline()
fs.readline()

# take the third line in the file (the column names), strip the new line character and split by tabs. Assign to a new variable "header"
header = fs.readline().rstrip("\n").split("\t") 

header = header[2:]


# print(header)

tissue_columns = {}

for tissue, samples in tissue_samples.items():      # the items function returns both the keys and values from the tissue_samples dictionary and stores the keys as the variable tissue and the values as the variable samples
    tissue_columns.setdefault(tissue, [])           # in the empty tissue_columns dictionary, any time a new tissue is encountered, give it a value of a blank list
    for sample in samples:                          # step through each sample in the list of samples associated with the given tissue
        if sample in header:
            position = header.index(sample)
            tissue_columns[tissue].append(position)

fs.close()

# print(tissue_columns)



# get target file name and open the file
filename = sys.argv[3]
fs = open(filename, mode = 'r')

fs.readline()
fs.readline()
fs.readline()

print("Expression", "Gene_ID", "Tissue")

i = 0
for line in fs:
    fields = line.rstrip("\n").split("\t")      # split line into fields
    gene_ID = fields[0]
    for gene, tissue in relevant_samples.keys():
        if gene_ID == gene:
            expression_values = []
            expression_values.append(fields[2:])
            expression_values_array = numpy.array(expression_values, dtype = float)
            sample_indices = tissue_columns.get(tissue)
            # print(tissue)
            # print(sample_indices)
            # print(expression_values_array)
            expression_val_interest = list(expression_values_array[0][sample_indices])
            # print(expression_val_interest)
            for index in expression_val_interest:
                print(expression_val_interest[i], gene, tissue)
                i += 1
            i = 0
                 
fs.close()




# Using Unix code, run this program and save the results as a tsv file, titled "genes_tissues_expression.tsv"






