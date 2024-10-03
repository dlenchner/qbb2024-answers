#!/usr/bin/env python3

# import package
import sys
import numpy

# open the first system argument (in this case, the vcf file) and save as the variable fs
fs = open(sys.argv[1], mode = 'r')

# create a new file "AF.txt" and put it in write mode
af_file = open('AF.txt', mode = 'w')

# create the header for the resulting AF.txt file and write it to the file
header = "allele_frequency"
af_file.write(header + "\n")                        # don't forget the new line character "\n"

# step through each line in the vcf file fs
for line in fs:
    if line.startswith('#'):                        # if the line is a header lines
        continue                                        # skip the line
    fields = line.rstrip("\n").split("\t")          # split line by tabs and remove new line character
    info_column = fields[7].split(";")              # isolate the info column and split the info column further by ";"
    # print(info_column)                    # test the above code to ensure that the correct column from the vcf file was taken and split
    af = info_column[3].lstrip("AF=")               # take the allele frequency AF from the third index in the info column and strip off the "AF="
    af_file.write(af + "\n")                        # write the allele frequency to the AF.txt file (don't forget the new line character "\n")


# close the files
af_file.close()                                    
fs.close()


# reopen the first system argument (in this case, the vcf file) and save as the variable fs
fs = open(sys.argv[1], mode = 'r')


# create a new file "DP.txt" and put it in write mode
dp_file = open('DP.txt', mode = 'w')

# create the header for the resulting DP.txt file and write it to the file
header = "read_depth"
dp_file.write(header + "\n")                        # don't forget the new line character "\n"

# step through each line in the vcf file fs
for line in fs:
    if line.startswith('#'):                        # if the line is a header lines
        continue                                        # skip the line
    fields = line.rstrip("\n").split("\t")          # split line by tabs and remove new line character
    for sample in fields[9:]:                       # step through the FORMAT field for each of the samples
        # print(sample)                         # test the above code to ensure that the FORMAT fields were taken for each sample of each variant
        attributes = sample.split(":")                  # split the FORMAT field for each sample by ":"
        # print(attributes)                     # test the above code
        dp_value = attributes[2]                        # isolate the read depth and save as a new variable dp_value
        # print(dp_value)                       # test the above code
        dp_file.write(dp_value + "\n")                  # write the read depth for each sample and each variant to the DP.txt file
    

# close the files
dp_file.close()
fs.close()
