#!/usr/bin/env python3

# import package
import sys

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

fs.close()
af_file.close()
