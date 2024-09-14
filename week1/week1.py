#!/usr/bin/env python3

import sys
import numpy

# Exercise 1
# Step 1.1

# assign variables to each of the factors needed to calculate the number of reads (genomoe size, the size of each read, and the coverage we're looking for)
genome_size = 1000000       # 1 Mbp = 1,000,000 bp
read_size = 100             # 100 bp per read
coverage = 3                # 3x coverage of the genome

# perform the calculation for the number of reads needed to cover a 1 Mbp
read_number = (genome_size * coverage) / read_size

# print(read_number)            # test the above code to calculate the number of reads needed

# Answer: 30,000 reads are needed to achieve 3x coverage of a 1 Mbp genome with 100 bp reads


# Exercise 1.2

genome_coverage = numpy.zeros(genome_size, int)

# print(genome_coverage)        # test the above code


# randomly simulate sequencing the genome and count the number of times each position is covered
# for each of the 30,000 iterations (reads) calculated above
for i in range(int(read_number)):
    start_position = numpy.random.randint(0, genome_size - read_size + 1)           # randomly assign a starting position that could fall anywhere from the 0th position to the 999,900th position of the genome_coverage array (because python is upper limit exclusionary) - genome position 1 to 999,901 
    end_position = start_position + read_size                                       # assign an ending position by adding the read size (100 bp) to the starting position assigned above
    genome_coverage[start_position:end_position] += 1                               # for every index of the genome_coverage array covered by the random read, add 1 to the count displayed in the array 
    

# print the number of times each position is covered by the random reads above
# create a header for the resulting read numbers to help when reading the file into R
print("number_reads")

# step through each position in the genome_coverage array
for position in genome_coverage:
    print(position)                                                                 # print the value at each position in the array, which represents the number of times that index in the genome was covered by a simulated sequence



# Using Unix code, run this program and save the results as a txt file


