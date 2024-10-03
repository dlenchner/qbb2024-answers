#!/usr/bin/env bash


### Question 1.1 ###



# view the first alignment and perform a character count on the longest line of data (which is the read length)
head -4 A01_09.fastq | wc -L -m


# Read lengths are 76 bp



### Question 1.2 ###

# determine the total number of lines in the document
less -s A01_09.fastq | wc -l

# Total number of lines = 2678192

# divide the total number of lines by 4 (since each read has four lines associated with it - sequence name, sequence, "+", quality score)
bc -e 2678192/4


# Number of reads in the file = 669548 reads



### Question 1.3 ###


# calculate the total number of base pairs covered in the fastq file
# NOTE: the '\' is required before the '*' so the program can identify it as a multiplicative character, rather than a wild card character
bc -e 669548\*76       


# Toal number of base pairs covered in the fastq file = 50885648
# Size of S. cerevisiae reference genome = 12.1 Mb

# calculate the approzimate coverage based on the number of base pairs covered in the fastq file and the known size of the S. cerevisiae genome
bc -S 3 -e 50885648/12.1e6


# Approximate coverage of S. cerevisiae genome = 4x coverage




### Question 1.4 ###

# determine the size of each file and sort based on file size
du *.fastq | sort

# Largest file = A01_62.fastq (304352 bytes)
# Smallest file = A01_27.fastq (224512 bytes)



### Question 1.5 ###

# run the FastQC program on the A01_09.fastq file 
fastqc A01_09.fastq

# open the HTML report and examine the results
open A01_09_fastqc.html


# Median base quality along the read: 35 or 36
# The high score along the read indicates that each given base is likely not an error (close to a 1 in 10^-4, or 0.01%, chance that each base call is an error)
# For this particular read, there is very little variation in quality with respect to position














