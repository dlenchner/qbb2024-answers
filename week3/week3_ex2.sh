#!/usr/bin/env bash


### STEP 2.1 ###

### Question 2.1 ###

# use grep to find each of the chromosome names
grep chr sacCer3.fa

# Total number of chromosomes in sacCer3 reference: 16 nulcear chromosomes plus 1 mitochrondrial chromosome 


### STEP 2.2 ###

bwa index sacCer3.fa


for my_sample in *.fastq
do
    my_sample=`basename ${my_sample} .fastq`
    bwa mem -R "@RG\tID:${my_sample}\tSM:${my_sample}" sacCer3.fa ${my_sample}.fastq > ${my_sample}.sam 
done



### STEP 2.3 ###


# take a look at the SAM file generated for the A01_09 sample
less -S A01_09.sam


### Question 2.2 ###


# use grep to find each of the sequence alignments and bypass the header lines. Pipe this to wc to count the number of lines/alignments
grep HWI A01_09.sam | wc -l

# Number of read alignments = 669548 alignments
    # NOTE: this is the same number of reads in the original A01_09.fastq file



### Question 2.3 ###

# use awk to isolate the column where chromosome number is listed. Pipe this to grep to find each alignment from chrIII. Pipe this to wc to count the number of chrII alignments
awk '{print$3}' A01_09.sam | grep chrIII | wc -l


# Number of alignments to loci on chromosome III = 17815 alignments



### STEP 2.4 ###

# adjusted for loop to include the formatting and indexing of the alignments and saving them as BAM files
for my_sample in *.fastq
do
    my_sample=`basename ${my_sample} .fastq`
    bwa mem -R "@RG\tID:${my_sample}\tSM:${my_sample}" sacCer3.fa ${my_sample}.fastq > ${my_sample}.sam 
    samtools sort -@ 4 -O bam -o ${my_sample}.bam ${my_sample}.sam
    samtools index ${my_sample}.bam
done


### STEP 2.5 ###


### Question 2.4 ###

# The depth of coverage appears to match relatively well with 
# the estimated 4x coverage. There are certainly areas with no 
# coverage, but I found areas with up to 10x coverage, so the 
# average coverage of the whole genome is likely around 4x.



### Question 2.5 ###

# There are two SNPs that I am very confident in, a T>C at position 
# 113,207 of chromosome I and an A>G at position 113,132 chromosome I. I am confident that these 
# are SNPs because all of the reads provided showcase the base change, 
# and there are 4-5 reads covering those particular segments in the genome.

# There is one SNP that I am confident in, but would prefer to see more 
# coverage before trusting it entirely. This is a C>G at position 113,326 chromosome I. 
# This base change is present in every read of the area, but the area only 
# has 2x coverage, so more eveidence would be nice to confirm this read 
# as a SNP.

# Beyond the SNPs mentioned above, there are 12 base changes observed 
# throughout the chromosomal region. However, these are likely not SNPs, 
# since they all only appear in one read of the many that cover those 
# particular regions.



### Question 2.6 ###

# The SNP in this region falls at position 825,834 of chromosome IV, 
# which falls between two genes (scc2 and sas4).







