#!/bin/bash


# Generate a new file and give it a header separated by tabs (\t)
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt 

# step through each of the MAF values
for i in 0.1 0.2 0.3 0.4 0.5
    do                                                                                      # initiate the for loop 
    bedtools coverage -a genome_chr1.bed -b chr1_snps_${i}.bed > snp_chr1_coverage.txt      # determine the coverage of the SNP file corresponding to the MAF value i relative to chromosome 1 in the chr1 bed file. Assign to a new temporary file "snp_chr1_coverage.txt"
    snp_chr1_count=$(awk '{s+=$4}END{print s}' snp_chr1_coverage.txt)                       # in the new file generated from the coverage function in line 10, add the number of SNPs and save as the variable "snp_chr1_count"
    chr1_size=$(awk '{s+=$6}END{print s}' snp_chr1_coverage.txt)                            # also in the new file generated in line 10, determine the total size of chromosome 1 and save as the variable "chr1_size"
    background=$(bc -l -e "${snp_chr1_count} / ${chr1_size}")                               # calculate the ratio of snp_chr1_count relative to chr1_size and save as the variable "background"

    # within the first for loop, step through each of the genome features
    for j in exons introns cCREs other
        do                                                                                  # initiate the second for loop
        bedtools coverage -a ${j}_chr1.bed -b chr1_snps_${i}.bed > snp_feat_coverage.txt    # determine the coverage of the SNP file corresponding to the MAF value i relative to the feature file corresponding to the feature j. Assign to a new temporary file "snp_feat_coverage.txt"
        snp_feat_count=$(awk '{s+=$4}END{print s}' snp_feat_coverage.txt)                   # in the new file generated from the coverage function in line 18, add the number of SNPs and save as the variable "snp_feat_count"
        feat_size=$(awk '{s+=$6}END{print s}' snp_feat_coverage.txt)                        # also in the new file generated in line 18, determine the toal size of the feature j and save as the variable "feat_size"
        snp_density=$(bc -l -e "${snp_feat_count} / ${feat_size}")                          # calculate the ratio of snp_feat_count relative to feat_size and save as the variable "snp_density"
        enrichment=$(bc -l -e "${snp_density} / ${background}")                             # calculate the enrichment value by dividing snp_density by the "background" calculated in the initial for loop. Save as the variable "enrichment" 
        echo -e "${i}\t${j}\t${enrichment}"                                                 # pass the information derived from this nested for loop (MAF value i, genome feature j, and enrichment value) to the next function (separate the information by tabs)
        done                                                                                # close the second for loop
    done >> snp_counts.txt                                                                  # append the information that was passed with the echo function (line 23) from each iteration into the file generated in line 5 (snp_counts.txt)


# remove the temporary files generated within the for loop from my directory (NOTE: only the final iteration of the for loop was saved to each file name)
rm snp_chr1_coverage.txt                                                                    
rm snp_feat_coverage.txt                                                      