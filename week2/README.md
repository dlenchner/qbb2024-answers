
Exercise 1.4

    The following command lines demonstrate the arguments I entered to sort (bedtools sort -i), merge (beftools merge), and save (> <feature>_chr.bed) each of the feature files

        bedtools sort -i genes.bed | bedtools merge > genes_chr1.bed

        bedtools sort -i exons.bed | bedtools merge > exons_chr1.bed

        bedtools sort -i cCREs.bed | bedtools merge > cCREs_chr1.bed


Exercise 1.5 

    The following command line demonstrates the arguments I entered to subtract the labeled regions in the exons feature file from the labeled regions in the genes feature file to obtain an introns feature file

        bedtools subtract -a genes_chr1.bed -b exons_chr1.bed > introns_chr1.bed


Exercise 1.6

    The following command line demonstrates the arguments I entered to subtract the regions from each of the feature files from the genome file to obtain a feature file that highlights regions in chromosome 1 that are not listed as exons, introns, or regulatory regions

        bedtools subtract -a genome_chr1.bed -b exons_chr1.bed -b introns_chr1.bed -b cCREs_chr1.bed > other_chr1.bed



Exercise 2.3

    1. Based on the SNP enrichment data, exons appear to be under the most stringent purifying selection. This is indicated by the considerably lower SNP enrichments observed across all MAF values in the exons relative to the introns, cCREs, and other genome features.

    2. The general decrease in enrichment observed as MAF value increases is likely due to purifying selection and genomic maintenance mechanisms. MAF reflects the frequency with which an SNP/alternative allele is present in a population. One would expect SNPs with high MAFs to be less frequent/less enriched in the genome given that cellular machinery and selective pressures work to prevent mutations from taking hold in a population.

    3. As mentioned in question 1, exons display the lowest SNP enrichment of any genomic feature. Among the other genomic features assessed, introns appear to have the highest SNP enrichment, though the SNP enrichments in cCREs and other features are very similar. This is an interesting finding as it suggests that introns and regulatory elements in particular are more prone to nucleotide changes than exons. Upon further thought, this trend makes sense, as SNPs in the coding exonic regions of the genome are more likely to alter the functional protein products of a gene (and so are more likely to be selected against and are less enriched), whereas SNPs in the non-coding intronic and regulatory regions of the genome are less likely to impact protein function (and so face less selective pressure and are more enriched).    