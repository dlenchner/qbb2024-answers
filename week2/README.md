
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

        bedtools subtract -a genome_chr1.bed -b exons_chr1.bed -c introns_chr1.bed -d cCREs_chr1.bed > other_chr1.bed


