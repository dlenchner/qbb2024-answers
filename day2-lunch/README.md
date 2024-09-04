# day2-lunch answers

## Exercise 1-1

- `cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c`
- There are 19,618 `protein_coding` genes

- I would like to learn more about the miRNA biotype. I studied the influence of microRNAs in C. elegans stress response, so I have experience and a lot of interest in miRNA genetics.

## Exercise 1-2

- `cut -f 1 hg38-gene-metadata-go.tsv | sort | uniq -c | sort -n`
- The gene ENSG00000168036 has the greatest number of `go_ids` (273)

- `grep -w "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -k3 -f | less -S > gene_idENSG00000168036-go_id-sorted.txt`
- After examining some of the gene ontology descriptions, it seems that gene ENSG00000168036 is involved in cell/tissue differentiation, morphogenesis, and adhesion across the human body 

## Exercise 2-1

- `grep "IG_" genes.gtf | grep -v "_pseudogene" | cut -f 1 | uniq -c`
- The following chromosomes contain IG Genes:
    - Chromosome 2: 52
    - Chromosome 14: 91
    - Chromosome 15: 16
    - Chromosome 16: 6
    - Chromosome 21: 1


- `grep "IG_" genes.gtf | grep "_pseudogene" | cut -f 1 | uniq -c`
- The following chromosomes contain IG Pseudogenes:
    - Chromosome 1: 1
    - Chromosome 2: 45
    - Chromosome 8: 1
    - Chromosome 9: 5
    - Chromosome 10: 1
    - Chromosome 14: 84

- When comparing distribution between IG Genes and IG Pseudogenes, there is tremendous overlap. Most of the IG Genes are on chromosomes 2, 14, and 15, and most of the IG Pseudogenes are on chromosomes 2 and 14


## Exercise 2-2

- Using `grep pseudogene genes.gtf` doesn't work to identify all pseudogenes because this technique does not specify to search for "pseudogene" within the "gene_type" descriptor. Instead, it searches for the character string p-s-e-u-d-o-g-e-n-e anywhere in the file. This is demonstrated with the transcripts the are tagged as overlapping a pseudogene (they are not psuedogenes, but they overlap psuedogenes)
- A better method to search for pseudogenes would be to grep for the phrase "processed_pseudogene". This would provide a list of genes where the gene_type is described as either processed_pseudogene OR unprocessed_pseudogene (since the "processed_pseudogene" search term is not restricted to search for the whole word/phrase)


## Exercise 2-3

- `cut -f 1,4-5,14 genes-tabs.gtf > genes-tabs.bed`

