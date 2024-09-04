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

