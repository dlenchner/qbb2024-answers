# WEEK 5 HOMEWORK #


## Exercise 1.1 Question ##

- None of the metrics in the fastqc reports seem particularly troubling to me. Per base sequence content is abnormal (meaning, not the expected 25% per nucleotide) within the first few basepairs of each read, but that is typical given the lack of fidelity that's often observed in sequencing polymerases at the start of a read. Additionally, the distribution of mean GC content has a higher peak than the theoretical distribution, indicating that GC contents at either extreme (close to 100% or close to 0%) are virtually non-existent, but that makes sense given that a relatively equal distribution of As, Ts, C, and Gs allows for the greatest genomic diversity.


## Exercise 1.2 Question ##

- The most overrepresented sequence in the sample originiates in Drosophila and maps to a handful of serine proteases. This makes sense given that proteases are often included in lysis buffers used for nucleic acid purification (to inhibit nucleases).


***


## Exercise 2 Questions ##

- If I were to reject samples with less than 45% unique reads, I would have zero samples to examine (the highest proportion of unique reads is 35.2% in sample SRX298288_2).


- In the DESeq2 map of sample-to-sample distances, the triplicates are clearly distinguishable, meaning each set of replicates is highly consistent.


***



## Exercise 3.3 Questions ##

- No, the third replicate from the LFC-Fe tissue is clustered with the second and third replicates from the Fe tissue, while the first replicate from the Fe tissue is clustered with the first and second replicates from the LFC-Fe tissue.

- The PCA plot suggests that the labels for two of the replicates were swapped (Replicate 1 of Fe tissue and Replicate 3 of LFC-Fe tissue). For an explanation of my answer, see above.


## Exercise 3.6 Question ##

- The results from the GO analysis are't what I initially expected. I expected high expression of genes related to digestion (and there are a few catabolic enzymes expressed), but most of the highly expressed genes are involved in the nervous system. Another large chunk of the highly expressing genes are involved in morphogenesis and development. While these results aren't what I initially expected, they do make sense. The digestive system is heavily intertwined with the nervous system, explaining the neuron-related gene expression. Additionally, this paper was particualrly interested in midgut stem cells in flies, which accounts for the morphogenesis/development-related gene expression.




