Exercise 1

Step 1.1

    genome_size = 1000000       # 1 Mbp = 1,000,000 bp
    read_size = 100             # 100 bp per read
    coverage = 3                # 3x coverage of the genome

    perform the calculation for the number of reads needed to cover a 1 Mbp

        read_number = (genome_size * coverage) / read_size

    Answer: 30,000 reads are needed to achieve 3x coverage of a 1 Mbp genome with 100 bp reads



Step 1.2

    See python code in week1.py



Step 1.3

    See R code in week1.R



Step 1.4

    1.  Using the grep function in Unix on the genome_coverage_3x.txt file I made, I found that 49,450 bp were not covered by any reads (roughly 5% of the genome)
    2.  The number of positions lacking coverage (~50,000) matches the expectation based on the Poisson distribution plotted with the histogram.
        For the most part, the normal distribution fits the data really well (the curve for the normal distribution overlaps the histogram bars almost perfectly). However, the normal distribution overestimates the number of samples with zero coverage to be ~70,000.



Step 1.5

    Number of reads: 100,000 reads are needed to achieve 10x coverage of a 1 Mbp genome with 100 bp reads

    1.  Using the grep function in Unix on the genome_coverage_10x.txt file I made, I found that only 23 bp were not covered by any reads
    2.  Again, the number of positions lacking coverage matches the expectation based on the Poisson distribution pretty well, though it is not quite as exact as   was observed in the 3x coverage simulation.
        In this instance, the normal distribution matches the Poisson distribution a bit more closely than was observed previously. However, both distributions stray from the data more than in the 3x coverage example.


Step 1.6

    Number of reads: 300,000 reads are needed to achieve 30x coverage of a 1 Mbp genome with 100 bp reads

    1.  Using the grep function in Unix on the genome_coverage_30x.txt file I made, I found that only 4 bp were not covered by any reads.
    2.  This matches the expectation based on the Poisson distribution, but only because the Poisson distribution loses so much confidence as it tends toward x = 0.
        In this case, the normal and Poisson distributions are almost perfectly aligned and both match the data relatively well (though not as tightly as was seen in the previous exercises).









