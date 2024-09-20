#!/usr/bin/env python3

# import packages
import sys
import numpy

# create a list of genome reads
reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

# set variable k = 3 (each k-mer has a length of 3 nucleotides)
k = 3

# create an empty DICTIONARY 
graph = {}

# generate the edges for the de Bruijn graph and add them to the empty graph list
# step through each of the 5-mer genome reads in the reads list
for read in reads:
    for i in range(len(read) - k):                              # within each read, iterate at indices 0 and 1 -> range(len(read) - k) = range(len(read) - 3) = range(5-3) = range(2) = [0, 1]
        kmer1 = read[i : i + k]                                 # assign kmer1 to the three nucleotides from index i to index i + k (because python is upper limit exclusionary, adding k = 3 to i will generate a kmer1 that is 3 nt long)
        kmer2 = read[i + 1 : i + 1 + k]                         # assign kmer2 to the three nucletoides shifted over 1 index from those in kmer1
        # print(kmer1 + " -> " + kmer2)            # test the variable assignments above
        graph.setdefault(kmer1 + " -> " + kmer2, 0)             # initialize the graph dictionary with each new edge as the key and "0" as the value
        graph[kmer1 + " -> " + kmer2] += 1                      # add 1 to the value associated with the edge/key
        
# print(graph)              # test the above code


# generate a header for the resulting file                # edited to make the file compatible with graphviz
header = "digraph {"
print(header)


# step through each edge in the graph list
for edge in graph:
    print(edge)                                             # print the edge

# generate a footer for the resulting file
footer = "}"
print(footer)

# Using Unix code, run this program and save the results as a dot file, "week1_deBruijn_edges.dot"







