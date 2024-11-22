#!/usr/bin/env python

import numpy
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly
import plotly.express as px


genes = ['APEX1', 'PIM2', 'POLR2B', 'SRSF1']
fields = ['field0', 'field1']
channels = ['DAPI', 'PCNA', 'nascentRNA']

img = imageio.v3.imread("APEX1_field0_DAPI.tif")

images = []

for gene in genes:
    for field in fields:
        array = numpy.zeros((img.shape[0], img.shape[1], 3), numpy.uint16)
        for i, channel in enumerate(channels):
            array[:, :, i] = imageio.v3.imread(f"{gene}_{field}_{channel}.tif")
        images.append(array)

plt.imshow(array)
plt.show()


mask = []
for i in range(len(images[:])):
    mask.append(images[i][:, :, 0] >= numpy.mean(images[i][:, :, 0]))


label_array = []

# Let's use a function to find nuclei bodies
def find_labels(mask):
    # Set initial label
    l = 0
    # Create array to hold labels
    labels = numpy.zeros(mask.shape, numpy.int32)
    # Create list to keep track of label associations
    equivalence = [0]
    # Check upper-left corner
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    # For each non-zero column in row 0, check back pixel label
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If back pixel has a label, use same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    # For each non-zero row
    for x in range(1, mask.shape[0]):
        # Check left-most column, up  and up-right pixels
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If up pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If up-right pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        # For each non-zero column except last in nonzero rows, check up, up-right, up-right, up-left, left pixels
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If up pixel has label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If not up but up-right pixel has label, need to update equivalence table
                    if mask[x - 1, y - 1]:
                        # If up-left pixel has label, relabel up-right equivalence, up-left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If left pixel has label, relabel up-right equivalence, left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither up-left or left pixels are labeled, use up-right equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If not up, or up-right pixels have labels but up-left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If not up, up-right, or up-left pixels have labels but left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, add new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        # Check last pixel in row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # if up pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # if not up but up-left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # if not up or up-left but left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = numpy.array(equivalence)
    # Go backwards through all labels
    for i in range(1, len(equivalence))[::-1]:
        # Convert labels to the lowest value in the set associated with a single object
        labels[numpy.where(labels == i)] = equivalence[i]
    # Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels


label_list = []
for i in range(len(mask[:])):
    label = find_labels(mask[i])
    label_list.append(label)


def filter_by_size(labels, minsize, maxsize):
    # Find label sizes
    sizes = numpy.bincount(labels.ravel())
    # Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        # If the number of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find all pixels for label
            where = numpy.where(labels == i)
            labels[where] = 0
    # Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels


label_array = numpy.array(label_list)
for i in range(len(label_array[:])):
    label_array[i] = filter_by_size(label_array[i],100,100000)


for i in range(len(label_array[:])):
    sizes = numpy.bincount(label_array[i].ravel())
    mean = numpy.mean(sizes[1:])
    sd = numpy.std(sizes[1:])
    label_array[i] = filter_by_size(label_array[i],mean-sd,mean+sd)

print("Gene_ID" + "\t" + "nascentRNA" + "\t" + "PCNA" + "\t" + "log2ratio")

for i in range(len(images)):
    image = images[i]
    DAPI = image[:,:,0]
    nascentRNA = image[:,:,1]
    PCNA = image[:,:,2]
    nuc_counts = numpy.amax(label_array)
    nuc_counts += 1
    for j in range(1, nuc_counts):
        where = numpy.where(label_array[i] == j)
        nascentRNA_signal = numpy.mean(nascentRNA[where])
        PCNA_signal = numpy.mean(PCNA[where])
        log2ratio = numpy.log2(nascentRNA_signal / PCNA_signal)
        if i in [0, 1]: 
            Gene = "APEX1"
        if i in [2, 3]:
            Gene = "PIM2"
        if i in [4,5]:
            Gene = "POLR2B"
        if i in [6,7]: 
            Gene = "SRSF1"
        print(Gene, nascentRNA_signal, PCNA_signal, log2ratio, sep = "\t")

# run the following code in terminal to generate a tsv file with the information gathered from the images
# ./week10_homework.py > week10_signals.tsv




