import os, sys, getopt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import Bio
from Bio import SeqIO

def inputFasta(records_file_location):
    return list(SeqIO.parse(records_file_location, "fasta"))

def align(grna, target):
    max = 0
    #For every starting position of target
    for k in range(len(target)):
      count = 0
      #For every nucleotide of the gRNA sequence 
      for l in range(min(len(target[k:]),len(grna))):
        if target[k+l] == grna[l]:
          count += 1
      if count >= max:
        max = count
        position = k
    return position

def hammingDistance(records_gRNA,records_targets,plot=True):
    #For each target/gRNA pair, shuffle starting point of alignment and compute the hamming distance, take the minimum hamming distance
    hamming_distances = []

    #For every target
    for i in range(len(records_targets)):
        hamming_distance_to_certain_grna = []
        #For every gRNA
        for j in range(len(records_gRNA)):
            max = 0
            #For every starting position of target
            for k in range(len(records_targets[i].seq)):
                count = 0
                #For every nucleotide of the gRNA sequence 
                for l in range(min(len(records_targets[i].seq[k:]),len(records_gRNA[j].seq))):
                    if records_targets[i].seq[k+l] == records_gRNA[j].seq[l]:
                        count += 1
                if count >= max:
                    max = count
            hamming_distance_to_certain_grna.append(20-max)

        hamming_distances.append(hamming_distance_to_certain_grna)
    hamming_distances = np.asarray(hamming_distances)

    if plot:
        plt.figure(figsize=(16, 16))
        plt.imshow(hamming_distances, cmap='bwr', interpolation='nearest')
        cbar = plt.colorbar(shrink = 0.75)
        cbar.set_label('Hamming distance', rotation=270, fontsize = 15, labelpad = 20)

        plt.xlabel('gRNA', fontsize = 15);plt.ylabel('target', fontsize = 15)
        plt.show()

        hamming_distances_1d = hamming_distances.ravel()
        plt.hist(hamming_distances_1d, bins = 15);plt.ylabel('Count');plt.xlabel('Hamming distance')
        plt.show()
    
    return hamming_distances

def one_hot_encode(seq):
    mapping = dict(zip("acgt", range(4)))    
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

def one_hot_encode_duplets(seq):
    mapping = dict(zip(('aa', 'ac', 'ag', 'at', 'ca', 'cc', 'cg', 'ct', 'ga', 'gc', 'gg', 'gt', 'ta', 'tc', 'tg', 'tt'), range(16)))
    seq2 = []
    for i in range(len(seq)-1):
      seq2.append(mapping[seq[i:i+2]])
    return np.eye(16)[seq2]

def calcGC(seq):
  gc_content = [] 
  count = 0
  for l in range(len(seq)):
    if seq[l] == 'g' or seq[l] == 'c':
      count += 1
  return count/len(seq)

def countMotifs(seq, motifs):
  counts = []
  for motif in motifs:
    count = 0
    for i in range(len(seq) - len(motif)):
      if seq[i:i+len(motif)] == motif:
        count += 1
    counts.append(count/(len(seq) - len(motif)))
  return counts

def createMotifs():
    bases = ['a','c','g','t']
    motifs_to_look_for = []
    for i in range(len(bases)):
        motifs_to_look_for.append(bases[i])
        for j in range(len(bases)):
            motifs_to_look_for.append(bases[i]+bases[j])
            for k in range(len(bases)):
                motifs_to_look_for.append(bases[i]+bases[j]+bases[k])
                # for l in range(len(bases)):
                #   motifs_to_look_for.append(bases[i]+bases[j]+bases[k]+bases[l])
    motifs_to_look_for.append('aaaa')
    motifs_to_look_for.append('cccc')
    motifs_to_look_for.append('gggg')
    motifs_to_look_for.append('tttt')