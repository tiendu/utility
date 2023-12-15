import argparse
import statistics
import copy
import re
import random
import numpy as np
from itertools import groupby
import networkx as nx
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='k-nucleotide clustering.')
parser.add_argument('-i', type = str, required = True)
parser.add_argument('-k', type = int, required = True)
parser.add_argument('-p', type = int, required = True)
args = parser.parse_args()

def fasta_iter(filename):
    fin = open(filename, 'rb')
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerstr = str(header.__next__(), 'utf-8')
        long_name = headerstr.strip().replace('>', '')
        name = long_name.split()[0]
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield name, seq

class IdSequence:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence
    def base_freq(self):
        bases = {i:0 for i in ["A", "T", "G", "C"]}
        for i in bases.keys():
            bases[i] = self.sequence.count(i) / len(self.sequence)
        return bases
    def kmer_freq(self, k):
        kmers = {i:0 for i in generate_kmer(k)}
        subseq = [self.sequence[i:i + k:1] for i in range(0, len(self.sequence) - k + 1, 1)]
        for i in kmers.keys():
            kmers[i] = sum(map(lambda k: k == i, subseq)) / len(subseq)
        return kmers
    def norm_kmer_freq(self, k):
        base_freq = self.base_freq()
        kmer_freq = self.kmer_freq(k)
        norm_freq = {}
        for i in kmer_freq.keys():
            norm_freq[i] = kmer_freq[i]
            nu = [*i]
            for j in nu:
                if base_freq[j] != 0:
                    norm_freq[i] = norm_freq[i] / base_freq[j]
        return norm_freq
    
def generate_kmer(k):
    bases_1 = ["A", "T", "G", "C"]
    bases_2 = bases_1.copy()
    i = 0
    while i < k - 1:
        i += 1
        temp = []
        for m in bases_1:
            for n in bases_2:
                temp.append(m + n)
        bases_2 = None
        bases_2 = temp
    return bases_2

def remove_tandem_repeat(sequence, k, c):
    for i in generate_kmer(k):
        sequence = re.sub("(" + i + ")" + "{" + str(c) + ",}", "", sequence)
        return sequence
    
def create_reverse_complement(sequence):
    rev = ""
    for i in range(len(sequence) - 1, -1, -1):
        nu = [*sequence]
        if nu[i] == "A":
            rev += "T"
        elif nu[i] == "T":
            rev += "A"
        elif nu[i] == "G":
            rev += "C"
        elif nu[i] == "C":
            rev += "G"
    return rev

labels = []
features = []
for i in fasta_iter(args.i):
    n = IdSequence(i[0], remove_tandem_repeat(i[1], args.k, 3))
    labels.append(n.id)
    if args.p == 1:
        palindromes = []
        for j in n.norm_kmer_freq(args.k).keys():
            if j == create_reverse_complement(j):
                palindromes.append(n.norm_kmer_freq(args.k)[j])
        features.append(palindromes)
    else:
        features.append(list(n.norm_kmer_freq(args.k).values()))
features = np.array(features)
labels = np.array(labels)
labels = labels.reshape(-1, 1)

G = nx.empty_graph()

for i in range(features.shape[0]):
    for j in range(i + 1, features.shape[0]):
        dist = np.linalg.norm(features[i] - features[j])
        G.add_edge(labels[i][0], labels[j][0], weight = dist)

threshold = statistics.median([i[2] for i in G.edges.data("weight")])
filtered = list(filter(lambda e: e[2] > threshold, (e for e in G.edges.data("weight"))))

sorted_filtered = sorted(filtered, key = lambda x: x[2], reverse = True)
for i in sorted_filtered:
    if G.degree(i[0]) > 1 and G.degree(i[1]) > 1:
        G.remove_edge(i[0], i[1])

# credit to https://github.com/53RT/Highly-Connected-Subgraphs-Clustering-HCS for making this possible
def highly_connected(G, E):
    """
    Checks if the graph G is highly connected
    Highly connected means, that splitting the graph G into subgraphs needs more than 0.5*|V| edge deletions
    :param G: Graph G
    :param E: Edges needed for splitting G
    :return: True if G is highly connected, otherwise False
    """
    return len(E) > len(G.nodes) / 2

def remove_edges(G, E):
    """
    Removes all edges E from G
    Iterates over all edges in E and removes them from G
    :param G: Graph to remove edges from
    :param E: One or multiple Edges
    :return: Graph with edges removed
    """
    for edge in E:
        G.remove_edge(*edge)
    return G

def HCS(G):
    """
    Basic HCS Algorithm
    cluster labels, removed edges are stored in global variables
    :param G: Input graph
    :return: Either the input Graph if it is highly connected, otherwise a Graph composed of
    Subgraphs that build clusters
    """
    E = nx.algorithms.connectivity.cuts.minimum_edge_cut(G)
    if not highly_connected(G, E):
        G = remove_edges(G, E)
        sub_graphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        if len(sub_graphs) == 2:
            H = HCS(sub_graphs[0])
            _H = HCS(sub_graphs[1])
            G = nx.compose(H, _H)
    return G

def improved_HCS(G):
    """
    Implements improvements mentioned in the paper
    1. Iterated HCS
    2. Singleton adoption
    3. Removing Low Degree Vertices
    """
    pass

def labelled_HCS(G):
    """
    Runs basic HCS and returns Cluster Labels
    :param G: Input graph
    :return: List of cluster assignments for the single vertices
    """
    _G = HCS(G)
    sub_graphs = (G.subgraph(c).copy() for c in nx.connected_components(_G))
    labels = [list(g) for g in sub_graphs]
    return labels

#G = HCS(G)

#edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
#nx.draw(G, edgelist = edges, edge_color = weights, width = 1, edge_cmap = plt.cm.Reds, with_labels = True)
#plt.show()

for i in labelled_HCS(G):
    print(len(i), i)
