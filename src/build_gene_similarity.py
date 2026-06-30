from pathlib import Path
import numpy as np
from itertools import combinations
from Bio import SeqIO

GENE_DIR = Path("data/genes_cds")


# ----------------------------
# SIMPLE KMER FUNCTION
# ----------------------------

def get_kmers(seq, k=4):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]


def kmer_profile(seq, k=4):
    kmers = get_kmers(seq, k)
    return set(kmers)


# ----------------------------
# JACCARD SIMILARITY
# ----------------------------

def similarity(set1, set2):
    if len(set1) == 0 or len(set2) == 0:
        return 0
    return len(set1 & set2) / len(set1 | set2)


# ----------------------------
# LOAD GENES
# ----------------------------

def load_genes():
    genes = {}

    for file in GENE_DIR.glob("*.fasta"):

        record = list(SeqIO.parse(file, "fasta"))[0]
        seq = str(record.seq).upper()

        genes[file.name] = kmer_profile(seq)

    return genes


# ----------------------------
# BUILD GRAPH EDGES
# ----------------------------

def main():

    print("\n=== Building Gene Similarity Network ===\n")

    genes = load_genes()

    edges = []

    for g1, g2 in combinations(genes.keys(), 2):

        sim = similarity(genes[g1], genes[g2])

        if sim > 0.3:   # threshold for potential HGT signal
            edges.append((g1, g2, sim))

    print(f"Total genes: {len(genes)}")
    print(f"Potential HGT edges: {len(edges)}")

    for e in edges[:10]:
        print(e)


if __name__ == "__main__":
    main()