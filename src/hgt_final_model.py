from pathlib import Path
from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx
import math

GENE_DIR = Path("data/genes_cds")


# ----------------------------
# KMER MODEL
# ----------------------------

def kmer_counts(seq, k=8):
    return Counter(seq[i:i+k] for i in range(len(seq)-k+1))


def cosine_similarity(c1, c2):
    intersection = set(c1) & set(c2)

    dot = sum(c1[k] * c2[k] for k in intersection)

    norm1 = math.sqrt(sum(v*v for v in c1.values()))
    norm2 = math.sqrt(sum(v*v for v in c2.values()))

    if norm1 == 0 or norm2 == 0:
        return 0

    return dot / (norm1 * norm2)


# ----------------------------
# LOAD GENES
# ----------------------------

def load_genes():
    genes = {}

    for file in GENE_DIR.glob("*.fasta"):
        record = list(SeqIO.parse(file, "fasta"))[0]
        seq = str(record.seq).upper()

        genes[file.stem] = {
            "seq": seq,
            "kmers": kmer_counts(seq, k=8),
            "species": file.stem.split("_")[0]
        }

    return genes


# ----------------------------
# BUILD GRAPH
# ----------------------------

def build_graph(genes):

    G = nx.Graph()

    for g in genes:
        G.add_node(g, species=genes[g]["species"])

    gene_list = list(genes.keys())

    for i in range(len(gene_list)):
        for j in range(i + 1, len(gene_list)):

            g1 = gene_list[i]
            g2 = gene_list[j]

            sim = cosine_similarity(
                genes[g1]["kmers"],
                genes[g2]["kmers"]
            )

            if sim > 0.2:
                G.add_edge(g1, g2, weight=sim)

    return G


# ----------------------------
# HGT SCORING ENGINE
# ----------------------------

def compute_hgt_scores(G, genes):

    print("\n=== HGT SCORING RESULTS ===\n")

    # gene occurrence frequency (rarity)
    gene_freq = Counter([g.split("_cds_")[0] for g in G.nodes()])

    scored_edges = []

    for u, v, data in G.edges(data=True):

        sim = data["weight"]

        species_u = genes[u]["species"]
        species_v = genes[v]["species"]

        # cross-species bonus
        cross_species = 1 if species_u != species_v else 0.2

        # rarity penalty (rare genes = higher score)
        rarity = 1 / (
            gene_freq[u.split("_cds_")[0]] +
            gene_freq[v.split("_cds_")[0]]
        )

        # FINAL HGT SCORE
        hgt_score = sim * cross_species * rarity

        scored_edges.append((u, v, hgt_score))

    scored_edges.sort(key=lambda x: x[2], reverse=True)

    return scored_edges


# ----------------------------
# MAIN
# ----------------------------

def main():

    print("\n=== PHASE 10: FINAL HGT MODEL ===\n")

    genes = load_genes()

    G = build_graph(genes)

    print("Nodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    scored_edges = compute_hgt_scores(G, genes)

    print("\n=== TOP HGT CANDIDATES (FINAL SCORE) ===\n")

    for e in scored_edges[:20]:
        print(e)


if __name__ == "__main__":
    main()