from pathlib import Path
from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx
import math
import matplotlib.pyplot as plt

GENE_DIR = Path("data/genes_cds")


# ----------------------------
# KMER FEATURES
# ----------------------------

def kmer_counts(seq, k=8):
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    return Counter(kmers)


def cosine_similarity(c1, c2):
    intersection = set(c1.keys()) & set(c2.keys())

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
# BUILD GRAPH + HGT SCORE
# ----------------------------

def main():

    print("\n=== Phase 8: HGT Graph + Scoring + Visualisation ===\n")

    genes = load_genes()

    G = nx.Graph()

    # species tracking for ubiquity penalty
    gene_species_map = defaultdict(set)

    for g in genes:
        G.add_node(g, species=genes[g]["species"])
        gene_species_map[g.split("_cds_")[0]].add(genes[g]["species"])

    edges_added = 0

    gene_list = list(genes.keys())

    for i in range(len(gene_list)):
        for j in range(i + 1, len(gene_list)):

            g1 = gene_list[i]
            g2 = gene_list[j]

            sim = cosine_similarity(
                genes[g1]["kmers"],
                genes[g2]["kmers"]
            )

            # remove noise
            if sim > 0.2:

                G.add_edge(g1, g2, weight=sim)
                edges_added += 1

    print("Nodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    # ----------------------------
    # HGT SCORING
    # ----------------------------

    hgt_edges = []

    for u, v, data in G.edges(data=True):

        if genes[u]["species"] != genes[v]["species"]:

            ubiquity_penalty = 1 / (
                len(gene_species_map[u.split("_cds_")[0]]) +
                len(gene_species_map[v.split("_cds_")[0]])
            )

            hgt_score = data["weight"] * ubiquity_penalty

            hgt_edges.append((u, v, hgt_score))

    hgt_edges.sort(key=lambda x: x[2], reverse=True)

    print("\n=== Top HGT Candidates ===\n")

    for e in hgt_edges[:15]:
        print(e)

    # ----------------------------
    # VISUALISATION
    # ----------------------------

    print("\nGenerating graph plot...")

    plt.figure(figsize=(10, 7))

    pos = nx.spring_layout(G, seed=42)

    # node colors by species
    species_list = list(set(nx.get_node_attributes(G, "species").values()))
    color_map = {s: i for i, s in enumerate(species_list)}

    node_colors = [
        color_map[G.nodes[n]["species"]] for n in G.nodes()
    ]

    nx.draw(
        G,
        pos,
        with_labels=True,
        node_color=node_colors,
        cmap=plt.cm.Set2,
        node_size=800,
        font_size=8
    )

    plt.title("HGT Gene Similarity Network")
    plt.show()

    # ----------------------------
    # EXPORT FOR BIOINFO TOOLS
    # ----------------------------

    nx.write_graphml(G, "data/processed/hgt_graph.graphml")

    print("\nSaved graph to: data/processed/hgt_graph.graphml")


if __name__ == "__main__":
    main()