from pathlib import Path
from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx
import math
import matplotlib.pyplot as plt

# ----------------------------
# CONFIG
# ----------------------------

GENE_DIR = Path("data/genes_cds")
OUTPUT_GRAPH = "data/processed/hgt_graph.graphml"
OUTPUT_FIG = "data/processed/hgt_graph_publication.png"

# ----------------------------
# HELPERS
# ----------------------------

def get_species(gene_id):
    return gene_id.split("_")[0]


# ----------------------------
# KMER FEATURES
# ----------------------------

def kmer_counts(seq, k=8):
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    c = Counter(kmers)
    total = sum(c.values())

    # normalize (important for cosine stability)
    return {k: v / total for k, v in c.items()}


def cosine_similarity(c1, c2):
    intersection = set(c1.keys()) & set(c2.keys())

    dot = sum(c1[k] * c2[k] for k in intersection)

    norm1 = math.sqrt(sum(v * v for v in c1.values()))
    norm2 = math.sqrt(sum(v * v for v in c2.values()))

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

        gene_id = file.stem

        genes[gene_id] = {
            "seq": seq,
            "kmers": kmer_counts(seq, k=8),
            "species": get_species(gene_id)
        }

    return genes


# ----------------------------
# MAIN PIPELINE
# ----------------------------

def main():

    print("\n=== Phase: HGT Graph + Publication Analysis ===\n")

    genes = load_genes()

    G = nx.Graph()

    # track gene family presence across species
    gene_species_map = defaultdict(set)

    # add nodes
    for g, data in genes.items():
        G.add_node(g, species=data["species"])
        gene_species_map[g.split("_cds_")[0]].add(data["species"])

    gene_list = list(genes.keys())

    # ----------------------------
    # BUILD GRAPH
    # ----------------------------

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

    print("Nodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    # ----------------------------
    # HGT SCORING
    # ----------------------------

    hgt_edges = []

    for u, v, data in G.edges(data=True):

        if genes[u]["species"] != genes[v]["species"]:

            u_family = u.split("_cds_")[0]
            v_family = v.split("_cds_")[0]

            ubiquity_penalty = 1 / (
                len(gene_species_map[u_family]) +
                len(gene_species_map[v_family])
            )

            hgt_score = data["weight"] * ubiquity_penalty

            hgt_edges.append((u, v, hgt_score))

    hgt_edges.sort(key=lambda x: x[2], reverse=True)

    print("\n=== Top HGT Candidates ===\n")
    for e in hgt_edges[:15]:
        print(e)

    # ----------------------------
    # VISUALISATION (PUBLICATION STYLE)
    # ----------------------------

    print("\nGenerating publication-quality figure...")

    plt.figure(figsize=(14, 10), dpi=300)

    pos = nx.spring_layout(G, seed=42, k=0.8)

    species_list = list(set(nx.get_node_attributes(G, "species").values()))
    color_map = {s: i for i, s in enumerate(species_list)}

    node_colors = [
        color_map[G.nodes[n]["species"]] for n in G.nodes()
    ]

    # edges
    nx.draw_networkx_edges(
        G,
        pos,
        alpha=0.25,
        width=0.8
    )

    # nodes
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        cmap=plt.cm.Set2,
        node_size=900,
        edgecolors="black",
        linewidths=0.5
    )

    # labels (FIXED)
    labels = {n: n for n in G.nodes()}

    nx.draw_networkx_labels(
        G,
        pos,
        labels=labels,
        font_size=7,
        font_color="black",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.7,
            boxstyle="round,pad=0.2"
        )
    )

    plt.title("HGT Gene Similarity Network", fontsize=14)
    plt.axis("off")
    plt.tight_layout()

    plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
    plt.show()

    # ----------------------------
    # EXPORT GRAPH
    # ----------------------------

    nx.write_graphml(G, OUTPUT_GRAPH)
    print(f"\nSaved graph to: {OUTPUT_GRAPH}")
    print(f"Saved figure to: {OUTPUT_FIG}")


if __name__ == "__main__":
    main()