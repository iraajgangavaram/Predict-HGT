from pathlib import Path
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

GENE_DIR = Path("data/genes_cds")
OUTPUT_FIG = Path("data/processed/hgt_graph_publication.png")


def get_kmers(seq, k=4):
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))


def similarity(a, b):
    if len(a) == 0 or len(b) == 0:
        return 0
    return len(a & b) / len(a | b)


def load_genes():
    genes = {}
    for file in GENE_DIR.glob("*.fasta"):
        record = list(SeqIO.parse(file, "fasta"))[0]
        seq = str(record.seq).upper()
        genes[file.stem] = {
            "seq": seq,
            "kmers": get_kmers(seq),
            "species": file.stem.split("_")[0]
        }
    return genes


def build_graph(genes):
    G = nx.Graph()
    for name, meta in genes.items():
        G.add_node(name, species=meta["species"])

    gene_list = list(genes.keys())
    for i in range(len(gene_list)):
        for j in range(i + 1, len(gene_list)):
            g1 = gene_list[i]
            g2 = gene_list[j]
            sim = similarity(genes[g1]["kmers"], genes[g2]["kmers"])
            if sim > 0.2:
                G.add_edge(g1, g2, weight=sim)
    return G


def plot_graph(G):
    OUTPUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    print("\nGenerating publication-quality figure...")
    plt.figure(figsize=(14, 10), dpi=300)
    pos = nx.spring_layout(G, seed=42, k=0.8)

    species_list = sorted(list(set(nx.get_node_attributes(G, "species").values())))
    color_map = {s: i for i, s in enumerate(species_list)}

    node_colors = [
        color_map[G.nodes[n]["species"]]
        for n in G.nodes()
    ]

    edge_widths = [
        G[u][v]["weight"] * 4
        for u, v in G.edges()
    ]

    nx.draw_networkx_edges(
        G,
        pos,
        width=edge_widths,
        alpha=0.35,
        edge_color="gray"
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        cmap=plt.cm.Set2,
        node_size=1100,
        edgecolors="black",
        linewidths=0.8
    )

    labels = {}
    for node in G.nodes():
        species = G.nodes[node]["species"].capitalize()
        cds = node.split("_cds_")[-1]
        labels[node] = f"{species}\nCDS {cds}"

    nx.draw_networkx_labels(
        G,
        pos,
        labels=labels,
        font_size=8,
        font_weight="bold",
        font_color="black"
    )

    legend_elements = []
    for species in species_list:
        colour = plt.cm.Set2(
            color_map[species] /
            max(1, len(species_list) - 1)
        )
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker='o',
                color='w',
                label=species.capitalize(),
                markerfacecolor=colour,
                markeredgecolor='black',
                markersize=10
            )
        )

    plt.legend(
        handles=legend_elements,
        title="Bacterial Species",
        loc="upper left",
        fontsize=9,
        title_fontsize=10,
        frameon=True
    )

    plt.title(
        "Figure 1. Gene Similarity Network of Candidate Horizontal Gene Transfer Events",
        fontsize=16,
        fontweight="bold",
        pad=20
    )

    plt.figtext(
        0.5,
        0.01,
        "Nodes represent coding DNA sequences (CDS) extracted from five bacterial genomes. "
        "Edges indicate cosine similarity between 8-mer frequency profiles. "
        "Node colours denote bacterial species and edge thickness is proportional to sequence similarity.",
        ha="center",
        fontsize=9,
        wrap=True
    )

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
    plt.savefig(str(OUTPUT_FIG.with_suffix(".pdf")), bbox_inches="tight")
    plt.close()


def main():
    genes = load_genes()
    G = build_graph(genes)
    print("Nodes (genes):", G.number_of_nodes())
    print("Edges:", G.number_of_edges())
    plot_graph(G)


if __name__ == "__main__":
    main()
