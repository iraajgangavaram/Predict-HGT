from pathlib import Path
import networkx as nx
from networkx.algorithms import community
from collections import Counter

GRAPH_PATH = "data/processed/hgt_graph.graphml"


# ----------------------------
# LOAD GRAPH
# ----------------------------

def load_graph():
    print("\n=== Loading HGT Graph ===\n")
    return nx.read_graphml(GRAPH_PATH)


# ----------------------------
# HUB GENES
# ----------------------------

def find_hub_genes(G):
    print("\n=== Hub Genes (Most Connected) ===\n")

    degree_dict = dict(G.degree())

    sorted_genes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)

    for gene, deg in sorted_genes[:10]:
        print(gene, "degree:", deg)


# ----------------------------
# COMMUNITY DETECTION
# ----------------------------

def detect_communities(G):

    print("\n=== HGT Communities ===\n")

    communities = community.greedy_modularity_communities(G)

    for i, c in enumerate(communities):

        print(f"\nCommunity {i+1} (size={len(c)}):")

        for gene in list(c)[:10]:
            print(" ", gene)


# ----------------------------
# SPECIES MIX ANALYSIS
# ----------------------------

def species_distribution(G):

    print("\n=== Species Distribution in Graph ===\n")

    species = [G.nodes[n]["species"] for n in G.nodes()]

    counts = Counter(species)

    for s, c in counts.items():
        print(s, ":", c)


# ----------------------------
# MAIN
# ----------------------------

def main():

    G = load_graph()

    print("\nNodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    species_distribution(G)
    find_hub_genes(G)
    detect_communities(G)


if __name__ == "__main__":
    main()