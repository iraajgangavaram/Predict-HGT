from pathlib import Path
from itertools import combinations
import pandas as pd
from Bio import SeqIO

GENE_DIR = Path("data/genes_cds")


# ----------------------------
# FEATURES
# ----------------------------

def get_kmers(seq, k=4):
    return set(seq[i:i+k] for i in range(len(seq)-k+1))


def similarity(a, b):
    if len(a) == 0 or len(b) == 0:
        return 0
    return len(a & b) / len(a | b)


def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq) if len(seq) > 0 else 0


def length(seq):
    return len(seq)


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
            "kmers": get_kmers(seq),
            "gc": gc_content(seq),
            "len": length(seq),
        }

    return genes


# ----------------------------
# BUILD DATASET
# ----------------------------

def main():

    print("\n=== Phase 6: Building Improved HGT Dataset ===\n")

    genes = load_genes()

    rows = []

    for g1, g2 in combinations(genes.keys(), 2):

        a = genes[g1]
        b = genes[g2]

        sim = similarity(a["kmers"], b["kmers"])
        gc_diff = abs(a["gc"] - b["gc"])
        len_diff = abs(a["len"] - b["len"])

        # improved heuristic label
        # (proxy: cross-species + moderate similarity = potential HGT)
        if g1.split("_")[0] != g2.split("_")[0] and 0.1 < sim < 0.6:
            label = 1
        else:
            label = 0

        rows.append([
            g1,
            g2,
            sim,
            gc_diff,
            len_diff,
            label
        ])

    df = pd.DataFrame(rows, columns=[
        "gene1",
        "gene2",
        "similarity",
        "gc_diff",
        "length_diff",
        "label"
    ])

    output = Path("data/processed/hgt_dataset_v2.csv")
    output.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output, index=False)

    print("Dataset shape:", df.shape)
    print("Saved to:", output)


if __name__ == "__main__":
    main()