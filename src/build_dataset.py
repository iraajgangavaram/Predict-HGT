from Bio import SeqIO
import gzip
import pandas as pd
from itertools import product
from collections import Counter
from pathlib import Path

# ----------------------------
# CONFIG
# ----------------------------

GENOMES_DIR = Path(__file__).resolve().parents[1] / "data" / "raw" / "genomes"
OUTPUT_FILE = "data/processed/genome_ml_dataset.csv"
K = 3

# ----------------------------
# KMER FUNCTIONS
# ----------------------------

def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def all_kmers(k):
    return [''.join(p) for p in product("ATGC", repeat=k)]


def compute_kmer_freq(seq, k):
    kmers = get_kmers(seq, k)
    counts = Counter(kmers)

    total = sum(counts.values())
    all_possible = all_kmers(k)

    return {kmer: counts.get(kmer, 0) / total for kmer in all_possible}


def compute_gc(seq):
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


# ----------------------------
# MAIN PIPELINE
# ----------------------------

def process_genome(file_path, label):
    records = list(SeqIO.parse(gzip.open(file_path, "rt"), "fasta"))
    seq = str(records[0].seq).upper()

    features = compute_kmer_freq(seq, K)
    features["gc_content"] = compute_gc(seq)
    features["label"] = label

    return features


def main():

    dataset = []

    print("\n=== Building Multi-Genome Dataset ===\n")

    for file in GENOMES_DIR.glob("*.fna.gz"):

        label = file.stem  # genome name

        print(f"Processing {label}")

        features = process_genome(file, label)
        dataset.append(features)

    df = pd.DataFrame(dataset)

    Path("data/processed").mkdir(parents=True, exist_ok=True)

    df.to_csv(OUTPUT_FILE, index=False)

    print("\nDONE")
    print(f"Dataset saved to: {OUTPUT_FILE}")
    print(f"Shape: {df.shape}")


if __name__ == "__main__":
    main()