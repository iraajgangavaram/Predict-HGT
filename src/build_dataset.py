from cProfile import label
from pyexpat import features

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
    try:
        handle = gzip.open(file_path, "rt")
        records = list(SeqIO.parse(handle, "fasta"))

        print(f"{label}: sequences found =", len(records))

        if len(records) == 0:
            print(f"[SKIP] {label} - no sequences")
            return None

        seq = str(records[0].seq).upper()

        features = compute_kmer_freq(seq, K)
        features["gc_content"] = compute_gc(seq)
        features["label"] = label

        return features

    except Exception as e:
        print(f"[ERROR] {label}: {e}")
        return None


def main():

    dataset = []

    print("\n=== Building Multi-Genome Dataset ===\n")

    for file in GENOMES_DIR.glob("*.fna.gz"):
        print("FOUND FILE:", file)

        label = file.stem  # genome name

        print(f"Processing {label}")

        print("PROCESSING:", label)
        features = process_genome(file, label)
        if features is None:
            print("FAILED:", label)
        else:
            dataset.append(features)
        print("ADDED:", label)

    df = pd.DataFrame(dataset)

    Path("data/processed").mkdir(parents=True, exist_ok=True)

    df.to_csv(OUTPUT_FILE, index=False)

    print("\nDONE")
    print(f"Dataset saved to: {OUTPUT_FILE}")
    print(f"Shape: {df.shape}")


if __name__ == "__main__":
    main()