from Bio import SeqIO
import gzip
import pandas as pd
from itertools import product
from collections import Counter
from pathlib import Path

# ----------------------------
# CONFIG
# ----------------------------

project_root = Path(__file__).resolve().parents[1]
file_path = project_root / "data" / "raw" / "genomes" / "ecoli.fna.gz"
K = 3  # k-mer size (start small)

# ----------------------------
# FUNCTIONS
# ----------------------------

def get_kmers(seq, k):
    """Extract k-mers from a DNA sequence."""
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def all_possible_kmers(k):
    """Generate all possible DNA k-mers."""
    return [''.join(p) for p in product("ATGC", repeat=k)]


def compute_kmer_freq(seq, k):
    """Compute normalized k-mer frequency vector."""
    kmers = get_kmers(seq, k)
    counts = Counter(kmers)

    total = sum(counts.values())
    all_kmers = all_possible_kmers(k)

    freq = {kmer: counts.get(kmer, 0) / total for kmer in all_kmers}

    return freq


def main():

    records = list(SeqIO.parse(gzip.open(file_path, "rt"), "fasta"))

    genome_seq = str(records[0].seq).upper()

    print("\n=== Feature Engineering ===\n")
    print(f"Genome length: {len(genome_seq)}")

    kmer_features = compute_kmer_freq(genome_seq, K)

    df = pd.DataFrame([kmer_features])

    output_path = project_root / "data" / "processed" / "ecoli_kmer_features.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    print(f"K-mer features saved to: {output_path}")
    print(f"Number of features: {len(df.columns)}")


if __name__ == "__main__":
    main()