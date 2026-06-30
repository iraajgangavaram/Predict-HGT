# Ensure Biopython is installed externally (do not run pip from within the script)
from Bio import SeqIO
import gzip
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
file_path = project_root / "data" / "raw" / "genomes" / "ecoli.fna.gz"

def compute_gc(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)

def main():
    with gzip.open(file_path, "rt") as handle:
        sequences = list(SeqIO.parse(handle, "fasta"))

    print("\n=== Genome Summary ===\n")
    print(f"Number of sequences: {len(sequences)}")

    total_length = 0
    total_gc = 0

    for record in sequences:
        seq = str(record.seq)
        total_length += len(seq)
        total_gc += len(seq) * compute_gc(seq)

    print(f"Total genome length: {total_length}")
    print(f"Average GC content: {total_gc / total_length:.4f}")

if __name__ == "__main__":
    main()