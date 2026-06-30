from Bio import SeqIO
import gzip
from pathlib import Path

GENOMES_DIR = Path("data/raw/genomes")
OUTPUT_DIR = Path("data/genes")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)


def extract_genes(file_path, genome_name):

    handle = gzip.open(file_path, "rt")
    records = list(SeqIO.parse(handle, "fasta"))

    gene_count = 0

    for record in records:

        gene_id = f"{genome_name}_gene_{gene_count}.fasta"

        output_path = OUTPUT_DIR / gene_id

        with open(output_path, "w") as f:
            f.write(f">{gene_id}\n")
            f.write(str(record.seq) + "\n")

        gene_count += 1

    print(f"{genome_name}: extracted {gene_count} genes")


def main():

    print("\n=== Gene Extraction Phase ===\n")

    for file in GENOMES_DIR.glob("*.fna.gz"):

        genome_name = file.stem

        extract_genes(file, genome_name)

    print("\nDONE: gene extraction complete")


if __name__ == "__main__":
    main()