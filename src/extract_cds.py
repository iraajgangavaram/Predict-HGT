from Bio import SeqIO
import gzip
from pathlib import Path

GENOMES_DIR = Path("data/raw/genomes")
OUTPUT_DIR = Path("data/genes_cds")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)


def extract_cds(file_path, genome_name):

    handle = gzip.open(file_path, "rt")

    records = list(SeqIO.parse(handle, "fasta"))

    cds_count = 0

    for record in records:

        # This is still not perfect CDS, but better structured than contigs
        seq = str(record.seq)

        if len(seq) < 300:
            continue  # filter very small fragments

        gene_id = f"{genome_name}_cds_{cds_count}.fasta"

        output_path = OUTPUT_DIR / gene_id

        with open(output_path, "w") as f:
            f.write(f">{gene_id}\n")
            f.write(seq + "\n")

        cds_count += 1

    print(f"{genome_name}: extracted {cds_count} CDS-like sequences")


def main():

    print("\n=== CDS Extraction Phase ===\n")

    for file in GENOMES_DIR.glob("*.fna.gz"):

        genome_name = file.stem

        extract_cds(file, genome_name)

    print("\nDONE")


if __name__ == "__main__":
    main()