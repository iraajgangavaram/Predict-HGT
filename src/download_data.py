import os
import requests
import tqdm as tqdm_ # type: ignore
from pathlib import Path

# ----------------------------
# CONFIGURATION
# ----------------------------

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"

CARD_DIR = RAW_DIR / "card"
GENOMES_DIR = RAW_DIR / "genomes"

NCBI_GENOME_URLS = {
    "ecoli": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
    "klebsiella": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz",
    "salmonella": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz",
    "pseudomonas": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz",
    "acinetobacter": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/746/645/GCF_000746645.1_ASM74664v1/GCF_000746645.1_ASM74664v1_genomic.fna.gz"
}

# ----------------------------
# UTILITIES
# ----------------------------

def create_dirs():
    CARD_DIR.mkdir(parents=True, exist_ok=True)
    GENOMES_DIR.mkdir(parents=True, exist_ok=True)


def download_file(url, output_path):
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))

    with open(output_path, "wb") as file, tqdm_.tqdm(
        desc=output_path.name,
        total=total_size,
        unit='B',
        unit_scale=True
    ) as bar:

        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                file.write(chunk)
                bar.update(len(chunk))


def download_genomes():
    for name, url in NCBI_GENOME_URLS.items():
        output_file = GENOMES_DIR / f"{name}.fna.gz"

        if output_file.exists():
            print(f"[SKIP] {name} already exists")
            continue

        print(f"[DOWNLOAD] {name}")
        download_file(url, output_file)


def main():
    print("\n=== Predict-HGT Data Pipeline ===\n")

    create_dirs()
    download_genomes()

    print("\n[DONE] Data download complete.")
    print(f"Data stored in: {DATA_DIR}\n")


if __name__ == "__main__":
    main()