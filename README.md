# Predict-HGT

Predict-HGT is a bioinformatics pipeline that investigates potential horizontal
gene transfer (HGT) between bacterial genomes using sequence similarity,
machine learning and graph-based network analysis.

This was one of my first bioinformatics projects. When I reviewed it in October
2026 I found that some early design choices do not measure what I intended, so
this README now says plainly what the code does and where it falls short (see
[Known problems](#known-problems)).

## What it does

1. Downloads five RefSeq reference genomes from NCBI.
2. Splits each genome file into its FASTA records. For these assemblies that
   means one record per replicon: the chromosome plus any plasmids. This gives
   12 sequences in total (Klebsiella has six plasmids, Salmonella has one).
3. Compares every pair of sequences by k-mer content and builds a similarity
   network with NetworkX.
4. Scores cross-species edges as possible HGT candidates, and runs community
   and hub analysis on the network.

| Species | RefSeq assembly |
|---|---|
| *Escherichia coli* | GCF_000005845.2 |
| *Klebsiella pneumoniae* | GCF_000240185.1 |
| *Salmonella enterica* | GCF_000006945.2 |
| *Pseudomonas aeruginosa* | GCF_000006765.1 |
| *Acinetobacter baumannii* | GCF_000746645.1 |

## Known problems

These are the reasons the results should not be read as HGT evidence.

- **The "genes" are whole replicons.** `extract_cds.py` writes out each FASTA
  record of at least 300 bp. It does not find coding sequences, so each node in
  the network is a chromosome or plasmid, not a gene. A real version would use
  the NCBI CDS files (`*_cds_from_genomic.fna.gz`) or call genes with Prodigal.
- **4-mer similarity saturates.** A sequence of more than a few kilobases
  contains nearly all 256 possible 4-mers, so the Jaccard similarity between
  any two of these sequences is close to 1. Across all 66 pairs the lowest
  value is 0.958, and 36 pairs score exactly 1.0. The 8-mer cosine score in
  `hgt_final_model.py` is less affected, but still compares whole replicons.
- **The labels are circular.** `build_hgt_dataset.py` labels a pair as HGT if
  it is cross-species with similarity between 0.1 and 0.6, and the classifier
  is then trained on that same similarity. Because of the saturation above, no
  pair falls in that range, so `hgt_dataset_v2.csv` has 0 positive labels out
  of 66 and the classifier has nothing to learn.
- **Small sample.** Five genomes and 12 sequences are too few for any
  statistical claim.

## Running it

```bash
git clone https://github.com/iraajgangavaram/Predict-HGT.git
cd Predict-HGT
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Run the scripts from the repository root, in this order:

```bash
python src/download_data.py       # fetch the five genomes into data/raw/genomes
python src/analyze_genome.py      # optional: length and GC summary for E. coli
python src/extract_cds.py         # split genomes into replicon FASTA files
python src/build_hgt_dataset.py   # pairwise features and heuristic labels
python src/train_hgt_v2.py        # random forest on those features
python src/build_hgt_graph.py     # similarity network, figure and .graphml
python src/analyze_hgt_graph.py   # communities and hub nodes
python src/hgt_final_model.py     # 8-mer cosine network and ranked candidates
```

## Scripts

| Script | Role |
|---|---|
| `download_data.py` | Downloads the genomes from NCBI |
| `analyze_genome.py` | Prints sequence count, length and GC for the E. coli genome |
| `extract_cds.py` | Writes each FASTA record of at least 300 bp to `data/genes_cds/` |
| `build_hgt_dataset.py` | Pairwise 4-mer Jaccard similarity, GC and length difference, heuristic label |
| `train_hgt_v2.py` | Random forest on the dataset above |
| `build_hgt_graph.py` | Builds the 4-mer similarity network, saves the figure and `hgt_graph.graphml` |
| `analyze_hgt_graph.py` | Community detection and hub nodes on the saved network |
| `hgt_final_model.py` | 8-mer cosine network with a cross-species and rarity score |

Earlier experiments, kept for reference but not part of the pipeline above:

| Script | What it was |
|---|---|
| `feature_engineering.py` | 3-mer frequencies for the E. coli genome |
| `build_dataset.py`, `train_model.py` | Genome-level 3-mer features and a classifier |
| `extract_genes.py` | Same as `extract_cds.py` without the length filter |
| `build_gene_similarity.py` | First version of the pairwise similarity step |
| `train_hgt_model.py` | First classifier; reads `hgt_dataset.csv`, which no current script writes |

## Outputs

`data/processed/` holds the generated tables, `hgt_graph.graphml` and the
network figure (`hgt_graph_publication.png` and `.pdf`). The saved network has
12 nodes and 46 edges.

## What I would change

- Use annotated CDS or Prodigal gene calls so the nodes are genes.
- Compare genes with alignment (BLAST or DIAMOND) or with k large enough not to
  saturate, rather than 4-mers.
- Label HGT from independent evidence, such as known mobile genetic elements or
  phylogenetic incongruence, instead of from the similarity score itself.
- Test the scoring on a set of documented transfer events before applying it
  to new genomes.

## Licence

MIT. See [LICENSE](LICENSE).
