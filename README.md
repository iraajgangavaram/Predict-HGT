# Predict-HGT

Predict-HGT is a bioinformatics pipeline that investigates potential horizontal
gene transfer (HGT) between bacterial genomes using sequence similarity,
machine learning and graph-based network analysis.

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

## Limitations

- Sequences are compared at the replicon level (chromosomes and plasmids).
  Gene-level analysis would use the NCBI CDS files
  (`*_cds_from_genomic.fna.gz`) or gene calls from Prodigal.
- 4-mer profiles of long sequences are very similar to one another (pairwise
  Jaccard similarity of 0.958 or higher across all 66 pairs), so longer k-mers,
  such as the 8-mer cosine score in `hgt_final_model.py`, separate sequences
  better.
- HGT labels in `build_hgt_dataset.py` are heuristic and based on similarity,
  rather than on experimentally validated transfer events.
- The dataset is small: five genomes and 12 sequences.

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
| `build_gene_similarity.py` | Earlier version of the pairwise similarity step |
| `train_hgt_model.py` | Earlier classifier trained on `hgt_dataset.csv` |

## Outputs

`data/processed/` holds the generated tables, `hgt_graph.graphml` and the
network figure (`hgt_graph_publication.png` and `.pdf`). The saved network has
12 nodes and 46 edges.

## Future work

- Use annotated CDS or Prodigal gene calls so the nodes are genes.
- Compare genes with alignment (BLAST or DIAMOND) or with k large enough not to
  saturate, rather than 4-mers.
- Label HGT from independent evidence, such as known mobile genetic elements or
  phylogenetic incongruence, instead of from the similarity score itself.
- Test the scoring on a set of documented transfer events before applying it
  to new genomes.

## Author

Iraaj Gangavaram

## Licence

MIT. See [LICENSE](LICENSE).
