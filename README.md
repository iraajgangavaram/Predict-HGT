## Overview

Project Overview

Predict-HGT is a bioinformatics pipeline that investigates potential Horizontal Gene Transfer (HGT) between bacterial genomes using sequence similarity, machine learning and graph-based network analysis.

Horizontal Gene Transfer is an important evolutionary process in bacteria that enables the movement of genetic material between organisms. It plays a major role in the spread of antibiotic resistance, virulence factors and metabolic adaptations. Detecting HGT computationally is challenging because sequence similarity alone cannot distinguish between shared ancestry and genuine transfer events.

This project implements an end-to-end computational workflow that begins with publicly available bacterial genomes and progresses through genome parsing, coding sequence (CDS) extraction, feature engineering, machine learning and network analysis. Each extracted gene is represented using k-mer frequency vectors, allowing pairwise sequence similarity to be calculated using cosine similarity. A weighted similarity graph is then constructed to identify highly connected genes, network communities and candidate cross-species relationships that may warrant further biological investigation.

The pipeline demonstrates how multiple computational biology techniques can be integrated into a single reproducible workflow using Python. While the project is intended as an exploratory proof-of-concept rather than a validated HGT detection tool, it illustrates many of the core methods used in modern comparative genomics and network biology.


## Workflow

```text
                 ┌──────────────────────────┐
                 │ Download Bacterial Genomes│
                 └─────────────┬────────────┘
                               │
                               ▼
                 ┌──────────────────────────┐
                 │ Parse FASTA Genome Files │
                 └─────────────┬────────────┘
                               │
                               ▼
                 ┌──────────────────────────┐
                 │ Extract CDS Sequences    │
                 └─────────────┬────────────┘
                               │
                               ▼
                 ┌──────────────────────────┐
                 │ Generate k-mer Features  │
                 └─────────────┬────────────┘
                               │
                               ▼
                 ┌──────────────────────────┐
                 │ Compute Cosine Similarity│
                 └─────────────┬────────────┘
                               │
                               ▼
                 ┌──────────────────────────┐
                 │ Build Gene Similarity    │
                 │ Graph (NetworkX)         │
                 └─────────────┬────────────┘
                               │
                 ┌─────────────┴────────────┐
                 ▼                          ▼
      Community Detection           Hub Gene Analysis
                 │                          │
                 └─────────────┬────────────┘
                               ▼
                 ┌──────────────────────────┐
                 │ HGT Candidate Scoring    │
                 └─────────────┬────────────┘
                               ▼
                 ┌──────────────────────────┐
                 │ Ranked HGT Candidates &  │
                 │ Network Visualisation    │
                 └──────────────────────────┘
```
## Installation

# Clone the repository

```bash
git clone https://github.com/iraajgangavaram/Predict-HGT.git
cd Predict-HGT
```

# Create a virtual environment

**macOS / Linux**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows**

```bash
python -m venv .venv
.venv\Scripts\activate
```

# Install dependencies

```bash
pip install -r requirements.txt
```

The project uses the following major Python libraries:

* BioPython
* pandas
* NumPy
* scikit-learn
* NetworkX
* Matplotlib
* tqdm

# Run the complete pipeline

```bash
python src/download_data.py
python src/analyse_genome.py
python src/process_genomes.py
python src/extract_cds.py
python src/build_hgt_dataset.py
python src/train_hgt_model.py
python src/build_hgt_graph.py
python src/analyze_hgt_graph.py
python src/hgt_final_model.py
```

## Objectives

- Extract biologically relevant features from genomic sequences
- Identify patterns associated with potential horizontal gene transfer
- Apply machine learning models to classify gene transfer signals
- Visualise results using publication-style scientific plots

## Data Sources

- Genomic sequence data from [NCBI / GTDB]
- Annotated gene regions used for feature extraction
- Processed using Python-based bioinformatics pipeline

## Project Structure

```
Predict-HGT/
│
├── data/
│   ├── raw/
│   │   └── genomes/                 # Downloaded bacterial genomes
│   │
│   ├── genes_cds/                   # Extracted CDS FASTA files
│   │
│   └── processed/
│       ├── genome_ml_dataset.csv
│       ├── hgt_dataset.csv
│       ├── hgt_graph.graphml
│       └── other generated outputs
│
├── src/
│   ├── download_data.py             # Download bacterial genomes
│   ├── analyse_genome.py            # Genome statistics
│   ├── process_genomes.py           # Feature extraction
│   ├── extract_cds.py               # CDS extraction
│   ├── build_hgt_dataset.py         # Construct HGT dataset
│   ├── train_model.py               # Baseline ML classifier
│   ├── train_hgt_model.py           # HGT prediction model
│   ├── build_hgt_graph.py           # Construct similarity network
│   ├── analyze_hgt_graph.py         # Community & hub analysis
│   └── hgt_final_model.py           # Final HGT scoring pipeline
│
├── README.md
├── requirements.txt
└── LICENSE
```
The HGT Similarity Network Graph is in the 'data' folder inside 'processed' folder and is named hgt_graph_publication.png.
# Directory Description

* **data/raw/** — Stores downloaded bacterial genome FASTA files.
* **data/genes_cds/** — Contains extracted coding DNA sequences used for downstream analysis.
* **data/processed/** — Stores generated datasets, graph files and analysis outputs.
* **src/** — Contains all Python scripts implementing the computational pipeline.
* **README.md** — Project documentation.
* **requirements.txt** — Python dependencies required to reproduce the project.

## Methodology

The project implements a computational pipeline for identifying candidate Horizontal Gene Transfer (HGT) events from bacterial genome sequences using sequence similarity and graph-based analysis.

Complete bacterial genomes from multiple species (Escherichia coli, Salmonella enterica, Klebsiella pneumoniae, Pseudomonas aeruginosa and Acinetobacter baumannii) were downloaded from public genomic databases. Genome FASTA files were parsed using BioPython, and coding DNA sequences (CDS) were extracted to create a gene-level dataset.

Each extracted sequence was represented using 8-mer frequency vectors. Cosine similarity was calculated between every pair of genes to quantify sequence similarity while reducing the effect of sequence length. These similarities were used to construct a weighted gene similarity network in which nodes represent genes and weighted edges represent pairwise sequence similarity.

To investigate potential HGT candidates, a scoring framework was developed that combines:

Cross-species sequence similarity
Graph connectivity
Gene rarity within the similarity network

The project also includes supervised machine learning experiments using Random Forest classifiers as an exploratory baseline for genome classification before transitioning to an unsupervised graph-based approach that is more appropriate for studying evolutionary relationships.

Finally, NetworkX was used to perform community detection and identify highly connected hub genes within the similarity network, providing an additional level of biological interpretation.

## Results

The complete analysis pipeline successfully processed genomes from five bacterial species and extracted twelve CDS sequences for downstream analysis.

Key outputs generated by the project include:

Automated genome downloading and parsing
CDS extraction pipeline
Genome feature dataset generation
Machine learning baseline classifier
Weighted gene similarity network
Community detection analysis
Hub gene identification
Ranked HGT candidate pairs

The final similarity network contained:

Metric	Value
Bacterial species analysed	5
Genes (nodes)	12
Similarity edges	46
Network communities detected	3

The graph-based HGT scoring system produced a ranked list of candidate cross-species similarities. The highest scoring relationships were observed between genes from Acinetobacter, Escherichia, Salmonella and Pseudomonas, indicating regions of high sequence conservation across species.

Community detection identified three distinct gene clusters, including one mixed-species community that may represent conserved evolutionary relationships worthy of further investigation.

## Biological Interpretation

The similarity network demonstrates that genes from different bacterial species frequently cluster together based on sequence composition. These relationships indicate shared evolutionary history and conserved genomic regions across multiple bacterial taxa.

Importantly, the results should not be interpreted as definitive evidence of horizontal gene transfer. Sequence similarity alone cannot distinguish between:

Conserved orthologous genes
Shared ancestry
Genuine horizontal gene transfer events

Instead, the project should be viewed as an exploratory framework that identifies candidate regions requiring further biological validation.

The graph analysis illustrates how network-based methods can be applied to comparative genomics by identifying highly connected hub genes and mixed-species communities that warrant further investigation. These computational approaches are commonly used as an initial screening step before integrating additional biological evidence such as phylogenetic analysis, mobile genetic element annotation or plasmid databases.

## Limitations

Several limitations should be considered when interpreting the results.

The dataset contains a relatively small number of genomes and extracted CDS sequences, limiting the statistical power of the analysis.
Gene extraction was simplified and does not distinguish between annotated genes, hypothetical proteins or non-coding regions.
Horizontal gene transfer labels were generated using heuristic rules rather than experimentally validated HGT events.
Sequence similarity is only one source of evidence and cannot independently confirm horizontal gene transfer.
No functional annotation, phylogenetic reconstruction, plasmid detection or mobile genetic element analysis was incorporated into the current pipeline.
The machine learning component should be regarded as a proof-of-concept rather than a predictive model suitable for biological inference.

These limitations reflect the scope of an undergraduate research project while highlighting opportunities for future development.

## Future Work

Several extensions could improve both the biological accuracy and computational capability of the pipeline.

Future work may include:

Analysing substantially larger collections of bacterial genomes from diverse taxonomic groups.
Integrating functional gene annotation using databases such as RefSeq or UniProt.
Incorporating known mobile genetic element and plasmid databases to improve HGT candidate prioritisation.
Comparing inferred HGT candidates with published benchmark datasets to evaluate prediction performance.
Applying phylogenetic incongruence methods to distinguish horizontal gene transfer from vertical inheritance.
Investigating Graph Neural Networks (GNNs) for learning directly from gene similarity networks.
Developing an interactive web application for visualising HGT networks and candidate transfer events.

These extensions would transform the current proof-of-concept pipeline into a more comprehensive comparative genomics framework suitable for larger-scale microbial evolutionary studies
