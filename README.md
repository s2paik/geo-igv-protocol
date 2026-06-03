# GEO-to-IGV Protocol: A Beginner-Friendly Workflow

A step-by-step guide for wet lab scientists to download public genomics data from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), visualize it in [IGV](https://igv.org/), and generate publication-quality figures — with no prior computational experience required.

---

## Dataset

This repository uses data from:

> **Zhu Y, Liu J, Wang J, et al.** "Integrative transcriptomic and epigenomic analysis identifies BCL6B as a novel regulator of human pluripotent stem cell to endothelial differentiation." *Protein & Cell* (2025). DOI: [10.1093/procel/pwaf039](https://doi.org/10.1093/procel/pwaf039) | PMID: [40318187](https://pubmed.ncbi.nlm.nih.gov/40318187/)
>
> **GEO Accession:** [GSE186755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186755)

The dataset contains ATAC-seq, CUT&Tag (H3K4me3, H3K27me3, H3K27ac), and SMART-seq2 RNA-seq across five stages of directed hPSC-to-endothelial cell differentiation: **hESC → VMC → EPC-1 → EPC-2 → EC**.

---

## Repository Structure
geo-igv-protocol/
│
├── README.md                           ← You are here
│
├── 01-reference-genome/                ← Download and prepare reference genome
│   └── (coming soon)
│       - Download Ensembl GRCh38 FASTA
│       - Filter GTF to MANE_Select transcripts only
│       - Fix chr prefix for IGV compatibility
│       - Sort and index with IGVTools
│
├── 02-geo-data-download/               ← Find and download GEO data
│   └── (coming soon)
│       - Navigate a GEO series page
│       - Identify file types (bigwig, FPKM, BAM)
│       - Download supplementary files
│
├── 03-expression-analysis/             ← Analyze RNA-seq expression data
│   └── plot_fpkm_annotated.py          ← Annotated Python script
│       - Reads FPKM table from GEO
│       - Generates heatmap + line plots
│       - Fully commented for beginners
│
├── 04-igv-visualization/               ← Load and interpret tracks in IGV
│   └── (coming soon)
│       - Load genome and annotation tracks
│       - Navigate to genes of interest
│       - Interpret ATAC-seq and histone marks
│
└── resources/                          ← Reference notes
└── (coming soon)
- Ensembl vs UCSC naming conventions
- FPKM normalization caveats
- Useful links
---

## Prerequisites

### Software
- **IGV** (Integrative Genomics Viewer) — [download](https://igv.org/doc/desktop/)
- **Python 3.x** — [download](https://www.python.org/downloads/)
- **WSL** (Windows Subsystem for Linux) — Windows users only

### Python packages
```bash
pip install pandas numpy matplotlib seaborn scipy
```

Or, using the requirements file:
```bash
pip install -r requirements.txt
```

---

## Quick Start

### 1. Clone this repository
```bash
git clone https://github.com/s2paik/geo-igv-protocol.git
cd geo-igv-protocol
```

### 2. Download the FPKM data
Go to [GSE186755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186755) and download `GSE186755_hPSC-EC_fpkm.txt.gz`. Decompress it.

### 3. Run the expression analysis
```bash
# Edit the FILE path in the script to point to your FPKM file
python 03-expression-analysis/plot_fpkm_annotated.py
```

This generates `figure_fpkm.pdf` and `figure_fpkm.png` with five panels:
- **A** — Heatmap (log2 fold-change vs hESC)
- **B** — Stemness markers (POU5F1, SOX2, NANOG)
- **C** — EC surface markers (PECAM1, CDH5, KDR, VWF)
- **D** — EC lineage TFs (ETV2, ERG, FLI1)
- **E** — IFN-γ pathway (RELA, IRF1, STAT1, CIITA)

---

## Key Lessons Documented

This repo was built while learning these workflows from scratch. Notes on common pitfalls are included throughout:

- **Chromosome naming mismatch**: Ensembl uses `1, 2, X` while UCSC/IGV uses `chr1, chr2, chrX`. Mixing them causes silent failures in IGV.
- **GTF isoform clutter**: Full Ensembl GTF shows dozens of isoforms per gene. Filtering to `MANE_Select` transcripts gives one clean canonical isoform per gene.
- **FPKM normalization**: FPKM normalizes within a sample but not across samples. Cross-stage comparisons should be interpreted as relative trends, not absolute quantification.
- **Z-score artifacts**: Z-scoring near-zero genes (like CIITA at ~0.09 FPKM) amplifies noise into apparent signal. Log2 fold-change vs a reference stage is more honest.
- **IGV: Load Genome vs Load Track**: The reference FASTA goes under Genomes → Load Genome. Annotation files (GTF, bigwig) go under File → Load from File. Mixing these up causes errors.

---

## License

This repository is for educational purposes. The original data is from Zhu et al. (2025) and is publicly available through NCBI GEO.
