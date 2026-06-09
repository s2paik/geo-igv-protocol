# GEO-to-IGV Protocol: A Beginner-Friendly Workflow

A step-by-step guide for wet lab scientists to download public genomics data from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), visualize it in [IGV](https://igv.org/), and generate publication-quality figures — with no prior computational experience required.

---

## Dataset

This repository uses data from:

> **Zhu Y, Liu J, Wang J, et al.** "Integrative transcriptomic and epigenomic analysis identifies BCL6B as a novel regulator of human pluripotent stem cell to endothelial differentiation." *Protein & Cell* (2025). DOI: [10.1093/procel/pwaf039](https://doi.org/10.1093/procel/pwaf039) | PMID: [40318187](https://pubmed.ncbi.nlm.nih.gov/40318187/)
>
> **GEO Accession:** [GSE186755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186755)

The dataset contains ATAC-seq, CUT&Tag (H3K4me3, H3K27me3, H3K27ac), and SMART-seq2 RNA-seq across five stages of directed hPSC-to-endothelial cell differentiation: **hESC -> VMC -> EPC-1 -> EPC-2 -> EC**.

---

## Repository Structure

| Folder | Contents | Status |
|--------|----------|--------|
| 01-reference-genome/ | Download Ensembl GRCh38 FASTA, filter GTF to MANE_Select, fix chr prefix, sort and index for IGV | Available |
| 02-geo-data-download/ | Navigate GEO, identify file types (bigwig, FPKM, BAM), download supplementary files | Available |
| 03-expression-analysis/ | Annotated Python script for FPKM heatmap and line plots | Available |
| 04-igv-visualization/ | Load genome and tracks, navigate to genes, interpret ATAC-seq and histone marks, CIITA case study | Available |
| 05-gene-explorer/ | Interactive gene search tool with alias support (CLI + Google Colab) | Available |
| resources/ | Ensembl vs UCSC naming, FPKM normalization caveats, useful links | Available |

---

## Prerequisites

### Software
- **IGV** (Integrative Genomics Viewer) - [download](https://igv.org/doc/desktop/)
- **Python 3.x** - [download](https://www.python.org/downloads/)
- **WSL** (Windows Subsystem for Linux) - Windows users only

### Python packages
