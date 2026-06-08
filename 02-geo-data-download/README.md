# 02 - Downloading Data from NCBI GEO

This guide walks through finding, understanding, and downloading public genomics data from the NCBI Gene Expression Omnibus (GEO).

---

## What is GEO?

GEO (Gene Expression Omnibus) is the main public repository where researchers deposit their sequencing data when they publish a paper. Every dataset gets a unique accession number starting with GSE (for a series/project) or GSM (for individual samples).

---

## Step 1: Navigate to a GEO series

Go to https://www.ncbi.nlm.nih.gov/geo/ and enter the accession number in the search box (e.g. GSE186755).

The series page shows the title, summary, linked publication, platform information, list of all samples (GSM numbers), and supplementary files available for download.

---

## Step 2: Understand the file types

GEO datasets typically contain several types of files:

| File type | Extension | What it contains | When you need it |
|-----------|-----------|-----------------|-----------------|
| FPKM/TPM table | .txt, .csv | Gene-level expression values | Quick gene lookups, plotting |
| Count matrix | .txt, .csv | Raw read counts (not normalized) | DESeq2/edgeR analysis |
| bigWig | .bw | Genome-wide signal track | Viewing in IGV |
| BAM | .bam | Individual aligned reads | Sashimi plots, read-level analysis |
| BED/bedGraph | .bed | Genomic intervals, peak calls | Viewing peaks in IGV |
| FASTQ | .fastq | Raw unaligned reads | Re-running alignment from scratch |

### Which files to download depends on your goal:

| Goal | Files needed |
|------|-------------|
| Look up gene expression across conditions | FPKM/TPM table only |
| Visualize chromatin/histone signals in IGV | bigWig files |
| Run proper differential expression analysis | Raw count matrix |
| Make sashimi plots of splice junctions | BAM files (from SRA) |
| Re-analyze everything from scratch | FASTQ files (from SRA) |

---

## Step 3: Find the supplementary files

On the GEO series page, scroll down to the section labeled Supplementary file.

For GSE186755, the supplementary files are:

| File | Size | Contents |
|------|------|----------|
| GSE186755_RAW.tar | 2.0 GB | All bigWig tracks (ATAC-seq, CUT and Tag) |
| GSE186755_hPSC-EC_fpkm.txt.gz | 582 KB | FPKM table for differentiation series |
| GSE186755_BK-EC_fpkm.txt.gz | 714 KB | FPKM table for BCL6B knockout experiment |

Click the (http) link next to each file to download directly in your browser.

---

## Step 4: Download and extract bigWig files

Click (http) next to GSE186755_RAW.tar on the series page. After downloading, extract in WSL with: tar -xf GSE186755_RAW.tar

Or download individual bigWig files by going to individual sample pages (GSM numbers).

---

## Step 5: Download expression data

Click (http) next to the .txt.gz file. Decompress after downloading. In WSL: gunzip GSE186755_hPSC-EC_fpkm.txt.gz

The resulting .txt file is tab-separated and can be opened in Excel, Python, R, or any text editor.

---

## Step 6: Check the genome build

The series page often does NOT state which genome assembly was used. To find it, click any individual sample (GSM number) and look for a line labeled Genome_build (e.g. hg38) or the Data processing section.

For GSE186755: data was aligned to hg38 (UCSC) using bowtie2. The bigWig files use chr-prefixed names (chr1, chr2, chrX) and are compatible with IGV built-in hg38 genome.

Loading hg38-aligned bigWigs against an Ensembl reference (which uses 1, 2, X) will show blank tracks with no error message. Always check the genome build first.

---

## Understanding sample naming conventions

GEO sample names encode the experiment. For GSE186755, a file like GSM5660478_ATAC_hESC_rep1.bw tells you the GEO accession (GSM5660478), the assay type (ATAC), the cell stage (hESC), and the replicate number (rep1).

The four assay types in this dataset:

| Assay | What it measures | Suggested IGV color |
|-------|-----------------|-------------------|
| ATAC-seq | Chromatin accessibility | Blue |
| H3K4me3 (CUT and Tag) | Active promoters | Red |
| H3K27me3 (CUT and Tag) | Repressed/silenced regions | Green |
| H3K27ac (CUT and Tag) | Active enhancers and promoters | Orange |

---

## Accessing raw data from SRA

If you need raw FASTQ or BAM files not provided as supplementary files, go to the SRA (Sequence Read Archive) link at the bottom of the GEO page. Install the SRA toolkit (sudo apt install sra-toolkit) and download with fastq-dump. This is only needed for advanced analyses like re-alignment or sashimi plots.

---

## Common pitfalls

| Problem | Symptom | Fix |
|---------|---------|-----|
| File downloads as a folder | Cannot open the file | The actual file is inside the folder |
| .gz file wont open in Excel | Garbled text | Decompress first with gunzip or 7-Zip |
| BigWig tracks blank in IGV | No signal, no error | Genome build mismatch (chr1 vs 1) |
| Cant match sample names | Confusing GSM numbers | Check the series page metadata table |
| No processed data available | Only raw data in SRA | You may need to process FASTQ yourself |
| Downloaded wrong FPKM file | BK-EC not hPSC-EC | BK-EC is knockout; hPSC-EC is differentiation |

---

## Quick reference: useful GEO URLs

| Resource | URL |
|----------|-----|
| GEO home | https://www.ncbi.nlm.nih.gov/geo/ |
| Search by accession | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=YOUR_ACCESSION |
| GEO DataSets (keyword search) | https://www.ncbi.nlm.nih.gov/gds/ |
| SRA (raw sequencing data) | https://www.ncbi.nlm.nih.gov/sra |
| GEO2R (online differential expression) | https://www.ncbi.nlm.nih.gov/geo/geo2r/ |

---

## Dataset reference

Zhu Y, Liu J, Wang J, et al. Integrative transcriptomic and epigenomic analysis identifies BCL6B as a novel regulator of human pluripotent stem cell to endothelial differentiation. Protein and Cell (2025). GEO: GSE186755
