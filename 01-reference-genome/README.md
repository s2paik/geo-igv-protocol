# 01 — Reference Genome Setup for IGV

This guide walks through downloading a reference genome from Ensembl, filtering
the gene annotation to show only one canonical isoform per gene, and loading
everything into IGV.

---

## Overview

| Step | What | Tool | Time |
|------|------|------|------|
| 1 | Download GRCh38 FASTA from Ensembl | Browser or curl | ~10 min |
| 2 | Download GTF annotation from Ensembl | Browser | ~1 min |
| 3 | Filter GTF to MANE_Select transcripts | grep (WSL) | ~2 min |
| 4 | Add chr prefix to GTF (for UCSC compatibility) | sed (WSL) | ~2 min |
| 5 | Sort and index GTF | IGVTools | ~5 min |
| 6 | Load into IGV | IGV GUI | ~1 min |

---

## Why not just use IGV's built-in hg38?

You can — and for many purposes it is the easiest option. IGV's built-in hg38
(Genomes > Load Genome from Server > Human hg38) works perfectly and uses
chr-prefixed naming compatible with most public datasets.

The reason to add a custom Ensembl GTF on top is the gene annotation. IGV's
default RefSeq track shows many isoforms per gene, making the display cluttered.
Filtering Ensembl's GTF to MANE_Select transcripts gives you one clean,
clinically agreed-upon canonical isoform per gene.

**Recommended setup:**
- Use IGV's built-in hg38 as the reference genome (no download needed)
- Load a filtered Ensembl GTF as a custom annotation track on top

---

## Step 1: Download the Ensembl GTF

Go to the Ensembl FTP site in your browser:
https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
Download this file (click it directly — it is about 50 MB):
Homo_sapiens.GRCh38.115.gtf.gz

Avoid files with "abinitio", "chr", or "chr_patch_hapl_scaff" in the name.

Decompress the file after downloading (right-click > Extract, or use 7-Zip).

---

## Step 2: Filter to MANE_Select transcripts

MANE_Select (Matched Annotation from NCBI and EMBL-EBI) gives you one
canonical transcript per gene — the gold standard agreed upon by both
NCBI and Ensembl.

### Why filter?

The full Ensembl GTF contains every annotated isoform, including low-confidence
and computationally predicted ones. For CIITA alone, there are 20+ isoforms.
After filtering to MANE_Select, you get exactly one.

### How to filter (in WSL terminal):

```bash
# Navigate to where your GTF file is
cd /mnt/c/Users/seung/Downloads

# Filter: keep only lines tagged as MANE_Select
# The input file may be inside a folder with the same name (common with browser downloads)
grep 'tag "MANE_Select"' Homo_sapiens.GRCh38.115.gtf > Homo_sapiens.GRCh38.115.MANE.gtf
```

Check the output:
```bash
ls -lh Homo_sapiens.GRCh38.115.MANE.gtf
# Should be ~250-300 MB (much smaller than the original ~1.5 GB)
# If it is 0 bytes, the input file path is wrong
```

### How grep works here:
```
grep 'tag "MANE_Select"' input.gtf > output.gtf
|         |                |            |
|         |                |            +-- write results to this new file
|         |                +-- the file to search through
|         +-- the text pattern to find (lines containing this exact string)
+-- the command: search for lines matching a pattern
```
Only lines containing the text `tag "MANE_Select"` are kept. Everything else
(non-MANE isoforms, predicted transcripts) is discarded.

---

## Step 3: Add chr prefix for UCSC/IGV compatibility

Ensembl uses numeric chromosome names (1, 2, X) while IGV's built-in hg38
uses UCSC-style names (chr1, chr2, chrX). If they don't match, tracks
silently fail to display — no error message, just blank tracks.

```bash
# Add 'chr' to the start of each line that begins with a chromosome number or X/Y/M
sed 's/^\([0-9XYM]\)/chr\1/' Homo_sapiens.GRCh38.115.MANE.gtf > Homo_sapiens.GRCh38.115.MANE.ucsc.gtf
```

Verify it worked:
```bash
head -1 Homo_sapiens.GRCh38.115.MANE.ucsc.gtf
# Should start with "chr1" not "1"
```

### How sed works here:
```
sed 's/^([0-9XYM])/chr\1/' input > output
|   |  |       |      |   |
|   |  |       |      |   +-- \1 = whatever was matched inside ( )
|   |  |       |      +-- replacement: "chr" followed by the original character
|   |  |       +-- match: any digit, X, Y, or M
|   |  +-- ^ means "at the start of the line only"
|   +-- s/find/replace/ = substitution command
+-- stream editor: processes text line by line
---
```
## Step 4: Sort and index with IGVTools

IGV needs sorted and indexed files for fast navigation.

1. Open **IGV**
2. Go to **Tools > Run igvtools**
3. Select command: **Sort**
   - Input file: `Homo_sapiens.GRCh38.115.MANE.ucsc.gtf`
   - Click **Run**
   - Output: `Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf`
4. Select command: **Index**
   - Input file: `Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf`
   - Click **Run**
   - Output: `Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf.idx`

**Important:** Always Sort first, then Index. The index only works on sorted files.

---

## Step 5: Load into IGV

1. Load the reference genome: **Genomes > Load Genome from Server > Human (hg38)**
2. Load your filtered annotation: **File > Load from File** > select `Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf`

You should now see one clean canonical isoform per gene.

**Do NOT load the GTF under "Genomes > Load Genome from File"** — that is for
FASTA files only. Loading a GTF as a genome causes errors.

---

## Common Pitfalls

| Problem | Symptom | Fix |
|---------|---------|-----|
| GTF downloaded as a folder | Permission denied errors when trying to read the file | The actual file is inside the folder. Use: `folder_name/file_name.gtf` |
| Chr naming mismatch | Tracks load but display nothing | Check with `head -1 file.gtf` — if it starts with "1" instead of "chr1", run the sed command |
| Forgot to Sort before Index | IGVTools index fails or IGV shows errors | Always Sort first, then Index |
| Loaded GTF as genome instead of track | "chrName is null" error, IGV crashes | Use Genomes menu for FASTA only. Use File menu for GTF. |
| File is 0 bytes after grep | Wrong input path | Check if your file is inside a subfolder with `ls -la` |

---

## Available tag filters

The GTF contains several useful tags beyond MANE_Select. You can swap the
grep filter depending on how many isoforms you want:

| Tag | What it gives you | Isoforms per gene |
|-----|-------------------|-------------------|
| `MANE_Select` | Gold standard canonical transcript | 1 |
| `Ensembl_canonical` | Ensembl's chosen representative | 1 |
| `MANE_Plus_Clinical` | MANE_Select + clinically relevant variants | 1-2 |
| `gencode_basic` | Well-supported transcripts only | 2-5 |
| (no filter) | Everything including predictions | 10-30+ |

To check what tags exist in your GTF:
```bash
grep -o 'tag "[^"]*"' your_file.gtf | sort -u
```

---

## File summary

After completing all steps, you should have these files:

| File | Size | Purpose |
|------|------|---------|
| `Homo_sapiens.GRCh38.115.gtf.gz` | ~50 MB | Original download (can delete) |
| `Homo_sapiens.GRCh38.115.gtf` | ~1.5 GB | Decompressed full GTF (can delete) |
| `Homo_sapiens.GRCh38.115.MANE.gtf` | ~270 MB | Filtered to MANE_Select (can delete) |
| `Homo_sapiens.GRCh38.115.MANE.ucsc.gtf` | ~270 MB | Chr prefix added (can delete) |
| `Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf` | ~270 MB | **Keep — load this in IGV** |
| `Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf.idx` | ~few KB | **Keep — IGV needs this alongside the sorted file** |
