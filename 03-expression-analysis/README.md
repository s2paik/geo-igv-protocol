# 03 — Expression Analysis: FPKM Visualization

This folder contains an annotated Python script that reads FPKM expression
data from GEO and generates publication-quality figures showing gene expression
changes across hPSC-to-endothelial cell differentiation.

---

## What the script produces

A single figure (PDF + PNG) with five panels:

| Panel | Contents | What to look for |
|-------|----------|-----------------|
| A | Heatmap — log2 fold-change vs hESC | Red = upregulated, blue = downregulated relative to undifferentiated state |
| B | Stemness markers (POU5F1, SOX2, NANOG) | Should decline sharply from hESC onward |
| C | EC surface markers (PECAM1, CDH5, KDR, VWF) | Should rise progressively from EPC to EC |
| D | EC lineage TFs (ETV2, ERG, FLI1) | ETV2 peaks transiently at EPC then drops; ERG and FLI1 rise and sustain |
| E | IFN-g pathway (RELA, IRF1, STAT1, CIITA) | RELA rises in EC; CIITA near-zero throughout (requires IFN-g to induce) |

---

## How to run

### Prerequisites

Python 3.x with the following packages:
pip install pandas numpy matplotlib seaborn scipy

pip install pandas numpy matplotlib seaborn scipy

source ~/bioenv/bin/activate
pip install pandas numpy matplotlib seaborn scipy
### Step 1: Download the data

Go to [GSE186755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186755)
and download `GSE186755_hPSC-EC_fpkm.txt.gz`. Decompress it.

### Step 2: Edit the file path

Open `plot_fpkm_annotated.py` and change these two lines near the top
to match where your file is saved:

```python
FILE = "/path/to/your/GSE186755_hPSC-EC_fpkm.txt"
OUT  = "/path/to/your/output/folder/"
```

For example on Windows via WSL:

```python
FILE = "/mnt/c/Users/YOUR_NAME/Downloads/GSE186755_hPSC-EC_fpkm.txt"
OUT  = "/mnt/c/Users/YOUR_NAME/Downloads/"
```

### Step 3: Run
python3 plot_fpkm_annotated.py
Output files appear in your specified output folder:
- `figure_fpkm.pdf` (vector — for publication, editable in Illustrator)
- `figure_fpkm.png` (raster — for quick preview)

---

## Heatmap normalization explained

The heatmap (Panel A) uses log2 fold-change relative to the mean hESC value
for each gene:
log2((FPKM_sample + 1) / (mean_FPKM_hESC + 1))
This means:
- hESC columns are near zero (white) since they are the reference
- Red = upregulated relative to hESC
- Blue = downregulated relative to hESC
- Each gene is normalized independently (row-wise)
- The pseudocount (+1) prevents division by zero for silent genes

### Why not z-score?

Z-scoring normalizes each gene to its own mean and standard deviation across
all samples. This works well for genes with real dynamic range, but creates
artifacts for near-zero genes. For example, CIITA has a maximum FPKM of 0.09
across all stages — biologically meaningless variation. Z-scoring amplifies
this noise into apparent signal, making it look like CIITA peaks at VMC.
The fold-change approach correctly shows CIITA as flat (near-white) because
all values divided by a near-zero reference remain near zero.

---

## FPKM normalization caveat

FPKM normalizes for sequencing depth and gene length within each library,
making genes comparable within a single sample. However, FPKM does not
properly normalize across samples because it assumes total RNA output per
cell is constant between conditions. During dramatic biological transitions
like hPSC-to-EC differentiation, this assumption may not hold.

Cross-stage comparisons in this figure should be interpreted as relative
expression trends, not absolute quantification. For rigorous differential
expression analysis with statistical testing, raw count data analyzed with
DESeq2 or edgeR would be more appropriate.

---

## How to customize

### Add or remove genes

Edit the gene lists in Section 3 of the script:

```python
GENES_STEM = ["POU5F1", "SOX2", "NANOG"]
GENES_ECMARKER = ["PECAM1", "CDH5", "KDR", "VWF"]
GENES_ECTF = ["ETV2", "ERG", "FLI1"]
GENES_IFN = ["RELA", "IRF1", "STAT1", "CIITA"]
```

Then add matching entries in the COLORS, LSTYLE, and MARKERS dictionaries
in Section 5. The script comments explain each option.

### Change colors

Edit hex codes in the COLORS dictionary. Use https://colorbrewer2.org
for colorblind-safe palettes.

### Change the reference stage

Replace "hESC" in the normalization section with any other stage:

```python
# Current: normalize against hESC
hesc_cols = [c for c in df.columns if c.startswith("hESC")]

# Alternative: normalize against EC
hesc_cols = [c for c in df.columns if c.startswith("EC")]
```

### Use a different dataset

The script works with any tab-separated expression table where:
- First column contains gene symbols
- Remaining columns contain expression values (FPKM, TPM, etc.)

You will need to update the STAGE_MAP dictionary to match your column names.

---

## Files in this folder

| File | Description |
|------|-------------|
| `plot_fpkm_annotated.py` | Main script — fully annotated with beginner-friendly comments |
| `README.md` | This file |

---

## Dataset reference

Zhu Y, Liu J, Wang J, et al. "Integrative transcriptomic and epigenomic
analysis identifies BCL6B as a novel regulator of human pluripotent stem
cell to endothelial differentiation." Protein & Cell (2025).
GEO: [GSE186755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186755)
