# 05 - Gene Expression Explorer

Interactive tools for searching and visualizing any gene in the GSE186755
FPKM dataset. Search by official gene symbol or common alias (CD31, OCT4,
VEGFR2, etc.). Genes are automatically categorized into biological groups.

---

## Two versions available

| Version | File | Best for |
|---------|------|----------|
| Python CLI | gene_explorer.py | Command-line users, scripting, batch processing |
| Google Colab | gene_explorer_colab.ipynb | Team sharing, no installation needed |

---

## Python CLI version

### Usage

    python3 gene_explorer.py --file GSE186755_hPSC-EC_fpkm.txt --genes CIITA ERG FLI1 RELA

    python3 gene_explorer.py --file GSE186755_hPSC-EC_fpkm.txt --genes CD31 OCT4 VEGFR2

    python3 gene_explorer.py --file GSE186755_hPSC-EC_fpkm.txt --group stemness ec_tfs

    python3 gene_explorer.py --file GSE186755_hPSC-EC_fpkm.txt --group ec_markers --genes CIITA

### Available preset groups

| Group name | Genes included |
|------------|---------------|
| stemness | POU5F1, SOX2, NANOG |
| ec_markers | PECAM1, CDH5, KDR, VWF |
| ec_tfs | ETV2, ERG, FLI1 |
| ifn_pathway | RELA, IRF1, STAT1, CIITA |
| hematopoietic | RUNX1, SPI1, CEBPB, GATA1 |
| smooth_muscle | ACTA2, TAGLN, CNN1, MYH11 |
| mesoderm | T, MESP1, MIXL1, EOMES |

### Output

Generates gene_explorer_output.pdf and gene_explorer_output.png with:
- Panel A: Heatmap (log2 fold-change vs hESC)
- Panels B+: Line plots grouped by biological category

Also prints a table of raw FPKM values to the terminal.

---

## Google Colab version

### How to use

1. Open gene_explorer_colab.ipynb in Google Colab
2. Run all cells (Runtime > Run all)
3. Upload your FPKM file when prompted
4. Type gene names in the search box and click Add
5. Use preset group buttons to add whole groups
6. Click Plot to generate the figure

### Sharing with your team

Upload the notebook to Google Drive and share the link. Anyone with access
can open it in Colab, upload the FPKM file, and search genes — no
installation required.

---

## Supported gene aliases

The tool recognizes common alternative names:

| You type | Finds |
|----------|-------|
| CD31 | PECAM1 |
| CD144, VE-cadherin | CDH5 |
| OCT4 | POU5F1 |
| VEGFR2, CD309 | KDR |
| TIE2, CD202B | TEK |
| CD105, endoglin | ENG |
| alpha-SMA, SMA | ACTA2 |
| PU.1 | SPI1 |
| ER71 | ETV2 |
| P65, NF-KB | RELA |
| MHC-II, MHC2TA | CIITA |
| HLA-DR | HLA-DRA |
| brachyury | T |

---

## Dataset reference

Zhu Y, Liu J, Wang J, et al. Integrative transcriptomic and epigenomic
analysis identifies BCL6B as a novel regulator of human pluripotent stem
cell to endothelial differentiation. Protein and Cell (2025).
GEO: GSE186755
