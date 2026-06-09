# 04 - IGV Visualization Guide

Step-by-step guide for loading and interpreting genomic tracks in IGV
(Integrative Genomics Viewer), specifically for ATAC-seq and histone
modification CUT&Tag data from the GSE186755 hPSC-to-EC differentiation dataset.

---

## What is IGV?

IGV is a desktop genome browser that lets you visually inspect sequencing
data at any genomic locus. You can see where chromatin is open (ATAC-seq),
where histone marks are deposited (CUT&Tag/ChIP-seq), and how gene
structures are organized - all by simply typing a gene name and scrolling.

Download: https://igv.org/doc/desktop/

---

## Step 1: Load the reference genome

The reference genome tells IGV what the DNA sequence looks like and where
genes are located. There are two ways to load it:

### Option A: Use IGV built-in hg38 (recommended)

1. Open IGV
2. Go to Genomes > Load Genome from Server
3. Select Human (hg38)

This is the easiest option and uses chr-prefixed naming (chr1, chr2, chrX)
which is compatible with the GSE186755 bigWig files.

### Option B: Load a custom Ensembl FASTA

Only needed if you want a custom reference. See Section 01 for details.

### CRITICAL: Genome vs Track

| Menu | What to load | File types |
|------|-------------|------------|
| Genomes > Load Genome | Reference genome only | .fa, .fasta |
| File > Load from File | Everything else | .gtf, .bw, .bam, .bed |

Loading a GTF under the Genomes menu causes crashes. Loading a FASTA
under the File menu causes errors. This is the single most common
beginner mistake in IGV.

---

## Step 2: Load a custom gene annotation (optional)

IGV comes with a built-in RefSeq gene track, but it shows many isoforms
per gene. For a cleaner view with one canonical isoform per gene, load
a MANE_Select-filtered Ensembl GTF (see Section 01 for how to create this).

1. Go to File > Load from File
2. Select your sorted GTF file (e.g. Homo_sapiens.GRCh38.115.MANE.ucsc.sorted.gtf)
3. The .idx index file must be in the same folder (IGV finds it automatically)

You will now have two gene tracks: the built-in RefSeq and your custom
Ensembl MANE track. You can right-click either track to change display:
- Expanded: shows all features with labels
- Collapsed: single line per gene (cleanest)
- Squished: compressed view, shows features but smaller

---

## Step 3: Load bigWig tracks

BigWig files contain the genome-wide signal from ATAC-seq and CUT&Tag.
Each file represents one sample + one assay type.

1. Go to File > Load from File
2. Select one or more .bw files
3. They appear as signal tracks in the main panel

### Organizing your tracks

For the GSE186755 dataset, organize tracks by assay type (not by stage).
This makes it easy to compare the same mark across differentiation stages:

    --- ATAC-seq ---
    ATAC_hESC_rep1.bw
    ATAC_VMC_rep1.bw
    ATAC_EPC-1_rep1.bw
    ATAC_EC_rep1.bw

    --- H3K4me3 ---
    H3K4me3_hESC_rep1.bw
    H3K4me3_VMC_rep1.bw
    H3K4me3_EPC-1_rep1.bw
    H3K4me3_EC_rep1.bw

    --- H3K27me3 ---
    (same pattern)

    --- H3K27ac ---
    (same pattern)

### Color coding (recommended convention)

Right-click each track > Change Track Color to maintain visual consistency:

| Assay | Suggested color | What it shows |
|-------|----------------|---------------|
| ATAC-seq | Blue | Open chromatin / accessible regions |
| H3K4me3 | Red | Active promoters |
| H3K27me3 | Green | Repressed / silenced regions |
| H3K27ac | Orange | Active enhancers and promoters |

### Adjusting track height and scale

- Right-click track > Set Data Range: match all tracks of the same assay
  to the same scale for fair comparison across stages
- Drag the bottom edge of a track to resize its height
- Right-click > Autoscale: lets IGV pick the range (quick but not
  comparable across tracks)

---

## Step 4: Navigate to a gene

### By gene name
Type a gene name in the search box at the top (e.g. CIITA, ERG, PECAM1)
and press Enter or click Go. IGV jumps to that locus.

### By coordinates
Type exact coordinates for precise navigation:
    chr16:10,866,000-10,920,000

This is useful when you want to zoom into a specific promoter region
rather than seeing the whole gene.

### Zoom controls
- Use the +/- buttons in the top right to zoom in and out
- Click and drag on the chromosome ideogram (top bar) to jump to a region
- Double-click the main track area to zoom in at that position

---

## Step 5: Interpreting the tracks

### What each mark tells you at a gene promoter

| Pattern at promoter | H3K4me3 | H3K27me3 | H3K27ac | ATAC | Interpretation |
|--------------------|---------|----------|---------|------|----------------|
| Active gene | High | Low | High | High | Actively transcribed |
| Poised/bivalent | Present | Present | Low | Variable | Ready to activate or silence |
| Silenced | Low | High | Low | Low | Repressed by Polycomb |
| Accessible but inactive | Low | Low | Low | High | Open but not transcribed |

### What to look for across differentiation stages

When comparing hESC vs VMC vs EPC vs EC at a gene of interest:

- ATAC peak appearing = chromatin opening at that stage
- ATAC peak disappearing = chromatin closing
- H3K4me3 appearing at TSS = promoter becoming active
- H3K27me3 disappearing = loss of repression
- H3K27ac appearing = enhancer or promoter activation
- H3K27me3 to H3K27ac switch = classic activation transition

### Important caveats

- Broad ATAC accessibility (as seen in hESCs) does not necessarily mean
  a gene is active. hPSCs have globally open chromatin at many loci that
  are not transcribed. This is sometimes called "phantom accessibility."

- The presence of H3K4me3 at a TSS without gene expression can indicate
  a poised or bivalent state (H3K4me3 + H3K27me3 together).

- Always cross-reference epigenetic tracks with expression data (FPKM).
  A gene showing open chromatin but near-zero FPKM is accessible but
  not transcribed - meaning something besides chromatin state is limiting
  expression (e.g. absence of required transcription factors).

---

## Step 6: Save and export views

### Save a session
File > Save Session saves all your loaded tracks and their settings.
You can reopen the exact same view later with File > Open Session.

### Export an image
File > Save Image (or the camera icon) saves a screenshot of the
current view as PNG or SVG.

For publication figures:
- Zoom to the exact region you want to show
- Set all tracks to the same data range for fair comparison
- Use SVG export for vector graphics (editable in Illustrator)

---

## Example: Viewing CIITA across differentiation

A practical walkthrough of what we analyzed in this project:

### 1. Navigate to CIITA
Type CIITA in the search box. The gene spans chr16:10,866,206-10,943,021.

### 2. What you see in hESC vs EC

At the CIITA promoter region:

- ATAC-seq: hESC shows a broad accessibility peak at the TSS region.
  EC shows comparable or slightly lower accessibility. This is the
  opposite of what you might expect if CIITA were simply "closed" in
  hESCs and "opened" in ECs.

- H3K4me3: May appear in later stages (VMC/EPC/EC) but not in hESC,
  suggesting the promoter is becoming transcriptionally primed during
  differentiation even though CIITA is not expressed under basal conditions.

- H3K27me3: Variable across stages. Check whether it resolves (decreases)
  at the CIITA locus during differentiation.

- H3K27ac: Generally low across all stages at CIITA under basal conditions,
  consistent with CIITA being IFN-gamma inducible rather than constitutively
  expressed.

### 3. Key biological insight

The ATAC-seq data shows that chromatin accessibility alone does not explain
why hESCs cannot activate CIITA in response to IFN-gamma while ECs can.
The answer lies in the co-activator environment: the presence of EC-specific
transcription factors (ERG, FLI1) and NF-kB (RELA) that are needed alongside
STAT1 to drive CIITA expression from its pIV promoter.

See Section 03 (Expression Analysis) for the FPKM data supporting this model.

---

## CIITA promoter structure

CIITA has four distinct promoters, each used by different cell types:

| Promoter | Cell type | Activation | Key TFs |
|----------|-----------|------------|---------|
| pI | Dendritic cells, macrophages | Constitutive | PU.1, IRF8 |
| pII | Minor/rare | Constitutive (weak) | Poorly defined |
| pIII | B cells, activated T cells | Constitutive | PAX5, B-lineage TFs |
| pIV | Non-hematopoietic (ECs, fibroblasts) | IFN-gamma inducible | STAT1, IRF1, RELA |

The four promoters span approximately 12-14 kb at the 5' end of the gene.
Each drives a unique first exon that splices to shared downstream exons,
producing different CIITA protein isoforms.

When looking at CIITA in IGV, be aware that the MANE_Select transcript
corresponds to one specific promoter. To see all promoter regions, use
the full (unfiltered) Ensembl or RefSeq annotation track.

---

## Common pitfalls

| Problem | Symptom | Fix |
|---------|---------|-----|
| Tracks show nothing | Loaded successfully but blank | Chromosome naming mismatch: check chr1 vs 1 |
| Tracks show nothing | Only at certain chromosomes | Some files may only contain data for specific chromosomes |
| All tracks look the same height | Cannot compare across stages | Set the same Data Range for all tracks of the same assay type |
| Gene track shows too many isoforms | Cluttered view | Right-click > Collapsed, or load MANE-filtered GTF |
| IGV crashes when loading GTF | chrName is null error | You loaded GTF as genome instead of as track |
| IGV is very slow | Laggy scrolling | Too many tracks loaded, or tracks are on a network drive. Copy to local disk. |
| Autoscale makes peaks look equal | Misleading visual comparison | Set manual data range instead of autoscale |
| Cannot find gene | Search returns nothing | Check spelling, try official symbol instead of alias |

---

## Dataset reference

Zhu Y, Liu J, Wang J, et al. Integrative transcriptomic and epigenomic
analysis identifies BCL6B as a novel regulator of human pluripotent stem
cell to endothelial differentiation. Protein and Cell (2025).
GEO: GSE186755
