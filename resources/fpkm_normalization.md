# FPKM Normalization: What It Does and Does Not Do

Understanding what FPKM normalizes for — and what it does not —
is essential for interpreting expression data correctly and knowing
when a more rigorous method is needed.

---

## What FPKM stands for

Fragments Per Kilobase of transcript per Million mapped reads.

Each word in the name describes a normalization step:

- Per Million mapped reads: corrects for how deeply the sample was
  sequenced. A sample with 50 million reads would otherwise show
  higher counts for every gene compared to a sample with 25 million reads.

- Per Kilobase of transcript: corrects for gene length. A 10 kb gene
  will naturally generate more reads than a 1 kb gene, even if both
  are expressed at the same level. Dividing by length makes genes of
  different sizes comparable within a sample.

---

## What FPKM normalizes (within-sample)

FPKM makes different genes comparable within a single sample:

    "Is gene A expressed more than gene B in this sample?"

This is a valid question to ask with FPKM. If gene A has FPKM=100
and gene B has FPKM=10 in the same sample, gene A is likely more
abundantly transcribed (after accounting for their length difference).

---

## What FPKM does NOT normalize (across samples)

FPKM does NOT reliably answer:

    "Is gene A expressed more in sample 1 than in sample 2?"

The problem is that FPKM is a relative measure. All FPKM values in a
sample are constrained to sum to a fixed total. This means:

- If a few genes are extremely highly expressed in one sample (e.g.
  pluripotency genes in hESCs), they mathematically suppress the FPKM
  of everything else in that sample

- A gene can appear to "go up" across conditions partly because the
  denominator changed, not because its actual transcription increased

- During dramatic biological transitions (like hESC-to-EC differentiation),
  the total RNA composition shifts substantially, violating the assumption
  that total RNA output is constant

---

## TPM vs FPKM

TPM (Transcripts Per Million) is a related measure that normalizes in
a slightly different order:

- FPKM: divide by depth first, then by length
- TPM: divide by length first, then by depth

TPM has the advantage that values across samples sum to the same total
(1 million), making cross-sample comparison slightly more consistent.
However, TPM shares the same fundamental limitation as FPKM — it is
still a relative/compositional measure.

For practical purposes, the difference between FPKM and TPM is minor.
Neither is suitable for rigorous cross-sample statistical testing.

---

## When FPKM is good enough

- Identifying which genes are expressed vs silent in a sample
- Qualitative trends: "gene X goes up during differentiation"
- Heatmaps showing relative patterns (especially with per-gene
  normalization like log2 fold-change vs a reference stage)
- Exploratory analysis before investing in rigorous methods

---

## When you need something better

- Formal differential expression testing (p-values, fold-changes
  with confidence intervals)
- Quantitative claims: "gene X is 50-fold higher in EC than hESC"
- Comparing expression across datasets or experiments
- Any analysis going into a publication's main statistical claims

---

## The gold standard: DESeq2 / edgeR

For rigorous cross-sample comparison, use raw read counts (not FPKM)
analyzed with DESeq2 or edgeR. These tools:

- Apply median-of-ratios (DESeq2) or TMM (edgeR) normalization, which
  explicitly corrects for the compositional bias that FPKM ignores
- Model biological variability between replicates
- Provide proper statistical tests with adjusted p-values
- Account for the mean-variance relationship in count data

These require raw count matrices, not FPKM tables. If the GEO dataset
only provides FPKM (as with GSE186755), you would need to download raw
FASTQ files from SRA and re-run the alignment and counting pipeline.

---

## Z-score artifacts on low-expression genes

A specific pitfall we encountered in this project:

Z-scoring normalizes each gene to its own mean and standard deviation
across all samples. This works well for genes with real dynamic range
(e.g. POU5F1 going from 1817 to 4 FPKM). But for near-zero genes like
CIITA (max 0.09 FPKM), z-scoring amplifies biological noise:

    CIITA FPKM values: 0.01, 0.01, 0.09, 0.08, 0.03, 0.03, 0.01, 0.01, 0.01, 0.07

    Z-scored: these tiny fluctuations become -0.8, -0.8, +1.9, +1.6, ...

The z-scored heatmap makes it look like CIITA peaks at VMC — but the
absolute difference is 0.07 FPKM, which is biologically meaningless.

### Solution: log2 fold-change vs a reference stage

Instead of z-scoring, normalize each gene against its own value in a
reference condition (e.g. hESC):

    log2((FPKM_sample + 1) / (mean_FPKM_hESC + 1))

The +1 pseudocount dampens near-zero divisions. CIITA correctly appears
flat (near-white) because all values divided by a near-zero reference
remain near zero. This is the normalization used in Section 03.

---

## Practical advice for this dataset

The GSE186755 FPKM data is suitable for:
- Identifying gene expression trends across differentiation (qualitative)
- Heatmaps with log2 fold-change normalization
- Supporting biological hypotheses with correlative evidence

Frame your claims as relative trends, not absolute quantification:
- Good: "ERG expression increases during endothelial differentiation"
- Avoid: "ERG is 137-fold higher in EC than hESC"

Note the limitation in figure legends:
"Expression shown as FPKM; cross-stage comparisons reflect relative
trends and were not re-normalized across libraries."
