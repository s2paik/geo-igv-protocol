# Ensembl vs UCSC Chromosome Naming

This is the single most common source of silent failures when working
with genomics data in IGV. Understanding it saves hours of debugging.

---

## The problem

Two major genome databases name chromosomes differently:

| Ensembl | UCSC | Same DNA? |
|---------|------|-----------|
| 1 | chr1 | Yes, identical sequence |
| 2 | chr2 | Yes |
| X | chrX | Yes |
| Y | chrY | Yes |
| MT | chrM | Yes (also different name) |

The underlying DNA sequence is identical. The only difference is whether
the chromosome name has "chr" in front of it. But this tiny difference
causes tools like IGV to completely fail — silently, with no error message.

---

## Why it matters in practice

When you load a bigWig track into IGV, IGV tries to match the chromosome
names in the track file to the chromosome names in the reference genome.

- Track says "chr1" + Reference says "chr1" = data displays correctly
- Track says "chr1" + Reference says "1" = blank track, no data, no error
- Track says "1" + Reference says "chr1" = blank track, no data, no error

This is especially frustrating because IGV gives you no warning. The track
loads without error, it just shows nothing. Many beginners assume their
data is empty or the download failed, when the actual problem is a three-
character naming mismatch.

---

## Which convention does each source use?

| Source | Convention | Example |
|--------|-----------|---------|
| Ensembl genome/GTF | Numeric (no chr) | 1, 2, X, MT |
| UCSC genome browser | chr-prefixed | chr1, chr2, chrX, chrM |
| IGV built-in hg38 | chr-prefixed | chr1, chr2, chrX, chrM |
| GENCODE GTF | chr-prefixed | chr1, chr2, chrX, chrM |
| Most GEO bigWig files | chr-prefixed (usually aligned to UCSC) | chr1, chr2 |
| RefSeq | chr-prefixed | chr1, chr2 |

---

## How to check which convention a file uses

### For a GTF file (in WSL terminal):

    head -1 your_file.gtf

If the line starts with "1" or "X", it is Ensembl-style (no chr prefix).
If it starts with "chr1" or "chrX", it is UCSC-style.

### For a bigWig file:

Load it in IGV. If tracks are blank against hg38 (chr-prefixed), the file
probably uses numeric naming. If tracks display correctly, naming matches.

### For a FASTA file:

    grep "^>" your_file.fa | head -5

Shows the sequence headers. Look for ">1" vs ">chr1".

---

## How to convert between conventions

### Add chr prefix (Ensembl to UCSC):

    sed 's/^\([0-9XYM]\)/chr\1/' input.gtf > output.gtf

What this does line by line:
- sed = stream editor, processes text line by line
- s/ = substitution command
- ^\([0-9XYM]\) = at the start of a line, match any digit, X, Y, or M
- /chr\1/ = replace with "chr" followed by whatever was matched
- So "1" at the start of a line becomes "chr1"

### Remove chr prefix (UCSC to Ensembl):

    sed 's/^chr//' input.gtf > output.gtf

Simply removes "chr" from the beginning of each line.

---

## The safe approach

Pick one convention and stick with it throughout your analysis:

1. Check what your data files use (most GEO data uses chr-prefixed)
2. Use a matching reference genome (IGV built-in hg38 is chr-prefixed)
3. If loading custom annotations (Ensembl GTF), convert them to match

For this repository, we use IGV built-in hg38 (chr-prefixed) as the
reference and convert Ensembl GTFs to chr-prefixed with sed before loading.

---

## Why do two conventions exist?

Historical reasons. UCSC Genome Browser was developed at UC Santa Cruz
and chose to add "chr" to chromosome names early on. Ensembl was developed
at the European Bioinformatics Institute (EBI) and used plain numbers.
Neither side changed because millions of existing files and pipelines
depend on each convention. The bioinformatics community has lived with
this inconsistency for over 20 years.

GENCODE, which provides the most comprehensive human gene annotation,
uses chr-prefixed naming despite being produced by the Ensembl team —
a pragmatic acknowledgment that most tools expect it.
