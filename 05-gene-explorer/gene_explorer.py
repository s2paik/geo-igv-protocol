#!/usr/bin/env python3
"""
=============================================================================
GENE EXPLORER — Interactive Gene Expression Lookup (CLI version)
=============================================================================

PURPOSE:
  Search any gene by official symbol or common alias (CD31, OCT4, VEGFR2)
  and generate a heatmap + grouped line plot showing expression across the
  hPSC-to-EC differentiation series (GSE186755).

USAGE EXAMPLES:
  # Search by official gene symbols
  python3 gene_explorer.py --file data.txt --genes CIITA ERG FLI1 RELA

  # Search by common aliases (automatically resolved)
  python3 gene_explorer.py --file data.txt --genes CD31 OCT4 VEGFR2

  # Use preset biological groups
  python3 gene_explorer.py --file data.txt --group stemness ec_tfs

  # Mix groups and individual genes
  python3 gene_explorer.py --file data.txt --group ec_markers --genes CIITA

  # Specify output directory
  python3 gene_explorer.py --file data.txt --genes ERG FLI1 --out ./figures

WHAT IT PRODUCES:
  - gene_explorer_output.pdf (vector, publication-ready)
  - gene_explorer_output.png (raster, quick preview)
  - FPKM values table printed to terminal

DEPENDENCIES:
  pip install pandas numpy matplotlib seaborn

DATASET:
  Zhu Y, Liu J, Wang J, et al. Protein & Cell (2025). GEO: GSE186755
=============================================================================
"""


# =========================================================================
# SECTION 1: IMPORT LIBRARIES
# =========================================================================
# argparse — parses command-line arguments (--file, --genes, --group)
#            this is what lets you run the script with different options
#            without editing the code each time
# sys      — system functions, used here for sys.exit() to stop on errors
# pandas   — data manipulation (reads FPKM table, filters, groups)
# numpy    — math operations (log2 transformation)
# matplotlib — core plotting engine
# seaborn  — built on matplotlib, makes heatmaps easier

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')    # 'Agg' backend writes directly to file
                         # without opening a window on screen.
                         # Required for WSL/server environments.
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import seaborn as sns


# =========================================================================
# SECTION 2: GENE ALIASES
# =========================================================================
# Many genes have multiple names — official symbols (PECAM1),
# CD numbers (CD31), protein names (VE-cadherin), etc.
#
# This dictionary maps common aliases to official gene symbols.
# When you search for "CD31", the tool looks it up here and finds
# "PECAM1", which is the name used in the FPKM file.
#
# Format:  "ALIAS": "OFFICIAL_SYMBOL"
#
# To add your own aliases, just add a new line like:
#   "MY_ALIAS": "OFFICIAL_GENE_SYMBOL",
#
# All aliases are case-insensitive (CD31 = cd31 = Cd31).

GENE_ALIASES = {
    # Endothelial markers
    "CD31": "PECAM1",           # platelet endothelial cell adhesion molecule
    "CD144": "CDH5",            # VE-cadherin
    "VE-CADHERIN": "CDH5",
    "VEGFR2": "KDR",            # VEGF receptor 2
    "FLK1": "KDR",              # mouse name for KDR
    "VEGFR1": "FLT1",           # VEGF receptor 1
    "TIE2": "TEK",              # angiopoietin receptor
    "ENDOGLIN": "ENG",          # CD105
    "CD105": "ENG",
    "ENOS": "NOS3",             # endothelial nitric oxide synthase
    "CD146": "MCAM",
    "CD309": "KDR",             # another CD number for VEGFR2
    "CD202B": "TEK",            # another CD number for TIE2
    "CD141": "THBD",            # thrombomodulin

    # Pluripotency
    "OCT4": "POU5F1",           # the most common alias in stem cell biology
    "OCT3": "POU5F1",

    # Transcription factors
    "ER71": "ETV2",             # early endothelial TF
    "PU.1": "SPI1",             # myeloid master TF
    "SCL": "TAL1",              # hematopoietic TF
    "COUP-TFII": "NR2F2",       # venous identity TF

    # Immune / IFN-gamma pathway
    "P65": "RELA",              # NF-kB p65 subunit
    "NFKB": "RELA",
    "NF-KB": "RELA",
    "MHC-II": "CIITA",          # MHC class II transactivator
    "MHCII": "CIITA",
    "MHC2TA": "CIITA",
    "HLA-DR": "HLA-DRA",        # MHC class II alpha chain
    "B2-MICROGLOBULIN": "B2M",   # MHC class I component

    # Mesoderm
    "BRACHYURY": "T",           # primitive streak marker

    # Smooth muscle
    "SMA": "ACTA2",             # smooth muscle actin
    "ALPHA-SMA": "ACTA2",
    "SM22": "TAGLN",            # smooth muscle marker
    "DESMIN": "DES",

    # Hematopoietic
    "CD45": "PTPRC",            # pan-leukocyte marker
    "CD41": "ITGA2B",           # megakaryocyte marker

    # Other
    "VON WILLEBRAND": "VWF",
    "VEGF": "VEGFA",
}


# =========================================================================
# SECTION 3: PREDEFINED BIOLOGICAL GROUPS
# =========================================================================
# These are curated sets of genes that belong to the same biological
# pathway or cell-type signature. When you use --group stemness, the
# tool adds all genes in that group at once.
#
# Each group has:
#   - label:  display name (appears in the plot title)
#   - genes:  list of official gene symbols
#   - color:  hex color code for plotting (all genes in a group
#             share the same color, distinguished by line style)
#
# To add a new group, copy one of the existing entries and modify it.

KNOWN_GROUPS = {
    "stemness": {
        "label": "Stemness markers",
        "genes": ["POU5F1", "SOX2", "NANOG"],
        "color": "#7209B7",     # purple family
    },
    "ec_markers": {
        "label": "EC surface markers",
        "genes": ["PECAM1", "CDH5", "KDR", "VWF"],
        "color": "#185FA5",     # blue family
    },
    "ec_tfs": {
        "label": "EC lineage TFs",
        "genes": ["ETV2", "ERG", "FLI1"],
        "color": "#007F5F",     # green family
    },
    "ifn_pathway": {
        "label": "IFN-g pathway / immune",
        "genes": ["RELA", "IRF1", "STAT1", "CIITA"],
        "color": "#E63946",     # red family
    },
    "hematopoietic": {
        "label": "Hematopoietic",
        "genes": ["RUNX1", "SPI1", "CEBPB", "GATA1"],
        "color": "#F4A261",     # orange
    },
    "smooth_muscle": {
        "label": "Smooth muscle / pericyte",
        "genes": ["ACTA2", "TAGLN", "CNN1", "MYH11"],
        "color": "#80B918",     # yellow-green
    },
    "mesoderm": {
        "label": "Mesoderm",
        "genes": ["T", "MESP1", "MIXL1", "EOMES"],
        "color": "#D4A373",     # tan
    },
}


# =========================================================================
# SECTION 4: STAGE MAPPING
# =========================================================================
# The FPKM file has columns named "hESC.rep1", "VMC.rep2", etc.
# This dictionary tells the script which differentiation stage each
# column belongs to, so it can group replicates for averaging.
#
# If you use this tool with a DIFFERENT dataset, you must update this
# dictionary to match your column names. For example:
#   "Day0_rep1": "Day0", "Day0_rep2": "Day0",
#   "Day3_rep1": "Day3", "Day3_rep2": "Day3",

STAGE_MAP = {
    "hESC.rep1": "hESC",    "hESC.rep2": "hESC",     # undifferentiated
    "VMC.rep1": "VMC",      "VMC.rep2": "VMC",        # vascular mesoderm
    "EPC-1.rep1": "EPC-1",  "EPC-1.rep2": "EPC-1",    # early endo progenitor
    "EPC-2.rep1": "EPC-2",  "EPC-2.rep2": "EPC-2",    # late endo progenitor
    "EC.rep1": "EC",        "EC.rep2": "EC",           # mature endothelial
}

# The order stages appear on the x-axis (left = undifferentiated, right = mature)
STAGE_ORDER = ["hESC", "VMC", "EPC-1", "EPC-2", "EC"]

# All column names from the FPKM file (derived from STAGE_MAP)
STAGE_COLS = list(STAGE_MAP.keys())

# Human-readable labels for common genes (used in plot legends)
# Genes not listed here just use their official symbol as-is
GENE_LABELS = {
    "POU5F1": "POU5F1 (OCT4)",
    "PECAM1": "PECAM1 (CD31)",
    "CDH5": "CDH5 (CD144)",
    "KDR": "KDR (VEGFR2)",
    "RELA": "RELA (NF-kB)",
}


# =========================================================================
# SECTION 5: VISUAL STYLING
# =========================================================================
# Genes within the same group share a color but are distinguished by
# different line styles and marker shapes.
#
# LINE_STYLES: how the connecting line is drawn
#   "-"            = solid line
#   "--"           = dashed line
#   (0,(3,1,1,1))  = dash-dot pattern
#   (0,(5,2))      = long dashes
#
# MARKER_SHAPES: the dot shape at each data point
#   "o" = circle    "s" = square    "^" = triangle up
#   "D" = diamond   "v" = triangle down   "P" = plus

LINE_STYLES = ["-", "--", (0,(3,1,1,1)), (0,(5,2))]
MARKER_SHAPES = ["o", "s", "^", "D", "v", "P"]


# =========================================================================
# SECTION 6: HELPER FUNCTIONS
# =========================================================================

def resolve_gene(name, available_genes):
    """
    Resolve a gene name or alias to the official symbol in the dataset.

    Search order:
    1. Exact match in the dataset (case-sensitive)
    2. Case-insensitive match in the dataset
    3. Alias lookup (CD31 -> PECAM1)

    Parameters:
        name: string the user typed (could be "CD31", "pecam1", "PECAM1")
        available_genes: set of all gene names in the FPKM file

    Returns:
        The official gene symbol if found, or None if not found.
    """
    upper = name.upper()

    # 1. Exact match
    if name in available_genes:
        return name

    # 2. Case-insensitive match
    for g in available_genes:
        if g.upper() == upper:
            return g

    # 3. Alias lookup
    if upper in GENE_ALIASES:
        official = GENE_ALIASES[upper]
        if official in available_genes:
            return official

    return None


def classify_gene(gene):
    """
    Assign a gene to a biological group based on the predefined groups.

    If the gene is in any KNOWN_GROUPS list, returns that group's label
    and color. Otherwise returns "Other" with a gray color.

    This is how the tool automatically sorts genes into the right
    line plot panel without you having to specify.
    """
    for gname, ginfo in KNOWN_GROUPS.items():
        if gene in ginfo["genes"]:
            return ginfo["label"], ginfo["color"]
    return "Other", "#888780"


# =========================================================================
# SECTION 7: MAIN FUNCTION
# =========================================================================
# This is where the actual work happens. It:
# 1. Parses command-line arguments
# 2. Reads the FPKM file
# 3. Resolves gene names and aliases
# 4. Groups genes by biological category
# 5. Generates the heatmap and line plots
# 6. Saves the output figure

def main():
    # --- Parse command-line arguments ---
    # argparse handles the --file, --genes, --group, --out flags
    # When you type: python3 gene_explorer.py --genes CD31 OCT4
    # argparse puts ["CD31", "OCT4"] into args.genes
    parser = argparse.ArgumentParser(
        description="Gene Expression Explorer for GSE186755",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 gene_explorer.py --file data.txt --genes CIITA ERG FLI1 RELA
  python3 gene_explorer.py --file data.txt --genes CD31 OCT4 VEGFR2
  python3 gene_explorer.py --file data.txt --group stemness ec_tfs
  python3 gene_explorer.py --file data.txt --group ec_markers --genes CIITA

Available groups: stemness, ec_markers, ec_tfs, ifn_pathway,
                  hematopoietic, smooth_muscle, mesoderm
        """,
    )
    parser.add_argument("--file", required=True,
                        help="Path to FPKM text file (tab-separated)")
    parser.add_argument("--genes", nargs="+", default=[],
                        help="Gene names or aliases to plot")
    parser.add_argument("--group", nargs="+", default=[],
                        help="Preset group names to add")
    parser.add_argument("--out", default=".",
                        help="Output directory (default: current folder)")
    args = parser.parse_args()


    # --- Read the FPKM data ---
    # pd.read_csv reads the tab-separated file
    # sep="\t" tells pandas the columns are separated by tabs
    # index_col=0 means the first column (gene names) becomes the row index
    print(f"Reading {args.file}...")
    df = pd.read_csv(args.file, sep="\t", index_col=0)
    available = set(df.index)    # set of all ~36,000 gene symbols
    print(f"  {len(available):,} genes loaded")


    # --- Collect genes from groups and individual names ---
    # selected is a list of tuples: (gene_symbol, group_label, color)
    selected = []

    # Add genes from preset groups (--group stemness ec_tfs)
    for gname in args.group:
        if gname not in KNOWN_GROUPS:
            print(f"  WARNING: Unknown group '{gname}'.")
            print(f"  Available: {', '.join(KNOWN_GROUPS.keys())}")
            continue
        ginfo = KNOWN_GROUPS[gname]
        for g in ginfo["genes"]:
            # Only add if the gene exists in the dataset and isn't already added
            if g in available and g not in [s[0] for s in selected]:
                selected.append((g, ginfo["label"], ginfo["color"]))
                print(f"  + {g} ({ginfo['label']})")

    # Add individual genes (--genes CD31 OCT4 CIITA)
    for name in args.genes:
        resolved = resolve_gene(name, available)
        if resolved is None:
            print(f"  WARNING: '{name}' not found (checked aliases too)")
            continue
        if resolved in [s[0] for s in selected]:
            continue    # skip duplicates
        label, color = classify_gene(resolved)
        selected.append((resolved, label, color))
        # Show the user if an alias was resolved
        alias_note = f" (alias for {resolved})" if name.upper() != resolved.upper() else ""
        print(f"  + {resolved}{alias_note} [{label}]")

    # Stop if nothing was selected
    if not selected:
        print("No valid genes selected. Use --genes or --group.")
        sys.exit(1)

    print(f"\nPlotting {len(selected)} genes...")


    # --- Group genes by biological category for separate panels ---
    # groups = {"Stemness markers": {"genes": ["POU5F1","SOX2"], "color": "#7209B7"}, ...}
    groups = {}
    for gene, label, color in selected:
        if label not in groups:
            groups[label] = {"genes": [], "color": color}
        groups[label]["genes"].append(gene)


    # --- Prepare data subset ---
    gene_list = [s[0] for s in selected]
    sub = df.loc[df.index.isin(gene_list)].reindex(gene_list)

    # Calculate mean hESC FPKM for fold-change normalization
    hesc_cols = [c for c in sub.columns if c.startswith("hESC")]
    hesc_mean = sub[hesc_cols].mean(axis=1)


    # --- Plot settings (Nature/Cell journal style) ---
    plt.rcParams.update({
        "font.family": "DejaVu Sans",  # clean sans-serif font
        "font.size": 8,                # base font size
        "axes.linewidth": 0.6,         # thin axis borders
        "pdf.fonttype": 42,            # TrueType fonts in PDF (journal requirement)
        "ps.fonttype": 42,
    })


    # --- Create figure layout ---
    # The figure has 1 heatmap panel + 1 line plot per biological group
    n_groups = len(groups)
    fig = plt.figure(figsize=(7.5, 3.5 + n_groups * 2.5))
    gs = gridspec.GridSpec(
        1 + n_groups, 1,               # rows = 1 heatmap + N line plots
        figure=fig, hspace=0.55,       # vertical spacing between panels
        height_ratios=[1.3] + [1] * n_groups,  # heatmap slightly taller
        top=0.95, bottom=0.05,
        left=0.11, right=0.74,         # right margin leaves room for legends
    )


    # =====================================================================
    # PANEL A: HEATMAP (log2 fold-change vs hESC)
    # =====================================================================
    # For each gene and each sample, we compute:
    #   log2((FPKM_sample + 1) / (mean_hESC_FPKM + 1))
    #
    # This makes hESC the zero baseline (white in the heatmap).
    # Red = upregulated vs hESC, Blue = downregulated vs hESC.
    # The +1 pseudocount prevents log2(0) = -infinity for silent genes.

    ax1 = fig.add_subplot(gs[0])

    # Compute fold-change matrix
    fc_mat = pd.DataFrame(index=sub.index, columns=sub.columns)
    for col in sub.columns:
        fc_mat[col] = np.log2((sub[col] + 1) / (hesc_mean + 1))
    fc_mat = fc_mat.astype(float)

    # Order columns by differentiation stage
    col_order = [c for s in STAGE_ORDER for c in STAGE_MAP if STAGE_MAP[c] == s]
    fc_mat = fc_mat[col_order]

    # Two-line column labels: "hESC\nrep1"
    col_labels = [f"{STAGE_MAP[c]}\n{c.split('.')[1]}" for c in col_order]

    # Symmetric color scale capped at 10 (prevents extreme values from
    # washing out the color range for all other genes)
    vmax = min(np.abs(fc_mat.values).max(), 10)

    # Use human-readable labels where available
    hm_labels = [GENE_LABELS.get(g, g) for g in gene_list]

    # Draw the heatmap using seaborn
    sns.heatmap(fc_mat, ax=ax1, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                linewidths=0.4, linecolor="white",
                xticklabels=col_labels, yticklabels=hm_labels,
                cbar_kws={"shrink": 0.5, "label": "log2 FC vs hESC"})
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0, ha="center", fontsize=7)
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, fontsize=7)
    ax1.tick_params(left=False, bottom=False)
    ax1.set_title("A", loc="left", fontweight="bold", fontsize=10, pad=6)

    # Add bold stage labels above columns
    stage_positions, prev = [], None
    for i, c in enumerate(col_order):
        s = STAGE_MAP[c]
        if s != prev:
            stage_positions.append(i)
            prev = s
    stage_positions.append(len(col_order))
    for i in range(len(stage_positions) - 1):
        mid = (stage_positions[i] + stage_positions[i+1]) / 2
        ax1.text(mid, -0.8, STAGE_ORDER[i], ha="center", va="top",
                 fontsize=8, fontweight="bold", transform=ax1.get_xaxis_transform())
        if i > 0:
            ax1.axvline(stage_positions[i], color="white", lw=2)


    # =====================================================================
    # PREPARE LONG-FORMAT DATA FOR LINE PLOTS
    # =====================================================================
    # "Long format" means one row per gene per sample, instead of one row
    # per gene with all samples as columns. This is easier to work with
    # for plotting because each data point is its own row.
    #
    # Wide format (FPKM file):
    #   gene    | hESC.rep1 | hESC.rep2 | VMC.rep1 | ...
    #   POU5F1  | 1896      | 1738      | 192      | ...
    #
    # Long format (what we build):
    #   gene    | stage | log2fpkm
    #   POU5F1  | hESC  | 10.89
    #   POU5F1  | hESC  | 10.76
    #   POU5F1  | VMC   | 7.59

    long = []
    for col in sub.columns:
        stage = STAGE_MAP.get(col)
        if not stage:
            continue
        for gene in gene_list:
            if gene in sub.index:
                long.append({"gene": gene, "stage": stage,
                             "log2fpkm": np.log2(sub.loc[gene, col] + 1)})
    long = pd.DataFrame(long)
    long["stage"] = pd.Categorical(long["stage"], categories=STAGE_ORDER, ordered=True)

    # Mean per gene per stage (averaging the two replicates)
    means = long.groupby(["gene", "stage"], observed=True)["log2fpkm"].mean().reset_index()

    # Map stage names to numeric x-positions: hESC=0, VMC=1, EPC-1=2, etc.
    x_pos = {s: i for i, s in enumerate(STAGE_ORDER)}


    # =====================================================================
    # LINE PLOT PANELS (one per biological group)
    # =====================================================================
    # Each group gets its own subplot. Genes within a group share the same
    # base color but are distinguished by line style and marker shape.

    panel_letters = "BCDEFGHIJ"

    for pi, (group_label, group_info) in enumerate(groups.items()):
        ax = fig.add_subplot(gs[1 + pi])
        genes_in_group = group_info["genes"]
        base_color = group_info["color"]

        for gi, gene in enumerate(genes_in_group):
            # Get mean values for the connecting line
            gm = means[means.gene == gene].sort_values("stage")
            xs = [x_pos[s] for s in gm.stage]
            ys = gm.log2fpkm.values

            # Draw the connecting line (mean values)
            ax.plot(xs, ys, color=base_color, linewidth=1.6,
                    linestyle=LINE_STYLES[gi % len(LINE_STYLES)],
                    zorder=2)    # zorder=2: line draws behind dots

            # Overlay individual replicate data points
            # This shows the spread between replicates at each stage
            gr = long[long.gene == gene]
            for _, row in gr.iterrows():
                ax.scatter(x_pos[row.stage], row.log2fpkm,
                           color=base_color, s=30, zorder=3,
                           marker=MARKER_SHAPES[gi % len(MARKER_SHAPES)],
                           edgecolors="white", linewidths=0.5)

        # --- Format axes ---
        ax.set_xticks(range(len(STAGE_ORDER)))
        ax.set_xticklabels(STAGE_ORDER, fontsize=8)
        ax.set_ylabel("log2(FPKM + 1)", fontsize=8)
        ax.set_xlim(-0.4, len(STAGE_ORDER) - 0.6)
        ax.spines["top"].set_visible(False)     # remove top border
        ax.spines["right"].set_visible(False)    # remove right border

        # Panel letter (B, C, D, ...)
        if pi < len(panel_letters):
            ax.set_title(panel_letters[pi], loc="left",
                         fontweight="bold", fontsize=10, pad=6)

        # Italic group label centered above the plot
        ax.text(0.5, 1.02, group_label, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=9,
                fontstyle="italic", color="#5F5E5A")

        # --- Legend outside the plot ---
        # Each gene gets a legend entry showing its line style and marker
        labels = [GENE_LABELS.get(g, g) for g in genes_in_group]
        handles = [mlines.Line2D([], [], color=base_color, linewidth=1.6,
                                 linestyle=LINE_STYLES[i % len(LINE_STYLES)],
                                 marker=MARKER_SHAPES[i % len(MARKER_SHAPES)],
                                 markersize=5, markeredgecolor="white",
                                 markeredgewidth=0.5, label=lab)
                   for i, lab in enumerate(labels)]
        ax.legend(handles=handles, frameon=False, fontsize=7.5,
                  bbox_to_anchor=(1.02, 1), loc="upper left",
                  handlelength=2.2, handletextpad=0.6)


    # =====================================================================
    # SAVE THE FIGURE
    # =====================================================================
    # dpi=300 meets minimum resolution for most journals
    # bbox_inches="tight" crops whitespace around the figure

    out_pdf = f"{args.out}/gene_explorer_output.pdf"
    out_png = f"{args.out}/gene_explorer_output.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {out_pdf}")
    print(f"Saved: {out_png}")


    # =====================================================================
    # PRINT FPKM TABLE TO TERMINAL
    # =====================================================================
    # Shows the raw FPKM values so you can quickly check numbers
    # without opening the figure

    print(f"\n{'Gene':<12} {'Group':<25} ", end="")
    for s in STAGE_ORDER:
        print(f"{s:>8}", end="")
    print()
    print("-" * 72)

    for gene, label, color in selected:
        print(f"{gene:<12} {label:<25} ", end="")
        for s in STAGE_ORDER:
            cols = [c for c in STAGE_COLS if STAGE_MAP[c] == s]
            val = np.mean([sub.loc[gene, c] for c in cols if gene in sub.index])
            print(f"{val:>8.1f}", end="")
        print()


# =========================================================================
# SECTION 8: ENTRY POINT
# =========================================================================
# This is standard Python boilerplate. It means:
# "Only run main() if this script is executed directly,
#  not if it is imported as a module by another script."
#
# When you type: python3 gene_explorer.py --file data.txt --genes ERG
# Python sees __name__ == "__main__" and calls main().

if __name__ == "__main__":
    main()


# =========================================================================
# QUICK REFERENCE: COMMON MODIFICATIONS
# =========================================================================
#
# ADD A NEW GENE ALIAS:
#   Add to GENE_ALIASES dict in Section 2:
#   "MY_ALIAS": "OFFICIAL_SYMBOL",
#
# ADD A NEW BIOLOGICAL GROUP:
#   Add to KNOWN_GROUPS dict in Section 3:
#   "my_group": {
#       "label": "My Group Name",
#       "genes": ["GENE1", "GENE2"],
#       "color": "#HEX_COLOR",
#   },
#
# USE WITH A DIFFERENT DATASET:
#   Update STAGE_MAP in Section 4 to match your column names.
#   Update STAGE_ORDER to match your experimental stages.
#
# CHANGE COLORS:
#   Edit the "color" value in the relevant KNOWN_GROUPS entry.
#   Use https://colorbrewer2.org for colorblind-safe palettes.
