"""
=============================================================================
FPKM EXPRESSION ANALYSIS — hPSC-to-EC Differentiation (GSE186755)
=============================================================================

PURPOSE:
  Visualize gene expression changes across 5 stages of directed endothelial
  differentiation (hESC -> VMC -> EPC-1 -> EPC-2 -> EC) for four functional
  gene groups: stemness markers, EC surface markers, EC lineage TFs, and
  IFN-gamma pathway genes.

OUTPUTS:
  - figure_fpkm.pdf  (vector, publication-ready)
  - figure_fpkm.png  (raster, for quick preview)

DATASET:
  Zhu Y, Liu J, Wang J, et al. Protein & Cell (2025)
  GEO: GSE186755

HOW TO RUN:
  source ~/bioenv/bin/activate
  python3 plot_fpkm_annotated.py

DEPENDENCIES:
  pip install pandas numpy matplotlib seaborn scipy
=============================================================================
"""

# =========================================================================
# SECTION 1: IMPORT LIBRARIES
# =========================================================================
# pandas       — reads and manipulates the FPKM data table
# numpy        — math operations (log2 transformation)
# matplotlib   — core plotting engine
# seaborn      — makes heatmaps easier with better defaults
# scipy        — statistical functions (unused now but available)

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')        # Write to file without opening a window
                             # Essential for WSL/server environments
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import seaborn as sns


# =========================================================================
# SECTION 2: FILE PATHS — CHANGE THESE TO MATCH YOUR SETUP
# =========================================================================
# Windows via WSL example:
#   FILE = "/mnt/c/Users/YOUR_NAME/Downloads/GSE186755_hPSC-EC_fpkm.txt"
#   OUT  = "/mnt/c/Users/YOUR_NAME/Downloads/"
# Mac/Linux example:
#   FILE = "/home/YOUR_NAME/data/GSE186755_hPSC-EC_fpkm.txt"
#   OUT  = "/home/YOUR_NAME/data/"

FILE = "/path/to/your/GSE186755_hPSC-EC_fpkm.txt"
OUT  = "/path/to/your/output/folder/"


# =========================================================================
# SECTION 3: DEFINE GENE LISTS
# =========================================================================
# Genes grouped by biological function.
# Names must EXACTLY match the gene symbols in the FPKM file.
# To add/remove genes, edit these lists and add matching entries
# in COLORS, LSTYLE, and MARKERS below.

GENES_STEM = ["POU5F1", "SOX2", "NANOG"]
# POU5F1 = OCT4, master pluripotency transcription factor
# SOX2   = cooperates with OCT4 to maintain stemness
# NANOG  = reinforces pluripotency downstream of OCT4/SOX2

GENES_ECMARKER = ["PECAM1", "CDH5", "KDR", "VWF"]
# PECAM1 = CD31, pan-endothelial surface marker
# CDH5   = CD144, VE-cadherin, endothelial adherens junctions
# KDR    = VEGFR2, main VEGF receptor on endothelial cells
# VWF    = von Willebrand factor, mature EC marker

GENES_ECTF = ["ETV2", "ERG", "FLI1"]
# ETV2 = early EC pioneer TF, transiently expressed during commitment
# ERG  = sustained EC master TF, maintains EC identity long-term
# FLI1 = ETS family member, co-occupies EC enhancers with ERG

GENES_IFN = ["RELA", "IRF1", "STAT1", "CIITA"]
# RELA  = NF-kB p65, co-activator at CIITA pIV promoter
# IRF1  = interferon regulatory factor 1, induced by STAT1
# STAT1 = primary signal transducer downstream of IFN-gamma receptor
# CIITA = master regulator of MHC class II expression

ALL_GENES = GENES_STEM + GENES_ECMARKER + GENES_ECTF + GENES_IFN


# =========================================================================
# SECTION 4: SAMPLE-TO-STAGE MAPPING
# =========================================================================
# Maps each column name in the FPKM file to its differentiation stage.
# This allows grouping replicates for averaging.

STAGE_MAP = {
    "hESC.rep1":  "hESC",    "hESC.rep2":  "hESC",
    "VMC.rep1":   "VMC",     "VMC.rep2":   "VMC",
    "EPC-1.rep1": "EPC-1",   "EPC-1.rep2": "EPC-1",
    "EPC-2.rep1": "EPC-2",   "EPC-2.rep2": "EPC-2",
    "EC.rep1":    "EC",      "EC.rep2":    "EC",
}

# Order on x-axis (left to right = differentiation progression)
STAGE_ORDER = ["hESC", "VMC", "EPC-1", "EPC-2", "EC"]


# =========================================================================
# SECTION 5: VISUAL STYLING
# =========================================================================
# Color families:  Purple=stemness, Blue=EC markers,
#                  Green=EC TFs,    Red/orange=IFN-gamma
COLORS = {
    "POU5F1": "#3A0CA3", "SOX2": "#7209B7", "NANOG": "#C77DFF",
    "PECAM1": "#185FA5", "CDH5": "#378ADD", "KDR": "#85B7EB", "VWF": "#042C53",
    "ETV2": "#80B918", "ERG": "#007F5F", "FLI1": "#2DC653",
    "RELA": "#E63946", "IRF1": "#F4A261", "STAT1": "#E9C46A", "CIITA": "#444441",
}

# Line styles: vary within each group for distinction
#   "-"            = solid
#   "--"           = dashed
#   (0,(3,1,1,1))  = dash-dot
#   (0,(5,2))      = long dash
LSTYLE = {
    "POU5F1": "-", "SOX2": "--", "NANOG": (0,(3,1,1,1)),
    "PECAM1": "-", "CDH5": "--", "KDR": (0,(3,1,1,1)), "VWF": (0,(5,2)),
    "ETV2": "-", "ERG": "--", "FLI1": (0,(3,1,1,1)),
    "RELA": "-", "IRF1": "--", "STAT1": (0,(3,1,1,1)), "CIITA": (0,(5,2)),
}

# Marker shapes: "o"=circle, "s"=square, "^"=triangle, "D"=diamond
MARKERS = {
    "POU5F1": "o", "SOX2": "s", "NANOG": "^",
    "PECAM1": "o", "CDH5": "s", "KDR": "^", "VWF": "D",
    "ETV2": "o", "ERG": "s", "FLI1": "^",
    "RELA": "o", "IRF1": "s", "STAT1": "^", "CIITA": "D",
}

# Human-readable legend labels
GENE_LABELS = {
    "POU5F1": "POU5F1 (OCT4)",
    "PECAM1": "PECAM1 (CD31)",
    "CDH5": "CDH5 (CD144)",
    "KDR": "KDR (VEGFR2)",
    "RELA": "RELA (NF-\u03baB)",
}


# =========================================================================
# SECTION 6: READ AND FILTER DATA
# =========================================================================
# Read the tab-separated file. index_col=0 = first column is gene names.
df = pd.read_csv(FILE, sep="\t", index_col=0)

# Keep only genes of interest and order them
df = df.loc[df.index.isin(ALL_GENES)].reindex(ALL_GENES)

# Log2 transform: adding 1 prevents log2(0) = -infinity
log2df = np.log2(df + 1)


# =========================================================================
# SECTION 7: GLOBAL PLOT SETTINGS (Nature/Cell style)
# =========================================================================
plt.rcParams.update({
    "font.family": "DejaVu Sans",   # sans-serif (Arial equivalent)
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,              # TrueType fonts in PDF (journal req)
    "ps.fonttype": 42,
})


# =========================================================================
# SECTION 8: CREATE FIGURE LAYOUT
# =========================================================================
fig = plt.figure(figsize=(7.5, 14))  # 7.5" wide = Nature single column
gs = gridspec.GridSpec(
    5, 1,                            # 5 rows, 1 column
    figure=fig, hspace=0.55,
    height_ratios=[1.3, 1, 1, 1, 1],
    top=0.96, bottom=0.04,
    left=0.11, right=0.74,          # right=0.74 leaves room for legends
)


# =========================================================================
# SECTION 9: PANEL A — HEATMAP (log2 fold-change vs hESC)
# =========================================================================
# For each gene: log2((FPKM_sample + 1) / (mean_hESC_FPKM + 1))
# This sets hESC as the zero baseline. Red=up, blue=down vs hESC.

ax1 = fig.add_subplot(gs[0])

# Mean of hESC replicates as the reference denominator
hesc_cols = [c for c in df.columns if c.startswith("hESC")]
hesc_mean = df[hesc_cols].mean(axis=1)

# Compute fold-change for every sample
fc_mat = pd.DataFrame(index=df.index, columns=df.columns)
for col in df.columns:
    fc_mat[col] = np.log2((df[col] + 1) / (hesc_mean + 1))
fc_mat = fc_mat.astype(float)

# Order columns by differentiation stage
col_order = [c for s in STAGE_ORDER for c in STAGE_MAP if STAGE_MAP[c] == s]
fc_mat = fc_mat[col_order]
col_labels = [f"{STAGE_MAP[c]}\n{c.split('.')[1]}" for c in col_order]

# Symmetric color scale, capped at 10
vmax = min(np.abs(fc_mat.values).max(), 10)

# Row labels with common names
hm_labels = [GENE_LABELS.get(g, g) for g in ALL_GENES]

# Draw heatmap
sns.heatmap(
    fc_mat, ax=ax1, cmap="RdBu_r",
    center=0, vmin=-vmax, vmax=vmax,
    linewidths=0.4, linecolor="white",
    xticklabels=col_labels, yticklabels=hm_labels,
    cbar_kws={"shrink": 0.5, "label": "log2 FC vs hESC",
              "orientation": "vertical"},
)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0, ha="center", fontsize=7)
ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, fontsize=7)
ax1.tick_params(left=False, bottom=False)

# Stage labels and white divider lines
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
             fontsize=8, fontweight="bold",
             transform=ax1.get_xaxis_transform())
    if i > 0:
        ax1.axvline(stage_positions[i], color="white", lw=2)

# White lines between gene groups
group_sizes = [len(GENES_STEM), len(GENES_ECMARKER), len(GENES_ECTF), len(GENES_IFN)]
running = 0
for size in group_sizes[:-1]:
    running += size
    ax1.axhline(running, color="white", lw=2)

ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.set_title("A", loc="left", fontweight="bold", fontsize=10, pad=6)


# =========================================================================
# SECTION 10: PREPARE DATA FOR LINE PLOTS
# =========================================================================
# Convert wide table to long format for easier plotting.

long = []
for col in df.columns:
    stage = STAGE_MAP.get(col)
    if not stage:
        continue
    for gene in ALL_GENES:
        if gene in df.index:
            long.append({
                "gene": gene,
                "stage": stage,
                "log2fpkm": np.log2(df.loc[gene, col] + 1),
            })

long = pd.DataFrame(long)
long["stage"] = pd.Categorical(long["stage"],
                               categories=STAGE_ORDER, ordered=True)

# Mean per gene per stage (averaging replicates)
means = long.groupby(["gene", "stage"], observed=True)["log2fpkm"].mean().reset_index()

# Map stages to x-positions
x_pos = {s: i for i, s in enumerate(STAGE_ORDER)}


# =========================================================================
# SECTION 11: REUSABLE LINE PLOT FUNCTION
# =========================================================================
# Called once per panel (B-E). Draws lines + individual replicate dots.

def draw_panel(ax, gene_list, ylabel, title_letter, group_label):

    for gene in gene_list:
        gm = means[means.gene == gene].sort_values("stage")
        xs = [x_pos[s] for s in gm.stage]
        ys = gm.log2fpkm.values

        # Line connecting means
        ax.plot(xs, ys, color=COLORS[gene], linewidth=1.6,
                linestyle=LSTYLE[gene], zorder=2)

        # Individual replicate dots
        gr = long[long.gene == gene]
        for _, row in gr.iterrows():
            ax.scatter(x_pos[row.stage], row.log2fpkm,
                       color=COLORS[gene], s=30, zorder=3,
                       marker=MARKERS[gene],
                       edgecolors="white", linewidths=0.5)

    # Formatting
    ax.set_xticks(range(len(STAGE_ORDER)))
    ax.set_xticklabels(STAGE_ORDER, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlim(-0.4, len(STAGE_ORDER) - 0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ["left", "bottom"]:
        ax.spines[sp].set_linewidth(0.6)
    ax.set_title(title_letter, loc="left", fontweight="bold", fontsize=10, pad=6)

    # Group label above plot
    ax.text(0.5, 1.02, group_label, transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9,
            fontstyle="italic", color="#5F5E5A")

    # Legend outside plot
    labels = [GENE_LABELS.get(g, g) for g in gene_list]
    handles = [
        mlines.Line2D([], [], color=COLORS[g], linewidth=1.6,
                      linestyle=LSTYLE[g], marker=MARKERS[g],
                      markersize=5, markeredgecolor='white',
                      markeredgewidth=0.5, label=lab)
        for g, lab in zip(gene_list, labels)
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7.5,
              bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0, handlelength=2.2, handletextpad=0.6)


# =========================================================================
# SECTION 12: DRAW PANELS B-E
# =========================================================================

ax2 = fig.add_subplot(gs[1])
draw_panel(ax2, GENES_STEM, "log2(FPKM + 1)", "B", "Stemness markers")

ax3 = fig.add_subplot(gs[2])
draw_panel(ax3, GENES_ECMARKER, "log2(FPKM + 1)", "C",
           "Endothelial surface markers")

ax4 = fig.add_subplot(gs[3])
draw_panel(ax4, GENES_ECTF, "log2(FPKM + 1)", "D",
           "EC lineage transcription factors")

ax5 = fig.add_subplot(gs[4])
draw_panel(ax5, GENES_IFN, "log2(FPKM + 1)", "E",
           "IFN-\u03b3 pathway / immune response")


# =========================================================================
# SECTION 13: SAVE
# =========================================================================
fig.savefig(OUT + "figure_fpkm.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUT + "figure_fpkm.png", dpi=300, bbox_inches="tight")
print("Saved: figure_fpkm.pdf and figure_fpkm.png")
