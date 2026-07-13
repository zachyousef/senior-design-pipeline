#!/usr/bin/env python3
"""Figure 1 - pipeline overview and the reproducibility (accept/reject) framework.
Reproducible schematic for the manuscript. Run: python figures_build_pipeline_schematic.py
Output: manuscript/figures/Figure_1_pipeline_schematic.(png|pdf)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
FIGDIR = os.path.join(os.path.dirname(__file__), "manuscript", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# palette
C_DATA="#DCE6F1"; C_EMB="#CFE8DC"; C_BENCH="#FCE7C8"; C_TEST="#E7D9EF"
C_ACCEPT="#1B9E77"; C_REJECT="#D95F02"; C_APP="#F2F2F2"; EDGE="#444444"

fig, ax = plt.subplots(figsize=(9.2, 12.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 14); ax.axis("off")

def box(x, y, w, h, text, fc, fs=10, bold=False, tc="black", ec=EDGE, r=0.12):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.02,rounding_size={r}",
                 linewidth=1.2, edgecolor=ec, facecolor=fc, zorder=2))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=tc, zorder=3, wrap=True)

def arrow(x1, y1, x2, y2, color=EDGE, lw=1.6, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=16,
                 linewidth=lw, color=color, zorder=1, shrinkA=2, shrinkB=2))

def label(x, y, text, fs=10, color="#333333", bold=True, ha="left"):
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs, color=color,
            fontweight="bold" if bold else "normal", zorder=3)

# ---------- Section A: data -> embeddings -> benchmark ----------
label(0.2, 13.7, "A   Build embeddings and benchmark clustering", fs=12)
box(1.4, 12.5, 7.2, 0.85, "Six mouse-retina scRNA-seq datasets (P15–P91)\nQC + scDblFinder doublet removal + SCTransform + Harmony integration", C_DATA, fs=9.5)
arrow(5, 12.5, 5, 12.1)
box(1.0, 11.25, 3.7, 0.8, "PCA embedding\n(linear)", C_EMB, fs=9.5, bold=True)
box(5.3, 11.25, 3.7, 0.8, "scDHA embedding\n(nonlinear)", C_EMB, fs=9.5, bold=True)
ax.text(5, 11.05, "Harmony re-applied to each embedding", ha="center", va="top", fontsize=8, style="italic", color="#555")
arrow(2.85, 11.25, 4.3, 10.55); arrow(7.15, 11.25, 5.7, 10.55)
box(1.4, 9.75, 7.2, 0.75, "Benchmark 7 (method × k) configurations on 5 non-redundant metrics\nROGUE  ·  balance  ·  per-cluster Jaccard  ·  marker specificity  ·  silhouette", C_BENCH, fs=9)
arrow(5, 9.75, 5, 9.35)

# ---------- Section B: two stability tests + decision ----------
label(0.2, 9.05, "B   Test the selected partition for two kinds of reproducibility", fs=12)
box(0.8, 7.9, 3.9, 0.95, "Within-embedding stability\nbootstrap-resample cells,\nre-cluster (fixed embedding)", C_TEST, fs=9)
box(5.3, 7.9, 3.9, 0.95, "Across-embedding stability\nre-train scDHA (5 seeds),\nre-cluster", C_TEST, fs=9, bold=True)
arrow(2.75, 7.9, 4.2, 7.15); arrow(7.25, 7.9, 5.8, 7.15)
# decision box
box(2.4, 6.15, 5.2, 0.95, "As k increases, do across-embedding reproducibility\nAND cluster quality rise together?", "#FFFFFF", fs=9.5, bold=True, ec="#222222")
arrow(3.6, 6.15, 2.6, 5.5, color=C_ACCEPT); arrow(6.4, 6.15, 7.4, 5.5, color=C_REJECT)
ax.text(2.7, 5.75, "YES", fontsize=10, fontweight="bold", color=C_ACCEPT)
ax.text(7.0, 5.75, "NO", fontsize=10, fontweight="bold", color=C_REJECT)
box(0.7, 4.55, 3.9, 0.9, "ACCEPT\nreproducible sub-structure", C_ACCEPT, fs=10, bold=True, tc="white", ec=C_ACCEPT)
box(5.4, 4.55, 3.9, 0.9, "REJECT\nover-clustering (embedding-dependent)", C_REJECT, fs=9.5, bold=True, tc="white", ec=C_REJECT)

# ---------- Section C: applications / validation ----------
label(0.2, 4.05, "C   Applications and ground-truth validation", fs=12)
box(0.5, 2.7, 2.85, 1.05, "Phase 1\nWhole retina → recover\n8 cell types\nARI 0.91 (label-free)", C_APP, fs=8.6)
box(3.6, 2.7, 2.85, 1.05, "Phase 2\nBipolar cells → ACCEPT\ncross-seed ARI → 0.93 at ~15 types\n(known subtypes)", "#DDF0E8", fs=8.2, bold=False, ec=C_ACCEPT)
box(6.7, 2.7, 2.85, 1.05, "Phase 3\nRods → REJECT\ncross-seed ARI 0.69\n(one dominant state)", "#FBE5D6", fs=8.6, bold=False, ec=C_REJECT)
arrow(2.6, 4.55, 1.9, 3.75, color=C_ACCEPT); arrow(2.6, 4.55, 5.0, 3.75, color=C_ACCEPT)
arrow(7.35, 4.55, 8.1, 3.75, color=C_REJECT)
box(1.4, 1.15, 7.2, 1.05, "Synthetic ground-truth validation (Splatter)\nRecovers real structure at the true k — down to a 5% rare subpopulation and weak separation\n(ARI vs truth 0.97–1.00)  ·  rejects null data (ARI vs truth 0)", "#EFE7F5", fs=8.8, ec="#7E57C2")
arrow(5, 2.7, 5, 2.2, color="#7E57C2")

plt.tight_layout(pad=0.5)
for ext in ("png","pdf"):
    fig.savefig(os.path.join(FIGDIR, f"Figure_1_pipeline_schematic.{ext}"), dpi=300 if ext=="png" else None, bbox_inches="tight")
print("Saved Figure_1_pipeline_schematic.png/pdf")
