#!/usr/bin/env python3
"""
Publication-quality figure generation for GBS-Genomics-Pipeline.

All figures are derived exclusively from validated pipeline outputs.
If required data are unavailable, a report file is written instead of
generating placeholder or synthetic visualizations.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False


# ── Abricate TSV column indices (0-based) ──────────────────────────────
COL_GENE = 5
COL_PRODUCT = 13


def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def parse_abricate(tsv_path: Path) -> pd.DataFrame:
    """Parse an Abricate TSV into a DataFrame of gene hits."""
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return pd.DataFrame(columns=["sample", "gene", "product"])

    df = pd.read_csv(tsv_path, sep="\t", comment="#", header=0, dtype=str)
    if df.empty:
        return pd.DataFrame(columns=["sample", "gene", "product"])

    sample = tsv_path.stem.replace("_amr", "").replace("_vf", "")
    genes = df.iloc[:, COL_GENE].dropna().unique().tolist()
    products = df.iloc[:, COL_PRODUCT].dropna().unique().tolist() if df.shape[1] > COL_PRODUCT else [""] * len(genes)

    rows = []
    for g in genes:
        rows.append({"sample": sample, "gene": g, "product": ""})
    return pd.DataFrame(rows)


def load_gene_matrix(directory: Path, suffix: str) -> Tuple[pd.DataFrame, List[str]]:
    """Build a sample × gene presence/absence matrix from Abricate TSV files."""
    files = sorted(directory.glob(f"*{suffix}"))
    if not files:
        return pd.DataFrame(), []

    all_rows = []
    for f in files:
        df = parse_abricate(f)
        if not df.empty:
            all_rows.append(df)

    if not all_rows:
        return pd.DataFrame(), []

    combined = pd.concat(all_rows, ignore_index=True)
    samples = sorted(combined["sample"].unique())
    genes = sorted(combined["gene"].unique())

    matrix = pd.DataFrame(0, index=samples, columns=genes)
    for _, row in combined.iterrows():
        matrix.at[row["sample"], row["gene"]] = 1

    return matrix, genes


def save_figure(fig, outdir: Path, basename: str) -> None:
    """Save figure in PNG, SVG, and PDF formats."""
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"{basename}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)
    log(f"Saved {basename}.{{png,svg,pdf}}")


def write_skip_report(outdir: Path, figure_name: str, reason: str) -> None:
    report = outdir / "Reports" / f"{figure_name}_SKIPPED.txt"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        f"Figure: {figure_name}\n"
        f"Status: NOT GENERATED\n"
        f"Reason: {reason}\n"
        f"Timestamp: {datetime.now().isoformat()}\n"
    )
    log(f"Skipped {figure_name}: {reason}")


def plot_heatmap(matrix: pd.DataFrame, title: str, outdir: Path, basename: str,
                 cmap: str = "YlOrRd") -> bool:
    if matrix.empty or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        write_skip_report(outdir, basename, "No gene detection data available.")
        return False

    fig, ax = plt.subplots(figsize=(max(8, matrix.shape[1] * 0.4), max(4, matrix.shape[0] * 0.5)))
    sns.heatmap(matrix, cmap=cmap, cbar_kws={"label": "Present"}, linewidths=0.5,
                ax=ax, annot=matrix.shape[1] <= 20)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Sample")
    save_figure(fig, outdir, basename)
    matrix.to_csv(outdir / f"{basename}_matrix.tsv", sep="\t")
    return True


def plot_frequency(matrix: pd.DataFrame, title: str, outdir: Path, basename: str) -> bool:
    if matrix.empty or matrix.shape[1] == 0:
        write_skip_report(outdir, basename, "No gene detection data available.")
        return False

    freq = matrix.sum(axis=0).sort_values(ascending=False)
    freq_df = freq.reset_index()
    freq_df.columns = ["gene", "count"]
    freq_df.to_csv(outdir / f"{basename}_counts.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(max(8, len(freq) * 0.35), 5))
    freq.plot(kind="bar", ax=ax, color="steelblue", edgecolor="black")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Number of Samples")
    ax.tick_params(axis="x", rotation=45)
    save_figure(fig, outdir, basename)
    return True


def plot_cooccurrence(matrix: pd.DataFrame, title: str, outdir: Path, basename: str) -> bool:
    genes = [g for g in matrix.columns if matrix[g].sum() >= 2]
    if len(genes) < 2:
        write_skip_report(outdir, basename,
                          f"Need ≥2 genes detected in ≥2 samples; found {len(genes)} qualifying genes.")
        return False

    sub = matrix[genes]
    pairs = []
    for g1, g2 in combinations(genes, 2):
        co = ((sub[g1] == 1) & (sub[g2] == 1)).sum()
        if co > 0:
            pairs.append({"gene_1": g1, "gene_2": g2, "co_occurrence": co})

    if not pairs:
        write_skip_report(outdir, basename, "No co-occurring gene pairs found.")
        return False

    pair_df = pd.DataFrame(pairs)
    pair_df.to_csv(outdir / f"{basename}_pairs.tsv", sep="\t", index=False)

    n = len(genes)
    co_mat = pd.DataFrame(0, index=genes, columns=genes)
    for _, row in pair_df.iterrows():
        co_mat.at[row["gene_1"], row["gene_2"]] = row["co_occurrence"]
        co_mat.at[row["gene_2"], row["gene_1"]] = row["co_occurrence"]

    fig, ax = plt.subplots(figsize=(max(6, n * 0.5), max(6, n * 0.5)))
    sns.heatmap(co_mat, annot=True, fmt="d", cmap="Blues", ax=ax)
    ax.set_title(title, fontsize=14, fontweight="bold")
    save_figure(fig, outdir, basename)
    return True


def parse_mlst(directory: Path) -> pd.DataFrame:
    """Parse MLST TSV outputs into a structured DataFrame."""
    files = sorted(directory.glob("*_mlst.tsv"))
    rows = []
    for f in files:
        if f.stat().st_size == 0:
            continue
        sample = f.stem.replace("_mlst", "")
        try:
            content = f.read_text().strip()
            if not content:
                continue
            parts = content.split("\t")
            if len(parts) >= 3:
                rows.append({
                    "sample": sample,
                    "scheme": parts[1] if len(parts) > 1 else "unknown",
                    "ST": parts[2] if len(parts) > 2 else "unknown",
                    "alleles": "\t".join(parts[3:]) if len(parts) > 3 else "",
                })
        except Exception as e:
            log(f"Warning: could not parse MLST file {f}: {e}")

    return pd.DataFrame(rows)


def plot_mlst_distribution(mlst_df: pd.DataFrame, outdir: Path) -> bool:
    if mlst_df.empty:
        write_skip_report(outdir, "mlst_distribution", "No MLST results available.")
        return False

    mlst_df.to_csv(outdir / "mlst_summary.tsv", sep="\t", index=False)
    st_counts = mlst_df["ST"].value_counts()

    fig, ax = plt.subplots(figsize=(max(6, len(st_counts) * 0.6), 5))
    st_counts.plot(kind="bar", ax=ax, color="teal", edgecolor="black")
    ax.set_title("MLST Sequence Type Distribution", fontsize=14, fontweight="bold")
    ax.set_xlabel("Sequence Type (ST)")
    ax.set_ylabel("Number of Samples")
    save_figure(fig, outdir, "mlst_distribution")
    return True


def plot_gene_count_bars(matrix: pd.DataFrame, title: str, outdir: Path, basename: str) -> bool:
    if matrix.empty:
        write_skip_report(outdir, basename, "No gene detection data available.")
        return False

    counts = matrix.sum(axis=1).sort_values(ascending=False)
    counts_df = counts.reset_index()
    counts_df.columns = ["sample", "gene_count"]
    counts_df.to_csv(outdir / f"{basename}_counts.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(max(6, len(counts) * 0.5), 5))
    counts.plot(kind="bar", ax=ax, color="coral", edgecolor="black")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Gene Count")
    save_figure(fig, outdir, basename)
    return True


def plot_dashboard(amr_matrix, vf_matrix, mlst_df, outdir: Path) -> bool:
    if amr_matrix.empty and vf_matrix.empty and mlst_df.empty:
        write_skip_report(outdir, "sample_summary_dashboard", "No summary data available.")
        return False

    samples = sorted(set(
        list(amr_matrix.index if not amr_matrix.empty else []) +
        list(vf_matrix.index if not vf_matrix.empty else []) +
        list(mlst_df["sample"] if not mlst_df.empty else [])
    ))

    summary_rows = []
    for s in samples:
        summary_rows.append({
            "sample": s,
            "amr_genes": int(amr_matrix.loc[s].sum()) if s in amr_matrix.index else 0,
            "virulence_genes": int(vf_matrix.loc[s].sum()) if s in vf_matrix.index else 0,
            "ST": mlst_df.loc[mlst_df["sample"] == s, "ST"].values[0]
                  if not mlst_df.empty and s in mlst_df["sample"].values else "NA",
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "sample_summary_dashboard.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    summary.plot(x="sample", y="amr_genes", kind="bar", ax=axes[0], color="firebrick", legend=False)
    axes[0].set_title("AMR Genes per Sample")
    axes[0].tick_params(axis="x", rotation=45)

    summary.plot(x="sample", y="virulence_genes", kind="bar", ax=axes[1], color="darkorange", legend=False)
    axes[1].set_title("Virulence Genes per Sample")
    axes[1].tick_params(axis="x", rotation=45)

    if not mlst_df.empty:
        st_map = summary.set_index("sample")["ST"]
        axes[2].bar(range(len(samples)), [1] * len(samples),
                    color=[plt.cm.Set3(hash(str(st_map.get(s, "NA"))) % 12) for s in samples])
        axes[2].set_xticks(range(len(samples)))
        axes[2].set_xticklabels([f"{s}\n(ST:{st_map.get(s, 'NA')})" for s in samples], rotation=45, ha="right")
        axes[2].set_title("Sample ST Labels")
    else:
        axes[2].text(0.5, 0.5, "No MLST data", ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_title("MLST (unavailable)")

    fig.suptitle("Sample Summary Dashboard", fontsize=16, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, outdir, "sample_summary_dashboard")
    return True


def render_tree(tree_path: Path, annotation: dict, outdir: Path, basename: str,
                title: str, color_key: str = None) -> bool:
    if not HAS_ETE3:
        write_skip_report(outdir, basename, "ete3 is not installed.")
        return False
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        write_skip_report(outdir, basename, "Phylogenetic tree file is missing or empty.")
        return False

    try:
        t = Tree(str(tree_path), format=1)
    except Exception as e:
        write_skip_report(outdir, basename, f"Could not parse tree: {e}")
        return False

    for leaf in t.iter_leaves():
        clean = leaf.name.replace(".fna", "").replace(".fasta", "").replace(".fsa", "")
        nstyle = NodeStyle()
        nstyle["size"] = 8
        if color_key and clean in annotation:
            val = annotation[clean]
            nstyle["fgcolor"] = "red" if val and val != "None" and val != "NA" else "green"
        leaf.set_style(nstyle)
        if clean in annotation:
            label = annotation[clean]
            if label and label != "None":
                leaf.add_face(TextFace(f" [{label}]", fsize=8), column=0, position="branch-right")

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.title.add_face(TextFace(title, fsize=14), column=0)

    for ext in ("png", "svg", "pdf"):
        t.render(str(outdir / f"{basename}.{ext}"), tree_style=ts)
    log(f"Saved {basename}.{{png,svg,pdf}}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate publication figures from pipeline outputs.")
    parser.add_argument("--amr-dir", required=True, help="Directory containing *_amr.tsv files")
    parser.add_argument("--vf-dir", required=True, help="Directory containing *_vf.tsv files")
    parser.add_argument("--mlst-dir", required=True, help="Directory containing *_mlst.tsv files")
    parser.add_argument("--tree", required=True, help="Path to Newick tree file")
    parser.add_argument("--outdir", required=True, help="Output directory for figures")
    args = parser.parse_args()

    amr_dir = Path(args.amr_dir)
    vf_dir = Path(args.vf_dir)
    mlst_dir = Path(args.mlst_dir)
    tree_path = Path(args.tree)
    outdir = Path(args.outdir)
    reports_dir = outdir / "Reports"
    outdir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    log("Stage: Integration Visualization — loading pipeline outputs")
    log(f"  AMR directory:   {amr_dir}")
    log(f"  Virulence directory: {vf_dir}")
    log(f"  MLST directory:  {mlst_dir}")
    log(f"  Tree file:       {tree_path}")
    log(f"  Output directory: {outdir}")

    amr_matrix, _ = load_gene_matrix(amr_dir, "_amr.tsv")
    vf_matrix, _ = load_gene_matrix(vf_dir, "_vf.tsv")
    mlst_df = parse_mlst(mlst_dir)

    log(f"  Samples with AMR data: {len(amr_matrix)}")
    log(f"  Samples with VF data:  {len(vf_matrix)}")
    log(f"  Samples with MLST data: {len(mlst_df)}")

    generated = []
    skipped = []

    figs = [
        ("amr_heatmap", lambda: plot_heatmap(amr_matrix, "AMR Gene Presence/Absence Heatmap", outdir, "amr_heatmap")),
        ("virulence_heatmap", lambda: plot_heatmap(vf_matrix, "Virulence Factor Presence/Absence Heatmap", outdir, "virulence_heatmap", "PuBu")),
        ("combined_amr_virulence_heatmap", lambda: plot_heatmap(
            pd.concat([amr_matrix.add_prefix("AMR:"), vf_matrix.add_prefix("VF:")], axis=1).fillna(0).astype(int),
            "Combined AMR + Virulence Heatmap", outdir, "combined_amr_virulence_heatmap", "RdYlBu_r")),
        ("amr_frequency", lambda: plot_frequency(amr_matrix, "AMR Gene Frequency Across Samples", outdir, "amr_frequency")),
        ("virulence_frequency", lambda: plot_frequency(vf_matrix, "Virulence Factor Frequency Across Samples", outdir, "virulence_frequency")),
        ("amr_cooccurrence", lambda: plot_cooccurrence(amr_matrix, "AMR Gene Co-occurrence", outdir, "amr_cooccurrence")),
        ("virulence_cooccurrence", lambda: plot_cooccurrence(vf_matrix, "Virulence Factor Co-occurrence", outdir, "virulence_cooccurrence")),
        ("mlst_distribution", lambda: plot_mlst_distribution(mlst_df, outdir)),
        ("resistance_gene_counts", lambda: plot_gene_count_bars(amr_matrix, "Resistance Gene Counts per Sample", outdir, "resistance_gene_counts")),
        ("virulence_gene_counts", lambda: plot_gene_count_bars(vf_matrix, "Virulence Gene Counts per Sample", outdir, "virulence_gene_counts")),
        ("sample_summary_dashboard", lambda: plot_dashboard(amr_matrix, vf_matrix, mlst_df, outdir)),
    ]

    for name, fn in figs:
        try:
            if fn():
                generated.append(name)
            else:
                skipped.append(name)
        except Exception as e:
            write_skip_report(outdir, name, f"Unexpected error: {e}")
            skipped.append(name)

    st_map = {}
    if not mlst_df.empty:
        st_map = dict(zip(mlst_df["sample"], mlst_df["ST"]))

    amr_ann = {s: ",".join(amr_matrix.columns[amr_matrix.loc[s] == 1]) or "None"
               for s in amr_matrix.index} if not amr_matrix.empty else {}
    vf_ann = {s: ",".join(vf_matrix.columns[vf_matrix.loc[s] == 1]) or "None"
              for s in vf_matrix.index} if not vf_matrix.empty else {}

    tree_figs = [
        ("tree_with_st_labels", st_map, "Phylogenetic Tree with ST Labels", None),
        ("tree_with_amr_annotation", amr_ann, "Phylogenetic Tree with AMR Annotation", "amr"),
        ("tree_with_virulence_annotation", vf_ann, "Phylogenetic Tree with Virulence Annotation", "vf"),
    ]
    for basename, ann, title, _ in tree_figs:
        if ann:
            if render_tree(tree_path, ann, outdir, basename, title):
                generated.append(basename)
            else:
                skipped.append(basename)
        else:
            write_skip_report(outdir, basename, "No annotation data available for tree labeling.")
            skipped.append(basename)

    summary_path = reports_dir / "visualization_summary.txt"
    summary_path.write_text(
        f"GBS-Genomics-Pipeline Visualization Report\n"
        f"Generated: {datetime.now().isoformat()}\n\n"
        f"Figures generated ({len(generated)}):\n" +
        "\n".join(f"  - {g}" for g in generated) + "\n\n" +
        f"Figures skipped ({len(skipped)}):\n" +
        "\n".join(f"  - {s} (see Reports/{s}_SKIPPED.txt)" for s in skipped) + "\n"
    )
    log(f"Visualization complete: {len(generated)} generated, {len(skipped)} skipped.")
    log(f"Summary report: {summary_path}")


if __name__ == "__main__":
    main()
