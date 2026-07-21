/*
 * Integration & Visualization — final workflow stage.
 *
 * Consumes real outputs from AMR, Virulence, MLST, and Phylogeny stages.
 * Generates publication-quality figures via bin/generate_figures.py.
 * Headless Docker-safe ETE3 rendering.
 */

process INTEGRATION_VISUALIZATION {

    tag "visualization"

    label 'low_compute'

    container 'kizitodevbio/strepto-pipeline:latest'


    /*
     * Publish generated visualization outputs
     */

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.png"

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.svg"

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.pdf"

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.tsv"

    publishDir "${params.outdir}/Reports",
        mode: 'copy',
        pattern: "Reports/**"


    input:

    path tree
    path amr_files
    path vf_files
    path mlst_files


    output:

    path "*.png", emit: figures_png
    path "*.svg", emit: figures_svg
    path "*.pdf", emit: figures_pdf
    path "*.tsv", emit: tables
    path "Reports/**", emit: reports
    path "integration_viz.log", emit: log


    script:

    """

    #!/usr/bin/env bash

    set -euo pipefail


    LOG="integration_viz.log"

    START=\$(date +%s)


    # -----------------------------------
    # Headless visualization configuration
    # -----------------------------------

    export QT_QPA_PLATFORM=offscreen
    export MPLBACKEND=Agg
    export DISPLAY=:99


    {

        echo "========================================"
        echo "Stage: Integration & Visualization"
        echo "Started: \$(date)"
        echo "Mode: Docker headless rendering"
        echo "========================================"

    } > \$LOG



    mkdir -p amr_inputs
    mkdir -p vf_inputs
    mkdir -p mlst_inputs
    mkdir -p Reports



    echo "Collecting inputs..." >> \$LOG



    # -----------------------------------
    # AMR inputs
    # -----------------------------------

    for f in ${amr_files}; do

        cp "\$f" amr_inputs/

        echo "AMR: \$(basename "\$f")" >> \$LOG

    done



    # -----------------------------------
    # Virulence inputs
    # -----------------------------------

    for f in ${vf_files}; do

        cp "\$f" vf_inputs/

        echo "VF: \$(basename "\$f")" >> \$LOG

    done



    # -----------------------------------
    # MLST inputs
    # -----------------------------------

    for f in ${mlst_files}; do

        cp "\$f" mlst_inputs/

        echo "MLST: \$(basename "\$f")" >> \$LOG

    done



    # -----------------------------------
    # Phylogenetic tree
    # -----------------------------------

    cp ${tree} tree_input.nwk

    echo "Tree: ${tree}" >> \$LOG



    # -----------------------------------
    # Check visualization environment
    # -----------------------------------

    echo "Checking Python visualization stack..." >> \$LOG


    python3 - <<'PY' >> \$LOG 2>&1

import matplotlib
matplotlib.use("Agg")

from ete3 import Tree

print("ETE3 OK")

tree = Tree("tree_input.nwk")

print("Tree loaded:", len(tree))

PY



    echo "Running figure generation..." >> \$LOG



    python3 ${projectDir}/bin/generate_figures.py \
        --amr-dir amr_inputs \
        --vf-dir vf_inputs \
        --mlst-dir mlst_inputs \
        --tree tree_input.nwk \
        --outdir . >> \$LOG 2>&1



    # -----------------------------------
    # Summary report
    # -----------------------------------

    cat > Reports/visualization_summary.txt <<EOF

GBS-Genomics-Pipeline
Integration Visualization Summary
=================================

Status: SUCCESS

Generated components:
- Antimicrobial resistance profiles
- Virulence factor profiles
- MLST distribution
- Core genome phylogeny
- AMR/Virulence integrated visualizations

Generated files:
- PNG figures
- SVG figures
- PDF figures
- TSV matrices

EOF



    ELAPSED=\$(( \$(date +%s)-START ))


    echo "----------------------------------------" >> \$LOG
    echo "Completed: \$(date)" >> \$LOG
    echo "Elapsed: \${ELAPSED}s" >> \$LOG
    echo "Status: SUCCESS" >> \$LOG
    echo "========================================" >> \$LOG


    """

}

