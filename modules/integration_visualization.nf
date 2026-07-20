/*
 * Integration & Visualization — final workflow stage.
 *
 * Consumes real outputs from AMR, Virulence, MLST, and Phylogeny stages.
 * Generates publication-quality figures via bin/generate_figures.py.
 * Writes skip reports when data are insufficient — never fabricates results.
 */

process INTEGRATION_VISUALIZATION {
    tag "visualization"
    label 'low_compute'

    conda 'bioconda::ete3=3.1.3 conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/Figures", mode: 'copy', pattern: "*.{png,svg,pdf,tsv}"
    publishDir "${params.outdir}/Reports", mode: 'copy', pattern: "Reports/**"

    input:
    path tree
    path amr_files
    path vf_files
    path mlst_files

    output:
    path "Reports/visualization_summary.txt", emit: summary
    path "integration_viz.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="integration_viz.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Integration & Visualization"
        echo "Started:     \$TS"
        echo "Output dir:  ${params.outdir}/Figures, ${params.outdir}/Reports"
        echo "Software:    generate_figures.py (Python)"
        echo "----------------------------------------"
    } > "\$LOG"

    mkdir -p amr_inputs vf_inputs mlst_inputs Reports

    echo "Organizing input files..." >> "\$LOG"

    for f in ${amr_files}; do
        cp "\$f" "amr_inputs/\$(basename \$f)"
        echo "  AMR:   \$(basename \$f)" >> "\$LOG"
    done

    for f in ${vf_files}; do
        cp "\$f" "vf_inputs/\$(basename \$f)"
        echo "  VF:    \$(basename \$f)" >> "\$LOG"
    done

    for f in ${mlst_files}; do
        cp "\$f" "mlst_inputs/\$(basename \$f)"
        echo "  MLST:  \$(basename \$f)" >> "\$LOG"
    done

    cp "${tree}" tree_input.nwk
    echo "  Tree:  \$(basename ${tree})" >> "\$LOG"

    python3 ${projectDir}/bin/generate_figures.py \\
        --amr-dir amr_inputs \\
        --vf-dir vf_inputs \\
        --mlst-dir mlst_inputs \\
        --tree tree_input.nwk \\
        --outdir . >> "\$LOG" 2>&1

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
