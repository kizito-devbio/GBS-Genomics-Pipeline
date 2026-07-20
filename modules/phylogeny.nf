/*
 * Phylogenetic Tree Construction — maximum-likelihood tree via IQ-TREE.
 *
 * Input:  core gene alignment (FASTA/ALN)
 * Output: Newick tree file, log
 */

process PHYLOGENY {
    tag "phylogeny"
    label 'high_compute'

    conda 'bioconda::iqtree=2.2.0'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/Phylogeny", mode: 'copy'
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "phylogeny.log"

    input:
    path aligned

    output:
    path "gbs_phylogeny_tree.nwk", emit: tree
    path "phylogeny.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="phylogeny.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Phylogenetic Tree Construction"
        echo "Started:     \$TS"
        echo "Input:       ${aligned}"
        echo "Output dir:  ${params.outdir}/Phylogeny"
        echo "Software:    IQ-TREE \$(iqtree2 --version 2>&1 | head -1)"
        echo "Parameters:  model=MFP, bootstrap=1000"
        echo "----------------------------------------"
    } > "\$LOG"

    if [ ! -s "${aligned}" ] || grep -q "NO_ALIGNMENT" "${aligned}" 2>/dev/null; then
        echo "WARNING: No valid alignment — skipping tree construction." >> "\$LOG"
        echo "()" > gbs_phylogeny_tree.nwk
        ELAPSED=\$(( \$(date +%s) - START ))
        echo "Status:      SKIPPED (insufficient samples)" >> "\$LOG"
        echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
        echo "========================================" >> "\$LOG"
        exit 0
    fi

    iqtree2 -s "${aligned}" -m MFP -bb 1000 -nt AUTO --redo --prefix phylogeny_out >> "\$LOG" 2>&1

    if [ -f "phylogeny_out.treefile" ]; then
        cp "phylogeny_out.treefile" "gbs_phylogeny_tree.nwk"
        NUM_TIPS=\$(grep -oE '[A-Za-z0-9_.-]+:[0-9.eE+-]+' gbs_phylogeny_tree.nwk | wc -l || echo 0)
        echo "Tree tips:   \$NUM_TIPS" >> "\$LOG"
    else
        echo "ERROR: IQ-TREE did not produce a treefile." >> "\$LOG"
        exit 1
    fi

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
