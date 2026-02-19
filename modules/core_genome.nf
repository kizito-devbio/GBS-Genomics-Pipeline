process CORE_GENOME {
    tag { "core_genome" }
    label 'high_compute'

    // FIXED: Removed strict versioning on gclib and pinned python to 3.9 for stability
    conda 'conda-forge::python=3.9 bioconda::panaroo=1.5.0 bioconda::gclib'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/core_genome", mode: 'copy'

    input:
    path gff_files

    output:
    path "core_alignment.aln", emit: alignment
    path "core_genome.log",    emit: log
    path "panaroo_results",    emit: all_results

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="core_genome.log"
    TIMESTAMP=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$TIMESTAMP] Core Genome Extraction (Panaroo) started" > "\$LOG"

    # Count GFFs
    num_gffs=\$(ls *.gff 2>/dev/null | wc -l || echo 0)
    echo "[\$TIMESTAMP] Number of GFFs found: \$num_gffs" >> "\$LOG"

    if [ "\$num_gffs" -lt 2 ]; then
        echo "[\$TIMESTAMP] ERROR: Need 2+ GFFs. Found: \$num_gffs" >> "\$LOG"
        exit 1
    fi

    echo "[\$TIMESTAMP] Running Panaroo on all GFF files..." >> "\$LOG"

    # FIXED: Corrected backslash syntax for the panaroo command
    panaroo \
        -i *.gff \
        -o panaroo_results \
        --clean-mode strict \
        -a core \
        --threads ${task.cpus} >> "\$LOG" 2>&1

    # Copy the alignment to the top level for the PHYLOGENY process
    if [ -f "panaroo_results/core_gene_alignment.aln" ]; then
        cp "panaroo_results/core_gene_alignment.aln" "core_alignment.aln"
        echo "[\$TIMESTAMP] Alignment produced successfully." >> "\$LOG"
    else
        echo "[\$TIMESTAMP] ERROR: Alignment not produced" >> "\$LOG"
        exit 1
    fi
    """
}	

