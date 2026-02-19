process BACKGROUND_SELECTION {
    tag "background_selection"
    label 'med_compute'

    // Fix: wget is in conda-forge, not bioconda
    conda 'conda-forge::wget=1.21.4'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/background", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'step_log.txt'

    input:
    path genomes

    output:
    path "background_genomes/*.fna", emit: background_fasta
    path "step_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    TIMESTAMP=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$TIMESTAMP] Step: Preparing Streptococcus background" > step_log.txt

    mkdir -p background_genomes

    # Copy all input genomes
    echo "[\$TIMESTAMP] Adding sample genomes" >> step_log.txt
    for genome_file in ${genomes}; do
        cp "\$genome_file" background_genomes/\$(basename "\$genome_file")
    done

    # Download reference genome only if not already present
    REF_GENOME="background_genomes/S_agalactiae_ref.fna"
    if [ ! -f "\$REF_GENOME" ]; then
        echo "[\$TIMESTAMP] Downloading Streptococcus reference genome" >> step_log.txt
        wget -q "${params.background_ref_url}" -O "\$REF_GENOME.gz"
        gunzip -f "\$REF_GENOME.gz"
    else
        echo "[\$TIMESTAMP] Reference genome already present" >> step_log.txt
    fi

    echo "[\$TIMESTAMP] Background genome preparation completed" >> step_log.txt
    """
}
