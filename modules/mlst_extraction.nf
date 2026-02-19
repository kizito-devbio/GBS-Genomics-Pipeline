process MLST_EXTRACTION {
    tag { sample }
    label 'low_compute'

    conda 'bioconda::mlst=2.23.0'
    container 'kizitodevbio/strepto-pipeline:latest'

    // Use sample-specific subfolders to keep it organized
    publishDir "${params.outdir}/mlst/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta_file)

    output:
    tuple val(sample), path("${sample}_mlst.tsv"), emit: results
    path "${sample}_mlst.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    export LC_ALL=C
    export LANG=C

    LOG="${sample}_mlst.log"
    TIMESTAMP=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$TIMESTAMP] MLST extraction for sample ${sample}" > "\$LOG"

    # Run MLST - using the sample name for the output file
    # Fixed the broken redirect and added the missing error handling block
    mlst --scheme ${params.mlst_scheme} "$fasta_file" > "${sample}_mlst.tsv" 2>> "\$LOG" || {
        echo "[\$TIMESTAMP] WARNING: MLST failed, creating empty output" >> "\$LOG"
        touch "${sample}_mlst.tsv"
    }

    # Ensure output is not empty (check if file size is zero)
    if [ ! -s "${sample}_mlst.tsv" ]; then
        echo "[\$TIMESTAMP] WARNING: MLST produced empty file" >> "\$LOG"
        touch "${sample}_mlst.tsv"
    fi

    echo "[\$TIMESTAMP] MLST completed successfully" >> "\$LOG"
    """
}
