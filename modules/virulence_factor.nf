process VIRULENCE_FACTOR {
    tag { sample }
    label 'low_compute'

    // Open Source Compatibility: Pins version to match functional annotation
    conda 'bioconda::abricate=1.0.1'
    container 'kizitodevbio/strepto-pipeline:latest'

    // Results go into a dedicated virulence folder
    publishDir "${params.outdir}/virulence/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta_file)

    output:
    // Organized output tuple
    tuple val(sample), path("${sample}_vf.tsv"), emit: results
    path "${sample}_virulence.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Standardize locale to avoid grep/sort issues
    export LC_ALL=C
    export LANG=C

    LOG="${sample}_virulence.log"
    TIMESTAMP=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$TIMESTAMP] Virulence factor detection for ${sample}" > "\$LOG"

    # Run ABRicate using the VFDB (Virulence Factor Database)
    echo "[\$TIMESTAMP] Running ABRicate with vfdb..." >> "\$LOG"
    abricate --db vfdb "$fasta_file" > "${sample}_vf.tsv" 2>> "\$LOG" || {
        echo "[\$TIMESTAMP] WARNING: ABRicate failed, creating empty VF file" >> "\$LOG"
        # Ensure header exists even if tool fails
        echo "#FILE\\tSEQUENCE\\tSTART\\tEND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%COVERAGE\\t%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT\\tRESISTANCE" > "${sample}_vf.tsv"
    }

    # Ensure the file is not empty/has a header
    if [ ! -s "${sample}_vf.tsv" ]; then
        echo "[\$TIMESTAMP] WARNING: ABRicate produced empty file" >> "\$LOG"
        touch "${sample}_vf.tsv"
    fi

    # Count detected virulence factors (excluding header)
    NUM_VF=\$(grep -v "FILE" "${sample}_vf.tsv" | wc -l || echo 0)
    echo "[\$TIMESTAMP] Virulence factors detected: \$NUM_VF" >> "\$LOG"
    """
}
