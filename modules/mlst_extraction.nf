/*
 * MLST Extraction — Multi-Locus Sequence Typing using mlst (PubMLST schemes).
 *
 * Input:  genome FASTA
 * Output: MLST assignment TSV, log
 */

process MLST_EXTRACTION {
    tag { sample }
    label 'low_compute'

    conda 'bioconda::mlst=2.23.0'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/MLST", mode: 'copy', pattern: "${sample}_mlst.tsv"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_mlst.log"

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
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       MLST Extraction"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${fasta_file}"
        echo "Output dir:  ${params.outdir}/MLST"
        echo "Software:    mlst \$(mlst --version 2>&1 | head -1)"
        echo "Parameters:  scheme=${params.mlst_scheme}"
        echo "----------------------------------------"
    } > "\$LOG"

    mlst --scheme ${params.mlst_scheme} "${fasta_file}" > "${sample}_mlst.tsv" 2>> "\$LOG" || {
        echo "WARNING: MLST assignment failed for ${sample}" >> "\$LOG"
        echo "${sample}\\tunknown\\t-" > "${sample}_mlst.tsv"
    }

    if [ ! -s "${sample}_mlst.tsv" ]; then
        echo "${sample}\\tunknown\\t-" > "${sample}_mlst.tsv"
    fi

    echo "MLST result: \$(cat ${sample}_mlst.tsv)" >> "\$LOG"

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
