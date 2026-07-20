/*
 * Virulence Factor Detection — Abricate screening against VFDB.
 *
 * Input:  annotated genome FASTA (Prokka .fna)
 * Output: virulence factor TSV, log
 */

process VIRULENCE_FACTOR {
    tag { sample }
    label 'low_compute'

    conda 'bioconda::abricate=1.0.1'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/Virulence", mode: 'copy', pattern: "${sample}_vf.tsv"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_virulence.log"

    input:
    tuple val(sample), path(fasta_file)

    output:
    tuple val(sample), path("${sample}_vf.tsv"), emit: results
    path "${sample}_virulence.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    export LC_ALL=C
    export LANG=C

    LOG="${sample}_virulence.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Virulence Factor Detection"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${fasta_file}"
        echo "Output dir:  ${params.outdir}/Virulence"
        echo "Software:    Abricate \$(abricate --version 2>&1 | head -1)"
        echo "Parameters:  database=vfdb"
        echo "----------------------------------------"
    } > "\$LOG"

    abricate --db vfdb "${fasta_file}" > "${sample}_vf.tsv" 2>> "\$LOG" || {
        echo "WARNING: Abricate returned non-zero exit code" >> "\$LOG"
        echo "#FILE\\tSEQUENCE\\tSTART\\tEND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%COVERAGE\\t%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT\\tRESISTANCE" > "${sample}_vf.tsv"
    }

    if [ ! -s "${sample}_vf.tsv" ]; then
        echo "# No virulence factors detected" > "${sample}_vf.tsv"
    fi

    NUM_VF=\$(grep -vc "^#" "${sample}_vf.tsv" 2>/dev/null || echo 0)
    echo "VF detected: \$NUM_VF" >> "\$LOG"

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
