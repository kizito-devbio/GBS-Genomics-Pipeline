/*
 * Functional Annotation — Prokka gene prediction and Abricate AMR screening (CARD).
 *
 * Input:  genome FASTA
 * Output: Prokka annotation directory, AMR TSV, annotation log
 */

process FUNCTIONAL_ANNOTATION {
    tag { sample }
    label 'med_compute'

    conda 'bioconda::prokka=1.14.6 bioconda::abricate=1.0.1 conda-forge::sed'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/Annotation/${sample}", mode: 'copy', pattern: "prokka_${sample}/**"
    publishDir "${params.outdir}/AMR", mode: 'copy', pattern: "amr_results/${sample}_amr.tsv"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_annotation.log"

    input:
    tuple val(sample), path(fasta_file)

    output:
    tuple val(sample), path("amr_results/${sample}_amr.tsv"), path("prokka_${sample}"), emit: annotation_results
    path "${sample}_annotation.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_annotation.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")
    PROKKA_OUT="prokka_${sample}"
    AMR_OUT="amr_results"

    {
        echo "========================================"
        echo "Stage:       Functional Annotation"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${fasta_file}"
        echo "Output dir:  ${params.outdir}/Annotation/${sample}, ${params.outdir}/AMR"
        echo "Software:    Prokka \$(prokka --version 2>&1 | head -1), Abricate \$(abricate --version 2>&1 | head -1)"
        echo "Parameters:  kingdom=Bacteria, genus=Streptococcus, AMR db=CARD"
        echo "----------------------------------------"
    } > "\$LOG"

    mkdir -p "\$PROKKA_OUT" "\$AMR_OUT"

    prokka \\
        --outdir "\$PROKKA_OUT" \\
        --prefix "${sample}" \\
        --kingdom Bacteria \\
        --genus Streptococcus \\
        --force \\
        --cpus ${task.cpus} \\
        "${fasta_file}" >> "\$LOG" 2>&1

    PROKKA_FNA="\$PROKKA_OUT/${sample}.fna"
    if [ -f "\$PROKKA_FNA" ]; then
        abricate --db card "\$PROKKA_FNA" > "\$AMR_OUT/${sample}_amr.tsv" 2>> "\$LOG"
        AMR_COUNT=\$(grep -vc "^#" "\$AMR_OUT/${sample}_amr.tsv" || echo 0)
        echo "AMR genes:   \$AMR_COUNT" >> "\$LOG"
    else
        echo "ERROR: Prokka did not produce .fna file" >> "\$LOG"
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
