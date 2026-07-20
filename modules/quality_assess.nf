/*
 * Assembly Quality Assessment — QUAST metrics for assembled genomes.
 *
 * Input:  assembled FASTA
 * Output: QUAST report TSV, N50/total/contig metrics, step log
 */

process QUALITY_ASSESS {
    tag "$sample"
    label 'med_compute'

    publishDir "${params.outdir}/Assembly", mode: 'copy', pattern: "${sample}_quast_report.tsv"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_qa_log.txt"

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path(assembly), path("n50.txt"), path("total_len.txt"), path("contigs.txt"), emit: metrics
    path "${sample}_quast_report.tsv", emit: report
    path "${sample}_qa_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_qa_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Assembly Quality Assessment"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${assembly}"
        echo "Output dir:  ${params.outdir}/Assembly"
        echo "Software:    QUAST \$(quast.py --version 2>&1 | head -1)"
        echo "Parameters:  min-contig=500"
        echo "----------------------------------------"
    } > "\$LOG"

    quast.py ${assembly} -o quast_out --min-contig 500 --no-plots --no-html >> "\$LOG" 2>&1
    cp quast_out/report.tsv ${sample}_quast_report.tsv

    grep '^N50' ${sample}_quast_report.tsv | cut -f2 > n50.txt
    grep '^Total length' ${sample}_quast_report.tsv | cut -f2 > total_len.txt
    grep '^# contigs' ${sample}_quast_report.tsv | cut -f2 > contigs.txt

    echo "N50:         \$(cat n50.txt)" >> "\$LOG"
    echo "Total len:   \$(cat total_len.txt)" >> "\$LOG"
    echo "Contigs:     \$(cat contigs.txt)" >> "\$LOG"

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
