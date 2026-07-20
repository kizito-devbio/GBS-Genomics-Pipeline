/*
 * Quality Control — raw read trimming and FastQC assessment.
 *
 * Input:  paired-end FASTQ (sample_1.fastq, sample_2.fastq)
 * Output: trimmed gzipped FASTQ, FastQC reports, step log
 */

process QC {
    tag "$sample"
    label 'med_compute'

    publishDir "${params.outdir}/QC/${sample}", mode: 'copy', pattern: "*.{fastq.gz,html,zip}"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_qc_log.txt"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_trimmed_{1,2}.fastq.gz"), emit: trimmed
    path "${sample}_qc_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_qc_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Quality Control (QC)"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input R1:    ${reads[0]}"
        echo "Input R2:    ${reads[1]}"
        echo "Output dir:  ${params.outdir}/QC/${sample}"
        echo "Software:    fastp \$(fastp --version 2>&1 | head -1), FastQC \$(fastqc --version 2>&1 | head -1)"
        echo "Parameters:  default fastp quality trimming"
        echo "----------------------------------------"
    } > "\$LOG"

    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample}_trimmed_1.fastq.gz \\
        -O ${sample}_trimmed_2.fastq.gz \\
        --thread ${task.cpus} >> "\$LOG" 2>&1

    READS_AFTER=\$(( \$(zcat ${sample}_trimmed_1.fastq.gz | wc -l) / 4 ))
    echo "Reads after trimming: \$READS_AFTER" >> "\$LOG"

    fastqc ${sample}_trimmed_1.fastq.gz ${sample}_trimmed_2.fastq.gz >> "\$LOG" 2>&1

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
