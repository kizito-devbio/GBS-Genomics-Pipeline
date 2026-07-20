/*
 * Genome Assembly — de novo assembly of paired-end reads using SPAdes.
 *
 * Input:  trimmed paired-end FASTQ (from QC or direct input)
 * Output: assembled contigs FASTA, assembly log
 */

process ASSEMBLY {
    tag "$sample"
    label 'high_compute'

    publishDir "${params.outdir}/Assembly", mode: 'copy', pattern: "${sample}_assembly.fasta"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_assembly_log.txt"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_assembly.fasta"), emit: assembled
    path "${sample}_assembly_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_assembly_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Genome Assembly"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input R1:    ${reads[0]}"
        echo "Input R2:    ${reads[1]}"
        echo "Output dir:  ${params.outdir}/Assembly"
        echo "Software:    SPAdes \$(spades.py --version 2>&1 | head -1)"
        echo "Parameters:  threads=${task.cpus}, phred-offset=33"
        echo "----------------------------------------"
    } > "\$LOG"

    spades.py --threads ${task.cpus} --phred-offset 33 \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o spades_output >> "\$LOG" 2>&1

    if [ -f spades_output/contigs.fasta ]; then
        mv spades_output/contigs.fasta ${sample}_assembly.fasta
        NUM_CONTIGS=\$(grep -c '^>' ${sample}_assembly.fasta)
        TOTAL_BP=\$(grep -v '^>' ${sample}_assembly.fasta | tr -d '\\n' | wc -c)
        echo "Contigs:     \$NUM_CONTIGS" >> "\$LOG"
        echo "Total size:  \$TOTAL_BP bp" >> "\$LOG"
    else
        echo "ERROR: SPAdes did not produce contigs.fasta" >> "\$LOG"
        exit 1
    fi

    rm -rf spades_output

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
