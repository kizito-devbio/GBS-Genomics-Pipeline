/*
 * Taxonomic Confirmation — BLASTn against Streptococcus reference database.
 *
 * Input:  genome FASTA, BLAST database directory
 * Output: taxonomy assignment text file, BLAST results, log
 */

process TAXONOMY {
    tag "$sample"
    label 'low_compute'

    conda 'bioconda::blast=2.15.0'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/Annotation", mode: 'copy', pattern: "${sample}_taxonomy.txt"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "${sample}_taxonomy.log"

    input:
    tuple val(sample), path(genome)
    path blast_db_folder

    output:
    tuple val(sample), path(genome), path("${sample}_taxonomy.txt"), emit: results
    path("${sample}_taxonomy.log"), emit: log

    shell:
    '''
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="!{sample}_taxonomy.log"
    START=$(date +%s)
    DATE=$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Taxonomic Confirmation"
        echo "Sample:      !{sample}"
        echo "Started:     $DATE"
        echo "Input:       !{genome}"
        echo "Output dir:  !{params.outdir}/Annotation"
        echo "Software:    BLAST+ $(blastn -version 2>&1 | head -1)"
        echo "Parameters:  perc_identity=95, max_hsps=1"
        echo "----------------------------------------"
    } > "$LOG"

    DB_INDEX="!{blast_db_folder}/strep_db"

    if [ ! -f "${DB_INDEX}.nhr" ]; then
        echo "ERROR: BLAST index not found at ${DB_INDEX}" >> "$LOG"
        exit 1
    fi

    blastn \
        -query "!{genome}" \
        -db "$DB_INDEX" \
        -perc_identity 95 \
        -max_hsps 1 \
        -num_alignments 5 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart ssend evalue bitscore stitle" \
        -out "!{sample}_blast_raw.tsv" >> "$LOG" 2>&1

    if [ -s "!{sample}_blast_raw.tsv" ]; then
        SPECIES_NAME=$(sort -k12,12rn "!{sample}_blast_raw.tsv" | awk -F'\t' 'NR==1 {
            if ($13 == "" || $13 == "N/A") print $2; else print $13
        }')
        echo "$SPECIES_NAME" > "!{sample}_taxonomy.txt"
        echo "Taxonomy:    $SPECIES_NAME" >> "$LOG"
    else
        echo "Unknown Streptococcus" > "!{sample}_taxonomy.txt"
        echo "WARNING: No BLAST hits above 95% identity" >> "$LOG"
    fi

    ELAPSED=$(( $(date +%s) - START ))
    echo "----------------------------------------" >> "$LOG"
    echo "Completed:   $(date +"%Y-%m-%d %H:%M:%S")" >> "$LOG"
    echo "Elapsed:     ${ELAPSED}s" >> "$LOG"
    echo "Status:      SUCCESS" >> "$LOG"
    echo "========================================" >> "$LOG"
    '''
}
