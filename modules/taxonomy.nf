process TAXONOMY {
    tag "$sample"
    label 'low_compute'

    conda 'bioconda::blast'
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*_taxonomy.log'

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

    DATE=$(date +"%Y-%m-%d %H:%M:%S")
    LOG="!{sample}_taxonomy.log"

    echo "[$DATE] Processing sample: !{sample}" > "$LOG"

    # Define path to BLAST DB index
    DB_INDEX="!{blast_db_folder}/strep_db"

    if [ ! -f "${DB_INDEX}.nhr" ]; then
        echo "[$DATE] ERROR: BLAST index 'strep_db' not found in !{blast_db_folder}" >> "$LOG"
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
        # Parse species name from the 13th column (stitle) or use 2nd (sseqid)
        SPECIES_NAME=$(sort -k12,12rn "!{sample}_blast_raw.tsv" | awk -F'\t' 'NR==1 { 
            if ($13 == "" || $13 == "N/A") 
                print $2; 
            else 
                print $13 
        }')
        echo "$SPECIES_NAME" > "!{sample}_taxonomy.txt"
        echo "[$DATE] Taxonomy identified for !{sample}: $SPECIES_NAME" >> "$LOG"
    else
        echo "Unknown Streptococcus" > "!{sample}_taxonomy.txt"
        echo "[$DATE] WARNING: No hits found for !{sample} above 95% identity" >> "$LOG"
    fi
    '''
}
