process DOWNLOAD_BLAST_DB {
    label 'low_compute'
    tag "download_blast_db"

    conda "${launchDir}/envs/blast_db.yml"
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir params.blast_db_dir, mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'blast_setup_log.txt'
    
    output:
    path "blast_db", emit: blast_db
    path "blast_setup_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    timestamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$timestamp] Starting BLAST DB Setup." > blast_setup_log.txt

    mkdir -p blast_db

    LOCAL_REF_DIR="${projectDir}/references/streptococcus_genus_db/blast_db"

    if [ -f "\$LOCAL_REF_DIR/strep_ref.fna" ]; then
        echo "[\$timestamp] Local data found. Copying..." >> blast_setup_log.txt
        cp "\$LOCAL_REF_DIR/strep_ref.fna" blast_db/strep_ref.fna
    else
        echo "[\$timestamp] Downloading Reference Genome via CURL..." >> blast_setup_log.txt
        
        # Direct URL to Streptococcus agalactiae NEM316 (Complete Genome)
        # This bypasses the 'datasets' CLI tool version issues entirely
        curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/265/GCF_000007265.1_ASM726v1/GCF_000007265.1_ASM726v1_genomic.fna.gz" -o strep.fna.gz

        if [ ! -f strep.fna.gz ]; then
            echo "[\$timestamp] ERROR: Download failed." >> blast_setup_log.txt
            exit 1
        fi

        gunzip -c strep.fna.gz > blast_db/strep_ref.fna
        rm strep.fna.gz
    fi

    echo "[\$timestamp] Formatting BLAST database..." >> blast_setup_log.txt
    makeblastdb -in blast_db/strep_ref.fna -dbtype nucl -out blast_db/strep_db
    """
}
