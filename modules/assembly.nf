process ASSEMBLY {
    tag "$sample"
    label 'high_compute'

    publishDir "${params.outdir}/assembly", mode: 'copy', pattern: "${sample}_assembly.fasta"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "${sample}_assembly_log.txt"

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
    TIMESTAMP=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "[\$TIMESTAMP] Starting Genome Assembly for sample: ${sample}" > "\$LOG"

    # Running SPAdes
    # --phred-offset 33 is standard for Illumina data
    spades.py --threads ${task.cpus} --phred-offset 33 \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o spades_output >> "\$LOG" 2>&1

    # Check if assembly was successful
    if [ -f spades_output/contigs.fasta ]; then
        mv spades_output/contigs.fasta ${sample}_assembly.fasta
        
        # Calculate assembly stats for the log
        num_contigs=\$(grep -c '^>' ${sample}_assembly.fasta)
        total_bp=\$(grep -v '^>' ${sample}_assembly.fasta | tr -d '\\n' | wc -c)
        
        echo "[\$(date +'%Y-%m-%d %H:%M:%S')] Assembly successful." >> "\$LOG"
        echo "Total contigs: \$num_contigs" >> "\$LOG"
        echo "Total assembly size: \$total_bp bp" >> "\$LOG"
    else
        echo "ERROR: SPAdes failed to produce contigs.fasta. Check the log for details." >> "\$LOG"
        exit 1
    fi

    # Optional: cleanup the large spades_output directory to save space
    rm -rf spades_output
    """
}
