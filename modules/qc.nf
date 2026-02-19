process QC {
    tag "$sample"
    label 'med_compute'

    publishDir "${params.outdir}/qc",   mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'step_log.txt'

    input:
    tuple val(sample), path(reads)

    output:
    // This 'emit' name must match what is used in pipeline.nf
    tuple val(sample), path("${sample}_trimmed_{1,2}.fastq.gz"), emit: trimmed
    path "step_log.txt",                                       emit: log

    script:
    """
    timestamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$timestamp] Step: Quality Control for sample ${sample}" > step_log.txt

    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample}_trimmed_1.fastq.gz \
        -O ${sample}_trimmed_2.fastq.gz

    num_reads_after=\$(( \$(zcat ${sample}_trimmed_1.fastq.gz | wc -l) / 4 ))
    echo "[\$timestamp] Reads after trimming: \$num_reads_after" >> step_log.txt

    fastqc ${sample}_trimmed_1.fastq.gz ${sample}_trimmed_2.fastq.gz
    """
}
