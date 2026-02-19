process QUALITY_ASSESS {
    tag "$sample"
    label 'med_compute'
    publishDir "${params.outdir}/quality_assess", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'step_log.txt'

    input:
    tuple val(sample), path(assembly)

    output:
    // We emit the tuple. Note: we read the values from files created in the script
    tuple val(sample), path(assembly), path("n50.txt"), path("total_len.txt"), path("contigs.txt"), emit: metrics
    path "${sample}_quast_report.tsv", emit: report
    path "step_log.txt",               emit: log

    script:
    """
    timestamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$timestamp] Step: Assembly Quality Assessment for sample ${sample}." > step_log.txt
    
    quast.py $assembly -o quast_out --min-contig 500 --no-plots --no-html || { echo "QUAST failed"; exit 1; }
    cp quast_out/report.tsv ${sample}_quast_report.tsv
    
    # Extract values and save to files so Nextflow can grab them as 'path' outputs
    grep '^N50' ${sample}_quast_report.tsv | cut -f2 > n50.txt
    grep '^Total length' ${sample}_quast_report.tsv | cut -f2 > total_len.txt
    grep '^# contigs' ${sample}_quast_report.tsv | cut -f2 > contigs.txt
    
    echo "[\$timestamp] Assessment complete." >> step_log.txt
    """
}
