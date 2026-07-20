#!/usr/bin/env nextflow
/*
 * GBS-Genomics-Pipeline — Streptococcus agalactiae whole-genome analysis workflow.
 *
 * Input pathways:
 *   --curated_dir  Pre-assembled FASTA files (.fa, .fna, .fasta)
 *   --raw_dir      Paired-end FASTQ files (*_1.fastq, *_2.fastq)
 *
 * Raw-read pathway: QC → Assembly → Quality Assessment (no decontamination)
 * See modules/frozen/ for previously removed human decontamination modules.
 */

nextflow.enable.dsl = 2

// ── Module imports ──────────────────────────────────────────────────
include { DOWNLOAD_BLAST_DB }       from './modules/download_blast_db.nf'
include { QC }                      from './modules/qc.nf'
include { ASSEMBLY }                from './modules/assembly.nf'
include { QUALITY_ASSESS }          from './modules/quality_assess.nf'
include { TAXONOMY }                from './modules/taxonomy.nf'
include { BACKGROUND_SELECTION }    from './modules/background_selection.nf'
include { FUNCTIONAL_ANNOTATION }   from './modules/functional_annotation.nf'
include { VIRULENCE_FACTOR }        from './modules/virulence_factor.nf'
include { MLST_EXTRACTION }         from './modules/mlst_extraction.nf'
include { CORE_GENOME }             from './modules/core_genome.nf'
include { PHYLOGENY }               from './modules/phylogeny.nf'
include { INTEGRATION_VISUALIZATION } from './modules/integration_visualization.nf'

// ── Parameter validation ────────────────────────────────────────────
def printHelp() {
    log.info """
    GBS-Genomics-Pipeline — Parameter Reference
    ==========================================

    Required (one of):
      --raw_dir       <dir>   Paired-end FASTQ input directory
      --curated_dir   <dir>   Pre-assembled FASTA input directory

    Optional:
      --outdir        <dir>   Output directory [default: ./results]
      --mlst_scheme   <str>   MLST scheme [default: sagalactiae]
      --min_n50       <int>   Minimum N50 for assembly filter [default: 10000]
      --max_cpus      <int>   Maximum CPU cores [auto-detected]
      --max_memory    <str>   Maximum memory [default: 6 GB]

    Profiles:
      -profile docker        Run with Docker container
      -profile singularity    Run with Singularity/Apptainer
      -profile conda          Run with Conda environments
      -profile cluster        Submit to SLURM cluster
      -profile test           Reduced resources for testing

    Example:
      nextflow run pipeline.nf -profile docker --curated_dir data/genomes --outdir results
    """
}

def validateParams() {
    if (params.help) {
        printHelp()
        exit 0
    }
    if (!params.curated_dir && !params.raw_dir) {
        error "Provide exactly one input: --curated_dir <path> OR --raw_dir <path>"
    }
    if (params.curated_dir && params.raw_dir) {
        error "Provide only one input pathway: --curated_dir OR --raw_dir, not both"
    }
    if (params.curated_dir && !file(params.curated_dir).exists()) {
        error "curated_dir does not exist: ${params.curated_dir}"
    }
    if (params.raw_dir && !file(params.raw_dir).exists()) {
        error "raw_dir does not exist: ${params.raw_dir}"
    }
}

workflow {
    validateParams()

    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║  GBS-Genomics-Pipeline                                   ║
    ║  Streptococcus agalactiae Whole-Genome Analysis          ║
    ╚══════════════════════════════════════════════════════════╝
    Input:    ${params.curated_dir ? "curated (${params.curated_dir})" : "raw reads (${params.raw_dir})"}
    Output:   ${params.outdir}
    Profile:  ${workflow.profile ?: 'default'}
    """

    // ── 1. Reference database setup ───────────────────────────────
    blast_db_ch = DOWNLOAD_BLAST_DB()

    // ── 2. Input acquisition ──────────────────────────────────────
    if (params.curated_dir) {
        ch_genomes = Channel
            .fromPath("${params.curated_dir}/*.{fa,fna,fasta}")
            .map { f -> [ f.baseName, f ] }
    } else {
        // Raw pathway: QC → Assembly (decontamination removed — see modules/frozen/)
        ch_raw_reads = Channel.fromFilePairs("${params.raw_dir}/*_{1,2}.fastq", flat: true)
        qc_out       = QC(ch_raw_reads)
        assembly_out = ASSEMBLY(qc_out.trimmed)
        qa_out       = QUALITY_ASSESS(assembly_out.assembled)

        ch_genomes = qa_out.metrics
            .filter { meta, fasta, n50_file, total_len, contigs ->
                n50_file.text.trim().toInteger() >= params.min_n50
            }
            .map { meta, fasta, n50, total_len, contigs -> [ meta, fasta ] }
    }

    // ── 3. Taxonomic confirmation ─────────────────────────────────
    ch_tax_out = TAXONOMY(ch_genomes, blast_db_ch.blast_db.collect()).results

    // ── 4. Background genome selection ────────────────────────────
    BACKGROUND_SELECTION(ch_tax_out.map { it[1] }.collect())

    // ── 5. Functional annotation (Prokka + AMR) ───────────────────
    ch_annot_results = FUNCTIONAL_ANNOTATION(ch_genomes).annotation_results

    // ── 6. Virulence factors & MLST ───────────────────────────────
    ch_fna_for_tools = ch_annot_results.map { sample, amr, prokka_dir ->
        [ sample, file("${prokka_dir}/${sample}.fna") ]
    }

    ch_vf_out   = VIRULENCE_FACTOR(ch_fna_for_tools).results
    ch_mlst_out = MLST_EXTRACTION(ch_fna_for_tools).results

    // ── 7. Core genome & phylogeny ────────────────────────────────
    ch_gffs = ch_annot_results
        .map { sample, amr, prokka_dir -> file("${prokka_dir}/${sample}.gff") }
        .collect()

    ch_core = CORE_GENOME(ch_gffs)
    ch_tree = PHYLOGENY(ch_core.alignment)

    // ── 8. Integration & visualization (final stage) ──────────────
    INTEGRATION_VISUALIZATION(
        ch_tree.tree.collect(),
        ch_annot_results.map { s, amr, prokka_dir -> amr }.collect(),
        ch_vf_out.map { s, vf -> vf }.collect(),
        ch_mlst_out.map { s, mlst -> mlst }.collect()
    )
}

workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║  Pipeline finished: ${status}
    ║  Results: ${params.outdir}
    ║  Reports: ${params.outdir}/Reports/
    ║  Figures: ${params.outdir}/Figures/
    ╚══════════════════════════════════════════════════════════╝
    """
}
