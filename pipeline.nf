#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ──────────────────────────────────────────────
// MODULES
// ──────────────────────────────────────────────
include { DOWNLOAD_HUMAN_REF } from './modules/download_human_ref.nf'
include { DOWNLOAD_BLAST_DB }  from './modules/download_blast_db.nf'
include { QC }                 from './modules/qc.nf'
include { DECONTAM }           from './modules/decontam.nf'
include { ASSEMBLY }           from './modules/assembly.nf'
include { QUALITY_ASSESS }      from './modules/quality_assess.nf'
include { TAXONOMY }           from './modules/taxonomy.nf'
include { BACKGROUND_SELECTION } from './modules/background_selection.nf'
include { FUNCTIONAL_ANNOTATION } from './modules/functional_annotation.nf'
include { VIRULENCE_FACTOR }      from './modules/virulence_factor.nf'
include { MLST_EXTRACTION }       from './modules/mlst_extraction.nf'
include { CORE_GENOME }           from './modules/core_genome.nf'
include { PHYLOGENY }             from './modules/phylogeny.nf'
include { INTEGRATION_VISUALIZATION } from './modules/integration_visualization.nf'

workflow {

    log.info "Starting Strepto Pipeline"

    // ── 1. Reference Setup ──
    human_ref_ch = DOWNLOAD_HUMAN_REF(params.human_accession)
    blast_db_ch = DOWNLOAD_BLAST_DB()

    // ── 2. Input Acquisition ──
    if (params.curated_dir) {
        ch_genomes = Channel.fromPath("${params.curated_dir}/*.{fa,fna,fasta}")
            .map { file -> [ file.baseName, file ] }
    } else if (params.raw_dir) {
        ch_raw_reads = Channel.fromFilePairs("${params.raw_dir}/*_{1,2}.fastq", flat: true)
        qc_out = QC(ch_raw_reads)
        decontam_out = DECONTAM(qc_out.trimmed, human_ref_ch.human_ref, human_ref_ch.index_files)
        assembly_out = ASSEMBLY(decontam_out.cleaned)
        qa_out = QUALITY_ASSESS(assembly_out.assembled)
        ch_genomes = qa_out.metrics
            .filter { meta, fasta, n50_file -> n50_file.text.trim().toInteger() >= params.min_n50 }
            .map { meta, fasta, n50 -> [ meta, fasta ] }
    } else {
        error "Provide --raw_dir or --curated_dir"
    }

    // ── 3. Taxonomy ──
    ch_tax_out = TAXONOMY(ch_genomes, blast_db_ch.blast_db.collect()).results

    // ── 4. Background Selection ──
    BACKGROUND_SELECTION(ch_tax_out.map { it[1] }.collect())

    // ── 5. Functional Annotation ──
    ch_annot_results = FUNCTIONAL_ANNOTATION(ch_genomes).annotation_results

    // ── 6. Virulence & MLST (FIXED LOGIC) ──
    // This explicit mapping prevents the 'call' error on ArrayLists
    ch_fna_for_tools = ch_annot_results.map { sample, amr, prokka_dir ->
        def fna_path = "${prokka_dir}/${sample}.fna"
        return [ sample, file(fna_path) ]
    }

    // Explicitly pass the channel to the tools
    ch_vf_out = VIRULENCE_FACTOR(ch_fna_for_tools).results
    MLST_EXTRACTION(ch_fna_for_tools)

    // ── 7. Pangenomics (Core Genome) ──
    ch_gffs = ch_annot_results.map { sample, amr, prokka_dir ->
        file("${prokka_dir}/${sample}.gff")
    }.collect()

    CORE_GENOME(ch_gffs)
    
    // ── 8. Results Integration ──
    // Prepare metadata files for the visualization tree
    ch_amr_metadata = ch_annot_results.map { sample, amr, dir -> [ sample, amr ] }
    
    // Join AMR and Virulence results by sample ID
    ch_viz_files = ch_amr_metadata.join(ch_vf_out)
        .map { sample, amr, vf -> [ amr, vf ] }
        .flatten()
        .collect()

    // PHYLOGENY generates the tree
    PHYLOGENY(CORE_GENOME.out.alignment)

    // INTEGRATION pulls the tree and the metadata files
    INTEGRATION_VISUALIZATION(PHYLOGENY.out.tree.collect(), ch_viz_files)
}

workflow.onComplete {
    log.info "Pipeline execution finished"
    log.info "Status: ${ workflow.success ? 'OK' : 'Failed' }"
}
