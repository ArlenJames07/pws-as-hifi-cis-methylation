#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGN_HIFI } from './modules/local/alignment'
include { CALL_SMALL_VARIANTS; FILTER_SMALL_VARIANTS; DISCOVER_SV; CALL_SV } from './modules/local/variants'
include { PHASE_VARIANTS } from './modules/local/phasing'
include { CALL_CNV } from './modules/local/cnv'
include { CALL_METHYLATION } from './modules/local/methylation'
include { MAKE_FIGURE_1; MAKE_FIGURE_2; MAKE_FIGURE_3; MAKE_FIGURE_4; MAKE_FIGURE_5 } from './modules/local/figures'

def helpMessage() {
    return """
PWS/AS PacBio HiFi analysis pipeline

Required:
  --input             CSV sample sheet (columns: sample,bam)
  --reference         T2T-CHM13v2.0 FASTA
  --reference_fai     FASTA index (defaults to <reference>.fai)
  --tandem_repeats    Tandem-repeat BED for pbsv
  --cnv_exclude       Excluded-regions BED for HiFiCNV
  --cpg_model         pb-CpG-tools pileup model (.tflite)

Common options:
  --outdir            Output directory [results]
  --run_figures       Run Figures 1-3 after the core workflow [false]
  -profile conda      Create pinned Conda environments per process
  -resume             Reuse completed work

Example:
  nextflow run main.nf -profile conda \\
    --input assets/samplesheet.csv \\
    --reference /data/chm13v2.0.fa \\
    --tandem_repeats /data/chm13v2.0.trf.bed \\
    --cnv_exclude /data/cnv.exclude.bed \\
    --cpg_model /opt/pb-CpG-tools/models/pileup_calling_model.v1.tflite

See README.md and nextflow_schema.json for the complete interface.
""".stripIndent()
}

def requiredFile(name, value) {
    if (!value) {
        error "Missing required parameter --${name}. Run with --help for usage."
    }
    return file(value, checkIfExists: true)
}

workflow CORE {
    take:
    samples
    reference
    reference_fai
    tandem_repeats
    cnv_exclude
    cpg_model

    main:
    ALIGN_HIFI(samples, reference, reference_fai)

    CALL_SMALL_VARIANTS(ALIGN_HIFI.out.alignments, reference, reference_fai)
    FILTER_SMALL_VARIANTS(CALL_SMALL_VARIANTS.out.calls)

    DISCOVER_SV(ALIGN_HIFI.out.alignments, tandem_repeats)
    CALL_SV(DISCOVER_SV.out.signatures, reference, reference_fai)

    phasing_inputs = FILTER_SMALL_VARIANTS.out.filtered
        .join(CALL_SV.out.calls, by: 0)
        .map { sample, bam, bai, small_vcf, small_tbi, sv_vcf, sv_tbi ->
            tuple(sample, bam, bai, small_vcf, small_tbi, sv_vcf, sv_tbi)
        }

    PHASE_VARIANTS(phasing_inputs, reference, reference_fai)
    CALL_CNV(PHASE_VARIANTS.out.phased, reference, reference_fai, cnv_exclude)
    CALL_METHYLATION(PHASE_VARIANTS.out.phased, cpg_model)

    emit:
    alignments  = ALIGN_HIFI.out.alignments
    small_calls = FILTER_SMALL_VARIANTS.out.filtered
    sv_calls    = CALL_SV.out.calls
    phased      = PHASE_VARIANTS.out.phased
    cnv         = CALL_CNV.out.cnv
    methylation = CALL_METHYLATION.out.methylation
}

workflow FIGURES {
    take:
    phased
    cnv
    methylation
    gtf
    metadata
    icr_bed
    segdup_bed
    imprintome_bed
    repeats_bed

    main:
    phased_files = phased
        .map { _sample, bam, bai, small_vcf, small_tbi, sv_vcf, sv_tbi, stats, blocks, summary ->
            [bam, bai, small_vcf, small_tbi, sv_vcf, sv_tbi, stats, blocks, summary]
        }
        .collect()
    cnv_files = cnv.map { _sample, files -> files }.collect()
    methylation_files = methylation.map { _sample, files -> files }.collect()

    MAKE_FIGURE_1(phased_files, cnv_files, methylation_files, gtf, metadata)
    MAKE_FIGURE_2(cnv_files, methylation_files, gtf, metadata, icr_bed, segdup_bed)
    MAKE_FIGURE_3(methylation_files, MAKE_FIGURE_1.out.results, gtf, segdup_bed, imprintome_bed, icr_bed, repeats_bed)

    emit:
    figure_1 = MAKE_FIGURE_1.out.results
    figure_2 = MAKE_FIGURE_2.out.results
    figure_3 = MAKE_FIGURE_3.out.results
}

workflow {
    if (params.help) {
        log.info helpMessage()
        return
    }

    input_file = requiredFile('input', params.input)
    reference = requiredFile('reference', params.reference)
    reference_fai = requiredFile('reference_fai', params.reference_fai ?: "${params.reference}.fai")
    tandem_repeats = requiredFile('tandem_repeats', params.tandem_repeats)
    cnv_exclude = requiredFile('cnv_exclude', params.cnv_exclude)
    cpg_model = requiredFile('cpg_model', params.cpg_model)

    parsed_samples = channel
        .fromPath(input_file, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.bam) {
                error "Every sample-sheet row needs non-empty 'sample' and 'bam' values."
            }
            if (!(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error "Invalid sample ID '${row.sample}'. Use letters, numbers, dots, underscores, or hyphens."
            }
            def bam_path = row.bam.startsWith('/') ? row.bam : "${input_file.parent}/${row.bam}"
            tuple(row.sample, file(bam_path, checkIfExists: true))
        }

    samples = parsed_samples
        .groupTuple(by: 0)
        .map { sample, bams ->
            if (bams.size() != 1) {
                error "Sample ID '${sample}' occurs ${bams.size()} times; sample IDs must be unique."
            }
            tuple(sample, bams.first())
        }

    CORE(samples, reference, reference_fai, tandem_repeats, cnv_exclude, cpg_model)

    if (params.run_figures) {
        gtf = requiredFile('gtf', params.gtf)
        metadata = requiredFile('metadata', params.metadata)
        icr_bed = requiredFile('icr_bed', params.icr_bed)
        segdup_bed = requiredFile('segdup_bed', params.segdup_bed)
        imprintome_bed = requiredFile('imprintome_bed', params.imprintome_bed)
        repeats_bed = requiredFile('repeats_bed', params.repeats_bed)

        FIGURES(
            CORE.out.phased,
            CORE.out.cnv,
            CORE.out.methylation,
            gtf,
            metadata,
            icr_bed,
            segdup_bed,
            imprintome_bed,
            repeats_bed
        )
    }

    if (params.figure4_tables) {
        MAKE_FIGURE_4(requiredFile('figure4_tables', params.figure4_tables))
    }
    if (params.figure5_tables) {
        MAKE_FIGURE_5(
            requiredFile('figure5_tables', params.figure5_tables),
            reference,
            requiredFile('gtf', params.gtf)
        )
    }
}
