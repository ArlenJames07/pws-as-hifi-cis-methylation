#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CALL_METHYLATION } from './methylation_process'

/*
 * Standalone methylation entry point.
 *
 * This workflow deliberately reads already-phased BAMs from params.phased_dir.
 * It does not include or schedule alignment, variant calling, phasing, or CNV.
 */
workflow {
    if (!params.phased_dir) {
        error "Missing required parameter --phased_dir"
    }
    if (!params.cpg_model) {
        error "Missing required parameter --cpg_model"
    }
    if (!params.methylation_samples) {
        error "Missing required parameter --methylation_samples (comma-separated sample IDs)"
    }

    phased_dir = file(params.phased_dir, checkIfExists: true)
    cpg_model = file(params.cpg_model, checkIfExists: true)

    if (!java.nio.file.Files.isDirectory(phased_dir)) {
        error "--phased_dir is not a directory: ${phased_dir}"
    }

    sample_ids = params.methylation_samples
        .toString()
        .split(',')
        .collect { it.trim() }
        .findAll { it }

    if (sample_ids.size() != sample_ids.toSet().size()) {
        error "--methylation_samples contains duplicate sample IDs"
    }

    methylation_rows = sample_ids.collect { sample ->
        if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "Invalid sample ID '${sample}'"
        }

        def matches = []
        java.nio.file.Files.newDirectoryStream(phased_dir, "*_${sample}.bam").withCloseable { paths ->
            paths.each { matches << it }
        }

        if (matches.size() != 1) {
            error "Expected exactly one '*_${sample}.bam' in ${phased_dir}; found ${matches.size()}"
        }

        def bam = matches.first()
        def bai = bam.resolveSibling("${bam.fileName}.bai")
        if (!java.nio.file.Files.exists(bai)) {
            error "Missing BAM index for ${bam}: expected ${bai}"
        }

        tuple(sample, bam, bai)
    }

    methylation_inputs = channel.fromList(methylation_rows)
    CALL_METHYLATION(methylation_inputs, cpg_model)
}
