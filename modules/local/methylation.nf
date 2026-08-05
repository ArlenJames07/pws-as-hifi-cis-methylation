process CALL_METHYLATION {
    tag "${sample}"
    label 'process_high'

    publishDir { "${params.outdir}/06_methylation/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai), path(small_vcf), path(small_tbi), path(sv_vcf), path(sv_tbi), path(stats), path(blocks), path(summary)
    path cpg_model

    output:
    tuple val(sample), path("${sample}.cpg.*"), emit: methylation
    path "${sample}.methylation.versions.yml", emit: versions

    script:
    """
    ${params.pbcpgtools_bin} \
        --bam ${bam} \
        --output-prefix ${sample}.cpg \
        --model ${cpg_model} \
        --threads ${task.cpus}

    cat > ${sample}.methylation.versions.yml <<-END_VERSIONS
    "${task.process}":
      pb_cpg_tools: "${params.pbcpgtools_version}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.cpg.combined.bed ${sample}.cpg.hap1.bed ${sample}.cpg.hap2.bed
    echo '${task.process}: {stub: true}' > ${sample}.methylation.versions.yml
    """
}
