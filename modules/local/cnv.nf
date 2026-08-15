process CALL_CNV {
    tag "${sample}"
    label 'process_high'
    conda "${moduleDir}/../../envs/hificnv.yml"

    publishDir { "${params.outdir}/05_cnv/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai), path(small_vcf), path(small_tbi), path(sv_vcf), path(sv_tbi), path(stats), path(blocks), path(summary)
    path reference
    path reference_fai
    path cnv_exclude

    output:
    tuple val(sample), path("${sample}.hificnv.*"), emit: cnv
    path "${sample}.hificnv.versions.yml", emit: versions

    script:
    """
    ${params.hificnv_bin} \
        --bam ${bam} \
        --ref ${reference} \
        --maf ${small_vcf} \
        --exclude ${cnv_exclude} \
        --threads ${task.cpus} \
        --output-prefix ${sample}.hificnv

    cat > ${sample}.hificnv.versions.yml <<-END_VERSIONS
    "${task.process}":
      hificnv: "\$(${params.hificnv_bin} --version 2>&1 | head -n 1)"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.hificnv.log ${sample}.hificnv.cnv.bed
    echo '${task.process}: {stub: true}' > ${sample}.hificnv.versions.yml
    """
}
