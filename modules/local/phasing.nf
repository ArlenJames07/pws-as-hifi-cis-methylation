process PHASE_VARIANTS {
    tag "${sample}"
    label 'process_high'
    conda "${moduleDir}/../../envs/hiphase.yml"

    publishDir "${params.outdir}/04_phasing", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai), path(small_vcf), path(small_tbi), path(sv_vcf), path(sv_tbi)
    path reference
    path reference_fai

    output:
    tuple val(sample),
        path("${sample}.phased.bam"), path("${sample}.phased.bam.bai"),
        path("${sample}.small.phased.vcf.gz"), path("${sample}.small.phased.vcf.gz.tbi"),
        path("${sample}.sv.phased.vcf.gz"), path("${sample}.sv.phased.vcf.gz.tbi"),
        path("${sample}.stats.csv"), path("${sample}.blocks.tsv"), path("${sample}.summary.tsv"),
        emit: phased
    path "${sample}.phasing.versions.yml", emit: versions

    script:
    """
    ${params.hiphase_bin} \
        --threads ${task.cpus} \
        --reference ${reference} \
        --bam ${bam} \
        --output-bam ${sample}.phased.bam \
        --vcf ${small_vcf} \
        --output-vcf ${sample}.small.phased.vcf.gz \
        --vcf ${sv_vcf} \
        --output-vcf ${sample}.sv.phased.vcf.gz \
        --stats-file ${sample}.stats.csv \
        --blocks-file ${sample}.blocks.tsv \
        --summary-file ${sample}.summary.tsv

    ${params.samtools_bin} index -@ ${task.cpus} ${sample}.phased.bam
    ${params.bcftools_bin} index --force --threads ${task.cpus} --tbi ${sample}.small.phased.vcf.gz
    ${params.bcftools_bin} index --force --threads ${task.cpus} --tbi ${sample}.sv.phased.vcf.gz

    cat > ${sample}.phasing.versions.yml <<-END_VERSIONS
    "${task.process}":
      hiphase: "\$(${params.hiphase_bin} --version 2>&1 | head -n 1)"
      samtools: "\$(${params.samtools_bin} --version | head -n 1)"
      bcftools: "\$(${params.bcftools_bin} --version | head -n 1)"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.phased.bam ${sample}.phased.bam.bai
    touch ${sample}.small.phased.vcf.gz ${sample}.small.phased.vcf.gz.tbi
    touch ${sample}.sv.phased.vcf.gz ${sample}.sv.phased.vcf.gz.tbi
    touch ${sample}.stats.csv ${sample}.blocks.tsv ${sample}.summary.tsv
    echo '${task.process}: {stub: true}' > ${sample}.phasing.versions.yml
    """
}
