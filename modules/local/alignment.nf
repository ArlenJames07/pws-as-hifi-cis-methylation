process ALIGN_HIFI {
    tag "${sample}"
    label 'process_high'
    conda "${moduleDir}/../../envs/pbmm2.yml"

    publishDir "${params.outdir}/01_alignment", mode: params.publish_mode

    input:
    tuple val(sample), path(input_bam)
    path reference
    path reference_fai

    output:
    tuple val(sample), path("${sample}.aligned.bam"), path("${sample}.aligned.bam.bai"), emit: alignments
    path "${sample}.alignment.versions.yml", emit: versions

    script:
    """
    ${params.pbmm2_bin} align \
        ${reference} \
        ${input_bam} \
        ${sample}.aligned.bam \
        --preset HIFI \
        --sort \
        --log-level INFO \
        --num-threads ${task.cpus}

    ${params.samtools_bin} index -@ ${task.cpus} ${sample}.aligned.bam

    cat > ${sample}.alignment.versions.yml <<-END_VERSIONS
    "${task.process}":
      pbmm2: "\$(${params.pbmm2_bin} --version 2>&1 | head -n 1)"
      samtools: "\$(${params.samtools_bin} --version | head -n 1)"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.aligned.bam ${sample}.aligned.bam.bai
    echo '${task.process}: {stub: true}' > ${sample}.alignment.versions.yml
    """
}
