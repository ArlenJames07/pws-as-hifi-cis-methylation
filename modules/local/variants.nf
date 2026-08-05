process CALL_SMALL_VARIANTS {
    tag "${sample}"
    label 'process_high'

    publishDir "${params.outdir}/02_small_variants/raw", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    path reference
    path reference_fai

    output:
    tuple val(sample), path(bam), path(bai), path("${sample}.deepvariant.vcf.gz"), path("${sample}.deepvariant.vcf.gz.tbi"), emit: calls
    tuple val(sample), path("${sample}.deepvariant.g.vcf.gz"), path("${sample}.deepvariant.g.vcf.gz.tbi"), emit: gvcf
    path "${sample}.deepvariant.versions.yml", emit: versions

    script:
    """
    ${params.deepvariant_bin} \
        --model_type=PACBIO \
        --ref=${reference} \
        --reads=${bam} \
        --output_vcf=${sample}.deepvariant.vcf.gz \
        --output_gvcf=${sample}.deepvariant.g.vcf.gz \
        --num_shards=${task.cpus}

    cat > ${sample}.deepvariant.versions.yml <<-END_VERSIONS
    "${task.process}":
      deepvariant: "${params.deepvariant_version}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.deepvariant.vcf.gz ${sample}.deepvariant.vcf.gz.tbi
    touch ${sample}.deepvariant.g.vcf.gz ${sample}.deepvariant.g.vcf.gz.tbi
    echo '${task.process}: {stub: true}' > ${sample}.deepvariant.versions.yml
    """
}

process FILTER_SMALL_VARIANTS {
    tag "${sample}"
    label 'process_medium'
    conda "${moduleDir}/../../envs/htslib.yml"

    publishDir "${params.outdir}/02_small_variants/filtered", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), path(tbi)

    output:
    tuple val(sample), path(bam), path(bai), path("${sample}.small.pass.vcf.gz"), path("${sample}.small.pass.vcf.gz.tbi"), emit: filtered
    path "${sample}.small_filter.versions.yml", emit: versions

    script:
    """
    bcftools view --threads ${task.cpus} --apply-filters PASS --output-type z \
        --output ${sample}.small.pass.vcf.gz ${vcf}
    bcftools index --threads ${task.cpus} --tbi ${sample}.small.pass.vcf.gz

    cat > ${sample}.small_filter.versions.yml <<-END_VERSIONS
    "${task.process}":
      bcftools: "\$(bcftools --version | head -n 1)"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.small.pass.vcf.gz ${sample}.small.pass.vcf.gz.tbi
    echo '${task.process}: {stub: true}' > ${sample}.small_filter.versions.yml
    """
}

process DISCOVER_SV {
    tag "${sample}"
    label 'process_high'
    conda "${moduleDir}/../../envs/pbsv.yml"

    publishDir "${params.outdir}/03_structural_variants/signatures", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    path tandem_repeats

    output:
    tuple val(sample), path("${sample}.svsig.gz"), emit: signatures
    path "${sample}.pbsv_discover.versions.yml", emit: versions

    script:
    """
    pbsv discover --tandem-repeats ${tandem_repeats} ${bam} ${sample}.svsig.gz

    cat > ${sample}.pbsv_discover.versions.yml <<-END_VERSIONS
    "${task.process}":
      pbsv: "\$(pbsv --version 2>&1 | head -n 1)"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.svsig.gz
    echo '${task.process}: {stub: true}' > ${sample}.pbsv_discover.versions.yml
    """
}

process CALL_SV {
    tag "${sample}"
    label 'process_high'
    conda "${moduleDir}/../../envs/pbsv.yml"

    publishDir "${params.outdir}/03_structural_variants/calls", mode: params.publish_mode

    input:
    tuple val(sample), path(signature)
    path reference
    path reference_fai

    output:
    tuple val(sample), path("${sample}.sv.pass.vcf.gz"), path("${sample}.sv.pass.vcf.gz.tbi"), emit: calls
    path "${sample}.pbsv_call.versions.yml", emit: versions

    script:
    """
    pbsv call --num-threads ${task.cpus} ${reference} ${signature} ${sample}.sv.vcf
    bcftools view --apply-filters PASS \
        --include 'INFO/SVTYPE="BND" || abs(INFO/SVLEN)>=${params.min_sv_length}' \
        --output-type z --output ${sample}.sv.pass.vcf.gz ${sample}.sv.vcf
    bcftools index --tbi ${sample}.sv.pass.vcf.gz

    cat > ${sample}.pbsv_call.versions.yml <<-END_VERSIONS
    "${task.process}":
      pbsv: "\$(pbsv --version 2>&1 | head -n 1)"
      bcftools: "\$(bcftools --version | head -n 1)"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.sv.pass.vcf.gz ${sample}.sv.pass.vcf.gz.tbi
    echo '${task.process}: {stub: true}' > ${sample}.pbsv_call.versions.yml
    """
}
