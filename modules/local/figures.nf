process MAKE_FIGURE_1 {
    label 'process_figure'
    conda "${moduleDir}/../../envs/figures.yml"

    publishDir "${params.outdir}/07_figures", mode: params.publish_mode

    input:
    path phased_files
    path cnv_files
    path methylation_files
    path gtf, stageAs: 'reference/genes.gtf'
    path metadata, stageAs: 'reference/metadata.csv'

    output:
    path "figure_1", emit: results

    script:
    """
    mkdir -p phased cnv methylation figure_1
    for f in ${phased_files}; do ln -sf "\$(realpath \"\$f\")" phased/; done
    for f in ${cnv_files}; do ln -sf "\$(realpath \"\$f\")" cnv/; done
    for f in ${methylation_files}; do ln -sf "\$(realpath \"\$f\")" methylation/; done

    python ${moduleDir}/../../scripts/05_make_figures/FIGURE_1.py \
        --vcf-dir phased \
        --bam-dir phased \
        --methylation-dir methylation \
        --cnv-dir cnv \
        --gtf ${gtf} \
        --metadata ${metadata} \
        --outdir figure_1
    """

    stub:
    """
    mkdir -p figure_1/tables figure_1/figures
    touch figure_1/tables/Figure1C_parental_assignment.tsv figure_1/figures/Figure1.png
    """
}

process MAKE_FIGURE_2 {
    label 'process_figure'
    conda "${moduleDir}/../../envs/figures.yml"

    publishDir "${params.outdir}/07_figures", mode: params.publish_mode

    input:
    path cnv_files
    path methylation_files
    path gtf, stageAs: 'reference/genes.gtf'
    path metadata, stageAs: 'reference/metadata.csv'
    path icr_bed, stageAs: 'reference/icr.bed'
    path segdup_bed, stageAs: 'reference/segdup.bed'

    output:
    path "figure_2", emit: results

    script:
    """
    mkdir -p cnv methylation figure_2
    for f in ${cnv_files}; do ln -sf "\$(realpath \"\$f\")" cnv/; done
    for f in ${methylation_files}; do ln -sf "\$(realpath \"\$f\")" methylation/; done

    python ${moduleDir}/../../scripts/05_make_figures/FIGURE_2.py \
        --methylation-dir methylation \
        --metadata ${metadata} \
        --gtf ${gtf} \
        --cnv-dir cnv \
        --icr-bed ${icr_bed} \
        --segdup-bed ${segdup_bed} \
        --outdir figure_2
    """

    stub:
    """
    mkdir -p figure_2/figures
    touch figure_2/figures/Figure2_reciprocal_cis_architecture.png
    """
}

process MAKE_FIGURE_3 {
    label 'process_figure'
    conda "${moduleDir}/../../envs/figures.yml"

    publishDir "${params.outdir}/07_figures", mode: params.publish_mode

    input:
    path methylation_files
    path figure_1_results
    path gtf, stageAs: 'reference/genes.gtf'
    path segdup_bed, stageAs: 'reference/segdup.bed'
    path imprintome_bed, stageAs: 'reference/imprintome.bed'
    path icr_bed, stageAs: 'reference/icr.bed'
    path repeats_bed, stageAs: 'reference/repeats.bed'

    output:
    path "figure_3", emit: results

    script:
    """
    mkdir -p methylation figure_3
    for f in ${methylation_files}; do ln -sf "\$(realpath \"\$f\")" methylation/; done

    python ${moduleDir}/../../scripts/05_make_figures/FIGURE_3.py \
        --methylation-dir methylation \
        --assignment-table ${figure_1_results}/tables/Figure1C_parental_assignment.tsv \
        --gtf ${gtf} \
        --segdup ${segdup_bed} \
        --imprintome-bed ${imprintome_bed} \
        --court2014-bed ${icr_bed} \
        --icr-bed ${icr_bed} \
        --repeats ${repeats_bed} \
        --outdir figure_3
    """

    stub:
    """
    mkdir -p figure_3/figures
    touch figure_3/figures/Figure3_boundary_mapping_improved.png
    """
}

process MAKE_FIGURE_4 {
    label 'process_figure'
    conda "${moduleDir}/../../envs/figures.yml"

    publishDir "${params.outdir}/07_figures", mode: params.publish_mode

    input:
    path table_dir

    output:
    path "figure_4", emit: results

    script:
    """
    python ${moduleDir}/../../scripts/05_make_figures/FIGURE_4.py \
        --table-dir ${table_dir} \
        --outdir figure_4
    """

    stub:
    """
    mkdir -p figure_4/figures
    touch figure_4/figures/Figure4.png
    """
}

process MAKE_FIGURE_5 {
    label 'process_figure'
    conda "${moduleDir}/../../envs/figures.yml"

    publishDir "${params.outdir}/07_figures", mode: params.publish_mode

    input:
    path table_dir
    path reference, stageAs: 'reference/genome.fa'
    path gtf, stageAs: 'reference/genes.gtf'

    output:
    path "figure_5", emit: results

    script:
    """
    python ${moduleDir}/../../scripts/05_make_figures/FIGURE_5.py \
        --table-dir ${table_dir} \
        --fasta ${reference} \
        --gtf ${gtf} \
        --outdir figure_5
    """

    stub:
    """
    mkdir -p figure_5/figures
    touch figure_5/figures/Figure5_v7.png
    """
}
