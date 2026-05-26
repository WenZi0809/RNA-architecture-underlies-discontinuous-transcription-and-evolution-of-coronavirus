#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow {
    index_file = buildIndex(params.refseq)
    bam_file = alignAndFilter(index_file, params.fqdata, params.title)
    sam_file = bamToSam(bam_file, params.title)
    bed_file = samToBed(params.v_chr, sam_file, params.title)
    countGene(params.v_chr, bed_file, params.title)
}

process buildIndex {
    publishDir "${params.outdir}/reference", mode: 'move'

    input:
    path refseq

    output:
    path "${refseq.baseName}.mmi"

    script:
    """
    minimap2 -d ${refseq.baseName}.mmi ${refseq}
    """
}

process alignAndFilter {
    publishDir "${params.outdir}/map", mode: 'move'

    input:
    path index_file
    path fqdata
    val title

    output:
    path "${title}_mapped.bam"

    script:
    """
    minimap2 -ax splice ${index_file} ${fqdata} | \
        samtools view -hF4 -b | \
        samtools sort -o ${title}_mapped.bam
    """
}

process bamToSam {
    publishDir "${params.outdir}/map", mode: 'move'

    input:
    path bam_file
    val title

    output:
    path "${title}_mapped.sam"

    script:
    """
    samtools view -h -o ${title}_mapped.sam ${bam_file}
    """
}

process samToBed {
    publishDir "${params.outdir}/analysis", mode: 'move'

    input:
    val v_chr
    path sam_file
    val title

    output:
    path "${title}_mapped.bed"

    script:
    """
    python sam_to_bed_linux.py ${v_chr} ${sam_file} ${title}_mapped.bed
    """
}

process countGene {
    publishDir "${params.outdir}/analysis", mode: 'move'

    input:
    val v_chr
    path bed_file
    val title

    output:
    path "${title}-2.bed"
    path "${title}-3.bed"

    script:
    """
    python counts_gene_linux.py ${v_chr} ${bed_file} ${title}
    """
}