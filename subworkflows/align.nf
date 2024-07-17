process call_paftools {
    label "singlecell"
    cpus 1
    input:
        path "ref_genes.gtf"
    output:
        path "ref_genes.bed", emit: ref_genes_bed
    """
    paftools.js gff2bed -j ref_genes.gtf > ref_genes.bed
    """
}

process get_chrom_sizes{
    label "singlecell"
    cpus 1
    input:
        path "ref_genome.fai"
    output:
        path 'chr_sizes', emit: ref_chrom_sizes
    """
    cut -f1,2 ref_genome.fai | sort -V > chr_sizes
    """
}

process align_to_ref {
    label "singlecell"
    cpus params.resources_mm2_max_threads
    input:
        tuple val(sample_id),
              path("reads.fastq")
        path "ref_genome.fasta"
    output:
        tuple val(sample_id), 
            path("${sample_id}.sorted.bam"), 
            path("${sample_id}.sorted.bam.bai"), 
            emit: bam_sort
    """
    minimap2 -ax map-ont -k 17 -t ${task.cpus} -L -y --secondary=no --MD --cap-kalloc=1g -K 10g ref_genome.fasta reads.fastq* \
        | samtools sort -@ 2 -o ${sample_id}.sorted.bam
    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """
}

// workflow module
workflow align {
    take:
        stranded_fq
        ref_genome
    main:
        //stranded_fq.groupTuple().view()
        align_to_ref(
            stranded_fq.groupTuple(),
            ref_genome)
    emit:
        bam_sort = align_to_ref.out.bam_sort
}
