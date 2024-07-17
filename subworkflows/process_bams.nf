import java.util.ArrayList;

process get_contigs {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              path(bam),
              path(bai)

    output:
        path("${sample_id}_contigs"),
        emit: contigs
    """
    samtools idxstats ${bam} \
        | gawk -v var=${sample_id} '/^[^*]/{print var,\$1}' \
        | gawk NF > "${sample_id}_contigs"
    """
}

process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "singlecell"
    cpus 2
    input:
        tuple path(bam),
              path(bai),
              val(meta),
              val(chr)
        path barcodes_dir

    output:
        // TODO: Do not write bams. Write mapping of read_id to barcode
        tuple val(meta.sample_id), 
              path("${meta.sample_id}_${chr}.bc_extract.sorted.tsv"),
              val(chr),
              emit: bc_uncorr_tsv
        tuple val(meta.sample_id),
              path("${meta.sample_id}_${chr}.uncorrected_bc_counts.tsv"), emit: barcode_counts

    """
    export NUMBA_CACHE_DIR="/staging/leuven/stg_00104/projects/demeulemeester_multiome/atac_testruns/tmp"
    workflow-glue extract_barcode \
    ${bam} ${barcodes_dir}/${meta['bc_long_list']}\
    -t $task.cpus \
    --kit ${meta['kit_name']} \
    --adapter1_suff_length $params.barcode_adapter1_suff_length \
    --min_barcode_qv $params.barcode_min_quality \
    --barcode_length ${meta['barcode_length']} \
    --umi_length ${meta['umi_length']} \
    --output_read_tags "${meta.sample_id}_${chr}.bc_extract.sorted.tsv" \
    --output_barcode_counts "${meta.sample_id}_${chr}.uncorrected_bc_counts.tsv" \
    --contig ${chr}
    """
}

process generate_whitelist{
    label "singlecell"

    publishDir "${params.outdir}/barcodes", mode:'copy'

    cpus 1
    input:
        tuple path(counts),
              val(meta)
    output:
        tuple val(meta.sample_id), 
              path("${meta.sample_id}.whitelist.tsv"), 
              emit: whitelist
        tuple val(meta.sample_id), 
              path("${meta.sample_id}.kneeplot.png"), 
              emit: kneeplot
    """
    workflow-glue knee_plot \
        ${counts} \
        --exp_cells ${meta['exp_cells']} \
        --output_whitelist "${meta.sample_id}.whitelist.tsv" \
        --output_plot "${meta.sample_id}.kneeplot.png"
    """
}

process assign_barcodes{
    label "singlecell"
    cpus 1
    input:
         tuple val(sample_id),
               path(whitelist),
               path(extract_barcodes),
               val(chr)
    output:
        tuple val(sample_id),
              val(chr),
              path("${sample_id}_${chr}_bc_assign_counts.tsv"),
              emit: chrom_assigned_barcode_counts
        tuple val(sample_id),
              val(chr),
              path("${sample_id}_${chr}_extract_barcodes_with_bc.tsv"),
              emit: tags
    """
    workflow-glue assign_barcodes \
        --output_tags ${sample_id}_${chr}_extract_barcodes_with_bc.tsv \
        --output_counts ${sample_id}_${chr}_bc_assign_counts.tsv \
        --max_ed $params.barcode_max_ed \
        --min_ed_diff $params.barcode_min_ed_diff \
        --extract_barcode_tags ${extract_barcodes} \
        --whitelist ${whitelist}
    """
}

process tag_bams {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path(bam),
              path(bai),
              path(tags)
    output:
         tuple val(sample_id),
               val(chr),
               path("${sample_id}_${chr}.tagged.bam"),
               path("${sample_id}_${chr}.tagged.bam.bai"),
               emit: tagged_bam
    """
    workflow-glue tag_bam \
        --in_bam ${bam} \
        --tags ${tags} \
        --out_bam ${sample_id}_${chr}.tagged.bam

    samtools index ${sample_id}_${chr}.tagged.bam
    """
}
/*
process deduplicate_bams {
    // Deduplicate bamfiles using gatk (picard) MarkDuplicates
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path("${sample_id}.${chr}.tagged.bam"),
              path("${sample_id}.${chr}.tagged.bam.bai")
    output:
        tuple val(sample_id),
              path("${sample_id}.${chr}.tagged.sorted.dedup.bam"),
              path("${sample_id}.${chr}.tagged.sorted.dedup.bam.bai"),
              emit: dedup_bam
        tuple val(sample_id),
              path("${sample_id}.dedup_metrics.txt"),
              emit: dedup_metrics
    
    """
    gatk MarkDuplicates -I ${bam} -O ${sample_id}.tagged.sorted.dedup.bam -M ${sample_id}.dedup_metrics.txt \
    --BARCODE_TAG CB
    """
}
*/

process combine_chrom_bams {
    // Merge all chromosome bams by sample_id
    label "singlecell"
    cpus Math.min(8, params.max_threads)

    publishDir "${params.outdir}/bamfiles", mode:'copy'

    input:
        tuple val(sample_id),
              val(chr),
              path(chrom_bams)
    output:
        tuple val(sample_id), 
              path("${sample_id}_merged.dedup.bam"), 
              path("${sample_id}_merged.dedup.bam.bai"),
              emit: bam_fully_tagged
    """
    samtools merge -@ ${task.cpus} -o "${sample_id}_merged.dedup.bam" ${chrom_bams}; 
    samtools index -@ ${task.cpus} "${sample_id}_merged.dedup.bam";
    """
}

process deduplicate_bams {
    tag "MarkDuplicatesSpark on ${sample}"
    //label 'process_high'
    cpus 5
    container "https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0"
    if (!params.merge_bam) {
        publishDir "${params.outdir}/bamfiles", mode:'copy', pattern: '*.bam*'
    }
    publishDir "${params.outdir}/metrics", mode:'copy', pattern: '*.metrics'

    input:
        tuple val(sample),
              val(chr),
              path(bam),
              path(bai)
        val stringency

    output:
        tuple val(sample),
              val(chr),
              path("${sample}_${chr}.dedup.bam"),
              emit: dedup_bam
        tuple val(sample),
              val(chr),
              path("${sample}_${chr}.dedup.metrics"),
              emit: metrics

    script:
        def avail_mem = 24576 //TODO: Change this back to 3072 and set proper memory labels also in script section
        if (!task.memory) {
            log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        """
        gatk --java-options "-Xmx24576M -XX:-UsePerfData" \\
            MarkDuplicates \\
            --INPUT ${bam} \\
            --OUTPUT ${sample}_${chr}.dedup.bam \\
            --BARCODE_TAG CB \\
            --METRICS_FILE ${sample}_${chr}.dedup.metrics \\
            --VALIDATION_STRINGENCY ${stringency} \\
            --TMP_DIR . 
        """
}

workflow process_bams {
    take:
        bam
        meta
        bc_longlist_dir
   
    main:

        get_contigs(bam)
        
        contigs = get_contigs.out.contigs
            .splitCsv(sep: " ")

        extract_barcodes(
            bam
            .cross(
                meta
                .cross(contigs).map{it -> it.flatten()})
                .map{it -> it.flatten()[1, 2, 4, 6]},
            bc_longlist_dir)

        generate_whitelist(
            extract_barcodes.out.barcode_counts
            .collectFile()
            .map {it -> tuple(it.getSimpleName(), it)}
            .join(meta).map {it -> it.tail()}) // Remove sample_id

       assign_barcodes(
            generate_whitelist.out.whitelist
            .cross(extract_barcodes.out.bc_uncorr_tsv)
            .map {it -> it.flatten()[0, 1, 3, 4]})

        tag_bams(
            bam.cross(
                assign_barcodes.out.tags
            ).map {it ->it.flatten()[0, 4, 1, 2, 5]})

        dedup_bam = deduplicate_bams(tag_bams.out.tagged_bam, params.validation_stringency)

        if (params.merge_bam) {
            combine_chrom_bams(dedup_bam.dedup_bam
                .groupTuple())
            // [sample_id, bam]
            tagged_bams = combine_chrom_bams.out.bam_fully_tagged

        } else {
            tagged_bams = dedup_bam.dedup_bam
            // [sample_id, bam, bai]
            .map {it -> it[0, 1, 2]}
            .groupTuple()
        }

    emit:
        results = assign_barcodes.out.tags
            .join(tagged_bams)
            .map{it -> it.flatten()}
}
