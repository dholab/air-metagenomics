#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_reads = Channel
        .fromPath( params.samplesheet )
        .splitCsv( header: true )
        .map { row -> tuple( row.raw_read_label, row.sample_id ) }
	
	// Workflow steps 
    MERGE_FASTQS (
        ch_reads
    )

    SAMPLE_QC (
        MERGE_FASTQS.out
    )

    MAP_TO_REFSEQS (
        SAMPLE_QC.out
    )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.merged_fastqs = params.results + "/1_merged_fastqs"
params.filtered_fastqs = params.results + "/2_filtered_fastqs"
params.bams = params.results + "/3_alignment_maps"
// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process MERGE_FASTQS {
	
	/* 
    Here we concatenate all the FASTQ files for each barcode,
    and name them according to the sample ID specified in the 
    samplesheet.
    */ 
	
	tag "${sample_id}"
    publishDir params.merged_fastqs, mode: 'symlink'
	
	input:
	tuple val(label), val(sample_id)
	
	output:
	tuple path("*.fastq.gz"), val(sample_id)
	
	script:
	"""
	cat ${params.raw_reads}/${label}*.fastq.gz > ${sample_id}.fastq.gz
	"""
}


process SAMPLE_QC {
	
	/*
    Here we run some trimming and quality filtering with the bbmap
    script `reformat.sh`
    */
	
	tag "${sample_id}"
    publishDir params.results, mode: 'copy'
	
	input:
	tuple path(fastq), val(sample_id)
	
	output:
	tuple path("*.fastq.gz"), val(sample_id)
	
	script:
	"""
	reformat.sh in=${fastq} \
    out=${sample_id}_filtered.fastq.gz \
    forcetrimleft=30 forcetrimright2=30 \
    mincalledquality=7 minlength=200 qin=33
	"""
}


process MAP_TO_REFSEQS {
	
	/* 
    This process maps each sample ID's reads to a reference FASTA 
    of 835 human viruses. The resulting BAMs will serve as "hits"
    for which viruses were present in a given air sample.
    */
	
	tag "${sample_id}"
    publishDir params.results, mode: 'copy'
	
	input:
	tuple path(fastq), val(sample_id)
	
	output:
	tuple path("*.bam*"), val(sample_id)
	
	script:
	"""
    minimap2 \
    -ax map-ont \
    ${params.virus_ref} \
    ${fastq} \
    --eqx \
    | reformat.sh \
    in=stdin.sam \
    ref=${params.virus_ref} \
    out=${sample_id}_filtered.bam \
    mappedonly=t
	"""
}


// --------------------------------------------------------------- //