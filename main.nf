#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_reads = Channel
        .fromPath( params.samplesheet )
        .splitCsv( header: true )
        .map { row -> 
            def parent_dir = row.parent_dir?.trim() ?: '/scratch'
            tuple( row.raw_read_label, row.sample_id, parent_dir ) 
        }
    
    ch_ref_seqs = Channel
        .fromPath( params.virus_ref )
    
    ch_contaminants = Channel
        .fromPath( params.contaminants_tar )
        // .splitFasta( by: 10, file: true )
        // .collect()
	
	// Workflow steps 
    FIND_AND_MERGE_FASTQS (
        ch_reads
    )

    SAMPLE_QC (
        FIND_AND_MERGE_FASTQS.out
    )

    FIND_NTC (
        FIND_AND_MERGE_FASTQS.out
    )

    CONVERT_TO_FASTA (
        SAMPLE_QC.out
    )

    DECOMPRESS_CONTAMINANTS (
        ch_contaminants
    )

    REMOVE_CONTAMINANTS (
        CONVERT_TO_FASTA.out,
        DECOMPRESS_CONTAMINANTS.out
    )

    REMOVE_NTC (
        REMOVE_CONTAMINANTS.out,
        FIND_NTC.out
    )

    MAP_TO_REFSEQS (
        REMOVE_NTC.out,
        ch_ref_seqs
    )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// specifying whether to run in low disk mode
if( params.low_disk_mode == true ) {
	params.publishMode = 'symlink'
}
else {
	params.publishMode = 'copy'
}

// Resources subdirectories
params.contam_ref = params.resources + "/contam_ref"

// Results subdirectories
params.merged_fastqs = params.results + "/1_merged_fastqs"
params.filtered_fastqs = params.results + "/2_filtered_fastqs"
params.fasta_cleaning = params.results + "/3_cleaned_fastas"
params.bams = params.results + "/4_alignment_maps"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process FIND_AND_MERGE_FASTQS {
	
	/*
    Here we determine if the read labels are SRA accessions or 
    local file names. If they are SRA accessions, they will be 
    downloaded and merged automatically from SRA servers. If 
    they are local, they will be located and merged.
    */ 
	
	tag "${sample_id}"
    publishDir params.merged_fastqs, mode: params.publishMode, overwite: true

    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 2
	
	input:
	tuple val(label), val(sample_id), path(parent_dir)
	
	output:
	tuple path("${sample_id}.fastq.gz"), val(sample_id)
	
	script:
    if ( label.startsWith("SRR") )
        """
        prefetch ${label}
        fasterq-dump ${label}/${label}.sra \
        --concatenate-reads --skip-technical --quiet && \
        gzip --no-name ${label}.fastq
        mv ${label}.fastq.gz ${sample_id}.fastq.gz
        rm -rf ${label}/
        rm -rf fasterq.tmp.*
        """
    else
        """
        find `realpath ${parent_dir}` -type f -name ${label}*.fastq.gz > fastq_list.txt
        touch ${sample_id}.fastq
        touch merged_list.txt 
        for i in `cat fastq_list.txt`;
        do
            echo "\$i" >> merged_list.txt 
            zcat \$i >> ${sample_id}.fastq
        done
        if [[ `cat merged_list.txt | wc -l` -eq `cat fastq_list.txt | wc -l` ]]; then
            gzip --no-name ${sample_id}.fastq
        else
            echo "Merging failed."
            exit 1
        fi
        """

}


process SAMPLE_QC {
	
	/*
    Here we run some trimming and quality filtering with the bbmap
    script `reformat.sh` for PacBio or Illumina reads, or cutadapt
    for ONT reads.
    */
	
	tag "${sample_id}"
    publishDir params.filtered_fastqs, mode: 'copy'
    
    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }
    maxRetries 2
	
	input:
	tuple path(fasta), val(sample_id)
	
	output:
	tuple path("*.fastq.gz"), val(sample_id)

    when:
    !sample_id.contains("NTC_")
	
	script:
    if ( params.ont == true )
        """
        reformat.sh in=${fasta} \
        out=${sample_id}_filtered.fastq.gz \
        forcetrimleft=30 forcetrimright2=30 \
        mincalledquality=9 qin=33 minlength=200 
        """
    else if( params.pacbio == true )
        """
        reformat.sh in=${fasta} \
        out=${sample_id}_filtered.fastq.gz \
        forcetrimleft=30 forcetrimright2=30 \
        mincalledquality=9 qin=33 minlength=200 
        """
    else 
        """
        reformat.sh in=${fasta} \
        out=${sample_id}_filtered.fastq.gz \
        forcetrimleft=30 forcetrimright2=30 \
        mincalledquality=9 qin=33 minlength=200 
        """

}


process FIND_NTC {

    /*
    Here we use filenames to identify which of the samples from your
    sequencing run were negative, no-template controls. One of these
    is required for each sequencing run and the associated workflow run.
    We have written this workflow to be run once for each sequencing run,
    and as such, it will produce errors if you put multiple sequencing
    runs worth of samples and controls in the same samplesheet.
    */

    tag "${sample_id}"

    input:
    tuple path(fastq), val(sample_id)

    output:
    path "NTC_*.fasta.gz"

    when:
    sample_id.contains("NTC_")

    script:
    """
    fastq_to_fasta.py \
    && gzip --no-name ${sample_id}.fasta
    """

}


process CONVERT_TO_FASTA {

    /*
    To save space, FASTQs are converted to FASTAs after QC.
    */

    tag "${sample_id}"
    publishDir params.fasta_cleaning, mode: params.publishMode, overwrite: true

    input:
    tuple path(fastq), val(sample_id)
    
    output:
    tuple path("*.fasta.gz"), val(sample_id)

    script:
    """
    fastq_to_fasta.py \
    && gzip --no-name ${sample_id}_filtered.fasta
    """
}


process DECOMPRESS_CONTAMINANTS {

    /*
    A tarball of potential library contaminants, including PhiX adapters,
    host reads, metagenomic contaminants, etc. is included with this workflow.
    To access them, this process decompresses them and sends them to downsream
    processes.
    */

    publishDir params.resources, mode: 'copy', overwrite: false

    input:
    path tar

    output:
    path "${folder_name}"

    script:
    folder_name = tar.getSimpleName()
    """
    tar -xvf ${tar}
    """

}


process REMOVE_CONTAMINANTS {

    /*
    map reads to NVD contaminant databases
    return unmapped reads
    this removes kit-ome and human sequences
    */

    tag "${sample_id}"
    publishDir params.fasta_cleaning, mode: params.publishMode, overwrite: true

    cpus 8

    input:
    tuple path(fasta), val(sample_id)
    path contaminants

    output:
	tuple path("${sample_id}_contam_removed.fasta.gz"), val(sample_id)

    script:
    """
    mv ${fasta} tmp.fasta.gz
    ls `realpath ${contaminants}/*.fa.gz` > contaminant_file_paths.txt
    for i in `cat contaminant_file_paths.txt`; 
    do
        basename=`basename \$i`
        echo "Now mapping to " \$basename
        mv tmp.fasta.gz tmp_\$basename.fasta.gz
        minimap2 -ax map-ont --eqx --secondary=no -t ${task.cpus} \$i \
        tmp_\$basename.fasta.gz \
        | reformat.sh unmappedonly=t in=stdin.sam \
        ref=\$i \
        out=tmp.fasta.gz
    done && \
    mv tmp.fasta.gz ${sample_id}_contam_removed.fasta.gz
    """

}


process REMOVE_NTC {

    /*
    map reads from sample to no-template water control reads
    remove reads that map to sequences found in no template control
    */
    
    tag "${sample_id}"
    publishDir params.fasta_cleaning, mode: params.publishMode, overwrite: true

    tag "${sample_id}"

    cpus 8

    input:
    tuple path(fasta), val(sample_id)
    path ntc

    output:
	tuple path("${sample_id}_ntc.fasta.gz"), val(sample_id)

    script:
    """
    minimap2 -ax map-ont --eqx --secondary=no -t ${task.cpus} \
    ${ntc} \
    ${fasta} \
    | reformat.sh unmappedonly=t in=stdin.sam \
    ref=${ntc} \
    out=${sample_id}_ntc.fasta.gz
    """

}


process MAP_TO_REFSEQS {
	
	/* 
    This process maps each sample ID's reads to a reference FASTA 
    of 835 human viruses. The resulting BAMs will serve as "hits"
    for which viruses were present in a given air sample.
    */
	
	tag "${sample_id}"
    publishDir params.bams, mode: 'copy', overwite: true

    errorStrategy 'ignore'
    
    cpus 4
	
	input:
	tuple path(fasta), val(sample_id)
    each path(refseq)
	
	output:
	tuple path("*.bam*"), val(sample_id)
	
	script:
	"""
    minimap2 \
    -ax map-ont \
    ${refseq} \
    ${fasta} \
    --eqx \
    -t 3 \
    | reformat.sh \
    in=stdin.sam \
    ref=${refseq} \
    out=${sample_id}_filtered.bam \
    mappedonly=t
	"""
}


// --------------------------------------------------------------- //