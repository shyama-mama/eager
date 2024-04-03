#!/usr/bin/env nextflow

// Small nextflow workflow to run bwa aln on RONIN 

Channel.fromPath(params.input)
    .splitCsv(sep: '\t')
    .map {
        it ->
            def samplename = it[0]
            def libraryid = it[1]
            def seqtype = it[2]
            def shard_idx = it[3]
            def r1 = it[4]
            [ samplename, libraryid, seqtype, shard_idx, r1 ]
    }
    .set { ch_bwa_input } 
    
process bwa_ronin {
    tag "${libraryid}"

    input:
    tuple samplename, libraryid, seqtype, shard_idx, r1 from ch_bwa_input

    script:
    def fasta = "${params.ref_dir}/${params.fasta_name}"
    def output_bam = "${libraryid}_${shard_idx}.mapped.bam"
    """
    # get r1 name
    r1_name=\$(basename $r1)

    # Copy input file from S3
    aws s3 cp ${r1} .

    # Check if ref exists
    if [ ! -f ${fasta} ]; then
        mkdir -p ${params.ref_dir}
        aws s3 cp --recursive ${params.ref_s3_path}/ ${params.ref_dir}
    fi

    bwa aln -t ${task.cpus} ${fasta} \${r1_name} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -o ${params.bwaalno} -f ${libraryid}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${samplename}\\tPL:illumina\\tPU:ILLUMINA-${libraryid}-${seqtype}" $fasta ${libraryid}.sai \$r1_name | samtools sort -@ ${task.cpus - 1} -O bam - > ${output_bam}
    samtools index ${output_bam}

    # Copy output file to S3
    aws s3 cp ${output_bam} ${params.output_s3}/${libraryid}/
    aws s3 cp ${output_bam}.bai ${params.output_s3}/${libraryid}/

    # clean up folder by removing fastq, sai, bam and bai
    rm -rf *gz 
    rm -rf *sai
    rm -rf *bam* 
    """    
}