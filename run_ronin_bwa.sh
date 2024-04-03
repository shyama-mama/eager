#!/bin/bash

nextflow run ronin_bwa.nf \
    -resume \
    -with-singularity \
    -c ${CONFIG} \
    -profile singularity \
    --input ${INPUT_FILE} \
    --ref_dir ${TMP_REF_DIR} \
    --fasta_name ${FASTA_NAME} \
    --ref_s3_path ${REF_S3_PATH} \
    --bwaalnn 0.01 --bwaalno 2 --bwaalnl 1024 \
    --bwaalnk 2 --output_s3 ${OUTPUT_S3}