#!/bin/bash

CONFIG=ronin.config
INPUT_FILE=inputFile.tsv
TMP_REF_DIR=/tmp/refs
FASTA_NAME=GCF_002263795.2_ARS-UCD1.3_nuclear-GCF_000754665.1_Bison_UMD1.0_mt.fna
REF_S3_PATH=s3://project-blue-babe.store.uoa.cloud/refs/cattle_bison_ref
OUTPUT_S3=s3://project-blue-babe.store.uoa.cloud/ronin_bwa_outputs

nextflow run ronin_bwa.nf \
    -resume \
    -with-singularity \
    -c ${CONFIG} \
    --input ${INPUT_FILE} \
    --ref_dir ${TMP_REF_DIR} \
    --fasta_name ${FASTA_NAME} \
    --ref_s3_path ${REF_S3_PATH} \
    --bwaalnn 0.01 --bwaalno 2 --bwaalnl 1024 \
    --bwaalnk 2 --output_s3 ${OUTPUT_S3} -with-trace

