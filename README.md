# General Introduction 
This is a version of the EAGER pipeline that lets you parallelize alignment with `bwa aln` and deduplication.

Alignment for ancient samples take a long time due to misincorporations from DNA damage and various contaminating sequences. 
Some HPCs may have a time limit for how long a task can run. Samples with a lot of reads could go over the limit, sometimes taking more than 72 hours to finish alignment. 

To overcome this limitation you can specify `--shard_bwa` and optionally `--chunk_size $n` to subset the input FASTQ with `n` reads per chunk (chunk is used synonymously with file) before aligning. The default is 1,000,000 reads per chunk. 

Additionally, the deduplication step can also be sharded. It works for both `MarkDuplicates` and `dedup` tools. Running `--shard_deduplication` will split the input BAM by contig and deduplicate the contigs in parallel. Then, the metrics and results BAMs are merged into one. A `--contig_file ${file}`, that contains 1 contig per line or multiple contigs in a line space separated, can be provided to only deduplicate and keep the specified contigs. 


## Details and Caveats of BWA sharding
### Details 
Running with `--shard_bwa` runs an extra process called `shardfastqs` which uses `seqtk split2 -s ${chunk_size}` to subset the input FASTQ into `n` sequences per chunk. 

`bwa aln` is then run on each chunk of the FASTQ in parallel. For example, a FASTQ with 10,000,000 reads can take 40 real time hours to align. Splitting into 10 chunks and aligning the 10 chunks in parallel can drastically reduce the alignment time. Assuming it takes 4 hours to align 1,000,000 million reads and all chunks are being aligned at the same time in parallel, alignment can finish in 4 real time hours for that sample. 

Post `bwa aln` the BAMs are merged using `samtools merge`

### Caveats and known issues 
Very large samples can create thousands of shards, which when merging creates a very large BAM header. So far `qualimap`, `damageprofiler` and `mtnucratio` error with `ran out of heap size`. This can be simply fixed with reheadering the BAM file. This will be automatically implemented in future releases. 

## Details and Caveats of Deduplication sharding
### Details
BAM files are split either by `bamtools split -reference`, when no `--contig_file` is provided or with `samtools view` when contigs are provided.

### Caveats and known issues.
Using `--shard_deduplication` without `--contig_file` can be dangerous. Most reference genomes contain the main contigs, some decoy contigs and alternate contigs, leading to more than 10,000 contigs in some cases. This means there will be one deduplication job for each contig. To prevent this, it is recommended to provide a `--contig_file` with one main contig per line and alternate contigs can be ignored from the file in cases it will not be used or alternate contigs can be provided in one line space separated. Multiqc summary table has funky sample names. Aim to fix this in the next release.

# Tutorial only for Phoenix HPC and ACAD users 
```
EXAMPLE_DIR=~/eager_example/
mkdir -p ${EXAMPLE_DIR}
cd ${EXAMPLE_DIR}

# Start Screen session
screen -S Example_Eager

# load Nextflow module which is needed to run EAGER
module load Nextflow/21.03.0 

INPUT_FILE=/hpcfs/groups/acad_users/shyrav/resources/eager_tutorial/example_input.tsv

FASTA=/hpcfs/groups/acad_users/shyrav/reference_genomes/human_homo_sapiens/hg38/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna
FASTA_INDEX=/hpcfs/groups/acad_users/shyrav/reference_genomes/human_homo_sapiens/hg38/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna.fai
SEQ_DICT=/hpcfs/groups/acad_users/shyrav/reference_genomes/human_homo_sapiens/hg38/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.dict
BWA_INDEX=/hpcfs/groups/acad_users/shyrav/reference_genomes/human_homo_sapiens/hg38/data/GCF_000001405.26/

CONTIG_FILE=/hpcfs/groups/acad_users/shyrav/reference_genomes/human_homo_sapiens/hg38/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.contig.list



nextflow run /hpcfs/groups/acad_users/shyrav/resources/eager-v2.4.5-sharding -c /hpcfs/groups/acad_users/nextflow_repos/phoenix.config -outdir ./ \
	-with-singularity \
	--input ${INPUT_FILE} \
	--fasta ${FASTA} \
	--fasta_index ${FASTA_INDEX} \
	--bwa_index ${BWA_INDEX} \
    --seq_dict ${SEQ_DICT} \
	--mapper 'bwaaln' \
    --shard_bwa \
	--bwaalnn 0.01 \
    --bwaalno 2 \
   	--bwaalnl 1024 \
	--dedupper 'markduplicates' \
    --contig_file ${CONTIG_FILE} \
    --shard_deduplication \
	--qualitymax 64 

# To exit screen 
# Ctrl + A and then D
```
