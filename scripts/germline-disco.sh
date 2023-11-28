#!/usr/bin/env bash

####################### NOTE ########################
#   This is a general script for example purposes   #
# Please modify it to fit your experimental design! #
##################################################### 

# Collated from Mohamed's pipeline scripts
# Refer to GATK's docs: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels
# Arguments (positional):
#   1: Reference FASTA
#   2: Directory to be created for GenomicsDB
#   3+: List of raw unmapped read files (uBAM or FASTQ)
# Example usage:
#   ./germline-disco.sh reference.fa your_database sample1.fq sample2.fq sample3.fq

GATK=/programs/gatk4/gatk
if ! [ -f $GATK ]; then
  GATK=gatk
fi

BWA=/programs/bwa-mem2-2.2.1/bwa-mem2
if ! [ -f $BWA ]; then
  BWA=bwa-mem2
fi

REFERENCE_FASTA=$1
GENOMICS_DB_PATH=$2

$BWA index $REFERENCE_FASTA
samtools faidx $REFERENCE_FASTA
$GATK CreateSequenceDictionary -R $REFERENCE_FASTA

for SAMPLE_FILE in "${@:3}"
do
    SAMPLE_NAME=$(basename -- "$SAMPLE_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%.*}"

    $BWA mem $REFERENCE_FASTA $SAMPLE_FILE > $SAMPLE_NAME.sam

    $GATK AddOrReplaceReadGroups -I $SAMPLE_NAME.sam -O ${SAMPLE_NAME}.readgroup.sam -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM $SAMPLE_NAME

    samtools sort -O BAM ${SAMPLE_NAME}.readgroup.sam > $SAMPLE_NAME.bam

    $GATK MarkDuplicates --INPUT $SAMPLE_NAME.bam --METRICS_FILE ${SAMPLE_NAME}.dedup.metrics --OUTPUT ${SAMPLE_NAME}.dedup.bam 

    samtools index ${SAMPLE_NAME}.dedup.bam

    $GATK HaplotypeCaller -R $REFERENCE_FASTA -I ${SAMPLE_NAME}.dedup.bam -O ${SAMPLE_NAME}.gvcf.gz -ERC GVCF
done

# NOTE: GenomicsDBImport requires an interval; in our case we will
#       iterate through all contigs defined in the reference index (.fa.fai)
#       and import them individually - you may not want this!
for CONTIG in $(cut -f1 "${REFERENCE_FASTA}.fai")
do
    VARIANT_ARGS=$(awk 'BEGIN{ORS=" "}; {print "-V " $0}' <(find . -type f -name "*.gvcf.gz"))
    $GATK GenomicsDBImport $VARIANT_ARGS --genomicsdb-workspace-path $GENOMICS_DB_PATH -L $CONTIG
done

for SAMPLE_FILE in "${@:3}"
do
    SAMPLE_NAME=$(basename -- "$SAMPLE_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%.*}"

    $GATK GenotypeGVCFs -R $REFERENCE_FASTA -V gendb://$GENOMICS_DB_PATH -O ${SAMPLE_NAME}.vcf.gz
done
