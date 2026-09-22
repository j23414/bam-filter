#! /usr/bin/env bash
set -euv

export REF=reference.fasta
export PRIMERS=primers_wo_restriction_enzyme.fasta
# module load bwa-mem2/2.2.1
bwa-mem2 index ${REF}

# -T 10, -k 10 to reduce the threshold since primer sequences are short
# -a print all matches, not just the top for primers matching multiple segments
bwa-mem2 mem -T 10 -k 10 -a ${REF} ${PRIMERS} > alignment.sam

# module load samtools/1.22.1   
samtools view -bS alignment.sam | samtools sort -o alignment.bam

# module load bedtools/2.31.0
bedtools bamtobed -i alignment.bam -cigar > alignment.bed

# Manually check the Bedfile, delete obviously wrong coordinates