#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 13:08:02
# @DESCRIPTION:

# Number of input parameters
param=$#

conda activate clair3

rawsam=/mnt/isilon/xing_lab/liuc9/projdata/m6a/test4m_result/alignment/M6A/0/aligned.sorted.sam
genomeref=/mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa

targetdir=/scr1/users/liuc9/m6a/
prefix=test4m.sorted.md
# Step1: add MD tag to SAM file
samtools fillmd ${rawsam} ${genomeref} >${targetdir}/${prefix}.sam

# Step2: convert sam file to bam and index bam
samtools view -bS --threads 10 ${targetdir}/${prefix}.sam -o ${targetdir}/${prefix}.bam
samtools index -@10 ${targetdir}/${prefix}.bam
