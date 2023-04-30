#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-04-29 19:21:54
# @DESCRIPTION:

# Number of input parameters
param=$#

prefetch --max-size 50G SRR3147674 --output-directory /scr1/users/liuc9/iclip &
prefetch --max-size 50G SRR3147675 --output-directory /scr1/users/liuc9/iclip &
prefetch --max-size 50G SRR3147676 --output-directory /scr1/users/liuc9/iclip &
prefetch --max-size 50G SRR3147677 --output-directory /scr1/users/liuc9/iclip &

wait

/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147674.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &
/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147675.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &
/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147676.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &
/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147677.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &

#!/bin/bash

mkdir filtering
cd filtering
# Trimming of 3' linker sequences
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  cutadapt --times 1 -e 0 -O 1 --quality-cutoff 5 -m 24 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -o $f.trim.fastq.gz ../fastq/$f.fastq >$f.cutadpt.log
done

# Note 1: It is good to check the number of reads by running the following command:
#e.g., for raw reads
for f in $(ls ../fastq/*.fastq); do
  c=$(cat $f | wc -l)
  c=$((c / 4))
  echo $f $c
done

#e.g., for trimmed reads
for f in $(ls *.trim.fastq.gz); do
  c=$(cat $f | wc -l)
  c=$((c / 4))
  echo $f $c
done

# Collapse exact duplicates

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/fastq2collapse.pl $f.trim.fastq.gz - | gzip -c >$f.trim.c.fastq.gz &
done

for f in $(ls *.trim.c.fastq.gz); do
  {
    c=$(zcat $f | wc -l)
    c=$((c / 4))
    echo $f $c
  } &
done

# Strip random barcode (UMI)

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl ~/tools/anaconda3/envs/ctk/lib/ctk/stripBarcode.pl -format fastq -len 9 $f.trim.c.fastq.gz - | gzip -c >$f.trim.c.tag.fastq.gz &
done

#Get the distribution of tag lengths for diagnostic purposes using the following command.
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  zcat $f.trim.c.tag.fastq.gz | awk '{if(NR%4==2) {print length($0)}}' | sort -n | uniq -c | awk '{print $2"\t"$1}' >$f.trim.c.tag.seqlen.stat.txt &
done

# Read mapping & parsing
cd mapping
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  {
    bwa aln -t 10 -n 0.06 -q 20 /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa ../filtering/$f.trim.c.tag.fastq.gz >$f.sai
    bwa samse /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa $f.sai ../filtering/$f.trim.c.tag.fastq.gz | gzip -c >$f.sam.gz
  } &
done

# Parsing SAM file
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file $f.mutation.txt $f.sam.gz $f.tag.bed &
done

samtools view -bS $f.sam | samtools sort - $f.sorted
samtools fillmd $f.sorted.bam /genomes/hg19/bwa/hg19.fa | gzip -c >$f.sorted.md.sam.gz

# Collapse PCR duplicates

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2collapse.pl -big -v --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name $f.tag.norRNA.bed $f.tag.uniq.bed &
done

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  awk '{print $3-$2}' $f.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' >$f.tag.uniq.len.dist.txt &
done

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  python /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/joinWrapper.py $f.mutation.txt $f.tag.uniq.bed 4 4 N $f.tag.uniq.mutation.txt
done
