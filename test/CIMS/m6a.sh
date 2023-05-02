#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-01 22:08:24
# @DESCRIPTION:

# Number of input parameters

# after espresso-internal get the aligned.sorted.sam file
# use samtools fillmd to get the mutation information

# Step1: add MD tag to SAM file
samtools fillmd aligned.sorted.sam /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa | gzip -c >aligned.sorted.md.sam.gz

# Step2: Parsing SAM file
# Parsing SAM file
conda activate ctk

nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file aligned.sorted.mutation.txt aligned.sorted.md.sam.gz aligned.sorted.tag.bed >nohup.step2.out &

wc -l *.tag.bed

# # Remove tags from rRNA and other repetitive RNA (optional)
# perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tagoverlap.pl -big -region /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/annotation/genomes/hg38/annotation/rmsk.RNA.bed -ss --complete-overlap -r --keep-tag-name --keep-score -v aligned.sorted.tag.bed aligned.sorted.tag.norRNA.bed

# Step3: Getting unique CLIP tags by collapsing PCR duplicates
# Collapse PCR duplicates
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2collapse.pl -big -v --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name aligned.sorted.tag.bed aligned.sorted.tag.uniq.bed >nohup.step3.out &

# diagnostic step
awk '{print $3-$2}' aligned.sorted.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' >aligned.sorted.tag.uniq.len.dist.txt

# Get the mutations in unique tags
nohup python /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/joinWrapper.py aligned.sorted.mutation.txt aligned.sorted.tag.uniq.bed 4 4 N aligned.sorted.tag.uniq.mutation.txt >nohup.step3.2.out &

# Step4: Merging biological replicates
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2rgb.pl -v -col "128,0,0" aligned.sorted.tag.uniq.bed aligned.sorted.tag.uniq.rgb.bed >nohup.step4.out &

# Step5: Annotating and visualizing CLIP tags
# Get genomic distribution of CLIP tags
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary Fox.pool.tag.uniq.annot.summary.txt Fox.pool.tag.uniq.rgb.bed Fox.pool.tag.uniq.annot.txt >nohup.step5.out &

# Generate bedgraph for visualization in the genome browser
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2profile.pl -v -big -ss -exact -of bedgraph -n ″Unique Tag Profile″ Fox.pool.tag.uniq.rgb.bed Fox.pool.tag.uniq.bedgraph >nohup.step5.2.out &
