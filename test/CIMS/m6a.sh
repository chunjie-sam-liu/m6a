#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-01 22:08:24
# @DESCRIPTION:

# Number of input parameters
param=$#

# Parsing SAM file
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file aligned.sorted.mutation.txt aligned.sorted.sam aligned.sorted.tag.bed

# Remove tags from rRNA and other repetitive RNA (optional)
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tagoverlap.pl -big -region /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/annotation/genomes/hg38/annotation/rmsk.RNA.bed -ss --complete-overlap -r --keep-tag-name --keep-score -v aligned.sorted.tag.bed aligned.sorted.tag.norRNA.bed

# Collapse PCR duplicates
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2collapse.pl -big -v --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name aligned.sorted.tag.norRNA.bed aligned.sorted.tag.uniq.bed

awk '{print $3-$2}' aligned.sorted.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' >aligned.sorted.tag.uniq.len.dist.txt

python /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/joinWrapper.py aligned.sorted.mutation.txt aligned.sorted.tag.uniq.bed 4 4 N aligned.sorted.tag.uniq.mutation.txt

# Merging biological replicates
