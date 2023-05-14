#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 21:05:51
# @DESCRIPTION:

# Number of input parameters

# ! real data --------------------------------------------------------------------

conda activate annotator

annotator \
  --input /mnt/isilon/xing_lab/liuc9/projdata/m6a/flare/APOBEC1YTH.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed \
  --output /mnt/isilon/xing_lab/liuc9/projdata/m6a/flare/APOBEC1YTH.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno \
  --gtfdb /mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.sqlite3 \
  --format ensembl

# ! https://github.com/byee4/annotator --------------------------------------------------------------------
# /mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.sqlite3
conda activate annotator

annotator \
  --input /scr1/users/liuc9/m6a/test4m.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed \
  --output /scr1/users/liuc9/m6a/test4m.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno \
  --gtfdb /mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.sqlite3 \
  --format ensembl
