#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 21:05:51
# @DESCRIPTION:

# Number of input parameters

# ! https://github.com/byee4/annotator --------------------------------------------------------------------
# /mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.sqlite3
/scr1/users/liuc9/tools/annotator/annotator/annotator.py \
  --input BED6_FILE \
  --output OUTPUT_FILE \
  --gtfdb gencode.v19.annotation.gtf.db \
  --species hg19
