#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-15 14:30:38
# @DESCRIPTION:

# Number of input parameters

# ! generate peakcalling --------------------------------------------------------------------

cd /home/liuc9/data/refdata/flare

python /scr1/users/liuc9/tools/FLARE/workflow_peakcaller/scripts/generate_peakcalling_regions.py \
  /mnt/isilon/xing_lab/liuc9/refdata/flare/Homo_sapiens.GRCh38.107.gtf \
  Homo_sapiens.GRCh38.107.gtf_peakcalling_regions
