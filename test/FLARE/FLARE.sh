#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 13:46:48
# @DESCRIPTION:

# Number of input parameters
param=$#
# ? we include with this release of FLARE an update of SAILOR to enable detection of all edit types.

snakemake \
  --snakefile /scr1/users/liuc9/tools/FLARE/workflow_sailor/Snakefile \
  --configfile /home/liuc9/github/m6a/test/FLARE/sailor.json \
  --verbose \
  --use-singularity \
  \
  --singularity-args '--bind /scr1/users/liuc9 --bind /mnt/isilon/xing_lab/liuc9 --bind /mnt/isilon/xing_lab' \
  -j1 # --singularity-args '--bind /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/ --bind /scr1/users/liuc9/m6a --bind /scr1/users/liuc9/m6a --bind /mnt/isilon/xing_lab/liuc9/refdata/dbsnp/knownsnp/' \
