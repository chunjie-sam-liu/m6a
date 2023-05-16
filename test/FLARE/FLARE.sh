#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 13:46:48
# @DESCRIPTION:

# Number of input parameters
param=$#

# ! real data --------------------------------------------------------------------

# ! peakcalling --------------------------------------------------------------------
# nohup snakemake \
#   --profile /home/liuc9/github/m6a/test/FLARE/profiles/peakcaller/ \
#   --configfile /home/liuc9/github/m6a/test/FLARE/peakcalling_real.json &

nohup snakemake \
  --snakefile /scr1/users/liuc9/tools/FLARE/workflow_peakcaller/Snakefile \
  --configfile /home/liuc9/github/m6a/test/FLARE/peakcalling_real.json \
  --verbose \
  --use-singularity \
  --singularity-args '--bind /scr1/users/liuc9 --bind /mnt/isilon/xing_lab/liuc9 --bind /mnt/isilon/xing_lab' \
  --singularity-prefix /mnt/isilon/xing_lab/liuc9/projdata/m6a/flare_singularity \
  -j8 &

# ! sailor --------------------------------------------------------------------

# nohup snakemake \
#   --profile /home/liuc9/github/m6a/test/FLARE/profiles/sailor/ \
#   --configfile /home/liuc9/github/m6a/test/FLARE/sailor_real.json &

nohup snakemake \
  --snakefile /scr1/users/liuc9/tools/FLARE/workflow_sailor/Snakefile \
  --configfile /home/liuc9/github/m6a/test/FLARE/sailor_real.json \
  --verbose \
  --use-singularity \
  --singularity-args '--bind /scr1/users/liuc9 --bind /mnt/isilon/xing_lab/liuc9 --bind /mnt/isilon/xing_lab' \
  --singularity-prefix /mnt/isilon/xing_lab/liuc9/projdata/m6a/flare_singularity \
  -j1 &

# ! peakcalling --------------------------------------------------------------------

# ! test4m --------------------------------------------------------------------

# ? we include with this release of FLARE an update of SAILOR to enable detection of all edit types.

snakemake \
  --snakefile /scr1/users/liuc9/tools/FLARE/workflow_sailor/Snakefile \
  --configfile /home/liuc9/github/m6a/test/FLARE/sailor.json \
  --verbose \
  --use-singularity \
  --singularity-args '--bind /scr1/users/liuc9 --bind /mnt/isilon/xing_lab/liuc9 --bind /mnt/isilon/xing_lab' \
  -j1 # singularity-prefix /home/liuc9/data/projdata/m6a/FLARE_singularity
