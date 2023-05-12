#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 12:57:28
# @DESCRIPTION:

# Number of input parameters
# param=$#

conda activate clair3

MODEL_NAME="r941_prom_sup_g5014"
THREADS=30

genomeref=/mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa
genomedir=/mnt/isilon/xing_lab/liuc9/refdata/bwa_index/
targetdir=/scr1/users/liuc9/m6a/
prefix=test4m.sorted.md

# /home/liuc9/tools/Clair3/run_clair3.sh \
#   --bam_fn=${targetdir}/${prefix}.bam \
#   --ref_fn=${genomeref} \
#   --threads=${THREADS} \
#   --platform="ont" \
#   --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \
#   --output=${targetdir}

singularity exec \
  -B ${genomedir},${targetdir} \
  /home/liuc9/sif/clair3_latest.sif \
  /opt/bin/run_clair3.sh \
  --bam_fn=${targetdir}/${prefix}.bam \
  --ref_fn=${genomedir}/Homo_sapiens.GRCh38.104.fa \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${targetdir}
