#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 14:38:26
# @DESCRIPTION:

# Number of input parameters
param=$#

chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
outputdir=/home/liuc9/data/refdata/dbsnp/knownsnp

# ! download from known snp --------------------------------------------------------------------

for chrom in ${chroms[@]}; do
  echo "Downloading known SNPs for chromosome ${chrom}..."
  cmd="wget -O ${outputdir}/bed_chr_${chrom}.bed.gz https://ftp.ncbi.nih.gov/snp/organisms/human_9606/BED/bed_chr_${chrom}.bed.gz"
  # echo ${cmd}
  # eval ${cmd} &
done

# ! gunzip --------------------------------------------------------------------
# * gunzip --------------------------------------------------------------------
# ? gunzip --------------------------------------------------------------------
# TODO gunzip --------------------------------------------------------------------

for chrom in ${chroms[@]}; do
  echo "Unzipping known SNPs for chromosome ${chrom}..."
  cmd="gunzip ${outputdir}/bed_chr_${chrom}.bed.gz"
  # echo ${cmd}
  # eval ${cmd} &
done

# ! merge --------------------------------------------------------------------

for chrom in ${chroms[@]}; do
  b=${outputdir}/bed_chr_${chrom}.bed
  echo $b
  # tail -n+2 $b | cut -f1,2,3 | sort >>human_dbsnp_combined.bed3
  awk 'NR > 1 {sub(/^chr/, "", $1); print $1"\t"$2"\t"$3}' $b | sort -k2n >>${outputdir}/human_dbsnp_combined.bed3
done
