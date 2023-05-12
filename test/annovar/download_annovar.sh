#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 15:17:37
# @DESCRIPTION:

# Number of input parameters

# ! available databases --------------------------------------------------------------------

# perl /scr1/users/liuc9/tools/annovar/annotate_variation.pl -webfrom annovar -downdb avdblist -buildver hg38 /home/liuc9/data/refdata/annovar/humandb

# ! donwload selected databases --------------------------------------------------------------------

download='/scr1/users/liuc9/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar'
outputdir='/home/liuc9/data/refdata/annovar/humandb'
# dbs=('refGene' 'knownGene' 'ensGene' 'avsift' 'ljb26_all' 'dbnsfp42c' 'dbscsnv11' 'intervar_20180118' 'cosmic70' 'esp6500siv2_ea' 'esp6500siv2_aa' 'esp6500siv2_all' 'exac03' 'gnomad_exome' 'gnomad_genome' 'kaviar_20150923' 'hrcr1' '1000g2015aug' 'gme' 'mcap' 'revel' 'avsnp150' 'nci60' 'clinvar_20180603' 'regsnpintron')
dbs=("1000g2015aug" "avsnp150" "clinvar_20221231" "cosmic70" "dbnsfp42c" "dbscsnv11" "ensGene" "ensGeneMrna" "esp6500siv2_all" "exac03" "gene4denovo201907" "gnomad312_genome" "hrcr1" "icgc28" "knownGene" "knownGeneMrna" "ljb26_all" "nci60" "refGene" "refGeneMrna")

for db in ${dbs[@]}; do
  echo "Downloading ${db}..."
  cmd="${download} ${db} ${outputdir}"
  echo ${cmd}
  eval ${cmd} &
done
