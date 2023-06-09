#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-12 16:01:18
# @DESCRIPTION:

# Number of input parameters

# ! real data --------------------------------------------------------------------

perl /scr1/users/liuc9/tools/annovar/convert2annovar.pl \
  -format vcf4 \
  -filter pass /mnt/isilon/xing_lab/liuc9/projdata/m6a/clair3/merge_output.vcf.gz \
  -outfile /mnt/isilon/xing_lab/liuc9/projdata/m6a/clair3/merge_output.vcf.gz.avinput

perl /scr1/users/liuc9/tools/annovar/table_annovar.pl \
  /mnt/isilon/xing_lab/liuc9/projdata/m6a/clair3/merge_output.vcf.gz.avinput \
  /home/liuc9/data/refdata/annovar/humandb \
  -buildver hg38 \
  -out /mnt/isilon/xing_lab/liuc9/projdata/m6a/clair3/merge_output.vcf.gz.avinput \
  -remove \
  -protocol refGene,ensGene,knownGene,avsnp150 \
  -operation g,g,g,f \
  -nastring .

# ! convert vcf to annovar input --------------------------------------------------------------------

perl /scr1/users/liuc9/tools/annovar/convert2annovar.pl -format vcf4 -filter pass /scr1/users/liuc9/m6a/merge_output.vcf.gz -outfile /scr1/users/liuc9/m6a/merge_output.vcf.avinput

perl /scr1/users/liuc9/tools/annovar/table_annovar.pl \
  /scr1/users/liuc9/m6a/merge_output.vcf.avinput \
  /home/liuc9/data/refdata/annovar/humandb \
  -buildver hg38 \
  -out /home/liuc9/scratch/m6a/merge_output.vcf.avinput \
  -remove \
  -protocol refGene,ensGene,knownGene,avsnp150 \
  -operation g,g,g,f \
  -nastring .
