#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-04-29 19:21:54
# @DESCRIPTION:

# Number of input parameters
param=$#

prefetch --max-size 50G SRR3147674 --output-directory /scr1/users/liuc9/iclip &
prefetch --max-size 50G SRR3147675 --output-directory /scr1/users/liuc9/iclip &
prefetch --max-size 50G SRR3147676 --output-directory /scr1/users/liuc9/iclip &
prefetch --max-size 50G SRR3147677 --output-directory /scr1/users/liuc9/iclip &

wait

/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147674.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &
/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147675.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &
/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147676.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &
/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump SRR3147677.sra --temp /scr1/users/liuc9/tmp/fasterq_dump --mem 50G --threads 10 &

#!/bin/bash

mkdir filtering
cd filtering
# Trimming of 3' linker sequences
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  cutadapt --times 1 -e 0 -O 1 --quality-cutoff 5 -m 24 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -o $f.trim.fastq.gz ../fastq/$f.fastq >$f.cutadpt.log
done

# Note 1: It is good to check the number of reads by running the following command:
#e.g., for raw reads
for f in $(ls ../fastq/*.fastq); do
  c=$(cat $f | wc -l)
  c=$((c / 4))
  echo $f $c
done

#e.g., for trimmed reads
for f in $(ls *.trim.fastq.gz); do
  c=$(cat $f | wc -l)
  c=$((c / 4))
  echo $f $c
done

# Collapse exact duplicates

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/fastq2collapse.pl $f.trim.fastq.gz - | gzip -c >$f.trim.c.fastq.gz &
done

for f in $(ls *.trim.c.fastq.gz); do
  {
    c=$(zcat $f | wc -l)
    c=$((c / 4))
    echo $f $c
  } &
done

# Strip random barcode (UMI)

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl ~/tools/anaconda3/envs/ctk/lib/ctk/stripBarcode.pl -format fastq -len 9 $f.trim.c.fastq.gz - | gzip -c >$f.trim.c.tag.fastq.gz &
done

#Get the distribution of tag lengths for diagnostic purposes using the following command.
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  zcat $f.trim.c.tag.fastq.gz | awk '{if(NR%4==2) {print length($0)}}' | sort -n | uniq -c | awk '{print $2"\t"$1}' >$f.trim.c.tag.seqlen.stat.txt &
done

# Read mapping & parsing
cd mapping
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  {
    bwa aln -t 10 -n 0.06 -q 20 /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa ../filtering/$f.trim.c.tag.fastq.gz >$f.sai
    bwa samse /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa $f.sai ../filtering/$f.trim.c.tag.fastq.gz | gzip -c >$f.sam.gz
  } &
done

# Parsing SAM file
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file $f.mutation.txt $f.sam.gz $f.tag.bed &
done

# samtools view -bS $f.sam | samtools sort - $f.sorted
# samtools fillmd $f.sorted.bam /genomes/hg19/bwa/hg19.fa | gzip -c >$f.sorted.md.sam.gz

# Remove tags from rRNA and other repetitive RNA (optional)
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tagoverlap.pl -big -region /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/annotation/genomes/hg38/annotation/rmsk.RNA.bed -ss --complete-overlap -r --keep-tag-name --keep-score -v $f.tag.bed $f.tag.norRNA.bed &
done

# Collapse PCR duplicates
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2collapse.pl -big -v --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name $f.tag.norRNA.bed $f.tag.uniq.bed &
done

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  awk '{print $3-$2}' $f.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' >$f.tag.uniq.len.dist.txt
done

for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  python /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/joinWrapper.py $f.mutation.txt $f.tag.uniq.bed 4 4 N $f.tag.uniq.mutation.txt &
done

# Merging biological replicates
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2rgb.pl -v -col "188,0,0" RBFOX2.rep1.tag.uniq.bed RBFOX2.rep1.tag.uniq.rgb.bed &
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2rgb.pl -v -col "128,0,128" RBFOX2.rep2.tag.uniq.bed RBFOX2.rep2.tag.uniq.rgb.bed &
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2rgb.pl -v -col "0,128,0" RBFOX2.rep3.tag.uniq.bed RBFOX2.rep3.tag.uniq.rgb.bed &
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2rgb.pl -v -col "51,102,255" RBFOX2.rep4.tag.uniq.bed RBFOX2.rep4.tag.uniq.rgb.bed &

cat RBFOX2.rep1.tag.uniq.rgb.bed RBFOX2.rep2.tag.uniq.rgb.bed RBFOX2.rep3.tag.uniq.rgb.bed RBFOX2.rep4.tag.uniq.rgb.bed >RBFOX2.pool.tag.uniq.rgb.bed
cat RBFOX2.rep1.tag.uniq.mutation.txt RBFOX2.rep2.tag.uniq.mutation.txt RBFOX2.rep3.tag.uniq.mutation.txt RBFOX2.rep4.tag.uniq.mutation.txt >RBFOX2.pool.tag.uniq.mutation.txt

# Annotating and visualizing CLIP tags
# Get the genomic distribution of CLIP tags
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary RBFOX2.pool.tag.uniq.annot.summary.txt RBFOX2.pool.tag.uniq.rgb.bed RBFOX2.pool.tag.uniq.annot.txt

# Generate bedgraph for visualization in the genome browser
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2profile.pl -v -ss -exact -of bedgraph -n ″Unique Tag Profile″ RBFOX2.pool.tag.noRNA.uniq.rgb.bed RBFOX2.pool.tag.uniq.bedgraph

# Peak calling
mkdir cluster
cd cluster
ln -s ../mapping/RBFOX2.pool.tag.uniq.rgb.bed
# Mode 1: Peak calling with no statistical significance
perl /usr/local/CTK/tag2peak.pl -big -ss -v --valley-seeking --valley-depth 0.9 RBFOX2.pool.tag.uniq.rgb.bed RBFOX2.pool.tag.uniq.peak.bed --out-boundary RBFOX2.pool.tag.uniq.peak.boundary.bed --out-half-PH RBFOX2.pool.tag.uniq.peak.halfPH.bed
perl /usr/local/CTK/bed2annotation.pl -dbkey hg19 -ss -big -region -v -summary RBFOX2.pool.tag.uniq.peak.annot.summary.txt RBFOX2.pool.tag.uniq.peak.bed RBFOX2.pool.tag.uniq.peak.annot.txt
# Mode 2: Peak calling with statistical significance
perl /usr/local/CTK/tag2peak.pl -big -ss -v --valley-seeking -p 0.05 --valley-depth 0.9 --multi-test --dbkey hg19 RBFOX2.pool.tag.uniq.rgb.bed RBFOX2.pool.tag.uniq.peak.sig.bed --out-boundary RBFOX2.pool.tag.uniq.peak.sig.boundary.bed --out-half-PH RBFOX2.pool.tag.uniq.peak.sig.halfPH.bed
# To get the number of significant peaks:
wc -l RBFOX2.pool.tag.noRNA.uniq.peak.sig.bed
# 35602 RBFOX2.pool.tag.uniq.peak.sig.bed
#extract from -500 to +500
awk '{print $1"\t"int(($2+$3)/2)-500"\t"int(($2+$3)/2)+500"\t"$4"\t"$5"\t"$6}' RBFOX2.pool.tag.uniq.peak.sig.bed >RBFOX2.pool.tag.uniq.peak.sig.center.ext1k.bed
for f in RBFOX2.rep1 RBFOX2.rep2 RBFOX2.rep3 RBFOX2.rep4; do
  perl /usr/local/CTK/tag2profile.pl -ss -region RBFOX2.pool.tag.uniq.peak.sig.boundary.bed -of bed -v $f.tag.uniq.bed $f.tag.uniq.peak.sig.boundary.count.bed
done

# CIMS
mkdir CIMS
cd CIMS

ln -s ../mapping/RBFOX2.pool.tag.uniq.rgb.bed
ln -s ../mapping/RBFOX2.pool.tag.uniq.mutation.txt
# Get specific types of mutations
perl ~/czlab_src/CTK/getMutationType.pl -t del RBFOX2.pool.tag.uniq.mutation.txt RBFOX2.pool.tag.uniq.del.bed
#deletions
perl ~/czlab_src/CTK/getMutationType.pl -t ins RBFOX2.pool.tag.uniq.mutation.txt RBFOX2.pool.tag.uniq.ins.bed
#insertions
perl ~/czlab_src/CTK/getMutationType.pl -t sub --summary RBFOX2.pool.tag.uniq.mutation.stat.txt RBFOX2.pool.tag.uniq.mutation.txt RBFOX2.pool.tag.uniq.sub.bed
#substitutions, as well as summary statistics of different types of mutations
# Get CIMS
perl /usr/local/CTK/CIMS.pl -big -n 10 -p -outp RBFOX2.pool.tag.uniq.del.pos.stat.txt -v RBFOX2.pool.tag.uniq.rgb.bed RBFOX2.pool.tag.uniq.del.bed RBFOX2.pool.tag.uniq.del.CIMS.txt
awk '{if($9<=0.001) {print $0}}' RBFOX2.pool.tag.uniq.del.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n >RBFOX2.pool.tag.uniq.del.CIMS.s30.txt
cut -f 1-6 RBFOX2.pool.tag.uniq.del.CIMS.s30.txt >RBFOX2.pool.tag.uniq.del.CIMS.s30.bed
wc -l RBFOX2.pool.tag.uniq.del.CIMS.s30.bed
awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' RBFOX2.pool.tag.uniq.del.CIMS.s30.bed >RBFOX2.pool.tag.uniq.del.CIMS.s30.21nt.bed
#+/-10 around CIMS
