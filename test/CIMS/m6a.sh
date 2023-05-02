#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-01 22:08:24
# @DESCRIPTION:

# Number of input parameters

# after espresso-internal get the aligned.sorted.sam file
# use samtools fillmd to get the mutation information

# Step1: add MD tag to SAM file
samtools fillmd aligned.sorted.sam /mnt/isilon/xing_lab/liuc9/refdata/bwa_index/Homo_sapiens.GRCh38.104.fa | gzip -c >aligned.sorted.md.sam.gz

# Step2: Parsing SAM file
# Parsing SAM file
conda activate ctk

nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file aligned.sorted.mutation.txt aligned.sorted.md.sam.gz aligned.sorted.tag.bed >nohup.step2.out &

wc -l *.tag.bed

# # Remove tags from rRNA and other repetitive RNA (optional)
# perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tagoverlap.pl -big -region /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/annotation/genomes/hg38/annotation/rmsk.RNA.bed -ss --complete-overlap -r --keep-tag-name --keep-score -v aligned.sorted.tag.bed aligned.sorted.tag.norRNA.bed

# Step3: Getting unique CLIP tags by collapsing PCR duplicates
# Collapse PCR duplicates
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2collapse.pl -big -v --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name aligned.sorted.tag.bed aligned.sorted.tag.uniq.bed >nohup.step3.out &

# diagnostic step
awk '{print $3-$2}' aligned.sorted.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' >aligned.sorted.tag.uniq.len.dist.txt

# Get the mutations in unique tags
nohup python /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/joinWrapper.py aligned.sorted.mutation.txt aligned.sorted.tag.uniq.bed 4 4 N aligned.sorted.tag.uniq.mutation.txt >nohup.step3.2.out &

# Step4: Merging biological replicates
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2rgb.pl -v -col "128,0,0" aligned.sorted.tag.uniq.bed aligned.sorted.tag.uniq.rgb.bed >nohup.step4.out &

# Step5: Annotating and visualizing CLIP tags
# Get genomic distribution of CLIP tags
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary aligned.sorted.tag.uniq.annot.summary.txt aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.annot.txt >nohup.step5.out &

# Generate bedgraph for visualization in the genome browser
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2profile.pl -v -big -ss -exact -of bedgraph -n ″Unique Tag Profile″ aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.bedgraph >nohup.step5.2.out &

# Step6 Peak calling
mkdir cluster
cd cluster
ln -s ../mapping/aligned.sorted.tag.uniq.rgb.bed .

# Peak calling without statistical assessment
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2peak.pl -big -ss -v --valley-seeking --valley-depth 0.9 aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.peak.bed --out-boundary aligned.sorted.tag.uniq.peak.boundary.bed --out-half-PH aligned.sorted.tag.uniq.peak.halfPH.bed >nohup.step6.out &

nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary aligned.sorted.tag.uniq.peak.annot.summary.txt aligned.sorted.tag.uniq.peak.bed aligned.sorted.tag.uniq.peak.annot.txt >nohup.step6.2.out &

# Mode 2: Peak calling with statistical assessment
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/tag2peak.pl -big -ss -v --valley-seeking -p 0.05 --valley-depth 0.9 --dbkey hg38 --multi-test aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.peak.sig.bed --out-boundary aligned.sorted.tag.uniq.peak.sig.boundary.bed --out-half-PH aligned.sorted.tag.uniq.peak.sig.halfPH.bed >nohup.step6.3.out &

# Step7: CIMS analysis
mkdir CIMS
cd CIMS
ln -s ../mapping/aligned.sorted.tag.uniq.rgb.bed
ln -s ../mapping/aligned.sorted.tag.uniq.mutation.txt
# Get specific types of mutations
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/getMutationType.pl -t del aligned.sorted.tag.uniq.mutation.txt aligned.sorted.tag.uniq.del.bed >nohup.step7.out &
#deletions
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/getMutationType.pl -t ins aligned.sorted.tag.uniq.mutation.txt aligned.sorted.tag.uniq.ins.bed >nohup.step7.2.out &
#insertions
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/getMutationType.pl -t sub --summary aligned.sorted.tag.uniq.mutation.stat.txt aligned.sorted.tag.uniq.mutation.txt aligned.sorted.tag.uniq.sub.bed >nohup.step7.3.out &
#substitutions, as well as summary statistics of different types of mutations

# Get CIMS

#del example
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/CIMS.pl -big -n 10 -p -outp aligned.sorted.tag.uniq.del.pos.stat.txt -v -c cache_del aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.del.bed aligned.sorted.tag.uniq.del.CIMS.txt >nohup.step7.4.out &

awk '{if($9<=0.001) {print $0}}' aligned.sorted.tag.uniq.del.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n >aligned.sorted.tag.uniq.del.CIMS.s30.txt
cut -f 1-6 aligned.sorted.tag.uniq.del.CIMS.s30.txt >aligned.sorted.tag.uniq.del.CIMS.s30.bed
wc -l aligned.sorted.tag.uniq.del.CIMS.s30.bed
awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' aligned.sorted.tag.uniq.del.CIMS.s30.bed >aligned.sorted.tag.uniq.del.CIMS.s30.21nt.bed
#+/-10 around CIMS
# ins example
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/CIMS.pl -big -n 10 -p -outp aligned.sorted.tag.uniq.ins.pos.stat.txt -v -c cache_ins aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.ins.bed aligned.sorted.tag.uniq.ins.CIMS.txt >nohup.step7.5.out &

awk '{if($9<=0.001) {print $0}}' aligned.sorted.tag.uniq.ins.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n >aligned.sorted.tag.uniq.ins.CIMS.s30.txt
cut -f 1-6 aligned.sorted.tag.uniq.ins.CIMS.s30.txt >aligned.sorted.tag.uniq.ins.CIMS.s30.bed
wc -l aligned.sorted.tag.uniq.ins.CIMS.s30.bed
awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' aligned.sorted.tag.uniq.ins.CIMS.s30.bed >aligned.sorted.tag.uniq.ins.CIMS.s30.21nt.bed
#+/-10 around CIMS
#sub example
nohup perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/CIMS.pl -big -n 10 -p -outp aligned.sorted.tag.uniq.sub.pos.stat.txt -v -c cache_sub aligned.sorted.tag.uniq.rgb.bed aligned.sorted.tag.uniq.sub.bed aligned.sorted.tag.uniq.sub.CIMS.txt >nohup.step7.6.out &

awk '{if($9<=0.001) {print $0}}' aligned.sorted.tag.uniq.sub.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n >aligned.sorted.tag.uniq.sub.CIMS.s30.txt
cut -f 1-6 aligned.sorted.tag.uniq.sub.CIMS.s30.txt >aligned.sorted.tag.uniq.sub.CIMS.s30.bed
wc -l aligned.sorted.tag.uniq.sub.CIMS.s30.bed
awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' aligned.sorted.tag.uniq.sub.CIMS.s30.bed >aligned.sorted.tag.uniq.sub.CIMS.s30.21nt.bed
