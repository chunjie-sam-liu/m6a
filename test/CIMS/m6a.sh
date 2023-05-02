#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-01 22:08:24
# @DESCRIPTION:

# Number of input parameters
param=$#

# Parsing SAM file
perl /scr1/users/liuc9/tools/anaconda3/envs/ctk/lib/ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file aligned.sorted.mutation.txt aligned.sorted.sam aligned.sorted.tag.bed
