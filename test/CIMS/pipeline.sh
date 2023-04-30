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
