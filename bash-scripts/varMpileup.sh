#!/bin/bash

# Generates an mpileup string for supplied Bams
PATH=$PATH:/home/genmol2/chad/tools/samtools/

args=0
ref=$2
out=" "
if [ $# -eq $args ]
 then
	echo "VCFfile ref.fa Sire-Bam Dam-Bam Child-Bam"
	exit
fi
IFS=$'\n\r'
sites=(`cut -f 1,2 $1 | grep -v "#"`)

for i in "${sites[@]}"
do

target=`echo $i | awk -F"\t" '{print $1 ":" $2 "-" $2}'`
params="samtools mpileup -f $2 -r $target $3 $4 $5 >> ${1}.pileup"
eval $params

done

