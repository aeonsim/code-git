#!/bin/bash

files=(`ls *bam`)

for i in "${files[@]}"
do

samtools index $i

done