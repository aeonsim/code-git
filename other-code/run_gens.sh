#!/bin/bash

CNT=0

scala -cp commons-math3-3.0.jar gens ped.txt 770K-genotypes.vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt

CNT=$((CNT+1)) 

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
CNT=$((CNT+1))

scala -cp commons-math3-3.0.jar gens ped$((CNT)).txt vcfout$((CNT)).vcf ped$((CNT+1)).txt vcfout$((CNT+1)).vcf mating$((CNT)).txt
