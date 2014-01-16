#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=24 --mem-per-cpu=2G 
###--cpus-per-task=24
#SBATCH --mail-type=END
#SBATCH --array=0-8

SAMTOOLS=/home/aeonsim/scripts/apps-damona-Oct13/samtools/samtools
BWA=/home/aeonsim/scripts/apps-damona-Oct13/bwa/bwa
REF=/home/aeonsim/refs/bosTau6.fasta

echo "ARRAY JOB: ${SLURM_ARRAY_TASK_ID}"

export FLOW=`echo $1 | awk '{split($0,arra,"/"); print arra[7]}'`
echo "FLOWCELL IS: ${FLOW}"

R1=(`ls $1*R1*`)
R2=(`ls $1*R2*`)

echo ${R1[@]}

NAME=`echo ${R1[$SLURM_ARRAY_TASK_ID]} | awk '{split($0,arra,"/"); split(arra[8],brra,"_"); print brra[1]}'`

echo "SAMPLE IS: ${NAME}"
RG="@RG\tID:${NAME}_${FLOW}\tPL:ILLUMINA\tPU:${FLOW}\tSM:${NAME}"
echo "BAM HEADER IS: ${RG}"

echo $SLURM_ARRAY_TASK_ID

echo "$BWA mem -t $SLURM_JOB_CPUS_PER_NODE -M -R ${RG} ${REF} ${R1[$SLURM_ARRAY_TASK_ID]} ${R2[$SLURM_ARRAY_TASK_ID]} | $SAMTOOLS view -bS - > /scratch/aeonsim/${NAME}_${FLOW}.bam"

$BWA mem -t $SLURM_JOB_CPUS_PER_NODE -M -R ${RG} ${REF} ${R1[$SLURM_ARRAY_TASK_ID]} ${R2[$SLURM_ARRAY_TASK_ID]} | $SAMTOOLS view -bS - > /scratch/aeonsim/${NAME}_${FLOW}.bam


echo "$SAMTOOLS sort -@ $SLURM_JOB_CPUS_PER_NODE -m 1800M /scratch/aeonsim/${NAME}_${FLOW}.bam ${NAME}_${FLOW}_sorted"

$SAMTOOLS sort -@ $SLURM_JOB_CPUS_PER_NODE -m 1800M /scratch/aeonsim/${NAME}_${FLOW}.bam ${NAME}_${FLOW}_sorted
