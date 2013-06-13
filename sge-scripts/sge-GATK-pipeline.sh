#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=3G,mem_free=3G
#$ -pe make 10
#$ -N GATK-pp
#$ -q genmol.q
#$ -t 1-3

jobID=${SGE_TASK_ID}-1
input=(
maraudeur_CGATGT_L001.sort.dedup.bam
5061_AGTCAA_L002001.sort.dedup.bam
5062_AGTTCC_L003001.sort.dedup.bam
)

java -jar /home/genmol2/chad/tools/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -R /home/genmol2/chad/refs/bosTau6.fasta -T  RealignerTargetCreator \
 -nt 4 -o ${input[jobID]}.intervals -I ${input[jobID]} --known /home/genmol2/chad/VCFs/UG.Genome.vcf.gz

java -jar /home/genmol2/chad/tools/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -R /home/genmol2/chad/refs/bosTau6.fasta -T IndelRealigner \
 -targetIntervals ${input[jobID]}.intervals -I ${input[jobID]} --known /home/genmol2/chad/VCFs/UG.Genome.vcf.gz -o ${input[jobID]}.ir.bam -compress 0
 
java -jar /home/genmol2/chad/tools/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -R /home/genmol2/chad/refs/bosTau6.fasta -T BaseRecalibrator \
 -I ${input[jobID]}.ir.bam --known /home/genmol2/chad/VCFs/UG.Genome.vcf.gz -o ${input[jobID]}.grp -nct 10

java -jar /home/genmol2/chad/tools/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -R /home/genmol2/chad/refs/bosTau6.fasta -T PrintReads \
-I ${input[jobID]}.ir.bam -BQSR ${input[jobID]}.grp -nct 10 -compress 6 -o ${input[jobID]}.ir.bqsr.bam

java -jar /home/genmol2/chad/tools/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -R /home/genmol2/chad/refs/bosTau6.fasta -T ReduceReads \
-I ${input[jobID]}.ir.bqsr.bam -o ${input[jobID]}.ir.bqsr.reduced.bam