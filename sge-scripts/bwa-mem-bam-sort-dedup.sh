#$ -S /bin/bash
## Scripting Language to use
#$ -cwd
## run in directory where it was submitted
#$ -l h_vmem=1G,mem_free=1G
## Amount of memory to request, is using pe total memory = pe slots x mem request (so 20x1 = 20GB here)
#$ -pe make 20
## the number of Threads needed you don't need this line if your program is single threaded
## instead comment it out with ## or delete it.
#$ -N BWA-MEM Mapping
## The Name of the Job
#$ -q genmol.q
## The Queue you wish to use
#$ -r y
## Automatically rerun the job if the compute node is rebooted during the jobs run time
#$ -t 1-31
##The number of subjobs you need to use


## To map your data create a file called input-files.txt in the directory
## containing a list of the files that Make up one read of each Paired end dataset via the command
## ls *_R1_*.fastq

files=($(< input-files.txt))
## Creates an Array of Files to process
threads=$NSLOTS
## sets threads equal to the number of Threads the SGE has assigned
curJob=${SGE_TASK_ID}-1
## Creating a curJob variable one less than the SGE subjob id as Bash Arrays start at 0 but
## SGE Sub job IDs start at 1.

## The Script you want to run goes on from here:

## note $files is the array of filenames (full path) I want to work with.
## ${files[1]} requests the filename in the files array as position 1.

read2=`echo ${files[$curJob]} | sed -r "s/_R1_/_R2_/"`
rgheader=`echo ${files[$curJob]} | sed -r "s/\/.*\///" | sed -r "s/_R1_.*bz2//"`

out=`echo ${files[$curJob]} | sed -r "s/\/.*\///" | sed -r "s/_R1_//"`
name=`echo ${files[$curJob]} | sed -r "s/\/.*\///" | sed -r "s/_.*bz2//"`

header=`echo "@RG\\tID:${rgheader}\\tSM:${name}"`

## Creating Variables containing the strings I need for my command
## This script works by loading a list of files from input-files.txt, determining the pair for that file
## Extracting the Name and other details that are needed from the filename then processing the data
## to generate a sorted, PCR deduplicated BAM file for each pair of Fastq files, using samtools & BWA MEM



## Only needed for your FASTQ is Bzip2 compressed if Gzipped then can use directly
in1=`echo "<bzip2 -dc ${files[$curJob]}"`
in2=`echo "<bzip2 -dc ${read2}"`

#Bzip2 Version
echo Mapping $name to ${out}.sam
/home/genmol2/chad/tools/bwa/bwa mem -M -R $header -t $threads /home/genmol2/chad/refs/bosTau6.lic.fa "$in1" "$in2" > ${out}.sam

## Gzip version
## /home/genmol2/chad/tools/bwa/bwa mem -M -R $header -t $threads /home/genmol2/chad/refs/bosTau6.lic.fa ${files[$curJob]} $read2 > ${out}.sam

echo Mapping Complete Sorting, Starting Bam conversion
/home/genmol2/chad/tools/samtools/samtools view -buS ${out}.sam | /home/genmol2/chad/tools/samtools/samtools sort -m 500000000 - ${out}
rm ${out}.sam
echo Starting Dedup
/home/genmol2/chad/tools/samtools/samtools rmdup ${out}.bam ${out}-dedup.bam
rm ${out}.bam
echo Indexing
/home/genmol2/chad/tools/samtools/samtools index ${out}-dedup.bam &

# You will need to merge the individual BAM files together before processing with GATK.


## Below is a non SGE version of the Same script using a Bash for Loop instead of SGEs sub jobs system.

: <<NONSGEVERSION

for i in "${files[@]}"
do
read2=`echo $i | sed -r "s/_R1_/_R2_/"`
rgheader=`echo $i | sed -r "s/\/.*\///" | sed -r "s/_R1_.*bz2//"`
out=`echo $i | sed -r "s/\/.*\///" | sed -r "s/_R1_//"`
name=`echo $i | sed -r "s/\/.*\///" | sed -r "s/_.*bz2//"`
header=`echo "@RG\\tID:${rgheader}\\tSM:${name}"`
in1=`echo "<bzip2 -dc ${i}"`
in2=`echo "<bzip2 -dc ${read2}"`
echo Mapping $name to ${out}.sam
/home/genmol2/chad/tools/bwa/bwa mem -M -R $header -t $threads /home/genmol2/chad/refs/bosTau6.lic.fa "$in1" "$in2" > ${out}.sam
echo Mapping Complete Sorting, Starting Bam conversion
/home/genmol2/chad/tools/samtools/samtools view -buS ${out}.sam | /home/genmol2/chad/tools/samtools/samtools sort -m 500000000 - ${out}
rm ${out}.sam
echo Starting Dedup
/home/genmol2/chad/tools/samtools/samtools rmdup ${out}.bam ${out}-dedup.bam
rm ${out}.bam
echo Indexing
/home/genmol2/chad/tools/samtools/samtools index ${out}-dedup.bam &

#java -jar /home/genmol2/chad/tools/picard-tools-1.79/SortSam.jar VALIDATION_STRINGENCY=SILENT I=${out}.sam O=${out}.bam QUIET=true SO=coordinate COMPRESSION_LEVEL=0 CREATE_INDEX=true 
#MAX_RECORDS_IN_RAM=50000000 ; \
#java -jar /home/genmol2/chad/tools/picard-tools-1.79/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT I=${out}.bam O=${out}-dedup.bam COMPRESSION_LEVEL=6 CREATE_INDEX=true M=${out}.metrics 
#REMOVE_DUPLICATES=true  &
done

NONSGEVERSION
