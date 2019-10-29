#!/bin/bash

## submission properties

#PBS -S /bin/bash                       ## shell selection
#PBS -q batch                           ## queue selection
#PBS -N starr-seq_DNA                   ## job name
#PBS -l nodes=1:ppn=10                  ## ppn=threads per node, Needs to match the software argument integer
#PBS -l walltime=48:00:00               ## total time running limit
#PBS -l mem=20gb                        ## memory limit
#PBS -M marand@uga.edu                  ## send to address upon initialization and completion of script
#PBS -m ae                              ##

## change to directory
cd $PBS_O_WORKDIR

## load modules
module load Bowtie/1.2.2-foss-2016b
module load Trimmomatic/0.36-Java-1.8.0_144
module load SAMtools/1.6-foss-2016b

## variables
prefix=starr-high-input-1918_S8_L00
new=dna-starr-trim
output=dna-starr-high-v2

# lane
ln1=1
ln2=2
ln3=3
ln4=4

# read suffix
read1=_R1_001.fastq.gz
read2=_R2_001.fastq.gz

# trimmed read suffix
trim1=_R1_001.trim.fastq
trim2=_R2_001.trim.fastq

# misc
INDEX=/lustre1/apm25309/reference_genomes/Zmays/Zm
threads=10
threads2=5

##-------------##
## trimmomatic ##
##-------------##

## lane1
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads $threads -phred33 $prefix$ln1$read1 $prefix$ln1$read2 $new$ln1$trim1 $new$ln1.R1_UN.fastq $new$ln1$trim2 $new$ln1.R2_UN.fastq ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:3:20 LEADING:0 TRAILING:0 MINLEN:30

## lane2
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads $threads -phred33 $prefix$ln2$read1 $prefix$ln2$read2 $new$ln2$trim1 $new$ln2_R1.UN.fastq $new$ln2$trim2 $new$ln2.R2_UN.fastq ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:3:20 LEADING:0 TRAILING:0 MINLEN:30

## lane 3
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads $threads -phred33 $prefix$ln3$read1 $prefix$ln3$read2 $new$ln3$trim1 $new$ln3.R1_UN.fastq $new$ln3$trim2 $new$ln3.R2_UN.fastq ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:3:20 LEADING:0 TRAILING:0 MINLEN:30

## lane 4
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads $threads -phred33 $prefix$ln4$read1 $prefix$ln4$read2 $new$ln4$trim1 $new$ln4.R1_UN.fastq $new$ln4$trim2 $new$ln4.R2_UN.fastq ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:3:20 LEADING:0 TRAILING:0 MINLEN:30

##----------------##
## bowtie mapping ##
##----------------##

## lane 1
bowtie $INDEX -t -p $threads2 -v 1 --best --strata -m 1 -X 2000 -S -1 $new$ln1$trim1 -2 $new$ln1$trim2 | \
	samtools view -@ $threads2 -bSh -q 1 -F 4 - > $new$ln1.bam

## lane 2
bowtie $INDEX -t -p $threads2 -v 1 --best --strata -m 1 -X 2000 -S -1 $new$ln2$trim1 -2 $new$ln2$trim2 | \
	samtools view -@ $threads2 -bSh -q 1 -F 4 - > $new$ln2.bam

## lane 3
bowtie $INDEX -t -p $threads2 -v 1 --best --strata -m 1 -X 2000 -S -1 $new$ln3$trim1 -2 $new$ln3$trim2 | \
	samtools view -@ $threads2 -bSh -q 1 -F 4 - > $new$ln3.bam

## lane 4
bowtie $INDEX -t -p $threads2 -v 1 --best --strata -m 1 -X 2000 -S -1 $new$ln4$trim1 -2 $new$ln4$trim2 | \
	samtools view -@ $threads2 -bSh -q 1 -F 4 - > $new$ln4.bam

## merge and sort
samtools merge -@ $threads2 $new.merged.bam $new$ln1.bam $new$ln2.bam $new$ln3.bam $new$ln4.bam 

## clean up
rm *.fastq
