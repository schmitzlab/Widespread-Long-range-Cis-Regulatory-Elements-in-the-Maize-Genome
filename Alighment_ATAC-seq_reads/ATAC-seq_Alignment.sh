###required softwares
#bedtools/2.26.0
#bowtie/1.2.2 
#samtools/1.3.1
#UCSC tools version 359
#picard 2.16.0
#Trimmomatic 0.36
#Java 1.8.0_144

thread=             ####number of CPU thread
input=              ####reads prefix
name=               ####output file name
INDEX=              ####bowtie1 INDEX
chrom_info=         ####chromosome information


echo "trimmomatic  version 0.32" >&2
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads $thread -phred33 $input.R1.fastq.gz $input.R2.fastq.gz ${name}.L.trim.fastq ${name}.L.trimU.fastq ${name}.R.trim.fastq ${name}.R.trimU.fastq ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:3:20 LEADING:0 TRAILING:0 MINLEN:30
# bowtie mapping
echo "bowtie   version 12.8" >&2
bowtie $INDEX -t -p 4 -v 2 --best --strata -m 1 -X 1000 -S ${name}.sam -1 ${name}.L.trim.fastq -2 ${name}.R.trim.fastq

# sort sam to bam
echo "sam to bam   samtools 1.3.1" >&2
echo "sort  samtools 1.3.1" >&2
samtools sort -O 'bam' -o ${name}.sorted.bam -T tmp ${name}.sam

# remove clonal
echo "remove clonal picard 2.16.0 " >&2
java -Xmx20g -classpath /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144 -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates INPUT=${name}.sorted.bam OUTPUT=${name}.clean.bam METRICS_FILE=XXX.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

# bam to bed
echo "bam to bed  bedtools 2.26.0" >&2
bedtools bamtobed -i ${name}.clean.bam > ${name}.bed
# Genome_coverage
echo "bedtools 2.26.0" >&2
bedtools genomecov -i ${name}.bed -split -bg -g $chrom_info > ${name}.bg
wigToBigWig ${name}.bg $chrom_info ${name}.bw
# Tn5_coverage
echo "bedtools 2.26.0" >&2
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5; else print $1, $3 - 6, $3 - 5}' ${name}.bed > ${name}.Tn5.bed
rm ${name}.sam
rm *.bg *trim*.fastq* 
