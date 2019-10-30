#!/bin/bash -e

if [ $# -lt 6 ]; then
    echo "usage:sh *.sh <output_name> <bin_size> <step> <genome_size> <path_to_sample_file> <path_to_input_file>"
    echo "
    To run this scripts, MACS2 and BEDTools should be in the system path;
    <output_name>: this scripts will create a folder under current path;
    <bin_size>: the length to split each peak;
    <step>: the overlapped length between to windowns(bins);
    <genome_size>: the total length of the genome;
    <path_to_sample files>: a txt file including the path to the sample files, BED format;
    <path_to_input files>: a txt file including the path to the input files, BED format;
    all the path should be the full path"
    exit -1; 
    exit -1; 
fi

name=$1
bin=$2
step=$3
size=$4
sample_list=$5
input_list=$6

#### set up sample and input samples ###
sample=($(cat $sample_list))
input=($(cat $input_list))
a=${#sample[@]}
b=${#input[@]}
x=$((${a}-1))

mkdir $name
cd $name
####call peak with macs2 ###

macs2 callpeak -t ${sample[@]} -c ${input[@]} -g $size --keep-dup all --nomodel --extsize 147 -n $name
peaks=${name}_peaks.narrowPeak
#### change BED reads to Tn5 intergration sites ####

for((i=0;i<$a;i++))
do
reads=${sample[$i]}
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5; else print $1, $3 - 6, $3 - 5}' $reads > sample_rep$i.Tn5
done

for((i=0;i<$b;i++))
do
reads=${input[$i]}
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5; else print $1, $3 - 6, $3 - 5}' $reads > input_rep$i.Tn5
done

#### make windows from macs2 peak file and split them into windows  (max line number 100000) #####
bedtools makewindows -b $peaks -w $bin -s $step|bedtools sort -i - |awk '{print $1,$2,$3,"bin_"NR}' OFS="\t" > $name.windows0000.bed
l=$(cat $name.windows0000.bed|wc -l)
m=$(($l/100000))
split -l 100000 -d -a ${#m} $name.windows0000.bed ${name}_split_

#### calculate the Tn5 integration frequency of each window #####
for n in `seq -w 0 $m`
do
 i=0
 reads=sample_rep$i.Tn5
 total$i=$(cat $reads |wc -l)
 bedtools coverage -a ${name}_split_$n -b $reads -counts |awk '{print $1,$2,$3,$4,'$size'*$5/'$[total$i]'/($3-$2)}' OFS="\t" > ${name}_split_$n.$i.bgx
 for((i=1;i<$a;i++))
 do
 	j=$(($i-1))
  reads=sample_rep$i.Tn5
  total$i=$(cat $reads |wc -l)
  bedtools coverage -a ${name}_split_$n -b $reads -counts |awk '{print '$size'*$5/'$[total$i]'/($3-$2)}' OFS="\t" > ${name}_split_$n.$i.bg
  paste ${name}_split_$n.$j.bgx ${name}_split_$n.$i.bg > ${name}_split_$n.$i.bgx
  rm ${name}_split_$n.$j.bgx ${name}_split_$n.$i.bg
 done
 rm ${name}_split_$n
done
cat ${name}_split_*.$x.bgx |awk '{c=0;for(i=5;i<=NF;++i){c+=$i};print $0,c}' OFS="\t" > $name.mergedBg

awk '{print $1,$2,$3,$4,$NF}' OFS="\t" $name.mergedBg > $name.density

