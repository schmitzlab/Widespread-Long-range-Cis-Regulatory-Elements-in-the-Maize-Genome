#!/bin/bash -e

if [ $# -lt 5 ]; then
    echo "usage:sh prepareTag.sh <output_tag_name> <chrom_info> <contig_prefix> <genome_size> <path_to_Tn5_sites> "
    echo "
    <output_tag_name>: this scripts will create a folder under current path;
    <chrom_info>: a txt include the chromosome name;
    <contig_prefix>: the prefix name of the contigs. If there're no contig, set it as <no>;
    <genome_size>: the total length of the genome;
    <path_to_Tn5_sites>: a txt file including the path to the Tn5 integration files;
    all the path should be the full path"
    exit -1; 
fi

#### specie ####
name=$1
####chromosome info #####
genome=$2
#### contig prefix ######
ctg=$3
### genome size ###
size=$4
#### sample list ####
sample_list=$5


mkdir $name

chrom=($(sed "/$ctg/d" $genome| cut -f1))
num=${#chrom[@]}

sample=($(cat $sample_list))
a=${#sample[@]}


for ((i=0;i<$num;i++))
do
chr=${chrom[$i]}
echo -n > ${name}/$chr.1.bg
 for((x=0;x<$a;x++))
 do
  reads=${sample[$x]}
  total$x=$(cat $reads |wc -l)
  awk '$1=="'$chr'"' ${reads} |awk '{print $1"-"$2}' OFS="\t" | uniq -c |sed "s/\ /\t/g" |awk '{print $NF,$(NF-1)}' OFS="\t" |awk '{print $1,'$size'*$2*/'$[total$x]'}' OFS="\t" >> ${name}/$chr.1.bg
 done
sort -k1 $name/$chr.1.bg | awk '{seen[$1]+=$2}END{for (i in seen) print i, seen[i]}' OFS="\t" |sed "s/-/\t/g" | sort -n -k2 > $name/$chr.bg
rm ${name}/$chr.1.bg
done


if [ $ctg = "no" ]
then
echo "no contig"
else
 echo -n > $name/ctg.1.bg
 for((x=0;x<$a;x++))
 do
  reads=${sample[$x]}
  total$x=$(cat $reads |wc -l)
  grep -e "$ctg" ${reads} |awk '{print $1"-"$2}' OFS="\t" | uniq -c |sed "s/\ /\t/g" |awk '{print $NF,$(NF-1)}' OFS="\t" |awk '{print $1,'$size'*$2*/'$[total$x]'}' OFS="\t" >>  $name/ctg.1.bg
 done
sort -k1,1 -k2n $name/ctg.1.bg | awk '{seen[$1]+=$2}END{for (i in seen) print i, seen[i]}' OFS="\t" |sed "s/-/\t/g" | sort -k1,1 -k2n > $name/ctg.bg
rm ${name}/ctg.1.bg
fi
