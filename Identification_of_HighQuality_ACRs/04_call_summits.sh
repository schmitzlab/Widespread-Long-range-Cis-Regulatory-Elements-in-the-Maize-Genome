#!/bin/bash -e

if [ $# -lt  7]; then
    echo "usage:sh call_summits.sh <output_folder_name> <peaks> <density-file> <tags> <chrom_info> <contig_prefix>"
    echo "
    <output_folder_name>: this scripts will create an folder under current path;
    <peaks>: the output filtered peaks from choose_cutoff.sh, it should include three column, <chr> <start> <end> and it's 0 based;
    <density-file>: the ouput *.density file from the calculate_density step;
    <tag>: the output tag folder from prepareTag;
    <chrom_info>: a file include the chromosome name;
    <contig_prefix>: the name of the contig prefix. If there're no contig, set it as <no>;
    all the path should be the full path."
    exit -1; 
fi

#### specie ####
out=$1
#### peaks to define###
input=$2
#### PATH to density files ####
density=$3
#### path to tag###
tag=$4
####chromosome prefix #####
genome=$5
#### contig prefix ######
ctg=$6

module load parallel

mkdir $out
cd $out

### define summit calling function for chromosomes ####
function cal {
chr=$1
tag=$2
name=$3
echo -n >$name.$chr.sum.bed
mkdir /dev/shm/$name.$chr.aabb  ####store the density files in memory ###
cp $tag/$chr.bg /dev/shm/$name.$chr.aabb/$chr.bg
cat $name.$chr.tmp |while read LINE
do
 s=$(echo "$LINE"|awk '{print $1+1}')  ### bed file is 0 based ##
 t=$(echo "$LINE"|cut -f3)
 n=$(echo "$LINE"|cut -f4)
 value=$(awk '{if($2<'$s') next}{if($2>='$s'&&$2<'$t') print $2,$3; else if($2>='$t') exit}' /dev/shm/$name.$chr.aabb/$chr.bg |sort -n -r -k2 |head -1)
 printf "$chr\t$s\t$t\t$n\t$value\n" >> $name.$chr.sum.bed
done
rm /dev/shm/$name.$chr.aabb/$chr.bg
}
export -f cal

### define summit calling function for contigs ####

function ctgc {
tag=$1
name=$2
echo -n >$name.ctg.sum.bed
mkdir /dev/shm/$name.ctg.aabb  ####store the density files in memory ###
cp $tag/ctg.bg /dev/shm/$name.ctg.aabb/$name.ctg.bg
cat $name.ctg.tmp |while read LINE
do
 chr=$(echo "$LINE"|cut -f1)
 s=$(echo "$LINE"|awk '{print $1+1}') ### bed file is 0 based ##
 t=$(echo "$LINE"|cut -f3)
 n=$(echo "$LINE"|cut -f4)
 value=$(awk '{if($1!="'$chr'") next}{if($1=="'$chr'"&&$2<'$s') next}{if($1=="'$chr'"&&$2>='$s'&&$2<'$t') print $2,$3; else exit}' /dev/shm/$name.ctg.aabb/$name.ctg.bg |sort -n -r -k2 |head -1)
 printf "$chr\t$s\t$t\t$n\t$value\n" >> $name.ctg.sum.bed
done
rm /dev/shm/$name.ctg.aabb/$name.ctg.bg
}

export -f ctgc

##### find the most enriched bins ###

awk '{print $1,$2,$3,$NF}' OFS="\t" $density |bedtools intersect -a - -b $input -wa -wb |awk '{print $1,$2,$3,$5"-"$6"-"$7,$4}' OFS="\t"  |\
sort -n -r -k5 |awk '!seen[$4]++' |bedtools sort -i - |sed "s/-/\t/g" |\
awk '{if($2<$5) print $1,$5,$3,$4"-"$5"-"$6,$7; else if($3>$6) print $1,$2,$6,$4"-"$5"-"$6,$7; else print $1,$2,$3,$4"-"$5"-"$6,$7}' OFS="\t" > $out.highestBIN 


#### output parrellel command ###

echo -n > $out.sum_Commands.txt

chrom=($(sed "/$ctg/d" $genome|cut -f1))
num=${#chrom[@]}

for ((i=0;i<$num;i++))
do
chr=${chrom[$i]}
awk '{if($1=="'$chr'") print $1,$2,$3,$4}' OFS="\t"  $out.highestBIN  > $out.$chr.tmp
echo "cal $chr $tag $out" >> $out.sum_Commands.txt
done

grep -e "$ctg"  $out.highestBIN |cut -f1-4  > $out.ctg.tmp
echo "ctgc $tag $out " >> $out.sum_Commands.txt

parallel < $out.sum_Commands.txt

cat $out.*.sum.bed |cut -f4- |sed "s/-/\t/g" |sed "s/\ /\t/g" |bedtools sort -i - > $out.sum.bed

rm $out.sum_Commands.txt $out.*.tmp $out.*.sum.bed

