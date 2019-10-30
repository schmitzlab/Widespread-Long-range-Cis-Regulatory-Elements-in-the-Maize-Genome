#!/bin/bash -e

if [ $# -lt 7 ]; then
    echo "usage:sh cutoff.sh <density_file> <cutoff> <merge-gaps> <peaks-length> <reps_counts> <genome_fasta> <NCBI_MtPt_Db> <name> "
    echo "
    <density_file>: the output file from calculate_density step; 
    <cutoff>: the density cutoff to choose for each bin;
    <merge-gaps>: the maxium length allowed to merge nearby bins;
    <peaks-length>: the minium length of a peak;
    <reps_counts>: sample replication number;
    <genome_fasta>: PATH to the genome fasta file;
    <NCBI_MtPt_Db>: PATH to the NCBI blast+ database of plant mitochrodria and plastids sequences;
    <name>: the name of the output peak files
    the NCBI blast+ and BEDTools should be in the system path."
    exit -1; 
fi

density=$1
cutoff=$2
gap=$3
len=$4
reps=$5
genome=$6
path_MtPt_Db=$7

awk '$NF>('$reps'*'$cutoff')' $density |bedtools sort -i - |bedtools merge -d $gap -i - |awk '$3-$2>'$len'' >  $density.$cutoff.$gap.bedo

export PATH=$PATH:/home/zefulu/soft/ncbi-blast-2.8.1+/bin

bedtools getfasta -fi $genome -bed $input.$density.$gap.bedo -fo $input.$density.$gap.fa

blastn -query $input.$density.$gap.fa -out $name.mtpt -db path_MtPt_Db -outfmt 7

sed '/#/d' $name.mtpt | cut -f1 |sed "s/\:/\t/g" |sed "s/-/\t/g" |bedtools sort -i - | uniq > $name.black

bedtools intersect -a  $density.$cutoff.$gap.bedo -b $name.black -v |cut -f1-3 > $name

rm $input.$density.$gap.fa
rm $name.mtpt $density.$cutoff.$gap.bedo
