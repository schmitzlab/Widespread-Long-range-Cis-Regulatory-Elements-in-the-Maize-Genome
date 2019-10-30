### required softwares
# bedtools/2.26.0
# MACS2 
# UCSC tools version 359


name=          ### output file prefix
HQ_ACR=        ### Coordinates of High Quality ACRs
sample_list=   ### Plain txt containing the sample path
control_list=  ### Plain txt containing the control path

bg=50000       ### background regions 
fold=2         ### enriched fold cutoff
len=50         ### windows to split to
size=2.1e9     ### genome sizes

################## code start here #################################

sample=($(cat $sample_list))
input=($(cat $input_list))

####call raw diff peak with macs2 ###

macs2 callpeak -t ${sample[@]} -c ${input[@]} -g $size --keep-dup all --nomodel --extsize 147 -n $name
diff=${name}_peaks.narrowPeak

#### filter raw diff Peaks with HQ_ACRs
bedtools intersect -a $HQ_ACR -b $diff |cut -f1-3 |awk '$3-$2>'$len'' > $name.tmp000
#generate background BED 
awk '{if($2>='$bg') print $1,$2-'$bg',$3+'$bg'; else print $1,"0",$3+'$bg'}' OFS="\t" $name.tmp000 >$name.background.bed

a=${#sample[@]}
x=$((${input}-1))


i=0
reads=${sample[$i]}
bedtools coverage -a $name.tmp000 -b $reads -counts > $name.pks.$i.bg
bedtools coverage -a $name.background.bed -b $reads -counts | cut -f4> $name.background.$i.bg
paste $name.pks.$i.bg $name.background.$i.bg |awk '{print $1,$2,$3,(2*'$bg'+$3-$2)/($3-$2)*$4/$5}' OFS="\t" > $name.sample.$i.dsx
rm $name.pks.$i.bg $name.background.$i.bg
for ((i=1;i<$a;i++))
do
 j=$(($i-1))
 reads=${sample[$i]}
 bedtools coverage -a $name.tmp000 -b $reads -counts > $name.pks.$i.bg
 bedtools coverage -a $name.background.bed -b $reads -counts | cut -f4 > $name.background.$i.bg
 paste $name.pks.$i.bg $name.background.$i.bg |awk '{print (2*'$bg'+$3-$2)/($3-$2)*$4/$5}' > $name.sample.$i.ds
 paste $name.sample.$j.dsx $name.sample.$i.ds > $name.sample.$i.dsx
 rm $name.pks.$i.bg $name.background.$i.bg $name.sample.$j.dsx $name.sample.$i.ds
done

mv $name.sample.$x.dsx $name.sample.density
rm $name.sample.*.dsx

control=($(cat $control_list))
b=${#control[@]}
y=$((${input}-1))


i=0
reads=${control[$i]}
bedtools coverage -a $name.tmp000 -b $reads -counts > $name.pks.$i.bg
bedtools coverage -a $name.background.bed -b $reads -counts | cut -f4> $name.background.$i.bg
paste $name.pks.$i.bg $name.background.$i.bg |awk '{print $1,$2,$3,(2*'$bg'+$3-$2)/($3-$2)*$4/$5}' OFS="\t" > $name.control.$i.dsx
rm $name.pks.$i.bg $name.background.$i.bg
for ((i=1;i<$b;i++))
do
 j=$(($i-1))
 reads=${control[$i]}
 bedtools coverage -a $name.tmp000 -b $reads -counts > $name.pks.$i.bg
 bedtools coverage -a $name.background.bed -b $reads -counts | cut -f4 > $name.background.$i.bg
 paste $name.pks.$i.bg $name.background.$i.bg |awk '{print (2*'$bg'+$3-$2)/($3-$2)*$4/$5}' > $name.control.$i.ds
 paste $name.control.$j.dsx $name.control.$i.ds > $name.control.$i.dsx
 rm $name.pks.$i.bg $name.background.$i.bg $name.control.$j.dsx $name.control.$i.ds
done

mv $name.control.$y.dsx $name.control.density
rm $name.control.*.dsx


awk '{for(i=4;i<=NF;i++) t+=$i; print $0,t; t=0}' $name.sample.density | awk '{print $1,$2,$3,$NF/(NF-4)}' OFS="\t"  > $name.sample.densityx
awk '{for(i=4;i<=NF;i++) t+=$i; print $0,t; t=0}' $name.control.density | awk '{print $1,$2,$3,$NF/(NF-4)}' OFS="\t"  > $name.control.densityx

paste $name.sample.densityx $name.control.densityx | cut -f1-4,8 > $name.densityx
rm $name.sample.densityx $name.control.densityx
awk '{if($5==0) print $0,"NA"; else if($4/$5>='$fold') print $0,$4/$5}' OFS="\t" $name.densityx > $name.diff.peaks.bed
