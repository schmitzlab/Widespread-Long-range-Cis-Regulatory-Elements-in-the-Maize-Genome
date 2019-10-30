### required softwares
### HiC-PRO 2.8.0
### java jdk1.8.0_131 
### cuda version 7.0.28
### GPU core was required


allValidPairs=         #### HiC pro output clean pairs 
chrom_info=            #### chromosome information
MboI_Frag=             #### Restriction enzyme digetion fragment files 
output=                #### output folder name

##### hicpro2juicebox.sh comes from HiC pro packages

sh hicpro2juicebox.sh -i $allValidPairs -g $chrom_info -j juicer_tools_0.7.0.jar -r $MboI_Frag -t . -o .

java -jar juicer_tools_0.7.0.jar hiccups -m 512 -r 1000,2000,5000,10000 -f 0.5,0.5,0.5,0.5 -p 10,8,4,2 -i 10,8,7,5 -d 5000,10000,20000,20000 --ignore_sparsity $allValidPairs.hic $output

