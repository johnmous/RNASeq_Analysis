#!/bin/bash
# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title: Convert SAM to BAM, sort and index

list=(1 2 3 7 8 12 13)

for i in ${list[*]}
do
  prefix=$(printf "S%02d" $i)
  
 echo -e "\n ####################### &&&&&&&&&&&&&&&&&& #############################"
 echo -e '\nWorking on sample '$prefix
 
 samtools view -hSb -o bamFiles/$prefix"_onBsubtilis.bam" $prefix"_onBsubtilis.sam" 
 samtools sort bamFiles/$prefix"_onBsubtilis.bam" bamFiles/$prefix"_onBsubtilis.sorted"
 samtools index bamFiles/$prefix"_onBsubtilis.sorted.bam"
 rm bamFiles/$prefix"_onBsubtilis.bam"
 
done
