#!/bin/bash
# Align reads on dastabases as indicated by Adrienne workflow
# Convert sam to bam to save space and allow sorting => load to Gene browser
# Output idx statistics for each alignment in a seperate file
# Author: Ioannis Moustakas, i.moustakas@uva.nl

readPath=/zfs/datastore0/group_root/MAD-RBAB/03_RBAB-internal/MAD1000-P013-cfRNA_mouse_tumors/MAD1000-P013-E002_2015_setup_inventory_rdekker2/Data-QC/raw/
databasePath=/zfs/datastore0/group_root/MAD-RBAB/03_RBAB-internal/MAD1000-P013-cfRNA_mouse_tumors/MAD1000-P013-E002_2015_setup_inventory_rdekker2/Scratch/databases

# list with the filenames prefixes 
readFileList=(miRNAseq_CS miRNAseq_FS RNAseq_CS RNAseq_FS)

# major output file with alignment statitstics 
bowtieReport=BowtieReport_rRNAMerged.txt

# remove the file where the output will be saved, if it already exists. 
# Otherwise, the script will start appending lines at the end of it
if [ -f $bowtieReport ]; then
    rm $bowtieReport
fi

# Path to save the intermediate output files
outputPath=PipelineOut_rRNAMerged

# create the output path if it doesnt exist
if [ ! -d "$outputPath" ]; then
  mkdir $outputPath
fi

for readsFile in ${readFileList[@]};
do
  ###### map on spikes
  databaseName=Spikes

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 -L 6 -i S,0,0.5 --norc --score-min L,-1,-0.6 \
  -x /zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/RBAB/spikes-smallRNAseq/spikes \
  -U $readPath$readsFile".fastq" -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport
		
  echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut

  ###### on map miRNA
  databaseName=miRNA

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 -L 6 -i S,0,0.5 --norc --score-min L,-1,-0.6 \
  -x $databasePath/miRNA/miRNA_mmu \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

    echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut


  ###### map on si/piRNA
  databaseName=si-piRNA

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 -L 6 -i S,0,0.5 --score-min L,-1,-0.6 \
  -x $databasePath/si-piRNA/si-piRNA_mmu \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

    echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut


  ###### map on sno/snRNA
  databaseName=sno-snRNA

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 -L 6 -i S,0,0.5 --score-min L,-1,-0.6 \
  -x $databasePath/sn-snoRNA/sn-snoRNA_mmu \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

    echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut


  ###### map on tRNA
  databaseName=tRNA

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 --norc \
  -x $databasePath/tRNA/tRNA_mmu \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

  echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut


  ###### map on rRNAmerged
  databaseName=rRNAmerged

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 --norc \
  -x $databasePath/rDNA/rRNAmerged \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

  echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut
  
  ###### map on lncRNADB
  databaseName=lncRNADB

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 --norc  \
  -x $databasePath/lncRNA/lncrnaDB \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

  echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 --norc \
  -x $databasePath/lncRNA/lncNONCODE \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

  echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut
  
  ###### map on cDNA
  databaseName=cDNA

  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################ \n" >> $bowtieReport
  echo -e "Aligning "$readsFile".fastq on "$databaseName"\n" >> $bowtieReport
  bowtie2 -q -p 20 --norc \
  -x /zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/mmu/DNA/cDNA/Mus_musculus.GRCm38.75.cdna.all \
  -U $fastqOut -S $outputPath/$readsFile"On"$databaseName".sam" &>> $bowtieReport

  echo -e " \n ##### convert sam output to bam, then sort, index and output idxstats"
  samtools view -S -b -o $outputPath/$readsFile"On"$databaseName".bam"  $outputPath/$readsFile"On"$databaseName".sam" 
  samtools sort $outputPath/$readsFile"On"$databaseName".bam" $outputPath/$readsFile"On"$databaseName".sorted"
  samtools index $outputPath/$readsFile"On"$databaseName".sorted.bam"	
  samtools idxstats $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/$readsFile"_On"$databaseName"ReadCount.txt"
  # remove the sam + unsorted bam file to save space
  rm $outputPath/$readsFile"On"$databaseName".sam"
  rm $outputPath/$readsFile"On"$databaseName".bam"
		
  echo -e " \n ##### Create a fastq file with reads that did not map on "$databaseName"\n" 
  # Extract unmapped reads with samtools
  samtools view -b -f 4 $outputPath/$readsFile"On"$databaseName".sorted.bam" > $outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam"
  fastqOut=$outputPath/unmapped_$readsFile"_On_"$databaseName".fastq"

  # picard tools sam => fastq file to feed to bowtie onthe next step
  java -jar /zfs/datastore0/software/src/picard-tools/SamToFastq.jar \
  INPUT=$outputPath/unmapped_$readsFile"On"$databaseName".sorted.bam" \
  FASTQ=$fastqOut
  
done
