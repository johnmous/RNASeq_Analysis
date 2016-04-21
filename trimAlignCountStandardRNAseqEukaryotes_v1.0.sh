#!/bin/bash
# Author: Ioannis Moustakas, i.moustakas@gmail.com
# Title: Standard script for RNAseq for Eucaryotes (makes use of Tophat)
# Description: 	1. Trim the fastq files using trimmomatric. Keep reads with length>25
#		2. Align on spikes/ERCCs 
#		3. Collect the reads that did not map on spikes/ERCCs
#		4. Map on genome with Tophat. Use the settings provided in the protocol
#		5. Go through alignment file and only keep reads with a number of mismatches less than 10% of their length.
#		   Create a new alignment file (.sam) to save reads that have #mismatches<0.10*read_length
#		6. Count the reads mapped on each gene with htseq-count. This results in a count table 
#		7. Count the reads that were discaded in step 5 and the reads left unmapped by tophat.
#		   Save these numbers at the bottom of the count table



echo -e "Standard script for RNAseq for Eucaryotes (makes use of Tophat) \n " > pipeline.out

rawData=/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Data-QC/basicQC/raw
outputDir=/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophat
db=/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/sly/DNA/genome/S_lycopersicum_chromosomes.2.50
countTableMAPQZeroDir=$outputDir/MAPQZero
countTableMAPQTenDir=$outputDir/MAPQTen

tomGenome="/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/sly/DNA/genome/S_lycopersicum_chromosomes.2.50"
tomatoGFF="/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/sly/DNA/genome/ITAG2.4_gene_models.gff3"

# if no path for the organism's transcriptome exist, create the transcriptome
# This can be omited, and run for every sample by tophat, but it saves time if it is done before hand 
transcrIndexPath=$outputDir/transcriptomeTomato 
if [ ! -d $transcrIndexPath ]; then
  tophat2 -G $tomatoGFF --transcriptome-index=$transcrIndexPath/S_lycopersicum $tomGenome
fi

if [ ! -d $outputDir/ERCCsAlignment ]; then
  mkdir $outputDir/ERCCsAlignment
fi

if [ ! -d $countTableMAPQZeroDir ]; then
  mkdir $countTableMAPQZeroDir
fi

if [ ! -d $countTableMAPQTenDir ]; then
  mkdir $countTableMAPQTenDir
fi
  
if [ ! -d $outputDir/trimmedFastq ]; then
  mkdir $outputDir/trimmedFastq
fi  

if [ ! -d $outputDir/TophatAlign/ ]; then
  mkdir $outputDir/TophatAlign/
fi 

# ls all fastq files in the basicQC/raw path
fastqFiles=$(ls $rawData)

for fastqFile in ${fastqFiles[@]};
do
  fastqFilePrefix=${fastqFile%.*}
  echo $fastqFilePrefix
  # Filtering reads
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Trimming "$fastqFilePrefix >> pipeline.out
  java -jar /mad/software/src/Trimmomatic-0.30/trimmomatic-0.30.jar SE -threads 20 -phred33 -trimlog log.out \
  $rawData/$fastqFilePrefix.fastq $outputDir/trimmedFastq/$fastqFilePrefix"_Trimmed.fastq"  MINLEN:25 &>> pipeline.out
  
  # align on ERCCs
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Aligning "$fastqFilePrefix" on ERCCs" >> pipeline.out
  bowtie2 -q -p 15 --norc -x /zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/ERCC/ERCC92 \
  -U $outputDir/trimmedFastq/$fastqFilePrefix"_Trimmed.fastq"  \
  -S $outputDir/ERCCsAlignment/$fastqFilePrefix"_onERCC92.sam" &>> pipeline.out
  
  # Collect the reads that did NOT map on ERCCs
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Collecting UnMapped reads from "$fastqFilePrefix" and puting them in a Fastq file" >> pipeline.out
  cat $outputDir/ERCCsAlignment/$fastqFilePrefix"_onERCC92.sam" \
  | awk '$3=="*" {printf("@%s\n%s\n+%s\n%s\n",$1,$10,$1,$11)}' > $outputDir/ERCCsAlignment/$fastqFilePrefix"_notMappedOnERCC92.fastq"

  # Run Tophat. Each sample will be saved in a different directory
  outputPath=$outputDir/TophatAlign/$fastqFilePrefix"Tophat"
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Aligning "$fastqFilePrefix" on genome using tophat" >> pipeline.out
  tophat2 -N 30 --read-gap-length 30 --read-edit-dist 30 -o $outputPath -p 20 --no-coverage-search --library-type fr-secondstrand \
  --no-novel-juncs --transcriptome-index=$transcrIndexPath/S_lycopersicum --transcriptome-only --no-novel-indels \
  $tomGenome $outputDir/ERCCsAlignment/$fastqFilePrefix"_notMappedOnERCC92.fastq" &>> pipeline.out
  
  
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Filtering on the alignments with number of mismatches > 10% for: "$fastqFilePrefix >> pipeline.out
  # Save bam file header
  samtools view -H $outputPath/accepted_hits.bam > $outputPath/acceptedHitsFiltered.sam
  
  # only save reads with #MisMaches/ReadLength < 10%
  samtools view $outputPath/accepted_hits.bam |\
  awk '{l=length($10); nm=substr($17, 6); if (nm/l<0.1) {print $0}}' \
  > $outputPath/acceptedHitsFiltered.sam
  
  # Calculate the number of reads that were discarded in the previous step
  discardedReadsMM=$(samtools view $outputPath/accepted_hits.bam |\
  awk '{l=length($10); nm=substr($17, 6); if (nm/l>=0.1) {countDiscarded+=1}} END \
  {print countDiscarded}')
 
  # Count the Number of reads left unmapped by Tophat
  tophatUnmappedReads=$(samtools view $outputPath/unmapped.bam | wc -l)
  
  # First count the reads with MAPQ zero. Then one can estimate the rRNA content
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Counting with minimum read MAPQ=0: "$fastqFilePrefix >> pipeline.out
  htseq-count -a 0 $outputPath/acceptedHitsFiltered.sam $tomatoGFF > $countTableMAPQZeroDir/countMinQualZero_$fastqFilePrefix.txt
  
  # Save the statistics that were previously calculated in the count table
  echo -e "__Reads Discarded (MisMatches>10%)\t" $discardedReadsMM >> $countTableMAPQZeroDir/countMinQualZero_$fastqFilePrefix.txt
  echo -e "__Reads UnMapped Tophat\t" $tophatUnmappedReads >> $countTableMAPQZeroDir/countMinQualZero_$fastqFilePrefix.txt
  
  # Now count the reads with MAPQ Ten. This will be reported to the client
  echo -e "\n ####################### &&&&&&&&&&&&&&&&&& ############################# \n" >> pipeline.out
  echo ">>>> Counting with minimum read MAPQ=10: "$fastqFilePrefix >> pipeline.out
  htseq-count -a 10 $outputPath/acceptedHitsFiltered.sam $tomatoGFF > $countTableMAPQTenDir/countMinQualTen_$fastqFilePrefix.txt
  
  # Save the statistics that were previously calculated in the count table
  echo -e "__Reads Discarded (MisMatches>10%)\t" $discardedReadsMM >> $countTableMAPQTenDir/countMinQualTen_$fastqFilePrefix.txt
  echo -e "__Reads UnMapped Tophat\t" $tophatUnmappedReads >> $countTableMAPQTenDir/countMinQualTen_$fastqFilePrefix.txt
done
