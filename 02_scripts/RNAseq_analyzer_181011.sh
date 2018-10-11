#!/usr/bin/env bash

################################################
# Program:
# RNAseq_analyzer_181011.sh
#
# Description:
# This is a very basic RNA-seq pipeline that I use for analyzing fastq reads. Step1 is a simple 
# wrapper that performs quality control, genome alignment, basic format conversions, and htseq-count tabulation 
# for paired-end RNA-seq samples using a specified genome. Step2 is a clean up program that removes unnecessary 
# files and compresses files to save space.
#
# Author:
# Erin Osborne Nishimura
#
# Date initiated:
# October 11, 2018
#
# Dependencies:
# Requires the installation of the follwing software: trimmomatic, fastqc, hisat2, featureCounts, samtools, deep-tools
# Requires access to the Nishimura Lab installed software sources on SUMMIT.
#
# Requires:
#    INPUT: .fastq.gz files. For each sample, paired forward and reverse sequencing files are required. 
#                            These should be placed in an input directory.
#
#    INPUT: _metadata.txt file: A metadata file with two columns. The first two columns are fastq.gz 
#                            file names. The third column is a "nickname" of each sample. Later columns 
#                            can be included with other metadata information. Metadata file should be 
#                            placed within the inputdir directory. Example of a metadata file:
#
#    BUILD: .bt2 files for the genome. These are produced using hisat2-build. For instructions 
#                            see https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer
#
#    GENOME: .fa file for the genome. This is the sequence of the genome.
#
#    GENOME: .gtf file for the genome. This is a genome annotation file of gene features. Version 
#                            and coordinates must match the genome sequence (.fa above)
#
# Executed with:
# $ bash RNAseq_analyzer_181011.sh metadata.txt 
# OR
# $ bash RRNAseq_analyzer_181011.sh metadata.txt 2>&1 | tee 181011_output.txt
#
################################################


echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"


####### MODIFY THIS SECTION #############

#The input samples (metadata file and _fastq.gz files) live in directory:
inputdir="../01_input/"

#This is where the bt2 files live:
hisat2path="/projects/erinnish@colostate.edu/genomes/mm10/from_ucsc/mm10"

#This is where the genome sequence lives:
genomefa="/projects/erinnish\@colostate.edu/genomes/mm10/from_ucsc/chromFa.tar.gz"

#This is where the gtf file lives:
gtffile="/projects/erinnish@colostate.edu/genomes/mm10/from_ensembl/gtf/Mus_musculus_GRCm38_2UCSC.gtf"
    
#Number of threads to use:
pthread=$2


########## DONE MODIFYING ###############





#This is the output_directory:
DATE=`date +%Y-%m-%d`
outputdir="../03_output/"$DATE"_output/"


echo -e ">>> MAKING output directory"
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir



####### META DATA #############



#These are the sample names, R1:
samples1=( $(cut -f 1 --output-delimiter=' ' $1) )

#These are the sample names, R2:
samples2=( $(cut -f 2 --output-delimiter=' ' $1) )

#These are the nicknames I want to give the files:
names=( $(cut -f 3 --output-delimiter=' ' $1) )



####### PIPELINE ##############

# Report back to the user which files will be processed and which names they'll be given:
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$1"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    echo -e "\t${samples1[$counter]}\t${samples2[$counter]}\t${names[$counter]}"
done



# UNZIP all the input files: 
echo -e "\n>>> GUNZIP: unzipping _R1.fastq.gz files:"
for fastqfile in ${samples1[@]}
do
   echo -e "\tgunzipping $fastqfile"
   gunzip -q $inputdir$fastqfile
done

echo -e ">>> GUNZIP: unzipping _R2.fastq.gz files:"
for fastqfile in ${samples2[@]}
do
   echo -e "\tgunzipping $fastqfile"
   gunzip -q $inputdir$fastqfile
done




# TRIMMOMATIC to remove unwanted sequences
# FASTQC to determine quality
echo -e "\n>>> TRIMMOMATIC: Trimming excess sequences from .fastq file"
mkdir -p $outputdir"01_trimmomatic"
#for fastqfile in ${samples1[@]}
#do
  
    #ENTER TRIMMOMATIC STUFF HERE
    
#done


# FASTQC to determine quality
echo -e "\n>>> FASTQC: analyzing quality of each .fastq file"
mkdir -p $outputdir"02_fastqc"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}
	
	#Chop off the .gz of each name:
	unzippedfile1=${sample1//.gz/}
	unzippedfile2=${sample2//.gz/}

	#Make output directories
   	mkdir -p $outputdir"02_fastqc/"$samplename"_1"
   	mkdir -p $outputdir"02_fastqc/"$samplename"_2"	
   	
   	# execute fastqc
   	cmd1="fastqc -o $outputdir"02_fastqc/"$samplename"_1" -t 20 $inputdir$unzippedfile1"
   	echo -e "\t$ ${cmd1}"
   	time eval $cmd1
   	
   	cmd2="fastqc -o $outputdir"02_fastqc/"$samplename"_2" -t 20 $inputdir$unzippedfile2"
   	echo -e "\t$ ${cmd2}"
   	time eval $cmd2
	
done



# HISAT2 to align to the genome
echo -e "\n>>> HISAT2: aligning each sample to the genome"
outhisat2=$outputdir"03_hisat2/"
mkdir -p $outhisat2

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}
	
	#Chop off the .gz of each name:
	unzippedfile1=${sample1//.gz/}
	unzippedfile2=${sample2//.gz/}

   	# execute hisat2
   	cmd3="hisat2 -x $hisat2path -1 ${inputdir}${unzippedfile1} -2 ${inputdir}${unzippedfile2} -S ${outhisat2}${samplename}.sam --summary-file ${outhisat2}${samplename}_summary.txt --no-unal -p $pthread"
   	echo -e "\t$ $cmd3"
	time eval $cmd3

done



# FEATURECOUNTS to tabulate reads per gene:
echo -e "\n>>> HTSEQ: Run HTSeq-counts on all files to tabulate read counts per gene"
outhtseq=$outputdir"04_feature/"
mkdir -p $outhtseq

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samfile=${names[$counter]}.sam
    
    cmd4="htseq-count --stranded=yes --minaqual=20 -r pos -q ${outhisat2}${samfile} $gtffile > ${outhtseq}${names[$counter]}_counts.txt"
    echo -e "\t$ $cmd4"
    time eval $cmd4

done



# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"05_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
    # echo
    echo -e "\tSamtools and BamCoverage convert: ${seqname}"
    
    # Samtools: compress .sam -> .bam
    cmd5="samtools view --threads $pthread -bS ${outhisat2}${seqname}.sam > ${samout}${seqname}.bam"
	echo -e "\t$ ${cmd5}"
	time eval $cmd5
	
    
    # Samtools: sort .bam -> _sort.bam
    cmd6="samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
    echo -e "\t$ ${cmd6}"
    time eval $cmd6
    
    
    # Samtools: index _sort.bam -> _sort.bam.bai
    cmd7="samtools index ${samout}${seqname}_sort.bam"
    echo -e "\t$ ${cmd7}"
    time eval $cmd7
    
    # bamCoverage: 
    cmd8="bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1"
    echo -e "\t$ ${cmd8}"
    time eval $cmd8
    
done




######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASEQC VERSION:"
fastqc --version
echo -e "\n>>> HISAT2 VERSION:"
hisat2 --version
echo -e "\n>>> SAMTOOLS VERSION:"
samtools --version
echo -e "\n>>> HTSEQ-COUNT VERSION:"
htseq-count --help | grep "version"
echo -e "\n>>> BAMCOVERAGE VERSION:"
export PYTHONPATH=/projects/dcking@colostate.edu/lib/python3.5/site-packages/
bamCoverage --version
echo -e ">>> END: Analayzer complete."
