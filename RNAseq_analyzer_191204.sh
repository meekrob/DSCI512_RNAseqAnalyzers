#!/usr/bin/env bash

################################################
# PROGRAM:
# RNAseq_analyzer_191204.sh
#
# DESCRIPTION:
# This is a very basic RNA-seq pipeline that I use for analyzing fastq reads. Step1 is a
# simple wrapper that performs quality control, genome alignment, basic format
# conversions, and htseq-count tabulation for paired-end RNA-seq samples using a specified
# genome. Step2 is a clean up program that removes unnecessary files and compresses files
# to save space.
#
# AUTHOR:
# Erin Osborne Nishimura
#
# START DATE:
# December 4, 2019
#
# DEPENDENCIES:
# 	Requires the installation of the follwing software: 
#		fastp
#		hisat2
#		featureCounts
#		samtools
#		deep-tools
#
# 	Requires access to the Nishimura Lab installed software sources on SUMMIT.
#
# REQUIRES:
#    INPUT: .fastq files.    For each sample, paired forward and reverse sequencing files
#								are required. These should be placed in an input
#								directory.
#
#    INPUT: _metadata.txt file: A metadata file with two columns. The first two columns
#								are fastq file names. The third column is a "nickname"
#								of each sample. Later columns can be included with other
#								metadata information. Metadata file should be placed
#								within the inputdir directory.
#
#
#    HISAT2 INDEXES: .ht2 files for the genome. These are produced using hisat2-build. For
#								instructions see
#	           https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer
#
#    GENOME SEQUENCE: .fa  or .tar.gz file for the genome. This is the sequence of the 
#                                genome.
#
#    GENOME ANNOTATION: .gtf file for the genome. This is a genome annotation file of gene
#								features. Version and coordinates must match the genome
#								sequence (.fa above).
#
# USAGE:
# $ bash RNAseq_analyzer_191204.sh <metadata.txt> <number of threads>
#
# OUTPUT:
#
# KNOWN BUGS:
#
# THINGS TO IMPROVE:
#
################################################


echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"


####### MODIFY THIS SECTION #############

#The input samples (metadata file and _fastq.gz files) live in directory:
inputdir="../../01_input/"

#This is where the ht2 files live:
hisat2path="/scratch/summit/erinnish@colostate.edu/DSCI512_RNAseq/PROJ01_ce11Build/ce11"

#This is where the genome sequence lives:
genomefa="/scratch/summit/erinnish@colostate.edu/DSCI512_RNAseq/PROJ01_ce11Build/chromFa.tar.gz"

#This is where the gtf file lives:
gtffile="../../01_input/ce11_annotation_ensembl_to_ucsc.gtf.gz"
    
#Number of threads to use:
pthread=$2


########## DONE MODIFYING ###############


########## SOURCE SOFTWARE ###############
module load singularity
module list

fastp='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif fastp'
hisat2='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif hisat2'
featureCounts='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif featureCounts'
samtools='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif samtools'
bamCoverage='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif bamCoverage'

########## BEGIN CODE ###############

#This is the output_directory:

DATE=`date +%Y-%m-%d`
#DATE='2018-10-16'
outputdir="../../03_output/"$DATE"_output/"


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




# FASTP to remove unwanted sequences
# FASTP to determine quality
echo -e "\n>>> FASTP: Trimming excess and low-quality sequences from .fastq file; generating quality report"
mkdir -p $outputdir"02_fastp"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}
	
	#Chop off the .gz of each name:
	unzippedfile1=${sample1}
	unzippedfile2=${sample2}

	#Make output directories
   	mkdir -p $outputdir"02_fastp/"$samplename
   	
   	## execute fastp
   	cmd1="$fastp -i $inputdir$unzippedfile1 -I $inputdir$unzippedfile2 -o $outputdir"02_fastp/"$samplename/$samplename"_trim_1.fastq" -O $outputdir"02_fastp/"$samplename/$samplename"_trim_2.fastq" -h $outputdir"02_fastp/"$samplename/$samplename"_report.html" -j $outputdir"02_fastp/"$samplename/$samplename"_report.json" --thread $pthread -x"   
    
   	echo -e "\t$ ${cmd1}"
   	time eval $cmd1
	
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
	unzippedfile1=${sample1}
	unzippedfile2=${sample2}

   	## execute hisat2
   	cmd3="$hisat2 -x $hisat2path -1 $outputdir"02_fastp/"$samplename/$samplename"_trim_1.fastq" -2 $outputdir"02_fastp/"$samplename/$samplename"_trim_2.fastq" -S ${outhisat2}${samplename}.sam --summary-file ${outhisat2}${samplename}_summary.txt --no-unal -p $pthread"
   	echo -e "\t$ $cmd3"
	time eval $cmd3

done



# FEATURECOUNTS to tabulate reads per gene:
echo -e "\n>>> FEATURECOUNTS: Run featureCounts on all files to tabulate read counts per gene"
outfeature=$outputdir"04_feature/"
mkdir -p $outfeature

# Acquire a list of .sam names
samfilePath=()
for (( counter=0; counter < ${#names[@]}; counter++ ))
do
    samfile=${names[$counter]}.sam
    samfilePath+=(${outhisat2}${samfile})

done

# Execute featureCounts
cmd4="$featureCounts -p -T $pthread -a $gtffile -o ${outfeature}counts.txt ${samfilePath[*]}"
echo -e "\t$ $cmd4"
time eval $cmd4


# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"05_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
    # echo
    echo -e "\tSamtools and BamCoverage convert: ${seqname}"
    
    # Samtools: compress .sam -> .bam
    cmd5="$samtools view --threads $pthread -bS ${outhisat2}${seqname}.sam > ${samout}${seqname}.bam"
	echo -e "\t$ ${cmd5}"
	time eval $cmd5
	
    
    # Samtools: sort .bam -> _sort.bam
    cmd6="$samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
    echo -e "\t$ ${cmd6}"
    time eval $cmd6
    
    
    # Samtools: index _sort.bam -> _sort.bam.bai
    cmd7="$samtools index ${samout}${seqname}_sort.bam"
    echo -e "\t$ ${cmd7}"
    time eval $cmd7
    
    
    # bamCoverage: 
    cmd8="$bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1"
    echo -e "\t$ ${cmd8}"
    time eval $cmd8
    
done




######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASTP VERSION:"
$fastp --version
echo -e "\n>>> HISAT2 VERSION:"
$hisat2 --version
echo -e "\n>>> SAMTOOLS VERSION:"
$samtools --version
echo -e "\n>>> FEATURECOUNTS VERSION:"
$featureCounts -v
echo -e "\n>>> BAMCOVERAGE VERSION:"
$bamCoverage --version
echo -e ">>> END: Analayzer complete."
