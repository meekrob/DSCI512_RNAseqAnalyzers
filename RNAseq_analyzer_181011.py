#!/usr/bin/env python3
################################################
# Program:
# RNAseq_analyzer_181011.sh
#
# Description: Python version by David C. King of a script by Erin Osborne Nishimura

import sys, os
import subprocess
import time
from datetime import date

def run_timed(*args):
    # run commands in this function to get timed output
    # may need to handle **kwargs
    return subprocess.run(["time"] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

progname=sys.argv[0]
args=" ".join(sys.argv[1:])
print("INITIATING analyzer with command:\n", progname, args, sep="\t")

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
pthread=sys.argv[2]


########## DONE MODIFYING ###############

timestamp=time.time()
DATE=date.fromtimestamp(timestamp)

#This is the output_directory:
outputdir="../03_output/%s_output/" % DATE

print(">>> MAKING output directory")
print("mkdir", outputdir)
os.makedirs(outputdir, exist_ok=True)

####### META DATA #############


meta_data_file=sys.argv[1]
meta_data_fh=open(meta_data_file)
samples1 = [] # sample names, R1
samples2 = [] # sample names, R2
names = [] # nicknames I want to give to files

for line in meta_data_file:
    R1,R2,name = line.strip().split()
    sample1.append(R1)
    sample2.append(R2)
    names.append(name)

meta_data_fh.close() 

####### PIPELINE ##############

# Report back to the user which files will be processed and which names they'll be given:
print(">>> INPUT: This script will process files from the metafile:\n\t%s" % meta_data_file)
print(">>> PLAN: This script will process the sample files into the following names: ")

print("\tSAMPLE1\tSAMPLE2\tNAMES")
for sample1,sample2,samplename in zip(samples1,samples2,names):
    print(sample1, sample2, name, sep="\t")


# UNZIP all the input files: 
print("\n>>> GUNZIP: unzipping _R1.fastq.gz files:")
for fastqfile in samples1:
    print("\tgunzipping", fastqfile)
    print( run_timed( [ "gunzip", "-q", fastqfile] ) )

print("\n>>> GUNZIP: unzipping _R2.fastq.gz files:")
for fastqfile in samples2:
    print("\tgunzipping", fastqfile)
    print( run_timed( [ "gunzip", "-q", fastqfile] ) )

# strip the .gz here?
unzipped1 = [f.replace('.gz', '') for f in samples1]
unzipped2 = [f.replace('.gz', '') for f in samples2]

# TRIMMOMATIC to remove unwanted sequences
# FASTQC to determine quality
print("\n>>> TRIMMOMATIC: Trimming excess sequences from .fastq file")
os.make_dirs("%s01_trimmomatic" % outputdir, exit_ok=True)

# FASTQC to determine quality
print("\n>>> FASTQC: analyzing quality of each .fastq file")
os.make_dirs( os.path.join(outputdir , "02_fastqc" ), exist_ok=True)

for samplename, unzippedfile1, unzippedfile2 in zip(names, unzipped1, unzipped2):

    #Make output directories
    sample_1_outputdir = os.path.join( outputdir , "02_fastqc/" , samplename + "_1" )
    sample_2_outputdir = os.path.join( outputdir , "02_fastqc/" , samplename + "_2" )

    inputfile = os.path.join(inputdir, unzippedfile1)

    # R1
    cmd1_args = [ 'fastqc', '-o', sample_1_outputdir, "-t 20", inputfile ]
    print("\t$", " ".join(cmd1_args))
    print( run_timed(cmd1_args) )
    # R2
    cmd2_args = [ 'fastqc', '-o', sample_2_outputdir, "-t 20", inputfile ]
    print("\t$", " ".join(cmd2_args))
    print( run_timed(cmd2_args) )

# HISAT2 to align to the genome
print("\n>>> HISAT2: aligning each sample to the genome")
outhisat2 = os.path.join( outputdir , "03_hisat2/" )
os.make_dirs(outhisat2, exist_ok=True)

for samplename, unzippedfile1, unzippedfile2 in zip(names, unzipped1, unzipped2):
    # hisat2 command
    cmd3_args = [ "hisat2",
                 '-x', hisat2path, # indexes
                 '-1', os.path.join(inputdir, unzippedfile1), # left reads
                 '-2', os.path.join(inputdir, unzippedfile2), # right reads
                 '-S', os.path.join(outhisat2, samplename + '.sam'), # SAM output
                 '--summary-file', os.path.join(outhisat2, samplename + '_summary.txt'),
                '--no-unal', # no unaligned
                '-p', pthread ]

    print("\t$"," ".join(cmd3_args))
    print( run_timed(cmd3_args) )

# FEATURECOUNTS to tabulate reads per gene:
print("\n>>> HTSEQ: Run HTSeq-counts on all files to tabulate read counts per gene")
outhtseq=$outputdir + "04_feature/"
os.make_dirs(outhtseq, exist_ok=True)

for samplename in names:
    samfile = os.path.join(outhisat2, samplename + '.sam')
    
    cmd4_args = ["htseq-count", "--stranded=yes", "--minaqual=20", "-r", "pos", "-q", samfile, gtffile ]
    print("\t$" + " ".join( cmd4_args ))

    # separate stdout,err in the PIPEs this time
    output = run_timed(cmd4_args)
    outfh = open( os.path.join( outhtseq, samplename + '_counts.txt' )
    outfh.write(output.stdout)
    outfh.close()
    print(output.stderr) # should be output from "time"


# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
print("\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:")
samout=os.path.join(outputdir, "05_samtools/")
os.make_dirs(samout, exit_ok=True)

for seqname in names:

    print("\tSamtools and BamCoverage convert:", seqname)

    samfile = os.path.join(outhisat2, samplename + '.sam')
    bamfile = os.path.join(samout, samplename + '.bam')
    
    # Samtools: compress .sam -> .bam
    cmd5_args = [
            "samtools", "view",
            "--threads", pthread,
            "-bS", samfile] 

	print("\t$", " ".join(cmd5_args))
	output = subprocess.run(["time"] + cmd5_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # write binary output
    outfh = open(bamfile, "rb")
    outfh.write(output.stdout)
    outfh.close()
    print(output.stderr) # output from 'time'
	
    
    # Samtools: sort .bam -> _sort.bam
    sortedbamfile = os.path.join(samout, samplename + '_sort.bam')
    cmd6_args = [
        "samtools", "sort",
        "--threads", pthread,
        "-o", sortedbamfile,
        "--reference", genomefa,
        bamfile]

    print("\t$ ", " ".join(cmd6_args))
    print( run_timed(cmd6_args) )
    
    
    # Samtools: index _sort.bam -> _sort.bam.bai
    cmd7_args= [ "samtools", "index", sortedbamfile]
    print("\t$", " ".join(cmd7_args))
    print(run_timed( cmd7_args) )
    
    # bamCoverage: 
    bw_file = os.path.join(samout, samplename + '.bw')
    cmd8_args= [
         "bamCoverage",
            "-b", sortedbamfile,
            "-o", bwfile,
            "--outFileFormat", "bigwig",
            "-p", pthread,
            "--normalizeUsing", "CPM",
            "--binSize", "1" ]

    print("\t$", " ".join(cmd8_args))
    print(run_timed(cmd8_args))
    
######## VERSIONS #############
print("\n>>> VERSIONS:")

print("\n>>> FASEQC VERSION:")
print(subprocess.run(["fastqc", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE))

print("\n>>> HISAT2 VERSION:")
print(subprocess.run(["hisat2", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE))

print("\n>>> SAMTOOLS VERSION:")
print(subprocess.run(["samtools", "--version"],  stdout=subprocess.PIPE, stderr=subprocess.PIPE))

print("\n>>> HTSEQ-COUNT VERSION:")
output = subprocess.run(['htseq-count', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
for line in output.stdout:
    if line.find("version") >= 0:
        print(line.strip()
        break

print("\n>>> BAMCOVERAGE VERSION:")
print(subprocess(["bamCoverage", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE))
print(">>> END: Analayzer complete.")
