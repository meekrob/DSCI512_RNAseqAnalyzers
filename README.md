# DSCI512_RNAseqAnalyzers
A simple set of wrappers and tools for RNA-seq analysis. These tools were designed for the DSCI512 RNA-seq analysis class

Below is a tutorial on the use of these scripts:
----


## Let's download the script templates I've written on github.

You will be able to tailor these templates to your own purposes.

DSCI512_RNAseqAnalyzers

**Exercise**

  * Locate the green **clone or download** button on the top right of this page. Click it.
  * Click on the clipboard icon. This will save a github URL address to your clipboard.
  * Switch over to summit.
  * Navigate into your folder for `PROJ04_GomezOrte2/02_scripts` and use `git clone` to pull the information from github to your location on summit.
  
```bash
$ cd /scratch/summit/<eID>@colostate.edu    #Replace <eID> with your EID
$ cd PROJ04_GomezOrte2/02_scripts
$ git clone <paste path to github repository here>
```

**Explore what you obtained.**


Notice that instead of having a single script, you now have a few scripts. These will work in a **Two step** method for executing jobs on summit. The `execute` script calls the `analyze` script. 

To execute the pipeline, you would do the following:

```bash
$ sbatch execute_RNAseq_pipeline.sbatch
```

By doing this, the `execute` script will start the `analyze` script by calling the following lines of code:

```bash

## execute the RNA-seq_pipeline
bash RNAseq_analyzer_191204.sh ../../01_input/metadata_gomezOrte_subset.txt $SLURM_NTASKS
```

**Usage:** `bash RNAseq_analyzer_191204.sh <metadatafile.txt> <number of threads>`
   *  Make sure the metadata file is correct.
   *  $SLURM_NTASKS automatically pulls the number of threads you have requested in the #SBATCH header.

**Exercise**
  * Open execute_RNAseq_pipeline.sbatch in an editor window
  * Tailor your execute_RNAseq_pipeline.sbatch script to suit your own inputs
  * Check the sbatch preamble to make sure everything appears to be working correctly
  * Include the proper path to your metadata file

  
**Answer**

  * Your execute script should look like this: 

```bash
#!/usr/bin/env bash

#SBATCH --job-name=test_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=6      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas  #modify this to reflect which queue you want to use. Options are 'shas' and 'shas-testing'
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=2:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_RNAseq_pipe_%j.txt


  ######### CHANGE <eID> TO YOUR EID: ############

## execute the RNA-seq_pipeline
bash RNAseq_analyzer_191204.sh ../../01_input/metadata_gomezOrte_subset.txt $SLURM_NTASKS
   ######### MODIFY the SECOND argument to point to YOUR metadata.file ######### 

## clean up by zipping .fastq files and deleting extra files
#bash RNAseq_cleanup_191204.sh ../../01_input/metadata_gomezOrte_subset.txt 
   ######### modify the SECOND argument to point to YOUR metadata.file ######### 
```

-----

## Let's navigate the analyze script

**Exercise**
  * Open RNAseq_analyze_191204.sh in an editor window
  * Tailor your script to suit your own inputs
  * There is a section ####### MODIFY THIS SECTION ############# Modify this section to match your requests
  * Add a relative path to your input directory
  * Add an absolute path to your hisat2path (needs to be int he format /path/to/hisat2indices/prefix)
  * Add an absolute path to your full genome fasta file (chrom.fa.gz)
  * Add the relative path to your annotation file.
  * Don't change the default output directory
  
  * Now run the script with the following command:
  
```bash
$ sbatch --reservation=csuhpc_dec05 execute_RNAseq_pipeline.sbatch
```

  * Note, don't use the option --reservation=csuhpc_dec05 when you are doing homework or working outside classtime

  * Check on the script with
  
```bash
$ squeue -u $USER
```

**Answer**

  * It should look something like this with your eID in place of eID
```bash


####### MODIFY THIS SECTION #############

#The input samples (metadata file and _fastq.gz files) live in directory:
inputdir="../../01_input/"

#This is where the ht2 files live:
hisat2path="/scratch/summit/eID@colostate.edu/DSCI512_RNAseq/PROJ01_ce11Build/ce11"

#This is where the genome sequence lives:
genomefa="/scratch/summit/eID@colostate.edu/DSCI512_RNAseq/PROJ01_ce11Build/chromFa.tar.gz"

#This is where the gtf file lives:
gtffile="../../01_input/ce11_annotation_ensembl_to_ucsc.gtf.gz"

#This is the output_directory:
DATE=`date +%Y-%m-%d`
#OR
#DATE='2018-10-16'
outputdir="../../03_output/"$DATE"_output/"


```


The script typically will take this long:

18 samples, 16 tasks = 38 min
9 samples, 16 tasks = 24 min
6 samples, 6 tasks = 


-----

## It's cooking show time! Let's see what we baked!

```bash
$ tree
.
├── 02_fastp
│   ├── EG01
│   │   ├── EG01_report.html
│   │   ├── EG01_report.json
│   │   ├── EG01_trim_1.fastq
│   │   └── EG01_trim_2.fastq
│   ├── EG02
│   │   ├── EG02_report.html
│   │   ├── EG02_report.json
│   │   ├── EG02_trim_1.fastq
│   │   └── EG02_trim_2.fastq
│   ├── EG03
│   │   ├── EG03_report.html
│   │   ├── EG03_report.json
│   │   ├── EG03_trim_1.fastq
│   │   └── EG03_trim_2.fastq
│   ├── EG04
│   │   ├── EG04_report.html
│   │   ├── EG04_report.json
│   │   ├── EG04_trim_1.fastq
│   │   └── EG04_trim_2.fastq
│   ├── EG05
│   │   ├── EG05_report.html
│   │   ├── EG05_report.json
│   │   ├── EG05_trim_1.fastq
│   │   └── EG05_trim_2.fastq
│   └── EG06
│       ├── EG06_report.html
│       ├── EG06_report.json
│       ├── EG06_trim_1.fastq
│       └── EG06_trim_2.fastq
├── 03_hisat2
│   ├── EG01.sam
│   ├── EG01_summary.txt
│   ├── EG02.sam
│   ├── EG02_summary.txt
│   ├── EG03.sam
│   ├── EG03_summary.txt
│   ├── EG04.sam
│   ├── EG04_summary.txt
│   ├── EG05.sam
│   ├── EG05_summary.txt
│   ├── EG06.sam
│   └── EG06_summary.txt
├── 04_feature
│   ├── counts.txt
│   └── counts.txt.summary
└── 05_samtools
    ├── EG01.bam
    ├── EG01_sort.bam
    ├── EG01_sort.bam.bai
    ├── EG01_sort.bw
    ├── EG02.bam
    ├── EG02_sort.bam
    ├── EG02_sort.bam.bai
    ├── EG02_sort.bw
    ├── EG03.bam
    ├── EG03_sort.bam
    ├── EG03_sort.bam.bai
    ├── EG03_sort.bw
    ├── EG04.bam
    ├── EG04_sort.bam
    ├── EG04_sort.bam.bai
    ├── EG04_sort.bw
    ├── EG05.bam
    ├── EG05_sort.bam
    ├── EG05_sort.bam.bai
    ├── EG05_sort.bw
    ├── EG06.bam
    ├── EG06_sort.bam
    ├── EG06_sort.bam.bai
    └── EG06_sort.bw

10 directories, 62 files
```

[Visualizing data using IGV](http://rna.colostate.edu/2019/doku.php?id=wiki:igv_visualization)

