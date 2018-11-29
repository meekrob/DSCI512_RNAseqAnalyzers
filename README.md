# DSCI512_RNAseqAnalyzers
A simple set of wrappers and tools for RNA-seq analysis. These tools were designed for the DSCI512 RNA-seq analysis class


----


## Let's download the script templates I've written on github.

You will be able to tailor these templates to your own purposes.

DSCI512_RNAseqAnalyzers

**Exercise**

  * Go to the github location.
  * Copy the github link by clicking on the green clone or download button
  * Click on the clipboard icon.
  * Switch over to summit.
  * Navigate into /scratch/summit/<eID>@colostate.edu
  
```bash
$ cd /scratch/summit/<eID>@colostate.edu    #Replace <eID> with your EID
$ git clone <paste path to github repository here>
```

**Explore what you obtained.**

**Copy and paste the template codes over to your `02_scripts` directory**

```bash
$ cp *.sh ../PROJ06_yeastDemo2/02_scripts
```

Notice that instead of having a single script, you now have two scripts. The `execute` script calls the `analyze` script. 

To execute the pipeline, you would do the following:

```bash
$ sbatch execute_RNAseq_pipeline.sh
```

By doing this, the `execute` script would start the `analyze` script by calling the following lines of code:

```bash
##
#source /scratch/summit/erinnish@colostate.edu/activate.bashrc
source /projects/dcking@colostate.edu/paths.bashrc

## execute the RNA-seq_pipeline
bash RNAseq_analyzer_181117.sh ../01_input/metadata_aceticAcid_subset.txt $SLURM_NTASKS
```

**Usage:** `bash RNAseq_analyzer_181117.sh <metadatafile.txt> <number of threads>`
   *  Make sure the metadata file is correct.
   * `$SLURM_NTASKS` automatically pulls the number of threads you have requested in the #SBATCH header.

**Exercise**
  * Tailor your execute_RNAseq_pipeline.sh script if you need to.

-----

## Let's navigate the analyze script

There is a section ####### MODIFY THIS SECTION #############

**Exercise** 

  * Modify this section to match your requests
  
**Exercise**
  
  * Now run the script with the following command:
  
```bash
$ sbatch execute_RNAseq_pipeline.sh
```

  * Check on the script with
  
```bash
$ scheck   #This is an alias of squeue -u $USER
$ squeue -u $USER --start
```

If you use 24 cores on 1 node and run it on the testing queue, the script should take 27 minutes to run.

-----

## It's cooking show time! Let's see what we baked!

```bash
.
├── 02_fastqc
│   ├── sample01_1
│   │   ├── SRR3567551_1_fastqc.html
│   │   └── SRR3567551_1_fastqc.zip
│   ├── sample01_2
│   │   ├── SRR3567551_2_fastqc.html
│   │   └── SRR3567551_2_fastqc.zip
│   ├── sample02_1
│   │   ├── SRR3567552_1_fastqc.html
│   │   └── SRR3567552_1_fastqc.zip
│   ├── sample02_2
│   │   ├── SRR3567552_2_fastqc.html
│   │   └── SRR3567552_2_fastqc.zip
│   ├── sample03_1
│   │   ├── SRR3567554_1_fastqc.html
│   │   └── SRR3567554_1_fastqc.zip
│   ├── sample03_2
│   │   ├── SRR3567554_2_fastqc.html
│   │   └── SRR3567554_2_fastqc.zip
│   ├── sample04_1
│   │   ├── SRR3567555_1_fastqc.html
│   │   └── SRR3567555_1_fastqc.zip
│   ├── sample04_2
│   │   ├── SRR3567555_2_fastqc.html
│   │   └── SRR3567555_2_fastqc.zip
│   ├── sample09_1
│   │   ├── SRR3567674_1_fastqc.html
│   │   └── SRR3567674_1_fastqc.zip
│   ├── sample09_2
│   │   ├── SRR3567674_2_fastqc.html
│   │   └── SRR3567674_2_fastqc.zip
│   ├── sample10_1
│   │   ├── SRR3567676_1_fastqc.html
│   │   └── SRR3567676_1_fastqc.zip
│   ├── sample10_2
│   │   ├── SRR3567676_2_fastqc.html
│   │   └── SRR3567676_2_fastqc.zip
│   ├── sample11_1
│   │   ├── SRR3567677_1_fastqc.html
│   │   └── SRR3567677_1_fastqc.zip
│   ├── sample11_2
│   │   ├── SRR3567677_2_fastqc.html
│   │   └── SRR3567677_2_fastqc.zip
│   ├── sample12_1
│   │   ├── SRR3567679_1_fastqc.html
│   │   └── SRR3567679_1_fastqc.zip
│   └── sample12_2
│       ├── SRR3567679_2_fastqc.html
│       └── SRR3567679_2_fastqc.zip
├── 03_hisat2
│   ├── sample01.sam
│   ├── sample01_summary.txt
│   ├── sample02.sam
│   ├── sample02_summary.txt
│   ├── sample03.sam
│   ├── sample03_summary.txt
│   ├── sample04.sam
│   ├── sample04_summary.txt
│   ├── sample09.sam
│   ├── sample09_summary.txt
│   ├── sample10.sam
│   ├── sample10_summary.txt
│   ├── sample11.sam
│   ├── sample11_summary.txt
│   ├── sample12.sam
│   └── sample12_summary.txt
├── 04_feature
│   ├── counts.txt
│   └── counts.txt.summary
└── 05_samtools
    ├── sample01.bam
    ├── sample01_sort.bam
    ├── sample01_sort.bam.bai
    ├── sample01_sort.bw
    ├── sample02.bam
    ├── sample02_sort.bam
    ├── sample02_sort.bam.bai
    ├── sample02_sort.bw
    ├── sample03.bam
    ├── sample03_sort.bam
    ├── sample03_sort.bam.bai
    ├── sample03_sort.bw
    ├── sample04.bam
    ├── sample04_sort.bam
    ├── sample04_sort.bam.bai
    ├── sample04_sort.bw
    ├── sample09.bam
    ├── sample09_sort.bam
    ├── sample09_sort.bam.bai
    ├── sample09_sort.bw
    ├── sample10.bam
    ├── sample10_sort.bam
    ├── sample10_sort.bam.bai
    ├── sample10_sort.bw
    ├── sample11.bam
    ├── sample11_sort.bam
    ├── sample11_sort.bam.bai
    ├── sample11_sort.bw
    ├── sample12.bam
    ├── sample12_sort.bam
    ├── sample12_sort.bam.bai
    └── sample12_sort.bw
    
20 directories, 82 files
```



