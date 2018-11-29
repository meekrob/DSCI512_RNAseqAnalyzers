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

If you use 24 cores on 1 node and run it on the testing queue, the script should take 26 minutes to run.

-----

## It's cooking show time! Let's see what we baked!

