#!/usr/bin/env bash

#SBATCH --job-name=test_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=6       # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas # keep. No shas-testing with array jobs.
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=2:00:00   # REDUCE THE TIME ALLOTMENT FOR WHAT YOU WOULD EXPECT RUNNING A SINGLE LINE OF THE METADATA FILE
#SBATCH --output=log_RNAseq_pipe_%j_%a.txt

METAFILE=$1

## Use the following example to submit THIS script (array_RNAseq_pipeline.sbatch):
## sbatch --array 1-J array_RNAseq_pipeline.sbatch metadata.txt

# Uncomment this line to run the analysis
#bash RNAseq_analyzer_191204.sh <(head -n $SLURM_ARRAY_TASK_ID $METAFILE | tail -1) $SLURM_NTASKS

## Uncomment this line to clean up by zipping .fastq files and deleting extra files
#bash RNAseq_cleanup_191204.sh <(head -n $SLURM_ARRAY_TASK_ID $METAFILE | tail -1)
