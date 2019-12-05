#!/usr/bin/env bash

#SBATCH --job-name=test_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=16      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas  #modify this to reflect which queue you want to use. Options are 'shas' and 'shas-testing'
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=2:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_RNAseq_pipe_%j.txt


  ######### CHANGE <eID> TO YOUR EID: ############

## execute the RNA-seq_pipeline
#bash RNAseq_analyzer_191204.sh ../../01_input/metadata_gomezOrte.txt $SLURM_NTASKS
   ######### MODIFY the SECOND argument to point to YOUR metadata.file ######### 

## clean up by zipping .fastq files and deleting extra files
bash RNAseq_cleanup_191204.sh ../../01_input/metadata_gomezOrte.txt 
   ######### modify the SECOND argument to point to YOUR metadata.file ######### 
