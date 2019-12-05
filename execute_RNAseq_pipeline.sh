#!/usr/bin/env bash

#SBATCH --job-name=test_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=6      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas-testing  #modify this to reflect which queue you want to use. Options are 'shas' and 'shas-testing'
#SBATCH --qos=testing     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=0:29:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_RNAseq_pipe_%j.txt



## Source software
module load singularity
module list

fastp='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif fastp'
hisat2='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif hisat2'
featureCounts='singularity exec /projects/dcking@colostate.edu/containers/Summit_RNAseq_container.sif featureCounts'


  ######### CHANGE <eID> TO YOUR EID: ############

## execute the RNA-seq_pipeline
bash RNAseq_analyzer_181117.sh ../01_input/metadata_aceticAcid_subset.txt $SLURM_NTASKS
   ######### MODIFY the SECOND argument to point to YOUR metadata.file ######### 

## OR, you can use a python script
#python RNAseq_analyzer_181011.py ../01_input/metadata_aceticAcid_subset.txt $SLURM_NTASKS


## clean up by zipping .fastq files and deleting extra files
#bash RNAseq_cleanup_mouse_181011.sh ../04_testing/metadata_mouse.txt
   ######### modify the SECOND argument to point to YOUR metadata.file ######### 
