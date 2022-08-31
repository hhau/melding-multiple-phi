#!/bin/bash

## Example SLURM script for BSU skylake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J aam71-fullmeld

## Enter the wall-clock time limit for for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=34:00:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 32.
#SBATCH --cpus-per-task=30

## Each task is allocated 5980M (skylake) or 12030M (skylake-himem). 
## If this is insufficient, uncomment and edit this line.
## Maximum value 191360M or 384960M.
## #SBATCH --mem=24G

## The system can send emails when your job starts and stops.
## Values include BEGIN, END, ALL, and TIME_LIMIT_80 and TIME_LIMIT_90 
## (reaching 80% or 90% of time limit.) Specify ARRAY_TASKS to receive
## a separate mail for each task. Multiple values can be given, separated by a comma.
#SBATCH --mail-type=ALL

## The project account name.
## Use mrc-bsu-sl2-cpu for skylake and mrc-bsu-sl2-gpu for pascal
#SBATCH -A mrc-bsu-sl2-cpu

## The partition. Use skylake for normal jobs, or skylake-himem if needed.
#SBATCH -p skylake

## GPU jobs only:
## Uncomment and specify the number of GPUs required per node, maximum 4.
## Note that there is a maximum of 3 cores per GPU.
## #SBATCH --gres=gpu:4

## Array jobs:
## #SBATCH --array=1-100

##  - - - - - - - - - - - - - -

## Section 2: Modules

# All scripts should include the first three lines.

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

# Load the latest R version.
# Before running your code, you should run R and install any required packages.
module load r-3.6.0-gcc-5.4.0-bzuuksv
module load jags-4.3.0-gcc-5.4.0-4z5shby

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel7/default-gpu

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

# You can use any arbitrary set of Linux commands here

CMD="make"

# Or for example:
# CMD="Rscript myScript.R"


###############################################################
### You should not have to change anything below this line ####
###############################################################

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then
        echo "Running on nodes: $SLURM_JOB_NODELIST"
else
        echo "Running on node: `hostname`"
fi

echo "Current directory: `pwd`"
echo -e "\nNum tasks = $SLURM_NTASKS, Num nodes = $SLURM_JOB_NUM_NODES, OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
