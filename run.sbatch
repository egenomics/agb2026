#!/bin/bash
#SBATCH --job-name=abgtemplate_integ    # Name of your master job (Group A -> B integration)
#SBATCH --output=nextflow_master_%j.log # Standard output log (%j inserts job ID)
#SBATCH --error=nextflow_master_%j.err  # Error log file
#SBATCH --ntasks=1                      # Nextflow only needs 1 task to orchestrate
#SBATCH --cpus-per-task=2               # 2 CPUs for the Java engine is plenty
#SBATCH --mem=24G                       # RAM for the Nextflow orchestrator process
#SBATCH --time=04:00:00                 # full 626-sample run; projected ~20–60 min
#SBATCH --partition=std                 # infinite walltime, default partition
#SBATCH --qos=normal                    # normal QOS allows >2h (express QOS capped lower)
#SBATCH --account=upfagb_normal         # account that grants the normal QOS

bash run.sh
