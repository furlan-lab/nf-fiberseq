#!/bin/bash
#SBATCH --job-name=NF-main
#SBATCH --output=nf%j%.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=tho3@fredhutch.org

# Ensure that the script fails if any of the commands fail
set -euo pipefail

# Nextflow Version
NXF_VER=24.04.3

# Nextflow Configuration File
NXF_CONFIG=/fh/fast/furlan_s/user/tho/nextflow/nextflow.config
NXF_APPTAINER_CACHEDIR=/fh/fast/furlan_s/user/tho/nextflow/apptainer_cache

# Workflow to Run (e.g. GitHub repository)
# use with nextflow -r <branch> <repository>
# WORKFLOW_REPO=FredHutch/workflow-repository


ml Nextflow
ml Apptainer

# Make sure that the apptainer executables are in the PATH
export PATH=$APPTAINERROOT/bin/:$PATH

# Run the workflow
NXF_VER=$NXF_VER \
nextflow \
    run main.nf \
    -c ${NXF_CONFIG} -profile slurm \
    -with-report nextflow.report.html \
    -resume
