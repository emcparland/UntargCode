#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=meta
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=emcparland@whoi.edu
#SBATCH --ntasks=1
#SBATCH --mem=4000
#SBATCH --time=00:20:00
#SBATCH --output=meta.out
#SBATCH --error=meta.error

dir="$HOME/untarg_xcms"
cd "$dir"

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate untargmetab

Rscript --no-save --no-restore --verbose $dir/scripts_dir/create_metadata.R "pos" > $dir/$SLURM_JOB_NAME.Rout 2>&1 >$dir/$SLURM_JOB_NAME.log