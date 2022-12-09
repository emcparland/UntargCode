#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=meta
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=emcparland@whoi.edu
#SBATCH --ntasks=1
#SBATCH --mem=4000
#SBATCH --time=00:10:00
#SBATCH --output=meta.out
#SBATCH --error=meta.error

dir="$HOME/untarg_xcms/mzML_MS2021"
cd "$dir"

in_dir="/vortexfs1/omics/kujawinski/data/biosscope_untarg_mzml"

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate untargmetabR4

Rscript --no-save --no-restore --verbose $dir/scripts_dir/create_metadata.R "$in_dir" "pos" > $dir/$SLURM_JOB_NAME.Rout 2>&1 >$dir/$SLURM_JOB_NAME.log
