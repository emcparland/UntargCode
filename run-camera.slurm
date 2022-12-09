#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=camera
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emcparland@whoi.edu
#SBATCH --ntasks-per-node=12
#SBATCH --mem=20gb
#SBATCH --time=16:00:00
#SBATCH --output=camera.out
#SBATCH --error=camera.error


dir="$HOME/untarg_xcms/mzML_MS2021/output_dir/xcms2"
cd "$dir"

work_dir="/output_dir/xcms2"

echo `date`
echo This is job $SLURM_JOB_ID

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate untargmetabR4

Rscript --no-save --no-restore --verbose $HOME/untarg_xcms/mzML_MS2021/scripts_dir/CAMERA.R "$work_dir"

echo `date`
