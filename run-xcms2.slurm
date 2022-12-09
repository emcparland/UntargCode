#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --job-name=xcms2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emcparland@whoi.edu
#SBATCH --ntasks-per-node=36
#SBATCH --mem=750gb
#SBATCH --time=4:00:00
#SBATCH --output=xcms2.out
#SBATCH --error=xcms2.error

dir="$HOME/untarg_xcms/mzML_MSn2022"
cd "$dir"

in_dir="$HOME/untarg_xcms/mzML_MSn2022/output_dir/xcms1"
out_dir="$dir/output_dir/$SLURM_JOB_NAME"

if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

echo `date`
echo This is job $SLURM_JOB_ID

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate untargmetabR4

Rscript --no-save --no-restore --verbose $dir/scripts_dir/xcms2.R "$in_dir" "$out_dir" "neg" > "$out_dir/${SLURM_JOB_NAME}_neg.Rout" 2>&1 >$out_dir/${SLURM_JOB_NAME}_neg.log

