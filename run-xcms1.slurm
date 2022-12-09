#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=xcms1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emcparland@whoi.edu
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=12:00:00
#SBATCH --output=xcms1_%A-%a.out
#SBATCH --error=xcms1_%A-%a.error
#SBATCH --array=1-514%40

dir="$HOME/untarg_xcms/mzML_MSn2022"
cd "$dir"

in_dir="/vortexfs1/omics/kujawinski/data/biosscope_untarg_mzML_MSn"
out_dir="$dir/output_dir/$SLURM_JOB_NAME"
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

echo `date`
echo This is job $SLURM_JOB_ID
echo This is task $SLURM_ARRAY_TASK_ID

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate untargmetabR4

Rscript --no-save --no-restore --verbose $dir/scripts_dir/xcms1.R "$in_dir" "pos" "$SLURM_ARRAY_TASK_ID" "$out_dir" > "$out_dir/${SLURM_JOB_NAME}_pos_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout" 2>&1 > "$out_dir/${SLURM_JOB_NAME}_pos_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

echo `date`
