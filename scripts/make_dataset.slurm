#!/bin/bash
#SBATCH -J IDRBlocks
#SBATCH --mail-user npradeep96
#SBATCH -p sched_mit_arupc
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu 4000
cd "$SLURM_SUBMIT_DIR"
echo $PWD

run_program()
{
    source activate IDRBlocks
    python make_dataset.py --i $idr_file --p $pr_file --o $out_file
    conda deactivate
}

run_program