#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J MultiSequenceAlignment
#SBATCH -o MultiSequenceAlignment.%J.out
#SBATCH -e MultiSequenceAlignment.%J.err
#SBATCH --time=48:00:00
#SBATCH --mem=123G
#SBATCH --constraint=[intel]
#SBATCH --cpus-per-task=16

#run the application:
module load anaconda3/4.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/

python3 alignment_pipeline.py