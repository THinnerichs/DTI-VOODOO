#!/bin/bash
#SBATCH -N 100
#SBATCH --partition=batch
#SBATCH -J MultiSequenceAlignment
#SBATCH -o MultiSequenceAlignment.%J.out
#SBATCH -e MultiSequenceAlignment.%J.err
#SBATCH --time=03-00:00:00
#SBATCH --mem=1T
#SBATCH --constraint=[intel]&[local_1T]
#SBATCH --cpus-per-task=160

#run the application:
module load anaconda3/4.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/

python3 alignment_pipeline.py
