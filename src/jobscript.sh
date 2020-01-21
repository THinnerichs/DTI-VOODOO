#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J MafftHmmBuild
#SBATCH -o jobscript_outputs/MafftHmmBuild.%J.out
#SBATCH -e jobscript_outputs/MafftHmmBuild.%J.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=120G
#SBATCH --constraint=[intel]
#SBATCH --cpus-per-task=20

#run the application:
module load anaconda3/4.4.0
module load hmmer/3.2.1
source /home/${USER}/.bashrc
module load gcc/6.4.0
conda activate ~/.conda/envs/dti/

python3 alignment_pipeline.py
