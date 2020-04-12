#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J drug_encoding
#SBATCH -o jobscript_outputs/drug_encoding.%J.out
#SBATCH -e jobscript_outputs/drug_encoding.%J.err
#SBATCH --time=0-06:00:00
#SBATCH --mem=300G
#SBATCH --gres=gpu:v100:2
#SBATCH --constraint=[gpu]
#SBATCH --cpus-per-task=24

#run the application:
module load anaconda3/4.4.0
# module load hmmer/3.2.1
source /home/${USER}/.bashrc
# module load gcc/6.4.0
conda activate ~/.conda/envs/dti/

# python3 alignment_pipeline.py
module load diamond/0.9.22
python3 molecular_predictor.py
