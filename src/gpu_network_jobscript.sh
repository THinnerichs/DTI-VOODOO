#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J MolPred
#SBATCH -o jobscript_outputs/MolPred.%J.out
#SBATCH -e jobscript_outputs/MolPred.%J.err
#SBATCH --time=1-00:00:00
#SBATCH --gres=gpu:v100:4
#SBATCH --mem=600G
#SBATCH --constraint=[gpu]
#SBATCH --cpus-per-gpu=6
# SBATCH --sockets-per-node=1
# SBATCH --gpus-per-socket=4

#run the application:
module load anaconda3/4.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/

module load cuda/10.0.130


# export CUDA_VISIBLE_DEVICES=0,1,2,3
# python3 torch_dti_predictor.py --num_proteins -1 --num_epochs=50 --batch_size=32 --num_folds 5

python3 molecular_predictor.py --batch_size 32384 --num_epochs 1 --lr 0.001 --model_id quick
