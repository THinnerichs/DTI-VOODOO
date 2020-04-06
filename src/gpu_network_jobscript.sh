#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J GCNNet
#SBATCH -o jobscript_outputs/GCNNet.%J.out
#SBATCH -e jobscript_outputs/GCNNet.%J.err
#SBATCH --time=3-00:00:00
#SBATCH --gres=gpu:v100:4
#SBATCH --mem=300G
#SBATCH --constraint=[gpu]

#run the application:
module load anaconda3/4.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/

module load cuda/10.0.130


export CUDA_VISIBLE_DEVICES=0,1,2,3
python3 torch_dti_predictor.py --num_proteins -1 --num_epochs=50 --batch_size=32 --num_folds 5
