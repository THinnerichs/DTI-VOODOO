#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J HPOPred
#SBATCH -o jobscript_outputs/HPOPred.%J.out
#SBATCH -e jobscript_outputs/HPOPred.%J.err
#SBATCH --time=4:00:00
#SBATCH --gres=gpu:v100:4
#SBATCH --mem=300G
#SBATCH --constraint=[gpu]
#SBATCH --cpus-per-task=24

#run the application:
module load anaconda3/4.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/

module load gcc/8.2.0
module load cuda/10.2.89

python3 HPO_predictor.py --num_epochs 100 --batch_size 1000000 --lr 0.0001 --fold 3 --biosnap_test
python3 torch_dti_predictor.py --split_mode quick_standard --arch GENConv --num_epochs 1000 --batch_size 250 --fold 3 --lr 0.0001 --biosnap_test


# export CUDA_VISIBLE_DEVICES=0,1,2,3
# python3 torch_dti_predictor.py --num_proteins -1 --num_epochs=50 --batch_size=32 --num_folds 5

# python3 molecular_predictor.py --batch_size 32384 --num_epochs 1 --lr 0.001 --model_id quick
