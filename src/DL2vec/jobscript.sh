#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J DL2vec
#SBATCH -o ../jobscript_outputs/DL2vec.%J.out
#SBATCH -e ../jobscript_outputs/DL2vec.%J.err
#SBATCH --time=0-24:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=64

#run the application:
source /home/hinnertr/.bashrc
module load gcc/6.4.0
conda activate ~/.conda/envs/dti/
module load groovy/3.0.6

python runDL2vec.py -embedsize 200 -ontology ../../data/PhenomeNET_data/phenomenet.owl -associations ../../data/PhenomeNET_data/drug_association_file -outfile ../../data/PhenomeNET_data/drug_embedding_model -entity_list ../../data/PhenomeNET_data/drug_entity_list -num_workers 64


