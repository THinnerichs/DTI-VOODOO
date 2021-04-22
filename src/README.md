## How to download the data
We provide a download script in `data/STITCH_data`,`data/STRING_data`, `data/PhenomeNET_data`. Run them in each directory with 
```
source download.sh
```
These files take about 30 GiB in total. Make sure you have enough disk space.

## How to run DTI-Voodoo data preparation pipeline
Please run the preprocessing pipeline with 
```
python DTI_data_preparation.py
```
The command will output 4 commands for DL2vec to execute. Please run each of them in `src/DL2vec/`. Each of them takes about 2 hours on 16 cores and 50GB RAM.

To run `DeepGOPlus`, switch to `results/protein_representation`, make sure to have `blastp` and `diamond` installed and run 
```
wget http://deepgoplus.bio2vec.net/data/data.tar.gz -P data/
```
and unpack the directory. If this is now available, use the more bloated `https://deepgo.cbrc.kaust.edu.sa/data/data-2016.tar.gz`
Then run
```
source predict.sh data/PPI_graph_protein_seqs.fasta results/output
```
Refer to [DeepGOPlus documentation](https://github.com/bio-ontology-research-group/deepgoplus) for more in-depth information.

Now run the following commands after installing the required packages to build the drug molecular features.
```
python molecular_predictor.py --compute_embeddings
```

## How to run DTI-Voodoo

To run DTI-Voodoo, we have to pretrain the `HPO_predictor` first. As each fold is taking significant time, we parallelized them by running each fold individualy in this 5-fold CV. 
To run on STITCH/Yamanishi/BioSnap substitute `{dataset}` with ` `(nothing)/`--yamanishi_test`/`--biosnap_test`
```
python HPO_predictor.py --num_epochs 100 --batch_size 10000 --lr 0.0001 --fold 3 --mode database {dataset}
```

To run full DTI-Voodoo after pretraining `HPO_predictor` with `GCNConv` run
```
python3 torch_dti_predictor.py --arch GCNConv --include_mol_features --num_epochs 200 --batch_size 50 --fold 3 --lr 0.0001 --PPI_min_score 700 --mode database {dataset}
```
or with `GENConv`
```
python3 torch_dti_predictor.py --arch GENConv --include_mol_features --num_epochs 200 --batch_size 50 --fold 3 --lr 0.0001 --PPI_min_score 700 --mode database {dataset}
```
with `{dataset}` set as described above.
