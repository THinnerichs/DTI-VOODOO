import pickle
import numpy as np
import random

import DTI_data_preparation

def write_predicted_DTIs(fold=3):
    filename = '../models/graph_models/PPI_network_model_with_mol_features_fold_'+str(fold)+'_predictions.pkl'

    drug_mapping = DTI_data_preparation.get_drug_to_name_mapping()
    protein_mapping = DTI_data_preparation.get_protein_to_EnsemblProtein_id()

    print('drug_dict.keys()', len(list(drug_mapping.keys())))
    print('protein_dict.keys()', len(list(protein_mapping.keys())))

    pred_list = None #[(drug, prot, label, pred)]
    with open(file=filename, mode='rb') as f:
        pred_list = pickle.load(f)

    drug_list = [tup[0] for tup in pred_list]
    protein_list = [tup[1] for tup in pred_list]

    print('drug_list', len(set(drug_list)))
    print('protein_list', len(set(protein_list)))
    print('intersection drug', len(set(drug_list) & set(drug_mapping.keys())))
    print('intersection protein', len(set(protein_list) & set(protein_mapping.keys())))

    # print(list(set(drug_list) - set(drug_mapping.keys()))[:100])

    # pos_sanity_list = [tup for tup in pred_list if round(tup[2] == 1)]
    # random.shuffle(pos_sanity_list)
    # print('sanity_check', pos_sanity_list[:200])


    pred_list = [tup for tup in pred_list if round(tup[2])==0 and tup[3] > 0.7]

    pred_list = sorted(pred_list, key=lambda tup: tup[3], reverse=True)
    print('len(pred_list):', len(pred_list))

    drug_list = list({drug for drug in drug_list if drug in drug_mapping.keys()})
    protein_list = list({protein for protein in protein_list if protein in protein_mapping.keys()})

    filename = '../data/PPI_data/protein_to_gene_dict.pkl'
    with open(file=filename, mode='rb') as f:
        protein_to_gene_mapping = pickle.load(f)

    given_gene_list = map(lambda prot: protein_to_gene_mapping[prot[5:]], protein_list)

    filename = '../data/PPI_data/all_samples_top10.txt'
    driver_gene_dict = {}
    with open(file=filename, mode='r') as f:
        f.readline()

        for line in f:
            split_line = line.strip().split('\t')
            gene = split_line[0]
            is_driver = bool(split_line[3])

            cancer_type, gene_name = split_line[1:3]

            if is_driver:
                driver_gene_dict[gene] = (cancer_type, gene_name)

    print('Num driver genes:', len(driver_gene_dict))

    protein_list = [prot for prot in protein_list if prot[5:] in protein_to_gene_mapping.keys() and
                                                     protein_to_gene_mapping[prot[5:]] in driver_gene_dict.keys()]

    print('Num proteins that are driver genes:', len(protein_list))

    filename = '../results/full_model_with_mol_feat_results/best_preds_false_negatives'
    with open(file=filename, mode='w') as f:
        print('drug\tprot\tconfidence\tdrug_alias\tprotein_alias\tcancer_type\tgene_name', file=f)
        for drug, protein, _, confidence in pred_list:
            if drug in drug_list and protein in protein_list:
                print('\t'.join([drug, protein, str(confidence), drug_mapping[drug], protein_mapping[protein]]),
                      str(driver_gene_dict[protein_to_gene_mapping[protein[5:]]])[1:-1], file=f)

if __name__ == '__main__':
    write_predicted_DTIs()




