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

    print(list(set(drug_list) - set(drug_mapping.keys()))[:100])

    raise Exception


    pos_sanity_list = [tup for tup in pred_list if round(tup[2] == 1)]
    random.shuffle(pos_sanity_list)
    print('sanity_check', pos_sanity_list[:200])


    pred_list = [tup for tup in pred_list if round(tup[2])==0 and tup[3] > 0.7]

    pred_list = sorted(pred_list, key=lambda tup: tup[3], reverse=True)
    print('len(pred_list):', len(pred_list))

    filename = '../results/full_model_with_mol_feat_results/best_preds_false_negatives'
    with open(file=filename, mode='w') as f:
        print('drug\tprot\tconfidence\tdrug_alias\tprotein_alias', file=f)
        for drug, protein, _, confidence in pred_list:
            if drug in drug_mapping.keys() and protein in protein_mapping.keys():
                print('\t'.join([drug, protein, str(confidence), drug_mapping[drug], protein_mapping[protein]]), file=f)





