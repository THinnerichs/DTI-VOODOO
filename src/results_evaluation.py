import pickle
import numpy as np
import random

import DTI_data_preparation
import DDI_utils
import dti_utils


import sklearn.metrics as metrics

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

    yamanishi_drug_mapping = DDI_utils.get_Yamanishi_db_to_PubChem_mapping_dict()

    print('yamanishi mapping length', len(yamanishi_drug_mapping))
    print(list(yamanishi_drug_mapping.items())[:10])
    drug_list = [drug for drug in drug_list if drug in yamanishi_drug_mapping.values()]

    print('drugs that are also present in yamanishi dataset:', len(drug_list))

    filename = '../results/full_model_with_mol_feat_results/best_preds_false_negatives_yama_subset'
    with open(file=filename, mode='w') as f:
        print('drug\tprot\tconfidence\tdrug_alias\tprotein_alias\tcancer_type', file=f)
        for drug, protein, _, confidence in pred_list:
            if drug in drug_list and protein in protein_list:
                print('\t'.join([drug, protein, str(confidence), drug_mapping[drug], protein_mapping[protein]]),
                      str(driver_gene_dict[protein_to_gene_mapping[protein[5:]]][0]), file=f)


def analyze_DTI_net_results():
    path = '../data/Yamanishi_data/'
    dti_matrix = np.loadtxt(path+'mat_drug_protein.txt')
    print(dti_matrix.shape)
    zscores = np.loadtxt(path+'Zscore.txt', delimiter=',')
    print(zscores.shape)

    print((dti_matrix.sum(axis=0)==0).sum())

    raise Exception

    print(zscores.max(), zscores.min())


    zscores = zscores>0.0
    ord = np.flip(np.argsort(zscores.flatten()))
    label = dti_matrix.flatten()[ord]

    p = label.sum()
    n = len(label)-p

    TP = np.cumsum(label)
    PP = np.arange(len(label)) + 1

    rocx = (PP - TP) / n
    rocy = TP / p

    print(metrics.auc(rocx, rocy))

    raise Exception

    noise = np.random.randn(708, 1512) * 0.001
    zscores = zscores + noise

    unique_values = np.unique(zscores[((0.0125<zscores) * (zscores<0.0135))])

    for t in range(0, 100, 5):
    # for i in range(2000, len(unique_values)-3000, 10):
        # threshold = unique_values[i]
        threshold = 0.8 - t/100
        y_pred = zscores>=threshold
        print(threshold, y_pred.sum(), 'auc:', dti_utils.dti_auroc(dti_matrix.flatten(), y_pred.flatten()))

def test_Yamanishi_AUC():
    path = '../data/Yamanishi_data/'
    # parse dti matrix provided by https://github.com/luoyunan/DTINet/
    dti_matrix = []
    with open(file=path + 'mat_drug_protein.txt', mode='r') as f:
        for line in f:
            dti_matrix.append(list(map(int, list(line.strip().replace(' ', '')))))
    dti_matrix = np.array(dti_matrix)
    print('dti_matrix.shape', dti_matrix.shape)
    print(dti_matrix[:, 0].sum())
    print(dti_matrix.sum())

    num_drugs, num_prots = dti_matrix.shape

    sum_vec = dti_matrix.sum(axis=0)

    print('sumvec:', sum_vec.shape, sum_vec[:20])

    idx = sum_vec.argsort()

    idx = np.flip(idx)

    dti_matrix = dti_matrix.flatten()

    for i in range(200, num_prots+1, 1):
        y_pred = np.zeros((num_drugs, num_prots))
        y_pred[:, idx[:i]] = 1
        y_pred = y_pred.flatten()

        print(f'i: {i}, ROCAUC: {dti_utils.dti_auroc(y_true=dti_matrix, y_pred=y_pred)}, '
              f'microAUC: {dti_utils.micro_AUC_per_prot(y_true=dti_matrix, y_pred=y_pred, num_drugs=num_drugs)}, '
              f'{dti_utils.micro_AUC_per_drug(dti_matrix, y_pred,num_drugs)}')


if __name__ == '__main__':
    analyze_DTI_net_results()
    # write_predicted_DTIs()




