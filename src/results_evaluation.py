import pickle
import numpy as np
import random

from tqdm import tqdm

import DTI_data_preparation
import DDI_utils
import dti_utils

import sklearn.metrics as metrics
from sklearn.model_selection import KFold



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

    print(zscores.max(), zscores.min())

    num_drugs, num_prots = dti_matrix.shape


    y_pred = zscores>=0.5
    print(y_pred.sum(), 'auc:', dti_utils.dti_auroc(dti_matrix.flatten(), y_pred.flatten()),
          'microauc:', dti_utils.micro_AUC_per_prot(dti_matrix.flatten(), y_pred.flatten(), num_drugs=num_drugs),
          'microauc:', dti_utils.micro_AUC_per_drug(dti_matrix.flatten(), y_pred.flatten(), num_drugs=num_drugs))


def test_Yamanishi_AUC():
    '''
    path = '../data/Yamanishi_data/'
    # parse dti matrix provided by https://github.com/luoyunan/DTINet/

    dti_matrix = []
    with open(file=path + 'mat_drug_protein.txt', mode='r') as f:
        for line in f:
            dti_matrix.append(list(map(int, list(line.strip().replace(' ', '')))))
    dti_matrix = np.array(dti_matrix)
    '''

    path = '../data/STITCH_data/'
    with open(path+'y_dti.npy', mode='rb') as f:
        dti_matrix = np.load(f)

    '''
    path = '../results/DTIGEMS/'
    filename = path + 'ic_admat_dgc.txt'

    dti_matrix = np.genfromtxt(filename, delimiter='\t', dtype=str)
    dti_matrix = dti_matrix[1:, :][:, 1:]
    dti_matrix = dti_matrix.astype(np.int)
    '''
    print('dti_matrix.shape', dti_matrix.shape)

    print(dti_matrix.sum())

    print((dti_matrix.sum(axis=1)==0).sum())

    num_drugs, num_prots = dti_matrix.shape

    sum_vec = dti_matrix.sum(axis=0)

    print('sumvec:', sum_vec.shape, sum_vec[:20])

    idx = sum_vec.argsort()

    idx = np.flip(idx)

    dti_matrix = dti_matrix.flatten()

    for i in range(1650, 1750, 10):
        y_pred = np.zeros((num_drugs, num_prots))
        y_pred[:, idx[:i]] = 1
        y_pred = y_pred.flatten()

        print(f'i: {i}, ROCAUC: {dti_utils.dti_auroc(y_true=dti_matrix, y_pred=y_pred)}, '
              f'microAUC: {dti_utils.micro_AUC_per_prot(y_true=dti_matrix, y_pred=y_pred, num_drugs=num_drugs)}, '
              f'{dti_utils.micro_AUC_per_drug(dti_matrix, y_pred,num_drugs)}')


    return

    kf = KFold(n_splits=5, random_state=1, shuffle=True)
    # X = np.zeros((num_drugs, 1))
    X = np.arange(num_drugs * num_prots)

    train_aucs = []
    train_micro_aucs = []
    test_aucs = []
    test_micro_aucs = []

    '''
    for train_drug_indices, test_drug_indices in kf.split(X):
        print(train_drug_indices.shape, test_drug_indices.shape)
        y_train = dti_matrix[train_drug_indices, :]
        y_test = dti_matrix[test_drug_indices, :]
    '''
    for train_indices, test_indices in kf.split(X):
        print(train_indices.shape, test_indices.shape)
        y_train = np.zeros(num_drugs * num_prots)
        y_train[train_indices] = dti_matrix.flatten()[train_indices]
        y_train = y_train.reshape((num_drugs, num_prots))

        y_test = np.zeros(num_drugs * num_prots)
        y_test[test_indices] = dti_matrix.flatten()[test_indices]


        sum_vec = y_train.sum(axis=0)

        print('sumvec:', sum_vec.shape, sum_vec[:20])
        idx = sum_vec.argsort()

        idx = np.flip(idx)
        y_train = y_train.flatten()
        y_test = y_test.flatten()

        max_train_auroc = max_train_micro_auc = 0
        max_test_auroc = max_test_micro_auc = 0
        for i in tqdm(range(0, int(num_prots / 2), 2)):
            y_pred = np.zeros((num_drugs, num_prots))
            y_pred[:, idx[:i]] = 1
            y_pred = y_pred.flatten()

            # help_val = dti_utils.dti_auroc(y_true=y_train.flatten(), y_pred=y_pred[train_drug_indices, :].flatten())
            help_val = dti_utils.dti_auroc(y_true=y_train[train_indices], y_pred=y_pred[train_indices])
            max_train_auroc = help_val if help_val > max_train_auroc else max_train_auroc

            # help_val = dti_utils.micro_AUC_per_prot(y_true=y_train, y_pred=y_pred[train_indices], num_drugs=len(train_drug_indices))
            # max_train_micro_auc = help_val if help_val > max_train_micro_auc else max_train_micro_auc

            help_val = dti_utils.dti_auroc(y_true=y_test[test_indices], y_pred=y_pred[test_indices].flatten())
            max_test_auroc = help_val if help_val > max_test_auroc else max_test_auroc

            # help_val = dti_utils.micro_AUC_per_prot(y_true=y_test, y_pred=y_pred[test_drug_indices, :].flatten(), num_drugs=len(test_drug_indices))
            # max_test_micro_auc = help_val if help_val > max_test_micro_auc else max_test_micro_auc

        print(f'Train: ROCAUC {max_train_auroc}, '
              # f'microAUC: {max_train_micro_auc}, '
              f'Test: ROCAUC {max_test_auroc}, '
              # f'microAUC: {max_test_micro_auc}'
              )

        train_aucs.append(max_train_auroc)
        # train_micro_aucs.append(max_train_micro_auc)
        test_aucs.append(max_test_auroc)
        # test_micro_aucs.append(max_test_micro_auc)

    for results in [train_aucs, train_micro_aucs, test_aucs, test_micro_aucs]:
        results = np.array(results)
        print(results.mean())

def parse_BioSnap():
    path = 'data/BioSnap_data/'

    train_file = 'train.csv'
    val_file = 'val.csv'
    test_file = 'test.csv'

    full_list = []
    for file in [train_file, val_file, test_file]:
        print(f'Parsing {path+file}')
        with open(file=path+file, mode='r') as f:
            f.readline()
            for i, line in enumerate(f):
                full_list.append(line.strip().split(','))
    full_table = np.array(full_list)

    '''
    train_table = np.genfromtxt(path+train_file, delimiter=',', dtype=np.str)[1:,:]
    val_table = np.loadtxt(path+val_file, delimiter=',', dtype=np.str)[1:,:]
    test_table = np.loadtxt(path+test_file, delimiter=',', dtype=np.str)[1,:]

    full_table = np.vstack((train_table, val_table, test_table))
    '''

    DB_drug_list = list(set(full_table[:, 3]))
    gene_list = list(set(full_table[:,4]))

    with open(file=path+'BioSnap_DB_drug_list', mode='w') as f:
        for drug in DB_drug_list:
            print(drug, file=f)
    with open(file=path+'BioSnap_gene_list', mode='w') as f:
        for gene in gene_list:
            print(gene, file=f)

    print('num drugs, num genes:', len(DB_drug_list), len(gene_list))

    dti_matrix = np.zeros((len(DB_drug_list), len(gene_list)))

    for i in tqdm(range(full_table.shape[0])):
        drug, gene, label = full_table[i,3:6]

        drug_index = DB_drug_list.index(drug)
        gene_index = gene_list.index(gene)
        label = float(label)

        dti_matrix[drug_index, gene_index] = label

    print(dti_matrix.sum(), dti_matrix.shape)

    print(f'Number of genes with at least one interactor {(dti_matrix.sum(axis=0)>0).sum()}')

    raise Exception

    naive_predictor_pair_split(dti_matrix=dti_matrix)

    # Naive predictor over DT pairs:
    num_drugs, num_prots = dti_matrix.shape

    sum_vec = dti_matrix.sum(axis=0)

    print('sumvec:', sum_vec.shape, sum_vec[:20])

    idx = sum_vec.argsort()

    idx = np.flip(idx)

    dti_matrix = dti_matrix.flatten()

    for i in range(0, num_prots, 10):
        y_pred = np.zeros((num_drugs, num_prots))
        y_pred[:, idx[:i]] = 1
        y_pred = y_pred.flatten()

        print(f'i: {i}, ROCAUC: {dti_utils.dti_auroc(y_true=dti_matrix, y_pred=y_pred)}, '
              f'microAUC: {dti_utils.micro_AUC_per_prot(y_true=dti_matrix, y_pred=y_pred, num_drugs=num_drugs)}, '
              f'{dti_utils.micro_AUC_per_drug(dti_matrix, y_pred, num_drugs)}')

def naive_predictor_pair_split(dti_matrix):
    # naive predictor on DT pair data
    # note that the specific micro_AUC computation needed, is extremely slow.

    kf = KFold(n_splits=5, random_state=42, shuffle=False)

    num_drugs, num_prots = dti_matrix.shape
    X = np.arange(num_drugs * num_prots)

    train_aucs = []
    train_micro_aucs = []
    test_aucs = []
    test_micro_aucs = []

    for train_indices, test_indices in kf.split(X):
        print(train_indices.shape, test_indices.shape)
        y_train = np.zeros(num_drugs * num_prots)
        y_train[train_indices] = dti_matrix.flatten()[train_indices]
        y_train = y_train.reshape((num_drugs, num_prots))

        y_test = np.zeros(num_drugs * num_prots)
        y_test[test_indices] = dti_matrix.flatten()[test_indices]

        sum_vec = y_train.sum(axis=0)

        idx = sum_vec.argsort()

        idx = np.flip(idx)
        y_train = y_train.flatten()
        y_test = y_test.flatten()

        max_train_auroc = max_train_micro_auc = 0
        max_test_auroc = max_test_micro_auc = 0
        for i in tqdm(range(0, int(num_prots / 2), 2)):
            y_pred = np.zeros((num_drugs, num_prots))
            y_pred[:, idx[:i]] = 1
            y_pred = y_pred.flatten()

            # help_val = dti_utils.dti_auroc(y_true=y_train.flatten(), y_pred=y_pred[train_drug_indices, :].flatten())
            help_val = dti_utils.dti_auroc(y_true=y_train[train_indices], y_pred=y_pred[train_indices])
            max_train_auroc = help_val if help_val > max_train_auroc else max_train_auroc

            help_val = dti_utils.micro_AUC_per_prot_DT_pairs(y_true=y_train, y_pred=y_pred, num_drugs=num_drugs, indices=train_indices)
            max_train_micro_auc = help_val if help_val > max_train_micro_auc else max_train_micro_auc

            help_val = dti_utils.dti_auroc(y_true=y_test[test_indices], y_pred=y_pred[test_indices].flatten())
            max_test_auroc = help_val if help_val > max_test_auroc else max_test_auroc

            help_val = dti_utils.micro_AUC_per_prot_DT_pairs(y_true=y_test, y_pred=y_pred, num_drugs=num_drugs, indices=test_indices)
            max_test_micro_auc = help_val if help_val > max_test_micro_auc else max_test_micro_auc

        print(f'Train: ROCAUC {max_train_auroc}, '
              f'microAUC: {max_train_micro_auc}, '
              f'Test: ROCAUC {max_test_auroc}, '
              f'microAUC: {max_test_micro_auc}'
              )

        train_aucs.append(max_train_auroc)
        train_micro_aucs.append(max_train_micro_auc)
        test_aucs.append(max_test_auroc)
        test_micro_aucs.append(max_test_micro_auc)

    for results in [train_aucs, train_micro_aucs, test_aucs, test_micro_aucs]:
        results = np.array(results)
        print(results.mean())


if __name__ == '__main__':
    # analyze_DTI_net_results()
    # test_Yamanishi_AUC()
    # write_predicted_DTIs()

    parse_BioSnap()




