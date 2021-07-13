import pickle
import numpy as np
import random

import matplotlib.pyplot as plt

from tqdm import tqdm

import DTI_data_preparation
import DDI_utils
import dti_utils

import sklearn.metrics as metrics
from sklearn.model_selection import KFold


def write_predicted_DTIs(fold=3):
    '''
    filename = '../models/graph_models/PPI_network_model_with_mol_features_fold_'+str(fold)+'_predictions.pkl'

    pred_list = None #[(drug, prot, label, pred)]
    with open(file=filename, mode='rb') as f:
        pred_list = pickle.load(f)
    '''

    pred_list = []
    for i in range(1,6):
        print(f'Parsing fold {str(i)}')
        path = '../results/full_model_STITCH_preds/all_'
        filename = path + 'results_fold_' + str(i)

        with open(file=filename, mode='r') as f:
            for line in f:
                drug, prot, val = line.strip().split('\t')
                pred_list.append((drug,prot,float(val)))

    print('num triples:', len(pred_list))

    drug_mapping = DTI_data_preparation.get_drug_to_name_mapping()
    protein_mapping = DTI_data_preparation.get_protein_to_EnsemblProtein_id()

    print('drug_dict.keys()', len(list(drug_mapping.keys())))
    print('protein_dict.keys()', len(list(protein_mapping.keys())))

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


    pred_list = [tup for tup in pred_list if tup[2] > 0.7]

    pred_list = sorted(pred_list, key=lambda tup: tup[2], reverse=True)
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

    drug_ATC_mapping = DTI_data_preparation.get_drug_ATC_classification_mapping()

    print('ATC mapping length', len(drug_ATC_mapping))
    print(list(drug_ATC_mapping.items())[:10])
    drug_list = [drug for drug in drug_list if drug in drug_ATC_mapping.keys()]

    print('drugs that are also present as ATC classes:', len(drug_list))

    filename = '../results/full_model_STITCH_preds/best_preds_false_negatives'
    with open(file=filename, mode='w') as f:
        print('drug\tprot\tconfidence\tdrug_alias\tprotein_alias\tcancer_type', file=f)
        for drug, protein, confidence in pred_list:
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

def analyze_ATC_Ipro_rel():
    # ipro_class_list = ['IPR013087', 'IPR000276', 'IPR000742', 'IPR003961', 'IPR007110', 'IPR000719', 'IPR002048', 'IPR001452', 'IPR000504', 'IPR001356']
    ipro_class_list = ['IPR013087', 'IPR017452', 'IPR000276', 'IPR001909', 'IPR001841', 'IPR003961', 'IPR001806',
                  'IPR001314', 'IPR001965', 'IPR001304', 'IPR002035', 'IPR000315', 'IPR001781', 'IPR000571']

    ATC_mapping = DTI_data_preparation.get_drug_ATC_classification_mapping()
    drug_list = None
    with open(file='../results/ATC_Interpro_results/drug_list', mode='rb') as f:
        drug_list = pickle.load(f)
    redundant_present_ATC_classes = np.array([ATC_mapping.get(drug, None) for drug in drug_list])
    present_ATC_classes = list(set(redundant_present_ATC_classes))
    present_ATC_classes.remove(None)

    res = [0] * len(present_ATC_classes)
    for drug in drug_list:
        query_res = ATC_mapping.get(drug, None)

        if query_res:
            res[present_ATC_classes.index(query_res)] += 1

    print(res[:10])
    sorted_classes = np.argsort(present_ATC_classes)

    '''
    plt.figure(figsize=(12, 8))
    plt.xticks(ticks=np.arange(len(present_ATC_classes)), labels=np.array(present_ATC_classes)[sorted_classes], rotation=90)
    plt.bar(np.array(present_ATC_classes)[sorted_classes], np.array(res)[sorted_classes])
    plt.tight_layout()
    plt.savefig('../results/ATC_Interpro_results/ATC_bar_distr.png', dpi=80)
    '''

    print('Number of present level 2 ATC classes:', len(present_ATC_classes))

    results_list = np.zeros((len(ipro_class_list), len(present_ATC_classes)))
    redundant_present_ATC_classes = np.array(redundant_present_ATC_classes)

    new_drug_target_list = []
    for ipro_index, ipro_class in enumerate(ipro_class_list):
        filename = f'../results/ATC_Interpro_results/{ipro_class}_'

        with open(file=filename+'results', mode='r') as f:
            # count lines of results file
            num_prots = sum(1 for _ in f)

        with open(file=filename + 'results', mode='r') as f:
            for line in f:
                prot, labels, predictions = line.strip().split('\t')
                labels = np.array([int(float(i)) for i in labels[1:-1].split(',')])
                # predictions = np.array([round(float(i)) for i in predictions[1:-1].split(',')])
                predictions = np.array([float(i) for i in predictions[1:-1].split(',')])

                mask = ((predictions-labels)>0.7).astype(int)
                for drug in np.array(drug_list)[mask==1]:
                    ret_val = ATC_mapping.get(drug, None)
                    if ret_val:
                        new_drug_target_list.append((drug, ATC_mapping[drug], prot, ipro_class))

                '''
                for ATC_index, ATC_class in enumerate(present_ATC_classes):
                    mask = redundant_present_ATC_classes==ATC_class
                    
                    results_list[ipro_index,ATC_index] = ((predictions[mask] - labels[mask])==1).astype(int).sum()
                '''
        # normalize by number of prots
        # results_list[ipro_index, :] /= num_prots

    print(len(new_drug_target_list))
    filename = ''
    with open(file='../results/ATC_Interpro_results/new_false_negatives_DTI_pairs.tsv', mode='w') as f:
        for d, c, t, family in new_drug_target_list:
            print(f'{d}\t{c}\t{t}\t{family}', file=f)

    raise Exception









    # Normalize with number of drugs in each ATC class
    for i, ATC_class in enumerate(present_ATC_classes):
        results_list[:, i] /= res[i]

    sorted_indices = np.argsort(results_list.flatten())
    highest_hits_list = []
    for i in sorted_indices[-100:]:
        ipro_index = i//71
        ATC_index = i%71

        highest_hits_list.append((ipro_class_list[ipro_index], present_ATC_classes[ATC_index], results_list.flatten()[i]))

    path = '../results/ATC_Interpro_results/'
    # np.savetxt(path+'ATC_IPRO_class_relations', results_list, fmt='%.4f', delimiter=' ', comments='# matrix (ipro x ATC classes) consisting of number of newly predicted interactions by DTI-Voodoo, normalized by the amount of proteins in that InterPro family; Class names present in ipro_classes and ATC_classes')
    # np.savetxt(path+'ipro_classes', np.array(ipro_class_list), fmt='%s', delimiter=',')
    # np.savetxt(path+'ATC_classes', np.array(present_ATC_classes), fmt='%s', delimiter=',')

    atc_sorted_indices = np.argsort(np.array(present_ATC_classes))
    present_ATC_classes = np.array(present_ATC_classes)
    plt.figure(figsize=(18, 8))
    plt.xticks(ticks=np.arange(len(present_ATC_classes)), labels=present_ATC_classes[atc_sorted_indices], rotation=90)
    plt.yticks(ticks=np.arange(len(ipro_class_list)), labels=ipro_class_list)
    hm = plt.imshow(results_list[:, atc_sorted_indices], interpolation="nearest")

    plt.colorbar(hm, shrink=0.5, aspect=20*0.5)

    plt.tight_layout()

    plt.savefig('../results/ATC_Interpro_results/ATC_Interpro_heatmap.png', dpi=80)

    # return results_list, ipro_class_list, present_ATC_classes
    return highest_hits_list


if __name__ == '__main__':
    # analyze_DTI_net_results()
    # test_Yamanishi_AUC()
    # write_predicted_DTIs()

    # parse_BioSnap()

    analyze_ATC_Ipro_rel()
