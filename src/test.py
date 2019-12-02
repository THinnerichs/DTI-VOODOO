import numpy as np

import os
import sys


def test_target_subset():
    path = "../data/drug_target_relations/"

    counter = 0
    counter1 = 0

    for filename in os.listdir(path):
        if 'm' not in filename:
            continue
        s_file = filename.replace('m','s')
        if os.path.exists(path+s_file):
            m_file_targets = set()
            with open(file=path+filename, mode='r') as m:
                for line in m:
                    m_file_targets.add(line.split('\t')[0])
            s_file_targets = set()
            with open(file=path+s_file, mode='r') as s:
                for line in s:
                    s_file_targets.add(line.split('\t')[0])

            if not s_file_targets.issubset(m_file_targets):
                counter += 1
                print(counter, filename)
                print(len(s_file_targets), len(m_file_targets))
                print('length', len(s_file_targets - m_file_targets), 'of', len(m_file_targets))
                continue

            counter1 += 1
            if counter1 % 100 == 0:
                # print("counter1", counter1)
                pass

def run_PPI_parallel():
    import queue
    import threading
    import subprocess

    q = queue.Queue()

    workers = 6

    for s,t in [(a.min(), a.max()) for a in np.array_split(np.arange(538), workers)]:
        q.put((s,t))

    def worker():
        while True:
            doublet = q.get()
            if doublet is None:  # EOF?
                return
            s,t = doublet
            command = "python3 PPI_utils.py "+str(s)+" "+str(t+1)
            subprocess.call(command, shell=True)

    threads = [threading.Thread(target=worker) for _i in range(workers)]
    for thread in threads:
        thread.start()
        q.put(None)  # one EOF marker for each thread


def merge_protein_to_subgraph_dicts(start):
    import pickle

    dict_path = "../data/PPI_data/"

    filename = 'protein_to_subgraph_dict_' + str(start) + '.pkl'

    super_dict = {}

    super_filename = dict_path+'protein_to_subgraph_dict.pkl'
    if os.path.exists(super_filename):
        with open(file=super_filename, mode='rb') as f:
            super_dict = pickle.load(f)

    with open(file=dict_path + filename, mode='rb') as f:
        batch_dict = pickle.load(f)
    print("batch_dict size:", len(batch_dict))

    super_dict = super_dict.update(batch_dict)

    print("Writing merged dict to disk ...")
    filename = "../data/PPI_data/protein_to_subgraph_dict"
    with open(file=filename, mode='wb') as f:
        pickle.dump(super_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished.\n")

    '''
    super_dict = {}
    print("Merging dicts ...")
    for filename in os.listdir(dict_path):
        if not 'protein_to_subgraph_dict_' in filename:
            continue
        if not any(char.isdigit() for char in filename):
            continue

        batch_dict = {}
        with open(file=dict_path+filename, mode='rb') as f:
            batch_dict = pickle.load(f)
        print("batch_dict size:", len(batch_dict))

        super_dict = {**super_dict, **batch_dict}

    print("Finished.\n")

    print("Writing merged dict to disk ...")
    filename = "../data/PPI_data/protein_to_subgraph_dict"
    with open(file=filename, mode='wb') as f:
        pickle.dump(super_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished.\n")
    '''

if __name__ == '__main__':
    # test_target_subset()

    # run_PPI_parallel()

    _, start = sys.argv()
    merge_protein_to_subgraph_dicts(start)

    pass


