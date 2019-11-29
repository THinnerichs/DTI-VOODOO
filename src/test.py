import numpy as np

import os


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



if __name__ == '__main__':
    test_target_subset()

