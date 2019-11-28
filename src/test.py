import numpy as np

import os


def test_target_subset():
    path = "../data/drug_target_relations/"

    counter = 0

    for filename in os.listdir(path):
        if 'm' not in filename:
            continue
        print("bumm")
        s_file = filename.replace('m','s')
        if os.path.exists(s_file):
            print(s_file)
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



if __name__ == '__main__':
    test_target_subset()


