import numpy as np

import similarity_measurement


def get_SIDER_to_HPO_mapping_dict():
    filename = '../data/SIDER_data/MedGen_HPO_Mapping.txt'

    mapping_dict = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()

        for line in f:
            split_line = line.strip().split('|')
            se_id, HPO_id = split_line[0], split_line[1]

            mapping_dict[se_id] = HPO_id

    return mapping_dict


