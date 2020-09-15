import numpy as np

import pickle as pkl

import DDI_utils


def write_PhenomeNET_files():
    path_prefix = "../data/PhenomeNET_data/"
    # parse drug HPO annotations and build updated drug list
    filename = "drug_SIDER_HPO_annotations.csv"
    print('Loading stereo to mono mapping...')
    drug_stereo_to_mono_mapping = DDI_utils.get_chemical_stereo_to_normal_mapping()

    drug_list = []
    drug_HPO_pairs = []
    print('Parsing SIDER associations...')
    with open(file=path_prefix + filename, mode='r') as f:
        # skip header
        f.readline()

        for line in f:
            split_line = line.split(',')
            drug, HPO_term = split_line[1].strip(), split_line[2].strip()
            drug = drug.replace('"', "")
            drug = drug.split("/")[-1]
            # map drugs from stereo to mono
            if drug.startswith('0'):
                drug = 'CIDm' + drug[1:]
            else:
                drug = drug_stereo_to_mono_mapping['CIDs' + drug[1:]]

            HPO_term = '<' + HPO_term[1:-1] + '>'

            print('drug, HPO_term:', drug, HPO_term)
            drug_list.append(drug)
            drug_HPO_pairs.append((drug, HPO_term))

            raise Exception

    # parse GO and MP annotations for proteins and update protein list
    onto_prefix = "<http://purl.obolibrary.org/obo/{entity}>"
    filename = "final_GO_ProteinID_human.txt"

    protein_list = []
    protein_GO_term_pairs = []
    print('Parsing protein-GO associations...')
    with open(file=path_prefix + filename, mode='r') as f:
        for line in f:
            GO_term, protein = line.strip().split('\t')
            GO_term = GO_term.replace(':', '_')
            GO_term = onto_prefix.format(entity=GO_term)

            protein_list.append(protein)
            protein_GO_term_pairs.append((protein, GO_term))

            raise Exception

    filename = "final_MP_ProteinID_human.txt"

    protein_MP_term_pairs = []
    print('Parsing protein-MP associations...')
    with open(file=path_prefix + filename, mode='r') as f:
        for line in f:
            MP_term, protein = line.strip().split('\t')
            MP_term = MP_term.replace(':', '_')
            MP_term = onto_prefix.format(entity=MP_term)

            protein_list.append(protein)
            protein_MP_term_pairs.append((protein, MP_term))

    # write drug and protein lists
    filename = 'PhenomeNET_drug_list'
    print('Writing drug and protein list...')
    with open(file=path_prefix+filename+'.pkl', mode='wb') as f:
        pkl.dump(drug_list, f, pkl.HIGHEST_PROTOCOL)

    filename = 'PhenomeNET_protein_list'
    with open(file=path_prefix + filename + '.pkl', mode='wb') as f:
        pkl.dump(protein_list, f, pkl.HIGHEST_PROTOCOL)

    print('drugs/proteins present:', len(drug_list), len(protein_list))

    # write drug-HPO-annotations and others to annotations file
    filename = "association_file"
    print('Writing association file...')
    with open(file=path_prefix+filename, mode='w') as f:
        for drug, HPO_term in drug_HPO_pairs:
            f.write(drug+' '+HPO_term+'\n')

        for protein, GO_term in protein_GO_term_pairs:
            f.write(protein+' '+GO_term+'\n')

        for protein, MP_term in protein_MP_term_pairs:
            f.write(protein+' '+MP_term+'\n')

    # write entity list
    filename = "entity_list"
    print('Writing entitiy list...')
    with open(file=path_prefix + filename, mode='w') as f:
        for drug in drug_list:
            f.write(drug+'\n')
        for protein in protein_list:
            f.write(protein+'\n')




if __name__ == '__main__':
    write_PhenomeNET_files()
