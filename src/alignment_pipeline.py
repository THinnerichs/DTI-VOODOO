import numpy as np

from Bio import SeqIO
import subprocess
import time

import pickle


def separate_prots_to_files(file_min_score=400,
                            min_score=400):

    # Create protein amino acid sequence from fasta
    print("Reading pruned dict...")
    dict_filename = "../data/prot_aa_seq_dict"
    protein_aa_seq_dict = None
    with open(dict_filename + '.pkl', 'rb') as f:
        protein_aa_seq_dict = pickle.load(f)
    print("Finished.")

    # process drug-protein-interaction file
    print("Processing drug protein data...")
    protein_chemical_links_filename = "../data/protein_chemical.links.min_score_"+str(file_min_score)+".tsv"


    counter = 0
    with open(file=protein_chemical_links_filename, mode='r') as f:
        f.readline()
        for line in f:
            counter += 1

            if counter%5000000==0:
                # print("Progress: %.2f %%\r"%(counter*100/total_lines))
                print("Processed lines:", counter)

            split_line = line.strip().split('\t')
            if int(split_line[10]) < min_score:
                continue

            drug = split_line[0]
            split_protein = split_line[1].strip().split('.')
            organism = split_protein.pop(0)
            protein = ".".join(split_protein)
            score = int(split_line[10])

            aa_seq = protein_aa_seq_dict[organism][protein]

            with open(file="../data/drug_target_relations/"+drug+"_targets", mode='a') as drug_handler:
                drug_handler.write(organism+'.'+protein + '\t' + aa_seq + '\t'+ str(score) + '\n')
    print("Finished.")

def create_fasta(drug_name,
                 min_score=700):
    
    para_filename = "../data/"+drug_name+"_targets"
    para_fasta_filename = "../data/"+mol_name+"_fasta_" + str(min_score) + "_min_score.fasta"
    print("Processing {} and writing {} ...".format(para_filename, para_fasta_filename))
    with open(file=para_filename, mode='r') as para_file, open(file=para_fasta_filename, mode='w') as fasta_file:
        for line in para_file:
            protein, aa_seq, score = line.split('\t')
            if int(score) < min_score:
                continue

            fasta_file.write(">"+protein+'\n')
            fasta_file.write(aa_seq+'\n')

    print("Finished.")

if __name__ == '__main__':
    separate_prots_to_files(file_min_score=400,
                            min_score=400)
