import numpy as np

from Bio import SeqIO
import subprocess
import time

import pickle


def separate_prots_to_files(file_min_score=400,
                            min_score=400):
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

            with open(file="../data/drug_target_relations/"+drug+"_targets", mode='a') as drug_handler:
                drug_handler.write(organism+'.'+protein + "\t" + str(score) + "\n")
    print("Finished.")


if __name__ == '__main__':
    separate_prots_to_files(file_min_score=400,
                            min_score=400)
