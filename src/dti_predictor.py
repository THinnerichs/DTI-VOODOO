import numpy as np

from Bio import SeqIO

import pickle

def write_pruned_SeqIO_fasta_dict():
    # Create protein amino acid sequence from fasta
    print("Reading fasta file to dict...")
    fasta_filename = "../data/protein.sequences.v10.fa"
    protein_aa_seq_dict = SeqIO.index(fasta_filename, 'fasta')
    print("Finished.")

    pruned_dict = {}
    for key, seqio_seq in protein_aa_seq_dict.items():
        split_list = key.split('.')

        organism = split_list.pop(0)
        protein = '.'.join(split_list)
        result = pruned_dict.get(organism, None)

        if not result:
            pruned_dict[organism] = {}

        pruned_dict[organism][protein] = str(seqio_seq.seq)
    print("Finished.")

    print("Writing pruned dict to disc")
    filename = "../data/prot_aa_seq_dict"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(pruned_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename)

def write_paracetamol_prots_to_file():

    # Create protein amino acid sequence from fasta
    print("Reading pruned dict...")
    dict_filename = "../data/prot_aa_seq_dict"
    protein_aa_seq_dict = None
    with open(dict_filename + '.pkl', 'rb') as f:
        protein_aa_seq_dict = pickle.load(f)
    print("Finished.")

    # process drug-protein-interaction file
    print("Processing drug protein data...")
    protein_chemical_links_filename = "../data/protein_chemical.links.transfer.v5.0.tsv"
    paracetamol_id = 'CIDs00001983'
    counter = 0
    total_lines = 7019540873
    paracetamol_prots_count = 0
    with open(file=protein_chemical_links_filename, mode='r') as f:
        for line in f:
            counter += 1
            paracetamol_prots_count += 1

            if counter%5000000==0:
                print("Progress: %.2f %%\r"%(counter*100/total_lines))

            if paracetamol_id not in line:
                continue

            split_line = line.strip().split('\t')
            if int(split_line[10]) < 400:
                continue

            drug = split_line[0]
            split_protein = split_line[1].strip().split('.')
            organism = split_protein.pop(0)
            protein = ".".join(split_protein)
            score = int(split_line[10])

            aa_seq = protein_aa_seq_dict[organism][protein]

            with open(file="../data/para_targets", mode='a') as para_handler:
                para_handler.write(organism+'.'+protein + "\t" + aa_seq+'\t' + str(score) + "\n")
    print("Finished.")

    print("Total line count:", counter)
    print("Paracetamol count:", paracetamol_prots_count)

def write_paracetamol_prots_to_fasta():

    para_filename = "../data/para_targets"
    para_fasta_filename = "../data/para_fasta.fasta"
    print("Processing {} and writing {} ...".format(para_filename, para_fasta_filename))
    with open(file=para_filename, mode='r') as para_file, open(file=para_fasta_filename, mode='w') as fasta_file:
        for line in para_file:
            protein, aa_seq = line.split('\t')

            fasta_file.write(">"+protein+'\n')
            fasta_file.write(aa_seq+'\n')

    print("Finished.")

def run_stitch_db_query():
    import requests
    import sys

    # build url
    paracetamol_id = "CID000001983"
    # paracetamol_id = "CID3676"
    required_score = 400
    url = "http://stitch.embl.de/api/tsv/interactors?identifier=" + paracetamol_id + \
          "&required_score=" + str(required_score) + "%0A"\
          "species=9606"
    url2 = "http://stitch.embl.de/api/tsv/resolve?identifier=" + paracetamol_id
    print(url2)

    # Run query with exception handling
    r = requests.get(url)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    seq = r.text

    print(seq)

def eliminate_para_target_duplicates():
    filename="../data/para_targets"
    protein_set = set()
    with open(file=filename, mode='r') as f:
        for line in f:
            org_prot, _, _ = line.split('\t')
            split_prot = org_prot.split('.')
            organism = split_prot.pop(0)
            protein = ".".join(split_prot)

            protein_set.add(protein)

    print("set size:", len(protein_set))

    with open(file="../data/setified_para_proteins", mode="w") as f:
        for ele in protein_set:
            f.write(ele+'\n')

def run_multi_sequence_alignment():
    pass







if __name__=='__main__':
    # write_pruned_SeqIO_fasta_dict()
    # write_paracetamol_prots_to_file()

    eliminate_para_target_duplicates()
    # run_stitch_db_query()