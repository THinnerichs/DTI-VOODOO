import numpy as np

from Bio import SeqIO
import subprocess
import time

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


def prune_drug_protein_db(min_score=700):

    filename = "../data/protein_chemical.links.transfer.v5.0.tsv"
    target_filename = "../data/protein_chemical.links.min_score_" + str(min_score) + ".tsv"

    print("Processing huge file ...")
    with open(file=filename, mode='r') as f, open(file=target_filename, mode='w') as targetfile:
        targetfile.write(f.readline())

        counter = 0

        for line in f:
            counter += 1
            if counter % 10000000 == 0:
                print("Processed lines:", counter)

            if int(line.strip().split('\t')[10]) < min_score:
                continue
            targetfile.write(line)
    print("Finished.")


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

def write_paracetamol_prots_to_fasta(min_score=400):

    para_filename = "../data/para_targets"
    para_fasta_filename = "../data/para_fasta_" + str(min_score) + "_min_score.fasta"
    print("Processing {} and writing {} ...".format(para_filename, para_fasta_filename))
    with open(file=para_filename, mode='r') as para_file, open(file=para_fasta_filename, mode='w') as fasta_file:
        for line in para_file:
            protein, aa_seq, score = line.split('\t')
            if int(score) < min_score:
                continue

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
            org_prot, _, score = line.split('\t')
            if int(score) < 700:
                continue
            split_prot = org_prot.split('.')
            organism = split_prot.pop(0)
            protein = ".".join(split_prot)

            protein_set.add(protein)

    print("set size:", len(protein_set))

    with open(file="../data/setified_para_proteins_700", mode="w") as f:
        for ele in protein_set:
            f.write(ele+'\n')

def run_para_multi_sequence_alignment(min_score=700,
                                      alignment_method='clustalo'):
    from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline, MafftCommandline, MSAProbsCommandline, TCoffeeCommandline

    in_file = "../data/para_fasta_" + str(min_score) + "_min_score.fasta"
    out_file = "../data/para_" + alignment_method + "aligned" + str(min_score) + "_min_score.fasta"
    start_time = time.time()

    command = None
    if alignment_method == 'clustalo':
        command = ClustalOmegaCommandline(infile=in_file,
                                          outfile=out_file,
                                          verbose=True,
                                          auto=True)
        command = "./" + str(command)
    elif alignment_method == 'muscle':
        command = MuscleCommandline(infile=in_file,
                                    outfile=out_file,
                                    verbose=True,
                                    auto=True)
        command = "./" + str(command)
    elif alignment_method == 'mafft':
        command = MafftCommandline(input=in_file)
        command = "./mafft-linux64/mafft.bat --anysymbol" + ' '.join(str(command).split(' ')[1:])

        print(command)

        print("Starting {} alignment ...".format(alignment_method))
        stdout= subprocess.check_output(str(command), shell=True)

        print("Finished in {} sec.".format(time.time()-start_time))

        with open(file=out_file, mode='w') as f:
            f.write(stdout)

    elif alignment_method == 'msaprobs':
        command = MSAProbsCommandline(infile=in_file,
                                      outfile=out_file,
                                      verbose=True,
                                      auto=True)
        command = "./MSAProbs-0.9.7/MSAProbs/msaprobs " + ' '.join(str(command).split(' ')[1:])
    '''
    elif alignment_method == 'tcoffee':
        command = TCoffeeCommandline(infile=in_file,
                                     outfile=out_file,
                                     verbose=True,
                                     auto=True)
    '''

    print(command)


    print("Starting {} alignment ...".format(alignment_method))
    subprocess.call(str(command), shell=True)
    print("Finished in {} sec.".format(time.time()-start_time))

def para_PWM_from_alignment(min_score=700,
                            alignment_method='clustalo'):
    from Bio.Alphabet import IUPAC, Gapped
    from Bio import AlignIO, Alphabet, motifs

    alphabet = Gapped(IUPAC.protein)

    filename = "../data/para_" + alignment_method + "aligned" + str(min_score) + "_min_score.fasta"

    alignment = AlignIO.read(filename, 'fasta', alphabet=alphabet)
    m = motifs.create([x.seq for x in alignment])

    print(m.consensus)
    print(m.counts)













if __name__=='__main__':
    # write_pruned_SeqIO_fasta_dict()
    # write_paracetamol_prots_to_file()

    # eliminate_para_target_duplicates()
    # run_stitch_db_query()
    # prune_drug_protein_db(min_score=400)

    '''
    min_score = 700

    write_paracetamol_prots_to_fasta(min_score=min_score)
    '''

    run_para_multi_sequence_alignment(min_score=700,
                                      alignment_method='mafft')

    # para_PWM_from_alignment()
