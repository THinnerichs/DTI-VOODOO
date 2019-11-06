import numpy as np

from Bio import SeqIO


def main():


    # Create protein amino acid sequence from fasta
    print("Reading fasta file to dict...")
    fasta_filename = "../data/protein.sequences.v10.fa"
    protein_aa_seq_dict = SeqIO.index(fasta_filename, 'fasta')
    print("Finished.")

    # process drug-protein-interaction file
    print("Processing drug protein data...")
    protein_chemical_links_filename = "../data/protein_chemical.links.transfer.v5.0.tsv"
    paracetamol_id = 'CID1983'
    counter = 0
    paracetamol_prots_count = 0
    with open(file=protein_chemical_links_filename, mode='r') as f:
        for line in f:
            counter += 1
            paracetamol_prots_count += 1

            if counter%100000==0:
                print("Processed lines: {}\r".format(counter))

            if paracetamol_id not in line:
                continue

            split_line = line.strip().split('\t')
            if int(split_line[10]) < 700:
                continue

            drug = split_line[0]
            organism, protein = split_line[1].strip().split('.')

            seqio_seq = protein_aa_seq_dict[organism+'.'+protein]

            with open(file="../data/para_targets", mode='a') as para_handler:
                para_handler.write(organism+'.'+protein + "\t" + seqio_seq.seq+"\n")
    print("Finished.")




    print("Total line count:", counter)
    print("Paracetamol count:", paracetamol_prots_count)




if __name__=='__main__':
    main()