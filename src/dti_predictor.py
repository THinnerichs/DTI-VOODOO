import numpy as np


def main():
    protein_chemical_links_filename = "../data/protein_chemical.links.transfer.v5.0.tsv"

    paracetamol_id = "CID1983"
    counter = 0
    with open(file=protein_chemical_links_filename, mode='r') as f:
        for line in f:
            counter += 1
            if counter%100000==0:
                print("Processed lines: {}\r".format(counter))

            '''
            if paracetamol_id not in line:
                continue
            split_line = line.strip().split('\t')
            drug = split_line[0]
            organism, protein = split_line[1].strip().split('.')
            '''


    print("Total line count:", counter)




if __name__=='__main__':
    main()