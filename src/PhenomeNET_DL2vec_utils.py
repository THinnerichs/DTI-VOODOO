import numpy as np

import pickle as pkl
import rdflib
from tqdm import tqdm

import DDI_utils


def write_PhenomeNET_files(mode='all'):
    path_prefix = "../data/PhenomeNET_data/"
    onto_prefix = "<http://purl.obolibrary.org/obo/{entity}>"

    # parse drug HPO annotations and build updated drug list
    # filename = "drug_SIDER_HPO_annotations.csv"
    filename = 'drugpheno_query_results.tsv'
    print('Loading stereo to mono mapping...')
    drug_stereo_to_mono_mapping = DDI_utils.get_chemical_stereo_to_normal_mapping()

    drug_list = []
    missed_drugs = []
    drug_HPO_pairs = []
    if mode=='drug' or mode=='all':
        print('Parsing SIDER associations...')
        missed_drugs_counter = 0
        with open(file=path_prefix + filename, mode='r') as f:
            # skip header
            f.readline()

            for line in f:
                drug, HPO_term = line.strip().split('\t')
                # map drugs from stereo to mono
                if drug.startswith('0'):
                    drug = 'CIDm' + drug[1:]
                else:
                    try:
                        drug = drug_stereo_to_mono_mapping['CIDs' + drug[1:]]
                    except:
                        missed_drugs_counter += 1
                        missed_drugs.append(drug)
                        continue

                HPO_term = onto_prefix.format(entity=HPO_term)

                drug_list.append(drug)
                drug_HPO_pairs.append((drug, HPO_term))
        print('Num drug-HPO-pairs:', len(drug_HPO_pairs))
        print('Num unique found drugs', len(set(drug_list)))
        print('Missed drugs from stereo/mono mapping:', missed_drugs_counter)

        return


    # parse GO and MP annotations for proteins and update protein list
    filename = "final_GO_ProteinID_human.txt"

    protein_GO_list = []
    protein_GO_term_pairs = []
    if mode=='GO' or mode=='all':
        print('Parsing protein-GO associations...')
        with open(file=path_prefix + filename, mode='r') as f:
            for line in f:
                GO_term, protein = line.strip().split('\t')
                GO_term = GO_term.replace(':', '_')
                GO_term = onto_prefix.format(entity=GO_term)

                protein_GO_list.append(protein)
                protein_GO_term_pairs.append((protein, GO_term))

        print('Num protein-GO-associations:', len(protein_GO_term_pairs))

    filename = "final_MP_ProteinID_human.txt"

    protein_MP_list = []
    protein_MP_term_pairs = []
    if mode=='MP' or mode=='all':
        print('Parsing protein-MP associations...')
        with open(file=path_prefix + filename, mode='r') as f:
            for line in f:
                MP_term, protein = line.strip().split('\t')
                MP_term = MP_term.replace(':', '_')
                MP_term = onto_prefix.format(entity=MP_term)

                protein_MP_list.append(protein)
                protein_MP_term_pairs.append((protein, MP_term))
        print('Num protein-MP-associations:', len(protein_MP_term_pairs))

    filename = "final_uberon_ProteinID_human.txt"

    protein_uberon_list = []
    protein_uberon_term_pairs = []
    if mode=='uberon' or mode=='all':
        print('Parsing protein-Uberon associations...')
        with open(file=path_prefix + filename, mode='r') as f:
            for line in f:
                uberon_term, protein = line.strip().split('\t')
                uberon_term = uberon_term.replace(':', '_')
                uberon_term = onto_prefix.format(entity=uberon_term)

                protein_uberon_list.append(protein)
                protein_uberon_term_pairs.append((protein, uberon_term))
        print('Num protein-uberon-associations:', len(protein_uberon_term_pairs))

    # build distinct drugs and proteins lists
    drug_list = sorted(list(set(drug_list)))
    protein_list = None
    if mode=='all':
        protein_list = sorted(list(set(protein_GO_list) & set(protein_MP_list) & set(protein_uberon_list)))
    else:
        protein_list = sorted(list(set(protein_GO_list) | set(protein_MP_list) | set(protein_uberon_list)))


    # write drug and protein lists
    filename = 'PhenomeNET_drug_list'
    if mode=='drug' or mode=='all':
        print('Writing drug and protein list...')
        with open(file=path_prefix+filename+'.pkl', mode='wb') as f:
            pkl.dump(drug_list, f, pkl.HIGHEST_PROTOCOL)

    filename = mode+'_protein_list' if mode in ['GO','MP','uberon','drug'] else 'PhenomeNET_protein_list'
    with open(file=path_prefix + filename + '.pkl', mode='wb') as f:
        pkl.dump(protein_list, f, pkl.HIGHEST_PROTOCOL)

    print('drugs/proteins present:', len(drug_list), len(protein_list))

    # write drug-HPO-annotations and others to annotations file
    filename = mode+'_association_file' if mode in ['GO','MP','uberon','drug'] else 'association_file'
    print('Writing association file...')
    with open(file=path_prefix+filename, mode='w') as f:
        if mode in ['drug','all']:
            for drug, HPO_term in drug_HPO_pairs:
                f.write(drug+' '+HPO_term+'\n')

        if mode in ['GO','all']:
            for protein, GO_term in protein_GO_term_pairs:
                f.write(protein+' '+GO_term+'\n')

        if mode in ['MP','all']:
            for protein, MP_term in protein_MP_term_pairs:
                f.write(protein+' '+MP_term+'\n')

        if mode in ['uberon','all']:
            for protein, uberon_term in protein_uberon_term_pairs:
                f.write(protein+' '+uberon_term+'\n')

    # write entity list
    filename = mode+'_entity_list' if mode in ['GO','MP','uberon','drug'] else 'entity_list'
    print('Writing entitiy list...')
    with open(file=path_prefix + filename, mode='w') as f:
        if mode =='drug' or mode=='all':
            for drug in drug_list:
                f.write(drug+'\n')
        if mode in ['GO','uberon', 'MP']:
            for protein in protein_list:
                f.write(protein+'\n')

def output_example_DL2vec_command(workers=48, embedsize=200):
    print('Running DL2vec embedding generator ...')
    path = '../../data/PhenomeNET_data/'
    onto = path + 'phenomenet.owl'
    asso = path + 'association_file'
    outfile = path + 'embedding_model'
    ents = path + 'entity_list'
    command = 'python runDL2vec.py -embedsize {embedsize} -ontology {onto} -associations {asso} -outfile {outfile} -entity_list {ents} -num_workers {num_workers}'.format(
        embedsize=embedsize,
        onto=onto,
        asso=asso,
        outfile=outfile,
        ents=ents,
        num_workers=workers)

    print('Command:', command)

def get_PhenomeNET_drug_list():
    filename = '../data/PhenomeNET_data/PhenomeNET_drug_list'
    with open(file=filename+'.pkl', mode='rb') as f:
        return pkl.load(f)

def get_PhenomeNET_protein_list(mode='all'):
    path = '../data/PhenomeNET_data/'
    filename = mode+'_protein_list' if mode in ['GO','MP','uberon'] else 'PhenomeNET_protein_list'
    with open(file=path+filename+'.pkl', mode='rb') as f:
        return pkl.load(f)

def write_query_drugpheno_rdf_graph():
    drugpheno_rdf_graph_filename = "../data/PhenomeNET_data/data-2020-12-07/drugphenotype.rdf"
    drugpheno_graph = rdflib.ConjunctiveGraph()
    print('Parsing drugpheno RDF graph...')
    result = drugpheno_graph.parse(drugpheno_rdf_graph_filename, format='xml')

    print('Query drugpheno RDF graph...')
    qres = drugpheno_graph.query(
        """SELECT ?concept ?phenotype
WHERE {
    ?association rdf:type rdf:Statement .
    ?association rdf:predicate obo:RO_0002200 .
    ?association rdf:object ?phenotype .
    ?association rdf:subject ?concept .
    ?concept rdf:type <http://phenomebrowser.net/Drug> .
    ?association dc:provenance ?prov .
    ?prov dc:creator ?creator .
    ?prov dcterms:source ?source .
}
""")

    outfile = '../data/PhenomeNET_data/drugpheno_query_results.tsv'

    print('Writing results...')

    results = []
    counter = 0
    for d,p in qres:
        results.append((d,p))
        counter += 1
        if counter % 1000 == 0:
            print(counter)
    print(len(results))
    with open(file=outfile, mode='w') as f:
        for drug, phenotype in tqdm(results):
            drug = drug.strip().split('/')[-1]
            phenotype = phenotype.strip().split('/')[-1]
            print(drug+'\t'+phenotype, file=f)
    print('Done.', len(results), 'drugpheno pairs.')


if __name__ == '__main__':
    write_PhenomeNET_files('drug')
    # write_query_drugpheno_rdf_graph()
