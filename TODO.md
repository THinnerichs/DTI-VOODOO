## TODO

1.1. Find DTI paper which uses drughub

Also write down failed approach?

Options:
- GCN vs GAT vs SAGE vs others
- Protein vs drug split

- molecular drug representation:
  - SMILES transformer?
  - proper GCN over molecular structure? (What node features? See other papers for that)
  - take model from them and enhance it with PPI graph?
- Model for drug ? For each protein?


Draw images for Networks
- PPI & DDI & Side effects data
- Molecular predictor
- Graph networks
- Images for splits

Robert talk:
  - search bioinformatics for drug repuposing
  - 10 years ago side effect similarity paper in Science  (0.7 AUC)
  - locality in PPI network coincides with phenotypes
    - See paper in Skype

  - function phenotypes

- What is usage of/hypothesis for 
  - GO / protein functions
  - anatomical site of expression
  - loss of function phenotypes
    - [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933875/)

- remove self loops? 
- Resnik Similarity drug-protein-similarity

- Remove self-loops from GCN
- use only boyce for DDIs


Data for experimental only:
- 300: num_nodes 1828, num_edges 4494
- 500: num_nodes 1725, num_edges 3620
- 700: num_nodes 1578, num_edges 3043

Robert talk 28.05.:
- side effects and phenotype similarity
  - see mouse model features 
  - use Siamese network

- normalize node features to Gaussian?
  - divide by mean and substract 1
  - plot histogram
- use whole PPI graph with less important edges and their scores too
  - quadratic increase in computational time

- semsim out of order? values at identity are not 1 
  - Is it algorithm or my dumbness?
  - Get maximum of each row and its index

Data for PPI graph reduction:
- no reduction: nodes: 16042, edges: 401752
- reduction: 
  - 100: nodes: 17982, edges: 1711038
  - 300: nodes: 15647, edges: 220167
  - 500: nodes: 12120, edges: 80325
  - 700: nodes: 8290, edges: 48558
  - 900: nodes: 5055, edges: 23601
