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
- try without DDI identitity
- use only boyce for DDIs
- use semsim-weighted PPI features 

- Transferred vs. normal
- take only experimental targets
  - Take proteins from nature paper (see Skype)
  - Take proteins and links from STITCH with appropriate score



