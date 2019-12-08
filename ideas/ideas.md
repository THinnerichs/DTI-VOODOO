## Important

- Extend VPN access (only valid until 30.12.19)
- Visit ithelpdesk on Sunday for VPN credentials

## TODO

- actually do some literature search!!!

- Put more comments to code

- how to map e values from [0,inf) with 0 as hit and inf as no hit to [0,1] where 0 no hit and 1 is hit

- Build alignment over all targets
  - build hmms over alignments
  - Query against only human target data
  - Build STRING graph only on human data (Y)
  - Use human interactions data to build training/testing data
    - Utilize non-binary score? Sigmoid over score?

- get phenotypes for drugs
  - download OMIM -> get mapping
  - Human phenotype ontology

- Finally execute hmm\_pipeline !!! 

- check if VPN is working

---------------------------------------------------
- read Roberts Hypothesis testing paper for title and outline

- What do links in PPI graph mean? (-> STRING data/documentation)
  - what does closeness in enriched PPI graph means?
- Why do i need the DDIs and the Side effect similarity

- Read PREDICT paper again

- are there papers on DTI using PPIs?

- test unsupervised node embdding learning
- Fix Hmm pipeline
- Hannes BA lesen
- Get other method to run on my dataset (literature!!)
  - HAve to get SMILES and stuff too -> meh
  - How can this be made proof? Model is tweaked on this dataset while foreign methods are not

- prepare presentation
- read PREDICT paper
- maybe try other stellargraph embeddings as well :)))
- run link prediction over whole graph


- run unsupervised methods now as more RAM is available


----------------------------------------------------
Talk with Robert on 03.12.:
- Test Graphsage
- Get Phenotypes from OMIM
- Hypothesis is about



## Ideas for presentation

- include numbers of interactions
  - histograms?
- Other approaches lack ability to generalize -> guilt by association

- histogram over scores of STRING and STITCH
- show image of network for STRING and STITCH
- Meaning of score 
- Functions to slides B)


## Other ideas

- Calculate similarity between drugs and metabolites...
	... over structural similarity
	... over alignment of e.g. SMILES

- get DD similarity trough ATC code/hierarchy?

## On performance of alignments

This [link](https://www.ebi.ac.uk/Tools/msa/)

## Other sources for MSAs

- [wiki article on MSAs](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software#Multiple_sequence_alignment)

## Questions



- What is the actual hypothesis? 
- What could be some meat to add here?
## Databases

- CTD

Tray Eidecke

