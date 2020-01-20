## Important

- Extend VPN access (only valid until 31.12.19)

## TODO

- actually do some literature search!!!

- Put more comments to code

- how to map e values from [0,inf) with 0 as hit and inf as no hit to [0,1] where 0 no hit and 1 is hit

---------------------------------------------------
- read Roberts Hypothesis testing paper for title and outline

- are there papers on DTI using PPIs?

- test unsupervised node embedding learning
- Get other method to run on my dataset (literature!!)
  - Have to get SMILES and stuff too -> meh
  - How can this be made proof? Model is tweaked on this dataset while foreign methods are not

- read PREDICT paper
- maybe try other stellargraph embeddings as well :)))

- run unsupervised methods now as more RAM is available

----------------------------------------------------

- extract predictions from search results

- Read paper how convolution over node features is done
  - Otherwise use same structure but with different node labels?

- update targets based on new alignments
  - update after new hmm\_search went trough :)

---------------------------------------------------
- DTI\_graph:
  - min\_score:
    - > 700: 18987
    - > 400: 33549
    - all: 171168

---------------------------------------------------

- Find SOTA for DTI
  - take intersection of data and perform on that
  - Maybe one top-down and one bottom up?

- Get Stellargraph models to run with my data -> Compare results
  - Rebuild their models (See model images)

- Finally attach filter data to graph!!

- Test different filter methods and their performance
  - Test different parameters for hmm\_build and search 

--------------------------------------------------

Talk with Robert on 09.01.:

Work until now:
  - performed mafft and FAMSA alignment with search -> poor results
    1. tweak this! because its new!
    2. take other bottom up approach to fix this
  - Tweaking GCN, GAT model -> not really sufficient, performs horrible 
    - heavy problems getting that to run -> decompose stellargraph methods
    - other library?! -> you had time over christmas
      - Top k pooling? pooling layers in general?
      - Pytorch geometric

  - artificially reduce PPI > 700 score

  - Take Proteins without links into model at the end, not for tuning
  
  - used DTIs accurate? currently just taking all with score aboth 700 (just captures cooccurence)
  - used PPIs accurate? Currently taking all of them

  - Deadline? Icsb is in June


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

