## TODO
- actually do some literature search!!!

- Put more comments to code

- run FAMSA on Ibex

- Run Alignment on min_score 700 to avoid that many 0 target drugs
- execute hmm_pipeline on MSA results


- parse similarity scores from Maxat script

- how to map e values from [0,inf) with 0 as hit and inf as no hit to [0,1] where 0 no hit and 1 is hit
- how to compare both similarity measures? (-> Maxat)

- calculate DDI of Boyce

- get STRING data

- get distibution over similarities and then take lowest, median and highest and fit above mapping to it

- iSCB conference deadline is 30. January

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

- What is the actual hypothesis? (Can drug-target interaction in other organisms can be used to determine DTIs in humans?) -> What is difference to GO, where all species are available as well
- My Fahrplan:
  - Proof similarity between semsim/jaccard similarity and found targets (most likely jaccard)
    - How to get correlation between two simlarities?
    - I only need alignments for drugs present in SIDER?
  - Remove orthologs -> query for humans (Orthologs are only the ones with same name?)
- Actually do some statistics on the removed orthologs
- What could be some meat to add here?
- What could be some backup plan here?

## Databases

- CTD

Tray Eidecke

