### ATC/InterPro results

This subdirectory contains nearly all results of our ATC/InterPro analysis presented in section 4.3.

We evaluated DTI-Voodoo over a 5-fold protein-split cross-validation within each InterPro protein family. The results are dumped into the `IPR..._results`
files. All used protein families are dumped into `ipro_classes`, while all associated ATC classes are available in `ATC_classes`

Each of these files contains the protein, the corresponding labels and DTI-Voodoos predictions. The drugs are ordered w.r.t. `./drug_list`, a dump of the corresponding numpy array.
Please use the python library `pickle` and `pkl.load(...)` to load `drug_list`. 

The actual novel predictions calculated, i.e. false positives with DTI-Voodoo confidence over the threshold of 0.7, are listed in in `new_false_positives_DTI_pairs.tsv`,
with the associated level 2 ATC group and fitting InterPro family for each drug and protein, respectively.
