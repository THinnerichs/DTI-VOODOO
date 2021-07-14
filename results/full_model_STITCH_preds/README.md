### full model STITCH preds

`best_preds_false_positives` contains the results of DTI-Voodoo trained and validated over the whole STITCH dataset, in contrast to the ATC/InterPro evaluation
where we evaluated for each InterPro protein family separately. 

To reduce the results to actually usable drugs, we filtered the drugs for existence of ATC annotations and the proteins for associations to cancer driver genes. 
For each drug we added the IUPAC notation for each drug and the Ensembl protein ID for each protein as aliases for better readibility. Further, we added the 
confidence provided by DTI-Voodoo for each interaction. Note that these confidence values are among the highest DTI-Voodoo will ever reaches. We cut these
values at the crisp threshold of 0.7. 
