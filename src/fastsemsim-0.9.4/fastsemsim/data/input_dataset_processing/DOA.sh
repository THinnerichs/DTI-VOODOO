
# Parse updated DOA from Peng et al, 2013, NAR
# RIF: Peng, K., Xu, W., Zheng, J., Huang, K., Wang, H., Tong, J., … Xia, T. (2013). The Disease and Gene Annotations (DGA): an annotation resource for human disease. Nucleic Acids Research, 41(Database issue), D553–60. doi:10.1093/nar/gks1244
# Download the rdf file from http://dga.nubic.northwestern.edu/pages/download.php and name it DO_Peng_NAR_2013_annotations.rdf
# Direct wget does not seem to work (ajax trick)
# Then...
OUTFILE=DO_human_Peng_NAR_`date "+%Y.%m.%d"`.txt
cat DO_Peng_NAR_2013_annotations.rdf | grep  "<cd:DOID" -A 2 > $OUTFILE
tr "\r" " " < $OUTFILE > temp.txt
tr "\n" " " < temp.txt > temp2.txt
tr "-" "\n" < temp2.txt > temp.txt
mv temp.txt $OUTFILE
cat $OUTFILE | grep  "DOID" > temp.txt
mv temp.txt $OUTFILE
cat $OUTFILE  | sed -e 's/  <cd:DOID>/DOID:/g' -e 's|</cd:DOID>||g' -e 's/<cd:GeneID>//g' -e 's|</cd:GeneID>||g' -e 's/<cd:PubMedID>//g' -e 's|</cd:PubMedID>||g' > temp.txt
mv temp.txt $OUTFILE
sed -e 's/ DOID/DOID/g' $OUTFILE > temp.txt
mv temp.txt $OUTFILE
cut -d " " -f 1,5,9  $OUTFILE > temp.txt
mv temp.txt $OUTFILE
cat $OUTFILE | gawk -F " " '{print $2"\t"$1"\t"$3}'  > temp.txt
mv temp.txt $OUTFILE
mv $OUTFILE ../ACs

# Parse DOA from Osborne work: Osborne, J. D., Flatow, J., Holko, M., Lin, S. M., Kibbe, W. a, Zhu, L. J., … Chisholm, R. L. (2009). Annotating the human genome with Disease Ontology. BMC Genomics, 10 Suppl 1, S6. doi:10.1186/1471-2164-10-S1-S6
# wget http://projects.bioinformatics.northwestern.edu/do_rif/generifs_basic
wget http://projects.bioinformatics.northwestern.edu/do_rif/do_rif.human.txt
OUTFILE=DO_human_Osborne_BMC_2010.08.24.txt # keep the data fixed, since it is the last update they have made
cp do_rif.human.txt $OUTFILE
cat $OUTFILE | gawk -F "\t" '{print $1"\t"$5"\t"$3}' > temp.txt
mv temp.txt $OUTFILE
mv $OUTFILE ../ACs


# clean up
rm do_rif.human.txt
rm temp*
rm DO_Peng_NAR_2013_annotations.rdf
