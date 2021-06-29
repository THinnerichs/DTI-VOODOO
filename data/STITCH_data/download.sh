#!/bin/bash

URL=http://stitch.embl.de/download/protein_chemical.links.transfer.v5.0/9606.protein_chemical.links.transfer.v5.0.tsv.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait

URL=http://stitch.embl.de/download/actions.v5.0/9606.actions.v5.0.tsv.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait


URL=http://stitch.embl.de/download/chemical.aliases.v5.0.tsv.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait

URL=http://stitch.embl.de/download/chemicals.v5.0.tsv.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait
