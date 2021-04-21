#!/bin/bash

URL=https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait

URL=https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait

URL=https://stringdb-static.org/download/protein.sequences.v11.0/9606.protein.sequences.v11.0.fa.gz
curl -O "$URL" && gunzip -f "${URL##*/}"
wait

