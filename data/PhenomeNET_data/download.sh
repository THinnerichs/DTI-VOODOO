#!/bin/bash

URL=http://aber-owl.net/media/ontologies/PhenomeNET/1/phenomenet.owl
curl -O "$URL" && gunzip -f "${URL##*/}"
wait
