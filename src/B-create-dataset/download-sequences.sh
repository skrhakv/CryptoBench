#!/bin/bash

OUTPUT_PATH="$1"

for acc in `cat $OUTPUT_PATH/uniprot_ids.txt` ; do curl -s "https://rest.uniprot.org/uniprotkb/$acc.fasta" ; done > $OUTPUT_PATH/uniprot_seqs.fasta 2>$OUTPUT_PATH/uniprot_download_error.txt
