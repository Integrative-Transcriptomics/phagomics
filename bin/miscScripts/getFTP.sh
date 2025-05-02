#!/bin/bash

while IFS="" read -r p || [ -n "$p" ]
do
    echo "${p}" |
    sed -E 's:_|\..+::g' |
    sed -E 's:.{3}:&/:g' |
    awk -v base="ftp://ftp.ncbi.nih.gov/genomes/all/" '{print base""$0}' >> ftp_links.txt
done < records.txt
