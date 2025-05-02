#!/bin/bash

while IFS="" read -r p || [ -n "$p" ]
do
    echo "${p}" |
    sed -E 's:_|\..+::g' |
    sed -E 's:.{3}:&/:g' |
    awk -v base="ftp://ftp.ncbi.nih.gov/genomes/all/" '{print base""$0}' |
    xargs curl -si |
    sed -n -e $(echo "${p}" |
    sed -E 's:.*\.(.*):\1p:g') |
    sed -E 's: +:\t:g' |
    cut -f9 >> full_acc.txt
done < records.txt
