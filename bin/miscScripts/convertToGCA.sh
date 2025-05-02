#!/bin/bash

while IFS="" read -r p || [ -n "$p" ]
do
    esearch -db assembly -query $p < /dev/null | \
    esummary | \
    xtract -pattern DocumentSummary -element AssemblyAccession >> records.txt
done < accession_list.txt