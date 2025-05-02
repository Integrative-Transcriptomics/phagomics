#!/bin/bash

paste -d@ ftp_links.txt full_acc.txt | while IFS="@" read -r f1 f2
do
    wget $f1/$f2/${f2}_genomic.gff.gz -P phage_data/$f2
    wget $f1/$f2/${f2}_genomic.fna.gz -P phage_data/$f2
    wget $f1/$f2/${f2}_translated_cds.faa.gz -P phage_data/$f2
done
