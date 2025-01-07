#!/usr/bin/env nextflow

// not used. Used for mmseqs cluster (not easy-cluster)
process createDB {
    debug true
    publishDir "$params.outDir/seq_cluster", mode: 'copy'
    
    input:
    path(path)

    output:
    path "protein_db", emit: db

    script:
    """
    mkdir -p protein_db
    mmseqs createdb $path protein_db/db
    """
}

process mmseqscluster {
    debug true
    publishDir "$params.outDir", mode: 'copy'

    input: 
    path(prot_db)

    output:
    path("protein_clusters/clu*.*")

    script:
    """
    mkdir -p protein_clusters/tmp
    mkdir -p protein_clusters/clu
    mmseqs easy-cluster $prot_db protein_clusters/clu protein_clusters/tmp
    """
}