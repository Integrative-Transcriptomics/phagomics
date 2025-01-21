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
    debug false
    publishDir "$params.outDir", mode: 'copy'

    input: 
    path(prot_db)

    output:
    path("protein_clusters/clu1_rep_seq.*"), emit: reps
    //path("protein_clusters/clu1*.*") // get rest of cluster files if needed (requires input change in refinement step)

    script:
    """
    mkdir -p protein_clusters/tmp
    mkdir -p protein_clusters/clu1
    mmseqs easy-cluster \
        $prot_db \
        protein_clusters/clu1 \
        protein_clusters/tmp \
        --max-seqs 50000 \
        -c 0.9 \
        --min-seq-id 0.9 \
        --threads 5
    """
}

process mmseqscluster_refine {
    debug false
    publishDir "$params.outDir", mode: 'copy'

    input: 
    path(prot_reps)

    output:
    path("protein_clusters/clu2_rep_seq.fasta"), emit: reps_refined
    path("protein_clusters/clu2_cluster.tsv")

    script:
    """
    mkdir -p protein_clusters/tmp
    mkdir -p protein_clusters/clu2
    mmseqs easy-cluster \
        $prot_reps \
        protein_clusters/clu2 \
        protein_clusters/tmp \
        --max-seqs 50000 \
        -c 0.75 \
        --min-seq-id 0.6 \
        --threads 5
    """
}

    