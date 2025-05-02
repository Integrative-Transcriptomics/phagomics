#!/usr/bin/env nextflow

process mmseqscluster {
    debug false
    publishDir "$params.outDir", mode: 'copy'
    containerOptions { "--rm" }

    input: 
    path(prot_db)

    output:
    path("protein_clusters/clu1_rep_seq.fasta"), emit: reps
    path("protein_clusters/clu1_cluster.tsv")

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
        --min-seq-id 0.9
    """
}

process mmseqsclusterRefine {
    debug false
    publishDir "$params.outDir", mode: 'copy'
    containerOptions { "--rm" }

    input: 
    path(prot_reps)

    output:
    path("protein_clusters/clu2_rep_seq.fasta"), emit: repsRefined
    path("protein_clusters/clu2_cluster.tsv"), emit: cluMembers

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
        --min-seq-id 0.6
    """
}