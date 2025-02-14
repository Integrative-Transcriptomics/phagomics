#!/usr/bin/env nextflow

include { mmseqscluster       as cluster}           from "../modules/cluster.nf"
include { mmseqsclusterRefine as clusterRefine}     from "../modules/cluster.nf"

workflow CLUSTER {
    take:
        allProteins
    
    main:
        cluster( allProteins )
        | clusterRefine

        // split into files of 500 sequences
        clusterRefine.out.repsRefined
        | splitFasta( by: 10, file: true )
        | take( 2 )
        | set{ splitClusterReps }

    emit:
        splitClusterReps
        clusterMembers = clusterRefine.out.cluMembers
}