#!/usr/bin/env nextflow

include { mmseqscluster       as cluster}           from "../modules/cluster.nf"
include { mmseqsclusterRefine as clusterRefine}     from "../modules/cluster.nf"

workflow CLUSTER {
    take:
        allProteins
    
    main:
        cluster( allProteins ).reps
        | clusterRefine

        // split into files of 200 sequences
        clusterRefine.out.repsRefined
        | splitFasta( by: 200, file: true )
        | set{ splitClusterReps }

    emit:
        splitClusterReps
        clusterMembers = clusterRefine.out.cluMembers
        allClusterReps = clusterRefine.out.repsRefined 
}