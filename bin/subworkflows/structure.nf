#!/usr/bin/env nextflow

include { colabfold_batch     } from "../modules/colabfold.nf"
include { colabfold_batch_wsl } from "../modules/colabfold.nf"

workflow STRUCTURE_PREDICITON {
    take:
        clusterReps

    main:
        if( params.wsl ) {
            // testsetup with cached results
            // Channel.fromPath(["./100test/*_rank_001*.pdb", "./100test/*_rank_001*.json"]) 
            // | map { it -> 
            //     tuple((it =~ /100test\/(.*?)_(unrelaxed|relaxed|scores)/)[0][1], it)
            // }
            // Output is colabofold .pdb and .json. The following extracts the id from filename and maps
            // the files according to id
            colabfold_batch_wsl( clusterReps )
            | flatten
            | map{ it -> 
                tuple((it =~ /colabfold\/(.*?)_(relaxed|scores)/)[0][1], it) }
            | groupTuple()
            | map{ id, paths ->
                [id:id, 
                pdb:  paths.find{ it -> it.toString().endsWith(".pdb") },
                json: paths.find{ it -> it.toString().endsWith(".json")}] }
            | set{ structures }
        } else {
            // Run normally
            colabfold_batch_wsl( clusterReps )
            | flatten
            | map{ it -> 
                tuple((it =~ /colabfold\/(.*?)_(relaxed|scores)/)[0][1], it) }
            | groupTuple()
            | map{ id, paths ->
                [id:id, 
                pdb:  paths.find{ it -> it.toString().endsWith(".pdb") },
                json: paths.find{ it -> it.toString().endsWith(".json")}] }
            | set{ structures }
        }

    emit:
        structures
}