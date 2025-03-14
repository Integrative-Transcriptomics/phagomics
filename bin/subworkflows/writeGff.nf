#!/usr/bin/env nextflow

include { writeGff } from "../modules/scripts.nf"

workflow GFF {
    take:
        report
    
    main:
        Channel.fromPath( params.features )
        | map{ it -> tuple(it.name[0..-14], it) }
        | concat( report )
        | groupTuple() // [id, [features.gff, report.json]]
        | map{ it -> tuple(it[0], it[1][0], it[1][1])}
        | filter{ it[2] != null } // during testing some .json reports were not there
        | writeGff



    
}