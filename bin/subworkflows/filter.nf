#!/usr/bin/env nextflow

workflow FILTER {
    take:
        proteins

    main: 
        // Split into known and unknown proteins
        proteins
        | branch { record ->
            unknown: record.desc =~ /\[protein=[^\]]*(hypothetical|putative|postulated)[^\]]*\]/
            known: true }
        | set{ ch_proteins }
        ///

        // Make seperate channels and files for sequences and descriptions (used in output report)
        ch_proteins.unknown
        | map{ it -> "${it.id}_unknown\t$it.desc" } 
        | set{ ch_proteins_unknown_descs }

        ch_proteins.unknown
        | map{ it -> ">${it.id}_unknown\n$it.seqString" } 
        | set{ ch_proteins_unknown_seqs }

        ch_proteins_unknown_seqs
        | collectFile( name:'unknownProteins.faa', newLine: true, storeDir: "$params.outDir/proteins" )
        | set{ unknownProteins }
        ///

        // Same for known proteins
        ch_proteins.known
        | map{ it -> "${it.id}_known\t$it.desc" } 
        | set{ ch_proteins_known_descs }

        ch_proteins.known
        | map{ it -> ">${it.id}_known\n$it.seqString" } 
        | set{ ch_proteins_known_seqs }

        ch_proteins_known_seqs
        | collectFile( name:'knownProteins.faa', newLine: true, storeDir: "$params.outDir/proteins" )
        ///

        // Collect all protein sequences into a file and channel
        ch_proteins_unknown_seqs
        | concat( ch_proteins_known_seqs )
        | collectFile( name:'allProteins.faa', newLine: true, storeDir: "$params.outDir/proteins" ) 
        | set{ allProteins }
        ///

        // Collect all protein descriptions into a file and channel
        ch_proteins_unknown_descs
        | concat( ch_proteins_known_descs )
        | collectFile( name:'proteinDescriptions.tsv', newLine: true, storeDir: "$params.outDir/proteins" ) 
        | set{ proteinDescriptions }
        ///

    emit:
        allProteins
        proteinDescriptions
        unknownProteins
}