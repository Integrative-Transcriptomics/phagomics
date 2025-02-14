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

        // Add _unknown tag to all sequences and
        // collect all unknown proteins into file
        ch_proteins.unknown
        | map{ it -> ">${it.id}_unknown \n $it.seqString" }
        | set{ ch_proteins_unknown }
        
        ch_proteins_unknown
        | collectFile( name:'unknownProteins.faa', newLine: true, storeDir: "$params.outDir/proteins" ) 

        // Same for known proteins
        ch_proteins.known
        | map{ it -> ">${it.id}_known \n $it.seqString" }
        | set{ ch_proteins_known }
        
        ch_proteins_known
        | collectFile( name:'knownProteins.faa', newLine: true, storeDir: "$params.outDir/proteins" )

        // Collect all proteins into a file and channel
        ch_proteins_unknown
        | concat( ch_proteins_known )
        | collectFile( name:'allProteins.faa', newLine: true, storeDir: "$params.outDir/proteins" ) 
        | set{ allProteins }

    emit:
        allProteins
}