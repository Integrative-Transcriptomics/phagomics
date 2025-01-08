#!/usr/bin/env nextflow

include { filterFunction }                      from "./bin/filter.nf"
include { annotateGenome as annotateGenome1 }   from "./bin/annotateGenomes.nf"
include { annotateGenome as annotateGenome2 }   from "./bin/annotateGenomes.nf"
include { mmseqscluster }                       from "./bin/cluster.nf"
include { mmseqscluster_refine }                from "./bin/cluster.nf"
include { colabfold }                           from "./bin/colabfold.nf"



workflow {

    // Import files and do some preliminary tagging
    // old combined channel of fastas + gffs
    // Channel.fromFilePairs( [params.genomes, params.features], checkIfExists: true )
    // | map { id, files ->
    //    [id:id, paths:files]
    // }
    // | view
    // | set{ ch_genomes }

    Channel.fromPath( params.features, checkIfExists: true) 
    | map{ path ->
        // could skip this because it gets undone at ch_filtered.functional... step
        meta = [id:path.getParent().getName()]
        [meta, path]
    }
    | set{ ch_annot }


    ///
    /// FILTERING
    ///
 
    // Call the filter Function. Emits 2 channels: functional, hypothetical
    filterFunction( ch_annot )
    | set{ ch_filtered }
    
    // Select the functional/hypothetical annotations and create channels
    // Add functional/hypothetical tag to gff files
    // Maybe outsource this to filter.nf (?)
    ch_filtered.functional
    | map { id, files ->
        [id.id, "functional", files]
    }
    | set{ ch_functional }

    // ch_filtered.hypothetical
    // | map { id, files ->
    //     [id.id, "hypothetical", files]
    // }
    // | set{ ch_hypothetical }


    ///
    /// ANNOTATION
    ///

    // Create channel for genome files
    Channel.fromPath( params.genomes, checkIfExists: true) 
    | map{ path ->
        // get id from parent folder
        [path.getParent().getName(), path]
    }
    | set{ ch_genomes }

    // Combine genome files and functional/hypothetical annotations channels,
    // Adds meta map [id:id, type:type]
    ch_functional
    | join( ch_genomes )
    | map{ id, type, path1, path2 -> 
        [[id:id, type:type], path1, path2]
    }
    | set{ ch_joinF }

    // ch_hypothetical
    // | join( ch_genomes )
    // | map{ id, type, path1, path2 -> 
    //     [[id:id, type:type], path1, path2]
    // }
    // | set{ ch_joinH }
    
    // Run annotations on both channels
    annotateGenome1( ch_joinF )
    // | view
    | set{ ch_known }

    // annotateGenome2( ch_joinH )
    // // | view


    ///
    /// CLUSTERING
    ///


    // potentially add flag to allow input of premade protein seqs db (.fasta format)

    // Concat all functional proteins of phages into prot_seqs and run clustering on them.
    prot_seqs = ch_known.collectFile( name: 'prots' ) { it[1] }

    cluster_reps = mmseqscluster( prot_seqs ) 
    cluster_reps_refined = mmseqscluster_refine( cluster_reps )


    ///
    /// STRUCTURE PREDICTION
    ///


    colabfold_search( cluster_reps_refined )
    colabfold_batch( cluster_reps_refined )

}   