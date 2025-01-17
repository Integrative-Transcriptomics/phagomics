#!/usr/bin/env nextflow

include { filterFunction }                      from "./bin/filter.nf"
include { annotateGenome as annotateGenome1 }   from "./bin/annotateGenomes.nf"
include { annotateGenome as annotateGenome2 }   from "./bin/annotateGenomes.nf"
include { addFlag as addFlag1 }                 from "./bin/annotateGenomes.nf"
include { addFlag as addFlag2 }                 from "./bin/annotateGenomes.nf"
include { mmseqscluster }                       from "./bin/cluster.nf"
include { mmseqscluster_refine }                from "./bin/cluster.nf"
include { colabfold_search }                    from "./bin/colabfold.nf"
include { colabfold_batch }                     from "./bin/colabfold.nf"
include { colabfold_batch_wsl }                 from "./bin/colabfold.nf"
include { foldseek }                            from "./bin/foldseek.nf"



workflow {

    // import feature data from .gff
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

    ch_filtered.hypothetical
    | map { id, files ->
        [id.id, "hypothetical", files]
    }
    | set{ ch_hypothetical }


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

    ch_hypothetical
    | join( ch_genomes )
    | map{ id, type, path1, path2 -> 
        [[id:id, type:type], path1, path2]
    }
    | set{ ch_joinH }
    
    // Run annotations on both channels
    annotateGenome1( ch_joinF )
    | addFlag1
    | set{ ch_known }

    annotateGenome2( ch_joinH )
    | addFlag2
    | set{ ch_unknown }


    ///
    /// CLUSTERING
    ///

    // potentially add flag to allow input of premade protein seqs db (.fasta format)

    // Concat all functional proteins of phages into prot_seqs_k/u and run clustering on them.
    // combine known and unkown here
    ch_known
    | concat(ch_unknown )
    | collectFile( name: 'prots' ) { it[1] }
    | set { prot_seqs }


    // run clustering
    mmseqscluster( prot_seqs )
    | mmseqscluster_refine
    | set { cluster_reps_refined }

    cluster_reps_refined.reps_refined
    /// use only 100 sequences for testing purposes
    | splitFasta()
    | take(100)
    | collectFile ( name: 'clust.fasta' )
    ///
    | set { cluster_reps_refined }


    ///
    /// STRUCTURE PREDICTION
    ///

    if (params.wsl) {
        colabfold_batch_wsl( cluster_reps_refined )
        | flatten
        | map{ it -> 
            // extract gene ID from colabfold output (=~ find operator with matching regex)
            [id:(it =~ /gene-(.*?)_unrelaxed/)[0][1], path:it]
        }
        | set{ ch_structures }
    } else {
        colabfold_batch( cluster_reps_refined )
        | flatten
        | map{ it -> 
            [id:(it =~ /gene-(.*?)_unrelaxed/)[0][1], path:it]
        }
        | set{ ch_structures }
    }


    ///
    /// FOLDSEEK SEARCH
    /// 
    
    foldseek( ch_structures)
}   