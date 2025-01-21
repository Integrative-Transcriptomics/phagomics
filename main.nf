#!/usr/bin/env nextflow

include { filterFunction }                      from "./bin/modules/filter.nf"
include { annotateGenome as annotateGenome1 }   from "./bin/modules/annotateGenomes.nf"
include { annotateGenome as annotateGenome2 }   from "./bin/modules/annotateGenomes.nf"
include { addFlag as addFlag1 }                 from "./bin/modules/annotateGenomes.nf"
include { addFlag as addFlag2 }                 from "./bin/modules/annotateGenomes.nf"
include { mmseqscluster }                       from "./bin/modules/cluster.nf"
include { mmseqscluster_refine }                from "./bin/modules/cluster.nf"
include { colabfold_search }                    from "./bin/modules/colabfold.nf"
include { colabfold_batch }                     from "./bin/modules/colabfold.nf"
include { colabfold_batch_wsl }                 from "./bin/modules/colabfold.nf"
include { foldseek }                            from "./bin/modules/foldseek.nf"

include { validateFoldseek }                    from "./bin/subworkflows/validate.nf" 

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
        //colabfold_batch_wsl( cluster_reps_refined )
        //| flatten
        Channel.fromPath("./100test/*_relaxed*_001_*.pdb") // take files from colabfold local run (not in docker)
        | map { it -> 
            // extract gene ID from colabfold output (=~ find operator with matching regex)
            [id:(it =~ /gene-(.*?)_rank/)[0][1], path:it]
        }
        | set{ ch_structures }
    } else {
        colabfold_batch( cluster_reps_refined )
        | flatten
        | map { it -> 
            [id:(it =~ /gene-(.*?)_rank/)[0][1], path:it]
        }
        | set{ ch_structures }
    }


    ///
    /// FOLDSEEK SEARCH
    /// 
    
    // Run foldseek on all predicted structures, known or unknown against given DB
    foldseek( ch_structures )
    | splitCsv( sep: "\t" )
    | filter { it ->                    // Filter to (from birth of proteins folds)
        it[1][14].toDouble() > 0.4 &&   // TM-score > 0.4 
        it[1][12].toDouble() < 0.001    // e-value  < 0.001
    } 
    // convoluted way to write a csv (from https://github.com/nf-core/sarek/blob/master/subworkflows/local/channel_variant_calling_create_csv/main.nf)
    | collectFile( keepHeader: true, skip: 1,sort: true, storeDir: "${params.outDir}/foldseek" ) { id, csv ->
        query      = csv[0]
        target     = csv[1] 
        evalue     = csv[12]
        alntmscore = csv[14]
        bits       = csv[13]
        ["result.csv", "query,target,evalue,TM-score,bits\n${query},${target},${evalue},${alntmscore},${bits}\n"] }
    
    // Validate known proteins by searching against BaselDB
    ch_structures
    | filter ( ~/.*_(known)_.*/ )
    | validateFoldseek


            // Branch known/unknown. DOn't know if needed
            // ch_structures
            // | branch { it ->
            //     known:      it =~ /.*_(known)_.*/
            //     unknown:    it =~ /.*_(unknown)_.*/
            // }
            // | set { struc }

            // struc.known
            // | view

            // struc.unknown
            // | view
}   