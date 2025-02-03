#!/usr/bin/env nextflow

include { filterFunction }                      from "./bin/modules/filter.nf"
include { annotateGenome as annotateGenome1 }   from "./bin/modules/annotateGenomes.nf"
include { annotateGenome as annotateGenome2 }   from "./bin/modules/annotateGenomes.nf"
include { addFlag as addFlag1 }                 from "./bin/modules/annotateGenomes.nf"
include { addFlag as addFlag2 }                 from "./bin/modules/annotateGenomes.nf"
include { mmseqscluster }                       from "./bin/modules/cluster.nf"
include { mmseqscluster_refine }                from "./bin/modules/cluster.nf"
include { colabfold_batch }                     from "./bin/modules/colabfold.nf"
include { colabfold_batch_wsl }                 from "./bin/modules/colabfold.nf"
include { foldseek }                            from "./bin/modules/foldseek.nf"
include { clusterMembers }                      from "./bin/modules/cluster.nf"

include { validateFoldseek }                    from "./bin/subworkflows/validate.nf"
include { outputReport }                        from "./bin/subworkflows/output.nf" 


workflow {


    // Import protein files and cap at length 1500 
    //Channel.fromPath( params.proteins, checkIfExists: true )
    Channel.fromPath( "./phage_data/NC_001604.1/*.faa", checkIfExists: true )
    | splitFasta( record: [id:true, desc:true, seqString:true] )
    | filter { record -> record.seqString.length() < 1500 }
    | set { ch_proteins_all }

    // Split into known and unknown proteins
    ch_proteins_all
    | branch { record ->
        unknown: record.desc =~ /\[protein=[^\]]*(hypothetical|putative|postulated)[^\]]*\]/
        known: true // the rest of the proteins 
    }
    | set{ ch_proteins }

    // Collect all unknown proteins into file
    ch_proteins.unknown
    | map { it -> ">${it.id}_[unknown] \n $it.seqString"}
    | set { ch_proteins_unknown }
    
    ch_proteins_unknown
    | collectFile( name:'unknown_prots.faa', newLine: true ) 
    | set{ proteins_unknown }

    // Collect all known proteins into file
    ch_proteins.known
    | map { it -> ">${it.id}_[known] \n $it.seqString"}
    | set { ch_proteins_known }
    
    ch_proteins_known
    | collectFile( name:'known_prots.faa', newLine: true ) 
    | set{ proteins_known }

    // Collect all proteins into file
    ch_proteins_unknown
    | concat( ch_proteins_known )
    | collectFile( name:'all_prots.faa', newLine: true ) 
    | set{ proteins_all }


    // ALL OLD
    // //old
    // // import feature data from .gff
    // Channel.fromPath( params.features, checkIfExists: true) 
    // | map{ path ->
    //     // could skip this because it gets undone at ch_filtered.functional... step
    //     meta = [id:path.getParent().getName()]
    //     [meta, path]
    // }
    // | set{ ch_annot }
    ///
    /// FILTERING
    ///
 
    // Call the filter Function. Emits 2 channels: functional, hypothetical
    // filterFunction( ch_proteins )
    // | set{ ch_filtered }
    
    // // Select the functional/hypothetical annotations and create channels
    // // Add functional/hypothetical tag to gff files
    // // Maybe outsource this to filter.nf (?)
    // ch_filtered.bedF
    // | map { id, files ->
    //     [id.id, "functional", files]
    // }
    // | set{ ch_functional }

    // ch_filtered.bedH
    // | map { id, files ->
    //     [id.id, "hypothetical", files]
    // }
    // | set{ ch_hypothetical }


    // ///
    // /// ANNOTATION
    // ///

    // // Create channel for genome files
    // Channel.fromPath( params.genomes, checkIfExists: true) 
    // | map{ path ->
    //     // get id from parent folder
    //     [path.getParent().getName(), path]
    // }
    // | set{ ch_genomes }

    // // Combine genome files and functional/hypothetical annotations channels,
    // // Adds meta map [id:id, type:type]
    // ch_functional
    // | join( ch_genomes )
    // | map{ id, type, path1, path2 -> 
    //     [[id:id, type:type], path1, path2]
    // }
    // | set{ ch_joinF }

    // ch_hypothetical
    // | join( ch_genomes )
    // | map{ id, type, path1, path2 -> 
    //     [[id:id, type:type], path1, path2]
    // }
    // | set{ ch_joinH }
    
    // // Run annotations on both channels
    // annotateGenome1( ch_joinF )
    // | addFlag1
    // | set{ ch_known }

    // annotateGenome2( ch_joinH )
    // | addFlag2
    // | set{ ch_unknown }
    
    // OLD ENDS HERE

    ///
    /// CLUSTERING
    ///

    //run clustering
    mmseqscluster( proteins_all )
    | mmseqscluster_refine
    | set { cluster_reps_refined }

    // split into files of 500 sequences
    cluster_reps_refined.reps_refined
    | splitFasta(by: 500, file: true)
    // | take (2)
    // | view
    | set { cluster_reps_refined_sample }


    ///
    /// STRUCTURE PREDICTION
    ///

    if (params.wsl) {
        // testsetup with cached results
        // Channel.fromPath(["./100test/*_rank_001*.pdb", "./100test/*_rank_001*.json"]) 
        // | map { it -> 
        //     tuple((it =~ /100test\/(.*?)_(unrelaxed|relaxed|scores)/)[0][1], it)
        // }
        // Output is colabofold .pdb and .json. The following extracts the id from filename and maps
        // the files according to id
        colabfold_batch_wsl( cluster_reps_refined_sample )
        | flatten
        | map { it -> 
            tuple((it =~ /colabfold\/(.*?)_(relaxed|scores)/)[0][1], it)
        }
        | groupTuple()
        | map { id, paths ->
            [id:id, 
            pdb:  paths.find{ it -> it.toString().endsWith(".pdb") },
            json: paths.find{ it -> it.toString().endsWith(".json")}] 
        }
        | set{ ch_structures }
    } else {
        // Run normally
        colabfold_batch( cluster_reps_refined_sample )
        | flatten
        | map { it -> 
            tuple((it =~ /colabfold\/(.*?)_(relaxed|scores)/)[0][1], it)
        }
        | groupTuple()
        | map { id, paths ->
            [id:id, 
            pdb:  paths.find{ it -> it.toString().endsWith(".pdb") },
            json: paths.find{ it -> it.toString().endsWith(".json")}] 
        }
        | set{ ch_structures }
    }


    ///
    /// FOLDSEEK SEARCH
    /// 

    //Branch known/unknown
    ch_structures
    | branch { it ->
        known:      it =~ /.*_[known]_.*/
        unknown:    it =~ /.*_[unknown]_.*/
    }
    | set { struc }

    // Run validation on known proteins
    struc.known
    | validateFoldseek
    
    // Run all proteins with foldseek
    // Alternative: run on only unknown proteins
    ch_structures
    // struc.unknown
    | foldseek
    | set{ out }

    // ch_structures
    // | view

    // Generate output report
    outputReport( out , cluster_reps_refined.clu_Members )

    
}   
