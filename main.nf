#!/usr/bin/env nextflow

include { FILTER                } from "./bin/subworkflows/filter.nf"
include { CLUSTER               } from "./bin/subworkflows/cluster.nf"
include { STRUCTURE_PREDICITON  } from "./bin/subworkflows/structure.nf"
include { SEARCH                } from "./bin/subworkflows/search.nf"
include { INTERPROSCAN          } from "./bin/subworkflows/interproscan.nf"
include { VALIDATE              } from "./bin/subworkflows/validate.nf"
include { REPORT_NEW            } from "./bin/subworkflows/output_new.nf"
include { GFF                   } from "./bin/subworkflows/writeGff.nf"


workflow {

    ///
    /// IMPORT & FILTERING
    ///

    // Import protein files and limit length to 1500 residues
    // Channel.fromPath( params.proteins, checkIfExists: true )
    // | splitFasta( record: [id:true, desc:true, seqString:true] )
    // | filter { record -> record.seqString.length() < params.maxProteinLength }
    // | map{ it -> [id:it.id.replace("lcl|", ""), desc:it.desc, seqString:it.seqString]} // remove lcl| from start of id
    // | set { ch_allProteins }

    // FILTER( ch_allProteins )


    ///
    /// CLUSTERING
    ///

    // CLUSTER( FILTER.out.allProteins )


    ///
    /// INTERPROSCAN
    ///

    // INTERPROSCAN( CLUSTER.out.splitClusterReps )

    //def my = file("resultsFull/protein_clusters/clu2_rep_seq.fasta")
    Channel.fromPath("resultsFull/protein_clusters/clu2_rep_seq.fasta")
    | set{ my }

    INTERPROSCAN( my )


    ///
    /// STRUCTURE PREDICTION
    ///

    // STRUCTURE_PREDICITON( CLUSTER.out.splitClusterReps )

    // Test setup
    Channel.fromPath(["./resultsFull/colabfold/*_rank_001*.pdb", "./resultsFull/colabfold/*_rank_001*.json"]) 
    | map { it -> 
        tuple((it =~ /colabfold\/(.*?)_(unrelaxed|relaxed|scores)/)[0][1], it)
    }
    | groupTuple()
    | map{ id, paths ->
        [id:id, 
        pdb:  paths.find{ it -> it.toString().endsWith(".pdb") },
        json: paths.find{ it -> it.toString().endsWith(".json")}] }
    | set{ temp }


    ///
    /// FOLDSEEK SEARCH
    /// 

    // STRUCTURE_PREDICITON.out
    temp
    | branch { it ->
        known:      it =~ /.*_known.*/
        unknown:    it =~ /.*_unknown.*/
    }
    | set { structures }

    SEARCH( structures.unknown )


    ///
    /// VALIDATION & OUTPUT
    /// 

    VALIDATE( structures.known )

    // REPORT_NEW ( 
    //     SEARCH.out, 
    //     CLUSTER.out.clusterMembers,
    //     FILTER.out.proteinDescriptions,
    //     CLUSTER.out.allClusterReps,
    //     INTERPROSCAN.out
    // )

    REPORT_NEW ( 
        SEARCH.out, 
        file("resultsFull/protein_clusters/clu2_cluster.tsv"),
        file("resultsFull/proteins/proteinDescriptions.tsv"),
        file("resultsFull/protein_clusters/clu2_rep_seq.fasta"),
        INTERPROSCAN.out
    )

    GFF ( REPORT_NEW.out )
}   

workflow.onComplete {
    workflow.workDir.deleteDir()
    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : Auto delete
        outDir      : ${params.outDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed:       ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}