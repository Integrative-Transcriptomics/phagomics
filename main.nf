#!/usr/bin/env nextflow

include { FILTER                } from "./bin/subworkflows/filter.nf"
include { CLUSTER               } from "./bin/subworkflows/cluster.nf"
include { STRUCTURE_PREDICITON  } from "./bin/subworkflows/structure.nf"
include { SEARCH                } from "./bin/subworkflows/search.nf"
include { INTERPROSCAN          } from "./bin/subworkflows/interproscan.nf"
include { VALIDATE              } from "./bin/subworkflows/validate.nf"
include { REPORT_NEW            } from "./bin/subworkflows/output_new.nf"
include { GFF                   } from "./bin/subworkflows/writeGff.nf"

def helpMsg() {
    println """
        Usage:
        The standard command to run the pipeline is as follows:
        nextflow run main.nf --foldseekdb=databaseDirectory/db/db_name --colabdb=databaseDirectory/colabfolddb --email=valid-E-mail-adress

        Mandatory arguments:
            --foldseekdb                Database used for Foldseek structural alignment (db_name is given during database creation)
            --email                     A valid E-mail is required for the InterProScan API 
            --colabdb                   Directory of ColabFold search databases. Ignored when --wsl is set
            OR
            --wsl                       Sets environmental flags to run the pipeline on Windows Subsystem for Linux. 
                                        Also runs structure prediction on public ColabFold servers, in case no colabdb is available [false]

        Optional arguments:
            --features                  Directory for feature GFF files [./phage_data/*/*.gff]
            --proteins                  Directory for protein FAA files [./phage_data/*/*.faa]
            --outDir                    Output directory [results]
            --maxProteinLength          Maximum length of proteins [1500]
            --validdb                   Internal validation database. If not set internal validation will be skipped.
            --workDir                   Nextflow work directory [work]
            --cleanup                   Automatically delete work directory after successful pipeline run [true]
            --help                      Show this message
    """
}

workflow {

    if (params.help) {
        helpMsg()
        exit 0
    }

    ///
    /// IMPORT & FILTERING
    ///

    // Import protein files and limit length to default:1500 residues
    Channel.fromPath( params.proteins, checkIfExists: true )
    | splitFasta( record: [id:true, desc:true, seqString:true] )
    | filter { record -> record.seqString.length() < params.maxProteinLength }
    | map{ it -> [id:it.id.replace("lcl|", ""), desc:it.desc, seqString:it.seqString]}  // remove lcl| from start of id
    | set { ch_allProteins }

    FILTER( ch_allProteins )


    ///
    /// CLUSTERING
    ///

    CLUSTER( FILTER.out.allProteins )


    ///
    /// INTERPROSCAN
    ///

    INTERPROSCAN( CLUSTER.out.splitClusterReps )


    ///
    /// STRUCTURE PREDICTION
    ///

    STRUCTURE_PREDICITON( CLUSTER.out.splitClusterReps )

    // Test setup for when structures were already predicted
    // Channel.fromPath(["./resultsT4/colabfold/*_rank_001*.pdb", "./resultsT4/colabfold/*_rank_001*.json"]) 
    // | map { it -> 
    //     tuple((it =~ /colabfold\/(.*?)_(unrelaxed|relaxed|scores)/)[0][1], it)
    // }
    // | groupTuple()
    // | map{ id, paths ->
    //     [id:id, 
    //     pdb:  paths.find{ it -> it.toString().endsWith(".pdb") },
    //     json: paths.find{ it -> it.toString().endsWith(".json")}] }
    // | set{ temp }


    ///
    /// FOLDSEEK SEARCH
    /// 

    STRUCTURE_PREDICITON.out
    // temp
    | branch { it ->
        known:      it =~ /.*_known.*/
        unknown:    it =~ /.*_unknown.*/
    }
    | set { structures }

    SEARCH( structures.unknown )


    ///
    /// VALIDATION & OUTPUT
    /// 


    if( params.validdb ) {
        VALIDATE( structures.known )
    }

    REPORT_NEW ( 
        SEARCH.out, 
        CLUSTER.out.clusterMembers,
        FILTER.out.proteinDescriptions,
        CLUSTER.out.allClusterReps,
        INTERPROSCAN.out
    )

    GFF ( REPORT_NEW.out )
}

workflow.onComplete {
    workflow.workDir.deleteDir()

    def message = """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : Auto delete
        outDir      : ${params.outDir}
        exit status : ${workflow.exitStatus}
        """

    println message
}