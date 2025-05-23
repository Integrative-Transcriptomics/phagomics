#!/usr/bin/env nextflow

// mount for database
if( (params.foldseekdb ==  ".") )
    throw new Exception("Missing foldseek database(s)!")
database = file(params.foldseekdb)

if( params.validdb ) {
    validateDB = file(params.validdb)
}

// Input: [gene-id, path] rank 1 predicted structures, .pdb format
// Output: Foldseek aln file
// aln = TSV query,target,fident,evalue,bits
// https://github.com/steineggerlab/foldseek?tab=readme-ov-file#alignment-mode <- OUTPUTFIELDS DESCRIBED HERE
process foldseek {
    //debug true
    publishDir "$params.outDir/foldseek", mode: 'copy'
    maxForks 6
    containerOptions { "--rm" }

    input:
    tuple val(id), path(path), path(json)

    output:
    tuple val(id), path("*.tsv"), path(json)

    script:
    """
    foldseek easy-search "$path" "$database" ${id}.tsv tmp --alignment-type 1 --tmscore-threshold 0.5 \
    --format-output "query,target,qstart,qend,tstart,tend,prob,alntmscore,evalue"
    """
}

/*
*   Searches structure representatives of known proteins against validationDB
*   validationDB = DB of known proteins (all, not just reps) 
*   and 50 each random proteins from model organisms (see DB readme).
*   This is to validate that cluster reps are processed correctly
*   and that the pipeline can find the correct proteins in a DB.
*   Returns channel of alignment files
*/  
process foldseekValidate {
    // Can't have TM-score (alntmscore) with current BaselDB setup
    //debug true
    publishDir "$params.outDir/foldseek/validate", mode: 'copy'
    maxForks 1 // to run locally with limited memory
    containerOptions { "--rm" }

    input:
    tuple val(id), path(structures), path(json)

    output:
    tuple val(id), path("*.tsv"), path(json)

    script:
    """
    foldseek easy-search "$structures" "$validateDB" ${id}.tsv tmp \
    --format-output "query,target,fident,evalue,bits" 
    """
}