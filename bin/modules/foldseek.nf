#!/usr/bin/env nextflow

// mount for database
if( (params.foldseekdb ==  ".") || (params.validdb ==  "."))
    throw new Exception("Missing foldseek database(s)!")
database = file(params.foldseekdb)
validateDB = file(params.validdb)

// Input: [gene-id, path] rank 1 predicted structures, .pdb format
// Output: Foldseek aln file
// aln = TSV query,target,fident,evalue,bits
// https://github.com/steineggerlab/foldseek?tab=readme-ov-file#alignment-mode <- OUTPUTFIELDS DESCRIBED HERE
process foldseek {
    //debug true
    publishDir "$params.outDir/foldseek", mode: 'copy', saveAs: { file -> file.endsWith("aln") ? "${id}_aln.tsv" : file }
    maxForks 4

    input:
    tuple val(id), path(path), path(json)

    output:
    tuple val(id), path("aln"), path(json)

    script:
    """
    foldseek easy-search "$path" "$database" aln tmp --alignment-type 1 \
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
    maxForks 4 // to run locally with limited memory

    input:
    tuple val(id), path(structures), path(json)

    output:
    tuple val(id), path(aln), path(json)

    script:
    """
    foldseek easy-search "$structures" "$validateDB" aln tmp \
    --format-output "query,target,fident,evalue,bits" 
    """
}