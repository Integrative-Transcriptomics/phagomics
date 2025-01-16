#!/usr/bin/env nextflow

// mount for database
database = file(params.db)

// Input: [gene-id, path] rank 1 predicted structures, .pdb format
// Output: Top 3 predicted functions found in AF/proteome [gene-id, alignment file]
// aln = TSV query,target,fident,alnlen,qlen,tlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore
// use Alphafold/SwissProt for now because its small to download
process foldseek {
    //debug true
    publishDir "$params.outDir/foldseek", mode: 'copy', saveAs: {aln -> "${id}_aln.tsv"}

    input:
    tuple val(id), path(path)

    output:
    tuple val(id), path(aln)

    script:
    """
    foldseek easy-search $path $database/db aln tmp
    """
    //cd BFVD
    //--format-output "query,target,fident,alnlen,qlen,tlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore"
}