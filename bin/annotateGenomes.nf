#!/usr/bin/env nextflow

process annotateGenome {
    // Generate protein sequences using gffread with a genome fasta and a features gff. 
    // Takes a stream of [[id, type], path/to/fasta, path/to/gff]
    // Return meta and annotated protein sequences
    debug true
    publishDir "$params.outDir/proteins", mode: 'copy'

    input:
    tuple val(meta), path(annotation), path(genome)

    output:
    tuple val(meta), path("*/*_*_proteins.fasta")
    // tuple val(meta), path("hypothetical/*__hypothetical_proteins.fasta"), emit: hypothetical_protein
   
    script:
    """
    mkdir -p functional
    mkdir -p hypothetical
    gffread -g ${genome} -y ${meta.type}/${meta.id}_${meta.type}_proteins.fasta ${annotation}
    """
}