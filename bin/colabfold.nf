#!/usr/bin/env nextflow

process colabofold_search {
    debug true
    publishDir "$params.outDir", mode: 'copy'

    input:
    path(path)

    output:
    path("colabfold_search/*")

    script:
    """
    echo "search!"
    """
}

process colabfold_batch {
    debug true
    publishDir "$params.outDir", mode: 'copy'

    input:
    path(path)

    output:
    path("colabfold/*")

    script:
    """
    mkdir colabfold
    echo "batch!"
    """
    // colabfold_batch $path colabfold/

    //colabfold_search --mmseqs /path/to/mmseqs /input.fasta /database msas > && colabfold_batch msas /predictions
}