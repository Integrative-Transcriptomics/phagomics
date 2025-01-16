#!/usr/bin/env nextflow

process colabfold_search {
    // Will be used for MSA generation
    debug true
    publishDir "$params.outDir", mode: 'copy'

    input:
    path(path)

    // output:
    // path("colabfold_search/*")

    script:
    """
    echo "search!"
    """
}

process colabfold_batch {
    // ONLY PASSES RANK 1 PREDICTIONS AT THE MOMENT!
    //debug true
    publishDir "$params.outDir", mode: 'copy'
    // docker flags, remove to run on CPU
    // containerOptions { "--runtime=nvidia --gpus 1" }

    input:
    val(path)

    output:
    path("colabfold/*_001_*.pdb")

    script:
    """
    colabfold_batch $path colabfold --msa-only
    colabfold_batch $path colabfold
    """
}

process colabfold_batch_wsl {
    // additional flags for WSL support
    //debug true
    publishDir "$params.outDir", mode: 'copy'
    containerOptions { "--runtime=nvidia --gpus 1" }

    input:
    val(path)

    output:
    path("colabfold/*_001_*.pdb")

    script:
    """
    export TF_FORCE_UNIFIED_MEMORY="1"
    export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
    export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
    export TF_FORCE_GPU_ALLOW_GROWTH="true"
    colabfold_batch $path colabfold --msa-only
    colabfold_batch $path colabfold --num-recycle 0
    """
}