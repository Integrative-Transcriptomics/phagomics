#!/usr/bin/env nextflow

// mount for database
if( params.colabdb ==  "." && !params.wsl)
    throw new Exception("Missing colabfold database!")
database = file(params.colabdb)

process colabfold_batch {
    // ONLY PASSES RANK 1 PREDICTIONS
    //debug true
    publishDir "$params.outDir", mode: 'copy'
    containerOptions { "--rm --runtime=nvidia --gpus all" }
    maxForks 1

    input:
    val(path)

    output:
    tuple path("colabfold/*_relaxed*_001_*.pdb"), path("colabfold/*_001_*.json")

    script:
    """
    colabfold_search $path $database msas
    colabfold_batch msas colabfold --amber --use-gpu-relax
    """
}

process colabfold_batch_wsl {
    // additional flags for WSL support and less recycles
    //debug true
    publishDir "$params.outDir", mode: 'copy'
    containerOptions { "--rm --runtime=nvidia --gpus 1" }
    maxForks 1

    input:
    val(path)

    output:
    tuple path("colabfold/*_relaxed*_001_*.pdb"), path("colabfold/*_001_*.json")

    script:
    """
    export TF_FORCE_UNIFIED_MEMORY="1"
    export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
    export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
    export TF_FORCE_GPU_ALLOW_GROWTH="true"
    colabfold_batch $path colabfold --msa-only
    colabfold_batch $path colabfold --num-recycle 0 --amber --use-gpu-relax
    """
}