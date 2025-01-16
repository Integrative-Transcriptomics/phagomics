#!/usr/bin/env nextflow

process annotateGenome {
    // Generate protein sequences using gffread with a genome fasta and a features gff. 
    // Takes a stream of [[id, type], path/to/fasta, path/to/gff]
    // Return meta and annotated protein sequences
    debug true
    //publishDir "$params.outDir/proteins", mode: 'copy'

    input:
    tuple val(meta), path(annotation), path(genome)

    output:
    tuple val(meta), path("*_*_proteins.fasta")
    // tuple val(meta), path("hypothetical/*__hypothetical_proteins.fasta"), emit: hypothetical_protein
   
    script:
    """
    gffread -g ${genome} -y ${meta.id}_${meta.type}_proteins.fasta ${annotation}
    """
}

process addFlag {
    publishDir "$params.outDir/proteins", mode: 'copy'

    input:
    tuple val(meta), path(path)

    output:
    tuple val(meta), path("*/*.fasta")

    script:
    """
    #!/usr/bin/env python

    import argparse, os

    os.makedirs("functional/", exist_ok=True)
    os.makedirs("hypothetical/", exist_ok=True)

    input = "${path}"
    output = "${meta.type}/${meta.id}_${meta.type}_proteins.fasta"
    

    with open(input, "r") as infile, open(output, "w") as outfile:
        tag = "_known" if ("${meta.type}" == "functional") else "_unknown"
        for line in infile:
            line = line.rstrip('\\n')
            if line.startswith('>'):
                line += tag
            
            outfile.write(line + '\\n')
    """
}