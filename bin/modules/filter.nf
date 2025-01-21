#!/usr/bin/env nextflow

process filterFunction {
    /*
    *   Filter proteins with known function and unknown function
    *   Creates directory structure of outDir/annotations_filtered/(functional|hypothetical)
    *   Takes a channel of [id, file]
    *   Emits 2 channels: functional and hypothetical
    */

    debug true
    publishDir "$params.outDir/annotations_filtered", mode: 'copy'

    input:
    tuple val(meta), path(path)

    output:
    tuple val(meta), path("functional/*_functional.gff"), emit: functional
    tuple val(meta), path("hypothetical/*_hypothetical.gff"), emit: hypothetical
   
    script:
    """
    #!/usr/bin/env python
    ### Filter out given words in filterWords from a gff file. Only consideres CDS regions.
    ### Takes a gff file
    ### Creates & writes new files ID_functional.gff and ID_hypothetical.gff
    ### in directories ./functional/ and ./hypothetical/

    
    import argparse, os

    os.makedirs("functional/", exist_ok=True)
    os.makedirs("hypothetical/", exist_ok=True)

    input = "${path}"
    outputFunction = "functional/${meta.id}_functional.gff"
    outputHypothetical = "hypothetical/${meta.id}_hypothetical.gff"
    filterWords = ["hypothetical", "putative"] # maybe add "postulated"

    with open(input, "r") as infile, open(outputFunction, "w") as outfileF, open(outputHypothetical, "w") as outfileH:

        for row in infile:
            if row.startswith("#"):
                outfileF.write(row)
                outfileH.write(row)
                continue

            columns = row.split("\t")

            # only take CDS segments, edit this to also collect rest
            if columns[2] == "CDS":
            
                if any(x in columns[8] for x in filterWords):
                    outfileH.write(row)
                else:
                    outfileF.write(row)

    """
}