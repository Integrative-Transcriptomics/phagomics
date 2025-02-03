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
    tuple val(meta), path("functional_beds/*.bed"), emit: bedF
    tuple val(meta), path("hypothetical_beds/*.bed"), emit: bedH
   
    script:
    """
    #!/usr/bin/env python
    ### Filter out given words in filterWords from a gff file. Only consideres CDS regions.
    ### Takes a gff file
    ### Creates & writes new files ID_functional.gff and ID_hypothetical.gff
    ### in directories ./functional/ and ./hypothetical/

import argparse, os, re

os.makedirs("functional/", exist_ok=True)
os.makedirs("hypothetical/", exist_ok=True)
os.makedirs("functional_beds/", exist_ok=True)
os.makedirs("hypothetical_beds/", exist_ok=True)

input = "${path}"
outputFunction = "functional/${meta.id}_functional.gff"
outputHypothetical = "hypothetical/${meta.id}_hypothetical.gff"
outputbedF = "functional_beds/${meta.id}_functional.bed"
outputbedH = "hypothetical_beds/${meta.id}_hypothetical.bed"
filterWords = ["hypothetical", "putative", "postulated"]

with open(input, "r") as infile, open(outputFunction, "w") as outfileF, \
open(outputHypothetical, "w") as outfileH, open(outputbedH, "w") as outfilebedH, open(outputbedF, "w") as outfilebedF:

    for row in infile:
        if row.startswith("#"):
            outfileF.write(row)
            outfileH.write(row)
            continue

        columns = row.split("\t")
        chrom = columns[0]
        start = int(columns[3]) - 1
        end = int(columns[4])
        attributes = columns[8]

        # only take CDS segments, edit this to also collect rest
        if columns[2] == "CDS":
            #Get ID from attributes
            match = re.search(r'ID=([^;]+)', attributes)

            if any(x in columns[8] for x in filterWords):
                outfileH.write(row)
                if match:
                    feature_id = match.group(1)
                    # Write to output file in BED format
                    outfilebedH.write(f"{chrom}\\t{start}\\t{end}\\t{feature_id}\\n")
            else:
                outfileF.write(row)
                if match:
                    feature_id = match.group(1)
                    # Write to output file in BED format
                    outfilebedF.write(f"{chrom}\\t{start}\\t{end}\\t{feature_id}\\n")

    """
}

process filterProteins {
    input:
    tuple val(meta), path(path)

    output:
    tuple val(meta), path("known/*_known.gff"), emit: known
    tuple val(meta), path("unknown/*_unknown.gff"), emit: unknown

    script:
    """
    #!/usr/bin/env python

    import argparse, os, re

    os.makedirs("known/", exist_ok=True)
    os.makedirs("unknown/", exist_ok=True)

    input = "${path}"
    outputFunction = "functional/${meta.id}_functional.faa"
    outputHypothetical = "hypothetical/${meta.id}_hypothetical.faa"
    filterWords = ["hypothetical", "putative", "postulated"]

    with open(input, "r") as infile, open(outputFunction, "w") as outfileF, \
    open(outputHypothetical, "w") as outfileH, open(outputbedH, "w") as outfilebedH, open(outputbedF, "w") as outfilebedF:

    for row in infile:
        if row.startswith("#"):
            outfileF.write(row)
            outfileH.write(row)
            continue

        columns = row.split("\t")
        chrom = columns[0]
        start = int(columns[3]) - 1
        end = int(columns[4])
        attributes = columns[8]

        # only take CDS segments, edit this to also collect rest
        if columns[2] == "CDS":
            #Get ID from attributes
            match = re.search(r'ID=([^;]+)', attributes)

            if any(x in columns[8] for x in filterWords):
                outfileH.write(row)
                if match:
                    feature_id = match.group(1)
                    # Write to output file in BED format
                    outfilebedH.write(f"{chrom}\\t{start}\\t{end}\\t{feature_id}\\n")
            else:
                outfileF.write(row)
                if match:
                    feature_id = match.group(1)
                    # Write to output file in BED format
                    outfilebedF.write(f"{chrom}\\t{start}\\t{end}\\t{feature_id}\\n")

    """
}