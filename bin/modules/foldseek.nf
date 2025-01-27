#!/usr/bin/env nextflow

// mount for database
database = file(params.db)
validateDB = file(params.validDB)

// Input: [gene-id, path] rank 1 predicted structures, .pdb format
// Output: Foldseek aln file
// aln = TSV query,target,fident,evalue,bits
process foldseek {
    //debug true
    publishDir "$params.outDir/foldseek", mode: 'copy', saveAs: {aln -> "${id}_aln.tsv"}
    maxForks 4

    input:
    tuple val(id), path(path), path(json)

    output:
    tuple val(id), path(aln), path(json)

    script:
    """
    foldseek easy-search $path $database/db aln tmp \
    --format-output "query,target,fident,evalue,bits"
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
    publishDir "$params.outDir/foldseek/validate", mode: 'copy', saveAs: {aln -> "${id}_aln.tsv"}
    maxForks 4 // to run locally with limited memory

    input:
    tuple val(id), path(structures), path(json)

    output:
    tuple val(id), path(aln), path(json)

    script:
    """
    foldseek easy-search $structures $validateDB aln tmp \
    --format-output "query,target,fident,evalue,bits" 
    """
}

/** 
*   Takes channel of foldseek aln files.
*   Writes entries with an findent > 0.8 (based on nothing). 
*   If there is no entry with fident > 0.8 it will write it still. 
*   Checks if target is same as query and writes to isQuery.
*   Returns a .tsv file with the fields: query,target,evalue,fident,bits,isQuery
*/ 
process generateValidationReport {
    debug true
    input:
    tuple val(id), path(file), path(json)

    output:
    path("validationReport.tsv"), emit: report

    script:
    """
    #!/usr/bin/env python

    ### Takes a .tsv file and keeps rows with a fident of > 0.8.
    ### Adds a column that is true/false depending if the query is the same as target

    import csv

    id = "${id}"
    input = "${file}"
    output = "validationReport.tsv"

    with open(input, newline='') as infile, open(output, "w") as outfile:
        reader = csv.reader(infile,  delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        writer.writerow(['query','target','evalue','fident','bits','isQuery'])
        
        for count, row in enumerate(reader):
            query = row[0]
            target = row[1]

            # make sure id = query
            if id not in query:
                raise Exception("Id is not query!")

            row[0] = id

            if float(row[2]) > 0.8:
                if target in query:
                    writer.writerow(row + ["true"])
                else:
                    writer.writerow(row + ["false"])
            else:
                # if there's no entry with fident > 0.8 write it anyways
                if count == 0:
                    if target in query:
                        writer.writerow(row + ["true"])
                    else:
                        writer.writerow(row + ["false"])
                break

        writer.writerow(" ")
    """
}

/** 
*   Takes channel of foldseek aln files and cluster Member file.
*   Writes top 3 entries with an evalue < 0.1. Looks for cluster member number 
*   and adds it to output fields.
*   Returns a .tsv file with the fields: query,clunum,target,evalue,fident,bits
*/ 
process generateFoldseekReport {
    debug true
    input:
    tuple val(id), path(alnfile), path(json)
    val(clusterMembers)

    output:
    path("reports.csv"), emit: reports

    script:
    """
    #!/usr/bin/env python

    import csv, json

    id = "${id}"
    input = "${alnfile}"
    output = "reports.csv"
    input_json = "${json}"
    clusters = "${clusterMembers}"

    # open all required files
    with open(input, "r")      as infile,  \
         open(output, "w")     as outfile, \
         open(clusters, "r")   as clufile, \
         open(input_json, "r") as jsonfile:

        # create readers/writers
        reader = csv.reader(infile,  delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        readerCluster = csv.reader(clufile, delimiter='\t')
        json_content = json.loads(jsonfile.read())
        clusterNumber = 0
        
        # Write header
        writer.writerow(['query','clunum','pTM','target','evalue','fident','bits'])

        for count, row in enumerate(reader):

            query = row[0]
            target = row[1]
            evalue = row[3]
            #alntmscore = row[5]

            # make sure id = query
            if id not in query:
                raise Exception("Id is not query!")

            if count == 3:
                break
                
            # Get number of cluster members
            for clu in readerCluster:
                if clu[0] in query:
                    clusterNumber = clu[1]

            # Filter on evalue
            if (float(evalue) < 0.1):
                writer.writerow([id, clusterNumber, json_content['ptm'], target, evalue, row[2], row[4]])
            else:
                count -= 1
                break
    """
}