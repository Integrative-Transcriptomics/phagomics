#!/usr/bin/env nextflow

// mount for database
database = file(params.db)
validateDB = file(params.validDB)

// Input: [gene-id, path] rank 1 predicted structures, .pdb format
// Output: Top 3 predicted functions found in AF/proteome [gene-id, alignment file]
// aln = TSV query,target,fident,alnlen,qlen,tlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore
// use Alphafold/SwissProt for now because its small to download
process foldseek {
    //debug true
    publishDir "$params.outDir/foldseek", mode: 'copy', saveAs: {aln -> "${id}_aln.tsv"}
    maxForks 1

    input:
    tuple val(id), path(path)

    output:
    tuple val(id), path(aln)

    script:
    """
    foldseek easy-search $path $database/db aln tmp \
    --format-output "query,target,fident,alnlen,qlen,tlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore"
    """
}

process foldseekValidate {
    //debug true
    publishDir "$params.outDir/foldseek/validate", mode: 'copy', saveAs: {aln -> "${id}_aln.tsv"}
    maxForks 1 // to run locally with limited memory

    input:
    tuple val(id), path(structures)

    output:
    path(aln), emit: file

    script:
    """
    foldseek easy-search $structures $validateDB aln tmp \
    --format-output "query,target,fident,evalue,bits" 
    """
}

process generateFoldseekReport {
    debug true
    input:
    path(file)

    output:
    path("report.tsv"), emit: report

    script:
    """
    #!/usr/bin/env python

    ### Takes a .tsv file and keeps rows with a fident of > 0.8.
    ### Adds a column that is true/false depending if the query is the same as target

    import csv

    input = "${file}"
    output = "report.tsv"

    with open(input, newline='') as infile, open(output, "w") as outfile:
        reader = csv.reader(infile,  delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        for row in reader:
            query = row[0]
            target = row[1]
            if float(row[2]) > 0.8:
                if target in query:
                    writer.writerow(row + ["true"])
                else:
                    writer.writerow(row + ["false"])
            else:
                break

        writer.writerow(" ")
    """

}