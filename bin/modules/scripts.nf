#!/usr/bin/env nextflow

///
/// All miscellaneous scripts go here
///

process clusterMembers {
    containerOptions { "--rm" }
    
    input:
    val(file)

    output:
    path("clusterNumbers.csv"), emit: result

    script:
    """
    #!/usr/bin/env python

    ### Takes a .tsv file and keeps rows with a evalue of < 1 and TM-score > 0.4
    ### csv fields: ...

    import csv
    from collections import defaultdict

    input = "${file}"
    output = "clusterNumbers.csv"

    cluster_count = defaultdict(int)

    with open(input, newline='') as infile, open(output, "w") as outfile:
        reader = csv.reader(infile,  delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        for row in reader:
            cluster_count[row[0]] += 1
        
        writer.writerow(["cluster", "members"])
        for cluster, count in cluster_count.items():
            writer.writerow([cluster, count])


    """
}

/** 
*   Takes channel of foldseek aln files and cluster Member file.
*   Writes top 3 entries sorted my TMscore. Looks for cluster member number 
*   and adds it to output fields.
*   Returns a .tsv file with the fields: 
*   "target", "function", "qstart", "qend", "tstart", "tend", "prob", "alntmscore", "evalue", "lddt"
*/ 
process generateFoldseekReport {
    debug true
    containerOptions { "--rm" }

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

/** 
*   Takes channel of foldseek aln files.
*   Writes entries with an findent > 0.8 (based on nothing). 
*   If there is no entry with fident > 0.8 it will write it still. 
*   Checks if target is same as query and writes to isQuery.
*   Returns a .tsv file with the fields: query,target,evalue,fident,bits,isQuery
*/ 
process generateValidationReport {
    debug true
    containerOptions { "--rm" }

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

// This returns a report per protein
process foldseekReport {
    debug true
    publishDir "$params.outDir/reports", mode: 'copy'
    errorStrategy 'retry'
    containerOptions { "--rm" }

    input:
    tuple val(id), path(file), path(json)

    when:
    file.size() > 0 //Ignote empty foldseek results

    output:
    path("*_fs.json")

    script:
    """
    #!/usr/bin/env python

    import requests, sys, re
    import pandas as pd
    import xml.etree.ElementTree as ET
    import json

    # generate foldseek part for outputReport JSON file.
    # Uses fields from alnfile aswell as getting protein function from uniprotAPI.
    # alnfile fields: query,target,qstart,qend,tstart,tend,qtmscore,prob,alntmscore,evalue,lddt
    def generateFoldseekReport(alnfile):
        df = pd.read_csv(alnfile, 
                        sep='\t', 
                        names=["query", "target", "qstart", "qend", "tstart", "tend", "prob", "alntmscore", "evalue"])
        
        # get top hit for each criteria
        topHits = [pd.DataFrame([df.loc[df[criteria].idxmax()]]) for criteria in ["evalue", "prob", "alntmscore"]]
        
        # applies regex to all dataframes in topHits, extracting the Alphafold-ID in "target"
        pattern = r"(?<=AF-)[A-Z0-9]+(?=-F)"
        topHitsId = map(lambda hit: re.search(pattern, hit.iat[0, 1]).group(), topHits)

        # API call for protein name
        protNames = []
        namespace = {'ns': 'https://uniprot.org/uniprot'}
        for elem in topHitsId:
            requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{elem}"
            r = requests.get(requestURL, headers={ "Accept" : "application/xml"})

            if not r.ok:
                name = None
            else:
                # create tree element from API response XML
                root = ET.fromstring(r.text)
                # UniProt entries can have a "recommendedName" and "submittedName"
                name = root.find(".//ns:protein/ns:recommendedName/ns:fullName", namespace)
                if name is None: name = root.find(".//ns:protein/ns:submittedName/ns:fullName", namespace)
                
            protNames.append(name.text) if name is not None else protNames.append("None")

        # Put top hits into new df and insert name column. Convert to JSON format
        newdf = pd.concat(topHits)
        newdf.insert(2, "name", protNames, True)
        newdf.insert(0, "criteria", ["best e-value", "best probability", "best aln-score"], True)

        # Format to json -> format to proper json format -> each json object (report) to own file
        objs = json.loads(newdf.to_json(orient='records', double_precision=2))
        id = re.match(r"^(.*?known)", objs[0]["query"]).group(1) # regex removes _relaxed... from coldabfold
        with open(f"{id}_fs.json", "w") as fh:
            fh.write('{\\n')
            fh.write('    "method":"foldseek",\\n')
            fh.write('    "hits": ')
            fh.write(json.dumps(objs, indent=4))
            fh.write('\\n}')

    generateFoldseekReport("$file")
    """
}

// This returns a report per protein
// Report can be with method "clusterMember" or "clusterRep"
// clusterMember: A non-rep member is known, function applies to all cluster members
// clusterRep: Rep is known, function applies to all members, 
//  but also creates seperate reports for each member (TODO: Fix this behaviour) 
process clusterReport {
    debug true
    publishDir "$params.outDir/reports", mode: 'copy'
    containerOptions { "--rm" }

    input:
    path(clusterfile)
    path(proteinDescriptionsfile)

    output:
    path("*_cl.json")

    script:
    """
    #!/usr/bin/env python

    import re, json
    import pandas as pd

    def generateClusterReport(clusterFile, descriptions):
        # Read cluster file and collect reps and their members
        clusterData = pd.read_csv(clusterFile, delimiter="\t", names=["id", "member"], index_col=False)
        clusterData = clusterData.groupby("id")["member"].apply(list).reset_index()

        # Read descriptions file
        descriptionData = pd.read_csv(descriptions, delimiter="\t", names=["id", "desc"], index_col=False)
        # Extract protein function from description
        f = lambda x: re.search(r"protein=([^\\]].[^\\]]*)", x).group(1)
        descriptionData["desc"] = descriptionData["desc"].apply(f)
        descriptionData.set_index('id', inplace=True)

        # Add function for each rep
        df = descriptionData.merge(clusterData, on="id", how="inner")

        # If any member is known and there are unknown members -> assign known function to unknown members
        outputDf = pd.DataFrame(columns=["method", "target", "clusterRep", "members", "function", "memberCount", "knownCount", "unknownCount"])
        for index, row in df.iterrows():
            rep     = row["id"]
            members = row["member"]
            knownCount, unknownCount = 0, 0
            for member in members:
                if member.endswith("_known"): knownCount += 1
                else: unknownCount += 1
        
            # Case singleton cluster
            if len(members) == 1:
                pass
            # Case non-Singleton cluster
            else:
                # If there are known and unknown members in the same cluster
                if any(member.endswith("_known") for member in members) and any(member.endswith("_unknown") for member in members):
                    # If the rep is the known member, assign function to other unknown members
                    if rep == members[0] and rep.endswith("_known"):
                        for member in members[1:]:
                            outputDf.loc[len(outputDf)] = [
                                "clusterRep", member, rep, members,
                                descriptionData.iloc[descriptionData.index.get_loc(rep), 0],
                                len(members), knownCount, unknownCount
                            ]
                    else:
                        # If rep is not the known member
                        known = None
                        # Get first known member
                        for member in members:
                            if member.endswith("_known"):
                                known = member
                                break
                        # If a member is unknown assign the known function
                        for member in members:
                            if member.endswith("_unknown"):
                                outputDf.loc[len(outputDf)] = [
                                    "clusterMember", member, rep, members,
                                    descriptionData.iloc[descriptionData.index.get_loc(known), 0],
                                    len(members), knownCount, unknownCount
                                ]

        # Format to json -> format to proper json format -> each json object (report) to own file
        objs = json.loads(outputDf.to_json(orient='records', double_precision=2))
        for obj in objs:
            id = obj["target"]
            with open(f"{id}_cl.json", "w") as fh:
                fh.write(json.dumps(obj, indent=4))


    generateClusterReport("$clusterfile", "$proteinDescriptionsfile")
    """
}

// Returns a report per protein
// function is the originally assigned annotation
// excludes "hypothetical protein" annotation
process postulatedReport {
    debug true
    publishDir "$params.outDir/reports", mode: 'copy'
    containerOptions { "--rm" }

    input:
    path(proteinDescriptionsfile)

    output:
    path("*_ps.json")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import re, json

    def postulatedReport(descriptions):
        # Read descriptions file
        descriptionData = pd.read_csv(descriptions, delimiter="\t", names=["id", "desc"], index_col=False)
        # Extract protein function from description
        f = lambda x: re.search(r"protein=([^\\]].[^\\]]*)", x).group(1)
        descriptionData["desc"] = descriptionData["desc"].apply(f)

        outputDf = pd.DataFrame(columns=["method", "target", "function"])
        for index, row in descriptionData.iterrows():
            if row["id"].endswith("_unknown") and row["desc"] != "hypothetical protein":
                outputDf.loc[len(outputDf)] = ["postulated function", row["id"], row["desc"]]

        # Format to json -> format to proper json format -> each json object (report) to own file
        objs = json.loads(outputDf.to_json(orient='records', double_precision=2))
        for obj in objs:
            id = obj["target"]
            with open(f"{id}_ps.json", "w") as fh:
                fh.write(json.dumps(obj, indent=4))

    postulatedReport("$proteinDescriptionsfile")
    """  
}

process mergeReports {
    publishDir "$params.outDir/reports/final", mode: 'copy'
    debug true
    containerOptions { "--rm" }

    input:
    tuple val(id), path(path)

    output:
    path("${id}_report.json")


    script:
    """
    #!/usr/bin/env python

    from collections import defaultdict
    import json

    resultReport = []
    proteinList = defaultdict(list)

    # Group same proteins together, so we can create a report per protein
    for item in "$path".split():
        key = item[0:-8]
        proteinList[key].append(item)

    # Write predictions (reports) for each protein to json object
    for sameProteinReports in proteinList:
        proteinEntry = {"protein": sameProteinReports, "predictions": []}

        for report in proteinList[sameProteinReports]:
            with open(report, "r") as infile:
                proteinEntry["predictions"].append(json.loads(infile.read()))

        resultReport.append(proteinEntry)
    
    with open("${id}_report.json", "w") as fh:
        fh.write(json.dumps(resultReport, indent=4))    
    """
}