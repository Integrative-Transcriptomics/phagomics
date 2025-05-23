#!/usr/bin/env nextflow

///
/// All miscellaneous scripts go here
///

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
    //publishDir "$params.outDir/reports/single", mode: 'copy'
    errorStrategy 'retry'
    containerOptions { "--rm" }

    input:
    tuple val(id), path(file), path(json)

    when:
    file.size() > 0 //Ignore empty foldseek results

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
        newdf.insert(2, "function", protNames, True)
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

// Known and unknown member in same cluster -> assign known funtion to unknown
// Currently does not confirm similarity of proteins in cluster so annotation could be wrong
process clusterReport {
    debug true
    //publishDir "$params.outDir/reports/single", mode: 'copy'
    containerOptions { "--rm" }

    input:
    path(clusterfile)
    path(proteinDescriptionsfile)

    output:
    path("*_cl.json"), optional: true

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

        df = descriptionData.merge(clusterData, on="id", how="inner")
        outputDf = pd.DataFrame(columns=["method", "target", "clusterRep", "members", "function", "memberCount", "knownCount", "unknownCount"])
        for index, row in df.iterrows():
            rep     = row["id"]
            members = row["member"]
            knownCount, unknownCount = 0, 0
            for member in members:
                if member.endswith("_known"): knownCount += 1
                else: unknownCount += 1
        
            # Case singleton cluster
            if len(members) == 1 and not member.endswith("_known"):
                outputDf.loc[len(outputDf)] = [
                    "singletonCluster", member, rep, members,
                    descriptionData.iloc[descriptionData.index.get_loc(members[0]), 0],
                    len(members), knownCount, unknownCount
                ]
            # Case non-Singleton cluster
            else:
                # Go through member list, add entry to outputdf for every unknown member
                if any(member.endswith("_known") for member in members) and any(member.endswith("_unknown") for member in members):
                    # Get first known member
                    known = None
                    for member in members:
                        if member.endswith("_known"):
                            known = member
                            break
                    # If a member is unknown, assign it the known function
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

// Returns a report per cluster Rep
// function is the originally assigned annotation
// excludes "hypothetical protein" annotation
process postulatedReport {
    debug true
    //publishDir "$params.outDir/reports/single", mode: 'copy'
    containerOptions { "--rm" }

    input:
    path(proteinDescriptionsfile)
    path(clusterReps)

    output:
    path("*_ps.json"), optional: true

    script:
    """
    #!/usr/bin/env python

import pandas as pd
import re, json

def postulatedReport(descriptions, reps):
    # get cluster reps
    repList = []
    for line in open(reps, "r"):
        if line.startswith(">"):
            repList.append(line[1:].strip())
    
    # Read descriptions file
    descriptionData = pd.read_csv(descriptions, delimiter="\t", names=["id", "desc"], index_col=False)
    
    # Extract protein function from description
    f = lambda x: re.search(r"protein=([^\\]].[^\\]]*)", x).group(1)
    descriptionData["desc"] = descriptionData["desc"].apply(f)
    # Filter out non-reps
    descriptionData = descriptionData[descriptionData["id"].isin(repList)]   
    
    # Get function from description
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

postulatedReport("$proteinDescriptionsfile", "$clusterReps")
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

process writeGff {
    publishDir "$params.outDir/reports/gff", mode: 'copy'
    debug true
    containerOptions { "--rm" }

    input:
    tuple val(id), path(gff), path(json)

    output:
    path("*.gff")

    script:
    """
    #!/usr/bin/env python
import pandas as pd

def safeGet(match):
    # Check for entries that are None
    # Try signature -> entry -> desc.
    # Try signature -> desc.
    # Try signature -> name
    # Else None

    try:
        domain = match.get("signature", {}).get("entry", {}).get("description")
    except:
        domain = None
    if domain is not None:
        return domain

    try:
        domain = match.get("signature", {}).get("description")
    except:
        domain = None
    if domain is not None:
        return domain

    try:
        domain = match.get("signature", {}).get("name")
    except:
        domain = None
    if domain is not None:
        return domain

    return None


def chooseFunction(input):
    with open(input, "r") as infile:
        df = pd.DataFrame(pd.read_json(infile))

    df["protein"] = df["protein"].str.extract(r"prot_(.*?)_\\d+_unknown")

    entryStrings = []

    for entries in df["predictions"]:
        entryString = ""

        for entry in entries:
            method = entry["method"]
            # Case: Has interproscan result -> append domains
            if method == "interproscan": 
                for match in entry["matches"]:
                   
                    # Check if match is using an accepted library
                   # libs = {"SUPERFAMILY", "PFAM", "PROSITE_PROFILES"}
                    #if match["signature"]["signatureLibraryRelease"]["library"] in libs:
                    
                    domain = safeGet(match)
                    if domain:
                        entryString += ";domain=" + domain
                        
            elif method != "postulated function":
                # Case foldseek
                if method == "foldseek":
                    function = ""
                    # if its an unacceptable result, fall back on other results
                    for i in [2, 0, 1]:  # = [best aln-score, e-value, prob] in foldseek report
                        function = entry["hits"][i]["function"]
                        if function not in ["Uncharacterized protein", "None"]:
                            break

                    entryString += ";foldseek=" + function

                # Case cluster
                if method == "clusterMember":
                    function = entry["function"]
                    rep = entry["clusterRep"]
                    entryString += ";cluster=" + function + ";clusterRep=" + rep
                
                # singleton cluster
                if method == "singletonCluster":
                    function = entry["function"]
                    rep = entry["clusterRep"]
                    entryString += ";cluster=singletonCluster"

        entryStrings.append(entryString)

    df["entryString"] = entryStrings

    return df


def writeGff(df, gfffile):
    # read original gff, and save headers
    headers = []
    cols = ["id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    with open(gfffile, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                headers.append(line.strip())
    headers.append("##predicted functions and domains via #github link#")

    gdf = pd.read_csv(gfffile, sep="\t", comment="#", names=cols)

    # go through gff file (gdf dataframe), search for protein id (from df) in the attributes field
    # if you find a match, meaning there is a potential function, append it to attributes string
    # and write the gff file.
    for index, row in gdf.iterrows():
        str = row["attributes"]
        match = df.loc[df["protein"].apply(lambda x: x in str), "entryString"]
        result = match.tolist()
        # ignore empty fields (putative funtion entries)
        if result and result[0] != "":
            gdf.at[index, "attributes"] += result[0]

    with open("${id}_features_predicted.gff", "w") as outfile:
        outfile.write("\\n".join(headers) + "\\n")  # write headers
        gdf.to_csv(outfile, sep="\t", index=False, header=False)


reports = "${json}"
gffFile = "${gff}"

writeGff(chooseFunction(reports), gffFile)
    """
}