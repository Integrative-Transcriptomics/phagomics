#!/usr/bin/env nextflow

process interproscan {
    //publishDir "$params.outDir/reports/single", mode: 'copy'
    maxForks 1 // make only 1 API request at a time
    errorStrategy 'retry'
    containerOptions { "--rm" }

    input:
    val(unknownProteins)

    output:
    path("*_ip.json"), optional: true

    script:
    """
#!/usr/bin/env python
import requests, json, time

# Interproscan API call
# Input is 30 fasta sequences
def interproscanAPI(input):
    REQUEST_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/"
    EMAIL = "${params.email}"

    def getPayload(input):
        seqs = {
            "email": EMAIL, 
            "sequence": input, 
            "goterms": True, 
            "pathways": True,
            "stype": "p"
        }
        return seqs
    
    def makeRequest(payload):
        headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
        'Accept': 'text/plain',
        }

        response = requests.post(f"{REQUEST_URL}run", headers=headers, data=payload)

        if response.status_code == 200:
            print(f"Job submitted successfully! Job ID: {response.text}")
            return response.text
        else:
            raise Exception("Job submission failed!")
        
    def checkStatus(jobID):
        maxTime = 60*5 # 5 mins
        timeout = time.time() + maxTime

        while time.time() < timeout:
            print("Task running...")
            response = requests.get(f"{REQUEST_URL}status/{jobID}", headers = {'Accept': 'text/plain'})
            time.sleep(10)

            if not response.text == "RUNNING":
                if response.text == "FINISHED":
                    print("Success!")
                    return True
                elif response.text == "FAILED":
                    raise Exception("Job failed!")


    def getResults(jobId):
        response = requests.get(f"{REQUEST_URL}result/{jobId}/json", headers = {'Accept': 'application/json'})
        
        if response.status_code == 200:
            print("Added result!")
            return response.text
        else:
            raise Exception("Job result failed!")
        
    # Run
    payload = getPayload(input)
    jobId   = makeRequest(payload)
    while not checkStatus(jobId):
        pass
    result = getResults(jobId)
    if result:
        return result
    
    return None


def filterMatches(contents):
    libs = {"SUPERFAMILY", "PFAM", "PROSITE_PROFILES"}
    for result in contents["results"]:
        result["matches"] = [
            match for match in result["matches"]
            if match["signature"]["signatureLibraryRelease"]["library"] in libs
        ]    

# input sets of 30 results
def interproscanReport(input):
    contents = json.loads(input)
    filterMatches(contents)
    for result in contents["results"]:
        if result["matches"]:
            id = result["xref"][0]["name"]
            with open(f"{id}_ip.json", "w") as fh:
                    result = {"method": "interproscan", **result}
                    fh.write(json.dumps(result, indent=4))

interproscanReport(interproscanAPI(\"""$unknownProteins\"""))
    """
}