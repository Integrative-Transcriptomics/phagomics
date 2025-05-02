from pathlib import Path
from Bio import SeqIO
import numpy as np
import re
import pandas as pd

# Get protein Ratios (Known/Unknown)
def proteinRatio(allProteins):
    df = pd.DataFrame(columns=["id", "known", "unknown"])
    with open(allProteins, "r") as infile:
        for row in infile:
            # get id
            if row.startswith(">"):
                row = row.rstrip()
                id = re.search(r">(.+?)_prot", row).group(1)

                # increment count
                if id in df["id"].values:
                    if row.endswith("_known"):
                        df.at[df.index[df["id"] == id][0], "known"] += 1
                    else:
                        df.at[df.index[df["id"] == id][0], "unknown"] += 1
                # add new entry
                else:
                    if row.endswith("_known"):
                        newRow = pd.DataFrame({"id": [id], "known": [1], "unknown": [0]})
                        df = pd.concat([df, newRow], ignore_index=True)
                    else:
                        newRow = pd.DataFrame({"id": [id], "known": [0], "unknown": [1]})
                        df = pd.concat([df, newRow], ignore_index=True)

    df.sort_values("id", inplace = True)
    print("Known/unknown count done!\n")
    return df


###
### Input
###
resultFolder = input("Enter name of result folder (results by default): ")
counts = proteinRatio(f"{resultFolder}/proteins/allProteins.faa")  # Fasta of all proteins
reportPath = Path(f"{resultFolder}/reports/final")  # Folder of all final reports
reps = f"{resultFolder}/protein_clusters/clu2_rep_seq.fasta"  # cluster reps from mmseqs2
descFile = f"{resultFolder}/proteins/proteinDescriptions.tsv"  # protein descriptions
 

counts["hits"] = 0
phageData = []
phageDatafiltered = []

# import phage reports
for file in reportPath.glob("*.json"):
    phageData.append(pd.read_json(file))

# load protein descriptions
descDf = pd.read_csv(descFile, sep="\t")

# Remove entries with only postulated function or only singletonCluster or only both
# This is to only get true Hits
def shouldRemove(entry):
    methods = [item.get("method") for item in entry]
    
    if len(methods) == 1 and methods[0] in {"singletonCluster", "postulated function"}:
        return True
    if len(methods) == 2 and set(methods) == {"singletonCluster", "postulated function"}:
        return True
    
    return False

for i, df in enumerate(phageData):
    phageDatafiltered.append(df[~df["predictions"].apply(shouldRemove)])

###############

### Calculate number of hits for each phage
### Hit = foldseek/clustermember/interpro entry in report
### No. of hits = entries in "protein" column of df

for i, df in enumerate(phageDatafiltered):
    totalHits = len(df["protein"])
    id = df.loc[0, "protein"]
    id = re.search(r"(.+?)_prot", id).group(1).strip()  # get protein id

    counts.loc[counts["id"] == id, "hits"] = totalHits

sumUk = counts['unknown'].sum()
sumHits = counts['hits'].sum()


###############

###
### Collects all the unknown sequence reps that have NO found hits
### elif: collects unknown reps that DO have a hit
###

noHit = []
ukRepsHit = []
avgLen = []
knownReps = []
phageDataAll = pd.concat(phageDatafiltered, ignore_index=True)
ids = phageDataAll["protein"].to_list()  # Proteins that have a hit

for record in SeqIO.parse(reps, "fasta"):
    if record.id not in ids and not record.id.endswith("_known"):
        noHit.append((record.id, record.seq))
        avgLen.append(len(record.seq))
    elif record.id in ids and not record.id.endswith("_known"):
        ukRepsHit.append((record.id, record.seq))
    elif record.id.endswith("_known"):
        knownReps.append((record.id, record.seq))


###############

###
### This is for: found function (foldseek) -> structural clustering -> function of clusters
###

temp = []
# Get entries with foldseek hit (foldseek will always be first in json report)
for i, df in enumerate(phageData):
    fsHit = df[(df["predictions"].apply(lambda x: x[0].get("method") == "foldseek"))]   
    temp.append(fsHit)

fsHitsdf = pd.concat(temp, ignore_index=True)
ids = fsHitsdf["protein"].to_list()

# Write Foldseek hits to fasta and also the function to fasta ID
# Optionally also include known proteins
fsHits = []
tempKnown = []  # only used for this function
output = "foldseekHits.fasta"
with open(output, "w") as outfile:
    for record in SeqIO.parse(reps, "fasta"):
        function = ""
        # get unknown protein foldseek hits
        if record.id in ids:
            for i in [2, 0, 1]:  # = [best aln-score, e-value, prob] in foldseek report
                # Awful term to get the protein function for each record 
                function = fsHitsdf.loc[fsHitsdf["protein"] == record.id, "predictions"].tolist()[0][0].get("hits")[i].get("function")
                if function not in ["Uncharacterized protein", "None"]:
                    break

            record.description = function
            fsHits.append((record.id, record.seq))
            SeqIO.write(record, outfile, "fasta")
        
        # get known proteins
        elif record.id in descDf.iloc[:, 0].values and record.id.endswith("_known"):  # Check if record.id is in column 1 of descDf
            function = descDf.loc[descDf.iloc[:, 0] == record.id, descDf.columns[1]].values[0]  # Get the description string from column 2
            proteinMatch = re.search(r"\[protein=([^\]]+)\]", function)  # Extract the [protein=...] part
            if proteinMatch:
                record.description = proteinMatch.group(1)  # Set description to extracted protein name
            
            SeqIO.write(record, outfile, "fasta")
        
print(f"\nNo of Foldseek hits: {len(fsHits)}")

###############


###
### Validation report eval
###

try: 
    validFile = f"{resultFolder}/validationReport.tsv" 
    validdf = pd.read_csv(validFile, sep="\t")
    sumTrue = validdf["isQuery"].sum()
    validationStatus = f"Successfully predicted and identified {sumTrue} protein structures."
except (all): 
    sumTrue = 0
    validationStatus = "Internal validation skipped!"

###
### Generate number file:
###

# Calculations
numberReps = sum(1 for _ in SeqIO.parse(reps, "fasta"))
noUnknownReps = numberReps-len(knownReps)

# Writing
stats = []
stats.append(f"Total phages: {counts['id'].count()}")
stats.append(f"Total proteins: {counts['known'].sum()+counts['unknown'].sum()}")
stats.append(f"Total known: {counts['known'].sum()}")
stats.append(f"Total unknown: {counts['unknown'].sum()}")
stats.append(f"Total hits: {sumHits}")
stats.append(f"% hits of total unknown proteins: {round(sumHits/sumUk*100, 2)}%")
stats.append(f"")
stats.append(f"Total Sequence representatives/clusters: {numberReps}")
stats.append(f"# of known reps: {len(knownReps)}")
stats.append(f"# of unknown reps: {noUnknownReps}")
stats.append(f"# of hits for unknown reps: {len(ukRepsHit)}")
stats.append(f"# of unknown reps with no hits: {len(noHit)}")
stats.append(f"% hits of unknown reps: {round(len(ukRepsHit)/noUnknownReps*100, 2)}%")
stats.append(f"")
stats.append(f"Mean protein length of reps with no hit: {round(np.mean(avgLen), 2)} residues")
stats.append(f"")
stats.append(f"# of structural (foldseek) hits: {len(fsHits)}")
stats.append(f"% structural hits of unknown reps: {round(len(fsHits)/noUnknownReps*100, 2)}%")
stats.append(f"")
stats.append(f"{validationStatus}")
stats.append(f"{round(sumTrue/len(knownReps)*100, 2)}% of known cluster reps")

with open(f"{resultFolder}Stats.txt", "w") as f:
    for item in stats:
        f.write(f"{item}\n")
        
print(f"Stats saved to: {resultFolder}Stats.txt")