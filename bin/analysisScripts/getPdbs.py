import os, shutil

###
### Collects all pds matching the ids in fastaFile from pdbFolder and copies them to outputFolder.
### Attaches the potential function from foldseek to the name, to make it visible when clustering with foldseek
###

fastaFile = "foldseekHits.fasta"
pdbFolder = "resultsFull/colabfold"
outputFolder = "out/fsHitsFull"

os.makedirs(outputFolder, exist_ok=True)

# Extract sequence IDs and descriptions from FASTA file
seqInfo = {} 
with open(fastaFile, "r") as fasta:
    for line in fasta:
        if line.startswith(">"):  
            parts = line[1:].strip().split(maxsplit=1)
            seqId = parts[0]  # First part is ID
            seqDesc = parts[1] if len(parts) > 1 else "unknown"  
            seqDesc = seqDesc.replace(" ", "_").replace("/", "-")  # Replace spaces and slashes
            seqInfo[seqId] = seqDesc

# process pdb files
for pdbFile in os.listdir(pdbFolder):

    pdbId, ext = os.path.splitext(pdbFile)
    
    # ignore non-pdb files
    if ext != ".pdb":  
        continue 

    # Find matching sequence ID
    matchingId = next((seqId for seqId in seqInfo if pdbId.startswith(seqId)), None)
    
    if matchingId:
        # name = phage_protein_seqDescription.pdb
        newName = f"{matchingId}_{seqInfo[matchingId]}.pdb"
        
        # Remove the last 49 characters (extra colabfold info) before renaming
        newName = pdbId[:-49] + "_" + seqInfo[matchingId] + ".pdb"

        oldPath = os.path.join(pdbFolder, pdbFile)
        newPath = os.path.join(outputFolder, newName)

        shutil.copy(oldPath, newPath)

print(f"Copied pdb files to {outputFolder}")