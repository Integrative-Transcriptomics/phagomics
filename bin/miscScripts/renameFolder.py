import os

def extractID(fnaPath):
    with open(fnaPath, "r") as f:
        firstLine = f.readline().strip()
        if firstLine.startswith(">"):
            return firstLine.split()[0][1:]  # Extract ID after '>' and before space
    return None

def renameFolders(parentDir):
    for folderName in os.listdir(parentDir):
        folderPath = os.path.join(parentDir, folderName)
        if os.path.isdir(folderPath):
            fnaFiles = [f for f in os.listdir(folderPath) if f.endswith(".fna")]
            if fnaFiles:
                fnaPath = os.path.join(folderPath, fnaFiles[0])
                fileId = extractID(fnaPath)
                if fileId:
                    newFolderPath = os.path.join(parentDir, fileId)
                    if folderPath != newFolderPath:
                        os.rename(folderPath, newFolderPath)
                        print(f"Folder {folderName} renamed to: {fileId}")
                    else:
                        print(f"Folder {folderName} name is already correct.")
                else:
                    print(f"Could not extract ID from {fnaFiles[0]}.")


parentDir = input("Enter the path to the parent directory: ")
renameFolders(parentDir)