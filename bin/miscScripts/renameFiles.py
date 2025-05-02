import os

def renameFiles(parent_directory):
    for folderName in os.listdir(parent_directory):
        folderPath = os.path.join(parent_directory, folderName)
        if os.path.isdir(folderPath):
            for fileName in os.listdir(folderPath):
                filePath = os.path.join(folderPath, fileName)
                if fileName.endswith(".fna"):
                    newName = f"{folderName}_genome.fna"
                elif fileName.endswith(".gff"):
                    newName = f"{folderName}_features.gff"
                elif fileName.endswith(".faa"):
                    newName = f"{folderName}_proteins.faa"
                else:
                    continue
                
                new_path = os.path.join(folderPath, newName)
                os.rename(filePath, new_path)
                print(f"Renamed {fileName} to {newName}")


parentDir = input("Enter the path to the parent directory: ")
renameFiles(parentDir)
