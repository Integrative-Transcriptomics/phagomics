import csv

file1_path = 'results/protein_clusters/clu2_cluster.tsv'  # File 1 contains the clusters
file2_path = 'results/proteins/proteinDescriptions.tsv'  # File 2 contains protein descriptions

# Open File 2 and create a dictionary mapping ID -> description
id_to_value = {}
with open(file2_path, mode='r', newline='') as file2:
    reader = csv.reader(file2, delimiter='\t')
    for row in reader:
        id_to_value[row[0]] = row[1]

# Open File 1, read its contents, and write the updated version
with open(file1_path, mode='r', newline='') as file1:
    reader = csv.reader(file1, delimiter='\t')
    rows = list(reader)

# Add the new column based on matching IDs
for row in rows:
    id_in_file1 = row[1]  # IDs are in the second column of File 1
    value_to_add = id_to_value.get(id_in_file1, 'None')  # Default to None if no match
    row.append(value_to_add)  # Append the value to the row

# Write the updated rows to a new file
with open('seqClustersDesc.tsv', mode='w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerows(rows)
