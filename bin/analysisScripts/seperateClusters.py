import csv

# Cluster tsv file from foldseek output.
input_file_path = 'res_cluster.tsv'

# Open the input file and read the data
with open(input_file_path, mode='r', newline='') as infile:
    reader = csv.reader(infile, delimiter='\t')
    rows = list(reader)

# Write the output with newlines between clusters
with open('cluster_seperate.tsv', mode='w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    
    # Track the previous representative
    previous_representative = None
    
    for row in rows:
        current_representative = row[0]  # Get the representative ID from first column
        
        # If the representative has changed, insert a newline
        if previous_representative and current_representative != previous_representative:
            outfile.write('\n')
        
        # Write the current row
        writer.writerow(row)
        
        # Update the previous representative to the current one
        previous_representative = current_representative
