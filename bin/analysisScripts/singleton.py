import pandas as pd

# Load the TSV file
file_path = "cluster_seperate.tsv"
df = pd.read_csv(file_path, sep='\t', header=None, names=["Representative", "Member"])

# Count number of members per representative to identify singleton clusters
cluster_sizes = df.groupby("Representative").size().reset_index(name="ClusterSize")

# Merge back into original dataframe to annotate singleton clusters
df = df.merge(cluster_sizes, on="Representative")

# Filter for singleton clusters only
singleton_clusters = df[df["ClusterSize"] == 1]

# Parse identifiers: phageID_proteinID_status_foundFunction
def parse_entry(entry):
    try:
        parts = entry.split("_")
        # phageID starts with either NC_ or MZ...
        if parts[0] == "NC":
            parts = [parts[0]+parts[1], parts[3]+parts[4]+"_"+parts[5], parts[6], "_".join(parts[7:])]
        else:
            parts = [parts[0], parts[2]+parts[3], parts[4], "_".join(parts[5:])]
        phage_id = parts[0]
        protein_id = parts[1]
        status = parts[2]
        found_function = parts[3]
        return pd.Series([phage_id, protein_id, status, found_function])
    except IndexError:
        return pd.Series([None, None, None, None])

# Apply parsing to each singleton
singleton_clusters[["PhageID", "ProteinID", "Status", "Function"]] = singleton_clusters["Representative"].apply(parse_entry)

# Filter for unknown proteins only
unknown_singletons = singleton_clusters[singleton_clusters["Status"] == "unknown"]

# Summary statistics
unknown_summary = unknown_singletons["Function"].value_counts().head(15)
phage_distribution = unknown_singletons["PhageID"].value_counts()

# Display output
print(unknown_singletons.head(), unknown_summary, phage_distribution)
unknown_singletons.to_csv(sep="\t", path_or_buf="singletons.tsv")
print(phage_distribution.to_string())