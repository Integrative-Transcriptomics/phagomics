import pandas as pd

# Load annotated cluster file
annotatedDf = pd.read_csv("cluster_with_descriptions.tsv", sep="\t", header=None, names=["cluster_rep", "member_with_desc"])

# map: cluster_rep => description
repDescMap = {}
for _, row in annotatedDf.iterrows():
    repId = row["cluster_rep"]
    member_with_desc = row["member_with_desc"]
    if repId not in repDescMap:
        parts = member_with_desc.split(" ", 1)
        if len(parts) == 2:
            repDescMap[repId] = parts[1]
        else:
            repDescMap[repId] = ""

# Load the second cluster file
secondClusterFile = "protein_clusters/clu1_cluster.tsv"
df2 = pd.read_csv(secondClusterFile, sep="\t", header=None, names=["cluster_rep", "member"])

# Annotate members and count successful annotations
descriptionCount = 0

def annotate_member(row):
    global descriptionCount
    repId = row["cluster_rep"]
    memberId = row["member"]
    desc = repDescMap.get(repId)
    if desc:
        descriptionCount += 1
        return f"{memberId} {desc}"
    else:
        return f"{memberId}"

df2["member_with_desc"] = df2.apply(annotate_member, axis=1)

# Output
output_file = "cluster2_with_descriptions.tsv"
df2[["cluster_rep", "member_with_desc"]].to_csv(output_file, sep="\t", index=False, header=False)

print(f"Added descriptions to {descriptionCount} proteins.")
