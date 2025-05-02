import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# Barplot showing known and unknown protein counts per phage 
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
                        # df.loc[df["id"] == id, "known"] += 1
                        df.at[df.index[df["id"] == id][0], "known"] += 1
                    else:
                        # df.loc[df["id"] == id, "unknown"] += 1
                        df.at[df.index[df["id"] == id][0], "unknown"] += 1
                # add new entry
                else:
                    if row.endswith("_known"):
                        newRow = pd.DataFrame({"id": [id], "known": [1], "unknown": [0]})
                        df = pd.concat([df, newRow], ignore_index=True)
                    else:
                        newRow = pd.DataFrame({"id": [id], "known": [0], "unknown": [1]})
                        df = pd.concat([df, newRow], ignore_index=True)


    # calculate average % of unknown proteins
    df['total'] = df['known'] + df['unknown']
    df['unknown_percent'] = (df['unknown'] / df['total']) * 100
    avgUk = df['unknown_percent'].mean()
    #print(avgUk)

    df.sort_values("id", inplace = True)

    plt.subplots(figsize=(16, 8))
    plt.bar(df['id'], df['known'], width=0.4, edgecolor='black', label='Known Proteins', color='#24b064')
    plt.bar(df['id'], df['unknown'], bottom=df['known'], width=0.4, edgecolor='black', label='Unknown Proteins', color='#F18046')

    plt.xlabel("Phage ID", fontsize=16)
    plt.ylabel("Count", fontsize=16)
    plt.xticks(df['id'], labels=range(1, 79), fontsize=13)
    plt.yticks(fontsize=13) 
    #plt.title("Known and Unknown Proteins per Phage", fontsize=10)
    plt.legend(fontsize=13)
    plt.tight_layout()
    #plt.show()


    ## Stacked plot

    half = len(df) // 2
    df1 = df.iloc[:half]
    df2 = df.iloc[half:]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), sharey=True)
    fonts = 18
    # Top half
    ax1.bar(df1['id'], df1['known'], width=0.4, edgecolor='black', label='Known Proteins', color='#24b064')
    ax1.bar(df1['id'], df1['unknown'], bottom=df1['known'], width=0.4, edgecolor='black', label='Unknown Proteins', color='#F18046')
    ax1.set_ylabel("Count", fontsize=fonts)
    ax1.set_xticks(df1['id'])
    ax1.set_xticklabels(range(1, half + 1), rotation=0, fontsize=fonts)
    ax1.legend(fontsize=fonts)
    ax1.tick_params(axis='y', labelsize=fonts)
    ax1.text(-0.05, 1.05, "A", transform=ax1.transAxes, fontsize=22, fontweight='bold', va='top', ha='left')
    #ax1.set_title("Phage ID (Phage 1–39)", fontsize=14)

    # Bottom half
    ax2.bar(df2['id'], df2['known'], width=0.4, edgecolor='black', label='Known Proteins', color='#24b064')
    ax2.bar(df2['id'], df2['unknown'], bottom=df2['known'], width=0.4, edgecolor='black', label='Unknown Proteins', color='#F18046')
    ax2.set_xlabel("Phage ID", fontsize=fonts)
    ax2.set_ylabel("Count", fontsize=fonts)
    ax2.set_xticks(df2['id'])
    ax2.set_xticklabels(range(half + 1, len(df) + 1), rotation=0, fontsize=fonts)
    ax2.tick_params(axis='y', labelsize=fonts)
    ax2.text(-0.05, 1.05, "B", transform=ax2.transAxes, fontsize=22, fontweight='bold', va='top', ha='left')
    #ax2.set_title("Phage ID (Phage 40–78)", fontsize=14)

    plt.tight_layout()
    plt.savefig("figures/phagesProteinCountsStack.svg", dpi=400, bbox_inches='tight')
    # plt.show()


def processClusterFile(input):
    infile = pd.read_csv(input, sep="\t", header=None, names=["cluster_rep", "member"])
    # group by cluster rep and collect members
    cluster_data = infile.groupby("cluster_rep")["member"].apply(list).reset_index()
    
    # get total member count, and counta of known and unknown proteins
    cluster_data["total_members"] = cluster_data["member"].apply(len)
    cluster_data["known_count"]   = cluster_data["member"].apply(lambda members: sum("_known" in m for m in members))
    cluster_data["unknown_count"] = cluster_data["member"].apply(lambda members: sum("_unknown" in m for m in members))

    return cluster_data

# Barplot showing cluster sizes and their makeup
def plot_cluster_data(clusterData):

    cluster_size_bins = np.sort(clusterData["total_members"].unique()) # get the cluster sizes
    known_counts = []
    unknown_counts = []

    # count known/unknown per size
    for size in cluster_size_bins:
        known_counts.append(clusterData[clusterData["total_members"] == size]["known_count"].sum())
        unknown_counts.append(clusterData[clusterData["total_members"] == size]["unknown_count"].sum())

    fonts=16
    plt.figure(figsize=(10, 6))
    plt.bar(cluster_size_bins, unknown_counts, bottom=known_counts, label="Unknown", color="#F18046", edgecolor="black")
    bars_known = plt.bar(cluster_size_bins, known_counts, label="Known", color="#24b064", edgecolor="black")

    # numbers above bars known/unknown
    for bar, known, unknown in zip(bars_known, known_counts, unknown_counts):
        plt.text(bar.get_x() + bar.get_width()/2, known+unknown + 1, f"{known}/{unknown}", 
                ha='center', va='bottom', fontsize=fonts, color='black')
        
    plt.xlabel("Cluster Size", fontsize=fonts)
    plt.ylabel("Frequency", fontsize=fonts)
    plt.tick_params(axis='y', labelsize=fonts-2)
    plt.tick_params(axis='x', labelsize=fonts-2)
    #plt.title("Total Members by Cluster Size (Functional vs Hypothetical)", fontsize=20)
    plt.legend(fontsize=fonts)
    plt.grid(axis="y", alpha=0.7)
    plt.tight_layout()
    plt.savefig("figures/clusterMakeup.svg", dpi=400, bbox_inches='tight')
    #plt.show()

# Scatterplot showing heterogeity of cluster (proportion of known members) per cluster
def clusterHet(clusterData):
    df = clusterData[clusterData["total_members"] > 1]
    df["known_proportion"] = df["known_count"] / df["total_members"]

    plt.figure(figsize=(10, 6))
    plt.scatter(df["cluster_rep"], df["known_proportion"], 
                s=df["total_members"] * 5, alpha=0.7, 
                c=df["known_proportion"], cmap="Spectral", edgecolors="black")
    plt.colorbar(label="Proportion of Known Members")
    plt.xticks([])
    plt.xlabel("Cluster")
    plt.ylabel("Proportion of Known Members")
    plt.title("Heterogeneity per Cluster")
    plt.savefig("figures/clusterHet.svg", dpi=400, bbox_inches='tight')
    #plt.show()
    
cluData = processClusterFile("results/protein_clusters/clu2_cluster.tsv")
#cluData = processClusterFile("")  # Structural clustering
proteinRatio("results/proteins/allProteins.faa")
plot_cluster_data(cluData)
clusterHet(cluData)