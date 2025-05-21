# Analysis Workflow using analysisScripts

This file explains what the scripts in bin/analysisScipts do, as well as the general analysis workflow used in the thesis.

0. Run the pipeline
1. analysis.py
    - This scripts generates run statistics. It promts for a results folder, which by default should be "results".
    - Generates "foldseekHits.fasta", which includes the proteins with putative functions assinged by Foldseek search.
2. getPdbs.py
    - Gets the .pdb files (protein structure files) from the colabfold subfolder, for structures that had a putative function assigned through Foldseek. These can then be clustered using Foldseek to identify similar predicted structures.
    - Files are saved to "out/fsHits"
3. Manual structural clustering
    - run foldseek easy-cluster out/fsHits res tmp
    - This generates files:
        - res_all_seqs.fasta: All protein sequences from fsHits
        - res_cluster.tsv: The structural clusters
        - res_rep_seq.fasta: The representative sequences
4. Cluster files (sequence and structure) were analysed using:
    - annotateClusters.py: Only sequence clusters. This script attaches the original protein descriptions to the cluster file, to analyse protein sequences that were clustered by sequence identity. Produces "seqClustersDesc.tsv"
    - seperateClusters.py: Used on sequence and structural cluster files. This script separates the clusters in the tsv file with newlines for easier reading. Produces "cluster_seperate.tsv"
    - singleton.py: Takes the separated structural cluster file and extracts singleton clusters, and saves them to "singletons.tsv". Also prints found functions with same names, and the amount of singleton clusters per phageID.