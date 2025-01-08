FROM ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2

WORKDIR /cache

# download model weights
RUN python -m colabfold.download

# download db setup script
RUN wget https://raw.githubusercontent.com/sokrypton/ColabFold/refs/heads/main/setup_databases.sh
RUN mkdir cache/databases
RUN cp ./setup_databases.sh cache/databases
RUN mmseqs -h

RUN apt-get update