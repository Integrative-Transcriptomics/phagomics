FROM ubuntu:22.04



RUN apt-get update && apt-get -y install wget
RUN wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz \
&&  tar xvfz foldseek-linux-gpu.tar.gz \
&&  rm foldseek-linux-gpu.tar.gz 


# temporary SP
#RUN mkdir db
# RUN foldseek/bin/foldseek databases Alphafold/Swiss-Prot db/AFSP tmp

# doesn't work ?
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PATH="/foldseek/bin/:$PATH"
