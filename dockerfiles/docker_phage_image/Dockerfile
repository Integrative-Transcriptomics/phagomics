FROM python:3.12-slim

RUN apt-get update && apt-get install -y wget bedtools

# install foldseek
RUN wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz \
&&  tar xvfz foldseek-linux-gpu.tar.gz \
&&  rm foldseek-linux-gpu.tar.gz 

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PATH="/foldseek/bin/:$PATH"
