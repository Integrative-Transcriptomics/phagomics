FROM python:3.12-slim

RUN apt-get update && apt-get install -y wget && pip install requests pandas

# Install foldseek
RUN wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz \
&&  tar xvfz foldseek-linux-gpu.tar.gz \
&&  rm foldseek-linux-gpu.tar.gz 

# Add foldseek to path
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PATH="/foldseek/bin/:$PATH"

# # Install InterProScan
# RUN mkdir interproscan \
# && cd interproscan \
# && wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/interproscan-5.73-104.0-64-bit.tar.gz \
# && wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/interproscan-5.73-104.0-64-bit.tar.gz.md5 \
# && md5sum -c interproscan-5.73-104.0-64-bit.tar.gz.md5 \
# && tar -pxvzf interproscan-5.73-104.0-*-bit.tar.gz

# # index hmm models for interproscan
# RUN cd interproscan-5.73-104.0/ \ 
# && python3 setup.py -f interproscan.properties \
# && cd .