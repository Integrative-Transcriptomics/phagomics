FROM ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2

WORKDIR /cache
# download model weights
RUN python -m colabfold.download

# Install NVIDIA GPU drivers
# RUN pip install --upgrade pip
# RUN pip install --upgrade "jax[cuda12_pip]<0.4.24"
#RUN pip install jax
#RUN pip install jaxlib

# Install NVIDIA Container Toolkit

RUN curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
    curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

RUN apt-get update
RUN apt-get install -y nvidia-container-toolkit



# download db setup script
# RUN   wget https://raw.githubusercontent.com/sokrypton/ColabFold/refs/heads/main/setup_databases.sh && \
# mkdir cache/databases && \
# cp    ./setup_databases.sh cache/databases
