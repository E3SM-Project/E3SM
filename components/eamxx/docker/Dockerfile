ARG CIME_VERSION=latest
FROM cime:${CIME_VERSION} AS base

ENV USER=root
ENV LOGNAME=root
ENV DEBIAN_FRONTEND="noninteractive"

SHELL ["/bin/bash", "-c"]

# Install additional packages
RUN mamba install --yes -c conda-forge \
            lapack \            
            blas && \
            rm -rf /opt/conda/pkgs/*

# install gcloud
RUN apt-get update && apt-get install -y  apt-transport-https ca-certificates gnupg curl gettext

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list &&\
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

RUN apt-get update && apt-get install -y google-cloud-sdk
RUN gcloud config set project vcm-ml