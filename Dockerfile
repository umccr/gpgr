FROM ubuntu:20.04
LABEL maintainer="https://github.com/pdiakumis"

ARG MINI_VERSION=4.12.0-0
ARG MINI_URL=https://github.com/conda-forge/miniforge/releases/download/${MINI_VERSION}/Mambaforge-${MINI_VERSION}-Linux-x86_64.sh
ARG MAMBA_PREFIX="/opt/mambaforge"

# install core pkgs, mambaforge
RUN apt-get update -qq && \
    apt-get install --yes --no-install-recommends -qq \
    git bash bzip2 curl ca-certificates && \
    apt-get clean && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    curl --silent -L "${MINI_URL}" -o "mambaforge.sh" && \
    /bin/bash mambaforge.sh -b -p "${MAMBA_PREFIX}" && \
    rm mambaforge.sh

# create conda env
ENV PATH="${MAMBA_PREFIX}/bin:$PATH"
ARG ENV_DIR=/home/gpgr_conda_env
ARG ENV_NAME="gpgr_env"
COPY ./conda/env/lock ${ENV_DIR}
RUN conda config --set always_yes yes && \
    mamba install -c conda-forge conda-lock==1.0.5
RUN conda-lock install --name ${ENV_NAME} ${ENV_DIR}/conda-lock.yml && \
    mamba clean --all --force-pkgs-dirs

ENV PATH="${MAMBA_PREFIX}/envs/${ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_PREFIX}/envs/${ENV_NAME}"

CMD [ "gpgr.R" ]
