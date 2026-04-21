FROM continuumio/miniconda3:latest

RUN conda create -n pymc-dev python=3.11 -y && \
    conda install -n pymc-dev -c bioconda -c conda-forge \
        umi_tools pyfastx pytensor numpy pandas scipy matplotlib compilers -y

ENV PATH="/opt/conda/envs/pymc-dev/bin:$PATH"
