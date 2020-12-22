FROM continuumio/miniconda3:4.9.2-alpine
LABEL authors="" \
      version="2.1.2" \
      description="Docker image for Pangolin"

# Install git for pangolin
RUN apk update && \
    apk add git bash

# Install the conda environment
COPY environment.yml /
COPY . /pangolin/
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/pangolin/bin:$PATH
WORKDIR /pangolin/
RUN pip install -e .
RUN pangolin --version &> pangolin-version.txt

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name pangolin > pangolin.yml
