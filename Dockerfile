FROM continuumio/miniconda3:23.5.2-0-alpine
LABEL version="2.1.2" \
      description="Docker image for Pangolin"

RUN conda install -n base -c conda-forge mamba
# Install git for pangolin
RUN apk update && \
    apk add git bash

COPY environment.yml /environment.yml
# Install the conda environment
RUN mamba env create --quiet -f /environment.yml && conda clean -a
# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH=/opt/conda/envs/pangolin/bin:$PATH

# Install Pangolin
COPY . /pangolin/
WORKDIR /pangolin/
RUN pip install . && rm -rf /root/.cache/pip
RUN pangolin --version &> /pangolin-version.txt

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name pangolin > /pangolin.yml
WORKDIR /tmp/
