FROM mambaorg/micromamba:latest

# Create env and install dependencies
RUN micromamba install -y \
    -n base \
    -c conda-forge \
    -c bioconda \
    pandas=1.4.3 bedtools=2.30.0 && \
    micromamba clean --all --yes

# Activate the environment to run mamba commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Set entrypoint in active environment
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
