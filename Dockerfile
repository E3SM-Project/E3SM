# Always use the latest Docker image of Miniconda,
# since we have to update it anyways.
# This is because there are/were many issues with older versions
# that prevented building from working properly.
FROM continuumio/miniconda:latest

LABEL maintainer="shaheen2@llnl.gov"

COPY conda/e3sm_diags_env.yml e3sm_diags_env.yml

# Set bash as the default shell for all RUN commands.
SHELL ["/bin/bash", "-c"]

# AVOID CREATING ANOTHER ANACONDA ENV IN THE CONTAINER!
# So we don't have:
#   RUN conda env create -f e3sm_diags_env.yml
# This caused many issues.
# We just update the root environment.
# We also need to update conda, to make
# sure that everything works.
RUN conda update -n base conda && \
        conda env update -n root --file e3sm_diags_env.yml && \
        conda clean --all -y

# Needs to be a list, otherwise arguments aren't 
# passed correctly when you run the container.
ENTRYPOINT ["e3sm_diags"]

# Where e3sm_diags will read the input files,
# and where it'll store the output *relative to this container*.
# When running the container, mount WORKDIR to a local directory
# on the machine running the container.
WORKDIR /e3sm_diags_container_io
