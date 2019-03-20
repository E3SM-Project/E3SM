# Specify the version of the Docker image of Miniconda.
# This is so we create the same image each time.
# Also, there are/were many issues with older versions
# that prevented the building process from working properly.
FROM continuumio/miniconda:4.5.4

LABEL maintainer="shaheen2@llnl.gov"
LABEL version="1.6.1"

# Copy the entire project dir because we'll install from source.
COPY . .

# Set bash as the default shell for all RUN commands.
SHELL ["/bin/bash", "-c"]

# AVOID CREATING ANOTHER ANACONDA ENV IN THE CONTAINER!
# So we don't have:
#   RUN conda env create -f e3sm_diags_env_dev.yml
# This caused many issues.
# We just update the base environment.
RUN conda env update -n base --file conda/e3sm_diags_env_dev.yml && \
        conda clean --all -y && \
        source activate base && \
        python setup.py install && \
        rm -r build/

# Needs to be a list, otherwise arguments aren't
# passed correctly when you run the container.
ENTRYPOINT ["e3sm_diags"]

# Set an environmental variable so e3sm_diags
# can know if it's running in a container.
ENV E3SM_DIAGS_CONTAINER true

# When running the container, you might need to mount WORKDIR
# to the cwd on the machine running the container, so that
# the input py and cfg files can be read properly.
WORKDIR /e3sm_diags_container_cwd
