SCREAM Docker
======

SCREAM Docker is built on top of CIME's Docker image. Make sure CIME image is built first, please refer to the CIME's [documentation](../../../cime/docker/README.md) for more information on building the base image.

## Building the container

From the top level of the SCREAM repository, run the following commands:
```bash
docker build -t cime:latest --target base cime/docker/
docker build --file components/eamxx/docker/Dockerfile -t scream:latest .
```

## Running the container

```bash
docker run -it --hostname docker -e CIME_MODEL=e3sm scream:latest bash
```

It is highly recommended to mount your current working scream repository to the container:

```bash
docker run -it -v $SCREAM_REPO:/src/E3SM --hostname docker -e CIME_MODEL=e3sm scream:latest bash
``` 
