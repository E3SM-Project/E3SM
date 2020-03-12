#!/bin/bash

set -e

source travis_ci/get_docs_version.bash

if [[ -z ${DOCS_VERSION} ]]; then
  exit 0
fi

conda create -y -n test python=3.7 sphinx mock sphinx_rtd_theme m2r
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate test


cd docs || exit 1
make html
cd ..
