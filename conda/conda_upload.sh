PKG_NAME=acme_diags
USER=UVCDAT

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`date +%Y.%m.%d`
conda build .
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l nightly $CONDA_BLD_PATH/$TRAVIS_OS_NAME/$PKG_NAME-$VERSION-0.tar.bz2 --force
