# Uploads the latest code from master as $VERSION on $OS as a main release
# usage:
# export OS=linux-64
# or 
# export OS=osx-64
# then
# export VERSION=v0.1a
PKG_NAME=acme_diags
USER=acme

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
conda build .
anaconda upload -u $USER $CONDA_BLD_PATH/$OS/$PKG_NAME-$VERSION-0.tar.bz2 --force
