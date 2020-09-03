#!/bin/bash
# based on https://stackoverflow.com/a/49516361/7728169

set -e

source travis_ci/get_docs_version.bash

if [[ -z ${DOCS_VERSION} || "$TRAVIS_PULL_REQUEST" != "false" ]]; then
  exit 0
fi

# We only make it this far when a PR gets merged or when there is a new
# MPAS-Model tag or release
echo "Docs version: $DOCS_VERSION"

DOCS_PATH="${DOCS_VERSION// /_}"

PUBLICATION_BRANCH=gh-pages
# Checkout the branch
REPO_PATH=$PWD
pushd $HOME || exit 1
git clone --branch=$PUBLICATION_BRANCH https://${GITHUB_TOKEN}@github.com/$TRAVIS_REPO_SLUG publish
cd publish || exit 1

# Update pages
if [[ -d "$DOCS_PATH" ]]; then
  git rm -rf "$DOCS_PATH" > /dev/null
fi
mkdir -p "$DOCS_PATH"
cp -r "$REPO_PATH"/docs/_build/html/* "$DOCS_PATH"
# Commit and push latest version
git add .
if git diff-index --quiet HEAD; then
  echo "No changes in the docs."
else
  git config user.name  "Travis"
  git config user.email "travis@travis-ci.org"
  git commit -m "Updated $DOCS_VERSION"
  git push -fq origin $PUBLICATION_BRANCH
fi
popd || exit 1
