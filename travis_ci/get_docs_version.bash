#!/bin/bash
# see if we should run at all
if [ -n "$TRAVIS_TAG" ]; then
  # this is a tag build
  export DOCS_VERSION="$TRAVIS_TAG"
elif [ -n "$TRAVIS_BRANCH" ]; then
  if [ "$TRAVIS_BRANCH" == "master" ]; then
    export DOCS_VERSION="stable"
  elif [ "$TRAVIS_BRANCH" == "develop" ]; then
    export DOCS_VERSION="latest"
  elif [ "$TRAVIS_BRANCH" == "ocean/develop" ]; then
    export DOCS_VERSION="latest ocean"
  elif [ "$TRAVIS_BRANCH" == "ocean/coastal" ]; then
    export DOCS_VERSION="latest coastal"
  elif [ "$TRAVIS_BRANCH" == "landice/develop" ]; then
    export DOCS_VERSION="latest land-ice"
  fi
fi
