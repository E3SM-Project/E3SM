#!/bin/bash

# Replace all ssh URLs to submodules with HTTP URLs
sed -i 's/git@github.com:/https:\/\/github.com\//' .gitmodules
git submodule update --init --recursive
