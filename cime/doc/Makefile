# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = on
SOURCEDIR     = source
BUILDDIR      = build
SPHINXAPI     = sphinx-apidoc
SCRIPTSDIR    = ../scripts

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# exit 0 is to prevent the api rule from being run by the catchall target
# and generating an error. 
api: CIME_api Tools_api Tools_user 
	exit 0

CIME_api: 
	@$(SPHINXAPI) --force -o $(SOURCEDIR)/$@ $(SCRIPTSDIR)/lib/CIME

Tools_api: 
	@$(SPHINXAPI) --force -o $(SOURCEDIR)/$@ $(SCRIPTSDIR)/Tools

Tools_user:
	rm -f $(SOURCEDIR)/$@/*.rst
	rm -f $(SOURCEDIR)/$@/temp_files/*
	./tools_autodoc.py

# It's too easy to forget to run 'make api' before running 'make html',
# so add a rule that ensures that the api documentation is regenerated
# whenever regenerating the html.
html: api

clean: clean_api

# Note that all of the files removed here are built using 'make api';
# these are not - or at least, should not be - files that exist in the
# repository.
clean_api:
	rm -f $(SOURCEDIR)/CIME_api/*.rst
	rm -f $(SOURCEDIR)/Tools_api/*.rst
	rm -f $(SOURCEDIR)/Tools_user/*.rst
	rm -f $(SOURCEDIR)/Tools_user/temp_files/*

.PHONY: help Makefile CIME_api Tools_api Tools_user clean_api

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
