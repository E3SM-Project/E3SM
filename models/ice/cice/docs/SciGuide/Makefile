#=======================================================================
# CVS: $Id: Makefile 5 2005-12-12 17:41:05Z mvr $
# CVS: $Source$
# CVS: $Name$
#=======================================================================
#--------------------------------------------------------
# This script will create dvi, ps, pdf and html
# versions of the CSIM4 Scientific Description from the
# LaTeX files
#--------------------------------------------------------

DOC   := ice_scidoc
VPATH := .
TEX   := $(wildcard $(addsuffix /*.tex , $(VPATH)))

all: $(DOC).dvi $(DOC).pdf $(DOC).ps 

dvi: $(DOC).dvi

pdf: $(DOC).pdf

ps: $(DOC).ps

#--------------------------------------------------------
#  (1) Build a .dvi file from LaTeX files 
#--------------------------------------------------------

$(DOC).dvi: $(TEX)
	latex $(DOC).tex
	bibtex $(DOC)
	latex $(DOC).tex
	latex $(DOC).tex

#--------------------------------------------------------
#  (2) Make .ps file for User's Guide
#--------------------------------------------------------

$(DOC).ps: $(DOC).dvi
	dvips -o $(DOC).ps $(DOC).dvi

#--------------------------------------------------------
#  (3) Make .pdf file for User's Guide
#--------------------------------------------------------

$(DOC).pdf: $(DOC).dvi
	dvipdf $(DOC).dvi $(DOC).pdf

clean:
	rm -f $(DOC).aux $(DOC).bbl $(DOC).blg $(DOC).dvi $(DOC).log $(DOC).toc \
              $(DOC).lof

realclean: clean
	rm -rf $(DOC).pdf $(DOC).ps

test:
	@echo "TEX := $(TEX)"


