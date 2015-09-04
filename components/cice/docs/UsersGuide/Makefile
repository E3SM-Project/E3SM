#=======================================================================
# SVN: $Id: Makefile 5 2005-12-12 17:41:05Z mvr $
#=======================================================================
#--------------------------------------------------------
# This script will create dvi, ps, pdf and html
# versions of the CICE4 User's Guide from the LaTeX files
#--------------------------------------------------------

DOC   := ice_usrdoc
VPATH := .
TEX   := $(wildcard $(addsuffix /*.tex , $(VPATH)))

all: $(DOC).dvi $(DOC).pdf $(DOC).ps $(DOC)/$(DOC).html

dvi: $(DOC).dvi

pdf: $(DOC).pdf

html: $(DOC)/$(DOC).html

ps: $(DOC).ps

#--------------------------------------------------------
#  (1) Build a .dvi file from LaTeX files 
#--------------------------------------------------------

$(DOC).dvi: $(TEX)
	latex $(DOC).tex
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
	pdflatex $(DOC).tex

#--------------------------------------------------------
#  (4) Produce an .html version
#      The html files will be in directory $(DOC)
#--------------------------------------------------------

$(DOC)/$(DOC).html: $(TEX)
	latex2html -white \
	           -toc_depth 7 \
	           -split 5 \
	           -show_section_numbers \
                   -address "<A HREF="mailto:csm@ucar.edu">csm@ucar.edu</A>" \
                   $(DOC).tex

clean:
	rm -f $(DOC).aux $(DOC).bbl $(DOC).blg $(DOC).dvi $(DOC).log $(DOC).toc

realclean: clean
	rm -rf $(DOC) $(DOC).pdf $(DOC).ps

test:
	@echo "TEX := $(TEX)"


