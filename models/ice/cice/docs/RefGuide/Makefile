#=======================================================================
# CVS: $Id: Makefile 5 2005-12-12 17:41:05Z mvr $
# CVS: $Source$
# CVS: $Name$
#=======================================================================
#--------------------------------------------------------
# This script will create dvi, ps, pdf and html
# versions of the CSIM4 Code Reference from the LaTeX files
#--------------------------------------------------------

#--------------------------------------------------------
# Info for creating the protex documentation
#--------------------------------------------------------
# JLS PRO_DOC := ice_protex
# JLS SRC_PATH := ../../src/source
# JLS SRC_CODE  := $(wildcard $(addsuffix /*.F , $(SRC_PATH)))

PRO_DOC := ice_protex
SRC_PATH := ../../src
SRC_CODE  := $(wildcard $(addsuffix */*.F90 , $(SRC_PATH)))

#--------------------------------------------------------
# Info for creating the code reference documentation
#--------------------------------------------------------
DOC   := ice_refdoc
VPATH := .
TEX   := $(wildcard $(addsuffix /*.tex , $(VPATH)))

all: $(DOC).dvi $(DOC).pdf $(DOC).ps $(DOC)/$(DOC).html
#all: $(PRO_DOC).tex $(DOC).dvi $(DOC).pdf $(DOC).ps $(DOC)/$(DOC).html

dvi: $(DOC).dvi

pdf: $(DOC).pdf

html: $(DOC)/$(DOC).html

ps: $(DOC).ps

# JLS protex: $(PRO_DOC).tex

#--------------------------------------------------------
#  (1) Build a .tex file from ProTeX in source code, then
#      build a .dvi file from LaTeX files 
#  NOTE:  Use the -b option for protex so that the pre-
#         amble is not written to the .tex file.  If the
#         protex stuff is the only .tex file, then don't
#         use the -b option.
#--------------------------------------------------------

# JLS $(PRO_DOC).tex: $(SRC_CODE) 
# JLS add tab here protex -b $(SRC_CODE) > $(PRO_DOC).tex

$(PRO_DOC).tex: $(SRC_CODE) 
	./protex -b $(SRC_CODE) > $(PRO_DOC).tex

$(DOC).dvi: $(TEX)
	latex $(DOC).tex
	bibtex $(DOC)
	latex $(DOC).tex
	latex $(DOC).tex

#--------------------------------------------------------
#  (2) Make .ps file for Code Reference
#--------------------------------------------------------

$(DOC).ps: $(DOC).dvi
	dvips -o $(DOC).ps $(DOC).dvi

#--------------------------------------------------------
#  (3) Make .pdf file for Code Reference
#--------------------------------------------------------

$(DOC).pdf: $(DOC).dvi
	pdflatex $(DOC).tex

#--------------------------------------------------------
#  (4) Produce an .html version
#      The html files will be in directory $(DOC)
#--------------------------------------------------------

$(DOC)/$(DOC).html: $(TEX)
	latex2html -white \
	           -toc_depth 5 \
	           -split 5 \
	           -show_section_numbers \
                   -address "<A HREF="mailto:csm@ucar.edu">csm@ucar.edu</A>" \
                   $(DOC).tex

clean:
	rm -f $(DOC).aux $(DOC).bbl $(DOC).blg $(DOC).dvi $(DOC).log $(DOC).toc \
              ice_protex.tex

realclean: clean
	rm -rf $(DOC) $(DOC).pdf $(DOC).ps ice_protex.tex

test:
	@echo "TEX := $(TEX)"


