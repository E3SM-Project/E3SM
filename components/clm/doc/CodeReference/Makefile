#
# Makefile to convert clm code documentation into html
# Currently using ProTex, eventually will be converted into
# Doxygen
#

# Set up special characters
null  :=
space := $(null) $(null)
comma := $(null),$(null)
dirs := . $(shell cat Filepath)

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir)))
VPATH    := $(subst $(space),:,$(VPATH))

# Get list of files
FIND_FILES = $(wildcard $(dir)/*.[Ff]90)
SOURCES    = $(foreach dir, $(dirs),$(FIND_FILES))

.SUFFIXES: .F90 .tex .html .dvi .pdf

HTMLCR  := index.html

all: $(HTMLCR)

clmcr.tex: $(SOURCES)
	protex -F $(SOURCES) > clmcr.tex

$(HTMLCR): clmcr.tex
	latex2html -noshow_section_numbers -nosubdir clmcr.tex

clean:
	rm -f $(HTMLCR)
realclean: clean
	rm -f clmcr.tex node*.html clmcr.css clmcr.html clmcr.aux clmcr.tex clmcr.pdf clmcr.toc internals.pl labels.pl WARNINGS
