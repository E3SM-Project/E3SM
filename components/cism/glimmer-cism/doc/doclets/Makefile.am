if DO_DOC_PDF
 pdf_doc = gc-branches.pdf
 dirs = figs
endif

doc_DATA = $(pdf_doc)
docdir = $(prefix)/doc
SUBDIRS = $(dirs)

EXTRA_DIST = gc-branches.tex

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc *.bbl *blg *.html 


.tex.dvi::
	$(LATEX) $<


%.pdf:	%.dvi
	$(DVIPDFT) -o $@ $<
