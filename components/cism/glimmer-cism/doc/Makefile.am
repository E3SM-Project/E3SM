if DO_DOC_PDF
 pdf_doc = glimmer.pdf
 dirs = dg ug num tut common doclets ext
endif
if DO_DOC_PS
 ps_doc = glimmer.ps
 dirs = dg ug num tut common doclets ext
endif
doc_DATA = $(ps_doc) $(pdf_doc)
docdir = $(prefix)/doc
SUBDIRS = $(dirs)

EXTRA_DIST = dg.tex dg_pdf.tex dg_ps.tex glimmer.tex glimmer_pdf.tex \
             glimmer_ps.tex num.tex num_pdf.tex num_ps.tex ug.tex ug_pdf.tex \
             ug_ps.tex tut_pdf.tex tut_ps.tex tut.tex \
             ext.tex ext_ps.tex ext_pdf.tex glimmer.bib

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc *.bbl *blg *.html \
             glimmer_html.4ct  glimmer_html.4tc  glimmer_html.css  glimmer_html.idv  \
             glimmer_html.lg  glimmer_html.tmp  glimmer_html.xref

editdoxy = sed -e 's|@top_srcdir[@]|$(top_srcdir)|g' \
               -e 's|@GLIMMER_VERSION[@]|$(GLIMMER_VERSION)|g' \
               -e 's|@PACKAGE_NAME[@]|$(PACKAGE_NAME)|g'
include $(top_srcdir)/m4macros/amdoxygen.am

glimmer-cism.cfg:	$(srcdir)/glimmer-cism.cfg.in Makefile
			$(editdoxy) '$(srcdir)/$@.in' >$@

glimmer_html::	glimmer.tex
	install -d www
	htlatex glimmer_html "html,3,info" "" " -d./www/ -m 644 "
	$(BIBTEX) glimmer_html
	htlatex glimmer_html "html,3,info" "" " -d./www/ -m 644 "
	find $(SUBDIRS) -name '*.png' -exec install -m 644 -D {} ./www/{} \;



WWW_DIR = docs
www::	doxygen-run glimmer_html doxygen-run
	rm -rf $(WWW_DIR)
	install -d $(WWW_DIR)/$(GLIMMER_VERSION)/API
	install -d $(WWW_DIR)/$(GLIMMER_VERSION)/manual
	cp -a doxygen/html/* $(WWW_DIR)/$(GLIMMER_VERSION)/API
	cp -a www/* $(WWW_DIR)/$(GLIMMER_VERSION)/manual
	cd $(WWW_DIR) && ln -s $(GLIMMER_VERSION) current
	tar cvzf $(PACKAGE)-doc.$(GLIMMER_VERSION).tar.gz $(WWW_DIR)
	rm -rf $(WWW_DIR)

.tex.dvi::
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<

%.ps:	%_ps.dvi
	dvips -o $@ $<

%.pdf:	%_pdf.dvi
	$(DVIPDFT) -o $@ $<
