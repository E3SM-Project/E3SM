## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005-2006, see the COPYRIGHT file for more information.
# This file builds the netcdf documentation.

# The "all" target does nothing in this directory. The "dist" target
# *does* build the documentation, because the distribution contains
# all of the built documentation.

# The "install" target installs the documentation in a directory
# called $(prefix)/doc/$(PACKAGE)-$(VERSION).

man_MANS = netcdf_fortran.3

# These are the source files for all the netcdf manuals.
# info_TEXINFOS = netcdf-f90.texi netcdf-f77.texi

# Get lists of the coresponding ps, info, and pdf files.
#html_docs = ${info_TEXINFOS:.texi=.html}
#html_dirs = ${info_TEXINFOS:.texi=}
#ps_docs = ${info_TEXINFOS:.texi=.ps}
#dvi_docs = ${info_TEXINFOS:.texi=.dvi}
#pdf_docs = ${info_TEXINFOS:.texi=.pdf}
#info_docs= ${info_TEXINFOS:.texi=.info}
#txt_docs= ${info_TEXINFOS:.texi=.txt}

# These are html versions of the man pages.
#html_mans = netcdf_fortran-man.html

# These files will be included with the dist.
EXTRA_DIST = netcdf.m4 $(man_MANS) \
CMakeLists.txt netcdf-50x50.png Doxyfile.in \
cfortran.doc

# All the documents depend on defines.texi for some common
# definitions.
# netcdf_fortran_TEXINFOS = defines.texi

# Get texi2dvi to shut the hell up, and work better for parallel
# builds. The --tidy option may be useful getting parallel builds to
# work someday.
# TEXI2DVI = texi2dvi -s

# Turn off parallel builds in this directory.
# .NOTPARALLEL:

# To build documentation you must have M4 in their path.
M4=m4

# Rule to generate the text copies of the documentation.
#.texi.txt:
#	makeinfo -o $*.txt --no-headers $(srcdir)/$<

# Generate the HTML version of the netcdf-3 F03 man page.
#netcdf_fortran-man.html: $(top_srcdir)/docs/netcdf_fortran.3
#	groff -man -Thtml $? > $@

#netcdf-fortran.html:
#	$(MAKEINFO) --html $(AM_MAKEINFOHTMLFLAGS) $(MAKEINFOFLAGS) -I$(srcdir) netcdf-fortran.texi
#	$(MAKEINFO) --no-split --html $(AM_MAKEINFOHTMLFLAGS) $(MAKEINFOFLAGS) -I$(srcdir) netcdf-fortran.texi

# This builds the docs from source, if necessary, and tars up
# everything needed for the website. Run this and copy the resulting
# tarball to the /contents/netcdf/docs directory to update the on-line
# docs.

#web-tarball: $(pdf_docs) $(html_mans) $(txt_docs) $(dvi_docs) $(ps_docs) $(info_docs) $(html_docs)
#	tar cf netcdf_info.tar *.info*
#	gzip -f netcdf_info.tar
#	zip netcdf_info *.info
#	tar cf netcdf_pdf.tar *.pdf
#	gzip -f netcdf_pdf.tar
#	zip netcdf_pdf *.pdf
#	tar cf netcdf_html.tar *.html $(html_dirs)
#	gzip -f netcdf_html.tar
#	zip netcdf_html *.html $(html_dirs)
#	tar cf netcdf_txt.tar *.txt
#	gzip -f netcdf_txt.tar
#	zip netcdf_txt *.txt
#	tar cf netcdf_man.tar netcdf_fortran.3
#	gzip -f netcdf_man.tar
#	zip netcdf_man netcdf_fortran.3
#	tar cf netcdf_docs.tar $(html_dirs) *.pdf *.txt *.html *.tar.gz *.zip
#	gzip -f netcdf_docs.tar

# Clean up all the junk we created!
MAINTAINERCLEANFILES = netcdf_fortran-man.html stamp-* version*

ARGS_MANPAGE = -DAPI=FORTRAN -DNETCDF4=TRUE
if BUILD_DAP
ARGS_MANPAGE += -DDAP=TRUE
endif
if BUILD_PARALLEL
ARGS_MANPAGE += -DPARALLEL_IO=TRUE
endif
manpages: $(top_srcdir)/docs/netcdf.m4
	m4 $(M4FLAGS) $(ARGS_MANPAGE) $? >netcdf_fortran.3 || rm netcdf_fortran.3
