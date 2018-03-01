SUBDIRS=IO
pypkglibdir = @pythondir@/PyGC
pypkglib_PYTHON = __init__.py GC_options.py GC_colourmap.py GC_autoconf.py
EXTRA_DIST = GC_autoconf.py.in
BUILT_SOURCES = GC_autoconf.py 

edit = sed -e 's|@pkgdatadir[@]|$(pkgdatadir)|g'

GC_autoconf.py:	Makefile GC_autoconf.py.in
		rm -f $@
		$(edit) '$(srcdir)/$@.in' > $@
