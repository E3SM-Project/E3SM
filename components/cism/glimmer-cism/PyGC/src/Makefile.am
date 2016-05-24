pypkglibdir = @pythondir@/PyGC
pypkglib_LTLIBRARIES  = _proj4.la

_proj4_la_SOURCES = projmodule.c
_proj4_la_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(PYTHON_SITE_PKG)/numpy/core/include -I$(PYTHON_SITE_PKG)/numpy/numarray/ 
_proj4_la_LDFLAGS = -module -shared -avoid-version 
_proj4_la_LIBADD = -lproj
