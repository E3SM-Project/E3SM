graphics = stress-strain.eps stress-strain-fit.eps

EXTRA_DIST = stress-strain.gp fit.gp stress-strain-fit.gp stress-strain.data 

BUILT_SOURCES = $(graphics)

stress-strain.gp:	fit.gp	stress-strain.data
stress-strain-fit.gp:	fit.gp	stress-strain.data

include $(top_srcdir)/doc/make_doc.rules
