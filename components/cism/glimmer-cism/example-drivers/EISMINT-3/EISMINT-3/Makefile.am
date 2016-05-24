data = e3gl_steady.nc
graphics = topo.eps

BUILT_SOURCES = 

EXTRA_DIST = e3gl_steady.config EISMINT3GL.input.nc

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/EISMINT-3/src/eismint3_glide

include $(top_srcdir)/extra_rules.am

topo.eps:	EISMINT3GL.input.nc
		$(plotGCvar) -v topg --legend=v $^ -o $@
