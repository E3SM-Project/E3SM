#if BUILD_GRAPHICS
#  graphics = e1-fm.1.ps e1-fm.2.ps e1-fm.3.ps e1-mm.1.ps e1-mm.2.ps e1-mm.3.ps
#endif
data = e1-fm.1.nc e1-fm.2.nc e1-fm.3.nc e1-mm.1.nc e1-mm.2.nc e1-mm.3.nc
BUILT_SOURCES =  

EXTRA_DIST = e1-fm.1.config  e1-fm.2.config  e1-fm.3.config  e1-mm.1.config  e1-mm.2.config  e1-mm.3.config

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/simple_glide/src/simple_glide

include $(top_srcdir)/extra_rules.am

.nc.ps:
	@GLIMMER_PREFIX@/bin/plotEISMINT.py -t 200 --pmt  -p ../PyCF/eismint.prof $< $@
