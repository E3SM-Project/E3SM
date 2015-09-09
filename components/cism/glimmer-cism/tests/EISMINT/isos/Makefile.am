#if BUILD_GRAPHICS
#  graphics = spots.ps profile.ps
#endif
data = isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
BUILT_SOURCES = 

EXTRA_DIST = isos.elfa.config  isos.elra.config  isos.llfa.config  isos.llra.config

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/simple_glide/src/simple_glide

spots.ps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
	@GLIMMER_PREFIX@/bin/plotSpot.py --ij 15 15 --ij 15 22 --ij 15 27 -vtopg $^ $@

profile.ps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
	@GLIMMER_PREFIX@/bin/plotProfile.py -t 18  -v topg -p ../PyCF/eismint.prof $^ $@

include $(top_srcdir)/extra_rules.am
