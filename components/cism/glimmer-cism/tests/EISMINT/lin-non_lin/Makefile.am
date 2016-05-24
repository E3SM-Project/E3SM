runs = l-5a.nc nl-5a.nc l-10a.nc nl-10a.nc l-20a.nc nl-20a.nc l-30a.nc nl-30a.nc

#if BUILD_GRAPHICS
#  graphics = e2-comp.ps e2-comp-t.ps e2-btemp.ps e2-bvel.ps
#endif
data = $(runs)
BUILT_SOURCES = 

EXTRA_DIST = l-10a.config  l-20a.config  l-5a.config  nl-10a.config  nl-20a.config  nl-5a.config l-30a.config  nl-30a.config move_data.sh

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/simple_glide/src/simple_glide

e2-comp.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e H $^ $@

e2-comp-t.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotCFstats.py -e -m -p ../PyCF/eismint.prof $^ $@

e2-btemp.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotCFvar.py --width 5 -T -1 -cthk -vtemp -l -1 $^ $@

e2-bvel.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotCFvar.py --width 5 -T -1 -cthk -vbvel -l -1 $^ $@

.config.nc:
	@GLIMMER_PREFIX@/bin/glide_launch.py --prefix @GLIMMER_PREFIX@ $<
	bash move_data.sh picard_info.data $*.piccard

include $(top_srcdir)/extra_rules.am
