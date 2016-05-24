#graphics = e2.a.eps e2.b.eps e2.c.eps e2.d.eps e2.f.eps e2.g.eps e2.h.eps e2.i.eps e2.j.eps e2.k.eps e2.l.eps e2.a-comp.eps e2.b-comp.eps e2.c-comp.eps e2.d-comp.eps e2.g-comp.eps e2.h-comp.eps
graphics = e2.a-comp.eps e2.b-comp.eps e2.c-comp.eps e2.d-comp.eps e2.g-comp.eps e2.h-comp.eps
data = e2.a.nc e2.b.nc e2.c.nc e2.d.nc e2.f.nc e2.g.nc e2.h.nc e2.i.nc e2.j.nc e2.k.nc e2.l.nc
BUILT_SOURCES = 

EXTRA_DIST = e2.a.config  e2.b.config  e2.c.config  e2.d.config  e2.f.config  e2.g.config  e2.h.config  e2.i.config  e2.j.config  e2.k.config  e2.l.config mound.nc trough.nc

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/simple_glide/src/simple_glide

e2.b.nc:        e2.a.nc
e2.c.nc:        e2.a.nc
e2.d.nc:        e2.a.nc
e2.j.nc:        e2.i.nc
e2.l.nc:        e2.k.nc

e2.a-comp.eps:	e2.a.nc
	$(plotEISMINT2stats) -e A $^ -o $@
e2.b-comp.eps:	e2.b.nc e2.a.nc
	$(plotEISMINT2stats) -e B $^ -o $@
e2.c-comp.eps:	e2.c.nc e2.a.nc
	$(plotEISMINT2stats) -e C $^ -o $@
e2.d-comp.eps:	e2.d.nc e2.a.nc
	$(plotEISMINT2stats) -e D $^ -o $@
e2.g-comp.eps:	e2.g.nc
	$(plotEISMINT2stats) -e G $^ -o $@
e2.h-comp.eps:	e2.h.nc
	$(plotEISMINT2stats) -e H $^ -o $@

include $(top_srcdir)/extra_rules.am

