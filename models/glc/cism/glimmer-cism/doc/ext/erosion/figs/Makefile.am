graphics = ice_velo.eps erosion_layers.eps pressure_stress_strain.eps plot_basal.eps
EXTRA_DIST = ice_velo.fig erosion_layers.fig pressure_stress_strain.fig plot_basal.py
BUILT_SOURCES = $(graphics)

plot_basal.eps:	plot_basal.py
	$(PYTHON) plot_basal.py $@

sediments.ps:	../../../tests/fenscan.nc
		$(PYTHON) ../../../python/plot_seds.py -T -1 -p fenscan.prof --not_p $< $@

sed_profile.ps:	../../../tests/fenscan.nc
		$(PYTHON) ../../../python/plot_seds_profile.py -t0. -e stages -p fenscan.prof --not_p $< $@

include $(top_srcdir)/doc/make_doc.rules
