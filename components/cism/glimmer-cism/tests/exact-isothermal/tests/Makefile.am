EXTRA_DIST = index.html template.config error.cpt

test_b = test_B_non-lin_50km_10.00a.nc test_B_lin_50km_10.00a.nc \
         test_B_non-lin_20km_2.00a.nc  test_B_lin_20km_2.00a.nc \
         test_B_non-lin_10km_0.50a.nc  test_B_lin_10km_0.50a.nc \
         test_B_non-lin_5km_0.10a.nc   test_B_lin_5km_0.10a.nc 
test_c = test_C_non-lin_50km_5.00a.nc test_C_lin_50km_5.00a.nc \
         test_C_non-lin_20km_1.00a.nc test_C_lin_20km_1.00a.nc \
         test_C_non-lin_10km_0.50a.nc test_C_lin_10km_0.50a.nc \
         test_C_non-lin_5km_0.10a.nc  test_C_lin_5km_0.10a.nc 
test_d = test_D_non-lin_50km_5.00a.nc test_D_lin_50km_5.00a.nc \
         test_D_non-lin_20km_1.00a.nc test_D_lin_20km_1.00a.nc

CLEANFILES = $(all_tests)

all_tests = $(test_b) $(test_c) $(test_d)

graphics:	vole_c.png vole_d.png errors.png time.png $(subst nc,png,$(all_tests))

errors.png:	$(all_tests)
		$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/plot_errors.py -o $@ $^

time.png:	results
		$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/plot_time.py -o $@ $^

vole_b.png:	$(test_b)
		$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/plot_vol.py -o $@ $^

vole_c.png:	$(test_c)
		$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/plot_vol.py -o $@ $^

vole_d.png:	$(test_d)
		$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/plot_vol.py -o $@ $^

results:	$(all_tests)

%.png:	%.nc
	$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/plot_diff.py -T-1 -o $@ $< 

%.nc:	$(top_srcdir)/tests/exact-isothermal/tests/template.config
	$(PYTHON) $(top_srcdir)/tests/exact-isothermal/scripts/run_verif.py -m ../src/verif_glide @SGE_OPTS@ -f $* $<


