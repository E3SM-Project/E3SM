TEST_CONFIGS = erosion0.config erosion1.config erosion2.config erosion3.config erosion4.config \
               er-transport0.config er-transport1.config er-transport2.config er-transport3.config er-transport4.config

EXTRA_DIST = topo.nc $(TEST_CONFIGS)

data = $(subst config,nc,$(TEST_CONFIGS))

MODEL_BINARY=$(top_builddir)/example-drivers/simple_glide/src/simple_glide

include $(top_srcdir)/extra_rules.am
