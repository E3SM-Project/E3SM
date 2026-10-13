include(ScreamUtils)

# wtrc_ratio and bulk_ocean_flux, incl. FPE safety (run in the fpe build too).
CreateUnitTest(water_tracers_functions
  SOURCES water_tracers_functions_tests.cpp
  LIBS water_tracers
  LABELS water_tracers physics
)

# Generic and property tracer ocean fluxes, and the invalid-ocean rule.
CreateUnitTest(water_tracers_ocean_flux
  SOURCES water_tracers_ocean_flux_tests.cpp
  LIBS water_tracers
  LABELS water_tracers physics
)

# tracers_file parsing and validation. The test writes its YAML fixtures
# into the working directory.
CreateUnitTest(water_tracers_registry
  SOURCES water_tracers_registry_tests.cpp
  LIBS water_tracers
  LABELS water_tracers physics
)
