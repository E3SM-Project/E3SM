include(ScreamUtils)

# Equilibrium water-isotope fractionation factor unit tests.
CreateUnitTest(water_isotopes_fractionation
  SOURCES eamxx_water_isotopes_fractionation_tests.cpp
  LIBS water_isotopes
  LABELS water_isotopes physics
)

# Runtime formulation selection unit tests.
CreateUnitTest(water_isotopes_runtime_options
  SOURCES water_isotopes_runtime_options_tests.cpp
  LIBS water_isotopes
  LABELS water_isotopes physics
)

# Isotope ocean evaporation flux: equilibrium invariants and H216O identity.
CreateUnitTest(water_isotopes_ocean_flux
  SOURCES water_isotopes_ocean_flux_tests.cpp
  LIBS water_isotopes
  LABELS water_isotopes physics
)
