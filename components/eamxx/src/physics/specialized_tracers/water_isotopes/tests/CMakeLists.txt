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
