# Append sources from this directory to the list.
list(APPEND circle_area_sources circle.F90)

# Export source list to parent (in CMake, "export" means something
# different).
sourcelist_to_parent(circle_area_sources)
