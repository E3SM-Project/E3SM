# TODO: Figure out how to generate pio.F90 automatically from pio.F90.in

list(APPEND pio_sources
  pio.F90
  )

sourcelist_to_parent(pio_sources)