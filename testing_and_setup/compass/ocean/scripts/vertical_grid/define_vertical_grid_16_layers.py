#!/usr/bin/env python

from make_vertical_grid import create_vertical_grid

# Standard 16 layer vertical grid
create_vertical_grid(
    num_vert_levels=16,
    bottom_depth=3000,
    min_layer_thickness=3.,
    max_layer_thickness=500.,
    plot_vertical_grid=True,
    outfile='vertical_grid.nc')
