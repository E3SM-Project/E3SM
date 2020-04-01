#!/usr/bin/env python

from make_vertical_grid import create_vertical_grid

# Standard 80 layer vertical grid
create_vertical_grid(
    num_vert_levels=80,
    bottom_depth=6500,
    min_layer_thickness=1.0,
    max_layer_thickness=200.,
    plot_vertical_grid=True,
    outfile='vertical_grid.nc')
