#!/usr/bin/env python

from make_vertical_grid import create_vertical_grid

# Standard 64 layer vertical grid
create_vertical_grid(
    num_vert_levels=64,
    bottom_depth=6000,
    min_layer_thickness=2.,
    max_layer_thickness=210.,
    plot_vertical_grid=True,
    outfile='vertical_grid.nc')
