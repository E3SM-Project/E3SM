#!/usr/bin/env python

from make_vertical_grid import create_vertical_grid

# Standard 64 layer vertical grid
create_vertical_grid(bottom_depth=5500, nz=64, dz1_in=2., dz2_in=200.,
                     plot_vertical_grid=True, outfile= 'vertical_grid.nc')
