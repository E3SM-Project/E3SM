#!/usr/bin/env python

import make_vertical_grid

# Standard 64 layer vertical grid
make_vertical_grid.create_vertical_grid(
     5500,            # bottom_depth
     64,              # num_vert_levels
     2,               # layer1_thickness
     200,             # maxLayer_thickness
     True,            # plot_vertical_grid
     maxit=1000,
     epsilon=1.0e-2,
     outFile= 'vertical_grid.nc')# output_file_name  
