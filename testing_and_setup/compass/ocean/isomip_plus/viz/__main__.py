#!/usr/bin/env python

import xarray
import argparse

from viz.streamfunction import compute_barotropic_streamfunction, \
    compute_overturning_streamfunction

from viz.plot import MoviePlotter, TimeSeriesPlotter

from viz.misomip import compute_misomip_interp_coeffs, interp_misomip
from viz.haney import compute_haney_number


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--folder", dest="folder",
                        help="Folder for plots", default='.')
    parser.add_argument("-e", "--expt", dest="expt",
                        help="Experiment number (0, 1 or 2)", default=0)
    parser.add_argument("--misomip", dest="misomip", action="store_true",
                        help="Indicates that standard MISOMIP output files "
                             "should be created")
    parser.add_argument("--streamfunctions", dest="streamfunctions", 
                        action="store_true",
                        help="Indicates that the barotropic and overturning "
                             "streamfunctions should be computed and plotted")
    parser.add_argument("--haney", dest="haney", action="store_true",
                        help="Indicates that the Haney number rx1 should be "
                             "computed and plotted")
    args = parser.parse_args()

    folder = args.folder
    expt = args.expt

    dsMesh = xarray.open_dataset('{}/init.nc'.format(folder))

    ds = xarray.open_mfdataset('{}/timeSeriesStatsMonthly*.nc'.format(folder),
                               concat_dim='Time')

    if args.streamfunctions:
        compute_barotropic_streamfunction(dsMesh, ds, folder)

        compute_overturning_streamfunction(dsMesh, ds, folder, dx=2e3, dz=5.)

    if args.haney:
        compute_haney_number(dsMesh, ds, folder)

    tsPlotter = TimeSeriesPlotter(inFolder=folder,
                                  outFolder='{}/plots'.format(folder),
                                  expt=expt)
    tsPlotter.plot_melt_time_series()
    tsPlotter = TimeSeriesPlotter(inFolder=folder,
                                  outFolder='{}/timeSeriesBelow300m'.format(folder),
                                  expt=expt)
    tsPlotter.plot_melt_time_series(sshMax=-300.)

    mPlotter = MoviePlotter(inFolder=folder,
                           outFolder='{}/plots'.format(folder),
                           expt=expt)

    if args.streamfunctions:
        mPlotter.plot_barotropic_streamfunction()
        mPlotter.plot_overturning_streamfunction()

    if args.haney:
        mPlotter.plot_haney_number()

    mPlotter.plot_melt_rates()
    mPlotter.plot_ice_shelf_boundary_variables()
    mPlotter.plot_temperature()
    mPlotter.plot_salinity()
    mPlotter.plot_potential_density()

    mPlotter.images_to_movies(outFolder='{}/movies'.format(folder),
                              framesPerSecond=30, extension='mp4')

    if args.misomip:
        compute_misomip_interp_coeffs(folder, expt)
        interp_misomip(folder, expt)


if __name__ == '__main__':
    main()
