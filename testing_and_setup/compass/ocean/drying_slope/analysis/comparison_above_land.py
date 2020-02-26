#!/usr/bin/env python
"""

Drying slope comparison betewen MPAS-O, analytical, and ROMS results from

Warner, J. C., Defne, Z., Haas, K., & Arango, H. G. (2013). A wetting and
drying scheme for ROMS. Computers & geosciences, 58, 54-61.

Phillip J. Wolfram and Zhendong Cao
04/30/2019

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

# render statically by default
plt.switch_backend('agg')

def setup_fig():
    fig, ax = plt.subplots(nrows=1,ncols=1, sharex=True, sharey=True)
    fig.text(0.04, 0.5, 'Channel depth (m)', va='center', rotation='vertical')
    fig.text(0.5, 0.02, 'Along channel distance (km)', ha='center')

def setup_subplot():
    plt.xlim(0,25)
    plt.ylim(-11, 11)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    x = np.linspace(0,25,100)
    y = -10 + 20.0/25.0*x
    plt.plot(x, y, 'k-', lw=3)

def plot_data(rval='0.0025', dtime='0.05', datatype='analytical', *args, **kwargs):
    datafile = 'data/r' + rval + 'd' + dtime + '-' + datatype + '.csv'
    data = pd.read_csv(datafile, header=None)
    measured=plt.scatter(data[0], data[1], *args, **kwargs)


def plot_datasets(rval, times, fileno, plotdata=True):
    for ii, dtime in enumerate(times):
        if plotdata:
            plot_data(rval=rval, dtime = dtime, datatype = 'analytical',
                      marker = '.', color = 'b', label='analytical')
            plot_data(rval=rval, dtime = dtime, datatype = 'roms',
                      marker = '.', color = 'g', label='ROMS')
        plot_MPASO([dtime], fileno, 'k-', lw=0.5, label='MPAS-O')

        if plotdata:
            if ii == 0:
                plt.legend(frameon=False, loc='lower left')
                place_time_labels(times)
                plt.text(0.5, 5, 'r = ' + str(rval))


def place_time_labels(times):
    locs = [9.3, 7.2, 4.2, 2.2, 1.2, 0.2]
    for atime, ay in zip(times, locs):
        plt.text(25.2, -10+20*ay, atime + ' days', size=8)

def plot_MPASO(times, fileno, *args, **kwargs):
    for atime in times:
        # factor of 1e-16 needed to account for annoying round-off issue to get right time slices
        plottime = int((float(atime)/0.2 + 1e-16)*24.0)
        ds = xr.open_dataset('output'+ fileno + '.nc')
        print('{} {} {}'.format(atime, plottime, ds.isel(Time=plottime).xtime.values))
        ds = ds.drop_vars(np.setdiff1d([i for i in ds.variables], ['yCell','ssh']))
        ymean = ds.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
        x = ymean.yCell.values/1000.0
        y = ymean.ssh.values
        #print('ymin={} ymax={}\n{}\n{}'.format(y.min(), y.max(),x, y))
        plt.plot(x, -y, *args, **kwargs)
        ds.close()


def plot_tidal_forcing_comparison():
    # data from MPAS-O on boundary
    for fileno in ['']:
        ds = xr.open_dataset('output' + fileno +'.nc')
        ympas = ds.ssh.where(ds.tidalInputMask).mean('nCells').values
        x = np.linspace(0, 1.0, len(ds.xtime))*12.0
        plt.plot(x, ympas, marker='o', label='MPAS-O forward' + fileno)

    # analytical case
    x = np.linspace(0,12.0,100)
    y = 20.0*np.sin(x*np.pi/12.0) - 10.0
    plt.plot(x, y, lw=3, color='black', label='analytical')

    plt.legend(frameon=False)
    plt.ylabel('Tidal amplitude (m)')
    plt.xlabel('Time (hrs)')
    plt.suptitle('Tidal amplitude forcing (right side) for MPAS-O and analytical')
    plt.savefig('tidalcomparison.png')


def main():
    ################ plot tidal forcing comparison ###############
    plot_tidal_forcing_comparison()
    ##############################################################

    ################ plot drying front comparison ###############
    setup_fig()
    times = ['0.50', '0.05', '0.40', '0.15', '0.30', '0.25']

    # subplot r = 0.01
    setup_subplot()
    plot_datasets(rval='0.01', times=times, fileno='', plotdata=False)

    plt.suptitle('Drying slope for land topography')
    for outtype in ['.png','.pdf']:
        plt.savefig('dryingslopecomparison' + outtype)
    ##############################################################


if __name__ == "__main__":
    main()
