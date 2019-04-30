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
from subprocess import call

# render statically by default
plt.switch_backend('agg')

def setup_fig():
    fig, ax = plt.subplots(nrows=2,ncols=1, sharex=True, sharey=True)
    fig.text(0.04, 0.5, 'Channel depth (m)', va='center', rotation='vertical')
    fig.text(0.5, 0.02, 'Along channel distance (km)', ha='center')

def setup_subplot():
    plt.xlim(0,25)
    plt.ylim(-1, 11)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    x = np.linspace(0,25,100)
    y = 10.0/25.0*x
    plt.plot(x, y, 'k-', lw=3)


def upper_plot():
    plt.subplot(2,1,1)
    plt.gca().set_xticklabels([])
    setup_subplot()


def lower_plot():
    plt.subplot(2,1,2)
    setup_subplot()


def plot_data(rval='0.0025', dtime='0.05', datatype='analytical', *args, **kwargs):
    datafile = 'data/r' + rval + 'd' + dtime + '-' + datatype + '.csv'
    data = pd.read_csv(datafile, header=None)
    measured=plt.scatter(data[0], data[1], *args, **kwargs)


def plot_datasets(rval, times, fileno, plotdata=True, frame=False):
    for ii, dtime in enumerate(times):
        if plotdata:
            plot_data(rval=rval, dtime = dtime, datatype = 'analytical',
                      marker = '.', color = 'b', label='analytical')
            plot_data(rval=rval, dtime = dtime, datatype = 'roms',
                      marker = '.', color = 'g', label='ROMS')
        if frame:
            timeslice = times[frame]
            plot_MPASO(timeslice, fileno, 'k-', lw=2, label='MPAS-O')
        else:
            plot_MPASO([dtime], fileno, 'k-', lw=0.5, label='MPAS-O')

        if plotdata:
            if ii == 0:
                plt.legend(frameon=False, loc='lower left')
                place_time_labels(times)
                plt.text(0.5, 5, 'r = ' + str(rval))


def place_time_labels(times):
    locs = [9.3, 7.2, 4.2, 2.2, 1.2, 0.2]
    for atime, ay in zip(times, locs):
        plt.text(25.2, ay, atime + ' days', size=8)

def plot_MPASO(times, fileno, *args, **kwargs):
    for atime in times:
        # factor of 1e-16 needed to account for annoying round-off issue to get right time slices
        plottime = int((float(atime)/0.2 + 1e-16)*24.0)
        ds = xr.open_mfdataset('output'+ fileno + '.nc')
        timestr = ds.isel(Time=plottime).xtime.values
        #print('{} {} {}'.format(atime, plottime, timestr))
        ds = ds.drop(np.setdiff1d([i for i in ds.variables], ['yCell','ssh']))
        ymean = ds.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
        x = ymean.yCell.values/1000.0
        y = ymean.ssh.values
        #print('ymin={} ymax={}\n{}\n{}'.format(y.min(), y.max(),x, y))
        plt.plot(x, -y, *args, **kwargs)
        ds.close()
        return timestr


def plot_tidal_forcing_comparison():
    # data from MPAS-O on boundary
    for fileno in ['1','2']:
        ds = xr.open_mfdataset('output' + fileno +'.nc')
        ympas = ds.ssh.where(ds.tidalInputMask).mean('nCells').values
        x = np.linspace(0, 1.0, len(ds.xtime))*12.0
        plt.plot(x, ympas, marker='o', label='MPAS-O forward' + fileno)

    # analytical case
    x = np.linspace(0,12.0,100)
    y = 10.0*np.sin(x*np.pi/12.0) - 10.0
    plt.plot(x, y, lw=3, color='black', label='analytical')

    plt.legend(frameon=False)
    plt.ylabel('Tidal amplitude (m)')
    plt.xlabel('Time (hrs)')
    plt.suptitle('Tidal amplitude forcing (right side) for MPAS-O and analytical')
    plt.savefig('tidalcomparison.png')


def make_movie(moviename):
   
    # assumes automatic rewrite of the file
    args = ['ffmpeg', '-y',
            '-framerate', '6',
            '-pattern_type', 'glob', 
            '-i', moviename + '*.png',
            '-c:v', 'libx264', 
            '-r', '30',
            '-profile:v', 'high',
            '-crf', '20', 
            '-pix_fmt', 'yuv420p', 
            moviename +  '.mp4']
    call(args)
    

def main():
    ################ plot tidal forcing comparison ###############
    plot_tidal_forcing_comparison()
    ##############################################################

    ################ plot drying front comparison ###############
    setup_fig()
    times = ['0.50', '0.05', '0.40', '0.15', '0.30', '0.25']

    # subplot r = 0.0025
    upper_plot()
    plot_datasets(rval='0.0025', times=times, fileno='1')

    # subplot r = 0.01
    lower_plot()
    plot_datasets(rval='0.01', times=times, fileno='2')

    plt.suptitle('Drying slope comparison between MPAS-O, analytical, and ROMS')
    for outtype in ['.png','.pdf']:
        plt.savefig('dryingslopecomparison' + outtype)

    ###################### make movie ##############################
    ii = 0
    frames = []
    moviename = 'dryingslopecomparisonmovie'
    for atime  in np.linspace(0,0.5,5*12+1):
        # plot snapshots
        print(atime)

        # initialize movie frame
        plt.close('all')
        setup_fig()
        # subplot r = 0.0025
        upper_plot()
        plot_datasets(rval='0.0025', times=times, fileno='1')

        # subplot r = 0.01
        lower_plot()
        plot_datasets(rval='0.01', times=times, fileno='2')

        # plot animated line on top image
        # subplot r = 0.0025
        upper_plot()
        timestr = plot_MPASO([atime], '1', 'k-', lw=2.0, label='MPAS-O')
        #plt.text(10, 10, 't = ' + str(timestr)[10:])
        #plt.text(10, 10, 't = {:1.4f} days'.format(atime))

        # subplot r = 0.01
        lower_plot()
        timestr = plot_MPASO([atime], '2', 'k-', lw=2.0, label='MPAS-O')
        #plt.text(10, 10, 't = ' + str(timestr)[10:])
        plt.text(10, 10, 't = {:1.2f} hrs'.format(atime*12))

        plt.suptitle('Drying slope comparison between MPAS-O, analytical, and ROMS')
        
        framename = '{}{:03d}.png'.format(moviename,ii)
        plt.savefig(framename)
        
        frames.append(framename)
        ii += 1

    make_movie(moviename)
    ##############################################################


if __name__ == "__main__":
    main()
