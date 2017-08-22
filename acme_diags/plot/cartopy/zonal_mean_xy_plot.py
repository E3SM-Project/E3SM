from __future__ import print_function

import os
import cdms2
import cdutil
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from acme_diags.metrics import rmse, corr
from acme_diags.driver.utils import get_output_dir, _chown

plotTitle = {'fontsize':12.5}
plotSideTitle = {'fontsize':11.5}

#panel = [(0.1691,0.4961,0.6465,0.2258),
#         (0.1691,0.2112,0.6465,0.2258),
#         ]
panel = [(0.1191,0.5461,0.6465*1.2,0.2258*1.3),
         (0.1191,0.1612,0.6465*1.2,0.2258*1.3),
         ]

def get_ax_size(fig,ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

def plot(reference, test, diff, metrics_dict, parameter):

    # Create figure, projection
    fig = plt.figure(figsize=[8.5, 11.0])
#    proj = ccrs.PlateCarree(central_longitude=180)
    proj = None
    ax = fig.add_axes(panel[0])
    ax.plot(test.getLatitude()[:],ma.squeeze( test.asma() ),'k',linewidth=2)
    
    ax.plot(reference.getLatitude()[:],ma.squeeze( reference.asma()),'r',linewidth=2)
    ax1 = fig.add_axes(panel[1])
    ax1.plot(diff.getLatitude()[:],ma.squeeze( diff.asma() ),'k', linewidth=2) 
    fig.text(panel[0][0],panel[0][2]+0.095,"Test: " + parameter.test_name,ha='left',fontdict=plotSideTitle,color = 'black')
    fig.text(panel[0][0],panel[0][2]+0.07,"Reference: " +parameter.reference_name,ha='left',fontdict=plotSideTitle,color='red')
    fig.text(panel[1][0],panel[0][2]-0.3,"Test-Reference",ha='left',fontdict=plotSideTitle)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])#crs=ccrs.PlateCarree())
    ax.set_xlim(-90,90)
    ax1.set_xticks([-90, -60, -30, 0, 30, 60, 90])#crs=ccrs.PlateCarree())
    ax1.set_xlim(-90,90)
##    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.0f')
    lat_formatter = LatitudeFormatter()
#    ax.xaxis.set_major_formatter(lon_formatter)
    #ax.xaxis.set_major_formatter(lat_formatter)
    #ax1.xaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=11.0, direction='out', pad=-2, width=1)
    ax1.tick_params(labelsize=11.0, direction='out', pad=-2, width=1)
    ax.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_ticks_position('bottom')
    ax.set_ylabel(test.long_name +' ('+ test.units+')')
    ax1.set_ylabel(test.long_name +' ('+ test.units+')')
#    ax.yaxis.set_ticks_position('left')
#
    fig.suptitle(parameter.main_title, x=0.5, y=0.95, fontsize=18)
#    plt.show()
#
    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter), parameter.output_file)
        plt.savefig(fnm + '.' + f)
        _chown(fnm + '.' + f, parameter.user)
        print('Plot saved in: ' + fnm + '.' + f)
