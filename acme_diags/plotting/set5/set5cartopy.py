import os
import cdms2
import cdutil
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from acme_diags.metrics import rmse, corr

plotTitle = {'fontsize':11.5}
plotSideTitle = {'fontsize':9.5}

panel = [(0.1691,0.6810,0.6465,0.2258),
         (0.1691,0.3961,0.6465,0.2258),
         (0.1691,0.1112,0.6465,0.2258),
         ]

def add_cyclic(var):
    lon = var.getLongitude() 
    return var(longitude=(lon[0],lon[0]+360.0,'coe'))

def get_ax_size(fig,ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

def plot_panel(n, fig, proj, var, clevels, cmap, title, stats=None):

    var_min = float(var.min())
    var_max = float(var.max())
    var_mean = cdutil.averager(var, axis='xy', weights='generate')
    var = add_cyclic(var)
    lon = var.getLongitude()
    lat = var.getLatitude()
    var = ma.squeeze( var.asma() )

    # Contour levels
    levels = None
    norm = None
    if len(clevels) > 0:
        levels = [-1.0e8] + clevels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

    # Contour plot
    ax = fig.add_axes(panel[n],projection=proj)
    ax.set_global()
    p1 = ax.contourf(lon, lat, var,
                     transform=ccrs.PlateCarree(), 
                     norm=norm,
                     levels=levels,
                     cmap=cmap,
                     extend='both',
                     )
    ax.set_aspect('auto')
    ax.coastlines(lw=0.3)
    if title[0] != None: ax.set_title(title[0], loc='left', fontdict=plotSideTitle)
    if title[1] != None: ax.set_title(title[1], fontdict=plotTitle)
    if title[2] != None: ax.set_title(title[2], loc='right', fontdict=plotSideTitle)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
    #ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction='out', pad=-2, width=1)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Color bar
    cbax = fig.add_axes((panel[n][0]+0.6635,panel[n][1]+0.0215,0.0326,0.1792))
    cbar = fig.colorbar(p1, cax=cbax)
    w,h = get_ax_size(fig,cbax)

    if levels == None:
        cbar.ax.tick_params(labelsize=9.0, length=w-2)

    else:
        cbar.set_ticks(levels[1:-1])
        labels = ["%4.1f" % l for l in levels[1:-1]]
        cbar.ax.set_yticklabels(labels,ha='right')
        cbar.ax.tick_params(labelsize=9.0, pad=25, length=w-2)

    # Min, Mean, Max
    fig.text(panel[n][0]+0.6635,panel[n][1]+0.2107,"Max\nMean\nMin",ha='left',fontdict=plotSideTitle)
    fig.text(panel[n][0]+0.7635,panel[n][1]+0.2107,"%.2f\n%.2f\n%.2f" % (var_max,var_mean,var_min),ha='right',fontdict=plotSideTitle)

    # RMSE, CORR
    if stats != None:
      fig.text(panel[n][0]+0.6635,panel[n][1]-0.0105,"RMSE\nCORR",ha='left',fontdict=plotSideTitle)
      fig.text(panel[n][0]+0.7635,panel[n][1]-0.0105,"%.2f\n%.2f" % stats,ha='right',fontdict=plotSideTitle)

def plot(reference, test, diff, metrics_dict, parameter):

    # Create figure, projection
    fig = plt.figure(figsize=[8.5, 11.0])
    proj = ccrs.PlateCarree(central_longitude=180)

    # First two panels
    plot_panel(0, fig, proj, test, parameter.contour_levels, 'viridis', (parameter.test_name,parameter.test_title,test.units))
    plot_panel(1, fig, proj, reference, parameter.contour_levels, 'viridis', (parameter.reference_name,parameter.reference_title,reference.units))

    # Third panel
    r = metrics_dict['misc']['rmse']
    c = metrics_dict['misc']['corr']
    plot_panel(2, fig, proj, diff, parameter.diff_levels, 'RdBu_r', (None,parameter.diff_title,None), stats=(r,c))

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # Save figure
    case_id = parameter.case_id
    if not os.path.exists(case_id):
        os.makedirs(case_id)
    
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        plt.savefig(case_id + '/' + parameter.output_file + '.' + f)

        print('Plot saved in: ' + case_id + '/' + parameter.output_file + '.' + f)
