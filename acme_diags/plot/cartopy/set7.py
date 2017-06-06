import os
import cdms2
import cdutil
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import matplotlib.path as mpath

from acme_diags.metrics import rmse, corr
from acme_diags.driver.utils import get_output_dir



plotTitle = {'fontsize':11.5}
plotSideTitle = {'fontsize':9.0}

panel = [(0.27,0.65,0.3235,0.25),
         (0.27,0.35,0.3235,0.25),
         (0.27,0.05,0.3235,0.25),
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

def plot_panel(n, fig, proj, pole, var, clevels, cmap, title, stats=None):

#    var_min = float(var.min())
#    var_max = float(var.max())
#    var_mean = cdutil.averager(var, axis='xy', weights='generate')

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

    ax.gridlines()
    if pole == 'N':
      ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
    elif pole == 'S':
      ax.set_extent([-180, 180, -55, -90], crs=ccrs.PlateCarree())


    p1 = ax.contourf(lon, lat, var,
                     transform=ccrs.PlateCarree(), 
                     norm=norm,
                     levels=levels,
                     cmap=cmap,
                     extend='both',
                     )
    ax.set_aspect('auto')
    ax.coastlines(lw=0.3)

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    if title[0] != None: ax.set_title(title[0], loc='left', fontdict=plotSideTitle)
    if title[1] != None: 
        t = ax.set_title(title[1], fontdict=plotTitle)
        if title[0] != None or title[2] != None: t.set_position([.5, 1.06])
    if title[2] != None: ax.set_title(title[2], loc='right', fontdict=plotSideTitle)

    # Color bar
    cbax = fig.add_axes((panel[n][0]+0.35,panel[n][1]+0.0354,0.0326,0.1792))
    cbar = fig.colorbar(p1, cax=cbax)
    w,h = get_ax_size(fig,cbax)

    if levels == None:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    else:
        cbar.set_ticks(levels[1:-1])
        labels = ["%4.1f" % l for l in levels[1:-1]]
        cbar.ax.set_yticklabels(labels,ha='right')
        cbar.ax.tick_params(labelsize=9.0, pad=25, length=0)

    # Min, Mean, Max
    fig.text(panel[n][0]+0.35,panel[n][1]+0.225,"Max\nMean\nMin",ha='left',fontdict=plotSideTitle)
    fig.text(panel[n][0]+0.45,panel[n][1]+0.225,"%.2f\n%.2f\n%.2f" % stats[0:3],ha='right',fontdict=plotSideTitle)

    # RMSE, CORR
    if len(stats) == 5:
      fig.text(panel[n][0]+0.35,panel[n][1]+0.,"RMSE\nCORR",ha='left',fontdict=plotSideTitle)
      fig.text(panel[n][0]+0.45,panel[n][1]+0.,"%.2f\n%.2f" % stats[3:5],ha='right',fontdict=plotSideTitle)

def plot(reference, test, diff, metrics_dict, parameter):

    # Create figure, projection
    fig = plt.figure(figsize=[8.5, 11.0])

    # Create projection
    print parameter.var_region
    if parameter.var_region.find('N') !=-1:
        pole = 'N'
        proj = ccrs.NorthPolarStereo(central_longitude=0)
    elif parameter.var_region.find('S') !=-1:
        pole = 'S'
        proj = ccrs.SouthPolarStereo(central_longitude=0)

    # First two panels
    min1  = metrics_dict['test']['min']
    mean1 = metrics_dict['test']['mean']
    max1  = metrics_dict['test']['max']
    if test.count() >1: 
        plot_panel(0, fig, proj, pole, test, parameter.contour_levels, 'viridis', (parameter.test_name,parameter.test_title,test.units),stats=(max1,mean1,min1))

    min2  = metrics_dict['ref']['min']
    mean2 = metrics_dict['ref']['mean']
    max2  = metrics_dict['ref']['max']

    if reference.count() >1: 
        plot_panel(1, fig, proj, pole, reference, parameter.contour_levels, 'viridis', (parameter.reference_name,parameter.reference_title,reference.units),stats=(max2,mean2,min2))

    # Third panel
    min3  = metrics_dict['diff']['min']
    mean3 = metrics_dict['diff']['mean']
    max3  = metrics_dict['diff']['max']
    r = metrics_dict['misc']['rmse']
    c = metrics_dict['misc']['corr']

    if diff.count() >1: 
        plot_panel(2, fig, proj, pole, diff, parameter.diff_levels, 'RdBu_r', (None,parameter.diff_title,None), stats=(max3,mean3,min3,r,c))

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=18)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir('7', parameter), parameter.output_file)
        plt.savefig(fnm + '.' + f)
        print('Plot saved in: ' + fnm + '.' + f)
