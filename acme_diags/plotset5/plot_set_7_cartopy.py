import os
import cdms2
import cdutil
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

from acme_diags.metrics import rmse, corr

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

    ax.gridlines()
    if pole == 'N':
      ax.set_extent([-180, 180, 48, 90], crs=ccrs.PlateCarree())
    elif pole == 'S':
      ax.set_extent([-180, 180, -48, -90], crs=ccrs.PlateCarree())

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
    if title[1] != None: 
        t = ax.set_title(title[1], fontdict=plotTitle)
        if title[0] != None or title[2] != None: t.set_position([.5, 1.06])
    if title[2] != None: ax.set_title(title[2], loc='right', fontdict=plotSideTitle)

    # Color bar
    cbax = fig.add_axes((panel[n][0]+0.35,panel[n][1]+0.0354,0.0326,0.1792))
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
    fig.text(panel[n][0]+0.35,panel[n][1]+0.225,"Max\nMean\nMin",ha='left',fontdict=plotSideTitle)
    fig.text(panel[n][0]+0.45,panel[n][1]+0.225,"%.2f\n%.2f\n%.2f" % (var_max,var_mean,var_min),ha='right',fontdict=plotSideTitle)

    # RMSE, CORR
    if stats != None:
      fig.text(panel[n][0]+0.35,panel[n][1]+0.,"RMSE\nCORR",ha='left',fontdict=plotSideTitle)
      fig.text(panel[n][0]+0.45,panel[n][1]+0.,"%.2f\n%.2f" % stats,ha='right',fontdict=plotSideTitle)

def plot(reference, test, reference_regrid, test_regrid, parameter, pole):

    # Create figure, projection
    fig = plt.figure(figsize=[8.5, 11.0])

    # Create projection
    if pole == 'N':
        proj = ccrs.NorthPolarStereo(central_longitude=0)
    elif pole == 'S':
       proj = ccrs.SouthPolarStereo(central_longitude=0)

    # First two panels
    plot_panel(0, fig, proj, pole, test, parameter.test_levels, 'viridis', (parameter.test_name,parameter.test_title,test.units))
    plot_panel(1, fig, proj, pole, reference, parameter.reference_levels, 'viridis', (parameter.reference_name,parameter.reference_title,reference.units))

    # Third panel
    r = rmse(reference_regrid, test_regrid)
    c = corr(reference_regrid, test_regrid)
    plot_panel(2, fig, proj, pole, test_regrid-reference_regrid, parameter.diff_levels, 'RdBu_r', (None,parameter.diff_title,None), stats=(r,c))

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=18)

    # Save figure
    case_id = parameter.case_id
    if not os.path.exists(case_id):
        os.makedirs(case_id)
    for f in parameter.format:
        plt.savefig(case_id + '/' + parameter.output_file + '_' + pole + '.' + f)
