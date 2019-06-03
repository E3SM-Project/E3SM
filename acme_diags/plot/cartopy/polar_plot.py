from __future__ import print_function

import os
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import matplotlib.path as mpath
from acme_diags.driver.utils.general import get_output_dir
from acme_diags.plot import get_colormap

plotTitle = {'fontsize': 11.5}
plotSideTitle = {'fontsize': 9.0}

# Position and sizes of subplot axes in page coordinates (0 to 1)
panel = [(0.27, 0.65, 0.3235, 0.25),
         (0.27, 0.35, 0.3235, 0.25),
         (0.27, 0.05, 0.3235, 0.25),
         ]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
border = (-0.02, -0.01, 0.14, 0.04)


def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0], lon[0] + 360.0, 'coe'))


def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height


def plot_panel(n, fig, proj, pole, var, clevels,
               cmap, title, parameters, stats=None):

    var = add_cyclic(var)
    lon = var.getLongitude()
    lat = var.getLatitude()
    var = ma.squeeze(var.asma())

    # Contour levels
    levels = None
    norm = None
    if len(clevels) > 0:
        levels = [-1.0e8] + clevels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

    # Contour plot
    ax = fig.add_axes(panel[n], projection=proj)
    ax.set_global()

    ax.gridlines()
    if pole == 'N':
        ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
    elif pole == 'S':
        ax.set_extent([-180, 180, -55, -90], crs=ccrs.PlateCarree())

    cmap = get_colormap(cmap, parameters)
    p1 = ax.contourf(lon, lat, var,
                     transform=ccrs.PlateCarree(),
                     norm=norm,
                     levels=levels,
                     cmap=cmap,
                     extend='both',
                     )
    ax.set_aspect('auto')
    ax.coastlines(lw=0.3)

    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # Plot titles
    for i in range(3):
        if title[i] is None:
            title[i] = ''

    t = ax.set_title('%s\n%s' %
                     (title[0], title[2]), loc='left', fontdict=plotSideTitle)

    t = ax.set_title(title[1], fontdict=plotTitle)
    if title[0] != '' or title[2] != '':
        t.set_position([.5, 1.06])

    # Color bar
    cbax = fig.add_axes(
        (panel[n][0] + 0.35, panel[n][1] + 0.0354, 0.0326, 0.1792))
    cbar = fig.colorbar(p1, cax=cbax)
    w, h = get_ax_size(fig, cbax)

    if levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    else:
        maxval = np.amax(np.absolute(levels[1:-1]))
        if maxval < 10.0:
            fmt = "%5.2f"
            pad = 25
        elif maxval < 100.0:
            fmt = "%5.1f"
            pad = 25
        else:
            fmt = "%6.1f"
            pad = 30
        cbar.set_ticks(levels[1:-1])
        labels = [fmt % l for l in levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha='right')
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)

    # Min, Mean, Max
    fig.text(panel[n][0] + 0.35, panel[n][1] + 0.225,
             "Max\nMean\nMin", ha='left', fontdict=plotSideTitle)
    fig.text(panel[n][0] + 0.45, panel[n][1] + 0.225, "%.2f\n%.2f\n%.2f" %
             stats[0:3], ha='right', fontdict=plotSideTitle)

    # RMSE, CORR
    if len(stats) == 5:
        fig.text(panel[n][0] + 0.35, panel[n][1] + 0.,
                 "RMSE\nCORR", ha='left', fontdict=plotSideTitle)
        fig.text(panel[n][0] + 0.45, panel[n][1] + 0., "%.2f\n%.2f" %
                 stats[3:5], ha='right', fontdict=plotSideTitle)


def plot(reference, test, diff, metrics_dict, parameter):

    # Create figure, projection
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    # Create projection
    if parameter.var_region.find('N') != -1:
        pole = 'N'
        proj = ccrs.NorthPolarStereo(central_longitude=0)
    elif parameter.var_region.find('S') != -1:
        pole = 'S'
        proj = ccrs.SouthPolarStereo(central_longitude=0)

    # First two panels
    min1 = metrics_dict['test']['min']
    mean1 = metrics_dict['test']['mean']
    max1 = metrics_dict['test']['max']
    if test.count() > 1:
        plot_panel(0, fig, proj, pole, test, parameter.contour_levels, parameter.test_colormap, [
                   parameter.test_name_yrs, parameter.test_title, test.units], parameter, stats=(max1, mean1, min1))

    min2 = metrics_dict['ref']['min']
    mean2 = metrics_dict['ref']['mean']
    max2 = metrics_dict['ref']['max']

    if reference.count() > 1:
        plot_panel(1, fig, proj, pole, reference, parameter.contour_levels, parameter.reference_colormap, [
                   parameter.ref_name_yrs, parameter.reference_title, reference.units], parameter, stats=(max2, mean2, min2))

    # Third panel
    min3 = metrics_dict['diff']['min']
    mean3 = metrics_dict['diff']['mean']
    max3 = metrics_dict['diff']['max']
    r = metrics_dict['misc']['rmse']
    c = metrics_dict['misc']['corr']

    if diff.count() > 1:
        plot_panel(2, fig, proj, pole, diff, parameter.diff_levels, parameter.diff_colormap, [
                   None, parameter.diff_title, None], parameter, stats=(max3, mean3, min3, r, c))

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=18)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir(parameter.current_set,
            parameter), parameter.output_file + '.' + f)
        plt.savefig(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
            ignore_container=True), parameter.output_file + '.' + f)
        print('Plot saved in: ' + fnm)

    # Save individual subplots
    for f in parameter.output_format_subplot:
        fnm = os.path.join(get_output_dir(
            parameter.current_set, parameter), parameter.output_file)
        page = fig.get_size_inches()
        i = 0
        for p in panel:
            # Extent of subplot
            subpage = np.array(p).reshape(2,2)
            subpage[1,:] = subpage[0,:] + subpage[1,:]
            subpage = subpage + np.array(border).reshape(2,2)
            subpage = list(((subpage)*page).flatten())
            extent = matplotlib.transforms.Bbox.from_extents(*subpage)
            # Save subplot
            fname = fnm + '.%i.' %(i) + f
            plt.savefig(fname, bbox_inches=extent)
            
            orig_fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
                ignore_container=True), parameter.output_file)
            fname = orig_fnm + '.%i.' %(i) + f
            print('Sub-plot saved in: ' + fname)
            
            i += 1

    plt.close()
