import os
import cdms2
import cdutil
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from acme_diags.metrics import rmse, corr
from acme_diags.driver.utils import get_output_dir

plotTitle = {'fontsize':11.5}
plotSideTitle = {'fontsize':9.5}

panel = [(0.1691,0.6810,0.6465,0.2258),
         (0.1691,0.3961,0.6465,0.2258),
         (0.1691,0.1112,0.6465,0.2258),
         ]


def get_ax_size(fig,ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

def plot_panel(n, fig, proj, var, clevels, cmap, title, stats=None):

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
    plt.pcolormesh(var2d)
    # Calculate 3 x 3 grids for cloud fraction for nine cloud class
    # Place cloud fraction of each cloud class in plot:
    cld_3x3= np.zeros((3,3))
    for j in range(3):
        for i in range(3):
            if j==2:
                cld_3x3[j,i]=var2d[4:7,2*i:2*i+2].sum()
                ax.text(i*2+1, j*2+1.5,'%.1f' %cld_3x3[j,i],horizontalalignment='center',verticalalignment='center',fontsize=25)
            else:
                cld_3x3[j,i]=var2d[2*j:2*j+2,2*i:2*i+2].sum()
                ax.text(i*2+1, j*2+1,'%.1f' %cld_3x3[j,i],horizontalalignment='center',verticalalignment='center',fontsize=25)
    
    # Place vertical/horizonal line to separate cloud class
    
    plt.axvline(x=2, linewidth=2, color='k')
    plt.axvline(x=4, linewidth=2, color='k')
    plt.axhline(y=2, linewidth=2, color='k')
    plt.axhline(y=4, linewidth=2, color='k')
    
    xlabels = ['%.1f' %i for i in x.getBounds()[:,0]]+['%.1f' %x.getBounds()[-1,-1]]
    ylabels = ['%.1f' %i for i in y.getBounds()[:,0]]+['%.1f' %y.getBounds()[-1,-1]]
    
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)
    #p1 = ax.contourf(lon, lat, var,
    #                 transform=ccrs.PlateCarree(), 
    #                 norm=norm,
    #                 levels=levels,
    #                 cmap=cmap,
    #                 extend='both',
    #                 )
    #ax.set_aspect('auto')
    #ax.coastlines(lw=0.3)
    #if title[0] != None: ax.set_title(title[0], loc='left', fontdict=plotSideTitle)
    #if title[1] != None: ax.set_title(title[1], fontdict=plotTitle)
    #if title[2] != None: ax.set_title(title[2], loc='right', fontdict=plotSideTitle)
    #ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
    ##ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    #ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    #lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.0f')
    #lat_formatter = LatitudeFormatter()
    #ax.xaxis.set_major_formatter(lon_formatter)
    #ax.yaxis.set_major_formatter(lat_formatter)
    #ax.tick_params(labelsize=8.0, direction='out', pad=-2, width=1)
    #ax.xaxis.set_ticks_position('bottom')
    #ax.yaxis.set_ticks_position('left')

    # Color bar
    cbax = fig.add_axes((panel[n][0]+0.6635,panel[n][1]+0.0215,0.0326,0.1792))
    cbar = fig.colorbar(p1, cax=cbax)
    w,h = get_ax_size(fig,cbax)

    if levels == None:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    else:
        cbar.set_ticks(levels[1:-1])
        labels = ["%4.1f" % l for l in levels[1:-1]]
        cbar.ax.set_yticklabels(labels,ha='right')
        #cbar.ax.set_yticklabels(labels,ha='right')
        cbar.ax.tick_params( labelsize=9.0,pad=25, length=0)


def plot(reference, test, diff, parameter):

    # Create figure, projection
    fig = plt.figure(figsize=[8.5, 11.0])

    plot_panel(0, fig, proj, test, parameter.contour_levels, 'viridis', (parameter.test_name,parameter.test_title,test.units))

#    min2  = metrics_dict['ref']['min']
#    mean2 = metrics_dict['ref']['mean']
#    max2  = metrics_dict['ref']['max']
#    plot_panel(1, fig, proj, reference, parameter.contour_levels, 'viridis', (parameter.reference_name,parameter.reference_title,reference.units),stats=(max2,mean2,min2))
#
#    # Third panel
#    min3  = metrics_dict['diff']['min']
#    mean3 = metrics_dict['diff']['mean']
#    max3  = metrics_dict['diff']['max']
#    r = metrics_dict['misc']['rmse']
#    c = metrics_dict['misc']['corr']
#    plot_panel(2, fig, proj, diff, parameter.diff_levels, 'RdBu_r', (None,parameter.diff_title,None), stats=(max3,mean3,min3,r,c))
#
#    # Figure title
#    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter), parameter.output_file)
        plt.savefig(fnm + '.' + f)
        print('Plot saved in: ' + fnm + '.' + f)
