from __future__ import print_function

import os
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LatitudeFormatter
from acme_diags.driver.utils.general import get_output_dir

plotTitle = {'fontsize': 12.5}
plotSideTitle = {'fontsize': 11.5}

# Position and sizes of subplot axes in page coordinates (0 to 1)
panel = [(0.1500, 0.5500, 0.7500, 0.3000),
         (0.1500, 0.1300, 0.7500, 0.3000),
         ]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
border = (-0.14, -0.06, 0.04, 0.08)

def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height


def plot(reference, test, diff, metrics_dict, parameter):

    # Create figure
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    # Top panel
    ax1 = fig.add_axes(panel[0])
    ax1.plot(test.getLatitude()[:], ma.squeeze(test.asma()), 'k', linewidth=2)
    ax1.plot(reference.getLatitude()[:], ma.squeeze(
        reference.asma()), 'r', linewidth=2)
    ax1.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax1.set_xlim(-90, 90)
    ax1.tick_params(labelsize=11.0, direction='out', width=1)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_ylabel(test.long_name + ' (' + test.units + ')')

    test_title = "Test" if parameter.test_title == '' else parameter.test_title
    test_title += ' : {}'.format(parameter.test_name_yrs)
    ref_title = "Reference" if parameter.reference_title == '' else parameter.reference_title
    ref_title += ' : {}'.format(parameter.ref_name_yrs)
    fig.text(panel[0][0], panel[0][1]+panel[0][3]+0.03, test_title,
             ha='left', fontdict=plotSideTitle, color='black')
    fig.text(panel[0][0], panel[0][1]+panel[0][3]+0.01, ref_title,
             ha='left', fontdict=plotSideTitle, color='red')

    # Bottom panel
    ax2 = fig.add_axes(panel[1])
    ax2.plot(diff.getLatitude()[:], ma.squeeze(diff.asma()), 'k', linewidth=2)
    ax2.axhline(y=0, color='0.5')
    ax2.set_title(parameter.diff_title, fontdict=plotSideTitle, loc='center') 
    ax2.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax2.set_xlim(-90, 90)
    ax2.tick_params(labelsize=11.0, direction='out', width=1)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_ylabel(test.long_name + ' (' + test.units + ')')

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.95, fontsize=18)

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
