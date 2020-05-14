from __future__ import print_function

import matplotlib
import numpy as np
import os
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from acme_diags.driver.utils.general import get_output_dir

panel = [(0.075, 0.70, 0.6, 0.225),
         (0.075, 0.425, 0.6, 0.225),
         (0.725, 0.425, 0.2, 0.5),
         (0.075, 0.075, 0.85, 0.275)]

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
border = (-0.06, -0.03, 0.13, 0.03)


def plot_panel(n, fig, plot_type, label_size, title, x, y, z=None, plot_colors=None, color_levels=None, color_ticks=None):
    # x,y,z should be of the form:
    # dict(axis_range=None, axis_scale=None, data=None, data_label=None, data2=None, data2_label=None, label=None)

    # Create new figure axis using dimensions from panel (hard coded)
    ax = fig.add_axes(panel[n])
    # Plot either a contourf or line plot
    if plot_type == 'contourf':
        if 'data' not in z:
            raise RuntimeError('Must set z["data"] to use plot_type={}.'.format(plot_type))
        p1 = ax.contourf(x['data'], y['data'], z['data'], color_levels, cmap=plot_colors)
        cbar = plt.colorbar(p1, ticks=color_ticks)
        cbar.ax.tick_params(labelsize=label_size)
    if plot_type == 'line':
        if 'data2' not in x or 'data2' not in y:
            raise RuntimeError('Must set data2 for both x and y to use plot_type={}.'.format(plot_type))
        elif 'data_label' not in x or 'data2_label' not in x:
            raise RuntimeError('Must set data_label and data2_label for x to use plot_type={}.'.format(plot_type))
        p1, = ax.plot(x['data'], y['data'], '-ok')
        p2, = ax.plot(x['data2'], y['data2'], '--or')
        plt.grid('on')
        ax.legend((p1, p2), (x['data_label'], x['data2_label']), loc='upper right', fontsize=label_size)
    ax.set_title(title, size=label_size, weight='demi')
    ax.set_xlabel(x['label'], size=label_size)
    ax.set_ylabel(y['label'], size=label_size)
    plt.yscale(y['axis_scale'])
    plt.ylim([y['axis_range'][0], y['axis_range'][1]])
    plt.yticks(size=label_size)
    plt.xscale(x['axis_scale'])
    plt.xlim([x['axis_range'][0], x['axis_range'][1]])
    plt.xticks(size=label_size)
    return ax


def plot(period_new, parameter, test, ref):
    label_size = 14

    fig = plt.figure(figsize=(14, 14))

    months = min(ref['qbo'].shape[0], test['qbo'].shape[0])
    x_test, y_test = np.meshgrid(np.arange(0, months), test['level'])
    x_ref, y_ref = np.meshgrid(np.arange(0, months), ref['level'])
    cmap2 = plt.cm.RdBu_r
    color_levels0 = np.arange(-50, 51, 100. / 20.)

    # Panel 0 (Top Left)
    x = dict(axis_range=[0, months], axis_scale='linear', data=x_test, label=' ')
    y = dict(axis_range=[100, 1], axis_scale='log', data=y_test, label='hPa')
    z = dict(data=test['qbo'].T)
    ax0 = plot_panel(0, fig, 'contourf', label_size, '{} U 5S-5N'.format(test['name']), x, y, z=z,
               plot_colors=cmap2, color_levels=color_levels0, color_ticks=[-50, -25, -5, 5, 25, 50])
    # Panel 1 (Middle Left)
    x = dict(axis_range=[0, months], axis_scale='linear', data=x_ref, label='month')
    y = dict(axis_range=[100, 1], axis_scale='log', data=y_ref, label='hPa')
    z = dict(data=ref['qbo'].T)
    ax1 = plot_panel(1, fig, 'contourf', label_size, '{} U 5S-5N'.format(ref['name']), x, y, z=z,
               plot_colors=cmap2, color_levels=color_levels0, color_ticks=[-50, -25, -5, 5, 25, 50])
    # Panel 2 (Top/Middle Right)
    # TODO: Are the test/ref values correct for data, data_label, data2, data2_label??
    x = dict(axis_range=[0, 30], axis_scale='linear', data=ref['amplitude'][:], data_label=ref['name'],
             data2=test['amplitude'][:], data2_label=test['name'], label='Amplitude (m/s)')
    y = dict(axis_range=[100, 1], axis_scale='log', data=ref['level'][:],
             data2=ref['level'][:], label='Pressure (hPa)')
    ax2 = plot_panel(2, fig, 'line', label_size, 'QBO Amplitude \n (period = 20-40 months)', x, y)
    # Panel 3 (Bottom)
    x = dict(axis_range=[0, 50], axis_scale='linear', data=period_new, data_label=ref['name'],
             data2=period_new, data2_label=test['name'], label='Period (months)')
    y = dict(axis_range=[-1, 25], axis_scale='linear', data=ref['amplitude_new'],
             data2=test['amplitude_new'], label='Amplitude (m/s)')
    ax3 = plot_panel(3, fig, 'line', label_size, 'QBO Spectral Density (Eq. 18-22 hPa zonal winds)', x, y)
    plt.tight_layout()

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.97, fontsize=15)

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/qbo/{parameter.case_id} 
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print('Output dir: {}'.format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/qbo/{parameter.case_id} 
    original_output_dir = get_output_dir(parameter.current_set, parameter, ignore_container=True)
    if parameter.print_statements:
        print('Original output dir: {}'.format(original_output_dir))
    # parameter.output_file is defined in acme_diags/driver/qbo_driver.py
    # {parameter.results_dir}/qbo/{parameter.case_id}/{parameter.output_file}
    file_path = os.path.join(output_dir, parameter.output_file)
    # {parameter.orig_results_dir}/qbo/{parameter.case_id}/{parameter.output_file}
    original_file_path = os.path.join(original_output_dir, parameter.output_file)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        plot_suffix = '.' + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        original_plot_file_path = original_file_path + plot_suffix
        print('Plot saved in: ' + original_plot_file_path)

    # TODO: The subplots don't come out so nicely. Same for ENSO Diags
    #  See Issue 294
    # Save individual subplots
    for f in parameter.output_format_subplot:
        # i = 0
        # for ax in [ax0, ax1, ax2, ax3]:
        #     # Extent of subplot
        #     # https://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib
        #     extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        #     # Save subplot
        #     subplot_suffix = ('.%i.' % i) + f
        #     subplot_file_path = file_path + subplot_suffix
        #     plt.savefig(subplot_file_path, bbox_inches=extent.expanded(1.1, 1.2))
        #     i += 1
        #     # Get the filename that the user has passed in and display that.
        #     # When running in a container, the paths are modified.
        #     original_subplot_file_path = original_file_path + subplot_suffix
        #     print('Sub-plot saved in: ' + original_subplot_file_path)
        page = fig.get_size_inches()
        i = 0
        for p in panel:
            # Extent of subplot
            subpage = np.array(p).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(border).reshape(2, 2)
            subpage = list(((subpage) * page).flatten())
            extent = matplotlib.transforms.Bbox.from_extents(*subpage)
            # Save subplot
            subplot_suffix = ('.%i.' % i) + f
            subplot_file_path = file_path + subplot_suffix
            plt.savefig(subplot_file_path, bbox_inches=extent)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            original_subplot_file_path = original_file_path + subplot_suffix
            print('Sub-plot saved in: ' + original_subplot_file_path)
            i += 1

    plt.close()
