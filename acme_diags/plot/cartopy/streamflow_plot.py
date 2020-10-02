from __future__ import print_function

import matplotlib
import matplotlib.lines as lines
import numpy as np
import os

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cdutil
from acme_diags.derivations.default_regions import regions_specs
from acme_diags.driver.utils.general import get_output_dir

plotTitle = {'fontsize': 11.5}
plotSideTitle = {'fontsize': 9.5}

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, width, height) in page coordinates
border = (-0.14, -0.06, 0.04, 0.08)


def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0], lon[0] + 360.0, 'coe'))


def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height


def determine_tick_step(degrees_covered):
    if degrees_covered > 180:
        return 60
    if degrees_covered > 60:
        return 30
    elif degrees_covered > 20:
        return 10
    else:
        return 1


def plot_panel_seasonality(plot_type, fig, proj, export, color_list, parameter):
    if plot_type == 'test':
        panel_index = 0
        seasonality_index_export_index = 5
        peak_month_export_index = 6
        title = (None, parameter.test_title, None)
    elif plot_type == 'ref':
        panel_index = 1
        seasonality_index_export_index = 3
        peak_month_export_index = 4
        title = (None, parameter.reference_title, None)
    else:
        raise Exception('Invalid plot_type={}'.format(plot_type))

    # Plot of streamflow gauges. Color -> peak month, marker size -> seasonality index.

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    panel = [(0.0900, 0.5500, 0.7200, 0.3000),
             (0.0900, 0.1300, 0.7200, 0.3000),
             ]
    ax = fig.add_axes(panel[panel_index], projection=proj)
    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    if 'domain' in region.keys():
        # Get domain to plot
        domain = region['domain']
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90., 90, 'ccb'))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if 'longitude' in kargs:
        lon_west, lon_east, _ = kargs['longitude']
    if 'latitude' in kargs:
        lat_south, lat_north, _ = kargs['latitude']
    lon_covered = lon_east - lon_west
    lon_step = determine_tick_step(lon_covered)
    xticks = np.arange(lon_west, lon_east, lon_step)
    # Subtract 0.50 to get 0 W to show up on the right side of the plot.
    # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the left side of the plot.
    # If a number is added, then the value won't show up at all.
    xticks = np.append(xticks, lon_east - 0.50)
    lat_covered = lat_north - lat_south
    lat_step = determine_tick_step(lat_covered)
    yticks = np.arange(lat_south, lat_north, lat_step)
    yticks = np.append(yticks, lat_north)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=proj)
    proj_function = ccrs.PlateCarree

    # Stream gauges
    si_2 = 2
    si_4 = 3
    si_6 = 4
    si_large = 5
    # `export` is the array of gauges. Each gauge has multiple fields -- e.g., lat is index 7
    for gauge in export:
        lat = gauge[7]
        lon = gauge[8]
        seasonality_index = gauge[seasonality_index_export_index]
        if seasonality_index < 2:
            markersize = si_2
        elif seasonality_index < 4:
            markersize = si_4
        elif seasonality_index < 6:
            markersize = si_6
        elif seasonality_index <= 12:
            markersize = si_large
        else:
            raise Exception('Invalid seasonality index={}'.format(seasonality_index))
        if seasonality_index == 1:
            color = 'black'
        else:
            peak_month = int(gauge[peak_month_export_index])
            color = color_list[peak_month]
        # https://scitools.org.uk/iris/docs/v1.9.2/examples/General/projections_and_annotations.html
        # Place a single marker point for each gauge.
        plt.plot(lon, lat, marker='o', color=color, markersize=markersize, transform=proj_function())
        # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
        # so for this one we transform the coordinates with a Cartopy call.
        at_x, at_y = ax.projection.transform_point(lon, lat, src_crs=proj_function())
    # https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/custom_legends.html
    legend_elements = [lines.Line2D([0], [0], marker='o', color='w', label='1 <= SI < 2',
                                    markerfacecolor='black', markersize=si_2),
                       lines.Line2D([0], [0], marker='o', color='w', label='2 <= SI < 4',
                                    markerfacecolor='black', markersize=si_4),
                       lines.Line2D([0], [0], marker='o', color='w', label='4 <= SI < 6',
                                    markerfacecolor='black', markersize=si_6),
                       lines.Line2D([0], [0], marker='o', color='w', label='6 <= SI <= 12',
                                    markerfacecolor='black', markersize=si_large)
                       ]
    seasonality_legend_title = 'Seasonality (SI)'
    plt.legend(handles=legend_elements, title=seasonality_legend_title, prop={'size': 8})

    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)
    if title[0] is not None:
        ax.set_title(title[0], loc='left', fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    if title[2] is not None:
        ax.set_title(title[2], loc='right', fontdict=plotSideTitle)
    ax.set_xticks(xticks, crs=proj_function())
    ax.set_yticks(yticks, crs=proj_function())
    lon_formatter = LongitudeFormatter(
        zero_direction_label=True, number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction='out', width=1)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Color bar
    cbax = fig.add_axes(
        (panel[panel_index][0] + 0.7535, panel[panel_index][1] + 0.0515, 0.0326, 0.1792))
    # https://matplotlib.org/tutorials/colors/colorbar_only.html
    num_colors = len(color_list)
    if parameter.print_statements:
        print('num_colors={}'.format(num_colors))
    cmap = colors.ListedColormap(color_list)
    cbar_label = 'Peak month'

    bounds = list(range(num_colors))
    # Set ticks to be in between the bounds
    ticks = list(map(lambda bound: bound + 0.5, bounds))
    # Add one more bound at the bottom of the colorbar.
    # `bounds` should be one longer than `ticks`.
    bounds += [bounds[-1] + 1]
    if parameter.print_statements:
        print('bounds={}'.format(bounds))
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cbax,
        boundaries=bounds,
        ticks=ticks,
        spacing='uniform',
        orientation='vertical',
        label=cbar_label,
    )
    # https://matplotlib.org/3.1.1/gallery/ticks_and_spines/colorbar_tick_labelling_demo.html
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    cbar.ax.set_yticklabels(months)
    cbar.ax.invert_yaxis()

    w, h = get_ax_size(fig, cbax)

    cbar.ax.tick_params(labelsize=9.0, length=0)


def plot_seasonality(export, parameter):
    if parameter.backend not in ['cartopy', 'mpl', 'matplotlib']:
        return

    # Create figure, projection
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=0)

    if parameter.test_title == '':
        parameter.test_title = "Test"
    if parameter.reference_title == '':
        parameter.reference_title = "Reference"
    # test and ref color lists
    color_list = [(0.12,0.96,0.87), (0.16,1.00,0.47), (0.22,1.00,0.19),
                  (0.70,1.00,0.20), (1.00,0.77,0.18), (1.00,0.20,0.11),
                  (1.00,0.05,0.30), (1.00,0.13,0.84), (0.61,0.15,1.00),
                  (0.07,0.14,0.99), (0.09,0.56,1.00), (0.16,0.91,0.99)]

    # First panel
    plot_panel_seasonality('test', fig, proj, export, color_list, parameter)

    # Second panel
    plot_panel_seasonality('ref', fig, proj, export, color_list, parameter)

    # Figure title
    fig.suptitle(parameter.main_title_seasonality, x=0.5, y=0.97, fontsize=15)

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/streamflow/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print('Output dir: {}'.format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/streamflow/{parameter.case_id}
    original_output_dir = get_output_dir(parameter.current_set, parameter, ignore_container=True)
    if parameter.print_statements:
        print('Original output dir: {}'.format(original_output_dir))
    # parameter.output_file_seasonality is defined in acme_diags/parameter/streamflow_parameter.py
    # {parameter.results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_seasonality}
    file_path = os.path.join(output_dir, parameter.output_file_seasonality)
    # {parameter.orig_results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_seasonality}
    original_file_path = os.path.join(original_output_dir, parameter.output_file_seasonality)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        plot_suffix = '.' + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        original_plot_file_path = original_file_path + plot_suffix
        # Always print, even without `parameter.print_statements`
        print('Plot saved in: ' + original_plot_file_path)

    # Save individual subplots
    for f in parameter.output_format_subplot:
        page = fig.get_size_inches()
        i = 0
        for p in panel:
            # Extent of subplot
            subpage = np.array(p).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(border).reshape(2, 2)
            subpage = list((subpage * page).flatten())
            extent = matplotlib.transforms.Bbox.from_extents(*subpage)
            # Save subplot
            subplot_suffix = '.%i.' % i + f
            subplot_file_path = file_path + subplot_suffix
            plt.savefig(subplot_file_path, bbox_inches=extent)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            original_subplot_file_path = original_file_path + subplot_suffix
            # Always print, even without `parameter.print_statements`
            print('Sub-plot saved in: ' + original_subplot_file_path)
            i += 1

    plt.close()


def plot_panel_bias(fig, proj, export, bias_array, parameter):
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    # panel = [(0.0900, 0.5500, 0.7200, 0.3000),
    #          (0.0900, 0.1300, 0.7200, 0.3000),
    #          ]
    panel = [(0.0900, 0.2000, 0.7200, 0.6000)]
    panel_index = 0
    title = (None, parameter.test_title + '\n' + parameter.reference_title, None)

    # Plot of streamflow gauges. Color -> peak month, marker size -> seasonality index.

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    ax = fig.add_axes(panel[panel_index], projection=proj)
    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    if 'domain' in region.keys():
        # Get domain to plot
        domain = region['domain']
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90., 90, 'ccb'))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if 'longitude' in kargs:
        lon_west, lon_east, _ = kargs['longitude']
    if 'latitude' in kargs:
        lat_south, lat_north, _ = kargs['latitude']
    lon_covered = lon_east - lon_west
    lon_step = determine_tick_step(lon_covered)
    xticks = np.arange(lon_west, lon_east, lon_step)
    # Subtract 0.50 to get 0 W to show up on the right side of the plot.
    # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the left side of the plot.
    # If a number is added, then the value won't show up at all.
    xticks = np.append(xticks, lon_east - 0.50)
    lat_covered = lat_north - lat_south
    lat_step = determine_tick_step(lat_covered)
    yticks = np.arange(lat_south, lat_north, lat_step)
    yticks = np.append(yticks, lat_north)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=proj)
    proj_function = ccrs.PlateCarree

    # Stream gauges
    discharge_10 = 2
    discharge_100 = 3
    discharge_1000 = 4
    discharge_10000 = 5
    discharge_large = 6
    # `export` is the array of gauges. Each gauge has multiple fields -- e.g., lat is index 7
    # Continuous colormap
    colormap = plt.get_cmap('jet')
    color_list = list(map(lambda index: colormap(index)[:3], range(colormap.N)))
    bias_min = -100
    bias_max = 100
    for gauge, i in zip(export, range(len(export))):
        discharge = gauge[1]
        bias = bias_array[i]
        if bias < bias_min:
            bias = bias_min
        if bias > bias_max:
            bias = bias_max
        if numpy.isnan(bias):
            continue
        lat = gauge[7]
        lon = gauge[8]
        if discharge < 10:
            markersize = discharge_10
        elif discharge < 100:
            markersize = discharge_100
        elif discharge < 1000:
            markersize = discharge_1000
        elif discharge < 10000:
            markersize = discharge_10000
        else:
            markersize = discharge_large
        # Rescale (min-max normalization)
        normalized_bias = (bias - bias_min) / (bias_max - bias_min)
        color = color_list[int(normalized_bias*(len(color_list)-1))]
        # https://scitools.org.uk/iris/docs/v1.9.2/examples/General/projections_and_annotations.html
        # Place a single marker point for each gauge.
        plt.plot(lon, lat, marker='o',  markersize=markersize, color=color, transform=proj_function())
        # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
        # so for this one we transform the coordinates with a Cartopy call.
        at_x, at_y = ax.projection.transform_point(lon, lat, src_crs=proj_function())

    # https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/custom_legends.html
    legend_elements = [
        lines.Line2D([0], [0], marker='o', color='w', label='<10', markerfacecolor='black', markersize=discharge_10),
        lines.Line2D([0], [0], marker='o', color='w', label='<100', markerfacecolor='black', markersize=discharge_100),
        lines.Line2D([0], [0], marker='o', color='w', label='<1000', markerfacecolor='black', markersize=discharge_1000),
        lines.Line2D([0], [0], marker='o', color='w', label='<10000', markerfacecolor='black', markersize=discharge_10000),
        lines.Line2D([0], [0], marker='o', color='w', label='>10000', markerfacecolor='black', markersize=discharge_large)
                       ]
    discharge_legend_title = 'Discharge $(m^3/s)$'
    plt.legend(handles=legend_elements, title=discharge_legend_title, prop={'size': 8})

    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)
    if title[0] is not None:
        ax.set_title(title[0], loc='left', fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    if title[2] is not None:
        ax.set_title(title[2], loc='right', fontdict=plotSideTitle)
    ax.set_xticks(xticks, crs=proj_function())
    ax.set_yticks(yticks, crs=proj_function())
    lon_formatter = LongitudeFormatter(
        zero_direction_label=True, number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction='out', width=1)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Color bar
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    cbax = fig.add_axes(
        (panel[panel_index][0] + 0.7535, panel[panel_index][1] + 0.0515, 0.0326, 0.1792*2))
    cmap = colors.ListedColormap(color_list)
    cbar_label = 'Bias of annual mean discharge'
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap),
        cax=cbax,
        label=cbar_label,
        extend='both'
    )
    w, h = get_ax_size(fig, cbax)
    step_size = (bias_max - bias_min)//5
    ticks = numpy.arange(int(bias_min), int(bias_max) + step_size, step_size)
    cbar.ax.tick_params(labelsize=9.0, length=0)
    cbar.ax.set_yticklabels(ticks)


def plot_bias(export, bias, parameter):
    if parameter.backend not in ['cartopy', 'mpl', 'matplotlib']:
        return

    # Create figure, projection
    figsize = [parameter.figsize[0], parameter.figsize[1]/2]
    fig = plt.figure(figsize=figsize, dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=0)

    if parameter.test_title == '':
        parameter.test_title = "Test"
    if parameter.reference_title == '':
        parameter.reference_title = "Reference"

    # Only panel
    plot_panel_bias(fig, proj, export, bias, parameter)

    # Figure title
    fig.suptitle(parameter.main_title_bias, x=0.5, y=0.97, fontsize=15)

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/streamflow/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print('Output dir: {}'.format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/streamflow/{parameter.case_id}
    original_output_dir = get_output_dir(parameter.current_set, parameter, ignore_container=True)
    if parameter.print_statements:
        print('Original output dir: {}'.format(original_output_dir))
    # parameter.output_file_bias is defined in acme_diags/parameter/streamflow_parameter.py
    # {parameter.results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_bias}
    file_path = os.path.join(output_dir, parameter.output_file_bias)
    # {parameter.orig_results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_bias}
    original_file_path = os.path.join(original_output_dir, parameter.output_file_bias)

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        plot_suffix = '.' + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        original_plot_file_path = original_file_path + plot_suffix
        # Always print, even without `parameter.print_statements`
        print('Plot saved in: ' + original_plot_file_path)

    # Save individual subplots
    for f in parameter.output_format_subplot:
        page = fig.get_size_inches()
        i = 0
        for p in panel:
            # Extent of subplot
            subpage = np.array(p).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(border).reshape(2, 2)
            subpage = list((subpage * page).flatten())
            extent = matplotlib.transforms.Bbox.from_extents(*subpage)
            # Save subplot
            subplot_suffix = '.%i.' % i + f
            subplot_file_path = file_path + subplot_suffix
            plt.savefig(subplot_file_path, bbox_inches=extent)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            original_subplot_file_path = original_file_path + subplot_suffix
            # Always print, even without `parameter.print_statements`
            print('Sub-plot saved in: ' + original_subplot_file_path)
            i += 1

    plt.close()
