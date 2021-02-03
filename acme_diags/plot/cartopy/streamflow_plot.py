from __future__ import print_function

import matplotlib
import matplotlib.lines as lines
import numpy as np
import os

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy
import scipy.stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cdutil
from acme_diags.derivations.default_regions import regions_specs
from acme_diags.driver.utils.general import get_output_dir

plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}

# Border padding relative to subplot axes for saving individual panels
# (left, bottom, width, height) in page coordinates
border = (-0.14, -0.06, 0.04, 0.08)


def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0], lon[0] + 360.0, "coe"))


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


def plot_panel_seasonality_map(
    plot_type, fig, proj, export, color_list, panel, parameter
):
    if plot_type == "test":
        panel_index = 0
        seasonality_index_export_index = 5
        peak_month_export_index = 6
        title = (None, parameter.test_title, None)
    elif plot_type == "ref":
        panel_index = 1
        seasonality_index_export_index = 3
        peak_month_export_index = 4
        title = (None, parameter.reference_title, None)
    else:
        raise Exception("Invalid plot_type={}".format(plot_type))

    # Plot of streamflow gauges. Color -> peak month, marker size -> seasonality index.
    ax = fig.add_axes(panel[panel_index], projection=proj)
    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    if "domain" in region.keys():
        # Get domain to plot
        domain = region["domain"]
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if "longitude" in kargs:
        lon_west, lon_east, _ = kargs["longitude"]
    if "latitude" in kargs:
        lat_south, lat_north, _ = kargs["latitude"]
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
            raise Exception(
                "Invalid seasonality index={}".format(seasonality_index)
            )
        if seasonality_index == 1:
            color = "black"
        else:
            peak_month = int(gauge[peak_month_export_index])
            color = color_list[peak_month]
        # https://scitools.org.uk/iris/docs/v1.9.2/examples/General/projections_and_annotations.html
        # Place a single marker point for each gauge.
        plt.plot(
            lon,
            lat,
            marker="o",
            color=color,
            markersize=markersize,
            transform=proj_function(),
        )
        # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
        # so for this one we transform the coordinates with a Cartopy call.
        at_x, at_y = ax.projection.transform_point(
            lon, lat, src_crs=proj_function()
        )
    # https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/custom_legends.html
    legend_elements = [
        lines.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="1 <= SI < 2",
            markerfacecolor="black",
            markersize=si_2,
        ),
        lines.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="2 <= SI < 4",
            markerfacecolor="black",
            markersize=si_4,
        ),
        lines.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="4 <= SI < 6",
            markerfacecolor="black",
            markersize=si_6,
        ),
        lines.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="6 <= SI <= 12",
            markerfacecolor="black",
            markersize=si_large,
        ),
    ]
    seasonality_legend_title = "Seasonality (SI)"
    plt.legend(
        handles=legend_elements,
        title=seasonality_legend_title,
        prop={"size": 8},
    )

    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    if title[2] is not None:
        ax.set_title(title[2], loc="right", fontdict=plotSideTitle)
    ax.set_xticks(xticks, crs=proj_function())
    ax.set_yticks(yticks, crs=proj_function())
    lon_formatter = LongitudeFormatter(
        zero_direction_label=True, number_format=".0f"
    )
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Color bar
    cbax = fig.add_axes(
        (
            panel[panel_index][0] + 0.7535,
            panel[panel_index][1] + 0.0515,
            0.0326,
            0.1792,
        )
    )
    # https://matplotlib.org/tutorials/colors/colorbar_only.html
    num_colors = len(color_list)
    if parameter.print_statements:
        print("num_colors={}".format(num_colors))
    cmap = colors.ListedColormap(color_list)
    cbar_label = "Peak month"

    bounds = list(range(num_colors))
    # Set ticks to be in between the bounds
    ticks = list(map(lambda bound: bound + 0.5, bounds))
    # Add one more bound at the bottom of the colorbar.
    # `bounds` should be one longer than `ticks`.
    bounds += [bounds[-1] + 1]
    if parameter.print_statements:
        print("bounds={}".format(bounds))
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cbax,
        boundaries=bounds,
        ticks=ticks,
        spacing="uniform",
        orientation="vertical",
        label=cbar_label,
    )
    # https://matplotlib.org/3.1.1/gallery/ticks_and_spines/colorbar_tick_labelling_demo.html
    months = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    cbar.ax.set_yticklabels(months)
    cbar.ax.invert_yaxis()

    w, h = get_ax_size(fig, cbax)

    cbar.ax.tick_params(labelsize=9.0, length=0)


def plot_seasonality_map(export, parameter):
    if parameter.backend not in ["cartopy", "mpl", "matplotlib"]:
        return

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    panel = [
        (0.0900, 0.5500, 0.7200, 0.3000),
        (0.0900, 0.1300, 0.7200, 0.3000),
    ]

    # Create figure, projection
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=0)

    # test and ref color lists
    # Selected from 'hsv' colormap:
    color_list = [
        (0.05, 0.00, 0.99),
        (0.03, 0.30, 0.98),
        (0.12, 0.92, 0.99),
        (0.13, 1.00, 0.65),
        (0.14, 1.00, 0.05),
        (0.98, 0.99, 0.04),
        (0.99, 0.67, 0.04),
        (0.99, 0.34, 0.03),
        (0.99, 0.07, 0.03),
        (0.99, 0.00, 0.53),
        (0.68, 0.00, 1.00),
        (0.29, 0.00, 1.00),
    ]

    # First panel
    plot_panel_seasonality_map(
        "test", fig, proj, export, color_list, panel, parameter
    )

    # Second panel
    plot_panel_seasonality_map(
        "ref", fig, proj, export, color_list, panel, parameter
    )

    # Figure title
    fig.suptitle(
        parameter.main_title_seasonality_map, x=0.5, y=0.97, fontsize=15
    )

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/streamflow/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print("Output dir: {}".format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/streamflow/{parameter.case_id}
    original_output_dir = get_output_dir(
        parameter.current_set, parameter, ignore_container=True
    )
    if parameter.print_statements:
        print("Original output dir: {}".format(original_output_dir))
    # parameter.output_file_seasonality_map is defined in acme_diags/parameter/streamflow_parameter.py
    # {parameter.results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_seasonality_map}
    file_path = os.path.join(output_dir, parameter.output_file_seasonality_map)
    # {parameter.orig_results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_seasonality_map}
    original_file_path = os.path.join(
        original_output_dir, parameter.output_file_seasonality_map
    )

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        plot_suffix = "." + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        original_plot_file_path = original_file_path + plot_suffix
        # Always print, even without `parameter.print_statements`
        print("Plot saved in: " + original_plot_file_path)

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
            subplot_suffix = ".%i." % i + f
            subplot_file_path = file_path + subplot_suffix
            plt.savefig(subplot_file_path, bbox_inches=extent)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            original_subplot_file_path = original_file_path + subplot_suffix
            # Always print, even without `parameter.print_statements`
            print("Sub-plot saved in: " + original_subplot_file_path)
            i += 1

    plt.close()


def plot_panel_annual_map(
    panel_index, fig, proj, export, bias_array, panel, parameter
):
    if panel_index == 0:
        panel_type = "test"
    elif panel_index == 1:
        panel_type = "ref"
    elif panel_index == 2:
        panel_type = "bias"
    else:
        raise Exception("Invalid panel_index={}".format(panel_index))

    # Plot of streamflow gauges. Color -> peak month, marker size -> seasonality index.

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    ax = fig.add_axes(panel[panel_index], projection=proj)
    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    if "domain" in region.keys():
        # Get domain to plot
        domain = region["domain"]
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if "longitude" in kargs:
        lon_west, lon_east, _ = kargs["longitude"]
    if "latitude" in kargs:
        lat_south, lat_north, _ = kargs["latitude"]
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
    # `export` is the array of gauges. Each gauge has multiple fields -- e.g., lat is index 7
    # Continuous colormap
    colormap = plt.get_cmap("jet_r")
    color_list = list(
        map(lambda index: colormap(index)[:3], range(colormap.N))
    )
    if panel_type in ["test", "ref"]:
        value_min, value_max = 1, 1e4
        # https://matplotlib.org/3.2.1/tutorials/colors/colormapnorms.html
        norm = matplotlib.colors.LogNorm(vmin=value_min, vmax=value_max)
    elif panel_type == "bias":
        if parameter.print_statements:
            value_min = numpy.floor(numpy.min(bias_array))
            value_max = numpy.ceil(numpy.max(bias_array))
            print(
                "Bias of mean annual discharge {} min={}, max={}".format(
                    panel_type, value_min, value_max
                )
            )

        value_min = -100
        value_max = 100
        norm = matplotlib.colors.Normalize()
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))
    for gauge, i in zip(export, range(len(export))):
        if panel_type == "test":
            # Test mean annual discharge
            value = gauge[1]
        elif panel_type == "ref":
            # Ref mean annual discharge
            value = gauge[0]
        elif panel_type == "bias":
            # Bias
            value = bias_array[i]
        else:
            raise Exception("Invalid panel_type={}".format(panel_type))
        if numpy.isnan(value):
            continue
        if value < value_min:
            value = value_min
        elif value > value_max:
            value = value_max
        if panel_type in ["test", "ref"]:
            # Logarithmic Rescale (min-max normalization) to [-1,1] range
            normalized_value = (
                numpy.log10(value) - numpy.log10(value_min)
            ) / (numpy.log10(value_max) - numpy.log10(value_min))
        elif panel_type == "bias":
            # Rescale (min-max normalization) to [-1,1] range
            normalized_value = (value - value_min) / (value_max - value_min)
        else:
            raise Exception("Invalid panel_type={}".format(panel_type))
        lat = gauge[7]
        lon = gauge[8]

        color = color_list[int(normalized_value * (len(color_list) - 1))]
        # https://scitools.org.uk/iris/docs/v1.9.2/examples/General/projections_and_annotations.html
        # Place a single marker point for each gauge.
        plt.plot(
            lon,
            lat,
            marker="o",
            markersize=2,
            color=color,
            transform=proj_function(),
        )
        # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
        # so for this one we transform the coordinates with a Cartopy call.
        at_x, at_y = ax.projection.transform_point(
            lon, lat, src_crs=proj_function()
        )

    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.RIVERS)
    if panel_type == "test":
        title = parameter.test_title
    elif panel_type == "ref":
        title = parameter.reference_title
    elif panel_type == "bias":
        title = "Relative Bias"
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))
    title = (None, title, None)
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    if title[2] is not None:
        ax.set_title(title[2], loc="right", fontdict=plotSideTitle)
    ax.set_xticks(xticks, crs=proj_function())
    ax.set_yticks(yticks, crs=proj_function())
    lon_formatter = LongitudeFormatter(
        zero_direction_label=True, number_format=".0f"
    )
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    # Color bar
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    cbax = fig.add_axes(
        (
            panel[panel_index][0] + 0.6635,
            panel[panel_index][1] + 0.0115,
            0.0326,
            0.1792,
        )
    )
    cmap = colors.ListedColormap(color_list)
    if panel_type in ["test", "ref"]:
        cbar_label = "Mean annual discharge ($m^3$/$s$)"
    elif panel_type == "bias":
        cbar_label = "Bias of mean annual discharge (%)\n(test-ref)/ref"
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        cax=cbax,
        label=cbar_label,
        extend="both",
    )
    w, h = get_ax_size(fig, cbax)
    if panel_type in ["test", "ref"]:
        pass
    elif panel_type == "bias":
        step_size = (value_max - value_min) // 5
        ticks = numpy.arange(
            int(value_min), int(value_max) + step_size, step_size
        )
        cbar.ax.tick_params(labelsize=9.0, length=0)
        cbar.ax.set_yticklabels(ticks)
    else:
        raise Exception("Invalid panel_type={}".format(panel_type))


def plot_annual_map(export, bias, parameter):
    if parameter.backend not in ["cartopy", "mpl", "matplotlib"]:
        return

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    panel = [
        (0.1691, 0.6810, 0.6465, 0.2258),
        (0.1691, 0.3961, 0.6465, 0.2258),
        (0.1691, 0.1112, 0.6465, 0.2258),
    ]

    # Create figure, projection
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=0)

    # First panel
    plot_panel_annual_map(0, fig, proj, export, bias, panel, parameter)

    # Second panel
    plot_panel_annual_map(1, fig, proj, export, bias, panel, parameter)

    # Third panel
    plot_panel_annual_map(2, fig, proj, export, bias, panel, parameter)

    # Figure title
    fig.suptitle(parameter.main_title_annual_map, x=0.5, y=0.97, fontsize=15)

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/streamflow/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print("Output dir: {}".format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/streamflow/{parameter.case_id}
    original_output_dir = get_output_dir(
        parameter.current_set, parameter, ignore_container=True
    )
    if parameter.print_statements:
        print("Original output dir: {}".format(original_output_dir))
    # parameter.output_file_annual_map is defined in acme_diags/parameter/streamflow_parameter.py
    # {parameter.results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_annual_map}
    file_path = os.path.join(output_dir, parameter.output_file_annual_map)
    # {parameter.orig_results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_annual_map}
    original_file_path = os.path.join(
        original_output_dir, parameter.output_file_annual_map
    )

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        plot_suffix = "." + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        original_plot_file_path = original_file_path + plot_suffix
        # Always print, even without `parameter.print_statements`
        print("Plot saved in: " + original_plot_file_path)

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
            subplot_suffix = ".%i." % i + f
            subplot_file_path = file_path + subplot_suffix
            plt.savefig(subplot_file_path, bbox_inches=extent)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            original_subplot_file_path = original_file_path + subplot_suffix
            # Always print, even without `parameter.print_statements`
            print("Sub-plot saved in: " + original_subplot_file_path)
            i += 1

    plt.close()


def plot_annual_scatter(xs, ys, zs, parameter):
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    panel = [(0.0900, 0.2000, 0.7200, 0.6000)]

    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    ax = fig.add_axes(panel[0])
    cmap = plt.get_cmap("jet")
    ax.scatter(xs, ys, label="Scatterplot", marker="o", s=10, c=zs, cmap=cmap)
    r, _ = scipy.stats.pearsonr(xs, ys)
    r2 = r * r
    r2_str = "{0:.2f}".format(r2)
    bounds = [0.01, 100000]
    ax.plot(bounds, bounds, color="red", linestyle="-")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(
        "{} streamflow ($m^3$/$s$)".format(parameter.reference_title),
        fontsize=12,
    )
    ax.set_ylabel(
        "{} streamflow ($m^3$/$s$)".format(parameter.test_title), fontsize=12
    )
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[0], bounds[1])
    ax.tick_params(axis="both", labelsize=12)

    # Color bar
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    cbax = fig.add_axes(
        (panel[0][0] + 0.7535, panel[0][1] + 0.0515, 0.0326, 0.1792 * 2)
    )
    cbar_label = "Drainage area bias (%)"
    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap), cax=cbax)
    cbar.ax.set_ylabel(cbar_label, fontsize=12)
    w, h = get_ax_size(fig, cbax)
    zs_max = numpy.ceil(numpy.max(zs))
    zs_min = numpy.floor(numpy.min(zs))
    step_size = (zs_max - zs_min) // 5
    try:
        ticks = numpy.arange(zs_min, zs_max + step_size, step_size)
        cbar.ax.set_yticklabels(ticks)
    except ValueError:
        # `zs` has invalid values (likely from no area_upstream being found).
        # Just use default colorbar.
        pass
    cbar.ax.tick_params(labelsize=12.0, length=0)

    # Figure title
    if parameter.main_title_annual_scatter == "":
        main_title_annual_scatter = "Annual mean streamflow\n{} vs {}".format(
            parameter.test_title, parameter.reference_title
        )
    else:
        main_title_annual_scatter = parameter.main_title_annual_scatter
    ax.set_title(main_title_annual_scatter, loc="center", y=1.05, fontsize=15)

    legend_title = "$R^2$={}, (n={})".format(r2_str, xs.shape[0])
    ax.legend(
        handles=[], title=legend_title, loc="upper left", prop={"size": 12}
    )

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/streamflow/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print("Output dir: {}".format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/streamflow/{parameter.case_id}
    original_output_dir = get_output_dir(
        parameter.current_set, parameter, ignore_container=True
    )
    if parameter.print_statements:
        print("Original output dir: {}".format(original_output_dir))
    # parameter.output_file_annual_scatter is defined in acme_diags/parameter/streamflow_parameter.py
    # {parameter.results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_annual_scatter}
    file_path = os.path.join(output_dir, parameter.output_file_annual_scatter)
    # {parameter.orig_results_dir}/streamflow/{parameter.case_id}/{parameter.output_file_annual_scatter}
    original_file_path = os.path.join(
        original_output_dir, parameter.output_file_annual_scatter
    )

    # Save figure
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        plot_suffix = "." + f
        plot_file_path = file_path + plot_suffix
        plt.savefig(plot_file_path)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        original_plot_file_path = original_file_path + plot_suffix
        print("Plot saved in: " + original_plot_file_path)

    plt.close()
