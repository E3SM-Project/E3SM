from __future__ import print_function

import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cdutil
import matplotlib
import numpy as np
import numpy.ma as ma
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.colors import hsv_to_rgb

from acme_diags.derivations.default_regions import regions_specs
from acme_diags.driver.utils.general import get_output_dir

matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402


plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}

# Position and sizes of subplot axes in page coordinates (0 to 1)
panel = [
    (0.1691, 0.55, 0.6465, 0.2758),
    (0.1691, 0.15, 0.6465, 0.2758),
]
# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
border = (-0.06, -0.03, 0.13, 0.03)


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
    if degrees_covered >= 270:
        return 60
    if degrees_covered >= 180:
        return 30
    if degrees_covered >= 90:
        return 25
    if degrees_covered >= 60:
        return 20
    elif degrees_covered >= 30:
        return 10
    elif degrees_covered >= 20:
        return 5
    else:
        return 1


def plot_panel(n, fig, proj, var, amp, amp_ref, title, parameter):

    normalize_test_amp = parameter.normalize_test_amp
    lat = var.getLatitude()
    var = ma.squeeze(var.asma())
    max_amp = amp.max()
    max_amp_ref = amp_ref.max()
    amp = ma.squeeze(amp.asma())
    amp_ref = ma.squeeze(amp_ref.asma())

    if normalize_test_amp:
        img = np.dstack(
            (
                (var / 24 - 0.5) % 1,
                (amp / max_amp * (max_amp / max_amp_ref)) ** 0.5,
                np.ones_like(amp),
            )
        )
        max_amp = max_amp_ref
    else:
        img = np.dstack(
            ((var / 24 - 0.5) % 1, (amp / max_amp) ** 0.5, np.ones_like(amp))
        )
    img = hsv_to_rgb(img)

    # imshow plot
    ax = fig.add_axes(panel[n], projection=proj)

    region_str = parameter.regions[0]
    region = regions_specs[region_str]
    global_domain = True
    full_lon = True
    if "domain" in region.keys():  # type: ignore
        # Get domain to plot
        domain = region["domain"]  # type: ignore
        global_domain = False
    else:
        # Assume global domain
        domain = cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
    kargs = domain.components()[0].kargs
    # lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    lon_west, lon_east, lat_south, lat_north = (-180, 180, -90, 90)
    if "longitude" in kargs:
        full_lon = False
        lon_west, lon_east, _ = kargs["longitude"]
        # Note cartopy Problem with gridlines across the dateline:https://github.com/SciTools/cartopy/issues/821. Region cross dateline is not supported yet.
        if lon_west > 180 and lon_east > 180:
            lon_west = lon_west - 360
            lon_east = lon_east - 360

    if "latitude" in kargs:
        lat_south, lat_north, _ = kargs["latitude"]
    lon_covered = lon_east - lon_west
    lon_step = determine_tick_step(lon_covered)
    xticks = np.arange(lon_west, lon_east, lon_step)
    # Subtract 0.50 to get 0 W to show up on the right side of the plot.
    # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the left side of the plot.
    # If a number is added, then the value won't show up at all.
    if global_domain or full_lon:
        xticks = [0, 60, 120, 180, 240, 300, 359.99]
    else:
        xticks = np.append(xticks, lon_east)
        proj = ccrs.PlateCarree()

    lat_covered = lat_north - lat_south
    lat_step = determine_tick_step(lat_covered)
    yticks = np.arange(lat_south, lat_north, lat_step)
    yticks = np.append(yticks, lat_north)

    ax.set_extent([lon_west, lon_east, lat_south, lat_north])  # , crs=proj)

    # Full world would be aspect 360/(2*180) = 1
    # ax.set_aspect((lon_east - lon_west)/(2*(lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    if title[0] is not None:
        ax.set_title(title[0], loc="left", fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    # ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    ax.set_ylim([yticks[0], yticks[-1]])

    # add the image. Because this image was a tif, the "origin" of the image is in the
    # upper left corner
    img_extent = [lon_west, lon_east, lat_south, lat_north]

    # When the requested region is inconsistent with what the data covered (`global` is secified but TRMM only has 50S-5N), set an arbitrary threhold.
    if global_domain and lat_covered - abs(lat[0] - lat[-1]) > 10:
        img_extent = [lon_west, lon_east, lat[0], lat[-1]]
    # ax.imshow(img, origin='lower', extent=img_extent, transform=ccrs.PlateCarree())
    ax.imshow(img, origin="lower", extent=img_extent, transform=proj)
    if region_str == "CONUS":
        ax.coastlines(resolution="50m", color="black", linewidth=0.3)
        state_borders = cfeature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lakes",
            scale="50m",
            facecolor="none",
        )
        ax.add_feature(state_borders, edgecolor="black")

    # Color bar
    bar_ax = fig.add_axes(
        (panel[n][0] + 0.67, panel[n][1] + 0.2, 0.07, 0.07), polar=True
    )
    theta, R = np.meshgrid(np.linspace(0, 2 * np.pi, 24), np.linspace(0, 1, 8))
    H, S = np.meshgrid(np.linspace(0, 1, 24), np.linspace(0, 1, 8))
    image = np.dstack(((H - 0.5) % 1, S ** 0.5, np.ones_like(S)))
    image = hsv_to_rgb(image)
    # bar_ax.set_theta_zero_location('N')
    bar_ax.set_theta_direction(-1)
    bar_ax.set_theta_offset(np.pi / 2)
    bar_ax.set_xticklabels(["0h", "3h", "6h", "9h", "12h", "15h", "18h", "21h"])
    bar_ax.set_yticklabels(["", "", "{:.2f}".format(max_amp)])
    bar_ax.set_rlabel_position(340)
    bar_ax.get_yticklabels()[-2].set_weight("bold")
    # We change the fontsize of minor ticks label
    bar_ax.tick_params(axis="both", labelsize=7, pad=0, length=0)
    bar_ax.text(
        0.2,
        -0.3,
        "Local Time",
        transform=bar_ax.transAxes,
        fontsize=7,
        verticalalignment="center",
    )
    bar_ax.text(
        -0.1,
        -0.9,
        "Max DC amp {:.2f}{}".format(amp.max(), "mm/d"),
        transform=bar_ax.transAxes,
        fontsize=7,
        fontweight="bold",
        verticalalignment="center",
    )
    bar_ax.text(
        -0.1,
        -0.5,
        "DC phase (Hue)",
        transform=bar_ax.transAxes,
        fontsize=7,
        verticalalignment="center",
    )
    bar_ax.text(
        -0.1,
        -0.7,
        "DC amplitude (Saturation)",
        transform=bar_ax.transAxes,
        fontsize=7,
        verticalalignment="center",
    )
    color = image.reshape((image.shape[0] * image.shape[1], image.shape[2]))
    pc = bar_ax.pcolormesh(theta, R, np.zeros_like(R), color=color, shading="auto")
    pc.set_array(None)


def plot(test_tmax, test_amp, ref_tmax, ref_amp, parameter):
    if parameter.backend not in ["cartopy", "mpl", "matplotlib"]:
        return

    # Create figure, projection
    fig = plt.figure(figsize=[8.5, 8.5], dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=180)

    # First panel
    plot_panel(
        0,
        fig,
        proj,
        test_tmax,
        test_amp,
        ref_amp,
        (parameter.test_name_yrs, parameter.test_title),
        parameter,
    )

    # Second panel
    plot_panel(
        1,
        fig,
        proj,
        ref_tmax,
        ref_amp,
        ref_amp,
        (parameter.ref_name_yrs, parameter.reference_title),
        parameter,
    )

    # Figure title
    fig.suptitle(parameter.main_title, x=0.5, y=0.9, fontsize=14)

    # Prepare to save figure
    # get_output_dir => {parameter.results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.results_dir}/enso_diags/{parameter.case_id}
    output_dir = get_output_dir(parameter.current_set, parameter)
    if parameter.print_statements:
        print("Output dir: {}".format(output_dir))
    # get_output_dir => {parameter.orig_results_dir}/{set_name}/{parameter.case_id}
    # => {parameter.orig_results_dir}/enso_diags/{parameter.case_id}
    original_output_dir = get_output_dir(
        parameter.current_set, parameter, ignore_container=True
    )
    if parameter.print_statements:
        print("Original output dir: {}".format(original_output_dir))
    # parameter.output_file is defined in acme_diags/driver/enso_diags_driver.py
    # {parameter.results_dir}/enso_diags/{parameter.case_id}/{parameter.output_file}
    file_path = os.path.join(output_dir, parameter.output_file)
    # {parameter.orig_results_dir}/enso_diags/{parameter.case_id}/{parameter.output_file}
    original_file_path = os.path.join(original_output_dir, parameter.output_file)

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

    # Save individual subplots
    for f in parameter.output_format_subplot:
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
            subplot_suffix = ".%i." % (i) + f
            subplot_file_path = file_path + subplot_suffix
            plt.savefig(subplot_file_path, bbox_inches=extent)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            original_subplot_file_path = original_file_path + subplot_suffix
            print("Sub-plot saved in: " + original_subplot_file_path)
            i += 1

    plt.close()
