import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from acme_diags.driver.utils.general import get_output_dir

matplotlib.use("agg")
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

plot_info = {}
# Each key gives a list with ax extent, x ticks , y ticks, title, clevs, reference
plot_info["aew"] = [
    [182, 359, 0, 35],
    [240, 300],
    [0, 15, 30],
    "African Easterly Wave Density",
    np.arange(0, 15.1, 1),
    "EAR5 (2000-2014)",
]
plot_info["cyclone"] = [
    [0, 359, -60, 60],
    [0, 60, 120, 180, 240, 300, 359.99],
    [-60, -30, 0, 30, 60],
    "TC Tracks Density",
    np.arange(0, 1, 0.2),
    "IBTrACS (1979-2018)",
]


def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0], lon[0] + 360.0, "coe"))


def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height


def plot_panel(n, fig, proj, var, var_num_years, region, title):

    ax = fig.add_axes(panel[n], projection=proj)
    ax.set_extent(plot_info[region][0], ccrs.PlateCarree())

    clevs = np.arange(0, 15.1, 1)
    clevs = plot_info[region][4]
    p1 = ax.contourf(
        var.getLongitude(),
        var.getLatitude(),
        var / var_num_years,
        transform=ccrs.PlateCarree(),
        levels=clevs,
        extend="both",
        cmap="jet",
    )
    ax.coastlines(lw=0.3)
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor="k")

    if title != "Observation":
        ax.set_title("{}".format(title), fontdict=plotTitle)
    else:
        ax.set_title("{}".format(plot_info[region][5]), fontdict=plotTitle)
    ax.set_xticks(plot_info[region][1], crs=ccrs.PlateCarree())
    ax.set_yticks(plot_info[region][2], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format=".0f")
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction="out", width=1)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    cbax = fig.add_axes((panel[n][0] + 0.6635, panel[n][1] + 0.0215, 0.0326, 0.1792))
    cbar = fig.colorbar(p1, cax=cbax)

    cbar.ax.tick_params(labelsize=9.0, length=0)
    return


def plot_map(test_data, ref_data, region, parameter):
    """Create figure, projection for maps"""

    test = test_data["{}_density".format(region)]
    test_num_years = test_data["{}_num_years".format(region)]

    ref = ref_data["{}_density".format(region)]
    ref_num_years = ref_data["{}_num_years".format(region)]

    fig = plt.figure(figsize=[8.5, 8.5], dpi=parameter.dpi)
    proj = ccrs.PlateCarree(central_longitude=180)

    # First panel
    plot_panel(
        0,
        fig,
        proj,
        test,
        test_num_years,
        region,
        parameter.test_title,
    )

    # Second panel
    plot_panel(
        1,
        fig,
        proj,
        ref,
        ref_num_years,
        region,
        parameter.ref_title,
    )

    # Figure title
    fig.suptitle(plot_info[region][3], x=0.5, y=0.9, fontsize=14)
    # plt.show()
    output_file_name = "{}-density-map".format(region)

    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()


def plot(test, ref, parameter, basin_dict):
    test_metrics = test["metrics"]
    ref_metrics = ref["metrics"]

    test_num_year = test_metrics["num_years"]
    ref_num_year = ref_metrics["num_years"]

    test_name = parameter.test_name
    ref_name = parameter.ref_name

    # TC intensity of each basins
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.15)
    axes = axes.ravel()

    ace_ref = []
    ace_test = []
    num_ref = []
    num_test = []
    for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
        ace_ref.append(ref_metrics[basin][0])
        ace_test.append(test_metrics[basin][0])
        num_ref.append(ref_metrics[basin][3])
        num_test.append(test_metrics[basin][3])
        ax = axes[ifig]
        ax.plot(
            np.arange(1, 7),
            test_metrics[basin][1] / test_num_year,
            "--k",
            linewidth=2,
            label="Test",
        )
        ax.plot(
            np.arange(1, 7),
            ref_metrics[basin][1] / ref_num_year,
            "k",
            linewidth=2,
            label="Ref",
        )
        ax.legend()
        ax.set_title(basin_info[0])
        ax.set_ylabel("TCs per year")
        plt.xticks([1, 2, 3, 4, 5, 6], ["Cat0", "Cat1", "Cat2", "Cat3", "Cat4", "Cat5"])
        ax.xaxis.set_tick_params(labelbottom=True)

    plt.suptitle(
        "Test: {}".format(test_name) + "\n" + "Ref: {}".format(ref_name),
        ha="left",
        x=0.1,
        y=0.99,
    )

    output_file_name = "tc-intensity"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    # TC frequency of each basins
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111)

    N = 6
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27

    ref_vals = num_ref / np.sum(num_ref)
    rects1 = ax.bar(ind + width / 2, ref_vals, width, color="darkgrey")
    test_vals = num_test / np.sum(num_test)
    rects2 = ax.bar(ind - width / 2, test_vals, width, color="black")
    print(
        "total number based on 6 basins",
        np.sum(num_test),
        num_test,
    )

    ax.set_xticks(ind)
    ax.set_xticklabels(
        (
            "North Atlantic",
            "Northwest Pacific",
            "Eastern Pacific",
            "North Indian",
            "South Indian",
            "South Pacific",
        )
    )

    ax.legend(
        (rects2[0], rects1[0]),
        (
            "{}: (~{})storms per year".format(test_name, int(np.sum(num_test))),
            "{}: (~{}) storms per year".format(ref_name, int(np.sum(num_ref))),
        ),
    )
    ax.set_ylabel("Fraction")
    ax.set_title("Relative frequency of TCs for each ocean basins")

    output_file_name = "tc-frequency"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    fig1 = plt.figure(figsize=(12, 6))
    ax = fig1.add_subplot(111)

    N = 6
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27

    ref_vals = ace_ref / np.sum(ace_ref)
    rects2 = ax.bar(ind - width / 2, ref_vals, width, color="black")
    test_vals = ace_test / np.sum(ace_test)
    rects1 = ax.bar(ind + width / 2, test_vals, width, color="darkgrey")

    ax.set_xticks(ind)
    ax.set_xticklabels(
        (
            "North Atlantic",
            "Northwest Pacific",
            "Eastern Pacific",
            "North Indian",
            "South Indian",
            "South Pacific",
        )
    )

    ax.legend((rects2[0], rects1[0]), (ref_name, test_name))
    ax.set_ylabel("Fraction")
    ax.set_title(
        "Distribution of accumulated cyclone energy (ACE) among various ocean basins"
    )
    output_file_name = "ace-distribution"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.4, wspace=0.15)
    axes = axes.ravel()

    for ifig, (basin, basin_info) in enumerate(basin_dict.items()):
        ax = axes[ifig]
        ax.plot(np.arange(1, 13), ref_metrics[basin][2], "k", linewidth=2, label="Test")
        ax.plot(
            np.arange(1, 13),
            test_metrics[basin][2],
            "--k",
            linewidth=2,
            label="Ref",
        )
        ax.legend()
        ax.set_title(basin_info[0])
        ax.set_ylabel("Fraction")
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"],
        )
        ax.xaxis.set_tick_params(labelbottom=True)

    plt.suptitle(
        "Test: {}".format(test_name) + "\n" + "Ref: {}".format(ref_name),
        ha="left",
        x=0.1,
        y=0.99,
    )

    output_file_name = "tc-frequency-annual-cycle"
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            output_file_name + "." + f,
        )
        plt.savefig(fnm)
        print("Plot saved in: " + fnm)
        plt.close()

    ##########################################################
    # Plot TC tracks density
    plot_map(test, ref, "aew", parameter)

    # Plot AEW density
    plot_map(test, ref, "cyclone", parameter)
