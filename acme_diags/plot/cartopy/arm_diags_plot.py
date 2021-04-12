import math
import os

import matplotlib
import numpy as np

from acme_diags.driver.utils.diurnal_cycle import fastAllGridFT
from acme_diags.driver.utils.general import get_output_dir

matplotlib.use("agg")

import matplotlib.pyplot as plt  # isort:skip  # noqa: E402


def plot_convection_onset_statistics(
    test_pr, test_prw, ref_pr, ref_prw, parameter, region
):
    # Original code: Kathleen Schiro, python version 22 Dec 2016, University of California Dept. of Atmospheric and Oceanic Sciences
    # Modifications: Baird Langenbrunner, Yi-Hung Kuo
    # Modifications: Jill Zhang, Cheng Tao
    # Scientific supervision: Prof. J David Neelin
    #
    # For related publications and research information see
    # the Neelin group webpage  http://www.atmos.ucla.edu/~csi/ #

    # Define precip threshold
    precip_threshold = 0.5  # default 0.5 (in mm/hr)

    # Define cwc bounds and bin_width for each site
    if region == "twpc1":  # twpc1
        cwv_max = 69
        cwv_min = 28
        bin_width = 1.5
        sitename = "Manus Island"
    if region == "twpc2":  # twpc2
        cwv_max = 70
        cwv_min = 28
        bin_width = 2.0
        sitename = "Nauru"
    if region == "twpc3":  # twpc3
        cwv_max = 85
        cwv_min = 28
        bin_width = 2.0
        sitename = "Darwin"
    if region == "sgp":  # sgp
        cwv_max = 75
        cwv_min = 20
        bin_width = 2.0
        sitename = "SGP"

    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    fig.subplots_adjust(wspace=0.3)
    title = ""
    for index in range(2):
        if index == 0:
            precip = test_pr
            cwv = test_prw
            data_name = "Test: " + parameter.test_name
            time_interval = 3
            line_color = ["black", "grey"]
        else:
            precip = ref_pr
            cwv = ref_prw
            data_name = "Ref: " + parameter.ref_name
            time_interval = 1
            line_color = ["blue", "steelblue"]
        var_time_absolute = cwv.getTime().asComponentTime()
        time_interval = int(var_time_absolute[1].hour - var_time_absolute[0].hour)

        number_of_bins = int(np.ceil((cwv_max - cwv_min) / bin_width))
        bin_center = np.arange(
            (cwv_min + (bin_width / 2)),
            (cwv_max - (bin_width / 2)) + bin_width,
            bin_width,
        )
        if len(bin_center) != number_of_bins:
            bin_center = np.arange(
                (cwv_min + (bin_width / 2)), (cwv_max - (bin_width / 2)), bin_width
            )

        # Define variables for binning
        bin_index = np.zeros([number_of_bins, cwv.size])
        precip_binned = np.empty([number_of_bins, cwv.size]) * np.nan
        precip_counts = np.zeros([number_of_bins, cwv.size])

        np.warnings.filterwarnings("ignore")
        # Bin the data by CWV value as specified above
        for i in range(0, number_of_bins):
            tmp1 = np.where(cwv > cwv_min + (i * bin_width))
            bin_index[i, tmp1] = 1
            tmp2 = np.where(cwv > cwv_min + (i * bin_width) + bin_width)
            bin_index[i, tmp2] = 0

        for i in range(0, number_of_bins):
            tmp1 = np.where(bin_index[i, :] == 1)
            precip_binned[i, tmp1] = precip[tmp1]
            tmp2 = np.where(bin_index[i, :] != 1)
            precip_binned[i, tmp2] = np.nan

        for i in range(0, number_of_bins):
            tmp1 = np.where(precip_binned[i, :] >= precip_threshold)
            precip_counts[i, tmp1] = 1
            for j in range(0, cwv.size):
                if np.isnan(precip_binned[i, j]):
                    precip_counts[i, j] = np.nan

        # Create binned arrays
        hist_cwv = np.empty([number_of_bins, 1]) * np.nan
        hist_precip_points = np.empty([number_of_bins, 1]) * np.nan
        pr_binned_mean = np.empty([number_of_bins, 1]) * np.nan
        pr_binned_std = np.empty([number_of_bins, 1]) * np.nan
        pr_probability = np.empty([number_of_bins, 1]) * np.nan
        errorbar_precip_points = np.empty([number_of_bins, 1]) * np.nan
        errorbar_precip = np.empty([number_of_bins, 1]) * np.nan

        ###
        errorbar_precip_binom = np.empty([number_of_bins, 2]) * np.nan

        # Fill binned arrays
        hist_cwv = bin_index.sum(axis=1)
        hist_cwv[hist_cwv <= 1] = 0
        hist_precip_points = np.nansum(precip_counts, axis=1)
        hist_precip_points[hist_precip_points <= 1] = 0
        pr_binned_mean = np.nanmean(precip_binned, axis=1)
        pr_binned_std = np.nanstd(precip_binned, axis=1)
        r = np.empty([1, number_of_bins]) * np.nan
        r = np.sum(~np.isnan(precip_counts), axis=1)
        pr_probability = np.nansum(precip_counts, axis=1) / r
        freq_cwv = (hist_cwv / bin_width) / np.nansum(hist_cwv)
        freq_precipitating_points = hist_precip_points / bin_width / np.nansum(hist_cwv)

        for i in range(0, number_of_bins):
            errorbar_precip[i] = pr_binned_std[i] / math.sqrt(hist_cwv[i])
            errorbar_precip_points[i] = (
                math.sqrt(hist_precip_points[i])
                / np.nansum(hist_cwv / bin_width)
                / bin_width
            )
            z = 0.675
            phat = hist_precip_points[i] / hist_cwv[i]
            errorbar_precip_binom[i, 0] = z * math.sqrt(phat * (1 - phat) / hist_cwv[i])
            errorbar_precip_binom[i, 1] = z * math.sqrt(phat * (1 - phat) / hist_cwv[i])
        axes_fontsize = 12  # size of font in all plots
        legend_fontsize = 9
        marker_size = 40  # size of markers in scatter plots
        xtick_pad = 10  # padding between x tick labels and actual plot
        bin_width = (np.max(bin_center) - np.min(bin_center)) / number_of_bins

        # create figure canvas

        # create figure 1
        ax1 = axes[0]
        xulim = 5 * np.ceil(np.max(np.round(bin_center + bin_width / 2)) / 5)
        xllim = 5 * np.floor(np.min(np.round(bin_center - bin_width / 2)) / 5)
        # ax1.set_xlim(xllim-10,xulim+15)
        # ax1.set_xticks(np.arange(np.ceil(xllim/10)*10-10,np.ceil(xulim/10)*10+15,15))
        # ax1.set_yticks(np.arange(0,5))
        ax1.tick_params(labelsize=axes_fontsize)
        ax1.tick_params(axis="x", pad=10)
        ax1.errorbar(
            bin_center,
            pr_binned_mean,
            xerr=0,
            yerr=errorbar_precip.squeeze(),
            ls="none",
            color="black",
        )
        # ax1.scatter(bin_center, pr_binned_mean, edgecolor='none', facecolor=scatter_colors, s=marker_size, clip_on=False, zorder=3)
        ax1.scatter(
            bin_center,
            pr_binned_mean,
            edgecolor="none",
            facecolor=line_color[0],
            s=marker_size,
            clip_on=True,
            zorder=3,
            label=data_name.split(":")[0],
        )
        ax1.set_xlim(xllim - 10, cwv_max)
        ax1.set_ylim(0, 3)
        ax1.set_ylabel("Precip (mm/hr)", fontsize=axes_fontsize)
        ax1.set_xlabel("CWV (mm)", fontsize=axes_fontsize)
        ax1.set_axisbelow(True)
        legend_handles, legend_labels = ax1.get_legend_handles_labels()
        ax1.legend(legend_handles, legend_labels, loc="upper left", frameon=False)

        # create figure 2 (probability pickup)
        ax2 = axes[1]
        xulim = 5 * np.ceil(np.max(np.round(bin_center + bin_width / 2)) / 5)
        xllim = 5 * np.floor(np.min(np.round(bin_center - bin_width / 2)) / 5)
        # ax2.set_xlim(xllim-10,xulim+15)
        ax2.tick_params(labelsize=axes_fontsize)
        ax2.errorbar(
            bin_center,
            pr_probability,
            xerr=0,
            yerr=errorbar_precip_binom.T,
            fmt="none",
            color="black",
        )
        ax2.tick_params(axis="x", pad=xtick_pad)
        ax2.scatter(
            bin_center,
            pr_probability,
            marker="d",
            s=marker_size,
            edgecolor="none",
            facecolor=line_color[0],
            zorder=3,
            label=data_name.split(":")[0],
        )
        ax2.set_xlim(xllim - 10, cwv_max)
        # ax2.set_xticks(np.arange(np.ceil(xllim/10)*10-10,np.ceil(xulim/10)*10+15,15))
        ax2.set_ylim(0, 1)
        ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax2.set_ylabel("Probability of Precip.", fontsize=axes_fontsize)
        ax2.set_xlabel("CWV (mm)", fontsize=axes_fontsize)
        # ax2.grid()
        ax2.set_axisbelow(True)
        legend_handles, legend_labels = ax2.get_legend_handles_labels()
        ax2.legend(legend_handles, legend_labels, loc="upper left", frameon=False)
        title = (
            title
            + data_name
            + ":  "
            + str(time_interval)
            + " hrly("
            + line_color[0]
            + ")\n"
        )

        # create figure 3 (non-normalized PDF)
        ax3 = axes[2]
        ax3.set_yscale("log")

        xulim = 5 * np.ceil(np.max(np.round(bin_center + bin_width / 2)) / 5)
        xllim = 5 * np.floor(np.min(np.round(bin_center - bin_width / 2)) / 5)
        # ax3.set_xlim(xllim-10,xulim+15)
        ax3.set_xlim(xllim - 10, cwv_max)
        ax3.set_xticks(
            np.arange(np.ceil(xllim / 10) * 10 - 10, np.ceil(xulim / 10) * 10 + 15, 15)
        )
        # low_lim = -6.0
        low_lim = -4.0
        ax3.set_ylim(10 ** low_lim, 100)
        ax3.set_yticks(10 ** np.arange(low_lim, 2, dtype="float64"))
        ax3.tick_params(labelsize=axes_fontsize)
        ax3.tick_params(axis="x", pad=xtick_pad)
        freq_precipitating_points[freq_precipitating_points == 0] = np.nan
        freq_cwv[freq_cwv == 0] = np.nan

        # ax3.errorbar(bin_center, freq_precipitating_points, xerr=0, yerr=errorbar_precip_points.squeeze(), ls='none', color='black')
        ax3.scatter(
            bin_center,
            freq_cwv,
            color=line_color[0],
            label=data_name.split(":")[0] + ": all",
        )
        ax3.scatter(
            bin_center,
            freq_precipitating_points,
            edgecolor="none",
            facecolor=line_color[1],
            s=marker_size,
            zorder=3,
            label=data_name.split(":")[0] + ": precip $>$ 0.5 mm/hr ",
        )
        ax3.set_ylabel("PDF", fontsize=axes_fontsize)
        ax3.set_xlabel("CWV (mm)", fontsize=axes_fontsize)
        ax3.set_axisbelow(True)

        # create legend
        legend_handles, legend_labels = ax3.get_legend_handles_labels()
        ax3.legend(
            legend_handles,
            legend_labels,
            loc="upper left",
            bbox_to_anchor=(0.1, 0.95),
            fontsize=legend_fontsize,
            scatterpoints=1,
            handlelength=0,
            labelspacing=0,
            borderpad=0,
            borderaxespad=0,
            frameon=False,
        )

    # set layout to tight (so that space between figures is minimized)
    # plt.tight_layout()
    plt.suptitle(
        "Convection Onset Metrics" + " at " + sitename, y=1.15, fontweight="bold"
    )
    plt.title(title, ha="left", x=-2, y=0.98)

    # save figure
    # mp.savefig(output_path +'/figures/conv_diagnostics_'+test+'_'+sites[0]+'.png', transparent=True, bbox_inches='tight')
    # Save the figure.
    output_file_name = parameter.output_file
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter), output_file_name + "." + f
        )
        plt.savefig(fnm, transparent=True, bbox_inches="tight")
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter, ignore_container=True),
            output_file_name + "." + f,
        )
        print("Plot saved in: " + fnm)

    plt.close()


def get_seasonal_mean(data):
    """
    Calculate annual mean and seasonal mean of input data (mean of 12 month)
    """
    ac = data.asma()
    ac_dec = np.concatenate((ac, ac))[11:23]
    season = np.nanmean(ac_dec.reshape(-1, 3), axis=1)
    ann = np.nanmean(season)
    return np.hstack([ann, season.data])


def plot_annual_cycle(var, vars_to_data, parameter):
    line_color = ["r", "b", "g", "m"]
    fig = plt.figure()  # Create figure

    ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8])  # Create axes
    xax = np.arange(1, 13, 1)

    refs = vars_to_data.refs
    test = vars_to_data.test
    ax1.plot(
        xax, test.asma(), "k", linewidth=2, label="Test: " + parameter.test_name
    )  # +' ({0:.1f})'.format(np.mean(test.asma())))
    test_season = get_seasonal_mean(test)
    for i_ref, ref in enumerate(refs):
        ref_season = get_seasonal_mean(ref)
        ax1.plot(
            xax,
            ref.asma(),
            line_color[i_ref],
            linewidth=2,
            label="Ref: " + parameter.ref_name,
        )  # +' ({0:.1f})'.format(np.mean(ref.asma())))
    my_xticks = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
    plt.xticks(xax, my_xticks)
    plt.xlim(1, 12)
    ymin, ymax = plt.gca().get_ylim()
    plt.ylim(0.8 * ymin, 1.2 * ymax)
    #     plt.ylim(ylim[va_ind])
    plt.xlabel("Month")
    plt.legend(loc="best", prop={"size": 10})

    if var == "PRECT":
        plt.ylabel("Total Precipitation Rate" + " (" + parameter.var_units + ")")
    else:
        plt.ylabel(parameter.var_name + " (" + parameter.var_units + ")")

    # Add a table at the bottom of the axes
    bias = test_season - ref_season
    cell_text = np.round(np.vstack((test_season, ref_season, bias)), 2)
    collabel = ("ANN", "DJF", "MAM", "JJA", "SON")
    rows = ("Test", "Ref", "Bias")
    plt.table(
        cellText=cell_text,
        rowLabels=rows,
        colLabels=collabel,
        alpha=0.8,
        bbox=[0.15, 0.0, 0.8, 0.2],
    )

    # Save the figure.
    output_file_name = parameter.output_file
    plt.title(output_file_name.replace("-", " "))

    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter), output_file_name + "." + f
        )
        plt.savefig(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter, ignore_container=True),
            output_file_name + "." + f,
        )
        print("Plot saved in: " + fnm)

    plt.close()


def plot_diurnal_cycle(var, vars_to_data, parameter):
    test = vars_to_data.test[0]
    ref = vars_to_data.refs[0][0]
    lst = vars_to_data.misc[0]
    t_conv = lst[0][0]

    output_file_name = parameter.output_file + "-" + "diurnal-cycle"

    fig = plt.figure()  # Create figure
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])  # Create axes

    for index in range(2):
        if index == 0:
            data = test
            line_c = "k"
            data_name = parameter.test_name
        else:
            data = ref
            line_c = "r"
            data_name = parameter.ref_name

        time_freq = len(data)
        res = int(24 / time_freq)
        c, maxvalue, tmax = fastAllGridFT(data, [0])
        xax = np.linspace(0, 48 - res, time_freq * 2)
        ax.plot(xax, np.concatenate((data, data)), "." + line_c, label=data_name)
        xax = np.linspace(0, 48 - res, time_freq * 2 * 3)
        w = 2.0 * np.pi / 24
        yax = (c + maxvalue[0] * np.sin(w * xax + np.pi / 2 - tmax[0] * w))[0]
        ax.plot(xax, yax, line_c, label="First harmonic")
        plt.xlim([24 - t_conv, 47 - t_conv + 1])
        plt.ylim([0, 5.5])
        # ymin, ymax = plt.gca().get_ylim()
        # plt.ylim(ymin, ymax)
        plt.xlabel("local solar time [hr]")
        # plt.ylabel(parameter.var_name + ' (' +parameter.var_units+ ')')
        plt.ylabel("Total Precipitation Rate" + " (" + parameter.var_units + ")")
        xax = np.arange(24 - t_conv, 47 - t_conv, 3)
        my_xticks = ["0h", "3h", "6h", "9h", "12h", "15h", "18h", "21h"]
        plt.xticks(xax, my_xticks)
        plt.legend(loc="upper right")
        plt.title(output_file_name.replace("-", " "))

    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter), output_file_name + "." + f
        )
        plt.savefig(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter, ignore_container=True),
            output_file_name + "." + f,
        )
        print("Plot saved in: " + fnm)

    plt.close()


def plot_diurnal_cycle_zt(var, vars_to_data, parameter):
    ref = vars_to_data.refs[0]
    test = vars_to_data.test
    lst = vars_to_data.misc
    month = [
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

    for index in range(2):
        fig, axs = plt.subplots(
            4,
            3,
            figsize=(15, 12),
            facecolor="w",
            edgecolor="k",
            sharex=True,
            sharey=True,
        )
        fig.subplots_adjust(hspace=0.4, wspace=0.1)
        axs = axs.ravel()
        t_conv = lst[0][0][0]
        for imon in range(12):
            if index == 0:
                title = parameter.ref_name
                data = ref
                data_name = "ref"
                if "armdiags" in title:
                    data = data[:, :, ::-1]

            else:
                title = parameter.test_name
                data = test
                data_name = "test"

            time_freq = data.shape[1]
            yy = np.linspace(0, 48, time_freq * 2)
            xx = np.linspace(100, 1000, 37)
            x, y = np.meshgrid(xx, yy)
            data_con = np.concatenate((data[imon, :, :], data[imon, :, :]), axis=0)
            im = axs[imon].pcolormesh(
                y, x, data_con[:, :], vmin=0, vmax=30, cmap="jet", shading="auto"
            )
            axs[imon].set_title(month[imon])
            plt.xlim([24 - t_conv, 47 - t_conv])
            xax = np.arange(24 - t_conv, 47 - t_conv, 3)
            my_xticks = ["0", "3", "6", "9", "12", "15", "18", "21"]
            plt.xticks(xax, my_xticks)
            axs[imon].xaxis.set_tick_params(labelbottom=True)

            # for ax in axs[9:12]:
            axs[imon].set_xlabel("Local time (hr)")
        for ax in axs[::3]:
            ax.set_ylabel("Pressure (mb)")
        axs[0].invert_yaxis()
        site = parameter.output_file.split("-")[-1]
        suptitle = "Cloud Fraction Monthly Diurnal Cycle " + site + "\n" + title
        plt.suptitle(suptitle, fontsize=20)
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        plt.title("cl (%)")

        output_file_name = parameter.output_file + "-" + data_name
        for f in parameter.output_format:
            f = f.lower().split(".")[-1]
            fnm = os.path.join(
                get_output_dir(parameter.current_set, parameter),
                output_file_name + "." + f,
            )
            plt.savefig(fnm)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            fnm = os.path.join(
                get_output_dir(parameter.current_set, parameter, ignore_container=True),
                output_file_name + "." + f,
            )
            print("Plot saved in: " + fnm)

        plt.close()
