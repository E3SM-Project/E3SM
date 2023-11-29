import os

import cartopy.crs as ccrs
import matplotlib
import numpy as np

from e3sm_diags.driver.utils.general import get_output_dir
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics import mean
from e3sm_diags.plot.cartopy.lat_lon_plot import plot_panel

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

plotTitle = {"fontsize": 11.5}
plotSideTitle = {"fontsize": 9.5}


def plot(test, test_site, ref_site, parameter):
    # Plot scatter plot
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # (left, bottom, width, height) in page coordinates
    panel = [
        (0.09, 0.40, 0.72, 0.30),
        (0.19, 0.2, 0.62, 0.30),
    ]
    # Border padding relative to subplot axes for saving individual panels
    # (left, bottom, right, top) in page coordinates
    border = (-0.06, -0.03, 0.13, 0.03)

    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.var_id, x=0.5, y=0.97)
    proj = ccrs.PlateCarree()
    max1 = test.max()
    min1 = test.min()
    mean1 = mean(test)
    plot_panel(
        0,
        fig,
        proj,
        test,
        parameter.contour_levels,
        parameter.test_colormap,
        (parameter.test_name_yrs, None, None),
        parameter,
        stats=(max1, mean1, min1),
    )

    ax = fig.add_axes(panel[1])
    ax.set_title(f"{parameter.var_id} from AERONET sites")

    # define 1:1 line, and x y axis limits

    if parameter.var_id == "AODVIS":
        x1 = np.arange(0.01, 3.0, 0.1)
        y1 = np.arange(0.01, 3.0, 0.1)
        plt.xlim(0.03, 1)
        plt.ylim(0.03, 1)
    else:
        x1 = np.arange(0.0001, 1.0, 0.01)
        y1 = np.arange(0.0001, 1.0, 0.01)
        plt.xlim(0.001, 0.3)
        plt.ylim(0.001, 0.3)

    plt.loglog(x1, y1, "-k", linewidth=0.5)
    plt.loglog(x1, y1 * 0.5, "--k", linewidth=0.5)
    plt.loglog(x1 * 0.5, y1, "--k", linewidth=0.5)

    corr = np.corrcoef(ref_site, test_site)
    xmean = np.mean(ref_site)
    ymean = np.mean(test_site)
    ax.text(
        0.3,
        0.9,
        f"Mean (test): {ymean:.3f} \n Mean (ref): {xmean:.3f}\n Corr: {corr[0, 1]:.2f}",
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
    )

    # axis ticks
    plt.tick_params(axis="both", which="major")
    plt.tick_params(axis="both", which="minor")

    # axis labels
    plt.xlabel(f"ref: {parameter.ref_name_yrs}")
    plt.ylabel(f"test: {parameter.test_name_yrs}")

    plt.loglog(ref_site, test_site, "kx", markersize=3.0, mfc="none")

    # legend
    plt.legend(frameon=False, prop={"size": 5})

    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            f"{parameter.output_file}" + "." + f,
        )
        plt.savefig(fnm)
        logger.info(f"Plot saved in: {fnm}")

    for f in parameter.output_format_subplot:
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            parameter.output_file,
        )
        page = fig.get_size_inches()
        i = 0
        for p in panel:
            # Extent of subplot
            subpage = np.array(p).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(border).reshape(2, 2)
            subpage = list(((subpage) * page).flatten())  # type: ignore
            extent = matplotlib.transforms.Bbox.from_extents(*subpage)
            # Save subplot
            fname = fnm + ".%i." % (i) + f
            plt.savefig(fname, bbox_inches=extent)

            orig_fnm = os.path.join(
                get_output_dir(parameter.current_set, parameter),
                parameter.output_file,
            )
            fname = orig_fnm + ".%i." % (i) + f
            logger.info(f"Sub-plot saved in: {fname}")

            i += 1
