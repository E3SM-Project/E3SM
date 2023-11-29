"""
Functionality to create the Taylor diagrams and
metrics table for the Latitude-Longitude set.
"""

import csv
import glob
import os

import matplotlib
import matplotlib.cbook as cbook
import numpy as np
import numpy.ma as ma
from bs4 import BeautifulSoup

import e3sm_diags
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot.cartopy.taylor_diagram import TaylorDiagram

from . import utils

logger = custom_logger(__name__)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

# Variables
VAR_DICT = [
    {
        "name": "Net TOA",
        "units": "W m$^{-2}$",
        "id": "RESTOM global ceres_ebaf_toa_v4.1",
        "exclude": (),
    },
    {
        "name": "SW CRE",
        "units": "W m$^{-2}$",
        "id": "SWCF global ceres_ebaf_toa_v4.1",
        "exclude": (),
    },
    {
        "name": "LW CRE",
        "units": "W m$^{-2}$",
        "id": "LWCF global ceres_ebaf_toa_v4.1",
        "exclude": (),
    },
    {
        "name": "prec",
        "units": "mm day$^{-1}$",
        "id": "PRECT global GPCP_v2.3",
        "exclude": ("CIESM",),
    },
    {"name": "tas land", "units": "K", "id": "TREFHT land ERA5", "exclude": ()},
    {"name": "SLP", "units": "hPa", "id": "PSL global ERA5", "exclude": ()},
    {
        "name": "u-200",
        "units": "m s$^{-1}$",
        "id": "U-200mb global ERA5",
        "exclude": (),
    },
    {
        "name": "u-850",
        "units": "m s$^{-1}$",
        "id": "U-850mb global ERA5",
        "exclude": (),
    },
    {
        "name": "Zg-500",
        "units": "hm",
        "id": "Z3-500mb global ERA5",
        "exclude": ("KIOST-ESM",),
    },
]

CMIP6_EXP = ["historical", "amip"]


def generate_lat_lon_metrics_table(
    lat_lon_table_info, seasons, viewer, root_dir, parameters
):
    """
    For each season in lat_lon_table_info, create a csv,
    convert it to an html and append that html to the viewer.
    """
    set_name = parameters[0].sets[0]
    table_name = ""
    if set_name == "lat_lon_land":
        table_name = "-land"

    if set_name == "lat_lon_river":
        table_name = "-river"

    table_dir = os.path.join(
        root_dir, f"table-data{table_name}"
    )  # output_dir/viewer/table-data

    if not os.path.exists(table_dir):
        os.mkdir(table_dir)

    for season in lat_lon_table_info:
        test_name = (
            parameters[0].short_test_name
            if parameters[0].short_test_name
            else parameters[0].test_name
        )
        if parameters[0].run_type == "model_vs_obs":
            ref_name = "Observation and Reanalysis"
        else:
            ref_name = (
                parameters[0].short_ref_name
                if parameters[0].short_ref_name
                else parameters[0].ref_name
            )
        csv_path = _create_csv_from_dict(
            lat_lon_table_info,
            table_dir,
            season,
            test_name,
            parameters[0].run_type,
        )
        html_path = _cvs_to_html(csv_path, season, test_name, ref_name)

        # Ex: change this: /Users/zshaheen/output_dir/viewer/table-data/ANN_metrics_table.html
        # to this: viewer/table-data/ANN_metrics_table.html
        html_path = "/".join(html_path.split("/")[-3:])

        lat_lon_table_info[season]["html_path"] = html_path

    url = _create_lat_lon_table_index(
        lat_lon_table_info, seasons, viewer, root_dir, table_name
    )
    utils.add_header(root_dir, os.path.join(root_dir, url), parameters)
    _edit_table_html(lat_lon_table_info, seasons, root_dir, table_name)

    return f"Table{table_name}", url


def _create_csv_from_dict(lat_lon_table_info, output_dir, season, test_name, run_type):
    """
    Create a csv for a season in lat_lon_table_info
    in output_dir and return the path to it.
    """
    table_path = os.path.join(output_dir, season + "_metrics_table.csv")

    col_names = [
        "Variables",
        "Unit",
        "Test_mean",
        "Ref._mean",
        "Mean_Bias",
        "Test_STD",
        "Ref._STD",
        "RMSE",
        "Correlation",
    ]

    with open(table_path, "w") as table_csv:
        writer = csv.writer(
            table_csv,
            delimiter=",",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
        )
        writer.writerow(col_names)
        for key, metrics_dic in list(lat_lon_table_info[season].items()):
            metrics = metrics_dic["metrics"]
            if run_type == "model_vs_model":
                key = key.split()[0] + " " + key.split()[1]

            if (
                metrics["test_regrid"]["mean"] == 999.999
                or metrics["ref_regrid"]["mean"] == 999.999
            ):
                mean_bias = 999.999
            else:
                mean_bias = (
                    metrics["test_regrid"]["mean"] - metrics["ref_regrid"]["mean"]
                )
            row = [
                key,
                metrics["unit"],
                round(metrics["test_regrid"]["mean"], 3),
                round(metrics["ref_regrid"]["mean"], 3),
                round(mean_bias, 3),
                round(metrics["test_regrid"]["std"], 3),
                round(metrics["ref_regrid"]["std"], 3),
                round(metrics["misc"]["rmse"], 3),
                round(metrics["misc"]["corr"], 3),
            ]
            writer.writerow(row)

    return table_path


def _cvs_to_html(csv_path, season, test_name, ref_name):
    """
    Convert the csv for a season located at csv_path
    to an HTML, returning the path to the HTML.
    """
    html_path = csv_path.replace("csv", "html")

    with open(html_path, "w") as htmlfile:
        htmlfile.write("<p><b>Test: {}</b><br>".format(test_name))
        htmlfile.write("<b>Reference: {}</b></p>".format(ref_name))
        htmlfile.write("<p><th><b>{} Mean </b></th></p>".format(season))
        htmlfile.write("<table>")

        with open(csv_path) as csv_file:
            read_csv = csv.reader(csv_file)

            # Generate the table's contents.
            for num, row in enumerate(read_csv):
                # Write the header row, assuming the first
                # row in csv contains the header.
                if num == 0:
                    htmlfile.write("<tr>")
                    for column in row:
                        htmlfile.write("<th>{}</th>".format(column))
                    htmlfile.write("</tr>")

                # Write all other rows.
                else:
                    htmlfile.write('<tr><div style="width: 50px">')
                    for column in row:
                        htmlfile.write("<td>{}</td>".format(column))
                    htmlfile.write("</div></tr>")

        htmlfile.write("</table>")

    return html_path


def _create_lat_lon_table_index(
    lat_lon_table_info, seasons, viewer, root_dir, table_name
):
    """
    Create an index in the viewer that links the
    individual htmls for the lat-lon table.
    """
    viewer.add_page(f"Table{table_name}", seasons)
    viewer.add_group("Summary Table")
    viewer.add_row("All variables")

    for s in seasons:
        if s in lat_lon_table_info:
            viewer.add_col(lat_lon_table_info[s]["html_path"], is_file=True, title=s)
        else:
            viewer.add_col("-----", is_file=True, title="-----")

    url = viewer.generate_page()
    return url


# --- Function to read CMIP6 model metrics  ---
def read_cmip6_metrics_from_csv(path, variables, seasons):
    models = []

    with open(path, "r") as fin:
        # skip 3 header lines and last 2 E3SMv2 composites
        rmse = fin.readlines()[3:-2]
        nmodels = len(rmse)
        nvariables = len(variables)
        nseasons = len(seasons)
        data = ma.array(np.zeros((nmodels, nvariables, nseasons)), mask=True)
        for imodel, line in enumerate(rmse):
            models.append(line.split(",")[0])
            model_line = line
            for ivariable in range(nvariables):
                rmse_seasons = model_line.split(",")[
                    1 + ivariable * 5 : 6 + ivariable * 5
                ]
                if any(x == "--" for x in rmse_seasons):
                    pass
                else:
                    data[imodel, ivariable, :] = [
                        float(x)
                        for x in model_line.split(",")[
                            1 + ivariable * 5 : 6 + ivariable * 5
                        ]
                    ]

    # Dictionary to hold data
    d = {}
    d["data"] = data.copy()
    d["models"] = models.copy()
    d["variables"] = variables.copy()
    d["seasons"] = seasons.copy()

    return d


# --- Function to read E3SM Diags metrics  ---
def read_e3sm_diags_metrics(path, variables, seasons, names=None):
    # List of available models
    models = []
    paths = []
    dirs = sorted(glob.glob(path + os.path.sep))
    for d in dirs:
        if not names:
            tmp = d.split(os.path.sep)
            # Note using tmp[-6] for model name only applys to specific e3sm_diags for CMIP6 metrics data structure: i.e.:/global/cfs/cdirs/e3sm/www/CMIP6_comparison_1985-2014_E3SMv2_golaz_etal_2022/*/historical/r1i1p1f1/viewer/table-data on Cori
            model = tmp[-6]
            models.append(model)
        paths.append(d)
    if names:
        models = names

    # Array to hold data
    nmodels = len(models)
    nvariables = len(variables)
    nseasons = len(seasons)
    data = ma.array(np.zeros((nmodels, nvariables, nseasons)), mask=True)

    # Fill data
    for imodel in range(nmodels):
        for iseason in range(nseasons):
            # Open metrics file
            fname = paths[imodel] + "/%s_metrics_table.csv" % (seasons[iseason])
            try:
                with open(fname, "r") as f:
                    content = f.readlines()
                    for ivariable in range(nvariables):
                        # Skip of model has been flagged for this variable
                        if models[imodel] in variables[ivariable]["exclude"]:
                            # print("Excluding: %s, %s, %s" % (models[imodel],variables[ivariable]['name'],seasons[iseason]) )
                            continue
                        lines = [
                            line
                            for line in content
                            if line.startswith(variables[ivariable]["id"])
                        ]
                        if len(lines) > 1:
                            logger.info("Found unexpected multiple entries")
                            pass
                        elif len(lines) == 1:
                            rmse = lines[0].split(",")[-2]
                            if rmse.upper() == "NAN":
                                logger.info(
                                    "NAN: %s, %s, %s"
                                    % (
                                        models[imodel],
                                        variables[ivariable]["name"],
                                        seasons[iseason],
                                    )
                                )
                            else:
                                data[imodel, ivariable, iseason] = float(rmse)
                        else:
                            logger.debug(
                                "Missing: %s, %s, %s"
                                % (
                                    models[imodel],
                                    variables[ivariable]["name"],
                                    seasons[iseason],
                                )
                            )
            except OSError as err:
                logger.debug(f"{err}")

    # Dictionary to hold data
    d = {}
    d["data"] = data.copy()
    d["models"] = models.copy()
    d["variables"] = variables.copy()
    d["seasons"] = seasons.copy()

    return d


def generate_lat_lon_cmip6_comparison(
    lat_lon_table_info, seasons, viewer, root_dir, parameters
):
    """
    For each season in lat_lon_table_info, create a csv, plot a figure for comparing rmse to those from CMIP6 models and append that html to the viewer.
    """
    cmip6_comparison_dir = os.path.join(root_dir, "cmip6-comparison-data")

    if not os.path.exists(cmip6_comparison_dir):
        os.mkdir(cmip6_comparison_dir)

    test_name = parameters[0].test_name_yrs

    # Create CMIP6 comparison figure

    # Seasons
    seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
    nseasons = len(seasons)

    # Read rmse for CMIP6 models (including e3smv1 and e3smv2) for historical r1i1pif1 ensembles averaging over 1985-2014.
    # example root_dir = "/Users/zhang40/Downloads/lat_lon_cmip6_test/viewer"
    test_path = root_dir + "/table-data"
    test_model = read_e3sm_diags_metrics(
        test_path,
        VAR_DICT,
        seasons,
        names=[
            test_name,
        ],
    )

    for exp in CMIP6_EXP:
        control_runs_path = os.path.join(
            e3sm_diags.INSTALL_PATH,
            "control_runs",
            f"cmip6_{exp}_seasonal_rmse_*.csv",
        )
        cmip6_csv_path = sorted(glob.glob(control_runs_path))[-1]
        cmip6_data_access = cmip6_csv_path.split("_")[-1][:6]

        cmip6 = read_cmip6_metrics_from_csv(cmip6_csv_path, VAR_DICT, seasons)

        # Create plot: compare test model with CMIP6, E3SMv1 and v2
        fig = plt.figure(figsize=[24, 18])
        nsx = 4
        nsy = 3
        for ivariable in range(len(VAR_DICT)):
            # CMIP6 data for box and whiskers
            data = []
            labels = []
            for iseason in range(nseasons):
                # Identify model with lowest RMSE
                # ibest = ma.argmin(cmip6["data"][:, ivariable, iseason].compressed())
                # print("Best model %s %s %s" % (variables[ivariable]['name'],seasons[iseason],cmip6['models'][ibest]))
                # Remove missing data using 'compressed()' function
                data.append(cmip6["data"][:, ivariable, iseason].compressed())
                labels.append(seasons[iseason])
            cmip6_stats = cbook.boxplot_stats(data, whis=[0, 100], labels=labels)

            # Plot panel
            ax = plt.subplot(nsy, nsx, ivariable + int(ivariable / 3) + 1)
            ax.set_box_aspect(1)

            # CMIP6 ensemble
            ax.bxp(cmip6_stats)

            # test model
            x = np.arange(nseasons) + 1.0
            ax.scatter(
                x,
                test_model["data"][0, ivariable, :],
                color="k",
                marker="o",
                label=test_name,
                s=60,
            )

            # E3SMv1
            x = np.arange(nseasons) + 0.8
            iE3SMv1 = cmip6["models"].index("E3SM-1-0")
            ax.scatter(
                x,
                cmip6["data"][iE3SMv1, ivariable, :],
                color="b",
                marker=">",
                label=f"E3SMv1 (0101), {exp} (1985-2014)",
                s=60,
            )

            # E3SMv2 (coupled)
            x = np.arange(nseasons) + 1.2
            iE3SMv2 = cmip6["models"].index("E3SM-2-0")
            ax.scatter(
                x,
                cmip6["data"][iE3SMv2, ivariable, :],
                color="r",
                marker="<",
                label=f"E3SMv2 (0101), {exp} (1985-2014)",
                s=60,
            )

            # Customize plot
            ax.set_title("(" + chr(97 + ivariable) + ")", loc="left")
            ax.set_title(
                f"{VAR_DICT[ivariable]['name']} ( {VAR_DICT[ivariable]['units']} )",
                loc="right",
            )
            ax.set_xlim([0.4, nseasons + 0.9])

        fig.subplots_adjust(wspace=0.3, hspace=0.3)

        # Legend base on last subplot
        handles, labels = ax.get_legend_handles_labels()
        ax.text(
            1.2,
            0.1,
            f"Comparison of RMSE (1985-2014) of an ensemble\nof CMIP6 models ({exp} r1i1p1f1 ensemble). \nBox and whiskers show 25th, 75th percentile, \nminimum and maximum RMSE of the ensemble. \nCMIP6 data access: {cmip6_data_access}",
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontsize=20,
        )
        fig.legend(handles, labels, loc=(0.65, 0.8))

        fig.savefig(cmip6_comparison_dir + f"/cmip6_{exp}.png", bbox_inches="tight")
        fig.savefig(cmip6_comparison_dir + f"/cmip6_{exp}.pdf", bbox_inches="tight")

    """
    Create an index in the viewer that links the
    individual htmls for the lat-lon table.
    """
    viewer.add_page("CMIP6 Comparison")
    viewer.add_group("Summary RMSE")
    for exp in CMIP6_EXP:
        viewer.add_row(
            f"RMSE from selected variables and all seasons, compare to CMIP6 {exp}(r1i1p1f1, 1985-2014 average) experiments"
        )
        # We need to make sure we have relative paths for viewers, and not absolute ones.
        pth = f"../viewer/cmip6-comparison-data/cmip6_{exp}.png"

        viewer.add_col(pth, is_file=True, title="output")

    url = viewer.generate_page()
    return "CMIP6 Comparison", url


def generate_lat_lon_taylor_diag(
    lat_lon_table_info, seasons, viewer, root_dir, parameters
):
    """
    For each season in lat_lon_table_info, create a csv, plot using
    taylor diagrams and append that html to the viewer.
    """
    taylor_diag_dir = os.path.join(
        root_dir, "taylor-diagram-data"
    )  # output_dir/viewer/taylor-diagram-data

    if not os.path.exists(taylor_diag_dir):
        os.mkdir(taylor_diag_dir)
    season_to_png = {}  # type: ignore
    for exp in CMIP6_EXP:
        season_to_png[exp] = {}
        for season in lat_lon_table_info:
            test_name = parameters[0].test_name_yrs
            if parameters[0].run_type == "model_vs_obs":
                ref_name = "Observation and Reanalysis"
            else:
                ref_name = (
                    parameters[0].short_ref_name
                    if parameters[0].short_ref_name
                    else parameters[0].ref_name
                )

            csv_path = _create_csv_from_dict_taylor_diag(
                lat_lon_table_info,
                taylor_diag_dir,
                season,
                test_name,
                parameters[0].run_type,
                ref_name,
                exp,
            )

            # Remove any reference to the results_dir when inserting the links into HTML pages.
            # This is because that folder can be renamed.
            csv_path = csv_path.split("viewer")[-1]
            csv_path = "viewer" + csv_path
            season_to_png[exp][season] = csv_path.replace("csv", "png")

    url = _create_taylor_index(seasons, viewer, root_dir, season_to_png)
    utils.add_header(root_dir, os.path.join(root_dir, url), parameters)

    return "Taylor Diagram", url


def _create_csv_from_dict_taylor_diag(
    lat_lon_table_info, output_dir, season, test_name, run_type, ref_name, exp
):
    """
    Create a csv for a season in lat_lon_table_info in
    output_dir and return the path to it.

    Since the Taylor diagram uses the same seasons
    as lat_lon_table_info, we can use that.
    """
    taylor_diag_path = os.path.join(
        output_dir, f"{season}_metrics_taylor_diag_{exp}.csv"
    )

    control_runs_path = os.path.join(
        e3sm_diags.INSTALL_PATH,
        "control_runs",
        f"{season}_metrics_table_taylor_diag_{exp}_1985-2014_E3SMv*.csv",
    )
    base_line_csv_paths = sorted(glob.glob(control_runs_path))

    col_names = ["Variables", "Test_STD", "Ref._STD", "Correlation"]

    with open(taylor_diag_path, "w") as table_csv:
        writer = csv.writer(
            table_csv,
            delimiter=",",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
        )
        writer.writerow(col_names)

        for key, metrics_dic in list(lat_lon_table_info[season].items()):
            # Only include variables from a certain list in the Taylor diagram.
            # An example of the VAR_DICT[ivariable]["id"]: 'RESTOM global ceres_ebaf_toa_v4.1'
            var_ids = []
            var_list = []
            for ivariable in range(len(VAR_DICT)):
                var_ids.append(VAR_DICT[ivariable]["id"].split(" ")[0])  # type: ignore
                var_list.append(VAR_DICT[ivariable]["id"])

            if run_type == "model_vs_obs":
                if key in var_list:
                    metrics = metrics_dic["metrics"]
                    row = [
                        key,
                        round(metrics["test_regrid"]["std"], 3),
                        round(metrics["ref_regrid"]["std"], 3),
                        round(metrics["misc"]["corr"], 3),
                    ]
                    writer.writerow(row)
            else:
                if key.split()[0] in var_ids:
                    metrics = metrics_dic["metrics"]
                    row = [
                        key,
                        round(metrics["test_regrid"]["std"], 3),
                        round(metrics["ref_regrid"]["std"], 3),
                        round(metrics["misc"]["corr"], 3),
                    ]
                    writer.writerow(row)

    with open(taylor_diag_path, "r") as taylor_csv:
        reader = csv.reader(taylor_csv, delimiter=",")
        data = list(reader)
        row_count = len(data)

    # Generate Taylor diagram plot if there is metrics
    # saved for any variable within the list.
    marker = ["o", "d", "+", "s", ">", "<", "v", "^", "x", "h", "X", "H"]
    color = ["k", "b", "r", "y", "m"]

    if row_count > 0:
        matplotlib.rcParams.update({"font.size": 20})
        fig = plt.figure(figsize=(9, 8))
        refstd = 1.0
        taylordiag = TaylorDiagram(refstd, fig=fig, rect=111, label="REF")
        ax = taylordiag._ax

        # Add samples to Taylor diagram.
        for irow in range(1, row_count):
            std_norm, correlation = (
                float(data[irow][1]) / float(data[irow][2]),
                float(data[irow][3]),
            )
            taylordiag.add_sample(
                std_norm,
                correlation,
                marker=marker[irow],
                c=color[0],
                ms=10,
                label=data[irow][0],
                markerfacecolor="None",
                markeredgecolor=color[0],
                linestyle="None",
            )

        # Add a legend to the figure.
        fig.legend(
            taylordiag.samplePoints,
            [p.get_label() for p in taylordiag.samplePoints],
            numpoints=1,
            loc="center right",
            bbox_to_anchor=(1.0, 0.5),
            prop={"size": 10},
        )

        # Add samples for baseline simulation.

        if run_type == "model_vs_obs":
            # Read the control run data.
            # Example base line csv file name: JJA_metrics_table_taylor_diag_historical_1985-2014_E3SMv1.csv
            for ibase, base_line_csv_path in enumerate(base_line_csv_paths):
                base_name = base_line_csv_path.split(".")[-2][-6:]  # ex E3SMv2
                long_base_name = f"{base_name} {exp} (1985-2014)"
                with open(base_line_csv_path, "r") as control_runs_taylor_csv:
                    reader = csv.reader(control_runs_taylor_csv, delimiter=",")
                    control_runs_data = list(reader)

                keys_control_runs = []
                for i in range(0, len(control_runs_data)):
                    keys_control_runs.append(control_runs_data[i][0])

                for irow in range(1, row_count):
                    if data[irow][0] in keys_control_runs:
                        control_irow = keys_control_runs.index(data[irow][0])
                        # std_norm = Test_STD/Ref._STD in the following order
                        # Variables,Unit,Test_mean,Ref._mean,Mean_Bias,Test_STD,Ref._STD,RMSE,Correlation
                        std_norm = float(control_runs_data[control_irow][5]) / float(
                            control_runs_data[control_irow][6]
                        )
                        correlation = float(control_runs_data[control_irow][8])
                        taylordiag.add_sample(
                            std_norm,
                            correlation,
                            marker=marker[irow],
                            c=color[1 + ibase],
                            ms=10,
                            label=data[irow][0] + long_base_name,
                            markerfacecolor="None",
                            markeredgecolor=color[1 + ibase],
                            linestyle="None",
                        )

                    baseline_text = long_base_name
                    ax.text(
                        0.7,
                        0.96 - ibase * 0.05,
                        baseline_text,
                        ha="left",
                        va="center",
                        transform=ax.transAxes,
                        color=color[1 + ibase],
                        fontsize=12,
                    )

        ax.text(
            0.7,
            1.01,
            test_name,
            ha="left",
            va="center",
            transform=ax.transAxes,
            color=color[0],
            fontsize=12,
        )
        if run_type == "model_vs_model":
            ax.text(
                0.6,
                0.95,
                "Ref. Model: " + ref_name,
                ha="left",
                va="center",
                transform=ax.transAxes,
                color="k",
                fontsize=12,
            )

        plt.title(season + ": Spatial Variability", y=1.08)
        fig.savefig(
            os.path.join(output_dir, season + f"_metrics_taylor_diag_{exp}.png")
        )

    return taylor_diag_path


def _create_taylor_index(seasons, viewer, root_dir, season_to_png):
    """
    Create an index in the viewer that links the
    individual htmls for the lat-lon table.
    """
    viewer.add_page("Taylor Diagram", seasons)
    viewer.add_group("Summary Taylor Diagrams")

    for exp in CMIP6_EXP:
        viewer.add_row(f"Selected variables for CMIP6 {exp} comparison")
        for s in seasons:
            if s in season_to_png[exp]:
                pth = os.path.join("..", season_to_png[exp][s])
                viewer.add_col(pth, is_file=True, title=s)
            else:
                viewer.add_col("-----", is_file=True, title="-----")

    url = viewer.generate_page()
    return url


def _edit_table_html(lat_lon_table_info, seasons, root_dir, table_name):
    """
    After the viewer is created, edit the table html to
    insert the custom htmls.
    """

    def _add_html_to_col(season, season_path, html_path):
        """
        Since the output viewer doesn't support html images, do this hack.
        For the col in the html at html_path, insert the link to col_path.
        """
        # Change:
        # <tr class="output-row">
        #  <!-- ... -->
        #  <td colspan="1">
        #   <!-- what needs to be changed -->
        #  </td>
        # <!-- ... -->
        # </tr>
        # to:
        # <tr class="output-row">
        #  <!-- ... -->
        #  <td colspan="1">
        #   <a href="{season_path}"> {season} </a> <!-- this was changed -->
        #  </td>
        # <!-- ... -->
        # </tr>

        soup = BeautifulSoup(open(html_path), "lxml")

        for tr in soup.find_all("tr", {"class": "output-row"}):
            lst = ["All variables"] + seasons
            index = lst.index(season)
            cols = tr.find_all(
                "td"
            )  # the cols are ['All variables', 'ANN', 'DJF', 'MAM', 'JJA', 'SON']
            td = cols[index]  # get the HTML element related to the season

            url = os.path.join("..", "..", season_path)
            a = soup.new_tag("a", href=url)
            a.append(season)

            td.string = ""
            td.append(a)

        html = soup.prettify("utf-8")
        with open(html_path, "wb") as f:
            f.write(html)

    for s in seasons:
        if s in lat_lon_table_info:
            _add_html_to_col(
                s,
                lat_lon_table_info[s]["html_path"],
                os.path.join(root_dir, f"table{table_name}", "index.html"),
            )
