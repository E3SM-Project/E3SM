#!/usr/bin/env python
# coding=utf-8
# Copyright (c) 2018-2022 UT-BATTELLE, LLC
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""The Kolmogorov-Smirnov Test:
This tests the null hypothesis that the reference (n) and modified (m) model
Short Independent Simulation Ensembles (SISE) represent the same climate
state, based on the equality of distribution of each variable's annual global
average in the standard monthly model output between the two simulations.

The (per variable) null hypothesis uses the non-parametric, two-sample (n and m)
Kolmogorov-Smirnov test as the univariate test of equality of distribution of
global means. The test statistic (t) is the number of variables that reject the
(per variable) null hypothesis of equality of distribution at a 95% confidence
level. The (overall) null hypothesis is rejected if t > α, where α is some
critical number of rejecting variables. The critical value, α, is obtained from
an empirically derived approximate null distribution of t using resampling
techniques.
"""

import argparse
import glob
import os
import re
from collections import OrderedDict
from functools import partial
from pathlib import Path
from pprint import pprint

import livvkit
import numpy as np
import pandas as pd
import six
from livvkit import elements as el
from livvkit.util import functions as fn
from livvkit.util.LIVVDict import LIVVDict
from scipy import stats
from statsmodels.stats import multitest as smm

from evv4esm import EVVException, human_color_names
# from evv4esm.ensembles import e3sm
from evv4esm.ensembles.tools import monthly_to_annual_avg, prob_plot
from evv4esm.utils import bib2html

from netCDF4 import Dataset


def component_file_instance(component, case_file):
    search_regex = r"{c}_[0-9]+".format(c=component)
    result = re.search(search_regex, case_file).group(0)
    return int(result.replace("{}_".format(component), ""))


def file_date_str(case_file, style="short", hist_name="h"):
    if style == "full":
        search_regex = r"{}\.[0-9]+-[0-9]+-[0-9]+-[0-9]+.nc".format(hist_name)
    elif style == "med":
        search_regex = r"{}\.[0-9]+-[0-9]+-[0-9]+.nc".format(hist_name)
    elif style == "short":
        search_regex = r"{}\.[0-9]+-[0-9]+.nc".format(hist_name)
    else:
        search_regex = r"{}\.[0-9]+-[0-9]+.nc".format(hist_name)

    result = re.search(search_regex, case_file).group(0)
    return result.replace("{}.".format(hist_name), "").replace(".nc", "")


def component_monthly_files(
    dir_, component, ninst, hist_name="h0", nmonth_max=12, date_style="short"
):
    if date_style == "full":
        date_search = "????-??-??-??"
    elif date_style == "med":
        date_search = "????-??-??"
    else:
        date_search = "????-??"

    base = "{d}/*{c}_????.{n}.{ds}.nc".format(
        d=dir_, c=component, n=hist_name, ds=date_search
    )
    search = os.path.normpath(base)
    result = sorted(glob.glob(search))

    instance_files = OrderedDict()
    _file_date_str = partial(file_date_str, style=date_style, hist_name=hist_name)
    for ii in range(1, ninst + 1):
        instance_files[ii] = sorted(
            filter(lambda x: component_file_instance(component, x) == ii, result),
            key=_file_date_str,
        )
        if len(instance_files[ii]) > nmonth_max:
            instance_files[ii] = instance_files[ii][-nmonth_max:]

    return instance_files


def get_variable_meta(dataset, var_name):
    try:
        _name = f": {dataset.variables[var_name].getncattr('long_name')}"
    except AttributeError:
        _name = ""
    try:
        _units = f" [{dataset.variables[var_name].getncattr('units')}]"
    except AttributeError:
        _units = ""
    return {"long_name": _name, "units": _units}


def gather_monthly_averages(ensemble_files, variable_set=None):
    monthly_avgs = []
    for case, inst_dict in six.iteritems(ensemble_files):
        for inst, i_files in six.iteritems(inst_dict):
            # Get monthly averages from files
            for file_ in i_files:
                date_str = file_date_str(file_)

                data = None
                try:
                    data = Dataset(file_)
                    if variable_set is None:
                        variable_set = set(data.variables.keys())
                except OSError as E:
                    six.raise_from(
                        BaseException(
                            "Could not open netCDF dataset: {}".format(file_)
                        ),
                        E,
                    )

                for var in data.variables.keys():
                    if var not in variable_set:
                        continue
                    if len(data.variables[var].shape) < 2 or var in [
                        "time_bnds",
                        "date_written",
                        "time_written",
                    ]:
                        continue
                    elif "ncol" not in data.variables[var].dimensions:
                        continue
                    else:
                        m = np.mean(data.variables[var][0, ...])

                    desc = "{long_name}{units}".format(**get_variable_meta(data, var))
                    monthly_avgs.append(
                        (case, var, "{:04}".format(inst), date_str, m, desc)
                    )

    monthly_avgs = pd.DataFrame(
        monthly_avgs,
        columns=("case", "variable", "instance", "date", "monthly_mean", "desc"),
    )
    return monthly_avgs


def load_mpas_climatology_ensemble(files, field_name, mask_value=None):
    # Get the first file to set up ensemble array output
    with Dataset(files[0], "r") as dset:
        _field = dset.variables[field_name][:].squeeze()
        var_desc = "{long_name}{units}".format(**get_variable_meta(dset, field_name))

    dims = _field.shape
    ens_out = np.ma.zeros([*dims, len(files)])
    ens_out[..., 0] = _field
    for idx, file_name in enumerate(files[1:]):
        with Dataset(file_name, "r") as dset:
            _field = dset.variables[field_name][:].squeeze()
            ens_out[..., idx + 1] = _field

    if mask_value:
        ens_out = np.ma.masked_less(ens_out, mask_value)

    return {"data": ens_out, "desc": var_desc}


def variable_set(name):
    var_sets = fn.read_json(os.path.join(os.path.dirname(__file__), "ksxx_vars.json"))
    try:
        the_set = var_sets[name.lower()]
        return set(the_set)
    except KeyError as e:
        six.raise_from(
            argparse.ArgumentTypeError(
                "Unknown variable set! Known sets are {}".format(var_sets.keys())
            ),
            e,
        )


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-c",
        "--config",
        type=fn.read_json,
        help="A JSON config file containing a `ks` dictionary defining "
        + "the options. NOTE: command line options will override file options.",
    )

    parser.add_argument("--test-case", default="default", help="Name of the test case.")

    parser.add_argument(
        "--test-dir",
        default=os.path.join(os.getcwd(), "archive"),
        help="Location of the test case run files.",
    )

    parser.add_argument(
        "--ref-case", default="fast", help="Name of the reference case."
    )

    parser.add_argument(
        "--ref-dir",
        default=os.path.join(os.getcwd(), "archive"),
        help="Location of the reference case run files.",
    )

    parser.add_argument(
        "--var-set",
        default="default",
        type=variable_set,
        help="Name of the variable set to analyze.",
    )

    parser.add_argument(
        "--ninst",
        default=30,
        type=int,
        help="The number of instances (should be the same for " "both cases).",
    )

    parser.add_argument(
        "--critical",
        default=13,
        type=int,
        help=(
            "The critical value (desired significance level) for rejecting the "
            "null hypothesis. Only used when --uncorrected / -u flag is passed. "
            "Otherwise the critical value is 1."
        ),
    )

    parser.add_argument("--img-dir", default=os.getcwd(), help="Image output location.")

    parser.add_argument(
        "--img-fmt", default="png", type=str, help="Format for output images"
    )

    parser.add_argument(
        "--component", default="eam", help="Model component name (e.g. eam, cam, ...)"
    )

    parser.add_argument(
        "--uncorrected",
        "-u",
        action="store_true",
        default=False,
        help="Do not use FDR correction to compute global pass / fail, will ",
    )

    parser.add_argument(
        "--alpha", default=0.05, type=float, help="Alpha threshold for pass / fail"
    )

    args, _ = parser.parse_known_args(args)

    # use config file arguments, but override with command line arguments
    if args.config:
        default_args = parser.parse_args([])

        for (
            key,
            val,
        ) in vars(args).items():
            if val != vars(default_args)[key]:
                args.config["ks"][key] = val

        config_arg_list = []
        _ = [
            config_arg_list.extend(["--" + key, str(val)])
            for key, val in args.config["ks"].items()
            if key != "config"
        ]
        args, _ = parser.parse_known_args(config_arg_list)

    return args


def col_fmt(dat):
    """Format results for table output."""
    if dat is None:
        _out = "-"
    elif isinstance(dat, (tuple, list)):
        try:
            _out = "{:.3e}, {:.3e}".format(*dat)
        except TypeError:
            _out = dat
    elif isinstance(dat, float):
        _out = f"{dat:.3e}"
    else:
        _out = "-"

    return _out


def run(name, config):
    """
    Runs the analysis.

    Args:
        name: The name of the test
        config: A dictionary representation of the configuration file

    Returns:
       The result of elements.page with the list of elements to display
    """

    config_arg_list = []
    [config_arg_list.extend(["--" + key, str(val)]) for key, val in config.items()]

    args = parse_args(config_arg_list)

    args.img_dir = os.path.join(livvkit.output_dir, "validation", "imgs", name)
    fn.mkdir_p(args.img_dir)

    details, img_gal = main(args)

    table_data = pd.DataFrame(details).T
    uc_rejections = (table_data["K-S test p-val"] < args.alpha).sum()
    _hdrs = [
        "h",
        "K-S test stat",
        "K-S test p-val",
        "K-S test p-val cor",
        "T test stat",
        "T test p-val",
        "mean test case",
        "mean ref. case",
        "std test case",
        "std ref. case",
    ]
    table_data = table_data[_hdrs]
    for _hdr in _hdrs[1:]:
        table_data[_hdr] = table_data[_hdr].apply(col_fmt)

    tables = [
        el.Table("Rejected", data=table_data[table_data["h"] == "reject"]),
        el.Table("Accepted", data=table_data[table_data["h"] == "accept"]),
        el.Table("Null", data=table_data[~table_data["h"].isin(["accept", "reject"])]),
    ]

    bib_html = bib2html(os.path.join(os.path.dirname(__file__), "ksxx.bib"))

    tabs = el.Tabs(
        {"Figures": img_gal, "Details": tables, "References": [el.RawHTML(bib_html)]}
    )
    rejects = [var for var, dat in details.items() if dat["h"] == "reject"]
    if args.uncorrected:
        critical = args.critical
    else:
        critical = 1

    results = el.Table(
        title="Results",
        data=OrderedDict(
            {
                # 'Test status': ['pass' if len(rejects) < args.critical else 'fail'],
                "Test status": ["pass" if len(rejects) < critical else "fail"],
                "Variables analyzed": [len(details.keys())],
                "Rejecting": [len(rejects)],
                "Critical value": [int(critical)],
                "Ensembles": [
                    "statistically identical"
                    if len(rejects) < critical
                    else "statistically different"
                ],
                "Un-corrected rejections": [uc_rejections],
            }
        ),
    )

    # FIXME: Put into a ___ function
    page = el.Page(name, __doc__.replace("\n\n", "<br><br>"), elements=[results, tabs])
    return page


def case_files(args):
    # ensure unique case names for the dictionary
    key1 = args.test_case
    key2 = args.ref_case
    if args.test_case == args.ref_case:
        key1 += "1"
        key2 += "2"

    f_sets = {
        # note in eamxx, we use the 'h' history file
        key1: component_monthly_files(args.test_dir, args.component, args.ninst, hist_name="h"),
        key2: component_monthly_files(args.ref_dir, args.component, args.ninst, hist_name="h"),
    }

    for key in f_sets:
        # Require case files for at least the last 12 months.
        if any(list(map(lambda x: x == [], f_sets[key].values()))[-12:]):
            raise EVVException(
                "Could not find all the required case files for case: {}".format(key)
            )

    return f_sets, key1, key2


def print_summary(summary):
    print("    Kolmogorov-Smirnov Test: {}".format(summary[""]["Case"]))
    print("      Variables analyzed: {}".format(summary[""]["Variables analyzed"]))
    print("      Rejecting: {}".format(summary[""]["Rejecting"]))
    print("      Critical value: {}".format(summary[""]["Critical value"]))
    print("      Ensembles: {}".format(summary[""]["Ensembles"]))
    print("      Test status: {}\n".format(summary[""]["Test status"]))


def print_details(details):
    for set_ in details:
        print("-" * 80)
        print(set_)
        print("-" * 80)
        pprint(details[set_])


def summarize_result(results_page):
    summary = {"Case": results_page.title}

    for elem in results_page.elements:
        if isinstance(elem, el.Table) and elem.title == "Results":
            summary["Test status"] = elem.data["Test status"][0]
            summary["Variables analyzed"] = elem.data["Variables analyzed"][0]
            summary["Rejecting"] = elem.data["Rejecting"][0]
            summary["Critical value"] = elem.data["Critical value"][0]
            summary["Ensembles"] = elem.data["Ensembles"][0]
            summary["Uncorrected Rejections"] = elem.data["Un-corrected rejections"][0]
            break

    return {"": summary}


def populate_metadata():
    """
    Generates the metadata responsible for telling the summary what
    is done by this module's run method
    """

    metadata = {
        "Type": "ValSummary",
        "Title": "Validation",
        "TableTitle": "Kolmogorov-Smirnov test",
        "Headers": [
            "Test status",
            "Variables analyzed",
            "Rejecting",
            "Critical value",
            "Ensembles",
        ],
    }
    return metadata


def compute_details(annual_avgs, common_vars, args):
    """Compute the detail table, perform a T Test and K-S test for each variable."""
    details = LIVVDict()
    for var in sorted(common_vars):
        annuals_1 = annual_avgs.query(
            "case == @args.test_case & variable == @var"
        ).monthly_mean.values
        annuals_2 = annual_avgs.query(
            "case == @args.ref_case & variable == @var"
        ).monthly_mean.values

        ttest_t, ttest_p = stats.ttest_ind(
            annuals_1, annuals_2, equal_var=False, nan_policy=str("omit")
        )
        ks_d, ks_p = stats.ks_2samp(annuals_1, annuals_2)

        if np.isnan([ttest_t, ttest_p]).any() or np.isinf([ttest_t, ttest_p]).any():
            ttest_t = None
            ttest_p = None

        details[var]["T test stat"] = ttest_t
        details[var]["T test p-val"] = ttest_p

        details[var]["K-S test stat"] = ks_d
        details[var]["K-S test p-val"] = ks_p

        details[var]["mean test case"] = annuals_1.mean()
        details[var]["mean ref. case"] = annuals_2.mean()

        details[var]["max test case"] = annuals_1.max()
        details[var]["max ref. case"] = annuals_2.max()

        details[var]["min test case"] = annuals_1.min()
        details[var]["min ref. case"] = annuals_2.min()

        details[var]["std test case"] = annuals_1.std()
        details[var]["std ref. case"] = annuals_2.std()

    # Now that the details have been computed, perform the FDR correction
    # Convert to a Dataframe, transposed so that the index is the variable name
    detail_df = pd.DataFrame(details).T
    # Create a null hypothesis rejection column for un-corrected p-values
    detail_df["h_uc"] = detail_df["K-S test p-val"] < args.alpha

    (detail_df["h_c"], detail_df["K-S test p-val cor"]) = smm.fdrcorrection(
        detail_df["K-S test p-val"], alpha=args.alpha, method="p", is_sorted=False
    )
    if args.uncorrected:
        _testkey = "h_uc"
    else:
        _testkey = "h_c"

    for var in common_vars:
        details[var]["K-S test p-val cor"] = detail_df.loc[var, "K-S test p-val cor"]

        if details[var]["T test stat"] is None:
            details[var]["h"] = "-"
        elif detail_df.loc[var, _testkey]:
            details[var]["h"] = "reject"
        else:
            details[var]["h"] = "accept"

    return details


def main(args):
    ens_files, key1, key2 = case_files(args)
    if args.test_case == args.ref_case:
        args.test_case = key1
        args.ref_case = key2

    monthly_avgs = gather_monthly_averages(ens_files, args.var_set)
    annual_avgs = (
        monthly_avgs.groupby(["case", "variable", "instance"])
        .monthly_mean.aggregate(monthly_to_annual_avg)
        .reset_index()
    )

    # now, we got the data, so let's get some stats
    test_set = set(monthly_avgs[monthly_avgs.case == args.test_case].variable.unique())
    ref_set = set(monthly_avgs[monthly_avgs.case == args.ref_case].variable.unique())
    common_vars = list(test_set & ref_set)
    if not common_vars:
        raise EVVException(
            "No common variables between {} and {} to analyze!".format(
                args.test_case, args.ref_case
            )
        )

    images = {"accept": [], "reject": [], "-": []}
    details = compute_details(annual_avgs, common_vars, args)

    for var in sorted(common_vars):
        annuals_1 = annual_avgs.query(
            "case == @args.test_case & variable == @var"
        ).monthly_mean.values
        annuals_2 = annual_avgs.query(
            "case == @args.ref_case & variable == @var"
        ).monthly_mean.values

        img_file = os.path.relpath(
            os.path.join(args.img_dir, f"{var}.{args.img_fmt}"), os.getcwd()
        )
        prob_plot(
            annuals_1,
            annuals_2,
            annuals_1.shape[0] // 2,
            img_file,
            test_name=args.test_case,
            ref_name=args.ref_case,
            pf=details[var]["h"],
            combine_hist=True,
        )
        _desc = monthly_avgs.query(
            "case == @args.test_case & variable == @var"
        ).desc.values[0]
        img_desc = (
            "Mean annual global average of {var}{desc} for <em>{testcase}</em> "
            "is {testmean:.3e} and for <em>{refcase}</em> is {refmean:.3e}. "
            "Pass (fail) is indicated by {cpass} ({cfail}) coloring of the "
            "plot markers and bars.".format(
                var=var,
                desc=_desc,
                testcase=args.test_case,
                testmean=details[var]["mean test case"],
                refcase=args.ref_case,
                refmean=details[var]["mean ref. case"],
                cfail=human_color_names["fail"][0],
                cpass=human_color_names["pass"][0],
            )
        )

        img_link = Path(*Path(args.img_dir).parts[-2:], Path(img_file).name)
        _img = el.Image(
            var, img_desc, img_link, relative_to="", group=details[var]["h"]
        )
        images[details[var]["h"]].append(_img)

    gals = []
    for group in ["reject", "accept", "-"]:
        _group_name = {
            "reject": "Failed variables",
            "accept": "Passed variables",
            "-": "Null variables",
        }
        if images[group]:
            gals.append(el.Gallery(_group_name[group], images[group]))

    return details, gals


if __name__ == "__main__":
    print_details(main(parse_args()))
