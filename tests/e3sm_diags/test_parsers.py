import argparse
import sys
from typing import List

import pytest

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parser.area_mean_time_series_parser import AreaMeanTimeSeriesParser
from e3sm_diags.parser.arm_diags_parser import ARMDiagsParser
from e3sm_diags.parser.core_parser import CoreParser
from e3sm_diags.parser.diurnal_cycle_parser import DiurnalCycleParser
from e3sm_diags.parser.enso_diags_parser import EnsoDiagsParser
from e3sm_diags.parser.meridional_mean_2d_parser import MeridionalMean2dParser
from e3sm_diags.parser.qbo_parser import QboParser
from e3sm_diags.parser.streamflow_parser import StreamflowParser
from e3sm_diags.parser.tc_analysis_parser import TCAnalysisParser
from e3sm_diags.parser.zonal_mean_2d_parser import ZonalMean2dParser
from e3sm_diags.parser.zonal_mean_2d_stratosphere_parser import (
    ZonalMean2dStratosphereParser,
)


class TestCoreParser:
    @pytest.fixture(autouse=True)
    def setup(self):
        # The base arguments added to `CoreParser.parser` through
        # `CoreParser.add_arguments()`.
        self.base_args = [
            "parameters",
            "other_parameters",
            "num_workers",
            "scheduler_addr",
            "granulate",
            "selectors",
            "set_name",
            "reference_data_set",
            "reference_data_path",
            "ref_timeseries_input",
            "ref_start_yr",
            "ref_end_yr",
            "ref_start_time_slice",
            "ref_end_time_slice",
            "ref_name",
            "ref_file",
            "test_data_set",
            "test_data_path",
            "test_timeseries_input",
            "test_start_yr",
            "test_end_yr",
            "test_start_time_slice",
            "test_end_time_slice",
            "test_file",
            "results_dir",
            "sets",
            "dataset",
            "run_type",
            "variables",
            "plevs",
            "plot_plevs",
            "plot_log_plevs",
            "seasons",
            "regions",
            "regrid_tool",
            "regrid_method",
            "case_id",
            "output_format",
            "output_format_subplot",
            "canvas_size_w",
            "canvas_size_h",
            "figsize",
            "dpi",
            "arrows",
            "logo",
            "contour_levels",
            "diff_levels",
            "reference_name",
            "test_name",
            "short_test_name",
            "diff_name",
            "main_title",
            "reference_title",
            "test_title",
            "diff_title",
            "reference_colormap",
            "test_colormap",
            "diff_colormap",
            "reference_units",
            "test_units",
            "diff_units",
            "backend",
            "multiprocessing",
            "save_netcdf",
            "no_viewer",
            "debug",
        ]

    def test__init__(self):
        CoreParser()

    def test_check_values_of_params_does_not_raise_error(self):
        param = CoreParameter()
        param.reference_data_path = "path"
        param.test_data_path = "path"
        param.results_dir = "path"

        # Should not raise RunTimeError
        CoreParser.check_values_of_params([param])

    def test_add_arguments_adds_base_args(self):
        parser = CoreParser()
        assert _is_args_added(parser, self.base_args)

    def test_parse_args_sets_cmd_used_and_returns_list_of_args(self):
        parser = CoreParser()

        result = parser.parse_args()
        result_args = list(vars(result).keys())

        assert result_args == self.base_args
        assert parser.cmd_used == sys.argv

    def test_view_args_returns_parser_namespace(self):
        parser = CoreParser()

        result = parser.view_args()
        assert isinstance(result, argparse.Namespace)

        result_args = list(vars(result).keys())
        assert result_args == self.base_args

    @pytest.mark.xfail
    def test_get_parameters_returns_cmdline_parameters(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_orig_parameters(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_other_parameters(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_default_vars(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_cmd_default_vars(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_parameter_from_defaults_of_the_command_line_arg(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_parameter_from_defaults_of_the_parameter_class(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_cartesian_product_of_granulate_attr(self):
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_returns_parameters_created_by_running_from_CLI(
        self,
    ):
        # FIXME: Should we deprecate this method for running `e3sm_diags`?
        # https://e3sm-project.github.io/e3sm_diags/_build/html/main/config-run.html#e3sm-diags-p-older-method
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_returns_parameters_created_by_cfg_file(self):
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_checks_values_in_cfg_file_and_returns_parameter(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_only_uses_argparse_values_and_returns_parameters(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_select_returns_cmdline_parameters_that_are_subset_of_the_main_parameters(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_select_returns_orig_parameters_that_are_subset_of_the_main_parameters(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_select_returns_other_parameters_that_are_subset_of_the_main_parameters(
        self,
    ):
        assert 0


def test_area_mean_time_series_parser_initializes():
    parser = AreaMeanTimeSeriesParser()
    custom_args = [
        "ref_names",
        "ref_timeseries_input",
        "test_timeseries_input",
        "start_yr",
        "end_yr",
    ]

    assert _is_args_added(parser, custom_args)


def test_arms_diags_parser_initializes():
    parser = ARMDiagsParser()
    custom_args = ["ref_names"]

    assert _is_args_added(parser, custom_args)


def test_diurnal_cycle_parser_initializes():
    parser = DiurnalCycleParser()
    custom_args = [
        "ref_timeseries_input",
        "test_timeseries_input",
        "start_yr",
        "end_yr",
    ]

    assert _is_args_added(parser, custom_args)


def test_enso_diags_parser_initializes():
    parser = EnsoDiagsParser()
    custom_args = [
        "ref_names",
        "ref_timeseries_input",
        "test_timeseries_input",
        "start_yr",
        "end_yr",
    ]

    assert _is_args_added(parser, custom_args)


def test_meridional_mean_2d_parser_initializes():
    parser = MeridionalMean2dParser()
    custom_args = [
        "plevs",
        "plot_plevs",
        "plot_log_plevs",
    ]

    assert _is_args_added(parser, custom_args)


def test_qbo_parser_initializes():
    parser = QboParser()
    custom_args = [
        "ref_timeseries_input",
        "test_timeseries_input",
        "start_yr",
        "end_yr",
    ]

    assert _is_args_added(parser, custom_args)


def test_streamflow_parser_initializes():
    parser = StreamflowParser()
    custom_args = [
        "gauges_path",
        "max_num_gauges",
        "print_statements",
        "ref_timeseries_input",
        "test_timeseries_input",
        "start_yr",
        "end_yr",
    ]

    assert _is_args_added(parser, custom_args)


def test_tc_analysis_parser_initializes():
    TCAnalysisParser()


def test_zonal_mean_2d_parser_initializes():
    parser = ZonalMean2dParser()
    custom_args = ["plevs", "plot_plevs", "plot_log_plevs"]

    assert _is_args_added(parser, custom_args)


def test_zonal_mean_2d_stratosphere_parser_initializes():
    parser = ZonalMean2dStratosphereParser()
    custom_args = ["plevs", "plot_plevs", "plot_log_plevs"]

    assert _is_args_added(parser, custom_args)


def _is_args_added(parser: CoreParser, custom_args: List[str]) -> bool:
    """Checks the parser's custom args are added to argparse.ArgumentParser.

    Parameters
    ----------
    parser : CoreParser
        The CoreParser-based object.
    custom_args : List[str]
        The list of custom arguments for this parser.

    Returns
    -------
    bool
        True if all custom arguments are added, else False.
    """
    namespace, _ = parser.parser.parse_known_args()
    namespace_args = vars(namespace).keys()

    for arg in custom_args:
        if arg not in namespace_args:
            return False

    return True
