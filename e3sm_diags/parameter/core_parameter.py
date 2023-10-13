import copy
import importlib
import sys
from typing import Any, Dict, List

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


class CoreParameter:
    def __init__(self):
        # File I/O
        # ------------------------
        # FIXME: (REQUIRED) attributes should be set through `__init__` since
        # they are required, but they are currently set after initializing the
        # CoreParameter object (e.g., param.results_dir = "dir"). The default
        # values are set to "", and `self.check_values()` will validate that
        # it is not "" since they must be set by the user.

        # (REQUIRED) Path to the reference (obs) data. If there are multiple
        # datasets in the path, use ref_name parameter to specify which dataset
        # should be used.
        self.reference_data_path: str = ""

        # (REQUIRED) Path to the test (model) data.
        self.test_data_path: str = ""

        # (REQUIRED) The name of the folder where all runs will be stored.
        self.results_dir: str = ""

        # The name of the folder where the results (plots and nc files) will be
        # stored for a single run
        self.case_id = ""

        # Set to True to not generate a Viewer for the result.
        self.no_viewer: bool = False

        # If True, stops running all of the diagnostics on the first failure.
        # If False (the default), all errors are caught and ignored. If there
        # was an error and a plot could not be created, there’s a ‘—’ for that
        # set of parameters in the viewer.
        self.debug: bool = False

        # Time series data
        # ------------------------
        # Set to True if the ref data is in timeseries format. Default False. If
        # True, both ref_start_yr and ref_end_yr must also be set.
        self.ref_timeseries_input: bool = False

        # Set to True if the test data is in timeseries format. Default False.
        # If True, both test_start_yr and # test_end_yr must also be set.
        self.test_timeseries_input: bool = False

        # Diagnostic run settings
        # ------------------------
        # The supported run type for the diagnostics. Possible options are:
        # 'model_vs_obs' (by default), 'model_vs_model', or 'obs_vs_obs'.
        self.run_type: str = "model_vs_obs"

        # A list of the sets to be run. Default is all sets
        self.sets: List[str] = [
            "zonal_mean_xy",
            "zonal_mean_2d",
            "zonal_mean_2d_stratosphere",
            "meridional_mean_2d",
            "lat_lon",
            "polar",
            "area_mean_time_series",
            "cosp_histogram",
            "enso_diags",
            "qbo",
            "streamflow",
            "diurnal_cycle",
            "arm_diags",
            "tc_analysis",
            "annual_cycle_zonal_mean",
            "lat_lon_land",
            "lat_lon_river",
            "aerosol_aeronet",
            "aerosol_budget",
        ]

        # The current set that is being ran when looping over sets in
        # `e3sm_diags_driver.run_diag()`.
        self.current_set: str = ""

        self.variables: List[str] = []
        self.seasons: List[str] = ["ANN", "DJF", "MAM", "JJA", "SON"]
        self.regions: List[str] = ["global"]

        self.regrid_tool: str = "esmf"
        self.regrid_method: str = "conservative"

        self.plevs: List[float] = []
        self.plot_log_plevs: bool = False
        self.plot_plevs: bool = False

        self.multiprocessing: bool = False
        self.num_workers: int = 4

        # Diagnostic plot settings
        # ------------------------
        self.main_title: str = ""
        self.backend: str = "mpl"
        self.save_netcdf: bool = False

        # Plot format settings
        self.output_file: str = ""
        self.output_format: List[str] = ["png"]
        self.output_format_subplot: List[str] = []
        self.canvas_size_w: int = 1212
        self.canvas_size_h: int = 1628
        self.figsize: List[float] = [8.5, 11.0]
        self.dpi: int = 150
        self.arrows: bool = True
        self.logo: bool = False
        self.contour_levels: List[str] = []

        # Test plot settings
        self.test_name: str = ""
        self.test_name_yrs: str = ""
        self.short_test_name: str = ""
        self.test_title: str = ""
        self.test_colormap: str = "cet_rainbow.rgb"
        self.test_units: str = ""

        # Reference plot settings
        self.ref_name: str = ""
        self.ref_name_yrs: str = ""
        self.reference_name: str = ""
        self.short_ref_name: str = ""
        self.reference_title: str = ""
        self.reference_colormap: str = "cet_rainbow.rgb"
        self.reference_units: str = ""

        # The reference file name based on the season and other parameters, for
        # climo files only.
        self.ref_file: str = ""

        # Variable ID is used as the reference file ID in `general.save_ncfiles()`
        self.var_id: str = ""

        # Variable region, used in polar plots.
        self.var_region: str = ""

        # Difference plot settings
        self.diff_name: str = ""
        self.diff_title: str = "Model - Observation"
        self.diff_colormap: str = "diverging_bwr.rgb"
        self.diff_levels: List[str] = []
        self.diff_units: str = ""
        self.diff_type: str = "absolute"

        # Other settings
        # ------------------------
        # TODO: Need documentation on these attributes here and
        # here: https://e3sm-project.github.io/e3sm_diags/_build/html/main/available-parameters.html
        self.dataset: str = ""
        self.granulate: List[str] = ["variables", "seasons", "plevs", "regions"]
        self.selectors: List[str] = ["sets", "seasons"]
        self.viewer_descr: Dict[str, str] = {}
        self.fail_on_incomplete: bool = False

        # List of user derived variables, set in `dataset.Dataset`.
        self.derived_variables: Dict[str, object] = {}

        # FIXME: This attribute is only used in `lat_lon_driver.py`
        self.model_only: bool = False

    def __add__(self, other):
        """Copies attributes from the other CoreParameter object to this one.

        This method overrides existing attributes if they are already set
        using the new values from the other object.

        It is based on cdp.cdp_parameter.CDPParameter.__add__().

        Parameters
        ----------
        other : CoreParameter
            The other CoreParameter (or sub-class) object.

        Returns
        -------
        CoreParameter
            A new instance of CoreParameter with new attributes.
        """
        new_obj = copy.deepcopy(self)
        new_attrs = other.__dict__

        for attr, value in new_attrs.items():
            setattr(new_obj, attr, value)

        return new_obj

    def check_values(self):
        # The default values for these attributes are set to "" in `__init__`.
        # FIXME: The user should pass these values into the object
        # initialization rather than having to check they are set here.
        must_have_params = ["reference_data_path", "test_data_path", "results_dir"]
        for param in must_have_params:
            if getattr(self, param) == "":
                msg = "You need to specify {p} in the parameters file or via the command line using --{p}".format(
                    p=param
                )
                raise RuntimeError(msg)

        if self.ref_timeseries_input and not (
            hasattr(self, "ref_start_yr") and hasattr(self, "ref_end_yr")
        ):
            msg = (
                "You need to define both the 'ref_start_yr' and 'ref_end_yr' parameter."
            )
            raise RuntimeError(msg)

        if self.test_timeseries_input and not (
            hasattr(self, "test_start_yr") and hasattr(self, "test_end_yr")
        ):
            msg = "You need to define both the 'test_start_yr' and 'test_end_yr' parameter."
            raise RuntimeError(msg)

    def _run_diag(self) -> List[Any]:
        """Run the diagnostics for each set in the parameter.

        Additional CoreParameter (or CoreParameter sub-class) objects are derived
        from the CoreParameter `sets` attribute, hence this function returns a
        list of CoreParameter objects.

        This method loops over the parameter's diagnostic sets and attempts to
        import and call the related `run_diags()` function.

        Returns
        -------
        List[Any]
            The list of CoreParameter objects with results from the diagnostic run.
            NOTE: `Self` type is not yet supported by mypy.
        """
        results = []

        for set_name in self.sets:
            self.current_set = set_name
            # FIXME: This is a shortcut to importing `run_diag`, but can break
            # easily because the module driver is statically imported via string.
            # Instead, the import should be done more progammatically via
            # direct Python import.
            mod_str = "e3sm_diags.driver.{}_driver".format(set_name)

            # Check if there is a matching driver module for the `set_name`.
            try:
                module = importlib.import_module(mod_str)
            except ModuleNotFoundError as e:
                logger.error(f"'Error with set name {set_name}'", exc_info=e)
                continue

            # If the module exists, call the driver module's `run_diag` function.
            try:
                single_result = module.run_diag(self)
                results.append(single_result)
            except Exception:
                logger.exception(f"Error in {mod_str}", exc_info=True)

                if self.debug:
                    sys.exit()

        return results
