"""
Multivariate test for climate reproducibility using the Kolmogrov-Smirnov (K-S)
test and based on The CESM/E3SM model's multi-instance capability is used to
conduct an ensemble of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

import os
import json
import logging
import shutil
import xml.etree.ElementTree as ET

from distutils import dir_util

import numpy as np
import netCDF4 as nc

import CIME.test_status
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.XML.machines import Machines

import evv4esm  # pylint: disable=import-error
from evv4esm.__main__ import main as evv  # pylint: disable=import-error

evv_lib_dir = os.path.abspath(os.path.dirname(evv4esm.__file__))
logger = logging.getLogger(__name__)
NINST = 30

PERT_FILE_TEMPLATE = "mpaso_oQ240_perturbed_inic_{ens:04d}.nc"

INIT_FILE_NAME = "20231204.GMPAS-NYF.T62_oQU240.v3beta.mpaso.rst.0401-01-01_00000.nc"

# No output is on by default for mpas-ocean
OCN_TEST_VARS = [
    "daysSinceStartOfSim",
    "ssh",
    "velocityZonal",
    "velocityMeridional",
    "activeTracers",
]

# daysSinceStartOfSim, icePresent, iceAreaCell, and iceVolumeCell are on by default
ICE_TEST_VARS = ["uVelocityGeo", "vVelocityGeo", "icePressure", "divergence", "shear"]


def perturb_init(infile, field_name, outfile, seed=None):
    """
    Create perturbed initial conditions file from another netCDF using random uniform.

    Parameters
    ----------
    infile : string
        Path to initial conditions input file which will be perturbed
    field_name : string
        Variable name in netCDF file to perturb
    outfile : string
        Path to output file
    seed : int
        Integer seed for the random number generator (optional)

    """

    init_f = nc.Dataset(infile, "r")
    field_temp = init_f.variables[field_name]

    # Initialize numpy random generator with a seed
    rng = np.random.default_rng(seed)

    # Perturbation between -1e-14 -- +1e-14, same size as the input field
    perturbation = rng.uniform(-1e-14, 1e-14, field_temp[:].shape)
    field = field_temp[:] + perturbation

    shutil.copy(infile, outfile)
    with nc.Dataset(outfile, "a") as out_f:

        field_out = out_f.variables[field_name]
        field_out[:] = field

    init_f.close()


def modify_stream(
    stream_file,
    input_file=None,
    var_names=None,
    output_stream_name="timeSeriesStatsClimatologyOutput",
):
    """Add variables to XML streams to specify output."""
    streams = ET.parse(stream_file)

    # Modify the input stream
    if input_file is not None:
        input_stream = list(
            filter(
                lambda x: x.get("name") == "input",
                streams.getroot().findall("immutable_stream"),
            )
        )
        input_stream[0].set("filename_template", input_file)

    # Add variables to a particular stream (Climatology output by default)
    clim_stream = list(
        filter(
            lambda x: x.get("name") == output_stream_name,
            streams.getroot().findall("stream"),
        )
    )
    clim_stream = clim_stream[0]
    if var_names is not None:
        # Convert the variable name into xml elements
        var_elements = [ET.Element("var", {"name": _var}) for _var in var_names]

        for clim_element in var_elements:
            clim_stream.append(clim_element)

    # Only output climatology files yearly
    clim_stream.set("output_interval", "01-00-00_00:00:00")

    if input_file is not None or var_elements is not None:
        # Only re-write the stream if changes were made
        streams.write(stream_file)


class MVKO(SystemTestsCommon):
    """MVK-Ocean SystemTest class."""

    def __init__(self, case):
        """
        initialize an object interface to the MVKO test
        """
        SystemTestsCommon.__init__(self, case)

        if self._case.get_value("MODEL") == "e3sm":
            self.ocn_component = "mpaso"
            self.ice_component = "mpassi"
        else:
            self.ocn_component = "pop"
            self.ice_component = "cice"

        if (
            self._case.get_value("RESUBMIT") == 0
            and self._case.get_value("GENERATE_BASELINE") is False
        ):
            self._case.set_value("COMPARE_BASELINE", True)
        else:
            self._case.set_value("COMPARE_BASELINE", False)

    def build_phase(self, sharedlib_only=False, model_only=False):
        """
        Create instance namelists, stream files, and initial conditions.
        """
        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        caseroot = self._case.get_value("CASEROOT")
        if not model_only:
            logging.warning("Starting to build multi-instance exe")
            for comp in self._case.get_values("COMP_CLASSES"):
                self._case.set_value(f"NTHRDS_{comp}", 1)

                ntasks = self._case.get_value(f"NTASKS_{comp}")

                self._case.set_value(f"NTASKS_{comp}", ntasks * NINST)
                if comp != "CPL":
                    self._case.set_value(f"NINST_{comp}", NINST)

            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)
        rundir = self._case.get_value("RUNDIR")

        init_file_path = os.path.join(
            self._case.get_value("DIN_LOC_ROOT"),
            "ocn/mpas-o/oQU240/inic/gmpas-nyf",
            INIT_FILE_NAME,
        )
        # Write yearly averages to custom output file
        tss_climatology_config = [
            "config_am_timeseriesstatsclimatology_enable = .true.\n",
            "config_am_timeseriesstatsclimatology_backward_output_offset = '01-00-00_00:00:00'\n",
            "config_am_timeseriesstatsclimatology_compute_interval = '00-00-00_01:00:00'\n",
            "config_am_timeseriesstatsclimatology_compute_on_startup = .false.\n",
            "config_am_timeseriesstatsclimatology_duration_intervals = '01-00-00_00:00'\n",
            "config_am_timeseriesstatsclimatology_operation = 'avg'\n",
            "config_am_timeseriesstatsclimatology_output_stream = 'timeSeriesStatsClimatologyOutput'\n",
            "config_am_timeseriesstatsclimatology_reference_times = '00-01-01_00:00:00'\n",
            "config_am_timeseriesstatsclimatology_repeat_intervals = '01-00-00_00:00:00'\n",
            "config_am_timeseriesstatsclimatology_reset_intervals = '0001-00-00_00:00:00'\n",
            "config_am_timeseriesstatsclimatology_restart_stream = 'timeSeriesStatsClimatologyRestart'\n",
            "config_am_timeseriesstatsclimatology_write_on_startup = .false.\n",
        ]

        for iinst in range(1, NINST + 1):
            pert_file_name = PERT_FILE_TEMPLATE.format(ens=iinst)
            pert_file = os.path.join(rundir, pert_file_name)
            if not os.path.exists(rundir):
                logging.warning("CREATE %s", rundir)
                os.mkdir(rundir)
            perturb_init(init_file_path, field_name="temperature", outfile=pert_file)

            # Set up ocean namelist to specify climatology output and reduce other output
            with open(
                f"user_nl_{self.ocn_component}_{iinst:04d}", "w", encoding="utf-8"
            ) as nl_ocn_file:

                for _config in tss_climatology_config:
                    nl_ocn_file.write(_config)

                # Disable ocean output we don't use for this test
                nl_ocn_file.write("config_am_highfrequencyoutput_enable = .false.\n")
                nl_ocn_file.write(
                    "config_am_timeseriesstatsmonthlymax_enable = .false.\n"
                )
                nl_ocn_file.write(
                    "config_am_timeseriesstatsmonthlymin_enable = .false.\n"
                )
                nl_ocn_file.write("config_am_timeseriesstatsmonthly_enable = .false.\n")
                nl_ocn_file.write("config_am_eddyproductvariables_enable = .false.\n")
                nl_ocn_file.write(
                    "config_am_meridionalheattransport_enable = .false.\n"
                )
                nl_ocn_file.write("config_am_oceanheatcontent_enable = .false.\n")
                # But we want to make sure globalStats output is on
                nl_ocn_file.write("config_am_globalstats_enable = .true.\n")

            # Set up sea ice namelist to reduce output
            with open(
                f"user_nl_{self.ice_component}_{iinst:04d}",
                "w",
                encoding="utf-8",
            ) as nl_ice_file:
                for _config in tss_climatology_config:
                    nl_ice_file.write(_config)
                nl_ice_file.write("config_am_timeseriesstatsdaily_enable = .false.\n")
                nl_ice_file.write("config_am_timeseriesstatsmonthly_enable = .false.\n")
                nl_ice_file.write("config_am_regionalstatistics_enable = .false.\n")

            if self.ocn_component == "mpaso":
                # Set up streams file to get perturbed init file and do yearly climatology
                ocn_stream_file = os.path.join(
                    caseroot,
                    "SourceMods",
                    "src.mpaso",
                    f"streams.ocean_{iinst:04d}",
                )
                ice_stream_file = os.path.join(
                    caseroot,
                    "SourceMods",
                    "src.mpassi",
                    f"streams.seaice_{iinst:04d}",
                )
                # buildnamelist creates the streams file for each instance, copy this to the
                # SourceMods directory to make changes to the correct version
                shutil.copy(f"{rundir}/streams.ocean_{iinst:04d}", ocn_stream_file)
                shutil.copy(f"{rundir}/streams.seaice_{iinst:04d}", ice_stream_file)

                modify_stream(
                    ocn_stream_file, input_file=pert_file, var_names=OCN_TEST_VARS
                )
                modify_stream(ice_stream_file, var_names=ICE_TEST_VARS)

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

    def run_phase(self):
        """Run the model."""
        self.run_indv()

    def _generate_baseline(self):
        """
        generate a new baseline case based on the current test
        """
        super()._generate_baseline()

        with CIME.utils.SharedArea():
            basegen_dir = os.path.join(
                self._case.get_value("BASELINE_ROOT"),
                self._case.get_value("BASEGEN_CASE"),
            )

            rundir = self._case.get_value("RUNDIR")
            ref_case = self._case.get_value("RUN_REFCASE")

            env_archive = self._case.get_env("archive")
            hists = env_archive.get_all_hist_files(
                self._case.get_value("CASE"),
                self.ocn_component,
                rundir,
                ref_case=ref_case,
            )
            logger.debug("MVKO additional baseline files: %s", hists)
            hists = [os.path.join(rundir, hist) for hist in hists]
            for hist in hists:
                basename = hist[hist.rfind(self.ocn_component) :]
                baseline = os.path.join(basegen_dir, basename)
                if os.path.exists(baseline):
                    os.remove(baseline)

                CIME.utils.safe_copy(hist, baseline, preserve_meta=False)

    def _compare_baseline(self):
        with self._test_status:
            if int(self._case.get_value("RESUBMIT")) > 0:
                # This is here because the comparison is run for each submission
                # and we only want to compare once the whole run is finished. We
                # need to return a pass here to continue the submission process.
                self._test_status.set_status(
                    CIME.test_status.BASELINE_PHASE, CIME.test_status.TEST_PEND_STATUS
                )
                return

            self._test_status.set_status(
                CIME.test_status.BASELINE_PHASE, CIME.test_status.TEST_FAIL_STATUS
            )

            run_dir = self._case.get_value("RUNDIR")
            case_name = self._case.get_value("CASE")
            base_dir = os.path.join(
                self._case.get_value("BASELINE_ROOT"),
                self._case.get_value("BASECMP_CASE"),
            )

            test_name = str(case_name.split(".")[-1])
            evv_config = {
                test_name: {
                    "module": os.path.join(evv_lib_dir, "extensions", "kso.py"),
                    "test-case": "Test",
                    "test-dir": run_dir,
                    "ref-case": "Baseline",
                    "ref-dir": base_dir,
                    "var-set": "default",
                    "ninst": NINST,
                    "critical": 0,
                    "component": self.ocn_component,
                    "alpha": 0.05,
                    "hist-name": "hist.am.timeSeriesStatsClimatology",
                }
            }

            json_file = os.path.join(run_dir, ".".join([case_name, "json"]))
            with open(json_file, "w", encoding="utf-8") as config_file:
                json.dump(evv_config, config_file, indent=4)

            evv_out_dir = os.path.join(run_dir, ".".join([case_name, "evv"]))
            evv(["-e", json_file, "-o", evv_out_dir])

            with open(
                os.path.join(evv_out_dir, "index.json"), encoding="utf-8"
            ) as evv_f:
                evv_status = json.load(evv_f)

            comments = ""
            for evv_ele in evv_status["Page"]["elements"]:
                if "Table" in evv_ele:
                    comments = "; ".join(
                        f"{key}: {val[0]}"
                        for key, val in evv_ele["Table"]["data"].items()
                    )
                    if evv_ele["Table"]["data"]["Test status"][0].lower() == "pass":
                        self._test_status.set_status(
                            CIME.test_status.BASELINE_PHASE,
                            CIME.test_status.TEST_PASS_STATUS,
                        )
                    break

            status = self._test_status.get_status(CIME.test_status.BASELINE_PHASE)
            mach_name = self._case.get_value("MACH")
            mach_obj = Machines(machine=mach_name)
            htmlroot = CIME.utils.get_htmlroot(mach_obj)
            urlroot = CIME.utils.get_urlroot(mach_obj)
            if htmlroot is not None:
                with CIME.utils.SharedArea():
                    dir_util.copy_tree(
                        evv_out_dir,
                        os.path.join(htmlroot, "evv", case_name),
                        preserve_mode=False,
                    )
                if urlroot is None:
                    urlroot = f"[{mach_name.capitalize()}_URL]"
                viewing = f"{urlroot}/evv/{case_name}/index.html"
            else:
                viewing = (
                    f"{evv_out_dir}\n"
                    "    EVV viewing instructions can be found at: "
                    "        https://github.com/E3SM-Project/E3SM/blob/master/cime/scripts/"
                    "climate_reproducibility/README.md#test-passfail-and-extended-output"
                )
            comments = (
                f"{CIME.test_status.BASELINE_PHASE} {status} for test '{test_name}'.\n"
                f"    {comments}\n"
                "    EVV results can be viewed at:\n"
                f"        {viewing}"
            )

            CIME.utils.append_testlog(comments, self._orig_caseroot)
