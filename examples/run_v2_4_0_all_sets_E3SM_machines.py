import os

from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.run import runner


def run_compy(html_prefix):
    # Run the following first:
    # srun --pty --nodes=1 --time=01:00:00 /bin/bash
    # source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
    ref_data_prefix = "/compyfs/e3sm_diags_data/obs_for_e3sm_diags"
    test_data_prefix = "/compyfs/e3sm_diags_data/test_model_data_for_acme_diags"

    test_data_prefix2 = "/compyfs/fors729"

    d = dict()

    d["obs_climo"] = os.path.join(ref_data_prefix, "climatology/")
    d["test_climo"] = os.path.join(test_data_prefix, "climatology/")

    d["obs_ts"] = os.path.join(ref_data_prefix, "time-series/")
    d["test_ts"] = os.path.join(test_data_prefix, "time-series/E3SM_v1/")

    d["dc_obs_climo"] = d["obs_climo"]
    d["dc_test_climo"] = os.path.join(
        test_data_prefix, "climatology/diurnal_cycle_climatology/"
    )

    d["streamflow_obs_ts"] = d["obs_ts"]
    d["streamflow_test_ts"] = os.path.join(
        test_data_prefix2, "time-series/streamflow_ts_postprocessed/"
    )

    d["arm_obs"] = os.path.join(ref_data_prefix, "arm-diags-data/")
    d["arm_test"] = os.path.join(test_data_prefix2, "20210122.F2010.armsites/")

    return run_all_sets(html_prefix, d)


def run_lcrc(html_prefix):
    # Run the following first:
    # srun --pty --nodes=1 --time=01:00:00 /bin/bash
    # source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
    # Or: source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh
    ref_data_prefix = "/lcrc/group/e3sm/public_html/diagnostics/observations/Atm"
    test_data_prefix = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data"

    ref_data_prefix2 = (
        "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/obs"
    )
    test_data_prefix2 = (
        "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/test"
    )

    d = dict()

    d["obs_climo"] = os.path.join(ref_data_prefix, "climatology/")
    d["test_climo"] = os.path.join(test_data_prefix, "climatology/")

    d["obs_ts"] = os.path.join(ref_data_prefix, "time-series/")
    d["test_ts"] = os.path.join(test_data_prefix, "time-series/E3SM_v1/")

    d["dc_obs_climo"] = os.path.join(ref_data_prefix2, "climatology/")
    d["dc_test_climo"] = os.path.join(
        test_data_prefix2, "climatology/dc_climo_postprocessed/"
    )

    d["streamflow_obs_ts"] = os.path.join(ref_data_prefix2, "time-series/")
    d["streamflow_test_ts"] = os.path.join(
        test_data_prefix2, "time-series/streamflow_ts_postprocessed/"
    )

    d["arm_obs"] = os.path.join(ref_data_prefix, "arm-diags-data/")
    d["arm_test"] = os.path.join(test_data_prefix, "20210122.F2010.armsites/")

    return run_all_sets(html_prefix, d)


def run_nersc(html_prefix):
    # Run the following first:
    # salloc --nodes=1 --partition=regular --time=01:00:00 -C haswell
    # source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
    ref_data_prefix = "/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags"
    test_data_prefix = (
        "/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags"
    )

    test_data_prefix2 = "/global/cfs/cdirs/e3sm/zhang40/postprocessing_for_e3sm_diags"

    d = dict()

    d["obs_climo"] = os.path.join(ref_data_prefix, "climatology/")
    d["test_climo"] = os.path.join(test_data_prefix, "climatology/")

    d["obs_ts"] = os.path.join(ref_data_prefix, "time-series/")
    d["test_ts"] = os.path.join(test_data_prefix, "time-series/E3SM_v1/")

    d["dc_obs_climo"] = d["obs_climo"]
    d["dc_test_climo"] = os.path.join(
        test_data_prefix2,
        "diurnal_climo/20180215.DECKv1b_H1.ne30_oEC.edison/1980-2014/rgr/",
    )

    d["streamflow_obs_ts"] = d["obs_ts"]
    d["streamflow_test_ts"] = os.path.join(
        test_data_prefix2,
        "monthly_ts/20180215.DECKv1b_H1.ne30_oEC.edison/1980-2014/rgr/",
    )

    d["arm_obs"] = os.path.join(ref_data_prefix, "arm-diags-data/")
    d["arm_test"] = os.path.join(test_data_prefix, "20210122.F2010.armsites/")

    return run_all_sets(html_prefix, d)


def run_all_sets(html_prefix, d):
    param = CoreParameter()

    param.reference_data_path = d["obs_climo"]
    param.test_data_path = d["test_climo"]
    param.test_name = "20161118.beta0.FC5COSP.ne30_ne30.edison"
    param.seasons = [
        "ANN",
        "JJA",
    ]  # Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

    param.results_dir = os.path.join(html_prefix, "v2_4_0_all_sets")
    param.multiprocessing = True
    param.num_workers = 30

    # Set specific parameters for new sets
    enso_param = EnsoDiagsParameter()
    enso_param.reference_data_path = d["obs_ts"]
    enso_param.test_data_path = d["test_ts"]
    enso_param.test_name = "e3sm_v1"
    enso_param.start_yr = "1990"
    enso_param.end_yr = "1999"

    qbo_param = QboParameter()
    qbo_param.reference_data_path = d["obs_ts"]
    qbo_param.test_data_path = d["test_ts"]
    qbo_param.test_name = "e3sm_v1"
    qbo_param.start_yr = "1990"
    qbo_param.end_yr = "1999"

    ts_param = AreaMeanTimeSeriesParameter()
    ts_param.reference_data_path = d["obs_ts"]
    ts_param.test_data_path = d["test_ts"]
    ts_param.test_name = "e3sm_v1"
    ts_param.start_yr = "1990"
    ts_param.end_yr = "1999"

    dc_param = DiurnalCycleParameter()
    dc_param.reference_data_path = d["dc_obs_climo"]
    dc_param.test_data_path = d["dc_test_climo"]
    dc_param.test_name = "20180215.DECKv1b_H1.ne30_oEC.edison"
    dc_param.short_test_name = "DECKv1b_H1.ne30_oEC"
    # Plotting diurnal cycle amplitude on different scales. Default is True
    dc_param.normalize_test_amp = False

    streamflow_param = StreamflowParameter()
    streamflow_param.reference_data_path = d["streamflow_obs_ts"]
    streamflow_param.test_data_path = d["streamflow_test_ts"]
    streamflow_param.test_name = "20180215.DECKv1b_H1.ne30_oEC.edison"
    streamflow_param.test_start_yr = "1980"
    streamflow_param.test_end_yr = "2014"
    # Streamflow gauge station data range from year 1986 to 1995
    streamflow_param.ref_start_yr = "1986"
    streamflow_param.ref_end_yr = "1995"

    arm_param = ARMDiagsParameter()
    arm_param.reference_data_path = d["arm_obs"]
    arm_param.ref_name = "armdiags"
    arm_param.test_data_path = d["arm_test"]
    arm_param.test_name = "20210122.F2010.armsites"
    arm_param.run_type = "model_vs_obs"
    arm_param.test_start_yr = "0001"
    arm_param.test_end_yr = "0001"
    arm_param.ref_start_yr = "0001"
    arm_param.ref_end_yr = "0001"

    runner.sets_to_run = [
        "lat_lon",
        "zonal_mean_xy",
        "zonal_mean_2d",
        "polar",
        "cosp_histogram",
        "meridional_mean_2d",
        "enso_diags",
        "qbo",
        "area_mean_time_series",
        "diurnal_cycle",
        "streamflow",
        "arm_diags",
    ]
    runner.run_diags(
        [param, enso_param, qbo_param, ts_param, dc_param, streamflow_param, arm_param]
    )

    return param.results_dir


if __name__ == "__main__":
    # Choose the `run` function based on what machine you're on.
    # Change <username>

    # Results will be at https://compy-dtn.pnl.gov/<username>/v2_4_0_all_sets/viewer/
    # run_compy('/compyfs/www/<username>/')

    # Results will be at https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/<username>/v2_4_0_all_sets/viewer/
    # run_lcrc('/lcrc/group/e3sm/public_html/diagnostic_output/<username>/')

    # Results will be at https://portal.nersc.gov/project/e3sm/<username>/v2_4_0_all_sets/viewer/
    run_nersc("/global/cfs/cdirs/e3sm/www/<username>/")
