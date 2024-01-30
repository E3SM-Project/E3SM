# Run E3SM Diags using zppy on NERSC


1. If simulation is archived with zstash and stored on NERSC HPSS, please first extract the data using zstash. An example usage is provided in following script ``zstash_extract_E3SM_data_from_NERSC_HPSS.sh``. To run the script:

bash zstash_extract_E3SM_data_from_NERSC_HPSS.sh

2. Generate a configuration file for zppy to run a complete set of E3SM Diags. An example is provided in ``zppy_config_for_complete_e3sm_diags_run_on_NERSC.cfg``. To run the script:

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
zppy -c zppy_config_for_complete_e3sm_diags_run_on_NERSC.cfg

