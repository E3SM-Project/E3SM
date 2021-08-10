**********************
Input Data Requirement
**********************


For **versions ealier than v1.6.0**, e3sm_diags runs with seasonal and annual mean climatology datasets generated using NCO. The E3SM output on native grid needs to be regridded/remapped first and split to climo files. (The program can not read h0 files directly) For instance, the usual way to generate climatology files (produce monthly + seasonal + annual climos from monthly input files) is to use ncclimo (a NCO operator) as follows:

 ``ncclimo -s start_yr -e end_yr -c run_id -i drc_in -o drc_out -r map_fl -O drc_rgr -a sdd --no_amwg_links`` for EAM output.

A short summary of the most common options is:

::

    -a: type of DJF average. Either -a scd (default) or -a sdd. scd means seasonally continuous December. The first month used will be Dec of the year before the start year you specify with -s. sdd means seasonally discontinuous December. The first month used will be Jan of the specified start year.
    -C: Climatology mode. Either "mth" (default, and for monthly input files) or "ann" (for annual input files). 
    -c: caseid, i.e., simulation name. For input files like famipc5_ne30_v0.3_00001.cam.h0.1980-01.nc, specify "-c famipc5_ne30_v0.3_00001". The ".cam." and ".h0." bits are added to the filenames internally by default, and can be modified via the "-m mdl_nm" and "-h hst_nm" switches if needed. See comments in ncclimo for documentation. 
    -e: end year (example: 2000). Unless the optional flag "-a sdd" is specified, the last month used will be Nov of the specified end year. If "-a sdd" is specified, the last month will be Dec of the specified end year.
    -h: history file volume that separates the model name from the date in the input file name. Default is "h0".  Other common values are "h1" and "h". 
    -i: directory containing all netcdf files to be used for input to this code.
    -m: model type. Default is "cam". Other options are "clm2", "ocn", "ice", "cism", "cice", "pop".
    -o: directory where computed native grid climo files will be placed. Regridded climos will also be placed here unless a separate directory for them is specified with -O (NB: capital "O") 
    -O: directory where regridded climo files will be placed.
    -s: start year (example: 1980). The first month used will be Dec of the year before the start year you specify (example Dec 1979 to allow for contiguous DJF climos). If "-a sdd" is specified, the first month used will be Jan of the specified start year.
    -v: variable list, e.g., FSNT,AODVIS,PREC.? (yes, regular expressions work so this expands to PRECC,PRECL,PRECSC,PRECSL)

Please see instructions on
`Generate, Regrid, and Split Climatologies (climo files) with ncclimo and ncremap <https://acme-climate.atlassian.net/wiki/spaces/SIM/pages/31129737/Generate+Regrid+and+Split+Climatologies+climo+files+with+ncclimo+and+ncremap>`_. The NCO documentation can be found `here <http://nco.sourceforge.net/nco.html#ncclimo>`_.

To be readable by e3sm_diags, the resulting climatology filenames should remain their NCO sytle, i.e.,starting with simulation name and the season name, e.g., ``20161118.beta0.F1850COSP.ne30_ne30.edison_ANN_climo.nc``. 

**Note**, the program searches files based on simulation name and season name, the input directory (test_data_path) should only contain the re-gridded/post_processed data for one run. If you have multiple climatology files generated with the same simulation (i.e., climo over different time period), please place each set of climatology files in individual folders. 


**Beginning version v1.6.0**, the software can also run with time-series datasets (i.e., multiple years included in one file for each variable). Before running the program, the raw model output should be ran through NCO, which created the time-series files. Below is an example to reshape a series of input files into outputs that are continuous timeseries of each variable taken from all input files.

::

    drc_in=/scratch2/scratchdirs/golaz/ACME_simulations/20161117.beta0.A_WCYCL1850S.ne30_oEC_ICG.edison/run
    map_fl=${DATA}/maps/map_ne30np4_to_fv129x256_aave.20150901.nc
    
    # Read list from file
    ls $drc_in/*cam.h0.0[012]??* > input_list
    ncclimo --dbg=0 --yr_srt=1 --yr_end=250 --var=FSNT,AODVIS --map=$map_fl --drc_out=$drc_out --drc_rgr=$drc_rgr < input_list
    # Pipe list to stdin
    cd $drc_in
    ls *cam.h0.0[012]??* | ncclimo --dbg=0 --yr_srt=1 --yr_end=250 --var=FSNT,AODVIS --map=$map_fl --drc_out=$drc_out --drc_rgr=$drc_rgr
    # List as positional arguments
    ncclimo --var=FSNT,AODVIS --yr_srt=1 --yr_end=250 --map=$map_fl --drc_out=$drc_out --drc_rgr=$drc_rgr $drc_in/*cam.h0.0[012]??*.nc
    # Read directory
    ncclimo --var=T,Q,RH --yr_srt=1 --yr_end=250 --drc_in=$drc_in --map=$map_fl --drc_out=$drc_out --drc_rgr=$drc_rgr

Note that ``map_ne30np4_to_fv129x256_aave.20150901.nc`` is an example mapping files we are using. For more mapping files available for e3sm, please refer to `<https://web.lcrc.anl.gov/public/e3sm/mapping/grids/>`_.

The output is a collection of per-variable timeseries such as ``FSNT_YYYYMM_YYYYMM.nc``, ``AODVIS_YYYYMM_YYYYMM.nc``, etc. The output is split into segments each containing no more than ypf_max (default 50) years-per-file, e.g., ``FSNT_000101_005012.nc``, ``FSNT_005101_009912.nc``, ``FSNT_010001_014912.nc``, etc. 

If you are using time series from CMIP style files, the model data file names must follow the naming conventions as follows, where you have
``<variable>_<start_yr>01_<end_yr>12.nc``. Ex: renaming ``tas_Amon_CESM1-CAM5_historical_r1i2p1_196001-201112.nc`` to ``tas_196001_201112.nc``.

All of the variables should be in the same directory. Please refer to the test data format available on NERSC (/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/) for examples.


Preprocessing Steps
^^^^^^^^^^^^^^^^^^^

Some diagnostic sets require preprocessing steps. See
`this script <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/tutorials/postprocessing_E3SM_data_for_e3sm_diags.sh>`_
for postprocessing E3SM data (i.e., preprocessing for E3SM Diags).

- QBO: monthly ``U`` (if different start_year and end_year are wanted for test/ref,
  the number of years need to remain the same)

- Streamflow: monthly ``RIVER_DISCHARGE_OVER_LAND_LIQ`` and ``areatotal2``.
  Example for splitting mosart timeseries files using latest NCO
  (NCO 4.9.6 includes the ``--var_xtr`` option, which can be used to specify the extra variables to include in
  timeseries that the splitter produces):

    ::

        cd ${drc_in}; ls *mosart.h0.20*nc | ncclimo --var_xtr=areatotal2 -v RIVER_DISCHARGE_OVER_LAND_LIQ --yr_srt=2000 --yr_end=2014 --drc_in=$drc_in --drc_out=$drc

- Diurnal Cycle: seasonal/annual mean diurnal cycle climo
  Example for Compy:
    ::

        #!/bin/bash

        source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

        # Low-res Compy simulations
        drc_in=/compyfs/zhen797/E3SM_simulations/20201027.alpha5_v1p-1.amip.ne30pg2_r05_oECv3.compy/archive/atm/hist

        drc_out=/qfs/people/zhan429/postprocess_scripts/diurnal_climo/native
        drc_rgr=/qfs/people/zhan429/postprocess_scripts/diurnal_climo/rgr
        caseid=20201027.alpha5_v1p-1.amip.ne30pg2_r05_oECv3.compy.eam.h4

        cd ${drc_in};ls ${caseid}*{1980..2014}*.nc | ncclimo --clm_md=hfc --caseid=$caseid -v PRECT --ypf=1 --yr_srt=1980 --yr_end=2014 --drc_out=${drc_out} -O $drc_rgr --map=/qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc

