More on Driver Namelists
=========================

There are a series of driver/coupler namelist input files created by the driver namelist generator ``$CIMEROOT/driver_cpl/cime_config/buildnml``. These are

- drv_in
- drv_flds_in
- cpl_modelio.nml, atm_modelio.nml, esp_modelio.nml, glc_modelio.nml, ice_modelio.nml, lnd_modelio.nmlo, ocn_modelio.nml, rof_modelio.nml, wav_modelio.nml
- seq_maps.rc

The ``*_modelio.nml`` files set the filename for the primary standard output file and also provide settings for the parallel IO library, PIO. 
The drv_in namelist file contains several different namelist groups associated with general options, time manager options, pe layout, timing output, and parallel IO settings.
The seq_maps.rc file specifies the mapping files for the configuration.
Note that seq_maps.rc is NOT a Fortran namelist file but the format should be relatively clear from the default settings.
