Setting up time-varying forcing for MPAS-O
==========================================
Implementation
---------------
 https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F 

API
---
 https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/framework/mpas_forcing.F 

Example datasets
----------------
 https://zenodo.org/record/1219191#.W9jWCKlRcUF with example from this dataset::

    > pwd 
    /Users/pwolfram/Downloads/MPAS-Seaice_test_dataset_V1/domain_QU240km

    > ncdump -h LYq_six_hourly.2000.nc
    netcdf LYq_six_hourly.2000 {
    dimensions:
        nCells = 1794 ;
        StrLen = 64 ;
         Time = UNLIMITED ; // (1460 currently)
    variables:
        char xtime(Time, StrLen) ;
        double airTemperature(Time, nCells) ;
        double airSpecificHumidity(Time, nCells) ;
        double uAirVelocity(Time, nCells) ;
        double vAirVelocity(Time, nCells) ;
    }

    > ncdump -v xtime LYq_six_hourly.2000.nc | tail
      "2000-12-29_24:00:00                                             ",
      "2000-12-30_06:00:00                                             ",
      "2000-12-30_12:00:00                                             ",
      "2000-12-30_18:00:00                                             ",
      "2000-12-30_24:00:00                                             ",
      "2000-12-31_06:00:00                                             ",
      "2000-12-31_12:00:00                                             ",
      "2000-12-31_18:00:00                                             ",
      "2000-12-31_24:00:00                                             " ;
    }

Data structure
--------------
1. Global variables needed: `type (MPAS_forcing_group_type), pointer :: seaiceForcingGroups`

2. Register against the `ForcingGroups` using the following initialization steps.  `ForcingGroups` is a pointer to a list of forcings, e.g., see  initialized group in https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L167 and members of the group https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L180 and https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L193 

3. Streams file defines these source inputs, use example list of steps to get the data into a registry defined variable field / array.

Initialization
--------------
Example of init steps (https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L232) 


1. Create forcing group via `MPAS_forcing_init_group` https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L167 

2. Initialize fields via `MPAS_forcing_init_field` https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L180 

3. Initialize the data via reading the file with `MPAS_forcing_init_field_data`: https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L232 

In time-step usage
------------------
Example of usage steps:

1. `MPAS_forcing_get_forcing` Loop over all the individual forcing fields in the forcing group and get the data and perform the time interpolation. https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L389 

2. `MPAS_forcing_get_forcing_time` Return the current forcing clock time for a forcing group. https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L517 

3. `MPAS_forcing_write_restart_times` Loop over the forcing groups in the forcing group object and write out the forcing clock times to registry variables that are included in the restart stream. 'forcingGroupHead' is the forcing group object https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/shared/mpas_seaice_forcing.F#L2080 

Registry changes needed (for restarts to work)
----------------------------------------------
Registry variables needed for restart:

1. `nForcingGroupMax` https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/Registry.xml#L287 

2. `nForcingGroupCounter` https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/Registry.xml#L3295 

3. `forcingGroupNames` https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/Registry.xml#L3297  

4. `forcingGroupRestartTimes` https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/Registry.xml#L3298 

Register to restart streams:

1. https://github.com/MPAS-Dev/MPAS-Model/blob/seaice/develop/src/core_seaice/Registry.xml#L1798
  
 | `<var name="forcingGroupNames"/>` 
 | `<var name="forcingGroupRestartTimes"/>`

Additional resources
--------------------
M Maltrud PRs using this information: https://github.com/MPAS-Dev/MPAS/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Aclosed+author%3Amaltrud
Old (depricated-- recommended to NOT USE) design doc: https://docs.google.com/document/d/1QtjmVCxLKS-S9Z_X8-WiHqXGn_yDu_VFHLA93ZZHg88/edit

Time estimate (MPAS / MPAS-O novice)
From start (this document) to merged PR:
2 solid weeks
4 weeks half-time
