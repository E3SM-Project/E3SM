#!/bin/csh -f

source /etc/profile.d/00softenv.csh
soft add +mpiwrapper-xl
soft add @ibm-compilers-2014-02
soft add +cmake
setenv PNETCDF /soft/libraries/pnetcdf/1.5.0/cnk-xl/V1R2M2-20140827
setenv NETCDF /soft/libraries/netcdf/4.3.0-f4.2/cnk-xl/V1R2M0-20131211
