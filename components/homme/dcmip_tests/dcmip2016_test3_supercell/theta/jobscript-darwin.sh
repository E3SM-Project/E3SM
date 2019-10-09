#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 3 on a mac running Darwin
#
# usage: ./jobscript-...
# jobscript-darwin          #run in nonhydrostatic mode
# jobscript-darwin true     #run in hydrostatic mode

hydrostatic="${1:-false}"
echo "theta_hydrostatic_mode=.${hydrostatic}."
# 4dg resolution
EXEC=../../../test_execs/theta-nlev40/theta-nlev40
#sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-r400.nl >& input.nl
sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-animation.nl >& input.nl

openmpiexec -n 6 $EXEC < ./input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl
