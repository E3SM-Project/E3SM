#! /bin/csh -f
#
# Create map grid files for all resolutions
#

set cdate="c"`date +%y%m%d`

# Currently ESMF regridding does NOT work for single-point files
# Although it does for regional files of at least 4 grid points
set resols = ( "0.9x1.25" "0.23x0.31" "1.9x2.5" "128x256" \
               "2.5x3.33" "32x64" "48x96" "4x5" "5x5_amazon" \
               "1x1_brazil" "1x1_camdenNJ" "1x1_urbanc_alpha" "1x1_mexicocityMEX" \
               "1x1_vancouverCAN" "1x1_tropicAtl" \
               "512x1024" "64x128" "8x16" "94x192" "10x15" \
             )
set query="../../bld/queryDefaultNamelist.pl -silent -onlyfiles "
set query="$query -justvalue -csmdata $CSMDATA -var fatmgrid "
set mask="nomask"
foreach res ( $resols )
   set griddata = `$query -res $res -options mask=$mask`
   cat >! namelist <<EOF
 &mkmapgrids_in
  fname_in  = '$griddata'
  fname_out = 'SCRIPgrid_${res}_${mask}_${cdate}.nc'
/
EOF
  echo "grid $res"
  cat mkmapgrids_in
  mkmapgrids < namelist || exit -1
end

