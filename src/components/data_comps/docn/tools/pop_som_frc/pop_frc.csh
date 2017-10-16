#/bin/csh -f

# This shell script creates a seasonal cycle climatology over a specified
# number of years. The pop history files are quite large, so each month is
# read for all years from the mass store, the climatology for that month
# is created and the history files are deleted. The coupler history files
# are much smaller and can all be read off the mass store at once.

set echo on

# Need to set the CASE name for the run to extract POP and CPL history
# files from and the start and end years of the period for the climatology.

setenv CASE b40.1850.track1.2deg.003

set BEGYR = 481
set ENDYR = 500

# Set the mass store userid where the files are located. If the run did
# not save the CPL history files, set CPLFILES to FALSE. If you have plenty
# of local space for the POP history files, set TONS_O_SPACE to TRUE.

set CCSMUSER = CCSM
set CPLFILES = FALSE
set TONS_O_SPACE = TRUE

set SCRIPT_HOME = ${HOME}/pop_frc
setenv PATH_MSS /${CCSMUSER}/csm/${CASE}/ocn/hist/

set DATE_FORMAT = 'yyyy-mm'

setenv WKDIR /biptmp/${USER}/${CASE}
setenv FILE_HEADER ${CASE}.pop.h.

mkdir -p $WKDIR
cd $WKDIR

if !( -e ${CASE}.pop.h.${BEGYR}-${ENDYR}.MAC.nc ) then

# If you have plenty of local disk space for the POP history files.
if ( $TONS_O_SPACE == 'TRUE' ) then
   ${SCRIPT_HOME}/read_from_mss.csh $DATE_FORMAT $CASE $BEGYR $ENDYR
endif

foreach MONTH ( 01 02 03 04 05 06 07 08 09 10 11 12 )

# Otherwise if you have limited local disk space, just read a single
# month from each year and delete.
if !( $TONS_O_SPACE == 'TRUE' ) then
   ${SCRIPT_HOME}/read_from_mss_month.csh $DATE_FORMAT $CASE $BEGYR $ENDYR $MONTH
endif

   ncra -O -vTEMP,SALT,UVEL,VVEL,SHF,QFLUX,MELTH_F,RESID_T,HBLT,REGION_MASK,TAREA,TLONG,TLAT,ANGLET *-${MONTH}.nc ${MONTH}.nc

if !( $TONS_O_SPACE == 'TRUE' ) then
   /bin/rm -f *-??.nc
endif

end

ncrcat -O ??.nc ${CASE}.pop.h.${BEGYR}-${ENDYR}.MAC.nc

/bin/rm -f ??.nc

endif

if ( $CPLFILES == TRUE ) then

setenv PATH_MSS /${CCSMUSER}/csm/${CASE}/cpl/hist/
setenv FILE_HEADER ${CASE}.cpl6.ha.

mkdir -p $WKDIR
cd $WKDIR

if !( -e ${CASE}.cpl6.ha.${BEGYR}-${ENDYR}.MAC.nc ) then

${SCRIPT_HOME}/read_from_mss.csh $DATE_FORMAT $CASE $BEGYR $ENDYR

ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-01.nc 01.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-02.nc 02.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-03.nc 03.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-04.nc 04.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-05.nc 05.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-06.nc 06.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-07.nc 07.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-08.nc 08.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-09.nc 09.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-10.nc 10.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-11.nc 11.nc
ncra -O -v avXc2i_i_So_u,avXc2i_i_So_v,avXc2i_i_So_dhdx,avXc2i_i_So_dhdy *-12.nc 12.nc

ncrcat -O ??.nc ${CASE}.cpl6.ha.${BEGYR}-${ENDYR}.MAC.nc
/bin/rm -f *-??.nc ??.nc

endif
endif

end
