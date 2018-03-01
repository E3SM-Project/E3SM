#!/bin/csh -f

#-------------------------------------------------------------------------------
# post process POP hist output, doing the following
#
# if ($OCN_TAVG_HIFREQ == TRUE)
#    1) Recreate monthly means for fields that have been moved to daily mean stream.
#    2) Remove fields from daily mean stream that are unnecessary.
#    For this processing, it is assumed
#       1) daily fields are bundled into monthly aggregates
#          (this is default behavior for nday1 output)
#       2) monthly aggregates files of daily fields contain all days of output
#       3) default OCN_TAVG_HIFREQ placement of variables into streams is being used
#
# Generally, it is assumed that
#    1) Model output resides on disk and that this script is called from the
#       directory where the history files to be processed reside.
#    2) $CASEROOT/Tools/ccsm_getenv has been called.
#    3) DOUT_L_HTAR is FALSE. (tarring is not handled at all)
#
# Author:
#   Keith Lindsay August 2010 
#
# Revision History:
#   Nancy Norton  31 August 2010 define and add my_pid to temp filenames to avoid possible collisions
#-------------------------------------------------------------------------------


# parse arguments passed to script
while ( 1 )
   if ( $#argv < 1 ) break
   switch ( $argv[1] )
      default:
         echo unknown argument $argv[1]
         echo "usage: $0"
         exit -1
   endsw
   shift argv
end

set my_pid = $$

if ($OCN_TAVG_HIFREQ == TRUE) then
   set missing_mon_fields = moc_components,EVAP_F,HBLT,IFRAC,IOFF_F,LWDN_F,LWUP_F,MELTH_F,MELT_F,MOC,PREC_F,QFLUX,QSW_HBL,QSW_HTP,RHO_VINT,ROFF_F,SALT_F,SENH_F,SFWF,SFWF_WRST,SHF,SHF_QSW,SNOW_F,SSH,TAUX,TAUX2,TAUY,TAUY2,UVEL,VISOP,VSUBM,VVEL,WISOP,WSUBM,WVEL,XBLT

   set day_fields_to_remove = VISOP,WISOP,VSUBM,WSUBM

   # determine which . seperated field from filename contains date string
   set case_fields = `echo ${CASE} | awk -F. '{print NF}'`
   @ nday1_date_field = ${case_fields} + 4

   foreach yyyy_mm_dd ( `ls | grep "$CASE.pop.h.nday1" | awk -F. '{print $'${nday1_date_field}'}'` )
      set yyyy_mm = `echo $yyyy_mm_dd | cut -d'-' -f1,2`

      set file_day = $CASE.pop.h.nday1.$yyyy_mm_dd.nc
      set file_mon = $CASE.pop.h.$yyyy_mm.nc

      date

      # Recreate monthly means for fields that have been moved to day
      # mean stream. The computed monthly means are initially stored
      # in an intermediate file. This is done so that global metadata
      # that appears in both the daily and monthly files, as well as the
      # time variable, are not overwritten in the monthly file.
      set file_mon_temp = $CASE.pop.h.$yyyy_mm.${my_pid}.temp.nc

      set field_name = `echo $missing_mon_fields | cut -f1 -d,`
      ncdump -h $file_mon | grep $field_name > /dev/null
      if ($status) then
         echo recreating monthly means in $file_mon from $file_day
         ncra -O -v $missing_mon_fields -h $file_day $file_mon_temp
         ncatted -a start_time,global,d,, \
                 -a nsteps_total,global,d,, \
                 -a tavg_sum,global,d,, \
                 -a tavg_sum_qflux,global,d,, \
                 -h $file_mon_temp
         ncks -A -v time -h $file_mon $file_mon_temp
         ncks -A -h $file_mon_temp $file_mon
         rm $file_mon_temp
      else
         echo first field of missing_mon_fields detected in $file_mon, skipping monthly mean recreation
      endif

      # remove fields from daily mean stream that are unnecessary
      set field_name = `echo $day_fields_to_remove | cut -f1 -d,`
      ncdump -h $file_day | grep $field_name > /dev/null
      if !($status) then
         echo removing unnecessary fields from $file_day
         ncks -O -v $day_fields_to_remove -x -h $file_day $file_day
      else
         echo first field of day_fields_to_remove not detected in $file_day, skipping daily mean field removal
      endif

   end
endif # $OCN_TAVG_HIFREQ == TRUE
