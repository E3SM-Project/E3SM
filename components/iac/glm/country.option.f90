program country

integer                                :: nc
character(len=60),dimension(192)       :: clist
character(len=1)                       :: smart_flow_option,converted_forest_land_option 
character(len=1)                       :: zdis_option,adjust_smart_flow_option,harvest_option
character(len=1)                       :: nodata_option


call getarg(1,smart_flow_option)
call getarg(2,converted_forest_land_option)
call getarg(3,zdis_option)
call getarg(4,adjust_smart_flow_option)
call getarg(5,harvest_option)
call getarg(6,nodata_option)


if ((nodata_option .eq. '5') .or. (nodata_option .eq. '6')) then
   nc=192
   open(14,file='clist192.in',status='old')
endif

open(16,file='clist.txt',status='unknown')

!* output order of options to clist.txt
! smart_flow_option, converted_forest_land_option, zdis_option, adjust_smart_flow_option,
! harvest_option

do i=1,nc

  read(14,'(a60)') clist(i)

  write(16,'(5(a1,1x))') smart_flow_option,converted_forest_land_option,zdis_option, &
                         adjust_smart_flow_option,harvest_option
 
end do
close(14)
close(16)

stop
end
