program regional

integer                                :: nr
character(len=60),dimension(25)        :: rlist
character(len=1)                       :: smart_flow_option,converted_forest_land_option 
character(len=1)                       :: zdis_option,adjust_smart_flow_option,harvest_option
character(len=2)                       :: num_reg


call getarg(1,smart_flow_option)
call getarg(2,converted_forest_land_option)
call getarg(3,zdis_option)
call getarg(4,adjust_smart_flow_option)
call getarg(5,harvest_option)
call getarg(6,num_reg)

if (num_reg .eq. '14') then
   nr=14
   open(14,file='rlist14.in',status='old')
else if (num_reg .eq. '24') then
   nr=24
   open(14,file='rlist24.in',status='old')
endif

open(16,file='rlist.txt',status='unknown')

!* output order of options to clist.txt
! smart_flow_option, converted_forest_land_option, zdis_option, adjust_smart_flow_option,
! harvest_option

do i=1,nr

  read(14,'(a60)') rlist(i)
  print *, rlist(i)
  write(16,'(5(a1,1x))') smart_flow_option,converted_forest_land_option,zdis_option, &
                         adjust_smart_flow_option,harvest_option
 
end do
close(14)
close(16)

stop
end
