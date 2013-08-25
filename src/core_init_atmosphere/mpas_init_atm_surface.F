! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_init_atm_surface
 use mpas_configure
 use mpas_grid_types
 use mpas_io_output
 use mpas_timekeeping
 use mpas_timer
 
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_read_met

 implicit none
 private
 public:: init_atm_case_sfc,interp_sfc_to_MPAS

 integer, parameter :: R4KIND = selected_real_kind(6)

 contains

!==================================================================================================
 subroutine init_atm_case_sfc(domain,dminfo,mesh,fg,state)
!==================================================================================================

!input arguments:
 type(domain_type), intent(inout):: domain
 type(dm_info), intent(in)       :: dminfo
 type(mesh_type), intent(inout)  :: mesh
 type(fg_type), intent(inout)    :: fg
 type (state_type), intent(inout):: state

!local variables:
 type(MPAS_Clock_type)       :: fg_clock
 type(MPAS_Time_type)        :: start_time,stop_time,curr_time
 type(MPAS_TimeInterval_type):: fg_interval

 type(io_output_object):: sfc_update_obj

 character(len=StrKIND) :: timeString

!==================================================================================================

!set up clock to step through all intermediate file dates to be processed:
 call mpas_set_time(start_time, dateTimeString=trim(config_start_time))
 call mpas_set_time(stop_time, dateTimeString=trim(config_stop_time))
 call mpas_set_timeInterval(fg_interval, S=config_fg_interval)
 call mpas_create_clock(fg_clock, start_time, fg_interval, stopTime=stop_time) 

!initialize the output file
 sfc_update_obj % time = 1
 sfc_update_obj % filename = trim(config_sfc_update_name)

 call mpas_output_state_init(sfc_update_obj, domain, "SFC")

!loop over all times:
 curr_time = mpas_get_clock_time(fg_clock, MPAS_NOW) 

 do while (curr_time <= stop_time)
    call mpas_get_time(curr_time, dateTimeString=timeString)
!   write(0,*) 'Processing ',trim(config_sfc_prefix)//':'//timeString(1:13)

    !read the sea-surface temperature and sea-ice data from the surface file, and interpolate the
    !data to the MPAS grid:
    call interp_sfc_to_MPAS(timeString(1:13),mesh,fg,dminfo)

    !write the interpolated SST/SKINTEMP field as a new time slice in the MPAS output file:
    call mpas_output_state_for_domain(sfc_update_obj, domain, sfc_update_obj % time)
    sfc_update_obj % time = sfc_update_obj % time + 1

    call mpas_advance_clock(fg_clock)
    curr_time = mpas_get_clock_time(fg_clock, MPAS_NOW) 

    call mpas_get_time(curr_time, dateTimeString=timeString)
    state % xtime % scalar = timeString

 enddo

 call mpas_output_state_finalize(sfc_update_obj, dminfo)
      
 end subroutine init_atm_case_sfc

!==================================================================================================
 subroutine interp_sfc_to_MPAS(timeString,mesh,fg,dminfo)
!==================================================================================================

!input arguments:
 character(len=*),intent(in):: timeString
 type(mesh_type), intent(in):: mesh
 type(dm_info),intent(in)   :: dminfo

!inout arguments:
 type(fg_type), intent(inout):: fg

!local variables:
 type(met_data) :: field !real*4 meteorological data.

 integer:: istatus
 integer:: masked
 integer,dimension(5):: interp_list
 integer,dimension(:),pointer:: mask_array
 logical :: have_landmask

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND),dimension(:,:),allocatable:: maskslab

 real(kind=RKIND), dimension(:), pointer:: destField1d

!==================================================================================================
 mask_array => mesh % landmask % array

!open intermediate file:
 call read_met_init(trim(config_sfc_prefix),.false.,timeString,istatus)
 if(istatus /= 0) then
    write(0,*) '********************************************************************************'
    write(0,*) 'Error opening surface file '//   &
                trim(config_sfc_prefix)//':'//timeString(1:13)
    write(0,*) '********************************************************************************'
    call mpas_dmpar_abort(dminfo)
 else
    write(0,*) 'Processing file ',trim(config_sfc_prefix)//':'//timeString(1:13)
 endif

!scan through all fields in the file, looking for the LANDSEA field:
 have_landmask = .false.
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
    if(index(field % field, 'LANDSEA') /= 0) then
       have_landmask = .true.
       if(.not.allocated(maskslab)) allocate(maskslab(-2:field % nx+3, field % ny))
       maskslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
       maskslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
       maskslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
       maskslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
       maskslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
       maskslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
       maskslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)
!      write(0,*) 'minval, maxval of LANDSEA = ', minval(maskslab), maxval(maskslab)
    endif
    deallocate(field % slab)
    call read_next_met_field(field,istatus)
 enddo
 call read_met_close()

!read sea-surface temperatures and seaice data. open intermediate file:
 call read_met_init(trim(config_sfc_prefix),.false.,timeString(1:13),istatus)
 if(istatus /= 0) then
    write(0,*) '********************************************************************************'
    write(0,*) 'Error opening surface file '//   &
                trim(config_sfc_prefix)//':'//timeString(1:13)
    write(0,*) '********************************************************************************'
    call mpas_dmpar_abort(dminfo)
 endif

 if(.not. have_landmask) then
    write(0,*) '********************************************************************************'
    write(0,*) 'Landsea mask not available from the surface file '
    write(0,*) '********************************************************************************'
 endif

!scan through all fields in the file, looking for the SST,SKINTEMP, or SEAICE field:
 call read_next_met_field(field,istatus)
 do while (istatus == 0)

    !sea-surface data:
    if(index(field % field, 'SKINTEMP') /= 0 .or. index(field % field, 'SST') /= 0) then
!      write(0,*) '... Processing SST:'
       fg % sst % array(1:mesh%nCells) = 0.0_RKIND
       destField1d => fg % sst % array

       !interpolation to the MPAS grid:
       interp_list(1) = FOUR_POINT
       interp_list(2) = SEARCH
       interp_list(3) = 0
       msgval  = -1.0e30_R4KIND !missing value
       masked  = -1
       maskval = -1.0_RKIND
       fillval =  0.0_RKIND
       if(have_landmask) then
          call interp_to_MPAS(mesh,field,destField1d,interp_list,msgval,masked,maskval,fillval, &
                              mask_array,maskslab)
       else
          call interp_to_MPAS(mesh,field,destField1d,interp_list,msgval,masked,maskval,fillval, &
                              mask_array)
       endif

       !field%slab was allocated in the subroutine read_next_met_field
       deallocate(field%slab)

    !sea-ice data:
    elseif(index(field % field, 'SEAICE') /= 0) then
!      write(0,*) '... Processing SEAICE:'
       fg % xice % array(1:mesh%nCells) = 0.0_RKIND
       destField1d => fg % xice % array

       !interpolation to the MPAS grid:
       interp_list(1) = FOUR_POINT
       interp_list(2) = W_AVERAGE4
       interp_list(3) = SEARCH
       interp_list(4) = 0
       msgval  = -1.0e30_R4KIND  !missing value
       masked  = 1
       maskval = 1.0_RKIND
       fillval = 0.0_RKIND
       if(have_landmask) then
          call interp_to_MPAS(mesh,field,destField1d,interp_list,msgval,masked,maskval,fillval, &
                              mask_array,maskslab)
       else
          call interp_to_MPAS(mesh,field,destField1d,interp_list,msgval,masked,maskval,fillval, &
                              mask_array)
       endif

       !field%slab was allocated in the subroutine read_next_met_field
       deallocate(field%slab)
        
    else
       deallocate(field%slab)

    endif

    call read_next_met_field(field,istatus)
 enddo

!close intermediate file:
 call read_met_close()
 if(allocated(maskslab)) deallocate(maskslab)

!freeze really cold oceans:
 where(fg%sst%array.lt.271.0_RKIND .and. mesh%landmask%array.eq.0) fg%xice%array = 1.0_RKIND

!limit XICE to values between 0 and 1. Although the input meteorological field is between 0. and 1.
!interpolation to the MPAS grid can yield values of XiCE less than 0. and greater than 1.:
 where (fg%xice%array < 0._RKIND) fg%xice%array = 0._RKIND
 where (fg%xice%array > 1._RKIND) fg%xice%array = 1._RKIND

 end subroutine interp_sfc_to_MPAS

!==================================================================================================
 subroutine interp_to_MPAS(mesh,field,destField1d,interp_list,msgval,masked,maskval,fillval, &
                           mask_array,maskslab)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(met_data),intent(in) :: field !real*4 meteorological data.

 integer,intent(in):: masked
 integer,dimension(5),intent(in):: interp_list
 integer,dimension(:),intent(in),pointer:: mask_array

 real(kind=RKIND),intent(in):: fillval,maskval,msgval
 real(kind=RKIND),intent(in),dimension(*),optional:: maskslab

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(:),pointer:: destField1d

!local variables:
 type(proj_info):: proj
 integer:: i,nInterpPoints
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:,:),allocatable:: rslab

 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints
 
!--------------------------------------------------------------------------------------------------

 call map_init(proj)   
 if(field % iproj == PROJ_LATLON) then
    call map_set(PROJ_LATLON, proj, &
                 latinc = real(field % deltalat,RKIND), &
                 loninc = real(field % deltalon,RKIND), &
                 knowni = 1.0_RKIND, &
                 knownj = 1.0_RKIND, &
                 lat1 = real(field % startlat,RKIND), &
                 lon1 = real(field % startlon,RKIND))
!   write(0,*) '--- The projection is PROJ_LATLON.'
 elseif(field % iproj == PROJ_GAUSS) then
    call map_set(PROJ_GAUSS, proj, &
                 nlat = nint(field % deltalat), &
                 loninc = real(field % deltalon,RKIND), &
                 lat1 = real(field % startlat,RKIND), &
                 lon1 = real(field % startlon,RKIND))
!   write(0,*) '--- The projection is PROJ_GAUSS.'
 elseif(field % iproj == PROJ_PS) then
    call map_set(PROJ_PS, proj, &
                 dx = real(field % dx,RKIND), &
                 truelat1 = real(field % truelat1,RKIND), &
                 stdlon = real(field % xlonc,RKIND), &
                 knowni = real(field % nx / 2.0,RKIND), &
                 knownj = real(field % ny / 2.0,RKIND), &
                 lat1 = real(field % startlat,RKIND), &
                 lon1 = real(field % startlon,RKIND))
!   write(0,*) '--- The projection is PROJ_PS.'
 endif

 nInterpPoints = mesh % nCells
 latPoints => mesh % latCell % array
 lonPoints => mesh % lonCell % array

 allocate(rslab(-2:field % nx+3, field % ny))
 rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
 rslab( 0, 1:field % ny) = field % slab(field % nx  , 1:field % ny)
 rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
 rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
 rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
 rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
 rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

 do i = 1,nInterpPoints
    if(mask_array(i) /= masked) then
       lat = latPoints(i) * DEG_PER_RAD
       lon = lonPoints(i) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(y < 0.5) then
          y = 1.0
       elseif(y >= real(field%ny)+0.5) then
          y = real(field % ny)
       endif
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       elseif (x >= real(field%nx)+0.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       endif
       if(present(maskslab)) then
          destField1d(i) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1, &
                           msgval,interp_list,1,maskval=maskval,mask_array=maskslab)
       else
          destField1d(i) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1, &
                           msgval,interp_list,1,maskval=maskval)
       endif
    else
       destField1d(i) = fillval
    endif
 enddo
 deallocate(rslab)

 end subroutine interp_to_MPAS

!==================================================================================================
 end module mpas_init_atm_surface
!==================================================================================================

