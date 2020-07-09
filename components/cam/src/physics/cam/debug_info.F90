module debug_info

  !-------------------------------------------------------------
  ! This is a light weight module to report parameters which
  ! might be of interest in case of a simulation crash.
  !
  ! Author: Balwinder Singh
  ! (Following suggestions from Peter Caldwell and Aaron Donahue)
  !-------------------------------------------------------------

  implicit none
  private

  !Variables updated during the simulation
  integer :: macmiciter, chunk, icol_p3
  !$OMP THREADPRIVATE(macmiciter,chunk, icol_p3 )

  public:: &
       get_debug_chunk,      &  !update chunk info from the simulation
       get_debug_column_id,  &  !update column id
       get_debug_macmiciter, &  !update macro/microphysics iteration count
       report_error_info        !report info on the current lat/lon and other simulation parameters

contains

  subroutine get_debug_chunk(chunk_in)

    !Obtain chunk # from tphycbc

    implicit none

    integer, intent(in) :: chunk_in

    chunk = chunk_in

  end subroutine get_debug_chunk

  subroutine get_debug_column_id(col_in)

    !Obtain column id from p3_main (micro_p3.F90)

    implicit none

    integer, intent(in) :: col_in

    icol_p3 = col_in

  end subroutine get_debug_column_id

  subroutine get_debug_macmiciter(iter)

    !Obtain iteration # of macro/microphysics sub-stepping iteration

    implicit none

    integer, intent(in) :: iter

    macmiciter = iter

  end subroutine get_debug_macmiciter


  subroutine report_error_info(cause, calling_func, icol, klev, prnt_macmic)

    !Reports simulation parameters
    use cam_logfile,    only: iulog
    use phys_grid,      only: get_rlat_p, get_rlon_p
    use time_manager,   only: get_nstep
    use physics_utils,  only: rtype
    use spmd_utils,     only: iam
    implicit none

    !intent-ins
    character(len=*),  intent(in) :: cause, calling_func

    !intent-in and optional
    integer, optional, intent(in) :: icol !column(lat/lon) index
    integer, optional, intent(in) :: klev !vertical index
    logical, optional, intent(in) :: prnt_macmic !whether print macmiciter or not

    !local variables
    real(rtype), parameter :: radian2deg = 180.0_rtype/(4.0_rtype*atan(1.0_rtype)) !180/pi = 57.296
    real(rtype)            :: lat, lon
    character(len=1000)    :: klevstr, macmicstr, icolstr

    !check if klev is present or not and populate klevstr accordingly
    if(present(klev)) then
       write(klevstr,*)'K level is:',klev
    else
       klevstr = 'K level not present'
    endif

    !check if prnt_macmic is present or not and populate macmicstr accordingly
    if(present(prnt_macmic) .and. prnt_macmic) then
       write(macmicstr,*)'macmic_it is:', macmiciter
    else
       macmicstr = 'prnt_macmic is not present or prnt_macmic is False '
    endif

    !update icolstr to add info about the column used and compute lat lon
    if(present(icol)) then
       icolstr = '*NOTE*: Using user specified column index for computing lat-lons'
       lat = get_rlat_p(chunk,icol)*radian2deg
       lon = get_rlon_p(chunk,icol)*radian2deg
    else
       icolstr = '*NOTE*: Using auto updated column index from P3 for computing lat-lons'
       lat = get_rlat_p(chunk,icol_p3)*radian2deg
       lon = get_rlon_p(chunk,icol_p3)*radian2deg
    endif

    !print proc info in the message and issue only one write statement
    write(iulog,'(a,/,a,i5,a,/,a,/,a,/,a,f12.8,a,f12.8,a,/,a,/,a,i5,/,a,/,a,i5,a,/,a)') &
         '', &
         '********************** Proc #', iam,' output *******************************************', &
         'Fail due to '//trim(adjustl(cause))//' in subroutine/function:'//trim(adjustl(calling_func)), &
         trim(adjustl(icolstr)), &
         'Column location in degrees(lat,lon) is:(',lat,',',lon,')', &
         trim(adjustl(macmicstr)), &
         'timestep:',get_nstep(), &
         trim(adjustl(klevstr)), &
         '********************** Proc #', iam,' output Ends **************************************', &
         ''
  end subroutine report_error_info


end module debug_info
