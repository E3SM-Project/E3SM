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
  integer :: macmiciter, chunk
  !$OMP THREADPRIVATE(macmiciter,chunk )

  public:: &
       get_debug_chunk, &       !update chunk info from the simulation
       get_debug_macmiciter, &  !update macro/microphysics iteration count
       report_error_info        !report info on the current lat/lon and other simulation parameters

contains

  subroutine get_debug_chunk(chunk_in)

    !Obtain chunk # from tphycbc

    implicit none

    integer :: chunk_in

    chunk = chunk_in

  end subroutine get_debug_chunk


  subroutine get_debug_macmiciter(iter)

    !Obtain iteration # of macro/microphysics sub-stepping iteration

    implicit none

    integer, intent(in) :: iter

    macmiciter = iter

  end subroutine get_debug_macmiciter


  subroutine report_error_info(cause, calling_func, icol, klev)

    !Reports simulation parameters
    use cam_logfile,    only: iulog
    use phys_grid,      only: get_rlat_p, get_rlon_p
    use time_manager,   only: get_nstep
    use physics_utils,  only: rtype
    implicit none

    integer, intent(in) :: icol !column(lat/lon) index
    integer, optional, intent(in) :: klev !vertical index
    character(len=*), intent(in) :: cause, calling_func

    !local variables
    real(rtype), parameter::  radian2deg = 180.0_rtype/(4.0_rtype*atan(1.0_rtype)) !180/pi = 57.296

    write(iulog,'(a)') '*****************************************************************'
    write(iulog,'(a)')'Fail due to '//trim(adjustl(cause))//' in subroutine/function:'//trim(adjustl(calling_func))
    write(iulog,'(a,f12.8,a,f12.8,a)')'Column location in degrees:(lat,lon) is:(',get_rlat_p(chunk,icol)*radian2deg,',', &
         get_rlon_p(chunk,icol)*radian2deg,')'
    write(iulog,'(a,i5)')'macmic_it is:', macmiciter
    write(iulog,'(a,i5)')'timestep=',get_nstep()

    if(present(klev)) then
       write(iulog,'(a,i5)') 'K level is:', klev
    endif

  end subroutine report_error_info


end module debug_info
