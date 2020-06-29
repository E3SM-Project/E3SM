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

    use cam_logfile,  only: iulog
    use phys_grid,    only: get_rlat_p, get_rlon_p
    use time_manager, only: get_nstep

    implicit none

    integer, intent(in) :: icol !column(lat/lon) index
    integer, intent(in) :: klev !vertical index
    character(len=*), intent(in) :: cause, calling_func

    !local variables
    write(iulog,*) 'Fail due to ',cause,' in ',calling_func, '.\n', &
         'At lat,lon,k,macmic_it, timestep=',get_rlat_p(chunk,icol),get_rlon_p(chunk,icol),klev,macmiciter,get_nstep()

  end subroutine report_error_info


end module debug_info
