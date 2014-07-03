module readinitial

implicit none

public :: read_initial

contains

subroutine read_initial(ncid)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Ensure that requisite netcdf variables are on the initial dataset.
    !          Set base day and date info using the "current" values from it.
    ! 
    ! Method: Issue proper netcdf wrapper calls.  Broadcast to slaves if SPMD
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plon
    use abortutils,     only: endrun
    use scamMod,        only: single_column   
    use cam_logfile,    only: iulog
    use pio, only : pio_inq_dimlen, pio_inq_varid,  pio_inq_dimid, &
         var_desc_t, pio_internal_error, pio_noerr, pio_bcast_error, &
         pio_seterrorhandling, file_desc_t

    ! Arguments
    type(file_desc_t), intent(inout) :: ncid  ! initial file

    ! Local variables
    integer :: lonid
    integer :: latid

    integer :: mlon             ! longitude dimension length from dataset
    integer :: morec            ! latitude dimension length from dataset

    integer :: ncollen
    integer :: ierr
    logical :: isncol = .false.

    character(len=*), parameter :: routine = 'read_initial'
    !-----------------------------------------------------------------------

    ! ****N.B.*** All this code belongs in the dyn_init methods.

    !  Tells pio to communicate any error to all tasks so that it can be handled locally
    !  otherwise and by default pio error are handled internally.

    call pio_seterrorhandling( ncid, PIO_BCAST_ERROR)
    ierr = PIO_inq_dimid( ncid, 'ncol', lonid)
    call pio_seterrorhandling( ncid, PIO_INTERNAL_ERROR)

    if(ierr==PIO_NOERR) then
       isncol=.true.
      ierr = PIO_inq_dimlen (ncid, lonid , ncollen)
    else

       !
       ! Get and check dimension/date info
       !
       ierr = PIO_inq_dimid (ncid, 'lon' , lonid)
       ierr = PIO_inq_dimid (ncid, 'lat' , latid)

       ierr = PIO_inq_dimlen (ncid, lonid , mlon)
       ierr = PIO_inq_dimlen (ncid, latid , morec)
    endif

    if ((.not. isncol) .and. (.not. single_column) .and. (mlon /= plon.or.morec /= plat)) then
       write(iulog,*) routine//': model parameters do not match initial dataset parameters'
       write(iulog,*)'Model Parameters:    plon = ',plon,' plat = ',plat
       write(iulog,*)'Dataset Parameters:  dlon = ',mlon,' dlat = ',morec
       call endrun(routine//': model parameters do not match initial dataset parameters')
    end if

  end subroutine read_initial

end module readinitial
