module topounit_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing topounit indices and associated variables and routines.
  !
  ! !USES:
  use ncdio_pio       , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
  use ncdio_pio       , only : ncd_io, check_var, ncd_inqfdims, check_dim, ncd_inqdid, ncd_inqdlen
  use clm_varctl      , only: fsurdat, iulog
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use pio
  use spmdMod
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !------------------------------------------------------------------
  ! Initialize topounit type constants
  !------------------------------------------------------------------
  
  !integer, parameter, public :: max_topounits  = 1 ! maximum number of topounits per gridcell
  integer, public :: max_topounits               ! maximum number of topounits per gridcell
  logical, public :: has_topounit                ! true => topounit dimension is on dataset
  
    !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: topounit_varcon_init  ! initialize constants in this module
  !-----------------------------------------------------------------------
  
  contains
  
  !-----------------------------------------------------------------------
  subroutine topounit_varcon_init(lfsurdat)
    !
    ! !DESCRIPTION:
    ! Initialize topounit parameters
    !
    ! !USES:
    use fileutils   , only : getfil    
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t)     :: ncid         ! netcdf id
    integer :: dimid
    integer :: topounits_size                      ! Size of topounit dimension in the surface data 
    character(len=256):: locfn                ! local file name
    character(len=*), parameter :: subname = 'topounit_varcon_init'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read topounit information from surface data .....'
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'
          !call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    endif
    
    ! Read surface data
    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    
     !Obtain the maximum number of topounits
    call ncd_inqdid(ncid, 'topounit', dimid, dimexist=has_topounit)
    if (.not. has_topounit) then
         max_topounits = 1
         if (masterproc) then
             write(iulog,*)'Surface dataset has no topounit dimention; max_topounits is set to 1'
         end if
    else
        call ncd_inqdlen(ncid, dimid, topounits_size)  ! Get the dimension size of topounit from file
        max_topounits = topounits_size
    end if

    call ncd_pio_closefile(ncid)

  end subroutine topounit_varcon_init
  
end module topounit_varcon
