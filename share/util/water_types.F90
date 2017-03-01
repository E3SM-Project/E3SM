
module water_types

!-----------------------------------------------------------------------
!
! Provide core functionality for types of condensed water to be used
! with the water vapor tracers.
!
! This module works in with "water_isotopes"  and "water_tracers".
!
! All interface routine are identified by wtype_*, etc.
!
! 5 types of water are available, three phases (vapor, cloud liquid
! and cloud ice) and precipitation (rain and snow).
!
! Author: Chuck Bardeen (2/4/2012)
!-----------------------------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8

  implicit none
  
  private
  save

!------------------------ Module Interfaces -----------------------------
!
! Public interfaces
!
  public :: wtype_init            ! initilize water types
  public :: wtype_get_itype       ! lookup a species index by name

  public :: wtype_get_alpha       ! isotope fractionation


!------------------- Module Variable Declarations -----------------------
!
! Water tracer type identifiers
  integer, parameter, public :: pwtype     = 7     ! number of water types

  integer, parameter, public :: iwtundef   = 0     ! not water type 
  integer, parameter, public :: iwtvap     = 1     ! water type is vapour
  integer, parameter, public :: iwtliq     = 2     ! water type is liquid
  integer, parameter, public :: iwtice     = 3     ! water type is ice
  integer, parameter, public :: iwtstrain  = 4     ! water type is stratiform rain
  integer, parameter, public :: iwtstsnow  = 5     ! water type is stratiform snow
  integer, parameter, public :: iwtcvrain  = 6     ! water type is convective rain
  integer, parameter, public :: iwtcvsnow  = 7     ! water type is convective snow

! Water type names
  character(len=8), dimension(pwtype), parameter, public :: & ! water type names 
      wtype_names      = (/ 'VAPOR   ', 'LIQUID  ', 'ICE     ', 'RAINS   ', 'SNOWS   ', 'RAINC   ', 'SNOWC   ' /)

! Water type Suffix
character(len=2), dimension(pwtype), parameter, public :: & ! suffix names 
      wtype_suffix     =   (/ '_v', '_l', '_i', '_R', '_S', '_r', '_s' /)


!
!-----------------------------------------------------------------------
contains

!=======================================================================
  subroutine wtype_init
!-----------------------------------------------------------------------
! Purpose: Initialize module internal data arrays
!-----------------------------------------------------------------------
    write(6,*) 'WTYPE_INIT: Initializing water types.'
    return
  end subroutine wtype_init


!=======================================================================
  function wtype_get_itype(name)
!-----------------------------------------------------------------------
! Purpose: Retrieve type index, based on type name
! Author: Chuck Bardeen
!-----------------------------------------------------------------------
    character(len=*),  intent(in)  :: name              ! water type name
    integer                        :: wtype_get_itype   ! return species index
!-----------------------------------------------------------------------
    do wtype_get_itype = 1, pwtype
      if (name == wtype_names(wtype_get_itype)) then
        return
      end if
    end do
    
    wtype_get_itype = iwtundef
    
    return
  end function wtype_get_itype
  
!=========================================================================


!=======================================================================
  function wtype_get_alpha(ispec, isrctype, idsttype, tk, rh, do_kinetic)
!-----------------------------------------------------------------------
! Purpose: Retrieve the fractionation for a process that goes from
!          the source water type to the destination water type.
!
! Author: Chuck Bardeen
!-----------------------------------------------------------------------
    use water_isotopes, only : wiso_alpl, wiso_alpi, wiso_akel, wiso_akci

    integer,  intent(in)            :: ispec           ! isotope species index
    integer,  intent(in)            :: isrctype        ! source water type index
    integer,  intent(in)            :: idsttype        ! destination water type index
    real(r8), intent(in)            :: tk              ! temperature (K)
    real(r8),  intent(in)           :: rh              ! relative humidity (fraction)
    logical, intent(in)             :: do_kinetic      ! use kinetic calculation
    real(r8)                        :: wtype_get_alpha ! return alpha
    
!-----------------------------------------------------------------------

    ! If their types are the same, then no fractionation occurs.
    wtype_get_alpha = 1._r8
    
    if (isrctype /= idsttype) then
      
      ! Is the source vapor?
      if (isrctype == iwtvap) then
       
        ! Is the destination a liquid?
        if ((idsttype == iwtliq) .or. (idsttype == iwtstrain) .or. (idsttype == iwtcvrain)) then
          wtype_get_alpha = wiso_alpl(ispec,tk)
          
          if (do_kinetic) then
            wtype_get_alpha = wiso_akel(ispec,tk,rh,wtype_get_alpha)
          end if
        else
          wtype_get_alpha = wiso_alpi(ispec,tk)
          
          if (do_kinetic) then
            wtype_get_alpha = wiso_akci(ispec,tk,wtype_get_alpha)
          end if
        end if
        
      ! Is the destination vapor?  
      else if (idsttype == iwtvap) then

        ! Is the source a liquid?
        if ((isrctype == iwtliq) .or. (isrctype == iwtstrain) .or. (isrctype == iwtcvrain)) then
          wtype_get_alpha = wiso_alpl(ispec,tk)
          
          if (do_kinetic) then
            wtype_get_alpha = wiso_akel(ispec,tk,rh,wtype_get_alpha)
          end if
          wtype_get_alpha = 1._r8 / wtype_get_alpha
        else
          wtype_get_alpha = 1._r8 !No fractionation occurs during sublimation
        end if
      end if
    end if
    
    return
  end function wtype_get_alpha
  
!=========================================================================

end module water_types
