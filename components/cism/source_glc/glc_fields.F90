!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module glc_fields

!BOP
! !MODULE: glc_fields

! !DESCRIPTION:
!  Holds coupling fields and other important data
!
! !REVISION HISTORY:
! 
!  Author: William Lipscomb, LANL
!
! !USES:

  use glc_kinds_mod
  use glad_main, only: glad_params

  implicit none
  save

! !PUBLIC MEMBER FUNCTIONS:

!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

  ! Fields received from CESM coupler
  
  real(r8),dimension(:,:), allocatable ::  & 
     tsfc        ,&! surface temperature (Celsius)
                   ! received from coupler in Kelvin, must be converted
     qsmb          ! flux of new glacier ice (kg/m^2/s)

  ! output to coupler

  ! the following need to be targets for the sake of the override code in glc_export
  real(r8),dimension(:,:), allocatable, target ::  &
     ice_covered,&         ! whether each grid cell is ice-covered [0,1]
     topo                  ! glacier surface elevation (m)

  real(r8),dimension(:,:), allocatable :: &
     rofi       ,&! ice runoff (calving) flux (kg/m^2/s)
     rofl       ,&! ice runoff (calving) flux (kg/m^2/s)
     hflx         ! heat flux from glacier interior, positive down (W/m^2)

  ! Note that there are two separate mask fields. Both of them provide information about
  ! where CISM is running. The difference is that ice_sheet_grid_mask includes icesheet
  ! areas that are diagtostic-only, whereas icemask_coupled_fluxes excludes icesheet
  ! areas where we are zeroing the fluxes sent to the coupler (thus, icesheets that are
  ! "diagnostic" in some sense). We need two separate maps, as opposed to a single map
  ! plus a scalar logical variable, in case we're running with multiple icesheet
  ! instances (e.g., Greenland & Antarctica), one of which is fully prognostic and one of
  ! which is diagnostic-only: in that case, ice_sheet_grid_mask would be non-zero over
  ! both Greenland and Antarctica, whereas icemask_coupled_fluxes would be non-zero over
  ! (e.g.) Greenland, but 0 over Antarctica.
  real(r8),dimension(:,:), allocatable :: &
     ice_sheet_grid_mask, &  ! mask of ice sheet grid coverage
     icemask_coupled_fluxes  ! mask of ice sheet grid coverage where we are potentially sending non-zero fluxes

  type(glad_params) :: ice_sheet   ! Parameters relevant to all model instances

!EOP
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_allocate_fields
! !INTERFACE:

 subroutine glc_allocate_fields (nx, ny)

! !DESCRIPTION:
!  Allocate fields declared here
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !USES:
   use glc_kinds_mod

! !INPUT/OUTPUT PARAMETERS:

   integer (i4), intent(in) :: &
        nx, ny           ! grid dimensions

!EOP
!BOC

   ! from coupler
   allocate(tsfc(nx,ny))
   allocate(qsmb(nx,ny))

   ! to coupler
   allocate(ice_covered(nx,ny))
   allocate(topo(nx,ny))
   allocate(rofi(nx,ny))
   allocate(rofl(nx,ny))
   allocate(hflx(nx,ny))
   allocate(ice_sheet_grid_mask(nx,ny))
   allocate(icemask_coupled_fluxes(nx,ny))
   
 end subroutine glc_allocate_fields

!***********************************************************************

!BOP
! !IROUTINE: glc_deallocate_fields
! !INTERFACE:

 subroutine glc_deallocate_fields

! !DESCRIPTION:
!  Deallocate global arrays on glc grid.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT/OUTPUT PARAMETERS:


!EOP
!BOC

   ! from coupler
   deallocate(tsfc)
   deallocate(qsmb)

   ! to coupler
   deallocate(ice_covered)
   deallocate(topo)
   deallocate(rofi)
   deallocate(rofl)
   deallocate(hflx)
   deallocate(ice_sheet_grid_mask)
   deallocate(icemask_coupled_fluxes)
  
   end subroutine glc_deallocate_fields

!***********************************************************************

 end module glc_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
