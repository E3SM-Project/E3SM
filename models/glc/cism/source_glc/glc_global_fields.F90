!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_global_fields

!BOP
! !MODULE: glc_global_fields

! !DESCRIPTION:
!  Holds fields and parameters related to the global grid
!
! !REVISION HISTORY:
! 
!  Author: William Lipscomb, LANL
!
! !USES:

   use glc_kinds_mod
   use glint_example_clim, only: glex_climate
   use glint_main, only: glint_params

   implicit none
   save

! !PUBLIC MEMBER FUNCTIONS:

!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

  ! grid dimensions
  ! NOTE: The glc grid is indexed from W to E and S to N.
  !       The glint grid is indexed from W to E and N to S.
  ! All fields declared below use the glint convention, so the latitude
  !  index must be reversed when exchanging with the coupler.

  ! Fields received from CESM coupler on glc grid
  ! NOTE: For SEB scheme, tsfc = ground surface temp
  !                       qsmb = accumulation minus ablation
  !       For PDD scheme, tsfc => 2m reference temp
  !                       qsmb => precipitation

  ! input from coupler (3rd dimension for elevation classes)

  real(r8),dimension(:,:,:), allocatable ::  & 
     tsfc        ,&! surface temperature (Celsius)
                   ! received from coupler in Kelvin, must be converted
     topo        ,&! surface elevation (m)
     qsmb          ! flux of new glacier ice (kg/m^2/s)

  ! output to coupler

  real(r8),dimension(:,:,:), allocatable ::  &   ! per elevation class
     gfrac       ,&! fractional glacier area [0,1] 
     gtopo       ,&! glacier surface elevation (m)
     grofi       ,&! ice runoff (calving) flux (kg/m^2/s)
     grofl       ,&! ice runoff (calving) flux (kg/m^2/s)
     ghflx         ! heat flux from glacier interior, positive down (W/m^2)

  type(glint_params) :: ice_sheet   ! Parameters relevant to all model instances, 
                                    ! i.e. pertaining to the global model 

  type(glex_climate) :: climate     ! Climate parameters and fields

  ! glint input fields
  real(r8),dimension(:,:), allocatable ::  & 
      temp        ,&! temperature     (deg C)
      precip      ,&! precipitation   (mm/s) 
      orog          ! orography       (m) 

!lipscomb - TO DO - Not all of the remaining fields are needed
!                   Some are commented out for now
 
  ! glint output fields (using glint N to S indexing convention)
  ! These are all on the global grid, except for orog_out
 
  real(r8),dimension(:,:), allocatable ::   &
!!      albedo      ,&! Fractional albedo
!!      orog_out    ,&! Output orography (m)
      ice_frac    ,&! Ice coverage fraction
      fw          ,&! Freshwater output flux (mm/s)
      fw_in         ! Freshwater input flux (mm/s)
 
  ! arrays which hold information about the global grid
 
!!  real(r8),dimension(:), allocatable ::  &
!!     lats_orog    ,&! Latitudes of global orography gridpoints
!!     lons_orog      ! Longitudes of global oropraphy gridpoints

  ! arrays which hold information about the ice model instances
 
!!  real(r8),dimension(:,:), allocatable ::   &
!!      coverage    ,&! Coverage map for normal global grid
!!      cov_orog      ! Coverage map for orography grid

!!  integer(i4) :: time           ! current time in hours

!EOP
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_allocate_global
! !INTERFACE:

 subroutine glc_allocate_global (nx, ny, glc_nec)

! !DESCRIPTION:
!  Allocate global arrays on glc grid.
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
      nx, ny           ! global grid dimensions

   integer (i4), intent(in), optional ::  &
      glc_nec              ! number of elevation classes
!EOP
!BOC

!!   integer (i4) :: nxo, nyo  ! not currently used

 ! from coupler
   allocate(tsfc(nx,ny,glc_nec))
   allocate(topo(nx,ny,glc_nec))
   allocate(qsmb(nx,ny,glc_nec))

 ! to coupler
   allocate(gfrac(nx,ny,glc_nec))
   allocate(gtopo(nx,ny,glc_nec))
   allocate(grofi(nx,ny,glc_nec))
   allocate(grofl(nx,ny,glc_nec))
   allocate(ghflx(nx,ny,glc_nec))

 ! Other fields

   allocate(temp(nx,ny),precip(nx,ny),orog(nx,ny))
   allocate(ice_frac(nx,ny),fw(nx,ny),fw_in(nx,ny))

!lipscomb - TO DO - I think these are not needed
!!   nxo = nx; nyo = ny
!!   allocate(orog_out(nxo,nyo),albedo(nx,ny))
!!   allocate(lats_orog(nyo),lons_orog(nxo))
!!   allocate(coverage(nx,ny),cov_orog(nxo,nyo))

   end subroutine glc_allocate_global

!***********************************************************************

!BOP
! !IROUTINE: glc_deallocate_global
! !INTERFACE:

 subroutine glc_deallocate_global

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
   deallocate(topo)
   deallocate(qsmb)

 ! to coupler
   deallocate(gfrac)
   deallocate(gtopo)
   deallocate(grofi)
   deallocate(grofl)
   deallocate(ghflx)

 ! Other fields
!lipscomb - TO DO - Some of these are not needed
  deallocate(temp)
  deallocate(orog)
  deallocate(precip)
!!  deallocate(orog_out)
!!  deallocate(albedo)
  deallocate(ice_frac)  
  deallocate(fw)  
  deallocate(fw_in)
!!  deallocate(lats_orog)
!!  deallocate(lons_orog)
!!  deallocate(coverage)
!!  deallocate(cov_orog)
  
   end subroutine glc_deallocate_global

!***********************************************************************

 end module glc_global_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
