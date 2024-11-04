module camsrfexch
!-----------------------------------------------------------------------
!
! Module to handle data that is exchanged between the CAM atmosphere
! model and the surface models (land, sea-ice, and ocean).
!
!-----------------------------------------------------------------------
!
! USES:
!
  use shr_kind_mod,  only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents,  only: pcnst
  use ppgrid,        only: pcols, begchunk, endchunk
  use phys_grid,     only: get_ncols_p, phys_grid_initialized
  use infnan,        only: posinf, assignment(=)
  use cam_abortutils,    only: endrun
  use cam_logfile,   only: iulog

  implicit none

!----------------------------------------------------------------------- 
! PRIVATE: Make default data and interfaces private
!----------------------------------------------------------------------- 
  private     ! By default all data is private to this module
!
! Public interfaces
!
  public atm2hub_alloc              ! Atmosphere to surface data allocation method
  public hub2atm_alloc              ! Merged hub surface to atmosphere data allocation method
  public atm2hub_deallocate
  public hub2atm_deallocate
  public cam_export
!
! Public data types
!
  public cam_out_t                  ! Data from atmosphere
  public cam_in_t                   ! Merged surface data

!---------------------------------------------------------------------------
! This is the data that is sent from the atmosphere to the surface models
!---------------------------------------------------------------------------

  type cam_out_t 
     integer  :: lchnk               ! chunk index
     integer  :: ncol                ! number of columns in chunk
     real(r8), allocatable :: tbot(:)     ! bot level temperature
     real(r8), allocatable :: zbot(:)     ! bot level height above surface
     real(r8), allocatable :: ubot(:)     ! bot level u wind
     real(r8), allocatable :: vbot(:)     ! bot level v wind
     real(r8), allocatable :: qbot(:,:)   ! bot level specific humidity
     real(r8), allocatable :: pbot(:)     ! bot level pressure
     real(r8), allocatable :: rho(:)      ! bot level density	
     real(r8), allocatable :: netsw(:)    !	
     real(r8), allocatable :: flwds(:)    ! 
     real(r8), allocatable :: precsc(:)   !
     real(r8), allocatable :: precsl(:)   !
     real(r8), allocatable :: precc(:)    ! 
     real(r8), allocatable :: precl(:)    ! 
     real(r8), allocatable :: soll(:)     ! 
     real(r8), allocatable :: sols(:)     ! 
     real(r8), allocatable :: solld(:)    !
     real(r8), allocatable :: solsd(:)    !
     real(r8), allocatable :: thbot(:)    ! 
     real(r8), allocatable :: co2prog(:)  ! prognostic co2
     real(r8), allocatable :: co2diag(:)  ! diagnostic co2
     real(r8), allocatable :: psl(:)
     real(r8), allocatable :: bcphiwet(:) ! wet deposition of hydrophilic black carbon
     real(r8), allocatable :: bcphidry(:) ! dry deposition of hydrophilic black carbon
     real(r8), allocatable :: bcphodry(:) ! dry deposition of hydrophobic black carbon
     real(r8), allocatable :: ocphiwet(:) ! wet deposition of hydrophilic organic carbon
     real(r8), allocatable :: ocphidry(:) ! dry deposition of hydrophilic organic carbon
     real(r8), allocatable :: ocphodry(:) ! dry deposition of hydrophobic organic carbon
     real(r8), allocatable :: dstwet1(:)  ! wet deposition of dust (bin1)
     real(r8), allocatable :: dstdry1(:)  ! dry deposition of dust (bin1)
     real(r8), allocatable :: dstwet2(:)  ! wet deposition of dust (bin2)
     real(r8), allocatable :: dstdry2(:)  ! dry deposition of dust (bin2)
     real(r8), allocatable :: dstwet3(:)  ! wet deposition of dust (bin3)
     real(r8), allocatable :: dstdry3(:)  ! dry deposition of dust (bin3)
     real(r8), allocatable :: dstwet4(:)  ! wet deposition of dust (bin4)
     real(r8), allocatable :: dstdry4(:)  ! dry deposition of dust (bin4)
     real(r8), allocatable :: wsresp(:)   ! first-order response of low-level wind to surface fluxes
     real(r8), allocatable :: tau_est(:)  ! stress estimated to be in equilibrium with ubot/vbot
     real(r8), allocatable :: ugust(:)    ! gustiness value
     real(r8), allocatable :: uovern(:)       ! ratio of wind speed/brunt vaisalla frequency  
  end type cam_out_t 

!---------------------------------------------------------------------------
! This is the merged state of sea-ice, land and ocean surface parameterizations
!---------------------------------------------------------------------------

  type cam_in_t    
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns
     real(r8), allocatable :: asdir(:)      ! albedo: shortwave, direct
     real(r8), allocatable :: asdif(:)      ! albedo: shortwave, diffuse
     real(r8), allocatable :: aldir(:)      ! albedo: longwave, direct
     real(r8), allocatable :: aldif(:)      ! albedo: longwave, diffuse
     real(r8), allocatable :: lwup(:)       ! longwave up radiative flux
     real(r8), allocatable :: lhf(:)        ! latent heat flux
     real(r8), allocatable :: shf(:)        ! sensible heat flux
     real(r8), allocatable :: h2otemp(:)    ! water temperature heat flux from ocean
     real(r8), allocatable :: wsx(:)        ! surface u-stress (N)
     real(r8), allocatable :: wsy(:)        ! surface v-stress (N)
     real(r8), allocatable :: tref(:)       ! ref height surface air temp
     real(r8), allocatable :: qref(:)       ! ref height specific humidity 
     real(r8), allocatable :: u10(:)        ! 10m wind speed
     real(r8), allocatable :: u10withgusts(:) ! 10m wind speed with gustiness
     real(r8), allocatable :: ts(:)         ! merged surface temp 
     real(r8), allocatable :: sst(:)        ! sea surface temp
     real(r8), allocatable :: snowhland(:)  ! snow depth (liquid water equivalent) over land
     real(r8), allocatable :: snowhice(:)   ! snow depth over ice
     real(r8), allocatable :: fco2_lnd(:)   ! co2 flux from lnd
     real(r8), allocatable :: fco2_ocn(:)   ! co2 flux from ocn
     real(r8), allocatable :: fdms(:)       ! dms flux
     real(r8), allocatable :: landfrac(:)   ! land area fraction
     real(r8), allocatable :: icefrac(:)    ! sea-ice areal fraction
     real(r8), allocatable :: ocnfrac(:)    ! ocean areal fraction
     real(r8), pointer, dimension(:) :: ram1       !aerodynamical resistance (s/m) (pcols)
     real(r8), pointer, dimension(:) :: fv         !friction velocity (m/s) (pcols)
     real(r8), pointer, dimension(:) :: soilw      !volumetric soil water (m3/m3)
     real(r8), allocatable :: cflx(:,:)     ! constituent flux (emissions)
     real(r8), allocatable :: ustar(:)      ! atm/ocn saved version of ustar
     real(r8), allocatable :: re(:)         ! atm/ocn saved version of re
     real(r8), allocatable :: ssq(:)        ! atm/ocn saved version of ssq
     real(r8), pointer, dimension(:,:) :: depvel   ! deposition velocities
     real(r8), pointer, dimension(:,:) :: dstflx   ! dust fluxes
     real(r8), pointer, dimension(:,:) :: meganflx ! MEGAN fluxes
  end type cam_in_t    

!===============================================================================
CONTAINS
!===============================================================================

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: hub2atm_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the surface to atmosphere data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine hub2atm_alloc( cam_in )
    use seq_drydep_mod,  only: lnd_drydep, n_drydep
    use cam_cpl_indices, only: index_x2a_Sl_ram1, index_x2a_Sl_fv, index_x2a_Sl_soilw, index_x2a_Fall_flxdst1
    use cam_cpl_indices, only: index_x2a_Fall_flxvoc
    use shr_megan_mod,   only: shr_megan_mechcomps_n

!
!!ARGUMENTS:
!
   type(cam_in_t), pointer ::  cam_in(:)     ! Merged surface state
!
!!LOCAL VARIABLES:
!
    integer :: c        ! chunk index
    integer :: ierror   ! Error code
!----------------------------------------------------------------------- 
! 
! EOP
!
    if ( .not. phys_grid_initialized() ) call endrun( "HUB2ATM_ALLOC error: phys_grid not called yet" )
    allocate (cam_in(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('HUB2ATM_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk
       nullify(cam_in(c)%ram1)
       nullify(cam_in(c)%fv)
       nullify(cam_in(c)%soilw)
       nullify(cam_in(c)%depvel)
       nullify(cam_in(c)%dstflx)
       nullify(cam_in(c)%meganflx)
    enddo  
    do c = begchunk,endchunk 
       allocate (cam_in(c)%asdir(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error asdir')

       allocate (cam_in(c)%asdif(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error asdif')

       allocate (cam_in(c)%aldir(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error aldir')

       allocate (cam_in(c)%aldif(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error aldif')

       allocate (cam_in(c)%lwup(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error lwup')

       allocate (cam_in(c)%lhf(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error lhf')

       allocate (cam_in(c)%shf(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error shf')

       allocate (cam_in(c)%h2otemp(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error h2otemp')

       allocate (cam_in(c)%wsx(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error wsx')

       allocate (cam_in(c)%wsy(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error wsy')

       allocate (cam_in(c)%tref(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error tref')

       allocate (cam_in(c)%qref(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error qref')

       allocate (cam_in(c)%u10(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error u10')

       allocate (cam_in(c)%u10withgusts(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error u10withgusts')

       allocate (cam_in(c)%ts(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ts')

       allocate (cam_in(c)%sst(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error sst')

       allocate (cam_in(c)%snowhland(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error snowhland')

       allocate (cam_in(c)%snowhice(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error snowhice')

       allocate (cam_in(c)%fco2_lnd(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fco2_lnd')

       allocate (cam_in(c)%fco2_ocn(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fco2_ocn')

       allocate (cam_in(c)%fdms(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fdms')

       allocate (cam_in(c)%landfrac(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error landfrac')

       allocate (cam_in(c)%icefrac(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error icefrac')

       allocate (cam_in(c)%ocnfrac(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ocnfrac')

       if (index_x2a_Sl_ram1>0) then
          allocate (cam_in(c)%ram1(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ram1')
       endif
       if (index_x2a_Sl_fv>0) then
          allocate (cam_in(c)%fv(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fv')
       endif
       if (index_x2a_Sl_soilw /= 0) then
          allocate (cam_in(c)%soilw(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error soilw')
       end if

       allocate (cam_in(c)%cflx(pcols,pcnst), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error cflx')

       allocate (cam_in(c)%ustar(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ustar')

       allocate (cam_in(c)%re(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error re')

       allocate (cam_in(c)%ssq(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ssq')

       if (index_x2a_Fall_flxdst1>0) then
          ! Assume 4 bins from surface model ....
          allocate (cam_in(c)%dstflx(pcols,4), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error dstflx')
       endif
       if ( index_x2a_Fall_flxvoc>0 .and. shr_megan_mechcomps_n>0 ) then
          allocate (cam_in(c)%meganflx(pcols,shr_megan_mechcomps_n), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error meganflx')
       endif
    end do

    if (lnd_drydep .and. n_drydep>0) then
       do c = begchunk,endchunk 
          allocate (cam_in(c)%depvel(pcols,n_drydep), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error depvel')
       end do
    endif

    do c = begchunk,endchunk
       cam_in(c)%lchnk = c
       cam_in(c)%ncol  = get_ncols_p(c)
       cam_in(c)%asdir    (:) = 0._r8
       cam_in(c)%asdif    (:) = 0._r8
       cam_in(c)%aldir    (:) = 0._r8
       cam_in(c)%aldif    (:) = 0._r8
       cam_in(c)%lwup     (:) = 0._r8
       cam_in(c)%lhf      (:) = 0._r8
       cam_in(c)%shf      (:) = 0._r8
       cam_in(c)%h2otemp  (:) = 0._r8
       cam_in(c)%wsx      (:) = 0._r8
       cam_in(c)%wsy      (:) = 0._r8
       cam_in(c)%tref     (:) = 0._r8
       cam_in(c)%qref     (:) = 0._r8
       cam_in(c)%u10      (:) = 0._r8
       cam_in(c)%u10withgusts(:) = 0._r8
       cam_in(c)%ts       (:) = 0._r8
       cam_in(c)%sst      (:) = 0._r8
       cam_in(c)%snowhland(:) = 0._r8
       cam_in(c)%snowhice (:) = 0._r8
       cam_in(c)%fco2_lnd (:) = 0._r8
       cam_in(c)%fco2_ocn (:) = 0._r8
       cam_in(c)%fdms     (:) = 0._r8
       cam_in(c)%landfrac (:) = posinf
       cam_in(c)%icefrac  (:) = posinf
       cam_in(c)%ocnfrac  (:) = posinf

       if (associated(cam_in(c)%ram1)) &
            cam_in(c)%ram1  (:) = 0.1_r8
       if (associated(cam_in(c)%fv)) &
            cam_in(c)%fv    (:) = 0.1_r8
       if (associated(cam_in(c)%soilw)) &
            cam_in(c)%soilw (:) = 0.0_r8
       if (associated(cam_in(c)%dstflx)) &
            cam_in(c)%dstflx(:,:) = 0.0_r8
       if (associated(cam_in(c)%meganflx)) &
            cam_in(c)%meganflx(:,:) = 0.0_r8

       cam_in(c)%cflx   (:,:) = 0._r8
       cam_in(c)%ustar    (:) = 0._r8
       cam_in(c)%re       (:) = 0._r8
       cam_in(c)%ssq      (:) = 0._r8
       if (lnd_drydep .and. n_drydep>0) then
          cam_in(c)%depvel (:,:) = 0._r8
       endif
    end do

  end subroutine hub2atm_alloc

!
!===============================================================================
!

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: atm2hub_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the atmosphere to surface data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine atm2hub_alloc( cam_out )
!
!!USES:
!
!
!!ARGUMENTS:
!
   type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
!
!!LOCAL VARIABLES:
!
    integer :: c            ! chunk index
    integer :: ierror       ! Error code
    !----------------------------------------------------------------------- 

    if ( .not. phys_grid_initialized() ) call endrun( "ATM2HUB_ALLOC error: phys_grid not called yet" )
    allocate (cam_out(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('ATM2HUB_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk 
       allocate (cam_out(c)%tbot(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error tbot')

       allocate (cam_out(c)%zbot(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error zbot')

       allocate (cam_out(c)%ubot(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error ubot')

       allocate (cam_out(c)%vbot(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error vbot')

       allocate (cam_out(c)%qbot(pcols,pcnst), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error qbot')

       allocate (cam_out(c)%pbot(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error pbot')

       allocate (cam_out(c)%rho(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error rho')

       allocate (cam_out(c)%netsw(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error netsw')

       allocate (cam_out(c)%flwds(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error flwds')

       allocate (cam_out(c)%precsc(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precsc')

       allocate (cam_out(c)%precsl(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precsl')

       allocate (cam_out(c)%precc(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precc')

       allocate (cam_out(c)%precl(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precl')

       allocate (cam_out(c)%soll(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error soll')

       allocate (cam_out(c)%sols(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error sols')

       allocate (cam_out(c)%solld(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error solld')

       allocate (cam_out(c)%solsd(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error solsd')

       allocate (cam_out(c)%thbot(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error thbot')

       allocate (cam_out(c)%co2prog(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error co2prog')

       allocate (cam_out(c)%co2diag(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error co2diag')

       allocate (cam_out(c)%psl(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB ALLOC error: allocation error psl')

       allocate (cam_out(c)%bcphiwet(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error bcphiwet')

       allocate (cam_out(c)%bcphidry(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error bcphidry')

       allocate (cam_out(c)%bcphodry(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error bcphodry')

       allocate (cam_out(c)%ocphiwet(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error ocphiwet')

       allocate (cam_out(c)%ocphidry(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error ocphidry')

       allocate (cam_out(c)%ocphodry(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error ocphodry')

       allocate (cam_out(c)%dstwet1(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstwet1')

       allocate (cam_out(c)%dstdry1(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstdry1')

       allocate (cam_out(c)%dstwet2(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstwet2')

       allocate (cam_out(c)%dstdry2(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstdry2')

       allocate (cam_out(c)%dstwet3(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstwet3')

       allocate (cam_out(c)%dstdry3(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstdry3')

       allocate (cam_out(c)%dstwet4(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstwet4')

       allocate (cam_out(c)%dstdry4(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error dstdry4')

       allocate (cam_out(c)%wsresp(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error wsresp')

       allocate (cam_out(c)%tau_est(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error tau_est')

       allocate (cam_out(c)%ugust(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error ugust')
       
       allocate (cam_out(c)%uovern(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error uovern')
    enddo  

    do c = begchunk,endchunk
       cam_out(c)%lchnk       = c
       cam_out(c)%ncol        = get_ncols_p(c)
       cam_out(c)%tbot(:)     = 0._r8
       cam_out(c)%zbot(:)     = 0._r8
       cam_out(c)%ubot(:)     = 0._r8
       cam_out(c)%vbot(:)     = 0._r8
       cam_out(c)%qbot(:,:)   = 0._r8
       cam_out(c)%pbot(:)     = 0._r8
       cam_out(c)%rho(:)      = 0._r8
       cam_out(c)%netsw(:)    = 0._r8
       cam_out(c)%flwds(:)    = 0._r8
       cam_out(c)%precsc(:)   = 0._r8
       cam_out(c)%precsl(:)   = 0._r8
       cam_out(c)%precc(:)    = 0._r8
       cam_out(c)%precl(:)    = 0._r8
       cam_out(c)%soll(:)     = 0._r8
       cam_out(c)%sols(:)     = 0._r8
       cam_out(c)%solld(:)    = 0._r8
       cam_out(c)%solsd(:)    = 0._r8
       cam_out(c)%thbot(:)    = 0._r8
       cam_out(c)%co2prog(:)  = 0._r8
       cam_out(c)%co2diag(:)  = 0._r8
       cam_out(c)%psl(:)      = 0._r8
       cam_out(c)%bcphidry(:) = 0._r8
       cam_out(c)%bcphodry(:) = 0._r8
       cam_out(c)%bcphiwet(:) = 0._r8
       cam_out(c)%ocphidry(:) = 0._r8
       cam_out(c)%ocphodry(:) = 0._r8
       cam_out(c)%ocphiwet(:) = 0._r8
       cam_out(c)%dstdry1(:)  = 0._r8
       cam_out(c)%dstwet1(:)  = 0._r8
       cam_out(c)%dstdry2(:)  = 0._r8
       cam_out(c)%dstwet2(:)  = 0._r8
       cam_out(c)%dstdry3(:)  = 0._r8
       cam_out(c)%dstwet3(:)  = 0._r8
       cam_out(c)%dstdry4(:)  = 0._r8
       cam_out(c)%dstwet4(:)  = 0._r8
       cam_out(c)%wsresp(:)   = 0._r8
       cam_out(c)%tau_est(:)  = 0._r8
       cam_out(c)%ugust(:)    = 0._r8
       cam_out(c)%uovern(:)   = 0._r8
    end do

  end subroutine atm2hub_alloc

  subroutine atm2hub_deallocate(cam_out)
    type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
    integer :: c

    if(associated(cam_out)) then
       do c = begchunk,endchunk
          deallocate(cam_out(c)%tbot)
          deallocate(cam_out(c)%zbot)
          deallocate(cam_out(c)%ubot)
          deallocate(cam_out(c)%vbot)
          deallocate(cam_out(c)%qbot)
          deallocate(cam_out(c)%pbot)
          deallocate(cam_out(c)%rho)
          deallocate(cam_out(c)%netsw)
          deallocate(cam_out(c)%flwds)
          deallocate(cam_out(c)%precsc)
          deallocate(cam_out(c)%precsl)
          deallocate(cam_out(c)%precc)
          deallocate(cam_out(c)%precl)
          deallocate(cam_out(c)%soll)
          deallocate(cam_out(c)%sols)
          deallocate(cam_out(c)%solld)
          deallocate(cam_out(c)%solsd)
          deallocate(cam_out(c)%thbot)
          deallocate(cam_out(c)%co2prog)
          deallocate(cam_out(c)%co2diag)
          deallocate(cam_out(c)%bcphiwet)
          deallocate(cam_out(c)%bcphidry)
          deallocate(cam_out(c)%bcphodry)
          deallocate(cam_out(c)%ocphiwet)
          deallocate(cam_out(c)%ocphidry)
          deallocate(cam_out(c)%ocphodry)
          deallocate(cam_out(c)%dstwet1)
          deallocate(cam_out(c)%dstdry1)
          deallocate(cam_out(c)%dstwet2)
          deallocate(cam_out(c)%dstdry2)
          deallocate(cam_out(c)%dstwet3)
          deallocate(cam_out(c)%dstdry3)
          deallocate(cam_out(c)%dstwet4)
          deallocate(cam_out(c)%dstdry4)
          deallocate(cam_out(c)%wsresp)
          deallocate(cam_out(c)%tau_est)
          deallocate(cam_out(c)%ugust)
          deallocate(cam_out(c)%uovern)
       enddo  

       deallocate(cam_out)
    end if
    nullify(cam_out)

  end subroutine atm2hub_deallocate
  subroutine hub2atm_deallocate(cam_in)
    type(cam_in_t), pointer :: cam_in(:)    ! Atmosphere to surface input
    integer :: c

    if(associated(cam_in)) then
       do c=begchunk,endchunk
          deallocate(cam_in(c)%asdir)
          deallocate(cam_in(c)%asdif)
          deallocate(cam_in(c)%aldir)
          deallocate(cam_in(c)%aldif)
          deallocate(cam_in(c)%lwup)
          deallocate(cam_in(c)%lhf)
          deallocate(cam_in(c)%shf)
          deallocate(cam_in(c)%h2otemp)
          deallocate(cam_in(c)%wsx)
          deallocate(cam_in(c)%wsy)
          deallocate(cam_in(c)%tref)
          deallocate(cam_in(c)%qref)
          deallocate(cam_in(c)%u10)
          deallocate(cam_in(c)%u10withgusts)
          deallocate(cam_in(c)%ts)
          deallocate(cam_in(c)%sst)
          deallocate(cam_in(c)%snowhland)
          deallocate(cam_in(c)%snowhice)
          deallocate(cam_in(c)%fco2_lnd)
          deallocate(cam_in(c)%fco2_ocn)
          deallocate(cam_in(c)%fdms)
          deallocate(cam_in(c)%landfrac)
          deallocate(cam_in(c)%icefrac)
          deallocate(cam_in(c)%ocnfrac)

          if(associated(cam_in(c)%ram1)) then
             deallocate(cam_in(c)%ram1)
             nullify(cam_in(c)%ram1)
          end if
          if(associated(cam_in(c)%fv)) then
             deallocate(cam_in(c)%fv)
             nullify(cam_in(c)%fv)
          end if
          if(associated(cam_in(c)%soilw)) then
             deallocate(cam_in(c)%soilw)
             nullify(cam_in(c)%soilw)
          end if

          deallocate(cam_in(c)%cflx)
          deallocate(cam_in(c)%ustar)
          deallocate(cam_in(c)%re)
          deallocate(cam_in(c)%ssq)

          if(associated(cam_in(c)%dstflx)) then
             deallocate(cam_in(c)%dstflx)
             nullify(cam_in(c)%dstflx)
          end if
          if(associated(cam_in(c)%meganflx)) then
             deallocate(cam_in(c)%meganflx)
             nullify(cam_in(c)%meganflx)
          end if
          if(associated(cam_in(c)%depvel)) then
             deallocate(cam_in(c)%depvel)
             nullify(cam_in(c)%depvel)
          end if
          
       enddo

       deallocate(cam_in)
    end if
    nullify(cam_in)

  end subroutine hub2atm_deallocate


!======================================================================

subroutine cam_export(state,cam_out,pbuf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transfer atmospheric fields into necessary surface data structures
! 
! Author: L. Bath  CMS Contact: M. Vertenstein
! 
!-----------------------------------------------------------------------
   use physics_types,    only: physics_state
   use ppgrid,           only: pver
   use cam_history,      only: outfld
   use comsrf,           only: psm1, srfrpdel, prcsnw
   use chem_surfvals,    only: chem_surfvals_get
   use co2_cycle,        only: co2_transport, c_i
   use physconst,        only: mwdry, mwco2
   use constituents,     only: pcnst
   use cam_control_mod,  only: rair
   use physics_buffer,   only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use phys_control,     only: phys_getopts
   implicit none

   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
   type(physics_state),  intent(in)    :: state
   type (cam_out_t),     intent(inout) :: cam_out
   type(physics_buffer_desc), pointer  :: pbuf(:)

   !
   !---------------------------Local variables-----------------------------
   !
   integer :: i              ! Longitude index
   integer :: m              ! constituent index
   integer :: lchnk          ! Chunk index
   integer :: ncol
   integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx,snow_sed_idx,prec_pcw_idx,snow_pcw_idx
   integer :: vmag_gust_idx, wsresp_idx, tau_est_idx
   real(r8) :: umb(pcols), vmb(pcols),vmag(pcols)
   logical :: linearize_pbl_winds ! Send wsresp and tau_est to coupler.

   real(r8), pointer :: prec_dp(:)                 ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)                 ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)                 ! snow from   Hack   convection
   real(r8), pointer :: prec_sed(:)                ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)                ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)                ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)                ! snow from Hack   convection
   real(r8), pointer :: vmag_gust(:)
   real(r8), pointer :: wsresp(:)                  ! First-order response of wind to surface stress
   real(r8), pointer :: tau_est(:)                 ! Estimated stress in equilibrium with ubot/vbot

   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   call phys_getopts(linearize_pbl_winds_out=linearize_pbl_winds)

   prec_dp_idx = pbuf_get_index('PREC_DP')
   snow_dp_idx = pbuf_get_index('SNOW_DP')
   prec_sh_idx = pbuf_get_index('PREC_SH')
   snow_sh_idx = pbuf_get_index('SNOW_SH')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')
   vmag_gust_idx = pbuf_get_index('vmag_gust')

   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
   call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)
   call pbuf_get_field(pbuf, vmag_gust_idx, vmag_gust)

   if (linearize_pbl_winds) then
      wsresp_idx = pbuf_get_index('wsresp')
      tau_est_idx = pbuf_get_index('tau_est')
      call pbuf_get_field(pbuf, wsresp_idx, wsresp)
      call pbuf_get_field(pbuf, tau_est_idx, tau_est)
   end if

!PMA adds gustiness to surface scheme c20181128

   do i=1,ncol
      cam_out%ubot(i)  = state%u(i,pver)
      cam_out%vbot(i)  = state%v(i,pver)
      cam_out%ugust(i) = vmag_gust(i)
      cam_out%tbot(i)  = state%t(i,pver)
      cam_out%thbot(i) = state%t(i,pver) * state%exner(i,pver)
      cam_out%zbot(i)  = state%zm(i,pver)
      cam_out%pbot(i)  = state%pmid(i,pver)
      cam_out%rho(i)   = cam_out%pbot(i)/(rair*cam_out%tbot(i))
      if (linearize_pbl_winds) then
         cam_out%wsresp(i)= max(wsresp(i), 0._r8)
         cam_out%tau_est(i)= tau_est(i)
      end if
      psm1(i,lchnk)    = state%ps(i)
      srfrpdel(i,lchnk)= state%rpdel(i,pver)
   end do
   do m = 1, pcnst
     do i = 1, ncol
        cam_out%qbot(i,m) = state%q(i,pver,m) 
     end do
   end do

   cam_out%co2diag(:ncol) = chem_surfvals_get('CO2VMR') * 1.0e+6_r8 
   if (co2_transport()) then
      do i=1,ncol
         cam_out%co2prog(i) = state%q(i,pver,c_i(4)) * 1.0e+6_r8 *mwdry/mwco2
      end do
   end if
   !
   ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
   ! Compute total convective and stratiform precipitation and snow rates
   !
   do i=1,ncol
      cam_out%precc (i) = prec_dp(i)  + prec_sh(i)
      cam_out%precl (i) = prec_sed(i) + prec_pcw(i)
      cam_out%precsc(i) = snow_dp(i)  + snow_sh(i)
      cam_out%precsl(i) = snow_sed(i) + snow_pcw(i)

      ! jrm These checks should not be necessary if they exist in the parameterizations
      if (cam_out%precc(i) .lt.0._r8) cam_out%precc(i)=0._r8
      if (cam_out%precl(i) .lt.0._r8) cam_out%precl(i)=0._r8
      if (cam_out%precsc(i).lt.0._r8) cam_out%precsc(i)=0._r8
      if (cam_out%precsl(i).lt.0._r8) cam_out%precsl(i)=0._r8
      if (cam_out%precsc(i).gt.cam_out%precc(i)) cam_out%precsc(i)=cam_out%precc(i)
      if (cam_out%precsl(i).gt.cam_out%precl(i)) cam_out%precsl(i)=cam_out%precl(i)
      ! end jrm
   end do

   ! total snowfall rate: needed by slab ocean model
   prcsnw(:ncol,lchnk) = cam_out%precsc(:ncol) + cam_out%precsl(:ncol)   

end subroutine cam_export

end module camsrfexch
