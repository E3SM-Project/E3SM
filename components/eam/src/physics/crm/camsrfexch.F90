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
  use phys_grid,     only: get_ncols_p, get_gcol_p,get_rlat_p,get_rlon_p, phys_grid_initialized
  use infnan,        only: posinf, assignment(=)
  use cam_abortutils,    only: endrun
  use cam_logfile,   only: iulog
  use seq_comm_mct, only : num_inst_atm
  !! Jungmin
  use spmd_utils,         only: masterproc
  !! Jungmin
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
     
     ! surface interface fields beteen ATM and LND/ICE components have an added
     ! dimension (num_inst_atm) for MMF-MAML configuration 
     real(r8), allocatable :: tbot_mi(:,:)     ! bot level temperature
     real(r8), allocatable :: thbot_mi(:,:)    ! bot level potential temperature
     real(r8), allocatable :: ubot_mi(:,:)     ! bot level u wind
     real(r8), allocatable :: vbot_mi(:,:)     ! bot level v wind
     real(r8), allocatable :: qbot_mi(:,:,:)   ! bot level specific humidity
     real(r8), allocatable :: rho_mi(:,:)      ! bot level density	
     real(r8), allocatable :: precsc_mi(:,:)   ! convective snow rate 
     real(r8), allocatable :: precsl_mi(:,:)   ! large scale snow rate
     real(r8), allocatable :: precc_mi(:,:)    ! convective precip. rate
     real(r8), allocatable :: precl_mi(:,:)    ! large-scale precip.rate
     real(r8), allocatable :: soll_mi(:,:)     ! near-IR, direct, down
     real(r8), allocatable :: sols_mi(:,:)     ! visible, direct, down
     real(r8), allocatable :: solld_mi(:,:)    ! near-IR, diffuse, down
     real(r8), allocatable :: solsd_mi(:,:)    ! visible, diffuse, down
     real(r8), allocatable :: flwds_mi(:,:)    ! longwave, down

     real(r8), allocatable :: tbot(:)     ! bot level temperature
     real(r8), allocatable :: thbot(:)    ! bot level potential temperature
     real(r8), allocatable :: ubot(:)     ! bot level u wind
     real(r8), allocatable :: vbot(:)     ! bot level v wind
     real(r8), allocatable :: qbot(:,:)   ! bot level specific humidity
     real(r8), allocatable :: rho(:)      ! bot level density	
     real(r8), allocatable :: zbot(:)      ! bot level height above surface	
     real(r8), allocatable :: pbot(:)      ! bot level pressure
     real(r8), allocatable :: precsc(:)   ! convective snow rate 
     real(r8), allocatable :: precsl(:)   ! large scale snow rate
     real(r8), allocatable :: precc(:)    ! convective precip. rate
     real(r8), allocatable :: precl(:)    ! large-scale precip.rate
     real(r8), allocatable :: soll(:)     ! near-IR, direct, down
     real(r8), allocatable :: sols(:)     ! visible, direct, down
     real(r8), allocatable :: solld(:)    ! near-IR, diffuse, down
     real(r8), allocatable :: solsd(:)    ! visible, diffuse, down
     real(r8), allocatable :: netsw(:)    ! surface solar, net absorbed	
     real(r8), allocatable :: flwds(:)    ! longwave, down
     
     real(r8), allocatable :: co2prog(:)    ! prognostic co2
     real(r8), allocatable :: co2diag(:)    ! diagnostic co2
     real(r8), allocatable :: psl(:)
     real(r8), allocatable :: bcphiwet(:)   ! wet deposition of hydrophilic black carbon
     real(r8), allocatable :: bcphidry(:)   ! dry deposition of hydrophilic black carbon
     real(r8), allocatable :: bcphodry(:)   ! dry deposition of hydrophobic black carbon
     real(r8), allocatable :: ocphiwet(:)   ! wet deposition of hydrophilic organic carbon
     real(r8), allocatable :: ocphidry(:)   ! dry deposition of hydrophilic organic carbon
     real(r8), allocatable :: ocphodry(:)   ! dry deposition of hydrophobic organic carbon
     real(r8), allocatable :: dstwet1(:)    ! wet deposition of dust (bin1)
     real(r8), allocatable :: dstdry1(:)    ! dry deposition of dust (bin1)
     real(r8), allocatable :: dstwet2(:)    ! wet deposition of dust (bin2)
     real(r8), allocatable :: dstdry2(:)    ! dry deposition of dust (bin2)
     real(r8), allocatable :: dstwet3(:)    ! wet deposition of dust (bin3)
     real(r8), allocatable :: dstdry3(:)    ! dry deposition of dust (bin3)
     real(r8), allocatable :: dstwet4(:)    ! wet deposition of dust (bin4)
     real(r8), allocatable :: dstdry4(:)    ! dry deposition of dust (bin4)
  end type cam_out_t 

!---------------------------------------------------------------------------
! This is the merged state of sea-ice, land and ocean surface parameterizations
!---------------------------------------------------------------------------

  type cam_in_t    
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns

     ! surface interface fields beteen ATM and LND/ICE components have an added
     ! dimension (num_inst_atm) for MMF-MAML configuration 
     real(r8), allocatable :: asdir_mi(:,:)      ! albedo: shortwave, direct
     real(r8), allocatable :: asdif_mi(:,:)      ! albedo: shortwave, diffuse
     real(r8), allocatable :: aldir_mi(:,:)      ! albedo: longwave, direct
     real(r8), allocatable :: aldif_mi(:,:)      ! albedo: longwave, diffuse
     real(r8), allocatable :: lwup_mi(:,:)       ! longwave up radiative flux
     real(r8), allocatable :: lhf_mi(:,:)        ! latent heat flux
     real(r8), allocatable :: cflx1_mi(:,:)        ! water vapor constituent flux 
     real(r8), allocatable :: shf_mi(:,:)        ! sensible heat flux
     real(r8), allocatable :: wsx_mi(:,:)        ! surface u-stress (N)
     real(r8), allocatable :: wsy_mi(:,:)        ! surface v-stress (N)
     real(r8), allocatable :: snowhland_mi(:,:)  ! snow depth (liquid water equivalent) over land
     real(r8), allocatable :: ts_mi(:,:)       ! surface temperature for MMF-MAML configuration
     
     ! Average across the instances for cam output
     real(r8), allocatable :: asdir(:)    ! albedo: shortwave, direct
     real(r8), allocatable :: asdif(:)      ! albedo: shortwave, diffuse
     real(r8), allocatable :: aldir(:)      ! albedo: longwave, direct
     real(r8), allocatable :: aldif(:)    ! albedo: longwave, diffuse
     real(r8), allocatable :: lwup(:)    ! longwave up radiative flux
     real(r8), allocatable :: lhf(:)        ! latent heat flux
     real(r8), allocatable :: shf(:)        ! sensible heat flux
     real(r8), allocatable :: wsx(:)        ! surface u-stress (N)
     real(r8), allocatable :: wsy(:)        ! surface v-stress (N)
     real(r8), allocatable :: snowhland(:)  ! snow depth (liquid water equivalent) over land
     
     real(r8), allocatable :: tref(:)         ! ref height surface air temp
     real(r8), allocatable :: qref(:)         ! ref height specific humidity 
     real(r8), allocatable :: u10(:)          ! 10m wind speed
     real(r8), allocatable :: ts(:)           ! merged surface temp 
     real(r8), allocatable :: sst(:)          ! sea surface temp
     real(r8), allocatable :: snowhice(:)     ! snow depth over ice
     real(r8), allocatable :: fco2_lnd(:)     ! co2 flux from lnd
     real(r8), allocatable :: fco2_ocn(:)     ! co2 flux from ocn
     real(r8), allocatable :: fdms(:)         ! dms flux
     real(r8), allocatable :: landfrac(:)     ! land area fraction
     real(r8), allocatable :: icefrac(:)      ! sea-ice areal fraction
     real(r8), allocatable :: ocnfrac(:)      ! ocean areal fraction
     real(r8), pointer, dimension(:) :: ram1  ! aerodynamical resistance (s/m) (pcols)
     real(r8), pointer, dimension(:) :: fv    ! friction velocity (m/s) (pcols)
     real(r8), pointer, dimension(:) :: soilw ! volumetric soil water (m3/m3)
     real(r8), allocatable :: cflx(:,:)       ! constituent flux (emissions)
     real(r8), allocatable :: ustar(:)        ! atm/ocn saved version of ustar
     real(r8), allocatable :: re(:)           ! atm/ocn saved version of re
     real(r8), allocatable :: ssq(:)          ! atm/ocn saved version of ssq
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
       allocate (cam_in(c)%asdir_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error asdir_mi')

       allocate (cam_in(c)%asdif_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error asdif_mi')

       allocate (cam_in(c)%aldir_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error aldir_mi')

       allocate (cam_in(c)%aldif_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error aldif_mi')

       allocate (cam_in(c)%lwup_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error lwup_mi')

       allocate (cam_in(c)%lhf_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error lhf_mi')
       
       allocate (cam_in(c)%cflx1_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error cflx1_mi')

       allocate (cam_in(c)%shf_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error shf_mi')

       allocate (cam_in(c)%wsx_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error wsx_mi')

       allocate (cam_in(c)%wsy_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error wsy_mi')

       allocate (cam_in(c)%snowhland_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error snowhland_mi')
       
       allocate (cam_in(c)%ts_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ts_mi')
       
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

       allocate (cam_in(c)%wsx(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error wsx')

       allocate (cam_in(c)%wsy(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error wsy')

       allocate (cam_in(c)%snowhland(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error snowhland')
       
       allocate (cam_in(c)%tref(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error tref')

       allocate (cam_in(c)%qref(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error qref')

       allocate (cam_in(c)%u10(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error u10')

       allocate (cam_in(c)%ts(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ts')

       allocate (cam_in(c)%sst(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error sst')

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
       cam_in(c)%asdir_mi    (:,:) = 0._r8
       cam_in(c)%asdif_mi    (:,:) = 0._r8
       cam_in(c)%aldir_mi    (:,:) = 0._r8
       cam_in(c)%aldif_mi    (:,:) = 0._r8
       cam_in(c)%lwup_mi     (:,:) = 0._r8
       cam_in(c)%lhf_mi      (:,:) = 0._r8
       cam_in(c)%cflx1_mi    (:,:) = 0._r8
       cam_in(c)%shf_mi      (:,:) = 0._r8
       cam_in(c)%wsx_mi      (:,:) = 0._r8
       cam_in(c)%wsy_mi      (:,:) = 0._r8
       cam_in(c)%snowhland_mi(:,:) = 0._r8
       cam_in(c)%ts_mi       (:,:) = 0._r8
       cam_in(c)%asdir    (:) = 0._r8
       cam_in(c)%asdif    (:) = 0._r8
       cam_in(c)%aldir    (:) = 0._r8
       cam_in(c)%aldif    (:) = 0._r8
       cam_in(c)%lwup     (:) = 0._r8
       cam_in(c)%lhf      (:) = 0._r8
       cam_in(c)%shf      (:) = 0._r8
       cam_in(c)%wsx      (:) = 0._r8
       cam_in(c)%wsy      (:) = 0._r8
       cam_in(c)%snowhland(:) = 0._r8
       cam_in(c)%tref     (:) = 0._r8
       cam_in(c)%qref     (:) = 0._r8
       cam_in(c)%u10      (:) = 0._r8
       cam_in(c)%ts       (:) = 0._r8
       cam_in(c)%sst      (:) = 0._r8
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
    !! Jungmin
    !if(masterproc) write(iulog,*) 'atm2hub_alloc: num_inst_atm:',num_inst_atm
    !! Jungmin
    if ( .not. phys_grid_initialized() ) call endrun( "ATM2HUB_ALLOC error: phys_grid not called yet" )
    allocate (cam_out(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('ATM2HUB_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk 
       allocate (cam_out(c)%tbot_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error tbot_mi')
       
       allocate (cam_out(c)%ubot_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error ubot_mi')

       allocate (cam_out(c)%vbot_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error vbot_mi')

       allocate (cam_out(c)%qbot_mi(pcols,pcnst,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error qbot_mi')

       allocate (cam_out(c)%rho_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error rho_mi')

       allocate (cam_out(c)%flwds_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error flwds_mi')

       allocate (cam_out(c)%precsc_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precsc_mi')

       allocate (cam_out(c)%precsl_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precsl_mi')

       allocate (cam_out(c)%precc_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precc_mi')

       allocate (cam_out(c)%precl_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error precl_mi')

       allocate (cam_out(c)%soll_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error soll_mi')

       allocate (cam_out(c)%sols_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error sols_mi')

       allocate (cam_out(c)%solld_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error solld_mi')

       allocate (cam_out(c)%solsd_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error solsd_mi')

       allocate (cam_out(c)%thbot_mi(pcols,num_inst_atm), stat=ierror)
       if ( ierror /= 0 ) call endrun('ATM2HUB_ALLOC error: allocation error thbot_mi')
       
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
    enddo  

    do c = begchunk,endchunk
       cam_out(c)%lchnk       = c
       cam_out(c)%ncol        = get_ncols_p(c)
       cam_out(c)%thbot_mi(:,:)    = 0._r8
       cam_out(c)%tbot_mi(:,:)     = 0._r8
       cam_out(c)%ubot_mi(:,:)     = 0._r8
       cam_out(c)%vbot_mi(:,:)     = 0._r8
       cam_out(c)%qbot_mi(:,:,:)   = 0._r8
       cam_out(c)%rho_mi(:,:)      = 0._r8
       cam_out(c)%precsc_mi(:,:)   = 0._r8
       cam_out(c)%precsl_mi(:,:)   = 0._r8
       cam_out(c)%precc_mi(:,:)    = 0._r8
       cam_out(c)%precl_mi(:,:)    = 0._r8
       cam_out(c)%soll_mi(:,:)     = 0._r8
       cam_out(c)%sols_mi(:,:)     = 0._r8
       cam_out(c)%solld_mi(:,:)    = 0._r8
       cam_out(c)%solsd_mi(:,:)    = 0._r8
       cam_out(c)%flwds_mi(:,:)    = 0._r8
       cam_out(c)%thbot(:)    = 0._r8
       cam_out(c)%tbot(:)     = 0._r8
       cam_out(c)%zbot(:)     = 0._r8
       cam_out(c)%ubot(:)     = 0._r8
       cam_out(c)%vbot(:)     = 0._r8
       cam_out(c)%qbot(:,:)   = 0._r8
       cam_out(c)%pbot(:)     = 0._r8
       cam_out(c)%rho(:)      = 0._r8
       cam_out(c)%precsc(:)   = 0._r8
       cam_out(c)%precsl(:)   = 0._r8
       cam_out(c)%precc(:)    = 0._r8
       cam_out(c)%precl(:)    = 0._r8
       cam_out(c)%soll(:)     = 0._r8
       cam_out(c)%sols(:)     = 0._r8
       cam_out(c)%solld(:)    = 0._r8
       cam_out(c)%solsd(:)    = 0._r8
       cam_out(c)%flwds(:)    = 0._r8
       cam_out(c)%netsw(:)    = 0._r8
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
    end do

  end subroutine atm2hub_alloc

  subroutine atm2hub_deallocate(cam_out)
    type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
    integer :: c

    if(associated(cam_out)) then
       do c = begchunk,endchunk
          deallocate(cam_out(c)%tbot_mi)
          deallocate(cam_out(c)%ubot_mi)
          deallocate(cam_out(c)%vbot_mi)
          deallocate(cam_out(c)%qbot_mi)
          deallocate(cam_out(c)%rho_mi)
          deallocate(cam_out(c)%flwds_mi)
          deallocate(cam_out(c)%precsc_mi)
          deallocate(cam_out(c)%precsl_mi)
          deallocate(cam_out(c)%precc_mi)
          deallocate(cam_out(c)%precl_mi)
          deallocate(cam_out(c)%soll_mi)
          deallocate(cam_out(c)%sols_mi)
          deallocate(cam_out(c)%solld_mi)
          deallocate(cam_out(c)%solsd_mi)
          deallocate(cam_out(c)%thbot_mi)
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
          deallocate(cam_in(c)%asdir_mi)
          deallocate(cam_in(c)%asdif_mi)
          deallocate(cam_in(c)%aldir_mi)
          deallocate(cam_in(c)%aldif_mi)
          deallocate(cam_in(c)%lwup_mi)
          deallocate(cam_in(c)%lhf_mi)
          deallocate(cam_in(c)%cflx1_mi)
          deallocate(cam_in(c)%shf_mi)
          deallocate(cam_in(c)%wsx_mi)
          deallocate(cam_in(c)%wsy_mi)
          deallocate(cam_in(c)%snowhland_mi)
          deallocate(cam_in(c)%ts_mi)
          deallocate(cam_in(c)%asdir)
          deallocate(cam_in(c)%asdif)
          deallocate(cam_in(c)%aldir)
          deallocate(cam_in(c)%aldif)
          deallocate(cam_in(c)%lwup)
          deallocate(cam_in(c)%lhf)
          deallocate(cam_in(c)%shf)
          deallocate(cam_in(c)%wsx)
          deallocate(cam_in(c)%wsy)
          deallocate(cam_in(c)%tref)
          deallocate(cam_in(c)%qref)
          deallocate(cam_in(c)%u10)
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

subroutine cam_export(state,cam_out,cam_in,pbuf)

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
   !! Jungmin
   use shr_const_mod,  only: shr_const_pi
   use cam_instance     , only: inst_index
   !! Jungmin
   !use maml_module,      only: cam_out_avg_mi
   implicit none

   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
   type(physics_state),  intent(in)    :: state
   type (cam_out_t),     intent(inout) :: cam_out
   type (cam_in_t),     intent(in) :: cam_in
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
   integer :: vmag_gust_idx
   real(r8) :: umb(pcols), vmb(pcols),vmag(pcols)
    

   real(r8), pointer :: prec_dp(:)                 ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)                 ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)                 ! snow from   Hack   convection
   real(r8), pointer :: prec_sed(:)                ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)                ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)                ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)                ! snow from Hack   convection
   real(r8), pointer :: vmag_gust(:)

   ! Variables defined below are for MMF-MAML
   integer :: j             ! Multi-instance index
   integer :: crm_t_idx
   integer :: crm_qv_idx
   integer :: crm_u_idx
   integer :: crm_v_idx
   integer :: crm_pcp_idx
   integer :: crm_snw_idx
   integer :: crm_angle_idx
   real(r8), pointer :: crm_t(:,:,:,:)
   real(r8), pointer :: crm_qv(:,:,:,:)
   real(r8), pointer :: crm_u(:,:,:,:)
   real(r8), pointer :: crm_v(:,:,:,:)
   real(r8), pointer :: crm_pcp(:,:,:)
   real(r8), pointer :: crm_snw(:,:,:)
   real(r8), pointer :: crm_angle(:)
   logical :: use_MAML
   real(r8) :: avgfac
   integer :: iter_max(10)

   !! Jungmin
   real(r8) :: clat,clon
   integer :: rcol
   !-----------------------------------------------------------------------
   iter_max = 10

   lchnk = state%lchnk
   ncol  = state%ncol

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

   call phys_getopts(use_MAML_out  = use_MAML)
   
   cam_out%pbot(1:ncol)         = state%pmid(1:ncol,pver)
   cam_out%zbot(1:ncol)         = state%zm  (1:ncol,pver)
   cam_out%qbot(1:ncol,2:pcnst) = state%q   (1:ncol,pver,2:pcnst)
   psm1        (1:ncol,lchnk)   = state%ps(1:ncol)
   srfrpdel    (1:ncol,lchnk)   = state%rpdel(1:ncol,pver)
   cam_out%co2diag(1:ncol)      = chem_surfvals_get('CO2VMR') * 1.0e+6_r8 
   if (co2_transport()) then
      cam_out%co2prog(1:ncol) = state%q(1:ncol,pver,c_i(4)) * 1.0e+6_r8 *mwdry/mwco2
   end if
   !! Jungmin
   !write(iulog,'("num_inst_atm:",I)') num_inst_atm
   !! Jungmin
   if(.not.use_MAML) then
      ! for default MMF, without MAML
   
      !PMA adds gustiness to surface scheme c20181128

      do i=1,ncol
         umb(i)           = state%u(i,pver)
         vmb(i)           = state%v(i,pver)
         vmag(i)          = max(1.e-5_r8,sqrt( umb(i)**2._r8 + vmb(i)**2._r8))
         cam_out%tbot(i)  = state%t(i,pver)
         cam_out%thbot(i) = state%t(i,pver) * state%exner(i,pver)
         cam_out%qbot(i,1)= state%q(i,pver,1) 
         cam_out%ubot(i)  = state%u(i,pver) * ((vmag_gust(i)+vmag(i))/vmag(i))
         cam_out%vbot(i)  = state%v(i,pver) * ((vmag_gust(i)+vmag(i))/vmag(i))
         cam_out%rho(i)   = cam_out%pbot(i)/(rair*cam_out%tbot(i))
      !
      ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
      ! Compute total convective and stratiform precipitation and snow rates
      !
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
      ! Even if use_MAML=false, cam_out%[X]_mi are still used as a interface between CRM atmosphere and land surface. 
      ! therefore, copy all manually
      do i = 1, ncol
         cam_out%tbot_mi(i,:) = cam_out%tbot(i)
         cam_out%thbot_mi(i,:) = cam_out%thbot(i)
         cam_out%ubot_mi(i,:) = cam_out%ubot(i)
         cam_out%vbot_mi(i,:) = cam_out%vbot(i)
         do m = 1,pcnst
            cam_out%qbot_mi(i,m,:) = cam_out%qbot(i,m)
         end do
         cam_out%rho_mi(i,:) = cam_out%rho(i)
         cam_out%precsc_mi(i,:) = cam_out%precsc(i)
         cam_out%precsl_mi(i,:) = cam_out%precsl(i)
         cam_out%precc_mi(i,:) = cam_out%precc(i)
         cam_out%precl_mi(i,:) = cam_out%precl(i)
         
         !! Jungmin
         !if(masterproc) then
         !   write(iulog,*) 'i:',i,'landfrac:',cam_in%landfrac(i)
         !   write(iulog,*),"cam_export:",i,' cam_out%tbot:',cam_out%tbot(i),'cam_out%tbot_mi:', (cam_out%tbot_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%ubot:',cam_out%ubot(i),'cam_out%ubot_mi:', (cam_out%ubot_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%vbot:',cam_out%vbot(i),'cam_out%vbot_mi:', (cam_out%vbot_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%rho:',cam_out%rho(i),'cam_out%rho_mi:', (cam_out%rho_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%precsc:',cam_out%precsc(i),'cam_out%precsc_mi:', (cam_out%precsc_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%precsl:',cam_out%precsl(i),'cam_out%precsl_mi:', (cam_out%precsl_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%precc:',cam_out%precc(i),'cam_out%precc_mi:', (cam_out%precc_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%precl:',cam_out%precl(i),'cam_out%precl_mi:', (cam_out%precl_mi(i,j),j=1,num_inst_atm)
         !   write(iulog,*),"cam_export:",i,' cam_out%qbot1:',cam_out%qbot(i,1),'cam_out%qbot1_mi:', (cam_out%qbot_mi(i,1,j),j=1,num_inst_atm)
         !end if
         !! Jungmin
      end do! i = 1, ncol   
   else   
      ! for MMF-MAML, surfaces are coupled to CRM atmosphere 
      crm_t_idx     = pbuf_get_index('CRM_T')
      crm_qv_idx    = pbuf_get_index('CRM_QV')
      crm_u_idx     = pbuf_get_index('CRM_U')
      crm_v_idx     = pbuf_get_index('CRM_V')
      crm_pcp_idx   = pbuf_get_index('CRM_PCP')
      crm_snw_idx   = pbuf_get_index('CRM_SNW')
      crm_angle_idx = pbuf_get_index('CRM_ANGLE')

      call pbuf_get_field(pbuf, crm_t_idx    , crm_t)
      call pbuf_get_field(pbuf, crm_qv_idx   , crm_qv)
      call pbuf_get_field(pbuf, crm_u_idx    , crm_u)
      call pbuf_get_field(pbuf, crm_v_idx    , crm_v)
      call pbuf_get_field(pbuf, crm_pcp_idx  , crm_pcp)
      call pbuf_get_field(pbuf, crm_snw_idx  , crm_snw)
      call pbuf_get_field(pbuf, crm_angle_idx, crm_angle)
      
      cam_out%precl_mi (1:ncol,1:num_inst_atm) = 0._r8     ! large-scale precip set to zero
      cam_out%precsl_mi(1:ncol,1:num_inst_atm) = 0._r8     ! large-scale precip set to zero

      do j=1,num_inst_atm
         cam_out%tbot_mi  (1:ncol,j)   = crm_t(1:ncol,j,1,1)
         cam_out%thbot_mi (1:ncol,j)   = cam_out%tbot_mi(1:ncol,j) * state%exner(1:ncol,pver) ! potential temperature
         cam_out%qbot_mi  (1:ncol,1,j) = crm_qv(1:ncol,j,1,1)
         ! u and v will use CRM value (must transform because of CRM orientation)
         cam_out%ubot_mi  (1:ncol,j)   = crm_u(1:ncol,j,1,1) * cos(crm_angle(1:ncol)) - crm_v(1:ncol,j,1,1) * sin(crm_angle(1:ncol))
         cam_out%vbot_mi  (1:ncol,j)   = crm_v(1:ncol,j,1,1) * cos(crm_angle(1:ncol)) + crm_u(1:ncol,j,1,1) * sin(crm_angle(1:ncol))
         cam_out%precc_mi (1:ncol,j)   = crm_pcp(1:ncol,j,1)  
         cam_out%precsc_mi(1:ncol,j)   = crm_snw(1:ncol,j,1)
         cam_out%rho_mi   (1:ncol,j)   = cam_out%pbot(1:ncol)/(rair*cam_out%tbot_mi(1:ncol,j))  ! air density
      end do 

      ! Average cam_out%[X]_mi across num_inst_atm and pass it to cam_out%[X]
      ! Data transfer between surface and CAM atmosphere will be done by cam_out%[X]
      ! Data transfer between surface and CRM atmosphere will be done by cam_out%[X]_mi
      ! TODO: call cam_out_avg_mi( cam_out )  
      avgfac = 1._r8/real(num_inst_atm,r8)
      ! Initialize again just to be sure
      cam_out%tbot = 0.
      cam_out%thbot = 0.
      cam_out%qbot(:,1) = 0.
      cam_out%precsc = 0.
      cam_out%precsl = 0.
      cam_out%precc = 0.
      cam_out%precl = 0.
      cam_out%ubot = 0.
      cam_out%vbot = 0.

      ! for non-MAML, num_inst_atm = avgfac = 1
      do i = 1,ncol
         do j = 1,num_inst_atm
            cam_out%tbot(i)   = cam_out%tbot(i)   + cam_out%tbot_mi(i,j)  *avgfac
            cam_out%thbot(i)  = cam_out%thbot(i)  + cam_out%thbot_mi(i,j) *avgfac
            cam_out%qbot(i,1) = cam_out%qbot(i,1) + cam_out%qbot_mi(i,1,j)*avgfac
            cam_out%ubot(i)   = cam_out%ubot(i)   + cam_out%ubot_mi(i,j)  *avgfac
            cam_out%vbot(i)   = cam_out%vbot(i)   + cam_out%vbot_mi(i,j)  *avgfac
            cam_out%rho(i)    = cam_out%rho(i)    + cam_out%rho_mi(i,j)   *avgfac
            cam_out%precsc(i) = cam_out%precsc(i) + cam_out%precsc_mi(i,j)*avgfac
            cam_out%precsl(i) = cam_out%precsl(i) + cam_out%precsl_mi(i,j)*avgfac
            cam_out%precc(i)  = cam_out%precc(i)  + cam_out%precc_mi(i,j) *avgfac
            cam_out%precl(i)  = cam_out%precl(i)  + cam_out%precl_mi(i,j) *avgfac
         end do
      end do! i = 1, ncol
      !! Jungmin
      do i = 1,ncol
         rcol = get_gcol_p(lchnk,i)         
         if(rcol.eq.223) then
            clat = get_rlat_p(lchnk,i)*180_r8/SHR_CONST_PI
            clon = get_rlon_p(lchnk,i)*180_r8/SHR_CONST_PI

            write(iulog,'("ATM_EXPORT: chunk_index=",I3," icol=",I3," rcol=",I5," inst_index=",I5, &
                          "lat=",F8.3," lon=",F8.3, &
                          " tbot=",F7.3)') &
                          lchnk,i,rcol,inst_index,&
                          clat,clon,&
                          cam_out%tbot(i)
            do j = 1, num_inst_atm
               write(iulog,'("      j=",I3," tbot_mi=",F7.3)')j,cam_out%tbot_mi(i,j)
            end do
         end if   
      end do
      !! Jungmin
   end if ! .not.use_MAML
   !
   ! total snowfall rate: needed by slab ocean model
   prcsnw(:ncol,lchnk) = cam_out%precsc(:ncol) + cam_out%precsl(:ncol)
end subroutine cam_export

end module camsrfexch
