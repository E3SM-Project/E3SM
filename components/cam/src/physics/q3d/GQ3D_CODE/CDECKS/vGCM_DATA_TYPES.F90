!#define single_run
! JUNG: comment out for use in CAM branch

module vGCM_data_types

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE parmsld,      only: nhalo_vgcm,nvgcm_seg,nlevel,ntracer,channel_seg_l

  implicit none
  private

! nhalo_vgcm: The size of vGCM halo (= 2 and declared at parmsld.F90)
! nvgcm_seg : The number of vGCM grids in a channel-segment
! nlevel    : The numer of vertical levels of vGCM (assumed to be same as nk1 for now)
! ntracer   : The numer of tracers

!***********************************************************************************************
! Define a type (1): vGCM variables declared for a channel-segment (covering 1/4 of a channel)
!***********************************************************************************************

  type, public :: vGCM_state_t
    integer :: nface   ! Number of face where the segment belongs
    real(kind=r8) :: alpha1,beta1  ! central angle coordinates of the 1st vGCM grid
                                   ! i.e., at lon(0,1) & lat(0,1)

    ! Top of atmosphere is index 1, nlayer is bottom (nlayer+1) for pint
    real(kind=r8), allocatable :: lat(:,:)      ! vGCM grid latitude
    real(kind=r8), allocatable :: lon(:,:)      ! vGCM grid longitude
    real(kind=r8), allocatable :: T(:,:,:)      ! temperature
    real(kind=r8), allocatable :: pint(:,:,:)   ! dry pressure at interfaces
    real(kind=r8), allocatable :: QV(:,:,:)     ! water vapor
    real(kind=r8), allocatable :: QC(:,:,:)     ! cloud liquid water
    real(kind=r8), allocatable :: QI(:,:,:)     ! cloud ice water
    real(kind=r8), allocatable :: QR(:,:,:)     ! cloud rain
    real(kind=r8), allocatable :: QS(:,:,:)     ! cloud snow
    real(kind=r8), allocatable :: QG(:,:,:)     ! cloud graupel
    real(kind=r8), allocatable :: QT(:,:,:,:)   ! tracer array

    real(kind=r8), allocatable :: U(:,:,:)      ! lon/lat zonal velocity component
    real(kind=r8), allocatable :: V(:,:,:)      ! lon/lat meridional velocity component

    real(kind=r8), allocatable :: omega(:,:,:)  ! vertical pressure velocity
    real(kind=r8), allocatable :: zm_int(:,:,:) ! geopotential height at interfaces
  end type vGCM_state_t
  
  type, public :: vGCM_map_t
    real(kind=r8), allocatable :: AM(:,:,:)     ! conversion matrix
    real(kind=r8), allocatable :: AMI(:,:,:)    ! inverse of conversion matrix
  end type vGCM_map_t

  type, public :: vGCM_tend_t
    real(kind=r8), allocatable :: dT(:,:)     ! temperature tendency
    real(kind=r8), allocatable :: dQV(:,:)    ! water vapor tendency
    real(kind=r8), allocatable :: dQC(:,:)    ! cloud liquid water tendency
    real(kind=r8), allocatable :: dQI(:,:)    ! cloud ice water tendency
    real(kind=r8), allocatable :: dQR(:,:)    ! cloud rain tendency
    real(kind=r8), allocatable :: dQS(:,:)    ! cloud snow tendency
    real(kind=r8), allocatable :: dQG(:,:)    ! cloud graupel tendency
    real(kind=r8), allocatable :: dQT(:,:,:)  ! tracer tendencies

    real(kind=r8), allocatable :: dU(:,:)     ! zonal velocity tendency
    real(kind=r8), allocatable :: dV(:,:)     ! meridional velocity tencency
  end type vGCM_tend_t

  type, public :: vGCM_out_t
    real(kind=r8), allocatable :: CLON(:)    ! Longitude at CRM q-point
    real(kind=r8), allocatable :: CLAT(:)    ! Latitude at CRM q-point

    ! diagnostic fields can be added here.
    real(kind=r8), allocatable :: SPREC(:)    ! surface precipitation
    real(kind=r8), allocatable :: WTH(:)      ! surface heat flux
    real(kind=r8), allocatable :: WQV(:)      ! surface moisture flux
    real(kind=r8), allocatable :: UW(:)       ! surface momentum flux (u"w")
    real(kind=r8), allocatable :: WV(:)       ! surface momentum flux (v"w")

    real(kind=r8), allocatable :: T(:,:)      ! temperature from the channel
    real(kind=r8), allocatable :: QV(:,:)     ! water vapor from the channel
    real(kind=r8), allocatable :: QC(:,:)     ! cloud liquid water from the channel
    real(kind=r8), allocatable :: QI(:,:)     ! cloud ice water from the channel
    real(kind=r8), allocatable :: QR(:,:)     ! cloud rain from the channel
    real(kind=r8), allocatable :: QS(:,:)     ! cloud snow from the channel
    real(kind=r8), allocatable :: QG(:,:)     ! cloud graupel from the channel
    real(kind=r8), allocatable :: QT(:,:,:)   ! tracer tendencies from the channel

    real(kind=r8), allocatable :: U(:,:)      ! zonal velocity from the channel
    real(kind=r8), allocatable :: V(:,:)      ! meridional velocity from the channel
    real(kind=r8), allocatable :: W(:,:)      ! vertical velocity from the channel
  end type vGCM_out_t

!***********************************************************************************************
! Define a type (2): Combination of 4 channel segments (a channel)
!***********************************************************************************************

  type, public :: channel_vGCM_t
        type(vGCM_state_t), public :: vGCM_state(4)
        type(vGCM_tend_t),  public :: vGCM_tend(4)
        type(vGCM_out_t),   public :: vGCM_out(4)
        type(vGCM_map_t),   public :: vGCM_map(4)
  end type  channel_vGCM_t

!******************************************************************************
! Define variables with derived type
!******************************************************************************

  type(channel_vGCM_t), public, pointer :: channel_vGCM(:) => NULL()


!***********************************************************************************************
! Define a type (3): Used in q3d_bg_module
!***********************************************************************************************

!------------------------------------------------------------------------------


!******************************************************************************
! Public interface subroutines
!******************************************************************************

  public :: vGCM_allocate_channel_data

!******************************************************************************

CONTAINS

  subroutine vGCM_allocate_channel_data(begchan, endchan)

#ifndef single_run
    use cam_abortutils,   only: endrun
#endif

    ! Dummy argument
    integer,              intent(in) :: begchan ! First channel to allocate
    integer,              intent(in) :: endchan ! Last channel to allocate

    ! Local data
    integer :: i,j

#ifndef single_run
    if (associated(channel_vGCM)) then
       call endrun("vGCM_allocate_channel_data: channel_vGCM already allocated")
    end if
#endif

    allocate(channel_vGCM(begchan:endchan))
    do j = begchan, endchan
       do i = 1, 4
          call allocate_channel_vGCM_section(channel_vGCM(j)%vGCM_state(i), &
                 channel_vGCM(j)%vGCM_tend(i), channel_vGCM(j)%vGCM_out(i), &
                 channel_vGCM(j)%vGCM_map(i))
       end do
    end do

  end subroutine vGCM_allocate_channel_data

  subroutine allocate_channel_vGCM_section (seg_state,seg_tend,seg_out,seg_map)
#ifndef single_run
     use physconst, only: pi
#endif
    ! Dummy argument
    type(vGCM_state_t), intent(inout) :: seg_state  ! vGCM state variables to be allocated
    type(vGCM_tend_t),  intent(inout) :: seg_tend   ! vGCM tendency variables to be allocated
    type(vGCM_out_t),   intent(inout) :: seg_out    ! vGCM output variables to be allocated
    type(vGCM_map_t),   intent(inout) :: seg_map    ! vGCM map variables to be allocated

    INTEGER numCRM
!-----------------------------------------------------------------------
!   State Variables
!-----------------------------------------------------------------------
    allocate(seg_state%lat(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM))
    allocate(seg_state%lon(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM))
#ifndef single_run
    ! Initialize lat and lon to impossible values to indicate unused columns
    seg_state%lat(:,:) = pi
    seg_state%lon(:,:) = 3.0_r8 * pi
#endif

    allocate(seg_state%T(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%pint(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel+1))
    allocate(seg_state%QV(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%QC(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%QI(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%QR(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%QS(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%QG(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%QT(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel,ntracer))

    allocate(seg_state%U(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%V(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))

    allocate(seg_state%omega(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%zm_int(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel+1))

!-----------------------------------------------------------------------
!   Map Variables
!-----------------------------------------------------------------------
    allocate(seg_map%AM(4,-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM))
    allocate(seg_map%AMI(4,-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM))
    
!-----------------------------------------------------------------------
!   Tendency Variables
!-----------------------------------------------------------------------
    allocate(seg_tend%dT(nVGCM_seg, nLevel))
    allocate(seg_tend%dQV(nVGCM_seg, nLevel))
    allocate(seg_tend%dQC(nVGCM_seg, nLevel))
    allocate(seg_tend%dQI(nVGCM_seg, nLevel))
    allocate(seg_tend%dQR(nVGCM_seg, nLevel))
    allocate(seg_tend%dQS(nVGCM_seg, nLevel))
    allocate(seg_tend%dQG(nVGCM_seg, nLevel))
    allocate(seg_tend%dQT(nVGCM_seg, nLevel, ntracer))

    allocate(seg_tend%dU(nVGCM_seg, nLevel))
    allocate(seg_tend%dV(nVGCM_seg, nLevel))

!-----------------------------------------------------------------------
!   Output Variables
!-----------------------------------------------------------------------
    allocate(seg_out%CLON(channel_seg_l))
    allocate(seg_out%CLAT(channel_seg_l))

    allocate(seg_out%SPREC(nVGCM_seg))
    allocate(seg_out%WTH(nVGCM_seg))
    allocate(seg_out%WQV(nVGCM_seg))
    allocate(seg_out%UW(nVGCM_seg))
    allocate(seg_out%WV(nVGCM_seg))

    allocate(seg_out%T(nVGCM_seg, nLevel))
    allocate(seg_out%QV(nVGCM_seg, nLevel))
    allocate(seg_out%QC(nVGCM_seg, nLevel))
    allocate(seg_out%QI(nVGCM_seg, nLevel))
    allocate(seg_out%QR(nVGCM_seg, nLevel))
    allocate(seg_out%QS(nVGCM_seg, nLevel))
    allocate(seg_out%QG(nVGCM_seg, nLevel))
    allocate(seg_out%QT(nVGCM_seg, nLevel, ntracer))

    allocate(seg_out%U(nVGCM_seg, nLevel))
    allocate(seg_out%V(nVGCM_seg, nLevel))
    allocate(seg_out%W(nVGCM_seg, nLevel))

  end subroutine allocate_channel_vGCM_section

end module vGCM_data_types
