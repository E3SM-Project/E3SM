module vGCM_data_types

  USE shr_kind_mod, only: dbl_kind => shr_kind_r8
  USE parmsld,      only: nhalo_vgcm,nvgcm_seg,nlevel,ntracer
  
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
    real(kind=dbl_kind) :: alpha1,beta1  ! central angle coordinates of the 1st vGCM grid 
                                         ! i.e., at lon(0,1) & lat(0,1)
    
    real(kind=dbl_kind), allocatable :: lat(:,:)    ! vGCM grid latitude
    real(kind=dbl_kind), allocatable :: lon(:,:)    ! vGCM grid longitude
    real(kind=dbl_kind), allocatable :: TH(:,:,:)   ! potential temperature
    real(kind=dbl_kind), allocatable :: QV(:,:,:)   ! water vapor
    real(kind=dbl_kind), allocatable :: QC(:,:,:)   ! cloud liquid water
    real(kind=dbl_kind), allocatable :: QI(:,:,:)   ! cloud ice water
    real(kind=dbl_kind), allocatable :: QR(:,:,:)   ! cloud rain
    real(kind=dbl_kind), allocatable :: QS(:,:,:)   ! cloud snow
    real(kind=dbl_kind), allocatable :: QG(:,:,:)   ! cloud groupel
    real(kind=dbl_kind), allocatable :: QT(:,:,:,:) ! tracer array

    real(kind=dbl_kind), allocatable :: Ucon(:,:,:) ! contravariant zonal velocity component
    real(kind=dbl_kind), allocatable :: Vcon(:,:,:) ! contravariant meridional velocity component
    real(kind=dbl_kind), allocatable :: Ucov(:,:,:) ! covariant zonal velocity component
    real(kind=dbl_kind), allocatable :: Vcov(:,:,:) ! covariant meridional velocity component
    real(kind=dbl_kind), allocatable :: W(:,:,:)    ! vertical velocity
  end type vGCM_state_t

  type, public :: vGCM_tend_t
    real(kind=dbl_kind), allocatable :: dTH(:,:)    ! potential temperature tendency
    real(kind=dbl_kind), allocatable :: dQV(:,:)    ! water vapor tendency
    real(kind=dbl_kind), allocatable :: dQC(:,:)    ! cloud liquid water tendency
    real(kind=dbl_kind), allocatable :: dQI(:,:)    ! cloud ice water tendency
    real(kind=dbl_kind), allocatable :: dQR(:,:)    ! cloud rain tendency
    real(kind=dbl_kind), allocatable :: dQS(:,:)    ! cloud snow tendency
    real(kind=dbl_kind), allocatable :: dQG(:,:)    ! cloud graupel tendency
    real(kind=dbl_kind), allocatable :: dQT(:,:,:)  ! tracer tendencies

    real(kind=dbl_kind), allocatable :: dU(:,:)     ! zonal velocity tendency
    real(kind=dbl_kind), allocatable :: dV(:,:)     ! meridional velocity tencency
  end type vGCM_tend_t

  type, public :: vGCM_out_t
    ! diagnostic fields can be added here.
    real(kind=dbl_kind), allocatable :: SPREC(:)    ! surface precipitation
    real(kind=dbl_kind), allocatable :: WTH(:)      ! surface heat flux
    real(kind=dbl_kind), allocatable :: WQV(:)      ! surface moisture flux
    real(kind=dbl_kind), allocatable :: UW(:)       ! surface momentum flux (u"w")
    real(kind=dbl_kind), allocatable :: WV(:)       ! surface momentum flux (v"w")
  end type vGCM_out_t

!***********************************************************************************************
! Define a type (2): Combination of 4 channel segments (a channel)
!***********************************************************************************************

  type, public :: channel_vGCM_t
        type(vGCM_state_t), public :: vGCM_state(4)
        type(vGCM_tend_t),  public :: vGCM_tend(4)
        type(vGCM_out_t),   public :: vGCM_out(4)
  end type  channel_vGCM_t

!***********************************************************************************************
! Define variables with derived type
!***********************************************************************************************

  type(channel_vGCM_t), public, pointer :: channel_vGCM(:) => NULL()

!***********************************************************************************************
  
  public :: allocate_channel_vGCM_data

CONTAINS

  subroutine allocate_channel_vGCM_data(channel_vGCM, nch)

    ! Dummy argument
    integer,              intent(in)    :: nch               ! Number of channels to be allocated  
    type(channel_vGCM_t), intent(inout) :: channel_vGCM(nch) ! vGCM variables to be allocated

    ! Local data
    integer :: i,j

    do j = 1, nch
    do i = 1, 4
      call allocate_channel_vGCM_section(channel_vGCM(j)%vGCM_state(i), &
                                         channel_vGCM(j)%vGCM_tend(i),  &
                                         channel_vGCM(j)%vGCM_out(i))
    enddo
    enddo

  end subroutine allocate_channel_vGCM_data

  subroutine allocate_channel_vGCM_section (seg_state,seg_tend,seg_out)

    ! Dummy argument
    type(vGCM_state_t), intent(inout) :: seg_state  ! vGCM state variables to be allocated 
    type(vGCM_tend_t),  intent(inout) :: seg_tend   ! vGCM tendency variables to be allocated
    type(vGCM_out_t),   intent(inout) :: seg_out    ! vGCM output variables to be allocated
!----------------------------------------------------------------------- 
!   State Variables
!-----------------------------------------------------------------------        
    allocate(seg_state%lat(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM))
    allocate(seg_state%lon(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM))

    allocate(seg_state%TH(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QV(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QC(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QI(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QR(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QS(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QG(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))   
    allocate(seg_state%QT(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel,ntracer)) 

    allocate(seg_state%Ucon(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%Vcon(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%Ucov(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%Vcov(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))
    allocate(seg_state%W(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel))    

!----------------------------------------------------------------------- 
!   Tendency Variables
!-----------------------------------------------------------------------        
    allocate(seg_tend%dTH(nVGCM_seg, nLevel))
    allocate(seg_tend%dQV(nVGCM_seg, nLevel))
    allocate(seg_tend%dQC(nVGCM_seg, nLevel))
    allocate(seg_tend%dQI(nVGCM_seg, nLevel))
    allocate(seg_tend%dQR(nVGCM_seg, nLevel))
    allocate(seg_tend%dQS(nVGCM_seg, nLevel))
    allocate(seg_tend%dQG(nVGCM_seg, nLevel))
    allocate(seg_tend%dQT(nVGCM_seg, nLevel))

    allocate(seg_tend%dU(nVGCM_seg, nLevel))
    allocate(seg_tend%dV(nVGCM_seg, nLevel))

!----------------------------------------------------------------------- 
!   Output Variables
!-----------------------------------------------------------------------        
    allocate(seg_out%SPREC(nVGCM_seg))
    allocate(seg_out%WTH(nVGCM_seg))
    allocate(seg_out%WQV(nVGCM_seg))
    allocate(seg_out%UW(nVGCM_seg))
    allocate(seg_out%WV(nVGCM_seg))

  end subroutine allocate_channel_vGCM_section
  
end module vGCM_data_types
