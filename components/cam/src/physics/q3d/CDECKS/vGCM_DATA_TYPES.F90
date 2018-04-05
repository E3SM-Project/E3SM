module vGCM_data_types

  USE shr_kind_mod, only: r8 => shr_kind_r8
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
    real(kind=r8) :: alpha1,beta1  ! central angle coordinates of the 1st vGCM grid
                                   ! i.e., at lon(0,1) & lat(0,1)

    real(kind=r8), allocatable :: lat(:,:)    ! vGCM grid latitude
    real(kind=r8), allocatable :: lon(:,:)    ! vGCM grid longitude
    real(kind=r8), allocatable :: TH(:,:,:)   ! potential temperature
    real(kind=r8), allocatable :: QV(:,:,:)   ! water vapor
    real(kind=r8), allocatable :: QC(:,:,:)   ! cloud liquid water
    real(kind=r8), allocatable :: QI(:,:,:)   ! cloud ice water
    real(kind=r8), allocatable :: QR(:,:,:)   ! cloud rain
    real(kind=r8), allocatable :: QS(:,:,:)   ! cloud snow
    real(kind=r8), allocatable :: QG(:,:,:)   ! cloud groupel
    real(kind=r8), allocatable :: QT(:,:,:,:) ! tracer array

    real(kind=r8), allocatable :: Ucon(:,:,:) ! contravariant zonal velocity component
    real(kind=r8), allocatable :: Vcon(:,:,:) ! contravariant meridional velocity component
    real(kind=r8), allocatable :: Ucov(:,:,:) ! covariant zonal velocity component
    real(kind=r8), allocatable :: Vcov(:,:,:) ! covariant meridional velocity component
    real(kind=r8), allocatable :: W(:,:,:)    ! vertical velocity
  end type vGCM_state_t

  type, public :: vGCM_tend_t
    real(kind=r8), allocatable :: dTH(:,:)    ! potential temperature tendency
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
    ! diagnostic fields can be added here.
    real(kind=r8), allocatable :: SPREC(:)    ! surface precipitation
    real(kind=r8), allocatable :: WTH(:)      ! surface heat flux
    real(kind=r8), allocatable :: WQV(:)      ! surface moisture flux
    real(kind=r8), allocatable :: UW(:)       ! surface momentum flux (u"w")
    real(kind=r8), allocatable :: WV(:)       ! surface momentum flux (v"w")
  end type vGCM_out_t

!***********************************************************************************************
! Define a type (2): Combination of 4 channel segments (a channel)
!***********************************************************************************************

  type, public :: channel_vGCM_t
        type(vGCM_state_t), public :: vGCM_state(4)
        type(vGCM_tend_t),  public :: vGCM_tend(4)
        type(vGCM_out_t),   public :: vGCM_out(4)
  end type  channel_vGCM_t

!******************************************************************************
! Define variables with derived type
!******************************************************************************

  type(channel_vGCM_t), public, pointer :: channel_vGCM(:) => NULL()

!******************************************************************************

!******************************************************************************
! Derived MPI types
!******************************************************************************

  integer, public, protected :: vGCM_state_mpi_type
  integer, public, protected :: vGCM_tend_mpi_type

!******************************************************************************

!******************************************************************************
! Public interface subroutines
!******************************************************************************

  public :: vGCM_allocate_channel_data
  public :: vGCM_create_MPI_types

!******************************************************************************

CONTAINS

  subroutine vGCM_allocate_channel_data(channel_vGCM, nch)

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

  end subroutine vGCM_allocate_channel_data

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
    allocate(seg_tend%dQT(nVGCM_seg, nLevel, ntracer))

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

  subroutine deallocate_channel_vGCM_section (seg_state,seg_tend,seg_out)

    ! Dummy argument
    type(vGCM_state_t), intent(inout) :: seg_state  ! vGCM state variables to be allocated
    type(vGCM_tend_t),  intent(inout) :: seg_tend   ! vGCM tendency variables to be allocated
    type(vGCM_out_t),   intent(inout) :: seg_out    ! vGCM output variables to be allocated
!-----------------------------------------------------------------------
!   State Variables
!-----------------------------------------------------------------------
    deallocate(seg_state%lat)
    deallocate(seg_state%lon)

    deallocate(seg_state%TH)
    deallocate(seg_state%QV)
    deallocate(seg_state%QC)
    deallocate(seg_state%QI)
    deallocate(seg_state%QR)
    deallocate(seg_state%QS)
    deallocate(seg_state%QG)
    deallocate(seg_state%QT)

    deallocate(seg_state%Ucon)
    deallocate(seg_state%Vcon)
    deallocate(seg_state%Ucov)
    deallocate(seg_state%Vcov)
    deallocate(seg_state%W)

!-----------------------------------------------------------------------
!   Tendency Variables
!-----------------------------------------------------------------------
    deallocate(seg_tend%dTH)
    deallocate(seg_tend%dQV)
    deallocate(seg_tend%dQC)
    deallocate(seg_tend%dQI)
    deallocate(seg_tend%dQR)
    deallocate(seg_tend%dQS)
    deallocate(seg_tend%dQG)
    deallocate(seg_tend%dQT)

    deallocate(seg_tend%dU)
    deallocate(seg_tend%dV)

!-----------------------------------------------------------------------
!   Output Variables
!-----------------------------------------------------------------------
    deallocate(seg_out%SPREC)
    deallocate(seg_out%WTH)
    deallocate(seg_out%WQV)
    deallocate(seg_out%UW)
    deallocate(seg_out%WV)

  end subroutine deallocate_channel_vGCM_section

  subroutine vGCM_create_MPI_types()
    use spmd_utils,      only: masterproc, masterprocid, iam, npes, mpicom
    use spmd_utils,      only: MPI_INTEGER, MPI_REAL8, MPI_CHARACTER
    use spmd_utils,      only: MPI_ADDRESS_KIND

    integer(kind=MPI_ADDRESS_KIND)   :: offsets(20)    ! For new MPI types
    integer                          :: origtypes(20)  ! For new MPI types
    integer                          :: lengths(20)    ! For new MPI types
    integer(kind=MPI_ADDRESS_KIND)   :: extent         ! For new MPI types
    integer                          :: num_fields     ! For new MPI types
    integer                          :: size_2d        ! Size of 2D arrays
    integer                          :: size_3d        ! Size of 3D arrays
    integer                          :: ierr           ! Error code
    integer                          :: index
    integer                          :: dummy_type     ! For new MPI types
    type(vGCM_state_t)               :: dummy_state(2) ! For new MPI types
    type(vGCM_tend_t)                :: dummy_tend(2)  ! For new MPI types
    type(vGCM_out_t)                 :: dummy_out(2)   ! For interface

    ! Define the MPI types needed to send state and tendencies between the
    ! vGCM and the  CRM
    call allocate_channel_vGCM_section(dummy_state(1), dummy_tend(1), dummy_out(1))
    call allocate_channel_vGCM_section(dummy_state(2), dummy_tend(2), dummy_out(2))

    ! type vGCM_state_t
    size_2d = ((2 * nhalo_vgcm) + 1) * (nvgcm_seg + (2 * nhalo_vgcm))
    size_3d = size_2d * nlevel
    num_fields = 0
    ! Most fields are MPI_REAL8
    origtypes(:) = MPI_REAL8
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%nface,  offsets(num_fields), ierr)
    lengths(num_fields) = 1
    origtypes(num_fields) = MPI_INTEGER
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%alpha1, offsets(num_fields), ierr)
    lengths(num_fields) = 1
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%beta1,  offsets(num_fields), ierr)
    lengths(num_fields) = 1
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%lat,    offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%lon,    offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%TH,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QV,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QC,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QI,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QR,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QS,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QG,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%QT,     offsets(num_fields), ierr)
    lengths(num_fields) = size_3d * ntracer
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%Ucon,   offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%Vcon,   offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%Ucov,   offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%Vcov,   offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_state(1)%W,      offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    ! Adjust offsets and create type
    do index = num_fields, 1, -1
      offsets(index) = offsets(index) - offsets(1)
    end do
    call MPI_type_create_struct(num_fields, lengths(1:num_fields),            &
         offsets(1:num_fields), origtypes(1:num_fields), dummy_type, ierr)
    ! Adjust for padding
    call MPI_Get_address(dummy_state(1)%nface, offsets(1), ierr)
    call MPI_Get_address(dummy_state(2)%nface, offsets(2), ierr)
    extent = offsets(2) - offsets(1)
    call MPI_type_create_resized(dummy_type, 0_MPI_ADDRESS_KIND, extent, &
         vGCM_state_mpi_type, ierr)
    call MPI_type_commit(vGCM_state_mpi_type, ierr)

    ! type vGCM_tend_t
    num_fields = 0
    size_2d = nvgcm_seg * nlevel
    size_3d = size_2d * ntracer
    ! All fields are MPI_REAL8
    origtypes(:) = MPI_REAL8
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dTH, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQV, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQC, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQI, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQR, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQS, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQG, offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dQT, offsets(num_fields), ierr)
    lengths(num_fields) = size_3d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dU,  offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    num_fields = num_fields + 1
    call MPI_Get_address(dummy_tend(1)%dV,  offsets(num_fields), ierr)
    lengths(num_fields) = size_2d
    do index = num_fields, 1, -1
      offsets(index) = offsets(index) - offsets(1)
    end do
    call MPI_type_create_struct(num_fields, lengths(1:num_fields),            &
         offsets(1:num_fields), origtypes(1:num_fields), dummy_type, ierr)
    ! Adjust for padding
    call MPI_Get_address(dummy_tend(1)%dtH, offsets(1), ierr)
    call MPI_Get_address(dummy_tend(2)%dtH, offsets(2), ierr)
    extent = offsets(2) - offsets(1)
    call MPI_type_create_resized(dummy_type, 0_MPI_ADDRESS_KIND, extent, &
         vGCM_tend_mpi_type, ierr)
    call MPI_type_commit(vGCM_tend_mpi_type, ierr)

    !  Cleanup
    call deallocate_channel_vGCM_section(dummy_state(1), dummy_tend(1), dummy_out(1))
    call deallocate_channel_vGCM_section(dummy_state(2), dummy_tend(2), dummy_out(2))

  end subroutine vGCM_create_MPI_types

end module vGCM_data_types
