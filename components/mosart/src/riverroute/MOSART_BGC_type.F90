
module MOSART_BGC_type

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: MOSART_BGC_type
!
! !DESCRIPTION:
! Module containing data structure for riverine biogeochemistry
!
!! !REVISION HISTORY:
! Hongyi Li: Created 03/2015
!USES:
  use mct_mod
  use RunoffMod     , only : Tctl, TUnit, rtmCTL
  use RtmSpmd       , only : masterproc
  use RtmVar        , only : iulog
  use RtmIO         , only : pio_subsystem, ncd_pio_openfile, ncd_pio_closefile
  use rof_cpl_indices, only : nt_rtm
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : denh2o => SHR_CONST_RHOFW, denice => SHR_CONST_RHOICE, grav => SHR_CONST_G, SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_flush, shr_sys_abort
  use MOSART_RES_type  , only : Tres_ctl, Tres_para, Tres
  use netcdf
  use pio


! !PUBLIC TYPES:
  implicit none
  real(r8), parameter :: TINYVALUE_s = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits

  private

! control information for subbasin-based representation

  ! sediment parameters
  type Tparameter_sedi
     character(len=200) :: paraFile             ! the path of the parameter files
     character(len=100) :: inputPath            ! the path of the input files

     real(r8), pointer :: d50(:)                ! median bed-material sediment particle size [m]
     real(r8), pointer :: u_f(:)                ! formative shear bed shear velocity [m/s]
  end type Tparameter_sedi 

    ! sediment storage and flux variables
    ! note most of the storage/flux items are included in runoff_flow structure in the RunoffMod.F90
    ! here we just add some unique items for sediment such as shear stress etc.
    !public :: TstatusFlux_sedi
    type TstatusFlux_sedi

        ! subnetwork channel
        real(r8), pointer :: Ssal_t(:)      ! sand-sediment storage in active layer, [kg]
        real(r8), pointer :: Smal_t(:)      ! mudd-sediment storage in active layer, [kg]
        real(r8), pointer :: ers_t(:)       ! sand-sediment erosion, [kg/s]
        real(r8), pointer :: ersal_t(:)     ! sand-sediment erosion @ active layer, [kg/s]
        real(r8), pointer :: ersb_t(:)      ! sand-sediment erosion @ bank, [kg/s]
        real(r8), pointer :: ses_t(:)       ! sand-sediment settled to bed zone, [kg/s]
        real(r8), pointer :: seout_t(:)     ! sand-sediment discharge out of local channel, [kg/s]
        real(r8), pointer :: erm_t(:)       ! mud-sediment erosion, [kg/s]
        real(r8), pointer :: ermal_t(:)     ! mud-sediment erosion @ active layer, [kg/s]
        real(r8), pointer :: ermb_t(:)      ! mud-sediment erosion @ bank, [kg/s]
        real(r8), pointer :: sem_t(:)       ! mud-sediment settled to bed zone, [kg/s]
        real(r8), pointer :: meout_t(:)     ! mud-sediment discharge out of local channel, [kg/s]
        real(r8), pointer :: erb_t(:)       ! total-sediment erosion @ bank, [kg/s]

        ! main channel
        real(r8), pointer :: Ssal_r(:)      ! sand-sediment storage in active layer, [kg]
        real(r8), pointer :: Smal_r(:)      ! mudd-sediment storage in active layer, [kg]
        real(r8), pointer :: ers_r(:)       ! sand-sediment erosion, [kg/s]
        real(r8), pointer :: ersal_r(:)     ! sand-sediment erosion @ active layer, [kg/s]
        real(r8), pointer :: ersb_r(:)      ! sand-sediment erosion @ bank, [kg/s]
        real(r8), pointer :: ses_r(:)       ! sand-sediment settled to bed zone, [kg/s]
        real(r8), pointer :: seout_r(:)     ! sand-sediment discharge out of local channel, [kg/s]
        real(r8), pointer :: erm_r(:)       ! mud-sediment erosion, [kg/s]
        real(r8), pointer :: ermal_r(:)     ! mud-sediment erosion @ active layer, [kg/s]
        real(r8), pointer :: ermb_r(:)      ! mud-sediment erosion @ bank, [kg/s]
        real(r8), pointer :: sem_r(:)       ! mud-sediment settled to bed zone, [kg/s]
        real(r8), pointer :: meout_r(:)     ! mud-sediment discharge out of local channel, [kg/s]
        real(r8), pointer :: erb_r(:)       ! total-sediment erosion @ bank, [kg/s]

    end type TstatusFlux_sedi

    type (TstatusFlux_sedi)     ,   public :: TSedi
    type (Tparameter_sedi)      ,   public :: TSedi_para


  public :: MOSART_sediment_init

  contains    

  subroutine MOSART_sediment_init(begr, endr, numr)

! !DESCRIPTION:
! initialize MOSART-sediment variables
! 
! !USES:
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: ier, ios
    character(len=1000) :: fname
    type(file_desc_t)  :: ncid       ! pio file desc
    type(var_desc_t)   :: vardesc    ! pio variable desc 
    type(io_desc_t)    :: iodesc_dbl ! pio io desc
    type(io_desc_t)    :: iodesc_int ! pio io desc
    integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
    integer :: dids(2)               ! variable dimension ids 
    integer :: dsizes(2)             ! variable dimension lengths
    integer :: begr, endr, numr, iunit, nn, n, cnt, nr, nt
    integer :: numDT_r, numDT_t
    integer :: lsize, gsize
    integer :: igrow, igcol, iwgt
    integer :: tcnt
    character(len=16384) :: rList             ! list of fields for SM multiply
    character(len=*),parameter :: FORMI = '(2A,2i10)'
    character(len=*),parameter :: FORMR = '(2A,2g15.7)'
    character(len=*),parameter :: subname='(MOSART_sediment_para)'

    allocate(TSedi%Ssal_t(begr:endr),       &
             TSedi%Smal_t(begr:endr),       &
             TSedi%ers_t(begr:endr),       &
             TSedi%ersal_t(begr:endr),     &
             TSedi%ersb_t(begr:endr),      &
             TSedi%ses_t(begr:endr),       &
             TSedi%erm_t(begr:endr),       &
             TSedi%ermal_t(begr:endr),     &
             TSedi%ermb_t(begr:endr),      &
             TSedi%sem_t(begr:endr),       &
             TSedi%Ssal_r(begr:endr),       &
             TSedi%Smal_r(begr:endr),       &
             TSedi%ers_r(begr:endr),       &
             TSedi%ersal_r(begr:endr),     &
             TSedi%ersb_r(begr:endr),      &
             TSedi%ses_r(begr:endr),       &
             TSedi%erm_r(begr:endr),       &
             TSedi%ermal_r(begr:endr),     &
             TSedi%ermb_r(begr:endr),      &
             TSedi%sem_r(begr:endr),       &
             TSedi_para%u_f(begr:endr),    &
             !TSedi_para%d50(begr:endr),    &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of sediment local flux/status arrays'
       call shr_sys_abort()
    end if

    TSedi%Ssal_t = 0._r8
    TSedi%Smal_t = 0._r8
    TSedi%ers_t = 0._r8
    TSedi%ersal_t = 0._r8
    TSedi%ersb_t = 0._r8
    TSedi%ses_t = 0._r8
    TSedi%erm_t = 0._r8
    TSedi%ermal_t = 0._r8
    TSedi%ermb_t = 0._r8
    TSedi%sem_t = 0._r8
    TSedi%Ssal_r = 0._r8
    TSedi%Smal_r = 0._r8
    TSedi%ers_r = 0._r8
    TSedi%ersal_r = 0._r8
    TSedi%ersb_r = 0._r8
    TSedi%ses_r = 0._r8
    TSedi%erm_r = 0._r8
    TSedi%ermal_r = 0._r8
    TSedi%ermb_r = 0._r8
    TSedi%sem_r = 0._r8
    !TSedi_para%d50 = 0._r8
    TSedi_para%u_f = 0._r8
    do iunit = begr, endr
        TSedi_para%u_f(iunit) = CRUF(TUnit%rdepth(iunit), TUnit%rslp(iunit))
    end do

  end subroutine MOSART_sediment_init  

!
!
!EOP
!-----------------------------------------------------------------------
    function CRUF(h_,slope_) result(uf_)
      ! !DESCRIPTION: calculate formative shear bed shear velocity following Lamb and Venditti (2016)
      implicit none
      real(r8), intent(in) :: h_,slope_   !
      real(r8)             :: uf_       ! formative shear bed shear velocity, [m/s]

      real(r8) :: vs   !
      real(r8) :: f_, gamma_ !
      real(r8) :: D90_to_D50  !

      f_ = 1.5_r8
      gamma_ = 0.9_r8
      D90_to_D50 = 3.0_r8

      vs = sqrt(grav*h_*slope_)
      uf_ = vs*sqrt(f_*D90_to_D50**(1-gamma_))

      if(abs(uf_) < TINYVALUE_s) then
          uf_ = 0._r8
      endif
      return
    end function CRUF


end module MOSART_BGC_type 