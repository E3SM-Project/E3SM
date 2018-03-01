module BeTR_nitrogenfluxType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_nitrogenflux_type
    real(r8), pointer :: nflx_input_litr_met_vr_col(:,:) => null() ! metabolic litter input, gC/m3/s
    real(r8), pointer :: nflx_input_litr_cel_vr_col(:,:) => null() ! cellulose litter input, gC/m3/s
    real(r8), pointer :: nflx_input_litr_lig_vr_col(:,:) => null() ! lignin litter input, gC/m3/s
    real(r8), pointer :: nflx_input_litr_cwd_vr_col(:,:) => null() ! coarse woody debris input, gC/m3/s
    real(r8), pointer :: nflx_input_litr_fwd_vr_col(:,:) => null() ! coarse woody debris input, gC/m3/s
    real(r8), pointer :: nflx_input_litr_lwd_vr_col(:,:) => null() ! coarse woody debris input, gC/m3/s

    !The only loss is through fire and no som is lost through burning
    real(r8), pointer :: nflx_output_litr_met_vr_col(:,:) => null() ! metabolic litter input, gN/m3/s
    real(r8), pointer :: nflx_output_litr_cel_vr_col(:,:) => null() ! cellulose litter input, gN/m3/s
    real(r8), pointer :: nflx_output_litr_lig_vr_col(:,:) => null() ! lignin litter input   , gN/m3/s
    real(r8), pointer :: nflx_output_litr_cwd_vr_col(:,:) => null() ! coarse woody debris input, gN/m3/s
    real(r8), pointer :: nflx_output_litr_fwd_vr_col(:,:) => null() ! coarse woody debris input, gN/m3/s
    real(r8), pointer :: nflx_output_litr_lwd_vr_col(:,:) => null() ! coarse woody debris input, gN/m3/s

    real(r8), pointer :: nflx_minn_input_nh4_vr_col(:,:)  => null() !mineral nh4 input through deposition & fertilization, gN/m3/s
    real(r8), pointer :: nflx_minn_input_no3_vr_col(:,:)  => null() !mineral no3 input through deposition & fertilization, gN/m3/s
    real(r8), pointer :: nflx_minn_nh4_fix_nomic_vr_col(:,:) => null()    !nitrogen fixation from non-microbe explicit calculation, gN/m3/s
  contains
    procedure, public  :: Init
    procedure, public  :: reset
    procedure, private :: InitAllocate
  end type betr_nitrogenflux_type

  contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
  implicit none
  class(betr_nitrogenflux_type)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  implicit none
  class(betr_nitrogenflux_type)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%nflx_input_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%nflx_input_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%nflx_input_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%nflx_input_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%nflx_input_litr_fwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%nflx_input_litr_lwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input

  allocate(this%nflx_output_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%nflx_output_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%nflx_output_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%nflx_output_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%nflx_output_litr_fwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%nflx_output_litr_lwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input

  allocate(this%nflx_minn_input_nh4_vr_col(begc:endc,lbj:ubj)) !mineral nh4 input through deposition & fertilization
  allocate(this%nflx_minn_input_no3_vr_col(begc:endc,lbj:ubj)) !mineral no3 input through deposition & fertilization
  allocate(this%nflx_minn_nh4_fix_nomic_vr_col(begc:endc,lbj:ubj))   !nh4 from fixation


  end subroutine InitAllocate


  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_nitrogenflux_type)  :: this
  real(r8), intent(in) :: value_column

  this%nflx_input_litr_met_vr_col(:,:) = value_column
  this%nflx_input_litr_cel_vr_col(:,:) = value_column
  this%nflx_input_litr_lig_vr_col(:,:)= value_column
  this%nflx_input_litr_cwd_vr_col(:,:)= value_column
  this%nflx_input_litr_fwd_vr_col(:,:)= value_column
  this%nflx_input_litr_lwd_vr_col(:,:)= value_column

  this%nflx_output_litr_met_vr_col(:,:)= value_column
  this%nflx_output_litr_cel_vr_col(:,:)= value_column
  this%nflx_output_litr_lig_vr_col(:,:)= value_column
  this%nflx_output_litr_cwd_vr_col(:,:)= value_column
  this%nflx_output_litr_fwd_vr_col(:,:)= value_column
  this%nflx_output_litr_lwd_vr_col(:,:)= value_column

  this%nflx_minn_nh4_fix_nomic_vr_col(:,:)= value_column
  this%nflx_minn_input_nh4_vr_col(:,:)= value_column
  this%nflx_minn_input_no3_vr_col(:,:)= value_column
  end subroutine reset

end module BeTR_nitrogenfluxType
