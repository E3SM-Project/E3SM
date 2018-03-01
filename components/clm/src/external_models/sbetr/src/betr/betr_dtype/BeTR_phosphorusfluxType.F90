module BeTR_phosphorusfluxType

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_phosphorusflux_type

    real(r8), pointer :: pflx_input_litr_met_vr_col(:,:) => null() ! metabolic litter input, gP/m3/s
    real(r8), pointer :: pflx_input_litr_cel_vr_col(:,:) => null() ! cellulose litter input, gP/m3/s
    real(r8), pointer :: pflx_input_litr_lig_vr_col(:,:) => null() ! lignin litter input, gP/m3/s
    real(r8), pointer :: pflx_input_litr_cwd_vr_col(:,:) => null() ! coarse woody debris input, gP/m3/s
    real(r8), pointer :: pflx_input_litr_fwd_vr_col(:,:) => null() ! coarse woody debris input, gP/m3/s
    real(r8), pointer :: pflx_input_litr_lwd_vr_col(:,:) => null() ! coarse woody debris input, gP/m3/s

    !The only loss is through fire and no som is lost through burning
    real(r8), pointer :: pflx_output_litr_met_vr_col(:,:) => null() ! metabolic litter input, gP/m3/s
    real(r8), pointer :: pflx_output_litr_cel_vr_col(:,:) => null() ! cellulose litter input, gP/m3/s
    real(r8), pointer :: pflx_output_litr_lig_vr_col(:,:) => null() ! lignin litter input, gP/m3/s
    real(r8), pointer :: pflx_output_litr_cwd_vr_col(:,:) => null() ! coarse woody debris input, gP/m3/s
    real(r8), pointer :: pflx_output_litr_fwd_vr_col(:,:) => null() ! coarse woody debris input, gP/m3/s
    real(r8), pointer :: pflx_output_litr_lwd_vr_col(:,:) => null() ! coarse woody debris input, gP/m3/s

    real(r8), pointer :: pflx_minp_input_po4_vr_col(:,:)  => null() !mineral phosphorus input through deposition & fertilization, gN/m3/s
    real(r8), pointer :: pflx_minp_weathering_po4_vr_col(:,:)  => null() !mineral phosphorus input through weathering, gN/m3/s

  contains
    procedure, public  :: Init
    procedure, public  :: reset
    procedure, private :: InitAllocate

  end type betr_phosphorusflux_type

  contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
  implicit none
  class(betr_phosphorusflux_type)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  implicit none
  class(betr_phosphorusflux_type)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%pflx_input_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%pflx_input_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%pflx_input_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%pflx_input_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%pflx_input_litr_fwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%pflx_input_litr_lwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input

  allocate(this%pflx_output_litr_met_vr_col(begc:endc,lbj:ubj))
  allocate(this%pflx_output_litr_cel_vr_col(begc:endc,lbj:ubj)) ! cellulose litter input
  allocate(this%pflx_output_litr_lig_vr_col(begc:endc,lbj:ubj)) ! lignin litter input
  allocate(this%pflx_output_litr_cwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%pflx_output_litr_fwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input
  allocate(this%pflx_output_litr_lwd_vr_col(begc:endc,lbj:ubj)) ! coarse woody debries input


  allocate(this%pflx_minp_input_po4_vr_col(begc:endc,lbj:ubj)) !mineral phosphorus input through weathering, deposition & fertilization

  allocate(this%pflx_minp_weathering_po4_vr_col(begc:endc,lbj:ubj)) !p from weathering


  end subroutine InitAllocate


  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_phosphorusflux_type)  :: this
  real(r8), intent(in) :: value_column

  this%pflx_input_litr_met_vr_col(:,:) = value_column
  this%pflx_input_litr_cel_vr_col(:,:) = value_column
  this%pflx_input_litr_lig_vr_col(:,:)= value_column
  this%pflx_input_litr_cwd_vr_col(:,:)= value_column
  this%pflx_input_litr_fwd_vr_col(:,:)= value_column
  this%pflx_input_litr_lwd_vr_col(:,:)= value_column

  this%pflx_output_litr_met_vr_col(:,:)= value_column
  this%pflx_output_litr_cel_vr_col(:,:)= value_column
  this%pflx_output_litr_lig_vr_col(:,:)= value_column
  this%pflx_output_litr_cwd_vr_col(:,:)= value_column
  this%pflx_output_litr_fwd_vr_col(:,:)= value_column
  this%pflx_output_litr_lwd_vr_col(:,:)= value_column


  this%pflx_minp_input_po4_vr_col(:,:)= value_column

  this%pflx_minp_weathering_po4_vr_col(:,:) = value_column
  end subroutine reset

end module BeTR_phosphorusfluxType
