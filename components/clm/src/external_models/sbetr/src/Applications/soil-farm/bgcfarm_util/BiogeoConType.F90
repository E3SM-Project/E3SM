module BiogeoConType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
implicit none

 private
  character(len=*), private, parameter :: filename = &
       __FILE__

 type, public :: BiogeoCon_type

  !decomposition
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding

  real(r8) :: rf_l1s1_bgc
  real(r8) :: rf_l2s1_bgc
  real(r8) :: rf_l3s2_bgc
  real(r8) :: rf_s2s1_bgc
  real(r8) :: rf_s3s1_bgc
  real(r8) :: cwd_fcel_bgc
  real(r8) :: cwd_flig_bgc
  real(r8) :: lwd_fcel_bgc
  real(r8) :: lwd_flig_bgc
  real(r8) :: fwd_fcel_bgc
  real(r8) :: fwd_flig_bgc
  real(r8) :: k_decay_lit1
  real(r8) :: k_decay_lit2
  real(r8) :: k_decay_lit3
  real(r8) :: k_decay_som1
  real(r8) :: k_decay_som2
  real(r8) :: k_decay_som3
  real(r8) :: k_decay_cwd
  real(r8) :: k_decay_lwd
  real(r8) :: k_decay_fwd
  !nitrification-denitrification
  real(r8) :: nitrif_n2o_loss_frac
  real(r8) :: organic_max
  real(r8) :: rij_kro_a
  real(r8) :: rij_kro_alpha
  real(r8) :: rij_kro_beta
  real(r8) :: rij_kro_gamma
  real(r8) :: rij_kro_delta
  real(r8) :: surface_tension_water

  real(r8) :: init_cn_met
  real(r8) :: init_cn_cel
  real(r8) :: init_cn_lig
  real(r8) :: init_cn_cwd
  real(r8) :: init_cn_lwd
  real(r8) :: init_cn_fwd
  real(r8) :: init_cn_som1
  real(r8) :: init_cn_som2
  real(r8) :: init_cn_som3

  real(r8) :: init_cp_met
  real(r8) :: init_cp_cel
  real(r8) :: init_cp_lig
  real(r8) :: init_cp_cwd
  real(r8) :: init_cp_lwd
  real(r8) :: init_cp_fwd
  real(r8) :: init_cp_som1
  real(r8) :: init_cp_som2
  real(r8) :: init_cp_som3

  real(r8) :: init_cc13_met
  real(r8) :: init_cc13_cel
  real(r8) :: init_cc13_lig
  real(r8) :: init_cc13_cwd
  real(r8) :: init_cc13_lwd
  real(r8) :: init_cc13_fwd
  real(r8) :: init_cc13_som1
  real(r8) :: init_cc13_som2
  real(r8) :: init_cc13_som3

  real(r8) :: init_cc14_met
  real(r8) :: init_cc14_cel
  real(r8) :: init_cc14_lig
  real(r8) :: init_cc14_cwd
  real(r8) :: init_cc14_lwd
  real(r8) :: init_cc14_fwd
  real(r8) :: init_cc14_som1
  real(r8) :: init_cc14_som2
  real(r8) :: init_cc14_som3

  logical :: use_c13
  logical :: use_c14

  !ECA nutrient competition
  real(r8), pointer :: vmax_minp_secondary_to_occlude(:)  => null() !maximum conversion rate of secondary P into occluded P
  real(r8), pointer :: vmax_minp_soluble_to_secondary(:)  => null() !maximum conversion rate of soluble P into secondary P

  !inorganic phosphorus cycling
  real(r8) :: frac_p_sec_to_sol                !fraction of released secondary phosphorus that goes into soluble form
  real(r8), pointer :: minp_secondary_decay(:) => null() !decay rate of secondary phosphorus
 contains
   procedure, public :: Init
   procedure, private :: InitAllocate
   procedure, private :: set_defpar_default
   procedure, private :: ReadNamelist
 end type BiogeoCon_type

 type(BiogeoCon_type), public :: bgc_con_eca
contains
  !--------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, bstatus)
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  character(len=betr_namelist_buffer_size_ext) , intent(in)    :: namelist_buffer
  type(betr_status_type)                   , intent(out) :: bstatus

  call bstatus%reset()

  call this%InitAllocate()
  call this%set_defpar_default()

  !update parameter from namelist
  call this%ReadNamelist(namelist_buffer, bstatus)

  end subroutine Init
  !--------------------------------------------------------------------
  subroutine InitAllocate(this)
  use betr_varcon, only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(BiogeoCon_type), intent(inout) :: this


  allocate(this%vmax_minp_secondary_to_occlude(betr_max_soilorder))
  allocate(this%minp_secondary_decay(betr_max_soilorder))
  allocate(this%vmax_minp_soluble_to_secondary(betr_max_soilorder))

  !the following will be actually calculated from CNP bgc
  end subroutine InitAllocate
  !--------------------------------------------------------------------
  subroutine set_defpar_default(this)
  implicit none
  class(BiogeoCon_type), intent(inout) :: this

  !decomposition
  this%Q10                   = 2._r8
  this%froz_q10              = 10._r8
  this%decomp_depth_efolding = 1._r8

  !following is based on Table 15.4 in CLM4.5 tech note
  this%rf_l1s1_bgc           = 0.55_r8
  this%rf_l2s1_bgc           = 0.5_r8
  this%rf_l3s2_bgc           = 0.5_r8
  this%rf_s2s1_bgc           = 0.55_r8
  this%rf_s3s1_bgc           = 0.55_r8
  this%cwd_fcel_bgc          = 0.76_r8
  this%cwd_flig_bgc          = 0.24_r8

  !following is based on Table 15.3 in CLM4.5 tech note
  this%k_decay_lit1          = 1._r8/(0.066_r8*86400._r8*365._r8)   !1/second
  this%k_decay_lit2          = 1._r8/(0.25_r8*86400._r8*365._r8)    !1/second
  this%k_decay_lit3          = 1._r8/(0.25_r8*86400._r8*365._r8)    !1/second
  this%k_decay_som1          = 1._r8/(0.17_r8*86400._r8*365._r8)    !1/second
  this%k_decay_som2          = 1._r8/(6.1_r8*86400._r8*365._r8)     !1/second
  this%k_decay_som3          = 1._r8/(270._r8*86400._r8*365._r8)    !1/second
  this%k_decay_cwd           = 1._r8/(4.1_r8*86400._r8*365._r8)     !1/second

  this%init_cn_met  = 90._r8  !mass based
  this%init_cn_cel  = 90._r8  !mass based
  this%init_cn_lig  = 90._r8  !mass based
  this%init_cn_cwd  = 90._r8  !mass based
  this%init_cn_som1 = 8._r8   !mass based
  this%init_cn_som2 = 11._r8  !mass based
  this%init_cn_som3 = 11._r8  !mass based

  this%init_cp_met  = 1600._r8
  this%init_cp_cel  = 2000._r8
  this%init_cp_lig  = 2500._r8
  this%init_cp_cwd  = 4500._r8
  this%init_cp_som1 = 110._r8 !mass based
  this%init_cp_som2 = 320._r8 !mass based
  this%init_cp_som3 = 114._r8 !mass based

  !nitrification-denitrification
  this%nitrif_n2o_loss_frac  = 1.e-4_r8   !Arah and Vinten, 1995
  this%organic_max           = 160._r8    !organic matter content (kg/m3) where soil is assumed to act like peat
  this%rij_kro_a             = 1.5e-10_r8 ! Arah and Vinten, 1995
  this%rij_kro_alpha         = 1.26_r8    ! Arah and Vinten, 1995
  this%rij_kro_beta          = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_gamma         = 0.6_r8     ! Arah and Vinten, 1995
  this%rij_kro_delta         = 0.85_r8    ! Arah and Vinten, 1995
  this%surface_tension_water = 73.e-3_r8  ! (J/m^2), Arah and Vinten, 1995

  !ECA nutrient competition
  this%vmax_minp_secondary_to_occlude(:) = 1.e-5_r8
  this%vmax_minp_soluble_to_secondary(:) = 1.e-5_r8
  !inorganic phosphorus cycling
  this%frac_p_sec_to_sol                 = 0.2_r8 !fraction of released secondary phosphorus that goes into soluble form
  this%minp_secondary_decay(:)           = 1.e-5_r8 !decay rate of secondary phosphorus

  this%use_c13 = .false.
  this%use_c14 = .false.

  this%init_cc13_met = 0._r8
  this%init_cc13_cel = 0._r8
  this%init_cc13_lig = 0._r8
  this%init_cc13_cwd = 0._r8
  this%init_cc13_som1= 0._r8
  this%init_cc13_som2= 0._r8
  this%init_cc13_som3= 0._r8

  this%init_cc14_met = 0._r8
  this%init_cc14_cel = 0._r8
  this%init_cc14_lig = 0._r8
  this%init_cc14_cwd = 0._r8
  this%init_cc14_som1= 0._r8
  this%init_cc14_som2= 0._r8
  this%init_cc14_som3= 0._r8
  end subroutine set_defpar_default

  !--------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer, bstatus)
  !
  ! DESCRIPTION
  ! reading bgc parameters
  ! will be updated later
  use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  use betr_ctrl      , only : iulog => biulog
  use bshr_log_mod   , only : errMsg => shr_log_errMsg
  use tracer_varcon, only : use_c13_betr, use_c14_betr
  implicit none
  class(BiogeoCon_type), intent(inout) :: this
  character(len=betr_namelist_buffer_size_ext) , intent(in)    :: namelist_buffer
  type(betr_status_type), intent(out) :: bstatus


  !
  ! !LOCAL VARIABLES:
  integer                                :: nml_error
  character(len=betr_string_length_long) :: ioerror_msg
  real(r8) :: tau_decay_lit1
  real(r8) :: tau_decay_lit2
  real(r8) :: tau_decay_lit3
  real(r8) :: tau_decay_som1
  real(r8) :: tau_decay_som2
  real(r8) :: tau_decay_som3
  real(r8) :: tau_decay_cwd

  real(r8), parameter :: year_sec=86400._r8*365._r8
  namelist / soibgc_ecaparam /                  &
       tau_decay_lit1, tau_decay_lit2, tau_decay_lit3, &
       tau_decay_som1, tau_decay_som2, tau_decay_som3


  call bstatus%reset()

  !years
  tau_decay_lit1          = 0.066_r8
  tau_decay_lit2          = 0.25_r8
  tau_decay_lit3          = 0.25_r8
  tau_decay_som1          = 0.17_r8
  tau_decay_som2          = 6.1_r8
  tau_decay_som3          = 270._r8
  tau_decay_cwd           = 4.1_r8


  if ( .false. )then
     ioerror_msg=''
     read(namelist_buffer, nml=soibgc_ecaparam, iostat=nml_error, iomsg=ioerror_msg)
     if (nml_error /= 0) then
        call bstatus%set_msg(msg="ERROR reading soibgc_ecaparam namelist "//errmsg(filename, __LINE__),err=-1)
        return
     end if
  end if

  this%use_c13 = use_c13_betr
  this%use_c14 = use_c14_betr

  this%k_decay_lit1          = 1._r8/(tau_decay_lit1*year_sec)    !1/second
  this%k_decay_lit2          = 1._r8/(tau_decay_lit2*year_sec)    !1/second
  this%k_decay_lit3          = 1._r8/(tau_decay_lit3*year_sec)    !1/second
  this%k_decay_som1          = 1._r8/(tau_decay_som1*year_sec)    !1/second
  this%k_decay_som2          = 1._r8/(tau_decay_som2*year_sec)    !1/second
  this%k_decay_som3          = 1._r8/(tau_decay_som3*year_sec)    !1/second
  this%k_decay_cwd           = 1._r8/(tau_decay_cwd*year_sec)     !1/second

  end subroutine ReadNamelist
end module BiogeoConType
