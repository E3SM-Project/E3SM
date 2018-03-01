module BgcCentCnpForcType

  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: centuryeca_forc_type
    logical :: debug
    real(r8), pointer :: ystates(:)
    !input
    real(r8) :: cflx_input_litr_met   !g C/m2/s
    real(r8) :: cflx_input_litr_cel   !g C/m2/s
    real(r8) :: cflx_input_litr_lig   !g C/m2/s
    real(r8) :: cflx_input_litr_cwd   !g C/m2/s
    real(r8) :: cflx_input_litr_fwd   !g C/m2/s
    real(r8) :: cflx_input_litr_lwd   !g C/m2/s

    real(r8) :: nflx_input_litr_met   !g N/m2/s
    real(r8) :: nflx_input_litr_cel
    real(r8) :: nflx_input_litr_lig
    real(r8) :: nflx_input_litr_cwd
    real(r8) :: nflx_input_litr_fwd
    real(r8) :: nflx_input_litr_lwd

    real(r8) :: pflx_input_litr_met
    real(r8) :: pflx_input_litr_cel
    real(r8) :: pflx_input_litr_lig
    real(r8) :: pflx_input_litr_cwd
    real(r8) :: pflx_input_litr_fwd
    real(r8) :: pflx_input_litr_lwd

    !loss
    !real(r8) :: cflx_output_litr_met
    !real(r8) :: cflx_output_litr_cel
    !real(r8) :: cflx_output_litr_lig
    !real(r8) :: cflx_output_litr_cwd
    !real(r8) :: cflx_output_litr_fwd
    !real(r8) :: cflx_output_litr_lwd

    !real(r8) :: nflx_output_litr_met
    !real(r8) :: nflx_output_litr_cel
    !real(r8) :: nflx_output_litr_lig
    !real(r8) :: nflx_output_litr_cwd
    !real(r8) :: nflx_output_litr_fwd
    !real(r8) :: nflx_output_litr_lwd

    !real(r8) :: pflx_output_litr_met
    !real(r8) :: pflx_output_litr_cel
    !real(r8) :: pflx_output_litr_lig
    !real(r8) :: pflx_output_litr_cwd
    !real(r8) :: pflx_output_litr_fwd
    !real(r8) :: pflx_output_litr_lwd

    real(r8) :: cflx_input_litr_met_c13   !g C/m2/s
    real(r8) :: cflx_input_litr_cel_c13   !g C/m2/s
    real(r8) :: cflx_input_litr_lig_c13   !g C/m2/s
    real(r8) :: cflx_input_litr_cwd_c13   !g C/m2/s
    real(r8) :: cflx_input_litr_fwd_c13   !g C/m2/s
    real(r8) :: cflx_input_litr_lwd_c13   !g C/m2/s

    real(r8) :: nflx_input_litr_met_c13   !g N/m2/s
    real(r8) :: nflx_input_litr_cel_c13
    real(r8) :: nflx_input_litr_lig_c13
    real(r8) :: nflx_input_litr_cwd_c13
    real(r8) :: nflx_input_litr_fwd_c13
    real(r8) :: nflx_input_litr_lwd_c13

    real(r8) :: pflx_input_litr_met_c13
    real(r8) :: pflx_input_litr_cel_c13
    real(r8) :: pflx_input_litr_lig_c13
    real(r8) :: pflx_input_litr_cwd_c13
    real(r8) :: pflx_input_litr_fwd_c13
    real(r8) :: pflx_input_litr_lwd_c13

    !loss
    !real(r8) :: cflx_output_litr_met_c13
    !real(r8) :: cflx_output_litr_cel_c13
    !real(r8) :: cflx_output_litr_lig_c13
    !real(r8) :: cflx_output_litr_cwd_c13
    !real(r8) :: cflx_output_litr_fwd_c13
    !real(r8) :: cflx_output_litr_lwd_c13

    !real(r8) :: nflx_output_litr_met_c13
    !real(r8) :: nflx_output_litr_cel_c13
    !real(r8) :: nflx_output_litr_lig_c13
    !real(r8) :: nflx_output_litr_cwd_c13
    !real(r8) :: nflx_output_litr_fwd_c13
    !real(r8) :: nflx_output_litr_lwd_c13

    !real(r8) :: pflx_output_litr_met_c13
    !real(r8) :: pflx_output_litr_cel_c13
    !real(r8) :: pflx_output_litr_lig_c13
    !real(r8) :: pflx_output_litr_cwd_c13
    !real(r8) :: pflx_output_litr_fwd_c13
    !real(r8) :: pflx_output_litr_lwd_c13

    real(r8) :: cflx_input_litr_met_c14   !g C/m2/s
    real(r8) :: cflx_input_litr_cel_c14   !g C/m2/s
    real(r8) :: cflx_input_litr_lig_c14   !g C/m2/s
    real(r8) :: cflx_input_litr_cwd_c14   !g C/m2/s
    real(r8) :: cflx_input_litr_fwd_c14   !g C/m2/s
    real(r8) :: cflx_input_litr_lwd_c14   !g C/m2/s

    real(r8) :: nflx_input_litr_met_c14   !g N/m2/s
    real(r8) :: nflx_input_litr_cel_c14
    real(r8) :: nflx_input_litr_lig_c14
    real(r8) :: nflx_input_litr_cwd_c14
    real(r8) :: nflx_input_litr_fwd_c14
    real(r8) :: nflx_input_litr_lwd_c14

    real(r8) :: pflx_input_litr_met_c14
    real(r8) :: pflx_input_litr_cel_c14
    real(r8) :: pflx_input_litr_lig_c14
    real(r8) :: pflx_input_litr_cwd_c14
    real(r8) :: pflx_input_litr_fwd_c14
    real(r8) :: pflx_input_litr_lwd_c14

    real(r8) :: frac_loss_lit_to_fire
    real(r8) :: frac_loss_cwd_to_fire
    !loss
    !real(r8) :: cflx_output_litr_met_c14
    !real(r8) :: cflx_output_litr_cel_c14
    !real(r8) :: cflx_output_litr_lig_c14
    !real(r8) :: cflx_output_litr_cwd_c14
    !real(r8) :: cflx_output_litr_fwd_c14
    !real(r8) :: cflx_output_litr_lwd_c14

    !real(r8) :: nflx_output_litr_met_c14
    !real(r8) :: nflx_output_litr_cel_c14
    !real(r8) :: nflx_output_litr_lig_c14
    !real(r8) :: nflx_output_litr_cwd_c14
    !real(r8) :: nflx_output_litr_fwd_c14
    !real(r8) :: nflx_output_litr_lwd_c14

    !real(r8) :: pflx_output_litr_met_c14
    !real(r8) :: pflx_output_litr_cel_c14
    !real(r8) :: pflx_output_litr_lig_c14
    !real(r8) :: pflx_output_litr_cwd_c14
    !real(r8) :: pflx_output_litr_fwd_c14
    !real(r8) :: pflx_output_litr_lwd_c14

    real(r8) :: sflx_minn_input_nh4       !nh4 from deposition and fertilization
    real(r8) :: sflx_minn_nh4_fix_nomic         !nh4 from fixation
    real(r8) :: sflx_minp_input_po4       !inorganic P from deposition and fertilization
    real(r8) :: sflx_minp_weathering_po4

    real(r8) :: rt_ar              !root autotrophic respiration, mol CO2/m3/s
    real(r8) :: rt_ar_c13
    real(r8) :: rt_ar_c14
    real(r8) :: temp               !temperature
    real(r8) :: depz               !depth of the soil
    real(r8) :: dzsoi              !soil thickness
    real(r8) :: sucsat             ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
    real(r8) :: soilpsi            ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
    real(r8) :: bsw
    real(r8) :: bd                 !bulk density
    real(r8) :: pct_sand
    real(r8) :: pct_clay
    real(r8) :: h2osoi_vol
    real(r8) :: h2osoi_liq
    real(r8) :: air_vol
    real(r8) :: finundated
    real(r8) :: watsat
    real(r8) :: watfc
    real(r8) :: cellorg
    real(r8) :: pH

    real(r8) :: aren_cond_n2
    real(r8) :: aren_cond_o2
    real(r8) :: aren_cond_n2o
    real(r8) :: aren_cond_co2
    real(r8) :: aren_cond_co2_c13
    real(r8) :: aren_cond_co2_c14
    real(r8) :: aren_cond_ar
    real(r8) :: aren_cond_ch4

    real(r8) :: ch4_g2b
    real(r8) :: co2_g2b
    real(r8) :: o2_g2b
    real(r8) :: o2_w2b             !conversion parameter for o2 from aqueous to bulk conc
    real(r8) :: n2_g2b
    real(r8) :: ar_g2b
    real(r8) :: n2o_g2b

    real(r8) :: conc_atm_n2   !n2 concentration in atmosphere, mol n2/m3
    real(r8) :: conc_atm_n2o
    real(r8) :: conc_atm_o2
    real(r8) :: conc_atm_ar
    real(r8) :: conc_atm_co2
    real(r8) :: conc_atm_co2_c13
    real(r8) :: conc_atm_co2_c14
    real(r8) :: conc_atm_ch4

    real(r8),pointer :: plant_froot_nn(:)
    real(r8),pointer :: plant_froot_np(:)
    integer, pointer :: plant_vtype(:)
    integer :: plant_ntypes
    integer :: soilorder
    real(r8):: msurf_nh4
    real(r8):: msurf_minp
  contains
    procedure, public :: init
    procedure, private:: initAllocate
    procedure, public :: set_defpar
  end type centuryeca_forc_type

  public :: create_century_forc_type
contains

  function create_century_forc_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(centuryeca_forc_type), pointer :: create_century_forc_type
    class(centuryeca_forc_type), pointer :: forc

    allocate(forc)
    create_century_forc_type => forc

  end function create_century_forc_type
  !--------------------------------------------------------------------
  subroutine Init(this, centurybgc_index)
  use BgcCentCnpIndexType       , only : centurybgc_index_type
  implicit none
  class(centuryeca_forc_type) , intent(inout) :: this
  type(centurybgc_index_type) , intent(in) :: centurybgc_index

  call this%InitAllocate(centurybgc_index)

  call this%set_defpar()
  this%debug =.false.
  end subroutine Init

  !--------------------------------------------------------------------

  subroutine InitAllocate(this, centurybgc_index)
  use BgcCentCnpIndexType , only : centurybgc_index_type
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(centuryeca_forc_type) , intent(inout) :: this
  type(centurybgc_index_type) , intent(in) :: centurybgc_index


  allocate(this%ystates(centurybgc_index%nstvars))
  allocate(this%plant_froot_nn(betr_maxpatch_pft))
  allocate(this%plant_froot_np(betr_maxpatch_pft))
  allocate(this%plant_vtype(betr_maxpatch_pft)); this%plant_vtype(:) = 0
  end subroutine InitAllocate
  !--------------------------------------------------------------------
  subroutine set_defpar(this)
  !
  !DESCRIPTION
  !set default parameter values
  !this is designed for offline testing
  implicit none
  class(centuryeca_forc_type), intent(inout) :: this

  !input
  this%cflx_input_litr_met =0._r8   !g C/m2/s
  this%cflx_input_litr_cel =0._r8   !g C/m2/s
  this%cflx_input_litr_lig =0._r8   !g C/m2/s
  this%cflx_input_litr_cwd =0._r8
  this%cflx_input_litr_fwd =0._r8
  this%cflx_input_litr_lwd =0._r8

  this%nflx_input_litr_met =0._r8
  this%nflx_input_litr_cel =0._r8
  this%nflx_input_litr_lig =0._r8
  this%nflx_input_litr_cwd =0._r8
  this%nflx_input_litr_fwd =0._r8
  this%nflx_input_litr_lwd =0._r8

  this%pflx_input_litr_met =0._r8
  this%pflx_input_litr_cel =0._r8
  this%pflx_input_litr_lig =0._r8
  this%pflx_input_litr_cwd =0._r8
  this%pflx_input_litr_fwd =0._r8
  this%pflx_input_litr_lwd =0._r8

    !loss
!  this%cflx_output_litr_met =0._r8
!  this%cflx_output_litr_cel =0._r8
!  this%cflx_output_litr_lig =0._r8
!  this%cflx_output_litr_cwd =0._r8
!  this%cflx_output_litr_fwd =0._r8
!  this%cflx_output_litr_lwd =0._r8

!  this%nflx_output_litr_met =0._r8
!  this%nflx_output_litr_cel =0._r8
!  this%nflx_output_litr_lig =0._r8
!  this%nflx_output_litr_cwd =0._r8
!  this%nflx_output_litr_fwd =0._r8
!  this%nflx_output_litr_lwd =0._r8

!  this%pflx_output_litr_met =0._r8
!  this%pflx_output_litr_cel =0._r8
!  this%pflx_output_litr_lig =0._r8
!  this%pflx_output_litr_cwd =0._r8
!  this%pflx_output_litr_fwd =0._r8
!  this%pflx_output_litr_lwd =0._r8

  !input
  this%cflx_input_litr_met_c13 =0._r8   !g C/m2/s
  this%cflx_input_litr_cel_c13 =0._r8   !g C/m2/s
  this%cflx_input_litr_lig_c13 =0._r8   !g C/m2/s
  this%cflx_input_litr_cwd_c13 =0._r8
  this%cflx_input_litr_fwd_c13 =0._r8
  this%cflx_input_litr_lwd_c13 =0._r8

  this%nflx_input_litr_met_c13 =0._r8
  this%nflx_input_litr_cel_c13 =0._r8
  this%nflx_input_litr_lig_c13 =0._r8
  this%nflx_input_litr_cwd_c13 =0._r8
  this%nflx_input_litr_fwd_c13 =0._r8
  this%nflx_input_litr_lwd_c13 =0._r8

  this%pflx_input_litr_met_c13 =0._r8
  this%pflx_input_litr_cel_c13 =0._r8
  this%pflx_input_litr_lig_c13 =0._r8
  this%pflx_input_litr_cwd_c13 =0._r8
  this%pflx_input_litr_fwd_c13 =0._r8
  this%pflx_input_litr_lwd_c13 =0._r8

    !loss
!  this%cflx_output_litr_met_c13 =0._r8
!  this%cflx_output_litr_cel_c13 =0._r8
!  this%cflx_output_litr_lig_c13 =0._r8
!  this%cflx_output_litr_cwd_c13 =0._r8
!  this%cflx_output_litr_fwd_c13 =0._r8
!  this%cflx_output_litr_lwd_c13 =0._r8

!  this%nflx_output_litr_met_c13 =0._r8
!  this%nflx_output_litr_cel_c13 =0._r8
!  this%nflx_output_litr_lig_c13 =0._r8
!  this%nflx_output_litr_cwd_c13 =0._r8
!  this%nflx_output_litr_fwd_c13 =0._r8
!  this%nflx_output_litr_lwd_c13 =0._r8

!  this%pflx_output_litr_met_c13 =0._r8
!  this%pflx_output_litr_cel_c13 =0._r8
!  this%pflx_output_litr_lig_c13 =0._r8
!  this%pflx_output_litr_cwd_c13 =0._r8
!  this%pflx_output_litr_fwd_c13 =0._r8
!  this%pflx_output_litr_lwd_c13 =0._r8

  !input
  this%cflx_input_litr_met_c14 =0._r8   !g C/m2/s
  this%cflx_input_litr_cel_c14 =0._r8   !g C/m2/s
  this%cflx_input_litr_lig_c14 =0._r8   !g C/m2/s
  this%cflx_input_litr_cwd_c14 =0._r8
  this%cflx_input_litr_fwd_c14 =0._r8
  this%cflx_input_litr_lwd_c14 =0._r8

  this%nflx_input_litr_met_c14 =0._r8
  this%nflx_input_litr_cel_c14 =0._r8
  this%nflx_input_litr_lig_c14 =0._r8
  this%nflx_input_litr_cwd_c14 =0._r8
  this%nflx_input_litr_fwd_c14 =0._r8
  this%nflx_input_litr_lwd_c14 =0._r8

  this%pflx_input_litr_met_c14 =0._r8
  this%pflx_input_litr_cel_c14 =0._r8
  this%pflx_input_litr_lig_c14 =0._r8
  this%pflx_input_litr_cwd_c14 =0._r8
  this%pflx_input_litr_fwd_c14 =0._r8
  this%pflx_input_litr_lwd_c14 =0._r8

    !loss
!  this%cflx_output_litr_met_c14 =0._r8
!  this%cflx_output_litr_cel_c14 =0._r8
!  this%cflx_output_litr_lig_c14 =0._r8
!  this%cflx_output_litr_cwd_c14 =0._r8
!  this%cflx_output_litr_fwd_c14 =0._r8
!  this%cflx_output_litr_lwd_c14 =0._r8

!  this%nflx_output_litr_met_c14 =0._r8
!  this%nflx_output_litr_cel_c14 =0._r8
!  this%nflx_output_litr_lig_c14 =0._r8
!  this%nflx_output_litr_cwd_c14 =0._r8
!  this%nflx_output_litr_fwd_c14 =0._r8
!  this%nflx_output_litr_lwd_c14 =0._r8

!  this%pflx_output_litr_met_c14 =0._r8
!  this%pflx_output_litr_cel_c14 =0._r8
!  this%pflx_output_litr_lig_c14 =0._r8
!  this%pflx_output_litr_cwd_c14 =0._r8
!  this%pflx_output_litr_fwd_c14 =0._r8
!  this%pflx_output_litr_lwd_c14 =0._r8

  this%sflx_minn_input_nh4      =0._r8       !nh4 from deposition and fertilization
  this%sflx_minn_nh4_fix_nomic        =0._r8      !nh4 from fixation
  this%sflx_minp_input_po4      =0._r8      !inorganic P from deposition and fertilization
  this%sflx_minp_weathering_po4 =0._r8

  this%rt_ar_c13  =0._r8             !root autotrophic respiration, mol CO2/m3/s
  this%rt_ar_c14  =0._r8             !root autotrophic respiration, mol CO2/m3/s
  this%rt_ar  =0._r8             !root autotrophic respiration, mol CO2/m3/s
  this%temp   =0._r8            !temperature
  this%depz   =0._r8             !depth of the soil
  this%dzsoi  =0._r8             !soil thickness
  this%sucsat =0._r8             ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
  this%soilpsi=0._r8             ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
  this%bsw    =0._r8
  this%bd     =0._r8             !bulk density
  this%pct_sand=0._r8
  this%h2osoi_vol=0._r8
  this%h2osoi_liq=0._r8
  this%air_vol=0._r8
  this%finundated=0._r8
  this%watsat=0._r8
  this%watfc=0._r8
  this%cellorg=0._r8
  this%pH=0._r8

  this%aren_cond_n2=0._r8
  this%aren_cond_o2=0._r8
  this%aren_cond_n2o=0._r8
  this%aren_cond_co2=0._r8
  this%aren_cond_co2_c13=0._r8
  this%aren_cond_co2_c14=0._r8
  this%aren_cond_ar=0._r8
  this%aren_cond_ch4=0._r8

  this%ch4_g2b=0._r8
  this%co2_g2b=0._r8
  this%o2_g2b =0._r8
  this%o2_w2b =0._r8       !conversion parameter for o2 from aqueous to bulk conc
  this%n2_g2b =0._r8
  this%ar_g2b =0._r8
  this%n2o_g2b=0._r8

  this%conc_atm_n2 =0._r8   !n2 concentration in atmosphere, mol n2/m3
  this%conc_atm_n2o=0._r8
  this%conc_atm_o2=0._r8
  this%conc_atm_ar=0._r8
  this%conc_atm_co2=0._r8
  this%conc_atm_co2_c13=0._r8
  this%conc_atm_co2_c14=0._r8
  this%conc_atm_ch4=0._r8

  this%plant_froot_nn(:)=0._r8
  this%plant_froot_np(:)=0._r8
  this%plant_vtype(:)=0
  this%plant_ntypes=0
  this%soilorder=1
  this%msurf_nh4=0._r8
  this%msurf_minp=0._r8

  this%ystates(:) = 0._r8
  end subroutine set_defpar
end module BgcCentCnpForcType
