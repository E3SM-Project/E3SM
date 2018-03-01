module BeTRTracerType
#include "bshr_assert.h"
  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! data type to configure betr simulations
  !
  ! !USES:
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use betr_constants  , only : betr_var_name_length
  use betr_ctrl       , only : do_betr_output,iulog => biulog
  use tracer_varcon   , only : sorp_isotherm_linear,  sorp_isotherm_langmuir
  !
  implicit none
  private
  character(len=*), parameter :: mod_filename = &
       __FILE__

  !----------------------------------------------------
  !betr tracer setup structure
  !----------------------------------------------------
  type, public :: BeTRtracer_type
   character(len=255) :: betr_simname                            ! name of the simulation
   integer :: nmem_max                                           ! maximum number of members in a transport group
   integer :: ntracers                                           ! total number of tracers, gas/aqueous tracers + solid tracers that undergo active mineral protection
   integer :: ngwmobile_tracers                                  ! total number of tracers potentially undergoing gas/aqueous movement
   integer :: nvolatile_tracers                                  ! number of volatile_tracers
   integer :: nsolid_equil_tracers                               ! number of tracers that undergo equilibrium adsorption in soil could include adsorbed doc, nh4(+)
   integer :: nsolid_passive_tracers                             ! number of tracers that undergo active mineral interaction
   integer :: nfrozen_tracers

   integer :: ntracer_groups                                     !
   integer :: ngwmobile_tracer_groups                            ! total number of groups for mobile tracers
   integer :: nvolatile_tracer_groups                            ! sub group within gwmobile group
   integer :: nsolid_equil_tracer_groups                         ! sub group in solid group
   integer :: nsolid_passive_tracer_groups                       ! sub group in solid group

   integer :: nh2o_tracers                                       ! number of h2o tracers, this will be used to compute vapor gradient and thermal gradient driven isotopic flow
   integer :: id_trc_n2, id_trc_beg_n2, id_trc_end_n2                           ! tag for n2
   integer :: id_trc_o2, id_trc_beg_o2, id_trc_end_o2                           ! tag for co2
   integer :: id_trc_ar, id_trc_beg_ar, id_trc_end_ar                           ! tag for ar
   integer :: id_trc_co2x, id_trc_beg_co2x, id_trc_end_co2x                       ! tag for co2 and its related species, co2x(CO2, H2CO3, HCO3(-), CO3(2-)),
   integer :: id_trc_ch4, id_trc_beg_ch4, id_trc_end_ch4                             ! tag for methane
   integer :: id_trc_Bm, id_trc_beg_Bm, id_trc_end_Bm
   integer :: id_trc_no, id_trc_beg_no, id_trc_end_no                               ! tag for no
   integer :: id_trc_n2o, id_trc_beg_n2o, id_trc_end_n2o                             ! tag for n2o
   integer :: id_trc_air_co2x                                    ! tag for atmospheric co2
   integer :: id_trc_arrt_co2x                                   ! tag for autotrophic co2
   integer :: id_trc_hrsoi_co2x                                  ! tag for heterotrophic co2
   integer :: id_trc_nh3x, id_trc_beg_nh3x, id_trc_end_nh3x                       ! tag for nh3 and its related species, nh3x(NH3, NH4OH,NH4(+))
   integer :: id_trc_no3x, id_trc_beg_no3x, id_trc_end_no3x                       ! tag for no3 and its related species, no3x(HNO3,NO3(-))
   integer :: id_trc_no2x, id_trc_beg_no2x, id_trc_end_no2x                       ! tag for no2 and its related species, no2x(HNO2,NO2(-))
   integer :: id_trc_dom, id_trc_beg_dom, id_trc_end_dom                         ! tag for generic dissolved organic matter
   integer :: id_trc_doc, id_trc_beg_doc, id_trc_end_doc                         ! tag for generic dissolved organic carbon, used for testing single carbon pool model
   integer :: id_trc_p_sol, id_trc_beg_p_sol, id_trc_end_p_sol                   ! tag for soluble inorganic P, this includes P in equilibrium adsorption

   integer :: id_trc_blk_h2o, id_trc_beg_blk_h2o, id_trc_end_blk_h2o             ! tag for bulk water, including all water and its isotopologues.
   integer :: id_trc_o18_h2o, id_trc_beg_o18_h2o, id_trc_end_o18_h2o                ! tag for H2O(18)
   integer :: id_trc_o17_h2o, id_trc_beg_o17_h2o, id_trc_end_o17_h2o                 ! tag for H2O(17)
   integer :: id_trc_o18_h2o_ice                                 ! tag for H2O(18) in ice
   integer :: id_trc_d_h2o, id_trc_beg_d_h2o, id_trc_end_d_h2o   ! tag for DHO
   integer :: id_trc_d_h2o_ice                                   ! tag for DHO in ice
   integer :: id_trc_c13_co2x, id_trc_beg_c13_co2x, id_trc_end_c13_co2x               ! tag for C(13)O2 and its related species
   integer :: id_trc_c14_co2x, id_trc_beg_c14_co2x, id_trc_end_c14_co2x               ! tag for C(14)O2 and its related species
   integer :: id_trc_o18_co2x, id_trc_beg_o18_co2x, id_trc_end_o18_co2x               ! tag for O(18)CO and its related species
   integer :: id_trc_o17_co2x, id_trc_beg_o17_co2x, id_trc_end_o17_co2x               ! tag for O(17)CO and its related species
   integer :: id_trc_litr, id_trc_beg_litr, id_trc_end_litr
   integer :: id_trc_wood, id_trc_beg_wood, id_trc_end_wood
   integer :: id_trc_som, id_trc_beg_som, id_trc_end_som
   integer :: id_trc_minp, id_trc_beg_minp, id_trc_end_minp

   integer :: id_trc_o18_o2, id_trc_beg_o18_o2, id_trc_end_o18_o2                                      ! tag for O(18)O and its related species
   integer :: id_trc_o17_o2, id_trc_beg_o17_o2, id_trc_end_o17_o2                                      ! tag for O(17)O and its related species
   integer, pointer :: id_trc_h2o_tags(:)                        !tagged h2o tracers

   logical, pointer :: is_volatile(:)    => null()                       !flag for volatile species,  true/false, (yes/no)
   logical, pointer :: is_diffusive(:)  => null()
   logical, pointer :: is_adsorb(:)     => null()                         !flag for adsorbable species, true/false (year/no), in equilibrium with aqueous phase and/or gaseous phase
   logical, pointer :: is_advective(:)  => null()                         !flag for advective species, some species, like non-dissolved som does not undergo advection, rather bioturbation is the major mechanism for vertical transport
   logical, pointer :: is_mobile(:)     => null()                         !flag indicating whether the tracer is mobile or inert, when it is innert, do not move it around
   logical, pointer :: is_h2o(:)        => null()                         !flag for water isotope
   logical, pointer :: is_co2tag(:)     => null()                         !tagged co2 tracer?
   logical, pointer :: is_dom(:)        => null()                         !true if it is a dom tracer, place holder for rtm bgc
   logical, pointer :: is_isotope(:)   => null()
   logical, pointer :: is_frozen(:)    => null()                          !true if the tracer could be frozen
   integer, pointer :: refisoid(:)     => null()                          !reference tracer for isotope calculation, this is setup only for non-h2o isotope now
   integer, pointer :: adsorbid(:)     => null()                          !which tracer is adsorbed
   integer, pointer :: volatileid(:)  => null()
   integer, pointer :: h2oid(:)    => null()
   integer, pointer :: adsorbgroupid(:)  => null()
   integer, pointer :: adsorb_isotherm(:) => null()
   integer, pointer :: volatilegroupid(:)  => null()                      !
   integer, pointer :: groupid(:) => null()
   integer, pointer :: frozenid(:) => null()

   logical :: is_tagged_h2o =.false.                             !no tagged h2o run by default
   real(r8),pointer :: tracer_solid_passive_diffus_scal_group(:) => null() !reference diffusivity for solid phase tracer, for modeling turbation
   real(r8),pointer :: tracer_solid_passive_diffus_thc_group(:) => null() !threshold diffusivity for solid phase tracer, for modeling turbation

   integer, pointer :: solid_passive_tracer_groupid(:,:)  => null()
   integer, pointer :: tracer_group_memid(:,:)            => null()       !grp, gmem
   character(len=betr_var_name_length),pointer :: tracernames(:)    => null()               !array with tracer names
   character(len=betr_var_name_length),pointer :: tracerfamilyname(:)=> null()
   character(len=betr_var_name_length),pointer :: units(:) => null()
   real(r8),pointer :: gram_mole_wt(:)      => null()                     !molecular weight of the master species, [g/mol]
   real(r8),pointer :: vtrans_scal(:)       => null()                     !scaling factor for plant tracer uptake through transpiration, for non-water neutral aqueous tracers
   logical :: debug
  contains
     procedure, public  :: Init
     procedure, public  :: init_scalars
     procedure, public  :: set_tracer
     procedure, private :: InitAllocate
     procedure, public  :: is_solidtransport
     procedure, public  :: add_tracer_group
  end type BeTRtracer_type



  contains


  subroutine Init(this)

    implicit none
    class(BeTRtracer_type), intent(inout) :: this

    this%ntracer_groups = this%nsolid_passive_tracer_groups + this%ngwmobile_tracer_groups
    !write(iulog,*)'BeTR: total tracer goups',this%ntracer_groups
    !write(iulog,*)'BeTR: total solid passive tracer groups',this%nsolid_passive_tracer_groups
    !write(iulog,*)'BeTR: total mobile tracer groups', this%ngwmobile_tracer_groups
    call this%InitAllocate()
    this%debug = .false.
  end subroutine Init
!--------------------------------------------------------------------------------
  subroutine init_scalars(this)

  ! !DESCRIPTION:
  ! initilaize scalar variables within the type

  implicit none
  class(BeTRtracer_type), intent(inout) :: this

  this%ntracers                     = 0      ! total number of tracers, gas/aqueous tracers + solid tracers that undergo active mineral protection
  this%ngwmobile_tracers            = 0      ! total number of tracers undergoing gas/aqueous movement
  this%nvolatile_tracers            = 0      ! number of volatile_tracers
  this%nsolid_equil_tracers         = 0      ! number of tracers that undergo equilibrium adsorption in soil could include adsorbed doc, nh4(+)
  this%nsolid_passive_tracers       = 0      ! number of tracers that undergo active mineral protection

  this%ntracer_groups               = 0
  this%ngwmobile_tracer_groups      = 0
  this%nvolatile_tracer_groups      = 0
  this%nsolid_equil_tracer_groups   = 0
  this%nsolid_passive_tracer_groups = 0
  this%nfrozen_tracers              = 0
  this%nh2o_tracers                 = 0      ! number of h2o tracers, this will be used to compute vapor gradient and thermal gradient driven isotopic flow

  this%id_trc_ch4          = 0; this%id_trc_beg_ch4 = 0;this%id_trc_end_ch4 = 0     ! tag for methane
  this%id_trc_o2           = 0; this%id_trc_beg_o2 = 0; this%id_trc_end_o2 = 0      ! tag for co2
  this%id_trc_n2           = 0; this%id_trc_beg_n2 = 0; this%id_trc_end_n2 = 0      ! tag for n2
  this%id_trc_no            = 0; this%id_trc_beg_no = 0; this%id_trc_end_no = 0      ! tag for no
  this%id_trc_n2o                   = 0; this%id_trc_beg_n2o= 0; this%id_trc_end_n2o= 0      ! tag for n2o
  this%id_trc_ar                    = 0; this%id_trc_beg_ar = 0      ! tag for ar
  this%id_trc_air_co2x              = 0      ! tag for atmospheric co2
  this%id_trc_arrt_co2x             = 0      ! tag for autotrophic co2
  this%id_trc_hrsoi_co2x            = 0      ! tag for heterotrophic co2

  this%id_trc_co2x                  = 0; this%id_trc_beg_co2x = 0; this%id_trc_end_co2x = 0      ! tag for co2 and its related species, co2x(CO2, H2CO3, HCO3(-), CO3(2-)),
  this%id_trc_nh3x                  = 0; this%id_trc_beg_nh3x = 0; this%id_trc_end_nh3x = 0      ! tag for nh3 and its related species, nh3x(NH3, NH4OH,NH4(+))
  this%id_trc_no3x                  = 0; this%id_trc_beg_no3x = 0; this%id_trc_end_no3x = 0      ! tag for no3 and its related species, no3x(HNO3,NO3(-))
  this%id_trc_no2x                  = 0; this%id_trc_beg_no2x = 0; this%id_trc_end_no2x = 0       ! tag for no2 and its related species, no2x(HNO2,NO2(-))
  this%id_trc_dom                   = 0; this%id_trc_beg_dom = 0;  this%id_trc_end_dom = 0      ! tag for generic dissolved organic matter
  this%id_trc_doc                   = 0; this%id_trc_beg_doc = 0; this%id_trc_end_doc = 0      ! tag for generic dissolved organic matter
  this%id_trc_p_sol                 = 0; this%id_trc_beg_p_sol = 0; this%id_trc_end_p_sol = 0
  this%id_trc_o18_h2o               = 0; this%id_trc_beg_o18_h2o=0; this%id_trc_end_o18_h2o=0      ! tag for H2O(18)
  this%id_trc_o17_h2o               = 0; this%id_trc_beg_o17_h2o=0; this%id_trc_end_o17_h2o=0      ! tag for H2O(17)
  this%id_trc_d_h2o                 = 0; this%id_trc_beg_d_h2o=0; this%id_trc_end_d_h2o=0      ! tag for DHO
  this%id_trc_c13_co2x              = 0; this%id_trc_beg_c13_co2x = 0; this%id_trc_end_c13_co2x = 0      ! tag for C(13)O2 and its related species
  this%id_trc_c14_co2x              = 0; this%id_trc_beg_c14_co2x = 0; this%id_trc_end_c14_co2x = 0      ! tag for C(14)O2 and its related species
  this%id_trc_o18_co2x              = 0; this%id_trc_beg_o18_co2x = 0; this%id_trc_end_o18_co2x = 0       ! tag for O(18)CO and its related species
  this%id_trc_o17_co2x              = 0; this%id_trc_beg_o17_co2x = 0; this%id_trc_end_o17_co2x = 0      ! tag for O(17)CO and its related species
  this%id_trc_o18_h2o_ice           = 0      ! tag for H2O(18) in ice
  this%id_trc_d_h2o_ice             = 0      ! tag for HDO in ice
  this%id_trc_o18_o2                = 0; this%id_trc_beg_o18_o2 = 0; this%id_trc_end_o18_o2 = 0      ! tag for O(18)O and its related species
  this%id_trc_o17_o2                = 0; this%id_trc_beg_o17_o2 = 0; this%id_trc_end_o17_o2 =0      ! tag for O(17)O and its related species

  this%id_trc_litr = 0; this%id_trc_beg_litr=0; this%id_trc_end_litr=0
  this%id_trc_wood = 0; this%id_trc_beg_wood=0; this%id_trc_end_wood=0
  this%id_trc_som = 0; this%id_trc_beg_som = 0; this%id_trc_end_som = 0
  this%id_trc_blk_h2o=0; this%id_trc_beg_blk_h2o=0; this%id_trc_end_blk_h2o=0
  this%id_trc_Bm = 0; this%id_trc_beg_Bm = 0; this%id_trc_end_Bm = 0

  this%betr_simname                 = ''
  end subroutine init_scalars


!--------------------------------------------------------------------------------
  subroutine InitAllocate(this)

  ! !DESCRIPTION:
  ! allocate memories for vectors

  implicit none
  class(BeTRtracer_type), intent(inout) :: this
  integer, parameter :: nanid=-1

  allocate(this%is_volatile        (this%ntracers));    this%is_volatile(:)     = .false.
  allocate(this%is_adsorb          (this%ntracers));    this%is_adsorb(:)       = .false.
  allocate(this%is_advective       (this%ntracers));             this%is_advective(:)    = .false.
  allocate(this%is_diffusive       (this%ntracers));             this%is_diffusive(:)    = .true.
  allocate(this%is_mobile          (this%ntracers));             this%is_mobile(:)       = .false.
  allocate(this%is_h2o             (this%ntracers));    this%is_h2o(:)          = .false.
  allocate(this%is_co2tag          (this%ntracers));    this%is_co2tag(:)       = .false.
  allocate(this%is_dom             (this%ntracers));    this%is_dom(:)          = .false.
  allocate(this%is_isotope         (this%ntracers));    this%is_isotope(:)      = .false.
  allocate(this%is_frozen          (this%ntracers));    this%is_frozen(:)       = .false.
  allocate(this%adsorbgroupid      (this%ntracers));    this%adsorbgroupid(:)   = nanid
  allocate(this%adsorb_isotherm    (this%ntracers));    this%adsorb_isotherm(:) = nanid
  allocate(this%adsorbid           (this%ntracers));    this%adsorbid(:)        = nanid

  allocate(this%volatileid         (this%ntracers));    this%volatileid(:)      = nanid
  allocate(this%volatilegroupid    (this%ntracers));    this%volatilegroupid(:) = nanid
  allocate(this%h2oid              (this%nh2o_tracers));         this%h2oid(:)           = nanid
  allocate(this%frozenid           (this%ntracers));    this%frozenid(:)        = nanid
  allocate(this%tracernames        (this%ntracers));             this%tracernames(:)     = ''
  allocate(this%tracerfamilyname   (this%ntracers));             this%tracerfamilyname(:)= ''
  allocate(this%units              (this%ntracers));             this%units(:)           = 'mol m-3'
  allocate(this%vtrans_scal        (this%ngwmobile_tracers));    this%vtrans_scal(:)     = 0._r8   !no transport through xylem transpiration

  allocate(this%tracer_solid_passive_diffus_scal_group(this%nsolid_passive_tracer_groups));
  this%tracer_solid_passive_diffus_scal_group(:) = 1._r8

  allocate(this%tracer_solid_passive_diffus_thc_group (this%nsolid_passive_tracer_groups));
  this%tracer_solid_passive_diffus_thc_group(:) = 1e-4_r8 / (86400._r8 * 365._r8) * 1.e-36_r8

  allocate(this%tracer_group_memid(this%ntracer_groups, this%nmem_max));
  this%tracer_group_memid(:,:) = nanid

  allocate(this%solid_passive_tracer_groupid(this%nsolid_passive_tracer_groups, 1:this%nmem_max));
  this%solid_passive_tracer_groupid(:,:) = nanid

  allocate(this%groupid(this%ntracers)); this%groupid(:) = nanid

  !write(iulog,*)'BeTR: total number of tracers',this%ntracers
  end subroutine InitAllocate

!--------------------------------------------------------------------------------

subroutine set_tracer(this, bstatus, trc_id, trc_name, is_trc_mobile, is_trc_advective, trc_group_id, &
   trc_group_mem, is_trc_diffusive, is_trc_volatile, trc_volatile_id, trc_volatile_group_id, &
   is_trc_h2o, trc_vtrans_scal, is_trc_adsorb, trc_adsorbid, trc_adsorbgroupid, trc_sorpisotherm, &
   is_trc_frozen, trc_frozenid, trc_family_name)

! !DESCRIPTION:
! set up tracer property based on input configurations
  use BetrStatusType  , only : betr_status_type
  implicit none
   ! !ARGUMENTS:
  class(BeTRtracer_type), intent(inout) :: this
  integer            , intent(in) :: trc_id
  character(len=*)   , intent(in) :: trc_name
  logical            , intent(in) :: is_trc_mobile
  logical            , intent(in) :: is_trc_advective
  integer            , intent(in) :: trc_group_id
  integer            , intent(in) :: trc_group_mem
  logical, optional  , intent(in) :: is_trc_diffusive
  logical, optional  , intent(in) :: is_trc_volatile
  integer, optional  , intent(in) :: trc_volatile_id
  integer, optional  , intent(in) :: trc_volatile_group_id
  real(r8),optional  , intent(in) :: trc_vtrans_scal
  logical ,optional  , intent(in) :: is_trc_h2o
  logical ,optional  , intent(in) :: is_trc_adsorb
  integer ,optional  , intent(in) :: trc_adsorbid
  integer ,optional  , intent(in) :: trc_adsorbgroupid
  character(len=*), optional, intent(in) :: trc_sorpisotherm
  logical ,optional  , intent(in) :: is_trc_frozen
  integer ,optional  , intent(in) :: trc_frozenid
  character(len=*),optional,intent(in) :: trc_family_name
  type(betr_status_type), intent(out)   :: bstatus

  call bstatus%reset()

  this%tracernames      (trc_id)=trim(trc_name)
  if(present(trc_family_name))then
    this%tracerfamilyname(trc_id)=trim(trc_family_name)
  else
    this%tracerfamilyname(trc_id)=trim(this%tracernames(trc_id))
  endif
  this%is_mobile        (trc_id)    = is_trc_mobile
  this%groupid          (trc_id)    = trc_group_id
  this%tracer_group_memid(trc_group_id,trc_group_mem) = trc_id
  this%is_advective     (trc_id)    = is_trc_advective

  if(present(is_trc_diffusive)) then
    this%is_diffusive (trc_id) = is_trc_diffusive
  endif
  if(present(is_trc_volatile))then
    this%is_volatile      (trc_id)    = is_trc_volatile
    if(this%is_volatile   (trc_id)) then
      if(.not.present(trc_volatile_id) .and. do_betr_output)then
        call bstatus%set_msg(msg='volatile tracer id is not provided for ' &
            //trim(trc_name)//errMsg(mod_filename, __LINE__),err=-1)
        return
      endif
      if(.not.present(trc_volatile_group_id) .and. do_betr_output)then
        call bstatus%set_msg(msg='volatile tracer group id is not provided for ' &
            //trim(trc_name)//errMsg(mod_filename, __LINE__), err=-1)
        return
      endif
      this%volatileid     (trc_id)    = trc_volatile_id
      this%volatilegroupid(trc_id)    = trc_volatile_group_id
      this%nvolatile_tracers          = this%nvolatile_tracers + 1
    endif
  endif

  if(present(trc_vtrans_scal))then
    this%vtrans_scal(trc_id) = trc_vtrans_scal
  endif

  if(present(is_trc_h2o))then
    this%is_h2o(trc_id) = is_trc_h2o
    if(is_trc_h2o)then
      this%nh2o_tracers = this%nh2o_tracers + 1
    endif
  endif

  if(present(is_trc_adsorb))then
    this%is_adsorb(trc_id) = is_trc_adsorb
    if(is_trc_adsorb)then
      if(.not.present(trc_adsorbid) .and. do_betr_output)then
        call bstatus%set_msg(msg='adsorb tracer id is not provided for ' &
            //trim(trc_name)//errMsg(mod_filename, __LINE__), err=-1)
        return
      endif
      if(.not.present(trc_adsorbgroupid) .and. do_betr_output)then
        call bstatus%set_msg(msg='adsorb tracer group id is not provided for ' &
           //trim(trc_name)//errMsg(mod_filename, __LINE__), err=-1)
        return
      endif
      if(.not.present(trc_sorpisotherm) .and. do_betr_output)then
        call bstatus%set_msg(msg='adsorb isotherm is not provided for ' &
           //trim(trc_name)//errMsg(mod_filename, __LINE__), err=-1)
        return
      endif
      this%adsorbid(trc_id) = trc_adsorbid
      this%adsorbgroupid(trc_id) = trc_adsorbgroupid
      this%nsolid_equil_tracers = this%nsolid_equil_tracers + 1
      if(trim(trc_sorpisotherm)=='LANGMUIR')then
         this%adsorb_isotherm(trc_id)=sorp_isotherm_langmuir
      elseif(trim(trc_sorpisotherm)=='LINEAR')then
         this%adsorb_isotherm(trc_id)=sorp_isotherm_linear
      endif
    endif
  endif

  if(present(is_trc_frozen))then
    this%is_frozen(trc_id) = is_trc_frozen
    if(is_trc_frozen)then
      if(.not. present(trc_frozenid) .and. do_betr_output)then
        call bstatus%set_msg(msg='frozen tracer id is not provided for ' &
            //trim(trc_name)//errMsg(mod_filename, __LINE__), err=-1)
        return
      endif
      this%frozenid(trc_id) = trc_frozenid
      this%nfrozen_tracers = this%nfrozen_tracers + 1
    endif
  endif
  end subroutine set_tracer

  !--------------------------------------------------------------------------------
  function is_solidtransport(this)result(yesno)

   ! !ARGUMENTS:
  implicit none
  class(BeTRtracer_type), intent(inout) :: this

  logical :: yesno

  yesno = (this%nsolid_passive_tracer_groups > 0)

  end function is_solidtransport

  !--------------------------------------------------------------------------------
  subroutine  add_tracer_group(this, trc_grp_cnt, mem, &
     trc_cnt, trc_grp, trc_grp_beg, trc_grp_end, is_trc_gw, is_trc_volatile,&
     is_trc_passive)
  implicit none
  class(BeTRtracer_type), intent(inout) :: this
  integer, intent(in)  :: trc_grp_cnt
  integer, intent(in)  :: mem
  integer, intent(inout):: trc_cnt
  integer, intent(out) :: trc_grp
  integer, intent(out):: trc_grp_beg
  integer, intent(out):: trc_grp_end
  logical, optional, intent(in)  :: is_trc_gw
  logical, optional, intent(in)  :: is_trc_volatile
  logical, optional, intent(in)  :: is_trc_passive

  trc_grp = trc_grp_cnt
  this%ntracers = this%ntracers + mem

  if(present(is_trc_gw))then
    if(is_trc_gw)then
      this%ngwmobile_tracers= this%ngwmobile_tracers + mem
      this%ngwmobile_tracer_groups = this%ngwmobile_tracer_groups + 1
    endif
  endif
  if(present(is_trc_volatile))then
    if(is_trc_volatile)then
      this%nvolatile_tracers            = this%nvolatile_tracers + mem                   !
      this%nvolatile_tracer_groups      = this%nvolatile_tracer_groups + 1               !
    endif
  endif

  if(present(is_trc_passive))then
    if(is_trc_passive)then
      this%nsolid_passive_tracer_groups = this%nsolid_passive_tracer_groups + 1
      this%nsolid_passive_tracers = this%nsolid_passive_tracers + mem
    endif
  endif
  trc_grp_beg = trc_cnt + 1
  trc_grp_end = trc_cnt + mem
  trc_cnt = trc_cnt + mem
  end subroutine add_tracer_group

end module BeTRTracerType
