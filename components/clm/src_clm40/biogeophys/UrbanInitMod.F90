module UrbanInitMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: UrbanInitMod
! 
! !DESCRIPTION: 
! Initialize urban data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun  
  use shr_sys_mod , only : shr_sys_flush 
  use clm_varctl  , only : iulog, use_vancouver, use_mexicocity
  use UrbanMod,     only : urban_traffic, urban_hac, urban_hac_off
!
! !PUBLIC TYPES:
  implicit none
  save

  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanInitTimeVar   ! Initialize urban time varying variables
  public :: UrbanInitTimeConst ! Initialize urban time constant variables
  public :: UrbanInitAero      ! Calculate urban landunit aerodynamic constants
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInitAero
!
! !INTERFACE:
  subroutine UrbanInitAero( )
!
! !DESCRIPTION: 
! Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
! Grimmond and Oke (1999)
!
! !USES:
    use clmtype 
    use clm_varcon, only : isturb, vkc
    use decompMod , only : get_proc_bounds
!
! !ARGUMENTS:
    implicit none
!
! local pointers to original implicit in arguments (urban clump)
!
    real(r8), pointer :: ht_roof(:)    ! height of urban roof (m)
    real(r8), pointer :: canyon_hwr(:) ! ratio of building height to street width (-)
    integer , pointer :: ltype(:)      ! landunit type
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: z_0_town(:)   ! urban landunit momentum roughness length (m)
    real(r8), pointer :: z_d_town(:)   ! urban landunit displacement height (m)
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Keith Oleson January 2005
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8), parameter :: alpha = 4.43_r8 ! coefficient used to calculate z_d_town
    real(r8), parameter :: beta = 1.0_r8   ! coefficient used to calculate z_d_town
    real(r8), parameter :: C_d = 1.2_r8    ! drag coefficient as used in Grimmond and Oke (1999)
    real(r8) :: plan_ai                    ! plan area index - ratio building area to plan area (-)
    real(r8) :: frontal_ai                 ! frontal area index of buildings (-)
    real(r8) :: build_lw_ratio             ! building short/long side ratio (-)
    integer  :: l,g                        ! indices
    integer  :: begp, endp                 ! clump beginning and ending pft indices
    integer  :: begc, endc                 ! clump beginning and ending column indices
    integer  :: begl, endl                 ! clump beginning and ending landunit indices
    integer  :: begg, endg                 ! clump beginning and ending gridcell indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit level)

    ltype      => lun%itype
    z_0_town   => lun%z_0_town
    z_d_town   => lun%z_d_town
    ht_roof    => lun%ht_roof
    canyon_hwr => lun%canyon_hwr

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do l = begl, endl 
      if (ltype(l) == isturb) then 

         ! Calculate plan area index 
         plan_ai = canyon_hwr(l)/(canyon_hwr(l) + 1._r8)

         ! Building shape shortside/longside ratio (e.g. 1 = square )
         ! This assumes the building occupies the entire canyon length
         build_lw_ratio = plan_ai

         ! Calculate frontal area index
         frontal_ai = (1._r8 - plan_ai) * canyon_hwr(l)

         ! Adjust frontal area index for different building configuration
         frontal_ai = frontal_ai * sqrt(1/build_lw_ratio) * sqrt(plan_ai)
         
         ! Calculate displacement height
         
         if (use_vancouver) then
            z_d_town(l) = 3.5_r8
         else if (use_mexicocity) then
            z_d_town(l) = 10.9_r8
         else
            z_d_town(l) = (1._r8 + alpha**(-plan_ai) * (plan_ai - 1._r8)) * ht_roof(l)
         end if

         ! Calculate the roughness length
         
         if (use_vancouver) then
            z_0_town(l) = 0.35_r8
         else if (use_mexicocity) then
            z_0_town(l) = 2.2_r8
         else
            z_0_town(l) = ht_roof(l) * (1._r8 - z_d_town(l) / ht_roof(l)) * &
                          exp(-1.0_r8 * (0.5_r8 * beta * C_d / vkc**2 * &
                          (1 - z_d_town(l) / ht_roof(l)) * frontal_ai)**(-0.5_r8))
         end if
      end if
   end do

 end subroutine UrbanInitAero

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInitTimeConst
!
! !INTERFACE:
  subroutine UrbanInitTimeConst()
!
! !DESCRIPTION: 
! Initialize urban time-constant variables
!
! !USES:
    use clmtype 
    use clm_varcon   , only : isturb, icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv, spval
    use decompMod    , only : get_proc_bounds, ldecomp
    use UrbanInputMod, only : urbinp
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: gdc(:)                 ! grid index for landunit 
    integer , pointer :: coli(:)                ! beginning column index for landunit 
    integer , pointer :: colf(:)                ! ending column index for landunit
    integer , pointer :: ctype(:)               ! column type
    integer , pointer :: ltype(:)               ! landunit type index
    integer , pointer :: lgridcell(:)           ! gridcell of corresponding landunit
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: canyon_hwr(:)          ! urban canyon height to width ratio
    real(r8), pointer :: emg(:)                 ! ground emissivity
    real(r8), pointer :: wtroad_perv(:)         ! weight of pervious column to total road
    real(r8), pointer :: ht_roof(:)             ! height of urban roof (m)
    real(r8), pointer :: wtlunit_roof(:)        ! weight of roof with respect to landunit
    real(r8), pointer :: wind_hgt_canyon(:)     ! height above road at which wind in canyon is to be computed (m)
    real(r8), pointer :: eflx_traffic_factor(:) ! multiplicative factor for sensible heat flux from urban traffic
    real(r8), pointer :: t_building_max(:)      ! maximum internal building temperature (K)
    real(r8), pointer :: t_building_min(:)      ! minimum internal building temperature (K)
    real(r8), pointer :: tk_wall(:,:)           ! thermal conductivity of urban wall (W/m/K)
    real(r8), pointer :: tk_roof(:,:)           ! thermal conductivity of urban roof (W/m/K)
    real(r8), pointer :: tk_improad(:,:)        ! thermal conductivity of urban impervious road (W/m/K)
    real(r8), pointer :: cv_wall(:,:)           ! thermal conductivity of urban wall (J/m^3/K)
    real(r8), pointer :: cv_roof(:,:)           ! thermal conductivity of urban roof (J/m^3/K)
    real(r8), pointer :: cv_improad(:,:)        ! thermal conductivity of urban impervious road (J/m^3/K)
    real(r8), pointer :: thick_wall(:)          ! thickness of urban wall (m)
    real(r8), pointer :: thick_roof(:)          ! thickness of urban roof (m)
    integer,  pointer :: nlev_improad(:)        ! number of impervious road layers (-)
!
!
! !OTHER LOCAL VARIABLES
!EOP
    integer  :: nc,fl,ib,l,c,p,g          ! indices
    integer  :: ier                       ! error status
    integer  :: begp, endp                ! clump beginning and ending pft indices
    integer  :: begc, endc                ! clump beginning and ending column indices
    integer  :: begl, endl                ! clump beginning and ending landunit indices
    integer  :: begg, endg                ! clump beginning and ending gridcell indices

    ! Assign local pointers to derived type members (landunit-level)

    ltype               => lun%itype
    lgridcell           => lun%gridcell
    coli                => lun%coli
    colf                => lun%colf
    canyon_hwr          => lun%canyon_hwr
    wtroad_perv         => lun%wtroad_perv 
    ht_roof             => lun%ht_roof
    wtlunit_roof        => lun%wtlunit_roof
    wind_hgt_canyon     => lun%wind_hgt_canyon
    eflx_traffic_factor => lef%eflx_traffic_factor
    t_building_max      => lps%t_building_max
    t_building_min      => lps%t_building_min
    canyon_hwr          => lun%canyon_hwr
    tk_wall             => lps%tk_wall
    tk_roof             => lps%tk_roof
    tk_improad          => lps%tk_improad
    cv_wall             => lps%cv_wall
    cv_roof             => lps%cv_roof
    cv_improad          => lps%cv_improad
    thick_wall          => lps%thick_wall
    thick_roof          => lps%thick_roof
    nlev_improad        => lps%nlev_improad

    ! Assign local pointers to derived type members (column-level)

    ctype               => col%itype
    emg                 => cps%emg
    
   ! Initialize time constant urban variables

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do l = begl, endl
       if (ltype(l) == isturb) then
          g = lun%gridcell(l)
          canyon_hwr(l)         = urbinp%canyon_hwr(g)
          wtroad_perv(l)        = urbinp%wtroad_perv(g)
          ht_roof(l)            = urbinp%ht_roof(g)
          wtlunit_roof(l)       = urbinp%wtlunit_roof(g)
          wind_hgt_canyon(l)    = urbinp%wind_hgt_canyon(g)
          tk_wall(l,:)          = urbinp%tk_wall(g,:)
          tk_roof(l,:)          = urbinp%tk_roof(g,:)
          tk_improad(l,:)       = urbinp%tk_improad(g,:)
          cv_wall(l,:)          = urbinp%cv_wall(g,:)
          cv_roof(l,:)          = urbinp%cv_roof(g,:)
          cv_improad(l,:)       = urbinp%cv_improad(g,:)
          thick_wall(l)         = urbinp%thick_wall(g)
          thick_roof(l)         = urbinp%thick_roof(g)
          nlev_improad(l)       = urbinp%nlev_improad(g)
          t_building_min(l)     = urbinp%t_building_min(g)
          t_building_max(l)     = urbinp%t_building_max(g)

          do c = coli(l),colf(l)
             if (ctype(c) == icol_roof       ) emg(c) = urbinp%em_roof(g)
             if (ctype(c) == icol_sunwall    ) emg(c) = urbinp%em_wall(g)
             if (ctype(c) == icol_shadewall  ) emg(c) = urbinp%em_wall(g)
             if (ctype(c) == icol_road_imperv) emg(c) = urbinp%em_improad(g)
             if (ctype(c) == icol_road_perv  ) emg(c) = urbinp%em_perroad(g)
          end do

          ! Inferred from Sailor and Lu 2004
          if (urban_traffic) then
             eflx_traffic_factor(l) = 3.6_r8 * (canyon_hwr(l)-0.5_r8) + 1.0_r8
          else
             eflx_traffic_factor(l) = 0.0_r8
          end if

          if (use_vancouver .or. use_mexicocity) then 
             ! Freely evolving
             t_building_max(l) = 380.00_r8
             t_building_min(l) = 200.00_r8
          else
             if (urban_hac == urban_hac_off) then
                ! Overwrite values read in from urbinp by freely evolving values
                t_building_max(l) = 380.00_r8
                t_building_min(l) = 200.00_r8
             end if
          end if
       else
          eflx_traffic_factor(l) = 0.0_r8
          t_building_max(l) = 0.0_r8
          t_building_min(l) = 0.0_r8
       end if
    end do

  end subroutine UrbanInitTimeConst

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInitTimeVar
!
! !INTERFACE:
  subroutine UrbanInitTimeVar( )
!
! !DESCRIPTION: 
! Initialize urban time-varying variables
!
! !USES:
    use clmtype
    use clm_varcon, only : isturb, spval, icol_road_perv
    use decompMod , only : get_proc_bounds
!
! !ARGUMENTS:
    implicit none
!
! local pointers to original implicit in arguments (urban clump)
!
    integer , pointer :: ltype(:)      ! landunit type
    integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
    integer , pointer :: clandunit(:)  ! landunit index of corresponding column
    integer , pointer :: plandunit(:)  ! landunit index of corresponding pft
    integer , pointer :: ctype(:)      ! column type
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: taf(:)                ! urban canopy air temperature (K)
    real(r8), pointer :: qaf(:)                ! urban canopy air specific humidity (kg/kg)
    real(r8), pointer :: eflx_building_heat(:) ! heat flux from urban building interior to walls, roof (W/m**2)
    real(r8), pointer :: eflx_urban_ac(:)      ! urban air conditioning flux (W/m**2)
    real(r8), pointer :: eflx_urban_heat(:)    ! urban heating flux (W/m**2)
    real(r8), pointer :: fcov(:)               ! fractional impermeable area
    real(r8), pointer :: fsat(:)               ! fractional area with water table at surface
    real(r8), pointer :: qcharge(:)            ! aquifer recharge rate (mm/s)
    real(r8), pointer :: t_building(:)         ! internal building temperature (K)
    real(r8), pointer :: eflx_traffic(:)       ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat(:)     ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_wasteheat_pft(:) ! sensible heat flux from urban heating/cooling sources of waste heat at pft level (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac_pft(:) ! sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: eflx_traffic_pft(:)   ! sensible heat flux from traffic (W/m**2)
    real(r8), pointer :: eflx_anthro(:)        ! total anthropogenic heat flux (W/m**2)
    real(r8), pointer :: t_ref2m_u(:)          ! Urban 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_min_u(:)      ! Urban daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_u(:)      ! Urban daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: rh_ref2m_u(:)         ! Urban 2 m height surface relative humidity (%)
    real(r8), pointer :: t_grnd_u(:)           ! Urban ground temperature (Kelvin)
    real(r8), pointer :: qflx_runoff_u(:)      ! Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
    real(r8), pointer :: fsa_u(:)              ! Urban absorbed solar radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net_u(:)   ! Urban net longwave radiation (W/m**2)
    real(r8), pointer :: eflx_lh_tot_u(:)      ! Urban latent heat flux (W/m**2)
    real(r8), pointer :: eflx_sh_tot_u(:)      ! Urban sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_soil_grnd_u(:)   ! Urban ground heat flux (W/m**2)
    real(r8), pointer :: eflx_snomelt_u(:)     ! Urban snow melt heat flux (W/m**2)
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Keith Oleson February 2005
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: l,g,c,p       ! indices
    integer :: begp, endp    ! clump beginning and ending pft indices
    integer :: begc, endc    ! clump beginning and ending column indices
    integer :: begl, endl    ! clump beginning and ending landunit indices
    integer :: begg, endg    ! clump beginning and ending gridcell indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit level)

    taf                => lps%taf
    qaf                => lps%qaf
    ltype              => lun%itype
    lgridcell          => lun%gridcell
    t_building         => lps%t_building
    eflx_traffic       => lef%eflx_traffic
    eflx_wasteheat     => lef%eflx_wasteheat

    ! Assign local pointers to derived type members (column level)

    clandunit          => col%landunit
    eflx_building_heat => cef%eflx_building_heat
    eflx_urban_ac      => cef%eflx_urban_ac
    eflx_urban_heat    => cef%eflx_urban_heat
    fcov               => cws%fcov
    fsat               => cws%fsat
    qcharge            => cws%qcharge
    ctype              => col%itype
    t_grnd_u           => ces%t_grnd_u
    qflx_runoff_u      => cwf%qflx_runoff_u
    eflx_snomelt_u     => cef%eflx_snomelt_u

    ! Assign local pointers to derived type members (pft level)

    t_ref2m_u          => pes%t_ref2m_u
    t_ref2m_min_u      => pes%t_ref2m_min_u
    t_ref2m_max_u      => pes%t_ref2m_max_u
    rh_ref2m_u         => pes%rh_ref2m_u
    plandunit          => pft%landunit
    eflx_wasteheat_pft => pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => pef%eflx_heat_from_ac_pft
    eflx_traffic_pft => pef%eflx_traffic_pft
    eflx_anthro => pef%eflx_anthro
    fsa_u              => pef%fsa_u
    eflx_lwrad_net_u   => pef%eflx_lwrad_net_u
    eflx_lh_tot_u      => pef%eflx_lh_tot_u
    eflx_sh_tot_u      => pef%eflx_sh_tot_u
    eflx_soil_grnd_u   => pef%eflx_soil_grnd_u

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do l = begl, endl 
       g = lgridcell(l)
       if (ltype(l) == isturb) then 
          if (use_vancouver) then
             taf(l) = 297.56_r8
             qaf(l) = 0.0111_r8
          else if (use_mexicocity) then
             taf(l) = 289.46_r8
             qaf(l) = 0.00248_r8
          else
             taf(l) = 283._r8
             ! Arbitrary set since forc_q is not yet available
             qaf(l) = 1.e-4_r8
          end if
       else
          t_building(l)     = 0.0_r8
          eflx_traffic(l)   = 0.0_r8
          eflx_wasteheat(l) = 0.0_r8
       end if
    end do

    do c = begc, endc 
       l = clandunit(c)
       if (ltype(l) == isturb) then 
          eflx_building_heat(c) = 0._r8
          eflx_urban_ac(c) = 0._r8
          eflx_urban_heat(c) = 0._r8
          !
          ! Set hydrology variables for urban to spvalue -- as only valid for pervious road
          !
          if (ctype(c) /= icol_road_perv  )then
             fcov(c)    = spval
             fsat(c)    = spval
             qcharge(c) = spval
          end if
       else
          eflx_building_heat(c) = 0._r8
          eflx_urban_ac(c) = 0._r8
          eflx_urban_heat(c) = 0.0_r8
          t_grnd_u(c) = 0.0_r8
          qflx_runoff_u(c) = 0.0_r8
          eflx_snomelt_u(c) = 0.0_r8
       end if
    end do

    do p = begp, endp 
       l = plandunit(p)
       if (ltype(l) /= isturb) then 
          t_ref2m_u(p)     = 0.0_r8
          t_ref2m_min_u(p) = 0.0_r8
          t_ref2m_max_u(p) = 0.0_r8
          rh_ref2m_u(p)     = 0.0_r8
          eflx_wasteheat_pft(p) = 0.0_r8
          eflx_heat_from_ac_pft(p) = 0.0_r8
          eflx_traffic_pft(p) = 0.0_r8
          eflx_anthro(p)    = 0.0_r8
          fsa_u(p)            = 0.0_r8
          eflx_lwrad_net_u(p) = 0.0_r8
          eflx_lh_tot_u(p)    = 0.0_r8
          eflx_sh_tot_u(p)    = 0.0_r8
          eflx_soil_grnd_u(p) = 0.0_r8
       else
          eflx_wasteheat_pft(p)    = 0.0_r8
          eflx_heat_from_ac_pft(p) = 0.0_r8
          eflx_traffic_pft(p)      = 0.0_r8
       end if
    end do
    
  end subroutine UrbanInitTimeVar
  
end module UrbanInitMod
