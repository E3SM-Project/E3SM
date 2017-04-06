module TracerParamsMod
#include "shr_assert.h"

  ! !DESCRIPTION:
  ! Module holding routines used to compute solubility, and phase conversion parameters
  ! to be used for BeTR 1D vertical tracer transport
  !
  !
  ! History
  ! Jinyun Tang created May 2014.
  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use decompMod             , only : bounds_type
  use clm_varpar            , only : nlevsoi
  use clm_varcon            , only : spval
  use PatchType             , only : pft
  use ColumnType            , only : col
  use tracer_varcon
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: tracer_param_init
  public :: set_multi_phase_diffusion
  public :: get_gas_diffusivity
  public :: get_diffusivity_ratio_gas2h2o, get_henrycef
  public :: convert_mobile2gas
  public :: set_phase_convert_coeff
  public :: calc_tracer_infiltration
  public :: pre_diagnose_soilcol_water_flux
  public :: diagnose_advect_water_flux
  public :: diagnose_drainage_water_flux
  public :: calc_smp_l
  public :: get_zwt
  public :: calc_aerecond
  public :: betr_annualupdate
  !parameters
  real(r8), parameter :: minval_diffus = 1.e-20_r8   !minimum diffusivity, m2/s
  real(r8), parameter :: minval_airvol = 1.e-10_r8   !minimum air-filled volume


  !declare a private tortuosity type
  type :: soil_tortuosity_type
     real(r8), pointer :: tau_gas(:,:)      !soil tortuosity for gaseous phase diffusion
     real(r8), pointer :: tau_liq(:,:)      !soil tortuosity for aqueous phase diffusion
  end type soil_tortuosity_type
  type(soil_tortuosity_type), target :: tau_soil
  real(r8), private, pointer :: h2osoi_liq_copy(:,:)
  !!
contains


  subroutine tracer_param_init(bounds)

    !
    ! !DESCRIPTION:
    !
    ! initialize the tracerParamsMod
    !
    ! !USES:
    use clm_varpar         , only : nlevtrc_soil

    implicit none
    type(bounds_type), intent(in) :: bounds   !bounds
    character(len=32)             :: subname ='tracer_param_init'

    allocate(tau_soil%tau_gas(bounds%begc:bounds%endc, 1 : nlevtrc_soil))
    tau_soil%tau_gas(:,:) = 0._r8
    allocate(tau_soil%tau_liq(bounds%begc:bounds%endc, 1 : nlevtrc_soil))
    tau_soil%tau_liq(:,:) = 0._r8

    allocate(h2osoi_liq_copy(bounds%begc:bounds%endc, 1:nlevsoi));  h2osoi_liq_copy(:, :) = spval
  end subroutine tracer_param_init

  !--------------------------------------------------------------------------------------------------------------
  subroutine Calc_gaseous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, soilstate_vars, waterstate_vars, tau_gas)
    !
    ! !DESCRIPTION:
    !
    ! compute soil tortuosity for gasesous diffusion

    ! !USES:
    use SoilStateType      , only : soilstate_type
    use WaterStateType     , only : Waterstate_Type

    implicit none
    !arguments
    type(bounds_type),        intent(in) :: bounds                             ! bounds
    integer,                  intent(in) :: num_soilc                          ! number of column soil points in column filter
    integer,                  intent(in) :: filter_soilc(:)                    ! column filter for soil points
    integer,                  intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0
    integer,                  intent(in) :: jtops(bounds%begc: )               ! top label of each column
    type(soilstate_type),     intent(in) :: soilstate_vars                     ! column soil physical state variables
    type(Waterstate_Type),    intent(in) :: waterstate_vars                    ! column soil water state variables
    real(r8),              intent(inout) :: tau_gas(bounds%begc: , lbj: )      !output variable

    !local variables
    integer :: n, fc, c     !indices
    character(len=255) :: subname = 'calc_gaseous_diffusion_soil_tortuosity'

    associate( &
         eff_porosity   => soilstate_vars%eff_porosity_col,                & !effective soil porosity
         bsw            => soilstate_vars%bsw_col         ,                & !clapp-hornber shape parameters
         air_vol        => waterstate_vars%air_vol_col                     & !volume possessed by air
         )

      SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
      SHR_ASSERT_ALL((ubound(tau_gas)         == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
      do n = lbj, ubj
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            if(n>=jtops(c))then
               tau_gas(c,n) = get_taugas(eff_porosity(c,n), air_vol(c,n),bsw(c,n))
            endif
         enddo
      enddo
    end associate

  end subroutine Calc_gaseous_diffusion_soil_tortuosity
  !--------------------------------------------------------------------------------------------------------------
  subroutine Calc_aqueous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, numf, filter, soilstate_vars, waterstate_vars, tau_liq)
    !
    ! DESCRIPTIONS
    ! compute soil tortuosity for aquesous diffusion
    !

    use SoilStateType      , only : soilstate_type
    use WaterStateType     , only : Waterstate_Type

    implicit none

    !arguments
    type(bounds_type),        intent(in) :: bounds                                ! bounds
    integer,                  intent(in) :: numf                                  ! number of columns in column filter
    integer,                  intent(in) :: filter(:)                             ! column filter
    integer,                  intent(in) :: lbj, ubj                              ! lower and upper bounds, make sure they are > 0
    integer,                  intent(in) :: jtops(bounds%begc: )                  ! top label of each column
    type(soilstate_type),     intent(in) :: soilstate_vars                        ! column soil physical state variables
    type(Waterstate_Type),    intent(in) :: waterstate_vars                       ! column soil water state variables
    real(r8),              intent(inout) :: tau_liq(bounds%begc: , lbj: )         !output variable

    !local variables
    integer :: n, fc, c     !indices
    character(len=255) :: subname = 'calc_aqueous_diffusion_soil_tortuosity'

    SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(tau_liq)         == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))

    associate( &
         eff_porosity   => soilstate_vars%eff_porosity_col,                & !effective soil porosity
         bsw            => soilstate_vars%bsw_col         ,                & !clapp-hornber shape parameters
         h2osoi_liqvol  => waterstate_vars%h2osoi_liqvol_col               & !soil volume possessed by liquid water
         )

      do n = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(n>=jtops(c))then
               tau_liq(c,n)=get_tauliq(eff_porosity(c,n), h2osoi_liqvol(c,n),bsw(c,n))
            endif
         enddo
      enddo
    end associate
  end subroutine calc_aqueous_diffusion_soil_tortuosity

  !--------------------------------------------------------------------------------------------------------------
  
  subroutine calc_bulk_diffusivity(bounds, lbj, ubj, jtops, numf, filter, bunsencef_col, &
       canopystate_vars, waterstate_vars, tau_soi, betrtracer_vars, t_soisno,  bulkdiffus)
    !
    ! !DESCRIPTION:
    ! compute the weighted bulk diffusivity in soil for dual-phase transport
    ! Reference: Tang and Riley, 2014, BG,  Simple formulations and solutions of &
    ! the dual-phase diffusive transport for biogeochemical modeling.
    !the formula for a volatile species is
    !D_bulk=(airvol*D_g*tau_g+bunsencef_col*h2osoi_liqvol*D_w*tau_w)

    ! !USES:
    use WaterStateType        , only : Waterstate_Type
    use BeTRTracerType        , only : betrtracer_type
    use CanopyStateType       , only : canopystate_type
    use clm_varcon            , only : zisoi

    type(bounds_type)          , intent(in) :: bounds                                  ! bounds
    integer                    , intent(in) :: numf                                    ! number of columns in column filter
    integer                    , intent(in) :: filter(:)                               ! column filter
    integer                    , intent(in) :: lbj, ubj                                ! lower and upper bounds, make sure they are > 0
    integer                    , intent(in) :: jtops(bounds%begc: )                    ! top label of each column
    real(r8)                   , intent(in) :: t_soisno(bounds%begc: ,  lbj: )         ! soil temperature
    real(r8)                   , intent(in) :: bunsencef_col(bounds%begc: ,lbj: ,1: )  ! bunsen coefficient for gaseous-aqueous conversion
    type(betrtracer_type)      , intent(in) :: betrtracer_vars                         ! betr configuration information
    type(Waterstate_Type)      , intent(in) :: waterstate_vars                         ! water state variables
    type(canopystate_type)     , intent(in) :: canopystate_vars                        ! canopy state variables
    type(soil_tortuosity_type) , intent(in) :: tau_soi                                 ! soil tortuosity
    real(r8)                   ,intent(out) :: bulkdiffus(bounds%begc: ,lbj: , 1: )    ! the returning variable

    !local variables
    real(r8) :: max_depth_cryoturb         = 3._r8  !m
    !parameters below will be encapsulated into a structure later
    real(r8) :: max_altdepth_cryoturbation = 1._r8  ! (m) maximum active layer thickness for cryoturbation to occur
    real(r8) :: cryoturb_diffusion_k       = 1e-4_r8 / (86400._r8 * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr = 1m^2/1000 yr
    real(r8) :: som_diffus                 = 5e-4_r8 / (86400._r8 * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr
    integer :: j, k, n, fc, c , trcid    !indices
    integer :: nsld
    real(r8) :: diffaqu, diffgas
    character(len=255) :: subname = 'calc_bulk_diffusivity'

    !array shape checking will be added later.
    SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno)        == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bunsencef_col, 1)== (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bunsencef_col, 2)== (/ubj/))        ,   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bunsencef_col, 3)== (/betrtracer_vars%nvolatile_tracer_groups/)),   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bulkdiffus, 1)   == (/bounds%endc/))                   ,   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bulkdiffus, 2)   == (/ubj/))                           ,   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bulkdiffus, 3)   == (/betrtracer_vars%ntracer_groups/)),   errMsg(__FILE__,__LINE__))

    associate(                                                                                              &
         ngwmobile_tracer_groups                => betrtracer_vars%ngwmobile_tracer_groups                , & ! Integer[intent(in)], number of dual phase (gw) tracers
         tracer_group_memid                     => betrtracer_vars%tracer_group_memid                     , & !
         ntracer_groups                         => betrtracer_vars%ntracer_groups                         , & ! Integer[intent(in)], total number of tracers
         is_volatile                            => betrtracer_vars%is_volatile                            , & ! logical[intent(in)], is a volatile tracer?
         is_h2o                                 => betrtracer_vars%is_h2o                                 , & ! logical[intent(in)], is a h2o tracer?
         volatilegroupid                        => betrtracer_vars%volatilegroupid                        , & ! integer[intent(in)], location in the volatile vector
         air_vol                                => waterstate_vars%air_vol_col                            , & ! volume possessed by air
         h2osoi_liqvol                          => waterstate_vars%h2osoi_liqvol_col                      , & ! soil volume possessed by liquid water
         altmax                                 => canopystate_vars%altmax_col                            , & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear                        => canopystate_vars%altmax_lastyear_col                   , & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth o
         tracer_solid_passive_diffus_scal_group => betrtracer_vars%tracer_solid_passive_diffus_scal_group , & !scaling factor for solid phase diffusivity
         tracer_solid_passive_diffus_thc_group  => betrtracer_vars%tracer_solid_passive_diffus_thc_group  , & !threshold for solid phase diffusivity
         tau_gas                                => tau_soi%tau_gas                                        , & ! real(r8)[intent(in)], gaseous tortuosity
         tau_liq                                => tau_soi%tau_liq                                          & ! real(r8)[intent(in)], aqueous tortuosity
         )

      bulkdiffus(:,:,:) = 1.e-40_r8                            !initialize to a very small number
      do j = 1, ngwmobile_tracer_groups
         trcid = tracer_group_memid(j,1)
         if(is_volatile(j))then
            !it is a volatile tracers
            k=volatilegroupid(trcid)

            do n=lbj, ubj
               do fc = 1, numf
                  c = filter(fc)
                  if(n>=jtops(c))then
                     !aqueous diffusivity
                     diffaqu=get_aqueous_diffusivity(trcid, t_soisno(c,n), betrtracer_vars)
                     !gaseous diffusivity
                     diffgas=get_gas_diffusivity(trcid, t_soisno(c,n), betrtracer_vars)

                     !bulk diffusivity
                     !the bulk diffusivity is calculated by assuming the diffusion equation is gas-primary
                     !accordingly the retardation factor is gas primary
                     if(is_h2o(trcid))then
                        !for water tracer, the aqueous phase is used as dominant species
                        bulkdiffus(c,n,j)=air_vol(c,n)*tau_gas(c,n)*diffgas/bunsencef_col(c,n,k)+ &
                             h2osoi_liqvol(c,n)*tau_liq(c,n)*diffaqu
                     else
                        bulkdiffus(c,n,j)=air_vol(c,n)*tau_gas(c,n)*diffgas+ &
                             h2osoi_liqvol(c,n)*tau_liq(c,n)*diffaqu*bunsencef_col(c,n,k)
                     endif
                     !to prevent division by zero
                     bulkdiffus(c,n,j)=max(bulkdiffus(c,n,j),minval_diffus)
                  endif
               enddo
            enddo

         else
            !it is not a volatile tracer
            do n = lbj, ubj
               do fc = 1, numf
                  c = filter(fc)
                  if(n>=jtops(c))then
                     !the retardation factor is 1.
                     diffaqu=get_aqueous_diffusivity(trcid, t_soisno(c,n), betrtracer_vars)
                     bulkdiffus(c,n,j)=diffaqu*h2osoi_liqvol(c,n)*tau_liq(c,n)
                     !to prevent division by zero
                     bulkdiffus(c,n,j)=max(bulkdiffus(c,n,j),minval_diffus) !avoid division by zero in following calculations
                  endif
               enddo
            enddo
         endif
      enddo

      !do solid phase passive tracers
      do j = ngwmobile_tracer_groups + 1, ntracer_groups
         nsld = j - ngwmobile_tracer_groups
         trcid = tracer_group_memid(j,1)
         do n = 1, ubj
            do fc = 1,numf
               c = filter(fc)

               if  ( ( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
                    ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
                  ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth

                  if ( zisoi(n) < max(altmax(c), altmax_lastyear(c)) ) then
                     bulkdiffus(c,n,j) = cryoturb_diffusion_k * tracer_solid_passive_diffus_scal_group(nsld)
                     bulkdiffus(c,n,j) = max(bulkdiffus(c,n,j), tracer_solid_passive_diffus_thc_group(nsld))
                  else
                     bulkdiffus(c,n,j) = max(cryoturb_diffusion_k * &
                          ( 1._r8 - ( zisoi(n) - max(altmax(c), altmax_lastyear(c)) ) / &
                          ( max_depth_cryoturb - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
                     bulkdiffus(c,n,j) = bulkdiffus(c,n,j) * tracer_solid_passive_diffus_scal_group(nsld)
                     bulkdiffus(c,n,j) = max(bulkdiffus(c,n,j), tracer_solid_passive_diffus_thc_group(nsld))
                  endif
               elseif (  max(altmax(c), altmax_lastyear(c)) > 0._r8 ) then
                  ! constant advection, constant diffusion
                  bulkdiffus(c,n,j) = som_diffus * tracer_solid_passive_diffus_scal_group(nsld)
                  bulkdiffus(c,n,j) = max(bulkdiffus(c,n,j), tracer_solid_passive_diffus_thc_group(nsld))
               else
                  ! completely frozen soils--no mixing
                  bulkdiffus(c,n,j) = 1e-4_r8 / (86400._r8 * 365._r8) * 1.e-36_r8  !set to very small number for numerical purpose
               endif
            enddo
         enddo
      enddo
    end associate
  end subroutine calc_bulk_diffusivity
!--------------------------------------------------------------------------------------------------------------


   subroutine calc_bulk_conductances(bounds, lbj, ubj, jtops, numf, filter, bulkdiffus, dz, betrtracer_vars,  hmconductance_col)
   !
   ! DESCRIPTIONS:
   ! Compute weighted conductances for diffusive/dispersive tracer transport
   ! dispersion is not modeled currently.
   ! The computation of diffusivity assumes to implement the Fick's law of diffusion,
   ! which is supposed to be of good accuracy when the chemical is in trace amount and
   ! total air pressure changes with a small amount. The Stefan-Maxwell relationship
   ! could be implemented, but it is too complicate to gain much given other sources
   ! of uncertainty.
   ! Reference: Tang and Riley, 2014, BG,  Simple formulations and solutions of &
   ! the dual-phase diffusive transport for biogeochemical modeling.

   ! jyt, Jan 6, 2014
   ! !USES:
   use transportmod       , only : calc_interface_conductance
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   !arguments
   type(bounds_type),      intent(in) :: bounds                      ! bounds
   integer,                intent(in) :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
   integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
   integer,                intent(in) :: numf                        ! number of columns in column filter
   integer,                intent(in) :: filter(:)                   ! column filter
   type(betrtracer_type),  intent(in) :: betrtracer_vars             ! betr configuration information
   real(r8),               intent(in) :: bulkdiffus(bounds%begc: ,lbj: ,1: )  !weighted bulk diffusivity for dual-phase diffusion
   real(r8),               intent(in) :: dz(bounds%begc: , lbj: )

   real(r8),            intent(inout) :: hmconductance_col(bounds%begc: , lbj: ,1: ) !weighted bulk conductance

   !local variables

   integer :: j, n, fc, c                 !indices
   character(len=255) :: subname = 'calc_bulk_conductances'

   SHR_ASSERT_ALL((ubound(jtops)             == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(dz)                == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bulkdiffus,1)      == (/bounds%endc/)),                     errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bulkdiffus,2)      == (/ubj/)),                             errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bulkdiffus,3)      == (/betrtracer_vars%ntracer_groups/)),  errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(hmconductance_col,1) == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(hmconductance_col,2) == (/ubj-1/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(hmconductance_col,3) == (/betrtracer_vars%ntracer_groups/)),   errMsg(__FILE__,__LINE__))

   associate( &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups   , & !Integer[intent(in)], number of gw tracers
    ntracer_groups             => betrtracer_vars%ntracer_groups            , & !Integer[intent(in)], total number of tracers
    is_volatile                => betrtracer_vars%is_volatile               , & !logical[intent(in)], is a volatile tracer?
    is_mobile                  => betrtracer_vars%is_mobile                 , &
    tracer_group_memid         => betrtracer_vars%tracer_group_memid        , &
    volatileid                 => betrtracer_vars%volatileid                  & !integer[intent(in)], location in the volatile vector
   )

!  compute the depth weighted diffusivities
   do j = 1, ntracer_groups
     if(.not. is_mobile(tracer_group_memid(j,1)))cycle
     call calc_interface_conductance(bounds, lbj, ubj, jtops, numf, filter , &
             bulkdiffus(bounds%begc:bounds%endc, lbj:ubj, j)               , &
             dz(bounds%begc:bounds%endc, lbj:ubj)                          , &
             hmconductance_col(bounds%begc:bounds%endc, lbj:ubj-1, j))
   enddo

   end associate
   end subroutine calc_bulk_conductances

!-------------------------------------------------------------------------------
   subroutine calc_henrys_coeff(bounds, lbj, ubj, jtops, numf, filter, t_soisno, soi_pH, &
       betrtracer_vars, aqu2neutralcef_col, henrycef_col)
   !
   ! DESCRIPTION
   ! compute henry's law constant for volatile tracers
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   !arguments
   type(bounds_type),      intent(in) :: bounds  ! bounds
   integer,                intent(in) :: lbj, ubj        ! lower and upper bounds, make sure they are > 0
   integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
   integer,                intent(in) :: numf                                          ! number of columns in column filter
   integer,                intent(in) :: filter(:)                                     ! column filter
   real(r8),               intent(in) :: t_soisno(bounds%begc: ,  lbj: )   !soil temperature
   real(r8),               intent(in) :: soi_pH(bounds%begc: , lbj: )      !pH profile
   type(betrtracer_type),  intent(in) :: betrtracer_vars        ! betr configuration information

   real(r8),   intent(inout):: aqu2neutralcef_col(bounds%begc: , lbj: , 1: ) !conversion parameter between bulk aqueous and neutral aqueous tracer
   real(r8),   intent(inout):: henrycef_col(bounds%begc: , lbj: ,  1: )       !henry's constant, mol/L/atm = M/atm

   !local variables
   integer :: j, k, n, fc, c, trcid   ! indices
   real(r8) :: scal
   character(len=255) :: subname='calc_henrys_coeff'

   SHR_ASSERT_ALL((ubound(jtops)             == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(t_soisno)          == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(soi_pH)            == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(aqu2neutralcef_col,1)== (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(aqu2neutralcef_col,2)== (/ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(aqu2neutralcef_col,3)== (/betrtracer_vars%ngwmobile_tracer_groups/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(henrycef_col,1)      == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(henrycef_col,2)      == (/ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(henrycef_col,3)      == (/betrtracer_vars%nvolatile_tracer_groups/)),   errMsg(__FILE__,__LINE__))


   associate( &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups , & !Integer[intent(in)], number of tracers
    is_volatile                => betrtracer_vars%is_volatile             , & !logical[intent(in)], is a volatile tracer?
    is_h2o                     => betrtracer_vars%is_h2o                  , & !logical[intent(in)], is a h2o tracer?
    tracer_group_memid         => betrtracer_vars%tracer_group_memid      , & !integer[intent(in)], tracer id
    volatilegroupid            => betrtracer_vars%volatilegroupid           & !integer[intent(in)], location in the volatile vector
   )

   do j = 1, ngwmobile_tracer_groups

     !for tagged co2 simulations, the henry's constants are assumed same for all co2 tracers
     trcid= tracer_group_memid(j, 1)
     if(is_volatile(trcid) .and. (.not. is_h2o(trcid)))then
       k = volatilegroupid(trcid)
       do n = lbj, ubj
         do fc = 1, numf
           c = filter(fc)
           if(n>=jtops(c))then
             !Henry's law constants
             henrycef_col(c,n,k)=get_henrycef(t_soisno(c,n), trcid, betrtracer_vars)
             scal = get_equilibrium_scal(t_soisno(c,n), soi_pH(c,n), trcid,betrtracer_vars)
             henrycef_col(c,n,k)=henrycef_col(c,n,k) * scal
             aqu2neutralcef_col(c,n,j)=1._r8/scal   !this will convert the bulk aqueous phase into neutral phase

           endif
         enddo
       enddo
     endif
   enddo
   end associate
   end subroutine calc_henrys_coeff
!-------------------------------------------------------------------------------
   subroutine calc_bunsen_coeff(bounds, lbj, ubj, jtops, numf, filter, henrycef_col, t_soisno, smp_l, betrtracer_vars, bunsencef_col)
   !
   ! DESCRIPTION
   ! compute Bunsen's coefficient
   !
   use clm_varcon        , only : denh2o
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   !arguments
   type(bounds_type),      intent(in) :: bounds                                    ! bounds
   integer,                intent(in) :: lbj, ubj                                  ! lower and upper bounds, make sure they are > 0
   integer,                intent(in) :: jtops(bounds%begc: )                      ! top label of each column
   integer,                intent(in) :: numf                                      ! number of columns in column filter
   integer,                intent(in) :: filter(:)                                 ! column filter
   real(r8),               intent(in) :: t_soisno(bounds%begc: ,  lbj: )           !soil temperature, K
   real(r8),               intent(in) :: smp_l(bounds%begc: , lbj: )               !soil matric pressure, mm
   type(betrtracer_type),  intent(in) :: betrtracer_vars                           ! betr configuration information
   real(r8),               intent(in) :: henrycef_col(bounds%begc: , lbj: ,  1: )  !henry's constant
   real(r8),            intent(inout) :: bunsencef_col(bounds%begc: , lbj: , 1: )  !returning variable

   !local variables
   integer :: j, k, n, fc, c , trcid       !indices
   real(r8) :: rho_vap(bounds%begc:bounds%endc, lbj:ubj)                           ! saturated vapor pressure for different layers

   character(len=255) :: subname = 'calc_bunsen_coeff'

   SHR_ASSERT_ALL((ubound(jtops)             == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(t_soisno)          == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(smp_l)             == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(henrycef_col,1)    == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(henrycef_col,2)    == (/ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(henrycef_col,3)    == (/betrtracer_vars%nvolatile_tracer_groups/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bunsencef_col,1)   == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bunsencef_col,2)   == (/ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bunsencef_col,3)   == (/betrtracer_vars%nvolatile_tracer_groups/)),   errMsg(__FILE__,__LINE__))

   associate( &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups      , & !Integer[intent(in)], number of tracers
    tracer_group_memid         => betrtracer_vars%tracer_group_memid           , &
    is_volatile                => betrtracer_vars%is_volatile                  , & !logical[intent(in)], is a volatile tracer?
    is_h2o                     => betrtracer_vars%is_h2o                       , & !logical[intent(in)], is a h2o tracer
    volatilegroupid            => betrtracer_vars%volatilegroupid                & !integer[intent(in)], location in the volatile vector
   )


   if(any(is_h2o))then
     call calc_rhovap(bounds, lbj, ubj, jtops, numf, filter, t_soisno, smp_l, rho_vap)
   endif

   do j = 1, ngwmobile_tracer_groups
     !for tagged co2 simulations, the henry's constant are assumed same for all co2 tracers
     trcid = tracer_group_memid(j, 1)
     if(is_volatile(trcid))then
       k = volatilegroupid(trcid)

       do n = lbj, ubj
         do fc = 1, numf
           c = filter(fc)
           if(n>=jtops(c))then
             bunsencef_col(c,n, k)= henrycef_col(c,n,k)*t_soisno(c,n)/12.2_r8
             !add the pH effect for tracers that can exist in multiple aqueous phases
             if(is_h2o(trcid))then
               !for water isotopes
               bunsencef_col(c,n,j) = get_equi_lv_h2oiso_fractionation(trcid, t_soisno(c,j), betrtracer_vars) * denh2o/rho_vap(c,n)
             endif
           endif
         enddo
       enddo
     endif
   enddo
   end associate
   end subroutine calc_bunsen_coeff

!-------------------------------------------------------------------------------

   subroutine calc_dual_phase_convert_coeff(bounds, lbj, ubj, jtops, numf, filter, &
     waterstate_vars, betrtracer_vars, tracercoeff_vars)

   !DESCRIPTIONS:
   !compute phase conversion coefficients between gaseous and aqueous phases
   ! The total aqueous and gases phase concentration is = theta*aqueous+epsilon*gaseous
   ! because aqueous = bunsen*gaseous, these coefficients are constant throughout the all period.

   !USES:
   use BeTRTracerType     , only : betrtracer_type
   use WaterStateType     , only : Waterstate_Type
   use TracerCoeffType    , only : tracercoeff_type
   use clm_varcon         , only : denh2o, denice
   implicit none
   !arguments
   type(bounds_type),      intent(in)    :: bounds                      ! bounds
   integer,                intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
   integer,                intent(in)    :: jtops(bounds%begc: )        ! top label of each column
   integer,                intent(in)    :: numf                        ! number of columns in column filter
   integer,                intent(in)    :: filter(:)                   ! column filter
   type(betrtracer_type),  intent(in)    :: betrtracer_vars             ! betr configuration information
   type(Waterstate_Type),  intent(in)    :: waterstate_vars             ! water state variables
   type(tracercoeff_type), intent(inout) :: tracercoeff_vars            ! structure containing tracer transport parameters

   !local variables
   integer :: j, n, k, fc, c , trcid  ! indices
   character(len=255) :: subname = 'calc_dual_phase_convert_coeff'

   SHR_ASSERT_ALL((ubound(jtops)                   == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))

   associate( &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups      , & !Input: [integer(:)], number of tracers
    tracer_group_memid         => betrtracer_vars%tracer_group_memid           , & !Input: [integer(:)], tracer id
    is_h2o                     => betrtracer_vars%is_h2o                       , & !Input: [logical(:)], is h2o tracer?
    is_volatile                => betrtracer_vars%is_volatile                  , & !Input: [logical(:)], is a volatile tracer?
    volatilegroupid            => betrtracer_vars%volatilegroupid              , & !Input: [logical(:)], location in the volatile vector
    adsorbgroupid              => betrtracer_vars%adsorbgroupid                , & !Input: [Integer(:)], tracer id
    h2osoi_liqvol              => waterstate_vars%h2osoi_liqvol_col            , & !Input: [real(r8)(:,:)], liquid h2o vol
    h2osoi_icevol              => waterstate_vars%h2osoi_icevol_col            , & !Input: [real(r8)(:,:)], ice h2o vol
    air_vol                    => waterstate_vars%air_vol_col                  , & !Input: [real(r8)(:,:)], air vol
    bunsencef_col              => tracercoeff_vars%bunsencef_col               , & !Input: [real(r8)(:,:)], bunsen coeff
    aqu2bulkcef_mobile         => tracercoeff_vars%aqu2bulkcef_mobile_col      , & !Output:[real(r8)(:,:)], phase conversion coeff
    aqu2equilsolidcef          => tracercoeff_vars%aqu2equilsolidcef_col       , & !Input: [real(r8)(:,:)], phase conversion coeff
    gas2bulkcef_mobile         => tracercoeff_vars%gas2bulkcef_mobile_col        & !Output:[real(r8)(:,:)], phase conversion coeff
   )


  do j = 1, ngwmobile_tracer_groups
    trcid = tracer_group_memid(j,1)
    if(is_volatile(trcid))then
      k = volatilegroupid(trcid)
      do n = lbj, ubj
        do fc = 1, numf
          c = filter(fc)
          if(n>=jtops(c))then
            !aqueous to bulk mobile phase
            aqu2bulkcef_mobile(c,n,j) = air_vol(c,n)/bunsencef_col(c,n,k)+h2osoi_liqvol(c,n)

            !gaseous to bulk mobile phase
            gas2bulkcef_mobile(c,n,k) = air_vol(c,n)+h2osoi_liqvol(c,n)*bunsencef_col(c,n,k)

            if(is_h2o(trcid))then
              !for water tracer, I assume the three phases are in equilibrium, such that
              aqu2bulkcef_mobile(c,n,j)= aqu2bulkcef_mobile(c,n,j) + aqu2equilsolidcef(c,n,adsorbgroupid(trcid))

              gas2bulkcef_mobile(c,n,k) = gas2bulkcef_mobile(c,n,k)+ aqu2equilsolidcef(c,n,adsorbgroupid(trcid)) * bunsencef_col(c,n,k)
            endif

            !correct for impermeable layer, to avoid division by zero in doing diffusive transport
            gas2bulkcef_mobile(c,n,k) = max(gas2bulkcef_mobile(c,n,k),air_vol(c,n),minval_airvol)
          endif
        enddo
      enddo
    else
      !when linear adsorption is used for some adsorptive aqueous tracers, the aqu2bulkcef will be the retardation factor
      !for the moment, it is set to one for all non-volatile tracers
      !It is assumed that ice have same equilibrium solublity as liquid water for soluable tracers
      do n = lbj, ubj
        do fc = 1, numf
          c = filter(fc)
          if(n>=jtops(c))then
            aqu2bulkcef_mobile(c, n, j) = h2osoi_liqvol(c,n)+denice/denh2o * h2osoi_icevol(c,n)
          endif
        enddo
      enddo
    endif
  enddo
  end associate
  end subroutine calc_dual_phase_convert_coeff


!-------------------------------------------------------------------------------
   function get_equilibrium_scal(temp, pH, tracer, betrtracer_vars)result(rscal)

   !DESCRIPTION:
   !
   !obtain the equilibrium scaling factor for species that
   !can exist in multiple aqueous forms through hydrolysis.
   !obtain mole fractions stored as carbonate and bicarbonate
   !dlogK/dT=dH/(2.303RT^2)
   !logK=logK(T_0)+dH*(1/(2.303R*T_0)-1/(2.303R*T))
   !
   ! through the scale, the bulk aqueous tracer is, rscal * aqueous(netural)
   ! !USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   real(r8), intent(in) :: temp, pH
   integer,  intent(in) :: tracer
   type(betrtracer_type), intent(in) :: betrtracer_vars                                          ! betr configuration information

   !local variables
   real(r8), parameter :: Tref=298.15 ! Kelvin
   real(r8), parameter :: co2reflogK1=-6.352_r8   !25 Celcius
   real(r8), parameter :: co2reflogK2=-10.33_r8   !25 Celcisum
   real(r8), parameter :: co2dH1=-2.0e3_r8 ! J/mol
   real(r8), parameter :: co2dH2=-3.5e3_r8 ! J/mol
   real(r8), parameter :: nh3logK=9.24_r8
   real(r8), parameter :: no3logK=1.30_r8
   real(r8), parameter :: R=8.3144
   real(r8) :: co2logK1, co2logK2, rscal
   character(len=255) :: subname ='get_equilibrium_scal'


   if(tracer==betrtracer_vars%id_trc_co2x)then
      !H2CO3  <--> H(+)+HCO3(-)    K1
      !HCO3(-)<--> H(+)+CO3(2-)    K2
      !I have to check why I need 1.e3_r8 for conversion
      co2logK1 = co2reflogK1+co2dH1*(1._r8/(2.303_r8*R*Tref)-1._r8/(2.303_r8*R*temp))
      co2logK2 = co2reflogK2+co2dH2*(1._r8/(2.303_r8*R*Tref)-1._r8/(2.303_r8*R*temp))
      rscal = 1._r8+10._r8**(co2logK1)*10._r8**(-pH)*(1._r8+10._r8**(co2logK2)*10._r8**(-pH))*1.e3_r8

   elseif(tracer==betrtracer_vars%id_trc_nh3x)then
      !NH3H2O <--> NH4(+) + OH(-)
      rscal = 1._r8+10._r8**(nh3logK)*10._r8**(-pH)

   elseif(tracer==betrtracer_vars%id_trc_no3x)then
      !HNO3 <--> NO3(-) + H(+)
      rscal = 1._r8+10._r8**(no3logK)*10._r8**(-pH)
   else
      rscal = 1._r8  ! no rescal for other tracers
   endif
   return
   end function get_equilibrium_scal
!-------------------------------------------------------------------------------
   function get_henrycef(temp, trcid, betrtracer_vars)result(henry)
   !
   ! !DESCRIPTION
   !Compute the henry's law coefficient
   !There are unconsidered isotopic effect on the henry's law constant.
   !Some theoretical comments on such effect can be found in Jancso, 2002
   !REFERENCE: compilation of henry's law constants for inorganic and organic
   ! species of potential importance in environmental chemistry, Rolf Sander

   !USES
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   real(r8), intent(in) :: temp
   integer, intent(in) :: trcid
   type(betrtracer_type), intent(in) :: betrtracer_vars                                          ! betr configuration information

   !local variable
   real(r8) :: henry    !unit[M/atm]=[mol_aq/dm3_aq]/[atm]
   real(r8) :: es
   character(len=255) :: subname ='get_henrycef'


   if(trcid == betrtracer_vars%id_trc_ar)then
      henry=1.4e-3_r8*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_o2)then
      henry=1.3e-3_r8*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid == betrtracer_vars%id_trc_nh3x)then
      henry=5.6e1_r8*exp(-4100._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_n2)then
      henry=6.1e-4_r8*exp(-1300._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_n2o)then
      henry=2.5e-2_r8*exp(-2600._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_no)then
      henry=1.9e-3*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_c13_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_c14_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_o17_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_o18_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_air_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid == betrtracer_vars%id_trc_ch4)then
      henry=1.3e-3_r8*exp(-1700._r8*(1._r8/temp-1._r8/298.15_r8))  !mol dm^-3 amt^-1
   endif
   return
   end function get_henrycef
!-------------------------------------------------------------------------------

   subroutine calc_rhovap(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, t_soisno, smp_l, rho_vap)
   !
   !DESCRIPTION
   !Compute actual vapor pressure inside the soil profile
   !
   ! uses
   use clm_varcon    , only : rwat, grav
   use QSatMod       , only : rhoSat

   implicit none
   !arguments
   type(bounds_type),      intent(in) :: bounds                           ! bounds
   integer,                intent(in) :: lbj, ubj                         ! lower and upper bounds, make sure they are > 0
   integer,                intent(in) :: jtops(bounds%begc: )             ! top label of each column
   integer,                intent(in) :: num_soilc                        ! number of soil filters
   integer,                intent(in) :: filter_soilc(:)                  ! filter
   real(r8),               intent(in) :: t_soisno(bounds%begc: , lbj: )   !soil temperature, K
   real(r8),               intent(in) :: smp_l(bounds%begc: , lbj: )      !liquid soil matric potential, mm
   real(r8),               intent(inout):: rho_vap(bounds%begc: , lbj: )  !actual vapor pressure, kg/m3

   !local variables
   real(r8) :: hh           !relative humidity
   real(r8) :: rho_sat      !saturated water vapor pressure
   integer :: c, fc, n
   character(len=255) :: subname ='calc_rhovap'

   SHR_ASSERT_ALL((ubound(jtops)        == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(t_soisno)     == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(smp_l)        == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(rho_vap)      == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))


   !be careful below, because snow and pure water has no definition of water matrix potential
   do n = lbj, ubj
     do fc = 1, num_soilc
       c = filter_soilc(fc)
       if(n>=jtops(c))then
         !calculate saturated vapor pressure
         if(n>=1)then
           hh = exp(smp_l(c,n)*1.e-3_r8*grav/(rwat*t_soisno(c,n)))          !relative humidity using Kelvin equation
         else
           hh = 1._r8                                                       !just in case it is snow layer
         endif
         call rhoSat(t_soisno(c,n), rho_sat)
         !I use max to avoid completely dry soil
         rho_vap(c,n) = rho_sat * max(hh,1.e-4_r8)                          !kg/m3
       endif
     enddo
   enddo
   end subroutine calc_rhovap

!-------------------------------------------------------------------------------
   function get_taugas(eff_por, airvol, bsw)result(taugas)
   !
   !Descriptions
   !compute the tortuosity for the gas diffusion
   !Reference:
   !Millington and Quirk, 1961; Maggi et al., 2008
   !Moldrup et al, 2003
   !
   ! USES
   !

   implicit none
   real(r8),           intent(in) :: eff_por   ! effective porosity
   real(r8),           intent(in) :: airvol    ! air filled volume
   real(r8), optional, intent(in) :: bsw       ! clapp-hornber shape parameter

   real(r8) :: taugas
   character(len=255) :: subname ='get_taugas'

   if(eff_por == 0._r8)then
      taugas = 0._r8
   else
      if(present(bsw))then
         !modified from Eq.(5) in Moldrup et al., 2003
         taugas = (airvol/eff_por)**(3._r8/bsw)*airvol
      else
         taugas= eff_por**(1._r8/3._r8)*(airvol/eff_por)**(7._r8/3._r8)
      endif
   endif
   end function get_taugas
!-------------------------------------------------------------------------------
   function get_tauliq(eff_por, liqvol, bsw)result(tauliq)
   !
   !DESCRIPTION:
   !compute tortuosity for solute diffusion
   !Reference:
   !Millington and Quirk, 1961; Maggi et al., 2008
   !Moldrup et al, 2003
   !
   ! USES
   !

   implicit none
   real(r8),           intent(in) :: eff_por   !effective porosity
   real(r8),           intent(in) :: liqvol    !liquid water filled volume
   real(r8), optional, intent(in) :: bsw       !clapp-hornberg shape parameter

   real(r8) :: tauliq
   character(len=255) :: subname ='get_tauliq'


   if(eff_por == 0._r8)then
      tauliq = 0._r8
   else
      if(present(bsw))then
         !Modified from Eq.(3) in Moldrup et al., 2003.
         tauliq=(min(liqvol/eff_por,1._r8))**(bsw/3._r8-1._r8)*liqvol
      else
         tauliq=eff_por**(1._r8/3._r8)*(min(liqvol/eff_por,1._r8))**(7._r8/3._r8)
      endif
   endif
   end function get_tauliq
!-------------------------------------------------------------------------------

   function get_gas_diffusivity(trcid, temp, betrtracer_vars)result(diff)
   !
   ! !DESCRIPTIONS
   !
   !compute gaseous diffusivity for volatile species
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   integer, intent(in)  :: trcid
   real(r8), intent(in) :: temp  !kelvin
   type(betrtracer_type), intent(in) :: betrtracer_vars                                          ! betr configuration information

   !local variable
   real(r8) :: diff
   character(len=255) :: subname = 'get_gas_diffusivity'

   if(trcid==betrtracer_vars%id_trc_n2)then
      diff=1.93e-5_r8*(temp/273.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_ar)then
      diff=1.61e-5_r8*(temp/273.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_o2)then
      diff=1.8e-5_r8*(temp/273.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_air_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_ch4)then
      diff=1.9e-5_r8*(temp/298.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_nh3x)then
      diff=0.211e-4_r8*(temp/273.15_r8)**1.75_r8
   elseif(trcid==betrtracer_vars%id_trc_no)then
      diff=0.199e-4_r8*(temp/293.15)**1.5_r8
   elseif(trcid==betrtracer_vars%id_trc_n2o)then
      diff=0.159e-4_r8*(temp/293.15)**1.5_r8
!isotopes
   !use the kinetic theory of gases, assuming the collison diameters are same
   !between the light and heavy isotopmer, assuming the bath gas is dry air
   elseif(trcid==betrtracer_vars%id_trc_o18_h2o)then
      diff = 0.9723_r8*0.226e-4_r8*(temp/273.15_r8)**1.75_r8     !from Merlivat, 1978
   elseif(trcid==betrtracer_vars%id_trc_d_h2o)then
      diff = 0.9755_r8*0.226e-4_r8*(temp/273.15_r8)**1.75_r8     !from Merlivat, 1978
   elseif(trcid==betrtracer_vars%id_trc_c13_co2x)then
      diff = 0.9958_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8      ! from Trudinger, 1997
   elseif(trcid==betrtracer_vars%id_trc_o18_co2x)then
      diff = 0.9913_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8      !e.g. Wingate et al., 2009
   elseif(trcid==betrtracer_vars%id_trc_c14_co2x)then
      diff = 0.9918_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8      ! from Trudinger, 1997
   elseif(trcid==betrtracer_vars%id_trc_o17_co2x)then
      diff = 0.9957_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   endif
   return
   end function get_gas_diffusivity
!-------------------------------------------------------------------------------

   function get_aqueous_diffusivity(trcid, temp, betrtracer_vars, is_dom)result(diff)
   !
   ! Descriptions:
   ! Compute aqueous diffusivity for volatile species
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   integer,   intent(in) :: trcid
   real(r8),  intent(in) :: temp
   type(betrtracer_type), intent(in) :: betrtracer_vars                                          ! betr configuration information
   logical, optional,  intent(in) :: is_dom

   !local variable
   real(r8) :: diff       ! m2/s
   character(len=255) :: subname ='get_aqueous_diffusivity'

   if(present(is_dom))then
     if(is_dom)then
       diff=2.6e-9_r8
       return
     endif
   endif
   if(trcid==betrtracer_vars%id_trc_n2)then
      diff=2.57e-9_r8*(temp/273.0_r8)
   elseif(trcid==betrtracer_vars%id_trc_o2)then
      diff=2.4e-9_r8*temp/298.0_r8
   elseif(trcid==betrtracer_vars%id_trc_ar)then
      diff=2.15e-9_r8*temp/298.0_r8
   elseif(trcid==betrtracer_vars%id_trc_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_air_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_ch4)then
      diff=1.5e-9_r8*temp/298.0_r8
   elseif(trcid==betrtracer_vars%id_trc_no3x)then
      diff=2.6e-9_r8*temp/298.15_r8
   elseif(trcid==betrtracer_vars%id_trc_no2x)then
      diff=2.6e-9_r8*temp/298.15_r8  !considers revision
   elseif(trcid==betrtracer_vars%id_trc_nh3x)then
      diff=1.64e-5_r8*temp/298.15_r8

!isotopes
   elseif(trcid==betrtracer_vars%id_trc_d_h2o)then
      diff=0.9833_r8*1e-9_r8*exp(-(535400._r8/temp-1393.3_r8)/temp+2.1876_r8)
   elseif(trcid==betrtracer_vars%id_trc_o18_h2o)then
      diff=0.9669_r8*1e-9_r8*exp(-(535400._r8/temp-1393.3_r8)/temp+2.1876_r8)
   elseif(trcid==betrtracer_vars%id_trc_o18_co2x)then
      !theoretical calculations based on molecular dynamics indicate the fractionation between carbonate and bicarbonate due to diffusion is less than 1 per mil.
      !so set it to the diffusivity of the base value
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_c13_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_c14_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   else
      diff=2.6e-9_r8
   endif

   end function get_aqueous_diffusivity


!-------------------------------------------------------------------------------
   function get_diffusivity_ratio_gas2h2o(trcid, temp, betrtracer_vars)result(ratio)
   !
   ! DESCRIPTIONS:
   ! Compute the ratio of gas phase diffusivities for different volatile species in air with respect to that of water vapor
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   integer,  intent(in) :: trcid
   real(r8), intent(in) :: temp
   type(betrtracer_type), intent(in) :: betrtracer_vars                                          ! betr configuration information

   real(r8) :: diff, ratio, diffh2o
   character(len=255) :: subname = 'get_diffusivity_ratio_gas2h2o'


   diffh2o=0.229e-4_r8*(temp/273.15_r8)**1.75_r8

   if(trcid==betrtracer_vars%id_trc_n2)then
      ratio=1.93e-5_r8*(temp/273.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_o2)then
      ratio=1.8e-5_r8*(temp/273.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_ar)then
      ratio=1.6e-5_r8*(temp/273.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_air_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_ch4)then
      ratio=1.9e-5_r8*(temp/298.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_nh3x)then
      ratio=0.211e-4_r8*(temp/273.15_r8)**1.75_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_no)then
      ratio=0.199e-4_r8*(temp/293.15)**1.5_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_n2o)then
      ratio=0.159e-4_r8*(temp/293.15)**1.5_r8/diffh2o

!isotopes
   elseif(trcid==betrtracer_vars%id_trc_o18_h2o)then
      ratio = 0.9723_r8            !from Merlivat, 1978
   elseif(trcid==betrtracer_vars%id_trc_d_h2o)then
      ratio = 0.9755_r8            !from Merlivat, 1978
   elseif(trcid==betrtracer_vars%id_trc_c13_co2x)then
      ratio = 0.9957_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_c14_co2x)then
      ratio = 0.9913_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_o18_co2x)then
      ratio = 0.9913_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_o17_co2x)then
      ratio = 0.9957_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   endif
   return
   end function get_diffusivity_ratio_gas2h2o

!-------------------------------------------------------------------------------
   subroutine convert_mobile2gas(bounds, lbj, ubj, jtops, numf, filter, do_forward, gas2bulkcef_mobile_col, betrtracer_vars, tracer_conc_mobile)
   !
   ! DESCRIPTIONS
   ! do conversion between bulk mobile phase and gaseous phase
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   !arguments
   type(bounds_type),      intent(in) :: bounds                  ! bounds
   integer,                intent(in) :: lbj, ubj                ! lower and upper bounds
   integer,                intent(in) :: jtops(bounds%begc: )    ! top label of each column
   integer,                intent(in) :: numf                    ! number of filters
   integer,                intent(in) :: filter(:)               ! filter
   logical,                intent(in) :: do_forward              ! true, dual_bulk => gaseous
   type(betrtracer_type),  intent(in) :: betrtracer_vars         ! betr configuration information
   real(r8),               intent(in) :: gas2bulkcef_mobile_col(bounds%begc: ,lbj: , 1: )  !conversion parameter
   real(r8),             intent(inout):: tracer_conc_mobile(bounds%begc: ,lbj: , 1: )  !bulk mobile tracer

   !local variables
   integer :: jj, kk, fc, c, j
   character(len=255) :: subname = 'convert_mobile2gas'

   SHR_ASSERT_ALL((ubound(jtops)                  == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(gas2bulkcef_mobile_col,1) == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(gas2bulkcef_mobile_col,2) == (/ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(gas2bulkcef_mobile_col,3) == (/betrtracer_vars%nvolatile_tracers/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(tracer_conc_mobile,1)     == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(tracer_conc_mobile,2)     == (/ubj/)),   errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(tracer_conc_mobile,3)     == (/betrtracer_vars%nvolatile_tracers/)),   errMsg(__FILE__,__LINE__))


   associate( &
    ngwmobile_tracers    => betrtracer_vars%ngwmobile_tracers            , & !Integer[intent(in)], number of tracers
    is_volatile => betrtracer_vars%is_volatile         , & !logical[intent(in)], is a volatile tracer?
    volatilegroupid  => betrtracer_vars%volatilegroupid            & !integer[intent(in)], location in the volatile vector
   )

   do jj = 1, ngwmobile_tracers
     if(is_volatile(jj))then
       kk = volatilegroupid(jj)
       if(do_forward)then
         do j = lbj, ubj
           do fc = 1, numf
             c = filter(fc)
             if(j>=jtops(c))then
               tracer_conc_mobile(c,j,jj) = tracer_conc_mobile(c,j,jj) / gas2bulkcef_mobile_col(c,j,kk)
             endif
           enddo
         enddo
       else
         do j = lbj, ubj
           do fc = 1, numf
             c = filter(fc)
             if(j>=jtops(c))then
               tracer_conc_mobile(c,j,jj) = tracer_conc_mobile(c,j,jj) * gas2bulkcef_mobile_col(c,j,kk)
             endif
           enddo
         enddo
       endif
     endif
   enddo
   end associate
   end subroutine convert_mobile2gas
!-------------------------------------------------------------------------------

   subroutine set_multi_phase_diffusion(bounds, lbj, ubj, jtops, numf, filter, soilstate_vars, waterstate_vars, &
      canopystate_vars, temperature_vars, chemstate_vars, betrtracer_vars, tracercoeff_vars)
   !
   ! DESCRIPTION
   ! set parameters for the dual phase diffusion
   !
   use TracerCoeffType    , only : tracercoeff_type
   use WaterStateType     , only : Waterstate_Type
   use SoilStateType      , only : soilstate_type
   use TemperatureType    , only : temperature_type
   use ChemStateType      , only : chemstate_type
   use BeTRTracerType     , only : betrtracer_type
   use CanopyStateType    , only : canopystate_type

   type(bounds_type)       , intent(in) :: bounds  ! bounds
   integer                 , intent(in) :: lbj, ubj             ! lower and upper bounds, make sure they are > 0
   integer                 , intent(in) :: jtops(bounds%begc: ) ! top label of each column
   integer                 , intent(in) :: numf                 ! number of columns in column filter
   integer                 , intent(in) :: filter(:)            ! column filter
   type(betrtracer_type)   , intent(in) :: betrtracer_vars      ! betr configuration information
   type(Waterstate_Type)   , intent(in) :: waterstate_vars          ! water state variables
   type(soilstate_type)    , intent(in) :: soilstate_vars          ! physical state variables
   type(temperature_type)  , intent(in) :: temperature_vars          ! energy state variable
   type(canopystate_type)  , intent(in) :: canopystate_vars
   type(chemstate_type)    , intent(in) :: chemstate_vars        ! chemistry state variable
   type(tracercoeff_type)  , intent(inout) :: tracercoeff_vars ! structure containing tracer transport parameters

   !
   real(r8) :: bulkdiffus(bounds%begc:bounds%endc,lbj:ubj,1:betrtracer_vars%ntracer_groups )  !weighted bulk diffusivity for dual-phase diffusion

   !maybe I should use tau_soil as a local variable, I will check this later
   character(len=255) :: subname='set_multi_phase_diffusion'


   SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))

   !compute tortuosity
   !gaseous phase
   call calc_gaseous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, numf, filter, soilstate_vars, waterstate_vars, &
      tau_soil%tau_gas(bounds%begc:bounds%endc,lbj:ubj))

   !aqueous phase
   call calc_aqueous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, numf, filter, soilstate_vars, waterstate_vars, &
      tau_soil%tau_liq(bounds%begc:bounds%endc,lbj:ubj))

   !compute bulk diffusivity
   call calc_bulk_diffusivity(bounds, lbj, ubj, jtops, numf, filter       , &
      tracercoeff_vars%bunsencef_col(bounds%begc:bounds%endc,lbj:ubj, : ) , &
      canopystate_vars, waterstate_vars, tau_soil, betrtracer_vars        , &
      temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj), bulkdiffus)

   !compute weigthed conductances
   call calc_bulk_conductances(bounds, lbj, ubj, jtops, numf, filter, bulkdiffus, &
      col%dz(bounds%begc:bounds%endc,lbj:ubj), betrtracer_vars,  &
      tracercoeff_vars%hmconductance_col(bounds%begc:bounds%endc, lbj:ubj-1, :))

   end subroutine set_multi_phase_diffusion


!--------------------------------------------------------------------------------
   subroutine set_phase_convert_coeff(bounds, lbj, ubj, jtops, numf, filter, dz, soilstate_vars, waterstate_vars, &
      temperature_vars, chemstate_vars, betrtracer_vars, tracercoeff_vars)
   !
   ! DESCRIPTION
   ! set parameters for phase conversion
   use TracerCoeffType    , only : tracercoeff_type
   use WaterStateType     , only : Waterstate_Type
   use SoilStateType      , only : soilstate_type
   use TemperatureType    , only : temperature_type
   use ChemStateType      , only : chemstate_type
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   type(bounds_type),      intent(in) :: bounds  ! bounds
   integer,                intent(in) :: lbj, ubj             ! lower and upper bounds, make sure they are > 0
   integer,                intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
   integer,                intent(in) :: numf                 ! number of columns in column filter
   integer,                intent(in) :: filter(:)            ! column filter
   real(r8),               intent(in) :: dz(bounds%begc: ,lbj: )
   type(betrtracer_type),  intent(in) :: betrtracer_vars             ! betr configuration information
   type(Waterstate_Type),  intent(in) :: waterstate_vars          ! water state variables
   type(soilstate_type),   intent(in) :: soilstate_vars          ! physical state variables
   type(temperature_type), intent(in) :: temperature_vars          ! energy state variable
   type(chemstate_type),   intent(in) :: chemstate_vars        ! chemical state variable
   type(tracercoeff_type), intent(inout) :: tracercoeff_vars ! structure containing tracer transport parameters
   character(len=255) :: subname = 'set_phase_convert_coeff'

   SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(dz)              == (/bounds%endc, ubj/)),        errMsg(__FILE__,__LINE__))

   !compute Henry's law constant
   call calc_henrys_coeff(bounds, lbj, ubj, jtops, numf, filter                    , &
       temperature_vars%t_soisno_col(bounds%begc:bounds%endc,lbj:ubj)              , &
       chemstate_vars%soil_pH(bounds%begc:bounds%endc, lbj:ubj),  betrtracer_vars  , &
       tracercoeff_vars%aqu2neutralcef_col(bounds%begc:bounds%endc,lbj:ubj, : )    , &
       tracercoeff_vars%henrycef_col(bounds%begc:bounds%endc, lbj:ubj, : ))

   !compute Bunsen's coefficients
   call calc_bunsen_coeff(bounds, lbj, ubj, jtops, numf, filter                    , &
        tracercoeff_vars%henrycef_col(bounds%begc:bounds%endc, lbj:ubj, : )        , &
        temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj)            , &
        waterstate_vars%smp_l_col    (bounds%begc:bounds%endc, lbj:ubj)            , &
        betrtracer_vars                                                            , &
        tracercoeff_vars%bunsencef_col(bounds%begc:bounds%endc, lbj:ubj, :))

   !compute equilibrium fraction to liquid phase conversion parameter
   if(betrtracer_vars%nsolid_equil_tracers>0)then
     call calc_equil_to_liquid_convert_coeff(bounds, lbj, ubj, jtops, numf, filter , &
        temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj)            , &
        waterstate_vars%h2osoi_ice_col(bounds%begc:bounds%endc,lbj:ubj)            , &
        dz(bounds%begc:bounds%endc, lbj:ubj)                                   , &
        betrtracer_vars, tracercoeff_vars%aqu2equilsolidcef_col(bounds%begc:bounds%endc, lbj:ubj,:))
   endif

   !compute phase conversion coefficients
   call calc_dual_phase_convert_coeff(bounds, lbj, ubj, jtops, numf, filter, waterstate_vars, betrtracer_vars, tracercoeff_vars)


   end subroutine set_phase_convert_coeff



  !------------------------------------------------------------------------
   subroutine calc_tracer_infiltration(bounds, lbj, ubj, jtops, numf, filter, bunsencef_topsoi, &
     betrtracer_vars, tracerboundarycond_vars, waterflux_vars, tracer_flx_infl)
   !
   ! DESCRIPTION
   ! calculate advection velocity for BeTR code
   ! this assumes the interfacial velocity qflx_adv (except infiltration) has been calcualted in doing vertical
   ! watermovement
   ! This assumes the advection solves the equation
   ! \frac{\parital m}{\partial t}+\frac{V*m}{\partial z} = 0
   !where m = C*vsm, therefore, V=ql/vsm
   !
   use WaterfluxType         , only : waterflux_type
   use TracerBoundaryCondType, only : tracerboundarycond_type
   use BeTRTracerType        , only : betrtracer_type
   use clm_varcon            , only : denh2o
   implicit none

   type(bounds_type)            , intent(in)    :: bounds
   integer                      , intent(in)    :: lbj, ubj
   integer                      , intent(in)    :: jtops(bounds%begc: )           ! top label of each column
   integer                      , intent(in)    :: numf                           ! number of columns in column filter
   integer                      , intent(in)    :: filter(:)                      ! column filter
   real(r8)                     , intent(in)    :: bunsencef_topsoi(bounds%begc: , 1: )
   type(BeTRTracer_Type)        , intent(in)    :: betrtracer_vars
   type(tracerboundarycond_type), intent(in)    :: tracerboundarycond_vars
   type(waterflux_type)         , intent(inout) :: waterflux_vars
   real(r8)                     , intent(inout) :: tracer_flx_infl(bounds%begc: , 1: )

   ! local variables
   integer :: fc, c, j

   associate(                                                                          &
    qflx_gross_infl_soil => waterflux_vars%qflx_gross_infl_soil_col                    & !real(r8) (:)  [intent(in)], infiltration, mm/s
   )

   SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bunsencef_topsoi) == (/bounds%endc, betrtracer_vars%nvolatile_tracers/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(tracer_flx_infl) == (/bounds%endc, betrtracer_vars%ngwmobile_tracers/)), errMsg(__FILE__,__LINE__))

   do j = 1, betrtracer_vars%ngwmobile_tracers

     !for a real mechanistic modeling, tracer_flx_infl should be derived from water flux coming from snow melt, surface ponding water,
     !and precipitation. I here just comparomise for a quick shot.

     if(j==betrtracer_vars%id_trc_o18_h2o)then
       do fc = 1, numf
         c = filter(fc)
         tracer_flx_infl(c,j) = qflx_gross_infl_soil(c)/denh2o    !kg m-2 s-1/ kg m-3 = m/s
       enddo
     else
       do fc = 1, numf
         c = filter(fc)

         if(betrtracer_vars%is_volatile(j) .and. betrtracer_vars%is_advective(j))then
           !for volatile non water tracer, infiltration is calculated based dissolution of the gas in the water, this may need
           !improvement when tracers are allowed to transport inside snow, such that the tracer infiltration is derived from mass balance in snow
           tracer_flx_infl(c,j) = bunsencef_topsoi(c,betrtracer_vars%volatilegroupid(j)) * tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1,j) * qflx_gross_infl_soil(c)*1.e-3_r8
         else
           tracer_flx_infl(c,j) = 0._r8
         endif
       enddo
     endif
   enddo
   end associate
   end subroutine calc_tracer_infiltration
   !------------------------------------------------------------------------
   subroutine pre_diagnose_soilcol_water_flux(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, h2osoi_liq)
   !
   ! DESCRIPTION
   ! pre diagnose advective water fluxes at different soil interfaces


   implicit none
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
   integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
   real(r8)                , intent(in)    :: h2osoi_liq(bounds%begc: , 1: )

   !local variables
   integer :: j, fc, c
   SHR_ASSERT_ALL((ubound(h2osoi_liq) == (/bounds%endc, nlevsoi/)), errMsg(__FILE__,__LINE__))

   do j = 1, nlevsoi
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       h2osoi_liq_copy(c,j) = h2osoi_liq(c,j)
     enddo
   enddo
   end subroutine pre_diagnose_soilcol_water_flux

   !------------------------------------------------------------------------
   subroutine diagnose_advect_water_flux(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, h2osoi_liq, qcharge, waterflux_vars)
   !
   ! DESCRIPTION
   ! diagnose advective water fluxes between different soil layers
   !

   use WaterFluxType        , only : waterflux_type
   use clm_time_manager     , only : get_step_size
   use clm_varcon           , only : denh2o
   implicit none
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
   integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
   type(waterflux_type)    , intent(inout) :: waterflux_vars
   real(r8)                , intent(in)    :: h2osoi_liq(bounds%begc: , 1: )  !mm H2O/m2 eqv. kg H2O/m2
   real(r8)                , intent(in)    :: qcharge(bounds%begc: )  ! mm H2O/s aquifer recharge rate

   !local variables
   integer :: j, fc, c
   real(r8):: dtime
   real(r8):: diff
   real(r8):: infl_tmp
   real(r8):: scal

   SHR_ASSERT_ALL((ubound(h2osoi_liq) == (/bounds%endc, nlevsoi/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(qcharge)    == (/bounds%endc/))         , errMsg(__FILE__,__LINE__))


   associate(                                                             & !
     qflx_rootsoi        =>    waterflux_vars%qflx_rootsoi_col          , & ! Iput  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/s) (+ = to atm)
     qflx_adv            =>    waterflux_vars%qflx_adv_col              , & ! Output: [real(r8) (:,:) ]  water flux at interfaces       (m H2O/s) (- = to atm)
     qflx_gross_infl_soil=>    waterflux_vars%qflx_gross_infl_soil_col  , & ! Output: [real(r8) (:)] gross infiltration (mm H2O/s)
     qflx_infl           =>    waterflux_vars%qflx_infl_col             , & ! Output: [real(r8) (:)] infiltration
     qflx_gross_evap_soil=>    waterflux_vars%qflx_gross_evap_soil_col    & ! Output: [real(r8) (:)] gross evaporation (mm H2O/s)
   )

   ! get time step
   dtime = get_step_size()
   !start from the bottom layer, because the water exchange between vadose zone soil and aquifer and plant root is known
   !the water flux at uppper surface can be inferred using the mass balance approach
   do j = nlevsoi, 1, -1
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(j==nlevsoi)then
         qflx_adv(c,j) = qcharge(c) * 1.e-3_r8
       else
         qflx_adv(c,j) = 1.e-3_r8 * (h2osoi_liq(c,j+1)-h2osoi_liq_copy(c,j+1))/dtime + qflx_adv(c,j+1) + qflx_rootsoi(c,j+1)
       endif

     enddo
   enddo

   ! correct gross infiltration and gross evaporation
   ! (h2osoi_liq(c,1)-h2osoi_liq_copy(c,1))/dtime=qflx_infl-q_out-qflx_rootsoi
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)

     !obtain the corrected infiltration
     qflx_infl(c) = (h2osoi_liq(c,1)-h2osoi_liq_copy(c,1))/dtime + (qflx_rootsoi(c,1)+qflx_adv(c,1))*1.e3_r8
     !the predicted net infiltration
     infl_tmp=qflx_gross_infl_soil(c)-qflx_gross_evap_soil(c)
     diff=qflx_infl(c)-infl_tmp
     if(abs(diff)/=0._r8)then
       if(infl_tmp==0._r8)then
         if(diff>0._r8)then
           qflx_gross_infl_soil(c)=diff
           qflx_gross_evap_soil(c)=0._r8
         else
           qflx_gross_infl_soil(c)=0._r8
           qflx_gross_evap_soil(c)=-diff
         endif
       else
         scal = (1._r8+diff/infl_tmp)
         qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c) * scal
         qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c) * scal
         if(qflx_gross_evap_soil(c)<0._r8)then
           !no negative evaporation allowed
           qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)-qflx_gross_evap_soil(c)
           qflx_gross_evap_soil(c) = 0._r8
         endif
         if(qflx_gross_infl_soil(c)<0._r8)then
           qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c)-qflx_gross_infl_soil(c)
           qflx_gross_infl_soil(c) = 0._r8
         endif
       endif
     endif

     qflx_adv(c,0) = qflx_gross_infl_soil(c) *.1e-3_r8  !surface infiltration

   enddo

   end associate
   end subroutine diagnose_advect_water_flux


   !------------------------------------------------------------------------
   subroutine diagnose_drainage_water_flux(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, h2osoi_liq,  waterflux_vars)
   !
   ! DESCRIPTION
   ! diagnose advective water fluxes between different soil layers
   !

   use WaterFluxType        , only : waterflux_type
   use clm_varcon           , only : denh2o
   use clm_time_manager     , only : get_step_size
   implicit none
   type(bounds_type)       , intent(in)    :: bounds               ! bounds
   integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
   integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
   type(waterflux_type)    , intent(inout) :: waterflux_vars
   real(r8)                , intent(in)    :: h2osoi_liq(bounds%begc: , 1: )  !mm H2O/m2 eqv. kg H2O/m2

   !local variables
   integer :: j, fc, c
   real(r8):: dtime

   SHR_ASSERT_ALL((ubound(h2osoi_liq) == (/bounds%endc, nlevsoi/)), errMsg(__FILE__,__LINE__))


   associate(                                                           & !
     qflx_drain_vr        =>    waterflux_vars%qflx_drain_vr_col        & ! Output  : [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/step) (to river +)
   )

   ! get time step
   dtime = get_step_size()
   !start from the bottom layer, because the water exchange between vadose zone soil and aquifer and plant root is known
   !the water flux at uppper surface can be inferred using the mass balance approach
   do j = nlevsoi, 1, -1
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qflx_drain_vr(c,j) = h2osoi_liq_copy(c,j)-h2osoi_liq(c,j)

       qflx_drain_vr(c,j) = min(qflx_drain_vr(c,j),h2osoi_liq_copy(c,j))/denh2o

     enddo
   enddo

   end associate
   end subroutine diagnose_drainage_water_flux

   !------------------------------------------------------------------------
   function get_equi_lv_h2oiso_fractionation(trcid, temp, betrtracer_vars)result(ans)
   !
   ! DESCRIPTION
   ! get equilibrium isotopic fractionation of liquid against gaseous phase
   !
   !
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   integer                   , intent(in) :: trcid
   real(r8)                  , intent(in) :: temp                !temperature
   type(betrtracer_type)     , intent(in) :: betrtracer_vars

   real(r8) :: ans



   !now it is set to one, for O18 and H/D, pleasure refer to Braud et al. (2005, J. Hydrology)
   ans = 1._r8
   return
   end function get_equi_lv_h2oiso_fractionation

   !------------------------------------------------------------------------
   function get_equi_sv_h2oiso_fractionation(trcid, temp, betrtracer_vars)result(ans)
   !
   ! DESCRIPTION
   ! get equilibrium isotopic fractionation of ice against vapor phase
   !
   !
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   integer                   , intent(in) :: trcid
   real(r8)                  , intent(in) :: temp                !temperature
   type(betrtracer_type)     , intent(in) :: betrtracer_vars

   real(r8) :: ans



   !now it is set to one, for O18, Roche (2013, GMD) gives some information
   ans = 1._r8
   return
   end function get_equi_sv_h2oiso_fractionation


   !------------------------------------------------------------------------
   function get_equi_sl_h2oiso_fractionation(trcid, temp, betrtracer_vars)result(ans)
   !
   ! DESCRIPTION
   ! get equilibrium isotopic fractionation of ice against liquid water
   !
   !
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   integer                   , intent(in) :: trcid
   real(r8)                  , intent(in) :: temp                !temperature
   type(betrtracer_type)     , intent(in) :: betrtracer_vars

   real(r8) :: ans



   !now it is set to one, it is equal to alpha_sv*alpha_vl
   ans = 1._r8
   return
   end function get_equi_sl_h2oiso_fractionation


   !------------------------------------------------------------------------
   subroutine calc_equil_to_liquid_convert_coeff(bounds, lbj, ubj, jtops, numf, filter, t_soisno, h2osoi_ice, dz, &
       betrtracer_vars, aqu2equilsolidcef_col)
   !
   ! DESCRIPTION
   ! calculate partition parameter between solid and aqueous phase tracers
   ! this could mean differnt things for different cases
   ! for water isotopes, this represents ice/liquid equilibrium partitioning, it could also
   ! mean linear isotherms of adsorption/desorption. Currently, it is only
   ! for water isotope partitioning between liquid water and ice
   !
   use clm_varcon            , only : denh2o, denice
   use BeTRTracerType        , only : betrtracer_type
   implicit none
   type(bounds_type)         , intent(in) :: bounds  ! bounds
   integer                   , intent(in) :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
   integer                   , intent(in) :: jtops(bounds%begc: )                              ! top label of each column
   integer                   , intent(in) :: numf                                              ! number of columns in column filter
   integer                   , intent(in) :: filter(:)                                         ! column filter
   real(r8)                  , intent(in) :: t_soisno(bounds%begc: , lbj: )
   type(betrtracer_type)     , intent(in) :: betrtracer_vars
   real(r8)                  , intent(in) :: h2osoi_ice(bounds%begc: , lbj: )
   real(r8)                  , intent(in) :: dz(bounds%begc: , lbj: )
   real(r8)                  , intent(inout) :: aqu2equilsolidcef_col(bounds%begc:bounds%endc, lbj:ubj, 1:betrtracer_vars%nsolid_equil_tracer_groups)

   real(r8) :: alpha_sl
   integer  :: fc, c, j

   SHR_ASSERT_ALL((ubound(t_soisno)          == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(h2osoi_ice)        == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(dz)                == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
   associate(                                           &
     is_h2o => betrtracer_vars%is_h2o                   &
   )

   if(any(is_h2o))then
     !doing a water isotope simulation
     if(betrtracer_vars%id_trc_o18_h2o>0)then
        do j = lbj, ubj
          do fc = 1, numf
            c = filter(fc)
            if(j>=jtops(c))then
              alpha_sl = get_equi_sl_h2oiso_fractionation(betrtracer_vars%id_trc_o18_h2o, t_soisno(c,j), betrtracer_vars)
              aqu2equilsolidcef_col(c,j, betrtracer_vars%id_trc_o18_h2o_ice) = alpha_sl * h2osoi_ice(c,j) / (denh2o * dz(c,j))
            endif
          enddo
        enddo
     endif
   endif

   end associate
   end subroutine calc_equil_to_liquid_convert_coeff

   !------------------------------------------------------------------------
   subroutine calc_smp_l(bounds, lbj, ubj, numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)

   use SoilStateType              , only : soilstate_type
   use WaterStateType             , only : waterstate_type
   use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
   use clm_varcon                 , only : grav,hfus,tfrz

   implicit none
   type(bounds_type)         , intent(in)    :: bounds  ! bounds
   integer                   , intent(in)    :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
   integer                   , intent(in)    :: numf                                              ! number of columns in column filter
   integer                   , intent(in)    :: filter(:)                                         ! column filter
   real(r8)                  , intent(in)    :: t_soisno(bounds%begc: , lbj: )                    ! soil temperature
   type(soilstate_type)      , intent(in)    :: soilstate_vars
   type(waterstate_type)     , intent(inout) :: waterstate_vars
   class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

   !local variables
   real(r8) :: s_node
   integer  :: fc, c, j

   SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(__FILE__,__LINE__))
   associate(                                                     & !
     h2osoi_vol        =>    waterstate_vars%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
     smp_l             =>    waterstate_vars%smp_l_col          , & ! Output: [real(r8) (:,:) ]  soil suction (mm)
     bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
     watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
     sucsat            =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
   )


   do j = lbj, ubj
     do fc = 1, numf
       c = filter(fc)
       if(j>=1)then

         if(t_soisno(c,j)<tfrz)then
           smp_l(c,j)= hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
         else
           s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
           call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j))
         endif

       endif
     enddo
   enddo
   end associate
   end subroutine calc_smp_l


  !-----------------------------------------------------------------------
  subroutine get_zwt (bounds, numf, filter, zi, &
       soilstate_vars, waterstate_vars, temperature_vars, zwt,jwt)
    !
    ! !DESCRIPTION:
    ! Finds the first unsaturated layer going up. Also allows a perched water table over ice.
    !
    use clm_varpar         , only : nlevsoi
    use clm_varcon         , only : tfrz
    use SoilStateType      , only : soilstate_type
    use WaterstateType     , only : waterstate_type
    use TemperatureType    , only : temperature_type
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)   :: bounds
    integer                , intent(in)   :: numf                ! number of column soil points in column filter
    integer                , intent(in)   :: filter(:)           ! column filter for soil points
    type(soilstate_type)   , intent(in)   :: soilstate_vars
    type(waterstate_type)  , intent(in)   :: waterstate_vars
    type(temperature_type) , intent(in)   :: temperature_vars
    real(r8)               , intent(in)   :: zi(bounds%begc: , 0: )
    real(r8)               , intent(inout):: zwt( bounds%begc: ) ! water table depth (-) [col]
    integer                , intent(inout):: jwt(bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: f_sat    ! volumetric soil water defining top of water table or where production is allowed
    integer  :: c,j,perch! indices
    integer  :: fc       ! filter column index

    !-----------------------------------------------------------------------
    f_sat = 0.95_r8   !a number borrowed from zack's ch4 code

    SHR_ASSERT_ALL((ubound(zwt) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(zi) == (/bounds%endc, nlevsoi/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(jwt) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))


    associate(                                          &
         watsat     => soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)  ] volumetric soil water at saturation (porosity)
         h2osoi_vol => waterstate_vars%h2osoi_vol_col , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         t_soisno   => temperature_vars%t_soisno_col    & ! Input:  [real(r8) (: ,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
         )


      ! The layer index of the first unsaturated layer, i.e., the layer right above
      ! the water table.
      ! ZS: Loop is currently not vectorized.
      do fc = 1, numf
         c = filter(fc)

         ! Check to see if any soil layers are frozen and saturated.  If so, start looking at the first layer above the top
         ! such layer.  This is potentially important for perched water tables in the Tundra.

         perch = nlevsoi
         do j = nlevsoi, 1, -1
            if (t_soisno(c,j) < tfrz .and. h2osoi_vol(c,j) > f_sat * watsat(c,j)) then
               ! strictly less than freezing because it could be permeable otherwise
               perch = j-1
            end if
         end do
         jwt(c) = perch

         do j = perch, 2, -1
            if(h2osoi_vol(c,j) > f_sat * watsat(c,j) .and. h2osoi_vol(c,j-1) < f_sat * watsat(c,j-1)) then
               jwt(c) = j-1
               zwt(c) = zi(c,jwt(c))
               exit
            end if
         enddo
         if (jwt(c) == perch .and. h2osoi_vol(c,1) > f_sat * watsat(c,1)) then ! missed that the top layer is saturated
            jwt(c) = 0
         endif

         zwt(c) = zi(c,jwt(c))
      end do

    end associate

  end subroutine get_zwt

  !-----------------------------------------------------------------------
  subroutine calc_aerecond(bounds, num_soilp, filter_soilp, jwt, rootfr, temperature_vars, betrtracer_vars, &
     canopystate_vars, carbonstate_vars, carbonflux_vars, tracercoeff_vars)
  !
  ! DESCRIPTION
  !
  ! calculate aerenchyma conductance (m/s)
  use clm_varcon            , only : tfrz, rpi
  use pftvarcon             , only : nc3_arctic_grass, crop, nc3_nonarctic_grass, nc4_grass, noveg
  use CNCarbonFluxType      , only : carbonflux_type
  use CNCarbonStateType     , only : carbonstate_type
  use CanopyStateType       , only : canopystate_type
  use BetrTracerType        , only : betrtracer_type
  use tracercoeffType       , only : tracercoeff_type
  use clm_varpar            , only : nlevsoi
  use TemperatureType       , only : temperature_type
  use MathfuncMod           , only : safe_div
  use clm_varctl            , only : use_cn
  use clm_time_manager      , only : get_step_size, get_nstep
  type(bounds_type)            , intent(in)   :: bounds
  integer                      , intent(in)   :: num_soilp                 ! number of column soil points in column filter
  integer                      , intent(in)   :: filter_soilp(:)           ! column filter for soil points
  integer                      , intent(in)   :: jwt(bounds%begc: )
  real(r8)                     , intent(in)   :: rootfr(bounds%begp: ,1: ) ! fraction of roots in each soil layer
  type(temperature_type)       , intent(in)   :: temperature_vars          ! energy state variable
  type(canopystate_type)       , intent(in)   :: canopystate_vars
  type(carbonstate_type)       , intent(in)   :: carbonstate_vars
  type(carbonflux_type)        , intent(in)   :: carbonflux_vars
  type(betrtracer_type)        , intent(in)   :: betrtracer_vars            ! betr configuration information
  type(tracercoeff_type)       , intent(inout) :: tracercoeff_vars

  real(r8) :: aerecond               !
  real(r8) :: nppratio
  real(r8) :: anpp
  real(r8) :: m_tiller
  real(r8) :: n_tiller
  real(r8) :: poros_tiller
  real(r8) :: area_tiller
  real(r8) :: lbl_rsc
  real(r8) :: porosmin = 0.05_r8                                     ! wait to be read in later
  real(r8) :: rob = 3._r8                                            ! ratio of root length to vertical depth ("obliquity"), wait to be read in later
  real(r8) :: nongrassporosratio = 0.33_r8                           ! non grass ratio
  real(r8) :: unsat_aere_ratio= 0.05_r8 / 0.3_r8
  logical  :: usefrootc = .false.                                    ! wait to be read in later
  integer  :: j, fp, p, c, g, kk, k, trcid

  SHR_ASSERT_ALL((ubound(jwt) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(rootfr) == (/bounds%endp, nlevsoi/)), errMsg(__FILE__, __LINE__))

  associate(                                                   & !
    z              =>    col%z                               , & ! Input:  [real(r8) (:,:)  ]  layer depth (m) (-nlevsno+1:nlevsoi)
    dz             =>    col%dz                              , & ! Input:  [real(r8) (:,:)  ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
    wtcol          =>    pft%wtcol                           , & ! Input:  [real(r8) (:)    ]  weight (relative to column)
    lbl_rsc_h2o    =>    canopystate_vars%lbl_rsc_h2o_patch  , & ! laminar layer resistance for h2o
    elai           =>    canopystate_vars%elai_patch         , &
    annsum_npp     =>    carbonflux_vars%annsum_npp_patch    , & ! Input:  [real(r8) (:) ]  annual sum NPP (gC/m2/yr)
    annavg_agnpp   =>    carbonflux_vars%annavg_agnpp_patch  , & ! Output: [real(r8) (:) ]  annual average above-ground NPP (gC/m2/s)
    annavg_bgnpp   =>    carbonflux_vars%annavg_bgnpp_patch  , & ! Output: [real(r8) (:) ]  annual average below-ground NPP (gC/m2/s)
    frootc         =>    carbonstate_vars%frootc_patch       , & ! Input:  [real(r8) (:)    ]  (gC/m2) fine root C
    is_volatile    =>    betrtracer_vars%is_volatile         , &
    volatilegroupid=>    betrtracer_vars%volatilegroupid     , &
    ngwmobile_tracer_groups=>  betrtracer_vars%ngwmobile_tracer_groups   , &
    tracer_group_memid => betrtracer_vars%tracer_group_memid , &
    t_veg          =>    temperature_vars%t_veg_patch        , &
    t_soisno       =>    temperature_vars%t_soisno_col       , &
    scal_aere_cond =>    tracercoeff_vars%scal_aere_cond_col , &
    tracer_diffusivity_air=> tracercoeff_vars%tracer_diffusivity_air_col, &
    aere_cond      =>    tracercoeff_vars%aere_cond_col        & !
  )


  do j=1,nlevsoi
    do fp = 1, num_soilp
      p = filter_soilp (fp)
      c = pft%column(p)
      g = col%gridcell(c)
      ! Calculate aerenchyma diffusion
      if (j > jwt(c) .and. t_soisno(c,j) > tfrz .and. pft%itype(p) /= noveg) then
        ! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
        ! the new components are in; if there is any tuning to be done to get a realistic global flux,
        ! this would probably be the place.  We will have to document clearly in the Tech Note
        ! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)
        if(use_cn)then
          anpp = annsum_npp(p) ! g C / m^2/yr
          anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

          if (annavg_agnpp(p) /= spval .and. annavg_bgnpp(p) /= spval .and. &
            annavg_agnpp(p) > 0._r8 .and. annavg_bgnpp(p) > 0._r8) then
            nppratio = annavg_bgnpp(p) / (annavg_agnpp(p) + annavg_bgnpp(p))
          else
            nppratio = 0.5_r8
          end if
        endif
        ! Estimate area of tillers (see Wania thesis)
        !m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
        !m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
        ! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
        ! done on a PFT-specific basis.

        if(.not. use_cn)then
          m_tiller = 0._r8        !this was set to zero purposely
        else
          if (usefrootc) then
            m_tiller = frootc(p) ! This will yield much smaller aere area.
          else
            m_tiller = anpp * nppratio * elai(p)
          end if
        endif
        n_tiller = m_tiller / 0.22_r8

        if (pft%itype(p) == nc3_arctic_grass .or. crop(pft%itype(p)) == 1 .or. &
          pft%itype(p) == nc3_nonarctic_grass .or. pft%itype(p) == nc4_grass) then
          poros_tiller = 0.3_r8  ! Colmer 2003
        else
          poros_tiller = 0.3_r8 * nongrassporosratio
        end if


        poros_tiller = poros_tiller * unsat_aere_ratio

        poros_tiller = max(poros_tiller, porosmin)

        area_tiller =  n_tiller * poros_tiller * rpi * 2.9e-3_r8**2._r8 ! (m2/m2)

        do k = 1, ngwmobile_tracer_groups
           trcid = tracer_group_memid(k, 1)
          if(is_volatile(trcid))then
            kk = volatilegroupid(k)
            tracer_diffusivity_air(c,kk) = get_gas_diffusivity(trcid,t_veg(p), betrtracer_vars)
            aerecond = scal_aere_cond(c, kk)*area_tiller * rootfr(p,j) * tracer_diffusivity_air(c,kk) / (z(c,j)*rob)
            ! Add in boundary layer resistance
            lbl_rsc = safe_div(lbl_rsc_h2o(p),  (get_diffusivity_ratio_gas2h2o(trcid, t_veg(p), betrtracer_vars))**(2._r8/3._r8))
            !laminar boundary resistance + resistance in the aerenchyma
            aerecond = safe_div(1._r8, (safe_div(1._r8,aerecond) + safe_div(1._r8,lbl_rsc)))
            aere_cond(c,kk) = aere_cond(c,kk) + wtcol(p) * aerecond
          endif
        enddo
      endif
    end do ! p filter
  end do ! over levels


  end associate

  end subroutine calc_aerecond


  !-----------------------------------------------------------------------
  subroutine betr_annualupdate(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonflux_vars, tracercoeff_vars)
    !
    ! !DESCRIPTION: Annual mean fields.
    !
    ! !USES:
    use clm_time_manager   , only : get_step_size, get_days_per_year, get_nstep
    use clm_varcon         , only : secspday
    use CNCarbonFluxType   , only : carbonflux_type
    use tracercoeffType    , only : tracercoeff_type
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                , intent(in)    :: num_soilp         ! number of soil points in pft filter
    integer                , intent(in)    :: filter_soilp(:)   ! patch filter for soil points
    type(Carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(tracercoeff_type) , intent(inout) :: tracercoeff_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p       ! indices
    integer :: fc        ! soil column filter indices
    integer :: fp        ! soil pft filter indices
    real(r8):: dt        ! time step (seconds)
    real(r8):: secsperyear
    logical :: newrun
    !-----------------------------------------------------------------------

    associate(                                                       &
         agnpp           =>    carbonflux_vars%agnpp_patch         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) aboveground NPP
         bgnpp           =>    carbonflux_vars%bgnpp_patch         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) belowground NPP
         tempavg_agnpp   =>    carbonflux_vars%tempavg_agnpp_patch , & ! Output: [real(r8) (:) ]  temporary average above-ground NPP (gC/m2/s)
         annavg_agnpp    =>    carbonflux_vars%annavg_agnpp_patch  , & ! Output: [real(r8) (:) ]  annual average above-ground NPP (gC/m2/s)
         tempavg_bgnpp   =>    carbonflux_vars%tempavg_bgnpp_patch , & ! Output: [real(r8) (:) ]  temporary average below-ground NPP (gC/m2/s)
         annavg_bgnpp    =>    carbonflux_vars%annavg_bgnpp_patch  , & ! Output: [real(r8) (:) ]  annual average below-ground NPP (gC/m2/s)
         

         annsum_counter  =>    tracercoeff_vars%annsum_counter_col   & ! Output: [real(r8) (:) ]  seconds since last annual accumulator turnover
         !finundated    =>    ch4_vars%finundated_col             , & ! Input:  [real(r8) (:) ]  fractional inundated area in soil column
         !tempavg_somhr =>    ch4_vars%tempavg_somhr_col          , & ! Output: [real(r8) (:) ]  temporary average SOM heterotrophic resp. (gC/m2/s)
         !annavg_somhr    =>    ch4_vars%annavg_somhr_col           , & ! Output: [real(r8) (:) ]  annual average SOM heterotrophic resp. (gC/m2/s)
         !tempavg_finrw  =>    ch4_vars%tempavg_finrw_col          , & ! Output: [real(r8) (:) ]  respiration-weighted annual average of finundated
         !annavg_finrw  =>    ch4_vars%annavg_finrw_col             & ! Output: [real(r8) (:) ]  respiration-weighted annual average of finundated
         )

      ! set time steps
      dt = real(get_step_size(), r8)
      secsperyear = real( get_days_per_year() * secspday, r8)

      newrun = .false.

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (annsum_counter(c) == spval) then
            ! These variables are now in restart files for completeness, but might not be in inicFile and are not.
            ! set for arbinit.
            newrun = .true.
            annsum_counter(c)    = 0._r8
            !tempavg_somhr(c)     = 0._r8
            !tempavg_finrw(c)     = 0._r8
         end if

         annsum_counter(c) = annsum_counter(c) + dt
      end do

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         if (newrun .or. tempavg_agnpp(p) == spval) then ! Extra check needed because for back-compatibility
            tempavg_agnpp(p) = 0._r8
            tempavg_bgnpp(p) = 0._r8
         end if
      end do


      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = pft%column(p)
         if (annsum_counter(c) >= secsperyear) then

            annavg_agnpp(p) = tempavg_agnpp(p)
            tempavg_agnpp(p) = 0._r8

            annavg_bgnpp(p) = tempavg_bgnpp(p)
            tempavg_bgnpp(p) = 0._r8

         else
            tempavg_agnpp(p) = tempavg_agnpp(p) + dt/secsperyear * agnpp(p)
            tempavg_bgnpp(p) = tempavg_bgnpp(p) + dt/secsperyear * bgnpp(p)
         end if
      end do

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         if (annsum_counter(c) >= secsperyear) annsum_counter(c) = 0._r8
      end do

    end associate

  end subroutine betr_annualupdate
end module TracerParamsMod
