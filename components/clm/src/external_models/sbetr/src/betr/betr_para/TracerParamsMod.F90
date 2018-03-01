module TracerParamsMod
#include "bshr_assert.h"

  ! !DESCRIPTION:
  ! Module holding routines used to compute solubility, and phase conversion parameters
  ! to be used for BeTR 1D vertical tracer transport
  !
  !
  ! History
  ! Jinyun Tang created May 2014.
  ! !USES:
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use bshr_log_mod             , only : errMsg => shr_log_errMsg
  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type
  use tracer_varcon            , only : nlevsoi  => betr_nlevsoi
  use betr_varcon              , only : spval => bspval
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  use TracerParamSetMod        , only : get_lgsorb_KL_Xsat, get_lnsorb_Kd
  use TracerParamSetMod        , only : get_henrycef, get_equilibrium_scal
  use TracerParamSetMod        , only : get_tauliq, get_taugas
  use TracerParamSetMod        , only : get_aqueous_diffusivity, get_gas_diffusivity
  use TracerParamSetMod        , only : get_diffusivity_ratio_gas2h2o, rhosat
  use TracerParamSetWatIsoMod  , only : get_equi_lv_h2oiso_fractionation
  use TracerParamSetWatIsoMod  , only : get_equi_sv_h2oiso_fractionation
  use TracerParamSetWatIsoMod  , only : get_equi_sl_h2oiso_fractionation
  implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: tracer_param_init
  public :: set_multi_phase_diffusion
  public :: convert_mobile2gas
  public :: set_phase_convert_coeff
  public :: calc_tracer_infiltration
  public :: get_zwt
  public :: calc_aerecond
  public :: betr_annualupdate
  !parameters
  character(len=*), parameter :: filename = &
       __FILE__
  real(r8), parameter         :: minval_diffus = 1.e-20_r8   !minimum diffusivity, m2/s
  real(r8), parameter         :: minval_airvol = 1.e-10_r8   !minimum air-filled volume

  !declare a private tortuosity type
  type :: soil_tortuosity_type
     real(r8), pointer :: tau_gas(:,:)      !soil tortuosity for gaseous phase diffusion
     real(r8), pointer :: tau_liq(:,:)      !soil tortuosity for aqueous phase diffusion
  end type soil_tortuosity_type
  type(soil_tortuosity_type), target :: tau_soil

  !!
contains


  subroutine tracer_param_init(bounds)

    !
    ! !DESCRIPTION:
    !
    ! initialize the tracerParamsMod
    !
    ! !USES:
    use tracer_varcon         , only : nlevtrc_soil => betr_nlevtrc_soil

    implicit none
    type(bounds_type), intent(in) :: bounds   !bounds
    character(len=32)             :: subname ='tracer_param_init'

    allocate(tau_soil%tau_gas(bounds%begc:bounds%endc, 1 : nlevtrc_soil))
    tau_soil%tau_gas(:,:) = 0._r8
    allocate(tau_soil%tau_liq(bounds%begc:bounds%endc, 1 : nlevtrc_soil))
    tau_soil%tau_liq(:,:) = 0._r8


  end subroutine tracer_param_init

  !--------------------------------------------------------------------------------------------------------------
  subroutine Calc_gaseous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       biophysforc, tau_gas, bstatus)
    !
    ! !DESCRIPTION:
    !
    ! compute soil tortuosity for gasesous diffusion

    ! !USES:
    use BetrStatusType         , only : betr_status_type
    implicit none
    !arguments
    type(bounds_type)                , intent(in)    :: bounds                             ! bounds
    integer                          , intent(in)    :: num_soilc                          ! number of column soil points in column filter
    integer                          , intent(in)    :: filter_soilc(:)                    ! column filter for soil points
    integer                          , intent(in)    :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0
    integer                          , intent(in)    :: jtops(bounds%begc: )               ! top label of each column
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    real(r8)                         , intent(inout) :: tau_gas(bounds%begc: , lbj: )      !output variable
    type(betr_status_type)           , intent(out)   :: bstatus
    !local variables
    integer :: n, fc, c     !indices
    character(len=255) :: subname = 'calc_gaseous_diffusion_soil_tortuosity'

    call bstatus%reset()
    associate(                                           &
         eff_porosity   => biophysforc%eff_porosity_col, & !effective soil porosity
         bsw            => biophysforc%bsw_col         , & !clapp-hornber shape parameters
         air_vol        => biophysforc%air_vol_col       & !volume possessed by air
         )

      SHR_ASSERT_ALL((ubound(jtops)   == (/bounds%endc/)),  errMsg(filename,__LINE__), bstatus)
      if(bstatus%check_status())return
      SHR_ASSERT_ALL((ubound(tau_gas) == (/bounds%endc, ubj/)),   errMsg(filename,__LINE__), bstatus)
      if(bstatus%check_status())return

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
  subroutine Calc_aqueous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, numf, filter, &
    biophysforc, tau_liq, betr_status)
    !
    ! DESCRIPTIONS
    ! compute soil tortuosity for aquesous diffusion
    !
    use BetrStatusType         , only : betr_status_type
    implicit none
    !arguments
    type(bounds_type)                , intent(in)    :: bounds                                ! bounds
    integer                          , intent(in)    :: numf                                  ! number of columns in column filter
    integer                          , intent(in)    :: filter(:)                             ! column filter
    integer                          , intent(in)    :: lbj, ubj                              ! lower and upper bounds, make sure they are > 0
    integer                          , intent(in)    :: jtops(bounds%begc: )                  ! top label of each column
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    real(r8)                         , intent(inout) :: tau_liq(bounds%begc: , lbj: )         !output variable
    type(betr_status_type)           , intent(out)   :: betr_status

    !local variables
    integer :: n, fc, c     !indices
    character(len=255) :: subname = 'calc_aqueous_diffusion_soil_tortuosity'

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops)    == (/bounds%endc/)),        errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(tau_liq)  == (/bounds%endc, ubj/)),   errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return

    associate(                                           &
         eff_porosity   => biophysforc%eff_porosity_col, & !effective soil porosity
         bsw            => biophysforc%bsw_col         , & !clapp-hornber shape parameters
         h2osoi_liqvol  => biophysforc%h2osoi_liqvol_col & !soil volume possessed by liquid water
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

  subroutine calc_bulk_diffusivity(bounds, col, lbj, ubj, jtops, numf, filter, bunsencef_col, &
       biophysforc, tau_soi, betrtracer_vars,  bulkdiffus, betr_status)
    !
    ! !DESCRIPTION:
    ! compute the weighted bulk diffusivity in soil for dual-phase transport
    ! Reference: Tang and Riley, 2014, BG,  Simple formulations and solutions of &
    ! the dual-phase diffusive transport for biogeochemical modeling.
    !the formula for a volatile species is
    !D_bulk=(airvol*D_g*tau_g+bunsencef_col*h2osoi_liqvol*D_w*tau_w)

    ! !USES:
    use BeTRTracerType        , only : betrtracer_type
    use BetrStatusType        , only : betr_status_type
    use betr_columnType       , only : betr_column_type
    implicit none
    type(bounds_type)                , intent(in)  :: bounds                                  ! bounds
    type(betr_column_type)           , intent(in)  :: col
    integer                          , intent(in)  :: numf                                    ! number of columns in column filter
    integer                          , intent(in)  :: filter(:)                               ! column filter
    integer                          , intent(in)  :: lbj, ubj                                ! lower and upper bounds, make sure they are > 0
    integer                          , intent(in)  :: jtops(bounds%begc: )                    ! top label of each column
    real(r8)                         , intent(in)  :: bunsencef_col(bounds%begc: ,lbj: ,1: )  ! bunsen coefficient for gaseous-aqueous conversion
    type(betrtracer_type)            , intent(in)  :: betrtracer_vars                         ! betr configuration information
    type(betr_biogeophys_input_type) , intent(in)  :: biophysforc
    type(soil_tortuosity_type)       , intent(in)  :: tau_soi                                 ! soil tortuosity
    real(r8)                         , intent(out) :: bulkdiffus(bounds%begc: ,lbj: , 1: )    ! the returning variable
    type(betr_status_type)           , intent(out) :: betr_status
    !local variables
    real(r8)           :: max_depth_cryoturb         = 3._r8  !m
    !parameters below will be encapsulated into a structure later
    real(r8)           :: max_altdepth_cryoturbation = 1._r8  ! (m) maximum active layer thickness for cryoturbation to occur
    real(r8)           :: cryoturb_diffusion_k       = 1e-4_r8 / (86400._r8 * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr = 1m^2/1000 yr
    real(r8)           :: som_diffus                 = 5e-4_r8 / (86400._r8 * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr
    integer            :: j, k, n, fc, c , trcid    !indices
    integer            :: nsld
    real(r8)           :: diffaqu, diffgas
    character(len=255) :: subname = 'calc_bulk_diffusivity'

    integer :: nvolatile_tracer_groups

    call betr_status%reset()
    nvolatile_tracer_groups = betrtracer_vars%nvolatile_tracer_groups

    !array shape checking will be added later.

    SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),   errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(jtops)           == (/bounds%endc/)),   errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(bunsencef_col,1) == bounds%endc), errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(bunsencef_col,2) == ubj), errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(bunsencef_col,3) == nvolatile_tracer_groups), errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(bulkdiffus,1)    == bounds%endc), errMsg(filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(bulkdiffus,2)    == ubj), errMsg(filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(bulkdiffus,3)    == betrtracer_vars%ntracer_groups), errMsg(filename,__LINE__),betr_status)
    if(betr_status%check_status())return

    associate(                                                                                              &
         ngwmobile_tracer_groups                => betrtracer_vars%ngwmobile_tracer_groups                , & ! Integer[intent(in)], number of dual phase (gw) tracers
         tracer_group_memid                     => betrtracer_vars%tracer_group_memid                     , & !
         ntracer_groups                         => betrtracer_vars%ntracer_groups                         , & ! Integer[intent(in)], total number of tracers
         is_volatile                            => betrtracer_vars%is_volatile                            , & ! logical[intent(in)], is a volatile tracer?
         is_h2o                                 => betrtracer_vars%is_h2o                                 , & ! logical[intent(in)], is a h2o tracer?
         volatilegroupid                        => betrtracer_vars%volatilegroupid                        , & ! integer[intent(in)], location in the volatile vector
         air_vol                                => biophysforc%air_vol_col                            ,     & ! volume possessed by air
         h2osoi_liqvol                          => biophysforc%h2osoi_liqvol_col                      ,     & ! soil volume possessed by liquid water
         altmax                                 => biophysforc%altmax_col                            ,      & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear                        => biophysforc%altmax_lastyear_col                   ,      & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth o
         tracer_solid_passive_diffus_scal_group => betrtracer_vars%tracer_solid_passive_diffus_scal_group , & !scaling factor for solid phase diffusivity
         tracer_solid_passive_diffus_thc_group  => betrtracer_vars%tracer_solid_passive_diffus_thc_group  , & !threshold for solid phase diffusivity
         tau_gas                                => tau_soi%tau_gas                                        , & ! real(r8)[intent(in)], gaseous tortuosity
         tau_liq                                => tau_soi%tau_liq                                        , & ! real(r8)[intent(in)], aqueous tortuosity
         zi                                     => col%zi                                                 , & ! real(r8)[intent(in)],
         t_soisno                               => biophysforc%t_soisno_col                                 & ! Input: [real(r8)(:,:)]
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
      !the following setup is adapted from CLM4.5
      do j = ngwmobile_tracer_groups + 1, ntracer_groups
         nsld = j - ngwmobile_tracer_groups
         trcid = tracer_group_memid(j,1)
         do n = 1, ubj
            do fc = 1,numf
               c = filter(fc)
               if  ( ( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
                  ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
                  ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
                  if ( zi(c,n) < max(altmax(c), altmax_lastyear(c)) ) then
                    bulkdiffus(c,n,j) = cryoturb_diffusion_k * tracer_solid_passive_diffus_scal_group(nsld)
                    bulkdiffus(c,n,j) = max(bulkdiffus(c,n,j), tracer_solid_passive_diffus_thc_group(nsld))
                  else
                    bulkdiffus(c,n,j) = max(cryoturb_diffusion_k * &
                          ( 1._r8 - ( zi(c,n) - max(altmax(c), altmax_lastyear(c)) ) / &
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


   subroutine calc_bulk_conductances(bounds, lbj, ubj, jtops, numf, filter, &
     bulkdiffus, dz, betrtracer_vars,  hmconductance_col, betr_status)
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
   use BetrStatusType     , only : betr_status_type
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
   type(betr_status_type), intent(out):: betr_status
   !local variables

   integer :: j, n, fc, c                 !indices
   character(len=255) :: subname = 'calc_bulk_conductances'
   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)             == (/bounds%endc/)),        errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(dz,1)       == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(dz,2)       == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bulkdiffus,1) == bounds%endc), errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bulkdiffus,2) == ubj), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bulkdiffus,3) == betrtracer_vars%ntracer_groups), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(hmconductance_col,1) == bounds%endc), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(hmconductance_col,2) == ubj-1), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(hmconductance_col,3) == betrtracer_vars%ntracer_groups), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return

   associate(                                                                 &
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
             hmconductance_col(bounds%begc:bounds%endc, lbj:ubj-1, j), betr_status)
     if(betr_status%check_status())return
   enddo

   end associate
   end subroutine calc_bulk_conductances

!-------------------------------------------------------------------------------
   subroutine calc_henrys_coeff(bounds, lbj, ubj, jtops, numf, filter, t_soisno, soi_pH, &
       betrtracer_vars, aqu2neutralcef_col, henrycef_col, betr_status)
   !
   ! DESCRIPTION
   ! compute henry's law constant for volatile tracers
   use BeTRTracerType     , only : betrtracer_type
   use BetrStatusType     , only : betr_status_type
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

   real(r8)            ,   intent(inout):: aqu2neutralcef_col(bounds%begc: , lbj: , 1: ) !conversion parameter between bulk aqueous and neutral aqueous tracer
   real(r8)            ,   intent(inout):: henrycef_col(bounds%begc: , lbj: ,  1: )       !henry's constant, mol/L/atm = M/atm
   type(betr_status_type), intent(out)  :: betr_status
   !local variables
   integer :: j, k, n, fc, c, trcid   ! indices
   real(r8) :: scal
   character(len=255) :: subname='calc_henrys_coeff'
   integer :: nvolatile_tracer_groups, ngwmobile_tracer_groups

   call betr_status%reset()
   ngwmobile_tracer_groups = betrtracer_vars%ngwmobile_tracer_groups
   nvolatile_tracer_groups = betrtracer_vars%nvolatile_tracer_groups

   SHR_ASSERT_ALL((ubound(jtops)             == (/bounds%endc/)),        errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(t_soisno,1)       ==  bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(t_soisno,2)       == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(soi_pH,1)         == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(soi_pH,2)         == ubj),  errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(aqu2neutralcef_col,1)== bounds%endc), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(aqu2neutralcef_col,2)== ubj), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(aqu2neutralcef_col,3)== ngwmobile_tracer_groups), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(henrycef_col,1)   == bounds%endc), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(henrycef_col,2)   == ubj), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(henrycef_col,3)   == nvolatile_tracer_groups), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return

   associate(                                                               &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups , & !Integer[intent(in)], number of tracers
    is_volatile                => betrtracer_vars%is_volatile             , & !logical[intent(in)], is a volatile tracer?
    is_h2o                     => betrtracer_vars%is_h2o                  , & !logical[intent(in)], is a h2o tracer?
    tracer_group_memid         => betrtracer_vars%tracer_group_memid      , & !integer[intent(in)], tracer id
    volatilegroupid            => betrtracer_vars%volatilegroupid         , & !integer[intent(in)], location in the volatile vector
    tracerfamilyname           => betrtracer_vars%tracerfamilyname          &
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
             scal = get_equilibrium_scal(t_soisno(c,n), soi_pH(c,n), tracerfamilyname(trcid),betrtracer_vars)
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
   subroutine calc_bunsen_coeff(bounds, lbj, ubj, jtops, numf, filter, &
        henrycef_col, t_soisno, smp_l, betrtracer_vars, bunsencef_col, betr_status)
   !
   ! DESCRIPTION
   ! compute Bunsen's coefficient
   !
   use betr_varcon        , only : denh2o  => bdenh2o
   use BeTRTracerType     , only : betrtracer_type
   use BetrStatusType     , only : betr_status_type
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
   type(betr_status_type) , intent(out)   :: betr_status
   !local variables
   integer            :: j, k, n, fc, c , trcid       !indices
   real(r8)           :: rho_vap(bounds%begc:bounds%endc, lbj:ubj)                           ! saturated vapor pressure for different layers
   character(len=255) :: subname = 'calc_bunsen_coeff'
   integer            :: nvolatile_tracer_groups

   call betr_status%reset()
   nvolatile_tracer_groups = betrtracer_vars%nvolatile_tracer_groups

   SHR_ASSERT_ALL((ubound(jtops)             == (/bounds%endc/)),        errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(t_soisno,1)     == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(t_soisno,2)     == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(smp_l,1)        == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(smp_l,2)        == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(henrycef_col,1)  == bounds%endc), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(henrycef_col,2)  == ubj), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(henrycef_col,3)  == nvolatile_tracer_groups), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bunsencef_col,1)  == bounds%endc), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bunsencef_col,2)  == ubj), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bunsencef_col,3)  == nvolatile_tracer_groups), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return

   associate(                                                                    &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups      , & !Integer[intent(in)], number of tracers
    tracer_group_memid         => betrtracer_vars%tracer_group_memid           , &
    is_volatile                => betrtracer_vars%is_volatile                  , & !logical[intent(in)], is a volatile tracer?
    is_h2o                     => betrtracer_vars%is_h2o                       , & !logical[intent(in)], is a h2o tracer
    volatilegroupid            => betrtracer_vars%volatilegroupid                & !integer[intent(in)], location in the volatile vector
   )
   if(any(is_h2o))then
     call calc_rhovap(bounds, lbj, ubj, jtops, numf, filter, t_soisno, smp_l, rho_vap, betr_status)
     if(betr_status%check_status())return
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
     biophysforc, betrtracer_vars, tracerstate_vars, tracercoeff_vars, betr_status)

   !DESCRIPTIONS:
   !compute phase conversion coefficients between gaseous and aqueous phases
   ! The total aqueous and gases phase concentration is = theta*aqueous+epsilon*gaseous
   ! because aqueous = bunsen*gaseous, these coefficients are constant throughout the all period.

   !USES:
   use BeTRTracerType     , only : betrtracer_type
   use TracerCoeffType    , only : tracercoeff_type
   use betr_varcon        , only : denh2o => bdenh2o, denice => bdenice
   use BetrStatusType     , only : betr_status_type
   use tracer_varcon      , only : sorp_isotherm_linear, sorp_isotherm_langmuir
   use tracerstatetype    , only : tracerstate_type
   implicit none
   !arguments
   type(bounds_type)                , intent(in)    :: bounds                      ! bounds
   integer                          , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
   integer                          , intent(in)    :: jtops(bounds%begc: )        ! top label of each column
   integer                          , intent(in)    :: numf                        ! number of columns in column filter
   integer                          , intent(in)    :: filter(:)                   ! column filter
   type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   type(tracerstate_type)           , intent(in)    :: tracerstate_vars            !
   type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars            ! structure containing tracer transport parameters
   type(betr_status_type)           , intent(out)   :: betr_status
   !local variables
   integer            :: j, n, k, fc, c , trcid  ! indices
   real(r8) :: Kd, KL, xs, Xsat, scal
   real(r8), parameter :: tiny_val = 1.e-12_r8
   character(len=255) :: subname = 'calc_dual_phase_convert_coeff'

   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)   == (/bounds%endc/)),        errMsg(filename,__LINE__),betr_status)
   associate(                                                                    &
    ngwmobile_tracer_groups    => betrtracer_vars%ngwmobile_tracer_groups      , & !Input: [integer(:)], number of tracers
    tracer_group_memid         => betrtracer_vars%tracer_group_memid           , & !Input: [integer(:)], tracer id
    is_h2o                     => betrtracer_vars%is_h2o                       , & !Input: [logical(:)], is h2o tracer?
    is_volatile                => betrtracer_vars%is_volatile                  , & !Input: [logical(:)], is a volatile tracer?
    volatilegroupid            => betrtracer_vars%volatilegroupid              , & !Input: [logical(:)], location in the volatile vector
    adsorbgroupid              => betrtracer_vars%adsorbgroupid                , & !Input: [Integer(:)], tracer id
    is_adsorb                  => betrtracer_vars%is_adsorb                    , & !Input: [logical(:)]
    adsorb_isotherm            => betrtracer_vars%adsorb_isotherm              , & !Input: [integer]
    tracerfamilyname           => betrtracer_vars%tracerfamilyname             , & !Input: [char(:)], tracer family name
    h2osoi_liqvol              => biophysforc%h2osoi_liqvol_col                , & !Input: [real(r8)(:,:)], liquid h2o vol
    h2osoi_icevol              => biophysforc%h2osoi_icevol_col                , & !Input: [real(r8)(:,:)], ice h2o vol
    air_vol                    => biophysforc%air_vol_col                      , & !Input: [real(r8)(:,:)], air vol
    watsat                     => biophysforc%watsat_col                       , & !Input: [real(r8)(:,:)]
    isoilorder                 => biophysforc%isoilorder                       , & !Input: [integer(:)]
    soil_pH                    => biophysforc%soil_pH                          , & !Input:
    t_soisno                   => biophysforc%t_soisno_col                     , & !Input:
    tracer_conc_mobile         =>  tracerstate_vars%tracer_conc_mobile_col     , & !Input
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
            if(is_h2o(j))then
              !this is a (bad) reverse hack because the hydrology code does not consider water vapor transport
              !jyt, Feb, 18, 2016, 1.e-12_r8 is a value for avoiding NaN

              aqu2bulkcef_mobile(c,n,j) = max(h2osoi_liqvol(c,n),tiny_val)
            else
              aqu2bulkcef_mobile(c,n,j) = air_vol(c,n)/bunsencef_col(c,n,k)+h2osoi_liqvol(c,n)
            endif
            !gaseous to bulk mobile phase
            gas2bulkcef_mobile(c,n,k) = air_vol(c,n)+h2osoi_liqvol(c,n)*bunsencef_col(c,n,k)
            !correct for impermeable layer, to avoid division by zero in doing diffusive transport
            gas2bulkcef_mobile(c,n,k) = max(gas2bulkcef_mobile(c,n,k),air_vol(c,n),minval_airvol)
          endif
        enddo
      enddo
      if(is_adsorb(trcid))then
        if(adsorb_isotherm(trcid)==sorp_isotherm_linear)then
          do n = lbj, ubj
            do fc = 1, numf
              c = filter(fc)
              if(n>=jtops(c))then
                scal = get_equilibrium_scal(t_soisno(c,n), soil_pH(c,n), tracerfamilyname(trcid),betrtracer_vars)
                Kd=get_lnsorb_Kd(tracerfamilyname(j))
                if(scal/=1._r8)then
                 !Because bunsen = bunsen0*scal
                 !
                  Kd = (scal-1._r8)/(1._r8+scal)*Kd
                endif
                aqu2bulkcef_mobile(c,n,j) = air_vol(c,n)/bunsencef_col(c,n,k) + (1._r8 +  Kd)*h2osoi_liqvol(c,n)
              endif
            enddo
          enddo
        elseif(adsorb_isotherm(trcid)==sorp_isotherm_langmuir)then
          call get_lgsorb_KL_Xsat(tracerfamilyname(j), isoilorder(c), KL, Xsat)
          do n = lbj, ubj
            do fc = 1, numf
              c = filter(fc)
              if(n>=jtops(c))then
                xs = Xsat/watsat(c,j) * h2osoi_liqvol(c,n)
                if(scal/=1._r8)then
                  KL=KL*scal/(scal-1._r8) * (air_vol(c,n)/bunsencef_col(c,n,k)+h2osoi_liqvol(c,n))
                else
                  KL=KL * (air_vol(c,n)/bunsencef_col(c,n,k)+h2osoi_liqvol(c,n))
                endif
                Kd = xs/(KL +tracer_conc_mobile(c,n,trcid))
                aqu2bulkcef_mobile(c,n,j) = (air_vol(c,n)/bunsencef_col(c,n,k)+h2osoi_liqvol(c,n)) * (1._r8+Kd)
                gas2bulkcef_mobile(c,n,k) = aqu2bulkcef_mobile(c,n,j)*bunsencef_col(c,n,k)
              endif
            enddo
          enddo
        endif
      endif
    else
      !when linear adsorption is used for some adsorptive aqueous tracers, the aqu2bulkcef will be the retardation factor
      !for the moment, it is set to one for all non-volatile tracers
      !It is assumed that ice have same equilibrium solubility as liquid water for soluble tracers
      do n = lbj, ubj
        do fc = 1, numf
          c = filter(fc)
          if(n>=jtops(c))then
            aqu2bulkcef_mobile(c, n, j) = max(h2osoi_liqvol(c,n),tiny_val)
          endif
        enddo
      enddo
      if(is_adsorb(trcid))then
        if(adsorb_isotherm(trcid)==sorp_isotherm_linear)then
          do n = lbj, ubj
            do fc = 1, numf
              c = filter(fc)
              if(n>=jtops(c))then
                Kd=get_lnsorb_Kd(tracerfamilyname(j))
                aqu2bulkcef_mobile(c, n, j) = aqu2bulkcef_mobile(c, n, j) * (1._r8+Kd)
              endif
            enddo
          enddo
        elseif(adsorb_isotherm(trcid)==sorp_isotherm_langmuir)then
          !the adsorption parameter should be a function of soil type, or soil order
          call get_lgsorb_KL_Xsat(tracerfamilyname(j), isoilorder(c), KL, Xsat)
          do n = lbj, ubj
            do fc = 1, numf
              c = filter(fc)
              if(n>=jtops(c))then
                KL= h2osoi_liqvol(c,n) * KL
                xs = Xsat/watsat(c,j) * h2osoi_liqvol(c,n)
                Kd = xs/(KL + tracer_conc_mobile(c,n,trcid))
                aqu2bulkcef_mobile(c, n, j) = aqu2bulkcef_mobile(c, n, j) * (1._r8+Kd)
              endif
            enddo
          enddo
        endif
      endif
    endif
  enddo
  end associate
  end subroutine calc_dual_phase_convert_coeff



!-------------------------------------------------------------------------------
   subroutine convert_mobile2gas(bounds, lbj, ubj, jtops, numf, filter, &
        do_forward, gas2bulkcef_mobile_col, betrtracer_vars, tracer_conc_mobile, &
        betr_status)
   !
   ! DESCRIPTIONS
   ! do conversion between bulk mobile phase and gaseous phase
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   use BetrStatusType     , only : betr_status_type
   implicit none
   !arguments
   type(bounds_type)     , intent(in)    :: bounds                  ! bounds
   integer               , intent(in)    :: lbj, ubj                ! lower and upper bounds
   integer               , intent(in)    :: jtops(bounds%begc: )    ! top label of each column
   integer               , intent(in)    :: numf                    ! number of filters
   integer               , intent(in)    :: filter(:)               ! filter
   logical               , intent(in)    :: do_forward              ! true, dual_bulk => gaseous
   type(betrtracer_type) , intent(in)    :: betrtracer_vars         ! betr configuration information
   real(r8)              , intent(in)    :: gas2bulkcef_mobile_col(bounds%begc: ,lbj: , 1: )  !conversion parameter
   real(r8)              , intent(inout) :: tracer_conc_mobile(bounds%begc: ,lbj: , 1: )  !bulk mobile tracer
   type(betr_status_type), intent(out)   :: betr_status
   !local variables
   integer :: jj, kk, fc, c, j
   character(len=255) :: subname = 'convert_mobile2gas'
   integer :: nvolatile_tracers

   call betr_status%reset()
   nvolatile_tracers = betrtracer_vars%nvolatile_tracers

   SHR_ASSERT_ALL((ubound(jtops)                  == (/bounds%endc/)),   errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(gas2bulkcef_mobile_col,1) == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(gas2bulkcef_mobile_col,2) == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(gas2bulkcef_mobile_col,3) == nvolatile_tracers),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(tracer_conc_mobile,1)     == bounds%endc),   errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(tracer_conc_mobile,2)     == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(tracer_conc_mobile,3)     == nvolatile_tracers),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return

   associate(                                                  &
    ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers   , & !Integer[intent(in)], number of tracers
    is_volatile       => betrtracer_vars%is_volatile         , & !logical[intent(in)], is a volatile tracer?
    volatilegroupid   => betrtracer_vars%volatilegroupid       & !integer[intent(in)], location in the volatile vector
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

   subroutine set_multi_phase_diffusion(bounds, col, lbj, ubj, jtops, numf, filter, &
      biophysforc, betrtracer_vars, tracercoeff_vars, betr_status)
   !
   ! DESCRIPTION
   ! set parameters for the dual phase diffusion
   !
   !USES
   use TracerCoeffType    , only : tracercoeff_type
   use BeTRTracerType     , only : betrtracer_type
   use BetrStatusType     , only : betr_status_type
   use betr_columnType    , only : betr_column_type
   implicit none
   !ARGUMENTS
   type(bounds_type)                , intent(in)    :: bounds  ! bounds
   type(betr_column_type)           , intent(in)    :: col
   integer                          , intent(in)    :: lbj, ubj             ! lower and upper bounds, make sure they are > 0
   integer                          , intent(in)    :: jtops(bounds%begc: ) ! top label of each column
   integer                          , intent(in)    :: numf                 ! number of columns in column filter
   integer                          , intent(in)    :: filter(:)            ! column filter
   type(betrtracer_type)            , intent(in)    :: betrtracer_vars      ! betr configuration information
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars ! structure containing tracer transport parameters
   type(betr_status_type)           , intent(out)   :: betr_status
   !
   real(r8) :: bulkdiffus(bounds%begc:bounds%endc,lbj:ubj,1:betrtracer_vars%ntracer_groups )  !weighted bulk diffusivity for dual-phase diffusion

   !maybe I should use tau_soil as a local variable, I will check this later
   character(len=255) :: subname='set_multi_phase_diffusion'

   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)   == (/bounds%endc/)), errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return

   !compute tortuosity
   !gaseous phase
   call calc_gaseous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, numf, filter, &
        biophysforc, tau_soil%tau_gas, betr_status)
   if(betr_status%check_status())return

   !aqueous phase
   call calc_aqueous_diffusion_soil_tortuosity(bounds, lbj, ubj, jtops, numf, filter, &
        biophysforc, tau_soil%tau_liq, betr_status)
   if(betr_status%check_status())return

   !compute bulk diffusivity
   call calc_bulk_diffusivity(bounds, col, lbj, ubj, jtops, numf, filter       , &
      tracercoeff_vars%bunsencef_col(bounds%begc:bounds%endc,lbj:ubj, : ) , &
      biophysforc, tau_soil, betrtracer_vars, bulkdiffus, betr_status)
   if(betr_status%check_status())return

   !compute weigthed conductances
   call calc_bulk_conductances(bounds, lbj, ubj, jtops, numf, filter, bulkdiffus, &
      col%dz(bounds%begc:bounds%endc,lbj:ubj), betrtracer_vars,  &
      tracercoeff_vars%hmconductance_col(bounds%begc:bounds%endc, lbj:ubj-1, :), betr_status)

   end subroutine set_multi_phase_diffusion


!--------------------------------------------------------------------------------
   subroutine set_phase_convert_coeff(bounds, lbj, ubj, jtops, numf, filter, &
        dz, biophysforc, betrtracer_vars, tracerstate_vars, tracercoeff_vars, betr_status)
   !
   ! DESCRIPTION
   ! set parameters for phase conversion
   use TracerCoeffType    , only : tracercoeff_type
   use BeTRTracerType     , only : betrtracer_type
   use BetrStatusType     , only : betr_status_type
   use tracerstatetype    , only : tracerstate_type
   implicit none
   type(bounds_type)                , intent(in)    :: bounds  ! bounds
   integer                          , intent(in)    :: lbj, ubj             ! lower and upper bounds, make sure they are > 0
   integer                          , intent(in)    :: jtops(bounds%begc:bounds%endc)        ! top label of each column
   integer                          , intent(in)    :: numf                 ! number of columns in column filter
   integer                          , intent(in)    :: filter(:)            ! column filter
   real(r8)                         , intent(in)    :: dz(bounds%begc: ,lbj: )
   type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   type(TracerState_type)           , intent(in)    :: tracerstate_vars
   type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars ! structure containing tracer transport parameters
   type(betr_status_type)           , intent(out)   :: betr_status

   character(len=255) :: subname = 'set_phase_convert_coeff'

   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)   == (/bounds%endc/)),        errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(dz)      == (/bounds%endc, ubj/)),   errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
    !compute Henry's law constant
   call calc_henrys_coeff(bounds, lbj, ubj, jtops, numf, filter                    , &
       biophysforc%t_soisno_col(bounds%begc:bounds%endc,lbj:ubj)              ,      &
       biophysforc%soil_pH(bounds%begc:bounds%endc, lbj:ubj),  betrtracer_vars  ,    &
       tracercoeff_vars%aqu2neutralcef_col(bounds%begc:bounds%endc,lbj:ubj, : )    , &
       tracercoeff_vars%henrycef_col(bounds%begc:bounds%endc, lbj:ubj, : ), betr_status)
    if(betr_status%check_status())return

   !compute Bunsen's coefficients
   call calc_bunsen_coeff(bounds, lbj, ubj, jtops, numf, filter                    , &
        tracercoeff_vars%henrycef_col(bounds%begc:bounds%endc, lbj:ubj, : )        , &
        biophysforc%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj)            ,      &
        biophysforc%smp_l_col    (bounds%begc:bounds%endc, lbj:ubj)            ,     &
        betrtracer_vars                                                            , &
        tracercoeff_vars%bunsencef_col(bounds%begc:bounds%endc, lbj:ubj, :), betr_status)
    if(betr_status%check_status())return

   !compute equilibrium fraction to liquid phase conversion parameter
   if(betrtracer_vars%nsolid_equil_tracers>0)then
     call calc_equil_to_liquid_convert_coeff(bounds, lbj, ubj, jtops, numf, filter , &
        biophysforc%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj)            ,      &
        biophysforc%h2osoi_ice_col(bounds%begc:bounds%endc,lbj:ubj)            ,     &
        dz(bounds%begc:bounds%endc, lbj:ubj)                                   ,     &
        betrtracer_vars, &
        tracercoeff_vars%aqu2equilsolidcef_col(bounds%begc:bounds%endc, lbj:ubj,:), betr_status)
     if(betr_status%check_status())return
   endif

   !compute phase conversion coefficients
   call calc_dual_phase_convert_coeff(bounds, lbj, ubj, jtops, numf, filter, &
      biophysforc, betrtracer_vars, tracerstate_vars, tracercoeff_vars, betr_status)

   end subroutine set_phase_convert_coeff



  !------------------------------------------------------------------------
   subroutine calc_tracer_infiltration(bounds, jtops, numf, filter, bunsencef_topsoi, &
     betrtracer_vars, tracerboundarycond_vars, biogeo_flux,  tracer_flx_infl, betr_status)
   !
   ! DESCRIPTION
   ! calculate advection velocity for BeTR code
   ! this assumes the interfacial velocity qflx_adv (except infiltration) has been calcualted in doing vertical
   ! watermovement
   ! This assumes the advection solves the equation
   ! \frac{\parital m}{\partial t}+\frac{V*m}{\partial z} = 0
   !where m = C*vsm, therefore, V=ql/vsm
   !
   !USES
   use TracerBoundaryCondType, only : tracerboundarycond_type
   use BeTRTracerType        , only : betrtracer_type
   use betr_varcon           , only : denh2o  => bdenh2o
   use BetrStatusType        , only : betr_status_type
   implicit none
   !ARGUMENTS
   type(bounds_type)            , intent(in)    :: bounds
   integer                      , intent(in)    :: jtops(bounds%begc: )           ! top label of each column
   integer                      , intent(in)    :: numf                           ! number of columns in column filter
   integer                      , intent(in)    :: filter(:)                      ! column filter
   real(r8)                     , intent(in)    :: bunsencef_topsoi(bounds%begc: , 1: )
   type(BeTRTracer_Type)        , intent(in)    :: betrtracer_vars
   type(tracerboundarycond_type), intent(in)    :: tracerboundarycond_vars
   type(betr_biogeo_flux_type)  , intent(in)    :: biogeo_flux
   real(r8)                     , intent(inout) :: tracer_flx_infl(bounds%begc: , 1: )
   type(betr_status_type)       , intent(out)   :: betr_status

   ! local variables
   integer :: fc, c, j

   call betr_status%reset()
    ! remove compiler warnings for unused dummy args
    !if (len(betrtracer_vars%betr_simname) > 0) continue

   associate(                                          &
    nvolatile_tracer_groups => betrtracer_vars%nvolatile_tracer_groups, &
    qflx_adv            =>    biogeo_flux%qflx_adv_col & !real(r8) (:)  [intent(in)], infiltration, mm/s
   )

   SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bunsencef_topsoi,1) == bounds%endc ), errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(bunsencef_topsoi,2) == nvolatile_tracer_groups), errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(tracer_flx_infl,1) == bounds%endc), errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(tracer_flx_infl,2) == betrtracer_vars%ngwmobile_tracers), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return

   do j = 1, betrtracer_vars%ngwmobile_tracers

     !for a real mechanistic modeling, tracer_flx_infl should be derived from water flux coming from snow melt, surface ponding water,
     !and precipitation. I here just comparomise for a quick shot.

     if(j==betrtracer_vars%id_trc_blk_h2o)then
       do fc = 1, numf
         c = filter(fc)
         tracer_flx_infl(c,j) = 1._r8 * qflx_adv(c,0) * denh2o
       enddo
     elseif(j==betrtracer_vars%id_trc_o18_h2o)then
         do fc = 1, numf
           c = filter(fc)
           tracer_flx_infl(c,j) = 1._r8 * qflx_adv(c,0) * denh2o
         enddo
     elseif(j==betrtracer_vars%id_trc_d_h2o)then
         do fc = 1, numf
           c = filter(fc)
           tracer_flx_infl(c,j) = 1._r8 * qflx_adv(c,0) * denh2o
         enddo
     else
       do fc = 1, numf
         c = filter(fc)
         if(betrtracer_vars%is_volatile(j) .and. betrtracer_vars%is_advective(j))then
           !for volatile non water tracer, infiltration is calculated based dissolution of the gas in the water, this may need
           !improvement when tracers are allowed to transport inside snow, such that the tracer infiltration is derived from mass balance in snow
           tracer_flx_infl(c,j) = bunsencef_topsoi(c,betrtracer_vars%volatilegroupid(j)) * &
                tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1,j) * qflx_adv(c,0)

         else
           tracer_flx_infl(c,j) = 0._r8
         endif

       enddo
     endif
   enddo
   end associate
   end subroutine calc_tracer_infiltration


   !------------------------------------------------------------------------
   subroutine calc_equil_to_liquid_convert_coeff(bounds, lbj, ubj, jtops, numf, filter,&
       t_soisno, h2osoi_ice, dz, betrtracer_vars, aqu2equilsolidcef_col, betr_status)
   !
   ! DESCRIPTION
   ! calculate partition parameter between solid and aqueous phase tracers
   ! this could mean differnt things for different cases
   ! for water isotopes, this represents ice/liquid equilibrium partitioning, it could also
   ! mean linear isotherms of adsorption/desorption. Currently, it is only
   ! for water isotope partitioning between liquid water and ice
   ! It is temporary not needed for water isotopes. Rather an
   ! explicit mass proportional partitioning during freeze-thaw is implemented.
   ! June 7, 2016. Jinyun Tang
   use betr_varcon           , only : denh2o => bdenh2o, denice => bdenice
   use BeTRTracerType        , only : betrtracer_type
   use BetrStatusType        , only : betr_status_type
   implicit none
   type(bounds_type)     , intent(in)    :: bounds  ! bounds
   integer               , intent(in)    :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
   integer               , intent(in)    :: jtops(bounds%begc: )                              ! top label of each column
   integer               , intent(in)    :: numf                                              ! number of columns in column filter
   integer               , intent(in)    :: filter(:)                                         ! column filter
   real(r8)              , intent(in)    :: t_soisno(bounds%begc: , lbj: )
   type(betrtracer_type) , intent(in)    :: betrtracer_vars
   real(r8)              , intent(in)    :: h2osoi_ice(bounds%begc: , lbj: )
   real(r8)              , intent(in)    :: dz(bounds%begc: , lbj: )
   real(r8)              , intent(inout) :: aqu2equilsolidcef_col(bounds%begc:bounds%endc, &
                                              lbj:ubj, 1:betrtracer_vars%nsolid_equil_tracer_groups)
   type(betr_status_type), intent(out)   :: betr_status

   !temporary variables
   real(r8) :: alpha_sl
   integer  :: fc, c, j

   call betr_status%reset()
   return
   SHR_ASSERT_ALL((ubound(t_soisno)   == (/bounds%endc, ubj/)), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(h2osoi_ice) == (/bounds%endc, ubj/)), errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(dz)         == (/bounds%endc, ubj/)), errMsg(filename,__LINE__),betr_status)

    ! remove unused dummy arg compiler warning
    if (numf > 0)                        continue
    if (size(filter) > 0)                continue
    if (size(jtops) > 0)                 continue
    if (size(aqu2equilsolidcef_col) > 0) continue

    associate(                        &
     is_h2o => betrtracer_vars%is_h2o &
   )

!the following code is now not used
!   if(any(is_h2o))then
     !doing a water isotope simulation
!     if(betrtracer_vars%id_trc_o18_h2o>0)then
!        do j = lbj, ubj
!          do fc = 1, numf
!            c = filter(fc)
!            if(j>=jtops(c))then
!              alpha_sl = get_equi_sl_h2oiso_fractionation(betrtracer_vars%id_trc_o18_h2o, t_soisno(c,j), betrtracer_vars)
!              aqu2equilsolidcef_col(c,j, betrtracer_vars%id_trc_o18_h2o_ice) = alpha_sl * h2osoi_ice(c,j) / (denh2o * dz(c,j))
!            endif
!          enddo
!        enddo
!     endif
!   endif

   end associate
   end subroutine calc_equil_to_liquid_convert_coeff


  !-----------------------------------------------------------------------
  subroutine get_zwt (bounds, numf, filter, zi, &
       biophysforc, zwt,jwt, betr_status)
    !
    ! !DESCRIPTION:
    ! Finds the first unsaturated layer going up. Also allows a perched water table over ice.
    !
    use tracer_varcon      , only : nlevsoi  => betr_nlevsoi
    use betr_varcon        , only : tfrz => btfrz
    use BetrStatusType     , only : betr_status_type
   implicit none
    ! !ARGUMENTS:
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: numf                ! number of column soil points in column filter
    integer                          , intent(in)    :: filter(:)           ! column filter for soil points
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    real(r8)                         , intent(in)    :: zi(bounds%begc: , 0: )
    real(r8)                         , intent(inout) :: zwt( bounds%begc: ) ! water table depth (-) [col]
    integer                          , intent(inout) :: jwt(bounds%begc: )
    type(betr_status_type)           , intent(out)   :: betr_status
    !
    ! !LOCAL VARIABLES:
    real(r8) :: f_sat    ! volumetric soil water defining top of water table or where production is allowed
    integer  :: c,j,perch! indices
    integer  :: fc       ! filter column index

    !-----------------------------------------------------------------------
    call betr_status%reset()

    f_sat = 0.95_r8   !a number borrowed from zack's ch4 code

    SHR_ASSERT_ALL((ubound(zwt) == (/bounds%endc/)), errMsg(filename, __LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(zi) == (/bounds%endc, nlevsoi/)), errMsg(filename, __LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(jwt) == (/bounds%endc/)), errMsg(filename, __LINE__),betr_status)
    if(betr_status%check_status())return

    associate(                                       &
         watsat     => biophysforc%watsat_col      , & ! Input:  [real(r8) (:,:)  ] volumetric soil water at saturation (porosity)
         h2osoi_vol => biophysforc%h2osoi_vol_col  , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         t_soisno   => biophysforc%t_soisno_col      & ! Input:  [real(r8) (: ,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
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
  subroutine calc_aerecond(bounds, col, pft, num_soilp, filter_soilp, jwt, &
      biophysforc,  betrtracer_vars,  betr_aerecond_vars, tracercoeff_vars, betr_status)
  !
  ! DESCRIPTION
  !
  ! calculate aerenchyma conductance (m/s)
  !USES
  use betr_varcon        , only : tfrz => btfrz, rpi => brpi
  use BeTR_pftvarconType , only : pftvarcon => betr_pftvarcon
  use BetrTracerType     , only : betrtracer_type
  use BeTR_aerocondType  , only : betr_aerecond_type
  use tracercoeffType    , only : tracercoeff_type
  use tracer_varcon      , only : nlevsoi  => betr_nlevsoi
  use MathfuncMod        , only : safe_div
  use betr_ctrl          , only : betr_use_cn
  use BetrStatusType     , only : betr_status_type
  use betr_columnType    , only : betr_column_type
  use BeTR_patchtype     , only : betr_patch_type
  implicit none
  type(bounds_type)                , intent(in)    :: bounds
  type(betr_column_type)           , intent(in)    :: col
  type(betr_patch_type)            , intent(in)    :: pft
  integer                          , intent(in)    :: num_soilp                 ! number of column soil points in column filter
  integer                          , intent(in)    :: filter_soilp(:)           ! column filter for soil points
  integer                          , intent(in)    :: jwt(bounds%begc: )
  type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
  type(betr_aerecond_type)         , intent(in)    :: betr_aerecond_vars
  type(betrtracer_type)            , intent(in)    :: betrtracer_vars            ! betr configuration information
  type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
  type(betr_status_type)           , intent(out)   :: betr_status

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

  call betr_status%reset()

  SHR_ASSERT_ALL((ubound(jwt) == (/bounds%endc/)), errMsg(filename, __LINE__),betr_status)
  if(betr_status%check_status())return

  associate(                                                                & !
    z                       =>    col%z                               ,     & ! Input:  [real(r8) (:,:)  ]  layer depth (m) (-nlevsno+1:nlevsoi)
    dz                      =>    col%dz                              ,     & ! Input:  [real(r8) (:,:)  ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
    wtcol                   =>    pft%wtcol                           ,     & ! Input:  [real(r8) (:)    ]  weight (relative to column)
    lbl_rsc_h2o             =>    biophysforc%lbl_rsc_h2o_patch  ,          & ! laminar layer resistance for h2o
    elai                    =>    biophysforc%elai_patch         ,          &
    annsum_npp              =>    biophysforc%annsum_npp_patch    ,         & ! Input:  [real(r8) (:) ]  annual sum NPP (gC/m2/yr)
    rootfr                  =>    biophysforc%rootfr_patch        ,         & ! Input: [real(r8) (:,:)]
    annavg_agnpp            =>    betr_aerecond_vars%annavg_agnpp_patch  ,  & ! Output: [real(r8) (:) ]  annual average above-ground NPP (gC/m2/s)
    annavg_bgnpp            =>    betr_aerecond_vars%annavg_bgnpp_patch  ,  & ! Output: [real(r8) (:) ]  annual average below-ground NPP (gC/m2/s)
    frootc                  =>    betr_aerecond_vars%plant_frootsc_patch ,  & ! Input:  [real(r8) (:)    ]  (gC/m2) fine root C
    is_volatile             =>    betrtracer_vars%is_volatile         ,     &
    volatilegroupid         =>    betrtracer_vars%volatilegroupid     ,     &
    ngwmobile_tracer_groups =>  betrtracer_vars%ngwmobile_tracer_groups   , &
    tracer_group_memid      => betrtracer_vars%tracer_group_memid ,         &
    t_veg                   =>    biophysforc%t_veg_patch        ,          &
    t_soisno                =>    biophysforc%t_soisno_col       ,          &
    scal_aere_cond          =>    tracercoeff_vars%scal_aere_cond_col ,     &
    tracer_diffusivity_air  => tracercoeff_vars%tracer_diffusivity_air_col, &
    aere_cond               =>    tracercoeff_vars%aere_cond_col            & !
  )

  do j=1,nlevsoi
    do fp = 1, num_soilp
      p = filter_soilp (fp)
      c = pft%column(p)
      ! Calculate aerenchyma diffusion
      if (j > jwt(c) .and. t_soisno(c,j) > tfrz .and. pftvarcon%is_grass_patch(pft%itype(p))) then
        ! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
        ! the new components are in; if there is any tuning to be done to get a realistic global flux,
        ! this would probably be the place.  We will have to document clearly in the Tech Note
        ! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)
        if(betr_use_cn)then
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

        if(.not. betr_use_cn)then
          m_tiller = 0._r8        !this was set to zero purposely
        else
          if (usefrootc) then
            m_tiller = frootc(p) ! This will yield much smaller aere area.
          else
            m_tiller = anpp * nppratio * elai(p)
          end if
        endif
        n_tiller = m_tiller / 0.22_r8

        if (pftvarcon%is_grass_patch(pft%itype(p)) .or. pft%crop(pft%itype(p)) == 1) then
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
  subroutine betr_annualupdate(betr_time, bounds, pft, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       biophysforc, betr_aerecond_vars, tracercoeff_vars, betr_status)
    !
    ! !DESCRIPTION: Annual mean fields.
    !
    ! !USES:
    use BeTR_TimeMod      , only : betr_time_type
    use betr_varcon       , only : secspday => bsecspday
    use tracercoeffType   , only : tracercoeff_type
    use BeTR_aerocondType , only : betr_aerecond_type
    use BetrStatusType    , only : betr_status_type
    use BeTR_patchtype    , only : betr_patch_type
    !
    ! !ARGUMENTS:
    type(betr_time_type)             , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds
    type(betr_patch_type)            , intent(in)    :: pft
    integer                          , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                          , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                          , intent(in)    :: num_soilp         ! number of soil points in pft filter
    integer                          , intent(in)    :: filter_soilp(:)   ! patch filter for soil points
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(betr_aerecond_type)         , intent(inout) :: betr_aerecond_vars
    type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
    type(betr_status_type)           , intent(out)   :: betr_status
    !
    ! !LOCAL VARIABLES:
    integer :: c,p       ! indices
    integer :: fc        ! soil column filter indices
    integer :: fp        ! soil pft filter indices
    real(r8):: dt        ! time step (seconds)
    real(r8):: secsperyear
    logical :: newrun
    !-----------------------------------------------------------------------

    ! remove unused dummy arg compiler warning
    if (bounds%begc > 0) continue
    call betr_status%reset()

    associate(                                                               &
         agnpp           =>    biophysforc%agnpp_patch                    ,  & ! Input:  [real(r8) (:) ]  (gC/m2/s) aboveground NPP
         bgnpp           =>    biophysforc%bgnpp_patch                    ,  & ! Input:  [real(r8) (:) ]  (gC/m2/s) belowground NPP
         tempavg_agnpp   =>    betr_aerecond_vars%tempavg_agnpp_patch      , & ! Output: [real(r8) (:) ]  temporary average above-ground NPP (gC/m2/s)
         annavg_agnpp    =>    betr_aerecond_vars%annavg_agnpp_patch       , & ! Output: [real(r8) (:) ]  annual average above-ground NPP (gC/m2/s)
         tempavg_bgnpp   =>    betr_aerecond_vars%tempavg_bgnpp_patch      , & ! Output: [real(r8) (:) ]  temporary average below-ground NPP (gC/m2/s)
         annavg_bgnpp    =>    betr_aerecond_vars%annavg_bgnpp_patch       , & ! Output: [real(r8) (:) ]  annual average below-ground NPP (gC/m2/s)
         annsum_counter  =>    tracercoeff_vars%annsum_counter_col           & ! Output: [real(r8) (:) ]  seconds since last annual accumulator turnover
         )

      ! set time steps
      dt = betr_time%get_step_size()
      secsperyear = real( betr_time%get_days_per_year() * secspday, r8)

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


!-------------------------------------------------------------------------------

   subroutine calc_rhovap(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
     t_soisno, smp_l, rho_vap, betr_status)
   !
   !DESCRIPTION
   !Compute actual vapor pressure inside the soil profile
   !
   ! uses
   use betr_varcon   , only : rwat => brwat, grav => bgrav
   use BetrStatusType, only : betr_status_type
   implicit none
   !arguments
   type(bounds_type) , intent(in)    :: bounds                           ! bounds
   integer           , intent(in)    :: lbj, ubj                         ! lower and upper bounds, make sure they are > 0
   integer           , intent(in)    :: jtops(bounds%begc: )             ! top label of each column
   integer           , intent(in)    :: num_soilc                        ! number of soil filters
   integer           , intent(in)    :: filter_soilc(:)                  ! filter
   real(r8)          , intent(in)    :: t_soisno(bounds%begc: , lbj: )   !soil temperature, K
   real(r8)          , intent(in)    :: smp_l(bounds%begc: , lbj: )      !liquid soil matric potential, mm
   real(r8)          , intent(inout) :: rho_vap(bounds%begc: , lbj: )  !actual vapor pressure, kg/m3
   type(betr_status_type),intent(out):: betr_status

   !local variables
   real(r8)           :: hh           !relative humidity
   real(r8)           :: rho_sat      !saturated water vapor pressure
   integer            :: c, fc, n
   character(len=255) :: subname ='calc_rhovap'
   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)        == (/bounds%endc/)),        errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(t_soisno,1)     == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(t_soisno,2)     == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(smp_l,1)        == bounds%endc),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(smp_l,2)        == ubj),   errMsg(filename,__LINE__),betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(rho_vap,1)      == bounds%begc),   errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   SHR_ASSERT_ALL((ubound(rho_vap,2)      == ubj),   errMsg(filename,__LINE__), betr_status)
   if(betr_status%check_status())return

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

end module TracerParamsMod
