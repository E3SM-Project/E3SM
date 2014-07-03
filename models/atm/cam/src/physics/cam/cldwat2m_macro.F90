
  module cldwat2m_macro

  !--------------------------------------------------- !
  ! Purpose     : CAM Interface for Cloud Macrophysics !
  ! Author      : Sungsu Park                          !
  ! Description : Park et al. 2010.                    !
  ! For questions, contact Sungsu Park                 !
  !                        e-mail : sungsup@ucar.edu   !
  !                        phone  : 303-497-1375       !
  !--------------------------------------------------- !

   use shr_kind_mod,     only: r8=>shr_kind_r8
   use spmd_utils,       only: masterproc
   use ppgrid,           only: pcols, pver, pverp
   use abortutils,       only: endrun
   use physconst,        only: cpair, latvap, latice, rh2o
   use wv_saturation,    only: qsat_water, svp_water, svp_ice
   use cam_history,      only: addfld, add_default, phys_decomp, outfld 
   use cam_logfile,      only: iulog
   use ref_pres,         only: top_lev=>trop_cloud_top_lev

   implicit none
   private
   save

   public :: ini_macro, mmacro_pcond, aist_vector, aist_single, astG_RHU_single, astG_PDF_single, CAMstfrac

   ! -------------- !
   ! Set Parameters !
   ! -------------- !

   ! -------------------------- !
   ! Parameters for Ice Stratus !
   ! -------------------------- !
   real(r8), public,  parameter :: rhmini       = 0.80_r8       ! Minimum rh for ice cloud fraction > 0.
   real(r8), public,  parameter :: rhmaxi       = 1.1_r8        ! rhi at which ice cloud fraction = 1.
   real(r8), parameter          :: qist_min     = 1.e-7_r8      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
   real(r8), parameter          :: qist_max     = 5.e-3_r8      ! Maximum in-stratus ice IWC constraint [ kg/kg ]

   ! ----------------------------- !
   ! Parameters for Liquid Stratus !
   ! ----------------------------- !

   logical,  parameter          :: CAMstfrac    = .false.       ! If .true. (.false.),
                                                                ! use Slingo (triangular PDF-based) liquid stratus fraction
   logical,  parameter          :: freeze_dry   = .false.       ! If .true., use 'freeze dry' in liquid stratus fraction formula
   real(r8), parameter          :: qlst_min     = 2.e-5_r8      ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: qlst_max     = 3.e-3_r8      ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: cc           = 0.1_r8        ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
   integer,  parameter          :: niter        = 2             ! For iterative computation of QQ with 'ramda' below.
   real(r8), parameter          :: ramda        = 0.5_r8        ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
   real(r8), private            :: rhminl                       ! Critical RH for low-level  liquid stratus clouds
   real(r8), private            :: rhminl_adj_land              ! rhminl adjustment for snowfree land
   real(r8), private            :: rhminh                       ! Critical RH for high-level liquid stratus clouds
   real(r8), private            :: premit                       ! Top    height for mid-level liquid stratus fraction
   real(r8), private            :: premib                       ! Bottom height for mid-level liquid stratus fraction
   integer,  private            :: iceopt                       ! option for ice cloud closure 
                                                                ! 1=wang & sassen 2=schiller (iciwc)  
                                                                ! 3=wood & field, 4=Wilson (based on smith)
                                                                ! 5=modified slingo (ssat & empyt cloud)        
   real(r8), private            :: icecrit                      ! Critical RH for ice clouds in Wilson & Ballard closure
                                                                ! ( smaller = more ice clouds )

   contains

   ! -------------- !
   ! Initialization !
   ! -------------- !

   subroutine ini_macro

   !--------------------------------------------------------------------- !
   !                                                                      ! 
   ! Purpose: Initialize constants for the liquid stratiform macrophysics !
   !          called from stratiform.F90.                                 !  
   ! Author:  Sungsu Park, Dec.01.2009.                                   !
   !                                                                      !
   !--------------------------------------------------------------------- !

   use cloud_fraction, only: cldfrc_getparams

   call cldfrc_getparams(rhminl_out = rhminl, rhminl_adj_land_out = rhminl_adj_land,    &
                         rhminh_out = rhminh, premit_out = premit,   premib_out = premib, &
                         iceopt_out = iceopt, icecrit_out = icecrit)

   if( masterproc ) then
       write(iulog,*) 'ini_macro: tuning parameters : rhminl = ', rhminl, &
                                                     'rhminl_adj_land = ', rhminl_adj_land, & 
                                                     'rhminh = ', rhminh, & 
                                                     'premit = ', premit, & 
                                                     'premib = ',  premib,  &
                                                     'iceopt = ',  iceopt,  &
                                                     'icecrit = ', icecrit
   end if

   return
   end subroutine ini_macro

   ! ------------------------------ !
   ! Stratiform Liquid Macrophysics !
   ! ------------------------------ !

   ! In the version, 'macro --> micro --> advective forcing --> macro...'
   ! A_...: only 'advective forcing' without 'microphysical tendency'
   ! C_...: only 'microphysical tendency'
   ! D_...: only 'detrainment of cumulus condensate'  
   ! So, 'A' and 'C' are exclusive. 

   subroutine mmacro_pcond( lchnk      , ncol       , dt         , p            , dp         ,              &
                            T0         , qv0        , ql0        , qi0          , nl0        , ni0        , &
                            A_T        , A_qv       , A_ql       , A_qi         , A_nl       , A_ni       , &
                            C_T        , C_qv       , C_ql       , C_qi         , C_nl       , C_ni       , C_qlst     , &
                            D_T        , D_qv       , D_ql       , D_qi         , D_nl       , D_ni       , &
                            a_cud      , a_cu0      , landfrac   , snowh        ,                           & 
                            s_tendout  , qv_tendout , ql_tendout , qi_tendout   , nl_tendout , ni_tendout , &
                            qme        , qvadj      , qladj      , qiadj        , qllim      , qilim      , &
                            cld        , al_st_star , ai_st_star , ql_st_star   , qi_st_star , do_cldice  )

   use constituents,     only : qmin, cnst_get_ind
   use time_manager,     only : is_first_step, get_nstep
   use wv_saturation,    only : findsp_vc

   implicit none

   integer   icol
   integer,  intent(in)    :: lchnk                        ! Chunk number
   integer,  intent(in)    :: ncol                         ! Number of active columns

   ! Input-Output variables

   real(r8), intent(inout) :: T0(pcols,pver)               ! Temperature [K]
   real(r8), intent(inout) :: qv0(pcols,pver)              ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(inout) :: ql0(pcols,pver)              ! Grid-mean liquid water content [kg/kg]
   real(r8), intent(inout) :: qi0(pcols,pver)              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(inout) :: nl0(pcols,pver)              ! Grid-mean number concentration of cloud liquid droplet [#/kg]
   real(r8), intent(inout) :: ni0(pcols,pver)              ! Grid-mean number concentration of cloud ice    droplet [#/kg]

   ! Input variables

   real(r8), intent(in)    :: dt                           ! Model integration time step [s]
   real(r8), intent(in)    :: p(pcols,pver)                ! Pressure at the layer mid-point [Pa]
   real(r8), intent(in)    :: dp(pcols,pver)               ! Pressure thickness [Pa] > 0

   real(r8), intent(in)    :: A_T(pcols,pver)              ! Non-microphysical advective external forcing of T  [K/s]
   real(r8), intent(in)    :: A_qv(pcols,pver)             ! Non-microphysical advective external forcing of qv [kg/kg/s]
   real(r8), intent(in)    :: A_ql(pcols,pver)             ! Non-microphysical advective external forcing of ql [kg/kg/s]
   real(r8), intent(in)    :: A_qi(pcols,pver)             ! Non-microphysical advective external forcing of qi [kg/kg/s]
   real(r8), intent(in)    :: A_nl(pcols,pver)             ! Non-microphysical advective external forcing of nl [#/kg/s]
   real(r8), intent(in)    :: A_ni(pcols,pver)             ! Non-microphysical advective external forcing of ni [#/kg/s] 

   real(r8), intent(in)    :: C_T(pcols,pver)              ! Microphysical advective external forcing of T  [K/s]
   real(r8), intent(in)    :: C_qv(pcols,pver)             ! Microphysical advective external forcing of qv [kg/kg/s]
   real(r8), intent(in)    :: C_ql(pcols,pver)             ! Microphysical advective external forcing of ql [kg/kg/s]
   real(r8), intent(in)    :: C_qi(pcols,pver)             ! Microphysical advective external forcing of qi [kg/kg/s]
   real(r8), intent(in)    :: C_nl(pcols,pver)             ! Microphysical advective external forcing of nl [#/kg/s]
   real(r8), intent(in)    :: C_ni(pcols,pver)             ! Microphysical advective external forcing of ni [#/kg/s] 
   real(r8), intent(in)    :: C_qlst(pcols,pver)           ! Microphysical advective external forcing of ql
                                                           ! within liquid stratus [kg/kg/s]

   real(r8), intent(in)    :: D_T(pcols,pver)              ! Cumulus detrainment external forcing of T  [K/s]
   real(r8), intent(in)    :: D_qv(pcols,pver)             ! Cumulus detrainment external forcing of qv [kg/kg/s]
   real(r8), intent(in)    :: D_ql(pcols,pver)             ! Cumulus detrainment external forcing of ql [kg/kg/s]
   real(r8), intent(in)    :: D_qi(pcols,pver)             ! Cumulus detrainment external forcing of qi [kg/kg/s]
   real(r8), intent(in)    :: D_nl(pcols,pver)             ! Cumulus detrainment external forcing of nl [#/kg/s]
   real(r8), intent(in)    :: D_ni(pcols,pver)             ! Cumulus detrainment external forcing of qi [#/kg/s] 

   real(r8), intent(in)    :: a_cud(pcols,pver)            ! Old cumulus fraction before update
   real(r8), intent(in)    :: a_cu0(pcols,pver)            ! New cumulus fraction after update

   real(r8), intent(in)    :: landfrac(pcols)              ! Land fraction
   real(r8), intent(in)    :: snowh(pcols)                 ! Snow depth (liquid water equivalent)
   logical, intent(in)     :: do_cldice                    ! Whether or not cldice should be prognosed

   ! Output variables

   real(r8), intent(out)   :: s_tendout(pcols,pver)        ! Net tendency of grid-mean s  from 'Micro+Macro' processes [J/kg/s]
   real(r8), intent(out)   :: qv_tendout(pcols,pver)       ! Net tendency of grid-mean qv from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: ql_tendout(pcols,pver)       ! Net tendency of grid-mean ql from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: qi_tendout(pcols,pver)       ! Net tendency of grid-mean qi from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: nl_tendout(pcols,pver)       ! Net tendency of grid-mean nl from 'Micro+Macro' processes [#/kg/s]
   real(r8), intent(out)   :: ni_tendout(pcols,pver)       ! Net tendency of grid-mean ni from 'Micro+Macro' processes [#/kg/s]

   real(r8), intent(out)   :: qme  (pcols,pver)            ! Net condensation rate [kg/kg/s]
   real(r8), intent(out)   :: qvadj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (vapor)
   real(r8), intent(out)   :: qladj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (liquid)
   real(r8), intent(out)   :: qiadj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (ice)
   real(r8), intent(out)   :: qllim(pcols,pver)            ! tendency from "instratus_condensate" call (liquid)
   real(r8), intent(out)   :: qilim(pcols,pver)            ! tendency from "instratus_condensate" call (ice)

   real(r8), intent(out)   :: cld(pcols,pver)              ! Net cloud fraction ( 0 <= cld <= 1 )
   real(r8), intent(out)   :: al_st_star(pcols,pver)       ! Physical liquid stratus fraction
   real(r8), intent(out)   :: ai_st_star(pcols,pver)       ! Physical ice stratus fraction
   real(r8), intent(out)   :: ql_st_star(pcols,pver)       ! In-stratus LWC [kg/kg] 
   real(r8), intent(out)   :: qi_st_star(pcols,pver)       ! In-stratus IWC [kg/kg] 

   ! --------------- !
   ! Local variables !
   ! --------------- !
   integer :: ixcldliq, ixcldice
 
   integer :: i, j, k, iter, ii, jj                        ! Loop indexes

   ! Thermodynamic state variables

   real(r8) T(pcols,pver)                                  ! Temperature of equilibrium reference state
                                                           ! from which 'Micro & Macro' are computed [K]
   real(r8) T1(pcols,pver)                                 ! Temperature after 'fice_force' on T01  
   real(r8) T_0(pcols,pver)                                ! Temperature after 'instratus_condensate' on T1
   real(r8) T_05(pcols,pver)                               ! Temperature after 'advection' on T_0 
   real(r8) T_prime0(pcols,pver)                           ! Temperature after 'Macrophysics (QQ)' on T_05star
   real(r8) T_dprime(pcols,pver)                           ! Temperature after 'fice_force' on T_prime
   real(r8) T_star(pcols,pver)                             ! Temperature after 'instratus_condensate' on T_dprime

   real(r8) qv(pcols,pver)                                 ! Grid-mean qv of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qv1(pcols,pver)                                ! Grid-mean qv after 'fice_force' on qv01  
   real(r8) qv_0(pcols,pver)                               ! Grid-mean qv after 'instratus_condensate' on qv1
   real(r8) qv_05(pcols,pver)                              ! Grid-mean qv after 'advection' on qv_0 
   real(r8) qv_prime0(pcols,pver)                          ! Grid-mean qv after 'Macrophysics (QQ)' on qv_05star
   real(r8) qv_dprime(pcols,pver)                          ! Grid-mean qv after 'fice_force' on qv_prime
   real(r8) qv_star(pcols,pver)                            ! Grid-mean qv after 'instratus_condensate' on qv_dprime

   real(r8) ql(pcols,pver)                                 ! Grid-mean ql of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ql1(pcols,pver)                                ! Grid-mean ql after 'fice_force' on ql01  
   real(r8) ql_0(pcols,pver)                               ! Grid-mean ql after 'instratus_condensate' on ql1
   real(r8) ql_05(pcols,pver)                              ! Grid-mean ql after 'advection' on ql_0 
   real(r8) ql_prime0(pcols,pver)                          ! Grid-mean ql after 'Macrophysics (QQ)' on ql_05star
   real(r8) ql_dprime(pcols,pver)                          ! Grid-mean ql after 'fice_force' on ql_prime
   real(r8) ql_star(pcols,pver)                            ! Grid-mean ql after 'instratus_condensate' on ql_dprime

   real(r8) qi(pcols,pver)                                 ! Grid-mean qi of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qi1(pcols,pver)                                ! Grid-mean qi after 'fice_force' on qi01  
   real(r8) qi_0(pcols,pver)                               ! Grid-mean qi after 'instratus_condensate' on qi1
   real(r8) qi_05(pcols,pver)                              ! Grid-mean qi after 'advection' on qi_0 
   real(r8) qi_prime0(pcols,pver)                          ! Grid-mean qi after 'Macrophysics (QQ)' on qi_05star
   real(r8) qi_dprime(pcols,pver)                          ! Grid-mean qi after 'fice_force' on qi_prime
   real(r8) qi_star(pcols,pver)                            ! Grid-mean qi after 'instratus_condensate' on qi_dprime

   real(r8) nl(pcols,pver)                                 ! Grid-mean nl of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) nl1(pcols,pver)                                ! Grid-mean nl after 'fice_force' on nl01  
   real(r8) nl_0(pcols,pver)                               ! Grid-mean nl after 'instratus_condensate' on nl1
   real(r8) nl_05(pcols,pver)                              ! Grid-mean nl after 'advection' on nl_0 
   real(r8) nl_prime0(pcols,pver)                          ! Grid-mean nl after 'Macrophysics (QQ)' on nl_05star
   real(r8) nl_dprime(pcols,pver)                          ! Grid-mean nl after 'fice_force' on nl_prime
   real(r8) nl_star(pcols,pver)                            ! Grid-mean nl after 'instratus_condensate' on nl_dprime

   real(r8) ni(pcols,pver)                                 ! Grid-mean ni of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ni1(pcols,pver)                                ! Grid-mean ni after 'fice_force' on ni01  
   real(r8) ni_0(pcols,pver)                               ! Grid-mean ni after 'instratus_condensate' on ni1
   real(r8) ni_05(pcols,pver)                              ! Grid-mean ni after 'advection' on ni_0 
   real(r8) ni_prime0(pcols,pver)                          ! Grid-mean ni after 'Macrophysics (QQ)' on ni_05star
   real(r8) ni_dprime(pcols,pver)                          ! Grid-mean ni after 'fice_force' on ni_prime
   real(r8) ni_star(pcols,pver)                            ! Grid-mean ni after 'instratus_condensate' on ni_dprime

   real(r8) a_st(pcols,pver)                               ! Stratus fraction of equilibrium reference state 
   real(r8) a_st_0(pcols,pver)                             ! Stratus fraction at '_0' state
   real(r8) a_st_star(pcols,pver)                          ! Stratus fraction at '_star' state

   real(r8) al_st(pcols,pver)                              ! Liquid stratus fraction of equilibrium reference state 
   real(r8) al_st_0(pcols,pver)                            ! Liquid stratus fraction at '_0' state
   real(r8) al_st_nc(pcols,pver)                           ! Non-physical liquid stratus fraction in the non-cumulus pixels

   real(r8) ai_st(pcols,pver)                              ! Ice stratus fraction of equilibrium reference state 
   real(r8) ai_st_0(pcols,pver)                            ! Ice stratus fraction at '_0' state
   real(r8) ai_st_nc(pcols,pver)                           ! Non-physical ice stratus fraction in the non-cumulus pixels

   real(r8) ql_st(pcols,pver)                              ! In-stratus LWC of equilibrium reference state [kg/kg] 
   real(r8) ql_st_0(pcols,pver)                            ! In-stratus LWC at '_0' state

   real(r8) qi_st(pcols,pver)                              ! In-stratus IWC of equilibrium reference state [kg/kg] 
   real(r8) qi_st_0(pcols,pver)                            ! In-stratus IWC at '_0' state

 ! Cumulus properties 

   real(r8) dacudt(pcols,pver)
   real(r8) a_cu(pcols,pver)

 ! Adjustment tendency in association with 'positive_moisture'

   real(r8) Tten_pwi1(pcols,pver)                          ! Pre-process T  tendency of input equilibrium state [K/s] 
   real(r8) qvten_pwi1(pcols,pver)                         ! Pre-process qv tendency of input equilibrium state [kg/kg/s]
   real(r8) qlten_pwi1(pcols,pver)                         ! Pre-process ql tendency of input equilibrium state [kg/kg/s]
   real(r8) qiten_pwi1(pcols,pver)                         ! Pre-process qi tendency of input equilibrium state [kg/kg/s]
   real(r8) nlten_pwi1(pcols,pver)                         ! Pre-process nl tendency of input equilibrium state [#/kg/s]
   real(r8) niten_pwi1(pcols,pver)                         ! Pre-process ni tendency of input equilibrium state [#/kg/s] 

   real(r8) Tten_pwi2(pcols,pver)                          ! Post-process T  tendency of provisional equilibrium state [K/s] 
   real(r8) qvten_pwi2(pcols,pver)                         ! Post-process qv tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qlten_pwi2(pcols,pver)                         ! Post-process ql tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qiten_pwi2(pcols,pver)                         ! Post-process qi tendency of provisional equilibrium state [kg/kg/s]
   real(r8) nlten_pwi2(pcols,pver)                         ! Post-process nl tendency of provisoonal equilibrium state [#/kg/s]
   real(r8) niten_pwi2(pcols,pver)                         ! Post-process ni tendency of provisional equilibrium state [#/kg/s] 

   real(r8) A_T_adj(pcols,pver)                            ! After applying external advective forcing [K/s]
   real(r8) A_qv_adj(pcols,pver)                           ! After applying external advective forcing [kg/kg/s]
   real(r8) A_ql_adj(pcols,pver)                           ! After applying external advective forcing [kg/kg/s]
   real(r8) A_qi_adj(pcols,pver)                           ! After applying external advective forcing [kg/kg/s]
   real(r8) A_nl_adj(pcols,pver)                           ! After applying external advective forcing [#/kg/s]
   real(r8) A_ni_adj(pcols,pver)                           ! After applying external advective forcing [#/kg/s]

 ! Adjustment tendency in association with 'instratus_condensate'

   real(r8) QQw1(pcols,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQi1(pcols,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQw2(pcols,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQi2(pcols,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]

   real(r8) QQnl1(pcols,pver)          ! Tendency of nl associated with QQw1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(r8) QQni1(pcols,pver)          ! Tendency of ni associated with QQi1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(r8) QQnl2(pcols,pver)          ! Tendency of nl associated with QQw2 only when QQw2<0 (net evaporation) [#/kg/s]
   real(r8) QQni2(pcols,pver)          ! Tendency of ni associated with QQi2 only when QQw2<0 (net evaporation) [#/kg/s]

 ! Macrophysical process tendency variables

   real(r8) QQ(pcols,pver)             ! Net condensation rate into water+ice           [kg/kg/s] 
   real(r8) QQw(pcols,pver)            ! Net condensation rate into water               [kg/kg/s] 
   real(r8) QQi(pcols,pver)            ! Net condensation rate into ice                 [kg/kg/s]
   real(r8) QQnl(pcols,pver)           ! Tendency of nl associated with QQw both for condensation and evaporation [#/kg/s]
   real(r8) QQni(pcols,pver)           ! Tendency of ni associated with QQi both for condensation and evaporation [#/kg/s]
   real(r8) ACnl(pcols,pver)           ! Cloud liquid droplet (nl) activation tendency [#/kg/s]
   real(r8) ACni(pcols,pver)           ! Cloud ice    droplet (ni) activation tendency [#/kg/s]

   real(r8) QQw_prev(pcols,pver)   
   real(r8) QQi_prev(pcols,pver)   
   real(r8) QQnl_prev(pcols,pver)  
   real(r8) QQni_prev(pcols,pver)  

   real(r8) QQw_prog(pcols,pver)   
   real(r8) QQi_prog(pcols,pver)   
   real(r8) QQnl_prog(pcols,pver)  
   real(r8) QQni_prog(pcols,pver)  

   real(r8) QQ_final(pcols,pver)                           
   real(r8) QQw_final(pcols,pver)                           
   real(r8) QQi_final(pcols,pver)                           
   real(r8) QQn_final(pcols,pver)                           
   real(r8) QQnl_final(pcols,pver)                          
   real(r8) QQni_final(pcols,pver)                          

   real(r8) QQ_all(pcols,pver)         ! QQw_all    + QQi_all
   real(r8) QQw_all(pcols,pver)        ! QQw_final  + QQw1  + QQw2  + qlten_pwi1 + qlten_pwi2 + A_ql_adj [kg/kg/s]
   real(r8) QQi_all(pcols,pver)        ! QQi_final  + QQi1  + QQi2  + qiten_pwi1 + qiten_pwi2 + A_qi_adj [kg/kg/s]
   real(r8) QQn_all(pcols,pver)        ! QQnl_all   + QQni_all
   real(r8) QQnl_all(pcols,pver)       ! QQnl_final + QQnl1 + QQnl2 + nlten_pwi1 + nlten_pwi2 + ACnl [#/kg/s]
   real(r8) QQni_all(pcols,pver)       ! QQni_final + QQni1 + QQni2 + niten_pwi1 + niten_pwi2 + ACni [#/kg/s]

 ! Coefficient for computing QQ and related processes

   real(r8) U(pcols,pver)                                  ! Grid-mean RH
   real(r8) U_nc(pcols,pver)                               ! Mean RH of non-cumulus pixels
   real(r8) G_nc(pcols,pver)                               ! d(U_nc)/d(a_st_nc)
   real(r8) F_nc(pcols,pver)                               ! A function of second parameter for a_st_nc
   real(r8) alpha                                          ! = 1/qs
   real(r8) beta                                           ! = (qv/qs**2)*dqsdT
   real(r8) betast                                         ! = alpha*dqsdT
   real(r8) gammal                                         ! = alpha + (latvap/cpair)*beta
   real(r8) gammai                                         ! = alpha + ((latvap+latice)/cpair)*beta
   real(r8) gammaQ                                         ! = alpha + (latvap/cpair)*beta
   real(r8) deltal                                         ! = 1 + a_st*(latvap/cpair)*(betast/alpha)
   real(r8) deltai                                         ! = 1 + a_st*((latvap+latice)/cpair)*(betast/alpha)
   real(r8) A_Tc                                           ! Advective external forcing of Tc [K/s]
   real(r8) A_qt                                           ! Advective external forcing of qt [kg/kg/s]
   real(r8) C_Tc                                           ! Microphysical forcing of Tc [K/s]
   real(r8) C_qt                                           ! Microphysical forcing of qt [kg/kg/s]
   real(r8) dTcdt                                          ! d(Tc)/dt      [K/s]
   real(r8) dqtdt                                          ! d(qt)/dt      [kg/kg/s]
   real(r8) dqtstldt                                       ! d(qt_alst)/dt [kg/kg/s]
   real(r8) dqidt                                          ! d(qi)/dt      [kg/kg/s]

   real(r8) dqlstdt                                        ! d(ql_st)/dt [kg/kg/s]
   real(r8) dalstdt                                        ! d(al_st)/dt  [1/s]
   real(r8) dastdt                                         ! d(a_st)/dt  [1/s]

   real(r8) anic                                           ! Fractional area of non-cumulus and non-ice stratus fraction
   real(r8) GG                                             ! G_nc(i,k)/anic

   real(r8) aa(2,2)
   real(r8) bb(2,1)

   real(r8) zeros(pcols,pver)

   real(r8) qmin1(pcols,pver)
   real(r8) qmin2(pcols,pver)
   real(r8) qmin3(pcols,pver)

   real(r8) esat_a(pcols)                                  ! Saturation water vapor pressure [Pa]
   real(r8) qsat_a(pcols)                                  ! Saturation water vapor specific humidity [kg/kg]
   real(r8) dqsdT_a(pcols)                                 ! dqsat/dT [kg/kg/K]
   real(r8) Twb_aw(pcols,pver)                             ! Wet-bulb temperature [K]
   real(r8) qvwb_aw(pcols,pver)                            ! Wet-bulb water vapor specific humidity [kg/kg]

   real(r8) esat_b(pcols)                                 
   real(r8) qsat_b(pcols)                                 
   real(r8) dqsdT_b(pcols)                                 

   real(r8) QQmax,QQmin,QQwmin,QQimin                      ! For limiting QQ
   real(r8) cone                                           ! Number close to but smaller than 1
   real(r8) qsmall                                         ! Smallest mixing ratio considered in the macrophysics

   cone            = 0.999_r8
   qsmall          = 1.e-18_r8
   zeros(:ncol,:)  = 0._r8

   ! ------------------------------------ !
   ! Global initialization of main output !
   ! ------------------------------------ !

     s_tendout(:ncol,:)     = 0._r8
     qv_tendout(:ncol,:)    = 0._r8
     ql_tendout(:ncol,:)    = 0._r8
     qi_tendout(:ncol,:)    = 0._r8
     nl_tendout(:ncol,:)    = 0._r8
     ni_tendout(:ncol,:)    = 0._r8

     qme(:ncol,:)           = 0._r8

     cld(:ncol,:)           = 0._r8
     al_st_star(:ncol,:)    = 0._r8
     ai_st_star(:ncol,:)    = 0._r8
     ql_st_star(:ncol,:)    = 0._r8
     qi_st_star(:ncol,:)    = 0._r8

   ! --------------------------------------- !
   ! Initialization of internal 2D variables !
   ! --------------------------------------- !

     T(:ncol,:)             = 0._r8
     T1(:ncol,:)            = 0._r8
     T_0(:ncol,:)           = 0._r8
     T_05(:ncol,:)          = 0._r8
     T_prime0(:ncol,:)      = 0._r8
     T_dprime(:ncol,:)      = 0._r8
     T_star(:ncol,:)        = 0._r8

     qv(:ncol,:)            = 0._r8
     qv1(:ncol,:)           = 0._r8
     qv_0(:ncol,:)          = 0._r8
     qv_05(:ncol,:)         = 0._r8
     qv_prime0(:ncol,:)     = 0._r8
     qv_dprime(:ncol,:)     = 0._r8
     qv_star(:ncol,:)       = 0._r8

     ql(:ncol,:)            = 0._r8
     ql1(:ncol,:)           = 0._r8
     ql_0(:ncol,:)          = 0._r8
     ql_05(:ncol,:)         = 0._r8
     ql_prime0(:ncol,:)     = 0._r8
     ql_dprime(:ncol,:)     = 0._r8
     ql_star(:ncol,:)       = 0._r8

     qi(:ncol,:)            = 0._r8
     qi1(:ncol,:)           = 0._r8
     qi_0(:ncol,:)          = 0._r8
     qi_05(:ncol,:)         = 0._r8
     qi_prime0(:ncol,:)     = 0._r8
     qi_dprime(:ncol,:)     = 0._r8
     qi_star(:ncol,:)       = 0._r8

     nl(:ncol,:)            = 0._r8
     nl1(:ncol,:)           = 0._r8
     nl_0(:ncol,:)          = 0._r8
     nl_05(:ncol,:)         = 0._r8
     nl_prime0(:ncol,:)     = 0._r8
     nl_dprime(:ncol,:)     = 0._r8
     nl_star(:ncol,:)       = 0._r8

     ni(:ncol,:)            = 0._r8
     ni1(:ncol,:)           = 0._r8
     ni_0(:ncol,:)          = 0._r8
     ni_05(:ncol,:)         = 0._r8
     ni_prime0(:ncol,:)     = 0._r8
     ni_dprime(:ncol,:)     = 0._r8
     ni_star(:ncol,:)       = 0._r8

     a_st(:ncol,:)          = 0._r8
     a_st_0(:ncol,:)        = 0._r8
     a_st_star(:ncol,:)     = 0._r8

     al_st(:ncol,:)         = 0._r8
     al_st_0(:ncol,:)       = 0._r8
     al_st_nc(:ncol,:)      = 0._r8

     ai_st(:ncol,:)         = 0._r8
     ai_st_0(:ncol,:)       = 0._r8
     ai_st_nc(:ncol,:)      = 0._r8

     ql_st(:ncol,:)         = 0._r8
     ql_st_0(:ncol,:)       = 0._r8

     qi_st(:ncol,:)         = 0._r8
     qi_st_0(:ncol,:)       = 0._r8

 ! Cumulus properties 

     dacudt(:ncol,:)        = 0._r8
     a_cu(:ncol,:)          = 0._r8

 ! Adjustment tendency in association with 'positive_moisture'

     Tten_pwi1(:ncol,:)     = 0._r8
     qvten_pwi1(:ncol,:)    = 0._r8
     qlten_pwi1(:ncol,:)    = 0._r8
     qiten_pwi1(:ncol,:)    = 0._r8
     nlten_pwi1(:ncol,:)    = 0._r8
     niten_pwi1(:ncol,:)    = 0._r8

     Tten_pwi2(:ncol,:)     = 0._r8
     qvten_pwi2(:ncol,:)    = 0._r8
     qlten_pwi2(:ncol,:)    = 0._r8
     qiten_pwi2(:ncol,:)    = 0._r8
     nlten_pwi2(:ncol,:)    = 0._r8
     niten_pwi2(:ncol,:)    = 0._r8

     A_T_adj(:ncol,:)       = 0._r8
     A_qv_adj(:ncol,:)      = 0._r8
     A_ql_adj(:ncol,:)      = 0._r8
     A_qi_adj(:ncol,:)      = 0._r8
     A_nl_adj(:ncol,:)      = 0._r8
     A_ni_adj(:ncol,:)      = 0._r8

     qvadj   (:ncol,:)      = 0._r8
     qladj   (:ncol,:)      = 0._r8
     qiadj   (:ncol,:)      = 0._r8

 ! Adjustment tendency in association with 'instratus_condensate'

     QQw1(:ncol,:)          = 0._r8
     QQi1(:ncol,:)          = 0._r8
     QQw2(:ncol,:)          = 0._r8
     QQi2(:ncol,:)          = 0._r8

     QQnl1(:ncol,:)         = 0._r8
     QQni1(:ncol,:)         = 0._r8
     QQnl2(:ncol,:)         = 0._r8
     QQni2(:ncol,:)         = 0._r8

     QQnl(:ncol,:)          = 0._r8
     QQni(:ncol,:)          = 0._r8

 ! Macrophysical process tendency variables

     QQ(:ncol,:)            = 0._r8
     QQw(:ncol,:)           = 0._r8
     QQi(:ncol,:)           = 0._r8
     QQnl(:ncol,:)          = 0._r8
     QQni(:ncol,:)          = 0._r8
     ACnl(:ncol,:)          = 0._r8
     ACni(:ncol,:)          = 0._r8

     QQw_prev(:ncol,:)      = 0._r8
     QQi_prev(:ncol,:)      = 0._r8
     QQnl_prev(:ncol,:)     = 0._r8
     QQni_prev(:ncol,:)     = 0._r8

     QQw_prog(:ncol,:)      = 0._r8
     QQi_prog(:ncol,:)      = 0._r8
     QQnl_prog(:ncol,:)     = 0._r8
     QQni_prog(:ncol,:)     = 0._r8

     QQ_final(:ncol,:)      = 0._r8                        
     QQw_final(:ncol,:)     = 0._r8                  
     QQi_final(:ncol,:)     = 0._r8           
     QQn_final(:ncol,:)     = 0._r8    
     QQnl_final(:ncol,:)    = 0._r8
     QQni_final(:ncol,:)    = 0._r8

     QQ_all(:ncol,:)        = 0._r8
     QQw_all(:ncol,:)       = 0._r8
     QQi_all(:ncol,:)       = 0._r8
     QQn_all(:ncol,:)       = 0._r8
     QQnl_all(:ncol,:)      = 0._r8
     QQni_all(:ncol,:)      = 0._r8

 ! Coefficient for computing QQ and related processes

     U(:ncol,:)             = 0._r8
     U_nc(:ncol,:)          = 0._r8
     G_nc(:ncol,:)          = 0._r8
     F_nc(:ncol,:)          = 0._r8

 ! Other

     qmin1(:ncol,:)         = 0._r8
     qmin2(:ncol,:)         = 0._r8
     qmin3(:ncol,:)         = 0._r8

     esat_b(:ncol)          = 0._r8     
     qsat_b(:ncol)          = 0._r8    
     dqsdT_b(:ncol)         = 0._r8 

     esat_a(:ncol)          = 0._r8     
     qsat_a(:ncol)          = 0._r8    
     dqsdT_a(:ncol)         = 0._r8 
     Twb_aw(:ncol,:)        = 0._r8
     qvwb_aw(:ncol,:)       = 0._r8

   ! ---------------- !
   ! Main computation ! 
   ! ---------------- !

   ! ---------------------------------- !
   ! Compute cumulus-related properties ! 
   ! ---------------------------------- !

   dacudt(:ncol,top_lev:pver) = &
        (a_cu0(:ncol,top_lev:pver) - a_cud(:ncol,top_lev:pver))/dt

   ! ---------------------------------------------------------------------- !
   ! set to zero for levels above
   ! ---------------------------------------------------------------------- !
   ql0(:ncol,:top_lev-1) = 0._r8
   qi0(:ncol,:top_lev-1) = 0._r8
   nl0(:ncol,:top_lev-1) = 0._r8
   ni0(:ncol,:top_lev-1) = 0._r8
   
   ! ---------------------------------------------------------------------- !
   ! Check if input non-cumulus pixels satisfie a non-negative constraint.  !
   ! If not, force all water vapor substances to be positive in all layers. !
   ! We should use 'old' cumulus properties for this routine.               !                
   ! ---------------------------------------------------------------------- !

   T1(:ncol,:)    =  T0(:ncol,:) 
   qv1(:ncol,:)   = qv0(:ncol,:) 
   ql1(:ncol,:)   = ql0(:ncol,:) 
   qi1(:ncol,:)   = qi0(:ncol,:) 
   nl1(:ncol,:)   = nl0(:ncol,:) 
   ni1(:ncol,:)   = ni0(:ncol,:) 

   
   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )


   qmin1(:ncol,:) = qmin(1)
   qmin2(:ncol,:) = qmin(ixcldliq)
   qmin3(:ncol,:) = qmin(ixcldice)

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv1, ql1, qi1, T1, qvten_pwi1, qlten_pwi1, &
                           qiten_pwi1, Tten_pwi1, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql1(i,k) .lt. qsmall ) then
          nlten_pwi1(i,k) = -nl1(i,k)/dt
          nl1(i,k)        = 0._r8
      endif 
      if( qi1(i,k) .lt. qsmall ) then
          niten_pwi1(i,k) = -ni1(i,k)/dt
          ni1(i,k)        = 0._r8
      endif 
   enddo
   enddo

   ! ------------------------------------------------------------- !
   ! Impose 'in-stratus condensate amount constraint'              !
   ! such that it is bounded by two limiting values.               !      
   ! This should also use 'old' cumulus properties since it is     !
   ! before applying external forcings.                            ! 
   ! Below 'QQw1,QQi1' are effective adjustive condensation        ! 
   ! Although this process also involves freezing of cloud         !
   ! liquid into ice, they can be and only can be expressed        !
   ! in terms of effective condensation.                           !
   ! ------------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( lchnk, ncol, k,                                   &
                                 p(:,k), T1(:,k), qv1(:,k), ql1(:,k), qi1(:,k),    &
                                 a_cud(:,k), zeros(:,k), zeros(:,k),               &
                                 zeros(:,k), zeros(:,k), zeros(:,k),               &
                                 landfrac, snowh,                                  &
                                 T_0(:,k), qv_0(:,k), ql_0(:,k), qi_0(:,k),        & 
                                 al_st_0(:,k), ai_st_0(:,k), ql_st_0(:,k), qi_st_0(:,k) )
      a_st_0(:ncol,k) = max(al_st_0(:ncol,k),ai_st_0(:ncol,k))
      QQw1(:ncol,k)   = (ql_0(:ncol,k) - ql1(:ncol,k))/dt
      QQi1(:ncol,k)   = (qi_0(:ncol,k) - qi1(:ncol,k))/dt
      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! Set a limit such that negative state not happens.  ! 
      ! -------------------------------------------------- !
      do i = 1, ncol
         if( QQw1(i,k) .le. 0._r8 ) then
             if( ql1(i,k) .gt. qsmall ) then
                 QQnl1(i,k) = QQw1(i,k)*nl1(i,k)/ql1(i,k)
                 QQnl1(i,k) = min(0._r8,cone*max(QQnl1(i,k),-nl1(i,k)/dt))
             else
                 QQnl1(i,k) = 0._r8
             endif  
         endif 
         if( QQi1(i,k) .le. 0._r8 ) then
             if( qi1(i,k) .gt. qsmall ) then
                 QQni1(i,k) = QQi1(i,k)*ni1(i,k)/qi1(i,k)
                 QQni1(i,k) = min(0._r8,cone*max(QQni1(i,k),-ni1(i,k)/dt))
             else
                 QQni1(i,k) = 0._r8
             endif  
         endif 
      enddo
   enddo
   nl_0(:ncol,top_lev:) = max(0._r8,nl1(:ncol,top_lev:)+QQnl1(:ncol,top_lev:)*dt) 
   ni_0(:ncol,top_lev:) = max(0._r8,ni1(:ncol,top_lev:)+QQni1(:ncol,top_lev:)*dt)

   ! ----------------------------------------------------------------------------- !
   ! Check if non-cumulus pixels of '_05' state satisfies non-negative constraint. !
   ! If not, force all water substances of '_05' state to be positive by imposing  !
   ! adjustive advection. We should use 'new' cumulus properties for this routine. !                
   ! ----------------------------------------------------------------------------- !

   T_05(:ncol,top_lev:)  =  T_0(:ncol,top_lev:) + (  A_T(:ncol,top_lev:) +  C_T(:ncol,top_lev:) ) * dt
   qv_05(:ncol,top_lev:) = qv_0(:ncol,top_lev:) + ( A_qv(:ncol,top_lev:) + C_qv(:ncol,top_lev:) ) * dt
   ql_05(:ncol,top_lev:) = ql_0(:ncol,top_lev:) + ( A_ql(:ncol,top_lev:) + C_ql(:ncol,top_lev:) ) * dt
   qi_05(:ncol,top_lev:) = qi_0(:ncol,top_lev:) + ( A_qi(:ncol,top_lev:) + C_qi(:ncol,top_lev:) ) * dt 
   nl_05(:ncol,top_lev:) = max(0._r8, nl_0(:ncol,top_lev:) + ( A_nl(:ncol,top_lev:) + C_nl(:ncol,top_lev:) ) * dt )
   ni_05(:ncol,top_lev:) = max(0._r8, ni_0(:ncol,top_lev:) + ( A_ni(:ncol,top_lev:) + C_ni(:ncol,top_lev:) ) * dt )

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv_05, ql_05, qi_05, T_05, A_qv_adj, &
                           A_ql_adj, A_qi_adj, A_T_adj, do_cldice)

   ! -------------------------------------------------------------- !
   ! Define reference state at the first iteration. This will be    !
   ! continuously updated within the iteration loop below.          !
   ! While equlibrium state properties are already output from the  !
   ! 'instratus_condensate', they will be re-computed within the    !
   ! each iteration process. At the first iteration, they will      !
   ! produce exactly identical results. Note that except at the     !
   ! very first iteration iter = 1, we must use updated cumulus     !
   ! properties at all the other iteration processes. Even at the   !
   ! first iteration, we should use updated cumulus properties      !
   ! when computing limiters for (Q,P,E).                           !
   ! -------------------------------------------------------------- !

   ! -------------------------------------------------------------- !
   ! Define variables at the reference state of the first iteration !
   ! -------------------------------------------------------------- !

   T(:ncol,top_lev:)     = T_0(:ncol,top_lev:)
   qv(:ncol,top_lev:)    = qv_0(:ncol,top_lev:)
   ql(:ncol,top_lev:)    = ql_0(:ncol,top_lev:)
   qi(:ncol,top_lev:)    = qi_0(:ncol,top_lev:)
   al_st(:ncol,top_lev:) = al_st_0(:ncol,top_lev:)
   ai_st(:ncol,top_lev:) = ai_st_0(:ncol,top_lev:)
   a_st(:ncol,top_lev:)  = a_st_0(:ncol,top_lev:)
   ql_st(:ncol,top_lev:) = ql_st_0(:ncol,top_lev:)
   qi_st(:ncol,top_lev:) = qi_st_0(:ncol,top_lev:)
   nl(:ncol,top_lev:)    = nl_0(:ncol,top_lev:)
   ni(:ncol,top_lev:)    = ni_0(:ncol,top_lev:)

   ! -------------------------- !
   ! Main iterative computation !
   ! -------------------------- !

   do iter = 1, niter

      ! ------------------------------------------ !
      ! Initialize array within the iteration loop !
      ! ------------------------------------------ !

      QQ(:,:)         = 0._r8
      QQw(:,:)        = 0._r8
      QQi(:,:)        = 0._r8
      QQnl(:,:)       = 0._r8
      QQni(:,:)       = 0._r8 
      QQw2(:,:)       = 0._r8
      QQi2(:,:)       = 0._r8
      QQnl2(:,:)      = 0._r8
      QQni2(:,:)      = 0._r8
      nlten_pwi2(:,:) = 0._r8
      niten_pwi2(:,:) = 0._r8
      ACnl(:,:)       = 0._r8
      ACni(:,:)       = 0._r8 
      aa(:,:)         = 0._r8
      bb(:,:)         = 0._r8

      do k = top_lev, pver

      ! "False" means that ice will not be taken into account.
      call findsp_vc(qv_05(:ncol,k),T_05(:ncol,k),p(:ncol,k), .false., &
           Twb_aw(:ncol,k),qvwb_aw(:ncol,k))

      call qsat_water(T_05(1:ncol,k), p(1:ncol,k), &
           esat_a(1:ncol), qsat_a(1:ncol), dqsdt=dqsdT_a(1:ncol))
      call qsat_water(T(1:ncol,k), p(1:ncol,k), &
           esat_b(1:ncol), qsat_b(1:ncol), dqsdt=dqsdT_b(1:ncol))

      if( iter .eq. 1 ) then
          a_cu(:ncol,k) = a_cud(:ncol,k)
      else
          a_cu(:ncol,k) = a_cu0(:ncol,k)
      endif
      do i = 1, ncol
         U(i,k)    =  qv(i,k)/qsat_b(i)
         U_nc(i,k) =  U(i,k)
      enddo
      if( CAMstfrac ) then
          call astG_RHU(U_nc(:,k),p(:,k),qv(:,k),landfrac(:),snowh(:),al_st_nc(:,k),G_nc(:,k),ncol)
      else
          call astG_PDF(U_nc(:,k),p(:,k),qv(:,k),landfrac(:),snowh(:),al_st_nc(:,k),G_nc(:,k),ncol)
      endif
      call aist_vector(qv(:,k),T(:,k),p(:,k),qi(:,k),landfrac(:),snowh(:),ai_st_nc(:,k),ncol)
      ai_st(:ncol,k)  =  (1._r8-a_cu(:ncol,k))*ai_st_nc(:ncol,k)
      al_st(:ncol,k)  =  (1._r8-a_cu(:ncol,k))*al_st_nc(:ncol,k)
      a_st(:ncol,k)   =  max(al_st(:ncol,k),ai_st(:ncol,k))  

      do i = 1, ncol

         ! -------------------------------------------------------- !
         ! Compute basic thermodynamic coefficients for computing Q !
         ! -------------------------------------------------------- !

         alpha  =  1._r8/qsat_b(i)
         beta   =  dqsdt_b(i)*(qv(i,k)/qsat_b(i)**2)
         betast =  alpha*dqsdT_b(i) 
         gammal =  alpha + (latvap/cpair)*beta
         gammai =  alpha + ((latvap+latice)/cpair)*beta
         gammaQ =  alpha + (latvap/cpair)*beta
         deltal =  1._r8 + a_st(i,k)*(latvap/cpair)*(betast/alpha)
         deltai =  1._r8 + a_st(i,k)*((latvap+latice)/cpair)*(betast/alpha)
         A_Tc   =  A_T(i,k)+A_T_adj(i,k)-(latvap/cpair)*(A_ql(i,k)+A_ql_adj(i,k))-((latvap+latice)/cpair)*(A_qi(i,k)+A_qi_adj(i,k))
         A_qt   =  A_qv(i,k) + A_qv_adj(i,k) + A_ql(i,k) + A_ql_adj(i,k) + A_qi(i,k) + A_qi_adj(i,k)
         C_Tc   =  C_T(i,k) - (latvap/cpair)*C_ql(i,k) - ((latvap+latice)/cpair)*C_qi(i,k)
         C_qt   =  C_qv(i,k) + C_ql(i,k) + C_qi(i,k)
         dTcdt  =  A_Tc + C_Tc
         dqtdt  =  A_qt + C_qt
       ! dqtstldt = A_qt + C_ql(i,k)/max(1.e-2_r8,al_st(i,k))                             ! Original  
       ! dqtstldt = A_qt - A_qi(i,k) - A_qi_adj(i,k) + C_ql(i,k)/max(1.e-2_r8,al_st(i,k)) ! New 1 on Dec.30.2009.
         dqtstldt = A_qt - A_qi(i,k) - A_qi_adj(i,k) + C_qlst(i,k)                        ! New 2 on Dec.30.2009.
       ! dqtstldt = A_qt + C_qt                                                           ! Original Conservative treatment
       ! dqtstldt = A_qt - A_qi(i,k) - A_qi_adj(i,k) + C_qt - C_qi(i,k)            ! New Conservative treatment on Dec.30.2009
         dqidt = A_qi(i,k) + A_qi_adj(i,k) + C_qi(i,k) 

         anic    = max(1.e-8_r8,(1._r8-a_cu(i,k)))
         GG      = G_nc(i,k)/anic
         aa(1,1) = gammal*al_st(i,k)
         aa(1,2) = GG + gammal*cc*ql_st(i,k)          
         aa(2,1) = alpha + (latvap/cpair)*betast*al_st(i,k)
         aa(2,2) = (latvap/cpair)*betast*cc*ql_st(i,k) 
         bb(1,1) = alpha*dqtdt - beta*dTcdt - gammai*dqidt - GG*al_st_nc(i,k)*dacudt(i,k) + F_nc(i,k) 
         bb(2,1) = alpha*dqtstldt - betast*(dTcdt + ((latvap+latice)/cpair)*dqidt) 
         call gaussj(aa(1:2,1:2),2,2,bb(1:2,1),1,1)
         dqlstdt = bb(1,1)
         dalstdt = bb(2,1)
         QQ(i,k) = al_st(i,k)*dqlstdt + cc*ql_st(i,k)*dalstdt - ( A_ql(i,k) + A_ql_adj(i,k) + C_ql(i,k) )

       ! ------------------------------------------------------------ !
       ! Limiter for QQ                                               !
       ! Here, 'fice' should be from the reference equilibrium state  !
       ! since QQ itself is computed from the reference state.        !
       ! From the assumption used for derivation of QQ(i), it must be !
       ! that QQw(i) = QQ(i)*(1._r8-fice(i)), QQi(i) = QQ(i)*fice(i)  !  
       ! ------------------------------------------------------------ !

         if( QQ(i,k) .ge. 0._r8 ) then
             QQmax    = (qv_05(i,k) - qmin(1))/dt ! For ghost cumulus & semi-ghost ice stratus
             QQmax    = max(0._r8,QQmax) 
             QQ(i,k)  = min(QQ(i,k),QQmax)
             QQw(i,k) = QQ(i,k)
             QQi(i,k) = 0._r8 
         else
             QQmin  = 0._r8
             if( qv_05(i,k) .lt. qsat_a(i) ) QQmin = min(0._r8,cone*(qv_05(i,k)-qvwb_aw(i,k))/dt)
             QQ(i,k)  = max(QQ(i,k),QQmin)
             QQw(i,k) = QQ(i,k)
             QQi(i,k) = 0._r8
             QQwmin   = min(0._r8,-cone*ql_05(i,k)/dt)
             QQimin   = min(0._r8,-cone*qi_05(i,k)/dt)
             QQw(i,k) = min(0._r8,max(QQw(i,k),QQwmin))
             QQi(i,k) = min(0._r8,max(QQi(i,k),QQimin))
         endif

       ! -------------------------------------------------- !
       ! Reduce droplet concentration if evaporation occurs !
       ! Note 'QQnl1,QQni1' are computed from the reference !
       ! equilibrium state but limiter is from 'nl_05'.     !
       ! -------------------------------------------------- !

         if( QQw(i,k) .lt. 0._r8 ) then
             if( ql_05(i,k) .gt. qsmall ) then
                 QQnl(i,k) = QQw(i,k)*nl_05(i,k)/ql_05(i,k)
                 QQnl(i,k) = min(0._r8,cone*max(QQnl(i,k),-nl_05(i,k)/dt))
             else
                 QQnl(i,k) = 0._r8
             endif  
         endif 

         if( QQi(i,k) .lt. 0._r8 ) then
             if( qi_05(i,k) .gt. qsmall ) then
                 QQni(i,k) = QQi(i,k)*ni_05(i,k)/qi_05(i,k)
                 QQni(i,k) = min(0._r8,cone*max(QQni(i,k),-ni_05(i,k)/dt))
             else
                 QQni(i,k) = 0._r8
             endif  
         endif 

      enddo
      enddo

    ! -------------------------------------------------------------------- !
    ! Until now, we have finished computing all necessary tendencies       ! 
    ! from the equilibrium input state (T_0).                              !
    ! If ramda = 0 : fully explicit scheme                                 !
    !    ramda = 1 : fully implicit scheme                                 !
    ! Note that 'ramda = 0.5 with niter = 2' can mimic                     !
    ! -------------------------------------------------------------------- !

      if( iter .eq. 1 ) then
          QQw_prev(:ncol,top_lev:)  = QQw(:ncol,top_lev:)       
          QQi_prev(:ncol,top_lev:)  = QQi(:ncol,top_lev:)   
          QQnl_prev(:ncol,top_lev:) = QQnl(:ncol,top_lev:)       
          QQni_prev(:ncol,top_lev:) = QQni(:ncol,top_lev:)   
      endif

      QQw_prog(:ncol,top_lev:)   = ramda*QQw(:ncol,top_lev:)   + (1._r8-ramda)*QQw_prev(:ncol,top_lev:)
      QQi_prog(:ncol,top_lev:)   = ramda*QQi(:ncol,top_lev:)   + (1._r8-ramda)*QQi_prev(:ncol,top_lev:)
      QQnl_prog(:ncol,top_lev:)  = ramda*QQnl(:ncol,top_lev:)  + (1._r8-ramda)*QQnl_prev(:ncol,top_lev:)
      QQni_prog(:ncol,top_lev:)  = ramda*QQni(:ncol,top_lev:)  + (1._r8-ramda)*QQni_prev(:ncol,top_lev:)

      QQw_prev(:ncol,top_lev:)   = QQw_prog(:ncol,top_lev:)
      QQi_prev(:ncol,top_lev:)   = QQi_prog(:ncol,top_lev:)
      QQnl_prev(:ncol,top_lev:)  = QQnl_prog(:ncol,top_lev:)
      QQni_prev(:ncol,top_lev:)  = QQni_prog(:ncol,top_lev:)

    ! -------------------------------------------------------- !
    ! Compute final prognostic state on which final diagnostic !
    ! in-stratus condensate adjustment is applied in the below.!
    ! Important : I must check whether there are any external  !  
    !             advective forcings of 'A_nl(i,k),A_ni(i,k)'. !
    !             Even they are (i.e., advection of aerosol),  !
    !             actual droplet activation will be performd   !
    !             in microphysics, so it will be completely    !
    !             reasonable to 'A_nl(i,k)=A_ni(i,k)=0'.       !
    ! -------------------------------------------------------- !

    do k = top_lev, pver
    do i = 1, ncol
       T_prime0(i,k)  = T_0(i,k)  + dt*( A_T(i,k)  +  A_T_adj(i,k) +  C_T(i,k) + &
            (latvap*QQw_prog(i,k)+(latvap+latice)*QQi_prog(i,k))/cpair )
       qv_prime0(i,k) = qv_0(i,k) + dt*( A_qv(i,k) + A_qv_adj(i,k) + C_qv(i,k) - QQw_prog(i,k) - QQi_prog(i,k) )
       ql_prime0(i,k) = ql_0(i,k) + dt*( A_ql(i,k) + A_ql_adj(i,k) + C_ql(i,k) + QQw_prog(i,k) )
       qi_prime0(i,k) = qi_0(i,k) + dt*( A_qi(i,k) + A_qi_adj(i,k) + C_qi(i,k) + QQi_prog(i,k) )
       nl_prime0(i,k) = max(0._r8,nl_0(i,k) + dt*( A_nl(i,k) + C_nl(i,k) + QQnl_prog(i,k) ))
       ni_prime0(i,k) = max(0._r8,ni_0(i,k) + dt*( A_ni(i,k) + C_ni(i,k) + QQni_prog(i,k) ))
       if( ql_prime0(i,k) .lt. qsmall ) nl_prime0(i,k) = 0._r8
       if( qi_prime0(i,k) .lt. qsmall ) ni_prime0(i,k) = 0._r8
    enddo
    enddo

   ! -------------------------------------------------- !
   ! Perform diagnostic 'positive_moisture' constraint. !
   ! -------------------------------------------------- !

   T_dprime(:ncol,top_lev:)  =  T_prime0(:ncol,top_lev:) 
   qv_dprime(:ncol,top_lev:) = qv_prime0(:ncol,top_lev:) 
   ql_dprime(:ncol,top_lev:) = ql_prime0(:ncol,top_lev:) 
   qi_dprime(:ncol,top_lev:) = qi_prime0(:ncol,top_lev:) 
   nl_dprime(:ncol,top_lev:) = nl_prime0(:ncol,top_lev:) 
   ni_dprime(:ncol,top_lev:) = ni_prime0(:ncol,top_lev:) 

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp,          & 
                           qv_dprime, ql_dprime, qi_dprime, T_dprime,  &
                           qvten_pwi2, qlten_pwi2, qiten_pwi2, Tten_pwi2, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql_dprime(i,k) .lt. qsmall ) then
          nlten_pwi2(i,k) = -nl_dprime(i,k)/dt
          nl_dprime(i,k)   = 0._r8
      endif 
      if( qi_dprime(i,k) .lt. qsmall ) then
          niten_pwi2(i,k) = -ni_dprime(i,k)/dt
          ni_dprime(i,k)   = 0._r8
      endif 
   enddo
   enddo

   ! -------------------------------------------------------------- !
   ! Add tendency associated with detrainment of cumulus condensate !
   ! This tendency is not used in computing Q                       !
   ! Since D_ql,D_qi,D_nl,D_ni > 0, don't need to worry about       !
   ! negative scalar.                                               !
   ! This tendency is not reflected into Fzs2, which is OK.         !
   ! -------------------------------------------------------------- !

   T_dprime(:ncol,top_lev:)   =  T_dprime(:ncol,top_lev:)  + D_T(:ncol,top_lev:) * dt 
   qv_dprime(:ncol,top_lev:)  = qv_dprime(:ncol,top_lev:) + D_qv(:ncol,top_lev:) * dt 
   ql_dprime(:ncol,top_lev:)  = ql_dprime(:ncol,top_lev:) + D_ql(:ncol,top_lev:) * dt
   qi_dprime(:ncol,top_lev:)  = qi_dprime(:ncol,top_lev:) + D_qi(:ncol,top_lev:) * dt
   nl_dprime(:ncol,top_lev:)  = nl_dprime(:ncol,top_lev:) + D_nl(:ncol,top_lev:) * dt 
   ni_dprime(:ncol,top_lev:)  = ni_dprime(:ncol,top_lev:) + D_ni(:ncol,top_lev:) * dt

   ! ---------------------------------------------------------- !
   ! Impose diagnostic upper and lower limits on the in-stratus !
   ! condensate amount. This produces a final equilibrium state !
   ! at the end of each iterative process.                      !
   ! ---------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( lchnk          , ncol           , k              , p(:,k)        , &
                                 T_dprime(:,k)  , qv_dprime(:,k) , ql_dprime(:,k) , qi_dprime(:,k), &
                                 a_cu0(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 & 
                                 zeros(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 &
                                 landfrac       , snowh          ,                                  &
                                 T_star(:,k)    , qv_star(:,k)   , ql_star(:,k)   , qi_star(:,k)  , & 
                                 al_st_star(:,k), ai_st_star(:,k), ql_st_star(:,k), qi_st_star(:,k) )
      a_st_star(:ncol,k)  = max(al_st_star(:ncol,k),ai_st_star(:ncol,k))
      QQw2(:ncol,k) = (ql_star(:ncol,k) - ql_dprime(:ncol,k))/dt
      QQi2(:ncol,k) = (qi_star(:ncol,k) - qi_dprime(:ncol,k))/dt
      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! -------------------------------------------------- !
      do i = 1, ncol
         if( QQw2(i,k) .le. 0._r8 ) then
             if( ql_dprime(i,k) .ge. qsmall ) then
                 QQnl2(i,k) = QQw2(i,k)*nl_dprime(i,k)/ql_dprime(i,k)
                 QQnl2(i,k) = min(0._r8,cone*max(QQnl2(i,k),-nl_dprime(i,k)/dt))
             else
                 QQnl2(i,k) = 0._r8
             endif  
         endif 
         if( QQi2(i,k) .le. 0._r8 ) then
             if( qi_dprime(i,k) .gt. qsmall ) then
                 QQni2(i,k) = QQi2(i,k)*ni_dprime(i,k)/qi_dprime(i,k)
                 QQni2(i,k) = min(0._r8,cone*max(QQni2(i,k),-ni_dprime(i,k)/dt))
             else
                 QQni2(i,k) = 0._r8
             endif  
         endif 
      enddo
   enddo
   nl_star(:ncol,top_lev:) = max(0._r8,nl_dprime(:ncol,top_lev:)+QQnl2(:ncol,top_lev:)*dt) 
   ni_star(:ncol,top_lev:) = max(0._r8,ni_dprime(:ncol,top_lev:)+QQni2(:ncol,top_lev:)*dt)

   ! ------------------------------------------ !
   ! Final adjustment of droplet concentration. !
   ! Set # to zero if there is no cloud.        !
   ! ------------------------------------------ !

   do k = top_lev, pver
   do i = 1, ncol 
      if( ql_star(i,k) .lt. qsmall ) then
          ACnl(i,k) = - nl_star(i,k)/dt
          nl_star(i,k) = 0._r8
      endif
      if( qi_star(i,k) .lt. qsmall ) then
          ACni(i,k) = - ni_star(i,k)/dt
          ni_star(i,k) = 0._r8
      endif
   enddo
   enddo

   ! ----------------------------------------------------- !
   ! Define equilibrium reference state for next iteration !
   ! ----------------------------------------------------- !

   T(:ncol,top_lev:)     = T_star(:ncol,top_lev:)
   qv(:ncol,top_lev:)    = qv_star(:ncol,top_lev:)
   ql(:ncol,top_lev:)    = ql_star(:ncol,top_lev:)
   qi(:ncol,top_lev:)    = qi_star(:ncol,top_lev:)
   al_st(:ncol,top_lev:) = al_st_star(:ncol,top_lev:)
   ai_st(:ncol,top_lev:) = ai_st_star(:ncol,top_lev:)
   a_st(:ncol,top_lev:)  = a_st_star(:ncol,top_lev:)
   ql_st(:ncol,top_lev:) = ql_st_star(:ncol,top_lev:)
   qi_st(:ncol,top_lev:) = qi_st_star(:ncol,top_lev:)
   nl(:ncol,top_lev:)    = nl_star(:ncol,top_lev:)
   ni(:ncol,top_lev:)    = ni_star(:ncol,top_lev:)

   enddo ! End of 'iter' prognostic iterative computation

   ! ------------------------------------------------------------------------ !
   ! Compute final tendencies of main output variables and diagnostic outputs !
   ! Note that the very input state [T0,qv0,ql0,qi0] are                      !
   ! marched to [T_star,qv_star,ql_star,qi_star] with equilibrium             !
   ! stratus informations of [a_st_star,ql_st_star,qi_st_star] by             !
   ! below final tendencies and [A_T,A_qv,A_ql,A_qi]                          !
   ! ------------------------------------------------------------------------ !

   ! ------------------ !
   ! Process tendencies !
   ! ------------------ !

   QQw_final(:ncol,top_lev:)  = QQw_prog(:ncol,top_lev:)
   QQi_final(:ncol,top_lev:)  = QQi_prog(:ncol,top_lev:)
   QQ_final(:ncol,top_lev:)   = QQw_final(:ncol,top_lev:) + QQi_final(:ncol,top_lev:)
   QQw_all(:ncol,top_lev:)    = QQw_prog(:ncol,top_lev:)  + QQw1(:ncol,top_lev:) + QQw2(:ncol,top_lev:) + &
        qlten_pwi1(:ncol,top_lev:) + qlten_pwi2(:ncol,top_lev:) + A_ql_adj(:ncol,top_lev:)
   QQi_all(:ncol,top_lev:)    = QQi_prog(:ncol,top_lev:)  + QQi1(:ncol,top_lev:) + QQi2(:ncol,top_lev:) + &
        qiten_pwi1(:ncol,top_lev:) + qiten_pwi2(:ncol,top_lev:) + A_qi_adj(:ncol,top_lev:)
   QQ_all(:ncol,top_lev:)     = QQw_all(:ncol,top_lev:)   + QQi_all(:ncol,top_lev:)
   QQnl_final(:ncol,top_lev:) = QQnl_prog(:ncol,top_lev:)
   QQni_final(:ncol,top_lev:) = QQni_prog(:ncol,top_lev:)
   QQn_final(:ncol,top_lev:)  = QQnl_final(:ncol,top_lev:) + QQni_final(:ncol,top_lev:)
   QQnl_all(:ncol,top_lev:)   = QQnl_prog(:ncol,top_lev:)  + QQnl1(:ncol,top_lev:) + QQnl2(:ncol,top_lev:) + &
        nlten_pwi1(:ncol,top_lev:) + nlten_pwi2(:ncol,top_lev:) + ACnl(:ncol,top_lev:) + A_nl_adj(:ncol,top_lev:)
   QQni_all(:ncol,top_lev:)   = QQni_prog(:ncol,top_lev:)  + QQni1(:ncol,top_lev:) + QQni2(:ncol,top_lev:) + &
        niten_pwi1(:ncol,top_lev:) + niten_pwi2(:ncol,top_lev:) + ACni(:ncol,top_lev:) + A_ni_adj(:ncol,top_lev:)
   QQn_all(:ncol,top_lev:)    = QQnl_all(:ncol,top_lev:)   + QQni_all(:ncol,top_lev:)
   qme(:ncol,top_lev:)        = QQ_final(:ncol,top_lev:)   
   qvadj(:ncol,top_lev:)      = qvten_pwi1(:ncol,top_lev:) + qvten_pwi2(:ncol,top_lev:) + A_qv_adj(:ncol,top_lev:)
   qladj(:ncol,top_lev:)      = qlten_pwi1(:ncol,top_lev:) + qlten_pwi2(:ncol,top_lev:) + A_ql_adj(:ncol,top_lev:)
   qiadj(:ncol,top_lev:)      = qiten_pwi1(:ncol,top_lev:) + qiten_pwi2(:ncol,top_lev:) + A_qi_adj(:ncol,top_lev:)
   qllim(:ncol,top_lev:)      = QQw1      (:ncol,top_lev:) + QQw2      (:ncol,top_lev:)
   qilim(:ncol,top_lev:)      = QQi1      (:ncol,top_lev:) + QQi2      (:ncol,top_lev:)

   ! ----------------- !
   ! Output tendencies !
   ! ----------------- !

   s_tendout(:ncol,top_lev:)  = cpair*( T_star(:ncol,top_lev:)  -  T0(:ncol,top_lev:) )/dt - &
        cpair*(A_T(:ncol,top_lev:)+C_T(:ncol,top_lev:))
   qv_tendout(:ncol,top_lev:) =    ( qv_star(:ncol,top_lev:) - qv0(:ncol,top_lev:) )/dt - &
        (A_qv(:ncol,top_lev:)+C_qv(:ncol,top_lev:))
   ql_tendout(:ncol,top_lev:) =    ( ql_star(:ncol,top_lev:) - ql0(:ncol,top_lev:) )/dt - &
        (A_ql(:ncol,top_lev:)+C_ql(:ncol,top_lev:))
   qi_tendout(:ncol,top_lev:) =    ( qi_star(:ncol,top_lev:) - qi0(:ncol,top_lev:) )/dt - &
        (A_qi(:ncol,top_lev:)+C_qi(:ncol,top_lev:))
   nl_tendout(:ncol,top_lev:) =    ( nl_star(:ncol,top_lev:) - nl0(:ncol,top_lev:) )/dt - &
        (A_nl(:ncol,top_lev:)+C_nl(:ncol,top_lev:))
   ni_tendout(:ncol,top_lev:) =    ( ni_star(:ncol,top_lev:) - ni0(:ncol,top_lev:) )/dt - &
        (A_ni(:ncol,top_lev:)+C_ni(:ncol,top_lev:))

   if (.not. do_cldice) then
      do k = top_lev, pver
         do i = 1, ncol

            ! Don't want either qi or ni tendencies, but the code above is somewhat convoluted and
            ! is trying to adjust both (small numbers). Just force it to zero here.
            qi_tendout(i,k) = 0._r8
            ni_tendout(i,k) = 0._r8
          end do
      end do
   end if

   ! ------------------ !
   ! Net cloud fraction !
   ! ------------------ !

   cld(:ncol,top_lev:) = a_st_star(:ncol,top_lev:) + a_cu0(:ncol,top_lev:)

   ! --------------------------------- !
   ! Updated grid-mean state variables !
   ! --------------------------------- !

   T0(:ncol,top_lev:)  = T_star(:ncol,top_lev:)
   qv0(:ncol,top_lev:) = qv_star(:ncol,top_lev:)
   ql0(:ncol,top_lev:) = ql_star(:ncol,top_lev:)
   qi0(:ncol,top_lev:) = qi_star(:ncol,top_lev:)
   nl0(:ncol,top_lev:) = nl_star(:ncol,top_lev:)
   ni0(:ncol,top_lev:) = ni_star(:ncol,top_lev:)

   return
   end subroutine mmacro_pcond

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine instratus_condensate( lchnk, ncol, k,                      &  
                                    p_in, T0_in, qv0_in, ql0_in, qi0_in, & 
                                    a_dc_in, ql_dc_in, qi_dc_in,         &
                                    a_sc_in, ql_sc_in, qi_sc_in,         & 
                                    landfrac, snowh,                     &
                                    T_out, qv_out, ql_out, qi_out,       &
                                    al_st_out, ai_st_out, ql_st_out, qi_st_out )

   ! ------------------------------------------------------- !
   ! Diagnostically force in-stratus condensate to be        ! 
   ! in the range of 'qlst_min < qc_st < qlst_max'           !
   ! whenever stratus exists in the equilibrium state        !
   ! ------------------------------------------------------- !

   use time_manager,  only: is_first_step, get_nstep

   implicit none

   integer,  intent(in)  :: lchnk                ! Chunk identifier
   integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
   integer,  intent(in)  :: k                    ! Layer index

   real(r8), intent(in)  :: p_in(pcols)          ! Pressure [Pa]
   real(r8), intent(in)  :: T0_in(pcols)         ! Temperature [K]
   real(r8), intent(in)  :: qv0_in(pcols)        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0_in(pcols)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0_in(pcols)        ! Grid-mean IWC [kg/kg]

   real(r8), intent(in)  :: a_dc_in(pcols)       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc_in(pcols)      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc_in(pcols)      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc_in(pcols)       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc_in(pcols)      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc_in(pcols)      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: landfrac(pcols)      ! Land fraction
   real(r8), intent(in)  :: snowh(pcols)         ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: T_out(pcols)         ! Temperature [K]
   real(r8), intent(out) :: qv_out(pcols)        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(out) :: ql_out(pcols)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(out) :: qi_out(pcols)        ! Grid-mean IWC [kg/kg]

   real(r8), intent(out) :: al_st_out(pcols)     ! Liquid stratus fraction
   real(r8), intent(out) :: ai_st_out(pcols)     ! Ice stratus fraction
   real(r8), intent(out) :: ql_st_out(pcols)     ! In-stratus LWC [kg/kg]
   real(r8), intent(out) :: qi_st_out(pcols)     ! In-stratus IWC [kg/kg]

   ! Local variables

   integer i                                     ! Column    index

   real(r8) p    
   real(r8) T0   
   real(r8) qv0    
   real(r8) ql0    
   real(r8) qi0    
   real(r8) a_dc   
   real(r8) ql_dc  
   real(r8) qi_dc  
   real(r8) a_sc   
   real(r8) ql_sc  
   real(r8) qi_sc  
   real(r8) esat0  
   real(r8) qsat0  
   real(r8) U0     
   real(r8) U0_nc  
   real(r8) G0_nc
   real(r8) al0_st_nc            
   real(r8) al0_st
   real(r8) ai0_st_nc            
   real(r8) ai0_st               
   real(r8) a0_st               
   real(r8) ql0_nc
   real(r8) qi0_nc
   real(r8) qc0_nc
   real(r8) ql0_st
   real(r8) qi0_st
   real(r8) qc0_st
   real(r8) T   
   real(r8) qv    
   real(r8) ql    
   real(r8) qi
   real(r8) ql_st
   real(r8) qi_st
   real(r8) es  
   real(r8) qs  
   real(r8) esat_in(pcols)  
   real(r8) qsat_in(pcols)  
   real(r8) U0_in(pcols)  
   real(r8) al0_st_nc_in(pcols)
   real(r8) ai0_st_nc_in(pcols)
   real(r8) G0_nc_in(pcols)
   integer  idxmod 
   real(r8) U
   real(r8) U_nc
   real(r8) al_st_nc
   real(r8) ai_st_nc
   real(r8) G_nc
   real(r8) a_st
   real(r8) al_st
   real(r8) ai_st
   real(r8) Tmin0
   real(r8) Tmax0
   real(r8) Tmin
   real(r8) Tmax
   integer caseid

   ! ---------------- !
   ! Main Computation ! 
   ! ---------------- !

   call qsat_water(T0_in(1:ncol), p_in(1:ncol), &
        esat_in(1:ncol), qsat_in(1:ncol))
   U0_in(:ncol) = qv0_in(:ncol)/qsat_in(:ncol)
   if( CAMstfrac ) then
       call astG_RHU(U0_in(:),p_in(:),qv0_in(:),landfrac(:),snowh(:),al0_st_nc_in(:),G0_nc_in(:),ncol)
   else
       call astG_PDF(U0_in(:),p_in(:),qv0_in(:),landfrac(:),snowh(:),al0_st_nc_in(:),G0_nc_in(:),ncol)
   endif
   call aist_vector(qv0_in(:),T0_in(:),p_in(:),qi0_in(:),landfrac(:),snowh(:),ai0_st_nc_in(:),ncol)

   do i = 1, ncol

      ! ---------------------- !
      ! Define local variables !
      ! ---------------------- !

      p   = p_in(i)

      T0  = T0_in(i)
      qv0 = qv0_in(i)
      ql0 = ql0_in(i)
      qi0 = qi0_in(i)

      a_dc  = a_dc_in(i)
      ql_dc = ql_dc_in(i)
      qi_dc = qi_dc_in(i)

      a_sc  = a_sc_in(i)
      ql_sc = ql_sc_in(i)
      qi_sc = qi_sc_in(i)

      ql_dc = 0._r8
      qi_dc = 0._r8
      ql_sc = 0._r8
      qi_sc = 0._r8

      es  = esat_in(i) 
      qs  = qsat_in(i) 

      idxmod = 0
      caseid = -1

      ! ------------------------------------------------------------ !
      ! Force the grid-mean RH to be smaller than 1 if oversaturated !
      ! In order to be compatible with reduced 3x3 QQ, condensation  !
      ! should occur only into the liquid in gridmean_RH.            !
      ! ------------------------------------------------------------ !

      if( qv0 .gt. qs ) then
          call gridmean_RH( lchnk, i, k, p, T0, qv0, ql0, qi0,      &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                            landfrac(i), snowh(i) )
          call qsat_water(T0, p, esat0, qsat0)
          U0      = (qv0/qsat0)
          U0_nc   =  U0 
          if( CAMstfrac ) then
              call astG_RHU_single(U0_nc,p,qv0,landfrac(i),snowh(i),al0_st_nc,G0_nc)
          else
              call astG_PDF_single(U0_nc,p,qv0,landfrac(i),snowh(i),al0_st_nc,G0_nc)
          endif
          call aist_single(qv0,T0,p,qi0,landfrac(i),snowh(i),ai0_st_nc)
          ai0_st  = (1._r8-a_dc-a_sc)*ai0_st_nc
          al0_st  = (1._r8-a_dc-a_sc)*al0_st_nc
          a0_st   = max(ai0_st,al0_st)         
          idxmod  = 1 
      else
          ai0_st  = (1._r8-a_dc-a_sc)*ai0_st_nc_in(i)
          al0_st  = (1._r8-a_dc-a_sc)*al0_st_nc_in(i)
      endif    
      a0_st   = max(ai0_st,al0_st)         

      ! ----------------------- ! 
      ! Handling of input state !
      ! ----------------------- !

      ql0_nc  = max(0._r8,ql0-a_dc*ql_dc-a_sc*ql_sc)
      qi0_nc  = max(0._r8,qi0-a_dc*qi_dc-a_sc*qi_sc)
      qc0_nc  = ql0_nc + qi0_nc 

      Tmin0 = T0 - (latvap/cpair)*ql0
      Tmax0 = T0 + ((latvap+latice)/cpair)*qv0

      ! ------------------------------------------------------------- !
      ! Do nothing and just exit if generalized in-stratus condensate !
      ! condition is satisfied. This includes the case I.             !
      ! For 4x4 liquid stratus, a0_st --> al0_st.                     ! 
      ! ------------------------------------------------------------- !
      if( ( ql0_nc .ge. qlst_min*al0_st ) .and. ( ql0_nc .le. qlst_max*al0_st ) ) then

          ! ------------------ !
          ! This is the case I !
          ! ------------------ ! 
             T = T0
             qv = qv0
             ql = ql0
             qi = qi0
             caseid = 0
             goto 10
      else
         ! ----------------------------- !
         ! This is case II : Dense Cloud !
         ! ----------------------------- !   
         if( al0_st .eq. 0._r8 .and. ql0_nc .gt. 0._r8 ) then
             ! ------------------------------------- !
             ! Compute hypothetical full evaporation !
             ! ------------------------------------- !
             T  = Tmin0
             qv = qv0 + ql0 
             call qsat_water(T, p, es, qs)
             U  = qv/qs
             U_nc = U  
             if( CAMstfrac ) then
                 call astG_RHU_single(U_nc,p,qv,landfrac(i),snowh(i),al_st_nc,G_nc)
             else
                 call astG_PDF_single(U_nc,p,qv,landfrac(i),snowh(i),al_st_nc,G_nc)
             endif
             al_st = (1._r8-a_dc-a_sc)*al_st_nc  
             caseid = 0

             if( al_st .eq. 0._r8 ) then
                 ql = 0._r8
                 qi = qi0
                 idxmod = 1
                 caseid = 1
                 goto 10
             else
                 ! ------------------------------------------- !
                 ! Evaporate until qc_st decreases to qlst_max !
                 ! ------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = T0 
                 call instratus_core( lchnk, i, k, p,                              &
                                      T0, qv0, ql0, 0._r8,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_max, Tmin, Tmax, landfrac(i), snowh(i), &
                                      T, qv, ql, qi )   
                 idxmod = 1
                 caseid = 2
                 goto 10
             endif
         ! ------------------------------ !
         ! This is case III : Empty Cloud !
         ! ------------------------------ !  
         elseif( al0_st .gt. 0._r8 .and. ql0_nc .eq. 0._r8 ) then
              ! ------------------------------------------ ! 
              ! Condense until qc_st increases to qlst_min !
              ! ------------------------------------------ !
              Tmin = Tmin0
              Tmax = Tmax0  
              call instratus_core( lchnk, i, k, p,                              &
                                   T0, qv0, ql0, 0._r8,                         &
                                   a_dc, ql_dc, qi_dc,                          &
                                   a_sc, ql_sc, qi_sc, ai0_st,                  &
                                   qlst_min, Tmin, Tmax, landfrac(i), snowh(i), &
                                   T, qv, ql, qi )   
              idxmod = 1 
              caseid = 3
              goto 10
         ! --------------- !
         ! This is case IV !
         ! --------------- !   
         elseif( al0_st .gt. 0._r8 .and. ql0_nc .gt. 0._r8 ) then

             if( ql0_nc .gt. qlst_max*al0_st ) then
                 ! --------------------------------------- !
                 ! Evaporate until qc_st drops to qlst_max !
                 ! --------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0
                 call instratus_core( lchnk, i, k, p,                              &
                                      T0, qv0, ql0, 0._r8,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_max, Tmin, Tmax, landfrac(i), snowh(i), &
                                      T, qv, ql, qi )   
                 idxmod = 1
                 caseid = 4
                 goto 10
             elseif( ql0_nc .lt. qlst_min*al0_st ) then
                 ! -------------------------------------------- !
                 ! Condensate until qc_st increases to qlst_min !
                 ! -------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0 
                 call instratus_core( lchnk, i, k, p,                              &
                                      T0, qv0, ql0, 0._r8,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_min, Tmin, Tmax, landfrac(i), snowh(i), & 
                                      T, qv, ql, qi )   
                 idxmod = 1
                 caseid = 5
                 goto 10
             else
                 ! ------------------------------------------------ !
                 ! This case should not happen. Issue error message !
                 ! ------------------------------------------------ !
                 write(iulog,*) 'Impossible case1 in instratus_condensate' 
                 call endrun
             endif
         ! ------------------------------------------------ !                   
         ! This case should not happen. Issue error message !
         ! ------------------------------------------------ !    
         else
             write(iulog,*) 'Impossible case2 in instratus_condensate' 
             write(iulog,*)  al0_st, a_sc, a_dc
             write(iulog,*)  1000*ql0_nc, 1000*(ql0+qi0)
             call endrun
         endif
      endif

10 continue   

   ! -------------------------------------------------- !
   ! Force final energy-moisture conserving consistency !
   ! -------------------------------------------------- !

     qi = qi0

     if( idxmod .eq. 1 ) then
         call aist_single(qv,T,p,qi,landfrac(i),snowh(i),ai_st_nc)
         ai_st = (1._r8-a_dc-a_sc)*ai_st_nc
         call qsat_water(T, p, es, qs)
         U     = (qv/qs)
         U_nc  =  U
         if( CAMstfrac ) then
             call astG_RHU_single(U_nc,p,qv,landfrac(i),snowh(i),al_st_nc,G_nc)
         else
             call astG_PDF_single(U_nc,p,qv,landfrac(i),snowh(i),al_st_nc,G_nc)
         endif
         al_st = (1._r8-a_dc-a_sc)*al_st_nc
     else
         ai_st  = (1._r8-a_dc-a_sc)*ai0_st_nc_in(i)
         al_st  = (1._r8-a_dc-a_sc)*al0_st_nc_in(i)
     endif

     a_st  = max(ai_st,al_st)

     if( al_st .eq. 0._r8 ) then
         ql_st = 0._r8
     else
         ql_st = ql/al_st
         ql_st = min(qlst_max,max(qlst_min,ql_st)) ! PJR
     endif
     if( ai_st .eq. 0._r8 ) then
         qi_st = 0._r8
     else
         qi_st = qi/ai_st
     endif

     qi    = ai_st*qi_st
     ql    = al_st*ql_st

     T     = T0 - (latvap/cpair)*(ql0-ql) - ((latvap+latice)/cpair)*(qi0-qi)
     qv    = qv0 + ql0 - ql + qi0 - qi

   ! -------------- !
   ! Send to output !
   ! -------------- !

   T_out(i)  = T
   qv_out(i) = qv
   ql_out(i) = ql
   qi_out(i) = qi
   al_st_out(i) = al_st
   ai_st_out(i) = ai_st
   ql_st_out(i) = ql_st
   qi_st_out(i) = qi_st

   enddo 

   return
   end subroutine instratus_condensate

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine instratus_core( lchnk, icol, k, p,                      &
                              T0, qv0, ql0, qi0,                      &
                              a_dc, ql_dc, qi_dc,                     & 
                              a_sc, ql_sc, qi_sc, ai_st,              &
                              qcst_crit, Tmin, Tmax, landfrac, snowh, &
                              T, qv, ql, qi )

   ! ------------------------------------------------------ !
   ! Subroutine to find saturation equilibrium state using  ! 
   ! a Newton iteration method, so that 'qc_st = qcst_crit' !
   ! is satisfied.                                          !
   ! ------------------------------------------------------ !

   use time_manager,  only: is_first_step, get_nstep

   implicit none

   integer,  intent(in)  :: lchnk      ! Chunk identifier
   integer,  intent(in)  :: icol       ! Number of atmospheric columns
   integer,  intent(in)  :: k          ! Layer index

   real(r8), intent(in)  :: p          ! Pressure [Pa]
   real(r8), intent(in)  :: T0         ! Temperature [K]
   real(r8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]

   real(r8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(r8), intent(in)  :: Tmin       ! Minimum temperature system can have [K]
   real(r8), intent(in)  :: Tmax       ! Maximum temperature system can have [K]
   real(r8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(r8), intent(in)  :: landfrac   ! Land fraction
   real(r8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: T          ! Temperature [K]
   real(r8), intent(out) :: qv         ! Grid-mean water vapor [kg/kg]
   real(r8), intent(out) :: ql         ! Grid-mean LWC [kg/kg]
   real(r8), intent(out) :: qi         ! Grid-mean IWC [kg/kg]

   ! Local variables

   integer i                           ! Iteration index

   real(r8) muQ0, muQ
   real(r8) ql_nc0, qi_nc0, qc_nc0, qc_nc    
   real(r8) fice0, fice    
   real(r8) ficeg0, ficeg   
   real(r8) esat0
   real(r8) qsat0
   real(r8) dqcncdt, dastdt, dUdt
   real(r8) alpha, beta
   real(r8) U, U_nc
   real(r8) al_st_nc, G_nc
   real(r8) al_st

   ! Variables for root-finding algorithm

   integer j                          
   real(r8)  x1, x2
   real(r8)  rtsafe
   real(r8)  df, dx, dxold, f, fh, fl, temp, xh, xl
   real(r8), parameter :: xacc = 1.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0._r8,ql0-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0._r8,qi0-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0._r8,ql0+qi0-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   fice0  = 0._r8
   ficeg0 = 0._r8
   muQ0   = 1._r8

   ! ------------ !
   ! Root finding !
   ! ------------ !

   x1 = Tmin
   x2 = Tmax
   call funcd_instratus( x1, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         fl, df, qc_nc, fice, al_st )
   call funcd_instratus( x2, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         fh, df, qc_nc, fice, al_st )
   if((fl > 0._r8 .and. fh > 0._r8) .or. (fl < 0._r8 .and. fh < 0._r8)) then
       call funcd_instratus( T0, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                             a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                             qcst_crit, landfrac, snowh,                    &
                             fl, df, qc_nc, fice, al_st )
       rtsafe = T0 
       goto 10       
   endif
   if( fl == 0._r8) then
           rtsafe = x1
           goto 10
   elseif( fh == 0._r8) then
           rtsafe = x2
           goto 10
   elseif( fl < 0._r8) then
           xl = x1
           xh = x2
   else
           xh = x1
           xl = x2
   end if
   rtsafe = 0.5_r8*(x1+x2)
   dxold = abs(x2-x1)
   dx = dxold
   call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                         qcst_crit, landfrac, snowh,                        &
                         f, df, qc_nc, fice, al_st )
   do j = 1, 20
      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0._r8 .or. abs(2.0_r8*f) > abs(dxold*df) ) then
           dxold = dx
           dx = 0.5_r8*(xh-xl)
           rtsafe = xl + dx
           if(xl == rtsafe) goto 10
      else
           dxold = dx
           dx = f/df
           temp = rtsafe
           rtsafe = rtsafe - dx
           if (temp == rtsafe) goto 10
      end if
    ! if(abs(dx) < xacc) goto 10
      call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                            qcst_crit, landfrac, snowh,                        &
                            f, df, qc_nc, fice, al_st )
    ! Sep.21.2010. Sungsu modified to enhance convergence and guarantee 'qlst_min <  qlst < qlst_max'.
      if( qcst_crit < 0.5_r8 * ( qlst_min + qlst_max ) ) then
          if( ( qc_nc*(1._r8-fice) .gt.          qlst_min*al_st .and. &
                qc_nc*(1._r8-fice) .lt. 1.1_r8 * qlst_min*al_st ) ) goto 10
      else
          if( ( qc_nc*(1._r8-fice) .gt. 0.9_r8 * qlst_max*al_st .and. &
                qc_nc*(1._r8-fice) .lt.          qlst_max*al_st ) ) goto 10
      endif
      if(f < 0._r8) then
          xl = rtsafe
      else
          xh = rtsafe
      endif

   enddo

10 continue

   ! ------------------------------------------- !
   ! Final safety check before sending to output !
   ! ------------------------------------------- !

   qc_nc = max(0._r8,qc_nc)

   T  = rtsafe
   ql = qc_nc*(1._r8-fice) + a_dc*ql_dc + a_sc*ql_sc
   qi = qc_nc*fice + a_dc*qi_dc + a_sc*qi_sc
   qv = qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc))
   qv = max(qv,1.e-12_r8) 

   return
   end subroutine instratus_core

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine funcd_instratus( T, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0,   &
                               a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,  &
                               qcst_crit, landfrac, snowh,                     &
                               f, fg, qc_nc, fice, al_st ) 

   ! --------------------------------------------------- !
   ! Subroutine to find function value and gradient at T !
   ! --------------------------------------------------- !

   implicit none

   real(r8), intent(in)  :: T          ! Iteration temperature [K]

   real(r8), intent(in)  :: p          ! Pressure [Pa]
   real(r8), intent(in)  :: T0         ! Initial temperature [K]
   real(r8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]
   real(r8), intent(in)  :: fice0      ! 
   real(r8), intent(in)  :: muQ0       ! 
   real(r8), intent(in)  :: qc_nc0     ! 

   real(r8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(r8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(r8), intent(in)  :: landfrac   ! Land fraction
   real(r8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: f          ! Value of minimization function at T
   real(r8), intent(out) :: fg         ! Gradient of minimization function 
   real(r8), intent(out) :: qc_nc      !
   real(r8), intent(out) :: al_st      !
   real(r8), intent(out) :: fice       !

   ! Local variables

   real(r8) es
   real(r8) qs
   real(r8) dqsdT
   real(r8) dqcncdt
   real(r8) alpha
   real(r8) beta
   real(r8) U
   real(r8) U_nc
   real(r8) al_st_nc
   real(r8) G_nc
   real(r8) dUdt
   real(r8) dalstdt
   real(r8) qv

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   call qsat_water(T, p, es, qs, dqsdt=dqsdT)

   fice    = fice0 
   qc_nc   = (cpair/latvap)*(T-T0)+muQ0*qc_nc0       
   dqcncdt = (cpair/latvap) 
   qv      = (qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc)))
   alpha   = (1._r8/qs)
   beta    = (qv/qs**2._r8)*dqsdT 

   U      =  (qv/qs)
   U_nc   =   U
   if( CAMstfrac ) then
       call astG_RHU_single(U_nc,p,qv,landfrac,snowh,al_st_nc,G_nc)
   else
       call astG_PDF_single(U_nc,p,qv,landfrac,snowh,al_st_nc,G_nc)
   endif
   al_st   =  (1._r8-a_dc-a_sc)*al_st_nc 
   dUdt    = -(alpha*dqcncdt+beta)
   dalstdt =  (1._r8/G_nc)*dUdt
   if( U_nc .eq. 1._r8 ) dalstdt = 0._r8

   f  = qc_nc   - qcst_crit*al_st
   fg = dqcncdt - qcst_crit*dalstdt

   return
   end subroutine funcd_instratus

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine gridmean_RH( lchnk, icol, k, p, T, qv, ql, qi,       &
                           a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                           landfrac, snowh )

   ! ------------------------------------------------------------- !
   ! Subroutine to force grid-mean RH = 1 when RH > 1              !
   ! This is condensation process similar to instratus_condensate. !
   ! During condensation, we assume 'fice' is maintained in this   !
   ! verison for MG not for RK.                                    !
   ! ------------------------------------------------------------- !

   use time_manager,  only: is_first_step, get_nstep

   implicit none

   integer,  intent(in)    :: lchnk      ! Chunk identifier
   integer,  intent(in)    :: icol       ! Number of atmospheric columns
   integer,  intent(in)    :: k          ! Layer index

   real(r8), intent(in)    :: p          ! Pressure [Pa]
   real(r8), intent(inout) :: T          ! Temperature [K]
   real(r8), intent(inout) :: qv         ! Grid-mean water vapor [kg/kg]
   real(r8), intent(inout) :: ql         ! Grid-mean LWC [kg/kg]
   real(r8), intent(inout) :: qi         ! Grid-mean IWC [kg/kg]

   real(r8), intent(in)    :: a_dc       ! Deep cumulus cloud fraction
   real(r8), intent(in)    :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)    :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)    :: a_sc       ! Shallow cumulus cloud fraction
   real(r8), intent(in)    :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)    :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)    :: landfrac   ! Land fraction
   real(r8), intent(in)    :: snowh      ! Snow depth (liquid water equivalent)

   ! Local variables

   integer m                             ! Iteration index

   real(r8)  ql_nc0, qi_nc0, qc_nc0
   real(r8)  Tscale
   real(r8)  Tc, qt, qc, dqcdt, qc_nc    
   real(r8)  es, qs, dqsdT
   real(r8)  al_st_nc, G_nc
   real(r8)  f, fg
   real(r8), parameter :: xacc = 1.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0._r8,ql-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0._r8,qi-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0._r8,ql+qi-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   Tc    = T - (latvap/cpair)*ql
   qt    = qv + ql

   do m = 1, 20
      call qsat_water(T, p, es, qs, dqsdt=dqsdT)
      Tscale = latvap/cpair
      qc     = (T-Tc)/Tscale
      dqcdt  = 1._r8/Tscale
      f      = qs + qc - qt 
      fg     = dqsdT + dqcdt
      fg     = sign(1._r8,fg)*max(1.e-10_r8,abs(fg))
    ! Sungsu modified convergence criteria to speed up convergence and guarantee RH <= 1.
      if( qc .ge. 0._r8 .and. ( qt - qc ) .ge. 0.999_r8*qs .and. ( qt - qc ) .le. 1._r8*qs ) then
          goto 10
      endif
      T = T - f/fg
   enddo
 ! write(iulog,*) 'Convergence in gridmean_RH is not reached. RH = ', ( qt - qc ) / qs
10 continue

   call qsat_water(T, p, es, qs)
 ! Sungsu modified 'qv = qs' in consistent with the modified convergence criteria above.
   qv = min(qt,qs) ! Modified
   ql = qt - qv
   T  = Tc + (latvap/cpair)*ql

   return
   end subroutine gridmean_RH

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine positive_moisture( ncol, dt, qvmin, qlmin, qimin, dp, &
                                 qv,   ql, qi,    t,     qvten, &
                                 qlten,    qiten, tten,  do_cldice)

   ! ------------------------------------------------------------------------------- !
   ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
   ! force them to be larger than minimum value by (1) condensating water vapor      !
   ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
   ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
   ! Update final state variables and tendencies associated with this correction.    !
   ! If any condensation happens, update (s,t) too.                                  !
   ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': top layer, 'pver' : near-surface layer         ! 
   ! ------------------------------------------------------------------------------- !

   implicit none
   integer,  intent(in)     :: ncol
   real(r8), intent(in)     :: dt
   real(r8), intent(in)     :: dp(pcols,pver), qvmin(pcols,pver), qlmin(pcols,pver), qimin(pcols,pver)
   real(r8), intent(inout)  :: qv(pcols,pver), ql(pcols,pver), qi(pcols,pver), t(pcols,pver)
   real(r8), intent(out)    :: qvten(pcols,pver), qlten(pcols,pver), qiten(pcols,pver), tten(pcols,pver)
   logical, intent(in)      :: do_cldice
   integer   i, k
   real(r8)  dql, dqi, dqv, sum, aa, dum 

   tten(:ncol,:pver)  = 0._r8
   qvten(:ncol,:pver) = 0._r8
   qlten(:ncol,:pver) = 0._r8
   qiten(:ncol,:pver) = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         if( qv(i,k) .lt. qvmin(i,k) .or. ql(i,k) .lt. qlmin(i,k) .or. qi(i,k) .lt. qimin(i,k) ) then
             goto 10
         endif
      enddo
      goto 11
   10 continue
      do k = top_lev, pver    ! From the top to the 1st (lowest) layer from the surface
         dql = max(0._r8,1._r8*qlmin(i,k)-ql(i,k))

         if (do_cldice) then
         dqi = max(0._r8,1._r8*qimin(i,k)-qi(i,k))
         else
           dqi = 0._r8
         end if

         qlten(i,k) = qlten(i,k) +  dql/dt
         qiten(i,k) = qiten(i,k) +  dqi/dt
         qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
         tten(i,k)  = tten(i,k)  + (latvap/cpair)*(dql/dt) + ((latvap+latice)/cpair)*(dqi/dt)
         ql(i,k)    = ql(i,k) + dql
         qi(i,k)    = qi(i,k) + dqi
         qv(i,k)    = qv(i,k) - dql - dqi
         t(i,k)     = t(i,k)  + (latvap * dql + (latvap+latice) * dqi)/cpair
         dqv        = max(0._r8,1._r8*qvmin(i,k)-qv(i,k))
         qvten(i,k) = qvten(i,k) + dqv/dt
         qv(i,k)    = qv(i,k)    + dqv
         if( k .ne. pver ) then 
             qv(i,k+1)    = qv(i,k+1)    - dqv*dp(i,k)/dp(i,k+1)
             qvten(i,k+1) = qvten(i,k+1) - dqv*dp(i,k)/dp(i,k+1)/dt
         endif
         qv(i,k) = max(qv(i,k),qvmin(i,k))
         ql(i,k) = max(ql(i,k),qlmin(i,k))
         qi(i,k) = max(qi(i,k),qimin(i,k))
      end do
      ! Extra moisture used to satisfy 'qv(i,pver)=qvmin' is proportionally 
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture. 
      if( dqv .gt. 1.e-20_r8 ) then
          sum = 0._r8
          do k = top_lev, pver
             if( qv(i,k) .gt. 2._r8*qvmin(i,k) ) sum = sum + qv(i,k)*dp(i,k)
          enddo
          aa = dqv*dp(i,pver)/max(1.e-20_r8,sum)
          if( aa .lt. 0.5_r8 ) then
              do k = top_lev, pver
                 if( qv(i,k) .gt. 2._r8*qvmin(i,k) ) then
                     dum        = aa*qv(i,k)
                     qv(i,k)    = qv(i,k) - dum
                     qvten(i,k) = qvten(i,k) - dum/dt
                 endif
              enddo 
          else 
              write(iulog,*) 'Full positive_moisture is impossible in Park Macro'
          endif
      endif 
11 continue
   enddo
   return

   end subroutine positive_moisture

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine astG_PDF_single( U, p, qv, landfrac, snowh, a, Ga, orhmin )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   implicit none

   real(r8), intent(in)  :: U                     ! Relative humidity
   real(r8), intent(in)  :: p                     ! Pressure [Pa]
   real(r8), intent(in)  :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac              ! Land fraction
   real(r8), intent(in)  :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a                     ! Stratus fraction
   real(r8), intent(out) :: Ga                    ! dU/da
   real(r8), optional, intent(out) :: orhmin      ! Critical RH

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhwght
                            
   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

       if( freeze_dry ) then
           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   endif

   if (present(orhmin)) orhmin = rhmin

   return
   end subroutine astG_PDF_single

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine astG_PDF( U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   implicit none

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: U_in(pcols)           ! Relative humidity
   real(r8), intent(in)  :: p_in(pcols)           ! Pressure [Pa]
   real(r8), intent(in)  :: qv_in(pcols)          ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols)    ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)       ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a_out(pcols)          ! Stratus fraction
   real(r8), intent(out) :: Ga_out(pcols)         ! dU/da

   real(r8)              :: U                     ! Relative humidity
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8)              :: landfrac              ! Land fraction
   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8)              :: a                     ! Stratus fraction
   real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhwght
                            
   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out(:)  = 0._r8
   Ga_out(:) = 0._r8

   do i = 1, ncol

   U        = U_in(i)      
   p        = p_in(i)        
   qv       = qv_in(i)       
   landfrac = landfrac_in(i) 
   snowh    = snowh_in(i)    

   if( p .ge. premib ) then

       if( land(i) .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

       if( freeze_dry ) then
           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   endif

   a_out(i)  = a
   Ga_out(i) = Ga 

   enddo

   return
   end subroutine astG_PDF

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine astG_RHU_single( U, p, qv, landfrac, snowh, a, Ga, orhmin )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! CAM35 cloud fraction formula.                             !
   ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core !  
   ! For the other cases, I should re-define 'rhminl,rhminh' & !
   ! 'premib,premit'.                                          !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.                          !
   ! --------------------------------------------------------- !

   implicit none

   real(r8), intent(in)  :: U               ! Relative humidity
   real(r8), intent(in)  :: p               ! Pressure [Pa]
   real(r8), intent(in)  :: qv              ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac        ! Land fraction
   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a               ! Stratus fraction
   real(r8), intent(out) :: Ga              ! dU/da
   real(r8), optional, intent(out) :: orhmin ! Critical RH

   ! Local variables
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhdif                                 ! Factor for stratus fraction
   real(r8) rhwght

   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   if( p .ge. premib ) then

       if( land .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0.0_r8))**2) 
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif
       if( freeze_dry ) then
           a  = a*max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e10_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   endif

   if (present(orhmin)) orhmin = rhmin

   return
   end subroutine astG_RHU_single

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine astG_RHU( U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! CAM35 cloud fraction formula.                             !
   ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core !  
   ! For the other cases, I should re-define 'rhminl,rhminh' & !
   ! 'premib,premit'.                                          !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.                          !
   ! --------------------------------------------------------- !

   implicit none

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: U_in(pcols)           ! Relative humidity
   real(r8), intent(in)  :: p_in(pcols)           ! Pressure [Pa]
   real(r8), intent(in)  :: qv_in(pcols)          ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols)    ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)       ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a_out(pcols)          ! Stratus fraction
   real(r8), intent(out) :: Ga_out(pcols)         ! dU/da

   real(r8)              :: U                     ! Relative humidity
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8)              :: landfrac              ! Land fraction
   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8)              :: a                     ! Stratus fraction
   real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer  i
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhdif                                 ! Factor for stratus fraction
   real(r8) rhwght

   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out(:) = 0._r8
   Ga_out(:) = 0._r8

   do i = 1, ncol

   U        = U_in(i)      
   p        = p_in(i)        
   qv       = qv_in(i)       
   landfrac = landfrac_in(i) 
   snowh    = snowh_in(i)    

   if( p .ge. premib ) then

       if( land(i) .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0.0_r8))**2) 
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif
       if( freeze_dry ) then
           a  = a*max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e20_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       rhdif = (U-rhmin)/(1.0_r8-rhmin)
       a  = min(1._r8,(max(rhdif,0._r8))**2)
       if( (U.ge.1._r8) .or. (U.le.rhmin) ) then
            Ga = 1.e10_r8
       else          
            Ga = 0.5_r8*(1._r8-rhmin)*((1._r8-rhmin)/(U-rhmin))
       endif

   endif

   a_out(i)  = a
   Ga_out(i) = Ga 

   enddo

   return
   end subroutine astG_RHU

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine aist_single( qv, T, p, qi, landfrac, snowh, aist )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   use physconst,     only: rair

   implicit none
  
   real(r8), intent(in)  :: qv              ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T               ! Temperature
   real(r8), intent(in)  :: p               ! Pressure [Pa]
   real(r8), intent(in)  :: qi              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac        ! Land fraction
   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: aist            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables
   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) es, qs

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme


   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8 
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     call qsat_water(T, p, es, qs)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qs * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

   return
   end subroutine aist_single

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine aist_vector( qv_in, T_in, p_in, qi_in, landfrac_in, snowh_in, aist_out, ncol )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   use physconst,     only: rair

   implicit none
  
   integer,  intent(in)  :: ncol 
   real(r8), intent(in)  :: qv_in(pcols)       ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T_in(pcols)        ! Temperature
   real(r8), intent(in)  :: p_in(pcols)        ! Pressure [Pa]
   real(r8), intent(in)  :: qi_in(pcols)       ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols) ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)    ! Snow depth (liquid water equivalent)
   real(r8), intent(out) :: aist_out(pcols)    ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables

   real(r8) qv                              ! Grid-mean water vapor[kg/kg]
   real(r8) T                               ! Temperature
   real(r8) p                               ! Pressure [Pa]
   real(r8) qi                              ! Grid-mean ice water content [kg/kg]
   real(r8) landfrac                        ! Land fraction
   real(r8) snowh                           ! Snow depth (liquid water equivalent)
   real(r8) aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) qs
   real(r8) esat_in(pcols)
   real(r8) qsat_in(pcols)

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme

   integer i


   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     aist_out(:) = 0._r8
     esat_in(:)  = 0._r8
     qsat_in(:)  = 0._r8

     call qsat_water(T_in(1:ncol), p_in(1:ncol), &
          esat_in(1:ncol), qsat_in(1:ncol))
     
     do i = 1, ncol

     landfrac = landfrac_in(i)     
     snowh = snowh_in(i)   
     T = T_in(i)
     qv = qv_in(i)
     p = p_in(i)
     qi = qi_in(i)
     qs = qsat_in(i)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land(i) .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qs * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

     aist_out(i) = aist

     enddo

   return
   end subroutine aist_vector

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      real(r8) a(np,np),b(np,mp)
      real(r8) aa(np,np),bb(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,ii,jj,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      real(r8) big,dum,pivinv

      aa(:,:) = a(:,:)
      bb(:,:) = b(:,:)

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0._r8
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(iulog,*) 'singular matrix in gaussj 1'
                do ii = 1, np
                do jj = 1, np
                   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
                end do
                end do   
                call endrun
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0._r8) then
            write(iulog,*) 'singular matrix in gaussj 2'
            do ii = 1, np
            do jj = 1, np
               write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
            end do
            end do   
            call endrun
        endif 
        pivinv=1._r8/a(icol,icol)
        a(icol,icol)=1._r8
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0._r8
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue

      return
      end subroutine gaussj

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

end module cldwat2m_macro
