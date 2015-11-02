module unicon

! ---------------------------------------------------------- !
!                                                            !
!                 UNIFIED CONVECTION SCHEME                  !
!                                                            !
!                         ( UNICON )                         !
!                                                            ! 
!                        Developed By                        ! 
!                                                            !
!            Sungsu Park, AMP/CGD/NCAR, Boulder.             !
!                                                            !
!                         Aug.2010.                          !
!                                                            !
!  I hope my scheme can provide many people with happiness.  ! 
!  However,                                                  ! 
!                                                            !
!                      <<< WARNING >>>                       ! 
!                                                            !
!            No one is allowed to use this UNICON            !
!  without an explicit personal permission from Sungsu Park  !
!         until it is released as a open source code         !
!                 as a part of CAM6/CESM2.                   !
!                                                            !
! ---------------------------------------------------------- !

use shr_kind_mod,    only : r8 => shr_kind_r8, i4 => shr_kind_i4
use cam_history,     only : outfld
use shr_spfn_mod,    only : erfc => shr_spfn_erfc
use time_manager,    only : get_nstep
use cam_abortutils,  only : endrun
use cam_logfile,     only : iulog
use constituents,    only : qmin, cnst_get_type_byind, cnst_get_ind, cnst_name
use wv_saturation,   only : qsat
use unicon_utils,    only : unicon_utils_init, exnf, conden, slope, area_overlap, &
                            envcon_flux, prod_prep_up, evap_prep_dn, progup_thlqt, &
                            progup_uv, progup_wu2, compute_dp, buosort_downdraft, &
                            compute_pdf, compute_epsdelnod, buosorts_uw, findsp_single
#ifdef MODAL_AERO
use modal_aero_data
#endif

implicit none
private
save

public :: &
   unicon_init,       &
   compute_unicon

real(r8) :: xlv       !  Latent heat of vaporization
real(r8) :: xlf       !  Latent heat of fusion
real(r8) :: xls       !  Latent heat of sublimation
real(r8) :: cp        !  Specific heat of dry air
real(r8) :: zvir      !  rh2o/rair - 1
real(r8) :: r         !  Gas constant for dry air
real(r8) :: g         !  Gravitational constant
real(r8) :: p00       !  Reference pressure for exner function
real(r8) :: rovcp     !  R/cp

! -------------------------- !
!                            !
! Define Dynamic Parameters  !
!                            !
! -------------------------- !

! ---------------------------------- !
! 3 key parameters : au_base, Ro, mu !
! ---------------------------------- !

integer,  parameter :: nseg               =  1           !  Number of updraft segments [ # ]
integer,  parameter :: inorm              =  2           !  Either 1 ( force 'au_base' is preserved ) or 2 (force 'cmfu_base' is preserved) for various nseg. 
! integer,  parameter :: iprpback           =  1         !  Either 1 ( backward differencing ) or 0 ( centered differencing ) or -1 ( previous semi-analytical ) 
!                                                        !  or    -2 ( correct full analytical ) only in 'evap_prep_dn'.  
integer,  parameter :: iprd_prep          =  5           !  -1 : Forward  Numerical,  0 : Centered  Numerical,                     1 : Backward  Numerical,
                                                         !  -5 : Forward Analytical, 10 : Centered Analytical (not available yet), 5 : Backward Analytical (Part.II Paper)
integer,  parameter :: ievp_prep          =  1           !  -1 : Forward  Numerical,  0 : Centered  Numerical,                     1 : Backward  Numerical (Part.II Paper)
                                                         !  -5 : Forward Analytical, 10 : Centered Analytical (not available yet), 5 : Backward Analytical (not available yet)
integer,  parameter :: nacc               =  1           !  Number of accretion iterations [ # ]. Can be any integer n_icc >=1. 
integer,  parameter :: niter              =  1           !  Number of whole iterations [ # ]. Should be 1 or 2.
logical,  parameter :: dbsort_con         = .false.      !  If .true. (.false.), do continuous (discontinuous) downdraft buoyancy sorting. .false. is the previous default.
                                                         !  This is applied only for downdraft buoyancy sorting, not the buoyancy sorting of convective updraft that generates
                                                         !  the mixing downdraft. That is, the generation (or source) of mixing downdraft is still done by using 
                                                         !  previous default method using 'ithv_minE, mu_mix, offset_minE' specified in the below parameter sentences.
integer,  parameter :: ithv_minE          =  1           !  If 1 ( -1 ), do buoynacy sorting ( both updraft and downdraft ) using 'thv_minE' ( thvl_minE ).
real(r8), parameter :: mu_mix             =  0.5_r8      !  Minimum downdraft = 0 <= mu <= 1 = Maximum Downdraft. Due to the displacement of interface, it may be set to 1. 
                                                         !  In order to deepen PBL, we should use the largest value 1.  
                                                         !  For mixing        downdraft.
real(r8), parameter :: mu_top             =  0.5_r8      !  For top           downdraft. 
real(r8), parameter :: mu_area            =  0.5_r8      !  For area-velocity downdraft.

real(r8), parameter :: offset_minE        =  0._r8       !  Final 'thv_minE = thv_minE from mu + offset_minE'. Same for thvl_minE. [K]. Set 'offset_minE < 0' to reduce downdraft.

real(r8), parameter :: epsz_dn            =  0.e-4_r8    !  Lateral entrainment rate of downdraft [ 1 / z ]. 5.e-5 = R:1000 [m], 1.e-4 = R:500 [m], 2.e-4 = R:250 [m]. 
real(r8), parameter :: delz_dn            =  0.e-4_r8    !  Lateral  derainment rate of downdraft [ 1 / z ]. 5.e-4 = R:100  [m], 1.e-3 = R:50  [m], 2.e-3 = R:25  [m]. 
                                                         !  Jul.10.2011. These two are very impportant in controlling the properties of detrained airs at surface.
                                                         !  From ARM95, too small ( even 5.e-4 ) value produced too dry and cold detrained airs at surface. It seems that 
                                                         !  I should at least use 1.e-3.

real(r8), parameter :: eps_wk             =  0.e-5_r8    !  Lateral entrainment rate from non-wake to     wake area within PBL [ 1 / s ] 
real(r8), parameter :: del_wk             =  0.e-5_r8    !  Lateral  derainment rate from     wake to non-wake area within PBL [ 1 / s ]. 2.32e-4 is 3hr when awk_PBL = 0.5.
                                                         !  Previous study used 0.
                                                         !  This is valid only for below 'int_del_wk = 0'.

integer,  parameter :: i_budget_coldpool  =  6           !  If 0 (1, 2), use budget-inconsistent A068j ( budget consistent with M^j_U=0, budget consistent with M^j_G=0 ) cold pool formula.
                                                         !  If 3, use clevely approximately budget consistent modified from '0' option ( ** This '3' is the best ** ) 
                                                         !  If 4, budget-inconsistent A068j but with M^j_G=0 ( ** This is for sensitivity simulation for writing paper ).
                                                         !  If 5, use clevely approximately budget consistent method with with M^j_G=0 ( ** This '5' is the another best ** ) 
                                                         !  If 6, use clevely approximately, budget consistent, raw cold pool formula modified from '3' option ( ** This '6' is the best ** ) 
integer,  parameter :: i_energy_coldpool  =  1           !  If 2 (1, 0), use the full energy-consistent (partially consistent, previous default) cold pool formula.
real(r8), parameter :: eps_wk0            =  0.e-5_r8    !  Lateral entrainment rate from non-wake to     wake area within PBL [ 1 / s ]. Used only when i_energy_coldpool  =  1. 
real(r8), parameter :: del_wk0            =  0.e-5_r8    !  Lateral  derainment rate from     wake to non-wake area within PBL [ 1 / s ]. Used only when i_energy_coldpool  =  1. 
                                                         !  This must be positive to generate organized flow. 
                                                         !   5 hr = 5.56e-5,  6 hr = 4.63e-5,  7 hr = 3.97e-5,  8 hr = 3.47e-5,  9 hr = 3.09e-5, 10 hr = 2.78e-5,
                                                         !  15 hr = 1.85e-5, 20 hr = 1.39e-5, 1 day = 1.16e-5
real(r8), parameter :: b1                 =  15.0_r8     !  Multiplication factor for decorrelation time scale of meso-scale TKE. No unit.
                                                         !  This 'b1' is used only when 'i_energy_coldpool = 2' is chosen. 

integer, parameter  :: int_del_wk         =  0           !  If 1 (0), use internally computed 'del_wk_eff = c_del_wk * tmp1 * awk_PBL * cmf_u(kpblhm)' (specified del_wk above ), 
                                                         !  But above eps_wk is still used in any cases. 
                                                         !  Note that '0 <= c_del_wk <= 1' which is specified below. 
                                                         !  Hopely, this will have an impact to retard the diurnal cycle of convective precipitation over land.
       							 !  CAUTION : Since analytical integration is performed, not only 'del_wk_eff' but also the
       			                                 !  format of 'taui,_orgforce' should be changed together. Thus, below option of 'int_del_wk .eq. 1'
       		                                         !  is incomplete at this stage. This should be refined later.

real(r8), parameter :: c_del_wk           =  0._r8       !  When 'int_del_wk = 1' above, 'del_wk_eff = c_del_wk * tmp1 * awk_PBL * cmf_u(kpblhm)'. 
                                                         !  Note that '0 <= c_del_wk <= 1'.
                                                         !  This is valid only when int_del_wk = 1.

integer,  parameter :: icudist_tail       =  0           !  If 0 ( 1 ), use whole ( tail ) distribution for cumulus updrafts. 
                                                         !  0 : Good to simulate CIN stabilizing effect. 1 : Simulate only strong updrafts and so conceptually attractable.
                                                         !  But '1' also produces negative updraft buoyancy near cloud base. If 1, we may need to use larger Ro if sigmaR = 0.
                                                         !  In general, '1' produces non-better results, such as not to good and unstable u,v profile. Thus, recommended to use 0.
                                                         !  Feb.08.2013. I should always choose 'icudist_tail = 0' because 'inorm = 2' is only supported for icudist_tail = 0.
 
real(r8), parameter :: au_base_min_ocn    =  0.045_r8    !  Updraft fractional area at the launching interface [ 0 - 0.5 ] when org = 0. 
real(r8), parameter :: au_base_max_ocn    =  0.045_r8    !  Updraft fractional area at the launching interface [ 0 - 0.5 ] when org = 1.
real(r8), parameter :: au_base_min_lnd    =  0.03_r8     !  Updraft fractional area at the launching interface [ 0 - 0.5 ] when org = 0. 
real(r8), parameter :: au_base_max_lnd    =  0.03_r8     !  Updraft fractional area at the launching interface [ 0 - 0.5 ] when org = 1.

integer,  parameter :: iau_base_ocn       =  1           !  If '0', use above 'au_base_min,au_base_max' in computing au_base.
                                                         !  If '1', use above 'au_base = au_base_min * ( 1._r8 - cuorg * awk_PBL_max )'. In this case only above au_base_min is used.
                                                         !  Aug.03.2012. In association with the internal computation of 'cdelta_s,cdelta_w', it is definitely reasonable to use
                                                         !               iau_base_ocn = 1 instead of 0.
integer,  parameter :: iau_base_lnd       =  1           !  If '0', use above 'au_base_min,au_base_max' in computing au_base.
                                                         !  If '1', use above 'au_base = au_base_min * ( 1._r8 - cuorg * awk_PBL_max )'. In this case only above au_base_min is used.
                                                         !  Aug.03.2012. In association with the internal computation of 'cdelta_s,cdelta_w', it is definitely reasonable to use
                                                         !               iau_base_lnd = 1 instead of 0.
real(r8), parameter :: cadj_area_ocn      =  3._r8       !  The multiplication factor of 'au_base_min_ocn' : overturning adjustment associated with convective organization occurs
                                                         !  over the 'cadj_area_ocn * au_base_ocn'. It must be '1 <= cadj_area_ocn <= ( 1. / au_base_min_ocn )'.
                                                         !  If cadj_area_ocn = ( 1. / au_base_min_ocn ), overturning adjustment occurs in 'a_U' as in the original formulation.
                                                         !  Dec.19.2012. Previous use of cdelta_s_ocn = 4 is equivalent to using cadj_area_ocn = 5 if au_base_min_ocn = 0.05.
                                                         !  Apr.05.2012. If this is set to be too small (e.g., 1 or 2), the model crashes as expected.
real(r8), parameter :: cadj_area_lnd      =  3._r8       !  The multiplication factor of 'au_base_min_lnd' : overturning adjustment associated with convective organization occurs
                                                         !  over the 'cadj_area_lnd * au_base_lnd'. It must be '1 <= cadj_area_lnd <= ( 1. / au_base_min_lnd )'.
                                                         !  If cadj_area_ocn = ( 1. / au_base_min_lnd ), overturning adjustment occurs in 'a_U' as in the original formulation.
                                                         !  Dec.19.2012. Previous use of cdelta_s_lnd = 8 is equivalent to using cadj_area_lnd = 5 if au_base_min_lnd = 0.025.
                                                         !  Apr.05.2012. If this is set to be too small (e.g., 1 or 2), the model crashes as expected.

integer,  parameter :: icridis            =  1           !  If '1', use internal 'cridis = rlc * cush', but if '0', use specified 'cridis = cridis_in'.
                                                         !  In order to impose positive feedback as in uwshcu, recommend to use '1'  with the corresponding setting of 'rlc'.
                                                         !  Mar.11.2011. Since we are using multiple plume, 'cush' is meaningless and so good to use 'icridis = 0'
real(r8), parameter :: rlc                =  0.15_r8     !  Critical distance for updraft buoyancy sorting [ 0 - 1 ]. Critical distance = rlc * cush. 
                                                         !  Active only when icridis = 1. If rlc = -1._r8, use 'cridis  =  dz_m'.
                                                         !  May.16.2011. From BOMEX L30/L80, this 'rlc = -1' was clearly shown to be the source of resolution sensitivity.
                                                         !               Especially, when 'cuorg' was very small non-zero, this 'rlc = -1' showed very strange sensitivity
                                                         !               in the BOMEX L80 simulation.  
                                                         !               Thus, I MUST NOT USE 'rlc = -1.' in any case. 
                                                         !               However, in order to impose a positive feedback from shallow to deep convection ( i.e., less mixing
                                                         !               for deep convective cases ), it is good not to use 'cridis_in = 1.e8' but to use 'icridis = 1' and
                                                         !               'rlc = 0.1' etc. Thus, from May.16.2011 today, I must use 'icridis = 1 & rlc = 0.1 etc'.   
real(r8), parameter :: cridis_in          =  1.e8_r8     !  Critical distance for updraft buoyancy sorting [ m ]. Default was 1.e8. 
                                                         !  Active only when icridis = 0.
integer,  parameter :: i_downloading      =  0           !  If '1' ( '0' ), include ( exclude ) precipitating-condensate-loading in computing the buoyancy of convective downdraft.
                                                         !  This will help to increase downdraft vertical velocity.
                                                         !  Oct.27.2011. Restored to zero following the '016' case. 
real(r8), parameter :: vfall_rain         =  1.e1_r8     !  Fall speed of rain droplet [ m/s ]. This is used for computing (1) rain mixing ratio and (2) evaporation time scale
                                                         !  of convective rain within environment in each layer for computing tkePBLorgEV.
real(r8), parameter :: vfall_snow         =  1.e1_r8     !  Fall speed of snow droplet [ m/s ]. This is used for computing (2) snow mixing ratio and (2) evaporation time scale
                                                         !  of convective snow within environment in each layer for computing tkePBLorgEV.   

real(r8), parameter :: prepminPBLH_org    =  0.0_r8      !  Minimum precipitation flux at the PBL top height interface required for initiating convective organization [ mm/day ].
                                                         !  This corresponds to 'individual' updraft segment, not the whole mean of various updrafts.
                                                         !  This should be carefully defined not to be sensitive to the choice of 'nseg'.
                                                         !  In order to completely remove the sensitivity to nseg, this should be set to zero in principle.
                                                         !  Setting this value to be positive (regardless of how small it is) will result in the sensitivity to nseg. 
                                                         !  Feb.06.2013. Now, without having sensitivity to nseg, we can use small positive value for this with one line modification
                                                         !               in the main program. See the main program part with 'prepminPBLH_org'. 

logical,  parameter :: iorg_adv           = .true.       !  If .true. (.false.), advect (do not advect) horizontal heterogeneity information associated with organization.

logical,  parameter :: orgfeedback_off    = .false.      !  If .true. (.false.), turn-off (turn-on) convective organization. Default is .false.
                                                         !  This switch with 'false' is used for doing sensitivity simulation without convective organization feedback.

real(r8), parameter :: norm_sgh           =  1.e3_r8     !  Normalization height for computing 'a_oro = sgh30 / norm_sgh' [ m ]
real(r8), parameter :: a_oro_max          =  0._r8       !  Maximally allowed forbidden area by surface orography [ fraction ]
                                                         !  Set this to be zero if I want to turn-off the orographic effect.
real(r8), parameter :: awk_PBL_min        =  0.05_r8     !  Minimum wake area used only for computing 'eps_wk_eff,del_wk_eff' to allow initial development of wake [fraction] > 0.
                                                         !  This should be set to be positive in order to allow wake development for a given set of non-zero 'eps_wk,del_wk'.
                                                         !  If we set 'eps_wk = del_wk = 0' in the parameter sentence, this does not do anything. 
                                                         !  Sep.16.2011. This should be appropriately set such that it prevents the development of weak wake but allows 
                                                         !               the development of reasonably strong wake. By doing this, we may be able to simulate diurnal cycle
                                                         !               very well - my impression is that we can control diurnal cycle by appropriately chosing 
                                                         !                   (1) 'awk_PBL_min', 
                                                         !                   (2) 'delta_thv_wc' 
                                                         !                   (3) 'sigmaR_min_lnd', 'sigmaR_max_lnd',  
                                                         !                   (4) 'au_base_max' , 'au_base_max'
                                                         !                   (5) 'iorg_detrain'
                                                         !               If awk_PBL = 0.05 ( 0.1, 0.5 ), then 'awk_PBL / ( 1 - awk_ PBL )' = 21, 11, 4. 
real(r8), parameter :: cdrag              =  1.5e-3_r8   !  Surface drag coefficient for computing damping time scale of wake within PBL [ no unit ]
                                                         !  Sep.16.2011. Since this should also reflect the neglected dffect of small entrainment flux at the PBL top
                                                         !               in the wake area due to enhanced stratification at the PBL top, we should use smaller value
                                                         !               than the typical allowed value of ~ 1.5e-3.
real(r8), parameter :: delta_thv_wc       = -0.1_r8      !  Critical thv difference between 'wake' and 'grid-mean' averaged over the PBL in order to be identified as the 'wake'. 
                                                         !  Must be negative value. [ K ]. This is a general wake selection parameter.
                                                         !  The 'wake' with thv offset less than this value ( e.g., -0.05 K ) will be identified as 'non-wake'. 
                                                         !  Sep.12.2011. Currently '-0.01' produces the best results. The '-0.05' produced too small 'cuorg' and unreasonably
                                                         !               large wake spreading velocity. Thus, let's use '-0.01'. 
real(r8), parameter :: kw_omega           =  1.414_r8    !  Control 'sigma_w = kw_omega * sqrt(tke_omega)'.
                                                         !  Sep.22.2011. Use 0.82 assuming isotropic meso-scale turbulence similar to kw.
real(r8), parameter :: kstar              =  0.1_r8      !  Compute tke_omega = kstar * tke_omega_max. In principle, 0 < kstar < 1.
                                                         !  By setting zero, we can turn-off density current parameterization on delta_w_PBL, i.e., delta_w_PBL = 0.
                                                         !  Oct.27.2011. Set to the kstar = 1 instead of 0.5.
                                                         !  Aug.03.2012. In association with the internally computed 'cdelta_s,cdelta_w', I can use more 
                                                         !               reasonably smaller values ( 0.05, 0.1 or 0.2 ) for kstar.  

integer,  parameter :: iorg_ent           =  1           !  If '1', use the detrained airs from previous time step as part of environmental airs for lateral mixing.
                                                         !  If '0', use the mean environmental airs at the current time step for lateral mixing.
integer,  parameter :: iorg_detrain       =  1           !  Choose the detrained airs that will be used for organized entrainment at the next time step. Default was 1.
                                                         !  This switch is active only when the above 'iorg_ent = 1'.
                                                         !  '1' : Detrained Updraft + Detrained Downdraft
                                                         !  '2' : Detrained Updraft
                                                         !  '3' : Detrained Downdraft
                                                         !  '4' : Updraft
                                                         !  '5' : Detrained Updraft + Detrained Downdraft + Updraft
                                                         !  '6' : Detrained Updraft +                       Updraft
                                                         !  '7' : Detrained Downdraft +                     Updraft

integer,  parameter :: i_detrain          =  1           !  If '0' : Convectively detrained air contains mixing environmental airs (previous old). This cause inconsistency between 
                                                         !           flux-convergence and subsidence-detrainment formula.
                                                         !  If '1' : Convectively detrained air is defined only using convective updraft (not with mixing environmental air) and convective downdraft airs 
                                                         !           This impose a full constency between 'flux-convergence' and 'subsidence-detrainment' formula.
                                                         !           It is highly recommended to use '1'. 

real(r8), parameter :: fac_org_ent        =  1._r8       !  Scale factor for computing 'org_ent = cuorg * fac_org_ent'. After this, we impose limits [0,1] on org_ent.
                                                         !  Active only when iorg_ent = '1'.
real(r8), parameter :: fac_org_rad        =  1._r8       !  Scale factor for computing 'org_rad = cuorg * fac_org_rad'. After this, we impose limits [0,1] on org_rad.

real(r8), parameter :: orp                =  1._r8       !  Power for the plume radius associated with organization. 0 < orp.
                                                         !  If orp = 1 (0.5, 2), R is a linear (square root, square) function of organization. In order to reduce SWCF over the
                                                         !  Panama, I strongly recommend to use 'orp' smaller than 1, e.g., 0.5. 
                                                         !  Apr.01.2013. If 'orp = -1', use 'sinusoidal' function.

real(r8), parameter :: Ro_min_ocn         =  100._r8     !  Minimum intercept updraft plume radius at surface at alpha = 0 [ m ]
                                                         !  This is a default value when org_Rad = 0. 
                                                         !  May.21.2011 : We may use slightly a large value in future.
real(r8), parameter :: Ro_max_ocn         = 4000._r8     !  Maximum intercept updraft plume radius at surface at alpha = 0 [ m ]
                                                         !  This is a default value when org_Rad = 1.
                                                         !  May.21.2011 : We may use slightly a large value in future.
real(r8), parameter :: sigmaR_min_ocn     =  100._r8     !  Standard deviation of updraft plume radius ( -infinity < alpha < infinity ) at surface [ m ]
                                                         !  This is a default value when org_Rad = 0.
                                                         !  May.21.2011 : We may use slightly a large value in future.
real(r8), parameter :: sigmaR_max_ocn     =  100._r8     !  Standard deviation of updraft plume radius ( -infinity < alpha < infinity ) at surface [ m ]
                                                         !  This is a default value when org_Rad = 1.

real(r8), parameter :: Ro_min_lnd         =  100._r8     !  Minimum intercept updraft plume radius at surface at alpha = 0 [ m ]
                                                         !  This is a default value when org_Rad = 0. 
                                                         !  May.21.2011 : We may use slightly a large value in future.
real(r8), parameter :: Ro_max_lnd         =10000._r8     !  Maximum intercept updraft plume radius at surface at alpha = 0 [ m ]
                                                         !  This is a default value when org_Rad = 1.
                                                         !  May.21.2011 : We may use slightly a large value in future.
real(r8), parameter :: sigmaR_min_lnd     =  100._r8     !  Standard deviation of updraft plume radius ( -infinity < alpha < infinity ) at surface [ m ]
                                                         !  This is a default value when org_Rad = 0.
                                                         !  May.21.2011 : We may use slightly a large value in future.
real(r8), parameter :: sigmaR_max_lnd     =  100._r8     !  Standard deviation of updraft plume radius ( -infinity < alpha < infinity ) at surface [ m ]
                                                         !  This is a default value when org_Rad = 1.
                                                         !  Aug.15.2011. Only this parameter is doubled compared to the ocean to compensate for the
                                                         !               relative small value of convective organization over the land. 

real(r8), parameter :: Ro_eps0            =  25._r8      !  Minimum updraft plume for computing mixing rate eps0 [ m ]. Originally 10 [ m ].
                                                         !  Aug.15.2011. Originally 10 but relaxed to 100 following 008a.
real(r8), parameter :: eta2               =  0.1_r8      !  Frac. of convective precipitation evaporated within downdraft (=ovc(a_p,a_d)/a_p) [ 0 (zero) - 1 (whole precipitation flux) ]
                                                         !  In order to allow evaporation within environment, good to use eta2 < 1. Recommend to use eta2 = 0.5 if specified.
                                                         !  Feb.08.2013. Since precipitation area is not necessarily coinciding with downdraft area (=ovc(a_p,a_d)/a_p), this value
                                                         !               should be smaller than 1. For example, eta2 = 0.5 seems to be more reasonable choice that 0.99.   
                                                         !  Apr.05.2013. Larger 'eta2' reduce MJO noise but increases PREH2O.
real(r8), parameter :: beta2              =  1.0_r8      !  Tilting parameter of convective updraft plume with height [ 0 (minimum tilting) - 1 (maximum tilting) ]
                                                         !  This control overlappings between 'evaporation ( or precipitation ) area' and 'wake area'.
                                                         !  Sep.13.2011. In principle, we can ( should ) set 'beta1 = beta2' since the same physical process is 
                                                         !               controlling both 'beta1' and 'beta2'.
                                                         !               In future, we can set this as a function of low-level wind shear.
                                                         !  Apr.25.2012. This is nothing to do with the choice of 'i_ovp' above. Only 'beta1' is related to 'i_ovp = 0'.

real(r8), parameter :: beta2_st           =  0.0_r8      !  Same as the above 'beta2' but for the stratiform precipitation part.
                                                         !  Aug.08.2013. Overlapping parameter between 'stratiform evaporation area' and 'cold-pool area'
                                                         !               0 (random overpap) <= beta2_st <= 1 (maximum overlap)

real(r8), parameter :: sigma_wo           =  0.0_r8      !  Background value as 'sigma_w = sigma_wo + kw * sqrt(tkes)'. Should be small positive value. > 0.
                                                         !  Oct.16.2011. This is designed to improve the timing of diurnal cycle of precipitation. 
                                                         !               For obtaining the desirable effect, this should be accompanied by the use of smaller kw.    


real(r8), parameter :: kw_min_ocn         =  0.35_r8     !  Control 'sigma_w = sigma_wo + kw * sqrt(tkes)' when org = 0 over ocean.
real(r8), parameter :: kw_max_ocn         =  0.35_r8     !  Control 'sigma_w = sigma_wo + kw * sqrt(tke1)' when org = 1 over ocean.
real(r8), parameter :: kw_min_lnd         =  0.35_r8     !  Control 'sigma_w = sigma_wo + kw * sqrt(tkes)' when org = 0 over land.
real(r8), parameter :: kw_max_lnd         =  0.35_r8     !  Control 'sigma_w = sigma_wo + kw * sqrt(tke1)' when org = 1 over land.

real(r8), parameter :: PGFc_up            =  0.9_r8      !  Effect of horizontal PGF on the vertical change of   updraft horizonal momentum [ 0 - 1 ]
real(r8), parameter :: PGFc_dn            =  0.9_r8      !  Effect of horizontal PGF on the vertical change of downdraft horizonal momentum [ 0 - 1 ]

integer,  parameter :: mclimit            =  1           !  If '1' ( '0' ), impose (not impose ) 'ql + qi > criqc' at the top interface after precipitation fall-out.
                                                         !  This is only valid when microcu = 1 above.

real(r8), parameter :: caer               =  0.15_r8     !  Wet scavenging efficiency of aerosols within convective updraft [ no unit or fraction ] ( 0 < caer < 1 ).
 
real(r8), parameter :: criqc_lnd          =  6.5e-4_r8   !  Critical in-cumulus LWC for the formation of precipitation over land  [ kg / kg ]. UWShCu used 7.e-4.    
                                                         !  Active when microcu = 0 or 1
                                                         !  Apr.11.2012. Land-Sea contrast is added to consider aerosol effects. 
                                                         !  From ARM97, this value over land is definitely smaller than over ocean. 
real(r8), parameter :: criqc_ocn          =  6.5e-4_r8   !  Critical in-cumulus LWC for the formation of precipitation over ocean [ kg / kg ]. UWShCu used 7.e-4.    
                                                         !  Active when microcu = 0 or 1
                                                         !  Apr.11.2012. Land-Sea contrast is added to consider aerosol effects. 
                                                         !  From TOGAII, this value over ocean is definitely larger than over land. 
real(r8), parameter :: c0_ac_lnd          =  1.0e-3_r8   !  Auto-conversion efficiency of cloud liquid ( ice ) to rain ( snow ) over land.
                                                         !  No unit between 0 and 1 (UWShCu) when microcu = 0, or '1 / m' when microcu = 1 ]
                                                         !  Active only when microcu = 0 or 1, both of them are vertical-resolution insensitive 
                                                         !  if microcu = 0 (1) is used with upward (centered) computation as in the current code.  
                                                         !  When microcu = 0, 0<= c0_ac <=1 where UWShCu uses 1. 
                                                         !  When microcu = 1, in CAM5.1, c0_ac = 0.0059 ( LAND ) and 0.0450 ( OCEAN ) 
  
real(r8), parameter :: c0_ac_ocn          =  1.0e-3_r8   !  Auto-conversion efficiency of cloud liquid ( ice ) to rain ( snow ) over ocean.
                                                         !  No unit between 0 and 1 (UWShCu) when microcu = 0, or '1 / m' when microcu = 1 ]
                                                         !  Active only when microcu = 0 or 1, both of them are vertical-resolution insensitive 
                                                         !  if microcu = 0 (1) is used with upward (centered) computation as in the current code.  
                                                         !  When microcu = 0, 0<= c0_ac <=1 where UWShCu uses 1.
                                                         !  When microcu = 1, in CAM5.1, c0_ac = 0.0059 ( LAND ) and 0.0450 ( OCEAN ) 

real(r8), parameter :: kevp_rain_dn_lnd   =  2.e-5_r8    !  Evaporation efficiency of rain flux within downdraft over land [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.
real(r8), parameter :: kevp_snow_dn_lnd   =  2.e-5_r8    !  Evaporation efficiency of snow flux within downdraft over land [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.
real(r8), parameter :: kevp_rain_lnd      =  2.e-5_r8    !  Evaporation efficiency of rain flux within environment over land [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.
real(r8), parameter :: kevp_snow_lnd      =  2.e-5_r8    !  Evaporation efficiency of snow flux within environment over land [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.

real(r8), parameter :: kevp_rain_dn_ocn   =  2.e-5_r8    !  Evaporation efficiency of rain flux within downdraft over ocean [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.
real(r8), parameter :: kevp_snow_dn_ocn   =  2.e-5_r8    !  Evaporation efficiency of snow flux within downdraft over ocean [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.
real(r8), parameter :: kevp_rain_ocn      =  2.e-5_r8    !  Evaporation efficiency of rain flux within environment over ocean [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.
real(r8), parameter :: kevp_snow_ocn      =  2.e-5_r8    !  Evaporation efficiency of snow flux within environment over ocean [ ( kg m^-2 s^-1 )^(-1/2) s^-1 ]. UWShCu used 1.e-6.

real(r8), parameter :: c0                 =  0.2_r8      !  For fractional mixing rate of eps0 when inverse-R formula is used [ no unit ]
integer,  parameter :: i_eps0             =  1           !  Parameter for eps0 : 
                                                         !    '0' : The original ramped eps0 as a function of ql_u + qi_u. 
                                                         !    '1' : Evaporative enhancement of mixing as a function of 'ql_u + qi_u' and 'rh_eg'.
                                                         !    '2' : The full eps0 = c0 / R / eeps following laboratory experiment.
real(r8), parameter :: cevpeps0           =  1.0_r8      !  Evaporative enhancement of mixing. This is valid if i_eps0 = 0 or 1 [ no unit ] 
                                                         !  If 'i_eps0 = 0' : 1 /  3 /  5 (  1 maybe optimum )
                                                         !  If 'i_eps0 = 1' : 5 / 10 / 15 ( 10 maybe optimum )  

integer,  parameter :: i_dnmixing         =  0           !  If '0' ( '1', '2', '3' ), convective downdraft is mixed with mean environmental air ( 
                                                         ! '1' : detrained airs at the previous time step, 
                                                         ! '2' : updraft at the same current time step, 
                                                         ! '3' : below PBL top, mix with 'wake' airs but above PBL top, mix with
                                                         !  the same mixing environmental airs as convective updraft ). 
                                                         !  Aug.31.2011. In association with wake parameterization, I added '3' option. 
                                                         !               However, it does not consider the amount of available mixing environmental
                                                         !               airs and 'wake' airs in computation. This is approximation but 
                                                         !               this '3' is probably the most realistic approach even though it involves
                                                         !               an inevitable approximation. 
                                                         !  Nov.15.2011. In addition, only some downdrafts falls into wake area. Thus, it is good to
                                                         !               use i_dnmixing = 0.
                                                         !  Sep.15.2011. Recommend to use '0', i.e., mean environmental airs for convective downdraft as opposed to convective updraft.

real(r8), parameter :: rbuoy_min          =  0.33_r8     !  Minimum Buoyancy coefficient [ no unit ]. Default was 0.37.
                                                         !  Aug.18.2011. Originally it was 0.33 but restored to 0.37 based on my previous LES analysis. 
real(r8), parameter :: rbuoy_max          =  1.00_r8     !  Maximum Buoyancy coefficient [ no unit ]. Default was 0.37.
                                                         !  Aug.18.2011. Originally it was 1 but restored to the value of 0.37 same as rbuoy_min. 
                                                         !  This is for enhancing vertical velocity perturbation when org is developed, so that deep convection is well developed.
                                                         !  It might be important to use the same value for rbuoy_min and rbuoy_max to obtain the effect of reasonable org.     
real(r8), parameter :: rdrag              =  2.0_r8      !  Drag coefficient [ no unit ]
real(r8), parameter :: rjet               =  0.0_r8      !  Jet coefficient associated with detrainment [ 0 - 1, no unit ]. 0 : No jet effect, 1 : Maximum jet effect. 
real(r8), parameter :: R_buo              = 100._r8      !  Scaling updraft plume radius where 0 < a(R)=0.5786 <=1.
                                                         !  Mar.27.2012. This is added on this day.
real(r8), parameter :: xc_min             =  0._r8       !  Minimum critical mixing fraction from the buoyancy sorting. 
                                                         !  Oct.16.2011. This must be set to be a very small positive value when i_eps0 = 1 is chosen. 1.e-3_r8. [ 0-1 ].
                                                         !  Aug.18.2011. Originally it was 0.15 but is reduced down to 0.1. In principle, I can use any positive value even small.
real(r8), parameter :: xc_max             =  1._r8       !  Maximum critical mixing fraction from the buoyancy sorting. Default = 1 [ 0-1 ].
                                                         !  Aug.18.2011. Originally 0.85 but restored to 1.0 to remove any unclear limitation.

real(r8), parameter :: droprad_liq        =  10.e-6_r8   !  Effectie droplet radius of detrained liquid [ m ]  
real(r8), parameter :: droprad_ice        =  85.e-6_r8   !  Effectie droplet radius of detrained    ice [ m ]

real(r8), parameter :: density_liq        =  997._r8     !  Density of cloud liquid droplets [ kg/m3 ]  
real(r8), parameter :: density_ice        =  500._r8     !  Density of cloud ice    crystals [ kg/m3 ]

real(r8), parameter :: droprad_rain       =  10.e-6_r8   !  Effectie droplet radius of rain [ m ]  
real(r8), parameter :: droprad_snow       =  85.e-6_r8   !  Effectie droplet radius of snow [ m ]

real(r8), parameter :: density_rain       =  1000._r8    !  Density of rain [ kg/m3 ]  
real(r8), parameter :: density_snow       =  250._r8     !  Density of snow [ kg/m3 ]

real(r8), parameter :: epsz0_max          =  0.05_r8     !  Maximum effective mixing rate allowed: (M(p+dp)/M(p)-1)/dp = (exp(eps0*dp)-1)/dp <= (epsz0_max/rho/g) [1/m].
                                                         !  Note that this has a unit of [1/m] not [1/Pa] since I should use 'densitiy' in the main body.                    
                                                         !  This is added on Nov.18.2013, replacing the above fmix_frac.
                                                         !  This value is computed by using the following formula of 
                                                         !             epsz0_max [1/m] ~ 0.5 / Ro_minimum [m] 
                                                         !  This may be potentially associated with 'Ro_eps0' but let's use 'epsz0_max' separately.   
                                                         !  This is consistently applied both for convective updraft and downdraft, preventing the mass flux from 
                                                         !  hugely increasing during vertical displacement in any layer.
integer,  parameter :: exp_cmf            =  1           !  If 1 (2,3), do original exponential ( simplified linear, combination ) computation of vertical variation of cmf_u and cmf_d.
                                                         !  If this is set to 2, above 'fmix_max, fmix_frac' are not used.
                                                         !  This linear option is critical important in the coarse resolution simulation to prevent unreasonably large
                                                         !  exponential growth of mass flux and so model crash. This inevitably induces inconsistency with the 
                                                         !  corresponding computation of 'thl,qt,u,v,q' and 'w' (so that affect vertical evolution of plume radius R),
                                                         !  but in order to use physically reasonable single 'eps0' in all the vertical prognostic equations with stable state,
                                                         !  the use of this linear option is very important. 
                                                         !  Jan.30.2013. With a new formulation on this day on imposing 'eps0(m), eps_dn, del_dn' with fmix_frac = ln(10) = 2.3026, 
                                                         !               I 'must' use exp_cmf = 1 ( full physical exponential function ) for full physical consistency of the model.
                                                         !  Feb.06.2013. Always choose 'exp_cmf = 1' and removes this option.
real(r8), parameter :: alpha_max          =  2._r8       !  Upper limit of mixing parameter of updraft mass flux PDF [ no unit ]
real(r8), parameter :: cmfmin             =  1.e-5_r8    !  Minimum updraft mass flux for identification as the non-detached updraft at the top interface [ kg / s / m^2 ]

real(r8), parameter :: au_max             =  1.e-1_r8    !  Maximum   updraft fractional area [ no unit ]
real(r8), parameter :: wumin              =  1.e-1_r8    !  Minimum   updraft vertical velocity > 0 [ m / s ]
real(r8), parameter :: wumax              =  20._r8      !  Maximum   updraft vertical velocity > 0 [ m / s ]

real(r8), parameter :: wdmin              =  1.e-1_r8    !  Minimum downdraft vertical velocity > 0 [ m / s ]. This also influences evaporation within downdraft. May increase to 0.3.

real(r8), parameter :: nonzero            =  1.e-20_r8   !  Non-zero minimal positive constant [ no unit ]
real(r8), parameter :: unity              =  0.9999_r8   !  Constant close to 1 but smaller than 1 [ no unit ]
real(r8), parameter :: thv_ref            =  300._r8     !  Reference virtual potential temperature for buoyancy computation [ K ]

integer,  parameter :: iup_par            =  4           !  Partitioning of convective surface updraft flux into other layers
                                                         !  1 : Lowest Layer, or equivalently, No Partitioning
                                                         !  2 : Minimum of 'PBL Top' and 'Cumulus Top' 
                                                         !  3 : PBL Layers ( this degrades diurnal cycle )
                                                         !  4 : Cumulus Layers ( ** the best option ** )
                                                         !  5 : Entire Layers
                                                         !  Added on Mar.20.2014 as the most general functionality.

integer,  parameter :: idn_par            =  1           !  Partitioning of convective surface downdraft flux into other layers
                                                         !  1 : Lowest Layer, or equivalently, No Partitioning
                                                         !  2 : Minimum of 'PBL Top' and 'Cumulus Top' ( ** the best option both for cold-pool and plumes in stable env. ** )  
                                                         !  3 : PBL Layers
                                                         !  4 : Cumulus Layers
                                                         !  5 : Entire Layers
                                                         !  Added on Mar.20.2014 as the most general functionality.

integer,  parameter :: islope_on_thlqttr  =  1           !  If 1 ( 0 ), turn-on ( off ) environmental profile reconstruction of 'thl, qt and tracers' in each layer. 
                                                         !  Strongly recommend to use 0 to obtain stable results and prevent unreasonable reversal of environmental buoyancy at the interface.
                                                         !  In the CAM5, simulation with '1' cuases model crash due to FINSDP_SINGLE in unicon.F90 --> But the cause was not due to this. 
                                                         !  Apr.03.2011 : I should do sensitivity test on this using BOMEX.
                                                         !  Sep.11.2011 : For full consistency with symmetric turbulence transport scheme, I must use 'islope_on_thlqttr = 1 ( not 0 )'
                                                         !                without no question in this asymmetric turbulence transport scheme.
integer,  parameter :: islope_on_uv       =  1           !  If 1 ( 0 ), turn-on ( off ) environmental profile reconstruction of 'u,v' in each layer. 
                                                         !  Recommend to use 0 to reasonably treat diabatic 'PGFc' effect similar to UW. 
                                                         !  Apr.03.2011 : I must use '1' for treating PGFc effect.
                                                         !  Sep.11.2011 : For full consistency with symmetric turbulence transport scheme, I must use 'islope_on_uv = 1 ( not 0 )'
                                                         !                without no question in this asymmetric turbulence transport scheme.
integer,  parameter :: iflux_env          =  1           !  Use the UW ( 0 ) or Park's reconstructed environmental profile ( 1 ) for flux computation.
                                                         !  Recommend to always use 1 since it is conceptually reasonable.
                                                         !  Sep.11.2011. For full conceptual consistency and consistent performance, I must absolutely use 'iflux_env = 1'.
integer,  parameter :: kiss               =  0           !  Launching interface : 0 ( 'surface' ) or 1 ( next interface ). In order to reduce model sensitivity to the thickness of 
                                                         !  the lowest layer, recommend to choose 0, which should be accompanied by sigfac = 1.
                                                         !  Sep.11.2011. For full conceptual consistency and coherence with 'ipartition = 1' choice above, 
                                                         !               I must absolutely use 'kiss = 0'.  
real(r8), parameter :: sigfac             =  1.0_r8      !  Reduction of surface flux from the surface to the next interface [ 0-1 ]
                                                         !  When kiss = 0, should be 1 but when kiss = 1, choose any reasonable value between [ 0 - 1 ].
                                                         !  Sep.11.2011. Needless to say, I must use 'sigfac = 1._r8' since I will always use kiss = 0.

!==================================================================================================
contains
!==================================================================================================

subroutine unicon_init(xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in)

   real(r8), intent(in) :: xlv_in     !  Latent heat of vaporization
   real(r8), intent(in) :: xlf_in     !  Latent heat of fusion
   real(r8), intent(in) :: cp_in      !  Specific heat of dry air
   real(r8), intent(in) :: zvir_in    !  rh2o/rair - 1
   real(r8), intent(in) :: r_in       !  Gas constant for dry air
   real(r8), intent(in) :: g_in       !  Gravitational constant

   call unicon_utils_init(&
      xlv_in, cp_in, xlf_in, zvir_in, r_in,                     &
      g_in, droprad_liq, droprad_ice, density_liq, density_ice, &
      mclimit)

   xlv   = xlv_in
   xlf   = xlf_in
   xls   = xlv + xlf
   cp    = cp_in
   zvir  = zvir_in
   r     = r_in
   g     = g_in
   p00   = 1.e5_r8
   rovcp = r/cp

end subroutine unicon_init

!==================================================================================================

subroutine compute_unicon( mix            , mkx           , iend          , ncnst         , dt            ,                 &
                           ps0_in         , zs0_in        , p0_in         , z0_in         , dp0_in        , dpdry0_in     , &
                           t0_in          , qv0_in        , ql0_in        , qi0_in        , tr0_in        ,                 & 
                           u0_in          , v0_in         , ast0_in       , tke0_in       , bprod0_in     ,                 &
                           kpblh_in       , pblh_in       , went_in       ,                                                 &
                           qflx_in        , shflx_in      , taux_in       , tauy_in       , aflx_in       ,                 & 
                           landfrac_in    , sgh30_in      ,                                                                 &
                           am_evp_st_in   , evprain_st_in , evpsnow_st_in ,                                                 &
                           cush_inout     , cushavg_inout , cuorg_inout   ,                                                 &
                           awk_PBL_inout                  , delta_thl_PBL_inout           , delta_qt_PBL_inout            , & 
                           delta_u_PBL_inout              , delta_v_PBL_inout             , delta_tr_PBL_inout            , &
                           cu_cmfum_out   , cu_cmfr_inout , cu_thlr_inout , cu_qtr_inout  , cu_ur_inout   , cu_vr_inout   , &
                           cu_qlr_inout   , cu_qir_inout  , cu_trr_inout  ,                                                 &
                           cu_cmfrd_out   , cu_thlrd_out  , cu_qtrd_out   , cu_urd_out    , cu_vrd_out    ,                 &
                           cu_qlrd_out    , cu_qird_out   , cu_trrd_out   ,                                                 &
                           am_u_out       , qlm_u_out     , qim_u_out     ,                                                 &
                           am_d_out       , qlm_d_out     , qim_d_out     ,                                                 &
                           cmf_u_out      , slflx_out     , qtflx_out     ,                                                 & 
                           qvten_out      , qlten_out     , qiten_out     , trten_out     ,                                 &
                           sten_out       , uten_out      , vten_out      ,                                                 &
                           qrten_out      , qsten_out     ,                                                                 & 
                           rqc_l_out      , rqc_i_out     , rqc_out       , rnc_l_out     , rnc_i_out     ,                 &
                           rliq_out       , precip_out    , snow_out      , evapc_out     ,                                 &
                           cnt_out        , cnb_out       , cmf_det_out   , ql_det_out    , qi_det_out    ,                 &
                           lchnk )

   !---------------------------------------------------!
   !                                                   ! 
   !     The Unified Convection Scheme - UNICON        !
   !     Developed by Sungsu Park                      !
   !     For detailed description, See JAS             !
   !     Copyright is owned by Sungsu Park             !
   !                                                   !
   !---------------------------------------------------!
 
   ! --------------------------- !
   ! Main Input-Output variables !
   ! --------------------------- !

   integer , intent(in)    :: lchnk
   integer , intent(in)    :: mix                                  !  Number of columns
   integer , intent(in)    :: iend                                 !  Number of columns 
   integer , intent(in)    :: mkx                                  !  k = 1 : Lowest layer, k = mkx : Top layer
   integer , intent(in)    :: ncnst                                !  Number of tracers
   real(r8), intent(in)    :: dt                                   !  Time step in seconds: 2*delta_t [s]

   real(r8), intent(in)    :: ps0_in(mix,0:mkx)                    !  Environmental pressure at the interface [Pa]
   real(r8), intent(in)    :: zs0_in(mix,0:mkx)                    !  Environmental height   at the interface [m]
   real(r8), intent(in)    :: p0_in(mix,mkx)                       !  Environmental pressure at the mid-point [m]
   real(r8), intent(in)    :: z0_in(mix,mkx)                       !  Environmental height   at the mid-point [m]
   real(r8), intent(in)    :: dp0_in(mix,mkx)                      !  Environmental layer pressure thickness  [Pa] > 0
   real(r8), intent(in)    :: dpdry0_in(mix,mkx)                   !  Environmental layer dry pressure thickness  [Pa] > 0
   real(r8), intent(in)    :: u0_in(mix,mkx)                       !  Environmental zonal wind [m/s]
   real(r8), intent(in)    :: v0_in(mix,mkx)                       !  Environmental meridional wind [m/s]
   real(r8), intent(in)    :: qv0_in(mix,mkx)                      !  Environmental water  vapor specific humidity [kg/kg]
   real(r8), intent(in)    :: ql0_in(mix,mkx)                      !  Environmental liquid water specific humidity [kg/kg]
   real(r8), intent(in)    :: qi0_in(mix,mkx)                      !  Environmental ice          specific humidity [kg/kg]
   real(r8), intent(in)    :: tr0_in(mix,mkx,ncnst)                !  Environmental tracers [ #/kg, kg/kg ]
   real(r8), intent(in)    :: t0_in(mix,mkx)                       !  Environmental temperature [K]
   real(r8), intent(in)    :: ast0_in(mix,mkx)                     !  Stratiform fractional area at the layer mid-point [ fraction ]
   real(r8), intent(in)    :: tke0_in(mix,0:mkx)                   !  TKE at the interface [ m2/s2 ]
   real(r8), intent(in)    :: bprod0_in(mix,0:mkx)                 !  Buoyancy production at the interface [ m2/s3 ]

   integer(i4), intent(in) :: kpblh_in(mix)                        !  Layer index with PBL top in it or at the base interface
   real(r8), intent(in)    :: pblh_in(mix)                         !  PBL top height [ m ]
   real(r8), intent(in)    :: went_in(mix)                         !  Entrainment rate at the PBL top interface directly from the UW PBL scheme [ m/s ]
   real(r8), intent(in)    :: qflx_in(mix)                         !  Upward water vapor flux into atmosphere at surface [ kg/m2/s ]
   real(r8), intent(in)    :: shflx_in(mix)                        !  Upward sensible heat flux into atmosphere at surface [ J/m2/s ]
   real(r8), intent(in)    :: taux_in(mix)                         !  Upward zonal      wind stress into atmosphere at surface [ kg m/s /m2/s ] 
   real(r8), intent(in)    :: tauy_in(mix)                         !  Upward meridional wind stress into atmosphere at surface [ kg m/s /m2/s ] 
   real(r8), intent(in)    :: aflx_in(mix,ncnst)                   !  Upward tracer fluxes          into atmosphere at surface [ #/m2/s, kg/m2/s ]

   real(r8), intent(in)    :: landfrac_in(mix)                     !  Land  Fraction [ fraction ]    

   real(r8), intent(in)    :: sgh30_in(mix)                        !  Standard deviation of subgrid topographic height at 30 s horizontal area [ meter ] 
                                                                   !  This 'sgh30' ( not sgh ) is used for the parameterization of tms. 
   ! Aug.08.2013. Evaporation of stratiform precipitation
   real(r8), intent(in)    :: am_evp_st_in(mix,mkx)                !  Evaporation area of stratiform precipitation [fraction]
   real(r8), intent(in)    :: evprain_st_in(mix,mkx)               !  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
   real(r8), intent(in)    :: evpsnow_st_in(mix,mkx)               !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.

   real(r8), intent(inout) :: cush_inout(mix)                      !  Cumulus top height [ m ]
   real(r8), intent(inout) :: cushavg_inout(mix)                   !  Mean cumulus top height weighted by updraft masss flux at surface [ m ]
   real(r8), intent(inout) :: cuorg_inout(mix)                     !  Covective organization parameter [ 0-1 ]
   real(r8), intent(inout) :: awk_PBL_inout(mix)                   !  Wake area within PBL [ 0 - 1 ]
   real(r8), intent(inout) :: delta_thl_PBL_inout(mix)             !  Difference of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
   real(r8), intent(inout) :: delta_qt_PBL_inout(mix)              !  Difference of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
   real(r8), intent(inout) :: delta_u_PBL_inout(mix)               !  Difference of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(r8), intent(inout) :: delta_v_PBL_inout(mix)               !  Difference of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(r8), intent(inout) :: delta_tr_PBL_inout(mix,ncnst)        !  Difference of tr  between off-wake region and grid-mean value averaged over the PBL [ kg/kg, #/kg ]

   real(r8), intent(out)   :: cu_cmfum_out(mix,mkx)                !  The amount of mass involved in the updraft buoyancy sorting at the previous time step [ kg/s/m2 ]
   real(r8), intent(inout) :: cu_cmfr_inout(mix,mkx)               !  The amount of detrained mass from convective updraft and downdraft at the previous time step [ kg/s/m2 ]
   real(r8), intent(inout) :: cu_thlr_inout(mix,mkx)               !  Mass-flux weighted mean 'thl' of detrained mass from convective updraft and downdraft at the previous time step [ K ]
   real(r8), intent(inout) :: cu_qtr_inout(mix,mkx)                !  Mass-flux weighted mean 'qt'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(r8), intent(inout) :: cu_ur_inout(mix,mkx)                 !  Mass-flux weighted mean 'u'   of detrained mass from convective updraft and downdraft at the previous time step [ m/s ]
   real(r8), intent(inout) :: cu_vr_inout(mix,mkx)                 !  Mass-flux weighted mean 'v'   of detrained mass from convective updraft and downdraft at the previous time step [ m/s ]
   real(r8), intent(inout) :: cu_qlr_inout(mix,mkx)                !  Mass-flux weighted mean 'ql'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(r8), intent(inout) :: cu_qir_inout(mix,mkx)                !  Mass-flux weighted mean 'qi'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(r8), intent(inout) :: cu_trr_inout(mix,mkx,ncnst)          !  Mass-flux weighted mean 'tr'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]

   real(r8)                :: cu_thvr_inout(mix,mkx)               !  Mass-flux weighted mean 'thv' of detrained mass from convective updraft and downdraft at the previous time step [ K ]
   real(r8)                :: cu_rhr_inout(mix,mkx)                !  Mass-flux weighted mean 'rh'  of detrained mass from convective updraft and downdraft at the previous time step [ ratio ]

   real(r8), intent(out)   :: cu_cmfrd_out(mix,mkx)                !  The amount of detrained mass from convective downdraft at the previous time step [ kg/s/m2 ]
   real(r8), intent(out)   :: cu_thlrd_out(mix,mkx)                !  Mass-flux weighted mean 'thl' of detrained mass from convective downdraft at the previous time step [ K ]
   real(r8), intent(out)   :: cu_qtrd_out(mix,mkx)                 !  Mass-flux weighted mean 'qt'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   real(r8), intent(out)   :: cu_urd_out(mix,mkx)                  !  Mass-flux weighted mean 'u'   of detrained mass from convective downdraft at the previous time step [ m/s ]
   real(r8), intent(out)   :: cu_vrd_out(mix,mkx)                  !  Mass-flux weighted mean 'v'   of detrained mass from convective downdraft at the previous time step [ m/s ]
   real(r8), intent(out)   :: cu_qlrd_out(mix,mkx)                 !  Mass-flux weighted mean 'ql'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   real(r8), intent(out)   :: cu_qird_out(mix,mkx)                 !  Mass-flux weighted mean 'qi'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   real(r8), intent(out)   :: cu_trrd_out(mix,mkx,ncnst)           !  Mass-flux weighted mean 'tr'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]

   ! Formal output variables

   real(r8), intent(out)   :: am_u_out(mix,mkx)                    !  Updraft fractional area [ fraction ] 
   real(r8), intent(out)   :: qlm_u_out(mix,mkx)                   !  Area-weighted in-cloud LWC within updraft fractional area [ kg / kg ]
   real(r8), intent(out)   :: qim_u_out(mix,mkx)                   !  Area-weighted in-cloud IWC within updraft fractional area [ kg / kg ]

   real(r8), intent(out)   :: am_d_out(mix,mkx)                    !  Downdraft fractional area [ fraction ] 
   real(r8), intent(out)   :: qlm_d_out(mix,mkx)                   !  Area-weighted in-cloud LWC within downdraft fractional area [ kg / kg ]
   real(r8), intent(out)   :: qim_d_out(mix,mkx)                   !  Area-weighted in-cloud IWC within downdraft fractional area [ kg / kg ]

   real(r8), intent(out)   :: cmf_u_out(mix,0:mkx)                 !  Upward convective mass flux at the interface [ kg / s / m2 ]
   real(r8), intent(out)   :: slflx_out(mix,0:mkx)                 !  Net upward convective flux of liquid static energy [ J / s / m2 ]
   real(r8), intent(out)   :: qtflx_out(mix,0:mkx)                 !  Net upward convective flux of total specific humidity [ kg / s / m2 ]
 
   real(r8), intent(out)   :: qvten_out(mix,mkx)                   !  Tendency of water vapor specific humidity [ kg / kg / s ]
   real(r8), intent(out)   :: qlten_out(mix,mkx)                   !  Tendency of liquid water mixing ratio [ kg / kg / s ]
   real(r8), intent(out)   :: qiten_out(mix,mkx)                   !  Tendency of ice mixing ratio [ kg / kg / s ]
   real(r8), intent(out)   :: sten_out(mix,mkx)                    !  Tendency of dry static energy [ J / kg / s ]
   real(r8), intent(out)   :: uten_out(mix,mkx)                    !  Tendency of zonal wind [ m / s / s ]
   real(r8), intent(out)   :: vten_out(mix,mkx)                    !  Tendency of meridional wind [ m / s / s ]
   real(r8), intent(out)   :: trten_out(mix,mkx,ncnst)             !  Tendency of tracers [ # / kg / s, kg / kg / s ]

   real(r8), intent(out)   :: qrten_out(mix,mkx)                   !  Production rate of rain by lateral expels of cumulus condensate [kg/kg/s]
   real(r8), intent(out)   :: qsten_out(mix,mkx)                   !  Production rate of snow by lateral expels of cumulus condensate [kg/kg/s]
   real(r8), intent(out)   :: precip_out(mix)                      !  Precipitation flux at surface in flux unit [ m / s ]
   real(r8), intent(out)   :: snow_out(mix)                        !  Snow flux at surface in flux unit [ m / s ]
   real(r8), intent(out)   :: evapc_out(mix,mkx)                   !  Evaporation rate of convective precipitation within environment [ kg/kg/s ]

   real(r8), intent(out)   :: rqc_out(mix,mkx)                     ! Production rate of raw detrained LWC+IWC  [kg/kg/s] > 0
   real(r8), intent(out)   :: rqc_l_out(mix,mkx)                   ! Production rate of raw detrained LWC      [kg/kg/s] > 0
   real(r8), intent(out)   :: rqc_i_out(mix,mkx)                   ! Production rate of raw detrained IWC      [kg/kg/s] > 0
   real(r8), intent(out)   :: rliq_out(mix)                        ! Vertical integral of 'rqc_out' in flux unit [m/s]
   real(r8), intent(out)   :: rnc_l_out(mix,mkx)                   ! Production rate of raw detrained droplet number of cloud liquid droplets [#/kg/s] > 0
   real(r8), intent(out)   :: rnc_i_out(mix,mkx)                   ! Production rate of raw detrained droplet number of cloud    ice droplets [#/kg/s] > 0

   real(r8), intent(out)   :: cnt_out(mix)                         ! Cloud top  interface index ( ki = kpen )
   real(r8), intent(out)   :: cnb_out(mix)                         ! Cloud base interface index ( ki = krel-1 )

   real(r8), intent(out)   :: cmf_det_out(mix,mkx)                 ! Detrained mass flux only from convective updraft (not from environmental air) and downdraft [ kg / s / m2 ] 
   real(r8), intent(out)   :: ql_det_out(mix,mkx)                  ! Detrained LWC without mixing with the environment ( flux-convergence & subsidence-detrainment consistent ) [ kg / kg ]
   real(r8), intent(out)   :: qi_det_out(mix,mkx)                  ! Detrained LWC without mixing with the environment ( flux-convergence & subsidence-detrainment consistent ) [ kg / kg ]

   ! ------------------------- ! 
   ! Internal output variables !
   ! ------------------------- !

   real(r8)                :: cmf_out(mix,0:mkx)                   !  Net upward convective mass flux at the interface [kg/s/m2]
   real(r8)                :: uflx_out(mix,0:mkx)                  !  Net upward convective flux of zonal momentum [m/s/s/m2]
   real(r8)                :: vflx_out(mix,0:mkx)                  !  Net upward convective flux of meridional momentum [m/s/s/m2]

   real(r8)                :: slflx_u_out(mix,0:mkx)               !  Upward convective flux of liquid static energy [J/s/m2]
   real(r8)                :: qtflx_u_out(mix,0:mkx)               !  Upward convective flux of total specific humidity [kg/s/m2]
   real(r8)                :: uflx_u_out(mix,0:mkx)                !  Upward convective flux of zonal momentum [kg m/s/s/m2]
   real(r8)                :: vflx_u_out(mix,0:mkx)                !  Upward convective flux of meridional momentum [kg m/s/s/m2]

   real(r8)                :: cmf_d_out(mix,0:mkx)                 !  Downward convective mass flux at the interface [kg/s/m2]
   real(r8)                :: slflx_d_out(mix,0:mkx)               !  Downward convective flux of liquid static energy [J/s/m2]
   real(r8)                :: qtflx_d_out(mix,0:mkx)               !  Downward convective flux of total specific humidity [kg/s/m2]
   real(r8)                :: uflx_d_out(mix,0:mkx)                !  Downward convective flux of zonal momentum [kg m/s/s/m2]
   real(r8)                :: vflx_d_out(mix,0:mkx)                !  Downward convective flux of meridional momentum [kg m/s/s/m2]

   real(r8)                :: thl_orgforce_out(mix)                !  Total organization forcing generating thl difference between 'non-wake' and 'grid-mean' areas [ K / s ] 
   real(r8)                :: qt_orgforce_out(mix)                 !  Total organization forcing generating qt  difference between 'non-wake' and 'grid-mean' areas [ kg / kg / s ]
   real(r8)                :: u_orgforce_out(mix)                  !  Total organization forcing generating u   difference between 'non-wake' and 'grid-mean' areas [ m / s / s ]
   real(r8)                :: v_orgforce_out(mix)                  !  Total organization forcing generating v   difference between 'non-wake' and 'grid-mean' areas [ m / s / s ]
   real(r8)                :: tr_orgforce_out(mix,ncnst)           !  Total organization forcing generating thv difference between 'non-wake' and 'grid-mean' areas [ kg / kg / s or # / kg / s ]
   real(r8)                :: awk_orgforce_out(mix)                !  Total organization forcing generating 'wake area' ( a_wk ) [ 1 / s ] 

   ! Below block is for detailed diagnostic output

   real(r8)                :: flxrain_out(mix,0:mkx)                  
   real(r8)                :: flxsnow_out(mix,0:mkx)                  

   real(r8)                :: thl_orgforce_flx_out(mix)            !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)                :: qt_orgforce_flx_out(mix)             !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)                :: u_orgforce_flx_out(mix)              !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)                :: v_orgforce_flx_out(mix)              !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)                :: awk_orgforce_flx_out(mix)            !  PBL top flux-related forcing for wake area [ 1 / s ]

   real(r8)                :: thl_orgforce_und_out(mix)            !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)                :: qt_orgforce_und_out(mix)             !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)                :: u_orgforce_und_out(mix)              !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)                :: v_orgforce_und_out(mix)              !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)                :: awk_orgforce_mix_out(mix)            !  Lateral-Mixing       forcing for wake area [ 1 / s ]

   real(r8)                :: thl_orgforce_env_out(mix)            !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)                :: qt_orgforce_env_out(mix)             !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)                :: u_orgforce_env_out(mix)              !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)                :: v_orgforce_env_out(mix)              !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)                :: cmf_d_org_pblh_out(mix)              !  Organization-inducing downdraft mass flux at the PBL top interface [ kg / m^2 / s ] 

   ! Above block is for detailed diagnostic output
       
   real(r8)                :: taui_thl_out(mix)                    !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' thl [ 1 / s ]
   real(r8)                :: taui_qt_out(mix)                     !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' qt [ 1 / s ]
   real(r8)                :: taui_u_out(mix)                      !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' u [ 1 / s ]
   real(r8)                :: taui_v_out(mix)                      !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' v [ 1 / s ]
   real(r8)                :: taui_tr_out(mix,ncnst)               !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' tracers [ 1 / s ]
   real(r8)                :: taui_awk_out(mix)                    !  Inverse of damping time scale of the wake area [ 1 / s ]

   real(r8)                :: del_org_out(mix)                     !  Detrainment rate of the cold pool [ 1 / s ]
   real(r8)                :: del0_org_out(mix)                    !  Effective detrainment rate of the cold pool [ 1 / s ]

   real(r8)                :: slten_u_out(mix,mkx)      
   real(r8)                :: qtten_u_out(mix,mkx)      
   real(r8)                :: uten_u_out(mix,mkx)       
   real(r8)                :: vten_u_out(mix,mkx)       
   real(r8)                :: sten_u_out(mix,mkx)       
   real(r8)                :: qvten_u_out(mix,mkx)      
   real(r8)                :: qlten_u_out(mix,mkx)      
   real(r8)                :: qiten_u_out(mix,mkx)      
   real(r8)                :: trten_u_out(mix,mkx,ncnst)      

   real(r8)                :: slten_d_out(mix,mkx)      
   real(r8)                :: qtten_d_out(mix,mkx)      
   real(r8)                :: uten_d_out(mix,mkx)       
   real(r8)                :: vten_d_out(mix,mkx)       
   real(r8)                :: sten_d_out(mix,mkx)       
   real(r8)                :: qvten_d_out(mix,mkx)      
   real(r8)                :: qlten_d_out(mix,mkx)      
   real(r8)                :: qiten_d_out(mix,mkx)      
   real(r8)                :: trten_d_out(mix,mkx,ncnst)      

   real(r8)                :: slten_evp_out(mix,mkx)    
   real(r8)                :: qtten_evp_out(mix,mkx)    
   real(r8)                :: uten_evp_out(mix,mkx)     
   real(r8)                :: vten_evp_out(mix,mkx)     
   real(r8)                :: sten_evp_out(mix,mkx)     
   real(r8)                :: qvten_evp_out(mix,mkx)    
   real(r8)                :: qlten_evp_out(mix,mkx)    
   real(r8)                :: qiten_evp_out(mix,mkx)    
   real(r8)                :: trten_evp_out(mix,mkx,ncnst)      

   real(r8)                :: slten_dis_out(mix,mkx)    
   real(r8)                :: qtten_dis_out(mix,mkx)    
   real(r8)                :: uten_dis_out(mix,mkx)     
   real(r8)                :: vten_dis_out(mix,mkx)     
   real(r8)                :: sten_dis_out(mix,mkx)     
   real(r8)                :: qvten_dis_out(mix,mkx)    
   real(r8)                :: qlten_dis_out(mix,mkx)    
   real(r8)                :: qiten_dis_out(mix,mkx)    
   real(r8)                :: trten_dis_out(mix,mkx,ncnst)      

   real(r8)                :: qlten_sub_out(mix,mkx)    
   real(r8)                :: qiten_sub_out(mix,mkx)    

   real(r8)                :: qlten_det_out(mix,mkx)    
   real(r8)                :: qiten_det_out(mix,mkx)    

   real(r8)                :: thl_u_out(mix,0:mkx)                 !  Mass-flux weighted updraft thl [ K ]
   real(r8)                :: qt_u_out(mix,0:mkx)                  !  Mass-flux weighted updraft qt [ kg / kg ]
   real(r8)                :: u_u_out(mix,0:mkx)                   !  Mass-flux weighted updraft u [ m / s ]
   real(r8)                :: v_u_out(mix,0:mkx)                   !  Mass-flux weighted updraft v [ m / s ]
   real(r8)                :: w_u_out(mix,0:mkx)                   !  Mass-flux weighted updraft w [ m / s ]
   real(r8)                :: ql_u_out(mix,0:mkx)                  !  Mass-flux weighted updraft in-cumulus ql [ kg / kg ]
   real(r8)                :: qi_u_out(mix,0:mkx)                  !  Mass-flux weighted updraft in-cumulus qi [ kg / kg ]
   real(r8)                :: tr_u_out(mix,0:mkx,ncnst)            !  Mass-flux weighted updraft tr [ # / kg, kg / kg ]
   real(r8)                :: wa_u_out(mix,0:mkx)                  !  Area weighted updraft w [ m / s ]
   real(r8)                :: qla_u_out(mix,0:mkx)                 !  Area weighted updraft in-cumulus ql [ kg / kg ]
   real(r8)                :: qia_u_out(mix,0:mkx)                 !  Area weighted updraft in-cumulus qi [ kg / kg ]
   real(r8)                :: a_u_out(mix,0:mkx)                   !  Updraft fractional area [ fraction ]
   real(r8)                :: rad_u_out(mix,0:mkx)                 !  Number weighted effective radius of updraft plumes [ m ]
   real(r8)                :: num_u_out(mix,0:mkx)                 !  Number concentration of updraft plumes [ # / m^2 ]
   real(r8)                :: gamw_u_out(mix,0:mkx)                !  Ratio of 'w_u_out / wa_u_out' [ no ]
   real(r8)                :: thva_u_out(mix,0:mkx)                !  Area weighted updraft thv [ K ]

   real(r8)                :: a_p_out(mix,0:mkx)                   !  Convective precipitation area [ fraction ]
   real(r8)                :: am_evp_out(mix,mkx)                  !  Convective evaporation area [ fraction ]
   real(r8)                :: am_pu_out(mix,mkx)                   !  Overlapping area between convective precipitation and saturated updraft [ fraction ]
   real(r8)                :: x_p_out(mix,0:mkx)                   !  Zonal displacement of the precipitation area from the surface [ m ]
   real(r8)                :: y_p_out(mix,0:mkx)                   !  Meridional displacement of the precipitation area from the surface [ m ]
   real(r8)                :: x_um_out(mix,mkx)                    !  Zonal displacement of the updraft area from the surface [ m ]
   real(r8)                :: y_um_out(mix,mkx)                    !  Meridional displacement of the updraft area from the surface [ m ]

   real(r8)                :: thl_d_out(mix,0:mkx)                 !  Mass-flux weighted downdraft thl [ K ]
   real(r8)                :: qt_d_out(mix,0:mkx)                  !  Mass-flux weighted downdraft qt [ kg / kg ]
   real(r8)                :: u_d_out(mix,0:mkx)                   !  Mass-flux weighted downdraft u [ m / s ]
   real(r8)                :: v_d_out(mix,0:mkx)                   !  Mass-flux weighted downdraft v [ m / s ]
   real(r8)                :: w_d_out(mix,0:mkx)                   !  Mass-flux weighted downdraft w [ m / s ]
   real(r8)                :: ql_d_out(mix,0:mkx)                  !  Mass-flux weighted downdraft in-cumulus ql [ kg / kg ]
   real(r8)                :: qi_d_out(mix,0:mkx)                  !  Mass-flux weighted downdraft in-cumulus qi [ kg / kg ]
   real(r8)                :: tr_d_out(mix,0:mkx,ncnst)            !  Mass-flux weighted downdraft tr [ # / kg, kg / kg ]
   real(r8)                :: wa_d_out(mix,0:mkx)                  !  Area weighted downdraft w [ m / s ]
   real(r8)                :: qla_d_out(mix,0:mkx)                 !  Area weighted downdraft in-cumulus ql [ kg / kg ]
   real(r8)                :: qia_d_out(mix,0:mkx)                 !  Area weighted downdraft in-cumulus qi [ kg / kg ]
   real(r8)                :: a_d_out(mix,0:mkx)                   !  Downdraft fractional area [ fraction ]

   real(r8)                :: thl_u_msfc_out(mix,0:mkx,nseg,niter)       !  Updraft             thl at the interface for each original updraft segment [ K ].
   real(r8)                :: qt_u_msfc_out(mix,0:mkx,nseg,niter)        !  Updraft              qt at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: u_u_msfc_out(mix,0:mkx,nseg,niter)         !  Updraft               u at the interface for each original updraft segment [ m / s ].
   real(r8)                :: v_u_msfc_out(mix,0:mkx,nseg,niter)         !  Updraft               v at the interface for each original updraft segment [ m / s ].
   real(r8)                :: w_u_msfc_out(mix,0:mkx,nseg,niter)         !  Updraft               w at the interface for each original updraft segment [ m / s ].
   real(r8)                :: ql_u_msfc_out(mix,0:mkx,nseg,niter)        !  Updraft              ql at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: qi_u_msfc_out(mix,0:mkx,nseg,niter)        !  Updraft              qi at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: tr_u_msfc_out(mix,0:mkx,nseg,ncnst,niter)  !  Updraft              tr at the interface for each original updraft segment [ # / kg, kg / kg ].
   real(r8)                :: cmf_u_msfc_out(mix,0:mkx,nseg,niter)       !  Updraft             cmf at the interface for each original updraft segment [ kg / s / m^2 ].
   real(r8)                :: a_u_msfc_out(mix,0:mkx,nseg,niter)         !  Updraft               a at the interface for each original updraft segment [ fraction ].
   real(r8)                :: num_u_msfc_out(mix,0:mkx,nseg,niter)       !  Updraft             num at the interface for each original updraft segment [ # / m^2 ].
   real(r8)                :: rad_u_msfc_out(mix,0:mkx,nseg,niter)       !  Updraft             rad at the interface for each original updraft segment [ m ].

   real(r8)                :: eps0_u_msfc_out(mix,0:mkx,nseg,niter)      !  Updraft            eps0 at the interface for each original updraft segment [ 1 / Pa ].
   real(r8)                :: eps_u_msfc_out(mix,0:mkx,nseg,niter)       !  Updraft             eps at the interface for each original updraft segment [ 1 / Pa ].
   real(r8)                :: del_u_msfc_out(mix,0:mkx,nseg,niter)       !  Updraft             del at the interface for each original updraft segment [ 1 / Pa ].
   real(r8)                :: eeps_u_msfc_out(mix,0:mkx,nseg,niter)      !  Updraft            eeps at the interface for each original updraft segment [ no ].
   real(r8)                :: ddel_u_msfc_out(mix,0:mkx,nseg,niter)      !  Updraft            ddel at the interface for each original updraft segment [ no ].
   real(r8)                :: xc_u_msfc_out(mix,0:mkx,nseg,niter)        !  Updraft              xc at the interface for each original updraft segment [ no ].
   real(r8)                :: xs_u_msfc_out(mix,0:mkx,nseg,niter)        !  Updraft              xs at the interface for each original updraft segment [ no ].
   real(r8)                :: xemin_u_msfc_out(mix,0:mkx,nseg,niter)     !  Updraft           xemin at the interface for each original updraft segment [ no ].
   real(r8)                :: xemax_u_msfc_out(mix,0:mkx,nseg,niter)     !  Updraft           xemax at the interface for each original updraft segment [ no ].
   real(r8)                :: cridis_u_msfc_out(mix,0:mkx,nseg,niter)    !  Updraft          cridis at the interface for each original updraft segment [ m ].
   real(r8)                :: thvcuenv_u_msfc_out(mix,0:mkx,nseg,niter)  !  Updraft        thvcuenv at the interface for each original updraft segment [ K ].
   real(r8)                :: thvegenv_u_msfc_out(mix,0:mkx,nseg,niter)  !  Updraft        thvegenv at the interface for each original updraft segment [ K ].
   real(r8)                :: thvxsenv_u_msfc_out(mix,0:mkx,nseg,niter)  !  Updraft        thvxsenv at the interface for each original updraft segment [ K ].
   real(r8)                :: fmix_u_msfc_out(mix,0:mkx,nseg,niter)      !  Updraft            fmix at the interface for each original updraft segment [ no ].
   real(r8)                :: cmfumix_u_msfc_out(mix,0:mkx,nseg,niter)   !  Updraft         cmfumix at the interface for each original updraft segment [ kg / s / m^2 ].

   real(r8)                :: thl_d_msfc_out(mix,0:mkx,nseg,niter)       !  Mass-flux weighted  mean downdraft      thl at the interface for each original updraft segment [ K ].
   real(r8)                :: qt_d_msfc_out(mix,0:mkx,nseg,niter)        !  Mass-flux weighted  mean downdraft       qt at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: u_d_msfc_out(mix,0:mkx,nseg,niter)         !  Mass-flux weighted  mean downdraft        u at the interface for each original updraft segment [ m / s ].
   real(r8)                :: v_d_msfc_out(mix,0:mkx,nseg,niter)         !  Mass-flux weighted  mean downdraft        v at the interface for each original updraft segment [ m / s ].
   real(r8)                :: w_d_msfc_out(mix,0:mkx,nseg,niter)         !  Mass-flux weighted  mean downdraft        w at the interface for each original updraft segment [ m / s ].
   real(r8)                :: ql_d_msfc_out(mix,0:mkx,nseg,niter)        !  Mass-flux weighted  mean downdraft       ql at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: qi_d_msfc_out(mix,0:mkx,nseg,niter)        !  Mass-flux weighted  mean downdraft       qi at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: tr_d_msfc_out(mix,0:mkx,nseg,ncnst,niter)  !  Mass-flux weighted  mean downdraft       tr at the interface for each original updraft segment [ # / kg, kg / kg ].
   real(r8)                :: wa_d_msfc_out(mix,0:mkx,nseg,niter)        !  Area-weighted       mean downdraft        w at the interface for each original updraft segment [ m / s ].
   real(r8)                :: qla_d_msfc_out(mix,0:mkx,nseg,niter)       !  Area-weighted       mean downdraft       ql at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: qia_d_msfc_out(mix,0:mkx,nseg,niter)       !  Area-weighted       mean downdraft       qi at the interface for each original updraft segment [ kg / kg ].
   real(r8)                :: cmf_d_msfc_out(mix,0:mkx,nseg,niter)       !  Net                      downdraft      cmf at the interface for each original updraft segment [ kg / s / m^2 ].
   real(r8)                :: a_d_msfc_out(mix,0:mkx,nseg,niter)         !  Net                      downdraft        a at the interface for each original updraft segment [ fraction ].

   real(r8)                :: ptop_msfc_out(mix,nseg,niter)              !  Updraft top height  of individual original updraft segment defined at surface [ Pa ]
   real(r8)                :: ztop_msfc_out(mix,nseg,niter)              !  Updraft top height  of individual original updraft segment defined at surface [ m ]

   real(r8)                :: thv_b_out(mix,0:mkx) 
   real(r8)                :: thv_t_out(mix,0:mkx) 
   real(r8)                :: thv_mt_out(mix,0:mkx) 
   real(r8)                :: thv_min_out(mix,0:mkx) 

   ! ----------------------------------------------------------- ! 
   ! One-dimensional local variables defined at each grid column !
   ! k = 0   : Surface interface                                 !
   ! k = 1   : Lowest layer                                      !
   ! k = mkx : Top layer                                         !
   ! ----------------------------------------------------------- !

   ! --------------- !
   ! Input variables !
   ! --------------- !

   real(r8)    ps0(0:mkx)                       !  Environmental pressure at the interface [Pa]
   real(r8)    zs0(0:mkx)                       !  Environmental height   at the interface [m]
   real(r8)    p0(mkx)                          !  Environmental pressure at the mid-point [Pa]
   real(r8)    z0(mkx)                          !  Environmental height   at the mid-point [m]
   real(r8)    dp0(mkx)                         !  Environmental layer pressure thickness  [Pa] ( > 0 )
   real(r8)    dpdry0(mkx)                      !  Environmental layer dry pressure thickness  [Pa] ( > 0 )

   real(r8)    u0(mkx)                          !  Environmental zonal wind [m/s]
   real(r8)    v0(mkx)                          !  Environmental meridional wind [m/s]
   real(r8)    qv0(mkx)                         !  Environmental water vapor specific humidity [kg/kg]
   real(r8)    ql0(mkx)                         !  Environmental liquid water mixing ratio [kg/kg]
   real(r8)    qi0(mkx)                         !  Environmental ice mixing ratio [kg/kg]
   real(r8)    tr0(mkx,ncnst)                   !  Environmental tracers [ #/kg, kg/kg ]
   real(r8)    t0(mkx)                          !  Environmental temperature [K]

   real(r8)    ast0(mkx)                        !  Stratiform fractional area at the layer mid-point [ fraction ]
   real(r8)    tke0(0:mkx)                      !  TKE [ m2/s2 ]
   real(r8)    bprod0(0:mkx)                    !  Buoyancy production [ m2/s3 ]

   ! -------------------------------------------------------- !
   ! Environmental variables derived from the input variables !
   ! -------------------------------------------------------- !

   real(r8)    dptr0(mkx,ncnst)                 !  Environmental layer pressure thickness for each dry or moist tracers [Pa] ( > 0 )
   real(r8)    dz0(mkx)                         !  Environmental layer thickness  [m] ( > 0 )
   real(r8)    dps0(0:mkx)                      !  Environmental interfacial layer pressure thickness  [Pa] ( > 0 )

   real(r8)    qt0(mkx)                         !  Environmental total specific humidity [kg/kg]
   real(r8)    thl0(mkx)                        !  Environmental liquid potential temperature [K]
   real(r8)    thv0(mkx)                        !  Environmental virtual potential temperature [K]
   real(r8)    rh0(mkx)                         !  Environmental mean rh [fraction]
   real(r8)    ssqt0(mkx)                       !  Vertical gradient of qt0 [kg/kg/Pa]
   real(r8)    ssthl0(mkx)                      !  Vertical gradient of thl0 [K/Pa]
   real(r8)    ssu0(mkx)                        !  Vertical gradient of u0 [m/s/Pa]
   real(r8)    ssv0(mkx)                        !  Vertical gradient of v0 [m/s/Pa]
   real(r8)    rho0(mkx)                        !  Environmental density [kg/m3]

   real(r8)    thl0bot(mkx)                     !  Environmental thl at the bottom interface [K]
   real(r8)    thl0top(mkx)                     !  Environmental thl at the top interface [K]
   real(r8)    qt0bot(mkx)                      !  Environmental qt at the bottom interface [kg/kg]
   real(r8)    qt0top(mkx)                      !  Environmental qt at the top interface [kg/kg]
   real(r8)    u0bot(mkx)                       !  Environmental u at the bottom interface [m/s]
   real(r8)    u0top(mkx)                       !  Environmental u at the top interface [m/s]
   real(r8)    v0bot(mkx)                       !  Environmental v at the bottom interface [m/s]
   real(r8)    v0top(mkx)                       !  Environmental v at the top interface [m/s]
   real(r8)    thv0bot(mkx)                     !  Environmental virtual potential temperature at the bottom interface [K]
   real(r8)    thv0top(mkx)                     !  Environmental virtual potential temperature at the top interface [K]
   real(r8)    thvl0bot(mkx)                    !  Environmental thvl at the bottom interface [K]
   real(r8)    thvl0top(mkx)                    !  Environmental thvl at the top interface [K]
   real(r8)    ql0bot(mkx)                      !  Environmental ql at the bottom interface [kg/kg]
   real(r8)    ql0top(mkx)                      !  Environmental ql at the top interface [kg/kg]
   real(r8)    qi0bot(mkx)                      !  Environmental qi at the bottom interface [kg/kg]
   real(r8)    qi0top(mkx)                      !  Environmental qi at the top interface [kg/kg]
   real(r8)    tr0bot(mkx,ncnst)                !  Environmental tracer at the bottom interface [#/kg, kg/kg]
   real(r8)    tr0top(mkx,ncnst)                !  Environmental tracer at the top    interface [#/kg, kg/kg]
   real(r8)    rho0bot(mkx)                     !  Environmental density at the bottom interface [kg/m3]
   real(r8)    rho0top(mkx)                     !  Environmental density at the top interface [kg/m3]
   real(r8)    rh0bot(mkx)                      !  Environmental RH at the bottom interface [0-1]
   real(r8)    exn0(mkx)                        !  Exner function at the mid-points
   real(r8)    exns0(0:mkx)                     !  Exner function at the interfaces
   real(r8)    ssql0(mkx)                       !  Vertical gradient of ql0   [kg/kg/Pa]
   real(r8)    ssqi0(mkx)                       !  Vertical gradient of qi0   [kg/kg/Pa]
   real(r8)    sstr0(mkx,ncnst)                 !  Vertical gradient of environmental tracers [ #/kg/Pa, kg/kg/Pa ]

   ! ----------------- !
   ! Cumulus variables !
   ! ----------------- !

   real(r8)    flxrain(0:mkx)                   !  Grid-mean convective rain flux after evaporation within downdraft and environment [kg/m2/s]
   real(r8)    flxsnow(0:mkx)                   !  Grid-mean convective snow flux after evaporation within downdraft and environment [kg/m2/s]
   real(r8)    flxtrrs(0:mkx,ncnst)             !  Grid-mean convective tracer flux after evaporation within downdraft and environment [kg(#)/m2/s]

   real(r8)    flxrain_msfc(0:mkx,nseg)         !  Grid-mean convective rain flux after evaporation within downdraft and environment for each original updraft segment [kg/m2/s]
   real(r8)    flxsnow_msfc(0:mkx,nseg)         !  Grid-mean convective snow flux after evaporation within downdraft and environment for each original updraft segment [kg/m2/s]
   real(r8)    flxtrrs_msfc(0:mkx,nseg,ncnst)   !  Grid-mean convective tracer flux after evaporation within downdraft and environment for each original updraft segment [kg(#)/m2/s]

   real(r8)    cmf_u_mix(mkx)                   !  Total amount of updraft mass flux involved in the buoyancy sorting [kg/m2/s]
   real(r8)    cmf_r(mkx)                       !  Total amount of detrained mass into the environment mass [kg/m2/s]
   real(r8)    thl_r(mkx)                       !  Mass flux weighted conservative scalar of detrained airs  
   real(r8)    qt_r(mkx)                        !  Same as above
   real(r8)    u_r(mkx)                         !  Same as above
   real(r8)    v_r(mkx)                         !  Same as above
   real(r8)    ql_r(mkx)                        !  Same as above
   real(r8)    qi_r(mkx)                        !  Same as above
   real(r8)    tr_r(mkx,ncnst)                  !  Same as above

   ! ------------------------------------------------------------------------------------------------------------------------ !
   ! Below '2' variables are same as the above, but only consider the detrained component purely from the convective updraft, !
   ! not from the environmental airs involved in the mixing. This approach is very important because this approach is         !
   ! fully consistently connecting the 'flux-convergence formula' to the 'subsidence-detrainment form' particularly           !
   ! for the budget of cloud condensate. This approach will be directly used for simulating the effect of convective          !
   ! detrainment on the critical relative humidity in the stratiform macrophysics. In addition, this approach may be          !
   ! used for defining the mixing environmental air, instead of the above approach.                                           !
   ! ------------------------------------------------------------------------------------------------------------------------ !  

   real(r8)    cmf_r2(mkx)                     
   real(r8)    thl_r2(mkx)                     
   real(r8)    qt_r2(mkx)                      
   real(r8)    u_r2(mkx)                       
   real(r8)    v_r2(mkx)                       
   real(r8)    ql_r2(mkx)                      
   real(r8)    qi_r2(mkx)                      
   real(r8)    tr_r2(mkx,ncnst)                

   real(r8)    cmf_u(0:mkx)                     !  Total updraft mass flux at the model interface [ kg / m2 / s ] 
   real(r8)    w_u(0:mkx)                       !  Mass-flux weighted updraft vertical velocity [ m / s ]
   real(r8)    wa_u(0:mkx)                      !  Area weighted updraft vertical velocity [ m / s ]
   real(r8)    a_u(0:mkx)                       !  Physical updraft fractional area [ fraction ]
   real(r8)    num_u(0:mkx)                     !  Number density of updraft plumes [ # / m^2 ]
   real(r8)    rad_u(0:mkx)                     !  Physical mean effective radius of updraft plumes [ m ] 
   real(r8)    thl_u(0:mkx)                     !  Mass-flux weighted updraft liquid potential temperature [ K ]
   real(r8)    qt_u(0:mkx)                      !  Mass-flux weighted updraft total specific humidity [ kg / kg ]
   real(r8)    u_u(0:mkx)                       !  Mass-flux weighted updraft zonal velocity [ m / s ]
   real(r8)    v_u(0:mkx)                       !  Mass-flux weighted updraft meridional velocity [ m / s ]
   real(r8)    ql_u(0:mkx)                      !  Mass-flux weighted in-cloud LWC     within convective   updraft [ kg / kg ]
   real(r8)    qi_u(0:mkx)                      !  Mass-flux weighted in-cloud IWC     within convective   updraft [ kg / kg ]
   real(r8)    qla_u(0:mkx)                     !  Area weighted in-cloud LWC within convective   updraft [ kg / kg ]
   real(r8)    qia_u(0:mkx)                     !  Area weighted in-cloud IWC within convective   updraft [ kg / kg ] 
   real(r8)    tr_u(0:mkx,ncnst)                !  Mass-flux weighted in-cloud tracers within convective   updraft [ # / kg, kg / kg ] 
   real(r8)    thva_u(0:mkx)                    !  Area weighted thv within updraft [ K ]

   real(r8)    cmf_u_dia(mkx)                   !  Total updraft mass flux at individual cloud tops [ kg / m2 / s ] 
   real(r8)    evp_thll_u(mkx)                  !  Mass-flux weighted diabatic change of updraft 'thl' at each cloud top due to evaporation of rain [ K ] <= 0.
   real(r8)    evp_qtl_u(mkx)                   !  Mass-flux weighted diabatic change of updraft 'qt'  at each cloud top due to evaporation of rain [ kg / kg ] >= 0.
   real(r8)    evp_thli_u(mkx)                  !  Mass-flux weighted diabatic change of updraft 'thl' at each cloud top due to evaporation of snow [ K ] <= 0.
   real(r8)    evp_qti_u(mkx)                   !  Mass-flux weighted diabatic change of updraft 'qt'  at each cloud top due to evaporation of snow [ kg / kg ] >= 0.
   real(r8)    evp_tr_u(mkx,ncnst)              !  Mass-flux weighted diabatic change of updraft  tracer  at each cloud top due to evaporation of rain + snow [ # / kg, kg / kg ]
   real(r8)    prep_thll_u(mkx)                 !  Mass-flux weighted diabatic change of updraft 'thl' at each cloud top due to production  of rain [ K ] >= 0.
   real(r8)    prep_qtl_u(mkx)                  !  Mass-flux weighted diabatic change of updraft 'qt'  at each cloud top due to production  of rain [ kg / kg ] <= 0.
   real(r8)    prep_thli_u(mkx)                 !  Mass-flux weighted diabatic change of updraft 'thl' at each cloud top due to production  of snow [ K ] >= 0.
   real(r8)    prep_qti_u(mkx)                  !  Mass-flux weighted diabatic change of updraft 'qt'  at each cloud top due to production  of snow [ kg / kg ] <= 0.
   real(r8)    prep_tr_u(mkx,ncnst)             !  Mass-flux weighted diabatic change of updraft  tracer  at each cloud top due to production  of rain + snow [ # / kg, kg / kg ]
   real(r8)    eff_ql_u(mkx)                    !  Mass-flux weighted diabatic change of updraft 'ql'  at each cloud top due to effective diabatic forcing on cloud condensate [ kg / kg ]
   real(r8)    eff_qi_u(mkx)                    !  Mass-flux weighted diabatic change of updraft 'qi'  at each cloud top due to effective diabatic forcing on cloud condensate [ kg / kg ]
   real(r8)    eff_tr_u(mkx,ncnst)              !  Mass-flux weighted diabatic change of updraft  tracer  at each cloud top due to effective diabatic forcing on the tracer [ # / kg, kg / kg ]
   real(r8)    PGF_u_u(mkx)                     !  Mass-flux weighted diabatic change of updraft 'u'   at each cloud top due to horizontal PGF forcing [ m / s ]
   real(r8)    PGF_v_u(mkx)                     !  Mass-flux weighted diabatic change of updraft 'v'   at each cloud top due to horizontal PGF forcing [ m / s ]

   real(r8)    f_srcd(mkx)                      !  Total source of downdraft generated from the updraft. f_srcd(k) = f_dd(k) + f_dud(k) + f_nud(k). [ ratio ] >= 0. 
   real(r8)    thl_srcd(mkx)                    !  Mass-flux weighted thl of net sources of downdraft [ K ]
   real(r8)    qt_srcd(mkx)                     !  Mass-flux weighted qt  of net sources of downdraft [ kg / kg ]
   real(r8)    u_srcd(mkx)                      !  Mass-flux weighted u   of net sources of downdraft [ m / s ]
   real(r8)    v_srcd(mkx)                      !  Mass-flux weighted v   of net sources of downdraft [ m / s ]
   real(r8)    tr_srcd(mkx,ncnst)               !  Mass-flux weighted tracer  of net sources of downdraft [ # / kg, kg / kg ]
   real(r8)    ql_srcd(mkx)                     !  Mass-flux weighted ql  of net sources of downdraft [ kg / kg ]
   real(r8)    qi_srcd(mkx)                     !  Mass-flux weighted qi  of net sources of downdraft [ kg / kg ]

   real(r8)    f_srcds(mkx,nseg,3)              !  Total source of downdraft generated from the updraft. f_srcd(k) = f_dd(k) + f_dud(k) + f_nud(k). [ ratio ] >= 0. 
   real(r8)    thl_srcds(mkx,nseg,3)            !  The thl of net sources of downdraft [ K ]
   real(r8)    qt_srcds(mkx,nseg,3)             !  The qt  of net sources of downdraft [ kg / kg ]
   real(r8)    u_srcds(mkx,nseg,3)              !  The u   of net sources of downdraft [ m / s ]
   real(r8)    v_srcds(mkx,nseg,3)              !  The v   of net sources of downdraft [ m / s ]
   real(r8)    tr_srcds(mkx,nseg,3,ncnst)       !  The tracer of net sources of downdraft [ # / kg, kg / kg ]
   real(r8)    ql_srcds(mkx,nseg,3)             !  The ql  of net sources of downdraft [ kg / kg ]
   real(r8)    qi_srcds(mkx,nseg,3)             !  The qi  of net sources of downdraft [ kg / kg ]

   real(r8)    f_srcr(mkx)                      !  Total source of remained airs generated from the updraft. f_srcr(k) = f_dr(k) + f_dur(k) + f_nur(k). [ ratio ] >= 0. 
   real(r8)    thl_srcr(mkx)                    !  Mass-flux weighted thl of net sources of remained airs [ K ]
   real(r8)    qt_srcr(mkx)                     !  Mass-flux weighted qt  of net sources of remained airs [ kg / kg ]
   real(r8)    u_srcr(mkx)                      !  Mass-flux weighted u   of net sources of remained airs [ m / s ]
   real(r8)    v_srcr(mkx)                      !  Mass-flux weighted v   of net sources of remained airs [ m / s ]
   real(r8)    tr_srcr(mkx,ncnst)               !  Mass-flux weighted tracer  of net sources of remained airs [ # / kg, kg / kg ]
   real(r8)    ql_srcr(mkx)                     !  Mass-flux weighted ql  of net sources of remained airs [ kg / kg ]
   real(r8)    qi_srcr(mkx)                     !  Mass-flux weighted qi  of net sources of remained airs [ kg / kg ]

   real(r8)    f_srcr2(mkx)                     
   real(r8)    thl_srcr2(mkx)                   
   real(r8)    qt_srcr2(mkx)                    
   real(r8)    u_srcr2(mkx)                     
   real(r8)    v_srcr2(mkx)                     
   real(r8)    tr_srcr2(mkx,ncnst)              
   real(r8)    ql_srcr2(mkx)                    
   real(r8)    qi_srcr2(mkx)                    

   real(r8)    f_srcrs(mkx,nseg,3)              !  Total source of remained airs generated from the updraft. f_srcr(k) = f_dr(k) + f_dur(k) + f_nur(k). [ ratio ] >= 0. 
   real(r8)    thl_srcrs(mkx,nseg,3)            !  The thl of net sources of remained airs [ K ]
   real(r8)    qt_srcrs(mkx,nseg,3)             !  The qt  of net sources of remained airs [ kg / kg ]
   real(r8)    u_srcrs(mkx,nseg,3)              !  The u   of net sources of remained airs [ m / s ]
   real(r8)    v_srcrs(mkx,nseg,3)              !  The v   of net sources of remained airs [ m / s ]
   real(r8)    tr_srcrs(mkx,nseg,3,ncnst)       !  The tracer  of net sources of remained airs [ # / kg , kg / kg ]
   real(r8)    ql_srcrs(mkx,nseg,3)             !  The ql  of net sources of remained airs [ kg / kg ]
   real(r8)    qi_srcrs(mkx,nseg,3)             !  The qi  of net sources of remained airs [ kg / kg ]

   real(r8)    f_srcrs2(mkx,nseg,3)            
   real(r8)    thl_srcrs2(mkx,nseg,3)          
   real(r8)    qt_srcrs2(mkx,nseg,3)           
   real(r8)    u_srcrs2(mkx,nseg,3)            
   real(r8)    v_srcrs2(mkx,nseg,3)            
   real(r8)    tr_srcrs2(mkx,nseg,3,ncnst)     
   real(r8)    ql_srcrs2(mkx,nseg,3)           
   real(r8)    qi_srcrs2(mkx,nseg,3)           

   real(r8)    cmf_ru(mkx)                      !  f_srcr(mkx) * cmf_u(km) [ kg / m^2 / s ]
   real(r8)    thl_ru(mkx)                      !  = thl_srcr(mkx)
   real(r8)    qt_ru(mkx)                       !  = qt_srcr(mkx)
   real(r8)    u_ru(mkx)                        !  = u_srcr(mkx)
   real(r8)    v_ru(mkx)                        !  = v_srcr(mkx)
   real(r8)    ql_ru(mkx)                       !  = ql_srcr(mkx)
   real(r8)    qi_ru(mkx)                       !  = qi_srcr(mkx)
   real(r8)    tr_ru(mkx,ncnst)                 !  = tr_srcr(mkx,ncnst)

   real(r8)    cmf_ru2(mkx)                     
   real(r8)    thl_ru2(mkx)                     
   real(r8)    qt_ru2(mkx)                      
   real(r8)    u_ru2(mkx)                       
   real(r8)    v_ru2(mkx)                       
   real(r8)    ql_ru2(mkx)                      
   real(r8)    qi_ru2(mkx)                      
   real(r8)    tr_ru2(mkx,ncnst)                

   real(r8)    cmf_ad(0:mkx,mkx,nseg,3)         !  Downdraft mass flux at the model interface     originated from the downdraft sources of the 'k' layer [ kg / m^2 / s ] >= 0.
   real(r8)    w_ad(0:mkx,mkx,nseg,3)           !  Mass-flux weighted downdraft vertical velocity originated from the downdraft sources of the 'k' layer [ m / s ]
   real(r8)    a_ad(0:mkx,mkx,nseg,3)           !  Physical downdraft fractional area             originated from the downdraft sources of the 'k' layer [ fraction ]
   real(r8)    thl_ad(0:mkx,mkx,nseg,3)         !  Mass-flux weighted downdraft thl               originated from the downdraft sources of the 'k' layer [ K ]
   real(r8)    qt_ad(0:mkx,mkx,nseg,3)          !  Mass-flux weighted downdraft qt                originated from the downdraft sources of the 'k' layer [ kg / kg ]
   real(r8)    u_ad(0:mkx,mkx,nseg,3)           !  Mass-flux weighted downdraft u                 originated from the downdraft sources of the 'k' layer [ m / s ]
   real(r8)    v_ad(0:mkx,mkx,nseg,3)           !  Mass-flux weighted downdraft v                 originated from the downdraft sources of the 'k' layer [ m / s ]
   real(r8)    ql_ad(0:mkx,mkx,nseg,3)          !  Mass-flux weighted downdraft ql                originated from the downdraft sources of the 'k' layer [ kg / kg ]
   real(r8)    qi_ad(0:mkx,mkx,nseg,3)          !  Mass-flux weighted downdraft qi                originated from the downdraft sources of the 'k' layer [ kg / kg ]
   real(r8)    tr_ad(0:mkx,mkx,nseg,3,ncnst)    !  Mass-flux weighted downdraft tracer            originated from the downdraft sources of the 'k' layer [ # / kg, kg / kg ]

   real(r8)    dpad(mkx,mkx,nseg,3)              

   real(r8)    cmf_d(0:mkx)                     !  Total downdraft mass flux at the model interface [ kg / m^2 / s ] >= 0.
   real(r8)    w_d(0:mkx)                       !  Mass-flux weighted downdraft vertical velocity   [ m / s ]
   real(r8)    wa_d(0:mkx)                      !  Area weighted downdraft vertical velocity        [ m / s ]
   real(r8)    a_d(0:mkx)                       !  Physical downdraft fractional area               [ fraction ]
   real(r8)    thl_d(0:mkx)                     !  Mass-flux weighted downdraft thl                 [ K ]
   real(r8)    qt_d(0:mkx)                      !  Mass-flux weighted downdraft qt                  [ kg / kg ]
   real(r8)    u_d(0:mkx)                       !  Mass-flux weighted downdraft u                   [ m / s ]
   real(r8)    v_d(0:mkx)                       !  Mass-flux weighted downdraft v                   [ m / s ]
   real(r8)    ql_d(0:mkx)                      !  Mass-flux weighted in-cloud LWC     within convective downdraft [ kg / kg ]
   real(r8)    qi_d(0:mkx)                      !  Mass-flux weighted in-cloud IWC     within convective downdraft [ kg / kg ] 
   real(r8)    qla_d(0:mkx)                     !  Area weighted in-cloud LWC within convective downdraft [ kg / kg ]
   real(r8)    qia_d(0:mkx)                     !  Area weighted in-cloud IWC within convective downdraft [ kg / kg ] 
   real(r8)    tr_d(0:mkx,ncnst)                !  Mass-flux weighted downdraft tracer              [ # / kg, kg / kg ]

   real(r8)    cmf_ar(mkx,mkx,nseg,3)           !  The mass flux of detrained airs into environment from the downdraft originated from the downdraft sources of the 'k' layer [ kg / m^2 / s ] >= 0.
   real(r8)    thl_ar(mkx,mkx,nseg,3)           !  Mass-flux weighted thl of detrained airs         from the downdraft originated from the downdraft sources of the 'k' layer [ K ]
   real(r8)    qt_ar(mkx,mkx,nseg,3)            !  Mass-flux weighted qt  of detrained airs         from the downdraft originated from the downdraft sources of the 'k' layer [ kg / kg ]
   real(r8)    u_ar(mkx,mkx,nseg,3)             !  Mass-flux weighted u   of detrained airs         from the downdraft originated from the downdraft sources of the 'k' layer [ m / s ]
   real(r8)    v_ar(mkx,mkx,nseg,3)             !  Mass-flux weighted v   of detrained airs         from the downdraft originated from the downdraft sources of the 'k' layer [ m / s ]
   real(r8)    tr_ar(mkx,mkx,nseg,3,ncnst)      !  Mass-flux weighted tracer  of detrained airs     from the downdraft originated from the downdraft sources of the 'k' layer [ # / kg, kg / kg ]
   real(r8)    ql_ar(mkx,mkx,nseg,3)            !  Mass-flux weighted ql  of detrained airs         from the downdraft originated from the downdraft sources of the 'k' layer [ kg / kg ]
   real(r8)    qi_ar(mkx,mkx,nseg,3)            !  Mass-flux weighted qi  of detrained airs         from the downdraft originated from the downdraft sources of the 'k' layer [ kg / kg ]

   real(r8)    cmf_rd(mkx)                      !  The mass flux of detrained airs into environment from the downdraft [ kg / m^2 / s ] >= 0.
   real(r8)    thl_rd(mkx)                      !  Mass-flux weighted thl     of detrained airs     from the downdraft [ K ]
   real(r8)    qt_rd(mkx)                       !  Mass-flux weighted qt      of detrained airs     from the downdraft [ kg / kg ]
   real(r8)    u_rd(mkx)                        !  Mass-flux weighted u       of detrained airs     from the downdraft [ m / s ]
   real(r8)    v_rd(mkx)                        !  Mass-flux weighted v       of detrained airs     from the downdraft [ m / s ]
   real(r8)    ql_rd(mkx)                       !  Mass-flux weighted ql      of detrained airs     from the downdraft [ kg / kg ]
   real(r8)    qi_rd(mkx)                       !  Mass-flux weighted qi      of detrained airs     from the downdraft [ kg / kg ]
   real(r8)    tr_rd(mkx,ncnst)                 !  Mass-flux weighted tracer  of detrained airs     from the downdraft [ # / kg, kg / kg ]

   real(r8)    cmf_ad_dia(mkx,mkx,nseg,3)       !  Downdraft mass flux at each downdraft base in each layer originated from the downdraft sources of the 'k' layer [ kg / m2 / s ] >= 0. 
   real(r8)    evp_thll_ad(mkx,mkx,nseg,3)      !  Diabatic change of downdraft 'thl'   at each downdraft base due to evaporation of rain [ K ] <= 0.
   real(r8)    evp_qtl_ad(mkx,mkx,nseg,3)       !  Diabatic change of downdraft 'qt'    at each downdraft base due to evaporation of rain [ kg / kg ] >= 0.
   real(r8)    evp_thli_ad(mkx,mkx,nseg,3)      !  Diabatic change of downdraft 'thl'   at each downdraft base due to evaporation of snow [ K ] <= 0.
   real(r8)    evp_qti_ad(mkx,mkx,nseg,3)       !  Diabatic change of downdraft 'qt'    at each downdraft base due to evaporation of snow [ kg / kg ] >= 0.
   real(r8)    evp_tr_ad(mkx,mkx,nseg,3,ncnst)  !  Diabatic change of downdraft  tracer at each downdraft base due to evaporation of rain + snow [ kg / kg ]
   real(r8)    prep_thll_ad(mkx,mkx,nseg,3)     !  Diabatic change of downdraft 'thl'   at each downdraft base due to production  of rain [ K ] >= 0.
   real(r8)    prep_qtl_ad(mkx,mkx,nseg,3)      !  Diabatic change of downdraft 'qt'    at each downdraft base due to production  of rain [ kg / kg ] <= 0.
   real(r8)    prep_thli_ad(mkx,mkx,nseg,3)     !  Diabatic change of downdraft 'thl'   at each downdraft base due to production  of snow [ K ] >= 0.
   real(r8)    prep_qti_ad(mkx,mkx,nseg,3)      !  Diabatic change of downdraft 'qt'    at each downdraft base due to production  of snow [ kg / kg ] <= 0.
   real(r8)    prep_tr_ad(mkx,mkx,nseg,3,ncnst) !  Diabatic change of downdraft  tracer at each downdraft base due to production  of rain + snow [ # / kg, kg / kg ]
   real(r8)    eff_ql_ad(mkx,mkx,nseg,3)        !  Diabatic change of downdraft 'ql'    at each downdraft base due to effective diabatic forcing on cloud condensate [ kg / kg ]
   real(r8)    eff_qi_ad(mkx,mkx,nseg,3)        !  Diabatic change of downdraft 'qi'    at each downdraft base due to effective diabatic forcing on cloud condensate [ kg / kg ]
   real(r8)    eff_tr_ad(mkx,mkx,nseg,3,ncnst)  !  Diabatic change of downdraft  tracer at each downdraft base due to effective diabatic forcing on tracer [ # / kg, kg / kg ]
   real(r8)    PGF_u_ad(mkx,mkx,nseg,3)         !  Diabatic change of downdraft 'u'     at each downdraft base due to horizontal PGF forcing [ m / s ]
   real(r8)    PGF_v_ad(mkx,mkx,nseg,3)         !  Diabatic change of downdraft 'v'     at each downdraft base due to horizontal PGF forcing [ m / s ]
   real(r8)    wdep_tr_ad(mkx,mkx,nseg,3,ncnst) !  Diabatic change of downdraft  tracer at each downdraft base due to wet deposition of aerosols within downdraft
                                                !                     including both the cloud-borne and interstitial aerosols within convective downdraft [ # / kg, kg / kg ]

   real(r8)    cmf_d_dia(mkx)                   !  Total downdraft mass flux at the downdraft base in each layer [ kg / m2 / s ] >= 0. 
   real(r8)    evp_thll_d(mkx)                  !  Mass-flux weighted diabatic change of downdraft 'thl'   at the downdraft base due to evaporation of rain [ K ] <= 0.
   real(r8)    evp_qtl_d(mkx)                   !  Mass-flux weighted diabatic change of downdraft 'qt'    at the downdraft base due to evaporation of rain [ kg / kg ] >= 0.
   real(r8)    evp_thli_d(mkx)                  !  Mass-flux weighted diabatic change of downdraft 'thl'   at the downdraft base due to evaporation of snow [ K ] <= 0.
   real(r8)    evp_qti_d(mkx)                   !  Mass-flux weighted diabatic change of downdraft 'qt'    at the downdraft base due to evaporation of snow [ kg / kg ] >= 0.
   real(r8)    evp_tr_d(mkx,ncnst)              !  Mass-flux weighted diabatic change of downdraft  tracer at the downdraft base due to evaporation of rain + snow [ # / kg, kg / kg ]
   real(r8)    prep_thll_d(mkx)                 !  Mass-flux weighted diabatic change of downdraft 'thl'   at the downdraft base due to production  of rain [ K ] >= 0.
   real(r8)    prep_qtl_d(mkx)                  !  Mass-flux weighted diabatic change of downdraft 'qt'    at the downdraft base due to production  of rain [ kg / kg ] <= 0.
   real(r8)    prep_thli_d(mkx)                 !  Mass-flux weighted diabatic change of downdraft 'thl'   at the downdraft base due to production  of snow [ K ] >= 0.
   real(r8)    prep_qti_d(mkx)                  !  Mass-flux weighted diabatic change of downdraft 'qt'    at the downdraft base due to production  of snow [ kg / kg ] <= 0.
   real(r8)    prep_tr_d(mkx,ncnst)             !  Mass-flux weighted diabatic change of downdraft  tracer at the downdraft base due to production  of rain + snow [ kg / kg ]
   real(r8)    eff_ql_d(mkx)                    !  Mass-flux weighted diabatic change of downdraft 'ql'    at the downdraft base due to effective diabatic forcing on cloud condensate [ kg / kg ]
   real(r8)    eff_qi_d(mkx)                    !  Mass-flux weighted diabatic change of downdraft 'qi'    at the downdraft base due to effective diabatic forcing on cloud condensate [ kg / kg ]
   real(r8)    eff_tr_d(mkx,ncnst)              !  Mass-flux weighted diabatic change of downdraft  tracer at the downdraft base due to effective diabatic forcing on tracer [ # / kg, kg / kg ]
   real(r8)    PGF_u_d(mkx)                     !  Mass-flux weighted diabatic change of downdraft 'u'     at the downdraft base due to horizontal PGF forcing [ m / s ]
   real(r8)    PGF_v_d(mkx)                     !  Mass-flux weighted diabatic change of downdraft 'v'     at the downdraft base due to horizontal PGF forcing [ m / s ]

   real(r8)    qlten_sub(mkx)                   !  Environmental tendency of ql due to compensating subsidence / upwelling [ kg / kg / s ]
   real(r8)    qiten_sub(mkx)                   !  Environmental tendency of qi due to compensating subsidence / upwelling [ kg / kg / s ] 

   real(r8)    rqc_l(mkx)                       !  Environmental tendency of ql      due to raw detrainment of updraft ( sum of 'cmf_dur(k)' and 'cmf_nur(k)' ) [ kg / kg / s ]
   real(r8)    rqc_i(mkx)                       !  Environmental tendency of qi      due to raw detrainment of updraft ( sum of 'cmf_dur(k)' and 'cmf_nur(k)' ) [ kg / kg / s ]
   real(r8)    rqc(mkx)                         !  Environmental tendency of ql + qi due to raw detrainment of updraft ( sum of 'cmf_dur(k)' and 'cmf_nur(k)' ) [ kg / kg / s ]
   real(r8)    rnc_l(mkx)                       !  Environmental tendency of nl      due to raw detrainment of updraft ( sum of 'cmf_dur(k)' and 'cmf_nur(k)' ) [  # / kg / s ]
   real(r8)    rnc_i(mkx)                       !  Environmental tendency of ni      due to raw detrainment of updraft ( sum of 'cmf_dur(k)' and 'cmf_nur(k)' ) [  # / kg / s ]

   real(r8)    qlten_det(mkx)                   !  Environmental tendency of ql      due to all detrainment of updraft and downdraft [ kg / kg / s ]
   real(r8)    qiten_det(mkx)                   !  Environmental tendency of qi      due to all detrainment of updraft and downdraft [ kg / kg / s ]

   real(r8)    am_u_msfc(mkx,nseg)              !  Updraft       fractional area at the layer mid-point for each original updraft segment [ fraction ].
   real(r8)    am_d_msfc(mkx,nseg)              !  Downdraft     fractional area at the layer mid-point for each original updraft segment [ fraction ].

   real(r8)    am_u(mkx)                        !  Updraft       fractional area at the layer mid-point [ fraction ]. 0 <= am_u(mkx) <= au_max <= 1.  
   real(r8)    am_d(mkx)                        !  Downdraft     fractional area at the layer mid-point [ fraction ]. 0 <= am_d(mkx) <= ad_max <= 1.  

   real(r8)    qlm_u_msfc(mkx,nseg)             !  Updraft LWC at the layer mid-point for each original updraft segment [ kg / kg ].
   real(r8)    qim_u_msfc(mkx,nseg)             !  Updraft IWC at the layer mid-point for each original updraft segment [ kg / kg ].
   real(r8)    thlm_u_msfc(mkx,nseg)            !  Updraft thl at the layer mid-point for each original updraft segment [ K ].
   real(r8)    qtm_u_msfc(mkx,nseg)             !  Updraft qt  at the layer mid-point for each original updraft segment [ kg / kg ].
   real(r8)    um_u_msfc(mkx,nseg)              !  Updraft u   at the layer mid-point for each original updraft segment [ m / s ].
   real(r8)    vm_u_msfc(mkx,nseg)              !  Updraft v   at the layer mid-point for each original updraft segment [ m / s ].
   real(r8)    trm_u_msfc(mkx,nseg,ncnst)       !  Updraft tr  at the layer mid-point for each original updraft segment [ kg / kg or # / kg ].
   
   real(r8)    qlm_u(mkx)                       !  Area-weighted updraft LWC at the layer mid-point [ kg / kg ].
   real(r8)    qim_u(mkx)                       !  Area-weighted updraft IWC at the layer mid-point [ kg / kg ].
   real(r8)    thlm_u(mkx)                      !  Area-weighted updraft thl at the layer mid-point [ K ].
   real(r8)    qtm_u(mkx)                       !  Area-weighted updraft qt  at the layer mid-point [ kg / kg ].
   real(r8)    um_u(mkx)                        !  Area-weighted updraft u   at the layer mid-point [ m / s ].
   real(r8)    vm_u(mkx)                        !  Area-weighted updraft v   at the layer mid-point [ m / s ].
   real(r8)    trm_u(mkx,ncnst)                 !  Area-weighted updraft tr  at the layer mid-point [ kg / kg or # / kg ].

   real(r8)    qlm_d_msfc(mkx,nseg)             !  Downdraft LWC at the layer mid-point for each original updraft segment [ kg / kg ].
   real(r8)    qim_d_msfc(mkx,nseg)             !  Downdraft IWC at the layer mid-point for each original updraft segment [ kg / kg ].

   real(r8)    qlm_d(mkx)                       !  Area-weighted downdraft LWC at the layer mid-point [ kg / kg ].
   real(r8)    qim_d(mkx)                       !  Area-weighted downdraft IWC at the layer mid-point [ kg / kg ].

   real(r8)    am_s(mkx)                        !  Stratiform    fractional area at the layer mid-point [ fraction ]. 0 <= am_s(mkx) <= 1.
   real(r8)    am_r(mkx)                        !  Clear-sky     fractional area at the layer mid-point [ fraction ]. 0 <= am_r(mkx)=1-am_u(mkx)-am_s(mkx)<= 1.  
   real(r8)    am_up(mkx)                       !  Precipitating updraft fractional area at the layer mid-point [ fraction ]. 0 <= am_up(mkx) <= am_u(mkx) <= 1.  
   real(r8)    am_us(mkx)                       !  Saturated     updraft fractional area at the layer mid-point [ fraction ]. 0 <= am_us(mkx) <= am_u(mkx) <= 1.  

   real(r8)    am_up_msfc(mkx,nseg)             !  Precipitating updraft fractional area at the layer mid-point for each original updraft segment [ fraction ].
   real(r8)    am_us_msfc(mkx,nseg)             !  Saturated     updraft fractional area at the layer mid-point for each original updraft segment [ fraction ].

   real(r8)    a_p(0:mkx)                       !  Physical convective precipitation area at the interface [ fraction ]
   real(r8)    a_p_msfc(0:mkx,nseg)             !  Physical convective precipitation area at the interface for each original updraft segment [ fraction ]
   real(r8)    a_pu                             !  Overlapping area between convective precipitation area at the top interface and   updraft area at the layer mid-point [ fraction ]
   real(r8)    a_pd                             !  Overlapping area between convective precipitation area at the top interface and downdraft area at the layer mid-point [ fraction ]
   real(r8)    a_ps                             !  Overlapping area between convective precipitation area at the top interface and   stratus area at the layer mid-point [ fraction ]
   real(r8)    a_pr                             !  Overlapping area between convective precipitation area at the top interface and     clear area at the layer mid-point [ fraction ]
   real(r8)    a_evp                            !  Fractional area where evaporation of precipitation occurs. 
   real(r8)    a_ovp                            ! 
   real(r8)    am_evp_msfc(mkx,nseg)            !  Same as 'a_evp' but into array variables.
   real(r8)    am_pu_msfc(mkx,nseg)
   real(r8)    am_pd_msfc(mkx,nseg)
   real(r8)    am_pr_msfc(mkx,nseg)
   real(r8)    am_ps_msfc(mkx,nseg)
   real(r8)    am_evp(mkx)         
   real(r8)    am_pu(mkx)
   real(r8)    am_pd(mkx)
   real(r8)    am_pr(mkx)
   real(r8)    am_ps(mkx)

   real(r8)    cmf_det(mkx)                     ! Detrained mass flux only from convective updraft (not from environmental air) and downdraft [ kg / s / m2 ] 
   real(r8)    ql_det(mkx)                      ! Detrained LWC without mixing with the environment ( flux-convergence & subsidence-detrainment consistent ) [ kg / kg ]
   real(r8)    qi_det(mkx)                      ! Detrained LWC without mixing with the environment ( flux-convergence & subsidence-detrainment consistent ) [ kg / kg ]

   real(r8)    am_evp_nw                        !  Overlapping area between 'evaporation area' and 'non-wake area' [ 0 - 1 ]
   real(r8)    am_p_nw                          !  Overlapping area between 'precipitation area' and 'non-wake area' [ 0 - 1 ]

   ! Aug.08.2013. Add stratiform part
   real(r8)    am_evp_nw_st                     !  Overlapping area between 'stratiform evaporation area' and 'non-wake area' [ 0 - 1 ]
   real(r8)    tmp2_st
   real(r8)    am_evp_st(mkx)
   real(r8)    evprain_st(mkx)
   real(r8)    evpsnow_st(mkx) 

   ! ---------------------------------------------------------------------- !
   ! Individual Updraft Segment Variables at 'msfc' index at each interface !
   ! ---------------------------------------------------------------------- !

   real(r8)    thl_u_msfc(0:mkx,nseg)           !  Updraft             thl at the interface for each original updraft segment [ K ].
   real(r8)    qt_u_msfc(0:mkx,nseg)            !  Updraft              qt at the interface for each original updraft segment [ kg / kg ].
   real(r8)    u_u_msfc(0:mkx,nseg)             !  Updraft               u at the interface for each original updraft segment [ m / s ].
   real(r8)    v_u_msfc(0:mkx,nseg)             !  Updraft               v at the interface for each original updraft segment [ m / s ].
   real(r8)    w_u_msfc(0:mkx,nseg)             !  Updraft               w at the interface for each original updraft segment [ m / s ].
   real(r8)    ql_u_msfc(0:mkx,nseg)            !  Updraft              ql at the interface for each original updraft segment [ kg / kg ].
   real(r8)    qi_u_msfc(0:mkx,nseg)            !  Updraft              qi at the interface for each original updraft segment [ kg / kg ].
   real(r8)    tr_u_msfc(0:mkx,nseg,ncnst)      !  Updraft              tr at the interface for each original updraft segment [ # / kg, kg / kg ].
   real(r8)    cmf_u_msfc(0:mkx,nseg)           !  Updraft             cmf at the interface for each original updraft segment [ kg / s / m^2 ].
   real(r8)    a_u_msfc(0:mkx,nseg)             !  Updraft               a at the interface for each original updraft segment [ fraction ].
   real(r8)    num_u_msfc(0:mkx,nseg)           !  Updraft             num at the interface for each original updraft segment [ # / m^2 ].
   real(r8)    rad_u_msfc(0:mkx,nseg)           !  Updraft             rad at the interface for each original updraft segment [ m ].

   real(r8)    eps0_u_msfc(0:mkx,nseg)          !  Updraft            eps0 at the interface for each original updraft segment [ 1 / Pa ].
   real(r8)    eps_u_msfc(0:mkx,nseg)           !  Updraft             eps at the interface for each original updraft segment [ 1 / Pa ].
   real(r8)    del_u_msfc(0:mkx,nseg)           !  Updraft             del at the interface for each original updraft segment [ 1 / Pa ].
   real(r8)    eeps_u_msfc(0:mkx,nseg)          !  Updraft            eeps at the interface for each original updraft segment [ no ].
   real(r8)    ddel_u_msfc(0:mkx,nseg)          !  Updraft            ddel at the interface for each original updraft segment [ no ].
   real(r8)    xc_u_msfc(0:mkx,nseg)            !  Updraft              xc at the interface for each original updraft segment [ no ].
   real(r8)    xs_u_msfc(0:mkx,nseg)            !  Updraft              xs at the interface for each original updraft segment [ no ].
   real(r8)    xemin_u_msfc(0:mkx,nseg)         !  Updraft           xemin at the interface for each original updraft segment [ no ].
   real(r8)    xemax_u_msfc(0:mkx,nseg)         !  Updraft           xemax at the interface for each original updraft segment [ no ].
   real(r8)    cridis_u_msfc(0:mkx,nseg)        !  Updraft          cridis at the interface for each original updraft segment [ m ].
   real(r8)    thvcuenv_u_msfc(0:mkx,nseg)      !  Updraft        thvcuenv at the interface for each original updraft segment [ K ].
   real(r8)    thvegenv_u_msfc(0:mkx,nseg)      !  Updraft        thvegenv at the interface for each original updraft segment [ K ].
   real(r8)    thvxsenv_u_msfc(0:mkx,nseg)      !  Updraft        thvxsenv at the interface for each original updraft segment [ K ].
   real(r8)    fmix_u_msfc(0:mkx,nseg)          !  Updraft            fmix at the interface for each original updraft segment [ no ].
   real(r8)    cmfumix_u_msfc(0:mkx,nseg)       !  Updraft         cmfumix at the interface for each original updraft segment [ kg / s / m^2 ]

   ! ------------------------------------------ !
   ! Variables associated with vertical overlap !
   ! ------------------------------------------ !

   real(r8)    x_um_msfc(mkx,nseg)              !  Location ( x ) of convective updraft at the layer mid-point relative to the convective updraft at surface [ m ]
   real(r8)    y_um_msfc(mkx,nseg)              !  Location ( y ) of convective updraft at the layer mid-point relative to the convective updraft at surface [ m ]
   real(r8)    x_p_msfc(0:mkx,nseg)             !  Location ( x ) of convective precipitation area at the interface relative to the convective updraft at surface [ m ]
   real(r8)    y_p_msfc(0:mkx,nseg)             !  Location ( y ) of convective precipitation area at the interface relative to the convective updraft at surface [ m ]

   ! -------------------------------------------------------------------------------------------------------- !
   ! Individual Mean Downdraft variables for each Updraft Segment Variables at 'msfc' index at each interface !
   ! -------------------------------------------------------------------------------------------------------- !

   real(r8)    thl_d_msfc(0:mkx,nseg)           !  Mass-flux weighted  mean downdraft      thl at the interface for each original updraft segment [ K ].
   real(r8)    qt_d_msfc(0:mkx,nseg)            !  Mass-flux weighted  mean downdraft       qt at the interface for each original updraft segment [ kg / kg ].
   real(r8)    u_d_msfc(0:mkx,nseg)             !  Mass-flux weighted  mean downdraft        u at the interface for each original updraft segment [ m / s ].
   real(r8)    v_d_msfc(0:mkx,nseg)             !  Mass-flux weighted  mean downdraft        v at the interface for each original updraft segment [ m / s ].
   real(r8)    w_d_msfc(0:mkx,nseg)             !  Mass-flux weighted  mean downdraft        w at the interface for each original updraft segment [ m / s ].
   real(r8)    ql_d_msfc(0:mkx,nseg)            !  Mass-flux weighted  mean downdraft       ql at the interface for each original updraft segment [ kg / kg ].
   real(r8)    qi_d_msfc(0:mkx,nseg)            !  Mass-flux weighted  mean downdraft       qi at the interface for each original updraft segment [ kg / kg ].
   real(r8)    tr_d_msfc(0:mkx,nseg,ncnst)      !  Mass-flux weighted  mean downdraft       tr at the interface for each original updraft segment [ # / kg, kg / kg ].
   real(r8)    wa_d_msfc(0:mkx,nseg)            !  Area-weighted       mean downdraft        w at the interface for each original updraft segment [ m / s ].
   real(r8)    qla_d_msfc(0:mkx,nseg)           !  Area-weighted       mean downdraft       ql at the interface for each original updraft segment [ kg / kg ].
   real(r8)    qia_d_msfc(0:mkx,nseg)           !  Area-weighted       mean downdraft       qi at the interface for each original updraft segment [ kg / kg ].
   real(r8)    cmf_d_msfc(0:mkx,nseg)           !  Net                      downdraft      cmf at the interface for each original updraft segment [ kg / s / m^2 ].
   real(r8)    a_d_msfc(0:mkx,nseg)             !  Net                      downdraft        a at the interface for each original updraft segment [ fraction ].

   ! ------------------------- !
   ! Updraft Segment Variables !
   ! ------------------------- !

   real(r8)    ytop(nseg)                       !  If 1 ( 0 ), updraft segment does ( not ) reach to the top interface [ no unit ]        
   real(r8)    xc(nseg)                         !  Critical   mixing fraction for buoyancy sorting [ fraction ]  0 <= xc <= 1.
   real(r8)    xs(nseg)                         !  Saturation mixing fraction for buoyancy sorting [ fraction ]  0 <= xs <= 1.
   real(r8)    eeps(nseg)                       !  Non-dimensional fractional entrainment rate [ no unit ]
   real(r8)    ddel(nseg)                       !  Non-dimensional fractional detrainment rate [ no unit ]
   real(r8)    eps(nseg)                        !  Fractional entrainment rate [ 1 / Pa ]
   real(r8)    del(nseg)                        !  Fractional detrainment rate [ 1 / Pa ]
   real(r8)    eps0(nseg)                       !  Fractional mixing rate [ 1 / Pa ]
   real(r8)    eps0org(nseg)                    !  Fractional mixing rate [ 1 / Pa ]
   real(r8)    xe_min(nseg)                     !  Minimum mixing fraction for conversion into downdraft [ fraction ]
   real(r8)    xe_max(nseg)                     !  Maximum mixing fraction for conversion into downdraft [ fraction ]
   real(r8)    dpa(nseg)                        !  Vertical distance that updraft segment can rise in each layer [ Pa ]. 0 <= dp <= dp0.
   real(r8)    dza(nseg)                        !  Vertical distance that updraft segment can rise in each layer [ m ].  0 <= dz <= dz0.
   real(r8)    ptop(nseg)                       !  Updraft top height [ Pa ]
   real(r8)    ztop(nseg)                       !  Updraft top height [ m ]
   real(r8)    ptops(mkx,nseg)                  !  Updraft top height for each segment in each layer [ Pa ]
   real(r8)    ztops(mkx,nseg)                  !  Updraft top height for each segment in each layer [ m ]
   integer     m_from_msfc(mkx,nseg)            !  Get 'm' index from 'msfc' index in each layer that updraft reaches [ no unit ]
   integer     msfc_from_m(mkx,nseg)            !  Get 'msfc' index from 'm' index in each layer that updraft reaches [ no unit ]
   integer     ktop_msfc(nseg)                  !  The top layer index of individual original updraft segment defined at surface [ no ]
   real(r8)    ptop_msfc(nseg)                  !  Updraft top height  of individual original updraft segment defined at surface [ Pa ]
   real(r8)    ztop_msfc(nseg)                  !  Updraft top height  of individual original updraft segment defined at surface [ m ]
   real(r8)    fmix(nseg)                       !  When multiplied by cmf_au(m), it becomes the amount of updraft mass involved in the buoyancy sorting mixing [ no unit ]
   real(r8)    f_wu(nseg)                       !  For updraft vertical velocity constraint [ no unit ]
   real(r8)    fmixd                            !  Same as 'fmix' but for downdraft.

   real(r8)    alpha(nseg)                      !  Mixing fraction within updraft at the k = 1 interface [ no unit ]
   real(r8)    Pmu(nseg)                        !  Normalized PDF of updraft mass flux            at the k = 1 interface [ 1 / d_alpha ]
   real(r8)    Pau(nseg)                        !  Normalized PDF of updraft fractional area      at the k = 1 interface [ 1 / d_alpha ]
   real(r8)    Pnu(nseg)                        !  Normalized PDF of updraft number concentration at the k = 1 interface [ 1 / d_alpha ]
   real(r8)    rnorm_a                          !  Area normalization constant to force the sum of updraft fraction area at surface to be the specified 'au_base' in all cases.
   real(r8)    rnorm_m                          !  Area normalization constant to force the sum of non-organized updraft mass flux at surface to be the specified 'cmfu_base' in all cases.
   real(r8)    cmfu_base                        !  Physical analytical non-organized updraft mass flux at surface for a given 'au_base'.

   real(r8)    cmf_au(nseg)                     !  Updraft mass flux            at the base interface [ kg / s / m^2 ]
   real(r8)    thl_au(nseg)                     !  Updraft thl                  at the base interface [ K ]
   real(r8)    qt_au(nseg)                      !  Updraft qt                   at the base interface [ kg / kg ] 
   real(r8)    u_au(nseg)                       !  Updraft u                    at the base interface [ m / s ]
   real(r8)    v_au(nseg)                       !  Updraft v                    at the base interface [ m / s ]
   real(r8)    w_au(nseg)                       !  Updraft vertical velocity    at the base interface [ m / s ]
   real(r8)    ql_au(nseg)                      !  Updraft ql                   at the base interface [ kg / kg ]
   real(r8)    qi_au(nseg)                      !  Updraft qi                   at the base interface [ kg / kg ]
   real(r8)    tr_au(nseg,ncnst)                !  Updraft tracers              at the base interface [ # / kg, kg / kg ]
   real(r8)    a_au(nseg)                       !  Updraft fractional area      at the base interface [ fraction ]
   real(r8)    num_au(nseg)                     !  Updraft number concentration at the base interface [ # / m^2 ]
   real(r8)    rad_au(nseg)                     !  Updraft radius               at the base interface [ m ]
   real(r8)    thv_au(nseg)                     !  Updraft thv                  at the base interface [ K ]
   real(r8)    S_b_ql_au(nseg)                  !  Updraft rain production rate at the base interface [ kg / kg / Pa ]
   real(r8)    S_b_qi_au(nseg)                  !  Updraft snow production rate at the base interface [ kg / kg / Pa ]

   real(r8)    cmf_aut(nseg)                    !  Updraft mass flux            at the cloud top or top interface [ kg / s / m^2 ]
   real(r8)    thl_aut(nseg)                    !  Updraft thl                  at the cloud top or top interface [ K ]
   real(r8)    qt_aut(nseg)                     !  Updraft qt                   at the cloud top or top interface [ kg / kg ] 
   real(r8)    u_aut(nseg)                      !  Updraft u                    at the cloud top or top interface [ m / s ]
   real(r8)    v_aut(nseg)                      !  Updraft v                    at the cloud top or top interface [ m / s ]
   real(r8)    w_aut(nseg)                      !  Updraft vertical velocity    at the cloud top or top interface [ m / s ]
   real(r8)    ql_aut(nseg)                     !  Updraft ql                   at the cloud top or top interface [ kg / kg ]
   real(r8)    qi_aut(nseg)                     !  Updraft qi                   at the cloud top or top interface [ kg / kg ]
   real(r8)    a_aut(nseg)                      !  Updraft fractional area      at the cloud top or top interface [ fraction ]
   real(r8)    num_aut(nseg)                    !  Updraft number concentration at the cloud top or top interface [ # / m^2 ]
   real(r8)    rad_aut(nseg)                    !  Updraft radius               at the cloud top or top interface [ m ]
   real(r8)    tr_aut(nseg,ncnst)               !  Updraft tracer               at the cloud top or top interface [ # / kg, kg / kg ] 
   real(r8)    thv_aut(nseg)                    !  Updraft thv                  at the cloud top or top interface [ K ]
   real(r8)    S_t_ql_au(nseg)                  !  Updraft rain production rate at the cloud top or top interface [ kg / kg / Pa ]
   real(r8)    S_t_qi_au(nseg)                  !  Updraft snow production rate at the cloud top or top interface [ kg / kg / Pa ]

   real(r8)    evp_thll_au(nseg)                !  Diabatic change of updraft 'thl' at each cloud top due to evaporation of rain [ K ] <= 0.
   real(r8)    evp_qtl_au(nseg)                 !  Diabatic change of updraft 'qt'  at each cloud top due to evaporation of rain [ kg / kg ] >= 0.
   real(r8)    evp_thli_au(nseg)                !  Diabatic change of updraft 'thl' at each cloud top due to evaporation of snow [ K ] <= 0.
   real(r8)    evp_qti_au(nseg)                 !  Diabatic change of updraft 'qt'  at each cloud top due to evaporation of snow [ kg / kg ] >= 0.
   real(r8)    evp_tr_au(nseg,ncnst)            !  Diabatic change of updraft  tracer  at each cloud top due to evaporation of rain + snow [ # / kg, kg / kg ]
   real(r8)    prep_thll_au(nseg)               !  Diabatic change of updraft 'thl' at each cloud top due to production  of rain [ K ] >= 0.
   real(r8)    prep_qtl_au(nseg)                !  Diabatic change of updraft 'qt'  at each cloud top due to production  of rain [ kg / kg ] <= 0.
   real(r8)    prep_thli_au(nseg)               !  Diabatic change of updraft 'thl' at each cloud top due to production  of snow [ K ] >= 0.
   real(r8)    prep_qti_au(nseg)                !  Diabatic change of updraft 'qt'  at each cloud top due to production  of snow [ kg / kg ] <= 0.
   real(r8)    prep_tr_au(nseg,ncnst)           !  Diabatic change of updraft  tracer  at each cloud top due to production  of rain + snow [ # / kg, kg / kg ]
   real(r8)    eff_ql_au(nseg)                  !  Diabatic change of updraft 'ql'  at each cloud top due to effective diabatic forcing on cloud condensate [ K ] >= 0.
   real(r8)    eff_qi_au(nseg)                  !  Diabatic change of updraft 'qi'  at each cloud top due to effective diabatic forcing on cloud condensate [ kg / kg ] 
   real(r8)    eff_tr_au(nseg,ncnst)            !  Diabatic change of updraft  tracer at each cloud top due to effective diabatic forcing on tracer [ # / kg, kg / kg ] 
   real(r8)    PGF_u_au(nseg)                   !  Diabatic change of updraft 'u'   at each cloud top due to horizontal PGF forcing [ m / s ]
   real(r8)    PGF_v_au(nseg)                   !  Diabatic change of updraft 'v'   at each cloud top due to horizontal PGF forcing [ m / s ]

   ! ---------------------------------------------------------------------------------------------------------------------------------------------- !
   ! Variables to compute effective diabatic forcings ( condensation-evaporation-freezing ) on cloud condensate and tracers and PGF on uv momentum. !
   ! ---------------------------------------------------------------------------------------------------------------------------------------------- !
    
   real(r8)    ql_aut_adi
   real(r8)    qi_aut_adi
   real(r8)    eff_ql
   real(r8)    eff_qi
   real(r8)    eff_tr(ncnst)
   real(r8)    u_aut_adi
   real(r8)    v_aut_adi

   ! ----------------------------------------------------------------------------- !
   ! Variables associated with the treatment of entrainment dilution of downdraft. !
   ! ----------------------------------------------------------------------------- !

   real(r8)    eps_dn, del_dn
   real(r8)    ql_db_adi, qi_db_adi, qv_db_adi, qv_db_adi_evp
   real(r8)    u_db_adi ,  v_db_adi
   real(r8)    ql_dt, qi_dt, qv_dt, th_dt 
   real(r8)    ql_db, qi_db, qv_db, qs_db, th_db
   real(r8)    wd2
   
   ! --------------------------------------------------------------------------------- !
   ! Variables related to the evaporation of precipitation within convective downdraft !
   ! --------------------------------------------------------------------------------- !

   real(r8)    evplflux                         !  Evaporation + Production of rain within downdraft in each layer [kg/s/m^2]
   real(r8)    evpiflux                         !  Evaporation + Production of snow within downdraft in each layer [kg/s/m^2]
   real(r8)    evptrflux(ncnst)                 !  Evaporation + Production + Wet Deposition of tracers within downdraft in each layer [kg(#)/s/m^2]
   real(r8)    fevp1_t_rate, fevp2_t_rate       !  Evaporation rate of precipitation within downdraft at the  top interface [(kg/kg)/Pa] >= 0.
   real(r8)    fevp1_b_rate, fevp2_b_rate       !  Evaporation rate of precipitation within downdraft at the base interface [(kg/kg)/Pa] >= 0.

   ! ------------------------ !
   ! Turbulent flux variables !
   ! ------------------------ !

   real(r8)    slflx_u(0:mkx)                   !  Convective updraft liquid static energy flux [J/s/m2]
   real(r8)    qtflx_u(0:mkx)                   !  Convective updraft total water flux [kg/s/m2]
   real(r8)    uflx_u(0:mkx)                    !  Convective updraft zonal momentum flux [kg m/s/s/m2]
   real(r8)    vflx_u(0:mkx)                    !  Convective updraft meridional momentum flux [kg m/s/s/m2]
   real(r8)    qlflx_u(0:mkx)                   !  Convective updraft ql flux [kg/s/m2]
   real(r8)    qiflx_u(0:mkx)                   !  Convective updraft qi flux [kg/s/m2]
   real(r8)    trflx_u(0:mkx,ncnst)             !  Convective updraft tracer flux [ #/s/m2, kg/s/m2 ]

   real(r8)    slflx_d(0:mkx)                   !  Convective downdraft liquid static energy flux [J/s/m2]
   real(r8)    qtflx_d(0:mkx)                   !  Convective downdraft total water flux [kg/s/m2]
   real(r8)    uflx_d(0:mkx)                    !  Convective downdraft zonal momentum flux [kg m/s/s/m2]
   real(r8)    vflx_d(0:mkx)                    !  Convective downdraft meridional momentum flux [kg m/s/s/m2]
   real(r8)    qlflx_d(0:mkx)                   !  Convective downdraft ql flux [kg/s/m2]
   real(r8)    qiflx_d(0:mkx)                   !  Convective downdraft qi flux [kg/s/m2]
   real(r8)    trflx_d(0:mkx,ncnst)             !  Convective downdraft tracer flux [ #/s/m2, kg/s/m2 ]

   real(r8)    uflx(0:mkx)                      !  Reconstructed convective zonal momentum flux [kg m/s/s/m2]
   real(r8)    vflx(0:mkx)                      !  Reconstructed convective meridional momentum flux [kg m/s/s/m2]

   real(r8)    thlflx_d_org_pblh                !  Adiabatic organization forcing. Convective downdraft flux of thl at the PBL top [kg*K/s/m2]
   real(r8)    qtflx_d_org_pblh                 !  Adiabatic organization forcing. Convective downdraft flux of qt  at the PBL top [kg*(kg/kg)/s/m2]
   real(r8)    uflx_d_org_pblh                  !  Adiabatic organization forcing. Convective downdraft flux of u   at the PBL top [kg*(m/s)/s/m2]
   real(r8)    vflx_d_org_pblh                  !  Adiabatic organization forcing. Convective downdraft flux of v   at the PBL top [kg*(m/s)/s/m2]
   real(r8)    trflx_d_org_pblh(ncnst)          !  Adiabatic organization forcing. Convective downdraft flux of tr  at the PBL top [kg*(kg/kg)/s/m2, kg*(#/kg)/s/m2]

   real(r8)    cmf_d_org_pblh                   !  Organization-inducing downdraft mass flux at the PBL top interface [kg/s/m2]               
   real(r8)    thl_d_org_pblh                
   real(r8)    qt_d_org_pblh                
   real(r8)    u_d_org_pblh                 
   real(r8)    v_d_org_pblh                 
   real(r8)    tr_d_org_pblh(ncnst)         

   real(r8)    thl_dia_d_org                    !  Diabatic organization forcing of thl      within convective downdraft [ K / s ]
   real(r8)    qt_dia_d_org                     !  Diabatic organization forcing of qt       within convective downdraft [ kg / kg / s ]
   real(r8)    tr_dia_d_org(ncnst)              !  Diabatic organization forcing of tracers  within convective downdraft [ kg / kg / s  or # / kg / s ]

   real(r8)    thl_dia_und_org                  !  Diabatic organization forcing of thl      within convective updraft and downdraft [ K / s ]
   real(r8)    qt_dia_und_org                   !  Diabatic organization forcing of qt       within convective updraft and downdraft [ kg / kg / s ]
   real(r8)    tr_dia_und_org(ncnst)            !  Diabatic organization forcing of tracers  within convective updraft and downdraft [ kg / kg / s  or # / kg / s ]

   real(r8)    thl_dia_env_org                  !  Diabatic organization forcing of thl      within environment [ K / s ]
   real(r8)    qt_dia_env_org                   !  Diabatic organization forcing of qt       within environment [ kg / kg / s ]
   real(r8)    tr_dia_env_org(ncnst)            !  Diabatic organization forcing of tracers  within environment [ kg / kg / s  or # / kg / s ]

 ! May.1.2014. Below '_d_orgU' and '_u_org' variables are added for 'budget consistent coldpool' treatment (i_budget_coldpool = 1,2 ).

   real(r8)    thlflx_d_orgU_pblh               !  Adiabatic organization forcing. Convective downdraft flux of thl at the PBL top [kg*K/s/m2]
   real(r8)    qtflx_d_orgU_pblh                !  Adiabatic organization forcing. Convective downdraft flux of qt  at the PBL top [kg*(kg/kg)/s/m2]
   real(r8)    uflx_d_orgU_pblh                 !  Adiabatic organization forcing. Convective downdraft flux of u   at the PBL top [kg*(m/s)/s/m2]
   real(r8)    vflx_d_orgU_pblh                 !  Adiabatic organization forcing. Convective downdraft flux of v   at the PBL top [kg*(m/s)/s/m2]
   real(r8)    trflx_d_orgU_pblh(ncnst)         !  Adiabatic organization forcing. Convective downdraft flux of tr  at the PBL top [kg*(kg/kg)/s/m2, kg*(#/kg)/s/m2]

   real(r8)    cmf_d_orgU_pblh                  !  Organization-inducing downdraft mass flux at the PBL top interface [kg/s/m2]               
   real(r8)    thl_d_orgU_pblh                
   real(r8)    qt_d_orgU_pblh                
   real(r8)    u_d_orgU_pblh                 
   real(r8)    v_d_orgU_pblh                 
   real(r8)    tr_d_orgU_pblh(ncnst)         

   real(r8)    thl_dia_d_orgU                   !  Diabatic organization forcing of thl      within convective downdraft [ K / s ]
   real(r8)    qt_dia_d_orgU                    !  Diabatic organization forcing of qt       within convective downdraft [ kg / kg / s ]
   real(r8)    tr_dia_d_orgU(ncnst)             !  Diabatic organization forcing of tracers  within convective downdraft [ kg / kg / s  or # / kg / s ]

   real(r8)    cmf_u_org_pblh                   !  = cmf_u(kpblhm)
   real(r8)    thlflx_u_org_pblh                !  Adiabatic organization forcing. Convective updraft flux of thl at the PBL top [kg*K/s/m2]
   real(r8)    qtflx_u_org_pblh                 !  Adiabatic organization forcing. Convective updraft flux of qt  at the PBL top [kg*(kg/kg)/s/m2]
   real(r8)    uflx_u_org_pblh                  !  Adiabatic organization forcing. Convective updraft flux of u   at the PBL top [kg*(m/s)/s/m2]
   real(r8)    vflx_u_org_pblh                  !  Adiabatic organization forcing. Convective updraft flux of v   at the PBL top [kg*(m/s)/s/m2]
   real(r8)    trflx_u_org_pblh(ncnst)          !  Adiabatic organization forcing. Convective updraft flux of tr  at the PBL top [kg*(kg/kg)/s/m2, kg*(#/kg)/s/m2]

   ! ------------------ !
   ! Tendency variables !
   ! ------------------ !

   real(r8)    qvten(mkx)                       !  Total tendency of qv [kg/kg/s]
   real(r8)    qlten(mkx)                       !  Total tendency of ql [kg/kg/s]
   real(r8)    qiten(mkx)                       !  Total tendency of qi [kg/kg/s]
   real(r8)    sten(mkx)                        !  Total tendency of s  [J/kg/s]
   real(r8)    uten(mkx)                        !  Total tendency of u  [m/s/s]
   real(r8)    vten(mkx)                        !  Total tendency of v  [m/s/s]
   real(r8)    trten(mkx,ncnst)                 !  Total tendency of tracers [ #/kg/s, kg/kg/s ]

   real(r8)    slten_u(mkx)                     !  Tendency of sl by updraft mass flux [J/kg/s]
   real(r8)    qtten_u(mkx)                     !  Tendency of qt by updraft mass flux [kg/kg/s]
   real(r8)    uten_u(mkx)                      !  Tendency of u  by updraft mass flux [m/s/s]
   real(r8)    vten_u(mkx)                      !  Tendency of v  by updraft mass flux [m/s/s]
   real(r8)    sten_u(mkx)                      !  Tendency of s  by updraft mass flux [J/kg/s]
   real(r8)    qvten_u(mkx)                     !  Tendency of qv by updraft mass flux [kg/kg/s]
   real(r8)    qlten_u(mkx)                     !  Tendency of ql by updraft mass flux [kg/kg/s]
   real(r8)    qiten_u(mkx)                     !  Tendency of qi by updraft mass flux [kg/kg/s]
   real(r8)    trten_u(mkx,ncnst)               !  Tendency of tracer by updraft mass flux [#/kg/s, kg/kg/s]

   real(r8)    slten_d(mkx)                     !  Tendency of sl by downdraft mass flux [J/kg/s]
   real(r8)    qtten_d(mkx)                     !  Tendency of qt by downdraft mass flux [kg/kg/s]
   real(r8)    uten_d(mkx)                      !  Tendency of u  by downdraft mass flux [m/s/s]
   real(r8)    vten_d(mkx)                      !  Tendency of v  by downdraft mass flux [m/s/s]
   real(r8)    sten_d(mkx)                      !  Tendency of s  by downdraft mass flux [J/kg/s]
   real(r8)    qvten_d(mkx)                     !  Tendency of qv by downdraft mass flux [kg/kg/s]
   real(r8)    qlten_d(mkx)                     !  Tendency of ql by downdraft mass flux [kg/kg/s]
   real(r8)    qiten_d(mkx)                     !  Tendency of qi by downdraft mass flux [kg/kg/s]
   real(r8)    trten_d(mkx,ncnst)               !  Tendency of tracer by downdraft mass flux [#/kg/s, kg/kg/s]

   real(r8)    slten_evp(mkx)                   !  Tendency of sl by convective precipitation and evaporation of convective precip. [J/kg/s]
   real(r8)    qtten_evp(mkx)                   !  Tendency of qt by convective precipitation and evaporation of convective precip. [kg/kg/s]
   real(r8)    uten_evp(mkx)                    !  Tendency of u  by convective precipitation and evaporation of convective precip. [m/s/s]
   real(r8)    vten_evp(mkx)                    !  Tendency of v  by convective precipitation and evaporation of convective precip. [m/s/s]
   real(r8)    sten_evp(mkx)                    !  Tendency of s  by convective precipitation and evaporation of convective precip. [J/kg/s]
   real(r8)    qvten_evp(mkx)                   !  Tendency of qv by convective precipitation and evaporation of convective precip. [kg/kg/s]
   real(r8)    qlten_evp(mkx)                   !  Tendency of ql by convective precipitation and evaporation of convective precip. [kg/kg/s]
   real(r8)    qiten_evp(mkx)                   !  Tendency of qi by convective precipitation and evaporation of convective precip. [kg/kg/s]
   real(r8)    trten_evp(mkx,ncnst)             !  Tendency of tracer by convective precipitation and evaporation of convective precip. [#/kg/s, kg/kg/s]
   real(r8)    trten_wdep(mkx,ncnst)            !  Tendency of tracer by wet deposition within environment by convective precip. [#/kg/s, kg/kg/s]

   real(r8)    slten_dis(mkx)                   !  Tendency of sl by dissipative heating of mean kinetic energy [J/kg/s]
   real(r8)    qtten_dis(mkx)                   !  Tendency of qt by dissipative heating of mean kinetic energy [kg/kg/s]
   real(r8)    uten_dis(mkx)                    !  Tendency of u  by dissipative heating of mean kinetic energy [m/s/s]
   real(r8)    vten_dis(mkx)                    !  Tendency of v  by dissipative heating of mean kinetic energy [m/s/s]
   real(r8)    sten_dis(mkx)                    !  Tendency of s  by dissipative heating of mean kinetic energy [J/kg/s]
   real(r8)    qvten_dis(mkx)                   !  Tendency of qv by dissipative heating of mean kinetic energy [kg/kg/s]
   real(r8)    qlten_dis(mkx)                   !  Tendency of ql by dissipative heating of mean kinetic energy [kg/kg/s]
   real(r8)    qiten_dis(mkx)                   !  Tendency of qi by dissipative heating of mean kinetic energy [kg/kg/s]
   real(r8)    trten_dis(mkx,ncnst)             !  Tendency of tracer by dissipative heating of mean kinetic energy [#/kg/s, kg/kg/s]

   real(r8)    slten_par(mkx)                   !  Tendency of sl by partitioning the tendency in the lowest layer in the layers within the PBL [J/kg/s]
   real(r8)    qtten_par(mkx)                   !  Tendency of qt by partitioning the tendency in the lowest layer in the layers within the PBL [kg/kg/s]
   real(r8)    qlten_par(mkx)                   !  Tendency of ql by partitioning the tendency in the lowest layer in the layers within the PBL [kg/kg/s]
   real(r8)    qiten_par(mkx)                   !  Tendency of qi by partitioning the tendency in the lowest layer in the layers within the PBL [kg/kg/s]
   real(r8)    uten_par(mkx)                    !  Tendency of u  by partitioning the tendency in the lowest layer in the layers within the PBL [m/s/s]
   real(r8)    vten_par(mkx)                    !  Tendency of v  by partitioning the tendency in the lowest layer in the layers within the PBL [m/s/s]
   real(r8)    trten_par(mkx,ncnst)             !  Tendency of tracer by partitioning the tendency in the lowest layer in the layers within the PBL [#/kg/s, kg/kg/s]

   real(r8)    slten_NUM(mkx)                   !  Final Numerical tendency of sl [J/kg/s]
   real(r8)    qtten_NUM(mkx)                   !  Final Numerical tendency of qt [kg/kg/s]
   real(r8)    uten_NUM(mkx)                    !  Final Numerical tendency of u  [m/s/s]
   real(r8)    vten_NUM(mkx)                    !  Final Numerical tendency of v  [m/s/s]
   real(r8)    qvten_NUM(mkx)                   !  Final Numerical tendency of qv [kg/kg/s]
   real(r8)    qlten_NUM(mkx)                   !  Final Numerical tendency of ql [kg/kg/s]
   real(r8)    qiten_NUM(mkx)                   !  Final Numerical tendency of qi [kg/kg/s]
   real(r8)    sten_NUM(mkx)                    !  Final Numerical tendency of s  [J/kg/s]
   real(r8)    trten_NUM(mkx,ncnst)             !  Final Numerical tendency of tracer [#/kg/s, kg/kg/s]

   real(r8)    qrten(mkx)                       !  Production rate of rain by the expels of in-cumulus excessive condensate [kg/kg/s]
   real(r8)    qsten(mkx)                       !  Production rate of snow by the expels of in-cumulus excessive condensate [kg/kg/s]

   real(r8)    qrten_u(mkx)                     !  Production rate of rain by the expels of in-cumulus excessive condensate within updraft [kg/kg/s]
   real(r8)    qsten_u(mkx)                     !  Production rate of snow by the expels of in-cumulus excessive condensate within updraft [kg/kg/s]

   real(r8)    qrten_u_msfc(mkx,nseg)           !  Production rate of rain by the expels of in-cumulus excessive condensate within updraft for each original updraft segment [kg/kg/s]
   real(r8)    qsten_u_msfc(mkx,nseg)           !  Production rate of snow by the expels of in-cumulus excessive condensate within updraft for each original updraft segment [kg/kg/s]
   real(r8)    trrsten_u_msfc(mkx,nseg,ncnst)   !  Production rate of precipitating tracers by the expels of in-cumulus tracers within updraft for each original updraft segment [kg/kg/s,#/kg/s]

   real(r8)    qrten_d(mkx)                     !  Production rate of rain by the expels of in-cumulus excessive condensate within downdraft [kg/kg/s]
   real(r8)    qsten_d(mkx)                     !  Production rate of snow by the expels of in-cumulus excessive condensate within downdraft [kg/kg/s]

   real(r8)    qrten_d_msfc(mkx,nseg)           !  Production rate of rain by the expels of in-cumulus excessive condensate within downdraft for each original updraft segment [kg/kg/s]
   real(r8)    qsten_d_msfc(mkx,nseg)           !  Production rate of snow by the expels of in-cumulus excessive condensate within downdraft for each original updraft segment [kg/kg/s]
   real(r8)    trrsten_d_msfc(mkx,nseg,ncnst)   !  Production rate of tracers by the expels of in-cumulus excessive tracers within downdraft for each original updraft segment [kg(#)/kg/s]

   real(r8)    snowmlt_e(mkx)                   !  Snow melting tendency within environment before evaporation within downdraft [kg/kg/s]
   real(r8)    snowmlt_e_msfc(mkx,nseg)         !  Snow melting tendency within environment before evaporation within downdraft for each original updraft segment [kg/kg/s]

   real(r8)    thlten_dia_u(mkx)                !  Diabatic tendency of thl within   updraft [K/s]
   real(r8)    thlten_dia_d(mkx)                !  Diabatic tendency of thl within downdraft [K/s]
   
   real(r8)    qtten_dia_u(mkx)                 !  Diabatic tendency of qt within   updraft [K/s]
   real(r8)    qtten_dia_d(mkx)                 !  Diabatic tendency of qt within downdraft [K/s]

   real(r8)    qlten_dia_u(mkx)                 !  Diabatic tendency of ql within   updraft ( exclude effective tendency ) [kg/kg/s]
   real(r8)    qlten_dia_d(mkx)                 !  Diabatic tendency of ql within downdraft ( exclude effective tendency ) [kg/kg/s]

   real(r8)    qiten_dia_u(mkx)                 !  Diabatic tendency of qi within   updraft ( exclude effective tendency ) [kg/kg/s]
   real(r8)    qiten_dia_d(mkx)                 !  Diabatic tendency of qi within downdraft ( exclude effective tendency ) [kg/kg/s]

   real(r8)    trten_dia_u(mkx,ncnst)           !  Diabatic tendency of tracers within   updraft ( exclude effective tendency ) [#/kg/s, kg/kg/s]
   real(r8)    trten_dia_d(mkx,ncnst)           !  Diabatic tendency of tracers within downdraft ( exclude effective tendency ) [#/kg/s, kg/kg/s]

   real(r8)    ntraprd(mkx)                     !  Net production rate of rain ( qrten(k) + snowmlt(k) - evprain(k) ) [kg/kg/s]
   real(r8)    ntsnprd(mkx)                     !  Net production rate of snow ( qsten(k) - snowmlt(k) - evpsnow(k) ) [kg/kg/s]
   real(r8)    nttrrsprd(mkx,ncnst)             !  Net production rate of tracer ( trrsten(k) - evptrrs + wdeptrrs ) [kg(#)/kg/s]

   real(r8)    ntraprd_msfc(mkx,nseg)           !  Net production rate of rain ( qrten(k) + snowmlt(k) - evprain ) for each original updraft segment [kg/kg/s]
   real(r8)    ntsnprd_msfc(mkx,nseg)           !  Net production rate of snow ( qsten(k) - snowmlt(k) - evpsnow ) for each original updraft segment [kg/kg/s]
   real(r8)    nttrrsprd_msfc(mkx,nseg,ncnst)   !  Net production rate of tracer ( trrsten(k) - evptrrs + wdeptrrs ) for each original updraft segment [kg(#)/kg/s]

   real(r8)    evprain_e(mkx)                   !  Evaporation rate of rain in the environment [kg/kg/s]
   real(r8)    evpsnow_e(mkx)                   !  Evaporation rate of snow in the environment [kg/kg/s]
   real(r8)    evptrrs_e(mkx,ncnst)             !  Evaporation rate of tracers in the environment [kg(#)/kg/s] > 0.
   real(r8)    wdeptrrs_e(mkx,ncnst)            !  Wet deposition rate of tracers in the environment [kg(#)/kg/s] > 0.

   real(r8)    evprain_d(mkx)                   !  Evaporation rate of rain in the downdraft [kg/kg/s]
   real(r8)    evpsnow_d(mkx)                   !  Evaporation rate of snow in the downdraft [kg/kg/s]
   real(r8)    evptrrs_d(mkx,ncnst)             !  Evaporation rate of tracers in the downdraft [kg(#)/kg/s] > 0.

   real(r8)    evprain_e_msfc(mkx,nseg)         !  Evaporation rate of rain for each original updraft segment in the environment [kg/kg/s]
   real(r8)    evpsnow_e_msfc(mkx,nseg)         !  Evaporation rate of snow for each original updraft segment in the environment [kg/kg/s]
   real(r8)    evptrrs_e_msfc(mkx,nseg,ncnst)   !  Evaporation rate of tracer for each original updraft segment in the environment [kg(#)/kg/s] > 0.
   real(r8)    wdeptrrs_e_msfc(mkx,nseg,ncnst)  !  Wet deposition rate of tracer for each original updraft segment in the environment [kg(#)/kg/s] > 0.

   real(r8)    evprain_d_msfc(mkx,nseg)         !  Evaporation rate of rain for each original updraft segment in the downdraft [kg/kg/s]
   real(r8)    evpsnow_d_msfc(mkx,nseg)         !  Evaporation rate of snow for each original updraft segment in the downdraft [kg/kg/s]
   real(r8)    evptrrs_d_msfc(mkx,nseg,ncnst)   !  Evaporation rate of tracer for each original updraft segment in the downdraft [kg(#)/kg/s] > 0.

   real(r8)    cvp_rainprd(mkx)                 !  Corrective production of rain from environmental qv0 [kg/kg/s]
   real(r8)    cvp_snowprd(mkx)                 !  Corrective production of snow from environmental qv0 [kg/kg/s]
   real(r8)    cvp_trrsprd(mkx,ncnst)           !  Corrective production of tracer from environmental tr0 [kg(#)/kg/s]

   real(r8)    cvp_rainprd_msfc(mkx,nseg)       !  Corrective production of rain from environmental qv0 for each original updraft segment [kg/kg/s]
   real(r8)    cvp_snowprd_msfc(mkx,nseg)       !  Corrective production of snow from environmental qv0 for each original updraft segment [kg/kg/s]
   real(r8)    cvp_trrsprd_msfc(mkx,nseg,ncnst) !  Corrective production of tracers from environmental tr0 for each original updraft segment [kg(#)/kg/s]

   real(r8)    qlten_eff_u(mkx)                 !  Effective diabatic tendency on the   updraft cloud condensate 'ql' [kg/kg/s]
   real(r8)    qiten_eff_u(mkx)                 !  Effective diabatic tendency on the   updraft cloud condensate 'qi' [kg/kg/s]

   real(r8)    qlten_eff_d(mkx)                 !  Effective diabatic tendency on the downdraft cloud condensate 'ql' [kg/kg/s]
   real(r8)    qiten_eff_d(mkx)                 !  Effective diabatic tendency on the downdraft cloud condensate 'qi' [kg/kg/s]

   real(r8)    trten_eff_u(mkx,ncnst)           !  Effective diabatic tendency on the   updraft tracers [#/kg/s, kg/kg/s]
   real(r8)    trten_eff_d(mkx,ncnst)           !  Effective diabatic tendency on the downdraft tracers [#/kg/s, kg/kg/s]

   real(r8)    uf(mkx)                          !  Provisional environmental      zonal wind
   real(r8)    vf(mkx)                          !  Provisional environmental meridional wind

   real(r8)    ql_env_ua(mkx)                   !  ql  of imported airs into the layer due to the compensating subsidence by   updraft mass flux at the top interface
   real(r8)    qi_env_ua(mkx)                   !  qi  of imported airs into the layer due to the compensating subsidence by   updraft mass flux at the top interface

   real(r8)    ql_env_da(mkx)                   !  ql  of imported airs into the layer due to the compensating subsidence by downdraft mass flux at the top interface
   real(r8)    qi_env_da(mkx)                   !  qi  of imported airs into the layer due to the compensating subsidence by downdraft mass flux at the top interface

   ! ------------------------------------------------------------------------------------ !
   ! Variables associated with mixing with multiple mixing environmental airs ( '_mxen' ) !
   ! ------------------------------------------------------------------------------------ !

   integer     ktop_mxen(niter)                 !  The top layer index of convective updraft

   real(r8)    cuorg_mxen                       !  Temporary convective organization used only for multiple mixing. Either 0 ( when iter = 1 ) or 1 ( when iter = 2 ) 
   real(r8)    cush_mxen(niter)                 !  Maximum updraft top height [ m ]
   real(r8)    cushavg_mxen(niter)              !  Updraft top height weighted by updraft mass flux at surface [ m ]

   real(r8)    cu_cmfum_mxen(mkx,niter)         !  Total amount of updraft mass flux involved in the buoyancy sorting [kg/m2/s]
   real(r8)    cu_cmfr_mxen(mkx,niter)          !  Total amount of detrained mass into the environment ( may also contain updraft properties ) [kg/m2/s]
   real(r8)    cu_thlr_mxen(mkx,niter)          !  Mass flux weighted anomalous conservative scalar of detrained airs
   real(r8)    cu_qtr_mxen(mkx,niter)       
   real(r8)    cu_ur_mxen(mkx,niter)         
   real(r8)    cu_vr_mxen(mkx,niter)         
   real(r8)    cu_qlr_mxen(mkx,niter)        
   real(r8)    cu_qir_mxen(mkx,niter)        
   real(r8)    cu_trr_mxen(mkx,ncnst,niter)       
   real(r8)    cu_cmfrd_mxen(mkx,niter)         !  Total amount of detrained mass into the environment from convective downdraft [kg/m2/s]
   real(r8)    cu_thlrd_mxen(mkx,niter)      
   real(r8)    cu_qtrd_mxen(mkx,niter)       
   real(r8)    cu_urd_mxen(mkx,niter)         
   real(r8)    cu_vrd_mxen(mkx,niter)         
   real(r8)    cu_qlrd_mxen(mkx,niter)        
   real(r8)    cu_qird_mxen(mkx,niter)        
   real(r8)    cu_trrd_mxen(mkx,ncnst,niter)

   real(r8)    cmf_u_mxen(0:mkx,niter)          
   real(r8)    cmf_d_mxen(0:mkx,niter)          
   real(r8)    slflx_u_mxen(0:mkx,niter)        
   real(r8)    slflx_d_mxen(0:mkx,niter)        
   real(r8)    qtflx_u_mxen(0:mkx,niter)        
   real(r8)    qtflx_d_mxen(0:mkx,niter)        
   real(r8)    uflx_u_mxen(0:mkx,niter)        
   real(r8)    uflx_d_mxen(0:mkx,niter)        
   real(r8)    vflx_u_mxen(0:mkx,niter)        
   real(r8)    vflx_d_mxen(0:mkx,niter)        

   real(r8)    flxrain_u_mxen(0:mkx,niter)        
   real(r8)    flxsnow_u_mxen(0:mkx,niter)        

   real(r8)    thl_orgforce_mxen(niter)          !  Total forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)    qt_orgforce_mxen(niter)           !  Total forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)    u_orgforce_mxen(niter)            !  Total forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)    v_orgforce_mxen(niter)            !  Total forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)    tr_orgforce_mxen(ncnst,niter)     !  Total forcing for organized difference between 'off-wake' and 'grid-mean' tracers [ kg / kg / s or # / kg / s ]    
   real(r8)    awk_orgforce_mxen(niter)          !  Total forcing for wake area [ 1 / s ]

  ! Below block is for detailed diagnostic output

   real(r8)    thl_orgforce_flx_mxen(niter)      !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)    qt_orgforce_flx_mxen(niter)       !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)    u_orgforce_flx_mxen(niter)        !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)    v_orgforce_flx_mxen(niter)        !  PBL top flux-related forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)    awk_orgforce_flx_mxen(niter)      !  PBL top flux-related forcing for wake area [ 1 / s ]

   real(r8)    thl_orgforce_und_mxen(niter)      !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)    qt_orgforce_und_mxen(niter)       !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)    u_orgforce_und_mxen(niter)        !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)    v_orgforce_und_mxen(niter)        !  Up-and-Down diabatic forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)    awk_orgforce_mix_mxen(niter)      !  Lateral-Mixing       forcing for wake area [ 1 / s ]

   real(r8)    thl_orgforce_env_mxen(niter)      !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' thl [ K / s ]
   real(r8)    qt_orgforce_env_mxen(niter)       !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' qt [ kg / kg / s ]
   real(r8)    u_orgforce_env_mxen(niter)        !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' u [ m / s / s ]
   real(r8)    v_orgforce_env_mxen(niter)        !  Environment diabatic forcing for organized difference between 'off-wake' and 'grid-mean' v [ m / s / s ]
   real(r8)    cmf_d_org_pblh_mxen(niter)        !  Organization-inducing downdraft mass flux at the PBL top interface [ kg / m^2 / s ] 

  ! Above block is for detailed diagnostic output

   real(r8)    taui_thl_mxen(niter)              !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' thl [ 1 / s ]
   real(r8)    taui_qt_mxen(niter)               !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' qt [ 1 / s ]
   real(r8)    taui_u_mxen(niter)                !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' u [ 1 / s ]
   real(r8)    taui_v_mxen(niter)                !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' v [ 1 / s ]
   real(r8)    taui_tr_mxen(ncnst,niter)         !  Inverse of damping time scale of the difference between 'off-wake' and 'grid-mean' tracers [ 1 / s ]
   real(r8)    taui_awk_mxen(niter)              !  Inverse of damping time scale of the wake area [ 1 / s ]

   real(r8)    del_org_mxen(niter)               !  Detrainment rate of the cold pool [ 1 / s ]
   real(r8)    del0_org_mxen(niter)              !  Effective detrainment rate of the cold pool [ 1 / s ]

   real(r8)    qvten_mxen(mkx,niter)          
   real(r8)    qlten_mxen(mkx,niter)          
   real(r8)    qiten_mxen(mkx,niter)          
   real(r8)    trten_mxen(mkx,ncnst,niter)          
   real(r8)    sten_mxen(mkx,niter)           
   real(r8)    uten_mxen(mkx,niter)           
   real(r8)    vten_mxen(mkx,niter)           
   real(r8)    qrten_mxen(mkx,niter)          
   real(r8)    qsten_mxen(mkx,niter)          

   real(r8)    rqc_l_mxen(mkx,niter)        
   real(r8)    rqc_i_mxen(mkx,niter)        
   real(r8)    rqc_mxen(mkx,niter)          
   real(r8)    rnc_l_mxen(mkx,niter)        
   real(r8)    rnc_i_mxen(mkx,niter)        

   real(r8)    cmf_det_mxen(mkx,niter)        
   real(r8)    ql_det_mxen(mkx,niter)        
   real(r8)    qi_det_mxen(mkx,niter)        

   real(r8)    evapc_mxen(mkx,niter)       

   real(r8)    am_u_mxen(mkx,niter)      
   real(r8)    qlm_u_mxen(mkx,niter)     
   real(r8)    qim_u_mxen(mkx,niter)     

   real(r8)    am_d_mxen(mkx,niter)      
   real(r8)    qlm_d_mxen(mkx,niter)     
   real(r8)    qim_d_mxen(mkx,niter)     

   real(r8)    rliq_mxen(niter)           
   real(r8)    precip_mxen(niter)         
   real(r8)    snow_mxen(niter)           

   real(r8)    cnt_mxen(niter)            
   real(r8)    cnb_mxen(niter)            

   real(r8)    slten_u_mxen(mkx,niter)    
   real(r8)    qtten_u_mxen(mkx,niter)    
   real(r8)    uten_u_mxen(mkx,niter)     
   real(r8)    vten_u_mxen(mkx,niter)     
   real(r8)    sten_u_mxen(mkx,niter)     
   real(r8)    qvten_u_mxen(mkx,niter)    
   real(r8)    qlten_u_mxen(mkx,niter)    
   real(r8)    qiten_u_mxen(mkx,niter)    
   real(r8)    trten_u_mxen(mkx,ncnst,niter) 

   real(r8)    slten_d_mxen(mkx,niter)     
   real(r8)    qtten_d_mxen(mkx,niter)     
   real(r8)    uten_d_mxen(mkx,niter)      
   real(r8)    vten_d_mxen(mkx,niter)      
   real(r8)    sten_d_mxen(mkx,niter)      
   real(r8)    qvten_d_mxen(mkx,niter)     
   real(r8)    qlten_d_mxen(mkx,niter)     
   real(r8)    qiten_d_mxen(mkx,niter)     
   real(r8)    trten_d_mxen(mkx,ncnst,niter)  

   real(r8)    slten_evp_mxen(mkx,niter)  
   real(r8)    qtten_evp_mxen(mkx,niter)  
   real(r8)    uten_evp_mxen(mkx,niter)   
   real(r8)    vten_evp_mxen(mkx,niter)   
   real(r8)    sten_evp_mxen(mkx,niter)   
   real(r8)    qvten_evp_mxen(mkx,niter)  
   real(r8)    qlten_evp_mxen(mkx,niter)  
   real(r8)    qiten_evp_mxen(mkx,niter)  
   real(r8)    trten_evp_mxen(mkx,ncnst,niter) 

   real(r8)    qlten_sub_mxen(mkx,niter)     
   real(r8)    qiten_sub_mxen(mkx,niter)     

   real(r8)    qlten_det_mxen(mkx,niter)     
   real(r8)    qiten_det_mxen(mkx,niter)     

   real(r8)    thl_u_mxen(0:mkx,niter)      
   real(r8)    qt_u_mxen(0:mkx,niter)       
   real(r8)    u_u_mxen(0:mkx,niter)        
   real(r8)    v_u_mxen(0:mkx,niter)        
   real(r8)    w_u_mxen(0:mkx,niter)        
   real(r8)    ql_u_mxen(0:mkx,niter)       
   real(r8)    qi_u_mxen(0:mkx,niter)       
   real(r8)    tr_u_mxen(0:mkx,ncnst,niter)    
   real(r8)    a_u_mxen(0:mkx,niter)      
   real(r8)    num_u_mxen(0:mkx,niter)    
   real(r8)    wa_u_mxen(0:mkx,niter)     
   real(r8)    qla_u_mxen(0:mkx,niter)    
   real(r8)    qia_u_mxen(0:mkx,niter)    
   real(r8)    rad_u_mxen(0:mkx,niter)    
   real(r8)    thva_u_mxen(0:mkx,niter)     

   real(r8)    a_p_mxen(0:mkx,niter)      
   real(r8)    am_evp_mxen(mkx,niter)      
   real(r8)    am_pu_mxen(mkx,niter)      
   real(r8)    x_p_mxen(0:mkx,niter)      
   real(r8)    y_p_mxen(0:mkx,niter)      
   real(r8)    x_um_mxen(mkx,niter)      
   real(r8)    y_um_mxen(mkx,niter)      

   real(r8)    thl_d_mxen(0:mkx,niter)     
   real(r8)    qt_d_mxen(0:mkx,niter)      
   real(r8)    u_d_mxen(0:mkx,niter)       
   real(r8)    v_d_mxen(0:mkx,niter)       
   real(r8)    w_d_mxen(0:mkx,niter)       
   real(r8)    ql_d_mxen(0:mkx,niter)      
   real(r8)    qi_d_mxen(0:mkx,niter)      
   real(r8)    tr_d_mxen(0:mkx,ncnst,niter)      
   real(r8)    a_d_mxen(0:mkx,niter)          
   real(r8)    wa_d_mxen(0:mkx,niter)         
   real(r8)    qla_d_mxen(0:mkx,niter)        
   real(r8)    qia_d_mxen(0:mkx,niter)        

   real(r8)    thl_u_msfc_mxen(0:mkx,nseg,niter) 
   real(r8)    qt_u_msfc_mxen(0:mkx,nseg,niter)  
   real(r8)    u_u_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    v_u_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    w_u_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    ql_u_msfc_mxen(0:mkx,nseg,niter)  
   real(r8)    qi_u_msfc_mxen(0:mkx,nseg,niter)  
   real(r8)    tr_u_msfc_mxen(0:mkx,nseg,ncnst,niter)
   real(r8)    cmf_u_msfc_mxen(0:mkx,nseg,niter)  
   real(r8)    a_u_msfc_mxen(0:mkx,nseg,niter)    
   real(r8)    num_u_msfc_mxen(0:mkx,nseg,niter)  
   real(r8)    rad_u_msfc_mxen(0:mkx,nseg,niter)  

   real(r8)    eps0_u_msfc_mxen(0:mkx,nseg,niter)    
   real(r8)    eps_u_msfc_mxen(0:mkx,nseg,niter)     
   real(r8)    del_u_msfc_mxen(0:mkx,nseg,niter)     
   real(r8)    eeps_u_msfc_mxen(0:mkx,nseg,niter)    
   real(r8)    ddel_u_msfc_mxen(0:mkx,nseg,niter)    
   real(r8)    xc_u_msfc_mxen(0:mkx,nseg,niter)      
   real(r8)    xs_u_msfc_mxen(0:mkx,nseg,niter)      
   real(r8)    xemin_u_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    xemax_u_msfc_mxen(0:mkx,nseg,niter)    
   real(r8)    cridis_u_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    thvcuenv_u_msfc_mxen(0:mkx,nseg,niter) 
   real(r8)    thvegenv_u_msfc_mxen(0:mkx,nseg,niter) 
   real(r8)    thvxsenv_u_msfc_mxen(0:mkx,nseg,niter) 
   real(r8)    fmix_u_msfc_mxen(0:mkx,nseg,niter)     
   real(r8)    cmfumix_u_msfc_mxen(0:mkx,nseg,niter)  

   real(r8)    thl_d_msfc_mxen(0:mkx,nseg,niter)       
   real(r8)    qt_d_msfc_mxen(0:mkx,nseg,niter)       
   real(r8)    u_d_msfc_mxen(0:mkx,nseg,niter)       
   real(r8)    v_d_msfc_mxen(0:mkx,nseg,niter)       
   real(r8)    w_d_msfc_mxen(0:mkx,nseg,niter)       
   real(r8)    ql_d_msfc_mxen(0:mkx,nseg,niter)      
   real(r8)    qi_d_msfc_mxen(0:mkx,nseg,niter)      
   real(r8)    tr_d_msfc_mxen(0:mkx,nseg,ncnst,niter) 
   real(r8)    cmf_d_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    a_d_msfc_mxen(0:mkx,nseg,niter)     
   real(r8)    wa_d_msfc_mxen(0:mkx,nseg,niter)    
   real(r8)    qla_d_msfc_mxen(0:mkx,nseg,niter)   
   real(r8)    qia_d_msfc_mxen(0:mkx,nseg,niter)   

   real(r8)    ptop_msfc_mxen(nseg,niter) 
   real(r8)    ztop_msfc_mxen(nseg,niter) 

   ! ---------------- !
   ! Single Variables !
   ! ---------------- !

   character(len=2)  numcha 
   integer     i, k, kv, ki, kvi, km, kp, ks, m, mm, ids, kc, mt, l, lspec, ixi, ixf
   integer     iter, iacc
   integer     nseg_det, nseg_nondet
   integer     ktop, ktop_up_par, ktop_dn_par, ks_top, ks_bot, msfc
   integer     ixcldliq, ixcldice, ixnumliq, ixnumice
   integer     i_awk, i_thl, i_qt, i_u, i_v

   integer     N_up(0:mkx)                      !  # of updraft segments at the base interface [ # ]

   integer     kpblh                            !  Layer index with PBL top in it or at the base interface 
   integer     kpblhm                           !  = kpblh - 1
   real(r8)    pblh                             !  PBL top height [ m ]
   real(r8)    pblhz                            !  Thickness of PBL depth in [ m ].  pblhz = zs0(kpblhm) - zs0(0).
   real(r8)    pblhp                            !  Thickness of PBL depth in [ Pa ]. pblhp = ps0(0) - ps0(kpblhm).
   real(r8)    went                             !  Entrainment rate at the PBL top interface directly from the UW PBL scheme [ m/s ]
   real(r8)    qflx                             !  Upward water vapor flux into atmosphere at surface [ kg/m2/s ]
   real(r8)    shflx                            !  Upward sensible heat flux into atmosphere at surface [ J/m2/s ]
   real(r8)    taux                             !  Upward zonal      wind stress into atmosphere at surface [ kg m/s /m2/s ]
   real(r8)    tauy                             !  Upward meridional wind stress into atmosphere at surface [ kg m/s /m2/s ] 
   real(r8)    aflx(ncnst)                      !  Upward tracer fluxes          into atmosphere at surface [ #/m2/s, kg/m2/s ]
   real(r8)    landfrac                         !  Land  Fraction [ fraction ]
   real(r8)    sgh30                            !  Standard deviation of subgrid topographic height at 30 s horizontal area [ meter ] 
   real(r8)    cush                             !  Input cumulus top height [ m ]
   real(r8)    cushavg                          !  Input mean cumulus top height weighted by updraft mass flux at surface [ m ]
   real(r8)    cuorg                            !  Input convective organization [ 0-1 ]
   real(r8)    awk_PBL_raw                      !  Wake area within PBL [ 0 - 1 ]
   real(r8)    delta_thl_PBL_raw                !  Difference of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
   real(r8)    delta_qt_PBL_raw                 !  Difference of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
   real(r8)    delta_u_PBL_raw                  !  Difference of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(r8)    delta_v_PBL_raw                  !  Difference of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(r8)    delta_thv_PBL_raw                !  Difference of thv between off-wake region and grid-mean value averaged over the PBL [ K ]
   real(r8)    delta_tr_PBL_raw(ncnst)          !  Difference of tr  between off-wake region and grid-mean value averaged over the PBL [ kg/kg, #/kg ]
   real(r8)    awk_PBL_max                      !  Maximum alloed wake area within PBL [ 0 - 1 ]
   real(r8)    awk_PBL                          !  Wake area within PBL [ 0 - 1 ]
   real(r8)    delta_thl_PBL                    !  Difference of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
   real(r8)    delta_qt_PBL                     !  Difference of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
   real(r8)    delta_u_PBL                      !  Difference of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(r8)    delta_v_PBL                      !  Difference of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(r8)    delta_thv_PBL                    !  Difference of thv between off-wake region and grid-mean value averaged over the PBL [ K ]
   real(r8)    delta_tr_PBL(ncnst)              !  Difference of tr  between off-wake region and grid-mean value averaged over the PBL [ kg/kg, #/kg ]
   real(r8)    delta_w_PBL                      !  Difference of w   between off-wake region and grid-mean value averaged over the PBL [ m/s ]

   real(r8)    cu_cmfr(mkx)      
   real(r8)    cu_thlr(mkx)      
   real(r8)    cu_qtr(mkx)       
   real(r8)    cu_ur(mkx)         
   real(r8)    cu_vr(mkx)         
   real(r8)    cu_qlr(mkx)        
   real(r8)    cu_qir(mkx)        
   real(r8)    cu_trr(mkx,ncnst)

   real(r8)    tke1                             !  TKE in the lowest model layer [ m2/s2 ] 
   real(r8)    wstar1                           !  wstar in the lowest model layer ( = ( 2.5 integral ( bprod * dz ) )^(1/3) ) [ m/s ]
   real(r8)    wstar2                           !  wstar2 = ( bprod_sfc * PBLH )^(1/3) [ m/s ] 
   real(r8)    tkePBL                           !  TKE within the PBL [ m2/s2 ] 
   real(r8)    wstarPBL                         !  wstar within the PBL [ m/s ]. In principle, should be the sams as the input wstar but not. 
   real(r8)    dpi, dzi
   real(r8)    qt0PBL, thl0PBL, u0PBL, v0PBL, tr0PBL(ncnst)
   real(r8)    qt0min_PBL, thl0min_PBL, tr0min_PBL(ncnst)
   real(r8)    eps_wk_eff, del_wk_eff
   real(r8)    a_oro                            !  Forbidden area from the subgrid scale variation of topography [ fraction ]
   real(r8)    a_forbid                         !  Total forbidden area ( awk_PBL + a_oro ) [ fraction ]
   real(r8)    c0_ac                            !  Auto-conversion efficiency of CAM5 deep convection scheme [ 1 / m ]
   real(r8)    criqc                            !  Critical condensate amount that updraft can hold [ kg / kg ] 
   real(r8)    au_base_max, au_base_min
   real(r8)    kevp_rain_dn, kevp_snow_dn
   real(r8)    kevp_rain, kevp_snow

   real(r8)    cnt, cnb                         !  Cloud top and base interface indices
   real(r8)    d_alpha
   real(r8)    z_b, z_t
   real(r8)    p_b, p_t
   real(r8)    dz_m, dp_m
   real(r8)    exn_b, exn_t
   real(r8)    thl_b
   real(r8)    qt_b
   real(r8)    tr_b(ncnst)
   real(r8)    u_b
   real(r8)    v_b
   real(r8)    ql_b
   real(r8)    qi_b
   real(r8)    thv_b, thv_t
   real(r8)    thv_mean_b, thv_mean_t
   real(r8)    thvbot, thvtop 
   real(r8)    plfc, plnb
   real(r8)    thvl_b, thvl_t
   real(r8)    rho_b, rho_m, rho_t
   real(r8)    thvl_minE, thv_minE
   real(r8)    dp, dz, rho
   real(r8)    thle_b(mkx)
   real(r8)    qte_b(mkx)
   real(r8)    tre_b(mkx,ncnst)
   real(r8)    ue_b(mkx)
   real(r8)    ve_b(mkx)
   real(r8)    we_b(mkx)
   real(r8)    we_t(mkx)
   real(r8)    qle_b(mkx)
   real(r8)    qie_b(mkx)
   real(r8)    ssthle(mkx)
   real(r8)    ssqte(mkx)
   real(r8)    ssue(mkx)
   real(r8)    ssve(mkx)
   real(r8)    ssqle(mkx)
   real(r8)    ssqie(mkx)
   real(r8)    sstre(mkx,ncnst)
   real(r8)    pe
   real(r8)    w_cu
   real(r8)    thl_cu, qt_cu
   real(r8)    ql_cu, qi_cu
   real(r8)    thl_eg, qt_eg, rh_eg
   real(r8)    thv_cu, thv_eg, thv_xs
   real(r8)    u_eg, v_eg, w_eg
   real(r8)    tr_eg(ncnst)
   real(r8)    thv_env
   real(r8)    cridis
   real(r8)    thl_cumC, thl_cumS
   real(r8)    qt_cumC, qt_cumS
   real(r8)    thv_cumC, thv_cumS, thv_cumE
   real(r8)    xdown_min, xdown_max
   real(r8)    zbar, zbar1, zbar2 
   real(r8)    zmass, zmass1, zmass2
   real(r8)    zmass_up, zmass_up1, zmass_up2
   real(r8)    th, qv, ql, qi, qse
   real(r8)    es, qs
   integer     id_check
   real(r8)    bogbot, bogtop, wu2 
   real(r8)    thv, thv_dt, thv_db
   real(r8)    thl_med, qt_med, ql_med, qi_med, qv_med, u_med, v_med
   real(r8)    thl_meu, qt_meu, ql_meu, qi_meu, u_meu, v_meu
   real(r8)    tr_med(ncnst), tr_meu(ncnst) 
   real(r8)    fac
   real(r8)    f_nu
   real(r8)    cmf_dt, thl_dt, qt_dt, u_dt, v_dt, w_dt
   real(r8)    tr_dt(ncnst)
   real(r8)    cmf_db, thl_db, qt_db, u_db, v_db, w_db
   real(r8)    tr_db(ncnst)
   real(r8)    qw_db, tw_db
   real(r8)    evp_thll, evp_qtl, evp_thli, evp_qti
   real(r8)    evp_qt
   real(r8)    evp_max
   real(r8)    evp_tr(ncnst)
   real(r8)    prep_thll, prep_qtl, prep_thli, prep_qti
   real(r8)    prep_tr(ncnst)
   real(r8)    exql, exqi, extr(ncnst)
   real(r8)    evpR, evpS, evpRStr(ncnst)
   real(r8)    u_grdPGF, v_grdPGF, PGF_u, PGF_v
   real(r8)    um, dm 
   real(r8)    thl_env_u, thl_env_d
   real(r8)    qt_env_u, qt_env_d
   real(r8)    u_env_u, u_env_d
   real(r8)    v_env_u, v_env_d
   real(r8)    ql_env_u, ql_env_d
   real(r8)    qi_env_u, qi_env_d
   real(r8)    tr_env_u, tr_env_d
   real(r8)    flxrain_top, flxsnow_top, flxrasn_top
   real(r8)    flxrain_top_in, flxsnow_top_in, flxrasn_top_in
   real(r8)    flxrain_bot_ee, flxsnow_bot_ee
   real(r8)    flxrain_bot_upee, flxsnow_bot_upee
   real(r8)    flxrain_bot_upeesm, flxsnow_bot_upeesm
   real(r8)    flxrain_bot, flxsnow_bot
   real(r8)    flxtrrs_top(ncnst)
   real(r8)    flxtrrs_bot_upee(ncnst)
   real(r8)    flxtrrs_bot_upeesm(ncnst)
   real(r8)    flxtrrs_bot(ncnst)
   real(r8)    tmp1, tmp2, tmp3, tmp4
   real(r8)    tmp_th, tmp_qv, tmp_ql, tmp_qi, tmp_qse
   real(r8)    thl_aut_tmp, qt_aut_tmp, tr_aut_tmp(ncnst)
   real(r8)    ql_aut_adi_prp, qi_aut_adi_prp 
   real(r8)    tmpx_bot, tmpx_top, tmpy_bot, tmpy_top, tmpw
   real(r8)    f_R, f_S
   real(r8)    eps_dia_L, eps_dia_I, eps_dia_V, srcg_V
   real(r8)    evprain_clr, evpsnow_clr, evplimit_clr_rain, evplimit_clr_snow, evplimit_clr
   real(r8)    subsat_clr, qw_clr, qv_clr
   real(r8)    tw
   real(r8)    precip               !  Convective rain+snow flux at surface [m/s]
   real(r8)    snow                 !  Convective snow flux at surface [m/s]
   real(r8)    evapc(mkx)           !  Evaporation rate of convective precipitation within environment [ kg/kg/s ]
   real(r8)    rliq                 !  Vertical integration of rqc in flux unit [m/s]
   real(r8)    sigma_w, sigma_thl, sigma_qt, sigma_u, sigma_v
   real(r8)    sigma_tr(ncnst)
   real(r8)    au_tent
   real(r8)    alpha_cri
   real(r8)    Ro, sigmaR
   real(r8)    Ro_min, Ro_max, sigmaR_min, sigmaR_max
   real(r8)    kw_min, kw_max
   real(r8)    cdelta_s, cdelta_w
   real(r8)    rbuoy_up, rbuoy_dn
   real(r8)    mu
   real(r8)    thv_max, thv_min, fddet  ! For continuous downdraft buoyancy sorting 

   ! ------------------------------------------------- !
   ! Variables associated with convective organization !
   ! ------------------------------------------------- !

   real(r8)    org_rad, org_ent, sum 
   real(r8)    kw, au_base, au_base_ocn, au_base_lnd, AOVTU 
   real(r8)    tmpm_array(mix,mkx), tmpi_array(mix,0:mkx)
   real(r8)    ws1, went_eff
   real(r8)    cd_thl, cd_qt, cd_u, cd_v, cd_tr(ncnst)

   ! --------------------------- !
   ! Diagnostic Output Variables !
   ! --------------------------- !

   real(r8)    kw_out(mix)

   real(r8)    sigma_w_out(mix)
   real(r8)    sigma_thl_out(mix)
   real(r8)    sigma_qt_out(mix)
   real(r8)    sigma_u_out(mix)
   real(r8)    sigma_v_out(mix)
   
   real(r8)    w_org_out(mix)
   real(r8)    thl_org_out(mix)
   real(r8)    qt_org_out(mix)
   real(r8)    u_org_out(mix)
   real(r8)    v_org_out(mix)
   
   real(r8)    tkes_out(mix)
   real(r8)    went_out(mix)
   real(r8)    went_eff_out(mix)

   ! ----------------------------------------------------------------------------------- !
   ! Rain and snow mixing ratios at the top and base interface of each layer.            !
   ! These are used for including condensate loading effect within convective downdraft. !
   ! ----------------------------------------------------------------------------------- !

   real(r8)    qrain_b, qsnow_b

   ! ------------------------------------------------------------------------------- !
   ! For mixing of convective downdraft with convective updraft or mean environment. !
   ! ------------------------------------------------------------------------------- !

   real(r8)    ssthl_tmp, ssqt_tmp, ssql_tmp, ssqi_tmp, ssqv_tmp, ssu_tmp, ssv_tmp
   real(r8)    sstr_tmp(ncnst)
   real(r8)    thl_tmp, qt_tmp, th_tmp, qv_tmp, ql_tmp, qi_tmp, qs_tmp, t_tmp, thv_tmp

   ! ----------------------------------------------------------------------------------------------- !
   ! Variables associated with unified treatment of various evaporation processes from top to bottom !
   ! ----------------------------------------------------------------------------------------------- !

   integer     ndb_evp     
   integer     ix_d_src(mkx,3)
   real(r8)    cmfdb_evp
   real(r8)    cmf_d_src(mkx,3), tr_d_src(mkx,3,ncnst) 
   real(r8)    thl_d_src(mkx,3), qt_d_src(mkx,3), ql_d_src(mkx,3), qi_d_src(mkx,3)
   real(r8)    u_d_src(mkx,3), v_d_src(mkx,3), w_d_src(mkx,3)
   real(r8)    fevp1_t_rate_src(mkx,3), fevp2_t_rate_src(mkx,3) 

   ! --------------------------- ! 
   ! For aerosol tendency output !
   ! --------------------------- !

   character(len=30)       :: varname

   !------------------------!
   !                        !
   ! Start Main Calculation !
   !                        !
   !------------------------!

   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )
   call cnst_get_ind( 'NUMLIQ', ixnumliq )
   call cnst_get_ind( 'NUMICE', ixnumice )

   ! ----------------------------------------------------------------------------- !
   ! Define index for multiple mixing environmental airs for nter = 1 and nter = 2 !
   ! ----------------------------------------------------------------------------- !
 
   if (niter .eq. 1) then
      ixi = 1
      ixf = 1
   else
      ixi = 1
      ixf = 2
   endif

   ! -------------------------------------------------------------- !
   ! Initialize formal output variables defined at all grid columns !
   ! -------------------------------------------------------------- !

   cmf_u_out(:iend,0:mkx)                            = 0._r8   
   slflx_out(:iend,0:mkx)                            = 0._r8
   qtflx_out(:iend,0:mkx)                            = 0._r8
   qvten_out(:iend,:mkx)                             = 0._r8
   qlten_out(:iend,:mkx)                             = 0._r8
   qiten_out(:iend,:mkx)                             = 0._r8
   sten_out(:iend,:mkx)                              = 0._r8
   uten_out(:iend,:mkx)                              = 0._r8
   vten_out(:iend,:mkx)                              = 0._r8 
   trten_out(:iend,:mkx,:ncnst)                      = 0._r8
   qrten_out(:iend,:mkx)                             = 0._r8
   qsten_out(:iend,:mkx)                             = 0._r8
   rqc_l_out(:iend,:mkx)                             = 0._r8
   rqc_i_out(:iend,:mkx)                             = 0._r8
   rqc_out(:iend,:mkx)                               = 0._r8
   rnc_l_out(:iend,:mkx)                             = 0._r8
   rnc_i_out(:iend,:mkx)                             = 0._r8
   rliq_out(:iend)                                   = 0._r8
   precip_out(:iend)                                 = 0._r8
   snow_out(:iend)                                   = 0._r8
   evapc_out(:iend,:mkx)                             = 0._r8
   cnt_out(:iend)                                    = real(mkx,r8)
   cnb_out(:iend)                                    = 0._r8
   am_u_out(:iend,:mkx)                              = 0._r8
   qlm_u_out(:iend,:mkx)                             = 0._r8
   qim_u_out(:iend,:mkx)                             = 0._r8
   am_d_out(:iend,:mkx)                              = 0._r8
   qlm_d_out(:iend,:mkx)                             = 0._r8
   qim_d_out(:iend,:mkx)                             = 0._r8

   cmf_det_out(:iend,:mkx)                           = 0._r8
   ql_det_out(:iend,:mkx)                            = 0._r8
   qi_det_out(:iend,:mkx)                            = 0._r8

   ! ---------------------------------------------------------------- !
   ! Initialize internal output variables defined at all grid columns !
   ! ---------------------------------------------------------------- !

   cmf_out(:iend,0:mkx)                              = 0._r8   
   uflx_out(:iend,0:mkx)                             = 0._r8 
   vflx_out(:iend,0:mkx)                             = 0._r8

   slflx_u_out(:iend,0:mkx)                          = 0._r8 
   qtflx_u_out(:iend,0:mkx)                          = 0._r8 
   uflx_u_out(:iend,0:mkx)                           = 0._r8 
   vflx_u_out(:iend,0:mkx)                           = 0._r8

   cmf_d_out(:iend,0:mkx)                            = 0._r8   
   slflx_d_out(:iend,0:mkx)                          = 0._r8 
   qtflx_d_out(:iend,0:mkx)                          = 0._r8 
   uflx_d_out(:iend,0:mkx)                           = 0._r8 
   vflx_d_out(:iend,0:mkx)                           = 0._r8

   thl_orgforce_out(:iend)                           = 0._r8 
   qt_orgforce_out(:iend)                            = 0._r8 
   u_orgforce_out(:iend)                             = 0._r8 
   v_orgforce_out(:iend)                             = 0._r8 
   tr_orgforce_out(:iend,:ncnst)                     = 0._r8 
   awk_orgforce_out(:iend)                           = 0._r8 

   ! Below block is for detailed diagnostic output

   flxrain_out(:iend,0:mkx)                          = 0._r8 
   flxsnow_out(:iend,0:mkx)                          = 0._r8 

   thl_orgforce_flx_out(:iend)                       = 0._r8 
   qt_orgforce_flx_out(:iend)                        = 0._r8 
   u_orgforce_flx_out(:iend)                         = 0._r8 
   v_orgforce_flx_out(:iend)                         = 0._r8 
   awk_orgforce_flx_out(:iend)                       = 0._r8 

   thl_orgforce_und_out(:iend)                       = 0._r8 
   qt_orgforce_und_out(:iend)                        = 0._r8 
   u_orgforce_und_out(:iend)                         = 0._r8 
   v_orgforce_und_out(:iend)                         = 0._r8 
   awk_orgforce_mix_out(:iend)                       = 0._r8 

   thl_orgforce_env_out(:iend)                       = 0._r8 
   qt_orgforce_env_out(:iend)                        = 0._r8 
   u_orgforce_env_out(:iend)                         = 0._r8 
   v_orgforce_env_out(:iend)                         = 0._r8 
   cmf_d_org_pblh_out(:iend)                         = 0._r8 

   ! Above block is for detailed diagnostic output

   taui_thl_out(:iend)                               = 0._r8 
   taui_qt_out(:iend)                                = 0._r8 
   taui_u_out(:iend)                                 = 0._r8 
   taui_v_out(:iend)                                 = 0._r8 
   taui_tr_out(:iend,:ncnst)                         = 0._r8 
   taui_awk_out(:iend)                               = 0._r8 

   del_org_out(:iend)                                = 0._r8 
   del0_org_out(:iend)                               = 0._r8 

   slten_u_out(:iend,:mkx)                           = 0._r8 
   qtten_u_out(:iend,:mkx)                           = 0._r8 
   uten_u_out(:iend,:mkx)                            = 0._r8 
   vten_u_out(:iend,:mkx)                            = 0._r8
   sten_u_out(:iend,:mkx)                            = 0._r8 
   qvten_u_out(:iend,:mkx)                           = 0._r8 
   qlten_u_out(:iend,:mkx)                           = 0._r8 
   qiten_u_out(:iend,:mkx)                           = 0._r8
   trten_u_out(:iend,:mkx,:ncnst)                    = 0._r8

   slten_d_out(:iend,:mkx)                           = 0._r8 
   qtten_d_out(:iend,:mkx)                           = 0._r8 
   uten_d_out(:iend,:mkx)                            = 0._r8 
   vten_d_out(:iend,:mkx)                            = 0._r8
   sten_d_out(:iend,:mkx)                            = 0._r8 
   qvten_d_out(:iend,:mkx)                           = 0._r8 
   qlten_d_out(:iend,:mkx)                           = 0._r8 
   qiten_d_out(:iend,:mkx)                           = 0._r8
   trten_d_out(:iend,:mkx,:ncnst)                    = 0._r8

   slten_evp_out(:iend,:mkx)                         = 0._r8 
   qtten_evp_out(:iend,:mkx)                         = 0._r8 
   uten_evp_out(:iend,:mkx)                          = 0._r8 
   vten_evp_out(:iend,:mkx)                          = 0._r8
   sten_evp_out(:iend,:mkx)                          = 0._r8 
   qvten_evp_out(:iend,:mkx)                         = 0._r8 
   qlten_evp_out(:iend,:mkx)                         = 0._r8 
   qiten_evp_out(:iend,:mkx)                         = 0._r8
   trten_evp_out(:iend,:mkx,:ncnst)                  = 0._r8

   slten_dis_out(:iend,:mkx)                         = 0._r8 
   qtten_dis_out(:iend,:mkx)                         = 0._r8 
   uten_dis_out(:iend,:mkx)                          = 0._r8 
   vten_dis_out(:iend,:mkx)                          = 0._r8
   sten_dis_out(:iend,:mkx)                          = 0._r8 
   qvten_dis_out(:iend,:mkx)                         = 0._r8 
   qlten_dis_out(:iend,:mkx)                         = 0._r8 
   qiten_dis_out(:iend,:mkx)                         = 0._r8
   trten_dis_out(:iend,:mkx,:ncnst)                  = 0._r8

   qlten_sub_out(:iend,:mkx)                         = 0._r8       
   qiten_sub_out(:iend,:mkx)                         = 0._r8       
   qlten_det_out(:iend,:mkx)                         = 0._r8       
   qiten_det_out(:iend,:mkx)                         = 0._r8       

   thl_u_out(:iend,0:mkx)                            = 0._r8
   qt_u_out(:iend,0:mkx)                             = 0._r8
   u_u_out(:iend,0:mkx)                              = 0._r8
   v_u_out(:iend,0:mkx)                              = 0._r8
   w_u_out(:iend,0:mkx)                              = 0._r8
   ql_u_out(:iend,0:mkx)                             = 0._r8
   qi_u_out(:iend,0:mkx)                             = 0._r8 
   tr_u_out(:iend,0:mkx,:ncnst)                      = 0._r8 
   wa_u_out(:iend,0:mkx)                             = 0._r8
   qla_u_out(:iend,0:mkx)                            = 0._r8
   qia_u_out(:iend,0:mkx)                            = 0._r8
   a_u_out(:iend,0:mkx)                              = 0._r8
   rad_u_out(:iend,0:mkx)                            = 0._r8
   num_u_out(:iend,0:mkx)                            = 0._r8
   gamw_u_out(:iend,0:mkx)                           = 0._r8
   thva_u_out(:iend,0:mkx)                           = 0._r8

   thl_d_out(:iend,0:mkx)                            = 0._r8
   qt_d_out(:iend,0:mkx)                             = 0._r8
   u_d_out(:iend,0:mkx)                              = 0._r8
   v_d_out(:iend,0:mkx)                              = 0._r8
   w_d_out(:iend,0:mkx)                              = 0._r8
   ql_d_out(:iend,0:mkx)                             = 0._r8
   qi_d_out(:iend,0:mkx)                             = 0._r8
   tr_d_out(:iend,0:mkx,:ncnst)                      = 0._r8 
   wa_d_out(:iend,0:mkx)                             = 0._r8
   qla_d_out(:iend,0:mkx)                            = 0._r8
   qia_d_out(:iend,0:mkx)                            = 0._r8
   a_d_out(:iend,0:mkx)                              = 0._r8

   a_p_out(:iend,0:mkx)                              = 0._r8
   am_evp_out(:iend,:mkx)                            = 0._r8
   am_pu_out(:iend,:mkx)                             = 0._r8
   x_p_out(:iend,0:mkx)                              = 0._r8
   y_p_out(:iend,0:mkx)                              = 0._r8
   x_um_out(:iend,:mkx)                              = 0._r8
   y_um_out(:iend,:mkx)                              = 0._r8

   thl_u_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   qt_u_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   u_u_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   v_u_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   w_u_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   ql_u_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   qi_u_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   tr_u_msfc_out(:iend,0:mkx,:nseg,:ncnst,:niter)    = 0._r8
   cmf_u_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   a_u_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   num_u_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   rad_u_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8

   eps0_u_msfc_out(:iend,0:mkx,:nseg,:niter)         = 0._r8
   eps_u_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   del_u_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   eeps_u_msfc_out(:iend,0:mkx,:nseg,:niter)         = 0._r8
   ddel_u_msfc_out(:iend,0:mkx,:nseg,:niter)         = 0._r8
   xc_u_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8 
   xs_u_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8 
   xemin_u_msfc_out(:iend,0:mkx,:nseg,:niter)        = 0._r8
   xemax_u_msfc_out(:iend,0:mkx,:nseg,:niter)        = 0._r8
   cridis_u_msfc_out(:iend,0:mkx,:nseg,:niter)       = 0._r8
   thvcuenv_u_msfc_out(:iend,0:mkx,:nseg,:niter)     = 0._r8
   thvegenv_u_msfc_out(:iend,0:mkx,:nseg,:niter)     = 0._r8
   thvxsenv_u_msfc_out(:iend,0:mkx,:nseg,:niter)     = 0._r8
   fmix_u_msfc_out(:iend,0:mkx,:nseg,:niter)         = 0._r8
   cmfumix_u_msfc_out(:iend,0:mkx,:nseg,:niter)      = 0._r8

   thl_d_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   qt_d_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   u_d_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   v_d_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   w_d_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8
   ql_d_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   qi_d_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   tr_d_msfc_out(:iend,0:mkx,:nseg,:ncnst,:niter)    = 0._r8
   wa_d_msfc_out(:iend,0:mkx,:nseg,:niter)           = 0._r8
   qla_d_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   qia_d_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   cmf_d_msfc_out(:iend,0:mkx,:nseg,:niter)          = 0._r8
   a_d_msfc_out(:iend,0:mkx,:nseg,:niter)            = 0._r8  
 
   ptop_msfc_out(:iend,:nseg,:niter)                 = 0._r8
   ztop_msfc_out(:iend,:nseg,:niter)                 = 0._r8

   thv_b_out(:iend,0:mkx)                            = 0._r8
   thv_t_out(:iend,0:mkx)                            = 0._r8
   thv_mt_out(:iend,0:mkx)                           = 0._r8
   thv_min_out(:iend,0:mkx)                          = 0._r8

   kw_out(:iend)                                     = 0._r8  

   sigma_w_out(:iend)                                = 0._r8  
   sigma_thl_out(:iend)                              = 0._r8
   sigma_qt_out(:iend)                               = 0._r8
   sigma_u_out(:iend)                                = 0._r8
   sigma_v_out(:iend)                                = 0._r8
   w_org_out(:iend)                                  = 0._r8
   thl_org_out(:iend)                                = 0._r8
   qt_org_out(:iend)                                 = 0._r8
   u_org_out(:iend)                                  = 0._r8
   v_org_out(:iend)                                  = 0._r8
   tkes_out(:iend)                                   = 0._r8
   went_out(:iend)                                   = 0._r8
   went_eff_out(:iend)                               = 0._r8

   ! ------------------------------------------------------ !
   ! Initialize other variables defined at all grid columns !
   ! ------------------------------------------------------ !

   !---------------------------------------------------------!
   !                                                         !
   ! Start the big i loop where i is a horozontal grid index !
   !                                                         !
   !---------------------------------------------------------!

   do i = 1, iend                                 

      ! ------------------------------------------------------------------------------------------------------- !
      ! Mar.27.2012. Initialize dissipation heating variables since dissipation heating is computed after doing ! 
      !              ensemble-mean average of iter = 1, niter = 2.                                              !
      ! ------------------------------------------------------------------------------------------------------- !

      slten_dis(:mkx)                                 = 0._r8
      qtten_dis(:mkx)                                 = 0._r8
      uten_dis(:mkx)                                  = 0._r8
      vten_dis(:mkx)                                  = 0._r8
      sten_dis(:mkx)                                  = 0._r8
      qvten_dis(:mkx)                                 = 0._r8
      qlten_dis(:mkx)                                 = 0._r8
      qiten_dis(:mkx)                                 = 0._r8
      trten_dis(:mkx,:ncnst)                          = 0._r8

      uf(:mkx)                                        = 0._r8
      vf(:mkx)                                        = 0._r8

      uflx(0:mkx)                                     = 0._r8
      vflx(0:mkx)                                     = 0._r8

      ! ----------------------------------------------------------------------------------------------- !
      ! Initialize variables associated with mixing with multiple mixing environmental airs ( '_mxen' ) !
      ! ----------------------------------------------------------------------------------------------- !

      ktop_mxen(:niter)                               = 0
      cush_mxen(:niter)                               = 0._r8
      cushavg_mxen(:niter)                            = 0._r8

      cu_cmfum_mxen(:mkx,:niter)                      = 0._r8
      cu_cmfr_mxen(:mkx,:niter)                       = 0._r8
      cu_thlr_mxen(:mkx,:niter)                       = 0._r8
      cu_qtr_mxen(:mkx,:niter)                        = 0._r8
      cu_ur_mxen(:mkx,:niter)                         = 0._r8
      cu_vr_mxen(:mkx,:niter)                         = 0._r8
      cu_qlr_mxen(:mkx,:niter)                        = 0._r8
      cu_qir_mxen(:mkx,:niter)                        = 0._r8
      cu_trr_mxen(:mkx,:ncnst,:niter)                 = 0._r8
      cu_cmfrd_mxen(:mkx,:niter)                      = 0._r8
      cu_thlrd_mxen(:mkx,:niter)                      = 0._r8
      cu_qtrd_mxen(:mkx,:niter)                       = 0._r8
      cu_urd_mxen(:mkx,:niter)                        = 0._r8  
      cu_vrd_mxen(:mkx,:niter)                        = 0._r8
      cu_qlrd_mxen(:mkx,:niter)                       = 0._r8
      cu_qird_mxen(:mkx,:niter)                       = 0._r8
      cu_trrd_mxen(:mkx,:ncnst,:niter)                = 0._r8

      cmf_u_mxen(0:mkx,:niter)                        = 0._r8
      cmf_d_mxen(0:mkx,:niter)                        = 0._r8
      slflx_u_mxen(0:mkx,:niter)                      = 0._r8
      slflx_d_mxen(0:mkx,:niter)                      = 0._r8
      qtflx_u_mxen(0:mkx,:niter)                      = 0._r8
      qtflx_d_mxen(0:mkx,:niter)                      = 0._r8
      uflx_u_mxen(0:mkx,:niter)                       = 0._r8
      uflx_d_mxen(0:mkx,:niter)                       = 0._r8
      vflx_u_mxen(0:mkx,:niter)                       = 0._r8
      vflx_d_mxen(0:mkx,:niter)                       = 0._r8

      flxrain_u_mxen(0:mkx,:niter)                    = 0._r8
      flxsnow_u_mxen(0:mkx,:niter)                    = 0._r8

      thl_orgforce_mxen(:niter)                       = 0._r8 
      qt_orgforce_mxen(:niter)                        = 0._r8 
      u_orgforce_mxen(:niter)                         = 0._r8 
      v_orgforce_mxen(:niter)                         = 0._r8 
      tr_orgforce_mxen(:ncnst,:niter)                 = 0._r8 
      awk_orgforce_mxen(:niter)                       = 0._r8 

      ! Below block is for detailed diagnostic output

      thl_orgforce_flx_mxen(:niter)                   = 0._r8 
      qt_orgforce_flx_mxen(:niter)                    = 0._r8 
      u_orgforce_flx_mxen(:niter)                     = 0._r8 
      v_orgforce_flx_mxen(:niter)                     = 0._r8 
      awk_orgforce_flx_mxen(:niter)                   = 0._r8 

      thl_orgforce_und_mxen(:niter)                   = 0._r8 
      qt_orgforce_und_mxen(:niter)                    = 0._r8 
      u_orgforce_und_mxen(:niter)                     = 0._r8 
      v_orgforce_und_mxen(:niter)                     = 0._r8 
      awk_orgforce_mix_mxen(:niter)                   = 0._r8 

      thl_orgforce_env_mxen(:niter)                   = 0._r8 
      qt_orgforce_env_mxen(:niter)                    = 0._r8 
      u_orgforce_env_mxen(:niter)                     = 0._r8 
      v_orgforce_env_mxen(:niter)                     = 0._r8 
      cmf_d_org_pblh_mxen(:niter)                     = 0._r8 

      ! Above block is for detailed diagnostic output

      taui_thl_mxen(:niter)                           = 0._r8 
      taui_qt_mxen(:niter)                            = 0._r8 
      taui_u_mxen(:niter)                             = 0._r8 
      taui_v_mxen(:niter)                             = 0._r8 
      taui_tr_mxen(:ncnst,:niter)                     = 0._r8 
      taui_awk_mxen(:niter)                           = 0._r8 

      del_org_mxen(:niter)                            = 0._r8 
      del0_org_mxen(:niter)                           = 0._r8 

      qvten_mxen(:mkx,:niter)                         = 0._r8
      qlten_mxen(:mkx,:niter)                         = 0._r8
      qiten_mxen(:mkx,:niter)                         = 0._r8
      trten_mxen(:mkx,:ncnst,:niter)                  = 0._r8 
      sten_mxen(:mkx,:niter)                          = 0._r8
      uten_mxen(:mkx,:niter)                          = 0._r8
      vten_mxen(:mkx,:niter)                          = 0._r8
      qrten_mxen(:mkx,:niter)                         = 0._r8
      qsten_mxen(:mkx,:niter)                         = 0._r8

      rqc_l_mxen(:mkx,:niter)                         = 0._r8
      rqc_i_mxen(:mkx,:niter)                         = 0._r8
      rqc_mxen(:mkx,:niter)                           = 0._r8
      rnc_l_mxen(:mkx,:niter)                         = 0._r8
      rnc_i_mxen(:mkx,:niter)                         = 0._r8

      cmf_det_mxen(:mkx,:niter)                       = 0._r8
      ql_det_mxen(:mkx,:niter)                        = 0._r8
      qi_det_mxen(:mkx,:niter)                        = 0._r8

      evapc_mxen(:mkx,:niter)                         = 0._r8

      am_u_mxen(:mkx,:niter)                          = 0._r8
      qlm_u_mxen(:mkx,:niter)                         = 0._r8
      qim_u_mxen(:mkx,:niter)                         = 0._r8

      am_d_mxen(:mkx,:niter)                          = 0._r8
      qlm_d_mxen(:mkx,:niter)                         = 0._r8
      qim_d_mxen(:mkx,:niter)                         = 0._r8

      rliq_mxen(:niter)                               = 0._r8
      precip_mxen(:niter)                             = 0._r8
      snow_mxen(:niter)                               = 0._r8

      cnt_mxen(:niter)                                = 0._r8
      cnb_mxen(:niter)                                = 0._r8

      slten_u_mxen(:mkx,:niter)                       = 0._r8
      qtten_u_mxen(:mkx,:niter)                       = 0._r8
      uten_u_mxen(:mkx,:niter)                        = 0._r8
      vten_u_mxen(:mkx,:niter)                        = 0._r8
      sten_u_mxen(:mkx,:niter)                        = 0._r8
      qvten_u_mxen(:mkx,:niter)                       = 0._r8
      qlten_u_mxen(:mkx,:niter)                       = 0._r8
      qiten_u_mxen(:mkx,:niter)                       = 0._r8
      trten_u_mxen(:mkx,:ncnst,:niter)                = 0._r8

      slten_d_mxen(:mkx,:niter)                       = 0._r8
      qtten_d_mxen(:mkx,:niter)                       = 0._r8
      uten_d_mxen(:mkx,:niter)                        = 0._r8
      vten_d_mxen(:mkx,:niter)                        = 0._r8
      sten_d_mxen(:mkx,:niter)                        = 0._r8
      qvten_d_mxen(:mkx,:niter)                       = 0._r8
      qlten_d_mxen(:mkx,:niter)                       = 0._r8
      qiten_d_mxen(:mkx,:niter)                       = 0._r8
      trten_d_mxen(:mkx,:ncnst,:niter)                = 0._r8

      slten_evp_mxen(:mkx,:niter)                     = 0._r8
      qtten_evp_mxen(:mkx,:niter)                     = 0._r8
      uten_evp_mxen(:mkx,:niter)                      = 0._r8
      vten_evp_mxen(:mkx,:niter)                      = 0._r8
      sten_evp_mxen(:mkx,:niter)                      = 0._r8
      qvten_evp_mxen(:mkx,:niter)                     = 0._r8
      qlten_evp_mxen(:mkx,:niter)                     = 0._r8
      qiten_evp_mxen(:mkx,:niter)                     = 0._r8
      trten_evp_mxen(:mkx,:ncnst,:niter)              = 0._r8

      qlten_sub_mxen(:mkx,:niter)                     = 0._r8
      qiten_sub_mxen(:mkx,:niter)                     = 0._r8

      qlten_det_mxen(:mkx,:niter)                     = 0._r8
      qiten_det_mxen(:mkx,:niter)                     = 0._r8

      thl_u_mxen(0:mkx,:niter)                        = 0._r8
      qt_u_mxen(0:mkx,:niter)                         = 0._r8
      u_u_mxen(0:mkx,:niter)                          = 0._r8
      v_u_mxen(0:mkx,:niter)                          = 0._r8
      w_u_mxen(0:mkx,:niter)                          = 0._r8
      ql_u_mxen(0:mkx,:niter)                         = 0._r8
      qi_u_mxen(0:mkx,:niter)                         = 0._r8
      tr_u_mxen(0:mkx,:ncnst,:niter)                  = 0._r8
      a_u_mxen(0:mkx,:niter)                          = 0._r8
      num_u_mxen(0:mkx,:niter)                        = 0._r8
      wa_u_mxen(0:mkx,:niter)                         = 0._r8
      qla_u_mxen(0:mkx,:niter)                        = 0._r8
      qia_u_mxen(0:mkx,:niter)                        = 0._r8
      rad_u_mxen(0:mkx,:niter)                        = 0._r8
      thva_u_mxen(0:mkx,:niter)                       = 0._r8

      a_p_mxen(0:mkx,:niter)                          = 0._r8
      am_evp_mxen(:mkx,:niter)                        = 0._r8
      am_pu_mxen(:mkx,:niter)                         = 0._r8
      x_p_mxen(0:mkx,:niter)                          = 0._r8
      y_p_mxen(0:mkx,:niter)                          = 0._r8
      x_um_mxen(:mkx,:niter)                          = 0._r8
      y_um_mxen(:mkx,:niter)                          = 0._r8

      thl_d_mxen(0:mkx,:niter)                        = 0._r8
      qt_d_mxen(0:mkx,:niter)                         = 0._r8
      u_d_mxen(0:mkx,:niter)                          = 0._r8
      v_d_mxen(0:mkx,:niter)                          = 0._r8
      w_d_mxen(0:mkx,:niter)                          = 0._r8
      ql_d_mxen(0:mkx,:niter)                         = 0._r8
      qi_d_mxen(0:mkx,:niter)                         = 0._r8
      tr_d_mxen(0:mkx,:ncnst,:niter)                  = 0._r8
      a_d_mxen(0:mkx,:niter)                          = 0._r8
      wa_d_mxen(0:mkx,:niter)                         = 0._r8
      qla_d_mxen(0:mkx,:niter)                        = 0._r8
      qia_d_mxen(0:mkx,:niter)                        = 0._r8

      thl_u_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      qt_u_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      u_u_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      v_u_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      w_u_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      ql_u_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      qi_u_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      tr_u_msfc_mxen(0:mkx,:nseg,:ncnst,:niter)       = 0._r8
      cmf_u_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      a_u_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      num_u_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      rad_u_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8

      eps0_u_msfc_mxen(0:mkx,:nseg,:niter)            = 0._r8
      eps_u_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      del_u_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      eeps_u_msfc_mxen(0:mkx,:nseg,:niter)            = 0._r8
      ddel_u_msfc_mxen(0:mkx,:nseg,:niter)            = 0._r8
      xc_u_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      xs_u_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      xemin_u_msfc_mxen(0:mkx,:nseg,:niter)           = 0._r8
      xemax_u_msfc_mxen(0:mkx,:nseg,:niter)           = 0._r8
      cridis_u_msfc_mxen(0:mkx,:nseg,:niter)          = 0._r8
      thvcuenv_u_msfc_mxen(0:mkx,:nseg,:niter)        = 0._r8
      thvegenv_u_msfc_mxen(0:mkx,:nseg,:niter)        = 0._r8
      thvxsenv_u_msfc_mxen(0:mkx,:nseg,:niter)        = 0._r8
      fmix_u_msfc_mxen(0:mkx,:nseg,:niter)            = 0._r8
      cmfumix_u_msfc_mxen(0:mkx,:nseg,:niter)         = 0._r8

      thl_d_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      qt_d_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      u_d_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      v_d_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      w_d_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      ql_d_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      qi_d_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      tr_d_msfc_mxen(0:mkx,:nseg,:ncnst,:niter)       = 0._r8
      cmf_d_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      a_d_msfc_mxen(0:mkx,:nseg,:niter)               = 0._r8
      wa_d_msfc_mxen(0:mkx,:nseg,:niter)              = 0._r8
      qla_d_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8
      qia_d_msfc_mxen(0:mkx,:nseg,:niter)             = 0._r8

      ptop_msfc_mxen(:nseg,:niter)                    = 0._r8
      ztop_msfc_mxen(:nseg,:niter)                    = 0._r8

      ! --------------------------------------------- !
      ! Define 1D input variables at each grid point  !
      ! Interface index from sfc : k = 0,1,2,...,mkx  !
      ! Mid-point index from sfc : k = 1,2,3,...,mkx  !
      ! --------------------------------------------- !

      kpblh                =                         kpblh_in(i)
      pblh                 =                          pblh_in(i)
      went                 =                          went_in(i)
      qflx                 =                          qflx_in(i)      
      shflx                =                         shflx_in(i)      
      taux                 =                          taux_in(i)      
      tauy                 =                          tauy_in(i)      
      aflx(:ncnst)         =                   aflx_in(i,:ncnst)
      landfrac             =                      landfrac_in(i)      
      sgh30                =                         sgh30_in(i)      

      ! ------------------------------------------------------------------------------------------------------ !
      ! Aug.31.2011. In order to simplify code later, impose the condition that kpblh >=2 and                  !
      !              corresponding pblh >= zs0_in(i,1). This condition is always satisfied when                !
      !              PBL is CL, but when PBL is STL, this is imposing conservative constraint.                 !
      !              since convective is likely not to very in-active in the STL, this is not an issue at all. !
      !              Even when convection is fired in STL, there is no problem at all.                         !
      ! Sep.09.2011. Add 'kpblhm', 'pblhz', 'pblhm' since these are used very frequently.                      !
      !              All the variables are replaced by these ( e.g., kpblh - 1 = kpblhm ).                     ! 
      !              Note that 'pblhz > 0' and 'pblhp > 0'.                                                    !
      ! ------------------------------------------------------------------------------------------------------ !

      kpblh               =                      max( kpblh, 2 )
      kpblhm              =                            kpblh - 1 
      pblh                =             max( pblh, zs0_in(i,1) )
      pblhz               =       zs0_in(i,kpblhm) - zs0_in(i,0)
      pblhp               =       ps0_in(i,0) - ps0_in(i,kpblhm)

      ! ---------------------------------------------------------------------------------------------- !    
      ! May.21.2011. Is it better to change the minimums of 'cush' and 'cushavg' from 1000 to 'pblh' ? !  
      !              I should think about this later.                                                  ! 
      ! May.21.2011. For full consistency with the other parts of the code, it is good to use 'pblh'   !
      !              instead of '1000._r8' as the minimum value of 'cush, cushavg'. So, I modified     !
      !              this.                                                                             !
      ! Aug.31.2011. Add delta_thl(qt,u,v,thv,tr)_PBL fields. I may need to impose a reasonable upper  !
      !              and lower limits on these excessive variables.                                    !  
      ! Sep.09.2011. I don't prognose 'thv' anymore - it is now diagnostically computed from the       !
      !              prognosed 'thl,qt' which imposes a full consistency as well as removing unfair    !
      !              assumption that both downdraft and environment are unsaturated for computing      !
      !              buoyancy flux.                                                                    !  
      ! ---------------------------------------------------------------------------------------------- !

      cush                     =                       cush_inout(i)
      cush                     =                  max( cush, pblhz )
      cushavg                  =                    cushavg_inout(i)
      cushavg                  =               max( cushavg, pblhz )
      cuorg                    =                      cuorg_inout(i)
      cuorg                    =   max( 0._r8, min( 1._r8, cuorg ) )

      ! ---------------------------------------------------------------------------------------- ! 
      ! Jun.07.2012. Include advection of horizontal heterogeneity associated with               !
      !              convective organization.                                                    !
      !              For the time being, these are treated only for                              !
      !              awk_PBL, delta_thl_PBL, delta_qt_PBL, delta_u_PBL, and delta_v_PBL.         !
      !              That is, delta_tr_PBL(:ncnst) are neglected to save computation time, which !
      !              should be included in future. Ideally, cush and cushavg should be treated   !
      !              in the same way.                                                            !  
      !              Note that a constant offset of -100. is extracted when retriving 'delta_xx' !
      !              values from tracer arrays.                                                  !
      !              Advection can cause 'awk_PBL_raw >  awk_PBL_max = 1._r8 - au_max'.          ! 
      !              In this case, re-set awk_PBL_raw to awk_PBL_max.                            !  
      ! ---------------------------------------------------------------------------------------- !

      if( iorg_adv ) then
         call cnst_get_ind( 'ORGawk', i_awk )
         call cnst_get_ind( 'ORGthl', i_thl )
         call cnst_get_ind( 'ORGqto',  i_qt )
         call cnst_get_ind( 'ORGuoo',   i_u )
         call cnst_get_ind( 'ORGvoo',   i_v )

         if( get_nstep() .eq. 0 ) then
            awk_PBL_raw       = 0._r8
            delta_thl_PBL_raw = 0._r8
            delta_qt_PBL_raw  = 0._r8
            delta_u_PBL_raw   = 0._r8
            delta_v_PBL_raw   = 0._r8
         else 
            tmp1                     = 0._r8
            awk_PBL_raw              = 0._r8
            delta_thl_PBL_raw        = 0._r8
            delta_qt_PBL_raw         = 0._r8
            delta_u_PBL_raw          = 0._r8
            delta_v_PBL_raw          = 0._r8
            do k = 1, kpblhm  ! Here, 'k' is a layer index.
               tmp1              =              tmp1 +                                                       dp0_in(i,k)
               awk_PBL_raw       =       awk_PBL_raw + min( tr0_in(i,k,i_awk), 1._r8 - au_max - 1.e-5_r8 ) * dp0_in(i,k)
! JHYoon : fixing advection problem of convective organization terms in SE
! (suggested by Sungsu Park)
!              delta_thl_PBL_raw = delta_thl_PBL_raw +    ( tr0_in(i,k,i_thl) - 100._r8 )                  * dp0_in(i,k)
!              delta_qt_PBL_raw  =  delta_qt_PBL_raw +    ( tr0_in(i,k, i_qt) - 100._r8 )                  * dp0_in(i,k)
!              delta_u_PBL_raw   =   delta_u_PBL_raw +    ( tr0_in(i,k,  i_u) - 100._r8 )                  * dp0_in(i,k)
!              delta_v_PBL_raw   =   delta_v_PBL_raw +    ( tr0_in(i,k,  i_v) - 100._r8 )                  * dp0_in(i,k)
               delta_thl_PBL_raw = delta_thl_PBL_raw +    ( tr0_in(i,k,i_thl) - 10._r8 )                  * dp0_in(i,k)
               delta_qt_PBL_raw  =  delta_qt_PBL_raw +    ( tr0_in(i,k, i_qt) - 0.01_r8 )                  * dp0_in(i,k)
               delta_u_PBL_raw   =   delta_u_PBL_raw +    ( tr0_in(i,k,  i_u) - 10._r8 )                  * dp0_in(i,k)
               delta_v_PBL_raw   =   delta_v_PBL_raw +    ( tr0_in(i,k,  i_v) - 10._r8 )                  * dp0_in(i,k)
! JHYoon
            enddo
            awk_PBL_raw       =       awk_PBL_raw / tmp1
            delta_thl_PBL_raw = delta_thl_PBL_raw / tmp1
            delta_qt_PBL_raw  =  delta_qt_PBL_raw / tmp1
            delta_u_PBL_raw   =   delta_u_PBL_raw / tmp1
            delta_v_PBL_raw   =   delta_v_PBL_raw / tmp1
         endif
      else
         awk_PBL_raw              =                    awk_PBL_inout(i)
         delta_thl_PBL_raw        =              delta_thl_PBL_inout(i)
         delta_qt_PBL_raw         =               delta_qt_PBL_inout(i)
         delta_u_PBL_raw          =                delta_u_PBL_inout(i)
         delta_v_PBL_raw          =                delta_v_PBL_inout(i)
      endif

      delta_tr_PBL_raw(:ncnst) =        delta_tr_PBL_inout(i,:ncnst)

      ! ----------------------------------------------------------------------------------- !
      ! Impose consistency between the input wake area and perturbations                    !
      ! This is a just minimal constraint not a sufficient one.                             !
      ! While this is a minimal constraint, we include additional consistncy constraint     !
      ! later in computing 'delta_thl_PBL' from  'delta_thl_PBL_raw' using tmp3 and tmp4    !
      ! later even though that is not perfrect either.                                      !
      ! If there is more consistent way, I should try to find that.                         !
      ! ----------------------------------------------------------------------------------- !

      if( awk_PBL_raw .lt. 1.e-5_r8 ) then 
         awk_PBL_raw              = 0._r8
         delta_thl_PBL_raw        = 0._r8
         delta_qt_PBL_raw         = 0._r8
         delta_u_PBL_raw          = 0._r8
         delta_v_PBL_raw          = 0._r8
         delta_tr_PBL_raw(:ncnst) = 0._r8
      endif

      cu_cmfr(:mkx)            =               cu_cmfr_inout(i,:mkx)
      cu_thlr(:mkx)            =               cu_thlr_inout(i,:mkx)
      cu_qtr(:mkx)             =                cu_qtr_inout(i,:mkx)
      cu_ur(:mkx)              =                 cu_ur_inout(i,:mkx)
      cu_vr(:mkx)              =                 cu_vr_inout(i,:mkx)
      cu_qlr(:mkx)             =                cu_qlr_inout(i,:mkx)
      cu_qir(:mkx)             =                cu_qir_inout(i,:mkx)
      cu_trr(:mkx,:ncnst)      =         cu_trr_inout(i,:mkx,:ncnst)

      ! ---------------------------------------- !
      ! Local environmental mean state variables !
      ! ---------------------------------------- !

      ps0(0:mkx)          =       ps0_in(i,0:mkx)
      zs0(0:mkx)          =       zs0_in(i,0:mkx)
      p0(:mkx)            =         p0_in(i,:mkx)
      z0(:mkx)            =         z0_in(i,:mkx)
      dp0(:mkx)           =        dp0_in(i,:mkx) 
      dpdry0(:mkx)        =     dpdry0_in(i,:mkx) 
      u0(:mkx)            =         u0_in(i,:mkx)
      v0(:mkx)            =         v0_in(i,:mkx)
      qv0(:mkx)           =        qv0_in(i,:mkx)
      ql0(:mkx)           =        ql0_in(i,:mkx)
      qi0(:mkx)           =        qi0_in(i,:mkx)
      do mt = 1, ncnst    
         tr0(:mkx,mt)     =     tr0_in(i,:mkx,mt)
      enddo
      t0(:mkx)            =         t0_in(i,:mkx)
      ast0(:mkx)          =       ast0_in(i,:mkx)
      tke0(0:mkx)         =      tke0_in(i,0:mkx)
      bprod0(0:mkx)       =    bprod0_in(i,0:mkx)

      ! Aug.08.2013. Evaporation of stratiform precipitation
      am_evp_st(:mkx)     =  am_evp_st_in(i,:mkx)
      evprain_st(:mkx)    = evprain_st_in(i,:mkx)
      evpsnow_st(:mkx)    = evpsnow_st_in(i,:mkx)

      ! --------------------------------------------------------- !
      ! Compute other basic thermodynamic variables directly from ! 
      ! the input variables at each grid point                    !
      ! --------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      ! Nov.30.2012. Layer thickness depending on dry or moist tracers. !
      ! --------------------------------------------------------------- !

      do mt = 1, ncnst
         if( cnst_get_type_byind(mt) .eq. 'wet' ) then
            dptr0(:mkx,mt) = dp0(:mkx)
         else
            dptr0(:mkx,mt) = dpdry0(:mkx)
         endif
      enddo

      ! --------------------------------------------- !
      ! Compute conservative scalars at the mid-point !
      ! --------------------------------------------- !     

      exn0(:mkx)   = ( p0(:mkx) / p00 )**rovcp
      exns0(0:mkx) = ( ps0(0:mkx) / p00 )**rovcp
      qt0(:mkx)    = ( qv0(:mkx) + ql0(:mkx) + qi0(:mkx) )
      thl0(:mkx)   = ( t0(:mkx) - xlv * ql0(:mkx) / cp - xls * qi0(:mkx) / cp ) / exn0(:mkx)
      rho0(:mkx)   = ( p0(:mkx) ) / ( r * t0(:mkx) * ( 1._r8 + zvir * qv0(:mkx) - ql0(:mkx) - qi0(:mkx) ) )
      do k = 1, mkx
         call conden( p0(k), thl0(k), qt0(k), th, qv, ql, qi, qse, id_check )
         thv0(k) = th * ( 1._r8 + zvir * qv - ql - qi )
         rh0(k)  = max( 0._r8, min( 1._r8, qv / max( nonzero, qse ) ) )
      enddo
      do k = 1, mkx
         dz0(k)  = zs0(k)  - zs0(k-1)
         if( k .eq. mkx ) then
            dps0(k) = p0(k) - ps0(k)
         else
            dps0(k) = p0(k)   - p0(k+1)
         endif
      end do
      dps0(0) = ps0(1) - p0(1)

      ! Dec.17.2012. Restore below block in addition to the computation of the other thermodynamic variables since
      !              these PBL-averaged variables are critically used for the revised 'bulk' computation of
      !              organized flux at the PBL top interface. 

      tmp1            = 0._r8
      qt0PBL          = 0._r8
      thl0PBL         = 0._r8
      u0PBL           = 0._r8
      v0PBL           = 0._r8
      tr0PBL(1:ncnst) = 0._r8 
      do k = 1, kpblhm  ! Here, 'k' is a layer index.
         dpi     = ps0(k-1) - ps0(k)
         tmp1    = tmp1     + dpi
         qt0PBL  = qt0PBL   + dpi*qt0(k)
         thl0PBL = thl0PBL  + dpi*thl0(k)
         u0PBL   = u0PBL    + dpi*u0(k)
         v0PBL   = v0PBL    + dpi*v0(k)
         do mt = 1, ncnst
            tr0PBL(mt) = tr0PBL(mt) + dpi*tr0(k,mt)
         enddo
      end do
      qt0PBL   = qt0PBL  / tmp1
      thl0PBL  = thl0PBL / tmp1
      u0PBL    = u0PBL   / tmp1
      v0PBL    = v0PBL   / tmp1
      do mt = 1, ncnst
         tr0PBL(mt) = tr0PBL(mt) / tmp1
      enddo

      ! ---------------------------------------------------------------------------- !
      ! Dec.20.2012. Compute reconstructed effective entrainment rate at the PBL top !
      !              interface, 'we_eff' [m/s] for use in computing damping time     ! 
      !              scale of organized flow within PBL.                             !
      !              Use the average of two constructed from 'qt0' and 'thl0'.       !
      !              Below computation may introduce a sensitivity to vertical       !
      !              resolution but hopely that effect is likely small.              ! 
      ! ---------------------------------------------------------------------------- !

      tmp1 = abs( thl0(kpblh) - thl0(kpblhm) ) / max( abs( thl0(kpblh) - thl0PBL ), 1._r8    )  + &
             abs(  qt0(kpblh) -  qt0(kpblhm) ) / max( abs(  qt0(kpblh) -  qt0PBL ), 1.e-3_r8 )
      went_eff = 0.5_r8 * tmp1 * went
      went_eff = min( max( 0._r8, went_eff ), 1._r8 )
      went_eff_out(i) = went_eff
      went_out(i)     = went

      ! ------------------------------------------------------------------------- !
      ! Compute in-layer slopes of conservative scalars                           !
      ! Dimension of slope is implicitly (1:mkx).                                 !
      ! Unit is [K/Pa] (for thl0) with a negative when thl increases vertically.  !
      ! It turns out that using non-zero slope induces inversion of 'thv0' across !
      ! the model interface, distorting model performance, including too much     !
      ! convective downdraft from updraft buoyancy sorting. Thus, it is much      !
      ! better to use the zero slope, which almost always guarantees stratified   !
      ! input environmental profile except in the lowest model layer.             ! 
      ! ------------------------------------------------------------------------- !

      ssthl0 = slope( mkx, thl0, p0 ) 
      ssqt0  = slope( mkx, qt0 , p0 )
      ssu0   = slope( mkx, u0  , p0 )
      ssv0   = slope( mkx, v0  , p0 )
      do mt = 1, ncnst
         sstr0(:mkx,mt) = slope( mkx, tr0(:mkx,mt), p0 )
      enddo
      if( islope_on_thlqttr .eq. 0 ) then           
         ssthl0(:mkx) = 0._r8
         ssqt0(:mkx)  = 0._r8
         do mt = 1, ncnst
            sstr0(:mkx,mt) = 0._r8
         enddo
      endif
      if( islope_on_uv .eq. 0 ) then           
         ssu0(:mkx)   = 0._r8
         ssv0(:mkx)   = 0._r8
      endif

      ! ------------------------------------------------------------- !
      ! Compute 'qt,thl,u,v,thv,thvl' at the top/bottom interfaces    !
      ! Note 'thv,thvl' are consistently computed from the top/bottom !
      ! interface values of 'thl,qt'.                                 !
      ! ------------------------------------------------------------- !

      do k = 1, mkx

         km = k - 1

         thl0bot(k)  = thl0(k) + ssthl0(k) * ( ps0(km) - p0(k) )
         qt0bot(k)   = qt0(k)  + ssqt0(k)  * ( ps0(km) - p0(k) )

         qt0bot(k)   = max( qt0bot(k), qmin(1) )

         u0bot(k)    = u0(k)   + ssu0(k)   * ( ps0(km) - p0(k) )
         v0bot(k)    = v0(k)   + ssv0(k)   * ( ps0(km) - p0(k) )
         do mt = 1, ncnst
            tr0bot(k,mt) = tr0(k,mt) + sstr0(k,mt) * ( ps0(km) - p0(k) )

            tr0bot(k,mt) = max( tr0bot(k,mt), qmin(mt) )

         enddo
         call conden( ps0(km), thl0bot(k), qt0bot(k), th, qv, ql, qi, qse, id_check )
         thvl0bot(k) = thl0bot(k) * ( 1._r8 + zvir * qt0bot(k) )
         thv0bot(k)  = th * ( 1._r8 + zvir * qv - ql - qi )
         ql0bot(k)   = ql
         qi0bot(k)   = qi
         if( islope_on_thlqttr .eq. 0 ) then
            thv0bot(k) = thv0(k)
            ql0bot(k)  = ql0(k)
            qi0bot(k)  = qi0(k)
         endif
         rho0bot(k)  = ps0(km) / ( r * thv0bot(k) * exns0(km) )
         rh0bot(k) = max( 0._r8, min( 1._r8, qv / max( nonzero, qse ) ) )

         thl0top(k)  = thl0(k) + ssthl0(k) * ( ps0(k) - p0(k) )
         qt0top(k)   = qt0(k)  + ssqt0(k)  * ( ps0(k) - p0(k) )

         qt0top(k)   = max( qt0top(k), qmin(1) )

         u0top(k)    = u0(k)   + ssu0(k)   * ( ps0(k) - p0(k) )
         v0top(k)    = v0(k)   + ssv0(k)   * ( ps0(k) - p0(k) )
         do mt = 1, ncnst
            tr0top(k,mt) = tr0(k,mt) + sstr0(k,mt) * ( ps0(k) - p0(k) )

            tr0top(k,mt) = max( tr0top(k,mt), qmin(mt) )

         enddo
         call conden( ps0(k), thl0top(k), qt0top(k), th, qv, ql, qi, qse, id_check )
         thvl0top(k) = thl0top(k) * ( 1._r8 + zvir * qt0top(k) )
         thv0top(k)  = th * ( 1._r8 + zvir * qv - ql - qi )
         ql0top(k)   = ql
         qi0top(k)   = qi
         if( islope_on_thlqttr .eq. 0 ) then
            thv0top(k) = thv0(k)
            ql0top(k)  = ql0(k)
            qi0top(k)  = qi0(k)
         endif
         rho0top(k)  = ps0(k) / ( r * thv0top(k) * exns0(k) )

         ssql0(k)    = ( ql0top(k)   -   ql0bot(k) ) / ( ps0(k) - ps0(km) )
         ssqi0(k)    = ( qi0top(k)   -   qi0bot(k) ) / ( ps0(k) - ps0(km) )

      end do   ! k = 1, mkx

      ! ---------------------------------------------------------------------------- !  
      ! Compute 'tke1' and 'wstar1' in the lowest model layer                        !
      ! Below compute average 'tke' and 'wstar' from the sfc                         !
      ! to the specified top interface, kc.                                          !
      ! Also define 'wstar2' only using 'surface buoyancy production' ( bprod0(0) )  !
      ! and 'PBLH'. This 'wstar2' seems to be the most reasonable choice of          !
      ! velocity scale to compute 'sigma_w'. This 'wstar2' is also independent of    !
      ! vertical resolution of GCM model grid.                                       !
      ! ---------------------------------------------------------------------------- !

      tmp1   = 0._r8
      tmp2   = 0._r8
      tke1   = 0._r8
      wstar1 = 0._r8
      kc = 1               ! Top interface of the average domain.
      do k = 0, kc         ! Here, 'k' is an interfacial layer index.
         if( k .eq. 0 ) then
            dpi = ps0(0) - p0(1)
            dzi = z0(1) - zs0(0)
         elseif( k .eq. kc ) then
            dpi = p0(kc) - ps0(kc)
            dzi = zs0(kc) - z0(kc)
         else
            dpi = p0(k) - p0(k+1)
            dzi = z0(k+1) - z0(k)
         endif
         tmp1 = tmp1 + dpi
         tmp2 = tmp2 + dzi
         tke1 = tke1 + dpi*tke0(k)
         wstar1 = wstar1 + dzi*bprod0(k)
      end do
      tke1   = tke1 / tmp1
      wstar1 = ( 2.5_r8 * max( 0._r8, wstar1 ) )**(1._r8/3._r8)
      tmp1 = ( qv0(1) + ql0(1) + qi0(1) )
      tmp2 = ( t0(1) - xlv * ql0(1) / cp - xls * qi0(1) / cp ) / ( ( p0(1) / p00 )**rovcp )
      call conden( p0(1), tmp2, tmp1, th, qv, ql, qi, qse, id_check )
      wstar2 = ( max( 0._r8, bprod0(0) ) * pblhz )**(1._r8/3._r8)
      wstar2 = max( 0._r8, min( wstar2, 10._r8 ) )

      ! OPTION
      ! Re-define 'tke1' as 'tke' at the 1st interface
      if( kiss .eq. 0 ) then
         ! Jun.28.2011. Directly use 'tkes' from the UW PBL scheme which does not include transport term.
         tke1   = tke0(0)
         ! tke1   = tkes 
      else
         tke1   = tke0(1)
      endif
      tke1   = max(1.e-5_r8,min(3._r8,tke1)) 
      ! OPTION

      tmp1     = 0._r8
      tmp2     = 0._r8
      tkePBL   = 0._r8
      wstarPBL = 0._r8
      kc       = kpblhm    ! Top interface of the average domain.
      do k = 0, kc         ! Here, 'k' is an interfacial layer index.
         if( k .eq. 0 ) then
            dpi = ps0(0) - p0(1)
            dzi = z0(1) - zs0(0)
         elseif( k .eq. kc ) then
            dpi = p0(kc) - ps0(kc)
            dzi = zs0(kc) - z0(kc)
         else
            dpi = p0(k) - p0(k+1)
            dzi = z0(k+1) - z0(k)
         endif
         tmp1   = tmp1   + dpi
         tmp2   = tmp2   + dzi
         tkePBL = tkePBL + dpi*tke0(k)
         wstarPBL = wstarPBL + dzi*bprod0(k)
      end do
      tkePBL   = tkePBL / tmp1
      tkePBL   = max(1.e-5_r8,min(3._r8,tkePBL)) 
      wstarPBL = ( 2.5_r8 * max( 0._r8, wstarPBL ) )**(1._r8/3._r8)

      tke1 = tkePBL

      ! ---------------------------------------------------------------------------- !
      ! Sep.16.2011. Compute forbidden area fraction by subgrid variation of surface !
      !              topograpgic height.                                             !
      !              Following 'tms', I am using 'sgh30' not 'sgh'.                  !
      !              Impose a upper limit of 'a_oro_max' on the 'a_oro'.             !   
      ! ---------------------------------------------------------------------------- !
  
      a_oro = sgh30 / norm_sgh
      a_oro = max( 0._r8, min( a_oro_max, a_oro ) )
  
      ! --------------------------------------------------------------------------------------------------------- !
      ! Sep.11.2011. Insert a condition preventing negative thl, qt , tracers both within wake and non-wake areas !
      !              before wake-inhomogeneity adjustment.                                                        !
      ! --------------------------------------------------------------------------------------------------------- !

!lim  In fact, below limiters before buoyancy adjustment are not necessary, since buoyancy adjument of cold pool
!lim  can be done regardless of what is the value of 'delta_xxx_PBL'. So, I commented out below limiters.     
!lim  Just for computing 'delta_thv_PBL_raw' below, impose a limiter for 'delta_thl_PBL_raw, delta_qt_PBL_raw'
!lim  Also add a limiter at the surface and top interface of the PBL for full consistency.

!lim  thl0min_PBL          = thl0bot(1)
      qt0min_PBL           = qt0bot(1)
      tr0min_PBL(:ncnst)   = tr0bot(1,:ncnst)
      do k = 1, kpblhm
!lim     thl0min_PBL       = min( thl0min_PBL, thl0(k) ) 
         qt0min_PBL        = min( qt0min_PBL,   qt0(k) ) 
         do mt = 1, ncnst
            tr0min_PBL(mt) = min( tr0min_PBL(mt), tr0(k,mt) )              
         enddo
      enddo
!lim  thl0min_PBL       = min( thl0min_PBL, thl0top(kpblhm) ) 
      qt0min_PBL        = min( qt0min_PBL,   qt0top(kpblhm) ) 
      do mt = 1, ncnst
         tr0min_PBL(mt) = min( tr0min_PBL(mt), tr0top(kpblhm,mt) )              
      enddo
             
!lim  delta_thl_PBL_raw = min( max( -thl0min_PBL, delta_thl_PBL_raw ), thl0min_PBL * awk_PBL_raw / ( 1._r8 - awk_PBL_raw ) ) 
      delta_qt_PBL_raw  = min( max( qmin(1) - qt0min_PBL, delta_qt_PBL_raw ), & 
                                  ( qt0min_PBL - qmin(1) ) * awk_PBL_raw / ( 1._r8 - awk_PBL_raw ) )
      do mt = 1, ncnst    
         delta_tr_PBL_raw(mt) = min( max( qmin(mt) - tr0min_PBL(mt), delta_tr_PBL_raw(mt) ), & 
                                        ( tr0min_PBL(mt) - qmin(mt) ) * awk_PBL_raw / ( 1._r8 - awk_PBL_raw ) )
      enddo

!lim2 It seems that below 4 limiters are too physically restricted without influencing model crash and 
!lim2 computation time. So, I removed them.  

      delta_thl_PBL_raw = max( - 2.0_r8, min( 2.0_r8, delta_thl_PBL_raw ) )
      delta_qt_PBL_raw  = max( - 0.2_r8*qt0min_PBL, min( 0.2_r8*qt0min_PBL, delta_qt_PBL_raw ) )
      delta_u_PBL_raw   = max( - 2.0_r8, min( 2.0_r8 , delta_u_PBL_raw ) )
      delta_v_PBL_raw   = max( - 2.0_r8, min( 2.0_r8 , delta_v_PBL_raw ) )

!lim  In fact, below computation from individual layer computation is also unnecessary.
!lim  We only need to do a bulk computation, which is also good for saving computation time.
!lim  However, in order to use 'conden' subroutine with p0(k), it seems that below computation is inevitable.
!lim  Instead of using the above block, I simply used 'max( qt0(k) + delta_qt_PBL_raw, qmin(1) )' as an argument below.

      delta_thv_PBL_raw = 0._r8
      do k = 1, kpblhm  ! Here, 'k' is a layer index.
         call conden( p0(k), thl0(k) + delta_thl_PBL_raw, qt0(k) + delta_qt_PBL_raw, & 
                      th, qv, ql, qi, qse, id_check )
         thv = th * ( 1._r8 + zvir * qv - ql - qi )
         delta_thv_PBL_raw = delta_thv_PBL_raw + ( thv - thv0(k) ) * dp0(k)
      enddo
      delta_thv_PBL_raw = delta_thv_PBL_raw / pblhp

      awk_PBL_max = 1._r8 - au_max
      if( awk_PBL_raw .lt. 0._r8 .or. awk_PBL_raw .gt. awk_PBL_max ) then
         write(iulog,*)
         write(iulog,*) 'UNICON : Unreasonable wake area awk_PBL_raw before inhomogeneity adjustment'
         write(iulog,*) 'awk_PBL_raw = ', awk_PBL_raw 
         call endrun('STOP : UNICON')
         write(iulog,*)  
      endif
      if( delta_thv_PBL_raw .lt. 0._r8 ) then
         awk_PBL                  = 0._r8
         delta_thl_PBL            = 0._r8
         delta_qt_PBL             = 0._r8
         delta_u_PBL              = 0._r8
         delta_v_PBL              = 0._r8
         delta_thv_PBL            = 0._r8
         do mt = 1, ncnst    
            delta_tr_PBL(mt)      = 0._r8
         enddo
         delta_w_PBL              = 0._r8
      else
         tmp1                     = delta_thv_wc * awk_PBL_raw / ( awk_PBL_raw - 1._r8 ) / max( nonzero, delta_thv_PBL_raw )
         tmp2                     = erfc( tmp1 / sqrt( 3.141592_r8 ) )
         tmp3                     = exp( - tmp1**2._r8 / 3.141592_r8 )
         tmp4                     = ( 1._r8 - awk_PBL_raw ) / ( 1._r8 - tmp2 * awk_PBL_raw ) 
         awk_PBL                  = tmp2 * awk_PBL_raw
         delta_thl_PBL            = tmp3 * tmp4 *   delta_thl_PBL_raw    
         delta_qt_PBL             = tmp3 * tmp4 *    delta_qt_PBL_raw
         delta_u_PBL              = tmp3 * tmp4 *     delta_u_PBL_raw
         delta_v_PBL              = tmp3 * tmp4 *     delta_v_PBL_raw
         delta_thv_PBL            = tmp3 * tmp4 *   delta_thv_PBL_raw 
         do mt = 1, ncnst    
            delta_tr_PBL(mt)      = tmp3 * tmp4 * delta_tr_PBL_raw(mt)
         enddo

    !lim Oct.4. Even though we apply non-negative limiter for 'delta_qt_PBL' and 'delta_tr_PBL' here, there is no way     
    !lim        to consistently change 'delta_thv_PBL' that influences the computation of 'delta_thv_PBL' below.
    !lim        Thus, the limiters for 'delta_qt_PBL' and 'delta_tr_PBL' are done later without any option.

         delta_w_PBL = kw_omega * sqrt( kstar * 0.5_r8 * ( g / thv_ref ) * pblhz * awk_PBL * max( 0._r8, delta_thv_PBL ) )
         if( i_energy_coldpool .eq. 2 ) then
             delta_w_PBL = sqrt( 2._r8 ) * ( g / thv_ref ) * b1 * pblhz**2._r8 * del_wk0 * awk_PBL * max( 0._r8, delta_thv_PBL ) 
             delta_w_PBL = delta_w_PBL**(1._r8/3._r8) * sqrt( awk_PBL ) / max( nonzero, ( 1._r8 - awk_PBL )**(1._r8/6._r8) )
         endif

         delta_w_PBL              = min( 10._r8, delta_w_PBL )
         if( awk_PBL .lt. 1.e-3_r8 ) delta_w_PBL = 0._r8
      endif
      if( delta_thv_PBL .lt. 0._r8 .or. awk_PBL .lt. 0._r8 .or. awk_PBL .gt. awk_PBL_max ) then
         write(iulog,*)
         write(iulog,*) 'UNICON : Unreasonable wake properties after inhomogeneity adjustment'
         write(iulog,*) 'delta_thv_PBL, awk_PBL = ', delta_thv_PBL, awk_PBL 
         call endrun('STOP : UNICON')  
         write(iulog,*)
      endif

      ! --------------------------------------------------------------------------------------------------------- !
      ! Sep.11.2011. Insert a condition preventing negative thl, qt , tracers both within wake and non-wake areas !
      !              after wake-inhomogeneity adjustment.                                                         !
      ! Aug.03.2012. I may need to below limiters in association with the use of internally computed              !
      !              cdelta_s and cdelta_w.                                                                       ! 
      ! Oct.04.2014. Consider 'tmp1 = 1/cdelta_s' too.                                                            !
      !              Also remove the limiter for 'delta_thl_PBL' since it is highly likely that it does not       !
      !              happens but takes computation time.                                                          ! 
      ! --------------------------------------------------------------------------------------------------------- !
             
      tmp1 = au_base_min_ocn * cadj_area_ocn + ( au_base_min_lnd * cadj_area_lnd - au_base_min_ocn * cadj_area_ocn ) * landfrac 

!lim  delta_thl_PBL = min( max(-tmp1 * thl0min_PBL, delta_thl_PBL ), thl0min_PBL * awk_PBL / ( 1._r8 - awk_PBL ) ) 
      delta_qt_PBL  = min( max( tmp1 * ( qmin(1) - qt0min_PBL ), delta_qt_PBL ), & 
                              ( qt0min_PBL - qmin(1) ) * awk_PBL / ( 1._r8 - awk_PBL ) )
      do mt = 1, ncnst    
         delta_tr_PBL(mt) = min( max( tmp1 * ( qmin(mt) - tr0min_PBL(mt) ), delta_tr_PBL(mt) ), & 
                                    ( tr0min_PBL(mt) - qmin(mt) ) * awk_PBL / ( 1._r8 - awk_PBL ) )
      enddo
      delta_tr_PBL(ixnumliq) = 0._r8
      delta_tr_PBL(ixnumice) = 0._r8

!lim4 I removed below 4 limiters because it seems to be too strict without any physical reason.

      delta_thl_PBL = max( - 2.0_r8, min( 2.0_r8, delta_thl_PBL ) )
      delta_qt_PBL  = max( - 0.2_r8*qt0min_PBL, min( 0.2_r8*qt0min_PBL, delta_qt_PBL ) )
      delta_u_PBL   = max( - 2.0_r8, min( 2.0_r8 , delta_u_PBL ) )
      delta_v_PBL   = max( - 2.0_r8, min( 2.0_r8 , delta_v_PBL ) )

!lim  Replace the above limiter as below.
!lim  Practically, only need to worry about tracer, and I should take into account of the effect of 'cdelta_s'.
!lim  Note that below 'tmp1 = 1 / cdelta_s' ( I should double check this ).

!lim  delta_qt_PBL = min( max( tmp1 * ( qmin(1) - qt0PBL ), delta_qt_PBL ), & 
!lim                    ( qt0PBL - qmin(1) ) * awk_PBL / ( 1._r8 - awk_PBL ) )
!lim  do mt = 1, ncnst    
!lim     delta_tr_PBL(mt) = min( max( tmp1 * ( qmin(mt) - tr0PBL(mt) ), delta_tr_PBL(mt) ), & 
!lim                           ( tr0PBL(mt) - qmin(mt) ) * awk_PBL / ( 1._r8 - awk_PBL ) )
!lim  enddo
!lim  ! Nov.28.2012. For consistent treatment in the mixing with organized detrained airs with condensate later, 
!lim  !              simply assume below two lines.
!lim  delta_tr_PBL(ixnumliq) = 0._r8
!lim  delta_tr_PBL(ixnumice) = 0._r8

!lim
!lim
!lim

      ! -------------------------------------------------------------------------------------------------------------------- !
      ! Sep.16.2011. Define effective 'eps_wk_eff, del_wk_eff' to allow the initial development of wake                      ! 
      !              from zero area by effectively reduce too large detrainment dilution effect of thermodynamic             !
      !              scalars, 'del_wk / awk_PBL * ( 1._r8 - awk_PBL) )' due to 'zero awk_PBL' at the beginning.              !
      ! -------------------------------------------------------------------------------------------------------------------- ! 

      if( i_energy_coldpool .eq. 1 .or. i_energy_coldpool .eq. 2 ) then
         eps_wk_eff  = eps_wk0 * awk_PBL                                                          
         del_wk_eff  = del_wk0 * awk_PBL
      endif

      ! ------------------------------------------------------------------------------------------------------------- !
      ! Compute convective organization parameter, 0 <= org_rad <= 1.                                                 ! 
      ! This  'org_rad' is used for computing updraft plume radius at surface                                         !
      ! Later 'org_ent' is defined for lateral mixing later as an option.                                             !
      ! The 'cuorg' can be computed in many different ways ( e.g., using convective precipitation flux at surface ).  !
      ! Here, I am using the ratio of 'downdraft mass flux' to 'updraft mass flux' at the previous time step because  !
      ! conceptually, strong downward motion is likely to generate meso-scale organized circulation near surface.     !
      ! ------------------------------------------------------------------------------------------------------------- !

      a_forbid = a_oro + awk_PBL
      a_forbid = max( 0._r8, min( awk_PBL_max, a_forbid ) ) 
      cuorg    = a_forbid / awk_PBL_max
      cuorg    = max( 0._r8, min( 1._r8, cuorg ) )

      if( orgfeedback_off ) then
         cuorg = 0._r8
      endif

      org_rad = fac_org_rad * cuorg                                   ! Mar.11.2011. Most Advanced.     
      org_rad = max( 0._r8, min( 1._r8, org_rad ) )                   ! Force to be [0,1] after multiplication of fac_org_rad.

      if( orp .ge. 0._r8 ) then
         tmp1 = org_rad**orp
      else
         tmp1 = 0.5_r8 * ( 1._r8 + sin( 3.141592_r8 * ( org_rad - 0.5_r8 ) ) )
      endif

      ! ------------------------------------------------------------------------------- !
      ! Compute Land-Ocean Avarage Parameters                                           ! 
      ! ------------------------------------------------------------------------------- !

      criqc               = criqc_ocn        + ( criqc_lnd        -        criqc_ocn ) * landfrac
      c0_ac               = c0_ac_ocn        + ( c0_ac_lnd        -        c0_ac_ocn ) * landfrac
      kevp_rain           = kevp_rain_ocn    + ( kevp_rain_lnd    -    kevp_rain_ocn ) * landfrac
      kevp_snow           = kevp_snow_ocn    + ( kevp_snow_lnd    -    kevp_snow_ocn ) * landfrac
      kevp_rain_dn        = kevp_rain_dn_ocn + ( kevp_rain_dn_lnd - kevp_rain_dn_ocn ) * landfrac
      kevp_snow_dn        = kevp_snow_dn_ocn + ( kevp_snow_dn_lnd - kevp_snow_dn_ocn ) * landfrac

      tmp2 = au_base_max_lnd
      if( iau_base_lnd .eq. 1 ) then
         tmp2 = au_base_min_lnd * ( 1._r8 - awk_PBL_max )
      endif
      tmp3 = au_base_max_ocn
      if( iau_base_ocn .eq. 1 ) then
         tmp3 = au_base_min_ocn * ( 1._r8 - awk_PBL_max )
      endif

      au_base_min = au_base_min_ocn + ( au_base_min_lnd - au_base_min_ocn ) * landfrac
      au_base_max = tmp3            + ( tmp2            -            tmp3 ) * landfrac
      Ro_min      = Ro_min_ocn      + ( Ro_min_lnd      -      Ro_min_ocn ) * landfrac
      Ro_max      = Ro_max_ocn      + ( Ro_max_lnd      -      Ro_max_ocn ) * landfrac
      sigmaR_min  = sigmaR_min_ocn  + ( sigmaR_min_lnd  -  sigmaR_min_ocn ) * landfrac
      sigmaR_max  = sigmaR_max_ocn  + ( sigmaR_max_lnd  -  sigmaR_max_ocn ) * landfrac
      kw_min      = kw_min_ocn      + ( kw_min_lnd      -      kw_min_ocn ) * landfrac
      kw_max      = kw_max_ocn      + ( kw_max_lnd      -      kw_max_ocn ) * landfrac

      Ro          = Ro_min          + ( Ro_max          -          Ro_min ) * tmp1
      sigmaR      = sigmaR_min      + ( sigmaR_max      -      sigmaR_min ) * tmp1
      kw          = kw_min          + ( kw_max          -          kw_min ) * tmp1

      au_base     = au_base_min     + ( au_base_max     -     au_base_min     ) * org_rad
      au_base_ocn = au_base_min_ocn + ( tmp3            -     au_base_min_ocn ) * org_rad
      au_base_lnd = au_base_min_lnd + ( tmp2            -     au_base_min_lnd ) * org_rad

      AOVTU = au_base_ocn * cadj_area_ocn + ( au_base_lnd * cadj_area_lnd - au_base_ocn * cadj_area_ocn ) * landfrac 
      cdelta_s = ( 1._r8 - awk_PBL ) / max( nonzero, AOVTU )
      cdelta_s = max( 1._r8, min( cdelta_s, 100._r8 ) ) 

      if( orgfeedback_off ) then
         cdelta_s = 0._r8
      endif

      cdelta_w = sqrt( cdelta_s ) 

      if( i_energy_coldpool .eq. 2 ) then
         cdelta_w = sqrt( cdelta_s ) / max( nonzero, ( AOVTU + awk_PBL )**(1._r8/6._r8) )   
      endif

      alpha_cri = 0._r8  

      ! --------------------------------------------- !
      !                                               !
      ! Main computation of cumulus updraft evolution !
      !                                               !
      ! --------------------------------------------- !

      do iter = 1, niter

         ! ------------------------------------------------------------------------------------------- !
         ! Aug.01.2011. Brian Juwon Park's Birthday. This 'do iter' routine can be used for explicit   !
         !              mixing of convective updraft and downdraft with different mixing environmental !
         !              airs because vertical evolution of convective updraft and downdraft is highly  !
         !              non-linear to the properties of mixing environmental airs. Thus, the simple    !
         !              use of mean environmental airs for mixing is not ideal. This explicit mixing   !
         !              process will inevitably increase computation time but seems to be important.   !
         !                 (1) iter = 1 : mixing with mean environmental airs at the current time step !
         !                                ( with the probability of '1._r8 - cuorg' )                  ! 
         !                 (2) iter = 2 : mixing with detrained + cumulus updraft airs at the previous !
         !                                time step ( with the probability of 'cuorg' )                !
         !              Note that this 'cuorg_mxen' is used only for 'org_ent', neither 'org_rad' nor  !
         !             'org_src' since 'cuorg_mxen' is designed to choose mixing environmental airs.   ! 
         ! ------------------------------------------------------------------------------------------- !

         if( niter .eq. 1 ) then
            cuorg_mxen = cuorg
         else
            if( iter .eq. 1 ) then 
               cuorg_mxen = 0._r8
            elseif( iter .eq. 2 ) then   
               cuorg_mxen = 1._r8
            endif
         endif

         ! ------------------------------------------------------------------------------------------------------- !
         !                                                                                                         !                
         ! Iteration for treating accretion should start here.                                                     ! 
         !                                                                                                         !
         ! The variables used for the next accretion iteration computations within 'subroutine prod_prep_up' are   !
         !                                                                                                         !
         !  1. flxrain_msfc(k,msfc), flxsnow_msfc(k,msfc)                                                          !
         !  2. a_p_msfc(k,msfc)                                                                                    ! 
         !  3. am_us_msfc(k,msfc)                                                                                  !
         !  4. am_pu_msfc(k,msfc)                                                                                  ! 
         !                                                                                                         !
         ! all of which will be computed before the end of accretion iteration loop.                               ! 
         !                                                                                                         !
         ! ------------------------------------------------------------------------------------------------------- !

         do iacc = 1, nacc

            if( iacc .eq. 1 ) then

               ! -------------------------------------------------------------------------------------- !
               ! Initialize below variables only at the first iteration, since at the second iteration, !
               ! non-zero values will be used for accretion.                                            ! 
               ! -------------------------------------------------------------------------------------- !

               flxrain_msfc(0:mkx,:nseg)                = 0._r8
               flxsnow_msfc(0:mkx,:nseg)                = 0._r8
               flxtrrs_msfc(0:mkx,:nseg,:ncnst)         = 0._r8

               a_p_msfc(0:mkx,:nseg)                    = 0._r8

               am_u_msfc(:mkx,:nseg)                    = 0._r8 
               am_up_msfc(:mkx,:nseg)                   = 0._r8 
               am_us_msfc(:mkx,:nseg)                   = 0._r8 

               am_evp_msfc(:mkx,:nseg)                  = 0._r8 
               am_pu_msfc(:mkx,:nseg)                   = 0._r8
               am_pd_msfc(:mkx,:nseg)                   = 0._r8
               am_pr_msfc(:mkx,:nseg)                   = 0._r8
               am_ps_msfc(:mkx,:nseg)                   = 0._r8

            endif

            ! --------------------------------------------------- !
            !                                                     !                
            !                                                     ! 
            ! Iteration for treating accretion should start here. ! 
            !                                                     !
            !                                                     !
            ! --------------------------------------------------- !

            ! --------------------------------------------------------------------------- !
            ! Initialization of ensemble-mean arrays in each layer inside the 'iter' loop !
            ! --------------------------------------------------------------------------- !

            flxrain(0:mkx)                           = 0._r8
            flxsnow(0:mkx)                           = 0._r8
            flxtrrs(0:mkx,:ncnst)                    = 0._r8

            cmf_u_mix(:mkx)                          = 0._r8
            cmf_r(:mkx)                              = 0._r8 
            thl_r(:mkx)                              = 0._r8
            qt_r(:mkx)                               = 0._r8
            u_r(:mkx)                                = 0._r8
            v_r(:mkx)                                = 0._r8
            ql_r(:mkx)                               = 0._r8
            qi_r(:mkx)                               = 0._r8
            tr_r(:mkx,:ncnst)                        = 0._r8

            cmf_r2(:mkx)                             = 0._r8 
            thl_r2(:mkx)                             = 0._r8
            qt_r2(:mkx)                              = 0._r8
            u_r2(:mkx)                               = 0._r8
            v_r2(:mkx)                               = 0._r8
            ql_r2(:mkx)                              = 0._r8
            qi_r2(:mkx)                              = 0._r8
            tr_r2(:mkx,:ncnst)                       = 0._r8

            thl_u(0:mkx)                             = 0._r8
            qt_u(0:mkx)                              = 0._r8
            u_u(0:mkx)                               = 0._r8
            v_u(0:mkx)                               = 0._r8
            cmf_u(0:mkx)                             = 0._r8
            w_u(0:mkx)                               = 0._r8
            wa_u(0:mkx)                              = 0._r8
            a_u(0:mkx)                               = 0._r8
            num_u(0:mkx)                             = 0._r8
            rad_u(0:mkx)                             = 0._r8  
            ql_u(0:mkx)                              = 0._r8   
            qi_u(0:mkx)                              = 0._r8   
            tr_u(0:mkx,:ncnst)                       = 0._r8   
            qla_u(0:mkx)                             = 0._r8   
            qia_u(0:mkx)                             = 0._r8   
            thva_u(0:mkx)                            = 0._r8   
            
            cmf_u_dia(:mkx)                          = 0._r8
            evp_thll_u(:mkx)                         = 0._r8
            evp_qtl_u(:mkx)                          = 0._r8
            evp_thli_u(:mkx)                         = 0._r8
            evp_qti_u(:mkx)                          = 0._r8
            prep_thll_u(:mkx)                        = 0._r8
            prep_qtl_u(:mkx)                         = 0._r8
            prep_thli_u(:mkx)                        = 0._r8
            prep_qti_u(:mkx)                         = 0._r8
            eff_ql_u(:mkx)                           = 0._r8
            eff_qi_u(:mkx)                           = 0._r8
            PGF_u_u(:mkx)                            = 0._r8
            PGF_v_u(:mkx)                            = 0._r8
            evp_tr_u(:mkx,:ncnst)                    = 0._r8
            prep_tr_u(:mkx,:ncnst)                   = 0._r8
            eff_tr_u(:mkx,:ncnst)                    = 0._r8

            thl_srcd(:mkx)                           = 0._r8
            qt_srcd(:mkx)                            = 0._r8
            u_srcd(:mkx)                             = 0._r8
            v_srcd(:mkx)                             = 0._r8
            tr_srcd(:mkx,:ncnst)                     = 0._r8
            f_srcd(:mkx)                             = 0._r8
            ql_srcd(:mkx)                            = 0._r8
            qi_srcd(:mkx)                            = 0._r8

            thl_srcds(:mkx,:nseg,1:3)                = 0._r8
            qt_srcds(:mkx,:nseg,1:3)                 = 0._r8
            u_srcds(:mkx,:nseg,1:3)                  = 0._r8
            v_srcds(:mkx,:nseg,1:3)                  = 0._r8
            tr_srcds(:mkx,:nseg,1:3,:ncnst)          = 0._r8
            f_srcds(:mkx,:nseg,1:3)                  = 0._r8
            ql_srcds(:mkx,:nseg,1:3)                 = 0._r8
            qi_srcds(:mkx,:nseg,1:3)                 = 0._r8

            thl_srcr(:mkx)                           = 0._r8
            qt_srcr(:mkx)                            = 0._r8
            u_srcr(:mkx)                             = 0._r8
            v_srcr(:mkx)                             = 0._r8
            tr_srcr(:mkx,:ncnst)                     = 0._r8
            f_srcr(:mkx)                             = 0._r8
            ql_srcr(:mkx)                            = 0._r8
            qi_srcr(:mkx)                            = 0._r8

            thl_srcr2(:mkx)                          = 0._r8
            qt_srcr2(:mkx)                           = 0._r8
            u_srcr2(:mkx)                            = 0._r8
            v_srcr2(:mkx)                            = 0._r8
            tr_srcr2(:mkx,:ncnst)                    = 0._r8
            f_srcr2(:mkx)                            = 0._r8
            ql_srcr2(:mkx)                           = 0._r8
            qi_srcr2(:mkx)                           = 0._r8

            thl_srcrs(:mkx,:nseg,1:3)                = 0._r8
            qt_srcrs(:mkx,:nseg,1:3)                 = 0._r8
            u_srcrs(:mkx,:nseg,1:3)                  = 0._r8
            v_srcrs(:mkx,:nseg,1:3)                  = 0._r8
            tr_srcrs(:mkx,:nseg,1:3,:ncnst)          = 0._r8
            f_srcrs(:mkx,:nseg,1:3)                  = 0._r8
            ql_srcrs(:mkx,:nseg,1:3)                 = 0._r8
            qi_srcrs(:mkx,:nseg,1:3)                 = 0._r8

            thl_srcrs2(:mkx,:nseg,1:3)               = 0._r8
            qt_srcrs2(:mkx,:nseg,1:3)                = 0._r8
            u_srcrs2(:mkx,:nseg,1:3)                 = 0._r8
            v_srcrs2(:mkx,:nseg,1:3)                 = 0._r8
            tr_srcrs2(:mkx,:nseg,1:3,:ncnst)         = 0._r8
            f_srcrs2(:mkx,:nseg,1:3)                 = 0._r8
            ql_srcrs2(:mkx,:nseg,1:3)                = 0._r8
            qi_srcrs2(:mkx,:nseg,1:3)                = 0._r8

            cmf_ru(:mkx)                             = 0._r8 
            thl_ru(:mkx)                             = 0._r8
            qt_ru(:mkx)                              = 0._r8
            u_ru(:mkx)                               = 0._r8
            v_ru(:mkx)                               = 0._r8
            ql_ru(:mkx)                              = 0._r8
            qi_ru(:mkx)                              = 0._r8
            tr_ru(:mkx,:ncnst)                       = 0._r8

            cmf_ru2(:mkx)                            = 0._r8 
            thl_ru2(:mkx)                            = 0._r8
            qt_ru2(:mkx)                             = 0._r8
            u_ru2(:mkx)                              = 0._r8
            v_ru2(:mkx)                              = 0._r8
            ql_ru2(:mkx)                             = 0._r8
            qi_ru2(:mkx)                             = 0._r8
            tr_ru2(:mkx,:ncnst)                      = 0._r8

            cmf_ad(0:mkx,:mkx,1:nseg,1:3)            = 0._r8
            thl_ad(0:mkx,:mkx,1:nseg,1:3)            = 0._r8
            qt_ad(0:mkx,:mkx,1:nseg,1:3)             = 0._r8
            u_ad(0:mkx,:mkx,1:nseg,1:3)              = 0._r8
            v_ad(0:mkx,:mkx,1:nseg,1:3)              = 0._r8
            w_ad(0:mkx,:mkx,1:nseg,1:3)              = 0._r8
            a_ad(0:mkx,:mkx,1:nseg,1:3)              = 0._r8
            ql_ad(0:mkx,:mkx,1:nseg,1:3)             = 0._r8
            qi_ad(0:mkx,:mkx,1:nseg,1:3)             = 0._r8
            tr_ad(0:mkx,:mkx,1:nseg,1:3,:ncnst)      = 0._r8

            dpad(:mkx,:mkx,1:nseg,1:3)               = 0._r8

            cmf_ar(:mkx,:mkx,1:nseg,1:3)             = 0._r8
            thl_ar(:mkx,:mkx,1:nseg,1:3)             = 0._r8
            qt_ar(:mkx,:mkx,1:nseg,1:3)              = 0._r8
            u_ar(:mkx,:mkx,1:nseg,1:3)               = 0._r8
            v_ar(:mkx,:mkx,1:nseg,1:3)               = 0._r8
            tr_ar(:mkx,:mkx,1:nseg,1:3,:ncnst)       = 0._r8
            ql_ar(:mkx,:mkx,1:nseg,1:3)              = 0._r8
            qi_ar(:mkx,:mkx,1:nseg,1:3)              = 0._r8

            cmf_ad_dia(:mkx,:mkx,1:nseg,1:3)         = 0._r8
            evp_thll_ad(:mkx,:mkx,1:nseg,1:3)        = 0._r8
            evp_qtl_ad(:mkx,:mkx,1:nseg,1:3)         = 0._r8
            evp_thli_ad(:mkx,:mkx,1:nseg,1:3)        = 0._r8
            evp_qti_ad(:mkx,:mkx,1:nseg,1:3)         = 0._r8
            prep_thll_ad(:mkx,:mkx,1:nseg,1:3)       = 0._r8
            prep_qtl_ad(:mkx,:mkx,1:nseg,1:3)        = 0._r8
            prep_thli_ad(:mkx,:mkx,1:nseg,1:3)       = 0._r8
            prep_qti_ad(:mkx,:mkx,1:nseg,1:3)        = 0._r8
            eff_ql_ad(:mkx,:mkx,1:nseg,1:3)          = 0._r8
            eff_qi_ad(:mkx,:mkx,1:nseg,1:3)          = 0._r8
            PGF_u_ad(:mkx,:mkx,1:nseg,1:3)           = 0._r8
            PGF_v_ad(:mkx,:mkx,1:nseg,1:3)           = 0._r8
            evp_tr_ad(:mkx,:mkx,1:nseg,1:3,:ncnst)   = 0._r8
            prep_tr_ad(:mkx,:mkx,1:nseg,1:3,:ncnst)  = 0._r8
            wdep_tr_ad(:mkx,:mkx,1:nseg,1:3,:ncnst)  = 0._r8
            eff_tr_ad(:mkx,:mkx,1:nseg,1:3,:ncnst)   = 0._r8

            cmf_d(0:mkx)                             = 0._r8
            thl_d(0:mkx)                             = 0._r8
            qt_d(0:mkx)                              = 0._r8
            u_d(0:mkx)                               = 0._r8
            v_d(0:mkx)                               = 0._r8
            w_d(0:mkx)                               = 0._r8
            a_d(0:mkx)                               = 0._r8
            wa_d(0:mkx)                              = 0._r8
            ql_d(0:mkx)                              = 0._r8   
            qi_d(0:mkx)                              = 0._r8   
            tr_d(0:mkx,:ncnst)                       = 0._r8   
            qla_d(0:mkx)                             = 0._r8   
            qia_d(0:mkx)                             = 0._r8   

            cmf_d_dia(:mkx)                          = 0._r8
            evp_thll_d(:mkx)                         = 0._r8
            evp_qtl_d(:mkx)                          = 0._r8
            evp_thli_d(:mkx)                         = 0._r8
            evp_qti_d(:mkx)                          = 0._r8
            prep_thll_d(:mkx)                        = 0._r8
            prep_qtl_d(:mkx)                         = 0._r8
            prep_thli_d(:mkx)                        = 0._r8
            prep_qti_d(:mkx)                         = 0._r8
            eff_ql_d(:mkx)                           = 0._r8
            eff_qi_d(:mkx)                           = 0._r8
            PGF_u_d(:mkx)                            = 0._r8
            PGF_v_d(:mkx)                            = 0._r8
            evp_tr_d(:mkx,:ncnst)                    = 0._r8
            prep_tr_d(:mkx,:ncnst)                   = 0._r8
            eff_tr_d(:mkx,:ncnst)                    = 0._r8

            cmf_rd(:mkx)                             = 0._r8
            thl_rd(:mkx)                             = 0._r8
            qt_rd(:mkx)                              = 0._r8
            u_rd(:mkx)                               = 0._r8
            v_rd(:mkx)                               = 0._r8
            ql_rd(:mkx)                              = 0._r8
            qi_rd(:mkx)                              = 0._r8
            tr_rd(:mkx,:ncnst)                       = 0._r8

            qlten_sub(:mkx)                          = 0._r8
            qiten_sub(:mkx)                          = 0._r8

            rqc_l(:mkx)                              = 0._r8
            rqc_i(:mkx)                              = 0._r8
            rqc(:mkx)                                = 0._r8
            rnc_l(:mkx)                              = 0._r8
            rnc_i(:mkx)                              = 0._r8

            cmf_det(:mkx)                            = 0._r8
            ql_det(:mkx)                             = 0._r8
            qi_det(:mkx)                             = 0._r8

            qlten_det(:mkx)                          = 0._r8 
            qiten_det(:mkx)                          = 0._r8 

            am_d_msfc(:mkx,:nseg)                    = 0._r8
            
            am_u(:mkx)                               = 0._r8
            am_d(:mkx)                               = 0._r8

            qlm_u_msfc(:mkx,:nseg)                   = 0._r8
            qim_u_msfc(:mkx,:nseg)                   = 0._r8
            thlm_u_msfc(:mkx,:nseg)                  = 0._r8
            qtm_u_msfc(:mkx,:nseg)                   = 0._r8
            um_u_msfc(:mkx,:nseg)                    = 0._r8
            vm_u_msfc(:mkx,:nseg)                    = 0._r8
            trm_u_msfc(:mkx,:nseg,:ncnst)            = 0._r8

            qlm_u(:mkx)                              = 0._r8
            qim_u(:mkx)                              = 0._r8
            thlm_u(:mkx)                             = 0._r8
            qtm_u(:mkx)                              = 0._r8
            um_u(:mkx)                               = 0._r8
            vm_u(:mkx)                               = 0._r8
            trm_u(:mkx,:ncnst)                       = 0._r8

            qlm_d_msfc(:mkx,:nseg)                   = 0._r8
            qim_d_msfc(:mkx,:nseg)                   = 0._r8

            qlm_d(:mkx)                              = 0._r8
            qim_d(:mkx)                              = 0._r8

            am_s(:mkx)                               = 0._r8
            am_r(:mkx)                               = 0._r8
            am_up(:mkx)                              = 0._r8
            am_us(:mkx)                              = 0._r8

            a_p(0:mkx)                               = 0._r8

            am_evp(:mkx)                             = 0._r8 
            am_pu(:mkx)                              = 0._r8
            am_pd(:mkx)                              = 0._r8
            am_pr(:mkx)                              = 0._r8
            am_ps(:mkx)                              = 0._r8

            thl_u_msfc(0:mkx,:nseg)                  = 0._r8
            qt_u_msfc(0:mkx,:nseg)                   = 0._r8
            u_u_msfc(0:mkx,:nseg)                    = 0._r8
            v_u_msfc(0:mkx,:nseg)                    = 0._r8
            w_u_msfc(0:mkx,:nseg)                    = 0._r8
            ql_u_msfc(0:mkx,:nseg)                   = 0._r8
            qi_u_msfc(0:mkx,:nseg)                   = 0._r8
            tr_u_msfc(0:mkx,:nseg,:ncnst)            = 0._r8
            cmf_u_msfc(0:mkx,:nseg)                  = 0._r8
            a_u_msfc(0:mkx,:nseg)                    = 0._r8
            num_u_msfc(0:mkx,:nseg)                  = 0._r8 
            rad_u_msfc(0:mkx,:nseg)                  = 0._r8

            eps0_u_msfc(0:mkx,:nseg)                 = 0._r8
            eps_u_msfc(0:mkx,:nseg)                  = 0._r8
            del_u_msfc(0:mkx,:nseg)                  = 0._r8
            eeps_u_msfc(0:mkx,:nseg)                 = 0._r8
            ddel_u_msfc(0:mkx,:nseg)                 = 0._r8
            xc_u_msfc(0:mkx,:nseg)                   = 0._r8
            xs_u_msfc(0:mkx,:nseg)                   = 0._r8
            xemin_u_msfc(0:mkx,:nseg)                = 0._r8
            xemax_u_msfc(0:mkx,:nseg)                = 0._r8
            cridis_u_msfc(0:mkx,:nseg)               = 0._r8
            thvcuenv_u_msfc(0:mkx,:nseg)             = 0._r8
            thvegenv_u_msfc(0:mkx,:nseg)             = 0._r8
            thvxsenv_u_msfc(0:mkx,:nseg)             = 0._r8
            fmix_u_msfc(0:mkx,:nseg)                 = 0._r8
            cmfumix_u_msfc(0:mkx,:nseg)              = 0._r8

            x_um_msfc(:mkx,:nseg)                    = 0._r8 
            y_um_msfc(:mkx,:nseg)                    = 0._r8
            x_p_msfc(0:mkx,:nseg)                    = 0._r8
            y_p_msfc(0:mkx,:nseg)                    = 0._r8

            thl_d_msfc(0:mkx,:nseg)                  = 0._r8
            qt_d_msfc(0:mkx,:nseg)                   = 0._r8
            u_d_msfc(0:mkx,:nseg)                    = 0._r8
            v_d_msfc(0:mkx,:nseg)                    = 0._r8
            w_d_msfc(0:mkx,:nseg)                    = 0._r8
            ql_d_msfc(0:mkx,:nseg)                   = 0._r8
            qi_d_msfc(0:mkx,:nseg)                   = 0._r8
            tr_d_msfc(0:mkx,:nseg,:ncnst)            = 0._r8
            wa_d_msfc(0:mkx,:nseg)                   = 0._r8
            qla_d_msfc(0:mkx,:nseg)                  = 0._r8
            qia_d_msfc(0:mkx,:nseg)                  = 0._r8
            cmf_d_msfc(0:mkx,:nseg)                  = 0._r8
            a_d_msfc(0:mkx,:nseg)                    = 0._r8

            slflx_u(0:mkx)                           = 0._r8
            qtflx_u(0:mkx)                           = 0._r8
            uflx_u(0:mkx)                            = 0._r8
            vflx_u(0:mkx)                            = 0._r8
            qlflx_u(0:mkx)                           = 0._r8
            qiflx_u(0:mkx)                           = 0._r8
            trflx_u(0:mkx,:ncnst)                    = 0._r8

            slten_u(:mkx)                            = 0._r8
            qtten_u(:mkx)                            = 0._r8
            uten_u(:mkx)                             = 0._r8
            vten_u(:mkx)                             = 0._r8
            sten_u(:mkx)                             = 0._r8
            qvten_u(:mkx)                            = 0._r8
            qlten_u(:mkx)                            = 0._r8
            qiten_u(:mkx)                            = 0._r8
            trten_u(:mkx,:ncnst)                     = 0._r8

            slflx_d(0:mkx)                           = 0._r8
            qtflx_d(0:mkx)                           = 0._r8
            uflx_d(0:mkx)                            = 0._r8
            vflx_d(0:mkx)                            = 0._r8
            qlflx_d(0:mkx)                           = 0._r8
            qiflx_d(0:mkx)                           = 0._r8
            trflx_d(0:mkx,:ncnst)                    = 0._r8

            slten_d(:mkx)                            = 0._r8
            qtten_d(:mkx)                            = 0._r8
            uten_d(:mkx)                             = 0._r8
            vten_d(:mkx)                             = 0._r8
            sten_d(:mkx)                             = 0._r8
            qvten_d(:mkx)                            = 0._r8
            qlten_d(:mkx)                            = 0._r8
            qiten_d(:mkx)                            = 0._r8
            trten_d(:mkx,:ncnst)                     = 0._r8

            slten_evp(:mkx)                          = 0._r8
            qtten_evp(:mkx)                          = 0._r8
            uten_evp(:mkx)                           = 0._r8
            vten_evp(:mkx)                           = 0._r8
            sten_evp(:mkx)                           = 0._r8
            qvten_evp(:mkx)                          = 0._r8
            qlten_evp(:mkx)                          = 0._r8
            qiten_evp(:mkx)                          = 0._r8
            trten_evp(:mkx,:ncnst)                   = 0._r8
            trten_wdep(:mkx,:ncnst)                  = 0._r8

            qrten(:mkx)                              = 0._r8
            qsten(:mkx)                              = 0._r8

            evapc(:mkx)                              = 0._r8

            qrten_u(:mkx)                            = 0._r8
            qsten_u(:mkx)                            = 0._r8

            qrten_u_msfc(:mkx,:nseg)                 = 0._r8
            qsten_u_msfc(:mkx,:nseg)                 = 0._r8
            trrsten_u_msfc(:mkx,:nseg,:ncnst)        = 0._r8

            qrten_d(:mkx)                            = 0._r8
            qsten_d(:mkx)                            = 0._r8

            qrten_d_msfc(:mkx,:nseg)                 = 0._r8
            qsten_d_msfc(:mkx,:nseg)                 = 0._r8
            trrsten_d_msfc(:mkx,:nseg,:ncnst)        = 0._r8

            snowmlt_e(:mkx)                          = 0._r8
            snowmlt_e_msfc(:mkx,:nseg)               = 0._r8

            thlten_dia_u(:mkx)                       = 0._r8       
            thlten_dia_d(:mkx)                       = 0._r8

            qtten_dia_u(:mkx)                        = 0._r8        
            qtten_dia_d(:mkx)                        = 0._r8

            qlten_dia_u(:mkx)                        = 0._r8        
            qlten_dia_d(:mkx)                        = 0._r8

            qiten_dia_u(:mkx)                        = 0._r8        
            qiten_dia_d(:mkx)                        = 0._r8

            trten_dia_u(:mkx,:ncnst)                 = 0._r8
            trten_dia_d(:mkx,:ncnst)                 = 0._r8

            ntraprd(:mkx)                            = 0._r8
            ntsnprd(:mkx)                            = 0._r8
            nttrrsprd(:mkx,:ncnst)                   = 0._r8

            ntraprd_msfc(:mkx,:nseg)                 = 0._r8
            ntsnprd_msfc(:mkx,:nseg)                 = 0._r8
            nttrrsprd_msfc(:mkx,:nseg,:ncnst)        = 0._r8

            evprain_e(:mkx)                          = 0._r8
            evpsnow_e(:mkx)                          = 0._r8
            evptrrs_e(:mkx,:ncnst)                   = 0._r8
            wdeptrrs_e(:mkx,:ncnst)                  = 0._r8

            evprain_d(:mkx)                          = 0._r8
            evpsnow_d(:mkx)                          = 0._r8
            evptrrs_d(:mkx,:ncnst)                   = 0._r8

            evprain_e_msfc(:mkx,:nseg)               = 0._r8
            evpsnow_e_msfc(:mkx,:nseg)               = 0._r8
            evptrrs_e_msfc(:mkx,:nseg,:ncnst)        = 0._r8
            wdeptrrs_e_msfc(:mkx,:nseg,:ncnst)       = 0._r8

            evprain_d_msfc(:mkx,:nseg)               = 0._r8
            evpsnow_d_msfc(:mkx,:nseg)               = 0._r8
            evptrrs_d_msfc(:mkx,:nseg,:ncnst)        = 0._r8

            cvp_rainprd(:mkx)                        = 0._r8
            cvp_snowprd(:mkx)                        = 0._r8
            cvp_trrsprd(:mkx,:ncnst)                 = 0._r8

            cvp_rainprd_msfc(:mkx,:nseg)             = 0._r8
            cvp_snowprd_msfc(:mkx,:nseg)             = 0._r8
            cvp_trrsprd_msfc(:mkx,:nseg,:ncnst)      = 0._r8

            qlten_eff_u(:mkx)                        = 0._r8
            qiten_eff_u(:mkx)                        = 0._r8

            qlten_eff_d(:mkx)                        = 0._r8
            qiten_eff_d(:mkx)                        = 0._r8

            trten_eff_u(:mkx,:ncnst)                 = 0._r8
            trten_eff_d(:mkx,:ncnst)                 = 0._r8

            qlten_par(:mkx)                          = 0._r8
            qiten_par(:mkx)                          = 0._r8
            qtten_par(:mkx)                          = 0._r8
            slten_par(:mkx)                          = 0._r8
            uten_par(:mkx)                           = 0._r8
            vten_par(:mkx)                           = 0._r8
            trten_par(:mkx,:ncnst)                   = 0._r8

            ql_env_ua(:mkx)                          = 0._r8
            qi_env_ua(:mkx)                          = 0._r8

            ql_env_da(:mkx)                          = 0._r8
            qi_env_da(:mkx)                          = 0._r8

            sten_NUM(:mkx)                           = 0._r8
            slten_NUM(:mkx)                          = 0._r8
            qtten_NUM(:mkx)                          = 0._r8
            uten_NUM(:mkx)                           = 0._r8
            vten_NUM(:mkx)                           = 0._r8
            qvten_NUM(:mkx)                          = 0._r8
            qlten_NUM(:mkx)                          = 0._r8
            qiten_NUM(:mkx)                          = 0._r8
            trten_NUM(:mkx,:ncnst)                   = 0._r8

            qlten(:mkx)                              = 0._r8
            qiten(:mkx)                              = 0._r8
            qvten(:mkx)                              = 0._r8
            sten(:mkx)                               = 0._r8
            uten(:mkx)                               = 0._r8
            vten(:mkx)                               = 0._r8
            trten(:mkx,:ncnst)                       = 0._r8

            N_up(0:mkx)                              = 0
            ptops(:mkx,:nseg)                        = 0._r8
            ztops(:mkx,:nseg)                        = 0._r8
            m_from_msfc(:mkx,:nseg)                  = 0 
            msfc_from_m(:mkx,:nseg)                  = 0      
            ktop_msfc(:nseg)                         = 0
            ptop_msfc(:nseg)                         = 0._r8
            ztop_msfc(:nseg)                         = 0._r8  

            thle_b(:mkx)                             = 0._r8
            qte_b(:mkx)                              = 0._r8
            tre_b(:mkx,:ncnst)                       = 0._r8
            ue_b(:mkx)                               = 0._r8
            ve_b(:mkx)                               = 0._r8
            we_b(:mkx)                               = 0._r8
            qle_b(:mkx)                              = 0._r8
            qie_b(:mkx)                              = 0._r8
            ssthle(:mkx)                             = 0._r8
            ssqte(:mkx)                              = 0._r8
            ssue(:mkx)                               = 0._r8 
            ssve(:mkx)                               = 0._r8
            ssqle(:mkx)                              = 0._r8
            ssqie(:mkx)                              = 0._r8
            sstre(:mkx,:ncnst)                       = 0._r8

            ! ---------------------------------------------------------------------------------------------------------------- !
            ! Nov.03.2012. REFINEMENT IS NECESSARY FOR TREATING ACCRETION OF CLOUD DROPLETS.                                   !
            !              The iteration loop for computing accretion of cloud droplet will start here.                        !
            !              However, correct initialization of several summed array variables may be necessary here again.      !
            !              I should be careful on this initialization.                                                         ! 
            ! ---------------------------------------------------------------------------------------------------------------- !

            ! --------------------------------------- !
            !                                         ! 
            ! Beginning of Updraft Computation Upward !
            !                                         ! 
            ! --------------------------------------- !

            do k = kiss + 1, mkx - 1    ! Here, 'k'  is a layer index.

               km = k - 1        ! Here, 'km' is a base interface index.

               ! ----------------------------------------- !
               ! Define environmental structure variables. ! 
               ! ----------------------------------------- !

               z_b         =      zs0(km)     
               z_t         =       zs0(k)

               p_b         =      ps0(km)
               p_t         =       ps0(k)

               dz_m        =       dz0(k)
               dp_m        =       dp0(k)
               
               exn_b       =    exns0(km)
               exn_t       =     exns0(k)          

               thl_b       =   thl0bot(k)

               qt_b        =    qt0bot(k)

               u_b         =     u0bot(k)

               v_b         =     v0bot(k)

               ql_b        =    ql0bot(k)

               qi_b        =    qi0bot(k)

               thv_b       =   thv0bot(k)
               thv_t       =   thv0top(k)         

               thvl_b      =  thvl0bot(k)
               thvl_t      =  thvl0top(k)         

               rho_b       =   rho0bot(k)
               rho_m       =      rho0(k)
               rho_t       =   rho0top(k)         

               do mt = 1, ncnst
                  tr_b(mt) = tr0bot(k,mt)
               enddo

               if( k .gt. 1 ) then
                  ! ------------------------------------- !
                  ! Below is applied for mixing downdraft !
                  ! ------------------------------------- !
                  mu        = mu_mix
                  tmp1      = ( 1._r8 - mu ) * thvl0top(km) + mu * thvl_b
                  tmp2      = ( 1._r8 - mu ) * thv0top(km)  + mu * thv_b
                  thvl_minE = min( min( thvl_b, thvl_t ), tmp1 ) 
                  thv_minE  = min( min( thv_b,  thv_t  ), tmp2 ) 
                ! Apr.15.2014. Modified formula for thv_minE to obtain a reasonable solution when the mean inversion exists.
                  thv_minE  = max( thv0top(km), tmp2 ) 
                ! Apr.15.2014. Modified formula for thv_minE
                  thvl_minE = thvl_minE + offset_minE
                  thv_minE  = thv_minE  + offset_minE
                  ! ----------------- !
                  ! Diagnostic Output ! 
                  ! ----------------- !
                  thv_b_out(i,km)   = thv_b
                  thv_t_out(i,km)   = thv_t    
                  thv_mt_out(i,km)  = tmp2 
                  thv_min_out(i,km) = thv_minE    
                  ! ----------------- !
                  ! Diagnostic Output ! 
                  ! ----------------- !
               else
                  thvl_minE  = -1.e8_r8 ! Always detrain downdraft in the lowest model layer after all diabatic forcings
                  thv_minE   = -1.e8_r8 ! Always detrain downdraft in the lowest model layer after all diabatic forcings
               endif

               ! -------------------------------------------------------- !
               !                                                          !
               ! Closure Conditions at the first interface k = 0 or k = 1 !
               ! Use 'rho_b,thl_b,...' not 'rho_m,thl_m'.                 !
               !                                                          !
               ! -------------------------------------------------------- !

               if( k .eq. kiss + 1 ) then    ! Here, 'k' is a layer-index.

                  d_alpha   = ( alpha_max - alpha_cri ) / nseg
                  sigma_w   =   sigma_wo + min( kw,   1.4142_r8 ) * sqrt( tke1 )
                  sigma_w   =   max( 1.e-1_r8, min( 1.0_r8, sigma_w ) ) 
                  sigma_thl =   sigfac * shflx / ( rho_b * cp * exn_b * sigma_w )
                  sigma_qt  =   sigfac * qflx  / ( rho_b *              sigma_w )
                  sigma_u   =   sigfac * taux  / ( rho_b *              sigma_w )
                  sigma_v   =   sigfac * tauy  / ( rho_b *              sigma_w )
                  do mt = 1, ncnst
                     sigma_tr(mt) = sigfac * aflx(mt)  / ( rho_b * sigma_w )
                  enddo
#ifdef MODAL_AERO
                  sigma_tr(lptr_dust_a_amode(modeptr_coarse)) = 0._r8
                  sigma_tr(lptr_dust_a_amode(modeptr_accum))  = 0._r8
                  sigma_tr(lptr_nacl_a_amode(modeptr_coarse)) = 0._r8
                  sigma_tr(lptr_nacl_a_amode(modeptr_accum))  = 0._r8
                  sigma_tr(lptr_nacl_a_amode(modeptr_aitken)) = 0._r8
#endif
                  
                  sigma_w   = sigma_wo + min( kw,   1.4142_r8 ) * sqrt( tke1 )
                  tmp1      = alpha_max 
                  sigma_w   = max( 1.e-1_r8, min( 1.0_r8, sigma_w ) ) 

            !lim  Oct.3. I commented out below limiter block. But instead, impose non-negative limiters later
            !lim         further below.  

            !lim  sigma_thl = max( - thl_b / tmp1, min( thl_b / tmp1 , sigma_thl ) )
                  sigma_thl = max( - 2.0_r8, min( 2.0_r8, sigma_thl ) )
            !lim  sigma_qt  = max( - qt_b / tmp1, min( qt_b / tmp1 , sigma_qt ) )
            !lim  ! sigma_qt  = max( - 5.e-4_r8, min( 1.e-3_r8, sigma_qt ) )
                  sigma_qt  = max( - 0.2_r8*qt_b, min( 0.2_r8*qt_b, sigma_qt ) )
                  sigma_u   = max( - 2.0_r8, min( 2.0_r8 , sigma_u ) ) 
                  sigma_v   = max( - 2.0_r8, min( 2.0_r8 , sigma_v ) )
            !lim  do mt = 1, ncnst
            !lim     sigma_tr(mt) = max( - tr_b(mt) / tmp1, min( tr_b(mt) / tmp1 , sigma_tr(mt) ) )
            !lim  enddo

                  cmfu_base = au_base * rho_b * ( sigma_w * sqrt( 2._r8 / 3.141592_r8 ) )
                  rnorm_a = 0._r8            
                  rnorm_m = 0._r8            
                  do m = 1, nseg
                     alpha(m)  = alpha_cri + d_alpha * ( m - 0.5_r8 )
                     if( nseg .eq. 1 ) alpha(m) = 1._r8
                     Pau(m)    = exp(-0.5_r8*alpha(m)**2._r8)/sqrt(2._r8*3.141592_r8)
                     Pau(m) = 2._r8 * au_base * Pau(m)
                     Pmu(m) = Pau(m) * rho_b * ( sigma_w * alpha(m) )

                     rnorm_a   = rnorm_a +  Pau(m) * d_alpha 
                     rnorm_m   = rnorm_m +  Pmu(m) * d_alpha
                  enddo

                  tmp1 = 0._r8
                  sum  = 0._r8            
                  do m = 1, nseg
                     alpha(m)  = alpha_cri + d_alpha * ( m - 0.5_r8 )
                     if( nseg .eq. 1 ) alpha(m) = 1._r8
                     Pau(m)    = exp(-0.5_r8*alpha(m)**2._r8)/sqrt(2._r8*3.141592_r8)
                     Pau(m) = 2._r8 * au_base * Pau(m)
                     if( inorm .eq. 1 ) then
                        Pau(m) = Pau(m) * ( au_base / rnorm_a )   
                     elseif( inorm .eq. 2 ) then
                        Pau(m) = Pau(m) * ( cmfu_base / rnorm_m )   
                     endif
                     Pmu(m)    = Pau(m) * rho_b * ( sigma_w * alpha(m) + delta_w_PBL * cdelta_w )
                     Pnu(m)    = Pau(m) / ( 3.141592_r8 * ( Ro + sigmaR * alpha(m) )**2._r8 )
                     a_au(m)   = Pau(m) * d_alpha                            ! Updraft fractional area          [ no unit ] 
                     cmf_au(m) = Pmu(m) * d_alpha                            ! Updraft mass flux                [ kg / s / m^2 ]
                     num_au(m) = Pnu(m) * d_alpha                            ! Number of updraft plumes         [ # / m^2 ]
                     rad_au(m) = Ro     + sigmaR   * alpha(m)                ! Physical radius of updraft plume [ m ]
                     sum       = sum    + cmf_au(m) 
                     tmp1      = tmp1   + a_au(m)
                  enddo

                  do m = 1, nseg

                     thl_au(m) = thl_b  + sigma_thl * alpha(m) + delta_thl_PBL * cdelta_s
                     qt_au(m)  = qt_b   + sigma_qt  * alpha(m) +  delta_qt_PBL * cdelta_s
                !lim
                     qt_au(m)  = max( qt_au(m), qmin(1) ) 
                !lim 
                     u_au(m)   = u_b    + sigma_u   * alpha(m) +   delta_u_PBL * cdelta_s
                     v_au(m)   = v_b    + sigma_v   * alpha(m) +   delta_v_PBL * cdelta_s
                     w_au(m)   =          sigma_w   * alpha(m) +   delta_w_PBL * cdelta_w
                     call conden( p_b, thl_au(m), qt_au(m), th, qv, ql, qi, qse, id_check )
                     ql_au(m)  = ql
                     qi_au(m)  = qi
                     thv_au(m) = th * ( 1._r8 + zvir * qv - ql - qi )
                     do mt = 1, ncnst
                        tr_au(m,mt) = tr_b(mt) + sigma_tr(mt) * alpha(m) + delta_tr_PBL(mt) * cdelta_s
                !lim 
                        tr_au(m,mt) = max( tr_au(m,mt), qmin(mt) )
                !lim 
                     enddo

                     ! Initialize 'S_b_ql_au(m) = S_b_qi_au(m) = 0' at surface.

                     S_b_ql_au(m) = 0._r8
                     S_b_qi_au(m) = 0._r8

                     ! -------------------------------------------------------------------------------------------------------------- !
                     ! Nov.28.2012. Imposing Consistency between 'droplet mass' and 'droplet number'                                  !
                     ! For the fully consistent treatment of microphysics process, it is important to                                 !
                     ! impose a full consistency between 'ql_au(m),tr_au(m,ixnumliq)' and 'qi_au(m),tr_au(m,ixnumice)'                !
                     ! here from the surface. Similar consistency will be imposed at all the interfaces above too as                  !
                     ! well as for convective downdraft.                                                                              !
                     ! Consistency should be imposed based on the 'ql_au(m),qi_au(m)' not on 'tr_au(m,ixnumliq),tr_au(m,ixnumice)'    !
                     ! since computation of 'droplet mass' is more reliable than 'droplet number' in the current situation.           !
                     ! If air is already saturated at surface with non-zero condensate, compute tr_au(m,ixnumliq) and                 ! 
                     ! tr_au(m,ixnumice) with a externally specified droplet radius. This can be done later using aerosol information !
                     ! or using stratiform cloud information at surface. However, this case will happen in a very rare way.           !
                     ! -------------------------------------------------------------------------------------------------------------- !

                     if( ql_au(m) .eq. 0._r8 ) then
                        tr_au(m,ixnumliq) = 0._r8
                     else
                        tr_au(m,ixnumliq) = ql_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq )
                     endif

                     if( qi_au(m) .eq. 0._r8 ) then
                        tr_au(m,ixnumice) = 0._r8
                     else
                        tr_au(m,ixnumice) = qi_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice )
                     endif

                     ! ----------------------------------------------------- !
                     ! Compute the original updraft segment index at surface !
                     ! ----------------------------------------------------- !

                     msfc_from_m(k,m) = m

                     thl_u(kiss)  = thl_u(kiss) +  thl_au(m) * cmf_au(m)
                     qt_u(kiss)   =  qt_u(kiss) +   qt_au(m) * cmf_au(m)
                     u_u(kiss)    =   u_u(kiss) +    u_au(m) * cmf_au(m)
                     v_u(kiss)    =   v_u(kiss) +    v_au(m) * cmf_au(m) 
                     w_u(kiss)    =   w_u(kiss) +    w_au(m) * cmf_au(m) 
                     wa_u(kiss)   =  wa_u(kiss) +    w_au(m) *   a_au(m)
                     ql_u(kiss)   =  ql_u(kiss) +   ql_au(m) * cmf_au(m)  
                     qi_u(kiss)   =  qi_u(kiss) +   qi_au(m) * cmf_au(m)  
                     do mt = 1, ncnst
                        tr_u(kiss,mt) = tr_u(kiss,mt) + tr_au(m,mt) * cmf_au(m)  
                     enddo
                     qla_u(kiss)  = qla_u(kiss) +   ql_au(m) *   a_au(m)  
                     qia_u(kiss)  = qia_u(kiss) +   qi_au(m) *   a_au(m)  
                     rad_u(kiss)  = rad_u(kiss) +  rad_au(m)**2._r8 * num_au(m) ! Effective plume radius [ m ]
                     cmf_u(kiss)  = cmf_u(kiss) +  cmf_au(m) 
                     a_u(kiss)    =   a_u(kiss) +    a_au(m)
                     num_u(kiss)  = num_u(kiss) +  num_au(m)
                     thva_u(kiss) =thva_u(kiss) +  thv_au(m) *   a_au(m)

                     ! ------------------------------------------------------------------------------------- !
                     ! Compute discretized surface flux being considered in the current convection scheme.   !
                     ! Below discrete computation is fully appropriate and should be used since the input    !
                     ! surface flux can be unreasonably large for some cases.                                !
                     ! Sep.12.2011. I double checked that my below 'ipartition' formula is perfect because   !
                     !              (1) it conserves column-integrated energy, and (2) it completely remove  !
                     !              the generation of unreasonable convective tendency in the lowest model   !
                     !              layer by convection, so that it grauantees computation of reasonable     ! 
                     !              surface heat, moisture, momentum, and tracer fluxes at surface in the    !  
                     !              following surface flux computation routine in the CAM. Also, I don't     !
                     !              need to modify any parts of CAM5 ( e.g., PBL scheme, surface flux        !
                     !              routine ), since all the required modifications are contained in the     !
                     !              UNICON in a fully reasonable way.                                        !
                     !              By using 'ipartition = 1' option, I don't need to combine 'symmetric     !
                     !              moist turbulence scheme' with the 'asymmetric moist turbulence scheme'   !
                     !              within the implicit iteration loop, so that I can save tremendous amount !    
                     !              of computation time.                                                     !
                     ! Mar.19.2014. I added 'qlflx_u(0),qiflx_u(0)' in the below lines for fully consistent  ! 
                     !              treatment in association with 'ipartition=1' and with the same tratment  !
                     !              for convective downdraft flux later.                                     !  
                     ! ------------------------------------------------------------------------------------- !

                        slflx_u(0) = slflx_u(0) + cp * exns0(0) * cmf_au(m) * ( thl_au(m)   -    thl_b )
                        qtflx_u(0) = qtflx_u(0) +                 cmf_au(m) * ( qt_au(m)    -     qt_b )
                        uflx_u(0)  =  uflx_u(0) +                 cmf_au(m) * ( u_au(m)     -      u_b )
                        vflx_u(0)  =  vflx_u(0) +                 cmf_au(m) * ( v_au(m)     -      v_b )
                        qlflx_u(0) = qlflx_u(0) +                 cmf_au(m) * ( ql_au(m)    -     ql_b )
                        qiflx_u(0) = qiflx_u(0) +                 cmf_au(m) * ( qi_au(m)    -     qi_b )
                        do mt = 1, ncnst
                           trflx_u(0,mt) = trflx_u(0,mt) +        cmf_au(m) * ( tr_au(m,mt) - tr_b(mt) )
                        enddo

                  enddo  ! do m = 1, nseg


                  ! Mean convective updraft values at the launching interface
                  ! By construction, 'cmf_u(kiss),a_u(kiss),num_u(kiss)' are non-zero.

                  thl_u(kiss)  = thl_u(kiss) / cmf_u(kiss)
                  qt_u(kiss)   =  qt_u(kiss) / cmf_u(kiss)
                  u_u(kiss)    =   u_u(kiss) / cmf_u(kiss)
                  v_u(kiss)    =   v_u(kiss) / cmf_u(kiss)
                  w_u(kiss)    =   w_u(kiss) / cmf_u(kiss)
                  wa_u(kiss)   =  wa_u(kiss) /   a_u(kiss)
                  ql_u(kiss)   =  ql_u(kiss) / cmf_u(kiss)
                  qi_u(kiss)   =  qi_u(kiss) / cmf_u(kiss)
                  do mt = 1, ncnst
                     tr_u(kiss,mt) = tr_u(kiss,mt) / cmf_u(kiss)
                  enddo
                  qla_u(kiss)  = qla_u(kiss) /   a_u(kiss)
                  qia_u(kiss)  = qia_u(kiss) /   a_u(kiss)
                  rad_u(kiss)  = sqrt( rad_u(kiss) / num_u(kiss) )    ! Effective plume radius [ m ]
                  N_up(kiss)   = nseg                                 ! Total number of updraft plumes at the launching interface
                  thva_u(kiss) =thva_u(kiss) /   a_u(kiss)            ! Mean potential temperature within convective updraft.

                  ! ---------------------------------------------------------- !
                  ! Diagnostic Output for Checking Final Source Air Properties !
                  ! ---------------------------------------------------------- !


                  kw_out(i)        = kw

                  sigma_w_out(i)   = sigma_w
                  sigma_thl_out(i) = sigma_thl
                  sigma_qt_out(i)  = sigma_qt * 1.e3_r8
                  sigma_u_out(i)   = sigma_u
                  sigma_v_out(i)   = sigma_v 

                  tkes_out(i)    = tke1
                  w_org_out(i)   = delta_w_PBL   * cdelta_w
                  thl_org_out(i) = delta_thl_PBL * cdelta_s
                  qt_org_out(i)  = delta_qt_PBL  * cdelta_s
                  u_org_out(i)   = delta_u_PBL   * cdelta_s
                  v_org_out(i)   = delta_v_PBL   * cdelta_s

                  ! ---------------------------------------------------------- !
                  ! Diagnostic Output for Checking Final Source Air Properties !
                  ! ---------------------------------------------------------- !             

               else      ! 'if( k .eq. 1 )'

                  N_up(km)     = nseg_nondet 

               endif     ! 'if( k .eq. 1 )'
 
               ! ---------------------------------------------------------------------------------------------------------- !
               ! Mar.07.2013.                                                                                               !
               ! Compute grid-mean (which includes both environment and convective updraft) virtual potential temperature   !
               ! at the base interface for use in computing vertical evolution of vertical velocity of convective           !
               ! updraft and downdraft for developing a full scale-adaptive parameterization.                               !
               ! Note that 'thv_mean_t' will be computed later in a separate way depending on 'thv_au(m)' for each updraft. !       
               ! ---------------------------------------------------------------------------------------------------------- !

               thv_mean_b = a_u(km) * thva_u(km) + ( 1._r8 - a_u(km) ) * thv_b

               ! --------------------------------------------------------------------------- !
               ! Environmental airs involved in buoyancy sorting via mesoscale organization. !
               ! This can be mixtures of ensemble-mean environmental air and detrained airs. !
               ! This is used only for the mixing purpose not for diabatic buoyancy forcing. !
               ! Assuming each updraft equally reaches to dp0(k), solve discrete ( not       !
               ! continuous ) mass flux equation to compute total amount of updraft mass that!
               ! is involved in the buoyancy sorting.                                        !   
               ! The available mass per unit time per unit area [ kg/s/m^2 ]                 !
               !     Environmental airs : ae(k)*dp_m/g/dt                                    !
               !     Detrained airs     : cmf_r_org(k)                                       !
               !     The updraft mass involved in the buoyancy sorting : sum                 !
               ! --------------------------------------------------------------------------- !

               ! --------------------------------------------------------------------------------------------- !
               ! Compute (1) fractional mixing rate, eps0 and (2) organization factor for lateral entrainment. !
               !   . sum :     The estimate of updraft mass flux that will be involved in                      ! 
               !               buoyancy sorting at the current time step. [ kg/m^2/s ]                         !
               !   . cu_cmfr : Total amount of detrained mass at the previous time [ kg/m^2/s ]                !
               ! --------------------------------------------------------------------------------------------- !

               sum = 0._r8
               do m = 1, N_up(km)
                  ! -------------------------------------------------------------------- !
                  ! Fractional Mixing Rates eps0 [ 1/Pa ] : This is the key of the model !
                  ! -------------------------------------------------------------------- !

                  tmp1 = sqrt( ( ( ql_au(m) + qi_au(m) ) / 1.e-3_r8 ) * ( 1._r8 - max( 0.0_r8, min( 1.00_r8, rh0bot(k) ) ) ) )
                  tmp2 = 1._r8 + cevpeps0 * tmp1
                  eps0(m) = tmp2 * c0 / max( rad_au(m), Ro_eps0 ) / ( rho_b * g )
                  eps0(m) = min( eps0(m), log( 1._r8 + epsz0_max * dz_m ) / dp_m )
                  sum     = sum + dp_m * eps0(m) * cmf_au(m) 
               enddo

               ! ----------------------------------------------------------------------------------------- !
               ! Aug.01.2011. Brian Juwon Park's 10th Birthday.                                            !
               !             'cuorg' is replaced by 'cuorg_mxen' in the below line computing 'org_ent' for ! 
               !              treating mixing with several mixing environmental airs.                      !
               !              cuorg_mxen = 0 ( iter = 1 ) or 1 ( iter = 2 ).                               !
               ! ----------------------------------------------------------------------------------------- !

               org_ent = min( fac_org_ent * cuorg_mxen , cu_cmfr(k) / max( sum, nonzero ) ) ! The second argument sets the upper limit of org.
               org_ent = max( 0._r8, min( 1._r8, org_ent ) ) 

               ! ----------------------------------------------- !
               ! Define environmental properties using 'org_ent' !
               ! ----------------------------------------------- !

               tmp1    = org_ent

               if( k .lt. kpblh ) then

                  thle_b(k)       =  max( 0._r8,    thl_b +               cdelta_s * delta_thl_PBL )
                     
                  qte_b(k)        =  max( qmin(1),  qt_b +                cdelta_s * delta_qt_PBL )

                  ue_b(k)         =     (             u_b +                 cdelta_s * delta_u_PBL )

                  ve_b(k)         =     (             v_b +                 cdelta_s * delta_v_PBL )

                  we_b(k)         =     (                                   cdelta_w * delta_w_PBL )

                  qle_b(k)        =  max( 0._r8,     ql_b +                                  0._r8 )

                  qie_b(k)        =  max( 0._r8,     qi_b +                                  0._r8 )

                  do mt = 1, ncnst
                     tre_b(k,mt)  =  max( qmin(mt), tr_b(mt) +         cdelta_s * delta_tr_PBL(mt) )
                  enddo

                  ssthle(k) =  ssthl0(k)
                  ssqte(k)  =   ssqt0(k)
                  ssue(k)   =    ssu0(k)
                  ssve(k)   =    ssv0(k)
                  ssqle(k)  =   ssql0(k)
                  ssqie(k)  =   ssqi0(k)
                  do mt = 1, ncnst
                     sstre(k,mt) = sstr0(k,mt)
                  enddo


               elseif( k .ge. kpblh ) then

                  thle_b(k)       =  max( 0._r8,    thl_b + tmp1 *   cu_thlr(k) )

                  qte_b(k)        =  max( qmin(1),   qt_b + tmp1 *    cu_qtr(k) )

                  ue_b(k)         =     (             u_b + tmp1 *     cu_ur(k) )

                  ve_b(k)         =     (             v_b + tmp1 *     cu_vr(k) )

                  we_b(k)         =  0._r8

                  qle_b(k)        =  max( 0._r8,     ql_b + tmp1 *    cu_qlr(k) )

                  qie_b(k)        =  max( 0._r8,     qi_b + tmp1 *    cu_qir(k) )

                  do mt = 1, ncnst
                     tre_b(k,mt)  =  max( qmin(mt), tr_b(mt) + tmp1 * cu_trr(k,mt) )
                  enddo

                  ssthle(k) =  ssthl0(k) * ( 1._r8 - tmp1 )
                  ssqte(k)  =   ssqt0(k) * ( 1._r8 - tmp1 )
                  ssue(k)   =    ssu0(k) * ( 1._r8 - tmp1 )
                  ssve(k)   =    ssv0(k) * ( 1._r8 - tmp1 )
                  ssqle(k)  =   ssql0(k) * ( 1._r8 - tmp1 )
                  ssqie(k)  =   ssqi0(k) * ( 1._r8 - tmp1 )
                  do mt = 1, ncnst
                     sstre(k,mt) = sstr0(k,mt) * ( 1._r8 - tmp1 ) 
                  enddo

               endif ! End of 'if( k .lt. kpblh )' and 'elseif( k .ge. kpblh)' blocks to define the properties of mixing environmental airs.

               ! ------------------------------------ !
               !                                      !
               ! Individual Updraft Plume Computation !
               !                                      ! 
               ! ------------------------------------ !

               ! ---------------------------------------- !
               ! Initialization of updraft segment arrays !
               ! ---------------------------------------- !

               ytop(:nseg)               = 0._r8       ! 0 : Not reach the top interface, 1 : Reach the top interface 
               xc(:nseg)                 = 0._r8       
               xs(:nseg)                 = 0._r8       
               eeps(:nseg)               = 0._r8
               ddel(:nseg)               = 0._r8
               eps(:nseg)                = 0._r8
               del(:nseg)                = 0._r8
               xe_min(:nseg)             = 0._r8
               xe_max(:nseg)             = 0._r8
               dpa(:nseg)                = 0._r8
               dza(:nseg)                = 0._r8
               ptop(:nseg)               = 0._r8
               ztop(:nseg)               = 0._r8
               fmix(:nseg)               = 0._r8
               f_wu(:nseg)               = 0._r8

               thl_aut(:nseg)            = 0._r8
               qt_aut(:nseg)             = 0._r8
               u_aut(:nseg)              = 0._r8
               v_aut(:nseg)              = 0._r8
               ql_aut(:nseg)             = 0._r8
               qi_aut(:nseg)             = 0._r8
               cmf_aut(:nseg)            = 0._r8
               w_aut(:nseg)              = 0._r8
               a_aut(:nseg)              = 0._r8
               num_aut(:nseg)            = 0._r8
               rad_aut(:nseg)            = 0._r8
               tr_aut(:nseg,:ncnst)      = 0._r8
               thv_aut(:nseg)            = 0._r8
               S_t_ql_au(:nseg)          = 0._r8 
               S_t_qi_au(:nseg)          = 0._r8 

               evp_thll_au(:nseg)        = 0._r8         
               evp_qtl_au(:nseg)         = 0._r8
               evp_thli_au(:nseg)        = 0._r8         
               evp_qti_au(:nseg)         = 0._r8
               prep_thll_au(:nseg)       = 0._r8         
               prep_qtl_au(:nseg)        = 0._r8
               prep_thli_au(:nseg)       = 0._r8         
               prep_qti_au(:nseg)        = 0._r8
               eff_ql_au(:nseg)          = 0._r8
               eff_qi_au(:nseg)          = 0._r8
               PGF_u_au(:nseg)           = 0._r8
               PGF_v_au(:nseg)           = 0._r8
               evp_tr_au(:nseg,:ncnst)   = 0._r8
               prep_tr_au(:nseg,:ncnst)  = 0._r8         
               eff_tr_au(:nseg,:ncnst)   = 0._r8

               do m = 1, N_up(km)

                  ! -------------------------------------------------- !
                  ! Compute R-dependent buoyancy coefficient, rbuoy_up !
                  ! -------------------------------------------------- !

                  rbuoy_up = rbuoy_min * ( 1._r8 + ( rbuoy_max / rbuoy_min - 1._r8 ) * exp( - rad_au(m) / R_buo ) )

                  ! --------------------------------------------------------------------------- !
                  ! Updraft Buoyancy Sorting                                                    !
                  ! --------------------------------------------------------------------------- !

                  pe        =  p_b
                  w_cu      =  w_au(m)
                  thl_cu    =  thl_au(m)
                  qt_cu     =  qt_au(m)
                  thl_eg    =  thle_b(k) 
                  qt_eg     =  qte_b(k)
                  u_eg      =  ue_b(k)        ! Not used here, but for consistent treatment of buoyancy sorting 
                  v_eg      =  ve_b(k)        ! with the detrained airs : 'thl_dd, etc.' later
                  w_eg      =  we_b(k) 
                  do mt = 1, ncnst
                     tr_eg(mt) = tre_b(k,mt)  ! Not used here in this buoyancy sorting but for defining tr_eg(mt). 
                  enddo
                  thv_env   =  thv_b
                  cridis    =  rlc * cushavg                       ! [ m ]. New for positive feedback. 
                  if( rlc .eq. -1._r8 ) cridis = dz_m 
                  call conden( pe, thl_cu, qt_cu, th, qv, ql, qi, qse, id_check )
                  ql_cu = ql
                  qi_cu = qi
                  thv  = th * ( 1._r8 + zvir * qv - ql - qi )

                  ! --------------------------------------------------------------------------- !
                  ! Mar.11.2013. Add 'w_eg' within PBL for treating the effect of organized flow!
                  !              on the buoyancy sorting.                                       !
                  !              Note that 'w_eg' is added only for 'buosorts_UW' subroutine,   !
                  !              so that the other buoyancy sorting subroutines should be       !
                  !              modified in future if I want to use them.                      ! 
                  ! --------------------------------------------------------------------------- ! 

                  call buosorts_UW(  rbuoy_up, pe, w_cu, thl_cu, qt_cu, w_eg, thl_eg, qt_eg, &
                     thv_env, cridis, xc(m), xs(m), thv_cu, thv_eg, thv_xs )

                  ! -------------------------------------------------------- !
                  ! We can impose a lower and upper limit on xc(m) as below. !
                  ! -------------------------------------------------------- !
                  xc(m) = max( xc_min, min( xc_max, xc(m) ) )
                  ! -------------------------------------------------------- !
                  ! Compute non-dimensional entrainment and detrainment rate !
                  ! -------------------------------------------------------- !
                  call compute_epsdelnod( 'PDFbsQ', xc(m), eeps(m), ddel(m) )
                  ! -------------------------------------------------------------------- !
                  !                                                                      !
                  ! Fractional Mixing Rates eps0 [ 1/Pa]  : This is the key of the model !
                  ! Aug.15.2011. In contrast to the computation in the above block for   !
                  !              limiting org_ent, precise value of eps0(m) should be    !
                  !              computed here.                                          ! 
                  !                                                                      !
                  ! -------------------------------------------------------------------- !  

                  call conden( pe, thl_eg, qt_eg, th, qv, ql, qi, qse, id_check )
                  rh_eg = max( 0._r8, min( 1._r8, qv / max( nonzero, qse ) ) )
                  tmp1 = sqrt( ( ( ql_au(m) + qi_au(m) ) / 1.e-3_r8 ) * ( 1._r8 - max( 0.0_r8, min( 1.00_r8, rh_eg ) ) ) )
                  tmp2 = 1._r8 + cevpeps0 * tmp1
                  eps0(m) = tmp2 * c0 / max( rad_au(m), Ro_eps0 ) / ( rho_b * g )

                  ! -------------------------------------------------------------------- !
                  !                                                                      !
                  ! Fractional Mixing Rates eps0 [ 1/Pa ] : This is the key of the model !
                  !                                                                      !
                  ! -------------------------------------------------------------------- !  

                  eps0(m)   =  max( 0._r8, eps0(m) )

                  ! ----------------------------------------------------------------------------------------------------------------------- !                    
                  ! CRITICALLY IMPORTANT : IMPOSING A REASONABLE LIMIT FOR eps0(m)                                                          !
                  ! ----------------------------------------------------------------------------------------------------------------------- !
                  
                  ! Feb.06.2013. Always choose physically reasonable exp_cmf = 1 as of this day.

                  eps0(m) = min( eps0(m), log( 1._r8 + epsz0_max * dz_m ) / dp_m )
                  eps(m)    = eps0(m) * eeps(m)
                  del(m)    = eps0(m) * ddel(m)

                  ! ----------------------------------------------------------------------------- !
                  ! Compute which mixtures can be the source of downdraft ( xe_min < x < xe_max ) !
                  ! ----------------------------------------------------------------------------- !
       
                  thl_cumC  = thl_cu  + xc(m) * ( thl_eg - thl_cu )
                  qt_cumC   = qt_cu   + xc(m) * ( qt_eg  -  qt_cu )
                  call conden( pe, thl_cumC, qt_cumC, th, qv, ql, qi, qse, id_check )
                  thv_cumC  = th *  ( 1._r8 + zvir * qv - ql - qi )

                  thl_cumS  = thl_cu  + xs(m) * ( thl_eg - thl_cu )
                  qt_cumS   = qt_cu   + xs(m) * ( qt_eg  -  qt_cu )
                  call conden( pe, thl_cumS, qt_cumS, th, qv, ql, qi, qse, id_check )
                  thv_cumS  = th *  ( 1._r8 + zvir * qv - ql - qi )

                  call conden( pe, thl_eg, qt_eg, th, qv, ql, qi, qse, id_check )
                  thv_cumE  = th *  ( 1._r8 + zvir * qv - ql - qi )

                  if( xc(m) .ge. xs(m) ) then
                     call buosort_downdraft( thv_cumC, thv_cumE, thv_minE, xdown_min, xdown_max )
                     xe_min(m) = xc(m) + ( 1._r8 - xc(m) ) * xdown_min
                     xe_max(m) = xc(m) + ( 1._r8 - xc(m) ) * xdown_max     
                  else
                     call buosort_downdraft( thv_cumC, thv_cumS, thv_minE, xdown_min, tmp1 )
                     xe_min(m) = xc(m) + ( xs(m) - xc(m) ) * xdown_min
                     call buosort_downdraft( thv_cumS, thv_cumE, thv_minE, tmp2, xdown_max )
                     xe_max(m) = xs(m) + ( 1._r8 - xs(m) ) * xdown_max
                  endif

                  ! ---------------------------------------------------------------------------------- !
                  ! Updraft Top Height & Vertical Velocity at the Top Interface                        ! 
                  ! ---------------------------------------------------------------------------------- !
                  
                  thvbot  = thv_au(m) 
                  bogbot  = rbuoy_up * ( thvbot / thv_mean_b - 1._r8 )

                  ! ------------------------------------------------------------------------- !
                  ! In this case, entrainment mixing occurs. So, simply use the previous code !
                  ! by assuming a simple linear profile of buoyancy from the base to the top  !
                  ! interface.                                                                !
                  ! Apr.17.2012. In order to remove ambiguity, I am using thl_meu, qt_meu     !
                  !              in the below block.                                          !
                  ! ------------------------------------------------------------------------- !
                  thl_meu = thle_b(k) + ssthle(k) * 0.5_r8 * ( p_t - p_b )
                  qt_meu  = qte_b(k)  + ssqte(k)  * 0.5_r8 * ( p_t - p_b )
                  ql_meu  = qle_b(k)  + ssqle(k)  * 0.5_r8 * ( p_t - p_b )
                  qi_meu  = qie_b(k)  + ssqie(k)  * 0.5_r8 * ( p_t - p_b )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, p_t, thl_meu, ssthle(k), thl_au(m), thl_aut_tmp )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, p_t,  qt_meu,  ssqte(k),  qt_au(m),  qt_aut_tmp )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, p_t,  ql_meu,  ssqle(k),  ql_au(m),  ql_aut_adi )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, p_t,  qi_meu,  ssqie(k),  qi_au(m),  qi_aut_adi )
                  call conden( p_t, thl_aut_tmp, qt_aut_tmp, th, qv, ql, qi, qse, id_check )
                  do mt = 1, ncnst
                     tr_meu(mt)  = tre_b(k,mt)  + sstre(k,mt)  * 0.5_r8 * ( p_t - p_b )
                     call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, p_t, tr_meu(mt), sstre(k,mt), tr_au(m,mt), tr_aut_tmp(mt) )
                  enddo

                  ! ------------------------------------------------------------------------------------ !
                  ! Feb.07.2013.                                                                         !
                  ! Compute precipitation production at the top interface.                               !
                  ! ------------------------------------------------------------------------------------ !

                  ! ------------------------------------------------------------------------- !
                  ! Compute 'exql,exqi' by solving analytical vertical integration of 'ql,qi' ! 
                  ! by including differential precipitation fall-out in the integration.      !
                  ! In the below simplified microphysics, 'eps_dia_L = eps_dia_I'. However,   !
                  ! in the future refined microphysics, they can differ. Thus, I am keeping   !
                  ! both 'eps_dia_L' and 'eps_dia_I'.                                         !   
                  ! ------------------------------------------------------------------------- !

                  if( iprd_prep .eq. 5 ) then  ! Backward Analytical Method

                     if( ( ql_cu + qi_cu ) .gt. criqc ) then
                        eps_dia_L = c0_ac * ( ( ql_cu + qi_cu ) - criqc ) / ( ql_cu + qi_cu )
                        eps_dia_I = c0_ac * ( ( ql_cu + qi_cu ) - criqc ) / ( ql_cu + qi_cu )
                     else
                        eps_dia_L = 0._r8
                        eps_dia_I = 0._r8
                     endif
                     call progup_thlqt( eps(m), eps_dia_L, 0._r8, p_b, p_t, ql_meu, ssqle(k), ql_au(m), ql_aut_adi_prp )
                     call progup_thlqt( eps(m), eps_dia_I, 0._r8, p_b, p_t, qi_meu, ssqie(k), qi_au(m), qi_aut_adi_prp )
                     exql = min( max( ql_aut_adi - ql_aut_adi_prp, 0._r8 ), 0.99_r8 * ql )  ! This should be guaranteed to be positive at this stage.
                     exqi = min( max( qi_aut_adi - qi_aut_adi_prp, 0._r8 ), 0.99_r8 * qi )  ! This should be guaranteed to be positive at this stage.
                     if( mclimit .eq. 1 ) then
                        tmp1 = exql + exqi
                        tmp2 = min( tmp1, max( ql + qi - criqc, 0._r8 ) )                 ! To impose a continuous variation across ql + qi = criqc.
                        exql = exql * ( tmp2 / max( tmp1, nonzero ) )
                        exqi = exqi * ( tmp2 / max( tmp1, nonzero ) )
                     endif
                     ! -------------------------------------------------- !
                     ! Evaporation within Updraft.                        !
                     ! Set it to be zero, but can be refined in future.   !
                     ! -------------------------------------------------- !
                     evpR = 0._r8
                     evpS = 0._r8

                  else ! The others of Backward Analytical Method

                     ! ------------------------------------------------------------------------- !
                     ! Below is previous block replaced by the above full                        !
                     ! analytical treatment of precipitation fall-out within convective updraft. !
                     ! ------------------------------------------------------------------------- !

                     msfc = msfc_from_m(k,m)
                     call prod_prep_up( z_b, z_t, p_b, p_t, exn_t, exn0(k),                         &
                        w_au(m), w_au(m),                                           &
                        thl_aut_tmp, qt_aut_tmp, ql, qi, tr_aut_tmp(1:ncnst),       &
                        S_b_ql_au(m), S_b_qi_au(m), iprd_prep,                      &
                        ql_cu, qi_cu, eps(m),                                       &
                        thl_meu, ssthle(k), thl_au(m), qt_meu, ssqte(k), qt_au(m),  & 
                        ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice, i, k, lchnk, &
                        flxrain_msfc(k,msfc), flxsnow_msfc(k,msfc),                 &     
                        a_p_msfc(k,msfc), am_u_msfc(k,msfc), am_pu_msfc(k,msfc),    & 
                        caer, criqc, c0_ac,                                         &
                        exql, exqi, extr(1:ncnst), S_t_ql_au(m), S_t_qi_au(m),      &
                        evpR, evpS, evpRStr(1:ncnst) )

                  endif ! End of Backward Analytical Method

                  ! ----------------------------------------------------------------------------------- !
                  ! Jun.16.2012. I should recompute the buoyancy using updated state variables as below !
                  !              similar to CAM5 shallow convection scheme. It turns out that           !
                  !              this update is important and has non-negligible impact on the          !
                  !              simulation. Comment-out above two lines.                               !
                  !              For simplicity, use exner function defined at the top interface since  !
                  !              p_t is defined at the top model interface.                             !
                  !              The correct cumulus top is computed later below.                       ! 
                  !              It turns out that computation of 'thvtop' has huge influences on the   !
                  !              CLDLOW and SWCF in the trade cumulus regime.                           !     
                  ! ----------------------------------------------------------------------------------- ! 
                  call conden( p_t, thl_aut_tmp + ( xlv / cp / exn_t ) * ( exql - evpR ) + ( xls / cp / exn_t ) * ( exqi - evpS ), &
                                  qt_aut_tmp - exql - exqi + evpR + evpS, &
                                  th, qv, ql, qi, qse, id_check )
                  thvtop  = th * ( 1._r8 + zvir * qv - ql - qi )
                  thv_mean_t = a_u(km) * thvtop + ( 1._r8 - a_u(km) ) * thv_t
                  bogtop  = rbuoy_up * ( thvtop / thv_mean_t - 1._r8 )

                  ! --------------------------------------------------------------------------------------------- !
                  ! Below block is generally formulated to treat the special case of 'bogbot < 0 and bogtop > 0'. !
                  ! So, even when 'OPTION.3' is selected above, below block itself can treat all the cases in a   !
                  ! most reasonable way. So, the CIN structure ( other than resolving LCL explicitly ) can be     !
                  ! resolved in a most reasonable way.                                                            !
                  ! Mar.11.2013. Add 'we_b(k)**2._r8' as the argument of 'progup_wu2' and 'compute_dp'.           !
                  ! --------------------------------------------------------------------------------------------- !

                  if( bogbot .lt. 0._r8 .and. bogtop .gt. 0._r8 ) then 
                     plfc = p_b - ( p_b - p_t ) * ( bogbot / ( bogbot - bogtop ) )
                     call progup_wu2( rdrag*eps(m) - rjet*del(m), rho_m, p_b, plfc, bogbot, 0._r8, w_au(m)**2._r8, &
                                      we_b(k)**2._r8, wu2 ) 
                     if( wu2 .ge. 0._r8 ) then
                        dpa(m) = p_b - p_t
                        call progup_wu2( rdrag*eps(m) - rjet*del(m), rho_m, p_b, p_t, bogbot, bogtop, w_au(m)**2._r8, &
                                         we_b(k)**2._r8, wu2 ) 
                        w_aut(m) = max( sqrt( max( wu2, nonzero ) ), wumin )       
                     else
                        dpa(m) = compute_dp( rdrag*eps(m) - rjet*del(m), rho_m, p_b, plfc, bogbot, 0._r8, w_au(m)**2._r8, &
                                 we_b(k)**2._r8 )  ! '0 <= dpa(m) <= p_b-plfc'
                        w_aut(m) = 0._r8
                     endif
                  else
                     call progup_wu2( rdrag*eps(m) - rjet*del(m), rho_m, p_b, p_t, bogbot, bogtop, w_au(m)**2._r8, &
                                      we_b(k)**2._r8, wu2 ) 
                     if( wu2 .ge. 0._r8 ) then
                        dpa(m) = p_b - p_t
                        w_aut(m) = max( sqrt( max( wu2, nonzero ) ), wumin )       
                     else
                        if( .not. ( ( bogbot .ge. 0._r8 .and. bogtop .lt. 0._r8 ) .or. &
                           ( bogbot .lt. 0._r8 .and. bogtop .le. 0._r8 ) ) ) then  
                           write(iulog,*) 'bogbot, bogtop = ', bogbot, bogtop
                           write(iulog,*) 'awk_PBL_raw, delta_thl_PBL_raw, delta_qt_PBL_raw, delta_u_PBL_raw, delta_v_PBL_raw = ', & 
                                         awk_PBL_raw, delta_thl_PBL_raw, delta_qt_PBL_raw, delta_u_PBL_raw, delta_v_PBL_raw
                           write(iulog,*) 'awk_PBL, delta_thl_PBL, delta_qt_PBL, delta_u_PBL, delta_v_PBL, delta_w_PBL = ', & 
                                         awk_PBL, delta_thl_PBL, delta_qt_PBL, delta_u_PBL, delta_v_PBL, delta_w_PBL
                           call endrun('UNICON : Impossible buoyancy case before compute_dp')  
                        endif
                        if( bogbot .ge. 0._r8 .and. bogtop .lt. 0._r8 ) then 
                           plnb = p_b - ( p_b - p_t ) * ( bogbot / ( bogbot - bogtop ) )
                           call progup_wu2( rdrag*eps(m) - rjet*del(m), rho_m, p_b, plnb, bogbot, 0._r8, w_au(m)**2._r8, &
                                            we_b(k)**2._r8, tmp1 )
                           dpa(m) = compute_dp( rdrag*eps(m) - rjet*del(m), rho_m, plnb, p_t, 0._r8, bogtop, tmp1, &
                                                we_b(k)**2._r8 )  ! '0 <= dpa(m) <= plnb-p_t'
                           dpa(m) = dpa(m) + ( p_b - plnb )
                        else
                           dpa(m) = compute_dp( rdrag*eps(m) - rjet*del(m), rho_m, p_b, p_t, bogbot, bogtop, w_au(m)**2._r8, &
                                                we_b(k)**2._r8 )  ! '0 <= dpa(m) <= p_b-p_t'
                        endif
                        w_aut(m) = 0._r8
                     endif
                  endif

                  dza(m)  = min( dz_m, max( 0._r8, dpa(m) / rho_m / g ) )
                  ptop(m) = p_b - dpa(m)
                  ztop(m) = z_b + dza(m)

            
                  tmp1 = cmf_au(m) * exp( dpa(m) * ( eps(m) - del(m) ) )

                  if( w_aut(m) .gt. nonzero .and. tmp1 .ge. cmfmin .and. k .ne. mkx ) then 
                     ytop(m)    = 1._r8
                     cmf_aut(m) = cmf_au(m) * exp( dp_m * ( eps(m) - del(m) ) )
                     a_aut(m)   = cmf_aut(m) / rho_t / w_aut(m)
                     num_aut(m) = num_au(m)
                     rad_aut(m) = sqrt( a_aut(m) / num_aut(m) /  3.141592_r8 )   ! Physical radius of updraft plume [ m ]

                  else

                     ytop(m)    = 0._r8
                     w_aut(m)   = 0._r8

                     cmf_aut(m) = cmf_au(m) * exp( dpa(m) * ( eps(m) - del(m) ) ) ! For the purpose of computing diabatic forcing later, retain this.
                     a_aut(m)   = a_au(m)
                     num_aut(m) = num_au(m)
                     rad_aut(m) = rad_au(m)

                  endif

                  ! ------------------------------------------------------------------ !
                  ! fmix(m) * dpa(m) * eps0(m) * cmf_au(m) :                           !
                  ! The amount of updraft mass involved in the buoyancy sorting mixing !
                  ! Oct.05.2010. This 'fmix' can be enormously large ( e.g., 500 ),    !
                  ! causing CAM5 crash. Thus, I used the discrete value of 1 by        !
                  ! commenting out the below 'if' block. This was the biggest bug.     !
                  ! Jan.30.2013. With a new constraint on eps0(m) on this day,         !
                  !              we don't need to use 'fmix_max' constraint in the     ! 
                  !              below case of 'exp_cmf .eq. 1'. So, I commented out   !
                  !              corresponding line.                                   !           
                  ! ------------------------------------------------------------------ !

                  if( dpa(m) * abs( eps(m) - del(m) ) .ge. 1.e-3_r8 ) then
                     fmix(m) = ( exp( dpa(m) * ( eps(m) - del(m) ) ) - 1._r8 ) / ( dpa(m) * ( eps(m) - del(m) ) )
                  else
                     fmix(m) = 1._r8
                  endif

                  ! ------------------------------------------------------------------------ !
                  ! Adiabatic Vertical Prognostic Equation                                   !
                  !  1. This should not include diabatic forcing.                            !
                  !  2. This should use the organized 'thle_m' not the 'thl_m'.              !
                  !  3. We are performing vertical integration upto the cloud top, not the   !
                  !     top interface. Thus, this is performed for all updraft plumes        !
                  !     including detached as well as non-detached updrafts.                 ! 
                  !  4. Also compute cloud condensate.                                       !
                  ! Apr.17.2012. In order to remove ambiguity due to variable use, I define  !
                  !             'thl_meu, qt_meu, u_meu, v_meu, ql_meu, qi_meu, tr_met(mt)'  !
                  !              similar to 'thl_med,...' for use in the below blocks.       !  
                  ! Nov.28.2012. In order to clarify adiabatic and diabatic processes,       !
                  !              any droplet activation due to lateral mixing SHOULD NOT be  !
                  !              treated in the 'progup_thlqt' below even though this        !
                  !              subroutine has a functionality to handle this by using      !
                  !              additional 'dia' input argument as the second argument.     !
                  !              In fact, we can handle diabatic process of 'nl,ni' in below !
                  !              but for that case, we should apply twice 'progup_thlqt' for ! 
                  !              each of the 'tr_au(m,ixnumliq)' and 'tr_au(m,ixnumice)',    !
                  !              first with dia = 0 and second with dia > 0 or < 0. Then by  !
                  !              subtracting the two, we should compute diabatic forcing     !
                  !              on each of 'nl,ni'.  This can be done in future.            !    
                  ! ------------------------------------------------------------------------ !

                  thl_meu = thle_b(k) + ssthle(k) * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  qt_meu  = qte_b(k)  + ssqte(k)  * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  u_meu   = ue_b(k)   + ssue(k)   * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  v_meu   = ve_b(k)   + ssve(k)   * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  ql_meu  = qle_b(k)  + ssqle(k)  * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  qi_meu  = qie_b(k)  + ssqie(k)  * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  do mt = 1, ncnst
                     tr_meu(mt)  = tre_b(k,mt) + sstre(k,mt) * ( 0.5_r8 * ( ptop(m) + p_b ) - p_b )
                  enddo

                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, ptop(m), thl_meu, ssthle(k), thl_au(m), thl_aut(m) )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, ptop(m), qt_meu, ssqte(k), qt_au(m), qt_aut(m) )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, ptop(m), ql_meu, ssqle(k), ql_au(m), ql_aut_adi )
                  call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, ptop(m), qi_meu, ssqie(k), qi_au(m), qi_aut_adi )
                  do mt = 1, ncnst
                     call progup_thlqt( eps(m), 0._r8, 0._r8, p_b, ptop(m), tr_meu(mt), sstre(k,mt), tr_au(m,mt), tr_aut(m,mt) )
                  enddo

                  ! Nov.28.2012. Impose consistency between droplet mass and droplet number.
                  !              Note that physically, droplet activation should not be performed here but
                  !              performed later in association with 'CEF' process.
                  !              In principle, I only need below two constraint lines.
                  !              However, in order to impose a constraint that in-cloud droplet radius is fixed
                  !              by the externally specified value for this version of code, I need two additional
                  !              lines further below.   

                  ! Below should be used for future generalized cloud microphysics.
                  ! if( ql_aut_adi .eq. 0._r8 ) tr_aut(m,ixnumliq) = 0._r8
                  ! if( qi_aut_adi .eq. 0._r8 ) tr_aut(m,ixnumice) = 0._r8           

                  ! Below is used instead of the above two lines to satisfy the constraint of constant in-cumulus
                  ! droplet radius. In future's generalized microphysics, above two lines should be used instead of below two lines. 
                  tr_aut(m,ixnumliq) = ql_aut_adi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq )
                  tr_aut(m,ixnumice) = qi_aut_adi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice )

                  ! CORRECTION
                  ! Note that the effect of PGFc is separately trested later. Thus, in order to
                  ! prevent double counting, we should set ssue = ssve = 0 below. 
                  ! Considering the case of non-constructed profile, I should definitely use the 
                  ! separate one at the top later below.   

                  u_grdPGF = ssu0(k)
                  v_grdPGF = ssv0(k)            

                  call progup_uv( eps(m), PGFc_up, p_b, ptop(m), u_meu, ssue(k), u_grdPGF, u_au(m), u_aut(m) )
                  call progup_uv( eps(m), PGFc_up, p_b, ptop(m), v_meu, ssve(k), v_grdPGF, v_au(m), v_aut(m) )
                  call progup_uv( eps(m), 0._r8,   p_b, ptop(m), u_meu, ssue(k), u_grdPGF, u_au(m), u_aut_adi )
                  call progup_uv( eps(m), 0._r8,   p_b, ptop(m), v_meu, ssve(k), v_grdPGF, v_au(m), v_aut_adi )
                  
                  ! -------------------------------------------------------------------- !
                  !                                                                      !
                  ! Treatment of Diabatic Forcings at the cloud top or the top interface !
                  !                                                                      !
                  ! -------------------------------------------------------------------- !

                  ! ------------------------------------------------------------------------- !
                  ! Apr.08.2013.                                                              !
                  ! Below block computing                                                     !
                  !    (1) 'exql, exqi, extr(mt)',                                            !
                  !    (2) 'evpR, evpS, evpRStr(mt)',                                         !
                  !    (3) 'eff_ql_au(m), eff_qi_au(m), eff_tr_au(m,mt)'                      ! 
                  ! are added on this day in association with fully analytical treatment      !
                  ! of precipitation fall-out within convective updraft.                      !
                  ! ------------------------------------------------------------------------- ! 

                  ! ------------------------------------------------------------------------- !
                  ! Compute 'exql,exqi' by solving analytical vertical integration of 'ql,qi' ! 
                  ! by including differential precipitation fall-out in the integration.      !
                  ! ------------------------------------------------------------------------- !

                  if( iprd_prep .eq. 5 ) then  ! Backward Analytical Method

                     call conden( ptop(m), thl_aut(m), qt_aut(m), th, qv, ql, qi, qse, id_check )
                     if( ( ql_cu + qi_cu ) .gt. criqc ) then
                        eps_dia_L = c0_ac * ( ( ql_cu + qi_cu ) - criqc ) / ( ql_cu + qi_cu )
                        eps_dia_I = c0_ac * ( ( ql_cu + qi_cu ) - criqc ) / ( ql_cu + qi_cu )
                     else
                        eps_dia_L = 0._r8
                        eps_dia_I = 0._r8
                     endif
                     call progup_thlqt( eps(m), eps_dia_L, 0._r8, p_b, ptop(m), ql_meu, ssqle(k), ql_au(m), ql_aut_adi_prp )
                     call progup_thlqt( eps(m), eps_dia_I, 0._r8, p_b, ptop(m), qi_meu, ssqie(k), qi_au(m), qi_aut_adi_prp )
                     exql = min( max( ql_aut_adi - ql_aut_adi_prp, 0._r8 ), 0.99_r8 * ql )  ! This should be guaranteed to be positive at this stage.
                     exqi = min( max( qi_aut_adi - qi_aut_adi_prp, 0._r8 ), 0.99_r8 * qi )  ! This should be guaranteed to be positive at this stage.
                     if( mclimit .eq. 1 ) then
                        tmp1 = exql + exqi
                        tmp2 = min( tmp1, max( ql + qi - criqc, 0._r8 ) )                 ! To impose a continuous variation across ql + qi = criqc.
                        exql = exql * ( tmp2 / max( tmp1, nonzero ) )
                        exqi = exqi * ( tmp2 / max( tmp1, nonzero ) )
                     endif
                     do mt = 1, ncnst
                        if( mt .eq. 1 ) then
                           extr(mt) = 0._r8
                        elseif( mt .eq. ixcldliq ) then
                           extr(mt) = exql
                        elseif( mt .eq. ixcldice ) then
                           extr(mt) = exqi
                        elseif( mt .eq. ixnumliq ) then
                           extr(mt) = exql * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq )
                        elseif( mt .eq. ixnumice ) then
                           extr(mt) = exqi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice )
                        else
                           ! ----------------------------------------------------------------------------------------- !
                           ! Wet deposition of aerosols (both interstitial and cloud-borne) within convective updarft. !
                           ! Below is a very simple treatment which should be refined in future.                       !
                           ! ----------------------------------------------------------------------------------------- !
                           extr(mt) = tr_aut(m,mt) * ( ( exql + exqi ) / max( qt_aut(m), nonzero ) )
                           ! Nov.29.2013. Following the reviewer's comments, use 'ql+qi' instead of 'qt_aut(m)' 
                           !              in computing 'extr(mt)' above.
                           extr(mt) = caer * tr_aut(m,mt) * min( 1._r8, ( ( exql + exqi ) / max( ql + qi, nonzero ) ) )
 
                        endif
                     enddo
                     ! -------------------------------------------------- !
                     ! Evaporation within Updraft.                        !
                     ! Set it to be zero, but can be refined in future.   !
                     ! -------------------------------------------------- !
                     evpR = 0._r8
                     evpS = 0._r8
                     do mt = 1, ncnst
                        evpRStr(mt) = 0._r8
                     enddo
                     ! ---------------------------------------------------------------------------------------------------------- !
                     ! Compute effective tendency of 'ql,qi' by 'CEF' process.                                                    !
                     ! In the below 'thl_tmp, qt_tmp' are the final correct updraft state variable at the top interface           !
                     ! that includes both 'mixing' and 'precipitation fall-out'.                                                  !
                     ! Note that I am using 'exn0(k)' instead of 'exnf(ptop(m))' in the below line to be consistent with          !
                     ! the other part of the code to satisfy energy-moisture conservation principle.                              ! 
                     ! IMPORTANT :                                                                                                !
                     !             In case of convective updraft when analytical integration is used as in the current code,      !
                     !             however, we should compute these 'effective CEF forcings' using the variables                  !
                     !             including 'precipitation fallout' both adiabatically and diabatically, since analytical        !
                     !             computation of precipitation fall-out will break the imposed fractional relationship           !
                     !             of 'ql / (ql + qi ) = f(T)', so that the restoration of this relationship using 'conden' and   !
                     !             associated heating-cooling process is eventually included as a part of 'CEF' forcing.          !
                     ! ---------------------------------------------------------------------------------------------------------- !
                     thl_tmp = thl_aut(m) + ( xlv / cp / exn0(k) ) * ( exql - evpR )  + ( xls / cp / exn0(k) ) * ( exqi - evpS )
                     qt_tmp  =  qt_aut(m) - exql - exqi + evpR + evpS
                     call conden( ptop(m), thl_tmp, qt_tmp, th_tmp, qv_tmp, ql_tmp, qi_tmp, qs_tmp, id_check )
                     eff_ql_au(m) = ql_tmp - ql_aut_adi_prp
                     eff_qi_au(m) = qi_tmp - qi_aut_adi_prp
                     do mt = 1, ncnst
                        if( mt .eq. ixnumliq ) then
                           ! Nov.28.2012. Below block is the new code. Note that in order to decide whether droplet activation
                           !              process should be performed or not, we should compare 'ql' with 'ql_u(m)' not
                           !              with 'ql_aut_adi'.
                           !              Below droplet activation form is not a general formula but a specific one designed 
                           !              to satisfy the constraint of a constant in-cloud droplet size specified externally.
                           !              In future's generalized cloud microphysics, more generalized droplet activation 
                           !              form should be used. The same is true for ice.   
                           ! if( ql_au(m) .eq. 0._r8 .and. ql .gt. 0._r8 ) then ! Droplet Activation   
                           !     eff_tr_au(m,mt) = eff_ql_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                           ! else
                           !   ! The second line assumes that only evaporation changes droplet number.
                           !   ! The choice of second line should be made in consistent with the one in the downdraft process.  
                           !     eff_tr_au(m,mt) =      eff_ql_au(m)          * ( tr_aut(m,mt) / max( ql_aut_adi, nonzero ) )
                           !   ! eff_tr_au(m,mt) = min( eff_ql_au(m), 0._r8 ) * ( tr_aut(m,mt) / max( ql_aut_adi, nonzero ) )
                           ! endif 
                           ! Below line is the old code.
                           eff_tr_au(m,mt) = eff_ql_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                        elseif( mt .eq. ixnumice ) then 
                           ! Nov.28.2012. Below block is the new code.
                           ! if( qi_au(m) .eq. 0._r8 .and. qi .gt. 0._r8 ) then ! Crystal Nucleation   
                           !     eff_tr_au(m,mt) = eff_qi_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 
                           ! else
                           !   ! The second line assumes that only evaporation changes droplet number.
                           !   ! The choice of second line should be made in consistent with the one in the downdraft process.  
                           !     eff_tr_au(m,mt) =      eff_qi_au(m)          * ( tr_aut(m,mt) / max( qi_aut_adi, nonzero ) )
                           !   ! eff_tr_au(m,mt) = min( eff_qi_au(m), 0._r8 ) * ( tr_aut(m,mt) / max( qi_aut_adi, nonzero ) )
                           ! endif 
                           ! Below line is the old code.
                           eff_tr_au(m,mt) = eff_qi_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 
                        else
                           eff_tr_au(m,mt) = 0._r8
                        endif
                        eff_tr_au(m,mt) = max( eff_tr_au(m,mt), qmin(mt) - tr_aut(m,mt) )   
                     enddo

                     ! ------------------------------------------------------------------------- !
                     ! Apr.08.2013.                                                              !
                     ! Above block computing                                                     !
                     !    (1) 'exql, exqi, extr(mt)',                                            !
                     !    (2) 'evpR, evpS, evpRStr(mt)',                                         !
                     !    (3) 'eff_ql_au(m), eff_qi_au(m), eff_tr_au(m,mt)'                      ! 
                     ! are added on this day in association with fully analytical treatment      !
                     ! of precipitation fall-out within convective updraft.                      !
                     ! Note that at this stage, none of (1),(2),(3) are used to update updraft   !
                     ! state variables at the top interface, which will be done further below.   !
                     ! ------------------------------------------------------------------------- ! 

                     ! ------------------------------------------------------------------------------------------------- !
                     ! TRACERS REFINEMENT NECESSARY : ADIABATIC CONDENSATION-EVAPORATION-FREEZING DURING VERTICAL MOTION !
                     ! Note that this CEF does not change 'thl_aut(m),qt_aut(m)'. So, here I only update 'tr_aut(m,mt)'. !
                     ! ------------------------------------------------------------------------------------------------- !
                     do mt = 1, ncnst
                        tr_aut(m,mt) = tr_aut(m,mt) + eff_tr_au(m,mt)
                     enddo
                     
                  else ! The others of Backward Analytical Method

                     ! ------------------------------------------------------------------------------- !
                     ! 1. Compute diabatic 'condensation-evaporation-freezing' on the cloud condensate !
                     !    In case of tracers, I temporarily set it to be zero. However, it should be   !
                     !    correctly computed later, depending on whether the tracers are cloud droplet !
                     !    number concentration or other aerosol tracers ( mass and number ).  If they  !
                     !    are droplet numbers, 'eff_tr_au' will be nonzero due to the evaporation and  !
                     !    freezing, but if they are other tracers, it is likely that 'eff_tr_au=0'.    !
                     !    Note that this 'CEF' process does not change 'thl,qt,u,v' which are          !
                     !    conserved scalars. Similarly, tracers other than 'nl,ni' are likely to be    !
                     !    invariant to this CEF process.                                               !
                     !    In the below, 'ql,qi' contains 'mixing' and 'CEF' already.                   !
                     ! Nov.28.2012. In order to correctly handle the changes of droplet numbers        !
                     !    associated with this 'CEF', we should be able to separately handle freezing  !
                     !    process. Computation of this separate freezing process is possible if we     !
                     !    think this CEF process as a seris of 'CE' --> 'F' process.                   !
                     !    In this case, eff_fz_au(m) = eff_qi_au(m). However, when there is only ice   !
                     !    not the liquid, this approach can be problematic, since we will              !
                     !    continuously generate ice crystal number by treating 'qi' increase as        !
                     !    freezing process.                                                            !    
                     !    This is a process reducing pressure without changing thl_aut(m) & qt_aut(m), !
                     !    that is, 'eff_ql_au(m) + eff_qi_au(m) > 0', so that condensate is always     !
                     !    generated. This property provides a clue for handling this process.          !
                     ! Apr.08.2013. Below block is commented-out with '!prp' since it is now being     !
                     !    computed above in association with analytical treatment of precipitation     !
                     !    fall-out within convective updraft above.                                    ! 
                     ! ------------------------------------------------------------------------------- !
                     call conden( ptop(m), thl_aut(m), qt_aut(m), th, qv, ql, qi, qse, id_check )
                     eff_ql_au(m) = ql - ql_aut_adi
                     eff_qi_au(m) = qi - qi_aut_adi
                     ! ------------------------------------------------------------------------------------------------- !
                     ! TRACERS REFINEMENT NECESSARY : ADIABATIC CONDENSATION-EVAPORATION-FREEZING DURING VERTICAL MOTION !
                     ! ------------------------------------------------------------------------------------------------- !
                     ! Nov.08.2011. Critical bug fix. I should correctly handle the droplet number concentration in order!
                     !              to prevent generating unreasonable source of droplet number which was the case in the!
                     !              previous wrong code.                                                                 ! 
                     ! Nov.28.2012.                                                                                      !
                     ! Below treatment implies that                                                                      !
                     !   (1) Activation occurs with the specified 'droprad_liq, droprad_ice'                             !
                     !   (2) Condensation occurs on the existing liquid and ice droplets,                                !
                     !   (3) Evaporation reduces the numbers of liquid and ice droplets.                                 !
                     ! where treatment of (2) allows the growth of cloud liquid and ice crystals, so that                !
                     ! future parameterization of precipitation production as a function of droplet size                 !
                     ! will be possible.                                                                                 !
                     ! The only caveat of this approach is that when freezing occurs so that some of                     !
                     ! the existing liquid droplets are converted into the ice crystals, corresponding                   !
                     ! increase of tr_aut(m,ixnumice) cannot be treated, even though decrease of tr_aut(m,ixnumliq)      !
                     ! is treated. However, some of these processes are treated by ice nucleation process, so that       !
                     ! below treatment is not so bad. However, more complete treatment should be made in future.         !
                     ! For the time being, for consistency with previous code, let's also change droplet number when     !
                     ! condensation occurs.                                                                              ! 
                     ! Apr.08.2013. Below block is commented-out with '!prp' since it is now being                       !
                     !              computed above in association with analytical treatment of precipitation             !
                     !              fall-out within convective updraft above.                                            ! 
                     ! ------------------------------------------------------------------------------------------------- !
                     do mt = 1, ncnst
                        if( mt .eq. ixnumliq ) then
                           ! Nov.28.2012. Below block is the new code. Note that in order to decide whether droplet activation
                           !              process should be performed or not, we should compare 'ql' with 'ql_u(m)' not
                           !              with 'ql_aut_adi'.
                           !              Below droplet activation form is not a general formula but a specific one designed 
                           !              to satisfy the constraint of a constant in-cloud droplet size specified externally.
                           !              In future's generalized cloud microphysics, more generalized droplet activation 
                           !              form should be used. The same is true for ice.   
                           ! if( ql_au(m) .eq. 0._r8 .and. ql .gt. 0._r8 ) then ! Droplet Activation   
                           !     eff_tr_au(m,mt) = eff_ql_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                           ! else
                           !   ! The second line assumes that only evaporation changes droplet number.
                           !   ! The choice of second line should be made in consistent with the one in the downdraft process.  
                           !     eff_tr_au(m,mt) =      eff_ql_au(m)          * ( tr_aut(m,mt) / max( ql_aut_adi, nonzero ) )
                           !   ! eff_tr_au(m,mt) = min( eff_ql_au(m), 0._r8 ) * ( tr_aut(m,mt) / max( ql_aut_adi, nonzero ) )
                           ! endif 
                           ! Below line is the old code.
                           eff_tr_au(m,mt) = eff_ql_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                        elseif( mt .eq. ixnumice ) then 
                           ! Nov.28.2012. Below block is the new code.
                           ! if( qi_au(m) .eq. 0._r8 .and. qi .gt. 0._r8 ) then ! Crystal Nucleation   
                           !     eff_tr_au(m,mt) = eff_qi_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 
                           ! else
                           !   ! The second line assumes that only evaporation changes droplet number.
                           !   ! The choice of second line should be made in consistent with the one in the downdraft process.  
                           !     eff_tr_au(m,mt) =      eff_qi_au(m)          * ( tr_aut(m,mt) / max( qi_aut_adi, nonzero ) )
                           !   ! eff_tr_au(m,mt) = min( eff_qi_au(m), 0._r8 ) * ( tr_aut(m,mt) / max( qi_aut_adi, nonzero ) )
                           ! endif 
                           ! Below line is the old code.
                           eff_tr_au(m,mt) = eff_qi_au(m) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 
                        else
                           eff_tr_au(m,mt) = 0._r8
                        endif
                        ! Nov.27.2012. Impose a constraint on the diabatic forcing to prevent the onset of negative tracer during vertical motion.
                        !              This constraint is imposed by using 'tr_aut(m,mt) + eff_tr_au(m,mt) > qmin(mt)' criteria,
                        !              where 'tr_aut(m,mt)' is a tracer at the base interface before adding diabatic forcing.
                        !              Similar constraint has also been imposed on the downdraft motion.
                        eff_tr_au(m,mt) = max( eff_tr_au(m,mt), qmin(mt) - tr_aut(m,mt) )   
                     enddo

                     ! ------------------------------------------------------------------------------------------------- !
                     ! TRACERS REFINEMENT NECESSARY : ADIABATIC CONDENSATION-EVAPORATION-FREEZING DURING VERTICAL MOTION !
                     ! Note that this CEF does not change 'thl_aut(m),qt_aut(m)'. So, here I only update 'tr_aut(m,mt)'. !
                     ! ------------------------------------------------------------------------------------------------- !
                     do mt = 1, ncnst
                        tr_aut(m,mt) = tr_aut(m,mt) + eff_tr_au(m,mt)
                     enddo

                     ! ------------------------------------------------------------------------- !
                     ! 2. Precipitation fallout at the cloud top                                 !
                     !    a. Future refinement of microphysics only requires recomputation of    !
                     !       production of precipitation, exql, exqi.                            !
                     !    b. I used 'exn_t' for consistency with later computation of 'slten_u'. !
                     !       Apr.21.2011. Use 'exn_top' instead of 'exn_t'.                      ! 
                     !    c. In case of tracers other than droplet number concentration,         !
                     !       I temporary set 'prep_tr_au' to be zero but correct value should be ! 
                     !       calculated by mimicking the in-cloud wet scavenging routine later.  !
                     ! ------------------------------------------------------------------------- !

                     ! -------------------------------------------------------------------------------------------------------------------- !
                     ! Feb.07.2013.                                                                                                         ! 
                     ! Compute precipitation production rate at the cumulus top during upward motion from 'p_b(z_b)' to 'ptop(ztop)'.       !
                     ! It is possible that w_aut(m) = 0 if 'ptop' is the final cumulus top instead of the top interface.                    !
                     ! Below subroutine 'prod_prep_up' computes the amount of precipitated condensate ( exql, exqi >= 0 in [kg/kg],         !
                     ! extr(1:ncnst) >= in [#(kg)/kg] ) during 'delta_t = ( ztop(m) - z_b ) / ( 0.5 * ( w_au(m) + w_aut(m) ) )'             !
                     ! at the top interface.                                                                                                !
                     ! The final layer-mean precipitation production rate will be computed later using 'exql, exqi, extr(1:ncnst)'.         !
                     ! Note that extr(1) = 0., extr(ixcldliq) = exql, extr(ixcldice) = exqi,                                                !
                     ! extr(ixnumliq) = exql * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ),                              !
                     ! extr(ixnumice) = exqi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ).                              ! 
                     ! Most importantly, 'extr(other indices)' is the wet deposition of aerosols (both interstiail and cloud-borne)         !
                     ! within cumulus updraft.                                                                                              !
                     ! Note since am_pu_msfc(k,msfc) is computed using am_us_msfc(k,msfc) not by using                                      !
                     ! am_u_msfc(k,msfc), I am using am_us_msfc(k,msfc) in the below subroutine as an input                                 !
                     ! argument. However, in association with the use of beta2 = 1 and evaporation within                                   !
                     ! PBL environment, I should come up with more satisfactory cloud overlapping structure                                 !
                     ! in future. This is always related with the treatment of wet deposition of aerosols                                   !
                     ! within unsaturated convective updraft. I should figure this out before AMWG.                                         !
                     ! Feb.09.2013.                                                                                                         !
                     ! Note that the input argument is 'am_u_msfc(k,msfc)' not 'am_us_msfc(k,msfc)' since evaporation within updraft is also! 
                     ! treated when updraft is unsaturated.                                                                                 !
                     ! Apr.08.2013. In association with fully analytical treatment precipitation fall-out within convective updraft above,  !
                     !              I commented-out below block with '!prp'. But it can be restored at any time in future.                  ! 
                     ! -------------------------------------------------------------------------------------------------------------------- ! 

                     msfc = msfc_from_m(k,m)
                     call prod_prep_up( z_b, ztop(m), p_b, ptop(m), (ptop(m)/p00)**rovcp, exn0(k),  &
                                     w_au(m), w_aut(m),                                          &
                                     thl_aut(m), qt_aut(m), ql, qi, tr_aut(m,1:ncnst),           &
                                     S_b_ql_au(m), S_b_qi_au(m), iprd_prep,                      & 
                                     ql_cu, qi_cu, eps(m),                                       &
                                     thl_meu, ssthle(k), thl_au(m), qt_meu, ssqte(k), qt_au(m),  & 
                                     ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice, i, k, lchnk, &
                                     flxrain_msfc(k,msfc), flxsnow_msfc(k,msfc),                 &     
                                     a_p_msfc(k,msfc), am_u_msfc(k,msfc), am_pu_msfc(k,msfc),    &
                                     caer, criqc, c0_ac,                                         &
                                     exql, exqi, extr(1:ncnst), S_t_ql_au(m), S_t_qi_au(m),      &
                                     evpR, evpS, evpRStr(1:ncnst) )

                  endif ! End of Backward Analytical Method


                  prep_qtl_au(m)  = - exql
                  prep_qti_au(m)  = - exqi

                  prep_thll_au(m) = - ( xlv / cp / exn0(k) ) * prep_qtl_au(m)
                  prep_thli_au(m) = - ( xls / cp / exn0(k) ) * prep_qti_au(m)    

                  ! ------------------------------------------------------------------------- !
                  ! TRACERS REFINEMENT NECESSARY : PRECIPITATION FALLOUT AT THE TOP INTERFACE !
                  ! ------------------------------------------------------------------------- !

                  do mt = 1, ncnst
                     prep_tr_au(m,mt) = - extr(mt)
                     prep_tr_au(m,mt) = max( prep_tr_au(m,mt), qmin(mt) - tr_aut(m,mt) )   
                  enddo
                  ! ------------------------------------------------------------------------- !
                  ! TRACERS REFINEMENT NECESSARY : PRECIPITATION FALLOUT AT THE TOP INTERFACE !
                  ! ------------------------------------------------------------------------- !
                  ! ------------------------------------------------------------------------------------------- !
                  ! Evaporation of Precipitation at the cloud top                                               !
                  ! ------------------------------------------------------------------------------------------- !    
                  evp_qtl_au(m)  =  evpR
                  evp_qti_au(m)  =  evpS
                  evp_thll_au(m) =  - ( xlv / cp / exn0(k) ) * evp_qtl_au(m)
                  evp_thli_au(m) =  - ( xls / cp / exn0(k) ) * evp_qti_au(m)
                  ! -------------------------------------------------------------------------------------- !
                  ! TRACERS REFINEMENT NECESSARY : EVAPORATION OF CONVECTIVE PRECIPITATION WITHIN UPDRAFTA !
                  ! -------------------------------------------------------------------------------------- !
                  do mt = 1, ncnst
                     evp_tr_au(m,mt)  =  evpRStr(mt)
                  enddo
                  ! ------------------------------------------------------------------------------------------- !
                  ! 0. Evaporation of Precipitation at the cloud top                                            !
                  !    a. Since updraft is unsaturated within the PBL, this evaporation should be treated here. !
                  !       This may be especially important for treating deep convection.                        !  
                  !    b. I should treat only evaporation of 'convective precipitation' not the 'stratiform     !
                  !       precipitation'. This is important to satisfy energy and moisture conservations.       !
                  !    c. Any diabatic forcings within convective updraft and downdraft will induce the         !
                  !       tendencies of environmental variables.                                                !
                  !    d. Simply neglect evaporation of precipitation within convctive updraft.                 !
                  ! ------------------------------------------------------------------------------------------- !
                  thl_aut(m)      = thl_aut(m) + prep_thll_au(m) + prep_thli_au(m) + evp_thll_au(m) + evp_thli_au(m) 
                  qt_aut(m)       = qt_aut(m)  +  prep_qtl_au(m) +  prep_qti_au(m) +  evp_qtl_au(m) +  evp_qti_au(m)
                  do mt = 1, ncnst
                     tr_aut(m,mt) = tr_aut(m,mt) + prep_tr_au(m,mt) + evp_tr_au(m,mt)
                  enddo
                  ! ----------------------------------------------------- !
                  ! 3. Horizontal PGF at the cloud top                    !
                  !    a. Instead of using the gradient within the layer, !
                  !       I should use the gradient between the layers.   !
                  ! ----------------------------------------------------- !
                  PGF_u_au(m) = u_aut(m) - u_aut_adi
                  PGF_v_au(m) = v_aut(m) - v_aut_adi

                  ! -------------------------------------------------------------------------- !
                  ! Computation of 'ql_aut(m)' and 'qi_aut(m)'                                 !
                  ! Note that this is done using a fully updated final 'thl_aut, qt_aut' which !
                  ! includes all the effect of 'mixing, precipitation fallout, evaporation.    !
                  ! -------------------------------------------------------------------------- !

                  call conden( ptop(m), thl_aut(m), qt_aut(m), th, qv, ql, qi, qse, id_check )
                  ql_aut(m)  = ql
                  qi_aut(m)  = qi
                  thv_aut(m) = th * ( 1._r8 + zvir * qv - ql - qi )
                  if( ql_aut(m) .eq. 0._r8 ) tr_aut(m,ixnumliq) = 0._r8
                  if( qi_aut(m) .eq. 0._r8 ) tr_aut(m,ixnumice) = 0._r8
 
                  ! ------------------------------------------------------------------------ !
                  ! Save cloud top height in each layer for each segment for future use.     !
                  ! Here, 'k' is defined as a layer mid-point index, not the base interface. !                
                  ! ------------------------------------------------------------------------ !

                  ptops(k,m) = ptop(m)
                  ztops(k,m) = ztop(m)

                  ! -------------------------------------------------------------------------------- !
                  ! Diagnostic Output : As of Jul.26.2011, these can be also used in the actual      !
                  ! as mentioned later below.                                                        !
                  ! Individual updraft segment ( 'msfc' ) properties                                 ! 
                  ! Mean downdraft properties for individual updraft segment will be computed later. !
                  ! Jul.26.2011. Below part of convective updraft properties at each interface are   !
                  !              moved into later part further below with additional inclusion of    !
                  !              non-zero cumulus updraft properties at the cumulus top.             !
                  ! -------------------------------------------------------------------------------- !
   
                  msfc = msfc_from_m(k,m)

                  eps0_u_msfc(km,msfc)     = eps0(m)
                  eps_u_msfc(km,msfc)      = eps(m)
                  del_u_msfc(km,msfc)      = del(m)                   
                  eeps_u_msfc(km,msfc)     = eeps(m)                   
                  ddel_u_msfc(km,msfc)     = ddel(m)
                  xc_u_msfc(km,msfc)       = xc(m)
                  xs_u_msfc(km,msfc)       = xs(m)
                  xemin_u_msfc(km,msfc)    = xe_min(m)
                  xemax_u_msfc(km,msfc)    = xe_max(m)                   
                  cridis_u_msfc(km,msfc)   = cridis                   
                  thvcuenv_u_msfc(km,msfc) = thv_cu - thv_env
                  thvegenv_u_msfc(km,msfc) = thv_eg - thv_env                   
                  thvxsenv_u_msfc(km,msfc) = thv_xs - thv_env
                  fmix_u_msfc(km,msfc)     = fmix(m)
                  cmfumix_u_msfc(km,msfc)  = fmix(m) * dpa(m) * eps0(m) * cmf_au(m)

                  ! ----------------- !
                  ! Diagnostic Output ! 
                  ! ----------------- !

               enddo  ! End of updraft segment 'm' loop

               ! ------------------------------------------------------------------------ !
               ! Number of detached and non-detached updraft segments                     ! 
               ! Also compute original ( at surface ) updraft segment index in each layer !
               ! Also compute the 'top layer index','top height in pressure and height'   !
               ! of the original updraft segment at surface.                              !  
               ! ------------------------------------------------------------------------ !

               nseg_nondet  = 0                      ! # of non-detached updraft segments
               do m = 1, N_up(km)
                  if( ytop(m) .gt. 0.5_r8 ) then
                     nseg_nondet = nseg_nondet + 1         
                     msfc_from_m(k+1,nseg_nondet) = msfc_from_m(k,m)
                  else 
                     ktop_msfc(msfc_from_m(k,m)) = k
                     ptop_msfc(msfc_from_m(k,m)) = ptops(k,m)
                     ztop_msfc(msfc_from_m(k,m)) = ztops(k,m)
                  endif
                  m_from_msfc(k,msfc_from_m(k,m)) = m
               enddo
               nseg_det = N_up(km) - nseg_nondet     ! # of detached updraft segments

               ! -------------------------------------------------------------------------- !
               ! Apply the constraint of updraft vertical velocity for each updraft segment !
               ! and compute the mass of the 3rd type of detrained airs.                    !
               ! Note that 'ybot(m)=1' always.                                              !
               ! -------------------------------------------------------------------------- !

               do m = 1, N_up(km)
                  if( ytop(m) .gt. 0.5_r8 .and. w_aut(m) .gt. wumax ) then
                     f_wu(m)        = ( 1._r8 - wumax / w_aut(m) )
                     f_srcds(k,m,3) = f_wu(m) * cmf_aut(m) / cmf_u(km)
                     f_srcrs(k,m,3) = 0._r8                     
                     f_srcrs2(k,m,3) = 0._r8                     
                     cmf_aut(m)     = cmf_aut(m) * ( 1._r8 - f_wu(m) )
                     w_aut(m)       =   w_aut(m) * ( 1._r8 - f_wu(m) )
                  endif
               enddo

               ! ------------------------------------------------------------------------ !
               ! Apply the constraint of updraft fractional area for each updraft segment !
               ! and compute the mass of the 3rd type of detrained airs.                  !
               ! Note that 'ybot(m)=1' always.                                            !
               ! ------------------------------------------------------------------------ !

               au_tent  = 0._r8 
               if( nseg_nondet .gt. 0.5_r8 ) then
                  do m = 1, N_up(km)
                     if( ytop(m) .gt. 0.5_r8 ) then
                        au_tent   = au_tent   +   a_aut(m)
                     endif
                  enddo
               endif

               f_nu = max( 0._r8, 1._r8 - au_max / max( nonzero, au_tent ) ) 

               if( au_tent .gt. au_max ) then
                  do m = 1, N_up(km)
                     if( ytop(m) .gt. 0.5_r8 ) then
                        f_srcds(k,m,3) = f_srcds(k,m,3) + f_nu * cmf_aut(m) / cmf_u(km)
                        cmf_aut(m) = cmf_aut(m) * ( 1._r8 - f_nu )
                        a_aut(m)   =   a_aut(m) * ( 1._r8 - f_nu )
                        rad_aut(m) = sqrt( a_aut(m) / num_aut(m) /  3.141592_r8 )   ! Physical radius of updraft plume [ m ]
                     endif
                  enddo
               endif

               ! ---------------------------------------------------------------------- !
               ! Mass-flux weighted mean or net updraft properties at the top interface !
               ! These are computed only using non-detached updrafts.                   !
               ! ---------------------------------------------------------------------- !

               if( nseg_nondet .gt. 0.5_r8 ) then

                  do m = 1, N_up(km)
                     if( ytop(m) .gt. 0.5_r8 ) then
                        cmf_u(k)  = cmf_u(k)  + cmf_aut(m)
                        num_u(k)  = num_u(k)  + num_aut(m)
                        a_u(k)    = a_u(k)    +   a_aut(m)
                        rad_u(k)  = rad_u(k)  + num_aut(m) * rad_aut(m)**2._r8
                        thl_u(k)  = thl_u(k)  + thl_aut(m) * cmf_aut(m)
                        qt_u(k)   = qt_u(k)   +  qt_aut(m) * cmf_aut(m)
                        u_u(k)    = u_u(k)    +   u_aut(m) * cmf_aut(m)
                        v_u(k)    = v_u(k)    +   v_aut(m) * cmf_aut(m)
                        w_u(k)    = w_u(k)    +   w_aut(m) * cmf_aut(m)
                        wa_u(k)   = wa_u(k)   +   w_aut(m) *   a_aut(m)
                        ql_u(k)   = ql_u(k)   +  ql_aut(m) * cmf_aut(m)  
                        qi_u(k)   = qi_u(k)   +  qi_aut(m) * cmf_aut(m)  
                        do mt = 1, ncnst
                           tr_u(k,mt) = tr_u(k,mt) + tr_aut(m,mt) * cmf_aut(m)  
                        enddo
                        qla_u(k)  = qla_u(k)  +  ql_aut(m) *   a_aut(m)  
                        qia_u(k)  = qia_u(k)  +  qi_aut(m) *   a_aut(m)  
                        thva_u(k) =thva_u(k)  + thv_aut(m) *   a_aut(m)  
                     endif
                  enddo
                  rad_u(k)  = sqrt( rad_u(k) / num_u(k) )                 ! Effective plume radius [ m ]
                  thl_u(k)  = thl_u(k) / cmf_u(k)
                  qt_u(k)   = qt_u(k)  / cmf_u(k)
                  u_u(k)    = u_u(k)   / cmf_u(k)
                  v_u(k)    = v_u(k)   / cmf_u(k)
                  w_u(k)    = w_u(k)   / cmf_u(k)
                  wa_u(k)   = wa_u(k)  / a_u(k)
                  ql_u(k)   = ql_u(k)  / cmf_u(k)                         ! Mass-flux weighted average of in-cloud liquid condensate   
                  qi_u(k)   = qi_u(k)  / cmf_u(k)                         ! Mass-flux weighted average of in-cloud liquid condensate
                  do mt = 1, ncnst
                     tr_u(k,mt) = tr_u(k,mt)  / cmf_u(k)             
                  enddo
                  qla_u(k)  = qla_u(k) / a_u(k)                           ! Area-weighting average of in-cloud liquid condensate   
                  qia_u(k)  = qia_u(k) / a_u(k)                           ! Area-weighting average of in-cloud liquid condensate
                  thva_u(k) =thva_u(k) / a_u(k)                           ! Area-weighting average of updraft buoyancy
               endif

               ! ----------------------------------------------------------------------------- !
               ! Mass-flux weighted diabatic change of conservative scalars at each cloud top. !
               ! This will be used for computing environmental tendency later.                 ! 
               ! We separately define cmf_u_dia(k) as well as cmf_u(k) to take into account of !
               ! non-zero updraft mass flux at just below the cloud top of each detached       !
               ! updraft for the purpose of treating diabatic forcing at the cloud top.        !
               ! ----------------------------------------------------------------------------- !

               do m = 1, N_up(km)
                  fac            = 0.5_r8 * ( cmf_au(m) + cmf_aut(m) ) * ( dpa(m) / dp_m )
                  if( ytop(m) .gt. 0.5_r8 ) fac = 0.5_r8 * ( cmf_au(m) + cmf_aut(m) / ( 1._r8 - f_nu ) / ( 1._r8 - f_wu(m) ) )
                  cmf_u_dia(k)   = cmf_u_dia(k)   +                    fac
                  evp_thll_u(k)  = evp_thll_u(k)  +   evp_thll_au(m) * fac
                  evp_qtl_u(k)   = evp_qtl_u(k)   +    evp_qtl_au(m) * fac
                  evp_thli_u(k)  = evp_thli_u(k)  +   evp_thli_au(m) * fac
                  evp_qti_u(k)   = evp_qti_u(k)   +    evp_qti_au(m) * fac
                  prep_thll_u(k) = prep_thll_u(k) +  prep_thll_au(m) * fac
                  prep_qtl_u(k)  = prep_qtl_u(k)  +   prep_qtl_au(m) * fac
                  prep_thli_u(k) = prep_thli_u(k) +  prep_thli_au(m) * fac
                  prep_qti_u(k)  = prep_qti_u(k)  +   prep_qti_au(m) * fac
                  PGF_u_u(k)     = PGF_u_u(k)     +      PGF_u_au(m) * fac
                  PGF_v_u(k)     = PGF_v_u(k)     +      PGF_v_au(m) * fac
                  eff_ql_u(k)    = eff_ql_u(k)    +     eff_ql_au(m) * fac
                  eff_qi_u(k)    = eff_qi_u(k)    +     eff_qi_au(m) * fac
                  do mt = 1, ncnst
                     evp_tr_u(k,mt)  =  evp_tr_u(k,mt) +   evp_tr_au(m,mt) * fac
                     prep_tr_u(k,mt) = prep_tr_u(k,mt) +  prep_tr_au(m,mt) * fac
                     eff_tr_u(k,mt)  =  eff_tr_u(k,mt) +   eff_tr_au(m,mt) * fac
                  enddo
                  cmf_u_mix(k)   = cmf_u_mix(k)   + fmix(m) * dpa(m) * eps0(m) * cmf_au(m)
                  ! ---------------------------------------------------------------------------- !
                  ! Compute individual segment's rain and snow production tendency by convective !
                  ! updrafts. Corresponding tendency by convective downdraft will be computed    !
                  ! later.                                                                       !  
                  ! ---------------------------------------------------------------------------- !
                  msfc = msfc_from_m(k,m)
                  qrten_u_msfc(k,msfc) = - ( g / dp0(k) ) * ( prep_qtl_au(m) + evp_qtl_au(m) ) * fac   ! >= 0
                  qsten_u_msfc(k,msfc) = - ( g / dp0(k) ) * ( prep_qti_au(m) + evp_qti_au(m) ) * fac   ! >= 0
                  ! --------------------------------------------------------------------------------------------------------- !
                  ! Nov.29.2012. I can compute corresponding 'trrsten_u_msfc(k,msfc,mt)' of tracers associated with           !
                  ! precipitation by adding ( prep_tr_au(m,mt) + evp_tr_au(m,mt) ) * fac'. This information will be used      !
                  ! to trace aerosol concentration within precipitation flux, which will be used to compute the increase of   !
                  ! aerosol concentration within convective downdraft and environment when precipitation is evaporated within !
                  ! convective downdraft and environment. Eventually, this will be used for computing aerosol wet scavenging  !
                  ! process associated with convective precipitation process instead of performing within wetdepa.            !
                  ! Note that I should use 'dptr0(k,mt)' instead of 'dp0(k)' for tracers.                                     !
                  ! Feb.06.2013. Note that wet deposition of aerosol within convetcive updraft (both cloud-borne and          !
                  !              interstitial) will be treated as a part of prep_tr_au(m,mt).                                 !
                  !              So, I don't need to use wdep_tr_au(m,mt) in the below trrsten_u_msfc().                      ! 
                  ! --------------------------------------------------------------------------------------------------------- !
                  do mt = 1, ncnst
                     if( mt .eq. ixcldliq ) then
                        trrsten_u_msfc(k,msfc,mt) = qrten_u_msfc(k,msfc)
                     elseif( mt .eq. ixcldice ) then
                        trrsten_u_msfc(k,msfc,mt) = qsten_u_msfc(k,msfc)
                     elseif( mt .eq. ixnumliq ) then
                        trrsten_u_msfc(k,msfc,mt) = qrten_u_msfc(k,msfc) * 3._r8 / &
                                                    ( 4._r8 * 3.141592_r8 * droprad_rain**3 * density_rain )                    
                     elseif( mt .eq. ixnumice ) then 
                        trrsten_u_msfc(k,msfc,mt) = qsten_u_msfc(k,msfc) * 3._r8 / &
                                                    ( 4._r8 * 3.141592_r8 * droprad_snow**3 * density_snow )                    
                     else
                        trrsten_u_msfc(k,msfc,mt) = - ( g / dptr0(k,mt) ) * ( prep_tr_au(m,mt) + evp_tr_au(m,mt) ) * fac   ! >= 0
                     endif
                  enddo
                  ! ---------------------------------------------------------------------------- !
                  !                                                                              !  
                  ! ---------------------------------------------------------------------------- !
               enddo
               if( cmf_u_dia(k) .gt. nonzero ) then
                  evp_thll_u(k)  = evp_thll_u(k)  / cmf_u_dia(k)
                  evp_qtl_u(k)   = evp_qtl_u(k)   / cmf_u_dia(k)
                  evp_thli_u(k)  = evp_thli_u(k)  / cmf_u_dia(k)
                  evp_qti_u(k)   = evp_qti_u(k)   / cmf_u_dia(k)
                  prep_thll_u(k) = prep_thll_u(k) / cmf_u_dia(k)
                  prep_qtl_u(k)  = prep_qtl_u(k)  / cmf_u_dia(k)
                  prep_thli_u(k) = prep_thli_u(k) / cmf_u_dia(k)
                  prep_qti_u(k)  = prep_qti_u(k)  / cmf_u_dia(k)
                  eff_ql_u(k)    = eff_ql_u(k)    / cmf_u_dia(k)
                  eff_qi_u(k)    = eff_qi_u(k)    / cmf_u_dia(k)
                  PGF_u_u(k)     = PGF_u_u(k)     / cmf_u_dia(k)
                  PGF_v_u(k)     = PGF_v_u(k)     / cmf_u_dia(k)
                  do mt = 1, ncnst
                     evp_tr_u(k,mt)  =  evp_tr_u(k,mt) / cmf_u_dia(k)
                     prep_tr_u(k,mt) = prep_tr_u(k,mt) / cmf_u_dia(k)
                     eff_tr_u(k,mt)  =  eff_tr_u(k,mt) / cmf_u_dia(k)
                  enddo
               else
                  cmf_u_dia(k)   = 0._r8
                  evp_thll_u(k)  = 0._r8
                  evp_qtl_u(k)   = 0._r8
                  evp_thli_u(k)  = 0._r8
                  evp_qti_u(k)   = 0._r8
                  prep_thll_u(k) = 0._r8
                  prep_qtl_u(k)  = 0._r8
                  prep_thli_u(k) = 0._r8
                  prep_qti_u(k)  = 0._r8
                  eff_ql_u(k)    = 0._r8
                  eff_qi_u(k)    = 0._r8
                  PGF_u_u(k)     = 0._r8
                  PGF_v_u(k)     = 0._r8
                  do mt = 1, ncnst
                     evp_tr_u(k,mt)  = 0._r8
                     prep_tr_u(k,mt) = 0._r8
                     eff_tr_u(k,mt)  = 0._r8
                  enddo
               endif

               ! ----------------------------------------------------------------------------- !
               !                                                                               !
               ! 3 sources of downdrafts and 1 source of detrained air from each segment level !
               !                                                                               !
               ! ----------------------------------------------------------------------------- !

               ! ------------------------------------------------------------------------------------- !                 
               ! 1. Mixing Downdraft + Detrained Updraft                                               !
               !    This is a lateral detrainment after updraft buoyancy sorting at the base interface !       
               ! ------------------------------------------------------------------------------------- !

               do m = 1, N_up(km)
                  ! ---------------------------------------------------------------------------- !
                  ! 1. Since updraft buoyancy sorting was performed at the base interface,       !
                  !    we should use 'thl_au(m)' not 'thl_aut(m)'.                               !
                  ! 2. For consistent treatment of buoyancy sorting with organized environmental ! 
                  !    airs, we should use 'thl_eg' not 'thl_b'.                                 !
                  ! 3. Since 'compute_PDF' sets 'zmass = 0' when 'xe_min(m) = xe_max(m)',        !  
                  !    below computation is correct in general case.                             !
                  ! ---------------------------------------------------------------------------- !
                  ! ---------------- !
                  ! Mixing Downdraft !
                  ! ---------------- !
                  ! ------------------------------------------------------------------ !
                  ! fmix(m) * dpa(m) * eps0(m) * cmf_au(m) :                           !
                  ! The amount of updraft mass involved in the buoyancy sorting mixing !
                  ! ------------------------------------------------------------------ !
                  
                  call compute_PDF( 'PDFbsQ', xe_min(m), xe_max(m), zbar, zmass, zmass_up )
                  f_srcds(k,m,1)   =    zmass * fmix(m) * dpa(m) * eps0(m) * 2._r8 * cmf_au(m) / cmf_u(km)
                  thl_srcds(k,m,1) = ( thl_au(m) + zbar * ( thl_eg - thl_au(m) ) )
                  qt_srcds(k,m,1)  = ( qt_au(m)  + zbar * ( qt_eg  -  qt_au(m) ) )
                  u_srcds(k,m,1)   = ( u_au(m)   + zbar * ( u_eg   -   u_au(m) ) )
                  v_srcds(k,m,1)   = ( v_au(m)   + zbar * ( v_eg   -   v_au(m) ) )
                  do mt = 1, ncnst
                     tr_srcds(k,m,1,mt) = ( tr_au(m,mt) + zbar * ( tr_eg(mt) - tr_au(m,mt) ) )
                  enddo
                  ! -------------------------------------------------------------------------------------------------- !
                  ! Nov.28.2012. Impose consistency between droplet mass and droplet number.                           !
                  !              Note that this should be computed here since it involves the use of 'zbar'.           ! 
                  !              Mixing downdraft is generated at the base interface which will be the source level of !
                  !              mixing downdraft layer. Thus, I am using 'ps0(km)' in the below 'conden' subroutine.  !
                  ! -------------------------------------------------------------------------------------------------- !
                  call conden( ps0(km), thl_srcds(k,m,1), qt_srcds(k,m,1), tmp_th, tmp_qv, tmp_ql, tmp_qi, tmp_qse, id_check )
                  ql_srcds(k,m,1)  = tmp_ql
                  qi_srcds(k,m,1)  = tmp_qi
                  tr_srcds(k,m,1,ixnumliq) = ql_srcds(k,m,1) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq )
                  tr_srcds(k,m,1,ixnumice) = qi_srcds(k,m,1) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice )

                  ! ----------------- !
                  ! Detrained Updraft !
                  ! ----------------- !
                  call compute_PDF( 'PDFbsQ', xc(m),    xe_min(m), zbar1, zmass1, zmass_up1 )
                  call compute_PDF( 'PDFbsQ', xe_max(m),    1._r8, zbar2, zmass2, zmass_up2 )
                  zmass_up = zmass_up1 + zmass_up2
                  zmass = zmass1 + zmass2
                  zbar  = 0._r8
                  if( zmass .gt. nonzero ) zbar = ( zbar1 * zmass1 + zbar2 * zmass2 ) / zmass
                  f_srcrs(k,m,1)   =    zmass * fmix(m) * dpa(m) * eps0(m) * 2._r8 * cmf_au(m) / cmf_u(km)
                  thl_srcrs(k,m,1) = ( thl_au(m) + zbar * ( thl_eg - thl_au(m) ) )
                  qt_srcrs(k,m,1)  = ( qt_au(m)  + zbar * ( qt_eg  -  qt_au(m) ) )
                  u_srcrs(k,m,1)   = ( u_au(m)   + zbar * ( u_eg   -   u_au(m) ) )
                  v_srcrs(k,m,1)   = ( v_au(m)   + zbar * ( v_eg   -   v_au(m) ) )
                  do mt = 1, ncnst
                     tr_srcrs(k,m,1,mt) = ( tr_au(m,mt) + zbar * ( tr_eg(mt) - tr_au(m,mt) ) )
                  enddo
                  ! -------------------------------------------------------------------------------------------------- !
                  ! Nov.28.2012. Impose consistency between droplet mass and droplet number.                           !
                  !              Note that this should be computed here since it involves the use of 'zbar'.           ! 
                  !              Detrained downdraft will be eventually detrained at the layer mid-point as has been   !
                  !              assumed in the previous code later.                                                   !
                  !              Thus, I am using 'p0(k)' in the below 'conden' subroutine.                            !
                  ! -------------------------------------------------------------------------------------------------- !
                  call conden( p0(k), thl_srcrs(k,m,1), qt_srcrs(k,m,1), tmp_th, tmp_qv, tmp_ql, tmp_qi, tmp_qse, id_check )
                  ql_srcrs(k,m,1)  = tmp_ql
                  qi_srcrs(k,m,1)  = tmp_qi
                  tr_srcrs(k,m,1,ixnumliq) = ql_srcrs(k,m,1) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq )
                  tr_srcrs(k,m,1,ixnumice) = qi_srcrs(k,m,1) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice )

                  ! -------------------------------------------------------------- !
                  ! Treatment of detrained airs purely from the convective updraft !
                  ! -------------------------------------------------------------- !
                  f_srcrs2(k,m,1)   =   zmass_up * fmix(m) * dpa(m) * eps0(m) * 2._r8 * cmf_au(m) / cmf_u(km)
                  thl_srcrs2(k,m,1) =   thl_au(m)
                  qt_srcrs2(k,m,1)  =   qt_au(m) 
                  u_srcrs2(k,m,1)   =   u_au(m)  
                  v_srcrs2(k,m,1)   =   v_au(m)  
                  do mt = 1, ncnst
                     tr_srcrs2(k,m,1,mt) = tr_au(m,mt)
                  enddo
                  ql_srcrs2(k,m,1)  = ql_au(m)
                  qi_srcrs2(k,m,1)  = qi_au(m)
                  tr_srcrs2(k,m,1,ixnumliq) = ql_srcrs2(k,m,1) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq )
                  tr_srcrs2(k,m,1,ixnumice) = qi_srcrs2(k,m,1) * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice )

               enddo

               ! ---------------- !
               ! 2. Top Downdraft !
               ! ---------------- !

               if( nseg_det .gt. 0.5_r8 ) then
                  do m = 1, N_up(km)
                     if( ytop(m) .lt. 0.5_r8 ) then
                        f_srcds(k,m,2)   = cmf_aut(m) / cmf_u(km)
                        thl_srcds(k,m,2) = thl_aut(m)
                        qt_srcds(k,m,2)  =  qt_aut(m)
                        u_srcds(k,m,2)   =   u_aut(m)
                        v_srcds(k,m,2)   =   v_aut(m)
                        do mt = 1, ncnst
                           tr_srcds(k,m,2,mt) = tr_aut(m,mt)
                        enddo
                        ql_srcds(k,m,2)  =  ql_aut(m)
                        qi_srcds(k,m,2)  =  qi_aut(m)
                     endif
                  enddo
               endif

               ! ----------------- !
               ! 3. Area downdraft !
               ! ----------------- !

               if( nseg_nondet .gt. 0.5_r8 ) then

                  do m = 1, N_up(km)
                     if( ytop(m) .gt. 0.5_r8 ) then
                        thl_srcds(k,m,3) = thl_aut(m)
                        qt_srcds(k,m,3)  =  qt_aut(m)
                        u_srcds(k,m,3)   =   u_aut(m)
                        v_srcds(k,m,3)   =   v_aut(m)
                        do mt = 1, ncnst
                           tr_srcds(k,m,3,mt) = tr_aut(m,mt)
                        enddo
                        ql_srcds(k,m,3)  =  ql_aut(m)
                        qi_srcds(k,m,3)  =  qi_aut(m)
                     endif
                  enddo
               endif

               ! ----------------------------------------------------------------------- !
               ! Mass-flux weighted average of 3 sources of downdraft and 1 detrained    !
               ! airs originated from the convective updrafts.                           !
               ! Jul.15.2010. Since I set ybot(m)=1, only ids=1 has non-zero mass flux   !
               ! of detrained remaining airs.                                            !
               ! ----------------------------------------------------------------------- !

               do m = 1, N_up(km)
                  do ids = 1, 3
                     if( f_srcds(k,m,ids) .gt. nonzero ) then
                        f_srcd(k)   =   f_srcd(k) + f_srcds(k,m,ids) 
                        thl_srcd(k) = thl_srcd(k) + f_srcds(k,m,ids) * thl_srcds(k,m,ids)
                        qt_srcd(k)  =  qt_srcd(k) + f_srcds(k,m,ids) *  qt_srcds(k,m,ids)
                        u_srcd(k)   =   u_srcd(k) + f_srcds(k,m,ids) *   u_srcds(k,m,ids)
                        v_srcd(k)   =   v_srcd(k) + f_srcds(k,m,ids) *   v_srcds(k,m,ids)
                        do mt = 1, ncnst
                           tr_srcd(k,mt) = tr_srcd(k,mt) + f_srcds(k,m,ids) * tr_srcds(k,m,ids,mt) 
                        enddo
                        ql_srcd(k)  =  ql_srcd(k) + f_srcds(k,m,ids) *  ql_srcds(k,m,ids)
                        qi_srcd(k)  =  qi_srcd(k) + f_srcds(k,m,ids) *  qi_srcds(k,m,ids)
                     endif
                     if( f_srcrs(k,m,ids) .gt. nonzero ) then
                        f_srcr(k)   =   f_srcr(k) + f_srcrs(k,m,ids) 
                        thl_srcr(k) = thl_srcr(k) + f_srcrs(k,m,ids) * thl_srcrs(k,m,ids)
                        qt_srcr(k)  =  qt_srcr(k) + f_srcrs(k,m,ids) *  qt_srcrs(k,m,ids)
                        u_srcr(k)   =   u_srcr(k) + f_srcrs(k,m,ids) *   u_srcrs(k,m,ids)
                        v_srcr(k)   =   v_srcr(k) + f_srcrs(k,m,ids) *   v_srcrs(k,m,ids)
                        ql_srcr(k)  =  ql_srcr(k) + f_srcrs(k,m,ids) * ql_srcrs(k,m,ids)
                        qi_srcr(k)  =  qi_srcr(k) + f_srcrs(k,m,ids) * qi_srcrs(k,m,ids)
                        do mt = 1, ncnst
                           tr_srcr(k,mt) = tr_srcr(k,mt) + f_srcrs(k,m,ids) * tr_srcrs(k,m,ids,mt)
                        enddo
                     endif
                     ! --------------------------------------------------------------- !
                     ! Treatment of detrained airs purely from the convective updrafts !
                     ! --------------------------------------------------------------- !
                     if( f_srcrs2(k,m,ids) .gt. nonzero ) then
                        f_srcr2(k)   =   f_srcr2(k) + f_srcrs2(k,m,ids) 
                        thl_srcr2(k) = thl_srcr2(k) + f_srcrs2(k,m,ids) * thl_srcrs2(k,m,ids)
                        qt_srcr2(k)  =  qt_srcr2(k) + f_srcrs2(k,m,ids) *  qt_srcrs2(k,m,ids)
                        u_srcr2(k)   =   u_srcr2(k) + f_srcrs2(k,m,ids) *   u_srcrs2(k,m,ids)
                        v_srcr2(k)   =   v_srcr2(k) + f_srcrs2(k,m,ids) *   v_srcrs2(k,m,ids)
                        ql_srcr2(k)  =  ql_srcr2(k) + f_srcrs2(k,m,ids) *  ql_srcrs2(k,m,ids)
                        qi_srcr2(k)  =  qi_srcr2(k) + f_srcrs2(k,m,ids) *  qi_srcrs2(k,m,ids)
                        do mt = 1, ncnst
                           tr_srcr2(k,mt) = tr_srcr2(k,mt) + f_srcrs2(k,m,ids) * tr_srcrs2(k,m,ids,mt)
                        enddo
                     endif
                  enddo

               enddo

               if( f_srcd(k) .gt. nonzero ) then  
                  thl_srcd(k) = thl_srcd(k) / f_srcd(k)
                  qt_srcd(k)  =  qt_srcd(k) / f_srcd(k)
                  u_srcd(k)   =   u_srcd(k) / f_srcd(k)
                  v_srcd(k)   =   v_srcd(k) / f_srcd(k)
                  do mt = 1, ncnst
                     tr_srcd(k,mt) = tr_srcd(k,mt) / f_srcd(k)
                  enddo
                  ql_srcd(k)  =  ql_srcd(k) / f_srcd(k)
                  qi_srcd(k)  =  qi_srcd(k) / f_srcd(k)
               else
                  f_srcd(k)   = 0._r8
                  thl_srcd(k) = 0._r8
                  qt_srcd(k)  = 0._r8
                  u_srcd(k)   = 0._r8
                  v_srcd(k)   = 0._r8
                  do mt = 1, ncnst
                     tr_srcd(k,mt) = 0._r8
                  enddo
                  ql_srcd(k)  = 0._r8
                  qi_srcd(k)  = 0._r8
               endif
               if( f_srcr(k) .gt. nonzero ) then  
                  thl_srcr(k) =  thl_srcr(k) / f_srcr(k)
                  qt_srcr(k)  =   qt_srcr(k) / f_srcr(k)
                  u_srcr(k)   =    u_srcr(k) / f_srcr(k)
                  v_srcr(k)   =    v_srcr(k) / f_srcr(k)
                  ql_srcr(k)  =   ql_srcr(k) / f_srcr(k)
                  qi_srcr(k)  =   qi_srcr(k) / f_srcr(k)
                  do mt = 1, ncnst
                     tr_srcr(k,mt) = tr_srcr(k,mt) / f_srcr(k)
                  enddo
               else
                  f_srcr(k)   = 0._r8
                  thl_srcr(k) = 0._r8
                  qt_srcr(k)  = 0._r8
                  u_srcr(k)   = 0._r8
                  v_srcr(k)   = 0._r8
                  ql_srcr(k)  = 0._r8
                  qi_srcr(k)  = 0._r8
                  do mt = 1, ncnst
                     tr_srcr(k,mt) = 0._r8
                  enddo
               endif
               ! ------------------------------------------------------------- !
               ! Treatment of detrained air purely from the convective updraft !
               ! ------------------------------------------------------------- !
               if( f_srcr2(k) .gt. nonzero ) then  
                  thl_srcr2(k) =  thl_srcr2(k) / f_srcr2(k)
                  qt_srcr2(k)  =   qt_srcr2(k) / f_srcr2(k)
                  u_srcr2(k)   =    u_srcr2(k) / f_srcr2(k)
                  v_srcr2(k)   =    v_srcr2(k) / f_srcr2(k)
                  ql_srcr2(k)  =   ql_srcr2(k) / f_srcr2(k)
                  qi_srcr2(k)  =   qi_srcr2(k) / f_srcr2(k)
                  do mt = 1, ncnst
                     tr_srcr2(k,mt) = tr_srcr2(k,mt) / f_srcr2(k)
                  enddo
               else
                  f_srcr2(k)   = 0._r8
                  thl_srcr2(k) = 0._r8
                  qt_srcr2(k)  = 0._r8
                  u_srcr2(k)   = 0._r8
                  v_srcr2(k)   = 0._r8
                  ql_srcr2(k)  = 0._r8
                  qi_srcr2(k)  = 0._r8
                  do mt = 1, ncnst
                     tr_srcr2(k,mt) = 0._r8
                  enddo
               endif

               ! --------------------------------------------------------------------------- !
               ! Allocation of detrained source airs from convective updraft into new arrays !
               ! --------------------------------------------------------------------------- !

               cmf_ru(k)  =     f_srcr(k) * cmf_u(km)
               thl_ru(k)  =   thl_srcr(k)
               qt_ru(k)   =    qt_srcr(k)
               u_ru(k)    =     u_srcr(k)
               v_ru(k)    =     v_srcr(k)
               ql_ru(k)   =    ql_srcr(k)  
               qi_ru(k)   =    qi_srcr(k)
               do mt = 1, ncnst
                  tr_ru(k,mt) = tr_srcr(k,mt)           
               enddo

               ! ------------------------------------------------------------- !
               ! Treatment of detrained air purely from the convective updraft !
               ! ------------------------------------------------------------- !

               cmf_ru2(k)  =     f_srcr2(k) * cmf_u(km)
               thl_ru2(k)  =   thl_srcr2(k)
               qt_ru2(k)   =    qt_srcr2(k)
               u_ru2(k)    =     u_srcr2(k)
               v_ru2(k)    =     v_srcr2(k)
               ql_ru2(k)   =    ql_srcr2(k)  
               qi_ru2(k)   =    qi_srcr2(k)
               do mt = 1, ncnst
                  tr_ru2(k,mt) = tr_srcr2(k,mt)           
               enddo
               
               ! ---------------------------------------------------------------------- !
               ! Compute 'cloud fraction' and 'in-cloud LWC,IWC' at the layer mid-point !
               ! for individual updraft segment.                                        ! 
               ! Note that 'ktop_msfc', 'ptop_msfc', 'ztop_msfc' should be printed out  !
               ! at the output side of vertical 'k' loop - that is, it should be        !
               ! printed-out where 'cushavg' is computed.                               !
               ! Sep.15.2011. All the other conservative scalars too for defining       !
               !              mixing environmental airs associated organization later.  !
               ! ---------------------------------------------------------------------- !
 
               do m = 1, N_up(km)
                  msfc = msfc_from_m(k,m)
                  am_u_msfc(k,msfc)         =  0.5_r8 * (     a_au(m) +     a_aut(m) ) * ( dpa(m) / dp_m )
                  qlm_u_msfc(k,msfc)        =  0.5_r8 * (    ql_au(m) +    ql_aut(m) ) * ( dpa(m) / dp_m )
                  qim_u_msfc(k,msfc)        =  0.5_r8 * (    qi_au(m) +    qi_aut(m) ) * ( dpa(m) / dp_m )
                  thlm_u_msfc(k,msfc)       =  0.5_r8 * (   thl_au(m) +   thl_aut(m) ) * ( dpa(m) / dp_m )
                  qtm_u_msfc(k,msfc)        =  0.5_r8 * (    qt_au(m) +    qt_aut(m) ) * ( dpa(m) / dp_m )
                  um_u_msfc(k,msfc)         =  0.5_r8 * (     u_au(m) +     u_aut(m) ) * ( dpa(m) / dp_m )
                  vm_u_msfc(k,msfc)         =  0.5_r8 * (     v_au(m) +     v_aut(m) ) * ( dpa(m) / dp_m )
                  do mt = 1, ncnst
                     trm_u_msfc(k,msfc,mt)  =  0.5_r8 * ( tr_au(m,mt) + tr_aut(m,mt) ) * ( dpa(m) / dp_m )
                  enddo
                  am_u(k)        =     am_u(k) +  am_u_msfc(k,msfc)
                  qlm_u(k)       =    qlm_u(k) +  am_u_msfc(k,msfc) *    qlm_u_msfc(k,msfc)
                  qim_u(k)       =    qim_u(k) +  am_u_msfc(k,msfc) *    qim_u_msfc(k,msfc)
                  thlm_u(k)      =   thlm_u(k) +  am_u_msfc(k,msfc) *   thlm_u_msfc(k,msfc)
                  qtm_u(k)       =    qtm_u(k) +  am_u_msfc(k,msfc) *    qtm_u_msfc(k,msfc)
                  um_u(k)        =     um_u(k) +  am_u_msfc(k,msfc) *     um_u_msfc(k,msfc)
                  vm_u(k)        =     vm_u(k) +  am_u_msfc(k,msfc) *     vm_u_msfc(k,msfc)
                  do mt = 1, ncnst
                     trm_u(k,mt) = trm_u(k,mt) +  am_u_msfc(k,msfc) * trm_u_msfc(k,msfc,mt)
                  enddo
               enddo
               if( am_u(k) .gt. nonzero ) then
                  qlm_u(k)       =    qlm_u(k) / am_u(k)
                  qim_u(k)       =    qim_u(k) / am_u(k)
                  thlm_u(k)      =   thlm_u(k) / am_u(k)
                  qtm_u(k)       =    qtm_u(k) / am_u(k)
                  um_u(k)        =     um_u(k) / am_u(k)
                  vm_u(k)        =     vm_u(k) / am_u(k)
                  do mt = 1, ncnst
                     trm_u(k,mt) = trm_u(k,mt) / am_u(k)
                  enddo
               else
                  ! Sep.16.2011. Below is not anomaly but total field. Thus, in order to reduce any 
                  !              potential bias grow later, I entered grid-mean value instead of zero.
                  qlm_u(k)       =      0._r8 
                  qim_u(k)       =      0._r8 
                  thlm_u(k)      =    thl0(k)
                  qtm_u(k)       =     qt0(k) 
                  um_u(k)        =      u0(k) 
                  vm_u(k)        =      v0(k)
                  do mt = 1, ncnst
                     trm_u(k,mt) =  tr0(k,mt)
                  enddo
               endif

               ! ----------------------------------------------------------------------------- !
               ! Jul.26.2011.                                                                  !
               ! Save all the convective updraft properties at each model interface and at the ! 
               ! cumulus top ( i.e., in the cumulus top layer in each 'm' segment, the value   !
               ! at the cumulus top not the value at the top interface of cumulus top layer is !
               ! saved ) for future use. For example, this can be used as alternative mixing   !
               ! environmental value for convective downdraft.                                 !
               ! Note that at the cumulus top, I set 'w_aut(m)=0, cmf_aut(m)>0' and            !
               ! 'a_aut(m) = a_au(m), num_aut(m) = num_au(m), rad_aut(m) = rad_au(m)'.         !  
               ! ----------------------------------------------------------------------------- !

               do m = 1, N_up(km)
                  msfc = msfc_from_m(k,m)
                  thl_u_msfc(k,msfc)           =   thl_aut(m)
                  qt_u_msfc(k,msfc)            =    qt_aut(m)       
                  u_u_msfc(k,msfc)             =     u_aut(m)       
                  v_u_msfc(k,msfc)             =     v_aut(m)       
                  w_u_msfc(k,msfc)             =     w_aut(m)       
                  ql_u_msfc(k,msfc)            =    ql_aut(m)
                  qi_u_msfc(k,msfc)            =    qi_aut(m)
                  do mt = 1, ncnst
                     tr_u_msfc(k,msfc,mt)      = tr_aut(m,mt)
                  enddo
                  cmf_u_msfc(k,msfc)           =   cmf_aut(m)
                  a_u_msfc(k,msfc)             =     a_aut(m)
                  num_u_msfc(k,msfc)           =   num_aut(m)
                  rad_u_msfc(k,msfc)           =   rad_aut(m)
                  if( k .eq. 1 ) then
                     thl_u_msfc(km,msfc)      =    thl_au(m)
                     qt_u_msfc(km,msfc)       =     qt_au(m)       
                     u_u_msfc(km,msfc)        =      u_au(m)       
                     v_u_msfc(km,msfc)        =      v_au(m)       
                     w_u_msfc(km,msfc)        =      w_au(m)       
                     ql_u_msfc(km,msfc)       =     ql_au(m)
                     qi_u_msfc(km,msfc)       =     qi_au(m)
                     do mt = 1, ncnst
                        tr_u_msfc(km,msfc,mt) =  tr_au(m,mt)
                     enddo
                     cmf_u_msfc(km,msfc)      =    cmf_au(m)
                     a_u_msfc(km,msfc)        =      a_au(m)
                     num_u_msfc(km,msfc)      =    num_au(m)
                     rad_u_msfc(km,msfc)      =    rad_au(m)
                  endif
               enddo

               ! -------------------------------------------------------------------- !
               ! Re-allocate updraft segment values for computation in the next layer !
               ! Assign the array only to the non-detached updraft.                   !
               ! Mar.12.2013. Add 'S_b_ql_au(mm),S_b_qi_au(mm)' parts.                ! 
               ! -------------------------------------------------------------------- !

               mm = 0 
               do m = 1, N_up(km)
                  if( ytop(m) .gt. 0.5_r8 ) then
                     mm = mm + 1
                     cmf_au(mm) = cmf_aut(m)
                     a_au(mm)   =   a_aut(m)
                     num_au(mm) = num_aut(m)
                     rad_au(mm) = rad_aut(m)
                     thl_au(mm) = thl_aut(m)
                     qt_au(mm)  =  qt_aut(m)
                     u_au(mm)   =   u_aut(m)
                     v_au(mm)   =   v_aut(m)
                     w_au(mm)   =   w_aut(m)
                     ql_au(mm)  =  ql_aut(m)
                     qi_au(mm)  =  qi_aut(m)
                     thv_au(mm) = thv_aut(m)
                     do mt = 1, ncnst
                        tr_au(mm,mt) = tr_aut(m,mt)    
                     enddo
                     S_b_ql_au(mm) =  S_t_ql_au(m)
                     S_b_qi_au(mm) =  S_t_qi_au(m)
                  endif
               enddo

               ! -------------------------- !
               ! Identify Cumulus Top Layer !
               ! -------------------------- !

               if( nseg_nondet .lt. 0.5_r8 ) then
                  ktop = k
                  cnt  = real(k,r8)
                  ! --------------------------------------------------------------------------------- !
                  ! Aug.02.2011. It seems to be more reasonable to set cnb  = real(0,r8) instead of 1 !
                  !              This should be done later.                                           !
                  ! --------------------------------------------------------------------------------- ! 
                  cnb  = real(1,r8)
                  cush_mxen(iter) = pblhz 
                  do mm = 1, N_up(ktop-1)
                     cush_mxen(iter) = max( cush_mxen(iter), ztops(ktop,mm) )
                  enddo
                  goto 50
               endif

            enddo ! k = 1, mkx - 1. Here, 'k' is a layer index.

50          continue

            ! ------------------------------------------------------ !
            ! Assign updraft top layer index to 'ktop_mxen' variable !
            ! ------------------------------------------------------ !

            ktop_mxen(iter) = ktop

            ! --------------------------------------------------------------------- !
            ! Compute mean-cumulus top height weighted by updraft mass flux         !
            ! at surface. This quantity will be used for computing cridis_in        !
            ! at the next time step instead of cush.                                !
            ! Sep.22.2011. Note that 'Pmu(m)' does not include the 'delta_w_PBL'.   !
            !              However, since 'delta_w_PBL' is added uniformly all over !
            !              the 'm' segments, below computation of 'cushavg_mxen' is !
            !              completely correct. Note that below block is the only    !
            !              part of the whole program to explicitly use the 'Pmu'    !
            !              except the computation of 'cmf_au(m)' at surface.        !  
            ! Sep.21.2011. Now, by re-defining 'Pmu(m)' above, below computation of !
            !              cushavg_mxen is perfectly correct                        !
            !              without any modification.                                !
            ! --------------------------------------------------------------------- !

            tmp1                  = 0._r8 
            cushavg_mxen(iter)    = 0._r8
            do m = 1, nseg
               tmp1               = tmp1               + Pmu(m)
               cushavg_mxen(iter) = cushavg_mxen(iter) + ztop_msfc(m) * Pmu(m)
            enddo
            cushavg_mxen(iter)    = max( cushavg_mxen(iter) / tmp1, pblhz )
                    
            ! ------------------------------------------- !
            !                                             !
            ! Vertical Evolution of Individual Downdrafts !
            !                                             !
            ! ------------------------------------------- !                

            rbuoy_dn = 0.5_r8 * ( rbuoy_min + rbuoy_max ) 

            do msfc = 1, nseg                         ! This 'msfc' is updraft segment index at surface.

               ! -------------------------------------------------------------------------- !
               ! Initialization of Downdraft Sources in Each Layer for Each Updraft Segment !
               ! Mar.11.2013. Add initialization of evaporation rate at the top interface.  ! 
               ! -------------------------------------------------------------------------- !
  
               cmf_d_src(1:mkx,1:3)        = 0._r8
               thl_d_src(1:mkx,1:3)        = 0._r8
               qt_d_src(1:mkx,1:3)         = 0._r8
               u_d_src(1:mkx,1:3)          = 0._r8
               v_d_src(1:mkx,1:3)          = 0._r8
               w_d_src(1:mkx,1:3)          = 0._r8
               tr_d_src(1:mkx,1:3,1:ncnst) = 0._r8
               ql_d_src(1:mkx,1:3)         = 0._r8 
               qi_d_src(1:mkx,1:3)         = 0._r8
               fevp1_t_rate_src(1:mkx,1:3) = 0._r8
               fevp2_t_rate_src(1:mkx,1:3) = 0._r8
               ix_d_src(1:mkx,1:3)         = 0

               ! ------------------------------------------------------------------------------------------------------- !
               ! Compute the location (x_um_msfc(k,msfc),y_um_msfc(k,msfc) in unit of [m]) of convective updraft center  !
               ! at the layer mid-point relative to the convective updraft center at surface                             !
               ! using the profiles of (u,v,w) of convective updraft.                                                    !
               ! Note that 'u_u_msfc, v_u_msfc, w_u_msfc' in the k = ktop_msfc is defined at the cumulus                 !
               ! top not at the top interface. Thus, my below computation is perfectly correct.                          !
               ! Note that this computation goes from the bottom to the cumulus top layer, so that I should use separate !
               ! vertical loop here from the lowest to the cumulus top layer.                                            !
               ! ------------------------------------------------------------------------------------------------------- !
               tmpx_bot = 0._r8
               tmpy_bot = 0._r8
               do k = 1, ktop_msfc(msfc)  ! This is a layer index
                  km = k - 1
                  if( k .eq. ktop_msfc(msfc) ) then
                     tmp3 = ztop_msfc(msfc) - zs0(km)
                  else
                     tmp3 = dz0(k)
                  endif
                  tmp1 = tmp3 / ( 0.5_r8 * ( w_u_msfc(k,msfc) + w_u_msfc(km,msfc) ) )
                  tmpx_top = tmpx_bot + 0.5_r8 * ( u_u_msfc(k,msfc) - u_u_msfc(km,msfc) ) * tmp1
                  tmpy_top = tmpy_bot + 0.5_r8 * ( v_u_msfc(k,msfc) - v_u_msfc(km,msfc) ) * tmp1
                  x_um_msfc(k,msfc) = 0.5_r8 * ( tmpx_bot + tmpx_top )
                  y_um_msfc(k,msfc) = 0.5_r8 * ( tmpy_bot + tmpy_top )
                  tmpx_bot = tmpx_top
                  tmpy_bot = tmpy_top
               end do

               ! --------------------------------------------------- !
               !                                                     ! 
               ! Start of Downdraft Vertical Evolution in Each Layer !
               !                                                     !
               ! --------------------------------------------------- !  

               do k   =  ktop_msfc(msfc), 1, -1                     ! This 'k' is a layer index where vertical evolution of downdraft is computed.

                  km  =  k - 1

                  ! ---------------------------------------------------------------------------------------------------- !
                  ! Define non-array 'rain/snow/tracer fluxes at the top interface for use in various computations below !
                  ! both for 'grid-mean' and 'in-precipitation-area' values.                                             !
                  ! ---------------------------------------------------------------------------------------------------- !

                  flxrain_top    = flxrain_msfc(k,msfc)
                  flxsnow_top    = flxsnow_msfc(k,msfc)
                  flxrasn_top    = flxrain_top + flxsnow_top

                  flxrain_top_in = flxrain_top / max( nonzero, a_p_msfc(k,msfc) )
                  flxsnow_top_in = flxsnow_top / max( nonzero, a_p_msfc(k,msfc) )
                  flxrasn_top_in = flxrain_top_in + flxsnow_top_in

                  do mt = 1, ncnst
                     if( mt .eq. ixcldliq ) then
                        flxtrrs_top(mt) = flxrain_top
                     elseif( mt .eq. ixcldice ) then
                        flxtrrs_top(mt) = flxsnow_top
                     elseif( mt .eq. ixnumliq ) then
                        flxtrrs_top(mt) = flxrain_top * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_rain**3 * density_rain )
                     elseif( mt .eq. ixnumice ) then
                        flxtrrs_top(mt) = flxsnow_top * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_snow**3 * density_snow )
                     else
                        flxtrrs_top(mt) = flxtrrs_msfc(k,msfc,mt)
                     endif
                  enddo

                  ! -------------------------------------------------------------------------------------------- !
                  ! Compute various areas to compute                                                             !
                  !  1. Production  of convective precipitation by accretion within convective updraft           !
                  !  2. Evaporation of convective precipitation within environment                               !
                  ! Below simply assumes that downdraft fractional area is zero for this purpose ( a_pd = 0. ).  !
                  ! In order to compute accretion rate, we should use                                            ! 
                  !       (a) a_pu : Overlapping area between precipitation area ( a_p_msfc(k,msfc) ) and        !
                  !                  saturated updraft fractional area ( am_us_msfc(k,msfc) ),                   !
                  !       (b) flxrain_top_in, flxsnow_top_in : Precipitation flux averaged over the              !
                  !                  precipitation area (not the grid-mean) at the top interface,                !
                  !       (c) qlm_u_msfc(k,msfc), qim_u_msfc(k,msfc), trm_u_msfc(k,msfc,mt) :                    !
                  !                  In-cumulus properties.                                                      !
                  ! Note 'a_p = a_pu + a_pd + a_pr + a_ps' and independent precipitation approximation is used,  !
                  ! so that in computing these overlapping areas, we really don't need to worry about the other  !
                  ! cumulus updraft segment' contribution.                                                       !     
                  !                                                                                              !
                  ! Feb.09.2013. In computing 'a_pu' below, change 'am_us_msfc' to 'am_u_msfc' since evaporation !
                  !              of precipitation is computed within the subroutine 'prod_prep_up' above.        !
                  !              Thus, here, we should only consider what is happening within environment and    !
                  !              downdraft, not within updraft regardless whether updraft is saturated or not.   !                                              
                  ! -------------------------------------------------------------------------------------------- !

                  am_s(k) = ast0(k)
                  if( ( am_u(k) + ast0(k) ) .gt. 1._r8 ) am_s(k) = 1._r8 - am_u(k)
                  am_r(k) = 1._r8 - am_u(k) - am_s(k)

                  am_us_msfc(k,msfc) = am_u_msfc(k,msfc)
                  if( ( qlm_u_msfc(k,msfc) + qim_u_msfc(k,msfc) ) .le. 0._r8 ) am_us_msfc(k,msfc) = 0._r8

                  a_pu  = area_overlap( x_p_msfc(k,msfc),  y_p_msfc(k,msfc),  a_p_msfc(k,msfc),   &
                     !j   x_um_msfc(k,msfc), y_um_msfc(k,msfc), am_u_msfc(k,msfc),  &
                          x_um_msfc(k,msfc), y_um_msfc(k,msfc), am_us_msfc(k,msfc), &
                          num_u_msfc(0,msfc) )
                  a_pd  = 0._r8
                  a_pr  = ( a_p_msfc(k,msfc) - a_pu - a_pd ) * am_r(k) / ( am_r(k) + am_s(k) )    ! Advanced.
                  a_ps  = ( a_p_msfc(k,msfc) - a_pu - a_pd ) * am_s(k) / ( am_r(k) + am_s(k) )    ! Advanced.
                  a_evp = a_pr

                  am_evp_msfc(k,msfc) = a_evp
                  am_pu_msfc(k,msfc)  = a_pu
                  am_pd_msfc(k,msfc)  = a_pd
                  am_pr_msfc(k,msfc)  = a_pr
                  am_ps_msfc(k,msfc)  = a_ps

                  ! ---------------------------------------------------------------------------------------------------------------------- !
                  ! Compute evaporation of rain/snow within environment ( evprain_e_msfc(k,msfc), evpsnow_e_msfc(k,msfc) >= 0. [kg/kg/s] ) !
                  ! ---------------------------------------------------------------------------------------------------------------------- !

                ! Apr.15.2014. Recalculate 'qv_clr' using the normalized stratus fraction instead of using
                !              the below original formula based on physical stratus fraction.
                  call qsat( t0(k), p0(k), es, qs )
                  tmp1         = am_s(k) / ( 1._r8 - am_u(k) )
                  qv_clr       = max( nonzero, qv0(k) - tmp1 * qs ) / max( nonzero, 1._r8 - tmp1 )
                ! qv_clr       = max( nonzero, qv0(k) - am_s(k) * qs ) / max( nonzero, 1._r8 - am_s(k) )
                  qv_clr       = min( min( qv_clr, qv0(k) ), qs )
                  subsat_clr   = min( 1._r8, max( 0._r8, 1._r8 - qv_clr / max( qs, nonzero ) ) )
                  call findsp_single( qv_clr, t0(k), p0(k), tw, qw_clr, i, k, lchnk )
                  evplimit_clr = max( 0._r8, ( qw_clr - qv_clr ) / dt )

                  evprain_clr  = kevp_rain * subsat_clr * sqrt( max( 0._r8, flxrain_top_in ) )
                  evpsnow_clr  = kevp_snow * subsat_clr * sqrt( max( 0._r8, flxsnow_top_in ) )         

                  evplimit_clr_rain = flxrain_top_in * g / dp0(k) ! New. Perfect.
                  evplimit_clr_snow = flxsnow_top_in * g / dp0(k) ! New. Perfect.

                  evprain_clr       = min(  evprain_clr, evplimit_clr_rain )
                  evpsnow_clr       = min(  evpsnow_clr, evplimit_clr_snow )
                  if( ( evprain_clr + evpsnow_clr ) .lt. ( - qv0(k) / dt / max( nonzero, a_evp ) ) ) then
                     call endrun('UNICON : Impossible correction of precipitation generation')
                  endif
                  if( ( evprain_clr + evpsnow_clr ) .gt. evplimit_clr ) then
                     if( evprain_clr .ge. 0._r8 .and. evpsnow_clr .ge. 0._r8 ) then
                        tmp1 = evprain_clr * evplimit_clr / ( evprain_clr + evpsnow_clr )
                        tmp2 = evpsnow_clr * evplimit_clr / ( evprain_clr + evpsnow_clr )
                        evprain_clr = tmp1
                        evpsnow_clr = tmp2
                     elseif( evprain_clr .lt. 0._r8 ) then
                        evpsnow_clr = evplimit_clr - evprain_clr
                     elseif( evpsnow_clr .lt. 0._r8 ) then
                        evprain_clr = evplimit_clr - evpsnow_clr
                     else
                        call endrun('UNICON : Impossible case in Limit 1a')
                     endif
                  endif

                  evprain_e_msfc(k,msfc)  = evprain_clr * a_evp
                  evpsnow_e_msfc(k,msfc)  = evpsnow_clr * a_evp
                  do mt = 1, ncnst
                     if( mt .eq. ixcldliq ) then
                        evptrrs_e_msfc(k,msfc,mt) = - evprain_e_msfc(k,msfc)
                     elseif( mt .eq. ixcldice ) then
                        evptrrs_e_msfc(k,msfc,mt) = - evpsnow_e_msfc(k,msfc)
                     elseif( mt .eq. ixnumliq ) then
                        evptrrs_e_msfc(k,msfc,mt) = - evprain_e_msfc(k,msfc) * 3._r8 / &
                                                    ( 4._r8 * 3.141592_r8 * droprad_rain**3 * density_rain )
                     elseif( mt .eq. ixnumice ) then
                        evptrrs_e_msfc(k,msfc,mt) = - evpsnow_e_msfc(k,msfc) * 3._r8 / &
                                                    ( 4._r8 * 3.141592_r8 * droprad_snow**3 * density_snow )
                     else
                        evptrrs_e_msfc(k,msfc,mt) = - flxtrrs_top(mt) * ( ( evprain_e_msfc(k,msfc) + &
                                                    evpsnow_e_msfc(k,msfc) ) / max( flxrasn_top, nonzero ) )
                     endif
                  enddo

                  ! ------------------------------------------------------------------------------------ !
                  ! Compute wet deposition of aerosols within denvironment.                              ! 
                  ! Feb.05.2013. Below wet deposition tendency within updraft should be updated later    !
                  !              possibly, in combination with the treatment of accretion process.       !
                  !              Note that wet deposition only influences tracers not cloud condensate.  !
                  !              If I turn-on this in future, wet deposition by convective precipitation !
                  !              in the separate wet deposition routine should be turned-off.            !
                  ! ------------------------------------------------------------------------------------ !

                  do mt = 1, ncnst
                     if( mt .eq. 1 .or. mt .eq. ixcldliq .or. mt .eq. ixcldice .or. mt .eq. ixnumliq .or. mt .eq. ixnumice ) then
                        wdeptrrs_e_msfc(k,msfc,mt) = 0._r8
                     else
                        wdeptrrs_e_msfc(k,msfc,mt) = 0._r8
                     endif
                  enddo

                  ! --------------------------------------------------------------------------------- !
                  ! Compute                                                                           !
                  ! Precipitation flux at the base interface by adding evaporation within environment !
                  ! Note that until the full 2-moment microphysics are implemented, I will            !
                  ! assume a fixed droplet size of rain and snow.                                     !
                  ! Note that wet deposition does not affect condensate but only influences tracers.  !
                  ! Mar.05.2013. For computing the location of precipitation area (x_p, y_p) at the   ! 
                  !              base interface and for diagnostic purpose, compute                   !
                  !              the 'flxrain(snow)_bot_up' and 'flxrain(snow)_bot_ee'.               !   
                  !              Note that it is guaranteed that 'flxrain(snow)_bot_ee >= 0'          !
                  !              since 'evprain(snow)_e_msfc' were computed from 'flxrain_top_in'     !
                  !              which is perfectly good.                                             ! 
                  ! --------------------------------------------------------------------------------- !

                  flxrain_bot_ee   = flxrain_top - evprain_e_msfc(k,msfc) * ( dp0(k) / g )
                  flxsnow_bot_ee   = flxsnow_top - evpsnow_e_msfc(k,msfc) * ( dp0(k) / g )

                  flxrain_bot_upee = flxrain_top + ( qrten_u_msfc(k,msfc) - evprain_e_msfc(k,msfc) ) * ( dp0(k) / g )
                  flxsnow_bot_upee = flxsnow_top + ( qsten_u_msfc(k,msfc) - evpsnow_e_msfc(k,msfc) ) * ( dp0(k) / g )

                  do mt = 1, ncnst
                     if( mt .eq. ixcldliq ) then
                        flxtrrs_bot_upee(mt) = flxrain_bot_upee
                     elseif( mt .eq. ixcldice ) then
                        flxtrrs_bot_upee(mt) = flxsnow_bot_upee
                     elseif( mt .eq. ixnumliq ) then
                        flxtrrs_bot_upee(mt) = flxrain_bot_upee * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_rain**3 * density_rain )
                     elseif( mt .eq. ixnumice ) then
                        flxtrrs_bot_upee(mt) = flxsnow_bot_upee * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_snow**3 * density_snow )
                     else
                        flxtrrs_bot_upee(mt) = flxtrrs_top(mt) + ( trrsten_u_msfc(k,msfc,mt) + evptrrs_e_msfc(k,msfc,mt) + &
                                               wdeptrrs_e_msfc(k,msfc,mt) ) * ( dptr0(k,mt) / g )
                     endif
                  enddo

                  ! ------------------------------------------------------------------------------------------------------- !
                  ! Snow Melting at the Base Interface ( snowmlt_e_msfc(k,msfc) >= 0. [kg/kg/s] )                           !
                  ! Since snow melting changes droplet number within precipitation, I should perform updated computation of ! 
                  ! droplet numbers and mass of precipitation.                                                              !
                  ! Note that below assume that snow melting does not change the other tracers concentration as expected.   !
                  ! ------------------------------------------------------------------------------------------------------- !

                  if( t0(k) .ge. 273.15_r8 ) then
                     snowmlt_e_msfc(k,msfc) = max( 0._r8, unity * flxsnow_bot_upee * g / dp0(k) )
                  else
                     snowmlt_e_msfc(k,msfc) = 0._r8
                  endif
                  flxrain_bot_upeesm = max( 0._r8, flxrain_bot_upee + snowmlt_e_msfc(k,msfc) * dp0(k) / g )
                  flxsnow_bot_upeesm = max( 0._r8, flxsnow_bot_upee - snowmlt_e_msfc(k,msfc) * dp0(k) / g )
                  do mt = 1, ncnst
                     if( mt .eq. ixcldliq ) then
                        flxtrrs_bot_upeesm(mt) = flxrain_bot_upeesm
                     elseif( mt .eq. ixcldice ) then
                        flxtrrs_bot_upeesm(mt) = flxsnow_bot_upeesm
                     elseif( mt .eq. ixnumliq ) then
                        flxtrrs_bot_upeesm(mt) = flxrain_bot_upeesm * 3._r8 / ( 4._r8 * 3.141592_r8 * &
                                                 droprad_rain**3 * density_rain )
                     elseif( mt .eq. ixnumice ) then
                        flxtrrs_bot_upeesm(mt) = flxsnow_bot_upeesm * 3._r8 / ( 4._r8 * 3.141592_r8 * &
                                                 droprad_snow**3 * density_snow )
                     else
                        flxtrrs_bot_upeesm(mt) = flxtrrs_bot_upee(mt) 
                     endif
                  enddo

                  ! ---------------------------------------------------------------------------------------------------------------------- !
                  ! Compute precipitation area at the base interface ( a_p_msfc(km,msfc) )                                                 !
                  ! While 'the location of the center of precipitation area ( x_p_msfc(km,msfc), y_p_msfc(km,msfc) )' are computed later   !
                  ! after doing evaporation of precipitation within downdraft, I should compute 'a_p_msfc(km,msfc)' here before            !
                  ! computing evaporation of precipitation within downdraft since that computation requires                                !
                  ! the use of 'a_p_msfc(km,msfc)'. This is completely correct approach assuming that evaporation                          !
                  ! within downdraft at the base interface does not completely evaporate 'flxrain_bot_upeesm', which is                    !
                  ! grauanteed by setting 'eta2' smaller than 1.                                                                           !
                  ! ---------------------------------------------------------------------------------------------------------------------- !

                  am_up_msfc(k,msfc) = am_u_msfc(k,msfc)
                  if( ( qrten_u_msfc(k,msfc) + qsten_u_msfc(k,msfc) ) .le. 0._r8 ) am_up_msfc(k,msfc) = 0._r8
                  a_ovp = area_overlap( x_p_msfc(k,msfc),  y_p_msfc(k,msfc),  a_p_msfc(k,msfc),   &
                     x_um_msfc(k,msfc), y_um_msfc(k,msfc), am_up_msfc(k,msfc), &
                     num_u_msfc(0,msfc) )
                  if( ( flxrasn_top_in - nonzero ) .le. ( evprain_clr + evpsnow_clr ) * dp0(k) / g ) then
                     a_p_msfc(km,msfc) = a_p_msfc(k,msfc) - a_evp + am_up_msfc(k,msfc) - a_ovp
                  else
                     a_p_msfc(km,msfc) = a_p_msfc(k,msfc) + am_up_msfc(k,msfc) - a_ovp
                  endif
                  if( ( flxrain_bot_upeesm + flxsnow_bot_upeesm ) .le. 0._r8 ) a_p_msfc(km,msfc) = 0._r8
                  a_p_msfc(km,msfc) = max( 0._r8, min( 1._r8, a_p_msfc(km,msfc) ) )

                  ! ---------------------------------------------------------------------------------------------- !
                  ! Initialize                                                                                     !
                  !                                                                                                !
                  !    flxrain_bot = flxrain_bot_upeesm                                                            !   
                  !    flxsnow_bot = flxsnow_bot_upeesm                                                            !
                  !                                                                                                ! 
                  ! where the final precipitation fluxes ( flxrain_bot, flxsnow_bot ) will be continuously updated !
                  ! within the individual downdraft loop ( do ks = ktop, k, -1, do ids = 1, 3 ) below, eventually  !
                  ! getting the final precipitation flux at the bottom interface.                                  ! 
                  ! ---------------------------------------------------------------------------------------------- !

                  flxrain_bot = flxrain_bot_upeesm
                  flxsnow_bot = flxsnow_bot_upeesm
                  do mt = 1, ncnst
                     flxtrrs_bot(mt) = flxtrrs_bot_upeesm(mt) 
                  enddo

                  ! ------------------------------------------------------- !
                  ! Define downdraft sources generated in the current layer !
                  ! ------------------------------------------------------- !

                  do ids  = 1, 3   ! This 'ids' is the type of downdraft source ( 1 : Mixing downdraft, 2 : Top downdraft, 3 : Constraint downdraft )
                     ! ----------------------------------------------------------------------------------------------------------------- !
                     ! m_from_msfc(k,msfc) : Convert updraft segment index at surface into shortened-updraft segment index in each layer !
                     !                       Since computation is done from k = 'ktop = ktop_msfc(msfc)' to 'k = 1' for each 'msfc',     !
                     !                       it is always grauanteed that 'm' is a positive integer ( m > 0 ). Thus, below computation   !
                     !                       is perfectly OK - if not, it will stop due to indexing error.                               !
                     ! Be careful : I should use index 'k' not 'ks' in the below block since in this case index 'k' denotes source layer !
                     ! ----------------------------------------------------------------------------------------------------------------- !
                     m                = m_from_msfc(k,msfc)
                     ix_d_src(k,ids)  = 1  
                     cmf_d_src(k,ids) = f_srcds(k,m,ids) * cmf_u(km)
                     thl_d_src(k,ids) = thl_srcds(k,m,ids)
                     qt_d_src(k,ids)  = qt_srcds(k,m,ids) 
                     u_d_src(k,ids)   = u_srcds(k,m,ids) 
                     v_d_src(k,ids)   = v_srcds(k,m,ids) 
                     w_d_src(k,ids)   = 0._r8
                     do mt = 1, ncnst
                        tr_d_src(k,ids,mt) = tr_srcds(k,m,ids,mt) 
                     enddo
                     ql_d_src(k,ids)  = ql_srcds(k,m,ids) 
                     qi_d_src(k,ids)  = qi_srcds(k,m,ids) 
                     ! ---------------------------------------------------------------------------------------------- !
                     ! Mar.11.2013. Add initialization of evaporation rate at the top interface                       !
                     !              Since below two variables are already initialized to zero above,                  !
                     !              below initialization is redundant. However, for clearness of the model structure, ! 
                     !              let's initialize again - it is no harm at all.                                    !
                     ! ---------------------------------------------------------------------------------------------- !
                     fevp1_t_rate_src(k,ids) = 0._r8
                     fevp2_t_rate_src(k,ids) = 0._r8
                     ! Note that I am using 'nonzero=1.e-20' instead of 'cmfmin=1.e-5' below, since downdraft sources generated 
                     ! in the current layer should go through vertical downward evolution and detrainment at the base interface
                     ! by the detrainment criteria. If I use 'cmfmin' here, those downdraft source with 'cmf_d_src < cmfmin'
                     ! cannot be detrained. 
                     if( cmf_d_src(k,ids) .lt. nonzero ) then ! No downdraft source in this layer
                        ix_d_src(k,ids)  = 0  
                        cmf_d_src(k,ids) = 0._r8
                     endif
                  enddo

                  ! ---------------------------------------------------------------------------------------------------------------- !
                  ! Compute total number of downdrafts (ndb_evp) that will be used for computating                                   !
                  ! evaporation of convective precipitation at the base interface of current layer.                                  !
                  ! Although mixing downdraft generated in the current layer ( cmf_d_src(k,1) ) will not evaporate convective        !
                  ! precipitation due to zero displacement diatance (i.e., the source level of mixing downdraft in the current layer !
                  ! is pt = ps0(km) ), I am including this contribution in computing 'ndb_evp >= 0'.                                 !
                  ! Note that ndb_evp = 0 is also possible, for example in the lowest model layer. However, in that case,            !
                  ! below constraint of 'goto 20 if ix_d_src(ks,ids) .eq. 0' will prevent division by zero within the main           !
                  ! computation part later. Thus, it is completely OK.                                                               !
                  ! Mar.18.2013. Do not perform evaporation of precipitation within 'constraint' downdraft which is a numerical not  !
                  !              physical downdraft. This removal of 'constraint' downdraft in treating evaporation may also help to !
                  !              reduce sensitivity to vertical resolution. Also neglect evaporation within 'top downdraft' since    !
                  !              geometrically, convective precipitation is below top downdraft.                                     !
                  !              removing these 'constraint' and 'top' downdrafts in doing evaporation of precipitation will also    !
                  !              help to reduce PREH20, which is extremely good.                                                     !
                  ! Mar.18.2013. In order to remove the contribution of mixing downdraft generated in the current layer, I changed   !
                  !              I changed  'do ks   = ktop_msfc(msfc), k, -1' to ' do ks   = ktop_msfc(msfc), k + 1, -1', which     !
                  !              should not do anything when 'ktop_msfc(msfc) = k' by construction.                                  ! 
                  ! ---------------------------------------------------------------------------------------------------------------- !
    
                ! Mar.17.2014. Since 'cmf_d_src(ks,1)' is the mass flux at the top interface of current 'k' layer, while
                !              the limitation of of evaporation within downdraft with 'evap_prep_dn' is imposed at the
                !              base interface, I multiplied 'tmp3' to the 'cmf_d_src(ks,1)' to compute downdraft mass
                !              flux at the base interface of curent layer, 'k'.
                !              Also note that I decided to exclude the mixing downdraft generated in the current layer 'k'
                !              by using 'ks = ktop_msfc(msfc), k + 1, -1' instead of 'ks = ktop_msfc(msfc), k, -1' since
                !              evaporation does not occur in the mixing downdraft generated in the 'k' layer. Note that
                !              this 'cmfdb_evp, ndb_evp' are used as ONLY THE LIMITER not for computing the actual 
                !              evaporation rate within the downdraft. Thus below my computation is completely correct.
                !              Note that we should only treat 'mixing downdraft', not 'top, constrained' downdrafts for
                !              treating evaporation of precipitation within downdraft. If we want to include 'top, constrained
                !              downdrafts, some portion of the codes (i.e., dz and dp in the below, and 'evap_prep_dn') should
                !              be consistently modified. 

                  ndb_evp = 0
                  cmfdb_evp = 0._r8
                  rho = 0.5_r8 * ( rho0top(k) + rho0bot(k) )
                  dz  = zs0(k) - zs0(km)
                  dp  = ps0(km) - ps0(k)
                  eps_dn = epsz_dn / ( rho * g )
                  del_dn = delz_dn / ( rho * g )
                  if( exp_cmf .eq. 1 ) then
                      tmp1 = eps_dn - del_dn
                      tmp2 = min( tmp1, log( 1._r8 + epsz0_max * dz ) / max( dp, nonzero ) )
                      if( abs( tmp2 - tmp1 ) .ge. nonzero ) then ! Constraint is imposed.
                          eps_dn = eps_dn * ( tmp2 / tmp1 ) 
                          del_dn = del_dn * ( tmp2 / tmp1 ) 
                      endif
                      tmp3 = exp( dp * ( eps_dn - del_dn ) )
                  elseif( exp_cmf .eq. 2 ) then
                      tmp3 = max( 0._r8, 1._r8 + dp * ( eps_dn - del_dn ) )
                  endif     
                  do ks   = ktop_msfc(msfc), k + 1, -1   ! This 'ks'  is a layer index where downdraft sources are generated.
                     if( ix_d_src(ks,1) .eq. 1 ) then
                         ndb_evp = ndb_evp + 1
                         cmfdb_evp = cmfdb_evp + cmf_d_src(ks,1) * tmp3
                     endif
                  enddo 

                  ! ---------------------------------------------------------------------------------------------------- !
                  ! Compute vertical evolution of individual downdraft from the top interface (or top of each downdraft) !
                  ! to the base interface within the given layer.                                                        !
                  ! ---------------------------------------------------------------------------------------------------- !

                  do ks = ktop_msfc(msfc), k, -1        ! This 'ks'  is a layer index where downdraft sources are generated.

                     ! ------------------------------------------------------------------------------------------- !
                     ! Convert updraft segment index at surface into shortened-updraft segment index in each layer !
                     ! ------------------------------------------------------------------------------------------- !

                     m = m_from_msfc(ks,msfc)

                     do ids = 1, 3      ! This 'ids' is the type of downdraft source ( 1 : Mixing downdraft, 2 : Top downdraft, 3 : Constraint downdraft )                   

                        ! ----------------------------------------------------------------------------------------- !
                        ! Define downdraft properties at the downdraft 'top'.                                       !
                        ! Perform downward evolution only for the surviving (existing) downdraft                    ! 
                        ! in the current layer.                                                                     !
                        ! Be careful : I should use index 'ks' not 'k' in the below block since all the source      !
                        ! variables of downdraft ( cmf_d_src(ks,ids), thl_d_src(ks,ids), ... ) are defined with the !
                        ! layer index of 'origination layer' ( 'ks' ), whose values however are continuously        !
                        ! updated in each layer as downdraft moves down into the layers below.                      !
                        ! ----------------------------------------------------------------------------------------- !
                        if( ix_d_src(ks,ids) .eq. 1 ) then
                           cmf_dt = cmf_d_src(ks,ids)
                           thl_dt = thl_d_src(ks,ids)
                           qt_dt  =  qt_d_src(ks,ids)
                           u_dt   =   u_d_src(ks,ids)
                           v_dt   =   v_d_src(ks,ids)
                           w_dt   =   w_d_src(ks,ids)
                           do mt = 1, ncnst
                              tr_dt(mt) = tr_d_src(ks,ids,mt)
                           enddo
                           ql_dt  =  ql_d_src(ks,ids)
                           qi_dt  =  qi_d_src(ks,ids)
                           ! ------------------------------------------------------------------------ !
                           ! Mar.11.2013. Add initialization of evaporation rate at the top interface ! 
                           ! ------------------------------------------------------------------------ !
                           fevp1_t_rate = fevp1_t_rate_src(ks,ids)
                           fevp2_t_rate = fevp2_t_rate_src(ks,ids)
                        else
                           cmf_d_src(ks,ids) = 0._r8
                           ix_d_src(ks,ids)  = 0
                           goto 20
                        endif
                        ! ------------------------------------------------------------------------------ !
                        ! Define mean properties at the downdraft 'top'                                  !
                        ! Here, 'top' is defined as the level where downdraft starts its downward motion !
                        ! in each layer. When 'k = ks', 'top' is defined in a different way depending on !
                        ! whether the source of downdraft is 'mixing downdraft', 'top downdraft'         !
                        ! and 'area downdraft'.                                                          !
                        ! Mar.07.2013. Define 'thv_mean_t' and 'thv_mean_b' by including updraft as well !
                        !              as environment for use in computing vertical evolution of         !
                        !              downdraft vertical velocity, simular to convective updraft.       !
                        !              Note that we are simply assuming that downdraft fractional area   !
                        !              is zero, so that only updraft information is included in          !
                        !              computing 'thv_mean_t' and 'thv_mean_b' below.                    !       
                        !              This is for developing scale-adaptive parameterization even for   !
                        !              prognostic convection scheme in future.                           !
                        !              Note we don't compute 'rho_mean_t, rho_mean_b' for simplicity.    ! 
                        ! ------------------------------------------------------------------------------ ! 
                        p_t   =     ps0(k)
                        p_b   =    ps0(km)
                        z_t   =     zs0(k)     
                        z_b   =    zs0(km)
                        thv_t = thv0top(k)
                        thv_b = thv0bot(k)
                        rho_t = rho0top(k)
                        rho_b = rho0bot(k)
                        thv_mean_t =  a_u(k) *  thva_u(k) + ( 1._r8 -  a_u(k) ) * thv0top(k)
                        thv_mean_b = a_u(km) * thva_u(km) + ( 1._r8 - a_u(km) ) * thv0bot(k)
                        if( ks .eq. k ) then
                           if( ids .eq. 1 ) then      ! Mixing Downdraft
                              p_t = ps0(km)
                              z_t = zs0(km)
                           elseif( ids .eq. 2 ) then  ! Top Downdraft                 
                              p_t = ptops(k,m)
                              z_t = ztops(k,m)
                           else                       ! Constraint Downdraft
                              p_t = ps0(k)   
                              z_t = zs0(k)
                           endif
                           p_t = min( p_t, ps0(km) )
                           z_t = max( z_t, zs0(km) )
                           thv_t = thv0bot(k) + ( p_t - ps0(km) ) * ( thv0top(k) - thv0bot(k) ) / ( ps0(k) - ps0(km) )
                           rho_t = rho0bot(k) + ( p_t - ps0(km) ) * ( rho0top(k) - rho0bot(k) ) / ( ps0(k) - ps0(km) )
                           thv_mean_t = thv_mean_b + ( p_t - ps0(km) ) * ( thv_mean_t - thv_mean_b ) / ( ps0(k) - ps0(km) )
                           cmf_ad(k,ks,m,ids)      =    cmf_dt
                           thl_ad(k,ks,m,ids)      =    thl_dt
                           qt_ad(k,ks,m,ids)       =     qt_dt 
                           u_ad(k,ks,m,ids)        =      u_dt 
                           v_ad(k,ks,m,ids)        =      v_dt 
                           w_ad(k,ks,m,ids)        =      w_dt
                           a_ad(k,ks,m,ids)        =     0._r8
                           do mt = 1, ncnst
                              tr_ad(k,ks,m,ids,mt) = tr_dt(mt)
                           enddo
                           ql_ad(k,ks,m,ids)       =     ql_dt
                           qi_ad(k,ks,m,ids)       =     qi_dt
                        endif

                        call conden( p_t, thl_dt, qt_dt, th, qv, ql, qi, qse, id_check )
                        th_dt  = th
                        qv_dt  = qv
                        ql_dt  = ql
                        qi_dt  = qi
                        thv_dt = th * ( 1._r8 + zvir * qv - ql - qi )
                        bogtop = rbuoy_dn * ( 1._r8 - thv_dt / thv_mean_t )
                        dp         = p_b - p_t
                        dz         = z_t - z_b
                        dpad(k,ks,m,ids) = dp
                        rho        = 0.5_r8 * ( rho_t + rho_b )
                        exn_b      = exns0(km)
                        if( k .gt. 1 ) then
                           ! ----------------------------------- !
                           ! Use different mu for each downdraft !
                           ! ----------------------------------- !
                           if( ids .eq. 1 ) then
                              mu = mu_mix
                           elseif( ids .eq. 2 ) then                    
                              mu = mu_top
                           elseif( ids .eq. 3 ) then
                              mu = mu_area
                           endif
                           tmp1 = ( 1._r8 - mu ) * thvl0top(km) + mu * thvl0bot(k)
                           tmp2 = ( 1._r8 - mu ) * thv0top(km)  + mu * thv0bot(k)
                           thvl_minE  = min( min( thvl0bot(k), thvl0top(k) ), tmp1 ) 
                           thv_minE   = min( min( thv0bot(k),  thv0top(k)  ), tmp2 ) 
                         ! Apr.15.2014. Modified formula for thv_minE to obtain a reasonable solution when the mean inversion exists.
                           thv_minE  = max( thv0top(km), tmp2 ) 
                         ! Apr.15.2014. Modified formula for thv_minE
                           thvl_minE  = thvl_minE + offset_minE
                           thv_minE   = thv_minE  + offset_minE
                         ! Mar.15.2014. For continuous buoyancy sorting in order to impose a stability in the code 
                         !              and in order to minimize perturbation growth.
                           thv_max    = max( thv0top(km), thv0bot(k) )
                           thv_min    = min( thv0top(km), thv0bot(k) )
                        else
                           thvl_minE  = -1.e8_r8 ! Always detrain downdraft in the lowest model layer after all diabatic forcings
                           thv_minE   = -1.e8_r8 ! Always detrain downdraft in the lowest model layer after all diabatic forcings
                         ! Mar.15.2014. Below two lines are just a place holder.
                           thv_max    = -1.e8_r8
                           thv_min    = -1.e8_r8
                        endif
                        ! ------------------------------------------------------------------------------------ !
                        ! Convert the unit of downdraft entrainment and detrainment rates from [1/z] to [1/Pa] !
                        ! Jan.30.2013. Following the treatment of convective updraft, impose similar advanced  !
                        !              constraint on convective downdraft below.                               !
                        ! Feb.06.2013. Always choose physically reasonable exp_cmf = 1 as of this day.         !
                        ! Nov.18.2013. Impose perfectly physical consistent limit. Only impose a upper limit   !
                        !              since we don't need to sorry about the decrease of downdraft mass flux  !
                        !              as downdraft moves down : we only need to worry about extremely huge    !
                        !              increase when downdraft moves down.                                     !
                        ! ------------------------------------------------------------------------------------ !
                        eps_dn = epsz_dn / ( rho * g )
                        del_dn = delz_dn / ( rho * g )
                        tmp1 = eps_dn - del_dn
                        tmp2 = min( tmp1, log( 1._r8 + epsz0_max * dz ) / max( dp, nonzero ) )
                        if( abs( tmp2 - tmp1 ) .ge. nonzero ) then ! Constraint is imposed.
                           eps_dn = eps_dn * ( tmp2 / tmp1 ) 
                           del_dn = del_dn * ( tmp2 / tmp1 ) 
                        endif

                        cmf_db  = cmf_dt * exp( dp * ( eps_dn - del_dn ) )
                        ! ----------------------------------------------- !
                        ! Update 'u_db,v_db' by diabatic horizontal PGF   !
                        ! Instead of using the gradient within the layer, ! 
                        ! I should use the gradient between the layers.   !
                        ! ----------------------------------------------- !

                        ssu_tmp = ssu0(k)
                        ssv_tmp = ssv0(k)
                        u_med   = u0bot(k) + ssu_tmp * ( 0.5_r8 * ( p_t + p_b ) - p_b )
                        v_med   = v0bot(k) + ssv_tmp * ( 0.5_r8 * ( p_t + p_b ) - p_b )

                        u_grdPGF = ssu0(k)
                        v_grdPGF = ssv0(k)            
                        
                        call progup_uv( -eps_dn, PGFc_dn, p_t, p_b, u_med, ssu_tmp, u_grdPGF, u_dt, u_db )
                        call progup_uv( -eps_dn, PGFc_dn, p_t, p_b, v_med, ssv_tmp, v_grdPGF, v_dt, v_db )

                        call progup_uv( -eps_dn,   0._r8, p_t, p_b, u_med, ssu_tmp, u_grdPGF, u_dt, u_db_adi )
                        call progup_uv( -eps_dn,   0._r8, p_t, p_b, v_med, ssv_tmp, v_grdPGF, v_dt, v_db_adi )

                        PGF_u = u_db - u_db_adi 
                        PGF_v = v_db - v_db_adi 

                        ! -------------------------------------------------------------------------- !
                        !                                                                            !
                        ! Update 'thl_db,qt_db' by diabatic evaporation of precipitation at the base !
                        !                                                                            ! 
                        ! -------------------------------------------------------------------------- !  
                        ! ----------------------------------------------------------------------- !
                        ! 1. Evaporation of Precipitation                                         !
                        !    Also compute effective diabatic forcing before evaporation of precip !
                        !    The effective forcing for tracers should be refined later by         !
                        !    considering 'evaporation-melting' especially for droplet numbers.    !  
                        ! ----------------------------------------------------------------------- !  

                        ssthl_tmp = ssthl0(k)
                        ssqt_tmp  =  ssqt0(k)
                        ssql_tmp  =  ssql0(k)
                        ssqi_tmp  =  ssqi0(k)
                        do mt = 1, ncnst
                           sstr_tmp(mt) = sstr0(k,mt)
                        enddo
                        thl_med   = thl0bot(k) + ssthl_tmp * ( 0.5_r8 * ( p_t + p_b ) - p_b )
                        qt_med    = qt0bot(k)  +  ssqt_tmp * ( 0.5_r8 * ( p_t + p_b ) - p_b )
                        ql_med    = ql0bot(k)  +  ssql_tmp * ( 0.5_r8 * ( p_t + p_b ) - p_b )
                        qi_med    = qi0bot(k)  +  ssqi_tmp * ( 0.5_r8 * ( p_t + p_b ) - p_b )
                        do mt = 1, ncnst
                           tr_med(mt)  = tr0bot(k,mt) + sstr_tmp(mt) * ( 0.5_r8 * ( p_t + p_b ) - p_b )
                        enddo

                        call progup_thlqt( -eps_dn, 0._r8, 0._r8, p_t, p_b, thl_med, ssthl_tmp, thl_dt, thl_db    )
                        call progup_thlqt( -eps_dn, 0._r8, 0._r8, p_t, p_b, qt_med,  ssqt_tmp,  qt_dt,  qt_db     )
                        call progup_thlqt( -eps_dn, 0._r8, 0._r8, p_t, p_b, ql_med,  ssql_tmp,  ql_dt,  ql_db_adi )
                        call progup_thlqt( -eps_dn, 0._r8, 0._r8, p_t, p_b, qi_med,  ssqi_tmp,  qi_dt,  qi_db_adi )
                        do mt = 1, ncnst
                           call progup_thlqt( -eps_dn, 0._r8, 0._r8, p_t, p_b, tr_med(mt), sstr_tmp(mt), tr_dt(mt), tr_db(mt) )
                        enddo

                        tr_db(ixnumliq) = ql_db_adi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                        tr_db(ixnumice) = qi_db_adi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 

                        call conden( p_b, thl_db, qt_db, th, qv, ql, qi, qse, id_check )
                        th_db = th 
                        qs_db = qse
                        qv_db = qv
                        ql_db = ql
                        qi_db = qi
                        eff_ql = ql_db - ql_db_adi
                        eff_qi = qi_db - qi_db_adi 

                        ! ------------------------------------------------------------------------------------------------ !
                        ! TRACERS REFINEMENT NECESSARY : ADIABATIC CONDENSATION-EVAPORATION-MELTING DURING VERTICAL MOTION !
                        ! ------------------------------------------------------------------------------------------------ !
                        do mt = 1, ncnst
                           if( mt .eq. ixnumliq ) then
                              eff_tr(mt) = eff_ql * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                           elseif( mt .eq. ixnumice ) then 
                              eff_tr(mt) = eff_qi * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 
                           else
                              eff_tr(mt) = 0._r8
                           endif
                           eff_tr(mt) = max( eff_tr(mt), qmin(mt) - tr_db(mt) )   
                        enddo

                        ! ------------------------------------------------------------------------------------------------ !
                        ! TRACERS REFINEMENT NECESSARY : ADIABATIC CONDENSATION-EVAPORATION-MELTING DURING VERTICAL MOTION !
                        ! ------------------------------------------------------------------------------------------------ !

                        do mt = 1, ncnst
                           tr_db(mt) = tr_db(mt) + eff_tr(mt)               
                        enddo

                        ! -------------------------------------------------------------------------------------------------------------------- !
                        ! Apr.17.2012.                                                                                                         !
                        ! Computation of (0) discrete numerical computation at the base interface or                                           !
                        !                (1) continuous analytical computation from the top to the base of convective downdraft in each layer, !
                        ! of evaporation of convective precipitation within convective downdraft.                                              !
                        ! Note that the entire goal of below if block is to compute separate (numerical or analytical) 'evp_qtl' and 'evp_qti' !
                        ! which are the resulting change of qt at the base interface [kg/kg] by the evaporation of convective rain and snow    !
                        ! within convective downdraft.                                                                                         !
                        ! -------------------------------------------------------------------------------------------------------------------- !

                        ! ------------------------------------------------------------------------------------------ !
                        ! Mar.15.2013. All materials computing 'evp_qtl, evp_qtr' are incorporated into              !
                        !              the below subroutines of 'evap_prep_dn'.                                      !
                        ! Mar.18.2013. Perform evaporation only for 'mixing' downdraft.                              !
                        ! Apr.09.2013. Add a switch to do an analytical integration of evaporation within downdraft  !
                        !              with 'iprpback = -1'.                                                         !  
                        ! ------------------------------------------------------------------------------------------ !

                        if( ids .eq. 1 ) then

                           if( ievp_prep .ne. -5 ) then
                
                              call evap_prep_dn( z_b, z_t, p_b, p_t, w_dt, bogtop,                                             &
                                       th_db, qv_db, ql_db, qi_db, tr_db(1:ncnst), qmin(1:ncnst),                              &
                                       fevp1_t_rate, fevp2_t_rate, ievp_prep,                                                  &
                                       flxrain_bot_upeesm, flxsnow_bot_upeesm, flxtrrs_bot_upeesm(1:ncnst), a_p_msfc(km,msfc), &
                                       ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice, ndb_evp, cmfdb_evp, i, k, ks, lchnk,     &
                                       rho, thv_mean_b, cmf_db, eps_dn, del_dn,                                                &
                                       kevp_rain_dn, kevp_snow_dn, eta2, rbuoy_dn, rdrag, rjet, nonzero, wdmin,                &
                                       evp_qtl, evp_qti, evp_tr(1:ncnst), fevp1_b_rate, fevp2_b_rate, w_db )

                           else

                              ! ---------------------------------------------------------------------------------------------------------  !
                              ! Analytical treatment of evaporation within downdraft                                                       !
                              ! Note that in order to prevent ambiguity due to snow melting, below uses 'flxrain(snow)_bot_upeesm'         !
                              ! instead of 'flxrain(snow)_top' in computing 'f_R,f_S', although the use of 'flxrain(snow)_bot_upeesm'      !
                              ! is still possible since the limiter is imposed based on 'flxrain(snow)_bot_upeesm'.                        !
                              ! IMPORTANT : Due to the constraint of 'evp_max = max( qw_db - qv_db, 0._r8 )', evaporation does not occur   !
                              !             when 'ql_db,qi_qb > 0'. If downdraft is unsaturated at the base interface ( evp_max > 0 ),     !
                              !             evaporation occurs only until it is saturated without generating condensate.                   !
                              !             Thus, the imposed fractional relationship of 'ql / (ql + qi ) = f(T)' from 'conden' is not     !
                              !             broken by the evaporation within downdraft.  Thus, computating 'eff_ql = ql_db - ql_db_adi'    !
                              !             and 'eff_qi = qi_db - qi_db_adi' without considering evaporation within downdraft is correct.  !
                              !             In case of convective updraft when analytical integration is used as in the current code,      !
                              !             however, we should compute these 'effective CEF forcings' using the variables                  !
                              !             including 'precipitation fallout' both adiabatically and diabatically, since analytical        !
                              !             computation of precipitation fall-out will break the imposed fractional relationship           !
                              !             of 'ql / (ql + qi ) = f(T)', so that the restoration of this relationship using 'conden' and   !
                              !             associated heating-cooling process is eventually included as a part of 'CEF' forcing.          !
                              ! ---------------------------------------------------------------------------------------------------------- !

                              call qsat(th_dt*exnf(p_t), p_t, es, qs)
                              f_R       = kevp_rain_dn * sqrt( max( 0._r8, flxrain_bot_upeesm / max( a_p_msfc(km,msfc),nonzero)))
                              f_S       = kevp_snow_dn * sqrt( max( 0._r8, flxsnow_bot_upeesm / max( a_p_msfc(km,msfc),nonzero)))
                              srcg_V    = ( f_R + f_S ) / ( rho * g * max( w_dt, wdmin ) )
                              eps_dia_V = ( f_R + f_S ) / ( rho * g * max( w_dt, wdmin ) * max( qs, nonzero ) ) 
                              qv_med    = qt_med - ql_med - qi_med
                              ssqv_tmp  = ssqt_tmp - ssql_tmp - ssqi_tmp 
                              call progup_thlqt( -eps_dn, -eps_dia_V, srcg_V, p_t, p_b, qv_med, ssqv_tmp, qv_dt, qv_db_adi_evp )
                              qv_db_adi = qt_db - ql_db_adi - qi_db_adi 
                              evp_qt    = max( 0._r8, qv_db_adi_evp - qv_db_adi )
                              evp_qtl   = evp_qt * ( f_R / max( f_R + f_S, nonzero ) )
                              evp_qti   = evp_qt * ( f_S / max( f_R + f_S, nonzero ) )
                              ! ---------------------------------------------------------------------------------- !
                              ! Limiter should be imposed using the thermodynamic properties at the base interface !
                              ! with entrainment dilution but without evaporation effect.                          !
                              ! ---------------------------------------------------------------------------------- !
                            ! evp_qtl   = max( 0._r8, min( evp_qtl, eta2 * flxrain_bot_upeesm / max( nonzero, cmf_db ) / &
                            !             max( 1._r8, real(ndb_evp,r8) ) ) )
                            ! evp_qti   = max( 0._r8, min( evp_qti, eta2 * flxsnow_bot_upeesm / max( nonzero, cmf_db ) / &
                            !             max( 1._r8, real(ndb_evp,r8) ) ) )
                              evp_qtl   = max( 0._r8, min( evp_qtl, eta2 * flxrain_bot_upeesm / max( nonzero, cmfdb_evp ) ) )
                              evp_qti   = max( 0._r8, min( evp_qti, eta2 * flxsnow_bot_upeesm / max( nonzero, cmfdb_evp ) ) )
                              call findsp_single( qv_db, th_db*exnf(p_b), p_b, tw_db, qw_db, i, k, lchnk )
                              evp_max   = max( qw_db - qv_db, 0._r8 )
                              if( ( evp_qtl + evp_qti ) .gt. evp_max ) then
                                 tmp1 = evp_qtl * evp_max / ( evp_qtl + evp_qti )
                                 tmp2 = evp_qti * evp_max / ( evp_qtl + evp_qti )
                                 evp_qtl = tmp1
                                 evp_qti = tmp2
                              endif
                              ! ------------------------------ !
                              ! Treatment of In-Cumulus Tracer !
                              ! ------------------------------ !
                              do mt = 1, ncnst
                                 if( mt .eq. 1 ) then
                                    evp_tr(mt) = evp_qtl + evp_qti
                                 elseif( mt .eq. ixcldliq .or. mt .eq. ixcldice .or. mt .eq. ixnumliq .or. mt .eq. ixnumice ) then
                                    evp_tr(mt) = 0._r8
                                 else
                                    evp_tr(mt) = flxtrrs_bot_upeesm(mt) * ( evp_qtl + evp_qti ) / max( ( flxrain_bot_upeesm + &
                                                 flxsnow_bot_upeesm ) , nonzero )
                                 endif
                                 evp_tr(mt) = max( evp_tr(mt), qmin(mt) - tr_db(mt) )
                              enddo
                              ! ----------------------------------------------------------------------- !
                              ! Compute 'w_db' at the base interface                                    !   
                              ! In the below block, I simply neglect the reduction of 'qrain_b,qsnow_b' !
                              ! due to the evaporation of precipitation within downdraft.               !
                              ! ----------------------------------------------------------------------- !
                              t_tmp     = th_db*exnf(p_b) - ( xlv / cp ) * evp_qtl - ( xls / cp ) * evp_qti
                              qv_tmp    = qv_db + evp_qtl + evp_qti
                              th_tmp    = t_tmp / exnf(p_b)
                              thv_tmp   = th_tmp * ( 1._r8 + zvir * qv_tmp - ql_db - qi_db )
                              bogbot    = rbuoy_dn * ( 1._r8 - thv_tmp / thv_mean_b )
                              call progup_wu2( -( rdrag*eps_dn - rjet*del_dn ), rho, p_t, p_b, -bogtop, -bogbot, w_dt**2._r8, &
                                               0._r8, wd2 )
                              w_db      = max( wdmin, sqrt( max( wd2, nonzero ) ) )
                              ! ------------------------------------------ !
                              ! Define the values of other null variables. !
                              ! ------------------------------------------ !
                              fevp1_b_rate    = 0._r8 
                              fevp2_b_rate    = 0._r8

                           endif

                        else

                           evp_qtl         = 0._r8 
                           evp_qti         = 0._r8
                           evp_tr(1:ncnst) = 0._r8 
                           fevp1_b_rate    = 0._r8 
                           fevp2_b_rate    = 0._r8

                        endif



                        evp_thll       = - ( xlv / cp / exn0(k) ) * evp_qtl
                        evp_thli       = - ( xls / cp / exn0(k) ) * evp_qti
                        thl_db         =   thl_db + evp_thll + evp_thli
                        qt_db          =   qt_db  +  evp_qtl +  evp_qti
                        do mt = 1, ncnst
                           tr_db(mt)   =   tr_db(mt) + evp_tr(mt)
                        enddo

                        ! ----------------------------------------------------- !
                        ! 2. Precipitation Fallout : Neglected in the downdraft !
                        ! ----------------------------------------------------- !
                        prep_qtl  = 0._r8
                        prep_qti  = 0._r8
                        prep_thll = - ( xlv / cp / exn0(k) ) * prep_qtl
                        prep_thli = - ( xls / cp / exn0(k) ) * prep_qti

                        ! -------------------------------------------------------------------------- !
                        ! TRACERS REFINEMENT NECESSARY : PRECIPITATION FALLOUT AT THE BASE INTERFACE !
                        ! -------------------------------------------------------------------------- !

                        do mt = 1, ncnst
                           if( mt .eq. 1 ) then
                              prep_tr(mt) = 0._r8
                           elseif( mt .eq. ixcldliq ) then
                              prep_tr(mt) = prep_qtl
                           elseif( mt .eq. ixcldice ) then 
                              prep_tr(mt) = prep_qti
                           elseif( mt .eq. ixnumliq ) then
                              prep_tr(mt) = prep_qtl * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_liq**3 * density_liq ) 
                           elseif( mt .eq. ixnumice ) then 
                              prep_tr(mt) = prep_qti * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_ice**3 * density_ice ) 
                           else
                              prep_tr(mt) = tr_db(mt) * ( ( prep_qtl + prep_qti ) / max( qt_db, nonzero ) )
                           endif
                           prep_tr(mt) = max( prep_tr(mt), qmin(mt) - tr_db(mt) )   
                        enddo

                        ! -------------------------------------------------------------------------- !
                        ! TRACERS REFINEMENT NECESSARY : PRECIPITATION FALLOUT AT THE BASE INTERFACE !
                        ! -------------------------------------------------------------------------- !

                        thl_db    = thl_db + prep_thll + prep_thli
                        qt_db     =  qt_db +  prep_qtl +  prep_qti
                        do mt = 1, ncnst
                           tr_db(mt) = tr_db(mt) + prep_tr(mt)               
                        enddo
                        
                        call conden( p_b, thl_db, qt_db, th_db, qv_db, ql_db, qi_db, qs_db, id_check )
                        thv_db  =  th_db * ( 1._r8 + zvir * qv_db - ql_db - qi_db )
                        if( ids .ne. 1 ) then
                           bogbot  = rbuoy_dn * ( 1._r8 - thv_db / thv_mean_b )
                           call progup_wu2( -( rdrag*eps_dn - rjet*del_dn ), rho, p_t, p_b, -bogtop, -bogbot, w_dt**2._r8, &
                                            0._r8, wd2 )
                           w_db   = max( wdmin, sqrt( max( wd2, nonzero ) ) )             
                        endif

                        ! ---------------------------------------------------------------------- !
                        ! Dec.18.2012. End of 'w_db' iteration loop for computing evaporation of ! 
                        !              convective precipitation within convective downdraft.     !
                        ! ---------------------------------------------------------------------- ! 

                        ! ----------------------------------------------------------------------------- !
                        ! Nov.28.2012. Impose final consistency between droplet mass and droplet number !
                        !              before moving to the next layer below.                           !
                        ! Mar.11.2013. Comment-out below two lines computing 'ql_db,qi_db' since they   !
                        !              have already been computed within the above iteration loop.      !  
                        ! ----------------------------------------------------------------------------- !
                        if( ql_db .eq. 0._r8 ) tr_db(ixnumliq) = 0._r8
                        if( qi_db .eq. 0._r8 ) tr_db(ixnumice) = 0._r8

                        ! -------------------------------------------------------------------------- !
                        ! Allocate diabatic forcings at the base interface into the array variable.  !
                        ! Diabatic forcing is always computed using the mean mass flux within        !
                        ! each layer, and with non-zero mass flux at the base interface of           !
                        ! detrainment layer.                                                         !
                        ! (3) 0.5_r8 * ( cmf_dt + cmf_db ) where cmf_db > 0.                         !
                        ! -------------------------------------------------------------------------- !
                        cmf_ad_dia(k,ks,m,ids)       = 0.5_r8 * ( cmf_db + cmf_dt ) * ( dp / dp0(k) )
                        evp_thll_ad(k,ks,m,ids)      =     evp_thll            
                        evp_qtl_ad(k,ks,m,ids)       =      evp_qtl
                        evp_thli_ad(k,ks,m,ids)      =     evp_thli            
                        evp_qti_ad(k,ks,m,ids)       =      evp_qti
                        prep_thll_ad(k,ks,m,ids)     =    prep_thll
                        prep_qtl_ad(k,ks,m,ids)      =     prep_qtl
                        prep_thli_ad(k,ks,m,ids)     =    prep_thli
                        prep_qti_ad(k,ks,m,ids)      =     prep_qti
                        eff_ql_ad(k,ks,m,ids)        =       eff_ql
                        eff_qi_ad(k,ks,m,ids)        =       eff_qi
                        PGF_u_ad(k,ks,m,ids)         =        PGF_u
                        PGF_v_ad(k,ks,m,ids)         =        PGF_v
                        do mt = 1, ncnst
                           evp_tr_ad(k,ks,m,ids,mt)  =   evp_tr(mt)
                           prep_tr_ad(k,ks,m,ids,mt) =  prep_tr(mt)
                           ! wdep_tr_ad(k,ks,m,ids,mt) =  wdep_tr(mt)
                           eff_tr_ad(k,ks,m,ids,mt)  =   eff_tr(mt)
                        enddo

                        ! ----------------------------------------------------------------------- !
                        ! Compute 'net evaporation of precipitation' within each layer [kg/s/m^2] !
                        ! and update precipitation fluxes at all interfaces below.                ! 
                        ! ----------------------------------------------------------------------- !

                        evplflux = ( evp_qtl + prep_qtl ) * cmf_ad_dia(k,ks,m,ids)
                        evpiflux = ( evp_qti + prep_qti ) * cmf_ad_dia(k,ks,m,ids)

                        do mt = 1, ncnst
                           if( mt .eq. ixcldliq ) then
                              evptrflux(mt) = evplflux
                           elseif( mt .eq. ixcldice ) then
                              evptrflux(mt) = evpiflux
                           elseif( mt .eq. ixnumliq ) then
                              evptrflux(mt) = evplflux * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_rain**3 * density_rain )         
                           elseif( mt .eq. ixnumice ) then
                              evptrflux(mt) = evpiflux * 3._r8 / ( 4._r8 * 3.141592_r8 * droprad_snow**3 * density_snow )         
                           else 
                              evptrflux(mt) = ( evp_tr(mt) + prep_tr(mt) ) * cmf_ad_dia(k,ks,m,ids)
                              ! evptrflux(mt) = ( evp_tr(mt) + prep_tr(mt) + wdep_tr(mt) ) * cmf_ad_dia(k,ks,m,ids)
                           endif
                        enddo

                        flxrain_bot = flxrain_bot - evplflux
                        flxsnow_bot = flxsnow_bot - evpiflux 
                        do mt = 1, ncnst
                           flxtrrs_bot(mt) = flxtrrs_bot(mt) - evptrflux(mt)
                        enddo

                        ! ------------------------------------------------------------- !
                        ! Check whether this air will move down into the below layer.   !
                        ! Also, compute the properties of detrained airs for treating   !
                        ! organization.                                                 ! 
                        ! Detrainment of downdraft always occurs at the base interface. !
                        ! ------------------------------------------------------------- !  

                        ! ---------------------------------------------------------------------------------- !
                        ! Compute laterally-detrained airs from downdraft while it moves from 'p_t' to 'p_b' !
                        ! assuming that thermodynamic properties of detrained airs is the properties at the  !
                        ! p_t similar to the treatment of detrainment from convective updraft.               !
                        ! ---------------------------------------------------------------------------------- !  

                        if( dp * abs( eps_dn - del_dn ) .ge. 1.e-3_r8 ) then
                           fmixd = ( exp( dp * ( eps_dn - del_dn ) ) - 1._r8 ) / ( dp * ( eps_dn - del_dn ) )
                        else
                           fmixd = 1._r8
                        endif

                        cmf_ar(k,ks,m,ids) = fmixd * del_dn * dp * cmf_dt
                        thl_ar(k,ks,m,ids) = 0.5_r8 * ( thl_dt + thl_db )
                        qt_ar(k,ks,m,ids)  = 0.5_r8 * (  qt_dt +  qt_db )
                        u_ar(k,ks,m,ids)   = 0.5_r8 * (   u_dt +   u_db )
                        v_ar(k,ks,m,ids)   = 0.5_r8 * (   v_dt +   v_db )
                        do mt = 1, ncnst
                           tr_ar(k,ks,m,ids,mt) = 0.5_r8 * ( tr_dt(mt) + tr_db(mt) )
                        enddo
                        ql_ar(k,ks,m,ids)  = 0.5_r8 * (  ql_dt +  ql_db )
                        qi_ar(k,ks,m,ids)  = 0.5_r8 * (  qi_dt +  qi_db )

                        ! --------------------------------------------------------------- !
                        ! Downdraft buoyancy sorting                                      !
                        ! Check whether downdraft will be detrained at the base interface !
                        ! --------------------------------------------------------------- !

                        fddet = 0._r8
                        if( dbsort_con ) then
                            fddet = ( thv_db - thv_min ) / max( thv_max - thv_min , nonzero )
                            fddet = max( 0._r8, min( 1._r8, fddet ) )
                        else
                           tmp1 = thv_db 
                           tmp2 = thv_minE 
                           if( tmp1 .gt. tmp2 ) fddet = 1._r8
                        endif 
                        if( k .eq. 1 .or. cmf_db * ( 1._r8 - fddet ) .lt. cmfmin ) fddet = 1._r8                

                      ! if( fddet .gt. 0._r8 ) then
                            tmp3 = cmf_ar(k,ks,m,ids)
                            cmf_ar(k,ks,m,ids) = tmp3 + cmf_db * fddet
                            thl_ar(k,ks,m,ids) = ( thl_ar(k,ks,m,ids) * tmp3 + thl_db * cmf_db * fddet ) / &
                                                 max( nonzero, cmf_ar(k,ks,m,ids) )
                            qt_ar(k,ks,m,ids)  = (  qt_ar(k,ks,m,ids) * tmp3 +  qt_db * cmf_db * fddet ) / &
                                                 max( nonzero, cmf_ar(k,ks,m,ids) )
                            u_ar(k,ks,m,ids)   = (   u_ar(k,ks,m,ids) * tmp3 +   u_db * cmf_db * fddet ) / &
                                                 max( nonzero, cmf_ar(k,ks,m,ids) )
                            v_ar(k,ks,m,ids)   = (   v_ar(k,ks,m,ids) * tmp3 +   v_db * cmf_db * fddet ) / &
                                                 max( nonzero, cmf_ar(k,ks,m,ids) )
                            do mt = 1, ncnst
                               tr_ar(k,ks,m,ids,mt) = (  tr_ar(k,ks,m,ids,mt) * tmp3 +  tr_db(mt) * cmf_db * fddet ) / &
                                                      max( nonzero, cmf_ar(k,ks,m,ids) )
                            enddo
                          ! Nov.28.2012. Below two lines are added. 
                            ql_ar(k,ks,m,ids)  = (  ql_ar(k,ks,m,ids) * tmp3 +  ql_db * cmf_db * fddet ) / &
                                                 max( nonzero, cmf_ar(k,ks,m,ids) )
                            qi_ar(k,ks,m,ids)  = (  qi_ar(k,ks,m,ids) * tmp3 +  qi_db * cmf_db * fddet ) / &
                                                 max( nonzero, cmf_ar(k,ks,m,ids) )
                      ! endif

                           ! ------------------------------------------------------------------------------ !
                           ! For diagnostic purpose and organization parameterization using cmf_d(0),       ! 
                           ! I am explicitly saving 'cmf_ad(0,ks,m,ids), thl_ad(0,ks,m,ids), and etc' here. !
                           ! This does not influence tendency computation in the lowest model layer         !
                           ! by downdraft since downdraft fluxes at surface will be set to zero before      !
                           ! computing tendency as will be explained later.                                 !
                           ! Jun.29.2011. Note that this explicit computation at surface is also used for   !
                           !              computing downdraft TKE within PBL ( DKE ) and anomalies of       !
                           !              downdraft conservative scalars within PBL later for convective    !
                           !              organization. Thus, this explicit adding here is important.       ! 
                           ! Jul.13.2011. This explicit non-zero setting at surface is extremely useful in  !
                           !              computing 'tkePBLorg = wa_d(0,it=1)' for computing additional     !
                           !              background mean w associated with convective organization.        !
                           ! Mar.19.2014. Below explicit computation of downdraft properties at the surface !
                           !              is now extremely important (not just diagnostic purpose), since   !
                           !              it is critically importantly used in partitioning convective      !
                           !              tendency in the lowest model layer intp the entire cumulus layers !
                           !              or the PBL depth in association with the ipartition = 1.          !        
                           ! ------------------------------------------------------------------------------ !
                           if( k .eq. 1 ) then
                              cmf_ad(km,ks,m,ids) = cmf_db            
                              thl_ad(km,ks,m,ids) = thl_db
                              qt_ad(km,ks,m,ids)  =  qt_db
                              u_ad(km,ks,m,ids)   =   u_db
                              v_ad(km,ks,m,ids)   =   v_db
                              w_ad(km,ks,m,ids)   =   w_db
                              a_ad(km,ks,m,ids)   = cmf_db / rho_b / w_db
                              do mt = 1, ncnst 
                                 tr_ad(km,ks,m,ids,mt) = tr_db(mt)
                              enddo
                              ql_ad(km,ks,m,ids)  =  ql_db
                              qi_ad(km,ks,m,ids)  =  qi_db
                           endif
                           ! Mar.15.2014. For continuous buoyancy sorting in order to impose a stability in the code 
                           !              and in order to minimize perturbation growth.

                        if( fddet .eq. 1._r8 ) then                ! Complete detrainment 
                           cmf_db             =  0._r8
                           thl_db             =  0._r8
                           qt_db              =  0._r8
                           u_db               =  0._r8
                           v_db               =  0._r8
                           w_db               =  0._r8
                           do mt = 1, ncnst
                              tr_db(mt)       =  0._r8
                           enddo
                           ql_db              =  0._r8
                           qi_db              =  0._r8
                           ix_d_src(ks,ids)   =  0
                           cmf_d_src(ks,ids)  =  0._r8
                           goto 20
                        endif
                        ! -------------------------------------------------------------- !
                        ! Allocate downdraft values at the base interface into the array !
                        ! -------------------------------------------------------------- !  
                        cmf_ad(km,ks,m,ids) = cmf_db * ( 1._r8 - fddet )           
                        thl_ad(km,ks,m,ids) = thl_db
                        qt_ad(km,ks,m,ids)  =  qt_db
                        u_ad(km,ks,m,ids)   =   u_db
                        v_ad(km,ks,m,ids)   =   v_db
                        w_ad(km,ks,m,ids)   =   w_db
                        a_ad(km,ks,m,ids)   = cmf_db / rho_b / w_db
                        do mt = 1, ncnst 
                           tr_ad(km,ks,m,ids,mt) = tr_db(mt)
                        enddo
                        ql_ad(km,ks,m,ids)  =  ql_db
                        qi_ad(km,ks,m,ids)  =  qi_db
                        ! --------------------------------------- !
                        ! Initialization for the next computation !
                        ! --------------------------------------- !  
                        ix_d_src(ks,ids)  = 1
                        cmf_d_src(ks,ids) = cmf_db * ( 1._r8 - fddet )
                        thl_d_src(ks,ids) = thl_db
                        qt_d_src(ks,ids)  =  qt_db
                        u_d_src(ks,ids)   =   u_db
                        v_d_src(ks,ids)   =   v_db
                        w_d_src(ks,ids)   =   w_db
                        do mt = 1, ncnst
                           tr_d_src(ks,ids,mt) = tr_db(mt)
                        enddo
                        ql_d_src(ks,ids)  =  ql_db
                        qi_d_src(ks,ids)  =  qi_db
                        ! ------------------------------------------------------------------------ !
                        ! Mar.11.2013. Add initialization of evaporation rate at the top interface ! 
                        !              for the next layer below.                                   !
                        ! ------------------------------------------------------------------------ !
                        fevp1_t_rate_src(ks,ids) = fevp1_b_rate
                        fevp2_t_rate_src(ks,ids) = fevp2_b_rate
                        ! ------------------------------------------------ !
                        ! End of downdraft sorting of one downdraft source !
                        ! ------------------------------------------------ !  
                        
20                      continue 
                     enddo                     ! ids  = 1, 3.                     This 'ids'  is a type of downdraft source.
                  enddo                        ! ks   = ktop_msfc(msfc), k, -1.   This 'ks'   is a layer index of downdraft source.

                  ! -------------------------------------------------------------------------------------- !
                  ! Assign final precipitation flux at the bottom interface into the array                 !
                  ! Also, compute 'evprain_d(k,msfc), evpsnow_d(k,msfc), evptrrs_d(k,msfc,mt)' >= 0        !
                  ! by differencing two flux variables.                                                    !
                  ! For safety and imposing full consistency, impose the constraint that droplet           !
                  ! size of precipitation rain/snow are externally specified.                              !
                  ! Note that 'flxtrrs_bot(mt)' also contains the effect of wet deposition of aerosols     !
                  ! as well as evaporation ( & production ) of convective precipitation within             !
                  ! convective downdraft. Note also that wet deposition does not influences condensate but !
                  ! only affect tracers.                                                                   !
                  ! -------------------------------------------------------------------------------------- !

                  flxrain_msfc(km,msfc)  = flxrain_bot
                  flxsnow_msfc(km,msfc)  = flxsnow_bot
                  evprain_d_msfc(k,msfc) = ( flxrain_bot_upeesm - flxrain_bot ) * ( g / dp0(k) ) ! >= 0.
                  evpsnow_d_msfc(k,msfc) = ( flxsnow_bot_upeesm - flxsnow_bot ) * ( g / dp0(k) ) ! >= 0.
                  do mt = 1, ncnst
                     if( mt .eq. ixcldliq ) then
                        flxtrrs_msfc(km,msfc,mt) = flxrain_msfc(km,msfc)
                     elseif( mt .eq. ixcldice ) then
                        flxtrrs_msfc(km,msfc,mt) = flxsnow_msfc(km,msfc)
                     elseif( mt .eq. ixnumliq ) then
                        flxtrrs_msfc(km,msfc,mt) = flxrain_msfc(km,msfc) * 3._r8 / ( 4._r8 * 3.141592_r8 * &
                                                   droprad_rain**3 * density_rain )         
                     elseif( mt .eq. ixnumice ) then
                        flxtrrs_msfc(km,msfc,mt) = flxsnow_msfc(km,msfc) * 3._r8 / ( 4._r8 * 3.141592_r8 * &
                                                   droprad_snow**3 * density_snow )         
                     else 
                        flxtrrs_msfc(km,msfc,mt) = flxtrrs_bot(mt)
                     endif
                     evptrrs_d_msfc(k,msfc,mt) = ( flxtrrs_msfc(km,msfc,mt) - flxtrrs_bot_upeesm(mt) ) * ( g / dptr0(k,mt) )
                  enddo

                  ! Sanity Check : In the new code, this case must not happen. 
                  if( flxrain_msfc(km,msfc) .lt. -1.e-18_r8 .or. flxsnow_msfc(km,msfc) .lt. -1.e-18_r8 ) then 
                     write(iulog,*) 'UNICON : Stop - Negative precipitation flux after computing evaporation within environment'
                     write(iulog,*)  k, flxrain_msfc(k,msfc), flxsnow_msfc(k,msfc), flxrain_msfc(km,msfc), flxsnow_msfc(km,msfc)
                     call endrun('STOP : UNICON - Negative Precipitation Flux')
                     write(iulog,*) 
                  endif

                  ! --------------------------------------------------------------------------------------------------- !
                  ! Final Net Rain and Snow Production Tendency in Each Layer for Each Updraft Segment.                 !
                  !                                                                                                     !
                  ! A. For evaporation and snow melting, the sign of tracer tendencies have been already reversed,      !
                  !    so that simplying adding all the tracer tendencies produces correct results.                     ! 
                  !                                                                                                     ! 
                  ! B. Note that if we turn-off area and vertical velocity constraint of convective downdraft, then     !
                  !                                                                                                     ! 
                  !    1. qrten_d_msfc(k,msfc)   = - evprain_d_msfc(k,msfc),                                            !
                  !    2. qrten_d_msfc(k,msfc)   = - evprain_d_msfc(k,msfc),                                            ! 
                  !    3. trrsten_d_msfc(k,msfc) =   evptrrs_d_msfc(k,msfc)                                             !
                  !                                                                                                     !
                  !    where 3 variables on the LSH will be computed later after performing vertical velocity and       !
                  !    are constraints on the convective downdraft. In fact, this is the versy reason why these are     !
                  !    computed later instead of downdraft computation loop above.                                      !
                  !                                                                                                     !
                  ! C. With the current treatment of new unified treatment of evaporation process with                  !
                  !    the 'do msfc = 1, nseg' loop outside of 'do k = ktop, 1, -1', there is no way to perfectly treat !
                  !    the area and vertical velocity constraint of convective downdraft. Thus, I should turn-off       !
                  !    the area and vertical velocity constraint of convective downdraft, which is OK given the huge    !
                  !    benefit of physically reasonable computation of evaporation process in a consistent way.         !
                  ! D. Note also that 'trrsten_u_msfc(k,msfc,mt),evptrrs_d_msfc(k,msfc,mt)' already contains the effect !
                  ! of wet deposition of aerosols within convective updraft and downdrafts.                             !
                  !                                                                                                     !
                  ! MODIFICATION IS REQUIRED :                                                                          ! 
                  !              Below 'x_p_msfc(km,msfc), y_p_msfc(km,msfc)' should be re-computed.                    ! 
                  !                                                                                                     !
                  ! Mar.05.2013. Recompute 'x_p_msfc(km,msfc), y_p_msfc(km,msfc)' using grid-mean precipitation fluxes  !
                  !              both from the one coming from at the top interface with evaporation within environment ! 
                  !              ( flxrain_bot_ee + flxsnow_bot_ee ) and the flux generated in the current layer from   !
                  !              convective updraft ( ( qrten_u_msfc(k,msfc) + qsten_u_msfc(k,msfc) ) * ( dp0(k) / g )  !
                  !              as weighting factors. Note that since snow melting and evaporation within downdraft    !
                  !              occur both in the above two sources in the same way, we don't need to consider the     !
                  !              effect of snow melting and evaporation within downdraft in computing x_p_msfc(km,msfc) !
                  !              and y_p_msfc(km,msfc). Note also that the weighting factor is the product of 'in-prep  !
                  !              precipitation flux' and the 'area', i.e., 'area' is also used as a weighting factor,   !
                  !              which is also a perfectly reasonable choice.                                           !
                  !              In addition, this new computation saves computation time.                              !
                  !              The fundamental simplifying assumption of this approach is that both precipitation and !
                  !              precipitating updraft areas can be described as a single disk with a single (x,y,R) at !
                  !              each model interface.                                                                  !  
                  ! --------------------------------------------------------------------------------------------------- !

                  tmpw = ( flxrain_bot_ee + flxsnow_bot_ee ) + ( qrten_u_msfc(k,msfc) + qsten_u_msfc(k,msfc) ) * ( dp0(k) / g )
                  tmpw = max( tmpw, nonzero )
                  x_p_msfc(km,msfc) = (  x_p_msfc(k,msfc) * ( flxrain_bot_ee + flxsnow_bot_ee ) + &
                     x_um_msfc(k,msfc) * ( qrten_u_msfc(k,msfc) + qsten_u_msfc(k,msfc) ) * ( dp0(k) / g ) ) / tmpw
                  y_p_msfc(km,msfc) = (  y_p_msfc(k,msfc) * ( flxrain_bot_ee + flxsnow_bot_ee ) + &
                     y_um_msfc(k,msfc) * ( qrten_u_msfc(k,msfc) + qsten_u_msfc(k,msfc) ) * ( dp0(k) / g ) ) / tmpw

               enddo                        ! k    = ktop_msfc(msfc), 1, -1.        This 'k'    is a layer index where vertical evolution of downdraft is computed.
            enddo                        ! msfc = 1, nseg                        This 'msfc' is a number of updraft segment at surface.

            ! ------------------------------------------------------------------------------------------------------- !
            !                                                                                                         !                
            ! Iteration for treating accretion should end here.                                                       ! 
            !                                                                                                         !
            ! The variables used for the next accretion iteration computations within 'subroutine prod_prep_up' are   !
            !                                                                                                         !
            !  1. flxrain_msfc(k,msfc), flxsnow_msfc(k,msfc)                                                          !
            !  2. a_p_msfc(k,msfc)                                                                                    ! 
            !  3. am_us_msfc(k,msfc)                                                                                  !
            !  4. am_pu_msfc(k,msfc)                                                                                  ! 
            !                                                                                                         !
            ! all of which are already computed above. Thus, this is the perfect location of the end of accretion     !
            ! iteration loop.                                                                                         !
            ! ------------------------------------------------------------------------------------------------------- !

         enddo                        ! Enf of iacc = 1, nacc                 This 'iacc' is a number of accretion iteration.

         ! ---------------------------------------------------------------------- !
         ! Compute mass-flux weighted mean                                        !
         !  (1) downdraft properties at the interfaces                            !
         !  (2) detrained properties at the layer mid-point                       !
         ! Also compute mass-flux weighted diabatic contribution.                 !
         ! ---------------------------------------------------------------------- !

         cmf_d_org_pblh           = 0._r8
         thl_d_org_pblh           = 0._r8
         qt_d_org_pblh            = 0._r8
         u_d_org_pblh             = 0._r8
         v_d_org_pblh             = 0._r8
         tr_d_org_pblh(1:ncnst)   = 0._r8
         qt_dia_d_org             = 0._r8
         thl_dia_d_org            = 0._r8
         tr_dia_d_org(1:ncnst)    = 0._r8
       ! May.1.2014. For the budget consistent cold-pool treatment.
         cmf_d_orgU_pblh          = 0._r8
         thl_d_orgU_pblh          = 0._r8
         qt_d_orgU_pblh           = 0._r8
         u_d_orgU_pblh            = 0._r8
         v_d_orgU_pblh            = 0._r8
         tr_d_orgU_pblh(1:ncnst)  = 0._r8
         qt_dia_d_orgU            = 0._r8
         thl_dia_d_orgU           = 0._r8
         tr_dia_d_orgU(1:ncnst)   = 0._r8

         do msfc = 1, nseg                         ! This 'msfc' is updraft segment index at surface.
            do ids  = 1, 3                            ! This 'ids' is the type of downdraft source ( 1 : mixing downdraft, 2 : top downdraft, 3 : area downdraft )
               if( ids .eq. 1 ) then
                  ks_top = ktop_msfc(msfc)
                  ks_bot = 1
               elseif( ids .eq. 2 ) then
                  ks_top = ktop_msfc(msfc)
                  ks_bot = ks_top
               elseif( ids .eq. 3 ) then
                  ks_top = ktop_msfc(msfc) - 1
                  ks_bot = 1
               endif
               do ks   = ks_top, ks_bot, -1   ! This 'ks'   is a layer index where downdraft sources are generated.
                  ! ------------------------------------------------------------------------------------------- !
                  ! Convert updraft segment index at surface into shortened-updraft segment index in each layer !
                  ! ------------------------------------------------------------------------------------------- !
                  m = m_from_msfc(ks,msfc)
                  do k = ks, 1, -1            ! This 'k'   is a layer index from the source layer to the 1st (not 2nd) layer. 
                     km = k - 1               ! This 'km'  is a base interface index
                     ! ------------------------------------------------- !
                     ! Sum of downdraft properties at the base interface !
                     ! ------------------------------------------------- ! 
                     if( cmf_ad(km,ks,m,ids) .gt. nonzero ) then
                        cmf_d(km) = cmf_d(km) + cmf_ad(km,ks,m,ids)            
                        a_d(km)   =   a_d(km) +   a_ad(km,ks,m,ids)            
                        thl_d(km) = thl_d(km) + thl_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                        qt_d(km)  =  qt_d(km) +  qt_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                        u_d(km)   =   u_d(km) +   u_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                        v_d(km)   =   v_d(km) +   v_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                        w_d(km)   =   w_d(km) +   w_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                        wa_d(km)  =  wa_d(km) +   w_ad(km,ks,m,ids) *   a_ad(km,ks,m,ids)
                        ql_d(km)  =  ql_d(km) +  ql_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)  
                        qi_d(km)  =  qi_d(km) +  qi_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)  
                        qla_d(km) = qla_d(km) +  ql_ad(km,ks,m,ids) *   a_ad(km,ks,m,ids)  
                        qia_d(km) = qia_d(km) +  qi_ad(km,ks,m,ids) *   a_ad(km,ks,m,ids)  
                        do mt = 1, ncnst
                           tr_d(km,mt)  =  tr_d(km,mt) +  tr_ad(km,ks,m,ids,mt) * cmf_ad(km,ks,m,ids)
                        enddo
                        ! ---------------------------------------------------------------------------------------- !
                        ! Refinement of Organization                                                               !
                        ! In order to compute convective organization at surface, only sum the downdraft mass flux !
                        ! that (1) was originated above PBLH, and (2) reaches down to surface by the rigorous      !
                        ! downdraft buoyancy sorting in the lowest model layer ( k = 1 ).                          !
                        ! This is an attempt to sort out only the downdraft mass flux forced by evaporation of     !
                        ! convective precipitation.                                                                !
                        ! Mar.11.2013. Comment-out below block since it is not used any more.                      !
                        !              In addition, 'th,qv,ql,qi' are not computed correctly since above 'conden'  !
                        !              was commented out.                                                          !  
                        ! ---------------------------------------------------------------------------------------- !

                        if( ids .eq. 1               .and. &  
                           ks   .gt. kpblh           .and. &
                           8.64e7_r8*(flxrain_msfc(kpblhm,msfc)/1000._r8+flxsnow_msfc(kpblhm,msfc)/250._r8) .gt. prepminPBLH_org &
                           .and. cmf_ad(1,ks,m,ids) .gt. nonzero ) then 
                           ! ---------------------------------------------------------- !
                           ! Compute adiabatic downdraft fluxes of conservative scalars !
                           ! ---------------------------------------------------------- !
                           if( km .eq. kpblhm ) then
                              cmf_d_org_pblh        =     cmf_d_org_pblh +                          cmf_ad(km,ks,m,ids)
                              thl_d_org_pblh        =     thl_d_org_pblh +    thl_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                              qt_d_org_pblh         =      qt_d_org_pblh +     qt_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                              u_d_org_pblh          =       u_d_org_pblh +      u_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                              v_d_org_pblh          =       v_d_org_pblh +      v_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                              do mt = 1, ncnst
                                 tr_d_org_pblh(mt)  =  tr_d_org_pblh(mt) +  tr_ad(km,ks,m,ids,mt) * cmf_ad(km,ks,m,ids)
                              enddo
                           endif
                           ! ----------------------------------------------------------------------------------------------------------- !
                           ! Compute diabatic forcings integrated all over the layers within PBL.                                        !
                           ! Note that there is no diabatic forcing for 'u,v' since 'PGF' forcing is a simple conversion                 !
                           ! between environment and convective downdraft.                                                               !
                           ! For tracers, tendencies in each layer is computed by using 'pdelx = dpdry0(k)' not 'dp0(k)'.                !
                           ! I multiplied 'dp0(k)' for vertical integration. However, I should check this in future.                     !
                           ! Note that below computation of diabatic forcing within downdraft does not contain either 'snow melting'     !
                           ! or 'corrective flux', both of which will be treated later as a part of environmental diabatic forcing.      !
                           ! Thus, there will be no missing process or double counting.                                                  !
                           ! Sep.10.2011. Note tha below computation of diabatic forcing within downdraft within PBL should nbe done in  !
                           !              a fully cocnsitently way as the above computation of adiabatic forcing, i.e., using the        !
                           !              exactly same set of convective downdraft as in the current code.                               !
                           ! Aug.02.2012. From above, I only allowed evaporation of precipitation within mixing downdraft not within     !
                           !              top downdraft and area downdraft anymore.                                                      !
                           ! ----------------------------------------------------------------------------------------------------------- !
                           if( k .le. kpblhm ) then
                              qt_dia_d_org        =       qt_dia_d_org     + & 
                                 g * ( prep_qtl_ad(k,ks,m,ids)  + prep_qti_ad(k,ks,m,ids)  + evp_qtl_ad(k,ks,m,ids)  &
                                 + evp_qti_ad(k,ks,m,ids)  ) * cmf_ad_dia(k,ks,m,ids)
                              thl_dia_d_org       =       thl_dia_d_org    + & 
                                 g * ( prep_thll_ad(k,ks,m,ids) + prep_thli_ad(k,ks,m,ids) + evp_thll_ad(k,ks,m,ids) &
                                 + evp_thli_ad(k,ks,m,ids) ) * cmf_ad_dia(k,ks,m,ids)
                              do mt = 1, ncnst
                                 tmp1 = dp0(k) / dptr0(k,mt)
                                 tr_dia_d_org(mt)   =     tr_dia_d_org(mt) + tmp1 * &
                                    g * ( evp_tr_ad(k,ks,m,ids,mt) + prep_tr_ad(k,ks,m,ids,mt) + wdep_tr_ad(k,ks,m,ids,mt) &
                                    + eff_tr_ad(k,ks,m,ids,mt) ) * cmf_ad_dia(k,ks,m,ids)
                              enddo
                           endif
                        endif

                      ! -------------------------------------------------------------------------------------------------------------------- !
                      ! May.1.2014. For budget consistent cold pool treatment.                                                               ! 
                      ! In contrast to the above part computing 'cmf_d_org_pblh' that exculsively sinks down into 'awk_PBL',                 !
                      ! this part computes 'cmf_d_orgU_pblh' that exclusively sinks down into '1-awk_PBL'.                                   !
                      ! Note that if 'i_budget_coldpool .eq. 1', there is no downdraft exclusively sinking into '1-awk_PBL', so that         !
                      ! it becomes 'cmf_d_orgU_pblh=0, thl_d_org_pblh=0, qt_dia_d_org=0' which is already done at the initialization at the  !
                      ! beginning of downdraft computation routine. Thus, I need to compute only for 'i_budget_cold_pool .eq. 2' when        !
                      ! all the downdrafts other than 'cmf_d_org_pblh' computed above exclusively fall into 'awk_PBL'.                       !
                      ! Below 'if' constraints is exact opposite to the above 'if' constraint, currently, but it can be further              !
                      ! generalized by similarly defining 'clf_d_orgG_pblh' in future.                                                       !
                      ! -------------------------------------------------------------------------------------------------------------------- !

                        if( i_budget_coldpool .eq. 2 .or. i_budget_coldpool .eq. 4 .or. i_budget_coldpool .eq. 5 ) then
               
                        if( .not. ( ids                                                          .eq. 1               .and. &  
                                    ks                                                           .gt. kpblh           .and. &
!?                                  8.64e7_r8*(flxrain_ava_msfc(kpblhm,msfc)/1000._r8+flxsnow_ava_msfc(kpblhm,msfc)/250._r8)  .gt. prepminPBLH_org .and. &
!!                                  8.64e7_r8*(flxrain_msfc(kpblhm)/1000._r8+flxsnow_msfc(kpblhm)/250._r8)                    .gt. prepminPBLH_org .and. &
                            8.64e7_r8*(flxrain_msfc(kpblhm,msfc)/1000._r8+flxsnow_msfc(kpblhm,msfc)/250._r8) &
                                                                                                 .gt. prepminPBLH_org .and. &
                            cmf_ad(1,ks,m,ids)                                                   .gt. nonzero ) ) then 
                            ! ---------------------------------------------------------- !
                            ! Compute adiabatic downdraft fluxes of conservative scalars !
                            ! ---------------------------------------------------------- !
                            if( km .eq. kpblhm ) then
                                cmf_d_orgU_pblh        =     cmf_d_orgU_pblh +                          cmf_ad(km,ks,m,ids)
                                thl_d_orgU_pblh        =     thl_d_orgU_pblh +    thl_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                                qt_d_orgU_pblh         =      qt_d_orgU_pblh +     qt_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                                u_d_orgU_pblh          =       u_d_orgU_pblh +      u_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                                v_d_orgU_pblh          =       v_d_orgU_pblh +      v_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids)
                                do mt = 1, ncnst
                                   tr_d_orgU_pblh(mt)  =  tr_d_orgU_pblh(mt) +  tr_ad(km,ks,m,ids,mt) * cmf_ad(km,ks,m,ids)
                                enddo 
                            endif
                            if( k .le. kpblhm ) then
                                qt_dia_d_orgU          =    qt_dia_d_orgU    + & 
                                                            g * ( prep_qtl_ad(k,ks,m,ids)  + prep_qti_ad(k,ks,m,ids)  + &
                                                            evp_qtl_ad(k,ks,m,ids)  + evp_qti_ad(k,ks,m,ids)  ) * & 
                                                            cmf_ad_dia(k,ks,m,ids)
                                thl_dia_d_orgU         =    thl_dia_d_orgU   + & 
                                                            g * ( prep_thll_ad(k,ks,m,ids) + prep_thli_ad(k,ks,m,ids) + &
                                                            evp_thll_ad(k,ks,m,ids) + evp_thli_ad(k,ks,m,ids) ) * & 
                                                            cmf_ad_dia(k,ks,m,ids)
                                do mt = 1, ncnst
                                   tmp1 = dp0(k) / dptr0(k,mt)
                                   tr_dia_d_orgU(mt)   =    tr_dia_d_orgU(mt) + tmp1 * &
                                                            g * ( evp_tr_ad(k,ks,m,ids,mt) + prep_tr_ad(k,ks,m,ids,mt) + &
                                                            wdep_tr_ad(k,ks,m,ids,mt) + eff_tr_ad(k,ks,m,ids,mt) ) * & 
                                                            cmf_ad_dia(k,ks,m,ids)
                                enddo 
                            endif
                        endif

                        endif

                     endif
                     ! --------------------------------------------- !
                     ! Sum of diabatic forcing at the base interface !
                     ! --------------------------------------------- !
                     if( cmf_ad_dia(k,ks,m,ids) .gt. nonzero ) then
                        cmf_d_dia(k)       =   cmf_d_dia(k)   +                             cmf_ad_dia(k,ks,m,ids)            
                        evp_thll_d(k)      =  evp_thll_d(k)   +   evp_thll_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        evp_qtl_d(k)       =   evp_qtl_d(k)   +    evp_qtl_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        evp_thli_d(k)      =  evp_thli_d(k)   +   evp_thli_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        evp_qti_d(k)       =   evp_qti_d(k)   +    evp_qti_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        prep_thll_d(k)     = prep_thll_d(k)   +  prep_thll_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        prep_qtl_d(k)      =  prep_qtl_d(k)   +   prep_qtl_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        prep_thli_d(k)     = prep_thli_d(k)   +  prep_thli_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        prep_qti_d(k)      =  prep_qti_d(k)   +   prep_qti_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        eff_ql_d(k)        =    eff_ql_d(k)   +     eff_ql_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        eff_qi_d(k)        =    eff_qi_d(k)   +     eff_qi_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        PGF_u_d(k)         =     PGF_u_d(k)   +      PGF_u_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        PGF_v_d(k)         =     PGF_v_d(k)   +      PGF_v_ad(k,ks,m,ids) * cmf_ad_dia(k,ks,m,ids)
                        do mt = 1, ncnst
                           evp_tr_d(k,mt)  =  evp_tr_d(k,mt)  +  evp_tr_ad(k,ks,m,ids,mt) * cmf_ad_dia(k,ks,m,ids)
                           prep_tr_d(k,mt) = prep_tr_d(k,mt)  + prep_tr_ad(k,ks,m,ids,mt) * cmf_ad_dia(k,ks,m,ids)
                           ! wdep_tr_d(k,mt) = wdep_tr_d(k,mt)  + wdep_tr_ad(k,ks,m,ids,mt) * cmf_ad_dia(k,ks,m,ids)
                           eff_tr_d(k,mt)  =  eff_tr_d(k,mt)  +  eff_tr_ad(k,ks,m,ids,mt) * cmf_ad_dia(k,ks,m,ids)
                        enddo
                     endif
                     ! ----------------------------------------------------------------- !
                     ! Sum of detrained properties from downdraft at the layer mid-point !
                     ! ----------------------------------------------------------------- !
                     if( cmf_ar(k,ks,m,ids) .gt. nonzero ) then
                        cmf_rd(k)      =   cmf_rd(k) +                        cmf_ar(k,ks,m,ids)            
                        thl_rd(k)      =   thl_rd(k) +   thl_ar(k,ks,m,ids) * cmf_ar(k,ks,m,ids)
                        qt_rd(k)       =    qt_rd(k) +    qt_ar(k,ks,m,ids) * cmf_ar(k,ks,m,ids)
                        u_rd(k)        =     u_rd(k) +     u_ar(k,ks,m,ids) * cmf_ar(k,ks,m,ids)
                        v_rd(k)        =     v_rd(k) +     v_ar(k,ks,m,ids) * cmf_ar(k,ks,m,ids)
                        ql_rd(k)       =    ql_rd(k) +    ql_ar(k,ks,m,ids) * cmf_ar(k,ks,m,ids)
                        qi_rd(k)       =    qi_rd(k) +    qi_ar(k,ks,m,ids) * cmf_ar(k,ks,m,ids)
                        do mt = 1, ncnst
                           tr_rd(k,mt) = tr_rd(k,mt) + tr_ar(k,ks,m,ids,mt) * cmf_ar(k,ks,m,ids)
                        enddo
                     endif
                     ! ------------------------------------------------------------------------------------- !
                     ! Compute mean downdraft properties at the layer mid-point or at the interface for each ! 
                     ! original updraft segment 'msfc'. These are diagnostic quantities.                     !
                     ! ------------------------------------------------------------------------------------- ! 
                     tmp1 = 0.5_r8 * ( a_ad(km,ks,m,ids) +  a_ad(k,ks,m,ids) ) * ( dpad(k,ks,m,ids) / dp0(k) )
                     am_d_msfc(k,msfc)    =  am_d_msfc(k,msfc) + tmp1
                     qlm_d_msfc(k,msfc)   = qlm_d_msfc(k,msfc) + tmp1 * 0.5_r8 * ( ql_ad(km,ks,m,ids) + ql_ad(k,ks,m,ids) ) * &
                                            ( dpad(k,ks,m,ids) / dp0(k) )
                     qim_d_msfc(k,msfc)   = qim_d_msfc(k,msfc) + tmp1 * 0.5_r8 * ( qi_ad(km,ks,m,ids) + qi_ad(k,ks,m,ids) ) * &
                                            ( dpad(k,ks,m,ids) / dp0(k) )
                     qrten_d_msfc(k,msfc) = qrten_d_msfc(k,msfc) - ( g / dp0(k) ) * ( prep_qtl_ad(k,ks,m,ids) + &
                                            evp_qtl_ad(k,ks,m,ids) ) * cmf_ad_dia(k,ks,m,ids)   ! <= 0
                     qsten_d_msfc(k,msfc) = qsten_d_msfc(k,msfc) - ( g / dp0(k) ) * ( prep_qti_ad(k,ks,m,ids) + &
                                            evp_qti_ad(k,ks,m,ids) ) * cmf_ad_dia(k,ks,m,ids)   ! <= 0
                     do mt = 1, ncnst
                        if( mt .eq. ixcldliq ) then
                           trrsten_d_msfc(k,msfc,mt) = qrten_d_msfc(k,msfc)
                        elseif( mt .eq. ixcldice ) then
                           trrsten_d_msfc(k,msfc,mt) = qsten_d_msfc(k,msfc)
                        elseif( mt .eq. ixnumliq ) then
                           trrsten_d_msfc(k,msfc,mt) = qrten_d_msfc(k,msfc) * 3._r8 / ( 4._r8 * 3.141592_r8 * &
                                                       droprad_rain**3 * density_rain )        
                        elseif( mt .eq. ixnumice ) then
                           trrsten_d_msfc(k,msfc,mt) = qsten_d_msfc(k,msfc) * 3._r8 / ( 4._r8 * 3.141592_r8 * &
                                                       droprad_snow**3 * density_snow )         
                        else
                           trrsten_d_msfc(k,msfc,mt) = trrsten_d_msfc(k,msfc,mt) - ( g / dptr0(k,mt) ) * &
                                                       ( prep_tr_ad(k,ks,m,ids,mt) + evp_tr_ad(k,ks,m,ids,mt) + &
                                                       wdep_tr_ad(k,ks,m,ids,mt) ) * cmf_ad_dia(k,ks,m,ids)   ! <= 0
                        endif
                     enddo
                     if( cmf_ad(km,ks,m,ids) .gt. nonzero ) then
                        cmf_d_msfc(km,msfc) = cmf_d_msfc(km,msfc) + cmf_ad(km,ks,m,ids)            
                        a_d_msfc(km,msfc)   =   a_d_msfc(km,msfc) +   a_ad(km,ks,m,ids)            
                        thl_d_msfc(km,msfc) = thl_d_msfc(km,msfc) + thl_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        qt_d_msfc(km,msfc)  =  qt_d_msfc(km,msfc) +  qt_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        u_d_msfc(km,msfc)   =   u_d_msfc(km,msfc) +   u_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        v_d_msfc(km,msfc)   =   v_d_msfc(km,msfc) +   v_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        w_d_msfc(km,msfc)   =   w_d_msfc(km,msfc) +   w_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        ql_d_msfc(km,msfc)  =  ql_d_msfc(km,msfc) +  ql_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        qi_d_msfc(km,msfc)  =  qi_d_msfc(km,msfc) +  qi_ad(km,ks,m,ids) * cmf_ad(km,ks,m,ids) 
                        wa_d_msfc(km,msfc)  =  wa_d_msfc(km,msfc) +   w_ad(km,ks,m,ids) *   a_ad(km,ks,m,ids) 
                        qla_d_msfc(km,msfc) = qla_d_msfc(km,msfc) +  ql_ad(km,ks,m,ids) *   a_ad(km,ks,m,ids) 
                        qia_d_msfc(km,msfc) = qia_d_msfc(km,msfc) +  qi_ad(km,ks,m,ids) *   a_ad(km,ks,m,ids) 
                        do mt = 1, ncnst
                           tr_d_msfc(km,msfc,mt)  =  tr_d_msfc(km,msfc,mt) +  tr_ad(km,ks,m,ids,mt) * cmf_ad(km,ks,m,ids)
                        enddo
                     endif
                     ! ---------------------------------------------------------------------------- !
                     !                                                                              !  
                     ! ---------------------------------------------------------------------------- !
                  enddo                     ! k    = ks, 1, -1.          This 'k'    is a layer index.
               enddo                        ! ks   = ks_top, ks_bot, -1. This 'ks'   is a layer index of downdraft source.
            enddo                        ! ids  = 1, 3.               This 'ids'  is a type of downdraft source.
         enddo                        ! msfc = 1, nseg             This 'msfc' is a number of updraft segment at surface.

         ! ------------------------------------------------------------------------------------- !
         ! Jul.10.2011. In order to do in-downdraft perturbation, I should divide by 'tmp3'      !
         !              However, TKE is accumulated quantities.                                  !
         !              Note that final 'thlPBLorg,qtPBLorg' are mass-flux weighted quantities.  ! 
         ! Aug.03.2011. Below may be modified in future such that 'if cmfPBLorg .le. nonzero',   !
         !              I should set cmfPBLorg = 0.                                              !
         ! ------------------------------------------------------------------------------------- ! 
      
         ! ----------------------------------------------------------------------------------------------------------------- !
         ! Aug.30.2011. Mean downdraft properties at the PBL top interface ( kpblhm ) for organized convective downdraft.    !
         ! Sep.06.2011. Also compute grid-mean vertically-averaged diabatic forcing within convective downdraft within PBL.  !
         ! Sep.09.2011. Below block is commented-out since all the wake-related adiabatic and diabatic forcings are treated  !
         !              later in a collective way.                                                                           !
         ! Sep.11.2011. I restored below block for selective chooce of convective downdrafts.
         ! ----------------------------------------------------------------------------------------------------------------- ! 

         if( cmf_d_org_pblh .gt. nonzero ) then
            thl_d_org_pblh          =     thl_d_org_pblh / cmf_d_org_pblh
            qt_d_org_pblh           =      qt_d_org_pblh / cmf_d_org_pblh 
            u_d_org_pblh            =       u_d_org_pblh / cmf_d_org_pblh 
            v_d_org_pblh            =       v_d_org_pblh / cmf_d_org_pblh 
            do mt = 1, ncnst
               tr_d_org_pblh(mt)    =  tr_d_org_pblh(mt) / cmf_d_org_pblh
            enddo
         else
            cmf_d_org_pblh          = 0._r8
            thl_d_org_pblh          = 0._r8
            qt_d_org_pblh           = 0._r8
            u_d_org_pblh            = 0._r8
            v_d_org_pblh            = 0._r8
            tr_d_org_pblh(1:ncnst)  = 0._r8 
         endif

         qt_dia_d_org        =     qt_dia_d_org / pblhp
         thl_dia_d_org       =    thl_dia_d_org / pblhp
         do mt = 1, ncnst 
            tr_dia_d_org(mt) = tr_dia_d_org(mt) / pblhp
         enddo

         ! --------------------------------------------------------------------------------------------------------- !
         ! May.1.2014.                                                                                               !
         ! Add exactly same part as the above but for 'cmf_d_orgU_pblh' that exclusively sinks down into '1-awk_PBL' !
         ! instead of 'awk_PBL' for 'i_budget_coldpool = 1,2' treatment (i.e., budget consistent cold pool).         !
         ! --------------------------------------------------------------------------------------------------------- !

         if( cmf_d_orgU_pblh .gt. nonzero ) then
             thl_d_orgU_pblh          =     thl_d_orgU_pblh / cmf_d_orgU_pblh
             qt_d_orgU_pblh           =      qt_d_orgU_pblh / cmf_d_orgU_pblh 
             u_d_orgU_pblh            =       u_d_orgU_pblh / cmf_d_orgU_pblh 
             v_d_orgU_pblh            =       v_d_orgU_pblh / cmf_d_orgU_pblh 
             do mt = 1, ncnst
                tr_d_orgU_pblh(mt)    =  tr_d_orgU_pblh(mt) / cmf_d_orgU_pblh
             enddo 
         else
             cmf_d_orgU_pblh          = 0._r8
             thl_d_orgU_pblh          = 0._r8
             qt_d_orgU_pblh           = 0._r8
             u_d_orgU_pblh            = 0._r8
             v_d_orgU_pblh            = 0._r8
             tr_d_orgU_pblh(1:ncnst)  = 0._r8 
         endif

         qt_dia_d_orgU        =     qt_dia_d_orgU / pblhp
         thl_dia_d_orgU       =    thl_dia_d_orgU / pblhp
         do mt = 1, ncnst 
            tr_dia_d_orgU(mt) = tr_dia_d_orgU(mt) / pblhp
         enddo

         ! --------------------------------------------------------------------------------------------------------- ! 
         ! Apr.21.2011. Below block is to compute diagnostic 'am_d,qlm_d,qim_d' extending element computation above. !
         !              This computation is same as corresnding computation of updraft ( am_u,qlm_u,qim_u ).         !
         ! --------------------------------------------------------------------------------------------------------- !  

         do msfc = 1, nseg
            do k = 1, ktop_msfc(msfc)
               qlm_d_msfc(k,msfc)  = qlm_d_msfc(k,msfc) / max( am_d_msfc(k,msfc), nonzero )
               qim_d_msfc(k,msfc)  = qim_d_msfc(k,msfc) / max( am_d_msfc(k,msfc), nonzero )
               am_d(k)  =  am_d(k) +  am_d_msfc(k,msfc)
               qlm_d(k) = qlm_d(k) +  am_d_msfc(k,msfc) * qlm_d_msfc(k,msfc)
               qim_d(k) = qim_d(k) +  am_d_msfc(k,msfc) * qim_d_msfc(k,msfc)
            enddo
         enddo
         qlm_d(k) = qlm_d(k) / max( am_d(k), nonzero )
         qim_d(k) = qim_d(k) / max( am_d(k), nonzero )

         ! ------------------------------------------ !
         ! Compute grid mean properties in each layer !
         ! ------------------------------------------ !

         do k = 1, ktop             ! This 'k'   is a layer index

            km = k - 1              ! This 'km'  is a base interface index of 'k'

            ! --------------------------------------------------------------------- !
            ! Mean downdraft properties at the base interface for each 'msfc' index !
            ! These are diagnostic quantities.                                      !
            ! --------------------------------------------------------------------- ! 

            do msfc = 1, nseg
               if( cmf_d_msfc(km,msfc) .gt. nonzero ) then
                  thl_d_msfc(km,msfc) = thl_d_msfc(km,msfc) / cmf_d_msfc(km,msfc) 
                  qt_d_msfc(km,msfc)  =  qt_d_msfc(km,msfc) / cmf_d_msfc(km,msfc)
                  u_d_msfc(km,msfc)   =   u_d_msfc(km,msfc) / cmf_d_msfc(km,msfc)
                  v_d_msfc(km,msfc)   =   v_d_msfc(km,msfc) / cmf_d_msfc(km,msfc)
                  w_d_msfc(km,msfc)   =   w_d_msfc(km,msfc) / cmf_d_msfc(km,msfc)
                  ql_d_msfc(km,msfc)  =  ql_d_msfc(km,msfc) / cmf_d_msfc(km,msfc)
                  qi_d_msfc(km,msfc)  =  qi_d_msfc(km,msfc) / cmf_d_msfc(km,msfc)
                  wa_d_msfc(km,msfc)  =  wa_d_msfc(km,msfc) /   a_d_msfc(km,msfc)
                  qla_d_msfc(km,msfc) = qla_d_msfc(km,msfc) /   a_d_msfc(km,msfc)
                  qia_d_msfc(km,msfc) = qia_d_msfc(km,msfc) /   a_d_msfc(km,msfc)
                  do mt = 1, ncnst
                     tr_d_msfc(km,msfc,mt) = tr_d_msfc(km,msfc,mt) / cmf_d_msfc(km,msfc)
                  enddo
               else
                  thl_d_msfc(km,msfc) = 0._r8
                  qt_d_msfc(km,msfc)  = 0._r8
                  u_d_msfc(km,msfc)   = 0._r8
                  v_d_msfc(km,msfc)   = 0._r8
                  w_d_msfc(km,msfc)   = 0._r8
                  ql_d_msfc(km,msfc)  = 0._r8
                  qi_d_msfc(km,msfc)  = 0._r8
                  wa_d_msfc(km,msfc)  = 0._r8
                  qla_d_msfc(km,msfc) = 0._r8
                  qia_d_msfc(km,msfc) = 0._r8
                  do mt = 1, ncnst
                     tr_d_msfc(km,msfc,mt) = 0._r8
                  enddo
               endif
            enddo

            ! ----------------------------------------------- !
            ! Mean downdraft properties at the base interface !
            ! ----------------------------------------------- ! 

            if( cmf_d(km) .gt. nonzero ) then
               thl_d(km) = thl_d(km) / cmf_d(km)
               qt_d(km)  =  qt_d(km) / cmf_d(km) 
               u_d(km)   =   u_d(km) / cmf_d(km)
               v_d(km)   =   v_d(km) / cmf_d(km)
               w_d(km)   =   w_d(km) / cmf_d(km)
               wa_d(km)  =  wa_d(km) /   a_d(km)
               ql_d(km)  =  ql_d(km) / cmf_d(km)
               qi_d(km)  =  qi_d(km) / cmf_d(km)
               qla_d(km) = qla_d(km) /   a_d(km)
               qia_d(km) = qia_d(km) /   a_d(km)
               do mt = 1, ncnst
                  tr_d(km,mt) = tr_d(km,mt) / cmf_d(km)
               enddo
            else
               cmf_d(km) = 0._r8
               thl_d(km) = 0._r8
               qt_d(km)  = 0._r8
               u_d(km)   = 0._r8
               v_d(km)   = 0._r8
               w_d(km)   = 0._r8
               wa_d(km)  = 0._r8
               ql_d(km)  = 0._r8
               qi_d(km)  = 0._r8
               qla_d(km) = 0._r8
               qia_d(km) = 0._r8
               do mt = 1, ncnst
                  tr_d(km,mt) = 0._r8
               enddo
            endif

            ! ------------------------------------------------- !
            ! Mean diabatic forcings on downdraft in each layer !
            ! ------------------------------------------------- !

            if( cmf_d_dia(k) .gt. nonzero ) then
               evp_thll_d(k)  =  evp_thll_d(k) / cmf_d_dia(k)
               evp_qtl_d(k)   =   evp_qtl_d(k) / cmf_d_dia(k) 
               evp_thli_d(k)  =  evp_thli_d(k) / cmf_d_dia(k)
               evp_qti_d(k)   =   evp_qti_d(k) / cmf_d_dia(k) 
               prep_thll_d(k) = prep_thll_d(k) / cmf_d_dia(k)
               prep_qtl_d(k)  =  prep_qtl_d(k) / cmf_d_dia(k) 
               prep_thli_d(k) = prep_thli_d(k) / cmf_d_dia(k)
               prep_qti_d(k)  =  prep_qti_d(k) / cmf_d_dia(k) 
               eff_ql_d(k)    =    eff_ql_d(k) / cmf_d_dia(k) 
               eff_qi_d(k)    =    eff_qi_d(k) / cmf_d_dia(k) 
               PGF_u_d(k)     =     PGF_u_d(k) / cmf_d_dia(k)
               PGF_v_d(k)     =     PGF_v_d(k) / cmf_d_dia(k)
               do mt = 1, ncnst
                  evp_tr_d(k,mt)  =  evp_tr_d(k,mt) / cmf_d_dia(k) 
                  prep_tr_d(k,mt) = prep_tr_d(k,mt) / cmf_d_dia(k)
                  eff_tr_d(k,mt)  =  eff_tr_d(k,mt) / cmf_d_dia(k) 
               enddo
            else
               cmf_d_dia(k)   = 0._r8
               evp_thll_d(k)  = 0._r8
               evp_qtl_d(k)   = 0._r8
               evp_thli_d(k)  = 0._r8
               evp_qti_d(k)   = 0._r8
               prep_thll_d(k) = 0._r8
               prep_qtl_d(k)  = 0._r8
               prep_thli_d(k) = 0._r8
               prep_qti_d(k)  = 0._r8
               eff_ql_d(k)    = 0._r8
               eff_qi_d(k)    = 0._r8
               PGF_u_d(k)     = 0._r8
               PGF_v_d(k)     = 0._r8
               do mt = 1, ncnst
                  evp_tr_d(k,mt)  = 0._r8
                  prep_tr_d(k,mt) = 0._r8
                  eff_tr_d(k,mt)  = 0._r8
               enddo
            endif

            ! ------------------------------------------------------- !
            ! Mean detrained properties from downdraft in each layer. !
            ! While all detrainment occurs at the base interface,     !
            ! we define detrained properties in each layer.           !
            ! ------------------------------------------------------- !

            if( cmf_rd(k) .gt. nonzero ) then
               thl_rd(k)  =  thl_rd(k) / cmf_rd(k)
               qt_rd(k)   =   qt_rd(k) / cmf_rd(k) 
               u_rd(k)    =    u_rd(k) / cmf_rd(k)
               v_rd(k)    =    v_rd(k) / cmf_rd(k)
               ql_rd(k)   =   ql_rd(k) / cmf_rd(k)
               qi_rd(k)   =   qi_rd(k) / cmf_rd(k)
               do mt = 1, ncnst
                  tr_rd(k,mt) = tr_rd(k,mt) / cmf_rd(k) 
               enddo
            else
               cmf_rd(k)  = 0._r8
               thl_rd(k)  = 0._r8
               qt_rd(k)   = 0._r8
               u_rd(k)    = 0._r8
               v_rd(k)    = 0._r8
               ql_rd(k)   = 0._r8 
               qi_rd(k)   = 0._r8
               do mt = 1, ncnst
                  tr_rd(k,mt) = 0._r8
               enddo
            endif

            ! ---------------------------------------------------------------------- !
            ! Computation of mean properties of all detrained source airs by summing !
            ! 1 detrained air from updrafts and 3 detrained airs from downdrafts.    !                  
            ! ---------------------------------------------------------------------- !

            cmf_r(k) =                cmf_ru(k) +             cmf_rd(k)             
            thl_r(k) =    thl_ru(k) * cmf_ru(k) + thl_rd(k) * cmf_rd(k) 
            qt_r(k)  =     qt_ru(k) * cmf_ru(k) +  qt_rd(k) * cmf_rd(k) 
            u_r(k)   =      u_ru(k) * cmf_ru(k) +   u_rd(k) * cmf_rd(k)
            v_r(k)   =      v_ru(k) * cmf_ru(k) +   v_rd(k) * cmf_rd(k)
            ql_r(k)  =     ql_ru(k) * cmf_ru(k) +  ql_rd(k) * cmf_rd(k)
            qi_r(k)  =     qi_ru(k) * cmf_ru(k) +  qi_rd(k) * cmf_rd(k)
            do mt = 1, ncnst 
               tr_r(k,mt) = tr_ru(k,mt) * cmf_ru(k) + tr_rd(k,mt) * cmf_rd(k) 
            enddo
            if( cmf_r(k) .gt. nonzero ) then
               thl_r(k) = thl_r(k) / cmf_r(k)
               qt_r(k)  =  qt_r(k) / cmf_r(k) 
               u_r(k)   =   u_r(k) / cmf_r(k)
               v_r(k)   =   v_r(k) / cmf_r(k)
               ql_r(k)  =  ql_r(k) / cmf_r(k)
               qi_r(k)  =  qi_r(k) / cmf_r(k)
               do mt = 1, ncnst
                  tr_r(k,mt) = tr_r(k,mt) / cmf_r(k) 
               enddo
            else
               cmf_r(k) = 0._r8
               thl_r(k) = 0._r8
               qt_r(k)  = 0._r8
               u_r(k)   = 0._r8
               v_r(k)   = 0._r8
               ql_r(k)  = 0._r8
               qi_r(k)  = 0._r8
               do mt = 1, ncnst
                  tr_r(k,mt) = 0._r8
               enddo
            endif

            ! ------------------------------------------------------------- !
            ! Treatment of detrained air purely from the convective updraft !
            ! ------------------------------------------------------------- !

            cmf_r2(k) =                 cmf_ru2(k) +             cmf_rd(k)             
            thl_r2(k) =    thl_ru2(k) * cmf_ru2(k) + thl_rd(k) * cmf_rd(k) 
            qt_r2(k)  =     qt_ru2(k) * cmf_ru2(k) +  qt_rd(k) * cmf_rd(k) 
            u_r2(k)   =      u_ru2(k) * cmf_ru2(k) +   u_rd(k) * cmf_rd(k)
            v_r2(k)   =      v_ru2(k) * cmf_ru2(k) +   v_rd(k) * cmf_rd(k)
            ql_r2(k)  =     ql_ru2(k) * cmf_ru2(k) +  ql_rd(k) * cmf_rd(k)
            qi_r2(k)  =     qi_ru2(k) * cmf_ru2(k) +  qi_rd(k) * cmf_rd(k)
            do mt = 1, ncnst 
               tr_r2(k,mt) = tr_ru2(k,mt) * cmf_ru2(k) + tr_rd(k,mt) * cmf_rd(k) 
            enddo
            if( cmf_r2(k) .gt. nonzero ) then
               thl_r2(k) = thl_r2(k) / cmf_r2(k)
               qt_r2(k)  =  qt_r2(k) / cmf_r2(k) 
               u_r2(k)   =   u_r2(k) / cmf_r2(k)
               v_r2(k)   =   v_r2(k) / cmf_r2(k)
               ql_r2(k)  =  ql_r2(k) / cmf_r2(k)
               qi_r2(k)  =  qi_r2(k) / cmf_r2(k)
               do mt = 1, ncnst
                  tr_r2(k,mt) = tr_r2(k,mt) / cmf_r2(k) 
               enddo
            else
               cmf_r2(k) = 0._r8
               thl_r2(k) = 0._r8
               qt_r2(k)  = 0._r8
               u_r2(k)   = 0._r8
               v_r2(k)   = 0._r8
               ql_r2(k)  = 0._r8
               qi_r2(k)  = 0._r8
               do mt = 1, ncnst
                  tr_r2(k,mt) = 0._r8
               enddo
            endif

            ! 'Flux-convergence' and 'Subsidence-detrainment' consistent diagnostic output for use in the macrophysics.

            cmf_det(k) = cmf_r2(k)
            ql_det(k)  =  ql_r2(k)
            qi_det(k)  =  qi_r2(k)


            ! ---------------------------------------------------------------------------------------------------------------------- !
            ! For more clear treatment fully consistent with the governing tendency equations, I can use 'cmf_r2(k), thl_r2(k),....' !
            ! in the below 'if' block instead of 'cmf_r(k), thl_r(k),....'. This should be tested in a near future.                  !
            ! I should do this today, after verifying that the modified code exactly reproduces the previous results.                !
            ! Nov.14.2014. I added 'i_detrain' option: '0' is previous formula-inconsistent default,                                 !
            !                                          '1' is a new formula consistent one.                                          !
            ! ---------------------------------------------------------------------------------------------------------------------- ! 

            if( i_detrain .eq. 0 ) then

               if( cmf_r(k) .gt. nonzero ) then
                  cu_cmfr_mxen(k,iter)      =      cmf_r(k)
                  cu_thlr_mxen(k,iter)      =      thl_r(k) -   thl0(k)
                  cu_qtr_mxen(k,iter)       =       qt_r(k) -    qt0(k)
                  cu_ur_mxen(k,iter)        =        u_r(k) -     u0(k) 
                  cu_vr_mxen(k,iter)        =        v_r(k) -     v0(k)
                  cu_qlr_mxen(k,iter)       =       ql_r(k) -    ql0(k)
                  cu_qir_mxen(k,iter)       =       qi_r(k) -    qi0(k)
                  do mt = 1, ncnst
                     cu_trr_mxen(k,mt,iter) =    tr_r(k,mt) - tr0(k,mt)
                  enddo
               endif
               
            else

               if( cmf_r2(k) .gt. nonzero ) then
                  cu_cmfr_mxen(k,iter)      =      cmf_r2(k)
                  cu_thlr_mxen(k,iter)      =      thl_r2(k) -   thl0(k)
                  cu_qtr_mxen(k,iter)       =       qt_r2(k) -    qt0(k)
                  cu_ur_mxen(k,iter)        =        u_r2(k) -     u0(k) 
                  cu_vr_mxen(k,iter)        =        v_r2(k) -     v0(k)
                  cu_qlr_mxen(k,iter)       =       ql_r2(k) -    ql0(k)
                  cu_qir_mxen(k,iter)       =       qi_r2(k) -    qi0(k)
                  do mt = 1, ncnst
                     cu_trr_mxen(k,mt,iter) =    tr_r2(k,mt) - tr0(k,mt)
                  enddo
               endif

            endif

            cu_cmfum_mxen(k,iter)         =  cmf_u_mix(k)

            if( cmf_rd(k) .gt. nonzero ) then
               cu_cmfrd_mxen(k,iter)         =     cmf_rd(k)
               cu_thlrd_mxen(k,iter)         =     thl_rd(k) -   thl0(k) 
               cu_qtrd_mxen(k,iter)          =      qt_rd(k) -    qt0(k)
               cu_urd_mxen(k,iter)           =       u_rd(k) -     u0(k)
               cu_vrd_mxen(k,iter)           =       v_rd(k) -     v0(k)
               cu_qlrd_mxen(k,iter)          =      ql_rd(k) -    ql0(k) 
               cu_qird_mxen(k,iter)          =      qi_rd(k) -    qi0(k)
               do mt = 1, ncnst
                  cu_trrd_mxen(k,mt,iter)    =   tr_rd(k,mt) - tr0(k,mt)
               enddo
            endif

         enddo

         ! --------------------------------------------------------------------------------------- !
         ! Compute rain and snow production tendencies & effective tendencies of cloud condensate. !
         ! Important Sanity Check for Positive Precipitation Flux.                                 !
         ! --------------------------------------------------------------------------------------- !

         do k = ktop, 1, -1  ! This is a layer index 
            km = k - 1
            qlten_eff_u(k) = ( g / dp0(k) ) * eff_ql_u(k) * cmf_u_dia(k) 
            qiten_eff_u(k) = ( g / dp0(k) ) * eff_qi_u(k) * cmf_u_dia(k) 
            qlten_eff_d(k) = ( g / dp0(k) ) * eff_ql_d(k) * cmf_d_dia(k) 
            qiten_eff_d(k) = ( g / dp0(k) ) * eff_qi_d(k) * cmf_d_dia(k) 
            do mt = 1, ncnst
               trten_eff_u(k,mt) = ( g / dptr0(k,mt) ) * eff_tr_u(k,mt) * cmf_u_dia(k) 
               trten_eff_d(k,mt) = ( g / dptr0(k,mt) ) * eff_tr_d(k,mt) * cmf_d_dia(k) 
            enddo
            ! ------------------------------------------------------------------------------------- !
            ! Below considers 'production  of precipitation within   updraft ( prep_qtl_u < 0 ) and !
            !                 'evaporation of precipitation within downdraft (  evp_qtl_d > 0 ).    !
            !  1. By multiplying mass flux, the updraft/downdraft area information are included.    !
            !  2. In a certain layer, qrten, qsten can be negative due to the evaporation. However, !
            !     its downward integrated 'rainflx, snowflx' is always positive at all interfaces   !
            !     due to the how 'evp_qtl_d' is computed.                                           !
            !  3. Snow melting effect is included into 'qrten(k), qsten(k)'. Note that it should    !
            !     not be included into 'qrten_u(k), qsten_u(k)' - very important.                   ! 
            ! Nov.29.2012. Add tracer 'trrsten' associated with rain/snow components.               !
            ! Feb.06.2013. Note that wet deposition of aerosols (both interstitial and cloud-borne) !
            !              within convective updraft and downdraft will be treated as a part of     !
            !              prep_tr_u(k,mt) and prep_tr_d(k,mt). Thus, I don't need separate use of  !
            !              wdep_tr_u(k,mt) and wdep_tr_d(k,mt) in the below computation of          !
            !              trten_dia_u(k,mt) and trten_dia_d(k,mt).                                 !   
            ! ------------------------------------------------------------------------------------- !
            thlten_dia_u(k) = ( g / dp0(k) ) * ( prep_thll_u(k) + prep_thli_u(k) +  evp_thll_u(k) +  evp_thli_u(k) ) * cmf_u_dia(k)
            qtten_dia_u(k)  = ( g / dp0(k) ) * ( prep_qtl_u(k)  + prep_qti_u(k)  +  evp_qtl_u(k)  +  evp_qti_u(k)  ) * cmf_u_dia(k)
            qlten_dia_u(k)  = ( g / dp0(k) ) * ( prep_qtl_u(k)                                                     ) * cmf_u_dia(k)
            qiten_dia_u(k)  = ( g / dp0(k) ) * ( prep_qti_u(k)                                                     ) * cmf_u_dia(k)
            thlten_dia_d(k) = ( g / dp0(k) ) * ( prep_thll_d(k) + prep_thli_d(k) +  evp_thll_d(k) +  evp_thli_d(k) ) * cmf_d_dia(k)
            qtten_dia_d(k)  = ( g / dp0(k) ) * ( prep_qtl_d(k)  + prep_qti_d(k)  +  evp_qtl_d(k)  +  evp_qti_d(k)  ) * cmf_d_dia(k)
            qlten_dia_d(k)  = ( g / dp0(k) ) * ( prep_qtl_d(k)                                                     ) * cmf_d_dia(k)
            qiten_dia_d(k)  = ( g / dp0(k) ) * ( prep_qti_d(k)                                                     ) * cmf_d_dia(k)
            do mt = 1, ncnst
               trten_dia_u(k,mt) = ( g / dptr0(k,mt) ) * ( prep_tr_u(k,mt) + evp_tr_u(k,mt) ) * cmf_u_dia(k) 
               trten_dia_d(k,mt) = ( g / dptr0(k,mt) ) * ( prep_tr_d(k,mt) + evp_tr_d(k,mt) ) * cmf_d_dia(k) 
            enddo
            qrten_u(k) = - ( g / dp0(k) ) * ( prep_qtl_u(k) + evp_qtl_u(k) ) * cmf_u_dia(k)   ! >= 0
            qrten_d(k) = - ( g / dp0(k) ) * ( prep_qtl_d(k) + evp_qtl_d(k) ) * cmf_d_dia(k)   ! <= 0
            qsten_u(k) = - ( g / dp0(k) ) * ( prep_qti_u(k) + evp_qti_u(k) ) * cmf_u_dia(k)   ! >= 0
            qsten_d(k) = - ( g / dp0(k) ) * ( prep_qti_d(k) + evp_qti_d(k) ) * cmf_d_dia(k)   ! <= 0
            qrten(k)   =   qrten_u(k) + qrten_d(k)
            qsten(k)   =   qsten_u(k) + qsten_d(k)
         end do

         ! ------------------------------------------------------------------------ !
         ! Grid-mean tendencies of non-conservative scalars: 'ql,qi'                !
         ! This is a sum of compensating subsidence and condensate detrainment.     !
         ! Note that 'ql(i)flx_u(i)' at surface is set to zero even though cmf_d(0) !
         ! was explicitly computed above. As a result, the tendency computation is  !
         ! absolutely correct, which is good.                                       ! 
         ! Mar.19.2014. Note that tendency computation will be done further below   !
         !              after computing flux first. As discussed further below,     !
         !              should compute non-zero convective flux (both updraft and   !
         !              downdraft) for use with 'ipartition = 1' option.            ! 
         ! ------------------------------------------------------------------------ ! 
         
         do k = 1, ktop - 1    ! This is a top interface index or layer index
            km = k - 1
            kp = k + 1
            um = abs( g * cmf_u(k) * dt )
            dm = abs( g * cmf_d(k) * dt )
            call envcon_flux( k, mkx, um, dm, ql0(1:mkx), ssql0(1:mkx), ps0(0:mkx), ql_env_u, ql_env_d )
            call envcon_flux( k, mkx, um, dm, qi0(1:mkx), ssqi0(1:mkx), ps0(0:mkx), qi_env_u, qi_env_d )
            ql_env_ua(k)  = ql_env_u 
            qi_env_ua(k)  = qi_env_u
            ql_env_da(kp) = ql_env_d 
            qi_env_da(kp) = qi_env_d
            if( iflux_env .eq. 0 ) then
               ql_env_ua(k)  = ql0bot(kp)
               qi_env_ua(k)  = qi0bot(kp)
               ql_env_da(kp) = ql0top(k) 
               qi_env_da(kp) = qi0top(k)
            endif
            qlflx_u(k) =   cmf_u(k) * ( ql_u(k)  -   ql_env_ua(k) )
            qiflx_u(k) =   cmf_u(k) * ( qi_u(k)  -   qi_env_ua(k) )
            qlflx_d(k) = - cmf_d(k) * ( ql_d(k)  -  ql_env_da(kp) )
            qiflx_d(k) = - cmf_d(k) * ( qi_d(k)  -  qi_env_da(kp) )
            do mt = 1, ncnst
               call envcon_flux( k, mkx, um, dm, tr0(1:mkx,mt), sstr0(1:mkx,mt), ps0(0:mkx), tr_env_u, tr_env_d )
               if( iflux_env .eq. 0 ) then
                  tr_env_u = tr0bot(kp,mt)
                  tr_env_d = tr0top(k,mt)
               endif
               trflx_u(k,mt) =   cmf_u(k) * ( tr_u(k,mt)  -  tr_env_u )
               trflx_d(k,mt) = - cmf_d(k) * ( tr_d(k,mt)  -  tr_env_d )
            enddo
         enddo
         ! ------------------------------------------------------------------------------------------- !
         ! Mar.19.2014. Compute downdraft flux at the base interface.                                  !
         !              Similar to the treatment of the updraft flux with 'ipartition = 1' option,     !
         !              I should compute the downdraft flux at the surface interface, and partition    !
         !              the downdraft surface flux uniformly over the entire cumulus layer or the      !
         !              PBL depth. I should do this for the 'ql,qi' and 'tracers' too.                 ! 
         !              This treatment of the 'ql,qi,tracer' fluxes are done here before computing     !
         !              corresponding tendencies, using the exactly same method used here below for    !
         !              conservative scalars.                                                          !
         !              Note that at this stage, we have already computed non-zero cmf_d(0), thl_d(0), ! 
         !              qt_d(0), u_d(0), v_d(0), so that below treatment is perfectly correct.         !
         !              Similar computation of convective fluxes at the surface for 'slflx_d(0),       !
         !              qtflx_d(0), uflx_d(0), vflx_d(0)' will be done later. The same computation for !
         !              convective updraft fluxes at the surface for all the scalars wre already done  !
         !              at the beginning portion of this module. Thus below computation for 'ql,qi'    !
         !              and 'tracers' of convective downdraft fluxes at the surface is enough.         !
         ! ------------------------------------------------------------------------------------------- !
         qlflx_d(0)       = -cmf_d(0) * ( ql_d(0)  -  ql0bot(1) )
         qiflx_d(0)       = -cmf_d(0) * ( qi_d(0)  -  qi0bot(1) )
         do mt = 1, ncnst
            trflx_d(0,mt) = -cmf_d(0) * ( tr_d(0,mt) - tr0bot(1,mt) )
         enddo

         ! ----------------------------------------------------------------------------------- !
         ! Mar.19.2014. Computation of tendencies of 'ql,qi,tracers' using the fluxes.         !
         !              Note that we are using non-zero surface fluxes here for consistent use ! 
         !              with the other conservativee scalars with ipartition=1.                !
         ! ----------------------------------------------------------------------------------- !

         do k = 1, ktop
            km = k - 1
            qlten_u(k) = ( g / dp0(k) ) * ( qlflx_u(km) - qlflx_u(k) ) + qlten_dia_u(k) + qlten_eff_u(k)
            qiten_u(k) = ( g / dp0(k) ) * ( qiflx_u(km) - qiflx_u(k) ) + qiten_dia_u(k) + qiten_eff_u(k)
            qlten_d(k) = ( g / dp0(k) ) * ( qlflx_d(km) - qlflx_d(k) ) + qlten_dia_d(k) + qlten_eff_d(k)
            qiten_d(k) = ( g / dp0(k) ) * ( qiflx_d(km) - qiflx_d(k) ) + qiten_dia_d(k) + qiten_eff_d(k)
            do mt = 1, ncnst
               trten_u(k,mt) = ( g / dptr0(k,mt) ) * ( trflx_u(km,mt) - trflx_u(k,mt) ) + trten_dia_u(k,mt) + trten_eff_u(k,mt)
               trten_d(k,mt) = ( g / dptr0(k,mt) ) * ( trflx_d(km,mt) - trflx_d(k,mt) ) + trten_dia_d(k,mt) + trten_eff_d(k,mt)
            enddo

         end do

         ! ----------------------------------------------------- !
         ! Tendency due to Compensating Subsidence : [ kg/kg/s ] !
         ! ----------------------------------------------------- ! 

         do k = 1, ktop
            km = k - 1
            kp = k + 1
            if( k .eq. 1 ) then
               qlten_sub(k) = g * 0.5_r8 * ( cmf_u(k) + cmf_u(km) - cmf_d(k) - cmf_d(km) ) * &
                  ( ql0(k+1) -  ql0(1) ) / ( p0(1) - p0(k+1) ) 
               qiten_sub(k) = g * 0.5_r8 * ( cmf_u(k) + cmf_u(km) - cmf_d(k) - cmf_d(km) ) * &
                  ( qi0(k+1) -  qi0(1) ) / ( p0(1) - p0(k+1) ) 
            elseif( k .eq. ktop ) then     
               qlten_sub(k) = g * 0.5_r8 * ( cmf_u(k) + cmf_u(km) - cmf_d(k) - cmf_d(km) ) * &
                  ( ql0(ktop) - ql0(k-1) ) / ( p0(k-1) - p0(ktop) ) 
               qiten_sub(k) = g * 0.5_r8 * ( cmf_u(k) + cmf_u(km) - cmf_d(k) - cmf_d(km) ) * &
                  ( qi0(ktop) - qi0(k-1) ) / ( p0(k-1) - p0(ktop) ) 
            else                           
               qlten_sub(k) = g * 0.5_r8 * ( cmf_u(k) + cmf_u(km) - cmf_d(k) - cmf_d(km) ) * &
                  ( ql0(k+1) -  ql0(k-1) ) / ( p0(k-1) - p0(k+1) ) 
               qiten_sub(k) = g * 0.5_r8 * ( cmf_u(k) + cmf_u(km) - cmf_d(k) - cmf_d(km) ) * &
                  ( qi0(k+1)  - qi0(k-1) ) / ( p0(k-1) - p0(k+1) ) 
            endif
         enddo

         ! ---------------------------------------------------- !
         ! Tendency due to Condensate Detrainment : [ kg/kg/s ] !
         ! ---------------------------------------------------- ! 

         rliq = 0._r8
         do k = 1, ktop   ! This is a layer index
            km = k - 1    ! This is a base interface index
            ! ------------------------------------------------------------------------------------------------------ !
            ! Nov.04.2014. For the fully consistent treatment directly relating 'flux-convergence' formula to        !
            ! the 'subsidence-detrainment' formula, I should use 'cmf_r2(k), ql_r2(k), qi_r2(k)' in the below block. ! 
            ! Thus, I changed to the correct formula today.                                                          !
            ! ------------------------------------------------------------------------------------------------------ !
          ! qlten_det(k) = cmf_r(k) * ( ql_r(k) - ql0(k) ) * ( g / dp0(k) )
          ! qiten_det(k) = cmf_r(k) * ( qi_r(k) - qi0(k) ) * ( g / dp0(k) )
            qlten_det(k) = cmf_r2(k) * ( ql_r2(k) - ql0(k) ) * ( g / dp0(k) )
            qiten_det(k) = cmf_r2(k) * ( qi_r2(k) - qi0(k) ) * ( g / dp0(k) )
            rqc_l(k) = 0._r8
            rqc_i(k) = 0._r8
            rnc_l(k) = 0._r8
            rnc_i(k) = 0._r8
            rqc(k)   = rqc_l(k) + rqc_i(k)
            rliq     = rliq  + rqc(k) * dp0(k) / g / 1000._r8       ! [ liquid m/s ]
         enddo ! k = 1, ktop. Here, 'k' is a layer index.

         ! -------------------------------------------- !
         !                                              !
         ! Grid-Mean Tendencies of Conservative Scalars !
         !                                              !
         ! -------------------------------------------- !

         ! ----------------------------------------------------------------- !
         ! When convective updraft plumes launch from the surface ( k = 0 ), ! 
         ! also specify non-zero updraft flux value at surface.              !
         ! By defult, the flux at surface is set to be zero.                 ! 
         ! As expected, this 'kiss=0' causes conservation error.             !
         ! However, by setting I_cri = 0 when kiss = 0, we can avoid         !
         ! conservation error in a most reasonable way.                      ! 
         ! ----------------------------------------------------------------- !     

         ! ---------------------------------- ! 
         ! Tendency due to Convective Updraft !
         ! ---------------------------------- !

         do k = 1, ktop - 1    ! This is a layer index        
            um = abs( g * cmf_u(k) * dt )
            dm = 0._r8
            call envcon_flux( k, mkx, um, dm, thl0(1:mkx), ssthl0(1:mkx), ps0(0:mkx), thl_env_u, thl_env_d )
            call envcon_flux( k, mkx, um, dm, qt0(1:mkx), ssqt0(1:mkx), ps0(0:mkx), qt_env_u, qt_env_d )
            call envcon_flux( k, mkx, um, dm, u0(1:mkx), ssu0(1:mkx), ps0(0:mkx), u_env_u, u_env_d )
            call envcon_flux( k, mkx, um, dm, v0(1:mkx), ssv0(1:mkx), ps0(0:mkx), v_env_u, v_env_d )
            ! ---------------------------- !
            ! Use UW approach as an option !
            ! ---------------------------- !
            if( iflux_env .eq. 0 ) then
               thl_env_u = thl0bot(k+1)
               qt_env_u  = qt0bot(k+1)
               u_env_u   = u0bot(k+1)
               v_env_u   = v0bot(k+1)
            endif

            slflx_u(k) = cp * exns0(k) * cmf_u(k) * ( thl_u(k) - thl_env_u )
            qtflx_u(k) =                 cmf_u(k) * ( qt_u(k)  -  qt_env_u )
            uflx_u(k)  =                 cmf_u(k) * ( u_u(k)   -   u_env_u )
            vflx_u(k)  =                 cmf_u(k) * ( v_u(k)   -   v_env_u )

         end do

         do k = 1, ktop
            km = k - 1
            slten_u(k) = ( g / dp0(k) ) * ( slflx_u(km) - slflx_u(k) ) + cp * exn0(k) * thlten_dia_u(k)
            qtten_u(k) = ( g / dp0(k) ) * ( qtflx_u(km) - qtflx_u(k) ) + qtten_dia_u(k)
            uten_u(k)  = ( g / dp0(k) ) * ( uflx_u(km)  -  uflx_u(k) )
            vten_u(k)  = ( g / dp0(k) ) * ( vflx_u(km)  -  vflx_u(k) )
            sten_u(k)  = slten_u(k) + xlv * qlten_u(k) + xls * qiten_u(k)
            qvten_u(k) = qtten_u(k) - qlten_u(k) - qiten_u(k)
         end do

         ! ------------------------------------------------------------------------ ! 
         ! Tendency due to Convective Downdraft                                     !
         ! Note that the fluxes at surface are set to zero even though cmf_d(0)     !
         ! was explicitly computed above. As a result, the tendency computation is  !
         ! absolutely correct, which is good.                                       ! 
         ! ------------------------------------------------------------------------ ! 

         do k = 1, ktop - 1    ! This is a top interface index or layer index
            um = abs( g * cmf_u(k) * dt )
            dm = abs( g * cmf_d(k) * dt )
            call envcon_flux( k, mkx, um, dm, thl0(1:mkx), ssthl0(1:mkx), ps0(0:mkx), thl_env_u, thl_env_d )
            call envcon_flux( k, mkx, um, dm, qt0(1:mkx), ssqt0(1:mkx), ps0(0:mkx), qt_env_u, qt_env_d )
            call envcon_flux( k, mkx, um, dm, u0(1:mkx), ssu0(1:mkx), ps0(0:mkx), u_env_u, u_env_d )
            call envcon_flux( k, mkx, um, dm, v0(1:mkx), ssv0(1:mkx), ps0(0:mkx), v_env_u, v_env_d )

            ! ---------------------------- !
            ! Use UW approach as an option !
            ! ---------------------------- !
            if( iflux_env .eq. 0 ) then
               thl_env_d = thl0top(k)
               qt_env_d  = qt0top(k)
               u_env_d   = u0top(k)
               v_env_d   = v0top(k)
            endif

            slflx_d(k) = - cp * exns0(k) * cmf_d(k) * ( thl_d(k) - thl_env_d )
            qtflx_d(k) = -                 cmf_d(k) * ( qt_d(k)  -  qt_env_d )
            uflx_d(k)  = -                 cmf_d(k) * ( u_d(k)   -   u_env_d )
            vflx_d(k)  = -                 cmf_d(k) * ( v_d(k)   -   v_env_d )

         end do
         ! ------------------------------------------------------------------------------------------- !
         ! Mar.19.2014. Compute downdraft flux at the base interface.                                  !
         !              Similar to the treatment of the updraft flux with 'ipartition = 1' option,     !
         !              I should compute the downdraft flux at the surface interface, and partition    !
         !              the downdraft surface flux uniformly over the entire cumulus layer or the      !
         !              PBL depth. I should do this for the 'ql,qi' and 'tracers' too.                 ! 
         !              This treatment of the 'ql,qi,tracer' fluxes are done above before computing    !
         !              corresponding tendencies, using the exactly same method used here below for    !
         !              conservative scalars.                                                          !
         !              Note that at this stage, we have already computed non-zero cmf_d(0), thl_d(0), ! 
         !              qt_d(0), u_d(0), v_d(0), so that below treatment is perfectly correct.         !
         ! ------------------------------------------------------------------------------------------- !
         slflx_d(0) = - cp * exns0(0) * cmf_d(0) * ( thl_d(0) - thl0bot(1) )
         qtflx_d(0) = -                 cmf_d(0) * ( qt_d(0)  -  qt0bot(1) )
         uflx_d(0)  = -                 cmf_d(0) * ( u_d(0)   -   u0bot(1) )
         vflx_d(0)  = -                 cmf_d(0) * ( v_d(0)   -   v0bot(1) )

         ! ----------------------------------------------------------------------------------------------------- !
         ! Aug.30.2011. Compute downdraft flux associated with convective organization at the PBL top interface. !
         !              Also compute buoyancy flux to compute 'convective organization' and 'vertical velocity   !
         !              perturbation at surface associated with convective organization.                         !
         ! Sep.07.2011. I should carefully choose whether I want to impose specific choice of only positive      !
         !              buoyancy flux.                                                                           ! 
         ! Sep.09.2011. New computation is made for imposing a full consistency and a full generality.           !
         !              Now, the case of negative buoyancy flux at the PBL top as well as positive buoyancy flux !
         !              is completely and generally treated.                                                     !
         ! ----------------------------------------------------------------------------------------------------- !
    
         if( cmf_d_org_pblh .gt. nonzero ) then 
            thlflx_d_org_pblh = - cmf_d_org_pblh * ( thl_d_org_pblh - thl0PBL )
            qtflx_d_org_pblh  = - cmf_d_org_pblh * ( qt_d_org_pblh  -  qt0PBL )
            uflx_d_org_pblh   = - cmf_d_org_pblh * ( u_d_org_pblh   -   u0PBL )
            vflx_d_org_pblh   = - cmf_d_org_pblh * ( v_d_org_pblh   -   v0PBL )
            do mt = 1, ncnst
               trflx_d_org_pblh(mt) = - cmf_d_org_pblh * ( tr_d_org_pblh(mt)  -  tr0PBL(mt) )
            enddo
         else
            cmf_d_org_pblh            = 0._r8
            thlflx_d_org_pblh         = 0._r8
            qtflx_d_org_pblh          = 0._r8
            uflx_d_org_pblh           = 0._r8
            vflx_d_org_pblh           = 0._r8
            ! thvflx_d_org_pblh         = 0._r8
            trflx_d_org_pblh(1:ncnst) = 0._r8
         endif

         ! -------------------------------------------------------------------------------------------------------------------------- !
         ! May.1.2014.                                                                                                                !
         ! For the treatment of budget consistent coldpool (i_budget_coldpool = 1,2), add the computation of 'qtflx_d_orgU_pblh' from !
         ! the downdrafts that exclusively sinks down into '1-awk_PBL' instead of 'awk_PBL' as the above, and also 'qtflx_u_org_pblh' !
         ! that is the flux by convective updraft defalted from the PBL. Note that 'qtflx_u_org_pblh' is slightly different from the  !
         ! already computed 'qtflx_u(kpblhm)' in that sense that 'qtflx_u_org_pblh' is using 'qt0PBL'. Since we are considering the   !
         ! budget of 'bulk' budget for coldpool treatment, we should use 'qtflx_u_org_pblh' instead of 'qtflx_u(kpblhm)' for the      !
         ! fully internally-self consistent treatment of cold pool.                                                                   !
         ! -------------------------------------------------------------------------------------------------------------------------- !

         if( cmf_d_orgU_pblh .gt. nonzero ) then 
             thlflx_d_orgU_pblh = - cmf_d_orgU_pblh * ( thl_d_orgU_pblh - thl0PBL )
             qtflx_d_orgU_pblh  = - cmf_d_orgU_pblh * ( qt_d_orgU_pblh  -  qt0PBL )
             uflx_d_orgU_pblh   = - cmf_d_orgU_pblh * ( u_d_orgU_pblh   -   u0PBL )
             vflx_d_orgU_pblh   = - cmf_d_orgU_pblh * ( v_d_orgU_pblh   -   v0PBL )
           ! Aug.30.2011. Modification is required in the below block for cloud liquid and ice droplet number 
           !              following the previous treatment above since we are assuming a certain fixed droplet
           !              radius for convective updraft and downdraft. However, since convective downdraft is
           !              likely not to have any condensate, below is likely to be OK for the time being.
             do mt = 1, ncnst
                trflx_d_orgU_pblh(mt) = - cmf_d_orgU_pblh * ( tr_d_orgU_pblh(mt)  -  tr0PBL(mt) )
             enddo
         else
             cmf_d_orgU_pblh            = 0._r8
             thlflx_d_orgU_pblh         = 0._r8
             qtflx_d_orgU_pblh          = 0._r8
             uflx_d_orgU_pblh           = 0._r8
             vflx_d_orgU_pblh           = 0._r8
             trflx_d_orgU_pblh(1:ncnst) = 0._r8
         endif

         cmf_u_org_pblh = cmf_u(kpblhm)
    !lim if( cmf_u_org_pblh .gt. nonzero ) then 
             thlflx_u_org_pblh          = cmf_u_org_pblh * ( thl_u(kpblhm) - thl0PBL )
             qtflx_u_org_pblh           = cmf_u_org_pblh * ( qt_u(kpblhm)  -  qt0PBL )
             uflx_u_org_pblh            = cmf_u_org_pblh * ( u_u(kpblhm)   -   u0PBL )
             vflx_u_org_pblh            = cmf_u_org_pblh * ( v_u(kpblhm)   -   v0PBL )
             do mt = 1, ncnst
                trflx_u_org_pblh(mt)    = cmf_u_org_pblh * ( tr_u(kpblhm,mt)  -  tr0PBL(mt) )
             enddo
    !lim else
    !lim     cmf_u_org_pblh             = 0._r8
    !lim     thlflx_u_org_pblh          = 0._r8
    !lim     qtflx_u_org_pblh           = 0._r8
    !lim     uflx_u_org_pblh            = 0._r8
    !lim     vflx_u_org_pblh            = 0._r8
    !lim     trflx_u_org_pblh(1:ncnst)  = 0._r8
    !lim endif

         ! -------------------------------------------------- !
         ! Compute grid-mean tendency by convective downdraft !
         ! -------------------------------------------------- !

         do k = 1, ktop
            km = k - 1
            slten_d(k) = ( g / dp0(k) ) * ( slflx_d(km) - slflx_d(k) ) + cp * exn0(k) * thlten_dia_d(k)
            qtten_d(k) = ( g / dp0(k) ) * ( qtflx_d(km) - qtflx_d(k) ) + qtten_dia_d(k)
            uten_d(k)  = ( g / dp0(k) ) * (  uflx_d(km) -  uflx_d(k) )
            vten_d(k)  = ( g / dp0(k) ) * (  vflx_d(km) -  vflx_d(k) )
            sten_d(k)  = slten_d(k) + xlv * qlten_d(k) + xls * qiten_d(k)
            qvten_d(k) = qtten_d(k) - qlten_d(k) - qiten_d(k)
         end do

         ! ------------------------------------------------------------------------------------- !
         ! Final Grid-Mean Evaporation Tendency                                                  !
         ! Nov.29.2012. Add corresponding tracer part. Here, we use proportional relationship    !
         !              using the fluxes of 'precipitation' and 'tracers' at the top interface.  !
         !              This seems to be correct and most reasonable approach.                   !
         !              If flx_all --> 0, then evprain_msfc(k,msfc) + evpsnow_msfc(k,msfc) --> 0,!
         !              so that unreasonable large-value is automatically prohibited, which is   !
         !              very good.                                                               !
         !              The same approach is used for the evaporation within downdraft.          !
         ! Dec.01.2012. Do separate treatment for precipitating droplet numbers.                 ! 
         ! ------------------------------------------------------------------------------------- !

         ! ------------------------------------------------------------ !
         ! Compute grid-mean tendencies by averaging or summing all the !
         ! updraft segment tendencies.                                  ! 
         ! ------------------------------------------------------------ !

         do k = ktop, 1, -1  ! 'k' is a layer index : 'mkx'('1') is the top ('bottom') layer          
            ! ------------------------------------------------------------ !
            ! Save the results into the array of grid-mean state variables.! 
            ! For tracers, I temporarily set it to be zero but it should   !
            ! be refined later.                                            !
            ! Oct.12.2010. Note that I added 'snowmlt_e(k)' in 'slten_evp' !
            ! ------------------------------------------------------------ !
            do msfc = 1, nseg
               snowmlt_e(k)   = snowmlt_e(k)   +   snowmlt_e_msfc(k,msfc) 
               ntraprd(k)     = ntraprd(k)     +     ntraprd_msfc(k,msfc) 
               ntsnprd(k)     = ntsnprd(k)     +     ntsnprd_msfc(k,msfc)
               evprain_e(k)   = evprain_e(k)   +   evprain_e_msfc(k,msfc)
               evpsnow_e(k)   = evpsnow_e(k)   +   evpsnow_e_msfc(k,msfc)
               evprain_d(k)   = evprain_d(k)   +   evprain_d_msfc(k,msfc)
               evpsnow_d(k)   = evpsnow_d(k)   +   evpsnow_d_msfc(k,msfc)
               flxrain(k)     = flxrain(k)     +     flxrain_msfc(k,msfc)
               flxsnow(k)     = flxsnow(k)     +     flxsnow_msfc(k,msfc)
               cvp_rainprd(k) = cvp_rainprd(k) + cvp_rainprd_msfc(k,msfc)
               cvp_snowprd(k) = cvp_snowprd(k) + cvp_snowprd_msfc(k,msfc)
               a_p(k)         = a_p(k)         +         a_p_msfc(k,msfc)
               am_up(k)       = am_up(k)       +       am_up_msfc(k,msfc)
               am_us(k)       = am_us(k)       +       am_us_msfc(k,msfc)
               am_evp(k)      = am_evp(k)      +      am_evp_msfc(k,msfc)
               am_pu(k)       = am_pu(k)       +       am_pu_msfc(k,msfc)
               am_pd(k)       = am_pd(k)       +       am_pd_msfc(k,msfc)
               am_pr(k)       = am_pr(k)       +       am_pr_msfc(k,msfc)
               am_ps(k)       = am_ps(k)       +       am_ps_msfc(k,msfc)
               ! Nov.29.2012. Add tracer block.
               ! Dec.13.2012. Add the line of wet deposition within enironment. Note that wet deposition
               !              within convective updraft and downdraft have already been computed.
               do mt = 1, ncnst
                  nttrrsprd(k,mt)   = nttrrsprd(k,mt)   +   nttrrsprd_msfc(k,msfc,mt) 
                  evptrrs_e(k,mt)   = evptrrs_e(k,mt)   +   evptrrs_e_msfc(k,msfc,mt)
                  evptrrs_d(k,mt)   = evptrrs_d(k,mt)   +   evptrrs_d_msfc(k,msfc,mt)
                  wdeptrrs_e(k,mt)  = wdeptrrs_e(k,mt)  +  wdeptrrs_e_msfc(k,msfc,mt)
                  flxtrrs(k,mt)     = flxtrrs(k,mt)     +     flxtrrs_msfc(k,msfc,mt)
                  cvp_trrsprd(k,mt) = cvp_trrsprd(k,mt) + cvp_trrsprd_msfc(k,msfc,mt)
               end do
            end do
            ! Since evaporation within downdraft has already been treated as a part of downdraft tendency ( qtten_d(k) ), 
            ! below treatment of '_evp(k)' should only contain the processes occuring within the environment.
            qlten_evp(k)  = 0._r8
            qiten_evp(k)  = 0._r8
            qvten_evp(k)  = evprain_e(k) - cvp_rainprd(k) + evpsnow_e(k) - cvp_snowprd(k)
            qtten_evp(k)  = qlten_evp(k) + qiten_evp(k) + qvten_evp(k)
            slten_evp(k)  = -xlv*(evprain_e(k) - cvp_rainprd(k)) - xls*(evpsnow_e(k) - cvp_snowprd(k)) - (xls - xlv)*snowmlt_e(k)
            sten_evp(k)   = slten_evp(k) + xlv * qlten_evp(k) + xls * qiten_evp(k)
            uten_evp(k)   = 0._r8
            vten_evp(k)   = 0._r8
            ! ----------------------------------------------------------------------------------------- !
            ! TRACERS REFINEMENT NECESSARY : EVAPORATION OF CONVECTIVE PRECIPITATION WITHIN ENVIRONMENT !
            ! Nov.29.2012. Below is updated and correctly computed.                                     !
            ! Dec.13.2012. I should carefully check whether below is double-counting with the separate  !
            !              routine of wet deposition of aerosols.                                       !
            ! Dec.13.2012. Add wet deposition part.                                                     !
            !              Note that 'evptrrs > 0' and 'wdeptrrs > 0'. Note also that wet deposition    !
            !              component is only applied to tracers not to the other thermodynamic scalars. !
            !              For convenience and following previous code, 'cvp_trrsprd' is added into     !
            !              the 'trten_evp'.                                                             !   
            !              Note that (-) sign should be multiplied in front of 'wdeptrrs'.              !
            ! Dec.13.2012. Note that 'cvp_trrsprd > 0' contains corrective tendencies associated with   !
            !              evaporation both within downdraft and environment, so that it is a function  !
            !              of 'evptrrs'. This re-geration of aerosols by evaporation of precipitation   !
            !              i.e., 'evptrrs' is to some degree already treated in the separate wet        !
            !              deposition routine. Thus, in order to precent double counting in the current !
            !              CAM5 structure, I should in principle remove 'evptrrs' in the below formula  !
            !              of 'trten_evp'. If then, however, 'trten_evp << 0', so that resulting        !
            !              concentration of aerosols after convection may become negative. Of course,   !
            !              by using 'positive_tracer' subroutine, negative tracer will be eventually    !
            !              converted into positive tracer. So, in order to reduce global AOD,           !
            !              let's remove 'evptrrs' in the below computation of trten_evp.                !   
            ! Feb.05.2013. Note that below 'trten_evp(k,mt),trten_wdep(k,mt)' are tendencies only by    !
            !              the processes within environment, since corresponding tendencies occuring    !
            !              within convective updraft and downdraft have already been treated and saved  !
            !              into 'trten_dia_u(k,mt)' and 'trten_dia_d(k,mt)'.                            !
            !              All final tendencies of tracers associated with diabatic forcings within     !
            !              environment ( trten_evp, trten_wdep ) are totally handled in the below block.!
            !              Thus, for correct computation of tracer tendency ( trten not trrsten ) by    !
            !              evaporation and wet deposition within environment, I only need to modify     !
            !              below block, which is sufficient and necessary.                              !    
            ! ----------------------------------------------------------------------------------------- !
            do mt = 1, ncnst 
               if( mt .eq. 1 .or. mt .eq. ixcldliq .or. mt .eq. ixcldice .or. mt .eq. ixnumliq .or. mt .eq. ixnumice ) then
                  trten_evp(k,mt)  = 0._r8
                  trten_wdep(k,mt) = 0._r8
               else
                  !?           ! In principle, I should use below two lines with appropriate computations of evptrrs_e(k,mt)
                  !?           ! and wdeptrrs_e(k,mt) in the main program. However, in order to avoid double counting with
                  !?           ! the similar treatment in the separate wet deposition routine in CAM5, I temporary set
                  !?           ! these two tendencies to zero. In future, if I turn-off wet deposition treatment by 
                  !?           ! convective precipitation in the separate routine in CAM5 ( i.e., by setting input cumulus
                  !?           ! area and convective precipitation to be zero ), I can use my robust below formula.    
                  ! trten_evp(k,mt)  =   - evptrrs_e(k,mt) - cvp_trrsprd(k,mt)
                  ! trten_wdep(k,mt) =  - wdeptrrs_e(k,mt)             
                  trten_evp(k,mt)  = - cvp_trrsprd(k,mt)
                  trten_wdep(k,mt) =   0._r8
               endif
            enddo
            ! ----------------------------------------------------------------------------------------- !
            ! TRACERS REFINEMENT NECESSARY : EVAPORATION OF CONVECTIVE PRECIPITATION WITHIN ENVIRONMENT !
            ! ----------------------------------------------------------------------------------------- !
         end do

         ! ----------------------------- !
         ! Precipitation flux at surface !
         ! Nov.29.2012. Add tracer block !
         ! ----------------------------- !
         do msfc = 1, nseg
            flxrain(0)    = flxrain(0) + flxrain_msfc(0,msfc)
            flxsnow(0)    = flxsnow(0) + flxsnow_msfc(0,msfc)
            ! Nov.29.2012. Add tracer block. 
            do mt = 1, ncnst 
               flxtrrs(0,mt) = flxtrrs(0,mt) + flxtrrs_msfc(0,msfc,mt)
            end do
            a_p(0)        = a_p(0)     +     a_p_msfc(0,msfc)
         end do
         precip  = ( flxrain(0) + flxsnow(0) ) / 1000._r8
         snow    =   flxsnow(0) / 1000._r8       

         ! -------------------------------------------------------------------------------------------------------- !
         ! Vertically-integrated grid-mean differential evaporation rate of convective precipitation ( aw*(Qn-Qw) ) !
         ! averaged over the PBL depth. Note that this should only consider didbatic forcing within environment not !
         ! within convective downdraft since diabatic forcing within convective downdraft has already been computed !
         ! above.                                                                                                   !
         ! I should also include 'corrective flux' ( 'cev' ) and 'snow melting'.                                    ! 
         ! For simplicity, I will assume that 'corrective flux' is homogeneous across the grid and so does not      !
         ! contribute to the computation of aw*(Qn-Qw).                                                             ! 
         ! Snow melting is assumed to occur even in the convective updraft in contrast to the                       !
         ! evaporation of convective precipitation. In order to treat snow melting more rigorously, I shoud use the !
         ! precipitation area before treating evaporation within downdraft ( a_p_prevp(k) defined at the interface )! 
         ! since snow melting was treated before evaporation within downdraft.                                      !
         ! As a more rigorous choice, I will assume that 'evaporation area, a_evp' is either within                 !
         ! completely within 'non-wake' area or within 'wake' area whenever possible. This different geometrical    !
         ! structure might be roughly described by using the similar tilting parameter, 'beta2'.                    !
         ! Conceptually, it seems to be most reasonable to assume that whenever possible, evaporation area exists   !
         ! within wake area ( i.e., beta2 = 1 ).                                                                    !       
         ! Future works : Treatment of tracers should be refined in future. Note that snow melting does not affect  !
         !                tracer concentrations. So, I can use tmp2 not tmp3.                                       !
         ! Sep.09.2011. Below block is commented-out since it will be computed later below in a collectively way    !
         !              for the whole wake forcing computation.                                                     !
         ! -------------------------------------------------------------------------------------------------------- !

         ! --------------------------------------------------------------------------- !
         !                                                                             !    
         ! Compute Grid-Mean Tendencies without repartitioning and dissipation heating !
         !                                                                             !
         ! --------------------------------------------------------------------------- !

         ! ------------------------------------------------------------------- !
         ! Currently, no constraint is imposed on qvten(k), qlten(k), qiten(k).!
         ! But negative condensate will be treated in positive_moisture in     ! 
         ! the above subroutine.                                               !
         ! ------------------------------------------------------------------- ! 

         do k = 1, ktop

            slten_NUM(k) =   slten_u(k) +     slten_d(k) +   slten_evp(k)
            qtten_NUM(k) =   qtten_u(k) +     qtten_d(k) +   qtten_evp(k)                
            uten_NUM(k)  =    uten_u(k) +      uten_d(k) +    uten_evp(k)                              
            vten_NUM(k)  =    vten_u(k) +      vten_d(k) +    vten_evp(k)                              

            qlten_NUM(k) =   qlten_u(k) +     qlten_d(k) +   qlten_evp(k)  
            qiten_NUM(k) =   qiten_u(k) +     qiten_d(k) +   qiten_evp(k)
            do mt = 1, ncnst
               trten_NUM(k,mt) = trten_u(k,mt) + trten_d(k,mt) + trten_evp(k,mt) + trten_wdep(k,mt)    
            enddo

            qvten_NUM(k) = qtten_NUM(k) -       qlten_NUM(k) -       qiten_NUM(k)
            sten_NUM(k)  = slten_NUM(k) + xlv * qlten_NUM(k) + xls * qiten_NUM(k)

         enddo

         ! ---------------------------------------------------------------------------------------- !
         ! Repartition the tendency in the lowest model layer                                       !
         ! into all the layers within the PBL or top of convection or whole atmospheric layer.      !
         ! It seems to be most reasonable to distribute to the whole convection layer.              !
         ! Note that the partition is done only for surface updraft flux not downdraft flux.        !
         ! Sep.12.2011. I double checked that my below 'ipartition' formula is perfect because      !
         !              (1) it conserves column-integrated energy, and (2) it completely remove     !
         !              the generation of unreasonable convective tendency in the lowest model      !
         !              layer by convection, so that it grauantees computation of reasonable        ! 
         !              surface heat, moisture, momentum, and tracer fluxes at surface in the       !  
         !              following surface flux computation routine in the CAM. Also, I don't        !
         !              need to modify any parts of CAM5 ( e.g., PBL scheme, surface flux           !
         !              routine ), since all the required modifications are contained in the        !
         !              UNICON in a fully reasonable way.                                           !
         !              By using 'ipartition = 1' option, I don't need to combine 'symmetric        !
         !              moist turbulence scheme' with the 'asymmetric moist turbulence scheme'      !
         !              within the implicit iteration loop, so that I can save tremendous amount    !    
         !              of computation time.                                                        !
         ! Oct.24.2011. For correct consistent output, I should also correct flux at each model     !
         !              interface, i.e., from each flux interface, I should subtract linear flux    !
         !              profile that is 'slflx_u(0)' at surface but 'zero' at k = ktop interface.   !
         ! Mar.19.2014. I should compute a similar non-zero partitioning tendency of 'ql,qi' below, !
         !              which is done on this day. In addition, I added the similar portion of      ! 
         !              convective downdraft fluxes in computing below partitioning tendencies.     !              
         !              Note that 'kpblhm >= 1' and 'iopt_partition = 2' might be more conceptually !
         !              consistent with the cold-pool formulation.                                  ! 
         ! ---------------------------------------------------------------------------------------- !


         if( iup_par .eq. 1 ) then             ! Lowest Layer - No Partitioning
             ktop_up_par = 1
         elseif( iup_par .eq. 2 ) then         ! Minimum of PBL Top and Cumulus Top
             ktop_up_par = min( kpblhm, ktop )
         elseif( iup_par .eq. 3 ) then         ! PBL Layers   
             ktop_up_par = kpblhm
         elseif( iup_par .eq. 4 ) then         ! Cumulus Layers
             ktop_up_par = ktop
         elseif( iup_par .eq. 5 ) then         ! Entire Layers
             ktop_up_par = mkx
         endif 

         if( idn_par .eq. 1 ) then             ! Lowest Layer - No Partitioning
             ktop_dn_par = 1
         elseif( idn_par .eq. 2 ) then         ! Minimum of PBL Top and Cumulus Top
             ktop_dn_par = min( kpblhm, ktop )
         elseif( idn_par .eq. 3 ) then         ! PBL Layers   
             ktop_dn_par = kpblhm
         elseif( idn_par .eq. 4 ) then         ! Cumulus Layers
             ktop_dn_par = ktop
         elseif( idn_par .eq. 5 ) then         ! Entire Layers
             ktop_dn_par = mkx
         endif 

         do k = ktop_up_par, 1, -1
            tmp2 = g / ( ps0(0) - ps0( ktop_up_par ) )
            slten_par(k) = - slflx_u(0) * tmp2
            qtten_par(k) = - qtflx_u(0) * tmp2
            uten_par(k)  = -  uflx_u(0) * tmp2
            vten_par(k)  = -  vflx_u(0) * tmp2
            qlten_par(k) = - qlflx_u(0) * tmp2
            qiten_par(k) = - qiflx_u(0) * tmp2
            do mt = 1, ncnst
               trten_par(k,mt) = - trflx_u(0,mt) * tmp2
            enddo
            ! -------------------------------------------------------------------------- !
            ! Oct.24.2011. Added below block for diagnostic 'ipartition' flux output.    !
            !              Note that below correction does not change simulation output. !
            ! Mar.20.2014. Note that fluxes at the surface are printed-out as non-zero   !
            !              purely for the diagnostic purpose. However, in the numerical  !
            !              computation they are treated to be zero with partitioning.    !  
            ! -------------------------------------------------------------------------- ! 
            tmp1         =   ( ps0(k) - ps0( ktop_up_par ) ) / ( ps0(0) - ps0( ktop_up_par ) )
            slflx_u(k)   =   slflx_u(k) - slflx_u(0) * tmp1
            qtflx_u(k)   =   qtflx_u(k) - qtflx_u(0) * tmp1
            uflx_u(k)    =    uflx_u(k) -  uflx_u(0) * tmp1
            vflx_u(k)    =    vflx_u(k) -  vflx_u(0) * tmp1
            qlflx_u(k)   =   qlflx_u(k) - qlflx_u(0) * tmp1
            qiflx_u(k)   =   qiflx_u(k) - qiflx_u(0) * tmp1
            do mt = 1, ncnst
               trflx_u(k,mt) = trflx_u(k,mt) - trflx_u(0,mt) * tmp1
            enddo
         enddo

         do k = ktop_dn_par, 1, -1
            tmp2 = g / ( ps0(0) - ps0( ktop_dn_par ) )
            slten_par(k) = slten_par(k) - slflx_d(0) * tmp2
            qtten_par(k) = qtten_par(k) - qtflx_d(0) * tmp2
            uten_par(k)  =  uten_par(k) -  uflx_d(0) * tmp2
            vten_par(k)  =  vten_par(k) -  vflx_d(0) * tmp2
            qlten_par(k) = qlten_par(k) - qlflx_d(0) * tmp2
            qiten_par(k) = qiten_par(k) - qiflx_d(0) * tmp2
            do mt = 1, ncnst
               trten_par(k,mt) = trten_par(k,mt) - trflx_d(0,mt) * tmp2
            enddo
            ! -------------------------------------------------------------------------- !
            ! Oct.24.2011. Added below block for diagnostic 'ipartition' flux output.    !
            !              Note that below correction does not change simulation output. !
            ! Mar.20.2014. Note that fluxes at the surface are printed-out as non-zero   !
            !              purely for the diagnostic purpose. However, in the numerical  !
            !              computation they are treated to be zero with partitioning.    !  
            ! -------------------------------------------------------------------------- ! 
            tmp1         =   ( ps0(k) - ps0( ktop_dn_par ) ) / ( ps0(0) - ps0( ktop_dn_par ) )
            slflx_d(k)   =   slflx_d(k) - slflx_d(0) * tmp1
            qtflx_d(k)   =   qtflx_d(k) - qtflx_d(0) * tmp1
            uflx_d(k)    =    uflx_d(k) -  uflx_d(0) * tmp1
            vflx_d(k)    =    vflx_d(k) -  vflx_d(0) * tmp1
            qlflx_d(k)   =   qlflx_d(k) - qlflx_d(0) * tmp1
            qiflx_d(k)   =   qiflx_d(k) - qiflx_d(0) * tmp1
            do mt = 1, ncnst
               trflx_d(k,mt) = trflx_d(k,mt) - trflx_d(0,mt) * tmp1
            enddo
         enddo

         do k = max( ktop_up_par, ktop_dn_par ), 1, -1
            slten_NUM(k) =  slten_NUM(k) +   slten_par(k)
            qtten_NUM(k) =  qtten_NUM(k) +   qtten_par(k)
            uten_NUM(k)  =   uten_NUM(k) +    uten_par(k)
            vten_NUM(k)  =   vten_NUM(k) +    vten_par(k)
            qlten_NUM(k) =  qlten_NUM(k) +   qlten_par(k)
            qiten_NUM(k) =  qiten_NUM(k) +   qiten_par(k)
            qvten_NUM(k) =  qtten_NUM(k) -       qlten_NUM(k) -       qiten_NUM(k)
            sten_NUM(k)  =  slten_NUM(k) + xlv * qlten_NUM(k) + xls * qiten_NUM(k)
            do mt = 1, ncnst
               trten_NUM(k,mt) = trten_NUM(k,mt) + trten_par(k,mt)
            enddo        
         enddo

         ! -------------------------------------------------------------------------- !
         ! Reset surface flux to be zero for column energy conservation.              !
         ! I commented-out below lines for looking at the detailed diagnostic output. !
         ! Note that this resetting does not influence numerical computation at all.  !
         ! Oct.24.2011. I restored below block for correct diagnostic output.         !
         ! Mar.19.2014. I added downdraft portition with the ipartition=1 condition   !
         !              and also 'ql,qi' portions.                                    !
         !              However, I commented-out below block on this day to see the   !
         !              actual non-zero convective fluxes at the surface.             !
         ! -------------------------------------------------------------------------- !
    !par if( ipartition .eq. 1 .and. kiss .eq. 0 ) then 
    !par     slflx_u(0) =  0._r8
    !par     qtflx_u(0) =  0._r8
    !par     uflx_u(0)  =  0._r8
    !par     vflx_u(0)  =  0._r8
    !par     qlflx_u(0) =  0._r8
    !par     qiflx_u(0) =  0._r8
    !par     do mt = 1, ncnst
    !par        trflx_u(0,mt) = 0._r8
    !par     enddo
    !par     slflx_d(0) =  0._r8
    !par     qtflx_d(0) =  0._r8
    !par     uflx_d(0)  =  0._r8
    !par     vflx_d(0)  =  0._r8
    !par     qlflx_d(0) =  0._r8
    !par     qiflx_d(0) =  0._r8
    !par     do mt = 1, ncnst
    !par        trflx_d(0,mt) = 0._r8
    !par     enddo
    !par endif

         ! ---------------------------------------------------------------------- !
         ! Choose either 'NUMerial' or 'ANAlytical' tendency as a final tendency. !
         ! ---------------------------------------------------------------------- !

         uten(:mkx)  =  uten_NUM(:mkx)
         vten(:mkx)  =  vten_NUM(:mkx)
         qlten(:mkx) = qlten_NUM(:mkx)
         qiten(:mkx) = qiten_NUM(:mkx)
         qvten(:mkx) = qvten_NUM(:mkx)
         sten(:mkx)  =  sten_NUM(:mkx)
         do mt = 1, ncnst
            trten(:mkx,mt) = trten_NUM(:mkx,mt)
         enddo

         ! ----------------------------------------------------------------------------------- !
         ! Save '_mxen' variables associated with multiple mixing environmental airs.          !
         ! Aug.01.2011. Brian Juwon Park's 10th Birthday.                                      !
         !              The explicit ensemble mixing process 'iter' routine are included here. !
         !              I have a hunch that this is the remaining final process I should       !
         !              implement into the UNICON.                                             !
         ! ----------------------------------------------------------------------------------- !

         cmf_u_mxen(0:mkx,iter)                   =                   cmf_u(0:mkx)
         cmf_d_mxen(0:mkx,iter)                   =                   cmf_d(0:mkx)
         slflx_u_mxen(0:mkx,iter)                 =                 slflx_u(0:mkx)
         slflx_d_mxen(0:mkx,iter)                 =                 slflx_d(0:mkx)
         qtflx_u_mxen(0:mkx,iter)                 =                 qtflx_u(0:mkx)
         qtflx_d_mxen(0:mkx,iter)                 =                 qtflx_d(0:mkx)
         uflx_u_mxen(0:mkx,iter)                  =                  uflx_u(0:mkx)
         uflx_d_mxen(0:mkx,iter)                  =                  uflx_d(0:mkx)
         vflx_u_mxen(0:mkx,iter)                  =                  vflx_u(0:mkx)
         vflx_d_mxen(0:mkx,iter)                  =                  vflx_d(0:mkx)

         flxrain_u_mxen(0:mkx,iter)               =                 flxrain(0:mkx)
         flxsnow_u_mxen(0:mkx,iter)               =                 flxsnow(0:mkx)

         ! ------------------------------------------------------------------------------------------------------------------------------------------- !
         ! Aug.31.2011. Add downdraft flux at the PBL top interface associated with convective organization                                            !
         !              in order to parameterize density current and convective organization.                                                          !
         ! Sep.07.2011. Compute total organization densitiy current forcing of conservative scalars and thv                                            !
         !              for the difference between non-wake area and grid-mean.                                                                        ! 
         !              Include not only adiabatic forcing but also diabatic forcing both within convective                                            !
         !              downdraft and environment. Also change the variable name from 'thlflx_d_org_pblh_mxen'                                         !
         !              to 'thl_orgforce_mxen' to denote that it includes all adiabatic and diabatic forcings.                                         !
         !              CAREFUL CONSIDERATION : Also impose a constraint based on the sign of total buoyancy forcing, if necessary.                    !
         !              Below is overall forcing averaged over the PBL.                                                                                !    
         !              The sum of below 3 should be positive in order to generate negative 'thv' anomaly within wake area or                          !
         !              equivalently, positive 'thv' anomaly within the non-wake area.                                                                 !
         ! 1. Adiabatic forcing                             : thvflx_d_org_pblh * g / pblhp                                                 [ K / s ]  !
         ! 2. Diabatic  forcing within convective downdraft : thl_dia_d_org * ( 1._r8 + zvir * qt0PBL ) + zvir * thl0PBL * qt_dia_d_org     [ K / s ]  !
         ! 3. Diabatic  forcing within environment          : thl_dia_env_org * ( 1._r8 + zvir * qt0PBL ) + zvir * thl0PBL * qt_dia_env_org [ K / s ]  !
         !                                                                                                                                             !
         ! Sep.07.2011. I also included 'inverse tau' into the mxen array in the below block.                                                          !
         !              Here, 'cd_' is non-dimensional darg coefficient for each conservative scalar and 'ws' is wind speed in the lowest model layer. !
         !              Temporarily, these drag coefficients are set to zero but we can include this component later.                                  !  
         !              Note that 'awk_' should not have this drag component since vertical mass exchange does not occur between the surface and       !
         !              atmosphere.                                                                                                                    !
         !              ( Example ) Fs [J/s/m2] = rho * cp * cd_thl * ws1 * ( Ts - thl0(1) ). Thus, if we know Ts, we can back up 'cd_thl'.            !
         ! Sep.07.2011. Instead of setting to zero, let's use a certain characteristic value of cd = 1.5e-3 for all scalars.                           !
         !              The ;thv' is not prognosed anymore for full consistency.                                                                       !
         ! ------------------------------------------------------------------------------------------------------------------------------------------- !

         ! -------------------------------------------- !
         !                                              !
         ! Computation of all the wake-related forcings !
         !                                              ! 
         ! -------------------------------------------- ! 
   
         ! -------------------------------------------------------------------------------------------------------------------- !
         ! 2. Differential diabatic forcing between 'non-wake' and 'all over the grid'.                                         !                                      
         !    ONLY WITHIN convective updraft and downdraft ( ['a_nw*Q_nw'/a_nw - 'grid-mean Q']_only_within_updraft_downdraft ) !
         !    NOT  WITHIN environmental portions of 'non-wake' and 'all over the grid'.                                         !
         !    The 'und' denotes 'updraft and downdraft' in contrast to 'env' which denotes environment.                         ! 
         !    Note that we should not include 'momentum forcing' here since 'momentum C' is conversion not diabatic forcing.    !
         !    The resulting units are [ kg / kg / s ], [ K / s ], [ # / kg / s ] or [ kg / kg / s ].                            !
         !    Sep.10.2011. Below block is alternatively chosen for using selectively chosen convective downdraft instead of     !
         !                 all convective downdrafts.                                                                           !
         !    Sep.13.2011. Through rigorous methematical derivation, I double-checked that                                      !
         !                 Below is always valid regardless whether I choose all convective downdraft or selective downdraft.   !
         ! -------------------------------------------------------------------------------------------------------------------- !
     
         qt_dia_und_org           = 0._r8 
         thl_dia_und_org          = 0._r8
         tr_dia_und_org(1:ncnst)  = 0._r8  
         tmp2                     = awk_PBL / ( 1._r8 - awk_PBL )

         do k = 1, kpblhm
            qt_dia_und_org        =     qt_dia_und_org + dp0(k) * ( tmp2 *   qtten_dia_u(k) ) 
            thl_dia_und_org       =    thl_dia_und_org + dp0(k) * ( tmp2 *  thlten_dia_u(k) ) 
            do mt = 1, ncnst
               tr_dia_und_org(mt) = tr_dia_und_org(mt) + dp0(k) * ( tmp2 * ( trten_dia_u(k,mt) + trten_eff_u(k,mt) ) )
            enddo
         enddo
         qt_dia_und_org           =     qt_dia_und_org / pblhp  -     qt_dia_d_org
         thl_dia_und_org          =    thl_dia_und_org / pblhp  -    thl_dia_d_org
         do mt = 1, ncnst
            tr_dia_und_org(mt)    = tr_dia_und_org(mt) / pblhp  - tr_dia_d_org(mt)  
         enddo

         ! ------------------------------------------------------------------------------------------------------------- !
         ! 2. Differential diabatic forcing between 'non-wake' and 'all over the grid'                                   !               
         !    ONLY WITHIN environmental portions of 'non-wake' and 'all over the grid'.                                  !
         !    This is caused by two processes : (1) evaporation of precipitation within environment,                     !
         !                                      (2) snow melting                 within environment.                     !
         !    The corrective tendendies 'cev' is inevitable assumed to occur uniformly all over the grid                 ! 
         !    and so does not contribute here.                                                                           !
         !    Note als that I am assuming dissipative heating is uniformly distributed all over the grid                 !
         !    and so does not contribute here.                                                                           !
         !    The resulting units are [ kg / kg / s ], [ K / s ], [ # / kg / s ] or [ kg / kg / s ].                     !
         !    Below 'tmp2' is from the the assumption of homogeneous distribution.                                       !
         !    tmp2 = - ( awk_PBL / ( 1._r8 - awk_PBL ) ) * ( am_u(k) / ( 1._r8 - am_u(k) )                               !
         !    Here, 'a_evp_wk' is overlapping area between 'evaporation area' and 'wake area'. Similar to the            !
         !    overlapping treatment between precipitation area and updraft fractional area, I am using additional        !
         !    tilting parameter to control the overlap between evaporation area ( am_evp(k) ) and wake area ( awk_PBL ). !
         !    Note that 'beta1' is likely to be related to 'beta2'.                                                      !
         !    Rigorously speaking, since below formula is explicitly using a_p_prevp(k), below treatment of snow melting !
         !    is only valid if 'i_snowmlt = 0' and so all snow melting occurs before evaporation within downdraft (      !  
         !    thus 'snowmlt(k) = 0' but 'snowmlt_e(k) > 0 ' )                                                            !
         !    However, since I will choose 'i_snowmlt = 0' as a default forever, below treatment is completely correct.  !
         !    Sep.09.2011. I checked couple of times that below formula is completely correct.                           !
         !                 Note that snow melting does not influence the tracers.                                        ! 
         !    Oct.03.2011. Since 'am_p_wk' and 'tmp3' are used for partitioning snow melting and since snow melting was  !
         !                 computed before evaporation of precipitation within downdraft and environment, I should use   !
         !                 the 'a_p_prevp' not the 'a_p' in computing 'am_p_wk' and 'tmp3' below. Thus, my below code    !
         !                 is completely correct.                                                                        !                            
         ! ------------------------------------------------------------------------------------------------------------- !

         qt_dia_env_org           = 0._r8 
         thl_dia_env_org          = 0._r8
         tr_dia_env_org(1:ncnst)  = 0._r8  
         do k = 1, kpblhm
            am_evp_nw  = ( 1._r8 - beta2 ) * am_evp(k)      * ( 1._r8 - awk_PBL ) + beta2 * max( am_evp(k)    - awk_PBL, 0._r8 )
            am_p_nw    = am_pu(k) + (1._r8 - beta2)*(a_p(k) - am_pu(k))*(1._r8 - awk_PBL) + &
                                           beta2*max( a_p(k) - am_pu(k) - awk_PBL, 0._r8 )
            tmp2                  =  ( 1._r8 / ( 1._r8 - awk_PBL ) ) * ( am_evp_nw - am_evp(k) * ( 1._r8 - awk_PBL ) )
            tmp3                  =  ( 1._r8 / ( 1._r8 - awk_PBL ) ) * ( am_p_nw   - a_p(k)    * ( 1._r8 - awk_PBL ) )     
            qt_dia_env_org  =    qt_dia_env_org +  tmp2*( evprain_e(k) + evpsnow_e(k) ) / max( am_evp(k), nonzero ) * dp0(k)
            thl_dia_env_org =   thl_dia_env_org + (tmp2*(-xlv*evprain_e(k) - xls*evpsnow_e(k))/max(am_evp(k),nonzero) + & 
                                                   tmp3*(-(xls-xlv)*snowmlt_e(k))/max(a_p(k),nonzero))/(cp*exn0(k))*dp0(k)

            am_evp_nw_st = (1._r8 - beta2_st)*am_evp_st(k)*(1._r8 - awk_PBL) + beta2_st*max( am_evp_st(k) - awk_PBL, 0._r8 )
            tmp2_st      =  ( 1._r8 / ( 1._r8 - awk_PBL ) ) * ( am_evp_nw_st - am_evp_st(k) * ( 1._r8 - awk_PBL ) )
            qt_dia_env_org  = qt_dia_env_org + tmp2_st*(evprain_st(k) + evpsnow_st(k)) / max( am_evp_st(k), nonzero ) * dp0(k)
            thl_dia_env_org = thl_dia_env_org + ( tmp2_st * ( - xlv * evprain_st(k) - xls * evpsnow_st(k) ) / &
                                                  max( am_evp_st(k), nonzero ) ) / ( cp * exn0(k) ) * dp0(k)
            do mt = 1, ncnst
               tr_dia_env_org(mt) = tr_dia_env_org(mt) +   tmp2 * ( trten_evp(k,mt) + trten_wdep(k,mt) ) * dp0(k)
            enddo
         enddo
         qt_dia_env_org           =     qt_dia_env_org / pblhp 
         thl_dia_env_org          =    thl_dia_env_org / pblhp  
         do mt = 1, ncnst
            tr_dia_env_org(mt)    = tr_dia_env_org(mt) / pblhp
         enddo

         ! ---------------------------------------------------------------------------------- !
         ! Computation of total wake forcing and relaxation time scale                        !
         ! Sep.09.2011. Note that I don't need to impose any limitation based on the the sign !
         !              of buoyancy flux or forcing since now UNICON can generally handle all !
         !              of the cases of positive and negative buoyancy forcings.              ! 
         !              Note that I don't need to prognose 'thv' anymore since it will be     !
         !              computed diagnostically from the prognosed 'thl,qt' at the beginning  !
         !              of nex time step. This will impose a full consistency into the model. !   
         ! Sep.16.2011. Note that 'cdrag' should be further reduced in principle since the    !
         !              entrainment flux at the PBL top in the 'a_D' area will be reduced due !
         !              to enhanced stratification at the PBL top within 'a_D'.               !
         !              This reduced entrainment effect should be simulated by the 'cdrag'    !
         !              alone in the current UNICON. Thus, we should use smaller 'cdrag'      !
         !              than the common value of 1.5e-3. The use of small value will also     !
         !              improve the diurnal cycle.                                            !
         !              In addition, in order to reduce model sensitivity to the vertical     !
         !              resolution, we may dfine 'ws1' using the PBL-averaged wind instead of !
         !              the value in the lowest model layer. This should be considered in     !
         !              future.                                                               !        
         ! Dec.20.2012. Bug fix. In computing inverse relaxation time scale of conservative   !
         !              scalars (taui below) other than 'taui_awk_mxen', I should include     !
         !              the following below term which has a unit of [1/s] :                  !
         !               ( g / pblhp ) * max( 0._r8, cmf_d(kpblhm) - cmf_d_org_pblh ) : [1/s] !
         !              This correction will help to reduce simulated convective organization !
         !              due to the increase of damping effect. This is fully physically and   !
         !              conceptually consistent.                                              !
         ! Dec.20.2012. In addition, I also add entrainment rate at the PBL top interface,    !
         !              went_eff [m/s] as well as 'cmf_d(kpblhm) - cmf_d_org_pblh' in         !
         !              computing 'tmp3' below. Note that the unit of rho0(kpblhm) * went_eff !
         !              is [ kg / m2 / s ], same as the unit of the mass flux.                !
         ! ---------------------------------------------------------------------------------- !

         ! ------------------------ !
         ! 1. Relaxation Time Scale !
         ! ------------------------ !   

         tmp1           = ( g / pblhp )
       ! May.1.2014. Correct budget of cold pool.
         if( i_budget_coldpool .eq. 0 ) then
             tmp3           = tmp1 * ( max( 0._r8, cmf_d(kpblhm) - cmf_d_org_pblh ) - cmf_u(kpblhm) + rho0(kpblhm) * went_eff )
         elseif( i_budget_coldpool .eq. 3 ) then
             tmp3           = tmp1 * ( max( 0._r8, cmf_d(kpblhm) - cmf_d_org_pblh ) + rho0(kpblhm) * went_eff - &
                                     ( ( 1._r8 - awk_PBL * cdelta_s ) / ( 1._r8 - awk_PBL ) ) * cmf_u(kpblhm) )
         elseif( i_budget_coldpool .eq. 6 ) then
             tmp3           = tmp1 * ( max( 0._r8, cmf_d(kpblhm) - cmf_d_org_pblh ) + rho0(kpblhm) * went_eff - &
                                     ( ( 1._r8 - awk_PBL * cdelta_s ) / ( 1._r8 - awk_PBL ) ) * cmf_u_org_pblh )
         elseif( i_budget_coldpool .eq. 1 .or. i_budget_coldpool .eq. 2 ) then
             tmp3           = tmp1 * ( max( 0._r8, cmf_d(kpblhm) - cmf_d_org_pblh - cmf_d_orgU_pblh ) + & 
                                     ( cmf_d_orgU_pblh - cmf_u_org_pblh ) / ( 1._r8 - awk_PBL ) + rho0(kpblhm) * went_eff )
         elseif( i_budget_coldpool .eq. 4 ) then 
             tmp3 = tmp1 * ( max( 0._r8, cmf_d_orgU_pblh / ( 1._r8 - awk_PBL ) ) - cmf_u(kpblhm) + rho0(kpblhm) * went_eff )
         elseif( i_budget_coldpool .eq. 5 ) then 
             tmp3           = tmp1 * ( cmf_d_orgU_pblh / ( 1._r8 - awk_PBL ) + rho0(kpblhm) * went_eff - &
                                     ( ( 1._r8 - awk_PBL * cdelta_s ) / ( 1._r8 - awk_PBL ) ) * cmf_u(kpblhm) )
         endif
         tmp4           = awk_PBL / ( 1._r8 - awk_PBL )
         ws1            = sqrt( u0(1)**2._r8 + v0(1)**2._r8 )
         cd_thl         = cdrag
         cd_qt          = cdrag
         cd_u           = cdrag
         cd_v           = cdrag
         cd_tr(:ncnst)  = cdrag

         if( int_del_wk .eq. 1 ) then
            del_wk_eff =  c_del_wk * tmp1 * awk_PBL * cmf_u(kpblhm)
         endif

         taui_thl_mxen(iter) = del_wk_eff / max( nonzero, awk_PBL * ( 1._r8 - awk_PBL ) ) +  cd_thl * tmp1 * ws1 * rho0(1) + tmp3
         taui_qt_mxen(iter)  = del_wk_eff / max( nonzero, awk_PBL * ( 1._r8 - awk_PBL ) ) +   cd_qt * tmp1 * ws1 * rho0(1) + tmp3
         taui_u_mxen(iter)   = del_wk_eff / max( nonzero, awk_PBL * ( 1._r8 - awk_PBL ) ) +    cd_u * tmp1 * ws1 * rho0(1) + tmp3
         taui_v_mxen(iter)   = del_wk_eff / max( nonzero, awk_PBL * ( 1._r8 - awk_PBL ) ) +    cd_v * tmp1 * ws1 * rho0(1) + tmp3

         do mt = 1, ncnst
            taui_tr_mxen(mt,iter) = del_wk_eff/max(nonzero, awk_PBL*(1._r8 - awk_PBL)) + cd_tr(mt) * tmp1 * ws1 * rho0(1) + tmp3
         enddo

       ! May.1.2014. Correct budget of cold pool.
         if( i_budget_coldpool .eq. 0 .or. i_budget_coldpool .eq. 3 ) then
             taui_awk_mxen(iter)                  = tmp1 * ( cmf_d_org_pblh - cmf_u(kpblhm) )
         elseif( i_budget_coldpool .eq. 6 ) then
             taui_awk_mxen(iter)                  = tmp1 * ( cmf_d_org_pblh - cmf_u_org_pblh )
         else
             taui_awk_mxen(iter)                  = tmp1 * ( cmf_d_org_pblh + cmf_d_orgU_pblh - cmf_u_org_pblh )
         endif 

         if( i_energy_coldpool .eq. 1 .or. i_energy_coldpool .eq. 2 ) then  
 
            taui_thl_mxen(iter) = del_wk0 / max( nonzero, ( 1._r8 - awk_PBL ) ) +     cd_thl * tmp1 * ws1 * rho0(1) + tmp3
            taui_qt_mxen(iter)  = del_wk0 / max( nonzero, ( 1._r8 - awk_PBL ) ) +      cd_qt * tmp1 * ws1 * rho0(1) + tmp3
            taui_u_mxen(iter)   = del_wk0 / max( nonzero, ( 1._r8 - awk_PBL ) ) +       cd_u * tmp1 * ws1 * rho0(1) + tmp3
            taui_v_mxen(iter)   = del_wk0 / max( nonzero, ( 1._r8 - awk_PBL ) ) +       cd_v * tmp1 * ws1 * rho0(1) + tmp3
            do mt = 1, ncnst
               taui_tr_mxen(mt,iter) = del_wk0 / max( nonzero, ( 1._r8 - awk_PBL ) ) +  cd_tr(mt) * tmp1 * ws1 * rho0(1) + tmp3
            enddo
 
          ! May.1.2014. Correct budget of cold pool.
            if( i_budget_coldpool .eq. 0 .or. i_budget_coldpool .eq. 3 ) then
                taui_awk_mxen(iter)                  = tmp1 * ( cmf_d_org_pblh - cmf_u(kpblhm) ) + ( del_wk0 - eps_wk0 )
            elseif( i_budget_coldpool .eq. 6 ) then
                taui_awk_mxen(iter)                  = tmp1 * ( cmf_d_org_pblh - cmf_u_org_pblh ) + ( del_wk0 - eps_wk0 )
            else
                taui_awk_mxen(iter) = tmp1 * ( cmf_d_org_pblh + cmf_d_orgU_pblh - cmf_u_org_pblh ) + ( del_wk0 - eps_wk0 )
            endif 

         endif
 
         del_org_mxen(iter)  = awk_PBL * ( 1._r8 - awk_PBL ) * ( g / pblhp ) * cmf_u(kpblhm)
         del0_org_mxen(iter) =           ( 1._r8 - awk_PBL ) * ( g / pblhp ) * cmf_u(kpblhm)
 
         ! --------------------- !
         ! 2. Total Wake Forcing !
         ! --------------------- !   

       ! May.1.2014. Correct budget of cold pool. Note that 'awk_force_mxen' is identical to the previous case.
 
         if( i_budget_coldpool .eq. 0 .or. i_budget_coldpool .eq. 3 ) then

             thl_orgforce_mxen(iter)                  =       thlflx_d_org_pblh * tmp1 +    thl_dia_und_org +    thl_dia_env_org
             qt_orgforce_mxen(iter)                   =        qtflx_d_org_pblh * tmp1 +     qt_dia_und_org +     qt_dia_env_org
             u_orgforce_mxen(iter)                    =         uflx_d_org_pblh * tmp1
             v_orgforce_mxen(iter)                    =         vflx_d_org_pblh * tmp1
             do mt = 1, ncnst
                tr_orgforce_mxen(mt,iter)             =    trflx_d_org_pblh(mt) * tmp1 + tr_dia_und_org(mt) +  tr_dia_env_org(mt)
             enddo

         elseif( i_budget_coldpool .eq. 1 .or. i_budget_coldpool .eq. 2 ) then

             thl_orgforce_mxen(iter)  =  (    thlflx_d_org_pblh - tmp4 * ( thlflx_u_org_pblh + thlflx_d_orgU_pblh ) ) * tmp1 + & 
                                                                thl_dia_und_org + tmp4 * thl_dia_d_orgU + thl_dia_env_org
             qt_orgforce_mxen(iter)   =  (     qtflx_d_org_pblh - tmp4 * (  qtflx_u_org_pblh +  qtflx_d_orgU_pblh ) ) * tmp1 + &  
                                                                 qt_dia_und_org + tmp4 *  qt_dia_d_orgU +  qt_dia_env_org
             u_orgforce_mxen(iter)    =  (      uflx_d_org_pblh - tmp4 * (   uflx_u_org_pblh +   uflx_d_orgU_pblh ) ) * tmp1
             v_orgforce_mxen(iter)    =  (      vflx_d_org_pblh - tmp4 * (   vflx_u_org_pblh +   vflx_d_orgU_pblh ) ) * tmp1
             do mt = 1, ncnst
                tr_orgforce_mxen(mt,iter)=(trflx_d_org_pblh(mt) - tmp4*( trflx_u_org_pblh(mt) + trflx_d_orgU_pblh(mt) ) ) * tmp1 +& 
                                                             tr_dia_und_org(mt) + tmp4 * tr_dia_d_orgU(mt) + tr_dia_env_org(mt)    
             enddo

         elseif( i_budget_coldpool .eq. 4 .or. i_budget_coldpool .eq. 5 ) then

             thl_orgforce_mxen(iter)     =  (    thlflx_d_org_pblh - tmp4 * (        0._r8 + thlflx_d_orgU_pblh ) ) * tmp1 + & 
                                                                thl_dia_und_org + tmp4 * thl_dia_d_orgU + thl_dia_env_org
             qt_orgforce_mxen(iter)      =  (     qtflx_d_org_pblh - tmp4 * (        0._r8 +  qtflx_d_orgU_pblh ) ) * tmp1 + &  
                                                                 qt_dia_und_org + tmp4 *  qt_dia_d_orgU +  qt_dia_env_org
             u_orgforce_mxen(iter)       =  (      uflx_d_org_pblh - tmp4 * (        0._r8 +   uflx_d_orgU_pblh ) ) * tmp1
             v_orgforce_mxen(iter)       =  (      vflx_d_org_pblh - tmp4 * (        0._r8 +   vflx_d_orgU_pblh ) ) * tmp1
             do mt = 1, ncnst
                tr_orgforce_mxen(mt,iter)= ( trflx_d_org_pblh(mt) - tmp4 * (        0._r8 + trflx_d_orgU_pblh(mt) ) ) * tmp1 + & 
                                                             tr_dia_und_org(mt) + tmp4 * tr_dia_d_orgU(mt) + tr_dia_env_org(mt)    
             enddo

         elseif( i_budget_coldpool .eq. 6 ) then

             thl_orgforce_mxen(iter) =  (thlflx_d_org_pblh - tmp4 * ( thlflx_u_org_pblh - cmf_u_org_pblh * cdelta_s * &
                                         delta_thl_PBL ) ) * tmp1 + thl_dia_und_org + thl_dia_env_org
             qt_orgforce_mxen(iter)  =  ( qtflx_d_org_pblh - tmp4 * (  qtflx_u_org_pblh - cmf_u_org_pblh * cdelta_s * &
                                         delta_qt_PBL ) ) * tmp1 + &
                                                             qt_dia_und_org +  qt_dia_env_org
             u_orgforce_mxen(iter)   =  (  uflx_d_org_pblh - tmp4 * (   uflx_u_org_pblh - cmf_u_org_pblh * cdelta_s * &
                                          delta_u_PBL ) ) * tmp1
             v_orgforce_mxen(iter)   =  (  vflx_d_org_pblh - tmp4 * (   vflx_u_org_pblh - cmf_u_org_pblh * cdelta_s * &
                                           delta_v_PBL ) ) * tmp1 
             do mt = 1, ncnst
                tr_orgforce_mxen(mt,iter) =  ( trflx_d_org_pblh(mt) - tmp4 * ( trflx_u_org_pblh(mt) - cmf_u_org_pblh * &
                                               cdelta_s * delta_tr_PBL(mt) ) ) * tmp1 + &  
                                                             tr_dia_und_org(mt) + tr_dia_env_org(mt)    
             enddo

         endif

         awk_orgforce_mxen(iter)                  =          cmf_d_org_pblh * tmp1 +  eps_wk_eff - del_wk_eff 
         
         if( i_energy_coldpool .eq. 1 .or. i_energy_coldpool .eq. 2 ) then  
 
            awk_orgforce_mxen(iter)              =          cmf_d_org_pblh * tmp1
 
         endif
 
         ! -------------------------------------------------- !
         ! 2-1. Individual Wake Forcing for Diagnostic Output !
         ! -------------------------------------------------- !   

         thl_orgforce_flx_mxen(iter)              =       thlflx_d_org_pblh * tmp1
         thl_orgforce_und_mxen(iter)              =                thl_dia_und_org
         thl_orgforce_env_mxen(iter)              =                thl_dia_env_org
         
         qt_orgforce_flx_mxen(iter)               =        qtflx_d_org_pblh * tmp1
         qt_orgforce_und_mxen(iter)               =                 qt_dia_und_org
         qt_orgforce_env_mxen(iter)               =                 qt_dia_env_org

         u_orgforce_flx_mxen(iter)                =         uflx_d_org_pblh * tmp1
         u_orgforce_und_mxen(iter)                =                          0._r8
         u_orgforce_env_mxen(iter)                =                          0._r8

         v_orgforce_flx_mxen(iter)                =         vflx_d_org_pblh * tmp1
         v_orgforce_und_mxen(iter)                =                          0._r8
         v_orgforce_env_mxen(iter)                =                          0._r8

         awk_orgforce_flx_mxen(iter)              =          cmf_d_org_pblh * tmp1
         awk_orgforce_mix_mxen(iter)              =        eps_wk_eff - del_wk_eff 
         
         cmf_d_org_pblh_mxen(iter)                =                 cmf_d_org_pblh

         ! ---------------------------------------------------------------------------------------------------- !
         ! Sep.07.2011. I should carefully consider whether I want to include below buoyancy constraint or not. !
         !              It seems to be more transparent and safe to include below buoyancy constraint block.    !
         !              Note that I still need to compute damping time scale for this case.                     !
         !              In case of 'taui_awk_mxen' and 'awk_orgforce_mxen', both of them are functions of       !
         !              cmf_d_org_pblh. Thus, when below if conditions happens, both the 'taui_awk_mxen' and    !
         !             'awk_orgforce_mxen' should be consistently modified. However, we already imposed the     !
         !              adiabatic forcing constraint on mass when buoyancy flux is negative, i.e,, I already    !
         !              imposed the constraint of cmf_d_org_pblh = 0 before. So, in principle, I don't need to  !
         !              do anything here. However, for safety, I also imposed the constraints on  awk as below. !
         ! Sep.09.2011. I don't need below if constraint any more sine UNICON is quite generally formulated to  !
         !              handle both positive/negative buoyancy forcing in wake.                                 ! 
         ! ---------------------------------------------------------------------------------------------------- !

         ! ------------------------- !
         ! End of Organization Block !
         ! ------------------------- !

         qvten_mxen(:mkx,iter)                    =                    qvten(:mkx)
         qlten_mxen(:mkx,iter)                    =                    qlten(:mkx)
         qiten_mxen(:mkx,iter)                    =                    qiten(:mkx)
         do mt = 1, ncnst
            trten_mxen(:mkx,mt,iter)              =                 trten(:mkx,mt)
         enddo
         sten_mxen(:mkx,iter)                     =                     sten(:mkx)
         uten_mxen(:mkx,iter)                     =                     uten(:mkx)
         vten_mxen(:mkx,iter)                     =                     vten(:mkx) 
         qrten_mxen(:mkx,iter)                    =                    qrten(:mkx) 
         qsten_mxen(:mkx,iter)                    =                    qsten(:mkx) 

         rqc_l_mxen(:mkx,iter)                    =                    rqc_l(:mkx)
         rqc_i_mxen(:mkx,iter)                    =                    rqc_i(:mkx)
         rqc_mxen(:mkx,iter)                      =                      rqc(:mkx)
         rnc_l_mxen(:mkx,iter)                    =                    rnc_l(:mkx)
         rnc_i_mxen(:mkx,iter)                    =                    rnc_i(:mkx)

         cmf_det_mxen(:mkx,iter)                  =                  cmf_det(:mkx)
         ql_det_mxen(:mkx,iter)                   =                   ql_det(:mkx)
         qi_det_mxen(:mkx,iter)                   =                   qi_det(:mkx)

         evapc_mxen(:mkx,iter)                    =                    evapc(:mkx)

         am_u_mxen(:mkx,iter)                     =                     am_u(:mkx)
         qlm_u_mxen(:mkx,iter)                    =                    qlm_u(:mkx)
         qim_u_mxen(:mkx,iter)                    =                    qim_u(:mkx)

         am_d_mxen(:mkx,iter)                     =                     am_d(:mkx)
         qlm_d_mxen(:mkx,iter)                    =                    qlm_d(:mkx)
         qim_d_mxen(:mkx,iter)                    =                    qim_d(:mkx)

         rliq_mxen(iter)                          =                           rliq
         precip_mxen(iter)                        =                         precip
         snow_mxen(iter)                          =                           snow

         cnt_mxen(iter)                           =                            cnt
         cnb_mxen(iter)                           =                            cnb

         slten_u_mxen(:mkx,iter)                  =                  slten_u(:mkx)
         qtten_u_mxen(:mkx,iter)                  =                  qtten_u(:mkx)
         uten_u_mxen(:mkx,iter)                   =                   uten_u(:mkx)
         vten_u_mxen(:mkx,iter)                   =                   vten_u(:mkx)
         sten_u_mxen(:mkx,iter)                   =                   sten_u(:mkx)
         qvten_u_mxen(:mkx,iter)                  =                  qvten_u(:mkx)
         qlten_u_mxen(:mkx,iter)                  =                  qlten_u(:mkx)
         qiten_u_mxen(:mkx,iter)                  =                  qiten_u(:mkx)
         do mt = 1, ncnst
            trten_u_mxen(:mkx,mt,iter)            =               trten_u(:mkx,mt)
         enddo

         slten_d_mxen(:mkx,iter)                  =                  slten_d(:mkx)
         qtten_d_mxen(:mkx,iter)                  =                  qtten_d(:mkx)
         uten_d_mxen(:mkx,iter)                   =                   uten_d(:mkx)
         vten_d_mxen(:mkx,iter)                   =                   vten_d(:mkx)
         sten_d_mxen(:mkx,iter)                   =                   sten_d(:mkx)
         qvten_d_mxen(:mkx,iter)                  =                  qvten_d(:mkx)
         qlten_d_mxen(:mkx,iter)                  =                  qlten_d(:mkx)
         qiten_d_mxen(:mkx,iter)                  =                  qiten_d(:mkx)
         do mt = 1, ncnst
            trten_d_mxen(:mkx,mt,iter)            =               trten_d(:mkx,mt)
         enddo

         slten_evp_mxen(:mkx,iter)                =                slten_evp(:mkx)
         qtten_evp_mxen(:mkx,iter)                =                qtten_evp(:mkx)
         uten_evp_mxen(:mkx,iter)                 =                 uten_evp(:mkx)
         vten_evp_mxen(:mkx,iter)                 =                 vten_evp(:mkx)
         sten_evp_mxen(:mkx,iter)                 =                 sten_evp(:mkx)
         qvten_evp_mxen(:mkx,iter)                =                qvten_evp(:mkx)
         qlten_evp_mxen(:mkx,iter)                =                qlten_evp(:mkx)
         qiten_evp_mxen(:mkx,iter)                =                qiten_evp(:mkx)
         do mt = 1, ncnst
            trten_evp_mxen(:mkx,mt,iter)          =             trten_evp(:mkx,mt)
         enddo

         qlten_sub_mxen(:mkx,iter)                =                qlten_sub(:mkx)
         qiten_sub_mxen(:mkx,iter)                =                qiten_sub(:mkx)

       ! Apr.15.2014. Temporary Hack

       ! qlten_sub_mxen(:mkx,iter)                =  qrten_u(:mkx) + qsten_u(:mkx)
       ! qiten_sub_mxen(:mkx,iter)                =  qrten_d(:mkx) + qsten_d(:mkx)

         qlten_det_mxen(:mkx,iter)                =                qlten_det(:mkx)
         qiten_det_mxen(:mkx,iter)                =                qiten_det(:mkx)

         thl_u_mxen(0:mkx,iter)                   =                   thl_u(0:mkx)
         qt_u_mxen(0:mkx,iter)                    =                    qt_u(0:mkx)
         u_u_mxen(0:mkx,iter)                     =                     u_u(0:mkx)
         v_u_mxen(0:mkx,iter)                     =                     v_u(0:mkx)
         w_u_mxen(0:mkx,iter)                     =                     w_u(0:mkx)
         ql_u_mxen(0:mkx,iter)                    =                    ql_u(0:mkx)
         qi_u_mxen(0:mkx,iter)                    =                    qi_u(0:mkx)
         do mt = 1, ncnst
            tr_u_mxen(0:mkx,mt,iter)              =                 tr_u(0:mkx,mt)
         enddo
         a_u_mxen(0:mkx,iter)                     =                     a_u(0:mkx)
         num_u_mxen(0:mkx,iter)                   =                   num_u(0:mkx)
         wa_u_mxen(0:mkx,iter)                    =                    wa_u(0:mkx)
         qla_u_mxen(0:mkx,iter)                   =                   qla_u(0:mkx)
         qia_u_mxen(0:mkx,iter)                   =                   qia_u(0:mkx)
         rad_u_mxen(0:mkx,iter)                   =                   rad_u(0:mkx)
         thva_u_mxen(0:mkx,iter)                  =                  thva_u(0:mkx)

         a_p_mxen(0:mkx,iter)                     =                     a_p(0:mkx)
         am_evp_mxen(:mkx,iter)                   =                   am_evp(:mkx)
         am_pu_mxen(:mkx,iter)                    =                    am_pu(:mkx)
         x_p_mxen(0:mkx,iter)                     =              x_p_msfc(0:mkx,1)
         y_p_mxen(0:mkx,iter)                     =              y_p_msfc(0:mkx,1)
         x_um_mxen(:mkx,iter)                     =              x_um_msfc(:mkx,1)
         y_um_mxen(:mkx,iter)                     =              y_um_msfc(:mkx,1)

         thl_d_mxen(0:mkx,iter)                   =                   thl_d(0:mkx)
         qt_d_mxen(0:mkx,iter)                    =                    qt_d(0:mkx)
         u_d_mxen(0:mkx,iter)                     =                     u_d(0:mkx)
         v_d_mxen(0:mkx,iter)                     =                     v_d(0:mkx)
         w_d_mxen(0:mkx,iter)                     =                     w_d(0:mkx)
         ql_d_mxen(0:mkx,iter)                    =                    ql_d(0:mkx)
         qi_d_mxen(0:mkx,iter)                    =                    qi_d(0:mkx)
         do mt = 1, ncnst
            tr_d_mxen(0:mkx,mt,iter)              =                 tr_d(0:mkx,mt)
         enddo
         a_d_mxen(0:mkx,iter)                     =                     a_d(0:mkx)
         wa_d_mxen(0:mkx,iter)                    =                    wa_d(0:mkx)
         qla_d_mxen(0:mkx,iter)                   =                   qla_d(0:mkx)
         qia_d_mxen(0:mkx,iter)                   =                   qia_d(0:mkx)

         thl_u_msfc_mxen(0:mkx,:nseg,iter)        =        thl_u_msfc(0:mkx,:nseg)
         qt_u_msfc_mxen(0:mkx,:nseg,iter)         =         qt_u_msfc(0:mkx,:nseg)
         u_u_msfc_mxen(0:mkx,:nseg,iter)          =          u_u_msfc(0:mkx,:nseg)
         v_u_msfc_mxen(0:mkx,:nseg,iter)          =          v_u_msfc(0:mkx,:nseg)
         w_u_msfc_mxen(0:mkx,:nseg,iter)          =          w_u_msfc(0:mkx,:nseg)
         ql_u_msfc_mxen(0:mkx,:nseg,iter)         =         ql_u_msfc(0:mkx,:nseg)
         qi_u_msfc_mxen(0:mkx,:nseg,iter)         =         qi_u_msfc(0:mkx,:nseg)
         do mt = 1, ncnst
            tr_u_msfc_mxen(0:mkx,:nseg,mt,iter)   =      tr_u_msfc(0:mkx,:nseg,mt) 
         enddo
         cmf_u_msfc_mxen(0:mkx,:nseg,iter)        =        cmf_u_msfc(0:mkx,:nseg)
         a_u_msfc_mxen(0:mkx,:nseg,iter)          =          a_u_msfc(0:mkx,:nseg)
         num_u_msfc_mxen(0:mkx,:nseg,iter)        =        num_u_msfc(0:mkx,:nseg)
         rad_u_msfc_mxen(0:mkx,:nseg,iter)        =        rad_u_msfc(0:mkx,:nseg)

         eps0_u_msfc_mxen(0:mkx,:nseg,iter)       =       eps0_u_msfc(0:mkx,:nseg)
         eps_u_msfc_mxen(0:mkx,:nseg,iter)        =        eps_u_msfc(0:mkx,:nseg)
         del_u_msfc_mxen(0:mkx,:nseg,iter)        =        del_u_msfc(0:mkx,:nseg)
         eeps_u_msfc_mxen(0:mkx,:nseg,iter)       =       eeps_u_msfc(0:mkx,:nseg)
         ddel_u_msfc_mxen(0:mkx,:nseg,iter)       =       ddel_u_msfc(0:mkx,:nseg)
         xc_u_msfc_mxen(0:mkx,:nseg,iter)         =         xc_u_msfc(0:mkx,:nseg)
         xs_u_msfc_mxen(0:mkx,:nseg,iter)         =         xs_u_msfc(0:mkx,:nseg)
         xemin_u_msfc_mxen(0:mkx,:nseg,iter)      =      xemin_u_msfc(0:mkx,:nseg)
         xemax_u_msfc_mxen(0:mkx,:nseg,iter)      =      xemax_u_msfc(0:mkx,:nseg)
         cridis_u_msfc_mxen(0:mkx,:nseg,iter)     =     cridis_u_msfc(0:mkx,:nseg)
         thvcuenv_u_msfc_mxen(0:mkx,:nseg,iter)   =   thvcuenv_u_msfc(0:mkx,:nseg)
         thvegenv_u_msfc_mxen(0:mkx,:nseg,iter)   =   thvegenv_u_msfc(0:mkx,:nseg)
         thvxsenv_u_msfc_mxen(0:mkx,:nseg,iter)   =   thvxsenv_u_msfc(0:mkx,:nseg)
         fmix_u_msfc_mxen(0:mkx,:nseg,iter)       =       fmix_u_msfc(0:mkx,:nseg)
         cmfumix_u_msfc_mxen(0:mkx,:nseg,iter)    =    cmfumix_u_msfc(0:mkx,:nseg)

         thl_d_msfc_mxen(0:mkx,:nseg,iter)        =        thl_d_msfc(0:mkx,:nseg)
         qt_d_msfc_mxen(0:mkx,:nseg,iter)         =         qt_d_msfc(0:mkx,:nseg)
         u_d_msfc_mxen(0:mkx,:nseg,iter)          =          u_d_msfc(0:mkx,:nseg)
         v_d_msfc_mxen(0:mkx,:nseg,iter)          =          v_d_msfc(0:mkx,:nseg)
         w_d_msfc_mxen(0:mkx,:nseg,iter)          =          w_d_msfc(0:mkx,:nseg)
         ql_d_msfc_mxen(0:mkx,:nseg,iter)         =         ql_d_msfc(0:mkx,:nseg)
         qi_d_msfc_mxen(0:mkx,:nseg,iter)         =         qi_d_msfc(0:mkx,:nseg)
         do mt = 1, ncnst
            tr_d_msfc_mxen(0:mkx,:nseg,mt,iter)   =      tr_d_msfc(0:mkx,:nseg,mt) 
         enddo
         cmf_d_msfc_mxen(0:mkx,:nseg,iter)        =        cmf_d_msfc(0:mkx,:nseg)
         a_d_msfc_mxen(0:mkx,:nseg,iter)          =          a_d_msfc(0:mkx,:nseg)
         wa_d_msfc_mxen(0:mkx,:nseg,iter)         =         wa_d_msfc(0:mkx,:nseg)
         qla_d_msfc_mxen(0:mkx,:nseg,iter)        =        qla_d_msfc(0:mkx,:nseg)
         qia_d_msfc_mxen(0:mkx,:nseg,iter)        =        qia_d_msfc(0:mkx,:nseg)

         ptop_msfc_mxen(:nseg,iter)               =               ptop_msfc(:nseg)
         ztop_msfc_mxen(:nseg,iter)               =               ztop_msfc(:nseg)

         ! ----------------------------------------------------------------------------------- !
         ! End of saving '_mxen' variables associated with multiple mixing environmental airs. !
         ! ----------------------------------------------------------------------------------- !

      enddo ! End of iter = 1, niter. This is an iteration loop of whole vertical layer

      ! -------------------------------------------------------------------------------- !
      !                                                                                  !
      !          Print-Out Formal Output Variables other than 'inout' Variables          !
      !                                                                                  ! 
      ! Aug.01.2011. Brian Juwon Park's 10th Birthday.                                   !
      !              The average of explicit ensemble mixing process                     !
      !              are treated here.                                                   !
      !              Below formula only considers two-types of mixing                    !
      !              with iter = 1, 2 but multi-types of mixing can                      !
      !              be treated in future.                                               !
      !                 (1) iter = 1 : with mean environmental airs at the current time  !
      !                 (2) iter = 2 : with cumulus updraft + detrained airs at the      !
      !                                previous time step.                               !  
      !              Note that I should use 'cuorg' not 'cuorg_mxen' in the below lines. !
      !                                                                                  !  
      ! -------------------------------------------------------------------------------- !

      ! ------------------------- !
      ! 1. Flux at the Interfaces !
      ! ------------------------- !

      do ki = 0, max( ktop_mxen(ixi) - 1, ktop_mxen(ixf) - 1 )
         kvi = mkx - ki
         cmf_u_out(i,ki) = (1._r8 - cuorg)*   cmf_u_mxen(ki,ixi)                         + cuorg*cmf_u_mxen(ki,ixf)   
         slflx_out(i,ki) = (1._r8 - cuorg)*(slflx_u_mxen(ki,ixi) + slflx_d_mxen(ki,ixi)) + cuorg*(slflx_u_mxen(ki,ixf) + &
                           slflx_d_mxen(ki,ixf) ) 
         qtflx_out(i,ki) = (1._r8 - cuorg)*(qtflx_u_mxen(ki,ixi) + qtflx_d_mxen(ki,ixi)) + cuorg*(qtflx_u_mxen(ki,ixf) + &
                           qtflx_d_mxen(ki,ixf) ) 
      enddo
 
      ! ------------------------ !
      ! 2. Layer-Mean Tendencies !
      ! ------------------------ !

      do k = 1, max( ktop_mxen(ixi), ktop_mxen(ixf) )
         kv = mkx + 1 - k
         qvten_out(i,k)          = ( 1._r8 - cuorg ) *     qvten_mxen(k,ixi) + cuorg *     qvten_mxen(k,ixf)
         qlten_out(i,k)          = ( 1._r8 - cuorg ) *     qlten_mxen(k,ixi) + cuorg *     qlten_mxen(k,ixf)
         qiten_out(i,k)          = ( 1._r8 - cuorg ) *     qiten_mxen(k,ixi) + cuorg *     qiten_mxen(k,ixf)
         do mt = 1, ncnst
            trten_out(i,k,mt)    = ( 1._r8 - cuorg ) *  trten_mxen(k,mt,ixi) + cuorg *  trten_mxen(k,mt,ixf)
         enddo
         sten_out(i,k)           = ( 1._r8 - cuorg ) *      sten_mxen(k,ixi) + cuorg *      sten_mxen(k,ixf)
         uten_out(i,k)           = ( 1._r8 - cuorg ) *      uten_mxen(k,ixi) + cuorg *      uten_mxen(k,ixf)  
         vten_out(i,k)           = ( 1._r8 - cuorg ) *      vten_mxen(k,ixi) + cuorg *      vten_mxen(k,ixf)
         qrten_out(i,k)          = ( 1._r8 - cuorg ) *     qrten_mxen(k,ixi) + cuorg *     qrten_mxen(k,ixf)
         qsten_out(i,k)          = ( 1._r8 - cuorg ) *     qsten_mxen(k,ixi) + cuorg *     qsten_mxen(k,ixf)
      enddo

      ! --------------------------------------------------------------------------------------------------- !
      ! SPECIAL : Compute Dissipation Heating - This must be done here to prevent energy conservation error !
      ! --------------------------------------------------------------------------------------------------- !

      ! --------------------------------------------------------------------------- !
      ! Compute diabatic tendency associated with KE dissipative heating            !
      ! In contrast to local symmetric turbulence scheme, KE dissipative heating by !
      ! nonlocal asymmetric turbulence (i.e., convection scheme)                    !
      ! can be either positive or negative.                                         !
      ! In order to suppress energy conservation error, this should be done here at !
      ! the end after choosing 'NUM' or 'ANA'.                                      ! 
      ! --------------------------------------------------------------------------- !
      do k = 1, max( ktop_mxen(ixi), ktop_mxen(ixf) )
         uf(k) = u0(k) + uten_out(i,k) * dt
         vf(k) = v0(k) + vten_out(i,k) * dt
      enddo
      ! ------------------------------------------------------------------------- !
      ! Reconstruct momemtum flux from the reconstructed uten_u(k) and uten_d(k)  !
      ! and vten_u(k) and vten_d(k) associated with ipartition = 1.               !
      ! Sep.12.2011. I should re-check whether below formula corretly incorporate !
      !              the 'ipartition = 1' effect. However, probably below is      !
      !              correct since it seems that now I can trust myself with      !
      !              reasonable amount of confidence.                             !  
      ! ------------------------------------------------------------------------- !
      do k = 1, max( ktop_mxen(ixi), ktop_mxen(ixf) )
         km = k - 1
         uflx(k) = uflx(km) - uten_out(i,k) * ( dp0(k) / g )
         vflx(k) = vflx(km) - vten_out(i,k) * ( dp0(k) / g )
      end do
      ! ----------------------------------------------------- !
      ! Add dissipation heating to the final heating tendency !
      ! ----------------------------------------------------- !
      do k = 1, max( ktop_mxen(ixi), ktop_mxen(ixf) )
         kp = k + 1
         km = k - 1 
         if( k .eq. 1 ) then
            sten_dis(k) = - g / 4._r8 * ( &
                            uflx(k)   * ( uf(kp) - uf(k) + u0(kp) - u0(k) ) / dps0(k) + & 
                            vflx(k)   * ( vf(kp) - vf(k) + v0(kp) - v0(k) ) / dps0(k) )
         elseif( k .ge. 2 .and. k .le. max( ktop_mxen(ixi), ktop_mxen(ixf) ) - 1 ) then
            sten_dis(k) = - g / 4._r8 * ( &
                            uflx(k)   * ( uf(kp) -  uf(k) + u0(kp) -  u0(k) ) /  dps0(k) + &
                            uflx(km)  * (  uf(k) - uf(km) +  u0(k) - u0(km) ) / dps0(km) + &
                            vflx(k)   * ( vf(kp) -  vf(k) + v0(kp) -  v0(k) ) /  dps0(k) + &
                            vflx(km)  * (  vf(k) - vf(km) +  v0(k) - v0(km) ) / dps0(km) )
         elseif( k .eq. max( ktop_mxen(ixi), ktop_mxen(ixf) ) ) then
            sten_dis(k) = - g / 4._r8 * ( &
                            uflx(km)  * ( uf(k) - uf(km) + u0(k) - u0(km) ) / dps0(km) + &
                            vflx(km)  * ( vf(k) - vf(km) + v0(k) - v0(km) ) / dps0(km) )
         endif
         slten_dis(k) = sten_dis(k)
         qtten_dis(k) = 0._r8 
         uten_dis(k)  = 0._r8
         vten_dis(k)  = 0._r8
         qvten_dis(k) = 0._r8
         qlten_dis(k) = 0._r8
         qiten_dis(k) = 0._r8
         do mt = 1, ncnst
            trten_dis(k,mt) = 0._r8 
         enddo
         sten_out(i,k)  = sten_out(i,k)  +  sten_dis(k)
      enddo

      ! -------------------------- !
      ! 3. Other Layer-Mean Values !
      ! -------------------------- !

      do k = 1, max( ktop_mxen(ixi), ktop_mxen(ixf) )
         kv = mkx + 1 - k
         ! ----------------------------------------------- !
         ! rqc_l, rqc_i, rqc are set to zero in the UNICON !
         ! ----------------------------------------------- !
         rqc_l_out(i,k)          = ( 1._r8 - cuorg ) *     rqc_l_mxen(k,ixi) + cuorg *     rqc_l_mxen(k,ixf) 
         rqc_i_out(i,k)          = ( 1._r8 - cuorg ) *     rqc_i_mxen(k,ixi) + cuorg *     rqc_i_mxen(k,ixf)
         rqc_out(i,k)            = ( 1._r8 - cuorg ) *       rqc_mxen(k,ixi) + cuorg *       rqc_mxen(k,ixf)
         rnc_l_out(i,k)          = ( 1._r8 - cuorg ) *     rnc_l_mxen(k,ixi) + cuorg *     rnc_l_mxen(k,ixf) 
         rnc_i_out(i,k)          = ( 1._r8 - cuorg ) *     rnc_i_mxen(k,ixi) + cuorg *     rnc_i_mxen(k,ixf)
         ! -------------------------------------------------------------------------------------------------------------- !
         ! Detrained mass flux and condensate: consistent between 'flux-convergence' and 'subsidence-detrainment' formula ! 
         ! -------------------------------------------------------------------------------------------------------------- !
         cmf_det_out(i,k)        = ( 1._r8 - cuorg ) *   cmf_det_mxen(k,ixi) + cuorg *   cmf_det_mxen(k,ixf) 
         ql_det_out(i,k)         = ( 1._r8 - cuorg ) *    ql_det_mxen(k,ixi) + cuorg *    ql_det_mxen(k,ixf) 
         qi_det_out(i,k)         = ( 1._r8 - cuorg ) *    qi_det_mxen(k,ixi) + cuorg *    qi_det_mxen(k,ixf) 
         ! ----------------------------------------------------------------------- !
         ! evapc : Evaporation rate of convective precipitation within environment !
         ! ----------------------------------------------------------------------- !
         evapc_out(i,k)          = ( 1._r8 - cuorg ) *     evapc_mxen(k,ixi) + cuorg *     evapc_mxen(k,ixf)
         ! ------------------------------------------------------------------------ !
         ! am_u(d)  : Updraft ( Downdraft ) fractional area at the layer mid-point. !
         ! qlm_u(d) : In-Updraft ( Downdraft ) LWC obtained by area weighting       !
         ! ------------------------------------------------------------------------ !
         am_u_out(i,k)  = ( 1._r8 - cuorg )* am_u_mxen(k,ixi)                    + cuorg* am_u_mxen(k,ixf)
         qlm_u_out(i,k) = ( 1._r8 - cuorg )*qlm_u_mxen(k,ixi) * am_u_mxen(k,ixi) + cuorg*qlm_u_mxen(k,ixf) * am_u_mxen(k,ixf)
         qim_u_out(i,k) = ( 1._r8 - cuorg )*qim_u_mxen(k,ixi) * am_u_mxen(k,ixi) + cuorg*qim_u_mxen(k,ixf) * am_u_mxen(k,ixf)
         am_d_out(i,k)  = ( 1._r8 - cuorg )* am_d_mxen(k,ixi)                    + cuorg* am_d_mxen(k,ixf)
         qlm_d_out(i,k) = ( 1._r8 - cuorg )*qlm_d_mxen(k,ixi) * am_d_mxen(k,ixi) + cuorg*qlm_d_mxen(k,ixf) * am_d_mxen(k,ixf)
         qim_d_out(i,k) = ( 1._r8 - cuorg )*qim_d_mxen(k,ixi) * am_d_mxen(k,ixi) + cuorg*qim_d_mxen(k,ixf) * am_d_mxen(k,ixf)
         if( am_u_out(i,k) .gt. nonzero ) then
            qlm_u_out(i,k)      = qlm_u_out(i,k) / am_u_out(i,k)
            qim_u_out(i,k)      = qim_u_out(i,k) / am_u_out(i,k)
         else
            am_u_out(i,k)       = 0._r8
            qlm_u_out(i,k)      = 0._r8
            qim_u_out(i,k)      = 0._r8
         endif
         if( am_d_out(i,k) .gt. nonzero ) then
            qlm_d_out(i,k)      = qlm_d_out(i,k) / am_d_out(i,k)
            qim_d_out(i,k)      = qim_d_out(i,k) / am_d_out(i,k)
         else
            am_d_out(i,k)       = 0._r8
            qlm_d_out(i,k)      = 0._r8
            qim_d_out(i,k)      = 0._r8
         endif
      enddo

      ! ------------------------------- !
      ! 4. Single-Value for Each Column !
      ! ------------------------------- !

      ! --------------------------------- !
      ! rliq is set to zero in the UNICON !
      ! --------------------------------- !
      rliq_out(i)                = ( 1._r8 - cuorg ) *        rliq_mxen(ixi) + cuorg *        rliq_mxen(ixf)
      precip_out(i)              = ( 1._r8 - cuorg ) *      precip_mxen(ixi) + cuorg *      precip_mxen(ixf)
      snow_out(i)                = ( 1._r8 - cuorg ) *        snow_mxen(ixi) + cuorg *        snow_mxen(ixf) 
      ! --------------------------------------------------------------------------------------------- !
      ! cnt, cnb : Updraft 'top' and 'base' interface index. It is reasonable to take the maximum and !
      !            minimum respectively for final output instead of cuorg weighting average.          !
      ! --------------------------------------------------------------------------------------------- !
      cnt_out(i)                 = max( cnt_mxen(ixi), cnt_mxen(ixf) )
      cnb_out(i)                 = min( cnb_mxen(ixi), cnb_mxen(ixf) )

      ! --------------------------------------------------------------------- !
      !                                                                       !
      !                  Print-Out Internal Output Variables                  !
      !                                                                       !
      ! Vertical index should be reversed for these internal output variables ! 
      !                                                                       !
      ! --------------------------------------------------------------------- !

      ! ------------------------ !
      ! 1. Flux Interface Values !
      ! ------------------------ !

      do ki = 0, max( ktop_mxen(ixi) - 1, ktop_mxen(ixf) - 1 )
         kvi = mkx - ki
         cmf_out(i,kvi)      = (1._r8 - cuorg)*(  cmf_u_mxen(ki,ixi) -  cmf_d_mxen(ki,ixi)) + cuorg*( cmf_u_mxen(ki,ixf) - &
                                                                                                     cmf_d_mxen(ki,ixf) )   
         uflx_out(i,kvi)     = (1._r8 - cuorg)*( uflx_u_mxen(ki,ixi) + uflx_d_mxen(ki,ixi)) + cuorg*(uflx_u_mxen(ki,ixf) + &
                                                                                                     uflx_d_mxen(ki,ixf) )
         vflx_out(i,kvi)     = (1._r8 - cuorg)*( vflx_u_mxen(ki,ixi) + vflx_d_mxen(ki,ixi)) + cuorg*(vflx_u_mxen(ki,ixf) + &
                                                                                                     vflx_d_mxen(ki,ixf) ) 
         slflx_u_out(i,kvi)  = (1._r8 - cuorg)* slflx_u_mxen(ki,ixi)                         + cuorg *  slflx_u_mxen(ki,ixf) 
         qtflx_u_out(i,kvi)  = (1._r8 - cuorg)* qtflx_u_mxen(ki,ixi)                         + cuorg *  qtflx_u_mxen(ki,ixf)
         uflx_u_out(i,kvi)   = (1._r8 - cuorg)*  uflx_u_mxen(ki,ixi)                         + cuorg *   uflx_u_mxen(ki,ixf)
         vflx_u_out(i,kvi)   = (1._r8 - cuorg)*  vflx_u_mxen(ki,ixi)                         + cuorg *   vflx_u_mxen(ki,ixf)
         cmf_d_out(i,kvi)    = (1._r8 - cuorg)*   cmf_d_mxen(ki,ixi)                         + cuorg *    cmf_d_mxen(ki,ixf) 
         slflx_d_out(i,kvi)  = (1._r8 - cuorg)* slflx_d_mxen(ki,ixi)                         + cuorg *  slflx_d_mxen(ki,ixf)
         qtflx_d_out(i,kvi)  = (1._r8 - cuorg)* qtflx_d_mxen(ki,ixi)                         + cuorg *  qtflx_d_mxen(ki,ixf)
         uflx_d_out(i,kvi)   = (1._r8 - cuorg)*  uflx_d_mxen(ki,ixi)                         + cuorg *   uflx_d_mxen(ki,ixf)
         vflx_d_out(i,kvi)   = (1._r8 - cuorg)*  vflx_d_mxen(ki,ixi)                         + cuorg *   vflx_d_mxen(ki,ixf)

         flxrain_out(i,kvi)         = ( 1._r8 - cuorg ) *  flxrain_u_mxen(ki,ixi)  + cuorg *  flxrain_u_mxen(ki,ixf)
         flxsnow_out(i,kvi)         = ( 1._r8 - cuorg ) *  flxsnow_u_mxen(ki,ixi)  + cuorg *  flxsnow_u_mxen(ki,ixf)

      enddo

      ! ----------------------------- !
      ! 2. Layer-Mean Tendency Values !
      ! ----------------------------- !

      do k = 1, max( ktop_mxen(ixi), ktop_mxen(ixf) )
         kv = mkx + 1 - k
         slten_u_out(i,kv)          = ( 1._r8 - cuorg ) *      slten_u_mxen(k,ixi) + cuorg *      slten_u_mxen(k,ixf)
         qtten_u_out(i,kv)          = ( 1._r8 - cuorg ) *      qtten_u_mxen(k,ixi) + cuorg *      qtten_u_mxen(k,ixf)
         uten_u_out(i,kv)           = ( 1._r8 - cuorg ) *       uten_u_mxen(k,ixi) + cuorg *       uten_u_mxen(k,ixf) 
         vten_u_out(i,kv)           = ( 1._r8 - cuorg ) *       vten_u_mxen(k,ixi) + cuorg *       vten_u_mxen(k,ixf)
         sten_u_out(i,kv)           = ( 1._r8 - cuorg ) *       sten_u_mxen(k,ixi) + cuorg *       sten_u_mxen(k,ixf)
         qvten_u_out(i,kv)          = ( 1._r8 - cuorg ) *      qvten_u_mxen(k,ixi) + cuorg *      qvten_u_mxen(k,ixf)
         qlten_u_out(i,kv)          = ( 1._r8 - cuorg ) *      qlten_u_mxen(k,ixi) + cuorg *      qlten_u_mxen(k,ixf)
         qiten_u_out(i,kv)          = ( 1._r8 - cuorg ) *      qiten_u_mxen(k,ixi) + cuorg *      qiten_u_mxen(k,ixf)
         do mt = 1, ncnst
            trten_u_out(i,kv,mt)    = ( 1._r8 - cuorg ) *   trten_u_mxen(k,mt,ixi) + cuorg *   trten_u_mxen(k,mt,ixf)
         enddo
         slten_d_out(i,kv)          = ( 1._r8 - cuorg ) *      slten_d_mxen(k,ixi) + cuorg *      slten_d_mxen(k,ixf)
         qtten_d_out(i,kv)          = ( 1._r8 - cuorg ) *      qtten_d_mxen(k,ixi) + cuorg *      qtten_d_mxen(k,ixf)
         uten_d_out(i,kv)           = ( 1._r8 - cuorg ) *       uten_d_mxen(k,ixi) + cuorg *       uten_d_mxen(k,ixf)
         vten_d_out(i,kv)           = ( 1._r8 - cuorg ) *       vten_d_mxen(k,ixi) + cuorg *       vten_d_mxen(k,ixf)
         sten_d_out(i,kv)           = ( 1._r8 - cuorg ) *       sten_d_mxen(k,ixi) + cuorg *       sten_d_mxen(k,ixf)
         qvten_d_out(i,kv)          = ( 1._r8 - cuorg ) *      qvten_d_mxen(k,ixi) + cuorg *      qvten_d_mxen(k,ixf)
         qlten_d_out(i,kv)          = ( 1._r8 - cuorg ) *      qlten_d_mxen(k,ixi) + cuorg *      qlten_d_mxen(k,ixf)
         qiten_d_out(i,kv)          = ( 1._r8 - cuorg ) *      qiten_d_mxen(k,ixi) + cuorg *      qiten_d_mxen(k,ixf)
         do mt = 1, ncnst
            trten_d_out(i,kv,mt)    = ( 1._r8 - cuorg ) *   trten_d_mxen(k,mt,ixi) + cuorg *   trten_d_mxen(k,mt,ixf) 
         enddo
         slten_evp_out(i,kv)        = ( 1._r8 - cuorg ) *    slten_evp_mxen(k,ixi) + cuorg *    slten_evp_mxen(k,ixf)
         qtten_evp_out(i,kv)        = ( 1._r8 - cuorg ) *    qtten_evp_mxen(k,ixi) + cuorg *    qtten_evp_mxen(k,ixf)
         uten_evp_out(i,kv)         = ( 1._r8 - cuorg ) *     uten_evp_mxen(k,ixi) + cuorg *     uten_evp_mxen(k,ixf)
         vten_evp_out(i,kv)         = ( 1._r8 - cuorg ) *     vten_evp_mxen(k,ixi) + cuorg *     vten_evp_mxen(k,ixf)
         sten_evp_out(i,kv)         = ( 1._r8 - cuorg ) *     sten_evp_mxen(k,ixi) + cuorg *     sten_evp_mxen(k,ixf)
         qvten_evp_out(i,kv)        = ( 1._r8 - cuorg ) *    qvten_evp_mxen(k,ixi) + cuorg *    qvten_evp_mxen(k,ixf)
         qlten_evp_out(i,kv)        = ( 1._r8 - cuorg ) *    qlten_evp_mxen(k,ixi) + cuorg *    qlten_evp_mxen(k,ixf)
         qiten_evp_out(i,kv)        = ( 1._r8 - cuorg ) *    qiten_evp_mxen(k,ixi) + cuorg *    qiten_evp_mxen(k,ixf)
         do mt = 1, ncnst
            trten_evp_out(i,kv,mt)  = ( 1._r8 - cuorg ) *  trten_evp_mxen(k,mt,ixi) + cuorg *  trten_evp_mxen(k,mt,ixf)
         enddo
         slten_dis_out(i,kv)        =     slten_dis(k)
         qtten_dis_out(i,kv)        =     qtten_dis(k)
         uten_dis_out(i,kv)         =      uten_dis(k)
         vten_dis_out(i,kv)         =      vten_dis(k)
         sten_dis_out(i,kv)         =      sten_dis(k)
         qvten_dis_out(i,kv)        =     qvten_dis(k)
         qlten_dis_out(i,kv)        =     qlten_dis(k)
         qiten_dis_out(i,kv)        =     qiten_dis(k)
         do mt = 1, ncnst
            trten_dis_out(i,kv,mt)  =  trten_dis(k,mt) 
         enddo
         qlten_sub_out(i,kv)        = ( 1._r8 - cuorg ) *    qlten_sub_mxen(k,ixi) + cuorg *    qlten_sub_mxen(k,ixf)
         qiten_sub_out(i,kv)        = ( 1._r8 - cuorg ) *    qiten_sub_mxen(k,ixi) + cuorg *    qiten_sub_mxen(k,ixf)
         qlten_det_out(i,kv)        = ( 1._r8 - cuorg ) *    qlten_det_mxen(k,ixi) + cuorg *    qlten_det_mxen(k,ixf)
         qiten_det_out(i,kv)        = ( 1._r8 - cuorg ) *    qiten_det_mxen(k,ixi) + cuorg *    qiten_det_mxen(k,ixf)

       ! Diagnostic output for UNICON paper
         am_evp_out(i,kv)           = ( 1._r8 - cuorg ) *       am_evp_mxen(k,ixi) + cuorg *       am_evp_mxen(k,ixf)
         am_pu_out(i,kv)            = ( 1._r8 - cuorg ) *        am_pu_mxen(k,ixi) + cuorg *        am_pu_mxen(k,ixf) 
         x_um_out(i,kv)             = ( 1._r8 - cuorg ) *         x_um_mxen(k,ixi) + cuorg *         x_um_mxen(k,ixf) 
         y_um_out(i,kv)             = ( 1._r8 - cuorg ) *         y_um_mxen(k,ixi) + cuorg *         y_um_mxen(k,ixf) 
       ! Diagnostic output for UNICON paper

      enddo

      ! --------------------------------------------------------------- !
      ! 3. Convective Updraft and Downdraft Properties at Interfaces    !
      !    Appropriate mass-flux, area, number weightings are required. !
      ! --------------------------------------------------------------- !

      do ki = 0, max( ktop_mxen(ixi) - 1, ktop_mxen(ixf) - 1 )
         kvi = mkx - ki
         ! ----------------------------- !
         ! Convective Updraft Properties !
         ! ----------------------------- !
         thl_u_out(i,kvi)= (1._r8 - cuorg)*thl_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg*thl_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf) 
         qt_u_out(i,kvi) = (1._r8 - cuorg)* qt_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg* qt_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf) 
         u_u_out(i,kvi)  = (1._r8 - cuorg)*  u_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg*  u_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf)
         v_u_out(i,kvi)  = (1._r8 - cuorg)*  v_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg*  v_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf)
         w_u_out(i,kvi)  = (1._r8 - cuorg)*  w_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg*  w_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf)
         ql_u_out(i,kvi) = (1._r8 - cuorg)* ql_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg* ql_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf)
         qi_u_out(i,kvi) = (1._r8 - cuorg)* qi_u_mxen(ki,ixi)*cmf_u_mxen(ki,ixi) + cuorg* qi_u_mxen(ki,ixf) * cmf_u_mxen(ki,ixf)
         do mt = 1, ncnst
            tr_u_out(i,kvi,mt) = (1._r8 - cuorg)*tr_u_mxen(ki,mt,ixi)*cmf_u_mxen(ki,ixi) + &
                                          cuorg*tr_u_mxen(ki,mt,ixf)*cmf_u_mxen(ki,ixf)
         enddo
         ! ------------------------------------------------------------------------------------------------------- ! 
         ! Note that 'cmf_u_out(i,ki)' has already been computed above with OPPOSITE layer array index ki not kvi. !
         ! The way how 'rad_u_out' is computed is the same as the one in the main code with number weighting.      !
         ! ------------------------------------------------------------------------------------------------------- !
         a_u_out(i,kvi)    = (1._r8 - cuorg)*   a_u_mxen(ki,ixi)                    + cuorg * a_u_mxen(ki,ixf) 
         num_u_out(i,kvi)  = (1._r8 - cuorg)* num_u_mxen(ki,ixi)                    + cuorg * num_u_mxen(ki,ixf) 
         wa_u_out(i,kvi)   = (1._r8 - cuorg)*  wa_u_mxen(ki,ixi) * a_u_mxen(ki,ixi) + cuorg * wa_u_mxen(ki,ixf)*a_u_mxen(ki,ixf)
         qla_u_out(i,kvi)  = (1._r8 - cuorg)* qla_u_mxen(ki,ixi) * a_u_mxen(ki,ixi) + cuorg * qla_u_mxen(ki,ixf)*a_u_mxen(ki,ixf) 
         qia_u_out(i,kvi)  = (1._r8 - cuorg)* qia_u_mxen(ki,ixi) * a_u_mxen(ki,ixi) + cuorg * qia_u_mxen(ki,ixf)*a_u_mxen(ki,ixf)
         rad_u_out(i,kvi)  = (1._r8 - cuorg)* rad_u_mxen(ki,ixi)**2._r8 * num_u_mxen(ki,ixi) + cuorg*rad_u_mxen(ki,ixf)**2._r8 * &
                                                                                                    num_u_mxen(ki,ixf) 
         thva_u_out(i,kvi) = (1._r8 - cuorg)*thva_u_mxen(ki,ixi) * a_u_mxen(ki,ixi) + cuorg *thva_u_mxen(ki,ixf)*a_u_mxen(ki,ixf) 

        ! Diagnostic output for UNICON paper
          a_p_out(i,kvi)                     = ( 1._r8 - cuorg ) *     a_p_mxen(ki,ixi)    + cuorg *     a_p_mxen(ki,ixf) 
          x_p_out(i,kvi)                     = ( 1._r8 - cuorg ) *     x_p_mxen(ki,ixi)    + cuorg *     x_p_mxen(ki,ixf) 
          y_p_out(i,kvi)                     = ( 1._r8 - cuorg ) *     y_p_mxen(ki,ixi)    + cuorg *     y_p_mxen(ki,ixf) 
        ! Diagnostic output for UNICON paper

         ! ------------------------------------------------------------------------------- !
         ! Final normalization and averaging                                               !
         ! Be careful of the different index ki ( not kvi ) only in the mass flux.         !
         ! Note also that nonzero net mass flux, area, number concentration are guaranteed !
         ! at all the interfaces we are considering in this do block.                      !
         ! However, if cuorg = 0 or 1, ero value may come out due to cuorg weighting.      !
         ! Thus, we should carefully treat below using if block.                           !
         ! For full consistency, I use all the mass-flux, area, and number concentration   !
         ! consistencies at the same time.                                                 ! 
         ! ------------------------------------------------------------------------------- !
         if( cmf_u_out(i,ki) .gt. nonzero .and. a_u_out(i,kvi) .gt. nonzero .and. num_u_out(i,kvi) .gt. nonzero ) then
            thl_u_out(i,kvi)                 =       thl_u_out(i,kvi) /  cmf_u_out(i,ki)
            qt_u_out(i,kvi)                  =        qt_u_out(i,kvi) /  cmf_u_out(i,ki)
            u_u_out(i,kvi)                   =         u_u_out(i,kvi) /  cmf_u_out(i,ki)
            v_u_out(i,kvi)                   =         v_u_out(i,kvi) /  cmf_u_out(i,ki)
            w_u_out(i,kvi)                   =         w_u_out(i,kvi) /  cmf_u_out(i,ki)
            ql_u_out(i,kvi)                  =        ql_u_out(i,kvi) /  cmf_u_out(i,ki)
            qi_u_out(i,kvi)                  =        qi_u_out(i,kvi) /  cmf_u_out(i,ki)
            do mt = 1, ncnst
               tr_u_out(i,kvi,mt)            =     tr_u_out(i,kvi,mt) /  cmf_u_out(i,ki) 
            enddo
            wa_u_out(i,kvi)                  =        wa_u_out(i,kvi) /   a_u_out(i,kvi)
            qla_u_out(i,kvi)                 =       qla_u_out(i,kvi) /   a_u_out(i,kvi)
            qia_u_out(i,kvi)                 =       qia_u_out(i,kvi) /   a_u_out(i,kvi)
            rad_u_out(i,kvi)                 = sqrt( rad_u_out(i,kvi) / num_u_out(i,kvi) )
            thva_u_out(i,kvi)                =      thva_u_out(i,kvi) /   a_u_out(i,kvi)
         else
            cmf_u_out(i,ki)                  = 0._r8
            a_u_out(i,kvi)                   = 0._r8
            num_u_out(i,kvi)                 = 0._r8
            thl_u_out(i,kvi)                 = 0._r8
            qt_u_out(i,kvi)                  = 0._r8
            u_u_out(i,kvi)                   = 0._r8
            v_u_out(i,kvi)                   = 0._r8
            w_u_out(i,kvi)                   = 0._r8
            ql_u_out(i,kvi)                  = 0._r8
            qi_u_out(i,kvi)                  = 0._r8
            do mt = 1, ncnst
               tr_u_out(i,kvi,mt)            = 0._r8
            enddo
            wa_u_out(i,kvi)                  = 0._r8
            qla_u_out(i,kvi)                 = 0._r8
            qia_u_out(i,kvi)                 = 0._r8
            rad_u_out(i,kvi)                 = 0._r8
            thva_u_out(i,kvi)                = 0._r8
         endif
         gamw_u_out(i,kvi)                    = w_u_out(i,kvi) /  max( wa_u_out(i,kvi), nonzero )
         ! ------------------------------- !
         ! Convective Downdraft Properties !
         ! ------------------------------- !
         thl_d_out(i,kvi) = (1._r8 - cuorg)*thl_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg*thl_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         qt_d_out(i,kvi)  = (1._r8 - cuorg)* qt_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg* qt_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         u_d_out(i,kvi)   = (1._r8 - cuorg)*  u_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg*  u_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         v_d_out(i,kvi)   = (1._r8 - cuorg)*  v_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg*  v_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         w_d_out(i,kvi)   = (1._r8 - cuorg)*  w_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg*  w_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         ql_d_out(i,kvi)  = (1._r8 - cuorg)* ql_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg* ql_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         qi_d_out(i,kvi)  = (1._r8 - cuorg)* qi_d_mxen(ki,ixi)*cmf_d_mxen(ki,ixi) + cuorg* qi_d_mxen(ki,ixf)*cmf_d_mxen(ki,ixf)
         do mt = 1, ncnst
            tr_d_out(i,kvi,mt) = (1._r8 - cuorg)*tr_d_mxen(ki,mt,ixi)*cmf_d_mxen(ki,ixi) + &
                                          cuorg*tr_d_mxen(ki,mt,ixf)*cmf_d_mxen(ki,ixf)
         enddo
         ! ------------------------------------------------------------------------------------------------- ! 
         ! Note that 'cmf_d_out(i,kvi)' has already been computed above with the SAME layer array index kvi, !
         ! in contrast to 'cmf_u_out(i,ki)'                                                                  !
         ! ------------------------------------------------------------------------------------------------- !
         a_d_out(i,kvi)   = (1._r8 - cuorg)*  a_d_mxen(ki,ixi)                  + cuorg*  a_d_mxen(ki,ixf)
         wa_d_out(i,kvi)  = (1._r8 - cuorg)* wa_d_mxen(ki,ixi)*a_d_mxen(ki,ixi) + cuorg* wa_d_mxen(ki,ixf)*a_d_mxen(ki,ixf)
         qla_d_out(i,kvi) = (1._r8 - cuorg)*qla_d_mxen(ki,ixi)*a_d_mxen(ki,ixi) + cuorg*qla_d_mxen(ki,ixf)*a_d_mxen(ki,ixf) 
         qia_d_out(i,kvi) = (1._r8 - cuorg)*qia_d_mxen(ki,ixi)*a_d_mxen(ki,ixi) + cuorg*qia_d_mxen(ki,ixf)*a_d_mxen(ki,ixf)
         ! ------------------------------------------------------------------------------- !
         ! Final normalization and averaging                                               !
         ! Be careful of the same index ki ( not kvi ) for the mass flux too in contrast   !
         ! to convective updraft.                                                          !
         ! ------------------------------------------------------------------------------- !
         if( cmf_d_out(i,kvi) .gt. nonzero .and. a_d_out(i,kvi) .gt. nonzero ) then
            thl_d_out(i,kvi)                 =       thl_d_out(i,kvi) / cmf_d_out(i,kvi)
            qt_d_out(i,kvi)                  =        qt_d_out(i,kvi) / cmf_d_out(i,kvi)
            u_d_out(i,kvi)                   =         u_d_out(i,kvi) / cmf_d_out(i,kvi)
            v_d_out(i,kvi)                   =         v_d_out(i,kvi) / cmf_d_out(i,kvi)
            w_d_out(i,kvi)                   =         w_d_out(i,kvi) / cmf_d_out(i,kvi)
            ql_d_out(i,kvi)                  =        ql_d_out(i,kvi) / cmf_d_out(i,kvi)
            qi_d_out(i,kvi)                  =        qi_d_out(i,kvi) / cmf_d_out(i,kvi)
            do mt = 1, ncnst
               tr_d_out(i,kvi,mt)            =     tr_d_out(i,kvi,mt) / cmf_d_out(i,kvi) 
            enddo
            wa_d_out(i,kvi)                  =        wa_d_out(i,kvi) /   a_d_out(i,kvi)
            qla_d_out(i,kvi)                 =       qla_d_out(i,kvi) /   a_d_out(i,kvi)
            qia_d_out(i,kvi)                 =       qia_d_out(i,kvi) /   a_d_out(i,kvi)
         else
            cmf_d_out(i,kvi)                 =   0._r8 
            a_d_out(i,kvi)                   =   0._r8
            thl_d_out(i,kvi)                 =   0._r8
            qt_d_out(i,kvi)                  =   0._r8
            u_d_out(i,kvi)                   =   0._r8
            v_d_out(i,kvi)                   =   0._r8
            w_d_out(i,kvi)                   =   0._r8
            ql_d_out(i,kvi)                  =   0._r8
            qi_d_out(i,kvi)                  =   0._r8
            do mt = 1, ncnst
               tr_d_out(i,kvi,mt)            =   0._r8
            enddo
            wa_d_out(i,kvi)                  =   0._r8
            qla_d_out(i,kvi)                 =   0._r8
            qia_d_out(i,kvi)                 =   0._r8
         endif
      enddo

      ! ---------------------------------------------------------------------------------- !
      ! 4. Individual Segment Properties of Convective Updraft and Downdraft at Interfaces !
      !    In this case, instead of doing 'cuorg' average, just print out individual       ! 
      !    mixing air properties.                                                          ! 
      ! ---------------------------------------------------------------------------------- !

      do iter = 1, niter
         do msfc = 1, nseg
            do ki = 0, max( ktop_mxen(ixi) - 1, ktop_mxen(ixf) - 1 )

               kvi = mkx - ki
               thl_u_msfc_out(i,kvi,msfc,iter)        =         thl_u_msfc_mxen(ki,msfc,iter)
               qt_u_msfc_out(i,kvi,msfc,iter)         =          qt_u_msfc_mxen(ki,msfc,iter)
               u_u_msfc_out(i,kvi,msfc,iter)          =           u_u_msfc_mxen(ki,msfc,iter)
               v_u_msfc_out(i,kvi,msfc,iter)          =           v_u_msfc_mxen(ki,msfc,iter)
               w_u_msfc_out(i,kvi,msfc,iter)          =           w_u_msfc_mxen(ki,msfc,iter)
               ql_u_msfc_out(i,kvi,msfc,iter)         =          ql_u_msfc_mxen(ki,msfc,iter)
               qi_u_msfc_out(i,kvi,msfc,iter)         =          qi_u_msfc_mxen(ki,msfc,iter)
               do mt = 1, ncnst
                  tr_u_msfc_out(i,kvi,msfc,mt,iter)   =       tr_u_msfc_mxen(ki,msfc,mt,iter)
               enddo
               cmf_u_msfc_out(i,kvi,msfc,iter)        =         cmf_u_msfc_mxen(ki,msfc,iter)
               a_u_msfc_out(i,kvi,msfc,iter)          =           a_u_msfc_mxen(ki,msfc,iter)
               num_u_msfc_out(i,kvi,msfc,iter)        =         num_u_msfc_mxen(ki,msfc,iter)
               rad_u_msfc_out(i,kvi,msfc,iter)        =         rad_u_msfc_mxen(ki,msfc,iter) 

               eps0_u_msfc_out(i,kvi,msfc,iter)       =        eps0_u_msfc_mxen(ki,msfc,iter)
               eps_u_msfc_out(i,kvi,msfc,iter)        =         eps_u_msfc_mxen(ki,msfc,iter)
               del_u_msfc_out(i,kvi,msfc,iter)        =         del_u_msfc_mxen(ki,msfc,iter)
               eeps_u_msfc_out(i,kvi,msfc,iter)       =        eeps_u_msfc_mxen(ki,msfc,iter)
               ddel_u_msfc_out(i,kvi,msfc,iter)       =        ddel_u_msfc_mxen(ki,msfc,iter)
               xc_u_msfc_out(i,kvi,msfc,iter)         =          xc_u_msfc_mxen(ki,msfc,iter)
               xs_u_msfc_out(i,kvi,msfc,iter)         =          xs_u_msfc_mxen(ki,msfc,iter)
               xemin_u_msfc_out(i,kvi,msfc,iter)      =       xemin_u_msfc_mxen(ki,msfc,iter)
               xemax_u_msfc_out(i,kvi,msfc,iter)      =       xemax_u_msfc_mxen(ki,msfc,iter)
               cridis_u_msfc_out(i,kvi,msfc,iter)     =      cridis_u_msfc_mxen(ki,msfc,iter)
               thvcuenv_u_msfc_out(i,kvi,msfc,iter)   =    thvcuenv_u_msfc_mxen(ki,msfc,iter) 
               thvegenv_u_msfc_out(i,kvi,msfc,iter)   =    thvegenv_u_msfc_mxen(ki,msfc,iter)
               thvxsenv_u_msfc_out(i,kvi,msfc,iter)   =    thvxsenv_u_msfc_mxen(ki,msfc,iter) 
               fmix_u_msfc_out(i,kvi,msfc,iter)       =        fmix_u_msfc_mxen(ki,msfc,iter)
               cmfumix_u_msfc_out(i,kvi,msfc,iter)    =     cmfumix_u_msfc_mxen(ki,msfc,iter)
  
               thl_d_msfc_out(i,kvi,msfc,iter)        =         thl_d_msfc_mxen(ki,msfc,iter) 
               qt_d_msfc_out(i,kvi,msfc,iter)         =          qt_d_msfc_mxen(ki,msfc,iter)
               u_d_msfc_out(i,kvi,msfc,iter)          =           u_d_msfc_mxen(ki,msfc,iter)
               v_d_msfc_out(i,kvi,msfc,iter)          =           v_d_msfc_mxen(ki,msfc,iter)
               w_d_msfc_out(i,kvi,msfc,iter)          =           w_d_msfc_mxen(ki,msfc,iter)
               ql_d_msfc_out(i,kvi,msfc,iter)         =          ql_d_msfc_mxen(ki,msfc,iter)
               qi_d_msfc_out(i,kvi,msfc,iter)         =          qi_d_msfc_mxen(ki,msfc,iter)
               do mt = 1, ncnst
                  tr_d_msfc_out(i,kvi,msfc,mt,iter)   =       tr_d_msfc_mxen(ki,msfc,mt,iter)
               enddo
               cmf_d_msfc_out(i,kvi,msfc,iter)        =         cmf_d_msfc_mxen(ki,msfc,iter) 
               a_d_msfc_out(i,kvi,msfc,iter)          =           a_d_msfc_mxen(ki,msfc,iter)
               wa_d_msfc_out(i,kvi,msfc,iter)         =          wa_d_msfc_mxen(ki,msfc,iter)
               qla_d_msfc_out(i,kvi,msfc,iter)        =         qla_d_msfc_mxen(ki,msfc,iter)
               qia_d_msfc_out(i,kvi,msfc,iter)        =         qia_d_msfc_mxen(ki,msfc,iter)

            enddo

            ptop_msfc_out(i,msfc,iter)             =             ptop_msfc_mxen(msfc,iter) 
            ztop_msfc_out(i,msfc,iter)             =             ztop_msfc_mxen(msfc,iter) 

         enddo
      enddo

      ! ------------------------------------------------------------- !
      !                                                               !
      ! Save the 'inout' variables related to convective organization ! 
      !                                                               !  
      ! ------------------------------------------------------------- ! 

      ! ------------------------------------------------------------------------------------------ ! 
      ! 1. Top height and average height of convective updraft at the single previous time step.   !
      !    Similar to previous treatment, take the maximum value for 'cush' while perform          !
      !    surface mass-flux weighted average for cush_avg. Note that both mxen = 1 and 2 have     !
      !    identical surface mass flux.                                                            !
      !    Thus, it is more reasonable to use cushavg instead of cush for computation of           !
      !    critical distance ( rlc ) for convective updraft mixing as is done in the current code. !
      ! ------------------------------------------------------------------------------------------ !
   
      cush_inout(i)    = max( cush_mxen(ixi), cush_mxen(ixf) )
      cushavg_inout(i) = ( 1._r8 - cuorg ) * cushavg_mxen(ixi) + cuorg * cushavg_mxen(ixf) 

     ! TEST. Apr.15.2014.
     ! if( get_nstep() .eq. 786 .or. get_nstep() .eq. 787 ) then  
     !     write(iulog,*)
     !     write(iulog,*) 'SPARKCONV: cushavg, nstep = ', get_nstep()
     !     write(iulog,*) 'Updated cushavg_inout, pblh, pblhz, pblhp =', cushavg_inout(i), pblh, pblhz, pblhp
     !     write(iulog,*) 'kpblh, kpblhm, kpblh_in(i) = ', kpblh, kpblhm, kpblh_in(i)
     !     write(iulog,*) 'zs0_in(i,0), zs0_in(i,1) = ', zs0_in(i,0), zs0_in(i,1)
     !     write(iulog,*) 'zs0_in(i,kpblhm-1), zs0_in(i,kpblhm), zs0_in(i,kpblh) = ', zs0_in(i,kpblhm-1), zs0_in(i,kpblhm), zs0_in(i,kpblh)
     !     write(iulog,*) 'ps0_in(i,0), ps0_in(i,kpblhm) = ', ps0_in(i,0), ps0_in(i,kpblhm)
     !     write(iulog,*) 'eps_u_msfc_out(i,19,1,1), eps_u_msfc_out(i,20,1,1) = ', & 
     !                 eps_u_msfc_out(i,19,1,1), eps_u_msfc_out(i,20,1,1)
     !     write(iulog,*) 'ktop_mxen(ixi), ktop_mxen(ixf) = ', ktop_mxen(ixi), ktop_mxen(ixf)
     !     write(iulog,*) 'zs0_in(i,ktop_mxen-1), zs0_in(i,ktop_mxen), zs0_in(i,ktop_mxen+1) = ', &
     !                 zs0_in(i,ktop_mxen-1), zs0_in(i,ktop_mxen), zs0_in(i,ktop_mxen+1)
     !     write(iulog,*)  
     ! endif 
     ! TEST. Apr.15.2014.

      ! --------------------------------------------------------------------------------------------- !
      ! 4. The mass and scalar properties of detrained airs into each layer from                      !  
      !    individual convective updraft and downdraft at the previous time step.                     !
      !    Due to the tracer array, I cannot save multi-time step informations.                       !
      !    Aug.01.2011. Brian Juwon Park's 9th Birthday. I modified below part to correctly compute   !
      !                 mass-flux weighted average for mixing with multiple mixing environmental air. !
      !                 Below bloak is a new code on this day.                                        !
      !                 Note that for this 'inout' variables, we must use 'do k = 1, mkx' not         !
      !                'k = 1, max( ktop_mxen(1), ktop_mxen(2) )' to correctly put zero values in the !
      !                 layers above 'max( ktop_mxen(1), ktop_mxen(2) )'.                             !
      ! --------------------------------------------------------------------------------------------- !

      do k = 1, mkx
         cu_cmfum_out(i,k)          = ( 1._r8 - cuorg ) * cu_cmfum_mxen(k,ixi) + cuorg * cu_cmfum_mxen(k,ixf)
         ! ---------------------------- !
         ! Total average detrained airs ! 
         ! ---------------------------- ! ------------------------------------------------------------------------------------------------- !
         ! Aug.12.2011. For more geographically organized mixing, it might be better to include only                                        !
         !              cu_thlr_mxen(k,2). But may not.                                                                                     !
         ! Sep.15.2011. In order to reduce the numerical sensitivity due to on-and-off behavior of convection, define below 'cu_thlr_inout' !
         !              variables as the average of previous two time steps, not just at the previous time step.                            !  
         !              In defining individual cu_thlr_mxen, I used 'am_u' instead of 'cu_cmfr_mxen' for nonzero constraint above.          !
         !              However, below block indicates that the use of 'cu_cmfr_mxen' instead of 'am_u' is much more easy to use.           !
         !              Thus, I modified previous computation using the 'cu_cmfr_mxen'. Since 'cu_thlr_inout' is anomaly values,            !
         !              this use of 'cu_cmfr_mxen' seems to perfectly OK even though 'cu_thlr_mxen' were defined using 'am_u'.              !
         !              In fact, exact two time step average is impossible unless we save multi-time step fields, which is very hard due to !
         !              large number of tracers. Thus, I simply do mass weighting average by using the newly generated 'cu_thlr' at         !
         !              the current time step and the passed 'cu_thlr' from previous time steps. This seems to be most reasonable.          !
         !              However, due to the uncertainty on how the treat the amount of mass itself ( clearly, we cannot add them all !!! ), !
         !              we cannot easily apply this method. Thus, I am leaving as a future work.                                            !
         !              However, current one-time step saving seems to be perfectly OK, which is actually very well mimicking what is       ! 
         !              happening in real nature.                                                                                           !
         !              also, when consideing cumulus updraft trailing previous cumulus updraft, current one time step methods seems to be  !
         !              most perfectly conceptually reasonable. So, let's keep the current one-time-step method which is perfect.           !
         ! -------------------------------------------------------------------------------------------------------------------------------- !  
         cu_cmfr_inout(i,k) =(1._r8 - cuorg)*cu_cmfr_mxen(k,ixi)                     + cuorg*                    cu_cmfr_mxen(k,ixf)
         cu_thlr_inout(i,k) =(1._r8 - cuorg)*cu_thlr_mxen(k,ixi)*cu_cmfr_mxen(k,ixi) + cuorg*cu_thlr_mxen(k,ixf)*cu_cmfr_mxen(k,ixf)
         cu_qtr_inout(i,k)  =(1._r8 - cuorg)* cu_qtr_mxen(k,ixi)*cu_cmfr_mxen(k,ixi) + cuorg* cu_qtr_mxen(k,ixf)*cu_cmfr_mxen(k,ixf)
         cu_ur_inout(i,k)   =(1._r8 - cuorg)*  cu_ur_mxen(k,ixi)*cu_cmfr_mxen(k,ixi) + cuorg*  cu_ur_mxen(k,ixf)*cu_cmfr_mxen(k,ixf)
         cu_vr_inout(i,k)   =(1._r8 - cuorg)*  cu_vr_mxen(k,ixi)*cu_cmfr_mxen(k,ixi) + cuorg*  cu_vr_mxen(k,ixf)*cu_cmfr_mxen(k,ixf)
         cu_qlr_inout(i,k)  =(1._r8 - cuorg)* cu_qlr_mxen(k,ixi)*cu_cmfr_mxen(k,ixi) + cuorg* cu_qlr_mxen(k,ixf)*cu_cmfr_mxen(k,ixf)
         cu_qir_inout(i,k)  =(1._r8 - cuorg)* cu_qir_mxen(k,ixi)*cu_cmfr_mxen(k,ixi) + cuorg* cu_qir_mxen(k,ixf)*cu_cmfr_mxen(k,ixf)
         do mt = 1, ncnst
            cu_trr_inout(i,k,mt) = (1._r8 - cuorg)*cu_trr_mxen(k,mt,ixi)*cu_cmfr_mxen(k,ixi) + &
                                            cuorg*cu_trr_mxen(k,mt,ixf)*cu_cmfr_mxen(k,ixf)
         enddo
         if( cu_cmfr_inout(i,k) .gt. nonzero ) then
            cu_thlr_inout(i,k)       =   cu_thlr_inout(i,k) / cu_cmfr_inout(i,k)
            cu_qtr_inout(i,k)        =    cu_qtr_inout(i,k) / cu_cmfr_inout(i,k)
            cu_ur_inout(i,k)         =     cu_ur_inout(i,k) / cu_cmfr_inout(i,k)
            cu_vr_inout(i,k)         =     cu_vr_inout(i,k) / cu_cmfr_inout(i,k)
            cu_qlr_inout(i,k)        =    cu_qlr_inout(i,k) / cu_cmfr_inout(i,k)
            cu_qir_inout(i,k)        =    cu_qir_inout(i,k) / cu_cmfr_inout(i,k)
            do mt = 1, ncnst
               cu_trr_inout(i,k,mt)  = cu_trr_inout(i,k,mt) / cu_cmfr_inout(i,k)
            enddo
         else
            cu_cmfr_inout(i,k)       = 0._r8
            cu_thlr_inout(i,k)       = 0._r8
            cu_qtr_inout(i,k)        = 0._r8
            cu_ur_inout(i,k)         = 0._r8
            cu_vr_inout(i,k)         = 0._r8
            cu_qlr_inout(i,k)        = 0._r8
            cu_qir_inout(i,k)        = 0._r8
            do mt = 1, ncnst
               cu_trr_inout(i,k,mt)  = 0._r8
            enddo
         endif
         ! ------------------------------------------------------------------------------------------------------------------ !
         ! Aug.2.2011. For diagnostic purpose, print-out 'cu_thvr_inout' and 'cu_rhr_inout'. Note that condensate is computed !
         !             from 'thlr' and 'qtr' by direct variable conversion.                                                   !
         ! ------------------------------------------------------------------------------------------------------------------ !
         call conden( p0(k), cu_thlr_inout(i,k) + thl0(k), cu_qtr_inout(i,k) + qt0(k), th, qv, ql, qi, qse, id_check )
         cu_thvr_inout(i,k) = th * ( 1._r8 + zvir * qv - ql - qi ) - thv0(k)
         cu_rhr_inout(i,k)  = max( 0._r8, min( 1._r8, qv / max( nonzero, qse ) ) ) - rh0(k)
         ! --------------------------------------------- !
         ! Detrained airs only from convective downdraft ! 
         ! --------------------------------------------- !
         cu_cmfrd_out(i,k) = (1._r8 - cuorg)*cu_cmfrd_mxen(k,ixi)                      + &
                                        cuorg*cu_cmfrd_mxen(k,ixf)
         cu_thlrd_out(i,k) = (1._r8 - cuorg)*cu_thlrd_mxen(k,ixi)*cu_cmfrd_mxen(k,ixi) + &
                                        cuorg*cu_thlrd_mxen(k,ixf)*cu_cmfrd_mxen(k,ixf)
         cu_qtrd_out(i,k)  = (1._r8 - cuorg)* cu_qtrd_mxen(k,ixi)*cu_cmfrd_mxen(k,ixi) + &
                                        cuorg* cu_qtrd_mxen(k,ixf)*cu_cmfrd_mxen(k,ixf)
         cu_urd_out(i,k)   = (1._r8 - cuorg)*  cu_urd_mxen(k,ixi)*cu_cmfrd_mxen(k,ixi) + &
                                        cuorg*  cu_urd_mxen(k,ixf)*cu_cmfrd_mxen(k,ixf)
         cu_vrd_out(i,k)   = (1._r8 - cuorg)*  cu_vrd_mxen(k,ixi)*cu_cmfrd_mxen(k,ixi) + &
                                        cuorg*  cu_vrd_mxen(k,ixf)*cu_cmfrd_mxen(k,ixf)
         cu_qlrd_out(i,k)  = (1._r8 - cuorg)* cu_qlrd_mxen(k,ixi)*cu_cmfrd_mxen(k,ixi) + &
                                        cuorg* cu_qlrd_mxen(k,ixf)*cu_cmfrd_mxen(k,ixf)
         cu_qird_out(i,k)  = (1._r8 - cuorg)* cu_qird_mxen(k,ixi)*cu_cmfrd_mxen(k,ixi) + &
                                        cuorg* cu_qird_mxen(k,ixf)*cu_cmfrd_mxen(k,ixf)
         do mt = 1, ncnst
            cu_trrd_out(i,k,mt) = ( 1._r8 - cuorg ) * cu_trrd_mxen(k,mt,ixi) * cu_cmfrd_mxen(k,ixi) + &
                                               cuorg * cu_trrd_mxen(k,mt,ixf) * cu_cmfrd_mxen(k,ixf)
         enddo
         if( cu_cmfrd_out(i,k) .gt. nonzero ) then
            cu_thlrd_out(i,k)      =   cu_thlrd_out(i,k) / cu_cmfrd_out(i,k)
            cu_qtrd_out(i,k)       =    cu_qtrd_out(i,k) / cu_cmfrd_out(i,k)
            cu_urd_out(i,k)        =     cu_urd_out(i,k) / cu_cmfrd_out(i,k)
            cu_vrd_out(i,k)        =     cu_vrd_out(i,k) / cu_cmfrd_out(i,k)
            cu_qlrd_out(i,k)       =    cu_qlrd_out(i,k) / cu_cmfrd_out(i,k)
            cu_qird_out(i,k)       =    cu_qird_out(i,k) / cu_cmfrd_out(i,k)
            do mt = 1, ncnst
               cu_trrd_out(i,k,mt) = cu_trrd_out(i,k,mt) / cu_cmfrd_out(i,k)
            enddo
         else
            cu_cmfrd_out(i,k)      = 0._r8
            cu_thlrd_out(i,k)      = 0._r8
            cu_qtrd_out(i,k)       = 0._r8
            cu_urd_out(i,k)        = 0._r8
            cu_vrd_out(i,k)        = 0._r8
            cu_qlrd_out(i,k)       = 0._r8
            cu_qird_out(i,k)       = 0._r8
            do mt = 1, ncnst
               cu_trrd_out(i,k,mt) = 0._r8
            enddo
         endif
      enddo

      ! -------------------------------------------------------------------------------------------------------- ! 
      ! 5. Convective Organization Parameter                                                                     !
      !    Compute updated convective organization for use at the next time step's ( t + dt ) convection scheme. !
      !    Impose lower and upper limits of [ 0, 1 ].                                                            !
      !    Note that cuorg at the next time step ( cuorg_out ) is updated using grid-mean organization forcing.  !
      !    Since 'tau_org' is and should be defined over the grid-mean, it is not a 'mxen-arrayed' variable but  !
      !    a single value for the entire column.                                                                 !
      ! -------------------------------------------------------------------------------------------------------- !

      cuorg_inout(i) = cuorg

      ! -------------------------------------------------------------------------------------------------------- ! 
      ! 6. Convective Organization Parameter - Time evolution of the difference of PBL-mean conservative scalars !
      !    between off-wake region and grid-mean. At the beginning of next time time, convective organization    !
      !    parameter ( cuorg ) will be diagnostically computed using the PBL-mean vertial potential difference   !
      !    between the off-wake area and the grid-mean.                                                          !
      !    The reason to compute cuorg at the next time step diagnostically is to impose consistency in the      !
      !    computed cuorg at the right time step for correct diagnostic output.                                  ! 
      !    The lower and upper limits of [ 0, 1 ] on the cuorg is also imposed.                                  !
      !    Note that cuorg at the next time step ( cuorg_out ) is updated using grid-mean organization forcing.  !
      !    Since 'tau_org' is and should be defined over the grid-mean, it is not a 'mxen-arrayed' variable but  !
      !    a single value for the entire column.                                                                 !
      !    Note also that we are using the same damping time scale for all scalars since physically, damping of  !
      !    the difference between off-wake and grid-mean value occurs due to lateral mixing along the wall of    !
      !    of the wake from surface to the PBL top interface.                                                    !
      !    CAUTION : We may want to impose a certain on the prognosed 'delta_' variables in order to prevent     !
      !              the disruption by unreasonable values.                                                      !
      !                                                                                                          !
      !    Sep.03.2011. As a forcing of 'delta_thl_PBL' variables in the below, I should also include diabatic   !
      !                 forcing within wake area ( aw * Qw ) as well as adiabatic flux at the PBL top interface  !
      !                 from the several selected downdraft components. Since diabatic forcing ( evaporation of  !
      !                 convective precipitation ) is larger over the land than over the ocean, this will        !
      !                 help to increase convective organization over the land, if I define 'cuorg' as the       !
      !                 normalized 'delta_thv_PBL' as in the current code.                                       !     
      !    Sep.07.2011. Redefine total forcings including both adiabatic forcing and two diabatic forcings both  !
      !                 within convective downdraft and environment.                                             ! 
      !                 Also rename the variables.                                                               !
      !    Sep.07.2011. Also includes prognostic equation for the wake area.                                     !
      !                 I should be very careful when I define 'tau' for 'niter' ensemble.                       !
      !                 In addition, I should define 'tau' differently for individual conservative scalars and   !
      !                 wake area. This is because adjustment of wake by surface flux can generate diffferent    !
      !                 relaxation time scale for each conservative scalar.                                      !
      !                 Note that 'tau' for conservative scalar for iter = 1 is the same as the value for        !
      !                 iter = 2 even though we include surface flux adjustment. However, 'tau' for wake area    !
      !                 are different for iter = 1 and iter = 2. Regardless of whether it is conservative scalar !
      !                 of wak area, we should perform 'cuorg' weighted average using '1/tau(iter)' and then     !
      !                 the resulting '1/tau(avg)' should be inversed to obtain the final 'tau(avg)'. This is    !
      !                 because the prognostic organization equations below are 1st order LINEAR equation.       !
      !                 By doing this, we can compuetely remove the ambiguity on the treatment of 'tau' in the   !
      !                 ensemble mean both for conservative scalar and wake area.                                !           
      !                 In the below block, 'taui' denotes '1/tau', i.e., inverse tau.                           ! 
      !                 Note that 'tau' for 'awk' can be negative, which is completely OK.                       !
      !    Sep.09.2011. I don't prognose 'thv' any more - it will be computed diagnostically at the beginning of !
      !                 the next time step for full model consistency.                                           !     
      ! -------------------------------------------------------------------------------------------------------- !

      ! -------------------------------- !
      ! Inverse of Adjustment Time Scale !
      ! -------------------------------- ! 

      taui_thl_out(i)              = ( 1._r8 - cuorg ) *            taui_thl_mxen(ixi) + cuorg *            taui_thl_mxen(ixf)
      taui_qt_out(i)               = ( 1._r8 - cuorg ) *             taui_qt_mxen(ixi) + cuorg *             taui_qt_mxen(ixf)
      taui_u_out(i)                = ( 1._r8 - cuorg ) *              taui_u_mxen(ixi) + cuorg *              taui_u_mxen(ixf)
      taui_v_out(i)                = ( 1._r8 - cuorg ) *              taui_v_mxen(ixi) + cuorg *              taui_v_mxen(ixf)
      do mt = 1, ncnst
         taui_tr_out(i,mt)         = ( 1._r8 - cuorg ) *          taui_tr_mxen(mt,ixi) + cuorg *          taui_tr_mxen(mt,ixf)
      enddo
      taui_awk_out(i)              = ( 1._r8 - cuorg ) *            taui_awk_mxen(ixi) + cuorg *            taui_awk_mxen(ixf)

      del_org_out(i)               = ( 1._r8 - cuorg ) *              del_org_mxen(ixi) + cuorg *             del_org_mxen(ixf)
      del0_org_out(i)              = ( 1._r8 - cuorg ) *             del0_org_mxen(ixi) + cuorg *            del0_org_mxen(ixf)
 
      ! -------------------------- !
      ! Total Organization Forcing ! 
      ! -------------------------- !

      thl_orgforce_out(i)          = ( 1._r8 - cuorg ) *        thl_orgforce_mxen(ixi) + cuorg *        thl_orgforce_mxen(ixf)
      qt_orgforce_out(i)           = ( 1._r8 - cuorg ) *         qt_orgforce_mxen(ixi) + cuorg *         qt_orgforce_mxen(ixf)
      u_orgforce_out(i)            = ( 1._r8 - cuorg ) *          u_orgforce_mxen(ixi) + cuorg *          u_orgforce_mxen(ixf)
      v_orgforce_out(i)            = ( 1._r8 - cuorg ) *          v_orgforce_mxen(ixi) + cuorg *          v_orgforce_mxen(ixf)
      do mt = 1, ncnst
         tr_orgforce_out(i,mt)     = ( 1._r8 - cuorg ) *      tr_orgforce_mxen(mt,ixi) + cuorg *      tr_orgforce_mxen(mt,ixf)
      enddo
      awk_orgforce_out(i)          = ( 1._r8 - cuorg ) *        awk_orgforce_mxen(ixi) + cuorg *        awk_orgforce_mxen(ixf)

      ! ----------------------------------------------------- !
      ! Individual Organization Forcing for Diagnostic Output ! 
      ! ----------------------------------------------------- !

      thl_orgforce_flx_out(i)      = ( 1._r8 - cuorg ) *    thl_orgforce_flx_mxen(ixi) + cuorg *    thl_orgforce_flx_mxen(ixf)
      thl_orgforce_und_out(i)      = ( 1._r8 - cuorg ) *    thl_orgforce_und_mxen(ixi) + cuorg *    thl_orgforce_und_mxen(ixf)
      thl_orgforce_env_out(i)      = ( 1._r8 - cuorg ) *    thl_orgforce_env_mxen(ixi) + cuorg *    thl_orgforce_env_mxen(ixf)

      qt_orgforce_flx_out(i)       = ( 1._r8 - cuorg ) *     qt_orgforce_flx_mxen(ixi) + cuorg *     qt_orgforce_flx_mxen(ixf)
      qt_orgforce_und_out(i)       = ( 1._r8 - cuorg ) *     qt_orgforce_und_mxen(ixi) + cuorg *     qt_orgforce_und_mxen(ixf)
      qt_orgforce_env_out(i)       = ( 1._r8 - cuorg ) *     qt_orgforce_env_mxen(ixi) + cuorg *     qt_orgforce_env_mxen(ixf)

      u_orgforce_flx_out(i)        = ( 1._r8 - cuorg ) *      u_orgforce_flx_mxen(ixi) + cuorg *      u_orgforce_flx_mxen(ixf)
      u_orgforce_und_out(i)        = ( 1._r8 - cuorg ) *      u_orgforce_und_mxen(ixi) + cuorg *      u_orgforce_und_mxen(ixf)
      u_orgforce_env_out(i)        = ( 1._r8 - cuorg ) *      u_orgforce_env_mxen(ixi) + cuorg *      u_orgforce_env_mxen(ixf)

      v_orgforce_flx_out(i)        = ( 1._r8 - cuorg ) *      v_orgforce_flx_mxen(ixi) + cuorg *      v_orgforce_flx_mxen(ixf)
      v_orgforce_und_out(i)        = ( 1._r8 - cuorg ) *      v_orgforce_und_mxen(ixi) + cuorg *      v_orgforce_und_mxen(ixf)
      v_orgforce_env_out(i)        = ( 1._r8 - cuorg ) *      v_orgforce_env_mxen(ixi) + cuorg *      v_orgforce_env_mxen(ixf)

      awk_orgforce_flx_out(i)      = ( 1._r8 - cuorg ) *    awk_orgforce_flx_mxen(ixi) + cuorg *    awk_orgforce_flx_mxen(ixf)
      awk_orgforce_mix_out(i)      = ( 1._r8 - cuorg ) *    awk_orgforce_mix_mxen(ixi) + cuorg *    awk_orgforce_mix_mxen(ixf)

      cmf_d_org_pblh_out(i)        = ( 1._r8 - cuorg ) *      cmf_d_org_pblh_mxen(ixi) + cuorg *      cmf_d_org_pblh_mxen(ixf)

      ! -------------------------------------------------------------------------------------------------------------- !
      ! Sep.07.2011. If the prognosed wake area becomes larger than the available maximum value ( awk_PBL_max ),       !
      !              increae the detrainment rate of the wake 'del_wk', such that the prognozed wake area becomes      !
      !              identical to awk_PBL_max. Also correspondingly modify 'awk_orgforce_out(i)' and 'taui_thl_out(i)' ! 
      !              etc. which are explicit function of 'del_wk'.                                                     !
      !              Note that this 'awk_PBL_max' constraint is applied to the prognozed 'awk_PBL' but BEFORE applying !
      !              non-homogeneous distribution assumption of wake properties. When we apply non-homogeneous distr.  !
      !              assumption of wake properties at the next time step, we an further decrease the wake area. So,    !
      !              we have double safety which is good. Note also that due to the coding structure, it is very hard  !
      !              to apply the awk_PBL_max constraint after applying non-homogeneous distribution assumption of     !
      !              wake properties. Current UNICON code is reasonable good and perfect.                              !
      ! Sep.09.2011. I don't prognose 'thv' anymore.                                                                   !
      !              I should re-check whether below formula is exactly correct or not - should re-derive for check.   !
      ! Sep.11.2011. In the below block, 'tmp1' is an adjusted new 'del_wk'.                                           !
      !              Note that the initial 'del_wk' specified from the parameter sentence can be zero. So, in order to !
      !              prevent division by zero, I used max( del_wk, nonzero ) in the below block.                       !
      ! REQUIREMENT: I must re-check whether below formula is not correct or not.                                      !
      !              I have not done re-checking yet as of Sep.11.2011.                                                !
      ! Sep.11.2011. In case of 'awk_PBL', I should adjust 'awk_orgforce' which is a function of del_wk, while         !
      !              in case of 'thl,qt' , I should adjust 'taui'         which is a function of del_wk.               !
      !              In the below block, 'tmp1' is a newly adjusted 'del_wk'.                                          !
      !              I performed correct re-checking and corrected the previously wrong code.                          ! 
      ! -------------------------------------------------------------------------------------------------------------- !

      if( abs(taui_awk_out(i)) .gt. nonzero ) then
         awk_PBL_inout(i)         =  awk_PBL *    exp( - dt * taui_awk_out(i) ) + awk_orgforce_out(i) / taui_awk_out(i) * &
                                    ( 1._r8 - exp( - dt * taui_awk_out(i) ) )
      else
         awk_PBL_inout(i)         =  awk_PBL * ( 1._r8 - dt * taui_awk_out(i) ) + awk_orgforce_out(i) * dt
      endif

      if( awk_PBL_inout(i) .gt. awk_PBL_max ) then
         awk_PBL_inout(i) = awk_PBL_max
         tmp2 = awk_orgforce_out(i)
         if( abs(taui_awk_out(i)) .gt. nonzero ) then
            awk_orgforce_out(i) = ( awk_PBL_max - awk_PBL * exp( - dt * taui_awk_out(i) ) ) * taui_awk_out(i) / ( 1._r8 - &
                                 exp( - dt * taui_awk_out(i) ) ) 
         else
            awk_orgforce_out(i) = ( awk_PBL_max - awk_PBL * ( 1._r8 - dt * taui_awk_out(i) ) ) / dt
         endif
         tmp1 = del_wk_eff + tmp2 - awk_orgforce_out(i)
         tmp3 = 1._r8 / max( nonzero, awk_PBL * ( 1._r8 - awk_PBL ) )
         taui_thl_out(i)      =  taui_thl_out(i)   + tmp3 * ( tmp1 - del_wk_eff )
         taui_qt_out(i)       =  taui_qt_out(i)    + tmp3 * ( tmp1 - del_wk_eff )
         taui_u_out(i)        =  taui_u_out(i)     + tmp3 * ( tmp1 - del_wk_eff )
         taui_v_out(i)        =  taui_v_out(i)     + tmp3 * ( tmp1 - del_wk_eff )
         do mt = 1, ncnst
            taui_tr_out(i,mt) =  taui_tr_out(i,mt) + tmp3 * ( tmp1 - del_wk_eff )
         enddo
      endif

      ! ------------------------------------------------------ !
      ! Compute final prognosed state                          !
      ! I may need to impose a bound on the prognosed results. !
      ! Sep.09.2011. I don't prognose 'thv' anymore.           !  
      ! ------------------------------------------------------ !

      if( abs(taui_thl_out(i)) .gt. nonzero ) then
         delta_thl_PBL_inout(i)      =    delta_thl_PBL *    exp( - dt * taui_thl_out(i) ) + thl_orgforce_out(i) / &
                                          taui_thl_out(i) * ( 1._r8 - exp( - dt * taui_thl_out(i) ) )
      else
         delta_thl_PBL_inout(i)      =    delta_thl_PBL * ( 1._r8 - dt * taui_thl_out(i) ) + thl_orgforce_out(i) * dt
      endif
      if( abs(taui_qt_out(i)) .gt. nonzero ) then
         delta_qt_PBL_inout(i)       =     delta_qt_PBL *    exp( - dt * taui_qt_out(i) ) + qt_orgforce_out(i) / &
                                           taui_qt_out(i) * ( 1._r8 - exp( - dt * taui_qt_out(i) ) )
      else
         delta_qt_PBL_inout(i)       =     delta_qt_PBL * ( 1._r8 - dt * taui_qt_out(i) ) + qt_orgforce_out(i) * dt
      endif
      if( abs(taui_u_out(i)) .gt. nonzero ) then
         delta_u_PBL_inout(i)        =      delta_u_PBL *    exp( - dt * taui_u_out(i) ) + u_orgforce_out(i) / &
                                            taui_u_out(i) * ( 1._r8 - exp( - dt * taui_u_out(i) ) )
      else
         delta_u_PBL_inout(i)        =      delta_u_PBL * ( 1._r8 - dt * taui_u_out(i) ) + u_orgforce_out(i) * dt
      endif
      if( abs(taui_v_out(i)) .gt. nonzero ) then
         delta_v_PBL_inout(i)        =      delta_v_PBL *    exp( - dt * taui_v_out(i) ) + v_orgforce_out(i) / &
                                            taui_v_out(i) * ( 1._r8 - exp( - dt * taui_v_out(i) ) )
      else
         delta_v_PBL_inout(i)        =      delta_v_PBL * ( 1._r8 - dt * taui_v_out(i) ) + v_orgforce_out(i) * dt
      endif
      do mt = 1, ncnst
         if ( abs(taui_tr_out(i,mt)) .gt. nonzero) then
            delta_tr_PBL_inout(i,mt) = delta_tr_PBL(mt) *    exp( - dt * taui_tr_out(i,mt) ) + tr_orgforce_out(i,mt) / &
              taui_tr_out(i,mt) * ( 1._r8 - exp( - dt * taui_tr_out(i,mt) ) )
         else
            delta_tr_PBL_inout(i,mt) = delta_tr_PBL(mt) * ( 1._r8 - dt * taui_tr_out(i,mt) ) + tr_orgforce_out(i,mt) * dt
         endif
      enddo

      ! ------------------------------------------------------------------------------------------------- !
      ! Substitute organization tendency to the tracer array in all layers when organization is advected. !
      ! In future, I can include the heterogeneity of individual 25 tracer components for completeness    !
      ! which however will increase computation time.                                                     ! 
      ! Add an offset to make the output tracer to be positive, so that advection scheme can handle.      !
      ! Note that the same amount of offset should be subtracted when 'delta_xxx' is extracted from the   !
      ! input tr0_in at the beginning of this program.                                                    !
      ! Note that to ensure that 'delta_xxx' are not modified by wet or dry deposition and other physical !
      ! processes other than UNICON and horizontal advection, a guaranteed update of tracer arrays        !
      ! just before advection scheme is done within tphysac.F90, which in fact makes below block          !
      ! unnecessary. However, for clarify, let's keep below block. It is not a harm at all.               !
      ! ------------------------------------------------------------------------------------------------- !

      if( iorg_adv ) then
         do k = 1, mkx
            trten_out(i,k,i_awk) = (                   awk_PBL_inout(i)             - tr0_in(i,k,i_awk) ) / dt    
! JHYoon
!           trten_out(i,k,i_thl) = ( max( 0._r8, delta_thl_PBL_inout(i) + 100._r8 ) - tr0_in(i,k,i_thl) ) / dt    
!           trten_out(i,k,i_qt)  = ( max( 0._r8,  delta_qt_PBL_inout(i) + 100._r8 ) - tr0_in(i,k,i_qt)  ) / dt    
!           trten_out(i,k,i_u)   = ( max( 0._r8,   delta_u_PBL_inout(i) + 100._r8 ) - tr0_in(i,k,i_u)   ) / dt    
!           trten_out(i,k,i_v)   = ( max( 0._r8,   delta_v_PBL_inout(i) + 100._r8 ) - tr0_in(i,k,i_v)   ) / dt    
            trten_out(i,k,i_thl) = ( max( 0._r8, delta_thl_PBL_inout(i) + 10._r8 ) - tr0_in(i,k,i_thl) ) / dt
            trten_out(i,k,i_qt)  = ( max( 0._r8,  delta_qt_PBL_inout(i) + 0.01_r8 ) - tr0_in(i,k,i_qt)  ) / dt
            trten_out(i,k,i_u)   = ( max( 0._r8,   delta_u_PBL_inout(i) + 10._r8 ) - tr0_in(i,k,i_u)   ) / dt
            trten_out(i,k,i_v)   = ( max( 0._r8,   delta_v_PBL_inout(i) + 10._r8 ) - tr0_in(i,k,i_v)   ) / dt
! JHYoon
         enddo
      endif

   enddo ! End of do i = 1, iend.  This is a column loop.
         
   ! ----------------------- !
   ! End of Main Computation !
   ! ----------------------- !

   ! ---------------------------------------- !
   ! Writing main diagnostic output variables !
   ! ---------------------------------------- !

   call outfld('cmf_SP'  ,                             cmf_out, mix, lchnk) 
   tmpi_array(:,0:mkx) = slflx_out(:,mkx:0:-1) 
   call outfld('slflx_SP',                          tmpi_array, mix, lchnk) 
   tmpi_array(:,0:mkx) = qtflx_out(:,mkx:0:-1)
   call outfld('qtflx_SP',                          tmpi_array, mix, lchnk) 
   call outfld('uflx_SP' ,                            uflx_out, mix, lchnk) 
   call outfld('vflx_SP' ,                            vflx_out, mix, lchnk) 

   call outfld('flxrain_SP' ,                      flxrain_out, mix, lchnk) 
   call outfld('flxsnow_SP' ,                      flxsnow_out, mix, lchnk) 

   tmpi_array(:,0:mkx) = cmf_u_out(:,mkx:0:-1)
   call outfld('cmf_u_SP'  ,                        tmpi_array, mix, lchnk) 
   call outfld('slflx_u_SP',                       slflx_u_out, mix, lchnk) 
   call outfld('qtflx_u_SP',                       qtflx_u_out, mix, lchnk) 
   call outfld('uflx_u_SP' ,                        uflx_u_out, mix, lchnk) 
   call outfld('vflx_u_SP' ,                        vflx_u_out, mix, lchnk) 

   call outfld('cmf_d_SP'  ,                         cmf_d_out, mix, lchnk) 
   call outfld('slflx_d_SP',                       slflx_d_out, mix, lchnk) 
   call outfld('qtflx_d_SP',                       qtflx_d_out, mix, lchnk) 
   call outfld('uflx_d_SP' ,                        uflx_d_out, mix, lchnk) 
   call outfld('vflx_d_SP' ,                        vflx_d_out, mix, lchnk) 

   call outfld('cush_SP',                           cush_inout, mix, lchnk) 
   call outfld('cushavg_SP',                     cushavg_inout, mix, lchnk) 
   call outfld('cuorg_SP',                         cuorg_inout, mix, lchnk) 
   call outfld('Radius_SP',                   rad_u_out(:,mkx), mix, lchnk) 
   call outfld('sgh30_SP',                            sgh30_in, mix, lchnk) 

   call outfld('kw_SP',                                 kw_out, mix, lchnk)

   call outfld('sigma_w_SP',                       sigma_w_out, mix, lchnk)
   call outfld('sigma_thl_SP',                   sigma_thl_out, mix, lchnk)
   call outfld('sigma_qt_SP',                     sigma_qt_out, mix, lchnk) 
   call outfld('sigma_u_SP',                       sigma_u_out, mix, lchnk) 
   call outfld('sigma_v_SP',                       sigma_v_out, mix, lchnk) 

   call outfld('w_org_SP',                           w_org_out, mix, lchnk)
   call outfld('thl_org_SP',                       thl_org_out, mix, lchnk)
   call outfld('qt_org_SP',                         qt_org_out, mix, lchnk)
   call outfld('u_org_SP',                           u_org_out, mix, lchnk)
   call outfld('v_org_SP',                           v_org_out, mix, lchnk)

   call outfld('tkes_SP',                             tkes_out, mix, lchnk)
   call outfld('went_SP',                             went_out, mix, lchnk)
   call outfld('went_eff_SP',                     went_eff_out, mix, lchnk)

   tmpm_array(:,1:mkx) = am_u_out(:,mkx:1:-1)
   call outfld('am_u_SP',                           tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = qlm_u_out(:,mkx:1:-1)
   call outfld('qlm_u_SP',                          tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = qim_u_out(:,mkx:1:-1)
   call outfld('qim_u_SP',                          tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = am_d_out(:,mkx:1:-1)
   call outfld('am_d_SP',                           tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = qlm_d_out(:,mkx:1:-1)
   call outfld('qlm_d_SP',                          tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = qim_d_out(:,mkx:1:-1)
   call outfld('qim_d_SP',                          tmpm_array, mix, lchnk) 

   call outfld('slten_u_SP' ,                      slten_u_out, mix, lchnk) 
   call outfld('qtten_u_SP' ,                      qtten_u_out, mix, lchnk) 
   call outfld('uten_u_SP'  ,                       uten_u_out, mix, lchnk) 
   call outfld('vten_u_SP'  ,                       vten_u_out, mix, lchnk)
   call outfld('sten_u_SP'  ,                       sten_u_out, mix, lchnk) 
   call outfld('qvten_u_SP' ,                      qvten_u_out, mix, lchnk) 
   call outfld('qlten_u_SP' ,                      qlten_u_out, mix, lchnk) 
   call outfld('qiten_u_SP' ,                      qiten_u_out, mix, lchnk) 
   call outfld('nlten_u_SP' ,        trten_u_out(:,:,ixnumliq), mix, lchnk) 
   call outfld('niten_u_SP' ,        trten_u_out(:,:,ixnumice), mix, lchnk) 
#ifdef MODAL_AERO
   do m = 1, ntot_amode
      l = numptr_amode(m)
      varname = trim(cnst_name(l))//'_u_SP'
      call outfld( trim(varname),        trten_u_out(:,:,l),    mix, lchnk )
      do lspec = 1, nspec_amode(m)
         l = lmassptr_amode(lspec,m)
         varname = trim(cnst_name(l))//'_u_SP'
         call outfld( trim(varname),        trten_u_out(:,:,l), mix, lchnk )
      enddo
   enddo
#endif

   call outfld('slten_d_SP' , slten_d_out, mix, lchnk) 
   call outfld('qtten_d_SP' , qtten_d_out, mix, lchnk) 
   call outfld('uten_d_SP'  ,  uten_d_out, mix, lchnk) 
   call outfld('vten_d_SP'  ,  vten_d_out, mix, lchnk)
   call outfld('sten_d_SP'  ,  sten_d_out, mix, lchnk) 
   call outfld('qvten_d_SP' , qvten_d_out, mix, lchnk) 
   call outfld('qlten_d_SP' , qlten_d_out, mix, lchnk) 
   call outfld('qiten_d_SP' , qiten_d_out, mix, lchnk) 
   call outfld('nlten_d_SP' , trten_d_out(:,:,ixnumliq), mix, lchnk) 
   call outfld('niten_d_SP' , trten_d_out(:,:,ixnumice), mix, lchnk) 
#ifdef MODAL_AERO
   do m = 1, ntot_amode
      l = numptr_amode(m)
      varname = trim(cnst_name(l))//'_d_SP'
      call outfld( trim(varname),  trten_d_out(:,:,l),    mix, lchnk )
      do lspec = 1, nspec_amode(m)
         l = lmassptr_amode(lspec,m)
         varname = trim(cnst_name(l))//'_d_SP'
         call outfld( trim(varname),  trten_d_out(:,:,l), mix, lchnk )
      enddo
   enddo
#endif
   
   call outfld('slten_evp_SP' , slten_evp_out, mix, lchnk) 
   call outfld('qtten_evp_SP' , qtten_evp_out, mix, lchnk) 
   call outfld('uten_evp_SP'  ,  uten_evp_out, mix, lchnk) 
   call outfld('vten_evp_SP'  ,  vten_evp_out, mix, lchnk)
   call outfld('sten_evp_SP'  ,  sten_evp_out, mix, lchnk) 
   call outfld('qvten_evp_SP' , qvten_evp_out, mix, lchnk) 
   call outfld('qlten_evp_SP' , qlten_evp_out, mix, lchnk) 
   call outfld('qiten_evp_SP' , qiten_evp_out, mix, lchnk) 
   call outfld('nlten_evp_SP' , trten_evp_out(:,:,ixnumliq), mix, lchnk) 
   call outfld('niten_evp_SP' , trten_evp_out(:,:,ixnumice), mix, lchnk) 
#ifdef MODAL_AERO
   do m = 1, ntot_amode
      l = numptr_amode(m)
      varname = trim(cnst_name(l))//'_evp_SP'
      call outfld( trim(varname),  trten_evp_out(:,:,l),    mix, lchnk )
      do lspec = 1, nspec_amode(m)
         l = lmassptr_amode(lspec,m)
         varname = trim(cnst_name(l))//'_evp_SP'
         call outfld( trim(varname),  trten_evp_out(:,:,l), mix, lchnk )
      enddo
   enddo
#endif

   call outfld('slten_dis_SP' , slten_dis_out, mix, lchnk) 
   call outfld('qtten_dis_SP' , qtten_dis_out, mix, lchnk) 
   call outfld('uten_dis_SP'  ,  uten_dis_out, mix, lchnk) 
   call outfld('vten_dis_SP'  ,  vten_dis_out, mix, lchnk)
   call outfld('sten_dis_SP'  ,  sten_dis_out, mix, lchnk) 
   call outfld('qvten_dis_SP' , qvten_dis_out, mix, lchnk) 
   call outfld('qlten_dis_SP' , qlten_dis_out, mix, lchnk) 
   call outfld('qiten_dis_SP' , qiten_dis_out, mix, lchnk) 
   call outfld('nlten_dis_SP' , trten_dis_out(:,:,ixnumliq), mix, lchnk) 
   call outfld('niten_dis_SP' , trten_dis_out(:,:,ixnumice), mix, lchnk) 
#ifdef MODAL_AERO
   do m = 1, ntot_amode
      l = numptr_amode(m)
      varname = trim(cnst_name(l))//'_dis_SP'
      call outfld( trim(varname),  trten_dis_out(:,:,l),    mix, lchnk )
      do lspec = 1, nspec_amode(m)
         l = lmassptr_amode(lspec,m)
         varname = trim(cnst_name(l))//'_dis_SP'
         call outfld( trim(varname),  trten_dis_out(:,:,l), mix, lchnk )
      enddo
   enddo
#endif
   
   call outfld('qlten_sub_SP'   , qlten_sub_out,   mix, lchnk) 
   call outfld('qiten_sub_SP'   , qiten_sub_out,   mix, lchnk) 
  
   call outfld('qlten_det_SP'   , qlten_det_out,   mix, lchnk) 
   call outfld('qiten_det_SP'   , qiten_det_out,   mix, lchnk) 

   call outfld('thl_u_SP'       ,     thl_u_out,   mix, lchnk)
   call outfld('qt_u_SP'        ,      qt_u_out,   mix, lchnk)
   call outfld('u_u_SP'         ,       u_u_out,   mix, lchnk)
   call outfld('v_u_SP'         ,       v_u_out,   mix, lchnk)
   call outfld('w_u_SP'         ,       w_u_out,   mix, lchnk)
   call outfld('ql_u_SP'        ,      ql_u_out,   mix, lchnk)
   call outfld('qi_u_SP'        ,      qi_u_out,   mix, lchnk)
   call outfld('wa_u_SP'        ,      wa_u_out,   mix, lchnk)
   call outfld('qla_u_SP'       ,     qla_u_out,   mix, lchnk)
   call outfld('qia_u_SP'       ,     qia_u_out,   mix, lchnk)
   call outfld('a_u_SP'         ,       a_u_out,   mix, lchnk)
   call outfld('rad_u_SP'       ,     rad_u_out,   mix, lchnk)
   call outfld('num_u_SP'       ,     num_u_out,   mix, lchnk)
   call outfld('gamw_u_SP'      ,    gamw_u_out,   mix, lchnk)
   call outfld('nl_u_SP'        ,      tr_u_out(:,:,ixnumliq),   mix, lchnk)
   call outfld('ni_u_SP'        ,      tr_u_out(:,:,ixnumice),   mix, lchnk)
   call outfld('thva_u_SP'      ,    thva_u_out,   mix, lchnk)

   call outfld('a_p_SP'         ,       a_p_out,   mix, lchnk)
   call outfld('am_evp_SP'      ,    am_evp_out,   mix, lchnk)
   call outfld('am_pu_SP'       ,     am_pu_out,   mix, lchnk)
   call outfld('x_p_SP'         ,       x_p_out,   mix, lchnk)
   call outfld('y_p_SP'         ,       y_p_out,   mix, lchnk)
   call outfld('x_um_SP'        ,      x_um_out,   mix, lchnk)
   call outfld('y_um_SP'        ,      y_um_out,   mix, lchnk)

   call outfld('thl_d_SP'       ,     thl_d_out,   mix, lchnk)
   call outfld('qt_d_SP'        ,      qt_d_out,   mix, lchnk)
   call outfld('u_d_SP'         ,       u_d_out,   mix, lchnk)
   call outfld('v_d_SP'         ,       v_d_out,   mix, lchnk)
   call outfld('w_d_SP'         ,       w_d_out,   mix, lchnk)
   call outfld('ql_d_SP'        ,      ql_d_out,   mix, lchnk)
   call outfld('qi_d_SP'        ,      qi_d_out,   mix, lchnk)
   call outfld('wa_d_SP'        ,      wa_d_out,   mix, lchnk)
   call outfld('qla_d_SP'       ,     qla_d_out,   mix, lchnk)
   call outfld('qia_d_SP'       ,     qia_d_out,   mix, lchnk)
   call outfld('a_d_SP'         ,       a_d_out,   mix, lchnk)
   call outfld('nl_d_SP'        ,      tr_d_out(:,:,ixnumliq),   mix, lchnk)
   call outfld('ni_d_SP'        ,      tr_d_out(:,:,ixnumice),   mix, lchnk)

   tmpi_array(:,0:mkx) = thv_b_out(:,mkx:0:-1)
   call outfld('thv_b_SP',         tmpi_array,     mix, lchnk)
   tmpi_array(:,0:mkx) = thv_t_out(:,mkx:0:-1)
   call outfld('thv_t_SP',         tmpi_array,     mix, lchnk)
   tmpi_array(:,0:mkx) = thv_mt_out(:,mkx:0:-1)
   call outfld('thv_mt_SP',        tmpi_array,     mix, lchnk)
   tmpi_array(:,0:mkx) = thv_min_out(:,mkx:0:-1)
   call outfld('thv_min_SP',       tmpi_array,     mix, lchnk)

   tmpm_array(:,1:mkx) = cu_cmfr_inout(:,mkx:1:-1)
   call outfld('cu_cmfr_SP',                       tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_thlr_inout(:,mkx:1:-1)
   call outfld('cu_thlr_SP',                       tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_qtr_inout(:,mkx:1:-1)
   call outfld('cu_qtr_SP',                        tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_qlr_inout(:,mkx:1:-1)
   call outfld('cu_qlr_SP',                        tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_qir_inout(:,mkx:1:-1)
   call outfld('cu_qir_SP',                        tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_ur_inout(:,mkx:1:-1)
   call outfld('cu_ur_SP',                         tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_vr_inout(:,mkx:1:-1)
   call outfld('cu_vr_SP',                         tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_thvr_inout(:,mkx:1:-1)
   call outfld('cu_thvr_SP',                       tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_rhr_inout(:,mkx:1:-1)
   call outfld('cu_rhr_SP',                        tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_trr_inout(:,mkx:1:-1,ixnumliq)
   call outfld('cu_nlr_SP',                        tmpm_array, mix, lchnk) 
   tmpm_array(:,1:mkx) = cu_trr_inout(:,mkx:1:-1,ixnumice)
   call outfld('cu_nir_SP',                        tmpm_array, mix, lchnk) 

   do msfc = 1, nseg        
      write(numcha,'(i2.2)') msfc        

      call outfld( 'thl_u'//numcha//'_SP',              thl_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'qt_u'//numcha//'_SP',                qt_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'u_u'//numcha//'_SP',                  u_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'v_u'//numcha//'_SP',                  v_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'w_u'//numcha//'_SP',                  w_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'ql_u'//numcha//'_SP',                ql_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'qi_u'//numcha//'_SP',                qi_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'cmf_u'//numcha//'_SP',              cmf_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'a_u'//numcha//'_SP',                  a_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'num_u'//numcha//'_SP',              num_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'rad_u'//numcha//'_SP',              rad_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'nl_u'//numcha//'_SP',       tr_u_msfc_out(:,:,msfc,ixnumliq,1), mix, lchnk) 
      call outfld( 'ni_u'//numcha//'_SP',       tr_u_msfc_out(:,:,msfc,ixnumice,1), mix, lchnk) 

      call outfld( 'eps0_u'//numcha//'_SP',            eps0_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'eps_u'//numcha//'_SP',              eps_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'del_u'//numcha//'_SP',              del_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'eeps_u'//numcha//'_SP',            eeps_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'ddel_u'//numcha//'_SP',            ddel_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'xc_u'//numcha//'_SP',                xc_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'xs_u'//numcha//'_SP',                xs_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'xemin_u'//numcha//'_SP',          xemin_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'xemax_u'//numcha//'_SP',          xemax_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'cridis_u'//numcha//'_SP',        cridis_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'thvcuenv_u'//numcha//'_SP',    thvcuenv_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'thvegenv_u'//numcha//'_SP',    thvegenv_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'thvxsenv_u'//numcha//'_SP',    thvxsenv_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'fmix_u'//numcha//'_SP',            fmix_u_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'cmfumix_u'//numcha//'_SP',      cmfumix_u_msfc_out(:,:,msfc,1), mix, lchnk) 

      call outfld( 'ptop'//numcha//'_SP',                ptop_msfc_out(:,msfc,1),       mix, lchnk) 
      call outfld( 'ztop'//numcha//'_SP',                ztop_msfc_out(:,msfc,1),       mix, lchnk) 

      call outfld( 'thl_d'//numcha//'_SP',              thl_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'qt_d'//numcha//'_SP',                qt_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'u_d'//numcha//'_SP',                  u_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'v_d'//numcha//'_SP',                  v_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'w_d'//numcha//'_SP',                  w_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'ql_d'//numcha//'_SP',                ql_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'qi_d'//numcha//'_SP',                qi_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'wa_d'//numcha//'_SP',                wa_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'qla_d'//numcha//'_SP',              qla_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'qia_d'//numcha//'_SP',              qia_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'cmf_d'//numcha//'_SP',              cmf_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'a_d'//numcha//'_SP',                  a_d_msfc_out(:,:,msfc,1), mix, lchnk) 
      call outfld( 'nl_d'//numcha//'_SP',       tr_d_msfc_out(:,:,msfc,ixnumliq,1), mix, lchnk) 
      call outfld( 'ni_d'//numcha//'_SP',       tr_d_msfc_out(:,:,msfc,ixnumice,1), mix, lchnk) 

   enddo

   call outfld('thl_orgfce_SP',                                             thl_orgforce_out, mix, lchnk)
   call outfld('qt_orgfce_SP',                                               qt_orgforce_out, mix, lchnk)
   call outfld('u_orgfce_SP',                                                 u_orgforce_out, mix, lchnk)
   call outfld('v_orgfce_SP',                                                 v_orgforce_out, mix, lchnk)
   call outfld('awk_orgfce_SP',                                             awk_orgforce_out, mix, lchnk)

   call outfld('thl_orgfce_f_SP',                                       thl_orgforce_flx_out, mix, lchnk)
   call outfld('qt_orgfce_f_SP',                                         qt_orgforce_flx_out, mix, lchnk)
   call outfld('u_orgfce_f_SP',                                           u_orgforce_flx_out, mix, lchnk)
   call outfld('v_orgfce_f_SP',                                           v_orgforce_flx_out, mix, lchnk)
   call outfld('awk_orgfce_f_SP',                                       awk_orgforce_flx_out, mix, lchnk)

   call outfld('thl_orgfce_u_SP',                                       thl_orgforce_und_out, mix, lchnk)
   call outfld('qt_orgfce_u_SP',                                         qt_orgforce_und_out, mix, lchnk)
   call outfld('u_orgfce_u_SP',                                           u_orgforce_und_out, mix, lchnk)
   call outfld('v_orgfce_u_SP',                                           v_orgforce_und_out, mix, lchnk)
   call outfld('awk_orgfce_m_SP',                                       awk_orgforce_mix_out, mix, lchnk)

   call outfld('thl_orgfce_e_SP',                                       thl_orgforce_env_out, mix, lchnk)
   call outfld('qt_orgfce_e_SP',                                         qt_orgforce_env_out, mix, lchnk)
   call outfld('u_orgfce_e_SP',                                           u_orgforce_env_out, mix, lchnk)
   call outfld('v_orgfce_e_SP',                                           v_orgforce_env_out, mix, lchnk)
   call outfld('cmf_d_orgh_SP',                                           cmf_d_org_pblh_out, mix, lchnk)

   call outfld('taui_thl_SP',                                                   taui_thl_out, mix, lchnk)
   call outfld('taui_qt_SP',                                                     taui_qt_out, mix, lchnk)
   call outfld('taui_u_SP',                                                       taui_u_out, mix, lchnk)
   call outfld('taui_v_SP',                                                       taui_v_out, mix, lchnk)
   call outfld('taui_awk_SP',                                                   taui_awk_out, mix, lchnk)

   call outfld('del_org_SP',                                                     del_org_out, mix, lchnk)
   call outfld('del0_org_SP',                                                   del0_org_out, mix, lchnk)

   return

end subroutine compute_unicon

end module unicon
