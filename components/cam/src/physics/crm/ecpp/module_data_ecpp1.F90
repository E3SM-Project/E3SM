!   file module_data_ecpp1.F
!-----------------------------------------------------------------------

module module_data_ecpp1

   use shr_kind_mod, only: r8=>shr_kind_r8

   ! integer, parameter :: r4=4
   ! integer, parameter :: r8=8


   ! The following are used to dimension several arrays 
   ! declared in module_ecpp_ppdriver.F with "save"
   ! in mmf framework, these arrays will be subr parameters
   ! in wrf-chem framework, doing this is just too much trouble
   ! because of registry limitations
   integer, parameter :: its_ecpptmp=1
   integer, parameter :: ite_ecpptmp=1
   integer, parameter :: jts_ecpptmp=1
   integer, parameter :: jte_ecpptmp=1
   integer, parameter :: kts_ecpptmp=1
   integer, parameter :: kte_ecpptmp=51
   integer, parameter :: ktebnd_ecpptmp=kte_ecpptmp
   integer, parameter :: ktecen_ecpptmp=kte_ecpptmp-1
   integer, parameter :: num_chem_ecpptmp=101


   ! maximum number of ecpp transport classes, used for dimensioning various arrays
   integer, parameter :: maxcls_ecpp=3
   integer, parameter :: maxsub_ecpp=maxcls_ecpp

   ! maximum number of "precipitation types" for wetscav diagnostics
   ! currently this is 1
   ! the wetscav diagnostics are done for each subarea type, so
   ! have info on where (up, down, quiescent) the scavenging happens.
   ! however, they do no account for the fact that precip formed
   ! in updraft can fall (or shift) into quiescent, etc.
   ! eventually it might be 2, so would have diagnostics involving
   ! where precip is formed -- quiescent versus (convective) up/downdrafts)
   integer, parameter :: max_wetdiagtype = 1

   ! set this to .false. for cam3-mmf
   logical, parameter :: hostcode_is_wrfchem = .false.


   ! these are possible values for mtype_updnenv_ecpp_3d & ..._clm3d & ..._clm3d & ..._clm
   integer, parameter :: mtype_updraft_ecpp=1
   integer, parameter :: mtype_dndraft_ecpp=2
   integer, parameter :: mtype_quiescn_ecpp=3
   integer, parameter :: mtype_upempty_ecpp=-1
   integer, parameter :: mtype_dnempty_ecpp=-2
   integer, parameter :: mtype_quempty_ecpp=-3

   ! these are possible values for mtype_clrcldy_ecpp_3d & ..._clm3d & ..._clm3d & ..._clm
   integer, parameter :: mtype_iscloud_ecpp=11
   integer, parameter :: mtype_nocloud_ecpp=0

   ! these are possible values for mtype_precip_ecpp_3d & ..._clm3d & ..._clm3d & ..._clm
   integer, parameter :: mtype_isprecip_ecpp=21
   integer, parameter :: mtype_noprecip_ecpp=0


   ! this flag determines whether updraft & dndraft profiles are calculated
   ! using the "primed" mass fluxes or "full" mass fluxes
   integer, save :: ppopt_updn_prof_aa
   ! these are possible values for the flag
   integer, parameter :: ppopt_updn_prof_aa_wfull=2001
   integer, parameter :: ppopt_updn_prof_aa_wprime=2002


   ! this flag determines whether quiescent subarea mass fluxes are
   ! provided by the host or calculated in the ppm
   integer, save :: ppopt_quiescn_mf
   integer, parameter :: ppopt_quiescn_mf_byhost=2101
   ! these are possible values for the flag
   integer, parameter :: ppopt_quiescn_mf_byppmx1=2101


   ! this flag determines how the quiescent subarea mixing ratios
   ! are obtained for source-sink calculations
   integer, save :: ppopt_quiescn_sosi
   ! these are possible values for the flag
   ! 2201 -- qe = qbar
   integer, parameter :: ppopt_quiescn_sosi_x1=2201
   ! 2202 -- ae*qe = max( 0.0, (qbar-au*qu-ad*qd) )
   integer, parameter :: ppopt_quiescn_sosi_x2=2202


   ! this flag determines how the subgrid vertical fluxes (and the
   ! finite differencing for flux divergence) is calculated
   integer, save :: ppopt_chemtend_wq
   ! these are possible values for the flag
   ! 2301 -- vertflux = mu*qu + md*qd - (mu+md)*qbar; 
   !         upstream approach for qbar at layer boundaries
   integer, parameter :: ppopt_chemtend_wq_wfullx1=2301
   ! 2302 -- vertflux = mu'*qu + md'*qd - (mu'+md')*qbar; 
   !         upstream approach for qbar at layer boundaries
   integer, parameter :: ppopt_chemtend_wq_wprimex1=2302


   ! this flag determines how the sub-time-step for integrating the
   ! d(qbar)/dt equation is determined
   ! (use sub-timesteps to keep courant number < 1 and 
   !  avoid negative mixing ratios)
   integer, save :: ppopt_chemtend_dtsub
   ! these are possible values for the flag
   integer, parameter :: ppopt_chemtend_dtsub_x1=2401
   ! 2401 -- dumcournomax = max( dumcourentmax, dumcouroutbmax )
   integer, parameter :: ppopt_chemtend_dtsub_x2=2402
   ! 2402 -- dumcournomax = max( dumcourentmax, dumcouroutamax, 
   !                                            dumcouroutbmax )
   integer, parameter :: ppopt_chemtend_dtsub_x3=2403
   ! 2403 -- dtstep_sub = largest value that does not produce
   !         negative mixing ratios


   ! this flag determines how frequently xxx
   ! is called to calculate up & dndraft profiles and source/sinks
   integer, save :: ppopt_chemtend_updnfreq
   ! these are possible values for the flag
   integer, parameter :: ppopt_chemtend_updnfreq_x1=2501
   ! 2501 -- called just once, when istep_sub=1
   integer, parameter :: ppopt_chemtend_updnfreq_x2=2502
   ! 2502 -- called for each istep_sub


   integer, parameter :: lunout = 0


   ! index of quiescent transport class
   integer, parameter :: jcls_quiescn = 1
   integer, parameter :: jcls_qu = jcls_quiescn


   ! subarea-average vertical mass fluxes (kg/m2/s) smaller than this 
   !    are treated as zero
   ! largest expected flux is ~1 (rho=1, w=10, afrac=0.1)
   !    so could expect truncation errors between 1e-7 and 1e-6
   real(r8), parameter :: mf_smallaa = 1.0e-6


   !! subarea-average vertical mass fluxes (kg/m2/s) smaller than
   !! aw_draft_cut*rho are treated as zero
   !! note that with a*w = 1e-4 m/s, dz over 1 day = 8.6 m which is small
   ! real(r8), parameter :: aw_draft_cut = 1.0e-4_r8   ! m/s

   !! maximum expected updraft
   ! real(r8), parameter :: w_draft_max = 50.0_r8   ! m/s

   !! fractional areas below afrac_cut are ignored
   ! real(r8), parameter :: afrac_cut = aw_draft_cut/w_draft_max
   ! real(r8), parameter :: afrac_cut_bb  = afrac_cut*0.5_r8
   ! real(r8), parameter :: afrac_cut_0p5 = afrac_cut*0.5_r8
   ! real(r8), parameter :: afrac_cut_0p2 = afrac_cut*0.2_r8
   ! real(r8), parameter :: afrac_cut_0p1 = afrac_cut*0.1_r8

   real(r8), save :: aw_draft_cut = 1.0e-4_r8   ! m/s
   ! maximum expected updraft
   real(r8), save :: w_draft_max = 50.0_r8   ! m/s
   ! fractional areas below afrac_cut are ignored
   real(r8), save :: afrac_cut
   real(r8), save :: afrac_cut_bb, afrac_cut_0p5, afrac_cut_0p2, afrac_cut_0p1


   ! draft lifetime (s)
   real, save :: draft_lifetime

   ! activat_onoff_ecpp - if positive, do aerosol activation in ecpp
   ! (set to +1 for normal runs)
   integer, save :: activat_onoff_ecpp

   ! cldchem_onoff_ecpp - if positive, do aerosol activation in ecpp
   ! (set to +1 for normal runs)
   integer, save :: cldchem_onoff_ecpp

   ! rename_onoff_ecpp - if positive, do aerosol activation in ecpp
   ! (set to +1 for normal runs)
   integer, save :: rename_onoff_ecpp

   ! wetscav_onoff_ecpp - if positive, do aerosol activation in ecpp
   ! (set to +1 for normal runs)
   integer, save :: wetscav_onoff_ecpp

   ! iflag_ecpp_startup_acw_partition - when positive, do 
   ! "special partitioning" of cloudborne and interstitial aerosol to
   ! clear and cloudy subareas (cloudy gets less interstitial than clear)
   ! in subr parampollu_tdx_startup
   ! for normal runs, set this to +1
   integer, save :: iflag_ecpp_startup_acw_partition

   ! iflag_ecpp_startup_host_chemtend - when positive, apply
   ! host changes to chem mixing ratios (e.g., emissions, gas chem)
   ! in subr parampollu_tdx_startup
   ! for normal runs, set this to +1
   integer, save :: iflag_ecpp_startup_host_chemtend

   ! iflag_ecpp_test_bypass_1 used for early testing -
   ! when positive, bypass the parampollu_td--- routine
   ! for normal runs, set this to 0
   integer, save :: iflag_ecpp_test_bypass_1

   ! iflag_ecpp_test_fixed_fcloud used for (early) testing with various fixed cloud fracs
   ! for normal runs, set this to zero
   integer, save :: iflag_ecpp_test_fixed_fcloud

   ! "method" flag for parameterized-pollutants module
   ! (set to +2223 for normal runs and in mmf)
   integer, save :: parampollu_opt

   ! minimum fractional area for total quiescent class   
   real(r8), save :: a_quiescn_minaa = 0.60_r8  ! min area for initial total quiescent
   real(r8), save :: a_quiescn_minbb = 0.30_r8  ! min area for final   total quiescent


   integer, save :: num_moist_ecpp
   integer, save :: num_chem_ecpp, param_first_ecpp

   integer, save :: num_moist
   integer, save :: num_chem
   integer, save :: p_qc
   integer, save :: p_qv

   integer, save :: p_num_a01, p_num_cw01
   integer, save :: p_oin_a01, p_oin_cw01
   integer, save :: p_num_a03, p_num_cw03
   integer, save :: p_oin_a03, p_oin_cw03

   ! time step for the ECPP
   ! It is fixed to be 1800 s. The GCM time step can be less than 1800s. 
   ! For example, if GCM time step is 600s, ECPP will be called at every third GCM time step
   real(r8), parameter :: dtstep_pp_input =  1800.0_r8            

end module module_data_ecpp1

