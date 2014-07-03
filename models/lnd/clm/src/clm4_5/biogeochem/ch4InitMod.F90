module ch4InitMod

  !-----------------------------------------------------------------------
  ! !MODULE: initch4Mod
  !
  ! !DESCRIPTION:
  ! Contains cold start initial values and time constant (and flux / diagnostic vars) 
  ! for CH4 scheme.
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initColdCH4
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  subroutine initColdCH4( bounds )
    !
    ! !DESCRIPTION:
    ! - Sets cold start values for time varying values.
    ! Initializes the following time varying variables:
    ! conc_ch4_sat, conc_ch4_unsat, conc_o2_sat, conc_o2_unsat, 
    ! lake_soilc, o2stress, finunduated
    ! - Sets variables for ch4 code that will not be input 
    ! from restart/inic file. 
    ! - Sets values for inactive CH4 columns to spval so that they will 
    ! not be averaged in history file.
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype      , only : pps_a, cps, cws, cch4, col, lun
    use clm_varpar   , only : nlevsoi, nlevgrnd, nlevdecomp
    use clm_varcon   , only : istsoil, istdlak, spval, istcrop
    use clm_varctl   , only : iulog
    use ch4varcon    , only : allowlakeprod
    use decompMod    , only : bounds_type
    use spmdMod      , only : masterproc
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: j,l,c,p      ! indices
    !-----------------------------------------------------------------------

    !   cps%watsat                  Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
    !   cps%cellorg                 Input:  [real(r8) (:,:) ]  column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)
    !   pps_a%grnd_ch4_cond         Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]       
    !   pps_a%rootfr                Output: [real(r8) (:,:) ]  column-averaged root fraction                   
    !   cws%finundated              Output: [real(r8) (:)   ]  inundated gridcell fractional area (excluding dedicated wetland columns)
    !   cch4%finundated_lag         Output: [real(r8) (:)   ]  time-lagged fractional inundated area             
    !   cch4%ch4_surf_diff_sat      Output: [real(r8) (:)   ]  CH4 surface flux (mol/m2/s)                       
    !   cch4%ch4_surf_diff_unsat    Output: [real(r8) (:)   ]  CH4 surface flux (mol/m2/s)                       
    !   cch4%ch4_surf_diff_lake     Output: [real(r8) (:)   ]  CH4 surface flux (mol/m2/s)                       
    !   cch4%ch4_surf_aere_sat      Output: [real(r8) (:)   ]  Total column CH4 aerenchyma (mol/m2/s)            
    !   cch4%ch4_surf_aere_unsat    Output: [real(r8) (:)   ]  Total column CH4 aerenchyma (mol/m2/s)            
    !   cch4%ch4_surf_ebul_sat      Output: [real(r8) (:)   ]  CH4 ebullition to atmosphere (mol/m2/s)           
    !   cch4%ch4_surf_ebul_unsat    Output: [real(r8) (:)   ]  CH4 ebullition to atmosphere (mol/m2/s)           
    !   cch4%ch4_surf_ebul_lake     Output: [real(r8) (:)   ]  CH4 ebullition to atmosphere (mol/m2/s)           
    !   cch4%ch4_oxid_depth_sat     Output: [real(r8) (:,:) ]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_oxid_depth_unsat   Output: [real(r8) (:,:) ]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_oxid_depth_lake    Output: [real(r8) (:,:) ]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_prod_depth_sat     Output: [real(r8) (:,:) ]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%ch4_prod_depth_unsat   Output: [real(r8) (:,:) ]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%ch4_prod_depth_lake    Output: [real(r8) (:,:) ]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%ch4_ebul_total_sat     Output: [real(r8) (:)   ]  Total column CH4 ebullition (mol/m2/s)            
    !   cch4%ch4_ebul_total_unsat   Output: [real(r8) (:)   ]  Total column CH4 ebullition (mol/m2/s)            
    !   cch4%ch4_ebul_depth_sat     Output: [real(r8) (:,:) ]  CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_ebul_depth_unsat   Output: [real(r8) (:,:) ]  CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_aere_depth_sat     Output: [real(r8) (:,:) ]  CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_aere_depth_unsat   Output: [real(r8) (:,:) ]  CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_tran_depth_sat     Output: [real(r8) (:,:) ]  CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_tran_depth_unsat   Output: [real(r8) (:,:) ]  CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%co2_aere_depth_sat     Output: [real(r8) (:,:) ]  CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%co2_aere_depth_unsat   Output: [real(r8) (:,:) ]  CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%o2_oxid_depth_sat      Output: [real(r8) (:,:) ]  O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%o2_oxid_depth_unsat    Output: [real(r8) (:,:) ]  O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%o2_decomp_depth_sat    Output: [real(r8) (:,:) ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%o2_decomp_depth_unsat  Output: [real(r8) (:,:) ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%o2_aere_depth_sat      Output: [real(r8) (:,:) ]  O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%o2_aere_depth_unsat    Output: [real(r8) (:,:) ]  O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%co2_decomp_depth_sat   Output: [real(r8) (:,:) ]  CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%co2_decomp_depth_unsat Output: [real(r8) (:,:) ]  CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
    !   cch4%co2_oxid_depth_sat     Output: [real(r8) (:,:) ]  CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%co2_oxid_depth_unsat   Output: [real(r8) (:,:) ]  CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !   cch4%ch4_dfsat_flux         Output: [real(r8) (:)   ]  CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]
    !   cch4%zwt_ch4_unsat          Output: [real(r8) (:)   ]  depth of water table for unsaturated fraction (m) 
    !   cch4%conc_ch4_lake          Output: [real(r8) (:,:) ]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
    !   cch4%conc_o2_lake           Output: [real(r8) (:,:) ]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
    !   cch4%fphr                   Output: [real(r8) (:,:) ]  fraction of potential HR                        
    !   cch4%sif                    Output: [real(r8) (:)   ]  (unitless) ratio applied to sat. prod. to account for seasonal inundation
    !   cch4%o2stress_unsat         Output: [real(r8) (:,:) ]  Ratio of oxygen available to that demanded by roots, aerobesmethanotrophs (nlevsoi)
    !   cch4%o2stress_sat           Output: [real(r8) (:,:) ]  Ratio of oxygen available to that demanded by roots, aerobesmethanotrophs (nlevsoi)
    !   cch4%ch4stress_unsat        Output: [real(r8) (:,:) ]  Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
    !   cch4%ch4stress_sat          Output: [real(r8) (:,:) ]  Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
    !   cch4%totcolch4              Output: [real(r8) (:)   ]  total methane in soil column (g C / m^2)          
    !   cch4%conc_ch4_sat           Output: [real(r8) (:,:) ]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
    !   cch4%conc_ch4_unsat         Output: [real(r8) (:,:) ]  CH4 conc in each soil layer (mol/m3) (nlevsoi)  
    !   cch4%conc_o2_sat            Output: [real(r8) (:,:) ]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
    !   cch4%conc_o2_unsat          Output: [real(r8) (:,:) ]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
    !   cch4%lake_soilc             Output: [real(r8) (:,:) ]  total soil organic matter found in level (g C / m^3) (nlevsoi)
    !   cch4%qflx_surf_lag          Output: [real(r8) (:)   ]  time-lagged surface runoff (mm H2O /s)            
    !   cch4%layer_sat_lag          Output: [real(r8) (:,:) ]  Lagged saturation status of soil layer in the unsaturated zone (1 = sat)

    if ( masterproc ) write (iulog,*) 'Setting initial data to non-spun up values for CH4 Mod'

    do c = bounds%begc,bounds%endc

       l = col%landunit(c)

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          cch4%conc_ch4_sat(c,1:nlevsoi)   = 0._r8
          cch4%conc_ch4_unsat(c,1:nlevsoi) = 0._r8
          cch4%conc_o2_sat(c,1:nlevsoi)    = 0._r8
          cch4%conc_o2_unsat(c,1:nlevsoi)  = 0._r8
          cch4%o2stress_sat(c,1:nlevsoi)   = 1._r8
          cch4%o2stress_unsat(c,1:nlevsoi) = 1._r8
          cch4%layer_sat_lag(c,1:nlevsoi)  = 1._r8
          cch4%qflx_surf_lag(c)            = 0._r8
          cch4%finundated_lag(c)           = 0._r8
          cws%finundated(c)                = 0._r8
          ! finundated will be used to calculate soil decomposition if anoxia is used 
          ! Note that finundated will be overwritten with cch4%fsat_bef upon reading
          ! a restart file - either in a continuation, branch or startup spun-up case

       else if (lun%itype(l) == istdlak) then
          cch4%conc_ch4_sat(c,1:nlevsoi) = 0._r8
          cch4%conc_o2_sat(c,1:nlevsoi)  = 0._r8
          cch4%lake_soilc(c,1:nlevsoi)   = 580._r8 * cps%cellorg(c,1:nlevsoi)
          ! Need to convert from kg/m^3 organic matter to g C / m^3 (org matter is defined to be 58% C)
          ! Note that cps%cellorg is called in initTimeConstMod - so that needs to be called
          ! before this routine

       end if

       ! Set values for all columns equal to zero below nlevsoi
       cch4%conc_ch4_sat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%conc_ch4_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
       cch4%conc_o2_sat(c,nlevsoi+1:nlevgrnd)    = 0._r8
       cch4%conc_o2_unsat(c,nlevsoi+1:nlevgrnd)  = 0._r8
       cch4%lake_soilc(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%o2stress_sat(c,nlevsoi+1:nlevgrnd)   = 1._r8
       cch4%o2stress_unsat(c,nlevsoi+1:nlevgrnd) = 1._r8
       cch4%layer_sat_lag(c,nlevsoi+1:nlevgrnd)  = 1._r8

       ! Set levels from nlevsoi+1 to nlevgrnd = 0
       cch4%ch4_prod_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%ch4_prod_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%ch4_prod_depth_lake(c,nlevsoi+1:nlevgrnd)    = 0._r8
       cch4%ch4_oxid_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%ch4_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%ch4_oxid_depth_lake(c,nlevsoi+1:nlevgrnd)    = 0._r8
       cch4%o2_oxid_depth_sat(c,nlevsoi+1:nlevgrnd)      = 0._r8
       cch4%o2_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd)    = 0._r8
       cch4%o2_decomp_depth_sat(c,nlevsoi+1:nlevgrnd)    = 0._r8
       cch4%o2_decomp_depth_unsat(c,nlevsoi+1:nlevgrnd)  = 0._r8
       cch4%o2_aere_depth_sat(c,nlevsoi+1:nlevgrnd)      = 0._r8
       cch4%o2_aere_depth_unsat(c,nlevsoi+1:nlevgrnd)    = 0._r8
       cch4%co2_decomp_depth_sat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%co2_decomp_depth_unsat(c,nlevsoi+1:nlevgrnd) = 0._r8
       cch4%co2_oxid_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%co2_oxid_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%ch4_aere_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%ch4_aere_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%ch4_tran_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%ch4_tran_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%co2_aere_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%co2_aere_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%ch4_ebul_depth_sat(c,nlevsoi+1:nlevgrnd)     = 0._r8
       cch4%ch4_ebul_depth_unsat(c,nlevsoi+1:nlevgrnd)   = 0._r8
       cch4%conc_ch4_lake(c,nlevsoi+1:nlevgrnd)          = 0._r8
       cch4%conc_o2_lake(c,nlevsoi+1:nlevgrnd)           = 0._r8
       pps_a%rootfr(c,nlevsoi+1:nlevgrnd)                = 0._r8
       cch4%ch4stress_unsat(c,nlevsoi+1:nlevgrnd)        = 0._r8
       cch4%ch4stress_sat(c,nlevsoi+1:nlevgrnd)          = 0._r8
       cch4%fphr(c,nlevdecomp+1:nlevgrnd)                = 0._r8

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          cch4%conc_ch4_lake(c,:)       = spval
          cch4%conc_o2_lake(c,:)        = spval
          cch4%ch4_surf_diff_lake(c)    = spval
          cch4%ch4_surf_ebul_lake(c)    = spval
          cch4%ch4_prod_depth_lake(c,:) = spval
          cch4%ch4_oxid_depth_lake(c,:) = spval

       else if (lun%itype(l) == istdlak .and. allowlakeprod) then

          cch4%ch4_prod_depth_unsat(c,:)   = spval
          cch4%ch4_oxid_depth_unsat(c,:)   = spval
          cch4%o2_oxid_depth_unsat(c,:)    = spval
          cch4%o2_decomp_depth_unsat(c,:)  = spval
          cch4%o2_aere_depth_unsat(c,:)    = spval
          cch4%co2_decomp_depth_unsat(c,:) = spval
          cch4%co2_oxid_depth_unsat(c,:)   = spval
          cch4%ch4_aere_depth_unsat(c,:)   = spval
          cch4%ch4_tran_depth_unsat(c,:)   = spval
          cch4%co2_aere_depth_unsat(c,:)   = spval
          cch4%ch4_surf_aere_unsat(c)      = spval
          cch4%ch4_ebul_depth_unsat(c,:)   = spval
          cch4%ch4_ebul_total_unsat(c)     = spval
          cch4%ch4_surf_ebul_unsat(c)      = spval
          cch4%ch4_surf_diff_unsat(c)      = spval
          cch4%ch4_dfsat_flux(c)           = spval
          cch4%zwt_ch4_unsat(c)            = spval
          cch4%fphr(c,:)                   = spval
          cch4%sif(c)                      = spval
          cch4%o2stress_unsat(c,:)         = spval
          cch4%ch4stress_unsat(c,:)        = spval
          cws%finundated(c)                = spval
          pps_a%rootfr(c,:)                = spval

       else  ! Inactive CH4 columns

          cch4%ch4_prod_depth_sat(c,:)     = spval
          cch4%ch4_prod_depth_unsat(c,:)   = spval
          cch4%ch4_prod_depth_lake(c,:)    = spval
          cch4%ch4_oxid_depth_sat(c,:)     = spval
          cch4%ch4_oxid_depth_unsat(c,:)   = spval
          cch4%ch4_oxid_depth_lake(c,:)    = spval
          cch4%o2_oxid_depth_sat(c,:)      = spval
          cch4%o2_oxid_depth_unsat(c,:)    = spval
          cch4%o2_decomp_depth_sat(c,:)    = spval
          cch4%o2_decomp_depth_unsat(c,:)  = spval
          cch4%o2_aere_depth_sat(c,:)      = spval
          cch4%o2_aere_depth_unsat(c,:)    = spval
          cch4%co2_decomp_depth_sat(c,:)   = spval
          cch4%co2_decomp_depth_unsat(c,:) = spval
          cch4%co2_oxid_depth_sat(c,:)     = spval
          cch4%co2_oxid_depth_unsat(c,:)   = spval
          cch4%ch4_aere_depth_sat(c,:)     = spval
          cch4%ch4_aere_depth_unsat(c,:)   = spval
          cch4%ch4_tran_depth_sat(c,:)     = spval
          cch4%ch4_tran_depth_unsat(c,:)   = spval
          cch4%co2_aere_depth_sat(c,:)     = spval
          cch4%co2_aere_depth_unsat(c,:)   = spval
          cch4%ch4_surf_aere_sat(c)        = spval
          cch4%ch4_surf_aere_unsat(c)      = spval
          cch4%ch4_ebul_depth_sat(c,:)     = spval
          cch4%ch4_ebul_depth_unsat(c,:)   = spval
          cch4%ch4_ebul_total_sat(c)       = spval
          cch4%ch4_ebul_total_unsat(c)     = spval
          cch4%ch4_surf_ebul_sat(c)        = spval
          cch4%ch4_surf_ebul_unsat(c)      = spval
          cch4%ch4_surf_ebul_lake(c)       = spval
          cch4%ch4_surf_diff_sat(c)        = spval
          cch4%ch4_surf_diff_unsat(c)      = spval
          cch4%ch4_surf_diff_lake(c)       = spval
          cch4%ch4_dfsat_flux(c)           = spval
          cch4%zwt_ch4_unsat(c)            = spval
          cch4%conc_ch4_lake(c,:)          = spval
          cch4%conc_o2_lake(c,:)           = spval
          cch4%fphr(c,:)                   = spval
          cch4%sif(c)                      = spval
          cch4%o2stress_unsat(c,:)         = spval
          cch4%o2stress_sat(c,:)           = spval
          cch4%ch4stress_unsat(c,:)        = spval
          cch4%ch4stress_sat(c,:)          = spval
          cws%finundated(c)                = spval
          pps_a%rootfr(c,:)                = spval
          pps_a%grnd_ch4_cond(c)           = spval
          ! totcolch4 Set to zero for inactive columns so that this can be used
          ! as an appropriate area-weighted gridcell average soil methane content.
          cch4%totcolch4(c)                = 0._r8  

       end if

    end do

  end subroutine initColdCH4

end module ch4InitMod
