module EDGrowthFunctionsMod

  ! ============================================================================
  ! Functions that control the trajectory of plant growth. 
  ! Ideally these would all use parameters that are fed in from the parameter file. 
  ! At present, there is only a single allocation trajectory. 
  ! ============================================================================

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals     , only : fates_log
  use EDPftvarcon        , only : EDPftvarcon_inst
  use EDTypesMod       , only : ed_cohort_type, nlevleaf, dinc_ed
  use FatesConstantsMod        , only : itrue,ifalse

  implicit none
  private

  public ::  bleaf
  public ::  hite
  public ::  ddbhdbd
  public ::  ddbhdbl
  public ::  dhdbd
  public ::  dbh
  public ::  bdead
  public ::  tree_lai
  public ::  tree_sai
  public ::  c_area
  public ::  mortality_rates

  logical :: DEBUG_growth = .false.

  ! ============================================================================
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  real(r8) function Dbh( cohort_in )

    ! ============================================================================
    !  Creates diameter in cm as a function of height in m
    !  Height(m) diameter(cm) relationships. O'Brien et al  - for 56 patch at BCI                                    
    ! ============================================================================

    type(ed_cohort_type), intent(in) :: cohort_in

    !FIX(SPM,040214) - move to param file
    real(r8) :: m ! parameter of allometric equation
    real(r8) :: c ! parameter of allometric equation

    m = EDPftvarcon_inst%allom_d2h1(cohort_in%pft)
    c = EDPftvarcon_inst%allom_d2h2(cohort_in%pft)

    dbh = (10.0_r8**((log10(cohort_in%hite) - c)/m))

    return

  end function dbh

! ============================================================================

  real(r8) function Hite( cohort_in )

    ! ============================================================================
    !  Creates height in m as a function of diameter in cm. 
    !  Height(m) diameter(cm) relationships. O'Brien et al  - for 56 pft at BCI                                    
    ! ============================================================================

    type(ed_cohort_type), intent(inout) :: cohort_in

    real(r8) :: m
    real(r8) :: c
    real(r8) :: h

    m = EDPftvarcon_inst%allom_d2h1(cohort_in%pft)
    c = EDPftvarcon_inst%allom_d2h2(cohort_in%pft)

    if(cohort_in%dbh <= 0._r8)then
       write(fates_log(),*) 'ED: dbh less than zero problem!'
       cohort_in%dbh = 0.1_r8
    endif

    ! if the hite is larger than the maximum allowable height (set by dbhmax) then 
    ! set the height to the maximum value. 
    ! this could do with at least re-factoring and probably re-thinking. RF
    if(cohort_in%dbh <= EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft)) then
       h = (10.0_r8**(log10(cohort_in%dbh) * m + c))
    else 
       h = (10.0_r8**(log10(EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft))*m + c))
    endif
    Hite = h 

    return

  end function Hite

! ============================================================================

  real(r8) function Bleaf( cohort_in )

    ! ============================================================================
    !  Creates leaf biomass (kGC) as a function of tree diameter.  
    ! ============================================================================

    type(ed_cohort_type), intent(in) :: cohort_in       
    
    real(r8) :: dbh2bl_a 
    real(r8) :: dbh2bl_b
    real(r8) :: dbh2bl_c

    dbh2bl_a  =  EDPftvarcon_inst%allom_d2bl1(cohort_in%pft)
    dbh2bl_b  =  EDPftvarcon_inst%allom_d2bl2(cohort_in%pft)
    dbh2bl_c  =  EDPftvarcon_inst%allom_d2bl3(cohort_in%pft)
    
    if(cohort_in%dbh < 0._r8.or.cohort_in%pft == 0.or.cohort_in%dbh > 1000.0_r8)then
       write(fates_log(),*) 'problems in bleaf',cohort_in%dbh,cohort_in%pft
    endif

    if(cohort_in%dbh <= EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft))then
       bleaf = dbh2bl_a * (cohort_in%dbh**dbh2bl_b) * EDPftvarcon_inst%wood_density(cohort_in%pft)**dbh2bl_c 
    else  
       bleaf = dbh2bl_a * (EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft)**dbh2bl_b) * &
            EDPftvarcon_inst%wood_density(cohort_in%pft)**dbh2bl_c
    endif  

    !write(fates_log(),*) 'bleaf',bleaf, slascaler,cohort_in%pft
    
    !Adjust for canopies that have become so deep that their bottom layer is not producing any carbon... 
    !nb this will change the allometry and the effects of this remain untested. RF. April 2014  
    
     bleaf = bleaf * cohort_in%canopy_trim

    return
  end function Bleaf

! ============================================================================

  real(r8) function tree_lai( cohort_in )

    ! ============================================================================
    !  LAI of individual trees is a function of the total leaf area and the total canopy area.   
    ! ============================================================================

    type(ed_cohort_type), intent(inout) :: cohort_in       

    real(r8) :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: slat               ! the sla of the top leaf layer. m2/kgC

    if( cohort_in%bl  <  0._r8 .or. cohort_in%pft  ==  0 ) then
       write(fates_log(),*) 'problem in treelai',cohort_in%bl,cohort_in%pft
    endif

    if( cohort_in%status_coh  ==  2 ) then ! are the leaves on? 
       slat = 1000.0_r8 * EDPftvarcon_inst%slatop(cohort_in%pft) ! m2/g to m2/kg
       cohort_in%c_area = c_area(cohort_in) ! call the tree area
       leafc_per_unitarea = cohort_in%bl/(cohort_in%c_area/cohort_in%n) !KgC/m2
       if(leafc_per_unitarea > 0.0_r8)then
          tree_lai = leafc_per_unitarea * slat  !kg/m2 * m2/kg = unitless LAI 
       else
          tree_lai = 0.0_r8
       endif
    else
       tree_lai = 0.0_r8
    endif !status
    cohort_in%treelai = tree_lai

    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevleaf default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(cohort_in%treelai > nlevleaf*dinc_ed)then
       write(fates_log(),*) 'too much lai' , cohort_in%treelai , cohort_in%pft , nlevleaf * dinc_ed
    endif

    return

  end function tree_lai
  
  ! ============================================================================

  real(r8) function tree_sai( cohort_in )

    ! ============================================================================
    !  SAI of individual trees is a function of the total dead biomass per unit canopy area.   
    ! ============================================================================

    type(ed_cohort_type), intent(inout) :: cohort_in       

    real(r8) :: bdead_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: sai_scaler     

    sai_scaler = EDPftvarcon_inst%allom_sai_scaler(cohort_in%pft) 

    if( cohort_in%bdead  <  0._r8 .or. cohort_in%pft  ==  0 ) then
       write(fates_log(),*) 'problem in treesai',cohort_in%bdead,cohort_in%pft
    endif

    cohort_in%c_area = c_area(cohort_in) ! call the tree area 
    bdead_per_unitarea = cohort_in%bdead/(cohort_in%c_area/cohort_in%n) !KgC/m2
    tree_sai = bdead_per_unitarea * sai_scaler !kg/m2 * m2/kg = unitless LAI 
   
    cohort_in%treesai = tree_sai

    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevleaf default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(cohort_in%treesai > nlevleaf*dinc_ed)then
       write(fates_log(),*) 'too much sai' , cohort_in%treesai , cohort_in%pft , nlevleaf * dinc_ed
    endif

    return

  end function tree_sai
  

! ============================================================================

  real(r8) function c_area( cohort_in )

    ! ============================================================================
    ! Calculate area of ground covered by entire cohort. (m2)
    ! Function of DBH (cm) canopy spread (m/cm) and number of individuals. 
    ! ============================================================================

    use EDTypesMod               , only : nclmax

    type(ed_cohort_type), intent(in) :: cohort_in       

    real(r8) :: dbh ! Tree diameter at breat height. cm. 
    real(r8) :: crown_area_to_dbh_exponent
    real(r8) :: spreadterm

    ! default is to use the same exponent as the dbh to bleaf exponent so that per-plant canopy depth remains invariant during growth,
    ! but allowed to vary via the allom_blca_expnt_diff term (which has default value of zero)
    crown_area_to_dbh_exponent = EDPftvarcon_inst%allom_d2bl2(cohort_in%pft) + &
          EDPftvarcon_inst%allom_blca_expnt_diff(cohort_in%pft)
    
    if (DEBUG_growth) then
       write(fates_log(),*) 'z_area 1',cohort_in%dbh,cohort_in%pft
       write(fates_log(),*) 'z_area 2',EDPftvarcon_inst%allom_dbh_maxheight
       write(fates_log(),*) 'z_area 3',EDPftvarcon_inst%woody
       write(fates_log(),*) 'z_area 4',cohort_in%n
       write(fates_log(),*) 'z_area 5',cohort_in%siteptr%spread
       write(fates_log(),*) 'z_area 6',cohort_in%canopy_layer
    end if
    
    dbh = min(cohort_in%dbh,EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft))
    
    ! ----------------------------------------------------------------------------------
    ! The function c_area is called during the process of canopy position demotion
    ! and promotion. As such, some cohorts are temporarily elevated to canopy positions
    ! that are outside the number of alloted canopy spaces.  Ie, a two story canopy
    ! may have a third-story plant, if only for a moment.  However, these plants
    ! still need to generate a crown area to complete the promotion, demotion process.
    ! So we allow layer index exceedence here and force it down to max.
    ! (rgk/cdk 05/2017)
    ! ----------------------------------------------------------------------------------
    
    ! apply site-level spread elasticity to the cohort crown allometry term
    spreadterm = cohort_in%siteptr%spread * EDPftvarcon_inst%allom_d2ca_coefficient_max(cohort_in%pft) + &
         (1._r8 - cohort_in%siteptr%spread) * EDPftvarcon_inst%allom_d2ca_coefficient_min(cohort_in%pft)
    !
    c_area = cohort_in%n * spreadterm * dbh ** crown_area_to_dbh_exponent

  end function c_area

! ============================================================================

  real(r8) function Bdead( cohort_in )

    ! ============================================================================
    ! Calculate stem biomass from height(m) dbh(cm) and wood density(g/cm3)
    ! default params using allometry of J.G. Saldarriaga et al 1988 - Rio Negro                                  
    ! Journal of Ecology vol 76 p938-958  
    !
    ! NOTE (RGK 07-2017) Various other biomass allometries calculate above ground
    ! biomass, and it appear Saldariagga may be an outlier that calculates total
    ! biomass (these parameters will have to be a placeholder for both)
    !
    ! ============================================================================

    type(ed_cohort_type), intent(in) :: cohort_in       

   real(r8) :: dbh2bd_a
   real(r8) :: dbh2bd_b
   real(r8) :: dbh2bd_c
   real(r8) :: dbh2bd_d
   
   dbh2bd_a =  EDPftvarcon_inst%allom_agb1(cohort_in%pft)
   dbh2bd_b =  EDPftvarcon_inst%allom_agb2(cohort_in%pft)
   dbh2bd_c =  EDPftvarcon_inst%allom_agb3(cohort_in%pft)  
   dbh2bd_d =  EDPftvarcon_inst%allom_agb4(cohort_in%pft)

   bdead = dbh2bd_a*(cohort_in%hite**dbh2bd_b)*(cohort_in%dbh**dbh2bd_c)* &
        (EDPftvarcon_inst%wood_density(cohort_in%pft)** dbh2bd_d)  

  end function Bdead

! ============================================================================

  real(r8) function dHdBd( cohort_in )

    ! ============================================================================
    ! convert changes in structural biomass to changes in height                                        
    ! consistent with Bstem and h-dbh allometries                               
    ! ============================================================================

    type(ed_cohort_type), intent(in)  :: cohort_in

    real(r8) :: dbddh ! rate of change of dead biomass (KgC) per unit change of height (m) 
    real(r8) :: dbh2bd_a
    real(r8) :: dbh2bd_b
    real(r8) :: dbh2bd_c
    real(r8) :: dbh2bd_d
    
    dbh2bd_a =  EDPftvarcon_inst%allom_agb1(cohort_in%pft)
    dbh2bd_b =  EDPftvarcon_inst%allom_agb2(cohort_in%pft)
    dbh2bd_c =  EDPftvarcon_inst%allom_agb3(cohort_in%pft)  
    dbh2bd_d =  EDPftvarcon_inst%allom_agb4(cohort_in%pft)
    
    dbddh =  dbh2bd_a*dbh2bd_b*(cohort_in%hite**(dbh2bd_b-1.0_r8))*(cohort_in%dbh**dbh2bd_c)* &
         (EDPftvarcon_inst%wood_density(cohort_in%pft)**dbh2bd_d)
    dHdBd = 1.0_r8/dbddh !m/KgC 
    
    return

  end function dHdBd

! ============================================================================
  real(r8) function dDbhdBd( cohort_in )

    ! ============================================================================
    ! convert changes in structural biomass to changes in diameter                                        
    ! consistent with Bstem and h-dbh allometries                               
    ! ============================================================================

    type(ed_cohort_type), intent(in) :: cohort_in 

    real(r8) :: dBD_dDBH !Rate of change of dead biomass (KgC) with change in DBH (cm) 
    real(r8) :: dH_dDBH  !Rate of change of height (m) with change in DBH (cm) 
    real(r8) :: m
    real(r8) :: c
    real(r8) :: h
    real(r8) :: dbh2bd_a
    real(r8) :: dbh2bd_b
    real(r8) :: dbh2bd_c
    real(r8) :: dbh2bd_d

    m = EDPftvarcon_inst%allom_d2h1(cohort_in%pft)
    c = EDPftvarcon_inst%allom_d2h2(cohort_in%pft)
    
    dbh2bd_a =  EDPftvarcon_inst%allom_agb1(cohort_in%pft)
    dbh2bd_b =  EDPftvarcon_inst%allom_agb2(cohort_in%pft)
    dbh2bd_c =  EDPftvarcon_inst%allom_agb3(cohort_in%pft)  
    dbh2bd_d =  EDPftvarcon_inst%allom_agb4(cohort_in%pft)
    
    dBD_dDBH =  dbh2bd_c*dbh2bd_a*(cohort_in%hite**dbh2bd_b)*(cohort_in%dbh**(dbh2bd_c-1.0_r8))* &
         (EDPftvarcon_inst%wood_density(cohort_in%pft)**dbh2bd_d)  

    if(cohort_in%dbh < EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft))then
       dH_dDBH = (10.0_r8**c)*m*(cohort_in%dbh**(m-1.0_r8))          

       dBD_dDBH =  dBD_dDBH + dbh2bd_b*dbh2bd_a*(cohort_in%hite**(dbh2bd_b - 1.0_r8))* &
            (cohort_in%dbh**dbh2bd_c)*(EDPftvarcon_inst%wood_density(cohort_in%pft)**dbh2bd_d)*dH_dDBH 
    endif

    dDbhdBd = 1.0_r8/dBD_dDBH

    return

  end function dDbhdBd

! ============================================================================

  real(r8) function dDbhdBl( cohort_in )

    ! ============================================================================
    ! convert changes in leaf biomass (KgC) to changes in DBH (cm)                                
    ! ============================================================================

    type(ed_cohort_type), intent(in) :: cohort_in

    real(r8) :: dblddbh ! Rate of change of leaf biomass with change in DBH
    real(r8) :: dbh2bl_a
    real(r8) :: dbh2bl_b
    real(r8) :: dbh2bl_c
    
    dbh2bl_a =  EDPftvarcon_inst%allom_d2bl1(cohort_in%pft)
    dbh2bl_b =  EDPftvarcon_inst%allom_d2bl2(cohort_in%pft)
    dbh2bl_c =  EDPftvarcon_inst%allom_d2bl3(cohort_in%pft)


    dblddbh = dbh2bl_b*dbh2bl_a*(cohort_in%dbh**dbh2bl_b)*(EDPftvarcon_inst%wood_density(cohort_in%pft)**dbh2bl_c)
    dblddbh = dblddbh*cohort_in%canopy_trim

    if( cohort_in%dbh<EDPftvarcon_inst%allom_dbh_maxheight(cohort_in%pft) ) then
        dDbhdBl = 1.0_r8/dblddbh
    else
        dDbhdBl = 1.0d15  ! At maximum size, the leaf biomass is saturated, dbl=0
    endif

    return

  end function dDbhdBl

! ============================================================================

  subroutine mortality_rates( cohort_in,cmort,hmort,bmort )

    ! ============================================================================
    !  Calculate mortality rates as a function of carbon storage       
    ! ============================================================================

    use EDParamsMod,  only : ED_val_stress_mort
    use FatesInterfaceMod,  only : hlm_use_ed_prescribed_phys

    type (ed_cohort_type), intent(in) :: cohort_in
    real(r8),intent(out) :: bmort ! background mortality : Fraction per year
    real(r8),intent(out) :: cmort  ! carbon starvation mortality
    real(r8),intent(out) :: hmort  ! hydraulic failure mortality

    real(r8) :: frac  ! relativised stored carbohydrate

    real(r8) :: hf_sm_threshold    ! hydraulic failure soil moisture threshold 


    if (hlm_use_ed_prescribed_phys .eq. ifalse) then

    ! 'Background' mortality (can vary as a function of density as in ED1.0 and ED2.0, but doesn't here for tractability) 
    bmort = EDPftvarcon_inst%bmort(cohort_in%pft) 

    ! Proxy for hydraulic failure induced mortality. 
    hf_sm_threshold = EDPftvarcon_inst%hf_sm_threshold(cohort_in%pft)

    if(cohort_in%patchptr%btran_ft(cohort_in%pft) <= hf_sm_threshold)then 
       hmort = ED_val_stress_mort
     else
       hmort = 0.0_r8
     endif 
    
    ! Carbon Starvation induced mortality.
    if ( cohort_in%dbh  >  0._r8 ) then
       if(Bleaf(cohort_in) > 0._r8 .and. cohort_in%bstore <= Bleaf(cohort_in))then
          frac = cohort_in%bstore/(Bleaf(cohort_in))
          cmort = max(0.0_r8,ED_val_stress_mort*(1.0_r8 - frac))
        else
          cmort = 0.0_r8
       endif

    else
       write(fates_log(),*) 'dbh problem in mortality_rates', &
            cohort_in%dbh,cohort_in%pft,cohort_in%n,cohort_in%canopy_layer
    endif

    !mortality_rates = bmort + hmort + cmort

    else ! i.e. hlm_use_ed_prescribed_phys is true
       if ( cohort_in%canopy_layer .eq. 1) then
          bmort = EDPftvarcon_inst%prescribed_mortality_canopy(cohort_in%pft)
       else
          bmort = EDPftvarcon_inst%prescribed_mortality_understory(cohort_in%pft)
       endif
       cmort = 0._r8
       hmort = 0._r8
    endif

 end subroutine mortality_rates

! ============================================================================

end module EDGrowthFunctionsMod
