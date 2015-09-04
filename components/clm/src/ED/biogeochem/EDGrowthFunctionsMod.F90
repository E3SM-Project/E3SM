module EDGrowthFunctionsMod

  ! ============================================================================
  ! Functions that control the trajectory of plant growth. 
  ! Ideally these would all use parameters that are fed in from the parameter file. 
  ! At present, there is only a single allocation trajectory. 
  ! ============================================================================

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use clm_varctl       , only : iulog 

  use EcophysConType   , only : ecophyscon
  use EDEcophysContype , only : EDecophyscon
  use EDtypesMod       , only : cohort, nlevcan_ed, dinc_ed

  implicit none
  save
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
    !  Height(m) diameter(cm) relationships. O'Brien et al  - for 56 pft at BCI                                    
    ! ============================================================================

    use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none

    type(cohort), intent(in) :: cohort_in

    !FIX(SPM,040214) - move to param file
    real(r8) :: m !parameter of allometric equation (needs to not be hardwired...
    real(r8) :: c !parameter of allometric equation (needs to not be hardwired...

    m = 0.64_r8
    c = 0.37_r8

    dbh = (10.0_r8**((log10(cohort_in%hite) - c)/m))

    return

  end function dbh

! ============================================================================

  real(r8) function Hite( cohort_in )

    ! ============================================================================
    !  Creates height in m as a function of diameter in cm. 
    !  Height(m) diameter(cm) relationships. O'Brien et al  - for 56 pft at BCI                                    
    ! ============================================================================

    use shr_kind_mod, only : r8 => shr_kind_r8;
    use clm_varctl  , only : iulog

    implicit none

    type(cohort), intent(inout) :: cohort_in

    real(r8) :: m
    real(r8) :: c
    real(r8) :: h

    m = 0.64_r8
    c = 0.37_r8       

    if(cohort_in%dbh <= 0._r8)then
       write(iulog,*) 'ED: dbh less than zero problem!',cohort_in%indexnumber
       cohort_in%dbh = 0.1_r8
    endif

    ! if the hite is larger than the maximum allowable height (set by dbhmax) then 
    ! set the height to the maximum value. 
    ! this could do with at least re-factoring and probably re-thinking. RF
    if(cohort_in%dbh <= EDecophyscon%max_dbh(cohort_in%pft)) then
       h = (10.0_r8**(log10(cohort_in%dbh) * m + c))
    else 
       h = (10.0_r8**(log10(EDecophyscon%max_dbh(cohort_in%pft))*m + c))
    endif
    Hite = h 

    return

  end function Hite

! ============================================================================

  real(r8) function Bleaf( cohort_in )

    ! ============================================================================
    !  Creates leaf biomass (kGC) as a function of tree diameter.  
    ! ============================================================================


    use shr_kind_mod, only : r8 => shr_kind_r8;
    use clm_varctl  , only : iulog

    implicit none

    type(cohort), intent(in) :: cohort_in       

    if(cohort_in%dbh < 0._r8.or.cohort_in%pft == 0.or.cohort_in%dbh > 1000_r8)then
       write(iulog,*) 'problems in bleaf',cohort_in%dbh,cohort_in%pft
    endif

    if(cohort_in%dbh <= EDecophyscon%max_dbh(cohort_in%pft))then
       bleaf = 0.0419_r8 * (cohort_in%dbh**1.56) * EDecophyscon%wood_density(cohort_in%pft)**0.55_r8
    else  
       bleaf = 0.0419_r8 * (EDecophyscon%max_dbh(cohort_in%pft)**1.56) * EDecophyscon%wood_density(cohort_in%pft)**0.55_r8       
    endif

    !Adjust for canopies that have become so deep that their bottom layer is not producing any carbon... 
    !nb this will change the allometry and the effects of this remain untested. RF. April 2014  
    bleaf = bleaf*cohort_in%canopy_trim

    return
  end function Bleaf

! ============================================================================

  real(r8) function tree_lai( cohort_in )

    ! ============================================================================
    !  LAI of individual trees is a function of the total leaf area and the total canopy area.   
    ! ============================================================================

    use shr_kind_mod, only : r8 => shr_kind_r8;
    use clm_varctl  , only : iulog

    implicit none

    type(cohort), intent(inout) :: cohort_in       

    real(r8) :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: slat               ! the sla of the top leaf layer. m2/kgC

    if( cohort_in%bl  <  0._r8 .or. cohort_in%pft  ==  0 ) then
       write(iulog,*) 'problem in treelai',cohort_in%bl,cohort_in%pft
    endif

    if( cohort_in%status_coh  ==  2 ) then ! are the leaves on? 
       slat = 1000.0_r8 * ecophyscon%slatop(cohort_in%pft) ! m2/g to m2/kg
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
    ! at the moments nlevcan_ed default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(cohort_in%treelai > nlevcan_ed*dinc_ed)then
       write(iulog,*) 'too much lai' , cohort_in%treelai , cohort_in%pft , nlevcan_ed * dinc_ed
    endif

    return

  end function tree_lai
  
  ! ============================================================================

  real(r8) function tree_sai( cohort_in )

    ! ============================================================================
    !  SAI of individual trees is a function of the total dead biomass per unit canopy area.   
    ! ============================================================================

    use shr_kind_mod, only : r8 => shr_kind_r8;
    use clm_varctl  , only : iulog

    implicit none

    type(cohort), intent(inout) :: cohort_in       

    real(r8) :: bdead_per_unitarea ! KgC of leaf per m2 area of ground.
    real(r8) :: sai_scaler     ! This is hardwired, but should be made a parameter  - 
             ! I need to add a new parameter to the 'standard' parameter file but don't have permission... RF 2 july.    

    sai_scaler = 0.05_r8 ! here, a high biomass of 20KgC per m2 gives us a high SAI of 1.0. 

    if( cohort_in%bdead  <  0._r8 .or. cohort_in%pft  ==  0 ) then
       write(iulog,*) 'problem in treesai',cohort_in%bdead,cohort_in%pft
    endif

    cohort_in%c_area = c_area(cohort_in) ! call the tree area 
    bdead_per_unitarea = cohort_in%bdead/(cohort_in%c_area/cohort_in%n) !KgC/m2
    tree_sai = bdead_per_unitarea * sai_scaler !kg/m2 * m2/kg = unitless LAI 
   
    cohort_in%treesai = tree_sai

    ! here, if the LAI exceeeds the maximum size of the possible array, then we have no way of accomodating it
    ! at the moments nlevcan_ed default is 40, which is very large, so exceeding this would clearly illustrate a 
    ! huge error 
    if(cohort_in%treesai > nlevcan_ed*dinc_ed)then
       write(iulog,*) 'too much sai' , cohort_in%treesai , cohort_in%pft , nlevcan_ed * dinc_ed
    endif

    return

  end function tree_sai
  

! ============================================================================

  real(r8) function c_area( cohort_in )

    ! ============================================================================
    ! Calculate area of ground covered by entire cohort. (m2)
    ! Function of DBH (cm) canopy spread (m/cm) and number of individuals. 
    ! ============================================================================


    use shr_kind_mod             , only : r8 => shr_kind_r8;
    use clm_varctl               , only : iulog
    use EDParamsMod              , only : ED_val_grass_spread

    implicit none

    type(cohort), intent(in) :: cohort_in       

    real(r8) :: dbh ! Tree diameter at breat height. cm. 

    if (DEBUG_growth) then
       write(iulog,*) 'z_area 1',cohort_in%dbh,cohort_in%pft
       write(iulog,*) 'z_area 2',EDecophyscon%max_dbh
       write(iulog,*) 'z_area 3',ecophyscon%woody
       write(iulog,*) 'z_area 4',cohort_in%n
       write(iulog,*) 'z_area 5',cohort_in%patchptr%spread
       write(iulog,*) 'z_area 6',cohort_in%canopy_layer
       write(iulog,*) 'z_area 7',ED_val_grass_spread
    end if

    dbh = min(cohort_in%dbh,EDecophyscon%max_dbh(cohort_in%pft))
    if(ecophyscon%woody(cohort_in%pft) == 1)then 
       c_area = 3.142_r8 * cohort_in%n * &
            (cohort_in%patchptr%spread(cohort_in%canopy_layer)*dbh)**1.56_r8
    else
       c_area = 3.142_r8 * cohort_in%n * (ED_val_grass_spread*dbh)**1.56_r8      
    end if

  end function c_area

! ============================================================================

  real(r8) function Bdead( cohort_in )

    ! ============================================================================
    ! Calculate stem biomass from height(m) dbh(cm) and wood density(g/cm3)
    ! using allometry of J.G. Saldarriaga et al 1988 - Rio Negro                                  
    ! Journal of Ecology vol 76 p938-958                                       
    ! ============================================================================

    use shr_kind_mod, only: r8 => shr_kind_r8;

    implicit none     

    type(cohort), intent(in) :: cohort_in       

    bdead = 0.06896_r8*(cohort_in%hite**0.572_r8)*(cohort_in%dbh**1.94_r8)* &
         (EDecophyscon%wood_density(cohort_in%pft)**0.931_r8)

  end function Bdead

! ============================================================================

  real(r8) function dHdBd( cohort_in )

    ! ============================================================================
    ! convert changes in structural biomass to changes in height                                        
    ! consistent with Bstem and h-dbh allometries                               
    ! ============================================================================

    use shr_kind_mod, only: r8 => shr_kind_r8;

    implicit none

    type(cohort), intent(in)  :: cohort_in

    real(r8) :: dbddh ! rate of change of dead biomass (KgC) per unit change of height (m) 

    dbddh = 0.06896_r8*0.572_r8*(cohort_in%hite**(-0.428_r8))*(cohort_in%dbh**1.94_r8)* &
         (EDecophyscon%wood_density(cohort_in%pft)**0.931_r8)
    dHdBd = 1.0_r8/dbddh !m/KgC 

    return

  end function dHdBd

! ============================================================================
  real(r8) function dDbhdBd( cohort_in )

    ! ============================================================================
    ! convert changes in structural biomass to changes in diameter                                        
    ! consistent with Bstem and h-dbh allometries                               
    ! ============================================================================

    use shr_kind_mod, only : r8 => shr_kind_r8;

    implicit none

    type(cohort), intent(in) :: cohort_in

    real(r8) :: dBD_dDBH !Rate of change of dead biomass (KgC) with change in DBH (cm) 
    real(r8) :: dH_dDBH  !Rate of change of height (m) with change in DBH (cm) 

    dBD_dDBH = 1.94_r8*0.06896_r8*(cohort_in%hite**0.572_r8)*(cohort_in%dbh**0.94_r8)* &
         (EDecophyscon%wood_density(cohort_in%pft)**0.931_r8)
    if(cohort_in%dbh < EDecophyscon%max_dbh(cohort_in%pft))then
       dH_dDBH = 1.4976_r8*(cohort_in%dbh**(-0.36_r8))
       dBD_dDBH = dBD_dDBH + 0.572_r8*0.06896_r8*(cohort_in%hite**(0.572_r8 - 1.0_r8))* &
            (cohort_in%dbh**1.94_r8)*(EDecophyscon%wood_density(cohort_in%pft)**0.931_r8)*dH_dDBH
    endif

    dDbhdBd = 1.0/dBD_dDBH

    return

  end function dDbhdBd

! ============================================================================

  real(r8) function dDbhdBl( cohort_in )

    ! ============================================================================
    ! convert changes in leaf biomass (KgC) to changes in DBH (cm)                                
    ! ============================================================================

    use shr_kind_mod, only: r8 => shr_kind_r8;

    implicit none

    type(cohort), intent(in) :: cohort_in

    real(r8) :: dblddbh ! Rate of change of leaf biomass with change in DBH

    dblddbh = 1.56_r8*0.0419_r8*(cohort_in%dbh**0.56_r8)*(EDecophyscon%wood_density(cohort_in%pft)**0.55_r8)
    dblddbh = dblddbh*cohort_in%canopy_trim

    dDbhdBl = 1.0_r8/dblddbh

    return

  end function dDbhdBl

! ============================================================================

  real(r8) function mortality_rates( cohort_in )

    ! ============================================================================
    !  Calculate mortality rates as a function of carbon storage       
    ! ============================================================================

    use shr_kind_mod, only : r8 => shr_kind_r8
    use EDParamsMod,  only : ED_val_stress_mort

    implicit none    

    type (cohort), intent(in) :: cohort_in

    real(r8) :: frac  ! relativised stored carbohydrate
    real(r8) :: smort ! stress mortality     : Fraction per year
    real(r8) :: bmort ! background mortality : Fraction per year

    ! 'Background' mortality (can vary as a function of density as in ED1.0 and ED2.0, but doesn't here for tractability) 
    bmort = 0.014_r8 
   
    ! Proxy for hydraulic failure induced mortality. 
    smort = 0.0_r8
    if(cohort_in%patchptr%btran_ft(cohort_in%pft) <= 0.000001_r8)then 
       smort = smort + ED_val_stress_mort 
    endif
    
    ! Carbon Starvation induced mortality.
    if ( cohort_in%dbh  >  0._r8 ) then
       if(Bleaf(cohort_in) > 0._r8.and.cohort_in%bstore <= Bleaf(cohort_in))then
          frac = cohort_in%bstore/(Bleaf(cohort_in))
          smort = smort + max(0.0_r8,ED_val_stress_mort*(1.0_r8 - frac))
       endif
    else
       write(iulog,*) 'dbh problem in mortality_rates', &
            cohort_in%dbh,cohort_in%pft,cohort_in%n,cohort_in%canopy_layer,cohort_in%indexnumber
    endif

    mortality_rates = smort + bmort

  end function mortality_rates

! ============================================================================

end module EDGrowthFunctionsMod
