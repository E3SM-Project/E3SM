module CNDVEstablishmentMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates establishment of new pfts
  ! Called once per year
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils  , only: endrun
  use decompMod   , only : bounds_type
  use shr_log_mod , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Establishment
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Establishment(bounds)
    !
    ! !DESCRIPTION:
    ! Calculates establishment of new pfts
    ! Called once per year
    !
    ! !USES:
    use clmtype
    use clm_varpar   , only : numpft
    use clm_varcon   , only : istsoil
    use clm_varctl   , only : iulog
    use pftvarcon    , only : noveg, nc3_arctic_grass
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_PI, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,p,m                        ! indices
    integer  :: fn, filterg(bounds%begg-bounds%endg+1)         ! local gridcell filter for error check
    ! gridcell level variables
    integer  :: ngrass(bounds%begg:bounds%endg)                ! counter
    integer  :: npft_estab(bounds%begg:bounds%endg)            ! counter
    real(r8) :: fpc_tree_total(bounds%begg:bounds%endg)        ! total fractional cover of trees in vegetated portion of gridcell
    real(r8) :: fpc_total(bounds%begg:bounds%endg)             ! old-total fractional vegetated portion of gridcell (without bare ground)
    real(r8) :: fpc_total_new(bounds%begg:bounds%endg)         ! new-total fractional vegetated portion of gridcell (without bare ground)
    ! pft level variables
    logical  :: survive(bounds%begp:bounds%endp)               ! true=>pft survives
    logical  :: estab(bounds%begp:bounds%endp)                 ! true=>pft is established
    real(r8) :: dstemc(bounds%begp:bounds%endp)                ! local copy of deadstemc
    ! local and temporary variables or parameters
    real(r8) :: taper          ! ratio of height:radius_breast_height (tree allometry)
    real(r8) :: estab_rate     ! establishment rate
    real(r8) :: estab_grid     ! establishment rate on grid cell
    real(r8) :: fpcgridtemp    ! temporary
    real(r8) :: stemdiam       ! stem diameter
    real(r8) :: stocking       ! #stems / ha (stocking density)
    real(r8) :: lai_ind        ! LAI per individual
    real(r8) :: lm_ind         ! leaf carbon (gC/ind)
    real(r8) :: fpc_ind        ! individual foliage projective cover
    real(r8)::  bm_delta
    ! parameters
    real(r8), parameter :: ramp_agddtw = 300.0
    ! minimum individual density for persistence of PFT (indiv/m2)
    real(r8), parameter :: nind_min = 1.0e-10_r8
    ! minimum precip. for establishment (mm/s)
    real(r8), parameter :: prec_min_estab = 100._r8/(365._r8*SHR_CONST_CDAY)
    ! maximum sapling establishment rate (indiv/m2)
    real(r8), parameter :: estab_max = 0.24_r8
    !-----------------------------------------------------------------------


   associate(& 
   agdd20                              =>    pdgvs%agdd20                                , & ! Input:  [real(r8) (:)]  20-yr running mean of agdd                        
   tmomin20                            =>    pdgvs%tmomin20                              , & ! Input:  [real(r8) (:)]  20-yr running mean of tmomin                      
   present                             =>    pdgvs%present                               , & ! InOut:  [logical (:)]  true=> PFT present in patch                        
   nind                                =>    pdgvs%nind                                  , & ! InOut:  [real(r8) (:)]  number of individuals (#/m**2)                    
   fpcgrid                             =>    pdgvs%fpcgrid                               , & ! Output: [real(r8) (:)]  foliar projective cover on gridcell (fraction)    
   crownarea                           =>    pdgvs%crownarea                             , & ! Output: [real(r8) (:)]  area that each individual tree takes up (m^2)     
   greffic                             =>    pdgvs%greffic                               , & ! Output: [real(r8) (:)]  lpj's growth efficiency                           
   heatstress                          =>    pdgvs%heatstress                            , & ! Output: [real(r8) (:)]                                                    
   annsum_npp                          =>    pepv%annsum_npp                             , & ! Input:  [real(r8) (:)]  annual sum NPP (gC/m2/yr)                         
   annsum_litfall                      =>    pepv%annsum_litfall                         , & ! Input:  [real(r8) (:)]  annual sum litfall (gC/m2/yr)                     
   prec365                             =>    pdgvs%prec365                               , & ! Input:  [real(r8) (:)]  365-day running mean of tot. precipitation        
   agddtw                              =>    pdgvs%agddtw                                , & ! Input:  [real(r8) (:)]  accumulated growing degree days above twmax       
   pftmayexist                         =>    pdgvs%pftmayexist                           , & ! Input:  [logical (:)]  exclude seasonal decid pfts from tropics [1=true, 0=false]
   crownarea_max                       =>    dgv_pftcon%crownarea_max                    , & ! Input:  [real(r8) (:)]  ecophys const - tree maximum crown area [m2]      
   twmax                               =>    dgv_pftcon%twmax                            , & ! Input:  [real(r8) (:)]  ecophys const - upper limit of temperature of the warmest month
   reinickerp                          =>    dgv_pftcon%reinickerp                       , & ! Input:  [real(r8) (:)]  ecophys const - parameter in allometric equation  
   allom1                              =>    dgv_pftcon%allom1                           , & ! Input:  [real(r8) (:)]  ecophys const - parameter in allometric           
   tcmax                               =>    dgv_pftcon%tcmax                            , & ! Input:  [real(r8) (:)]  ecophys const - maximum coldest monthly mean temperature
   tcmin                               =>    dgv_pftcon%tcmin                            , & ! Input:  [real(r8) (:)]  ecophys const - minimum coldest monthly mean temperature
   gddmin                              =>    dgv_pftcon%gddmin                           , & ! Input:  [real(r8) (:)]  ecophys const - minimum growing degree days (at or above 5 C)
   leafcmax                            =>    pcs%leafcmax                                , & ! Input:  [real(r8) (:)]  (gC/m2) ann max leaf C                            
   deadstemc                           =>    pcs%deadstemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   slatop                              =>    pftcon%slatop                               , & ! Input:  [real(r8) (:)] specific leaf area at top of canopy, projected area basis [m^2/gC]
   dsladlai                            =>    pftcon%dsladlai                             , & ! Input:  [real(r8) (:)] dSLA/dLAI, projected area basis [m^2/gC]           
   dwood                               =>    pftcon%dwood                                , & ! Input:  [real(r8) (:)]  ecophys const - wood density (gC/m3)              
   woody                               =>    pftcon%woody                                  & ! Input:  [real(r8) (:)]  ecophys const - woody pft or not                  
   )

    ! **********************************************************************
    ! Slevis version of LPJ's subr. bioclim
    ! Limits based on 20-year running averages of coldest-month mean
    ! temperature and growing degree days (5 degree base).
    ! For SURVIVAL, coldest month temperature and GDD should be
    ! at least as high as PFT-specific limits.
    ! For REGENERATION, PFT must be able to survive AND coldest month
    ! temperature should be no higher than a PFT-specific limit.
    ! **********************************************************************

    taper = 200._r8 ! make a global constant as with dwood (lpj's wooddens)

    ! Initialize gridcell-level metrics

    do g = bounds%begg,bounds%endg
       ngrass(g) = 0
       npft_estab(g) = 0
       fpc_tree_total(g) = 0._r8
       fpc_total(g) = 0._r8
       fpc_total_new(g) = 0._r8
    end do

    do p = bounds%begp,bounds%endp
       ! Set the presence of pft for this gridcell

       if (nind(p) == 0._r8) present(p) = .false.
       if (.not. present(p)) then
          nind(p) = 0._r8
          fpcgrid(p) = 0._r8
       end if
       survive(p) = .false.
       estab(p)   = .false.
       dstemc(p)  = deadstemc(p)
    end do

    ! Must go thru all 16 pfts and decide which can/cannot establish or survive
    ! Determine present, survive, estab.  Note: Even if tmomin20>tcmax, crops
    ! and 2nd boreal summergreen tree cannot exist (see
    ! EcosystemDynini) because this model cannot simulate such pfts, yet.
    ! Note - agddtw is only defined at the pft level and has now been moved
    ! to an if-statement below to determine establishment of boreal trees

    do p = bounds%begp,bounds%endp
       if (tmomin20(p) >= tcmin(pft%itype(p)) + SHR_CONST_TKFRZ ) then
          if (tmomin20(p) <= tcmax(pft%itype(p)) + SHR_CONST_TKFRZ  .and. agdd20(p) >= gddmin(pft%itype(p))) then
             estab(p) = .true.
          end if
          survive(p) = .true.
          ! seasonal decid. pfts that would have occurred in regions without
          ! short winter day lengths (see CNPhenology)
          if (.not. pftmayexist(p)) then
             survive(p) = .false.
             estab(p) = .false.
             pftmayexist(p) = .true.
          end if
       end if
    end do

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)

       ! Case 1 -- pft ceases to exist -kill pfts not adapted to current climate

       if (present(p) .and. (.not. survive(p) .or. nind(p)<nind_min)) then
          present(p) = .false.
          fpcgrid(p) = 0._r8
          nind(p) = 0._r8
       end if

       ! Case 2 -- pft begins to exist - introduce newly "adapted" pfts

       if (lun%itype(l) == istsoil) then
          if (.not. present(p) .and. prec365(p) >= prec_min_estab .and. estab(p)) then
             if (twmax(pft%itype(p)) > 999._r8 .or. agddtw(p) == 0._r8) then

                present(p) = .true.
                nind(p) = 0._r8
                ! lpj starts with fpcgrid=0 and calculates
                ! seed fpcgrid from the carbon of saplings;
                ! with CN we need the seed fpcgrid up front
                ! to scale seed leafc to lm_ind to get fpcgrid;
                ! sounds circular; also seed fpcgrid depends on sla,
                ! so theoretically need diff value for each pft;slevis
                fpcgrid(p) = 0.000844_r8
                if (woody(pft%itype(p)) < 1._r8) then
                   fpcgrid(p) = 0.05_r8
                end if

                ! Seed carbon for newly established pfts
                ! Equiv. to pleaf=1 & pstor=1 set in subr pftwt_cnbal (slevis)
                ! ***Dangerous*** to hardwire leafcmax here; find alternative!
                ! Consider just assigning nind and fpcgrid for newly
                ! established pfts instead of entering the circular procedure
                ! outlined in the paragraph above
                leafcmax(p) = 1._r8
                if (dstemc(p) <= 0._r8) dstemc(p) = 0.1_r8

             end if   ! conditions required for establishment
          end if   ! conditions required for establishment
       end if   ! if soil

       ! Case 3 -- some pfts continue to exist (no change) and some pfts
       ! continue to not exist (no change). Do nothing for this case.

    end do

    ! Sapling and grass establishment
    ! Calculate total woody FPC, FPC increment and grass cover (= crown area)
    ! Calculate total woody FPC and number of woody PFTs present and able to establish

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       if (present(p)) then
          if (woody(pft%itype(p)) == 1._r8) then
             fpc_tree_total(g) = fpc_tree_total(g) + fpcgrid(p)
             if (estab(p)) npft_estab(g) = npft_estab(g) + 1
          else if (woody(pft%itype(p)) < 1._r8 .and. pft%itype(p) > noveg) then !grass
             ngrass(g) = ngrass(g) + 1
          end if
       end if
    end do

    ! Above grid-level establishment counters are required for the next steps.

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)

       if (present(p) .and. woody(pft%itype(p)) == 1._r8 .and. estab(p)) then

          ! Calculate establishment rate over available space, per tree PFT
          ! Max establishment rate reduced by shading as tree FPC approaches 1
          ! Total establishment rate partitioned equally among regenerating woody PFTs

          estab_rate = estab_max * (1._r8-exp(5._r8*(fpc_tree_total(g)-1._r8))) / real(npft_estab(g))

          ! Calculate grid-level establishment rate per woody PFT
          ! Space available for woody PFT establishment is fraction of grid cell
          ! not currently occupied by woody PFTs

          estab_grid = estab_rate * (1._r8-fpc_tree_total(g))

          ! Add new saplings to current population

          nind(p) = nind(p) + estab_grid

          !slevis: lpj's lm_ind was the max leaf mass for the year;
          !now lm_ind is the max leaf mass for the year calculated in CNFire
          !except when a pft is newly established (nind==0); then lm_ind
          !is assigned a leafcmax above
 
          lm_ind = leafcmax(p) * fpcgrid(p) / nind(p) ! nind>0 for sure
          if (fpcgrid(p) > 0._r8 .and. nind(p) > 0._r8) then
             stocking = nind(p)/fpcgrid(p) !#ind/m2 nat veg area -> #ind/m2 pft area
             ! stemdiam derived here from cn's formula for htop found in
             ! CNVegStructUpdate and cn's assumption stemdiam=2*htop/taper
             ! this derivation neglects upper htop limit enforced elsewhere
             stemdiam = (24._r8 * dstemc(p) / (SHR_CONST_PI * stocking * dwood(pft%itype(p)) * taper))**(1._r8/3._r8)
          else
             stemdiam = 0._r8
          end if
          ! Eqn D (now also in Light; need here for 1st yr when pfts haven't established, yet)
          crownarea(p) = min(crownarea_max(pft%itype(p)), allom1(pft%itype(p))*stemdiam**reinickerp(pft%itype(p)))

          ! Update LAI and FPC

          if (crownarea(p) > 0._r8) then
             if (dsladlai(pft%itype(p)) > 0._r8) then
                ! make lai_ind >= 0.001 to avoid killing plants at this stage
                lai_ind = max(0.001_r8,((exp(lm_ind*dsladlai(pft%itype(p)) + log(slatop(pft%itype(p)))) - &
                              slatop(pft%itype(p)))/dsladlai(pft%itype(p))) / crownarea(p))
             else ! currently redundant because dsladlai=0 for grasses only
                lai_ind = lm_ind * slatop(pft%itype(p)) / crownarea(p) ! lpj's formula
             end if
          else
             lai_ind = 0._r8
          end if

          fpc_ind = 1._r8 - exp(-0.5_r8*lai_ind)
          fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind

       end if   ! add new saplings block
       if (present(p) .and. woody(pft%itype(p)) == 1._r8) then
          fpc_total_new(g) = fpc_total_new(g) + fpcgrid(p)
       end if
    end do   ! close loop to update fpc_total_new

    ! Adjustments- don't allow trees to exceed 95% of vegetated landunit

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       if (fpc_total_new(g) > 0.95_r8) then
          if (woody(pft%itype(p)) == 1._r8 .and. present(p)) then
             nind(p) = nind(p) * 0.95_r8 / fpc_total_new(g)
             fpcgrid(p) = fpcgrid(p) * 0.95_r8 / fpc_total_new(g)
          end if
          fpc_total(g) = 0.95_r8

       else
          fpc_total(g) = fpc_total_new(g)
       end if
    end do

    ! Section for grasses. Grasses can establish in non-vegetated areas

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       if (present(p) .and. woody(pft%itype(p)) < 1._r8) then
          if (leafcmax(p) <= 0._r8 .or. fpcgrid(p) <= 0._r8 ) then
             present(p) = .false.
             nind(p) = 0._r8
          else
             nind(p) = 1._r8 ! in case these grasses just established
             crownarea(p) = 1._r8
             lm_ind = leafcmax(p) * fpcgrid(p) / nind(p)
             if (dsladlai(pft%itype(p)) > 0._r8) then
                lai_ind = max(0.001_r8,((exp(lm_ind*dsladlai(pft%itype(p)) + log(slatop(pft%itype(p)))) - &
                              slatop(pft%itype(p)))/dsladlai(pft%itype(p))) / crownarea(p))
             else ! 'if' is currently redundant b/c dsladlai=0 for grasses only
                lai_ind = lm_ind * slatop(pft%itype(p)) / crownarea(p)
             end if
             fpc_ind = 1._r8 - exp(-0.5_r8*lai_ind)
             fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind
             fpc_total(g) = fpc_total(g) + fpcgrid(p)
          end if
       end if
    end do   ! end of pft-loop

    ! Adjustment of fpc_total > 1 due to grasses (pft%itype >= nc3_arctic_grass)

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)

       if (fpc_total(g) > 1._r8) then
          if (pft%itype(p) >= nc3_arctic_grass .and. fpcgrid(p) > 0._r8) then
             fpcgridtemp = fpcgrid(p)
             fpcgrid(p) = max(0._r8, fpcgrid(p) - (fpc_total(g)-1._r8))
             fpc_total(g) = fpc_total(g) - fpcgridtemp + fpcgrid(p)
          end if
       end if

       ! Remove tiny fpcgrid amounts

       if (fpcgrid(p) < 1.e-15_r8) then
          fpc_total(g) = fpc_total(g) - fpcgrid(p)
          fpcgrid(p) = 0._r8
          present(p) = .false.
          nind(p) = 0._r8
       end if

       ! Set the fpcgrid for bare ground if there is bare ground in
       ! vegetated landunit and pft is bare ground so that everything
       ! can add up to one.

       if (fpc_total(g) < 1._r8 .and. pft%itype(p) == noveg) then
          fpcgrid(p) = 1._r8 - fpc_total(g)
          fpc_total(g) = fpc_total(g) + fpcgrid(p)
       end if

    end do

    ! Annual calculations used hourly in GapMortality
    ! Ultimately may wish to place in separate subroutine...

    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)

       ! Stress mortality from lpj's subr Mortality

       if (woody(pft%itype(p)) == 1._r8 .and. nind(p) > 0._r8 .and. &
           leafcmax(p) > 0._r8 .and. fpcgrid(p) > 0._r8) then

          if (twmax(pft%itype(p)) < 999._r8) then
             heatstress(p) = max(0._r8, min(1._r8, agddtw(p) / ramp_agddtw))
          else
             heatstress(p) = 0._r8
          end if

          ! Net individual living biomass increment
          ! NB: lpj's turnover not exactly same as cn's litfall:
          ! lpj's sap->heartwood turnover not included in litfall (slevis)

          bm_delta = max(0._r8, annsum_npp(p) - annsum_litfall(p))
          lm_ind = leafcmax(p) * fpcgrid(p) / nind(p)

          ! Growth efficiency (net biomass increment per unit leaf area)

          if (dsladlai(pft%itype(p)) > 0._r8) then
             greffic(p) = bm_delta / (max(0.001_r8,                     &
                 ( ( exp(lm_ind*dsladlai(pft%itype(p)) + log(slatop(pft%itype(p)))) &
                     - slatop(pft%itype(p)) ) / dsladlai(pft%itype(p)) )))
          else ! currently redundant because dsladlai=0 for grasses only
             greffic(p) = bm_delta / (lm_ind * slatop(pft%itype(p)))
          end if
       else
          greffic(p) = 0.
          heatstress(p) = 0.
       end if

    end do

    ! Check for error in establishment
    fn = 0
    do g = bounds%begg,bounds%endg
       if (abs(fpc_total(g) - 1._r8) > 1.e-6) then
          fn = fn + 1
          filterg(fn) = g
       end if
    end do
    ! Just print out the first error
    if (fn > 0) then
       g = filterg(1)
       write(iulog,*) 'Error in Establishment: fpc_total =',fpc_total(g), ' at gridcell ',g
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    end associate 
   end subroutine Establishment

end module CNDVEstablishmentMod
