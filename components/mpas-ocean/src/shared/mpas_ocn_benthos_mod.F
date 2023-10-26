! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_benthos_mod
!
!> \brief MPAS ocean benthos
!> \author Nicole Jeffery
!> \date   07/15/2020
!> \details
!------------------------------------------------------------------------------
!> This subroutine computes the biogeochemistry of a ~30 cm sea-floor sediment
!> layer resolved at the cm scale.  The 35 component coupled-biogeochemistry
!> is based on Reed et al 2011 and Krumins et al 2013.  The vertical transport scheme
!> is similar to Jeffery et al. 2011 for sea ice and modified for a  fixed depth
!> porous sediment layer.
!------------------------------------------------------------------------------

 module ocn_benthos_mod

   use ocn_benthos_parms

   implicit none

   private
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      benthos_init,           &
      benthos_flags_init,     &
      benthos_SurfaceFluxes,  &
      benthos_SourceSink

!-----------------------------------------------------------------------
!  module variables
!-----------------------------------------------------------------------

   integer (KIND=benthos_i4), public, parameter :: &
      benthos_tracer_cnt = 35, &
      element_tracer_cnt = 7

   integer (KIND=benthos_i4), public :: &
      nBenthicTracers, &
      nElements

   integer (KIND=benthos_i4), public :: &
      poca_ind, &
      pocb_ind, &
      pocc_ind, &
      pona_ind, &
      ponb_ind, &
      ponc_ind, &
      popa_ind, &
      popb_ind, &
      popc_ind, &
      o2_ind, &
      nh4_ind, &
      h2po4_ind, &
      co2_ind, &
      no3_ind, &
      mno2a_ind, &
      mno2b_ind, &
      mn_ind, &
      feoh3a_ind, &
      feoh3b_ind, &
      fe_ind, &
      fepa_ind, &
      fepb_ind, &
      so4_ind, &
      h2s_ind, &
      ch4_ind, &
      dic_ind, &
      alk_ind, &
      hco3_ind, &
      fes_ind, &
      fes2_ind, &
      s_ind, &
      caco3a_ind, &
      caco3b_ind, &
      co3_ind, &
      camgco3_ind, &
      totc_ind, &
      toto_ind, &
      totn_ind, &
      tots_ind, &
      totp_ind, &
      totmn_ind, &
      totfe_ind

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: benthos_init
! !INTERFACE:

 subroutine benthos_init(benthos_indices, benthos_element_indices, &
      nBenthicVertLevels, nBenthicTracers_in, nElements_in, benthosInterfaceLayerThickness, &
      benthosVertGridPoints)

! !DESCRIPTION:
!  Initialize tracegas tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
     nBenthicVertLevels, &
     nBenthicTracers_in, &
     nElements_in

! !INPUT/OUTPUT PARAMETERS:

  type(benthos_indices_type), intent(inout) :: benthos_indices
  type(benthos_element_indices_type), intent(inout) :: benthos_element_indices

! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(out) :: &
       benthosInterfaceLayerThickness, & ! (m)
       benthosVertGridPoints

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  !NJ-TEST
  logical(KIND=benthos_log) :: test_flag
  integer(KIND=benthos_i4) :: iBenthicTracers, iSecondaryReactions, iElements, iBenthicVertLevels
  !NJ-END
!-----------------------------------------------------------------------
!  Local copy of indices
!-----------------------------------------------------------------------

   poca_ind    = benthos_indices%poca_ind
   pocb_ind    = benthos_indices%pocb_ind
   pocc_ind    = benthos_indices%pocc_ind
   pona_ind    = benthos_indices%pona_ind
   ponb_ind    = benthos_indices%ponb_ind
   ponc_ind    = benthos_indices%ponc_ind
   popa_ind    = benthos_indices%popa_ind
   popb_ind    = benthos_indices%popb_ind
   popc_ind    = benthos_indices%popc_ind
   o2_ind      = benthos_indices%o2_ind
   nh4_ind     = benthos_indices%nh4_ind
   h2po4_ind   = benthos_indices%h2po4_ind
   co2_ind     = benthos_indices%co2_ind
   no3_ind     = benthos_indices%no3_ind
   mno2a_ind   = benthos_indices%mno2a_ind
   mno2b_ind   = benthos_indices%mno2b_ind
   mn_ind      = benthos_indices%mn_ind
   feoh3a_ind  = benthos_indices%feoh3a_ind
   feoh3b_ind  = benthos_indices%feoh3b_ind
   fe_ind      = benthos_indices%fe_ind
   fepa_ind    = benthos_indices%fepa_ind
   fepb_ind    = benthos_indices%fepb_ind
   so4_ind     = benthos_indices%so4_ind
   h2s_ind     = benthos_indices%h2s_ind
   ch4_ind     = benthos_indices%ch4_ind
   dic_ind     = benthos_indices%dic_ind
   alk_ind     = benthos_indices%alk_ind
   hco3_ind    = benthos_indices%hco3_ind
   fes_ind     = benthos_indices%fes_ind
   fes2_ind    = benthos_indices%fes2_ind
   s_ind       = benthos_indices%s_ind
   caco3a_ind  = benthos_indices%caco3a_ind
   caco3b_ind  = benthos_indices%caco3b_ind
   co3_ind     = benthos_indices%co3_ind
   camgco3_ind = benthos_indices%camgco3_ind

   ! element indices
   totc_ind = benthos_element_indices%carbon_ind
   toto_ind = benthos_element_indices%oxygen_ind
   totn_ind = benthos_element_indices%nitrogen_ind
   totp_ind = benthos_element_indices%phosphorus_ind
   tots_ind = benthos_element_indices%sulfur_ind
   totmn_ind = benthos_element_indices%manganese_ind
   totfe_ind = benthos_element_indices%iron_ind
!-------------------------------------------------
! initialize tracer and element numbers
!-------------------------------------------------

   nElements = nElements_in
   nBenthicTracers = nBenthicTracers_in

   if (nBenthicTracers.ne.benthos_tracer_cnt) then
      write(*,*) 'error: nBenthicTracers is inconsistent with assumed value'
      write(*,*) 'nBenthicTracers:', nBenthicTracers
      write(*,*) 'benthos_tracer_cnt:', benthos_tracer_cnt
   end if
!-----------------------------------------------------------------------
!  Define benthos tracers
!-----------------------------------------------------------------------

   benthos_indices%units(:) = 'mmol/m3'
   benthos_element_indices%units(:) = 'mmol/m3'

   benthos_indices%short_name(poca_ind)='POCa'
   benthos_indices%long_name(poca_ind)='Particulate labile organic carbon'
   benthos_indices%units(poca_ind)='mmol/kg'

   benthos_indices%short_name(pocb_ind)='POCb'
   benthos_indices%long_name(pocb_ind)='Particulate semi-labile organic carbon'
   benthos_indices%units(pocb_ind)='mmol/kg'
   benthos_indices%short_name(pocc_ind)='POCc'
   benthos_indices%long_name(pocc_ind)='Particulate refractory organic carbon'
   benthos_indices%units(pocc_ind)='mmol/kg'

   benthos_indices%short_name(pona_ind)='PONa'
   benthos_indices%long_name(pona_ind)='Particulate labile organic nitrogen'
   benthos_indices%units(pona_ind)='mmol/kg'

   benthos_indices%short_name(ponb_ind)='PONb'
   benthos_indices%long_name(ponb_ind)='Particulate semi-labile organic nitrogen'
   benthos_indices%units(ponb_ind)='mmol/kg'

   benthos_indices%short_name(ponc_ind)='PONc'
   benthos_indices%long_name(ponc_ind)='Particulate refractory organic nitrogen'
   benthos_indices%units(ponc_ind)='mmol/kg'

   benthos_indices%short_name(popa_ind)='POPa'
   benthos_indices%long_name(popa_ind)='Particulate labile organic phosphorus'
   benthos_indices%units(popa_ind)='mmol/kg'

   benthos_indices%short_name(popb_ind)='POPb'
   benthos_indices%long_name(popb_ind)='Particulate semi-labile organic phosphorus'
   benthos_indices%units(popb_ind)='mmol/kg'

   benthos_indices%short_name(popc_ind)='POPc'
   benthos_indices%long_name(popc_ind)='Particulate refractory organic phosphorus'
   benthos_indices%units(popc_ind)='mmol/kg'

   benthos_indices%short_name(o2_ind)='O2'
   benthos_indices%long_name(o2_ind)='Dissolved oxygen'

   benthos_indices%short_name(nh4_ind)='NH4'
   benthos_indices%long_name(nh4_ind)='Dissolved ammonium'

   benthos_indices%short_name(h2po4_ind)='H2PO4'
   benthos_indices%long_name(h2po4_ind)='Dissolved phosphate'

   benthos_indices%short_name(co2_ind)='CO2'
   benthos_indices%long_name(co2_ind)='Dissolved carbon dioxide'

   benthos_indices%short_name(no3_ind)='NO3'
   benthos_indices%long_name(no3_ind)='Dissolved nitrate'

   benthos_indices%short_name(mno2a_ind)='MnO2a'
   benthos_indices%long_name(mno2a_ind)='Particulate amorphous manganese oxyhydroxide'
   benthos_indices%units(mno2a_ind)='mmol/kg'

   benthos_indices%short_name(mno2b_ind)='MnO2b'
   benthos_indices%long_name(mno2b_ind)='Particulate crystalline manganese oxyhydroxide'
   benthos_indices%units(mno2b_ind)='mmol/kg'

   benthos_indices%short_name(mn_ind)='Mn'
   benthos_indices%long_name(mn_ind)='Dissolved manganese'

   benthos_indices%short_name(feoh3a_ind)='Fe(OH)3a'
   benthos_indices%long_name(feoh3a_ind)='Particulate amorphous iron oxyhydroxide'
   benthos_indices%units(feoh3a_ind)='mmol/kg'

   benthos_indices%short_name(feoh3b_ind)='Fe(OH)3b'
   benthos_indices%long_name(feoh3b_ind)='Particulate crystalline iron oxyhydroxide'
   benthos_indices%units(feoh3b_ind)='mmol/kg'

   benthos_indices%short_name(fe_ind)='Fe(II)'
   benthos_indices%long_name(fe_ind)='Dissolved iron'

   benthos_indices%short_name(fepa_ind)='[Fe-P]a'
   benthos_indices%long_name(fepa_ind)='particulate amorphous iron bound phosphate'
   benthos_indices%units(fepa_ind)='mmol/kg'

   benthos_indices%short_name(fepb_ind)='[Fe-P]b'
   benthos_indices%long_name(fepb_ind)='particulate crystalline iron bound phosphate'
   benthos_indices%units(fepb_ind)='mmol/kg'

   benthos_indices%short_name(so4_ind)='SO4'
   benthos_indices%long_name(so4_ind)='Dissolved sulfate'

   benthos_indices%short_name(h2s_ind)='H2S'
   benthos_indices%long_name(h2s_ind)='Dissolved hydrogen sulfide'

   benthos_indices%short_name(ch4_ind)='CH4'
   benthos_indices%long_name(ch4_ind)='Dissolved methane'

   benthos_indices%short_name(dic_ind)='DIC'
   benthos_indices%long_name(dic_ind)='Dissolved inorganic carbon'

   benthos_indices%short_name(alk_ind)='ALK'
   benthos_indices%long_name(alk_ind)='Total alkalinity'
   benthos_indices%units(fepa_ind)='meq/m3'

   benthos_indices%short_name(hco3_ind)='HCO3'
   benthos_indices%long_name(hco3_ind)='Dissolved bicarbonate'

   benthos_indices%short_name(fes_ind)='FeS'
   benthos_indices%long_name(fes_ind)='Partculate iron monosulfide'
   benthos_indices%units(fes_ind)='mmol/kg'

   benthos_indices%short_name(fes2_ind)='FeS2'
   benthos_indices%long_name(fes2_ind)='Particulate pyrite'
   benthos_indices%units(fes2_ind)='mmol/kg'

   benthos_indices%short_name(s_ind)='S'
   benthos_indices%long_name(s_ind)='Particulate elemental sulfur'
   benthos_indices%units(s_ind)='mmol/kg'

   benthos_indices%short_name(caco3a_ind)='CaCO3a'
   benthos_indices%long_name(caco3a_ind)='Particulate calcite'
   benthos_indices%units(caco3a_ind)='mmol/kg'

   benthos_indices%short_name(caco3b_ind)='CaCO3b'
   benthos_indices%long_name(caco3b_ind)='Particulate aragonite'
   benthos_indices%units(caco3b_ind)='mmol/kg'

   benthos_indices%short_name(co3_ind)='CO3'
   benthos_indices%long_name(co3_ind)='Dissolved carbonate'

   benthos_indices%short_name(camgco3_ind)='Ca.85Mg.15CO3'
   benthos_indices%long_name(camgco3_ind)='Particulate 15% magnesium calcite'
   benthos_indices%units(camgco3_ind)='mmol/kg'

   benthos_element_indices%short_name(totc_ind)='Carbon'
   benthos_element_indices%long_name(totc_ind)='Total carbon'

   benthos_element_indices%short_name(toto_ind)='Oxygen'
   benthos_element_indices%long_name(toto_ind)='Total oxygen'

   benthos_element_indices%short_name(totn_ind)='Nitrogen'
   benthos_element_indices%long_name(totn_ind)='Total nitrogen'

   benthos_element_indices%short_name(totp_ind)='Phosphorus'
   benthos_element_indices%long_name(totp_ind)='Total phosphorus'

   benthos_element_indices%short_name(tots_ind)='Sulfur'
   benthos_element_indices%long_name(tots_ind)='Total sulfur'

   benthos_element_indices%short_name(totmn_ind)='Manganese'
   benthos_element_indices%long_name(totmn_ind)='Total manganese'

   benthos_element_indices%short_name(totfe_ind)='Iron'
   benthos_element_indices%long_name(totfe_ind)='Total iron'

   ! initialize transport parameters

   call benthos_parameters_init (nBenthicVertLevels, nBenthicTracers, nElements) ! now variable grid

   do iBenthicVertLevels = 1, nBenthicVertLevels+1
      benthosInterfaceLayerThickness(iBenthicVertLevels) =  vertBenthosGridThickI(iBenthicVertLevels)*benthosDepth
      benthosVertGridPoints(iBenthicVertLevels) = benthosGridPointsI(iBenthicVertLevels)
   end do

   !NJ-TEST
   !write(*,*) 'after benthos_parameters_init'
   !NJ-END

   !NJ-TEST: VERIFIED
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'diagnostics: benthos_parameters_init'
      do iBenthicVertLevels = 1, nBenthicVertLevels+2
         if (iBenthicVertLevels .lt. nBenthicVertLevels+2) then

            write(*,*) 'benthos_parameters_init variables: iBenthicVertLevels:',iBenthicVertLevels
            write(*,*) 'benthosMidPointI(iBenthicVertLevels):',benthosMidPointI(iBenthicVertLevels)
            write(*,*) 'benthosGridPointsI(iBenthicVertLevels):',benthosGridPointsI(iBenthicVertLevels)
            write(*,*) 'benthosPorosityI(iBenthicVertLevels):',benthosPorosityI(iBenthicVertLevels)
            write(*,*) 'benthosSolidPorosityI(iBenthicVertLevels):',benthosSolidPorosityI(iBenthicVertLevels)
            write(*,*) 'benthosDiffusivityI(iBenthicVertLevels):',benthosDiffusivityI(iBenthicVertLevels)
            write(*,*) 'benthosSolidDiffusivityI(iBenthicVertLevels):',benthosSolidDiffusivityI(iBenthicVertLevels)

            if (iBenthicVertLevels .lt. nBenthicVertLevels+1) write(*,*) 'vertBenthosGridThickI(iBenthicVertLevels):',vertBenthosGridThickI(iBenthicVertLevels)

            write(*,*) ' Grid point variables '
            write(*,*) 'benthosGridPoints(iBenthicVertLevels):',benthosGridPoints(iBenthicVertLevels)
            write(*,*) 'benthosPorosity(iBenthicVertLevels):',benthosPorosity(iBenthicVertLevels)
            write(*,*) 'benthosSolidPorosity(iBenthicVertLevels):',benthosSolidPorosity(iBenthicVertLevels)
            write(*,*) 'benthosDiffusivity(iBenthicVertLevels):',benthosDiffusivity(iBenthicVertLevels)
            write(*,*) 'benthosSolidDiffusivity(iBenthicVertLevels):',benthosSolidDiffusivity(iBenthicVertLevels)
         else

            write(*,*) 'iBenthicVertLevels:',iBenthicVertLevels
            write(*,*) 'benthosGridPoints(iBenthicVertLevels):',benthosGridPoints(iBenthicVertLevels)
            write(*,*) 'benthosPorosity(iBenthicVertLevels):',benthosPorosity(iBenthicVertLevels)
            write(*,*) 'benthosSolidPorosity(iBenthicVertLevels):',benthosSolidPorosity(iBenthicVertLevels)
            write(*,*) 'benthosDiffusivity(iBenthicVertLevels):',benthosDiffusivity(iBenthicVertLevels)
            write(*,*) 'benthosSolidDiffusivity(iBenthicVertLevels):',benthosSolidDiffusivity(iBenthicVertLevels)
         end if
      end do
   end if
   !  NJ:VERIFIED
   test_flag = .false.
   if (test_flag) then
      do iElements = 1,nElements
         do iBenthicTracers = 1,nBenthicTracers
            write(*,*) 'iElements, iBenthicTracers:',iElements, iBenthicTracers
            write(*,*) 'elementRatio:', elementRatios(iBenthicTracers, iElements)
         end do
      end do

   end if
   !NJ-END

   ! define stoichiometry of secondary reactions

   call secondaryStoichMatrix

   !NJ-TEST : VERIFIED
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'diagnostics: secondaryStoichMatrix'
      do iBenthicTracers = 1,nBenthicTracers
         do iSecondaryReactions = 1,nSecondaryReactions
            write(*,*) 'iBenthicTracers, iSecondaryReactions):',iBenthicTracers,iSecondaryReactions
            write(*,*) 'secondarySourceStoich(iBenthicTracers,iSecondaryReactions):',secondarySourceStoich(iBenthicTracers,iSecondaryReactions)
            write(*,*) 'secondarySinkStoich(iBenthicTracers,iSecondaryReactions):',secondarySinkStoich(iBenthicTracers,iSecondaryReactions)
         end do
      end do
   end if
   !NJ-END

   ! define stoichiometry of carbonate reactions

   call carbonateStoichMatrix

   !NJ-TEST : VERIFIED
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'diagnostics: carbonatesStoichMatrix'
      do iBenthicTracers = 1,nBenthicTracers
         do iSecondaryReactions = 1,nCarbonateReactions
            write(*,*) 'iBenthicTracers, iSecondaryReactions:',iBenthicTracers,iSecondaryReactions
            write(*,*) 'carbonateSourceStoich(iBenthicTracers,iSecondaryReactions):',carbonateSourceStoich(iBenthicTracers,iSecondaryReactions)
            write(*,*) 'carbonateSinkStoich(iBenthicTracers,iSecondaryReactions):',carbonateSinkStoich(iBenthicTracers,iSecondaryReactions)
         end do
      end do
      write(*,*) 'end benthos_init'
   end if
!   write(*,*) 'end benthos_init'
   !NJ-END

!-----------------------------------------------------------------------
!EOC

 end subroutine benthos_init

!*****************************************************************************
!BOP
! !IROUTINE: benthos_parameters_init
! !INTERFACE:

 subroutine benthos_parameters_init (nBenthicVertLevels, nBenthicTracers, nElements)


! !DESCRIPTION:
!  Specify parameter profiles for diffusivity, porosity and the vertical grid
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:

   integer (KIND=benthos_i4), intent(in) :: nBenthicVertLevels, nBenthicTracers, nElements

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
   !-----------------------------------------------------------------------

   integer(KIND=benthos_i4) :: &
        k, &
        iElements, &
        iLevels, &
        iTracers

   real(KIND=benthos_r8) :: &
        dx    , &
        Dm    , & ! molecular diffusion in m2/s
        iDin_o    ! maximum biodiffusion in m2/s

   real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
        tortuosity, &
        dz, &
        bioGrid

   logical(KIND=benthos_log) :: &
        test_flag

   integer(KIND=benthos_i4) :: &
        n1,n2, &   ! grid levels where transition between high and low resolutions
        n3         ! begins and ends

   real(KIND=benthos_r8), parameter :: &
        z1 = 0.001_benthos_r8, &  ! (m) high resolution spacing (1 mm upper layers)
        z2 = 0.02_benthos_r8      ! (m) low resolution spacing (2 cm deep bottom layers)

   ! allocate parameter functions

   allocate(vertBenthosGridThickI(nBenthicVertLevels+1))
   allocate(vertBenthosGridThick(nBenthicVertLevels))
   allocate(benthosMidPointI(nBenthicVertLevels+1))
   allocate(benthosGridPointsI(nBenthicVertLevels+1))
   allocate(benthosGridPoints(nBenthicVertLevels+2))
   allocate(benthosPorosity(nBenthicVertLevels+2))
   allocate(benthosPorosityI(nBenthicVertLevels+1))
   allocate(benthosSolidPorosity(nBenthicVertLevels+2))
   allocate(benthosSolidPorosityI(nBenthicVertLevels+1))
   allocate(benthosDiffusivity(nBenthicVertLevels+2))
   allocate(benthosDiffusivityI(nBenthicVertLevels+1))
   allocate(benthosSolidDiffusivity(nBenthicVertLevels+2))
   allocate(benthosSolidDiffusivityI(nBenthicVertLevels+1))
   allocate(benthosMolecularDiffusivity(nBenthicTracers))
   allocate(benthosSoluteDiffusivityI(nBenthicTracers,nBenthicVertLevels+1))
   allocate(benthosSoluteDiffusivity(nBenthicTracers,nBenthicVertLevels+2))
   allocate(transport_index(nBenthicTracers))

   ! allocate primary reaction matrices

   allocate(primarySourceStoich(nBenthicTracers,nPrimaryReactions))
   allocate(primarySinkStoich(nBenthicTracers,nPrimaryReactions))

   ! allocate element ratios for tracer to net element conversion

   allocate(elementRatios(nBenthicTracers,nElements))

   elementRatios(:,:) = c0_benthos

   ! carbon (mmol/ m^3), solid units  (mmol/kg)
   elementRatios(poca_ind,1) = c1_benthos*sediment_density
   elementRatios(pocb_ind,1) = c1_benthos*sediment_density
   elementRatios(pocc_ind,1) = c1_benthos*sediment_density
   elementRatios(dic_ind,1) = c1_benthos
   !elementRatios(co3_ind,1) = c1_benthos
   !elementRatios(co2_ind,1) = c1_benthos
   !elementRatios(hco3_ind,1) = c1_benthos
   elementRatios(caco3a_ind,1) = c1_benthos*sediment_density
   elementRatios(caco3b_ind,1) = c1_benthos*sediment_density
   elementRatios(camgco3_ind,1) = c1_benthos*sediment_density
   elementRatios(ch4_ind,1) = c1_benthos

  ! oxygen (mmol/m^3)
   elementRatios(poca_ind,2) = c1_benthos*sediment_density
   elementRatios(pocb_ind,2) = c1_benthos*sediment_density
   elementRatios(pocc_ind,2) = c1_benthos*sediment_density
   elementRatios(popa_ind,2) = 4.0_benthos_r8*sediment_density
   elementRatios(popb_ind,2) = 4.0_benthos_r8*sediment_density
   elementRatios(popc_ind,2) = 4.0_benthos_r8*sediment_density
   elementRatios(o2_ind,2) = c2_benthos
   elementRatios(h2po4_ind,2) = 4.0_benthos_r8
   elementRatios(co2_ind,2) = c2_benthos
   elementRatios(no3_ind,2) = 3.0_benthos_r8
   elementRatios(mno2a_ind,2) = c2_benthos*sediment_density
   elementRatios(mno2b_ind,2) = c2_benthos*sediment_density
   elementRatios(feoh3a_ind,2) = 3.0_benthos_r8*sediment_density
   elementRatios(feoh3b_ind,2) = 3.0_benthos_r8*sediment_density
   elementRatios(so4_ind,2) = 4.0_benthos_r8
   elementRatios(hco3_ind,2) = 3.0_benthos_r8
   elementRatios(caco3a_ind,2) = 3.0*sediment_density
   elementRatios(caco3b_ind,2) = 3.0*sediment_density
   elementRatios(co3_ind,2) = 3.0_benthos_r8

  ! nitrogen (mmol/m^3)
   elementRatios(nh4_ind,3)  = c1_benthos
   elementRatios(pona_ind,3) = c1_benthos*sediment_density
   elementRatios(ponb_ind,3) = c1_benthos*sediment_density
   elementRatios(ponc_ind,3) = c1_benthos*sediment_density
   elementRatios(no3_ind,3)  = c1_benthos

! phosphorus (mmol/m^3)
   elementRatios(popa_ind,4)  = c1_benthos*sediment_density
   elementRatios(popb_ind,4)  = c1_benthos*sediment_density
   elementRatios(popc_ind,4)  = c1_benthos*sediment_density
   elementRatios(h2po4_ind,4) = c1_benthos
   elementRatios(fepa_ind,4)  = c1_benthos*sediment_density
   elementRatios(fepb_ind,4)  = c1_benthos*sediment_density

! sulfur (mmol/m3)
   elementRatios(so4_ind,5)  = c1_benthos
   elementRatios(h2s_ind,5)  = c1_benthos
   elementRatios(fes_ind,5)  = c1_benthos*sediment_density
   elementRatios(fes2_ind,5) = c2_benthos*sediment_density
   elementRatios(s_ind,5)    = c1_benthos*sediment_density

! manganese (mmol/m^3)
   elementRatios(mno2a_ind,6) = c1_benthos*sediment_density
   elementRatios(mno2b_ind,6) = c1_benthos*sediment_density
   elementRatios(mn_ind,6)    = c1_benthos

   ! iron (mmol/m3)
   elementRatios(feoh3a_ind,7) = c1_benthos*sediment_density
   elementRatios(feoh3b_ind,7) = c1_benthos*sediment_density
   elementRatios(fe_ind,7)     = c1_benthos
   elementRatios(fes_ind,7)    = c1_benthos*sediment_density
   elementRatios(fes2_ind,7)   = c1_benthos*sediment_density

   ! define the tracers which are solved in the vertical transport
   transport_index(:) = 0
   transport_index(1) = poca_ind
   transport_index(2) = pocb_ind
   transport_index(3) = pocc_ind
   transport_index(4) = pona_ind
   transport_index(5) = ponb_ind
   transport_index(6) = ponc_ind
   transport_index(7) = popa_ind
   transport_index(8) = popb_ind
   transport_index(9) = popc_ind
   transport_index(10) = o2_ind
   transport_index(11) = nh4_ind
   transport_index(12) = h2po4_ind
   transport_index(13) = no3_ind
   transport_index(14) = mno2a_ind
   transport_index(15) = mno2b_ind
   transport_index(16) = mn_ind
   transport_index(17) = feoh3a_ind
   transport_index(18) = feoh3b_ind
   transport_index(19) = fe_ind
   transport_index(20) = fepa_ind
   transport_index(21) = fepb_ind
   transport_index(22) = so4_ind
   transport_index(23) = h2s_ind
   transport_index(24) = ch4_ind
   transport_index(25) = dic_ind
   transport_index(26) = alk_ind
   transport_index(27) = fes_ind
   transport_index(28) = fes2_ind
   transport_index(29) = s_ind
   transport_index(30) = caco3a_ind
   transport_index(31) = caco3b_ind
   transport_index(32) = camgco3_ind

   benthosMolecularDiffusivity(:) = c0_benthos
   
   ! Define the benthos variable resolution vertical grid
   n1 = 5
   n2 = nBenthicVertLevels-n1
   n3 = n2-n1-1

   dx = (z2-z1)/real(n3,KIND=benthos_r8)
   dz(1) = z1
   vertBenthosGridThickI(:) = c0_benthos
   vertBenthosGridThick(:) = c0_benthos
   do iLevels = 1,nBenthicVertLevels+1
      if (iLevels .gt. 1) dz(iLevels) = dz(iLevels-1) + dx
      benthosDiffusivityI(iLevels)         = c0_benthos
      benthosSolidDiffusivityI(iLevels)    = c0_benthos
      benthosSolidDiffusivity(iLevels)     = c0_benthos
      benthosDiffusivity(iLevels)          = c0_benthos

      do iTracers = 1,nBenthicTracers
         benthosSoluteDiffusivityI(iTracers,iLevels) = c0_benthos
         benthosSoluteDiffusivity(iTracers,iLevels)  = c0_benthos
         benthosSoluteDiffusivity(iTracers,iLevels+1)= c0_benthos
      end do
   end do

   benthosSolidDiffusivity(nBenthicVertLevels+2)     = c0_benthos
   benthosDiffusivity(nBenthicVertLevels+2)          = c0_benthos

   !benthosStorageConc = linspace(1,nBenthicTracers,nBenthicTracers)*0   !!!  (mmol/m3)

   bioGrid(1) = c0_benthos

   do k = 2, n1
      bioGrid(k) = bioGrid(k-1) + z1
   end do
   do k = n1+1, n2
      bioGrid(k) = bioGrid(k-1) + dz(k-n1)
   end do
   n3 = nBenthicVertLevels-n2+1
   dx = (benthosDepth - bioGrid(n2))/real(n3,KIND=benthos_r8)

   do k = n2+1, nBenthicVertLevels
      bioGrid(k) = bioGrid(k-1) + dx
   end do

   bioGrid(nBenthicVertLevels+1) = benthosDepth
   benthosGridPointsI(nBenthicVertLevels+1) = c1_benthos

   do k = 1, nBenthicVertLevels
      benthosGridPointsI(k) = bioGrid(k)/benthosDepth
   end do

   benthosGridPoints(1) = c0_benthos
   benthosGridPoints(nBenthicVertLevels+2) = c1_benthos
   do k = 2, nBenthicVertLevels+1
      benthosGridPoints(k) = (benthosGridPointsI(k) + benthosGridPointsI(k-1))/c2_benthos
      vertBenthosGridThickI(k-1) = benthosGridPoints(k) - benthosGridPoints(k-1)
      vertBenthosGridThick(k-1) = benthosGridPointsI(k)-benthosGridPointsI(k-1)
   end do
   vertBenthosGridThickI(nBenthicVertLevels+1) = benthosGridPoints(nBenthicVertLevels+2) - benthosGridPoints(nBenthicVertLevels+1)

   do k = 1, nBenthicVertLevels+1
      benthosMidPointI(k) = benthosDepth*benthosGridPoints(k+1)
   enddo

   ! Initialize porosity profile (Reed et al 2012)

   if (useDepthDependentPorosity) then

      do k =1, nBenthicVertLevels+1
         benthosPorosityI(k) = porosity_inf + (porosity_o-porosity_inf)*exp(-benthosGridPointsI(k)*benthosDepth/porosity_gamma)
      enddo

      benthosPorosity(1) = benthosPorosityI(1)
      benthosPorosity(nBenthicVertLevels+2) = benthosPorosityI(nBenthicVertLevels+1)

      do k = 2, nBenthicVertLevels+1
         benthosPorosity(k) = (benthosPorosityI(k-1) + benthosPorosityI(k))/c2_benthos
      end do

   else
      benthosPorosityI(:) = p5_benthos
      benthosPorosity(:)  = p5_benthos
   endif

   ! Define  solid porosity

   do k = 1, nBenthicVertLevels+1
      benthosSolidPorosityI(k) = c1_benthos - benthosPorosityI(k)
      benthosSolidPorosity(k) = c1_benthos - benthosPorosity(k)

      ! Tortuosity

      tortuosity(k) = c1_benthos - c2_benthos*log(benthosPorosityI(k))
   enddo

   benthosSolidPorosity(1) = benthosSolidPorosityI(1)
   benthosSolidPorosity(nBenthicVertLevels+2) = benthosSolidPorosity(nBenthicVertLevels+1)

   ! Define diffusivity (m2/s)

   Dm = molecular_diff/sec_per_year
   iDin_o = max_bio_diff*m2percm2/sec_per_year

   if (useNonZeroDiffusivity) then
      if (useConstantDiffusivity) then
         benthosDiffusivityI(:)      = iDin_o
         benthosSolidDiffusivityI(:) = iDin_o
         benthosDiffusivity(:)       = iDin_o
         benthosSolidDiffusivity(:)  = iDin_o
      else
         benthosDiffusivityI(1) = benthosPorosityI(1)*(Dm/tortuosity(1)+iDin_o)/benthosDepth**2
         benthosSolidDiffusivityI(1) = benthosSolidPorosityI(1)*iDin_o/benthosDepth**2

         do k = 2,nBenthicVertLevels+1
            benthosDiffusivityI(k) = benthosPorosityI(k)*(Dm/tortuosity(k) + iDin_o *exp(-benthosGridPointsI(k)*benthosDepth/x_biodiffusion))/benthosDepth**2
            benthosSolidDiffusivityI(k) =benthosSolidPorosityI(k)*(iDin_o *exp(-benthosGridPointsI(k)*benthosDepth/x_biodiffusion))/benthosDepth**2
         enddo

         benthosDiffusivity(1) = (Dm+iDin_o)/benthosDepth**2
         benthosDiffusivity(2) = benthosDiffusivityI(1)
         benthosDiffusivity(nBenthicVertLevels+2) = benthosDiffusivityI(nBenthicVertLevels+1)
         benthosDiffusivity(nBenthicVertLevels+1) = benthosDiffusivityI(nBenthicVertLevels+1)

         benthosSolidDiffusivity(nBenthicVertLevels+1) = benthosSolidDiffusivityI(nBenthicVertLevels+1)
         benthosSolidDiffusivity(nBenthicVertLevels+2) = benthosSolidDiffusivityI(nBenthicVertLevels+1)

         benthosSolidDiffusivity(1) = c0_benthos !benthosSolidDiffusivityI(1)
         benthosSolidDiffusivity(2) = benthosSolidDiffusivityI(1)

         do k = 3,nBenthicVertLevels
            benthosDiffusivity(k) = (benthosDiffusivityI(k-1) + benthosDiffusivityI(k))/c2_benthos
            benthosSolidDiffusivity(k) = (benthosSolidDiffusivityI(k-1) + benthosSolidDiffusivityI(k))/c2_benthos
         enddo

     endif   ! constant diffusivity
  endif
  test_flag = .false.
  if (test_flag) then
     do k = 1,nBenthicVertLevels+1
        write(*,*) 'benthosDiffusivity(k):', benthosDiffusivity(k)
        write(*,*) 'benthosSolidDiffusivity(k):', benthosSolidDiffusivity(k)
        write(*,*) 'benthosPorosity(k):', benthosPorosity(k)
        write(*,*) 'benthosSolidPorosity(k):', benthosSolidPorosity(k)
     end do
  end if

!
!  Initialize the deep benthic storage concentrations

! benthosStorageConc(:) = 0

!EOC

 end subroutine benthos_parameters_init

!*****************************************************************************
!BOP
! !IROUTINE: benthos_parameters_dt
! !INTERFACE:

 subroutine benthos_parameters_dt (nBenthicVertLevels,oceanBottomTemperature, oceanBottomSalinity, &
      oceanBottomPressure, benthos_indices)

! !DESCRIPTION:
!  Update time-dependent parameters
!
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:

   integer (KIND=benthos_i4), intent(in) :: nBenthicVertLevels

   real(KIND=benthos_r8), intent(in) :: &
        oceanBottomTemperature, & ! oC
        oceanBottomSalinity, &
        oceanBottomPressure

   type(benthos_indices_type), intent(in) :: benthos_indices

   !INOUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer(KIND=benthos_i4) :: &
        iBenthicTracers, &
        iBenthicLevels, &
        k, &
        iElements

   real(KIND=benthos_r8) :: &
        dz    , &
        Dm    , & ! molecular diffusion in m2/s
        iDin_o    ! maximum biodiffusion in m2/s

   real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
        tortuosity

   logical(KIND=benthos_log) :: &
        test_flag

   do iBenthicLevels = 1, nBenthicVertLevels+1

      benthosDiffusivityI(iBenthicLevels)        = c0_benthos
      benthosSolidDiffusivityI(iBenthicLevels)   = c0_benthos
      benthosDiffusivity(iBenthicLevels)         = c0_benthos
      benthosSolidDiffusivity(iBenthicLevels)    = c0_benthos

      ! Tortuosity
      tortuosity(iBenthicLevels) = c1_benthos - c2_benthos*log(benthosPorosityI(iBenthicLevels))
   enddo
   benthosDiffusivity(nBenthicVertLevels+2)         = c0_benthos
   benthosSolidDiffusivity(nBenthicVertLevels+2)    = c0_benthos

   ! Define diffusivity (m2/s)  NJ-Validated
   call benthos_molecular_diff (oceanBottomTemperature,oceanBottomSalinity, &
        oceanBottomPressure, benthosMolecularDiffusivity, benthos_indices)

   ! constant molecular_diff for solutes (for testing)
   Dm = molecular_diff/sec_per_year
   iDin_o = max_bio_diff*m2percm2/sec_per_year*abs((oceanBottomTemperature+degreeC_to_K)/T_max_biodiffusion)

   if (useNonZeroDiffusivity) then
      if (useConstantDiffusivity) then

         do iBenthicLevels = 1, nBenthicVertLevels+1
            benthosDiffusivityI(iBenthicLevels)      = iDin_o
            benthosSolidDiffusivityI(iBenthicLevels) = iDin_o
            benthosDiffusivity(iBenthicLevels)       = iDin_o
            benthosSolidDiffusivity(iBenthicLevels)  = iDin_o
         end do
         benthosDiffusivity(nBenthicVertLevels+2)       = iDin_o
         benthosSolidDiffusivity(nBenthicVertLevels+2)  = iDin_o
      else

      do  iBenthicTracers = 1,nBenthicTracers
         benthosSoluteDiffusivityI(iBenthicTracers,1) = &
              benthosPorosityI(1)*(benthosMolecularDiffusivity(iBenthicTracers)/tortuosity(1)+iDin_o)/benthosDepth**2
         do iBenthicLevels = 2,nBenthicVertLevels+1

            benthosSoluteDiffusivityI(iBenthicTracers,iBenthicLevels) = &
                 benthosPorosityI(iBenthicLevels)*(benthosMolecularDiffusivity(iBenthicTracers)/&
                 tortuosity(iBenthicLevels) + iDin_o *exp(-benthosGridPointsI(iBenthicLevels)* &
                 benthosDepth/x_biodiffusion))/benthosDepth**2

         end do

         benthosSoluteDiffusivity(iBenthicTracers,1) =(benthosMolecularDiffusivity(iBenthicTracers)+ &
              iDin_o)/benthosDepth**2
         benthosSoluteDiffusivity(iBenthicTracers,nBenthicVertLevels+2) = &
              benthosSoluteDiffusivityI(iBenthicTracers,nBenthicVertLevels+1)
         benthosSoluteDiffusivity(iBenthicTracers,2) = benthosSoluteDiffusivityI(iBenthicTracers,1)
         benthosSoluteDiffusivity(iBenthicTracers,nBenthicVertLevels+1) = &
              benthosSoluteDiffusivityI(iBenthicTracers,nBenthicVertLevels+1)

         do iBenthicLevels = 3,nBenthicVertLevels
            benthosSoluteDiffusivity(iBenthicTracers,iBenthicLevels) = &
                 (benthosSoluteDiffusivityI(iBenthicTracers,iBenthicLevels-1) + &
                 benthosSoluteDiffusivityI(iBenthicTracers,iBenthicLevels))/c2_benthos
         end do
      end do

      benthosDiffusivityI(1) = benthosPorosityI(1)*(Dm/tortuosity(1)+iDin_o)/benthosDepth**2
      benthosSolidDiffusivityI(1) = benthosSolidPorosityI(1)*iDin_o/benthosDepth**2

      do iBenthicLevels = 2,nBenthicVertLevels+1
         benthosDiffusivityI(iBenthicLevels) = benthosPorosityI(iBenthicLevels)*(Dm/ &
            tortuosity(iBenthicLevels) + iDin_o *exp(-benthosGridPointsI(iBenthicLevels)* &
            benthosDepth/x_biodiffusion))/benthosDepth**2
         benthosSolidDiffusivityI(iBenthicLevels) =benthosSolidPorosityI(iBenthicLevels)* &
            (iDin_o *exp(-benthosGridPointsI(iBenthicLevels)*benthosDepth/x_biodiffusion))/benthosDepth**2
      enddo

      benthosDiffusivity(1) = (Dm+iDin_o)/benthosDepth**2
      benthosDiffusivity(2) = benthosDiffusivityI(1)
      benthosDiffusivity(nBenthicVertLevels+2) = benthosDiffusivityI(nBenthicVertLevels+1)
      benthosDiffusivity(nBenthicVertLevels+1) = benthosDiffusivityI(nBenthicVertLevels+1)

      benthosSolidDiffusivity(1) = c0_benthos
      benthosSolidDiffusivity(2) = benthosSolidDiffusivityI(1)
      benthosSolidDiffusivity(nBenthicVertLevels+2) = benthosSolidDiffusivityI(nBenthicVertLevels+1)
      benthosSolidDiffusivity(nBenthicVertLevels+1) = benthosSolidDiffusivityI(nBenthicVertLevels+1)

      do iBenthicLevels = 3,nBenthicVertLevels
         benthosDiffusivity(iBenthicLevels) = (benthosDiffusivityI(iBenthicLevels-1) + &
            benthosDiffusivityI(iBenthicLevels))/c2_benthos
         benthosSolidDiffusivity(iBenthicLevels) = (benthosSolidDiffusivityI(iBenthicLevels-1) + &
            benthosSolidDiffusivityI(iBenthicLevels))/c2_benthos
      enddo

     endif   ! constant diffusivity
  endif

  !EOC
  ! NJ-Verified
!  write(*,*)' iDin_o :', iDin_o
!  write(*,*)'max_bio_diff:',max_bio_diff
!  write(*,*)'oceanBottomTemperature:',oceanBottomTemperature*m2percm2/sec_per_year*abs(oceanBottomTemperature/T_max_biodiffusion)
!
!            write(*,*)'benthosDiffusivity(:):',benthosDiffusivity(:)
!            write(*,*)'benthosDiffusivityI(:):',benthosDiffusivityI(:)
!            write(*,*)'benthosSolidDiffusivity(:):',benthosSolidDiffusivity(:)
!            write(*,*)'benthosSolidDiffusivityI(:):',benthosSolidDiffusivityI(:)
!            write(*,*)'benthosSoluteDiffusivityI(9,:):',benthosSoluteDiffusivityI(9,:)
!            write(*,*)'benthosSoluteDiffusivityI(1,:):',benthosSoluteDiffusivityI(1,:)

 end subroutine benthos_parameters_dt

!*****************************************************************************
!BOP
! !IROUTINE: benthos_molecular_diff
! !INTERFACE:

 subroutine benthos_molecular_diff(oceanBottomTemperature, oceanBottomSalinity, oceanBottomPressure, &
      benthosMolecularDiffusivity, benthos_indices)

! !DESCRIPTION:
!  Calculate the molecular diffusivity as a function of temperature and viscosity for each
! Solute.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(benthos_indices_type), intent(in) :: benthos_indices

  real(KIND=benthos_r8), intent(in) :: &
       oceanBottomTemperature, & ! (oC)
       oceanBottomSalinity, &
       oceanBottomPressure    ! N/m2

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(inout) :: &
       benthosMolecularDiffusivity ! (m2/s)

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  !NJ-TEST
  logical(KIND=benthos_log) :: test_flag
  !NJ-END

  real(KIND=benthos_r8) :: &
       viscosity, & !
       viscosity_fw, & !
       pressure_bar, &
       temperature_K, &
       temperature, &
       Av, Bv, Sa   !

  real(KIND=benthos_r8), dimension(nBenthicTracers) :: &
       Vb, mo, m1, D_o  ! parameters for solute molecular diffusion

  integer(KIND=benthos_i4) :: iBenthicTracers

  real(KIND=benthos_r8), parameter :: &
       u_star = 0.001_benthos_r8, &  ! m/s sqrt(bottom stress/fluid density)  from ocean model (fixed for now) !!!
       mu_fw = 0.01002_benthos_r8, &  ! viscosity parameter 1
       S1 = 35.165041_benthos_r8/35.0_benthos_r8, & ! salinity scaling
       A1 = 1.541_benthos_r8, &
       A2 = 0.01998_benthos_r8, &
       A3 = 9.52e-5_benthos_r8, &
       B1 = 7.974_benthos_r8, &
       B2 = 0.07561_benthos_r8, &
       B3 = 4.724e-5_benthos_r8, &
       Beta_o = 0.0889_benthos_r8, &
       D1 = 4.72e-9_benthos_r8/10000.0_benthos_r8 ! cm2/s to m2/s

  ! parameters for viscosity as function of T and pressure
  real(KIND=benthos_r8), parameter :: &
       E = 4.753_benthos_r8, & ! kJ/mol
       eta_o = 2.4055e-5_benthos_r8, & ! Pa s
       theta = 139.7_benthos_r8, & ! K
       a_1 = 4.42e-4_benthos_r8, & ! 1/bar
       b_1 = 9.565e-4_benthos_r8, & ! kJ/(mol bar)
       c_1 = 1.24e-2_benthos_r8, & ! K/bar
       R = 8.314462618_benthos_r8, & ! kg m2 /s2/K/mol
       Nm2_to_bar = 1e-5_benthos_r8, &  ! converts N/m2 to bar
       Kelvin = 273.15_benthos_r8
   !--------------
   ! viscosity in poise
   !---------------

   pressure_bar = oceanBottomPressure*Nm2_to_bar
   temperature_K = oceanBottomTemperature + Kelvin !oceanBottomTemperature + Kelvin
   temperature = max(c0_benthos, oceanBottomTemperature)

   Av = A1 + A2 * temperature - A3*temperature**2
   Bv = B1 - B2 * temperature + B3*temperature**2
   Sa = S1  * oceanBottomSalinity
!   viscosity = mu_fw * (c1_benthos + Av * Sa + Bv * Sa)
   viscosity = 1.7910_benthos_r8 - (6.144e-2_benthos_r8)*temperature + &
        (1.4510e-3_benthos_r8)*temperature**2 - (1.6826e-5_benthos_r8)*temperature**3 &
        - (1.5290e-4_benthos_r8)*pressure_bar + (8.3885e-8_benthos_r8)*pressure_bar**2 + &
        (2.4727e-3_benthos_r8)*oceanBottomSalinity + temperature*&
        (6.0574e-6_benthos_r8*pressure_bar - 2.6760e-9*pressure_bar**2) + &
        oceanBottomSalinity*(4.8429e-5_benthos_r8*temperature  - &
        4.7172e-6_benthos_r8*temperature**2 + 7.5986e-8_benthos_r8*temperature**3);
    viscosity = viscosity/100.0_benthos_r8    ! poise
!   viscosity_fw = eta_o * exp(a_1 * pressure_bar + (E - b_1 * pressure_bar) &
!        /(R * (temperature_K - theta - c_1*pressure_bar))) * 1.0e3_benthos_r8
!   viscosity = viscosity_fw *(c1_benthos + Av * Sa + Bv * Sa)

   !-----------------------------------------------------------------------
   !  Local copy of indices
   !-----------------------------------------------------------------------

   poca_ind    = benthos_indices%poca_ind
   pocb_ind    = benthos_indices%pocb_ind
   pocc_ind    = benthos_indices%pocc_ind
   pona_ind    = benthos_indices%pona_ind
   ponb_ind    = benthos_indices%ponb_ind
   ponc_ind    = benthos_indices%ponc_ind
   popa_ind    = benthos_indices%popa_ind
   popb_ind    = benthos_indices%popb_ind
   popc_ind    = benthos_indices%popc_ind
   o2_ind      = benthos_indices%o2_ind
   nh4_ind     = benthos_indices%nh4_ind
   h2po4_ind   = benthos_indices%h2po4_ind
   co2_ind     = benthos_indices%co2_ind
   no3_ind     = benthos_indices%no3_ind
   mno2a_ind   = benthos_indices%mno2a_ind
   mno2b_ind   = benthos_indices%mno2b_ind
   mn_ind      = benthos_indices%mn_ind
   feoh3a_ind  = benthos_indices%feoh3a_ind
   feoh3b_ind  = benthos_indices%feoh3b_ind
   fe_ind      = benthos_indices%fe_ind
   fepa_ind    = benthos_indices%fepa_ind
   fepb_ind    = benthos_indices%fepb_ind
   so4_ind     = benthos_indices%so4_ind
   h2s_ind     = benthos_indices%h2s_ind
   ch4_ind     = benthos_indices%ch4_ind
   dic_ind     = benthos_indices%dic_ind
   alk_ind     = benthos_indices%alk_ind
   hco3_ind    = benthos_indices%hco3_ind
   fes_ind     = benthos_indices%fes_ind
   fes2_ind    = benthos_indices%fes2_ind
   s_ind       = benthos_indices%s_ind
   caco3a_ind  = benthos_indices%caco3a_ind
   caco3b_ind  = benthos_indices%caco3b_ind
   co3_ind     = benthos_indices%co3_ind
   camgco3_ind = benthos_indices%camgco3_ind

   !------------
   ! parameters for solute molecular diffusion
   !-------------
   do iBenthicTracers = 1, nBenthicTracers
      mo(iBenthicTracers) = c0_benthos
      m1(iBenthicTracers) = c0_benthos
      Vb(iBenthicTracers) = c0_benthos
      D_o(iBenthicTracers) = c0_benthos
   end do

   Vb(o2_ind) = 27.9_benthos_r8
   Vb(co2_ind) = 37.3_benthos_r8
   Vb(h2s_ind) = 35.2_benthos_r8
   Vb(ch4_ind) = 37.7_benthos_r8

   mo(nh4_ind) = 9.5_benthos_r8
   m1(nh4_ind) = 0.413_benthos_r8
   mo(fe_ind) = 3.31_benthos_r8
   m1(fe_ind) = 0.150_benthos_r8
   mo(mn_ind) = 3.18_benthos_r8
   m1(mn_ind) = 0.155_benthos_r8
   mo(hco3_ind) = 5.06_benthos_r8
   m1(hco3_ind) = 0.275_benthos_r8
   mo(h2po4_ind) = 4.02_benthos_r8
   m1(h2po4_ind) = 0.223_benthos_r8
   mo(no3_ind) = 9.50_benthos_r8
   m1(no3_ind) = 0.388_benthos_r8
   mo(co3_ind) = 4.33_benthos_r8
   m1(co3_ind) = 0.199_benthos_r8
   mo(so4_ind) = 4.88_benthos_r8
   m1(so4_ind) = 0.232_benthos_r8
! Table 4.3 Boudreau book
   D_o(o2_ind) = D1   ! molecular diff at 25oC ~ 2.29e-5 cm2/s
   D_o(co2_ind) = D1  ! ~ 1.92e-5 cm2/s
   D_o(h2s_ind) = D1  ! ~ 2.10e-5 cm2/s
   D_o(ch4_ind) = D1  ! ~ 1.67e-5 cm2/s

   Do iBenthicTracers = 1,nBenthicTracers
      if (D_o(iBenthicTracers) .gt. c0_benthos) then ! Hayduk and Laudie (1974)
         benthosMolecularDiffusivity(iBenthicTracers) = D_o(iBenthicTracers) * (temperature_K)/ &
           (viscosity *Vb(iBenthicTracers)**0.6_benthos_r8)  ! Boudreau et al
      elseif (mo(iBenthicTracers) .gt. c0_benthos) then
         !!!T is in oC and not K!!!
         benthosMolecularDiffusivity(iBenthicTracers) = (mo(iBenthicTracers) + m1(iBenthicTracers)*temperature) * 1.0e-9_benthos_r8
      endif
   enddo
   benthosMolecularDiffusivity(alk_ind) = (benthosMolecularDiffusivity(hco3_ind) + benthosMolecularDiffusivity(co3_ind))/c2_benthos

   benthosMolecularDiffusivity(dic_ind) = (benthosMolecularDiffusivity(hco3_ind) + benthosMolecularDiffusivity(co3_ind) + benthosMolecularDiffusivity(co2_ind))/3.0_benthos_r8

!   write(*,*)'in benthos_molecular_diff'
!   write(*,*)'viscosity:', viscosity
!   write(*,*)'viscosity_fw:', viscosity_fw
!   write(*,*)'oceanBottomTemperature:',oceanBottomTemperature
!   write(*,*)'temperature:',temperature
!   write(*,*)'oceanBottomPressure:',oceanBottomPressure
!   write(*,*)'pressure_bar:',pressure_bar
!   write(*,*)'benthosMolecularDiffusivity(no3_ind):',benthosMolecularDiffusivity(no3_ind)
!   write(*,*)'benthosMolecularDiffusivity(co3_ind):',benthosMolecularDiffusivity(co3_ind)
!   write(*,*)'benthosMolecularDiffusivity(o2_ind):',benthosMolecularDiffusivity(o2_ind)
!   write(*,*)'benthosMolecularDiffusivity(co2_ind):',benthosMolecularDiffusivity(co2_ind)
!   write(*,*)'benthosMolecularDiffusivity(h2s_ind):',benthosMolecularDiffusivity(h2s_ind)
!   write(*,*)'benthosMolecularDiffusivity(ch4_ind):',benthosMolecularDiffusivity(ch4_ind)
!   write(*,*)'benthosMolecularDiffusivity(nh4_ind):',benthosMolecularDiffusivity(nh4_ind)
!   write(*,*)'benthosMolecularDiffusivity(fe_ind):',benthosMolecularDiffusivity(fe_ind)
!   write(*,*)'benthosMolecularDiffusivity(mn_ind):',benthosMolecularDiffusivity(mn_ind)
!   write(*,*)'benthosMolecularDiffusivity(hco3_ind):',benthosMolecularDiffusivity(hco3_ind)
!   write(*,*)'benthosMolecularDiffusivity(h2po4_ind):',benthosMolecularDiffusivity(h2po4_ind)
!   write(*,*)'benthosMolecularDiffusivity(alk_ind):',benthosMolecularDiffusivity(alk_ind)
!   write(*,*)'benthosMolecularDiffusivity(dic_ind):',benthosMolecularDiffusivity(dic_ind)

 end subroutine benthos_molecular_diff

!*****************************************************************************
!BOP
! !IROUTINE: carbonateStoichMatrix
! !INTERFACE:

 subroutine carbonateStoichMatrix


! !DESCRIPTION:
!  Specify mole fractions for sources and sinks in the following carbonate reactions
!
!  carbonate reactions
! 3 total carbonate equations
!
!
! p1. CaCO3 (calcite) --> Ca + CO3
! p2. CaCO3 (aragonite) --> Ca + CO3
! p3. Ca.85Mg.15CO3 --> 0.85 Ca + 0.15 Mg + CO3
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:
!
!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  integer(KIND=benthos_i4) :: iTracers, iReactions

   ! allocate parameter functions

   allocate(carbonateSourceStoich(nBenthicTracers,nCarbonateReactions))
   allocate(carbonateSinkStoich(nBenthicTracers,nCarbonateReactions))

   do iReactions = 1,nCarbonateReactions
      do iTracers = 1,nBenthicTracers
         carbonateSourceStoich(iTracers,iReactions) = c0_benthos
         carbonateSinkStoich(iTracers,iReactions) = c0_benthos
      end do
   end do

   ! SINK TERMS: Rates are in mmol/m3/s
   carbonateSinkStoich(caco3a_ind,1)  = c1_benthos/sediment_density
   carbonateSinkStoich(caco3b_ind,2)  = c1_benthos/sediment_density
   carbonateSinkStoich(camgco3_ind,3) = c1_benthos/sediment_density

   ! SOURCE TERMS
   carbonateSourceStoich(co3_ind,:) = c1_benthos
   carbonateSourceStoich(dic_ind,:) = c1_benthos
   carbonateSourceStoich(alk_ind,:) = c2_benthos

 end subroutine carbonateStoichMatrix

!*****************************************************************************
!BOP
! !IROUTINE: secondaryStoichMatrix
! !INTERFACE:

 subroutine secondaryStoichMatrix


! !DESCRIPTION:
!  Specify mole fractions for sources and sinks in the following secondary reactions

!  secondary reactions
!
! s1. O2 + NH4/2 + HCO3 --> NO3/2 + CO2 + 3/2H2O
! s2. O2 + 2Mn + 4HCO3 --> 2MnO2^a + 4CO2 + 2H2O
! s3. O2 + 4Fe + 8HCO3 + 2 H20 + 4 gamma H2PO4 + 24 gamma H^+ -->
!          4Fe(OH)3^a + 4 gamma [Fe-P]^a + 8CO2 + 16 gamma H2O
! s4. O2 + FeS/2 --> SO4/2 + Fe/2
! s5. O2+ 2/7FeS2 +  2/7H2O --> 4/7SO4 +  2/7Fe+ 4/7H
! s6. O2 + H2S/2 + HCO3 -->  SO4/2 + CO2 + H2O
! s7. O2 + CH4 --> CO2 + 4H^+
! s8. MnO2a + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gamma H -->
!             2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s9)s8b. MnO2b + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gammaH -->
!             2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s10)s9. Mn02a + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s11)s9b. Mn02b + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s12)s10. Fe(OH)3a+   gamma[Fe-P]a + H2S/2+ 2CO2 + 4 gamma H2O-->
!             Fe +   gamma H2PO4 + S/2 + 2HCO3 + (2+6gamma)H
! (s13)s10b. Fe(OH)3b+   gamma[Fe-P]b + H2S/2+ 2CO2 + 4 gamma H2O-->
!             Fe +  gamma H2PO4+ S/2 + 2HCO3 + (2+6 gamma)H
! (s14)s11. Fe + H2S -->  FeS + 2H
! (s15)s12. SO4 +  CH4 +  CO2 -->  2HCO3 + H2S
! (s16)s13. S +  H2O -->  3/4H2S +  SO4/4 +  H/2
! (s17)s14. FeS +  S --> FeS2
! (s18)s15. Fe(OH)3a + gamma [Fe-P]a --> Fe(OH)3b + gamma [Fe-P]b
! (s19)s16 MnO2a --> MnO2b
!
! 19 total secondary equations with 3 repetitions
! gamma = "ironBoundPFraction" is the fraction of  co-precipitated P with iron
!
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:
!
!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  integer(KIND=benthos_i4) :: iTracers, iReactions

   ! allocate parameter functions

   allocate(secondarySourceStoich(nBenthicTracers,nSecondaryReactions))
   allocate(secondarySinkStoich(nBenthicTracers,nSecondaryReactions))

   do iReactions = 1,nSecondaryReactions
      do iTracers = 1,nBenthicTracers
         secondarySourceStoich(iTracers,iReactions) = c0_benthos
         secondarySinkStoich(iTracers,iReactions) = c0_benthos
      end do
   end do

   ! add units conversion for secondary rates which are in mmol/m3/s except for the last
   ! three reactions (s17-s19) which only involve solid phase material. These are in mmol/kg/s

   ! SINK TERMS
   ! s1
   secondarySinkStoich(o2_ind,1) = c1_benthos
   !secondarySinkStoich(hco3_ind,1) = c1_benthos
   secondarySinkStoich(dic_ind,1) = c1_benthos
   secondarySinkStoich(nh4_ind,1) = p5_benthos
   secondarySinkStoich(alk_ind,1) = c1_benthos

   !s2

   secondarySinkStoich(mn_ind,2) = c2_benthos
   !secondarySinkStoich(hco3_ind,2) = 4.0_benthos_r8
   secondarySinkStoich(dic_ind,2) = 4.0_benthos_r8
   secondarySinkStoich(o2_ind,2) = c1_benthos
   secondarySinkStoich(alk_ind,2) = 4.0_benthos_r8

   ! s3
   secondarySinkStoich(o2_ind,3) = c1_benthos
   !secondarySinkStoich(hco3_ind,3) = 8.0_benthos_r8
   secondarySinkStoich(fe_ind,3) = 4.0_benthos_r8
   secondarySinkStoich(h2po4_ind,3) = 4.0_benthos_r8*ironBoundPFraction
   secondarySinkStoich(alk_ind,3) = 8.0_benthos_r8

   !s4
   secondarySinkStoich(o2_ind,4) = c1_benthos
   secondarySinkStoich(fes_ind,4) = p5_benthos/sediment_density

   !s5
   secondarySinkStoich(o2_ind,5) = c1_benthos
   secondarySinkStoich(fes2_ind,5) = c2_benthos/7.0_benthos_r8/sediment_density
   secondarySinkStoich(alk_ind,5) = 4.0_benthos_r8/7.0_benthos_r8

   !s6
   secondarySinkStoich(o2_ind,6) = c1_benthos
   !secondarySinkStoich(hco3_ind,6) = c1_benthos
   secondarySinkStoich(dic_ind,6) = c1_benthos
   secondarySinkStoich(h2s_ind,6) = p5_benthos
   secondarySinkStoich(alk_ind,6) = c1_benthos

   !s7
   secondarySinkStoich(o2_ind,7) = c1_benthos
   secondarySinkStoich(ch4_ind,7) = c1_benthos

   !s8 and s9
   secondarySinkStoich(mno2a_ind,8) = c1_benthos/sediment_density
   secondarySinkStoich(fe_ind,8) = c2_benthos
   secondarySinkStoich(fe_ind,9) = c2_benthos
   secondarySinkStoich(h2po4_ind,8) = c2_benthos*ironBoundPFraction
   secondarySinkStoich(h2po4_ind,9) = c2_benthos*ironBoundPFraction
   !secondarySinkStoich(hco3_ind,8:9) = c2_benthos
   secondarySinkStoich(dic_ind,8) = c2_benthos
   secondarySinkStoich(dic_ind,9) = c2_benthos
   secondarySinkStoich(mno2b_ind,9) = c1_benthos/sediment_density
   secondarySinkStoich(alk_ind,8) = c2_benthos
   secondarySinkStoich(alk_ind,9) = c2_benthos

   !s10 and s11
   secondarySinkStoich(h2s_ind,10) = c1_benthos
   secondarySinkStoich(h2s_ind,11) = c1_benthos
   secondarySinkStoich(mno2a_ind,10) = c1_benthos/sediment_density
   secondarySinkStoich(mno2b_ind,11) = c1_benthos/sediment_density
   !secondarySinkStoich(co2_ind,10:11) = c2_benthos
   secondarySinkStoich(dic_ind,10) = c2_benthos
   secondarySinkStoich(dic_ind,11) = c2_benthos

   !s12 and s13
   secondarySinkStoich(h2s_ind,12) = p5_benthos
   secondarySinkStoich(h2s_ind,13) = p5_benthos
   !econdarySinkStoich(co2_ind,12:13) = c2_benthos
   secondarySinkStoich(dic_ind,12) = c2_benthos
   secondarySinkStoich(dic_ind,13) = c2_benthos
   secondarySinkStoich(feoh3a_ind,12) = c1_benthos/sediment_density
   secondarySinkStoich(feoh3b_ind,13) = c1_benthos/sediment_density
   secondarySinkStoich(fepa_ind,12) = c1_benthos*ironBoundPFraction/sediment_density
   secondarySinkStoich(fepb_ind,13) = c1_benthos*ironBoundPFraction/sediment_density
   secondarySinkStoich(alk_ind,12) = 6.0_benthos_r8*ironBoundPFraction
   secondarySinkStoich(alk_ind,13) = 6.0_benthos_r8*ironBoundPFraction

   !s14
   secondarySinkStoich(h2s_ind,14) = c1_benthos
   secondarySinkStoich(fe_ind,14) = c1_benthos
   secondarySinkStoich(alk_ind,14) = c2_benthos

   !s15
   secondarySinkStoich(ch4_ind,15) = c1_benthos
   !secondarySinkStoich(co2_ind,15) = c1_benthos
   secondarySinkStoich(so4_ind,15) = c1_benthos
   secondarySinkStoich(dic_ind,15) = c1_benthos

   !s16
   secondarySinkStoich(s_ind,16) = c1_benthos/sediment_density

   !s17
   secondarySinkStoich(fes_ind,17) = c1_benthos
   secondarySinkStoich(s_ind,17) = c1_benthos

   !s18
   secondarySinkStoich(feoh3a_ind,18) = c1_benthos
   secondarySinkStoich(fepa_ind,18) = c1_benthos*ironBoundPFraction

   !s19
   secondarySinkStoich(mno2a_ind,19) = c1_benthos

   ! SOURCE TERMS
   !s1
   secondarySourceStoich(no3_ind,1) = p5_benthos
   !secondarySourceStoich(co2_ind,1) = c1_benthos
   secondarySourceStoich(dic_ind,1) = c1_benthos

   !s2
   secondarySourceStoich(mno2a_ind,2) = c2_benthos/sediment_density
   !secondarySourceStoich(co2_ind,2) = 4.0_benthos_r8
   secondarySourceStoich(dic_ind,2) = 4.0_benthos_r8

   !s3
   secondarySourceStoich(feoh3a_ind,3) = 4.0_benthos_r8/sediment_density
   secondarySourceStoich(fepa_ind,3) = 4.0_benthos_r8/sediment_density*ironBoundPFraction
   !secondarySourceStoich(co2_ind,3) = 8.0_benthos_r8
   secondarySourceStoich(dic_ind,3) = 8.0_benthos_r8
   secondarySourceStoich(alk_ind,3) = 24.0_benthos_r8*ironBoundPFraction

   !s4
   secondarySourceStoich(so4_ind,4) = p5_benthos
   secondarySourceStoich(fe_ind,4) = p5_benthos

   !s5
   secondarySourceStoich(so4_ind,5) = 4.0_benthos_r8/7.0_benthos_r8
   secondarySourceStoich(fe_ind,5) = c2_benthos/7.0_benthos_r8

   !s6
   secondarySourceStoich(so4_ind,6) = p5_benthos
   !secondarySourceStoich(co2_ind,6) = c1_benthos
   secondarySourceStoich(dic_ind,6) = c1_benthos

   !s7
   !secondarySourceStoich(co2_ind,7) = c1_benthos
   secondarySourceStoich(dic_ind,7) = c1_benthos
   secondarySourceStoich(alk_ind,7) = 4.0_benthos_r8

   !s8
   secondarySourceStoich(feoh3a_ind,8) = c2_benthos/sediment_density
   secondarySourceStoich(fepa_ind,8) = c2_benthos/sediment_density*ironBoundPFraction
   secondarySourceStoich(mn_ind,8) = c1_benthos
   !secondarySourceStoich(co2_ind,8) = c2_benthos
   secondarySourceStoich(dic_ind,8) = c2_benthos
   secondarySourceStoich(alk_ind,8) = 12.0_benthos_r8*ironBoundPFraction

   !s9
   secondarySourceStoich(feoh3a_ind,9) = c2_benthos/sediment_density
   secondarySourceStoich(fepa_ind,9) = c2_benthos/sediment_density*ironBoundPFraction
   secondarySourceStoich(mn_ind,9) = c1_benthos
   !secondarySourceStoich(co2_ind,9) = c2_benthos
   secondarySourceStoich(dic_ind,9) = c2_benthos
   secondarySourceStoich(alk_ind,9) = 12.0_benthos_r8*ironBoundPFraction

   !s10
   secondarySourceStoich(mn_ind,10) = c1_benthos
   secondarySourceStoich(s_ind,10) = c1_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,10) = c2_benthos
   secondarySourceStoich(dic_ind,10) = c2_benthos
   secondarySourceStoich(alk_ind,10) = c2_benthos

   !s11
   secondarySourceStoich(mn_ind,11) = c1_benthos
   secondarySourceStoich(s_ind,11) = c1_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,11) = c2_benthos
   secondarySourceStoich(dic_ind,11) = c2_benthos
   secondarySourceStoich(alk_ind,11) = c2_benthos

   !s12
   secondarySourceStoich(fe_ind,12) = c1_benthos
   secondarySourceStoich(h2po4_ind,12) = ironBoundPFraction
   secondarySourceStoich(s_ind,12) = p5_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,12) = c2_benthos
   secondarySourceStoich(dic_ind,12) = c2_benthos

   !s13
   secondarySourceStoich(fe_ind,13) = c1_benthos
   secondarySourceStoich(h2po4_ind,13) = ironBoundPFraction
   secondarySourceStoich(s_ind,13) = p5_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,13) = c2_benthos
   secondarySourceStoich(dic_ind,13) = c2_benthos

   !s14
   secondarySourceStoich(fes_ind,14) = c1_benthos/sediment_density

   !s15
   !secondarySourceStoich(hco3_ind,15) = c2_benthos
   secondarySourceStoich(h2s_ind,15) = c1_benthos
   secondarySourceStoich(dic_ind,15) = c2_benthos

   !s16
   secondarySourceStoich(h2s_ind,16) = 3.0_benthos_r8/4.0_benthos_r8
   secondarySourceStoich(so4_ind,16) = 0.25_benthos_r8

   !s17
   secondarySourceStoich(fes2_ind,17) = c1_benthos

   !s18
   secondarySourceStoich(feoh3b_ind,18) = c1_benthos
   secondarySourceStoich(fepb_ind,18) = c1_benthos*ironBoundPFraction

   !s19
   secondarySourceStoich(mno2b_ind,19) = c1_benthos

 end subroutine secondaryStoichMatrix
 
!*****************************************************************************
!BOP
! !IROUTINE: computeNetElements
! !INTERFACE:

 subroutine computeNetElements (netElements, nBenthicVertLevels, benthosTracerBulk, &
      vertBenthosGridThickI)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels
  real(KIND=benthos_r8), dimension(:,:), intent(in) :: benthosTracerBulk
  real(KIND=benthos_r8), dimension(:), intent(in) :: vertBenthosGridThickI

! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nElements), intent(out) :: netElements

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
!   real(KIND=benthos_r8), dimension(:,:), allocatable :: &
!     elementRatios

   real(KIND=benthos_r8), dimension(nBenthicTracers) :: &
     columnBenthosTracerBulk

   integer (KIND=benthos_i4) :: &
     iTracers, &
     iElements

   netElements(:) = c0_benthos

   do iTracers = 1,nBenthicTracers
       columnBenthosTracerBulk(iTracers) = c0_benthos
      call sum_benthos_column(columnBenthosTracerBulk(iTracers), benthosTracerBulk(iTracers,:),nBenthicVertLevels,vertBenthosGridThickI)
   end do
   do iElements = 1,nElements
      do iTracers = 1,nBenthicTracers
           netElements(iElements) = netElements(iElements) + elementRatios(iTracers,iElements)*columnBenthosTracerBulk(iTracers)
      end do
   end do

 end subroutine computeNetElements

!*****************************************************************************
!BOP
! !IROUTINE: sum_benthos_column
! !INTERFACE:

 subroutine sum_benthos_column(columnBenthosTracerBulk, benthosTracerBulk, &
      nBenthicVertLevels,vertBenthosGridThickI)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels

  real(KIND=benthos_r8), dimension(:), intent(in):: &
       benthosTracerBulk, &
       vertBenthosGridThickI

! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: columnBenthosTracerBulk

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
   integer(KIND=benthos_i4) :: &
     iLevels

   columnBenthosTracerBulk = c0_benthos

   do iLevels = 1,nBenthicVertLevels+1
      columnBenthosTracerBulk = columnBenthosTracerBulk + benthosTracerBulk(iLevels)*vertBenthosGridThickI(iLevels)*benthosDepth
   end do

 end subroutine sum_benthos_column

!*****************************************************************************
!BOP
! !IROUTINE: updateBenthosBoundaryConditions
! !INTERFACE:

 subroutine updateBenthosBoundaryConditions  (oceanSedFlux, oceanPrecipFlux, totalOceanSolidFlux, &
      oceanTracerConc, deepBenthosBottomFlux, oceanBottomDensity, nBenthicVertLevels,FluxSed, &
      dhtop, vel, dt, oceanBottomPhosphate, oceanBottomAmmonium, oceanBottomNitrate, oceanBottomOxygen, &
      oceanBottomIron, oceanBottomDIC, oceanBottomAlkalinity, oceanPOCFlux, oceanPONFlux, &
      oceanPOPFlux, oceanParticulateIronFlux, oceanCalciteFlux,benthosSolidPorosityTop, &
      oceanBiogenicSilicaFlux, oceanDustFlux, oceanBottomSilicate, oceanBottomManganese, &
      benthos_output, column, config_use_benthos_monthly_forcing, oceanBottomDepth, &
      terraPOCFlux, lithogenicFlux)

! !DESCRIPTION:
!  compute benthos boundary fluxes and concentrations at the start of each timestep
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nBenthicVertLevels, &  ! number of vertical interfaces in the benthos
       column

  real(KIND=benthos_r8), intent(in) :: &
       dt, &
       oceanBottomDensity, &
       oceanBottomPhosphate, &
       oceanBottomAmmonium, &
       oceanBottomNitrate, &
       oceanBottomOxygen, &
       oceanBottomIron, &
       oceanBottomManganese, &
       oceanBottomDIC, &
       oceanBottomAlkalinity, &
       oceanPOCFlux, &
       oceanPOPFlux, &
       oceanPONFlux, &
       oceanParticulateIronFlux, &
       oceanCalciteFlux, &
       benthosSolidPorosityTop, &
       oceanBottomDepth, &
       oceanBiogenicSilicaFlux, &
       oceanDustFlux, &
       terraPOCFlux, &
       lithogenicFlux

  logical(KIND=benthos_log), intent(in) :: &
       config_use_benthos_monthly_forcing

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(inout) :: &
     oceanSedFlux,    & ! particulate flux of sinking ocean tracers (mmol/kg*m/s)
     oceanTracerConc, & ! bottom ocean concentration of solutes (mmol/m3)
     oceanPrecipFlux, & ! flux of precipitated material from ocean bottom (mmol/kg*m/s)
     totalOceanSolidFlux, & ! total solid flux divided by benthosDepth (mmol/kg/s)
     deepBenthosBottomFlux  ! prescribed deep storage flux (mmol/m2/s)

  real(KIND=benthos_r8),  intent(inout) :: &
     oceanBottomSilicate

  type(benthos_output_type), intent(inout) :: benthos_output

! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: &
     FluxSed, & ! kg/m2/s
     dhtop,   & ! (m) increase in benthos surface over the timestep
     vel        ! (1/s) normalized velocity of the surface and bottom boundaries

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  integer(KIND=benthos_i4) :: iTracers

  real (KIND=benthos_r8) :: &
       pocFluxRate, &  ! mmol (m/kg)/s
       ponFluxRate, &
       popFluxRate, &
       ocean_feoh3aFlux, &   ! mmol (m/kg)/s
       ocean_caco3Flux, &
       fluxSed_middelburg, &  ! Middelburg et al (1997) parameterization
       fluxSed_forcing  ! based on ocean BGC forcing and monthly climatology

!-----------------------------------------------------------------------
!EOC

  ! pocFluxRate = pocFluxRate_o/sediment_density

   if (useOceanBottomState .or. config_use_benthos_monthly_forcing) then
      pocFluxRate = (oceanPOCFlux+terraPOCFlux)/sediment_density/benthosSolidPorosityTop
      ponFluxRate = pocFluxRate/CtoP*NtoP
      popFluxRate = pocFluxRate/CtoP
      ocean_feoh3aFlux = oceanParticulateIronFlux/sediment_density/benthosSolidPorosityTop
      ocean_caco3Flux = oceanCalciteFlux/sediment_density/benthosSolidPorosityTop
      fluxSed_forcing = (oceanPOCFlux + terraPOCFlux) * kg_per_mg * gC_per_molC * mass_per_massPOC &
         + oceanBiogenicSilicaFlux * kg_per_mg * gSi_per_molSi * massOpal_per_massbSi &
         + oceanDustFlux + lithogenicFlux &
         + oceanCalciteFlux * kg_per_mg * mass_per_molPIC
   else
      pocFluxRate = 6000.0_benthos_r8/sediment_density/sec_per_year
      ponFluxRate = pocFluxRate/CtoP*NtoP
      popFluxRate = pocFluxRate/CtoP
      ocean_feoh3aFlux = 57.7_benthos_r8/sec_per_year/sediment_density
      ocean_caco3Flux = 1000.0_benthos_r8 / sec_per_year/sediment_density
      oceanBottomSilicate = 11.0_benthos_r8
      fluxSed_forcing = c0_benthos
   end if

   ! sediment sinking flux will eventually come from the ocean model
   !  sea bed accumulation rates (ob, Yenisey  0.2-0.7 cm/y

   if (useSedimentation) then
      if (useOceanBottomState .or. config_use_benthos_monthly_forcing) then
         fluxSed = fluxSed_forcing
      else
         fluxSed_middelburg = sedFlux3 * 10.0_benthos_r8**(-sedFlux1-sedFlux2*oceanBottomDepth)*sediment_density * benthosSolidPorosityTop
         if (useFastSedimentation) then
            fluxSed = 1000.0_benthos_r8*fluxSed_o
         else
            fluxSed = fluxSed_middelburg !fluxSed_o
         end if
      end if
      if (.not. usePositiveSedimentation) fluxSed = -fluxSed
   else
      fluxSed = c0_benthos
   end if

   dhtop = fluxSed/sediment_density/benthosSolidPorosityTop * dt
   vel = dhtop/benthosDepth/dt

   do iTracers = 1,nBenthicTracers
      oceanSedFlux(iTracers) = c0_benthos
      deepBenthosBottomFlux(iTracers) = c0_benthos
      oceanTracerConc(iTracers) = c0_benthos
      oceanPrecipFlux(iTracers) = c0_benthos
   end do

   if (useDeepSource) then
   ! only a couple have a source term  (H2S has an upward flux -50.0  mmol/m2/y
      deepBenthosBottomFlux(h2s_ind) = deepFlux_h2s
   end if

   if (useBGCSinkingFlux .or. config_use_benthos_monthly_forcing) then
      !!!! For testing only
      if (benthosTestCase .eq. 4) then
         oceanSedFlux(poca_ind) = dhtop / dt * 0.0158_benthos_r8
      else
         oceanSedFlux(poca_ind) = pocFluxRate*fracHighlyReactive   ! mmol (m/kg)/s  (49% highy reactive, Westrich and Berner 1984)
      end if
      oceanSedFlux(pocb_ind) = pocFluxRate*fracLessReactive ! 15% moderately reactive
      oceanSedFlux(pocc_ind) = pocFluxRate*fracRefractory   ! 34% nonreactive
      oceanSedFlux(popa_ind) = oceanSedFlux(poca_ind)/CtoP  !pocFluxRate/popFluxRate*fracHighlyReactive   ! mmol (m/kg)/s  (49% highy reactive, Westrich and Berner 1984)
      oceanSedFlux(popb_ind) = oceanSedFlux(pocb_ind)/CtoP_r !popFluxRate*fracLessReactive ! 15% moderately reactive
      oceanSedFlux(popc_ind) = oceanSedFlux(pocc_ind)/CtoP_r !popFluxRate*fracRefractory ! 34% nonreactive
      oceanSedFlux(pona_ind) = oceanSedFlux(popa_ind)*NtoP   !ponFluxRate*fracHighlyReactive   ! mmol (m/kg)/s  (49% highy reactive, Westrich and Berner 1984)
      oceanSedFlux(ponb_ind) = oceanSedFlux(popb_ind)*NtoP_r !ponFluxRate*fracLessReactive ! 15% moderately reactive
      oceanSedFlux(ponc_ind) = oceanSedFlux(popc_ind)*NtoP_r !ponFluxRate*fracRefractory ! 34% nonreactive
      oceanSedFlux(mno2a_ind) = 0.03_benthos_r8/sec_per_year/sediment_density !mmol/m2/y to mmol (m/kg)/s
      oceanSedFlux(mno2b_ind) = c0_benthos
      oceanSedFlux(feoh3a_ind) = ocean_feoh3aFlux
      oceanSedFlux(feoh3b_ind) = c0_benthos
      oceanSedFlux(fepa_ind) = oceanSedFlux(feoh3a_ind)*ironBoundPFraction
      oceanSedFlux(fepa_ind) = c0_benthos
      oceanSedFlux(caco3a_ind) = ocean_caco3Flux * fracCalcite
      oceanSedFlux(caco3b_ind) = ocean_caco3Flux * fracAragonite
      oceanSedFlux(camgco3_ind) = ocean_caco3Flux * fracMgCalcite
   end if

   ! initial Ocean tracer concentrations  : zero for solids
   ! non-zero for solutes  (some from ocean model)
   !
   ! concentrations in mmol/m3 (umol/L)

   if (useOceanConc) then
      if (useOceanBottomState .or. config_use_benthos_monthly_forcing) then
         oceanTracerConc(o2_ind) = oceanBottomOxygen !350.0_benthos_r8  !umol/L (mmol/m3) Reed et al, 2011
         oceanTracerConc(nh4_ind) = oceanBottomAmmonium !c2_benthos  !umol/L  (Reed)
         oceanTracerConc(h2po4_ind) = oceanBottomPhosphate !c1_benthos  !umol/L Reed et al. 2011
         oceanTracerConc(no3_ind) = oceanBottomNitrate !6.0_benthos_r8 ! umol/L Reed.
         oceanTracerConc(mn_ind) = oceanBottomManganese !c0_benthos ! umol/L
         oceanTracerConc(fe_ind) = oceanBottomIron !c1_benthos ! umol/L Reed
         oceanTracerConc(so4_ind) = 13.0_benthos_r8*mM_umolperL !12.0_benthos_r8*mM_umolperL !28000.0_benthos_r8 !higher? 28000?   mmol/L Reed  to mmol/m3
         oceanTracerConc(h2s_ind) = c0_benthos ! Reed
         oceanTracerConc(ch4_ind) = c0_benthos !
         ! these shouldn't matter
         oceanTracerConc(co2_ind) = c0_benthos !0.125_benthos_r8*mM_umolperL  !mM Krumins et al, 2013
         oceanTracerConc(hco3_ind) = c0_benthos !c2_benthos*mM_umolperL  !mM  Krumins et al  (essentially *DIC*) and TA?
         oceanTracerConc(co3_ind) = c0_benthos !70.0e-3_benthos_r8*oceanBottomDensity !bottom_water_density  ! mmol/m3
         !umol/kg*density of water 1.027 g/m3  (Sulpis et al 2018)
         oceanTracerConc(dic_ind) = oceanBottomDIC !oceanTracerConc(co3_ind) + oceanTracerConc(co2_ind) + oceanTracerConc(hco3_ind)
         oceanTracerConc(alk_ind) = oceanBottomAlkalinity !2306.0_benthos_r8  ! mmol/m3  Krumins et al.
      else
         oceanTracerConc(o2_ind) = 350.0_benthos_r8  !umol/L (mmol/m3) Reed et al, 2011
         oceanTracerConc(nh4_ind) = c2_benthos  !umol/L  (Reed)
         oceanTracerConc(h2po4_ind) = c1_benthos  !umol/L Reed et al. 2011
         oceanTracerConc(co2_ind) = 0.125_benthos_r8*mM_umolperL  !mM Krumins et al, 2013
         oceanTracerConc(no3_ind) = 6.0_benthos_r8 ! umol/L Reed.
         oceanTracerConc(mn_ind) = c0_benthos ! umol/L
         oceanTracerConc(fe_ind) = c0_benthos ! umol/L !c1_benthos ! umol/L Reed
         oceanTracerConc(so4_ind) = 12.0_benthos_r8*mM_umolperL ! mmol/L Reed  to mmol/m3 or  28000_benthos_r8 !
         oceanTracerConc(h2s_ind) = c0_benthos ! Reed
         oceanTracerConc(ch4_ind) = c0_benthos !
         oceanTracerConc(hco3_ind) = c2_benthos*mM_umolperL  !mM  Krumins et al  (essentially *DIC*)
         oceanTracerConc(co3_ind) = 70.0e-3_benthos_r8 !*oceanBottomDensity !bottom_water_density  ! mmol/m3
         !umol/kg*density of water 1.027 g/m3  (Sulpis et al 2018)
         oceanTracerConc(dic_ind) = 2.022_benthos_r8*mM_umolperL !oceanTracerConc(co3_ind) + oceanTracerConc(co2_ind) + oceanTracerConc(hco3_ind)
         oceanTracerConc(alk_ind) = 2306.0_benthos_r8  ! mmol/m3  Krumins et al.
      end if
   end if

   !  make sure dic is consistent
   !initialBenthicTracerBulk(dic_ind,:) = initialBenthicTracerBulk(hco3_ind,:) + initialBenthicTracerBulk(co2_ind,:)  +  initialBenthicTracerBulk(co3_ind,:)
   ! define ocean precipitation flux: ignore this right now
   if (useBGCSinkingFlux) then
      ! use secondary equation 2 in ocean bottom layer
      oceanPrecipFlux(mno2a_ind) = oceanTracerConc(mn_ind) * oceanTracerConc(o2_ind) &
         * k_s2 * c2_benthos/sediment_density/benthosSolidPorosityTop 
      !Old spin-up fixes the ocean bottom manganese concentration at 0.0096 mmol mn/m3, i.e.
      !0.56_benthos_r8*oceanTracerConc(o2_ind) * k_s2 *  c2_benthos/sediment_density*benthosDepth

      ! use secondary equation 3 in ocean bottom layer
      oceanPrecipFlux(feoh3a_ind) = oceanTracerConc(fe_ind) * oceanTracerConc(o2_ind) &
         * k_s3 * 4.0_benthos_r8/sediment_density/benthosSolidPorosityTop
      !Old spin-up fixes the ocean bottom iron concentration at 0.0017 mmol fe/m3 everywhere, i.e.
      !0.10_benthos_r8* oceanTracerConc(o2_ind)*k_s3 * 4.0_benthos_r8/sediment_density*benthosDepth
      oceanPrecipFlux(fepa_ind) = oceanPrecipFlux(feoh3a_ind)*ironBoundPFraction
   end if

   do iTracers = 1,nBenthicTracers
      totalOceanSolidFlux(iTracers) = (oceanPrecipFlux(iTracers) + oceanSedFlux(iTracers))/benthosDepth
   end do

   !  Assign bottom ocean fields and fluxes to benthos_output

   benthos_output % oceanTotalSedimentMassFlux(column) = FluxSed
   benthos_output % oceanBottomOxygen(column)   = oceanTracerConc(o2_ind)
   benthos_output % oceanBottomAmmonium(column) = oceanTracerConc(nh4_ind)
   benthos_output % oceanBottomPhosphate(column) = oceanTracerConc(h2po4_ind)
   benthos_output % oceanBottomCarbonDioxide(column) = oceanTracerConc(co2_ind)
   benthos_output % oceanBottomNitrate(column) = oceanTracerConc(no3_ind)
   benthos_output % oceanBottomManganese(column) = oceanTracerConc(mn_ind)
   benthos_output % oceanBottomIron(column) = oceanTracerConc(fe_ind)
   benthos_output % oceanBottomHydrogenSulfide(column) = oceanTracerConc(h2s_ind)
   benthos_output % oceanBottomMethane(column) = oceanTracerConc(ch4_ind)
   benthos_output % oceanBottomBicarbonate(column) = oceanTracerConc(hco3_ind)
   benthos_output % oceanBottomCarbonate(column) = oceanTracerConc(co3_ind)
   benthos_output % oceanBottomDIC(column) = oceanTracerConc(dic_ind)
   benthos_output % oceanBottomAlkalinity(column) = oceanTracerConc(alk_ind)
   benthos_output % oceanBottomSulfate(column) = oceanTracerConc(so4_ind)
   benthos_output % oceanBottomSilicate(column) = oceanBottomSilicate

   benthos_output % oceanPOCaFlux(column) =  oceanSedFlux(poca_ind)
   benthos_output % oceanPOCbFlux(column) =  oceanSedFlux(pocb_ind)
   benthos_output % oceanPOCcFlux(column) =  oceanSedFlux(pocc_ind)
   benthos_output % oceanPOPaFlux(column) =  oceanSedFlux(popa_ind)
   benthos_output % oceanPOPbFlux(column) =  oceanSedFlux(popb_ind)
   benthos_output % oceanPOPcFlux(column) =  oceanSedFlux(popc_ind)
   benthos_output % oceanPONaFlux(column) =  oceanSedFlux(pona_ind)
   benthos_output % oceanPONbFlux(column) =  oceanSedFlux(ponb_ind)
   benthos_output % oceanPONcFlux(column) =  oceanSedFlux(ponc_ind)
   benthos_output % oceanManganeseOxideaFlux(column) = oceanSedFlux(mno2a_ind)
   benthos_output % oceanManganeseOxidebFlux(column) = oceanSedFlux(mno2b_ind)
   benthos_output % oceanIronOxyhydroxideaFlux(column) = oceanSedFlux(feoh3a_ind)
   benthos_output % oceanIronOxyhydroxidebFlux(column) = oceanSedFlux(feoh3b_ind)
   benthos_output % oceanIronBoundPhosphorusaFlux(column) = oceanSedFlux(fepa_ind)
   benthos_output % oceanIronBoundPhosphorusbFlux(column) = oceanSedFlux(fepb_ind)
   benthos_output % oceanCalciteFlux(column) = oceanSedFlux(caco3a_ind)
   benthos_output % oceanAragoniteFlux(column) = oceanSedFlux(caco3b_ind)
   benthos_output % oceanMgCalciteFlux(column) = oceanSedFlux(camgco3_ind)

   ! define ocean precipitation flux: ignore this right now
   benthos_output % oceanPrecipManganeseOxidea = oceanPrecipFlux(mno2a_ind)
   benthos_output % oceanPrecipManganeseOxideb = oceanPrecipFlux(mno2b_ind)
   benthos_output % oceanPrecipIronOxyhydroxidea = oceanPrecipFlux(feoh3a_ind)
   benthos_output % oceanPrecipIronOxyhydroxideb = oceanPrecipFlux(feoh3b_ind)
   benthos_output % oceanPrecipIronBoundPhosphorusa = oceanPrecipFlux(fepa_ind)
   benthos_output % oceanPrecipIronBoundPhosphorusb = oceanPrecipFlux(fepb_ind)

 end subroutine updateBenthosBoundaryConditions

!*****************************************************************************
!BOP
! !IROUTINE: benthos_flags_init
! !INTERFACE:

 subroutine benthos_flags_init (setBenthosFlags, useCarbonateSaturation_in, useSecondaryReactions_in, &
            useDepthDependentPorosity_in, useNonZeroDiffusivity_in, useSedimentation_in, &
            useFastSedimentation_in, usePositiveSedimentation_in, useBgcSinkingFlux_in, &
            useStepInitialProfiles_in, useBenthicReactions_in, useFluxCorrection_in, &
            useOceanConc_in, useDeepSource_in, useConstantDiffusivity_in, useOceanBottomState_in, &
            useBenthosMonthlyForcing_in)

! !DESCRIPTION:
!  Initialize flags for testcases
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer (KIND=benthos_i4), intent(in) :: setBenthosFlags

  logical (KIND=benthos_log), intent(inout) :: &
       useCarbonateSaturation_in,    & ! true turns on carbonate chemistry
       useSecondaryReactions_in,     & ! true computes secondary reaction set
       useDepthDependentPorosity_in, & ! true uses a depth dependent porosity
       useNonZeroDiffusivity_in,     & ! false sets the diffusivity to zero
       useSedimentation_in,          & ! true allows sediment accumulation or loss at the benthos surface
       useFastSedimentation_in,      & ! true enhances sedimentation rates to test advection
       usePositiveSedimentation_in,  & ! false switches from a sediment flux to resuspension at the benthic surface
       useBgcSinkingFlux_in,         & ! true uses ocean sinking flux of particulate benthos tracers
       useStepInitialProfiles_in,    & ! true sets the initial benthos profiles as step functions
       useBenthicReactions_in,       & ! turns on and off reaction bio-chemistry
       useFluxCorrection_in,         & ! true tracks transport errors in the storage flux
       useOceanConc_in,              & ! false sets ocean concentrations of solutes to zero
       useDeepSource_in,             & ! false sets to zero all deep storage fluxes
       useConstantDiffusivity_in,    & ! true removes the depth dependence of the diffusivity
       useOceanBottomState_in,       & ! Use the ocean biogeochemistry to force the benthos
       useBenthosMonthlyForcing_in     ! Use data file to force the benthos

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  !NJ-TEST
  logical(KIND=benthos_log) :: test_flag
  integer(KIND=benthos_i4) :: iBenthicTracers, iSecondaryReactions, iElements, iBenthicLevels
  !NJ-END
!-----------------------------------------------------------------------
!EOC
!NJ-TEST
!  write(*,*) 'test benthos_flags_init'
!NJ-END

  benthosTestCase = setBenthosFlags
  useOceanBottomState = useOceanBottomState_in
  useBenthosMonthlyForcing = useBenthosMonthlyForcing_in

 SELECT CASE (setBenthosFlags)
    CASE(:-1)

       ! 'User Specified Test Case: Use Flags in Namelist'

       useCarbonateSaturation = useCarbonateSaturation_in
       useSecondaryReactions = useSecondaryReactions_in
       useDepthDependentPorosity = useDepthDependentPorosity_in
       useNonZeroDiffusivity = useNonZeroDiffusivity_in
       useSedimentation = useSedimentation_in
       usePositiveSedimentation = usePositiveSedimentation_in
       useFastSedimentation = useFastSedimentation_in
       useBgcSinkingFlux = useBgcSinkingFlux_in
       useStepInitialProfiles = useStepInitialProfiles_in
       useBenthicReactions = useBenthicReactions_in
       useFluxCorrection = useFluxCorrection_in
       useOceanConc = useOceanConc_in
       useDeepSource = useDeepSource_in
       useConstantDiffusivity = useConstantDiffusivity_in

    CASE(1)

       !  'Positive fast advection test case, fixed porosity'

       useDepthDependentPorosity = .false.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .false.
       useSecondaryReactions = .false.

    CASE(2)

       ! 'Negative (resuspension) fast advection test case'

       useDepthDependentPorosity = .false.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .false.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(3)

       ! 'Positive fast advection with variable porosity test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(4)

       ! 'Positive fast advection with surface bio source and  variable porosity test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(5)

       ! 'Positive fast advection with surface Bio source + reaction and variable porosity test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .true.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(6)

        ! 'Constant Diffusion (no bottom fluxes and no solid surface fluxes balanced by Q_top and Q_bot)  test case'

       useDepthDependentPorosity = .false.
       useNonZeroDiffusivity = .true.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .false.
       useConstantDiffusivity = .true.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

   CASE(7)

      ! 'Constant Diffusion with variable porosity non-zero ocean Conc for solutes  test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .false.
       useConstantDiffusivity = .true.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(8)

       ! 'Variable Diffusion and variable porosity non-zero ocean Conc for solutes  test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useConstantDiffusivity = .false.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(9)

       ! 'Variable Diffusion, variable porosity + velocities'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useConstantDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .false.
       useCarbonateSaturation = .false.
       useSecondaryReactions = .false.

    CASE(10)

       ! 'Reactions only with variable porosity'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useConstantDiffusivity = .false.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .false.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .false.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(11)

       ! 'Reactions only with variable porosity'

       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.
       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useConstantDiffusivity = .false.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .false.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .false.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.

    CASE DEFAULT    ! 0  change useOceanConc = false (uncouple with BEC)

       ! 'Default test case'

       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.
       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .false.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .false.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .true.
       useConstantDiffusivity = .false.

 END SELECT

 useCarbonateSaturation_in = useCarbonateSaturation
 useSecondaryReactions_in = useSecondaryReactions
 useDepthDependentPorosity_in = useDepthDependentPorosity
 useNonZeroDiffusivity_in = useNonZeroDiffusivity
 useSedimentation_in = useSedimentation
 usePositiveSedimentation_in = usePositiveSedimentation
 useFastSedimentation_in = useFastSedimentation
 useBgcSinkingFlux_in = useBgcSinkingFlux
 useStepInitialProfiles_in = useStepInitialProfiles
 useBenthicReactions_in = useBenthicReactions
 useFluxCorrection_in = useFluxCorrection
 useOceanConc_in = useOceanConc
 useDeepSource_in = useDeepSource
 useConstantDiffusivity_in = useConstantDiffusivity

 fluxCorrect = c0_benthos
 if (useFluxCorrection) fluxCorrect = c1_benthos

 end subroutine benthos_flags_init

!***********************************************************************
!BOP
! !IROUTINE: benthos_SourceSink
! !INTERFACE:

 subroutine benthos_SourceSink(benthos_input, benthos_forcing, benthos_indices,&
                           benthos_output, benthos_diagnostic_fields, &
                           nBenthicVertLevels,numBenthicColumnsMax, &
                           dt,config_use_benthos_monthly_forcing,err)

! !DESCRIPTION:
!  Compute time derivatives for benthos tracers
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(benthos_input_type),  intent(in) :: benthos_input

  type(benthos_forcing_type),     intent(in) :: benthos_forcing

  type(benthos_indices_type), intent(in) :: benthos_indices

  integer (KIND=benthos_i4), intent(in) :: nBenthicVertLevels, numBenthicColumnsMax

  real (KIND=benthos_r8), intent(in) :: dt

  logical (KIND=benthos_log), intent(in) :: config_use_benthos_monthly_forcing

! !INPUT/OUTPUT PARAMETERS:

  integer (KIND=benthos_i4), intent(inout) :: err
  type(benthos_output_type), intent(inout) :: benthos_output
  type(benthos_diagnostics_type), intent(inout) :: benthos_diagnostic_fields

! !!OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  integer (KIND=benthos_i4), parameter :: &
       nOMtype = 3, &
       column = 1

  logical(KIND=benthos_log) :: &
       l_stop, &
       tooLarge, &
       lcomp_co3_coeff_benthos, &
       lcalc_co2_terms, &
       test_flag

  real(KIND=benthos_r8) :: &
       dts,              &       ! timestep
       dts_o,            &       ! timestep if subcycling
       oceanBottomTemperature, & ! oC
       oceanBottomSalinity, & !  (g salt/kg seawater)
       oceanBottomPressure, & ! N/m2
       oceanBottomSilicate, & !
       oceanBottomDepth, &
       oceanBottomDensity, &
       oceanBottomPhosphate, &
       oceanBottomAmmonium, &
       oceanBottomNitrate, &
       oceanBottomOxygen, &
       oceanBottomManganese, &
       oceanBottomIron, &
       oceanBottomDIC, &
       oceanBottomAlkalinity, &
       oceanPOCFlux, &
       oceanPONFlux, &
       oceanPOPFlux, &
       oceanParticulateIronFlux, &
       oceanCalciteFlux, &
       oceanBiogenicSilicaFlux, &
       oceanDustFlux, &
       terraPOCFlux, &
       lithogenicFlux

  real(KIND=benthos_r8) :: &
       work1, work2, work3, work4, work5, & !
       FluxSed,      & ! sediment sinking flux (kg/m2/s)
       dhtop,          & ! change in benthos upper surface (m)
       vel,              &  ! velocity of upper surface (m/s)
       C_topm,      & ! surface boundary concentration for solute or solid
       C_botm,      & ! bottom boundary concentration for solute or solid
       C_topSolid, & !
       sum_old, &
       sum_new, &
       sum_tot, &
       sum_init, &
       totalOceanSolidFluxm, & !
       benthosBottomFluxm, &
       dt_old, &
       del_ph, &
       phlo_3d_init, &
       phhi_3d_init, &
       reactionTend, &
       reactionTendprimary, &
       reactionTendsecondary, &
       reactionTendcarbonate, &
       trueTend, &
       tempR, &
       tempSol, &
       tempS, &
       benthosTracerBeforeReaction

  real(KIND=benthos_r8), dimension(:), allocatable :: &
       oceanSedFluxdt, & ! particulate flux * dt of sinking ocean tracers (mmol/kg m/s dt)
       oceanDiffVelFluxdt, & ! diffusive/vel surface boundary flux (mmol/kg or mmol/m3) * m/s dt
       benthosSedFluxdt, & ! particulate flux to deep storage (mmol/kg or mmol/m3) * m/s dt
       benthosDiffVelFluxdt, & ! diffusive/vel bottom boundary flux (mmol/kg or mmol/m3) * m/s * dt
       diffError, & ! mmol/m3 or mmol/kg
       totFluxesdt, & ! sum of tracer boundary fluxes * dt
       oceanSedFlux,    & ! particulate flux of sinking ocean tracers (mmol/kg*m/s)
       oceanTracerConc, & ! bottom ocean concentration of solutes (mmol/m3)
       oceanPrecipFlux, & ! flux of precipitated material from ocean bottom (mmol/kg*m/s)
       totalOceanSolidFlux, & ! total solid flux divided by benthosDepth (mmol/kg/s)
       deepBenthosBottomFlux, &  ! prescribed deep storage flux (mmol/m2/s)
       netElementsStart, & ! initial column element concentration (mmol/m2)
       totalColumnTracersInitial, &  ! initial column integrated tracer (mmol/m2)
       totalColumnTracersFinal, &  ! final column integrated tracer (mmol/m2)
       bDiff, & ! fluid or solid Diffusivity on the benthos grid
       bpor, & ! fluid or solid porosity on the benthos grid
       iDiff, & ! fluid or solid Diffusivity on the interfaces
       ipor, &  ! fluid or solid porosities on the interfaces
       C_top_vel, & ! (mmol/m3) surface concentration for velocity boundary condition
       C_bot_vel, & ! (mmol/m3) bottom concentration for velocity boundary condition
       C_top, & ! (mmol/m3) surface concntration for diffusive boundary condition
       C_bot, & ! (mmol/m3) bottom concentration for diffusive boundary condition
       source, & !
       bottomsource, &
       C_tot, &
       rel_error, &
       rel_error_tot, &
       PH_PREV_3D, &
       netElementsRelError
!       netTransport_tend

  real(KIND=benthos_r8), dimension(nBenthicTracers) :: &
       Source_top, &
       Source_bot, &
       Sink_bot_d, &
       Sink_bot_s, &
       Sink_top_d, &
       Sink_top_s, &
       flux_bio, &
       benthosStorageConc

  real(KIND=benthos_r8), dimension(nElements) :: &
       oceanElementExchange, &
       burialElementExchange, &
       netElements_end_Transport, &
       oceanElementSedimentation, &
       deepBurial, &
       netElements_start_reactions, &
       netElements_after_ph_adjustment, &
       netElements_start_reactions_tmp, &
       netElements_end_reactions_tmp, &
       netElements_end_reactions

  real(KIND=benthos_r8), dimension(nCarbonateReactions) :: &
       carbonateRates

  real(KIND=benthos_r8), dimension(nSecondaryReactions) :: &
       secondaryRates

  real(KIND=benthos_r8), dimension(nOMtype,nPrimaryReactions) :: &
       primaryRates

  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       initcons, & ! array for initial tracer concentration
       biocons, & ! array for tracer concentration
       biomat_low, &
       D_spdiag, D_sbdiag, ML_diag, &
       ML_sbdiag, ML_spdiag, &
       rhs, spdiag, diag, sbdiag, &
       sat_calc, sat_arag, sat_mgcalc, &
       K_calc, K_arag, K_mgcalc

  real(KIND=benthos_r8), dimension(:,:), allocatable :: &
       benthosTracerBulk, &    ! local tracer concentration
       initBenthosTracerBulk, &
       initBenthosTracerBulk_Carbonate, &
       primarySourceTend, &
       primarySinkTend, &
       carbonateSourceTend, &
       carbonateSinkTend, &
       secondarySourceTend, &
       secondarySinkTend, &
       netReactionTend

  integer (KIND=benthos_i4) :: &
       iLevels, &
       iTracers, &
       iElements, &
       ii, &
       subtt, &
       subN, &
       subN_o, &
       cnt_sub

  real (KIND=benthos_r8) :: &
       co3_mol, &   ! comp_CO3terms
       hco3_mol, &
       co2_mol, &
       phlo, &
       phhi, &
       ph, &
       kw, &
       kb, &
       ks, &
       kf, &
       k1p, &
       k2p, &
       k3p, &
       ksi, &
       bt, &
       st, &
       ft, &
       dic, &
       ta, &
       pt, &
       sit, &
       error

!-----------------------------------------------------------------------
!EOC
!-----------------------------------------------------------------------
!  write(*,*)'begin benthos_SourceSink'

  ! initialize arrays

   allocate(deepBenthosBottomFlux(nBenthicTracers))
   allocate(oceanSedFlux(nBenthicTracers))
   allocate(oceanPrecipFlux(nBenthicTracers))
   allocate(totalOceanSolidFlux(nBenthicTracers))
   allocate(oceanTracerConc(nBenthicTracers))
   allocate(netElementsStart(nElements))
   allocate(netElementsRelError(nElements))
   allocate(benthosTracerBulk(nBenthicTracers, nBenthicVertLevels+1))
   allocate(initBenthosTracerBulk(nBenthicTracers, nBenthicVertLevels+1))
   allocate(initBenthosTracerBulk_Carbonate(nBenthicTracers, nBenthicVertLevels+1))
   allocate(totalColumnTracersInitial(nBenthicTracers))
   allocate(totalColumnTracersFinal(nBenthicTracers))
   allocate(bDiff(nBenthicVertLevels+2))
   allocate(iDiff(nBenthicVertLevels+1))
   allocate(bpor(nBenthicVertLevels+2))
   allocate(ipor(nBenthicVertLevels+1))
   allocate(C_top_vel(nBenthicTracers))
   allocate(C_bot_vel(nBenthicTracers))
   allocate(C_top(nBenthicTracers))
   allocate(C_bot(nBenthicTracers))
   allocate(source(nBenthicTracers))
   allocate(bottomsource(nBenthicTracers))
   allocate(C_tot(nBenthicTracers))
!   allocate(netTransport_tend(nBenthicTracers))
   allocate(oceanSedFluxdt(nBenthicTracers))
   allocate(oceanDiffVelFluxdt(nBenthicTracers))
   allocate(benthosSedFluxdt(nBenthicTracers))
   allocate(benthosDiffVelFluxdt(nBenthicTracers))
   allocate(diffError(nBenthicTracers))
   allocate(totFluxesdt(nBenthicTracers))
   allocate(primarySourceTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(primarySinkTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(carbonateSourceTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(carbonateSinkTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(secondarySourceTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(secondarySinkTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(netReactionTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(rel_error(nElements))
   allocate(rel_error_tot(nElements))
   allocate(PH_PREV_3D(nBenthicVertLevels+1))

   err = 0

   do iTracers = 1,nBenthicTracers
      deepBenthosBottomFlux(iTracers) = c0_benthos !   Bottom boundary condition  (mmol/m2/s)
      oceanSedFlux(iTracers)          = c0_benthos
      oceanPrecipFlux (iTracers)      = c0_benthos
      totalOceanSolidFlux(iTracers)   = c0_benthos
      oceanTracerConc(iTracers)       = c0_benthos
   end do
   netElementsRelError(:) = c0_benthos
   rel_error(:) = c0_benthos
   rel_error_tot(:) = c0_benthos

   !Initialize local variables
   oceanBottomPhosphate  = benthos_input%oceanBottomPhosphate(column)
   oceanBottomAmmonium  = benthos_input%oceanBottomAmmonium(column)
   oceanBottomNitrate  = benthos_input%oceanBottomNitrate(column)
   oceanBottomOxygen  = benthos_input%oceanBottomOxygen(column)
   oceanBottomIron  = benthos_input%oceanBottomIron(column)
   oceanBottomDIC  = benthos_input%oceanBottomDIC(column)
   oceanBottomManganese = benthos_input%oceanBottomManganese(column)
   oceanBottomAlkalinity  = benthos_input%oceanBottomAlkalinity(column)
   ! Particulate fluxes are in mmol/m3*m/s
   oceanPOCFlux  = benthos_input%oceanPOCFlux(column)
   oceanPONFlux  = benthos_input%oceanPOCFlux(column)/CtoP*NtoP
   oceanPOPFlux  = benthos_input%oceanPOCFlux(column)/CtoP
   oceanParticulateIronFlux  = benthos_input%oceanParticulateIronFlux(column)
   oceanCalciteFlux  = benthos_input%oceanCalciteFlux(column)
   oceanBiogenicSilicaFlux  = benthos_input%oceanBiogenicSilicaFlux(column)
   oceanDustFlux  = benthos_input%oceanDustFlux(column)
   oceanBottomSilicate = benthos_input%oceanBottomSilicate(column)
   terraPOCFlux = benthos_input%terraPOCFlux(column)
   lithogenicFlux = benthos_input%lithogenicFlux(column)

   if (useOceanBottomState .or. config_use_benthos_monthly_forcing) then
      oceanBottomDepth = benthos_input%oceanBottomDepth(column)
      oceanBottomTemperature = benthos_input%oceanBottomTemperature(column)
      oceanBottomSalinity = benthos_input%oceanBottomSalinity(column)
      oceanBottomPressure = benthos_input%oceanBottomPressure(column)
      oceanBottomDensity = benthos_input%oceanBottomDensity(column)

      !NJ-TEST: VERIFIED
      test_flag = .false.
      if (test_flag) then
         if (oceanBottomDepth .lt. 50.0_benthos_r8) then
            write(*,*) 'ocean_bottom_state: oceanBottomDepth:',oceanBottomDepth
            write(*,*) 'oceanBottomTemperature:',oceanBottomTemperature
            write(*,*) 'oceanBottomSalinity:',oceanBottomSalinity
            write(*,*) 'oceanBottomPressure:',oceanBottomPressure
            write(*,*) 'oceanBottomDensity:',oceanBottomDensity
         end if
      end if
   else
      oceanBottomDepth = 53.0_benthos_r8
      oceanBottomTemperature = c0_benthos
      oceanBottomSalinity = 32.0_benthos_r8     ! Chukchi
      oceanBottomPressure = 4.236e5_benthos_r8 ! N/m2
      oceanBottomDensity = 1027.0_benthos_r8    ! kg/m3
   end if

   benthos_output%oceanBottomTemperature(column) = oceanBottomTemperature
   benthos_output%oceanBottomDepth(column) = oceanBottomDepth
   benthos_output%oceanBottomSalinity(column) = oceanBottomSalinity
   benthos_output%oceanBottomPressure(column) = oceanBottomPressure
   benthos_output%oceanBottomDensity(column) = oceanBottomDensity

   do iTracers = 1, nBenthicTracers
      benthosStorageConc(iTracers) = benthos_input%deepStorage(column,iTracers)
      do iLevels = 1, nBenthicVertLevels+1
         ! initialize output tendencies
         benthos_output%benthosTendencies(iLevels,column,iTracers) = c0_benthos
         benthos_output%benthosTransportTendencies(iLevels,column,iTracers) = c0_benthos
         benthos_output%benthosReactionTendencies(iLevels,column,iTracers) = c0_benthos

         benthosTracerBulk(iTracers,iLevels) = benthos_input%benthosTracerBulk(iLevels,column,iTracers)
         if (benthosTracerBulk(iTracers,iLevels) .lt. c0_benthos) then
            write(*,*) 'benthosTracerBulk < 0 at start of SourceSink'
            write(*,*) 'benthosTracerBulk(iTracers,iLevels):',benthosTracerBulk(iTracers,iLevels)
            write(*,*) 'iTracers,iLevels,column:',iTracers,iLevels,column
            if (benthosTracerBulk(iTracers,iLevels) .gt. -puny) then
               write(*,*) 'small negative. setting tracer to zero'
               benthosTracerBulk(iTracers,iLevels) = c0_benthos
            else
               err = 1
               return
            end if
         end if
      end do
   end do

   ! update time-dependent parameters (temperature dependent diffusion)  NJ-Verified
   call benthos_parameters_dt (nBenthicVertLevels, oceanBottomTemperature, oceanBottomSalinity, &
        oceanBottomSalinity,benthos_indices)

   ! Initialize the bottom flux boundary conditions
   call updateBenthosBoundaryConditions (oceanSedFlux, oceanPrecipFlux, totalOceanSolidFlux, &
        oceanTracerConc, deepBenthosBottomFlux, oceanBottomDensity, nBenthicVertLevels, &
        FluxSed, dhtop, vel, dt, oceanBottomPhosphate, oceanBottomAmmonium, oceanBottomNitrate, &
        oceanBottomOxygen, oceanBottomIron, oceanBottomDIC, oceanBottomAlkalinity, oceanPOCFlux, &
        oceanPONFlux, oceanPOPFlux, oceanParticulateIronFlux, oceanCalciteFlux,benthosSolidPorosityI(1), &
        oceanBiogenicSilicaFlux, oceanDustFlux, oceanBottomSilicate, oceanBottomManganese, &
        benthos_output, column, config_use_benthos_monthly_forcing, oceanBottomDepth, &
        terraPOCFlux, lithogenicFlux)

   call computeNetElements (netElementsStart, nBenthicVertLevels, benthosTracerBulk,&
        vertBenthosGridThickI)

   ! NJ:Validated
   !write(*,*) 'After computeNetElements ... netElementsStart:',netElementsStart

   do ii = 1,nBenthicTracers-3
      ! co3, co2 and hco3 is transported as dic
      iTracers = transport_index(ii)
      call sum_benthos_column (totalColumnTracersInitial(iTracers),benthosTracerBulk(iTracers,:),nBenthicVertLevels, vertBenthosGridThickI)

      bDiff(nBenthicVertLevels+2) = benthosSoluteDiffusivity(iTracers,nBenthicVertLevels+2)*&
           real(benthos_input%tracerType(iTracers),KIND=benthos_r8) +  &
           benthosSolidDiffusivity(nBenthicVertLevels+2)*&
           real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)

      bpor(nBenthicVertLevels+2) = benthosPorosity(nBenthicVertLevels+2)*&
           real(benthos_input%tracerType(iTracers),KIND=benthos_r8) +  &
           benthosSolidPorosity(nBenthicVertLevels+2)*&
           real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)

      do iLevels = 1,nBenthicVertLevels+1

!         iDiff(iLevels) = benthosDiffusivityI(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
!              benthosSolidDiffusivityI(iLevels)*(c1_benthos - real(benthos_input%tracerType(iTracers),KIND=benthos_r8))
!         bDiff(iLevels) = benthosDiffusivity(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
!              benthosSolidDiffusivity(iLevels)*(c1_benthos-real(benthos_input%tracerType(iTracers),KIND=benthos_r8))


         iDiff(iLevels) = benthosSoluteDiffusivityI(iTracers,iLevels)* &
              real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidDiffusivityI(iLevels)*real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)
         bDiff(iLevels) = benthosSoluteDiffusivity(iTracers,iLevels)* &
              real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidDiffusivity(iLevels)*real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)

         ipor(iLevels) = benthosPorosityI(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidPorosityI(iLevels)* real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)
         bpor(iLevels) = benthosPorosity(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidPorosity(iLevels)*real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)

         initcons(iLevels) = benthosTracerBulk(iTracers,iLevels)
         biocons(iLevels) = benthosTracerBulk(iTracers,iLevels)
      end do

      C_top_vel(iTracers) = benthosTracerBulk(iTracers,1)
      C_bot_vel(iTracers) = benthosTracerBulk(iTracers,nBenthicVertLevels+1)
      C_top(iTracers) = oceanTracerConc(iTracers)
      C_bot(iTracers) = C_bot_vel(iTracers)
      C_botm = C_bot(iTracers)
      C_topm = C_top(iTracers)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + C_top_vel(iTracers)*real((1-benthos_input%tracerType(iTracers)),KIND=benthos_r8)

      totalOceanSolidFluxm = totalOceanSolidFlux(iTracers)
      benthosBottomFluxm = deepBenthosBottomFlux(iTracers)/benthosDepth  ! normalized   positive down
      source(iTracers) = totalOceanSolidFlux(iTracers)    ! from ocean model
      bottomsource(iTracers) = deepBenthosBottomFlux(iTracers)/benthosDepth   ! positive down  *Prescribed*

      !NJ-TEST: VERIFIED
      test_flag = .false.
      if (test_flag .and. iTracers .eq. 12) then
      write(*,*) 'iTracers:',iTracers
      write(*,*) 'Before compute_FCT_matrix, totalOceanSolidFluxm:', totalOceanSolidFluxm
      write(*,*) 'benthosBottomFluxm:', benthosBottomFluxm
      write(*,*) 'oceanTracerConc(iTracers):',oceanTracerConc(iTracers)
      write(*,*) 'source(iTracers):', source(iTracers)
      write(*,*) 'bottomsource(iTracers):', bottomsource(iTracers)
      write(*,*) 'C_top_vel(iTracers):', C_top_vel(iTracers)
      write(*,*) 'C_bot_vel(iTracers):', C_bot_vel(iTracers)
      write(*,*) 'C_top(iTracers):', C_top(iTracers)
      write(*,*) 'C_bot(iTracers):', C_bot(iTracers)
      write(*,*) 'C_topm:', C_topm
      write(*,*) 'C_botm:', C_botm
      write(*,*) 'benthos_input%tracerType(iTracers):',benthos_input%tracerType(iTracers)

      write(*,*) 'biocons:', biocons
      write(*,*) 'bDiff:', bDiff
      write(*,*) 'iDiff:', iDiff
      write(*,*) 'bpor:', bpor
      write(*,*) 'ipor:', ipor
      write(*,*) 'benthosSolidDiffusivity:', benthosSolidDiffusivity
      write(*,*) 'benthosSolidDiffusivityI:', benthosSolidDiffusivityI
      err = 1
      return
      end if
      call compute_FCT_matrix_CN (sbdiag, diag, spdiag, rhs, ML_diag, ML_sbdiag, ML_spdiag, &
           D_sbdiag, D_spdiag, &
           Source_top(iTracers), Source_bot(iTracers), Sink_bot_d(iTracers), &
           Sink_bot_s(iTracers), &
           Sink_top_d(iTracers), Sink_top_s(iTracers),  &
           biocons, dt, nBenthicVertLevels,  dhtop,  ipor,  iDiff, benthosDepth,    &
           totalOceanSolidFluxm, bpor,   C_topm, C_botm,   C_top_vel(iTracers), &
           C_bot_vel(iTracers),  benthosBottomFluxm, vel, bDiff,vertBenthosGridThickI, &
           vertBenthosGridThick)

      call tridiag_solverz (biocons, nBenthicVertLevels+1, sbdiag,  diag, spdiag, rhs)

      test_flag = .false.
      if (test_flag .and. iTracers .eq. 12) then !biocons(1) .gt. initcons(1) .and. iTracers .eq. 12) then
      write(*,*) 'iTracers:',iTracers
      write(*,*) 'after tridiag, totalOceanSolidFluxm:', totalOceanSolidFluxm
      write(*,*) 'benthosBottomFluxm:', benthosBottomFluxm
      write(*,*) 'sbdiag:', sbdiag
      write(*,*) 'spdiag:', spdiag
      write(*,*) 'diag:', diag
      write(*,*) 'ML_sbdiag:', ML_sbdiag
      write(*,*) 'ML_spdiag:', ML_spdiag
      write(*,*) 'ML_diag:', ML_diag
      write(*,*) 'rhs:', rhs
      write(*,*) 'dhtop:', dhtop
      write(*,*) 'Sink_top_d(iTracers):', Sink_top_d(iTracers)
      write(*,*) 'Sink_top_s(iTracers):', Sink_top_s(iTracers)
      write(*,*) 'oceanTracerConc(iTracers):',oceanTracerConc(iTracers)
      write(*,*) 'source(iTracers):', source(iTracers)
      write(*,*) 'bottomsource(iTracers):', bottomsource(iTracers)
      write(*,*) 'C_top_vel(iTracers):', C_top_vel(iTracers)
      write(*,*) 'C_bot_vel(iTracers):', C_bot_vel(iTracers)
      write(*,*) 'C_top(iTracers):', C_top(iTracers)
      write(*,*) 'C_bot(iTracers):', C_bot(iTracers)
      write(*,*) 'C_topm:', C_topm
      write(*,*) 'C_botm:', C_botm
      write(*,*) 'benthos_input%tracerType(iTracers):',benthos_input%tracerType(iTracers)
      write(*,*) 'biocons:', biocons
      write(*,*) 'initcons:', initcons
      write(*,*) 'benthosMolecularDiffusivity(iTracers):', benthosMolecularDiffusivity(iTracers)
      write(*,*) 'bDiff:', bDiff
      write(*,*) 'iDiff:', iDiff
      write(*,*) 'bpor:', bpor
      write(*,*) 'ipor:', ipor
      write(*,*) 'benthosSolidDiffusivity:', benthosSolidDiffusivity
      write(*,*) 'benthosSolidDiffusivityI:', benthosSolidDiffusivityI
      err = 1
      return
      end if
      ! Better
      ! write(*,*) 'after tridiag_solverz, biocons:', biocons

      call check_conservation_FCT (flux_bio(iTracers), l_stop, benthosStorageConc(iTracers), &
           oceanSedFluxdt(iTracers),oceanDiffVelFluxdt(iTracers),benthosSedFluxdt(iTracers), &
           benthosDiffVelFluxdt(iTracers),diffError(iTracers),totFluxesdt(iTracers), initcons,&
           biocons, Source_top(iTracers), Source_bot(iTracers), Sink_bot_d(iTracers), &
           Sink_bot_s(iTracers), Sink_top_d(iTracers), Sink_top_s(iTracers), dt, &
           nBenthicVertLevels, bottomsource(iTracers), vel, benthosDepth,        &
           iTracers,benthosBottomFluxm,totalOceanSolidFluxm,vertBenthosGridThickI, iDiff, oceanBottomTemperature)

      if (l_stop) then
         err = 1
         return
      end if
      do iLevels = 1, nBenthicVertLevels+1
         biomat_low(iLevels) = biocons(iLevels)
      end do

 !     write(*,*) 'before compute_FCT_corr'

      call compute_FCT_corr (initcons,  &
                                 biocons, dt, nBenthicVertLevels, &
                                 D_sbdiag, D_spdiag, ML_diag, ML_spdiag, &
                                 ML_sbdiag, vertBenthosGridThickI)

      test_flag = .false.
      if (test_flag .and. biocons(1) .gt. initcons(1) .and. iTracers .eq. 12) then
      write(*,*) 'iTracers:',iTracers
      write(*,*) 'after FCT_corr, initcons:',initcons
      write(*,*) 'biocons:', biocons
      end if
      do iLevels = 1, nBenthicVertLevels+1
         benthosTracerBulk(iTracers,iLevels) = max(c0_benthos,biocons(iLevels))
      end do

      sum_old = c0_benthos
      sum_new = c0_benthos
      sum_tot = c0_benthos
      sum_init = c0_benthos

      call sum_benthos_column &
          (sum_old,biomat_low,nBenthicVertLevels,vertBenthosGridThickI)

      call sum_benthos_column &
          (sum_new,biocons,nBenthicVertLevels,vertBenthosGridThickI)

      call sum_benthos_column &
          (sum_tot,benthosTracerBulk(iTracers,:),nBenthicVertLevels,vertBenthosGridThickI)

      call sum_benthos_column &
          (sum_init,initcons,nBenthicVertLevels,vertBenthosGridThickI)

      error = max(maxval(biomat_low),maxval(biocons))
      error = max(error, max(sum_new,sum_old))*puny

      if ((abs(sum_new-sum_old) .gt. error) .or. &
          (minval(biocons) .lt. c0_benthos)   .or.  l_stop .or. (sum_new .gt. 1e12_benthos_r8)) then
         write(*,*) 'Benthic FCT tracer solution failed,iTracers,iLevels', iTracers,iLevels
         write(*,*)'or Benthic tracer blowing up after FCT,iTracers', iTracers
         write(*,*) 'sum_new,sum_old:',sum_new,sum_old
         write(*,*) 'l_stop:',l_stop
         write(*,*) 'error:',error
         write(*,*) 'maxval(biocons):',maxval(biocons)
         write(*,*) 'maxval(biomat_low):',maxval(biomat_low)
         write(*,*) 'puny*max(sum_old,sum_new):',puny*max(sum_old,sum_new)
         write(*,*) 'abs(sum_old-sum_new):',abs(sum_old-sum_new)
         write(*,*) 'iTracers,biocons(:):',iTracers,biocons(:)
         write(*,*) 'biomat_low:',biomat_low
         write(*,*) 'iDiff(:):',iDiff(:)
         write(*,*)  'initcons(:):',initcons(:)
         write(*,*) 'benthosTracerBulk(iTracers,:):',benthosTracerBulk(iTracers,:)
         write(*,*) 'dhtop',dhtop
         write(*,*)'source(iTracers):', source(iTracers)
         write(*,*)'bottomsource(iTracers):', bottomsource(iTracers)
         !stop_label = 'zbgc FCT tracer solution failed'
         err = 1
         return
       end if

       ! don't include transport tendencies for individual CO2 speciation

       do iLevels = 1,nBenthicVertLevels+1
          benthos_output%benthosTransportTendencies(iLevels,column,iTracers) = &
               (benthosTracerBulk(iTracers,iLevels) - benthos_input%benthosTracerBulk(iLevels,column,iTracers))/dt

          benthos_output%benthosTendencies(iLevels,column,iTracers) = &
               benthos_output%benthosTransportTendencies(iLevels,column,iTracers)
      end do  ! iLevels

      call sum_benthos_column (C_tot(iTracers),benthosTracerBulk(iTracers,:),nBenthicVertLevels,vertBenthosGridThickI)
      call sum_benthos_column (benthos_diagnostic_fields % diag_netBenthosTransportTend(column,iTracers), &
           benthos_output%benthosTransportTendencies(:,column,iTracers),nBenthicVertLevels, &
           vertBenthosGridThickI)

      benthos_diagnostic_fields % relErrorTransport(column,iTracers) = diffError(iTracers)/&
           (C_tot(iTracers) + c1_benthos)

      do iElements = 1,nElements
         netElementsRelError(iElements) = netElementsRelError(iElements) + &
              (elementRatios(iTracers,iElements)*diffError(iTracers))**2
      end do
   end  do !iTracers

   call computeElementFluxes (oceanElementExchange, burialElementExchange, oceanSedFluxdt, &
        oceanDiffVelFluxdt, benthosSedFluxdt, benthosDiffVelFluxdt,benthosDepth)

   call computeNetElements (netElements_end_Transport, nBenthicVertLevels, benthosTracerBulk,vertBenthosGridThickI)

   call netElementExchange (oceanElementSedimentation, deepBurial, totalOceanSolidFlux, &
        benthosStorageConc,dt,benthosDepth)

   !NJ-TEST:
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'computeElementFluxes, oceanElementExchange:',oceanElementExchange
      write(*,*) 'burialElementExchange:',burialElementExchange
      write(*,*) 'netElements_end_Transport:',netElements_end_Transport
      write(*,*) 'oceanElementSedimentation:',oceanElementSedimentation
   end if

   !
   ! Reaction terms : all the tracers have been updated for
   !  advection and transport.  Now define the reaction tendencies
   !

   ! initialize
   dts_o = 300.0_benthos_r8
   dt_old = dt
   subN_o = NINT(dt/dts_o,KIND=benthos_i4) !p1_benthos
   tooLarge = .true.

   !  make sure there is consistency to start:
   do iLevels = 1, nBenthicVertLevels+1
      !initialize ph for reactions
      PH_PREV_3D(iLevels) = benthos_input%PH_PREV_3D(iLevels, column)

      do iTracers = 1, nBenthicTracers
         initBenthosTracerBulk(iTracers,iLevels) = benthosTracerBulk(iTracers,iLevels)
         initBenthosTracerBulk_Carbonate(iTracers,iLevels) = benthosTracerBulk(iTracers,iLevels)
      end do
   end do

   !
   !  Save initial mmols of each ion
   ! Count the total mols of the following ions:
   ! In the order: C, O, N, P, S, Mn, Fe

   call computeNetElements (netElements_start_reactions,nBenthicVertLevels,&
        benthos_input%benthosTracerBulk(:,column,:),vertBenthosGridThickI)

   if (useBenthicReactions) then

      !-------------------------------
      !  Carbonate chemistry
      !-------------------------------

      ! equilibrium concentrations
      lcomp_co3_coeff_benthos = .true.
      lcalc_co2_terms = .true.
      del_ph = 0.2_benthos_r8
      phlo_3d_init = 6.0_benthos_r8
      phhi_3d_init = 9.0_benthos_r8

      !sumC_initial = benthicTracerBulk(dic_ind,iLevels) %benthicTracerBulk(co3_ind,iLevels) + benthicTracerBulk(co2_ind,iLevels) + benthicTracerBulk(hco3_ind,iLevels)
      !sumco3_initial = benthicTracerBulk(co3_ind,iLevels)
      !sumco2_initial = benthicTracerBulk(co2_ind,iLevels)
      !sumhco3_initial = benthicTracerBulk(hco3_ind,iLevels)

      if (useCarbonateSaturation) then

         do iLevels = 1, nBenthicVertLevels+1
            work3 = c0_benthos
            work4 = c0_benthos

            if (lcalc_co2_terms) then
                if (PH_PREV_3D(iLevels) .ne. c0_benthos) then
                   work1 = PH_PREV_3D(iLevels) - del_ph
                   work2 = PH_PREV_3D(iLevels) + del_ph
                else
                   work1 = phlo_3d_init
                   work2 = phhi_3d_init
                end if

                work5 = oceanBottomDepth + benthosMidPointI(iLevels)   ! m
                ! confirmed
                call comp_CO3terms_benthos (co3_mol,hco3_mol,co2_mol,ph,kw,kb, &
                   ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,dic,ta,pt,sit, oceanBottomTemperature,  &
                   oceanBottomSalinity,work5, &
                   lcomp_co3_coeff_benthos, oceanBottomSilicate, work1, work2, &
                   benthosTracerBulk(dic_ind,iLevels), benthosTracerBulk(alk_ind,iLevels), &
                   benthosTracerBulk(h2po4_ind,iLevels), iLevels, oceanBottomDensity, &
                   benthosPorosityI(iLevels))

                work3 = ph   ! do i need this?
                PH_PREV_3D(iLevels) = ph
            else
                co2_mol  = c0_benthos !H2CO3 = 0
                hco3_mol = c0_benthos !HCO3  = 0
                co3_mol  = c0_benthos !CO3   = 0
                PH_PREV_3D(iLevels) = 8.0_benthos_r8
            end if

            benthosTracerBulk(alk_ind,iLevels) = ta
!            benthosTracerBulk(h2po4_ind,iLevels) = pt
!            oceanBottomSilicate = sit

            benthosTracerBulk(co3_ind,iLevels) = co3_mol
            benthosTracerBulk(co2_ind,iLevels) = co2_mol
            benthosTracerBulk(hco3_ind,iLevels) = hco3_mol
            benthosTracerBulk(dic_ind,iLevels) = dic

            ! saturation state

            call  comp_co3_sat_benthos(sat_calc(iLevels), sat_arag(iLevels), sat_mgcalc(iLevels),&
                K_arag(iLevels),K_calc(iLevels),K_mgcalc(iLevels),oceanBottomTemperature,oceanBottomSalinity,&
                work5, oceanMagnesium,co3_mol,oceanBottomDensity,benthosPorosityI(iLevels))

            do iTracers = 1, nBenthicTracers
               initBenthosTracerBulk_Carbonate(iTracers,iLevels) = max(c0_benthos,benthosTracerBulk(iTracers,iLevels))
            end do
         end do

         call computeNetElements (netElements_after_ph_adjustment,nBenthicVertLevels,&
          initBenthosTracerBulk_Carbonate(:,:),vertBenthosGridThickI)

      end if ! useCarbonateSaturation
      cnt_sub = 0
      do while (tooLarge)
         cnt_sub = cnt_sub+1
         if (cnt_sub == 1) then
            subN = 1
            dts = dt
         else
            subN = subN_o
            dts = dts_o
         end if
         tooLarge = .false.
         do iTracers = 1, nBenthicTracers
            do iLevels = 1, nBenthicVertLevels+1
               netReactionTend(iTracers,iLevels) = c0_benthos
               benthosTracerBulk(iTracers,iLevels) = initBenthosTracerBulk_Carbonate(iTracers,iLevels)
            end do
         end do
         rel_error(:) = c0_benthos
         rel_error_tot(:) = c0_benthos
         do subtt = 1,subN
            do iLevels = 1,nBenthicVertLevels+1

               call computeNetElements (netElements_start_reactions_tmp, nBenthicVertLevels, benthosTracerBulk,vertBenthosGridThickI)

               !  Primary reactions
               do iTracers = 1,nBenthicTracers
                  primarySourceTend(iTracers,iLevels) = c0_benthos
                  primarySinkTend(iTracers,iLevels) = c0_benthos
                  carbonateSourceTend(iTracers,iLevels) = c0_benthos
                  carbonateSinkTend(iTracers,iLevels) = c0_benthos
                  secondarySourceTend(iTracers,iLevels) = c0_benthos
                  secondarySinkTend(iTracers,iLevels) = c0_benthos
                  if (subtt == 1) netReactionTend(iTracers,iLevels) = c0_benthos
               end do

               call primaryStoichMatrix (primarySourceStoich, primarySinkStoich, &
                    benthosTracerBulk(:,iLevels))

               call primaryRateConstants (primaryRates,benthosTracerBulk(:,iLevels),nOMtype)

               call primaryBenthosReactions (primarySourceTend(:,iLevels),primarySinkTend(:,iLevels),&
                    benthosTracerBulk(:,iLevels),primarySourceStoich,&
                    primarySinkStoich,primaryRates,dts,nOMtype)

               ! carbonate reactions

               if (useCarbonateSaturation) then

                 call carbonateRateConstants (carbonateRates, &
                   benthosTracerBulk(:,iLevels),dts,sat_calc(iLevels), sat_arag(iLevels),&
                   sat_mgcalc(iLevels),oceanBottomDensity)

                  call carbonateBenthosReactions (carbonateSourceTend(:,iLevels),carbonateSinkTend(:,iLevels),&
                    carbonateSourceStoich,carbonateSinkStoich,carbonateRates,dts)

               end if

               ! Secondary Reactions

               if (useSecondaryReactions) then

                  call secondaryRateConstants(secondaryRates, benthosTracerBulk(:,iLevels))

                  call secondaryBenthosReactions(secondarySourceTend(:,iLevels),secondarySinkTend(:,iLevels), &
                       secondarySourceStoich,secondarySinkStoich,secondaryRates,dts)

               end if

               do iTracers = 1,nBenthicTracers

                  reactionTend = c0_benthos
                  trueTend = c0_benthos
                  tempR = c0_benthos
                  tempSol = c0_benthos
                  tempS = c0_benthos

                  if (benthosTracerBulk(iTracers,iLevels) .ge. puny) then
                     tempR = (primarySinkTend(iTracers,iLevels) + secondarySinkTend(iTracers,iLevels)+&
                          carbonateSinkTend(iTracers,iLevels))/(benthosTracerBulk(iTracers,iLevels) )
                  else
                     tempR = c0_benthos
                  end if

                  tempS = (primarySourceTend(iTracers,iLevels) + &
                       secondarySourceTend(iTracers,iLevels)+ &
                       carbonateSourceTend(iTracers,iLevels))

                  if (tempR .lt. puny) then
                     tempSol = max(c0_benthos,benthosTracerBulk(iTracers,iLevels)  + tempS* (c1_benthos - tempR/c2_benthos))
                  else
                     tempSol = max(c0_benthos,benthosTracerBulk(iTracers,iLevels)*exp(-tempR) + tempS*(c1_benthos-exp(-tempR))/tempR)
                  end if

                  trueTend =  (tempSol-benthosTracerBulk(iTracers,iLevels))

                  reactionTendprimary = (primarySourceTend(iTracers,iLevels)-primarySinkTend(iTracers,iLevels))
                  reactionTendsecondary = (secondarySourceTend(iTracers,iLevels)-secondarySinkTend(iTracers,iLevels))
                  reactionTendcarbonate =  (carbonateSourceTend(iTracers,iLevels)-carbonateSinkTend(iTracers,iLevels))

                  !benthos_output % benthosReactionTendencies(iLevels, column, iTracers) = trueTend

                  ! For testing  elemental conservation of reaction set at each iLevels
                  ! cumulative

                  netReactionTend(iTracers,iLevels) = netReactionTend(iTracers,iLevels) + trueTend !reactionTendsecondary + reactionTendprimary+ reactionTendcarbonate ! sub over subcycles

                  benthosTracerBeforeReaction = benthosTracerBulk(iTracers,iLevels)
                  benthosTracerBulk(iTracers,iLevels) = benthosTracerBulk(iTracers,iLevels) + trueTend

                  if (iTracers .eq. nBenthicTracers .and. iLevels .eq. nBenthicVertLevels+1) then
                     call computeNetElements (netElements_end_reactions_tmp,nBenthicVertLevels,benthosTracerBulk,vertBenthosGridThickI)
                     rel_error(totc_ind) = netElements_end_reactions_tmp(totc_ind) - netElements_start_reactions_tmp(totc_ind)
                     rel_error_tot(totc_ind) = rel_error_tot(totc_ind) + rel_error(totc_ind)
                     rel_error(totc_ind) = rel_error(totc_ind)/(netElements_end_reactions_tmp(totc_ind) + c1_benthos)

                     do iElements = totp_ind,nElements
                        rel_error(iElements) = netElements_end_reactions_tmp(iElements) - netElements_start_reactions_tmp(iElements)
                        rel_error_tot(iElements) = rel_error_tot(iElements) + rel_error(iElements)
                        rel_error(iElements) = rel_error(iElements)/(netElements_end_reactions_tmp(iElements) + c1_benthos)
                     end do
                     if ( abs(maxval(rel_error)) .gt. 1.0e-11_benthos_r8 .and. &
                          .not. tooLarge .and.  (subN .lt. subN_o)) then
                        tooLarge = .true. ! subcycle
                        rel_error_tot(:) = c0_benthos
                        exit
                     end if
                  else
                     rel_error(:) = c0_benthos
                  end if

                  if (benthosTracerBulk(iTracers,iLevels) .lt. c0_benthos) then  !-puny) then

                     write(*,*)'Negative benthic tracers after  reactions in benthos_mod.F90'
                     write(*,*)'tempR,tempS, tempSol, benthosTracerBeforeReaction:',tempR,tempS, tempSol, benthosTracerBeforeReaction
                     write(*,*)'tracer index (iTracers,iLevels):',iTracers,iLevels
                     write(*,*)'benthosTracerBulk(iTracers,iLevels):',benthosTracerBulk(iTracers,iLevels)
                     err = 1

                  else if (benthosTracerBulk(iTracers,iLevels) .gt. 1.0e20_benthos_r8) then

                     write(*,*)'benthic tracers blowing up after  reactions'
                     write(*,*)'tempR,tempS, tempSol, benthosTracerBeforeReaction:',tempR,tempS, tempSol, benthosTracerBeforeReaction
                     write(*,*)'tracer index (iTracers,iLevels):',iTracers,iLevels
                     write(*,*)'benthosTracerBulk(iTracers,iLevels):',benthosTracerBulk(iTracers,iLevels)
                     err = 1

                  end if
               end do ! iTracers
            end do  ! iLevels
         end do !subtt  subgrid timestep
      end do !  while

   end if ! if useBenthicReactions

   call computeNetElements (netElements_end_reactions, nBenthicVertLevels, benthosTracerBulk,vertBenthosGridThickI)

   do iTracers = 1, nElements
      benthos_diagnostic_fields % relErrorElements(column,iTracers) = sqrt(netElementsRelError(iTracers) + &
           rel_error_tot(iTracers)**2)/(netElements_end_reactions(iTracers)+c1_benthos)
      benthos_diagnostic_fields % relErrorElementsReactions(column,iTracers) = rel_error_tot(iTracers)/&
           (netElements_end_reactions(iTracers)+c1_benthos)
      benthos_diagnostic_fields % diag_netElements(column,iTracers) = &
           netElements_end_reactions(iTracers)
   end do

   do iTracers = 1, nBenthicTracers

      do iLevels = 1,nBenthicVertLevels+1

      ! if (benthosTracerBulk(iTracers,iLevels) .lt. puny) then
      !     benthosTracerBulk(iTracers,iLevels) = c0_benthos
      ! end if

       !  benthos_output%benthosReactionTendencies(iLevels,column,iTracers) = &
       !       (benthosTracerBulk(iTracers,iLevels) - benthos_input%benthosTracerBulk(iLevels,column,iTracers))/dt
         benthos_output%benthosReactionTendencies(iLevels,column,iTracers) = &
              (benthosTracerBulk(iTracers,iLevels) - initBenthosTracerBulk(iTracers,iLevels))/dt

         benthos_output%benthosTendencies(iLevels,column,iTracers) = &
              benthos_output%benthosTendencies(iLevels,column,iTracers) + &
              benthos_output%benthosReactionTendencies(iLevels,column,iTracers)

      end do  ! iLevels

      call sum_benthos_column(benthos_diagnostic_fields % diag_netBenthosReactionTend(column,iTracers), &
         benthos_output % benthosReactionTendencies(:, column, iTracers),nBenthicVertLevels,vertBenthosGridThickI)

      benthos_diagnostic_fields % diag_netBenthosTend(column,iTracers) = &
      benthos_diagnostic_fields % diag_netBenthosTransportTend(column,iTracers) + &
      benthos_diagnostic_fields % diag_netBenthosReactionTend(column,iTracers)

      benthos_output%deepStorage(column,iTracers) = benthosStorageConc(iTracers)
      !do iLevels = 1, nBenthicVertLevels+1
      !   benthos_input%benthosTracerBulk(iLevels,column,iTracers) = benthosTracerBulk(iTracers,iLevels)
      !end do

   end do

   do iLevels = 1, nBenthicVertLevels+1
      benthos_output%PH_PREV_3D(iLevels,column) = PH_PREV_3D(iLevels)
   end do
   !
   ! call compute conservation and error diagnostics
   !
   deallocate(deepBenthosBottomFlux)
   deallocate(oceanSedFlux)
   deallocate(oceanPrecipFlux)
   deallocate(totalOceanSolidFlux)
   deallocate(oceanTracerConc)
   deallocate(netElementsStart)
   deallocate(netElementsRelError)
   deallocate(benthosTracerBulk)
   deallocate(initBenthosTracerBulk)
   deallocate(initBenthosTracerBulk_Carbonate)
   deallocate(totalColumnTracersInitial)
   deallocate(totalColumnTracersFinal)
   deallocate(bDiff)
   deallocate(iDiff)
   deallocate(bpor)
   deallocate(ipor)
   deallocate(C_top_vel)
   deallocate(C_bot_vel)
   deallocate(C_top)
   deallocate(C_bot)
   deallocate(source)
   deallocate(bottomsource)
   deallocate(C_tot)
!   deallocate(netTransport_tend)
   deallocate(oceanSedFluxdt)
   deallocate(oceanDiffVelFluxdt)
   deallocate(benthosSedFluxdt)
   deallocate(benthosDiffVelFluxdt)
   deallocate(diffError)
   deallocate(totFluxesdt)
   deallocate(primarySourceTend)
   deallocate(primarySinkTend)
   deallocate(carbonateSourceTend)
   deallocate(carbonateSinkTend)
   deallocate(secondarySourceTend)
   deallocate(secondarySinkTend)
   deallocate(netReactionTend)
   deallocate(rel_error)
   deallocate(rel_error_tot)
   deallocate(PH_PREV_3D)

 end subroutine benthos_SourceSink

!*****************************************************************************
!BOP
! !IROUTINE: compute_FCT_matrix_CN
! !INTERFACE:

 subroutine compute_FCT_matrix_CN (sbdiag, diag, spdiag, rhs, ML_diag, ML_sbdiag, ML_spdiag, &
      D_sbdiag, D_spdiag, &
      Qtop_out, Qbot_out, Sink_bot_diag,Sink_bot_sbdiag, Sink_top_diag,Sink_top_spdiag, &
      bioncons,  dt,  nBenthicVertLevels,   dhtop, &
      ipor, iDiff, benthosDepth,  totalOceanSolidFluxm, bpor,  &
      C_top, C_bot, C_top_vel, C_bot_vel,benthosBottomFluxm, vel, bDiff,zspace, dzspace)

! !DESCRIPTION:
!  compute the matrix elements for the flux-corrected conservative transport scheme
! Crank-nicholson version:
! [Mjk - dt * gamma *(Kjk + Sjk)]C'k = Mjk Ck + dt qk + (1-gamma)*(Kjk + Sjk)Ck
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels

   real(KIND=benthos_r8), dimension(:), intent(in) :: &
        bpor,  & ! grid fluid or solid porosity
        ipor,  & ! interface fluid or solid porosity
        iDiff, & !  (1/s) normalized  interface diffusivity
        bioncons, & ! mmol/m3 or mmol/kg tracer concentration
        bDiff,    & ! (1/s) normalized vertical grid diffusivity
        zspace,   & ! vertBenthosGridThickI (m) grid thickness of interface layers
        dzspace     ! vertBenthosGridThick  (m) grid thickness of layers

   real(KIND=benthos_r8), intent(in) :: &
        vel,                               & ! per s (normalized rate of surface growth or loss)
        benthosBottomFluxm, & ! (mmol/m3/s)  normalized positive down prescribed bottom flux
        C_bot_vel,                   & ! bottom concentration for velocity fluxes (mmol/m3 or mmol/kg)
        C_top_vel,                   & ! top concentration for velocity fluxes (mmol/m3 or mmol/kg)
        C_bot,                          & ! bottom concentration for diffusion(mmol/m3)
        C_top,                          & ! top concentration for diffusion (mmol/m3)
        totalOceanSolidFluxm, & ! precipitation and sinking flux of solids
        benthosDepth,              & ! (m) depth of the active benthos layer
        dhtop,                            & ! (m) change in the surface boundary in dt
        dt                                      ! (s) time-step

  ! !INPUT/OUTPUT PARAMETERS:

  ! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: &
       Sink_top_diag,   & ! (m/s * dt) diagonal matrix contribution to surface flux
       Sink_top_spdiag, & ! (m/s * dt) off-diagonal matrix contribution to surface flux
       Sink_bot_sbdiag, & ! (m/s * dt) off-diagonal matrix contribution to bottom flux
       Sink_bot_diag,   & ! (m/s * dt) diagonal matrix contribution to bottom flux
       Qbot_out,        & ! (mmol/m2/s * dt) rhs contribution to bottom flux
       Qtop_out           ! (mmol/m2/s * dt) ths contribution to surface flux

  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1), intent(out) :: &
       D_spdiag,             & !  matrix elements
       D_sbdiag,             & !
       ML_diag,              & !
       ML_sbdiag,            & !
       ML_spdiag,            & !
       rhs,                  & ! right hand side
       spdiag,               & !
       diag,                 &
       sbdiag

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
  !-----------------------------------------------------------------------------
  real(KIND=benthos_r8) :: &
       Sink_bot, &
       Sink_top, &
       dphi_dx, &
       dphi_dx1, &
       dphi_dxN, &
       gamma, &
       Sink_top_diag_tmp,   & ! (m/s * dt) diagonal matrix contribution to surface flux
       Sink_top_spdiag_tmp, & ! (m/s * dt) off-diagonal matrix contribution to surface flux
       Sink_bot_sbdiag_tmp, & ! (m/s * dt) off-diagonal matrix contribution to bottom flux
       Sink_bot_diag_tmp      ! (m/s * dt) diagonal matrix contribution to bottom flux

  real (KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       K_diag, &
       D_diag, &
       K_spdiag, &
       K_sbdiag, &
       S_diag, &
       S_spdiag, &
       S_sbdiag, &
       Q_bot, &
       Q_top, &
       diagTemp, &
       ML_diag_CN, &
       iDin, &
       topz, &
       interiorz

  integer (KIND=benthos_i4) :: &
       iLevels

!---------------------------------------------------------------------
!  Diag (jj) solve for j = 1:nBenthicVertLevels+1
!  spdiag(j) == (j,j+1) solve for j = 1:nBenthicVertLevels otherwise 0
!  sbdiag(j) == (j,j-1) solve for j = 2:nBenthicVertLevels+1 otherwise 0
!---------------------------------------------------------------------
     Qbot_out = c0_benthos
     Qtop_out = c0_benthos
     Sink_bot = c0_benthos
     Sink_top = c0_benthos

     K_spdiag(:) = c0_benthos
     K_sbdiag(:) = c0_benthos
     K_diag(:) = c0_benthos

     S_spdiag(:) = c0_benthos
     S_sbdiag(:) = c0_benthos
     S_diag(:) = c0_benthos

     D_spdiag(:) = c0_benthos
     D_sbdiag(:) = c0_benthos
     D_diag(:) = c0_benthos
     ML_diag(:) = c0_benthos
     spdiag(:) = c0_benthos
     sbdiag(:) = c0_benthos
     diag(:) = c0_benthos
     rhs(:) = c0_benthos
     topz(:) = c0_benthos
     interiorz(:) = c1_benthos

     ! define lumped mass matrix
     ML_spdiag(:) = c0_benthos
     ML_sbdiag(:) = c0_benthos

     ML_spdiag(1) = (zspace(1) + zspace(2))/12.0_benthos_r8
     ML_sbdiag(nBenthicVertLevels+1) = (zspace(nBenthicVertLevels)+zspace(nBenthicVertLevels+1))/&
          12.0_benthos_r8
     ML_diag(1) = zspace(1)
     ML_diag(nBenthicVertLevels+1) = zspace(nBenthicVertLevels+1)
     do iLevels = 2,nBenthicVertLevels
        ML_diag(iLevels) = zspace(iLevels)
        ML_spdiag(iLevels) = (zspace(iLevels) + zspace(iLevels+1))/12.0_benthos_r8
        ML_sbdiag(iLevels) = (zspace(iLevels-1) + zspace(iLevels))/12.0_benthos_r8
     end do

    ! compute matrix K: K_diag , K_sbdiag, K_spdiag
    ! compute matrix S: S_diag, S_sbdiag, S_spdiag

    topz(1) = c1_benthos
    interiorz(1) = c0_benthos

    do iLevels = 1 , nBenthicVertLevels+1
       iDin(iLevels) = iDiff(iLevels) !topz(iLevels) * bDiff(iLevels) + interiorz(iLevels) * iDiff(iLevels)
    end do

    !iLevels == 1
    iLevels = 1
    !dphi_dx = (ipor(iLevels+1)-bpor(iLevels))/&
    !     (c2_benthos*(benthosGridPointsI(iLevels+1) - benthosGridPointsI(iLevels)))
    dphi_dx = c0_benthos
    dphi_dx1 = dphi_dx

    K_sbdiag(iLevels)= c0_benthos
    K_spdiag(iLevels)= -p5_benthos*(min(c0_benthos,vel) +  &
       bDiff(iLevels)/(ipor(iLevels)*ipor(iLevels+1))*dphi_dx)
    K_diag(iLevels) =   -max(c0_benthos,vel)- (bDiff(iLevels))/&
       (dzspace(iLevels)*ipor(iLevels)) &
       + p5_benthos*(vel + bDiff(iLevels)/ipor(iLevels)**2*dphi_dx)

    S_sbdiag(iLevels) = c0_benthos
    S_diag(iLevels) = -bDiff(iLevels+1)/dzspace(iLevels)/ipor(iLevels)
    S_spdiag(iLevels) = bDiff(iLevels+1)/dzspace(iLevels)/ipor(iLevels+1)

    !(iLevels == nBenthicVertLevels+1
    iLevels = nBenthicVertLevels+1

    !dphi_dx = (bpor(iLevels+1)-ipor(iLevels-1))/abs(c1_benthos-benthosGridPointsI(iLevels-1))
    dphi_dx = c0_benthos
    dphi_dxN = dphi_dx

    K_spdiag(iLevels) = c0_benthos
    K_sbdiag(iLevels)= p5_benthos*(max(c0_benthos,vel) +  &
       iDin(iLevels)/(ipor(iLevels)*ipor(iLevels-1))*dphi_dx)
    K_diag(iLevels) =   - (bDiff(iLevels)) / (dzspace(iLevels-1)*ipor(iLevels)) - &
       p5_benthos*(vel + iDin(iLevels)/ipor(iLevels)**2*dphi_dx )+min(c0_benthos,vel)

    S_spdiag(iLevels) = c0_benthos
    S_sbdiag(iLevels) =  bDiff(iLevels)/ipor(iLevels-1)/dzspace(iLevels-1)
    S_diag(iLevels) = -bDiff(iLevels)/dzspace(iLevels-1)/ipor(iLevels)

    do iLevels = 2 , nBenthicVertLevels

       ! dphi_dx = (ipor(iLevels+1)-ipor(iLevels-1))/abs(benthosGridPointsI(iLevels+1)-benthosGridPointsI(iLevels-1))
       dphi_dx = c0_benthos

       K_sbdiag(iLevels)= p5_benthos*(max(c0_benthos,vel))  + &
          p5_benthos*( iDin(iLevels)/(ipor(iLevels)*ipor(iLevels-1))*dphi_dx)
       K_spdiag(iLevels)= -p5_benthos*(min(c0_benthos,vel)) - &
          p5_benthos*( iDin(iLevels)/(ipor(iLevels)*ipor(iLevels+1))*dphi_dx)
       K_diag(iLevels)  = -K_sbdiag(iLevels)-K_spdiag(iLevels)

       S_diag(iLevels)  = -(bDiff(iLevels)/dzspace(iLevels-1)/ipor(iLevels) + bDiff(iLevels+1)/dzspace(iLevels)/ipor(iLevels))
       S_sbdiag(iLevels)= bDiff(iLevels)/ipor(iLevels-1)/dzspace(iLevels-1)
       S_spdiag(iLevels) = bDiff(iLevels+1)/ipor(iLevels+1)/dzspace(iLevels)

    end do

    ! compute matrix artificial D: D_spdiag, D_diag  (D_spdiag(iLevels) = D_sbdiag(iLevels+1))

     do iLevels = 1,nBenthicVertLevels
        D_spdiag(iLevels)    = max(c0_benthos,max(-K_spdiag(iLevels), -K_sbdiag(iLevels+1)))
        D_sbdiag(iLevels+1)  = D_spdiag(iLevels)
      end do
      do  iLevels = 1,nBenthicVertLevels+1
         D_diag(iLevels) =  D_diag(iLevels) - D_spdiag(iLevels) - D_sbdiag(iLevels)
         Q_top(iLevels) = c0_benthos
         Q_bot(iLevels) = c0_benthos
      end do

      ! compute Q_top and Q_bot: top and bottom sources

      Q_top(1) =  (bDiff(1))*C_top/(c2_benthos*zspace(1)*bpor(1))+ &
           max(c0_benthos,totalOceanSolidFluxm)*p5_benthos

      Qtop_out = Q_top(1)

      Q_bot(nBenthicVertLevels+1) =  (iDin(nBenthicVertLevels+1))*C_bot &
           /(c2_benthos*zspace(nBenthicVertLevels+1)*bpor(nBenthicVertLevels+2)) - &
           min(c0_benthos,benthosBottomFluxm) - min(c0_benthos,vel)*C_bot*p5_benthos

       Qbot_out = Q_bot(nBenthicVertLevels+1)

       if (vel .LE. c0_benthos) then

           Sink_top_diag = K_diag(1)
           Sink_top_spdiag = -p5_benthos*(bDiff(1)/(ipor(1)*ipor(2))*dphi_dx1)

           Sink_bot_diag =   - (iDin(nBenthicVertLevels+1)) / &
                (c2_benthos*zspace(nBenthicVertLevels+1)*&
                ipor(nBenthicVertLevels+1)) - &
                p5_benthos*(iDin(nBenthicVertLevels+1)/&
                ipor(nBenthicVertLevels+1)**2*dphi_dxN )
           Sink_bot_sbdiag = K_sbdiag(nBenthicVertLevels+1)

        else   ! velocity advects downward away from the surface boundary

           Sink_top_diag = -(bDiff(1))/(c2_benthos*zspace(1)*ipor(1)) &
                + p5_benthos*(bDiff(1)/ipor(1)**2*dphi_dx1)
           Sink_top_spdiag = K_spdiag(1)

           Sink_bot_diag  = -iDin(nBenthicVertLevels+1)/&
                (c2_benthos*zspace(nBenthicVertLevels+1)*&
                ipor(nBenthicVertLevels+1))- &
                p5_benthos*(vel+iDin(nBenthicVertLevels+1)/ &
                ipor(nBenthicVertLevels+1)**2*dphi_dxN)

           Sink_bot_sbdiag =  p5_benthos*( iDin(nBenthicVertLevels+1)/&
                (ipor(nBenthicVertLevels+1)&
                *ipor(nBenthicVertLevels))*dphi_dxN )

        end if
        do iLevels = 1, nBenthicVertLevels+1
           spdiag(iLevels) = -dt *(D_spdiag(iLevels) + K_spdiag(iLevels) + S_spdiag(iLevels))
           sbdiag(iLevels) = -dt * (D_sbdiag(iLevels) + K_sbdiag(iLevels) + S_sbdiag(iLevels))
           diagTemp(iLevels) = -dt *(D_diag(iLevels) + K_diag(iLevels) + S_diag(iLevels))
           diag = ML_diag(iLevels) + diagTemp(iLevels)
           ML_diag_CN((iLevels)) = c0_benthos
        end do

        ML_diag_CN(1) = spdiag(1) * bioncons(2) +diagTemp(1)* bioncons(1)
        ML_diag_CN(nBenthicVertLevels+1) = sbdiag(nBenthicVertLevels+1) * &
             bioncons(nBenthicVertLevels)  + diagTemp(nBenthicVertLevels+1) * &
             bioncons(nBenthicVertLevels+1)

        do  iLevels = 2,nBenthicVertLevels
           ML_diag_CN(iLevels) = -spdiag(iLevels) * bioncons(iLevels+1) - &
                sbdiag(iLevels) * bioncons(iLevels-1) - diagTemp(iLevels)*bioncons(iLevels)
        end do

       !
       ! Crank-Nicholson Formulation
       !

        gamma = c1_benthos

        do iLevels = 1, nBenthicVertLevels+1
           spdiag(iLevels) = -dt *gamma*(D_spdiag(iLevels) + K_spdiag(iLevels) + S_spdiag(iLevels))
           sbdiag(iLevels) = -dt * gamma*(D_sbdiag(iLevels) + K_sbdiag(iLevels) + S_sbdiag(iLevels))
           diag(iLevels) = ML_diag(iLevels) - dt *gamma* (D_diag(iLevels) + K_diag(iLevels) + &
               S_diag(iLevels))

          rhs(iLevels) = ML_diag(iLevels) * bioncons(iLevels) + dt * Q_top(iLevels) + &
               dt* Q_bot(iLevels) + (c1_benthos-gamma)*ML_diag_CN(iLevels)*dt
       end do

     end subroutine compute_FCT_matrix_CN

!*****************************************************************************
!BOP
! !IROUTINE: tridiag_solverz
! !INTERFACE:

subroutine tridiag_solverz(xout, nmat, sbdiagz, diagz, spdiagz, rhsz)

! !DESCRIPTION:
!  tridiagonal matrix solver
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nmat

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       diagz,     &
       spdiagz, &
       sbdiagz, &
       rhsz

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8),dimension(:), intent(inout) :: xout

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  integer(KIND=benthos_i4) :: k

  real(KIND=benthos_r8) :: wbeta

  real(KIND=benthos_r8), dimension(nmat) :: &
       wgamma

   wbeta = diagz(1)
   xout(1) = rhsz(1) / wbeta

   do k = 2, nmat
      wgamma(k) = spdiagz(k-1) / wbeta
      wbeta = diagz(k) - sbdiagz(k)*wgamma(k)
      xout(k) = (rhsz(k) - sbdiagz(k)*xout(k-1)) / wbeta
   end do

    do k = nmat-1, 1, -1
         xout(k) = xout(k) - wgamma(k+1)*xout(k+1)
    end do

  end subroutine tridiag_solverz

!*****************************************************************************
!BOP
! !IROUTINE: check_conservation_FCT
! !INTERFACE:

  subroutine check_conservation_FCT(fluxbio, l_stop,deep_storage,oceanSedFluxdt, &
       oceanDiffVelFluxdt, benthosSedFluxdt, benthosDiffVelFluxdt, diffError, sources, &
       C_init, C_new, S_top,  S_bot, L_bot_d, L_bot_s,L_top_d,L_top_s, dt,     &
       nBenthicVertLevels, bottomsource, vel,benthosDepth, iTracers, &
       benthosBottomFluxm, totalOceanSolidFluxm,vertBenthosGridThickI, iDiff, oceanBottomTemperature)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nBenthicVertLevels, &
       iTracers

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       C_init, &
       C_new,  &
       vertBenthosGridThickI, &
       iDiff

  real(KIND=benthos_r8), intent(in) :: &
       S_top, &
       S_bot, &
       L_bot_d, &
       L_bot_s, &
       L_top_s, &
       L_top_d, &
       dt, &
       bottomsource, &
       vel, &
       benthosDepth, &
       benthosBottomFluxm, &
       totalOceanSolidFluxm, &
       oceanBottomTemperature

  ! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(inout) :: &
       deep_storage, & !
       oceanSedFluxdt, & !
       benthosSedFluxdt, &
       sources

  ! !OUTPUT PARAMETERS:

  logical(KIND=benthos_log), intent(out) :: l_stop

  real(KIND=benthos_r8), intent(out) :: &
       fluxbio, & !
       oceanDiffVelFluxdt, & !
       benthosDiffVelFluxdt, &
       diffError

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       dC

  real(KIND=benthos_r8) :: &
       C_init_tot, &
       C_new_tot, &
       C_max, &
       dC_tot, &
       sources_o, &
       sources_top, &
       sources_bot, &
       deep_storage_tmp, &
       diffCnewCinit, &
       diffError_o, &
       sources_top_corr, &
       sources_bot_corr, &
       sources_corr, &
       diffError_corr, &
       accuracy

  integer(KIND=benthos_i4) :: &
       iLevels

  l_stop = .false.

!-------------------------------------
!  Ocean flux: positive into the ocean
!-------------------------------------

   do iLevels = 1, nBenthicVertLevels + 1
      dC(iLevels) = C_new(iLevels) - C_init(iLevels)
   end do

   call sum_benthos_column (C_init_tot, C_init, nBenthicVertLevels,vertBenthosGridThickI)
   call sum_benthos_column (C_new_tot, C_new,nBenthicVertLevels,vertBenthosGridThickI)
   call sum_benthos_column (dC_tot, dC,nBenthicVertLevels,vertBenthosGridThickI)

    ! This is the estimate of the source and source from the top and bottom
    !       sources = (S_top+S_bot+L_bot_d*C_new(nBenthicVertLevels+1)+
    !  L_top_d*C_new(1)+L_bot_s*C_new(nBenthicVertLevels)+L_top_s*C_new(2))*dt*benthosDepth

   sources_o = (S_top+S_bot+L_bot_d*C_init(nBenthicVertLevels+1)+L_top_d*C_init(1)+ &
        L_bot_s*C_init(nBenthicVertLevels)+L_top_s*C_init(2))*dt*benthosDepth

   sources_top = (S_top+L_top_d*C_new(1)+L_top_s*C_new(2))*dt*benthosDepth

   sources_bot = (S_bot+L_bot_d*C_new(nBenthicVertLevels+1)+ &
        L_bot_s*C_new(nBenthicVertLevels))*dt*benthosDepth   ! mmol/kg * m or mmol/m2
              ! positive implies positive increase in benthic  so flux into the Benthic from below

   sources = sources_top + sources_bot

   deep_storage_tmp = deep_storage-sources_bot
   oceanSedFluxdt = totalOceanSolidFluxm*dt*benthosDepth
   benthosSedFluxdt = benthosBottomFluxm*dt*benthosDepth
   oceanDiffVelFluxdt = sources_top - oceanSedFluxdt
   benthosDiffVelFluxdt = -(sources_bot + benthosSedFluxdt)

   diffCnewCinit = (C_new_tot-C_init_tot)      ! mmol/m2
   diffError = (diffCnewCinit - sources)
   diffError_o = (diffCnewCinit-sources_o)    ! estimate of sources from old C_init

   benthosDiffVelFluxdt = benthosDiffVelFluxdt - fluxCorrect*diffError
   deep_storage = deep_storage_tmp - fluxCorrect*diffError

  ! if (useFluxCorrection) then
  !    benthosDiffVelFluxdt = benthosDiffVelFluxdt - diffError
  !    deep_storage = deep_storage_tmp - diffError
  ! else
  !    deep_storage = deep_storage_tmp
  ! end if

   sources_top_corr = oceanDiffVelFluxdt + oceanSedFluxdt
   sources_bot_corr = -benthosDiffVelFluxdt - benthosSedFluxdt
   sources_corr = sources_top_corr + sources_bot_corr
   diffError_corr = (diffCnewCinit - sources_corr)

   fluxbio = oceanDiffVelFluxdt/dt    !passed to the ocean bottom layers

! implicit time integration,  errors are linear in dt  and decrease with the grid spacing zspace*h.
   C_max = max(maxval(C_new),maxval(C_init))
   accuracy = puny*C_max

   if (minval(C_new) < c0_benthos .or.abs(diffError_corr) .gt. accuracy) then
      write(*,*)'Positivity of zbgc low order solution failed: C_new:',C_new
      write(*,*)'OR Conservation of zbgc low order solution failed: diffError (mmol/m2):', diffError
      write(*,*)'C_init:',C_init
      write(*,*)'S_top, S_bot, L_bot_d, L_bot_s, L_top_s, L_top_d:'
      write(*,*) S_top, S_bot, L_bot_d, L_bot_s, L_top_s, L_top_d
      write(*,*) 'oceanBottomTemperature:', oceanBottomTemperature
      write(*,*) 'vertBenthosGridThickI = zspace : ', vertBenthosGridThickI
      write(*,*) 'iDiff : ', iDiff
      write(*,*)'bottomsource, vel, benthosDepth, benthosBottomFluxm, totalOceanSolidFluxm'
      write(*,*) bottomsource, vel, benthosDepth, benthosBottomFluxm, totalOceanSolidFluxm
      write(*,*)'dt,sources_top_corr, sources_bot_corr, sources_corr, diffError_corr:'
      write(*,*) dt,sources_top_corr, sources_bot_corr, sources_corr, diffError_corr
      write(*,*)'sources', sources
      write(*,*)'sources_o',sources_o
      write(*,*)'iTracers,fluxbio:', iTracers,fluxbio
      write(*,*)'diffError_o', diffError_o
      write(*,*)'diffError_corr', diffError_corr
      write(*,*)'accuracy:',accuracy
      write(*,*)'Start deep_storage',deep_storage
      write(*,*)'intermediate deep_storage',deep_storage_tmp
      write(*,*)'Final deep_storage',deep_storage
      write(*,*)'diffCnewCinit (mmol/m2):',diffCnewCinit
      write(*,*)'oceanSedFluxdt',oceanSedFluxdt
      write(*,*)'oceanDiffVelFluxdt',oceanDiffVelFluxdt
      write(*,*)'fluxbio:',fluxbio
      write(*,*)'bottom final tracer', C_new(nBenthicVertLevels+1)
      write(*,*)'top final tracer', C_new(1)
      write(*,*)'bottom init tracer', C_init(nBenthicVertLevels+1)
      write(*,*)'top init tracer', C_init(1)
      write(*,*)'vel:',vel
      write(*,*)'Total initial tracer', C_init_tot
      write(*,*)'Total final1  tracer', C_new_tot
      l_stop = .true.
   end if

   diffError = diffError_corr !/(max(C_new_tot, C_init_tot)+c1_benthos)

   end subroutine check_conservation_FCT

!*****************************************************************************
!BOP
! !IROUTINE: computeElementFluxes
! !INTERFACE:

   subroutine computeElementFluxes (oceanElementExchange, burialElementExchange, oceanSedFluxdt, &
        oceanDiffVelFluxdt, benthosSedFluxdt, benthosDiffVelFluxdt,benthosDepth)

! !DESCRIPTION:
!  calculates net boundary exchange fluxes for each of 7 elements: C, O, N, P, S, Mn, Fe
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       benthosDepth

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosDiffVelFluxdt, &
       benthosSedFluxdt, &
       oceanDiffVelFluxdt, &
       oceanSedFluxdt

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:
  real(KIND=benthos_r8), dimension(nElements), intent(out) :: &
       burialElementExchange, &
       oceanElementExchange

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
   real (KIND=benthos_r8), dimension(nBenthicTracers) :: &
       topExchange, & !
       botExchange

   integer (KIND=benthos_i4) :: &
        iTracers, &
        iElements

   do iTracers = 1,nBenthicTracers
      topExchange(iTracers) = oceanDiffVelFluxdt(iTracers) + oceanSedFluxdt(iTracers)      ! positive down
      botExchange(iTracers) = benthosSedFluxdt(iTracers) + benthosDiffVelFluxdt(iTracers)   ! positive down
   end do

   do iElements = 1,nElements
      burialElementExchange(iElements) = c0_benthos
      oceanElementExchange(iElements) = c0_benthos
      do iTracers = 1,nBenthicTracers
         oceanElementExchange(iElements) = oceanElementExchange(iElements) + &
            elementRatios(iTracers,iElements)*topExchange(iTracers)
         burialElementExchange(iElements) = burialElementExchange(iElements) + &
            elementRatios(iTracers,iElements)*botExchange(iTracers)
      end do
   end do

   end subroutine computeElementFluxes

!*****************************************************************************
!BOP
! !IROUTINE: netElementExchange
! !INTERFACE:

   subroutine netElementExchange (oceanElementSedimentation, deepBurial, &
        totalOceanSolidFlux,benthosStorageConc,dt,benthosDepth)

! !DESCRIPTION:
!  calculates net boundary exchange fluxes for each of 7 elements: C, O, N, P, S, Mn, Fe
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       benthosDepth, &
       dt

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosStorageConc, &
       totalOceanSolidFlux

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:
  real(KIND=benthos_r8), dimension(nElements), intent(out) :: &
       oceanElementSedimentation, &
       deepBurial

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

   integer (KIND=benthos_i4) :: &
        iTracers, &
        iElements

   do iElements = 1,nElements
      do iTracers = 1,nBenthicTracers
         oceanElementSedimentation(iElements) = oceanElementSedimentation(iElements) + &
              elementRatios(iTracers,iElements)*totalOceanSolidFlux(iTracers)*benthosDepth*dt
         deepBurial(iElements) = deepBurial(iElements) + elementRatios(iTracers,iElements)*benthosStorageConc(iTracers)
      end do
   end do

   end subroutine netElementExchange

!*****************************************************************************
!BOP
! !IROUTINE: compute_FCT_corr
! !INTERFACE:

subroutine compute_FCT_corr (C_in,  &
                             C_low, dt, nBenthicVertLevels, &
                             D_sbdiag, D_spdiag, ML, ML_spdiag, &
                             ML_sbdiag, zspace)

! !DESCRIPTION:
!  computes the correction to the initial solver
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels

  real(KIND=benthos_r8), intent(in) ::  dt

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       C_in, &
       D_sbdiag, &
       D_spdiag, &
       ML, &
       ML_spdiag, &
       ML_sbdiag, &
       zspace

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(inout) :: &
       C_low

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       F_diag, &
       F_spdiag, &
       F_sbdiag, &
       a_spdiag, &
       a_sbdiag, &
       Pplus, &
       Pminus, &
       Qplus, &
       Qminus, &
       Rplus, &
       Rminus

  integer(KIND=benthos_i4) :: &
       iLevels

  F_diag(:) = c0_benthos
  F_spdiag(:) = c0_benthos
  F_sbdiag(:) = c0_benthos
  a_spdiag(:) = c0_benthos
  a_sbdiag(:) = c0_benthos
! C_low_new(nBenthicVertLevels+1) = C_low(nBenthicVertLevels+1)

  do iLevels = 1,nBenthicVertLevels

 !    C_low_new(iLevels) = C_low(iLevels)
     F_spdiag(iLevels) = ML_spdiag(iLevels)*(C_low(iLevels)-C_in(iLevels) - C_low(iLevels+1)+ &
          C_in(iLevels+1))/dt + D_spdiag(iLevels)*(C_low(iLevels)-C_low(iLevels+1))
     F_sbdiag(iLevels+1) =  ML_sbdiag(iLevels+1)*(C_low(iLevels+1)-C_in(iLevels+1) - &
          C_low(iLevels)+ &
          C_in(iLevels))/dt + D_sbdiag(iLevels+1)*(C_low(iLevels+1)-C_low(iLevels))

     if (F_spdiag(iLevels)*(C_low(iLevels) - C_low(iLevels+1)) > c0_benthos) &
          F_spdiag(iLevels) = c0_benthos

     if (F_sbdiag(iLevels+1)*(C_low(iLevels+1) - C_low(iLevels)) > c0_benthos) &
          F_sbdiag(iLevels+1) = c0_benthos

  end do

  if (maxval(abs(F_spdiag)) > c0_benthos) then

     ! compute the weighting factors: a_spdiag, a_sbdiag

     Pplus(1)  = max(c0_benthos, F_spdiag(1))
     Pminus(1) = min(c0_benthos, F_spdiag(1))
     Pplus(nBenthicVertLevels+1)  = max(c0_benthos, F_sbdiag(nBenthicVertLevels+1))
     Pminus(nBenthicVertLevels+1) = min(c0_benthos, F_sbdiag(nBenthicVertLevels+1))
     Qplus(1) = max(c0_benthos,C_low(2)-C_low(1))
     Qminus(1)= min(c0_benthos,C_low(2)-C_low(1))
     Qplus(nBenthicVertLevels+1) = max(c0_benthos,C_low(nBenthicVertLevels)- &
          C_low(nBenthicVertLevels+1))
     Qminus(nBenthicVertLevels+1)= min(c0_benthos,C_low(nBenthicVertLevels)- &
          C_low(nBenthicVertLevels+1))
     Rplus(1)  = min(c1_benthos, ML(1)*Qplus(1)/dt/(Pplus(1)+puny16))
     Rminus(1) = min(c1_benthos, ML(1)*Qminus(1)/dt/(Pminus(1)-puny16))
     Rplus(nBenthicVertLevels+1)  = min(c1_benthos, ML(nBenthicVertLevels+1)* &
          Qplus(nBenthicVertLevels+1)&
          /dt/(Pplus(nBenthicVertLevels+1)+puny16))
     Rminus(nBenthicVertLevels+1) = min(c1_benthos, ML(nBenthicVertLevels+1)* &
          Qminus(nBenthicVertLevels+1) &
          /dt/(Pminus(nBenthicVertLevels+1)-puny16))

     do iLevels = 2,nBenthicVertLevels
        Pplus(iLevels)  = max(c0_benthos,F_spdiag(iLevels)) + max(c0_benthos,F_sbdiag(iLevels))
        Pminus(iLevels) = min(c0_benthos,F_spdiag(iLevels)) + min(c0_benthos,F_sbdiag(iLevels))

        Qplus(iLevels)  = max(c0_benthos, max(C_low(iLevels+1)-C_low(iLevels), &
             C_low(iLevels-1)-C_low(iLevels)))
        Qminus(iLevels) = min(c0_benthos, min(C_low(iLevels+1)-C_low(iLevels), &
             C_low(iLevels-1)-C_low(iLevels)))

        Rplus(iLevels)  = min(c1_benthos, ML(iLevels)*Qplus(iLevels)/dt/(Pplus(iLevels)+puny16))
        Rminus(iLevels) = min(c1_benthos, ML(iLevels)*Qminus(iLevels)/dt/(Pminus(iLevels)-puny16))

     end do

     do iLevels = 1,nBenthicVertLevels

        a_spdiag(iLevels) = min(Rminus(iLevels),Rplus(iLevels+1))
        if (F_spdiag(iLevels) > c0_benthos) &
             a_spdiag(iLevels) = min(Rplus(iLevels),Rminus(iLevels+1))

        a_sbdiag(iLevels+1) = min(Rminus(iLevels+1),Rplus(iLevels))

        if (F_sbdiag(iLevels+1) > c0_benthos) &
             a_sbdiag(iLevels+1) = min(Rplus(iLevels+1),Rminus(iLevels))

     end do

     ! compute F_diag

     F_diag(1) = a_spdiag(1)*F_spdiag(1)
     F_diag(nBenthicVertLevels+1) = a_sbdiag(nBenthicVertLevels+1)* &
          F_sbdiag(nBenthicVertLevels+1)

     C_low(1) = C_low(1) + dt*F_diag(1)/ML(1)
     C_low(nBenthicVertLevels+1) = C_low(nBenthicVertLevels+1) + &
          dt*F_diag(nBenthicVertLevels+1)/ &
          ML(nBenthicVertLevels+1)

     do iLevels = 2,nBenthicVertLevels

        F_diag(iLevels) = a_spdiag(iLevels)*F_spdiag(iLevels) + &
             a_sbdiag(iLevels)*F_sbdiag(iLevels)
        C_low(iLevels) = C_low(iLevels) + dt*F_diag(iLevels)/ML(iLevels)

     end do

  end if  !F_spdiag is nonzero

end subroutine compute_FCT_corr

!*****************************************************************************
!BOP
! !IROUTINE: primaryStoichMatrix
! !INTERFACE:

subroutine primaryStoichMatrix (primarySourceStoich, primarySinkStoich, &
                    benthosTracerBulkLevel)

! !DESCRIPTION:
!
!  computes the stoichiometry matrix for the primary reactions
!
!  primary reactions OM has the form (CH2O)a(NH3)b(H3PO4)c or POC_a PON_b POP_c
!
! P1.  OM/a+ O2 + (b-c)/aH^+ --> CO2 + b/aNH4 + c/aH2PO4 + H2O
! P2.  OM/a + 4/5 NO3 + (4/5 + (b-c)/a)H^+ ---> CO2 + 2/5 N2 + b/aNH4 + c/aH2PO4 + 7/5H2O
! P3.  OM/a + 2MnO2 + 4H + (b-c)/aH --> CO2 + b/aNH4 + c/aH2PO4 + 2Mn + 3H2O
! P4.  OM/a + 4Fe(OH)3 + 4(gamma)[Fe-P] + (12 + (b-c)/a-4gamma)H^+ --> CO2 + b/aNH4 +  (c/a+4gamma) H2PO4 + 4Fe + 13H2O
! P5.  OM/a +  1/2 SO4 + (1+(b-c)/a)H^+ --> CO2 + b/aNH4 + c/aH2PO4 + 1/2H2S + H2O
! P6.  OM/a + (b-c)/aH^+ --> 1/2CO2 + 1/2 CH4 + b/aNH4 + c/aH2PO4 + H2O
!  gamma is the fraction of iron-bound phosphorus
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulkLevel

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:,:), intent(out) :: &
       primarySourceStoich, &
       primarySinkStoich

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  integer(KIND=benthos_i4) :: iTracers, iReactions

  real(KIND=benthos_r8), dimension(3) :: &
       a_ratio, &
       b_ratio

  real(KIND=benthos_r8) :: &
       c_div_a, &
       b_div_a

   a_ratio(:) = CtoP
   b_ratio(:) = NtoP

   if (benthosTracerBulkLevel(popa_ind) .gt. puny .and. benthosTracerBulkLevel(poca_ind) .gt. puny .and. benthosTracerBulkLevel(pona_ind) .gt. puny) then
      a_ratio(1) = benthosTracerBulkLevel(poca_ind)/benthosTracerBulkLevel(popa_ind)
      b_ratio(1) = benthosTracerBulkLevel(pona_ind)/benthosTracerBulkLevel(popa_ind)
   end if

   if (benthosTracerBulkLevel(popb_ind) .gt. puny .and. benthosTracerBulkLevel(pocb_ind) .gt. puny .and. benthosTracerBulkLevel(ponb_ind) .gt. puny) then
      a_ratio(2) = benthosTracerBulkLevel(pocb_ind)/benthosTracerBulkLevel(popb_ind)
      b_ratio(2) = benthosTracerBulkLevel(ponb_ind)/benthosTracerBulkLevel(popb_ind)
   end if

   if (benthosTracerBulkLevel(popc_ind) .gt. puny .and. benthosTracerBulkLevel(pocc_ind) .gt. puny .and. benthosTracerBulkLevel(ponc_ind) .gt. puny) then
      a_ratio(3) = benthosTracerBulkLevel(pocc_ind)/benthosTracerBulkLevel(popc_ind)
      b_ratio(3) = benthosTracerBulkLevel(ponc_ind)/benthosTracerBulkLevel(popc_ind)
   end if

   c_div_a = (c1_benthos/a_ratio(1) + c1_benthos/a_ratio(2) + c1_benthos/a_ratio(3))/3.0_benthos_r8
   b_div_a = (b_ratio(1)/a_ratio(1) + b_ratio(2)/a_ratio(2) + b_ratio(3)/a_ratio(3))/3.0_benthos_r8

   do iReactions = 1,nPrimaryReactions
      do iTracers = 1,nBenthicTracers
         primarySourceStoich(iTracers,iReactions) = c0_benthos
         primarySinkStoich(iTracers,iReactions) = c0_benthos
      end do

   !define ratios based on concentrations
   !
   ! if (benthosTracerBulkLevel(popa_ind) .gt.  puny)
   !    a_ratio = benthosTracerBulkLevel(poca_ind)/benthosTracerBulkLevel(popa_ind)
   !    b_ratio = benthosTracerBulkLevel(pona_ind)/benthosTracerBulkLevel(popa_ind)
   ! end
      ! Rates are in mmol/m3/s.  Add conversion for solids to mmol/kg/s
      ! epsilon_Phosphate is accelerated POP remineralization (Reed et al 2015)
      ! epsilon_Phosphate = 2 works okay for Arkona Basin but may not work in other places

      primarySinkStoich(poca_ind,iReactions) = c1_benthos/sediment_density
      primarySinkStoich(pocb_ind,iReactions) = c1_benthos/sediment_density
      primarySinkStoich(pocc_ind,iReactions) = c1_benthos/sediment_density
      primarySinkStoich(pona_ind,iReactions) = b_ratio(1)/a_ratio(1)/sediment_density
      primarySinkStoich(ponb_ind,iReactions) = b_ratio(2)/a_ratio(2)/sediment_density
      primarySinkStoich(ponc_ind,iReactions) = b_ratio(3)/a_ratio(3)/sediment_density
      primarySinkStoich(popa_ind,iReactions) = c1_benthos/a_ratio(1)/sediment_density*epsilon_Phosphate
      primarySinkStoich(popb_ind,iReactions) = c1_benthos/a_ratio(2)/sediment_density*epsilon_Phosphate
      primarySinkStoich(popc_ind,iReactions) = c1_benthos/a_ratio(3)/sediment_density*epsilon_Phosphate
      primarySinkStoich(alk_ind,iReactions) = c_div_a

      primarySourceStoich(dic_ind,iReactions) = c1_benthos
      primarySourceStoich(co2_ind,iReactions) = c1_benthos
      primarySourceStoich(nh4_ind,iReactions) = b_div_a
      primarySourceStoich(h2po4_ind,iReactions) = c_div_a*epsilon_Phosphate

   end do

   primarySinkStoich(o2_ind,1) = c1_benthos
   primarySinkStoich(alk_ind,4) = (c_div_a*epsilon_Phosphate+4.0*ironBoundPFraction)
   primarySinkStoich(alk_ind,5) = c_div_a*epsilon_Phosphate
   primarySinkStoich(alk_ind,6) = c_div_a*epsilon_Phosphate
   primarySinkStoich(no3_ind,2) = 4.0_benthos_r8/5.0_benthos_r8
   primarySinkStoich(mno2a_ind,3) = c2_benthos/sediment_density
   primarySinkStoich(feoh3a_ind,4) = 4.0_benthos_r8/sediment_density
   primarySinkStoich(fepa_ind,4) = 4.0_benthos_r8*ironBoundPFraction/sediment_density
   primarySinkStoich(so4_ind,5) = p5_benthos

   primarySourceStoich(dic_ind,6) = p5_benthos
   primarySourceStoich(co2_ind,6) = p5_benthos
   primarySourceStoich(alk_ind,1) = b_div_a
   primarySourceStoich(alk_ind,2) = (4.0_benthos_r8/5.0_benthos_r8 + b_div_a)
   primarySourceStoich(alk_ind,3) = (4.0_benthos_r8+b_div_a)
   primarySourceStoich(alk_ind,4) = (12.0_benthos_r8+b_div_a)
   primarySourceStoich(alk_ind,5) = (c1_benthos+b_div_a)
   primarySourceStoich(alk_ind,6) = b_div_a
   primarySourceStoich(h2po4_ind,4) =  (c_div_a*epsilon_Phosphate + 4.0_benthos_r8*ironBoundPFraction)
   primarySourceStoich(mn_ind,3) = c2_benthos
   primarySourceStoich(fe_ind,4) = 4.0_benthos_r8
   primarySourceStoich(h2s_ind,5) = p5_benthos
   primarySourceStoich(ch4_ind,6) = p5_benthos

   end subroutine primaryStoichMatrix

!*****************************************************************************
!BOP
! !IROUTINE: primaryRateConstants
! !INTERFACE:

   subroutine primaryRateConstants (primaryRates,benthosTracerBulkLevel,nOMtype)

! !DESCRIPTION:
! compute the rates for the primary redox reactions
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nOMtype

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulkLevel

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nOMtype,nPrimaryReactions), intent(out) :: &
       primaryRates   ! mmol/m3/s

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables

  real(KIND=benthos_r8), dimension(nOMtype) :: &
       k_OM

  real(KIND=benthos_r8), dimension(nOMtype) :: &
       OM

  integer (KIND=benthos_i4) :: &
       iType, iReactions

   k_OM = (/k_alpha, k_beta, k_ref/)
   OM = (/ benthosTracerBulkLevel(poca_ind), benthosTracerBulkLevel(pocb_ind), &
        benthosTracerBulkLevel(pocc_ind) /)*sediment_density

   do iType = 1,nOMtype

      do iReactions = 1,nPrimaryReactions
            primaryRates(iType,iReactions) = c0_benthos
      end do

      primaryRates(iType,1) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(o2_ind)/(k_o2 +  benthosTracerBulkLevel(o2_ind)))

      primaryRates(iType,2) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(no3_ind)/(k_no3 +  benthosTracerBulkLevel(no3_ind)))* &
           (k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))

      primaryRates(iType,3) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(mno2a_ind)/(k_mno2a + &
           benthosTracerBulkLevel(mno2a_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))*&
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))

      primaryRates(iType,4) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(feoh3a_ind)/(k_feoh3a + &
           benthosTracerBulkLevel(feoh3a_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))*&
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))*&
           (k_mno2a/(benthosTracerBulkLevel(mno2a_ind)+k_mno2a))

      primaryRates(iType,5) = &
           reductionFractionSO4*k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(so4_ind)/(k_so4 + &
           benthosTracerBulkLevel(so4_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))* &
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))*&
           (k_mno2a/(benthosTracerBulkLevel(mno2a_ind)+k_mno2a))* &
           (k_feoh3a/(benthosTracerBulkLevel(feoh3a_ind)+k_feoh3a))

      primaryRates(iType,6) = &
           reductionFractionSO4*k_OM(iType)*OM(iType)*(k_so4/(k_so4 + &
           benthosTracerBulkLevel(so4_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))*&
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))*&
           (k_mno2a/(benthosTracerBulkLevel(mno2a_ind)+k_mno2a))* &
           (k_feoh3a/(benthosTracerBulkLevel(feoh3a_ind)+k_feoh3a))

   end do

 end subroutine primaryRateConstants

!*****************************************************************************
!BOP
! !IROUTINE: primaryBenthosReactions
! !INTERFACE:

 subroutine primaryBenthosReactions (primarySourceTendk,primarySinkTendk,&
                    benthosTracerBulkLevels,primarySourceStoich,&
                    primarySinkStoich,primaryRates,dt,nOMtype)

! !DESCRIPTION:
!  compute tendencies for the primary redox reactions
!
!  primary reactions OM has the form (CH2O)a(NH3)b(H3PO4)c or POC_a PON_b POP_c
!
! P1.  OM/a+ O2 + (b-c)/aH^+ --> CO2 + b/aNH4 + c/aH2PO4 + H2O
! P2.  OM/a + 4/5 NO3 + (4/5 + (b-c)/a)H^+ ---> CO2 + 2/5 N2 + b/aNH4 + c/aH2PO4 + 7/5H2O
! P3.  OM/a + 2MnO2 + 4H + (b-c)/aH --> CO2 + b/aNH4 + c/aH2PO4 + 2Mn + 3H2O
! P4.  OM/a + 4Fe(OH)3 + 5(gamma)[Fe-P] + (12 + (b-c)/a-4gamma)H^+ --> CO2 + b/aNH4 +  (c/a+4gamma) H2PO4 + 4Fe + 13H2O
! P5.  OM/a +  1/2 SO4 + (1+(b-c)/a)H^+ --> CO2 + b/aNH4 + c/aH2PO4 + 1/2H2S + H2O
! P6.  OM/a + (b-c)/aH^+ --> 1/2CO2 + 1/2 CH4 + b/aNH4 + c/aH2PO4 + H2O
!  gamma is the fraction of iron-bound phosphorus

!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nOMtype

  real(KIND=benthos_r8), intent(in) :: &
       dt

  real(KIND=benthos_r8),dimension(:,:), intent(in) :: &
       primarySourceStoich, &
       primarySinkStoich, &
       primaryRates

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulklevels

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:),  intent(inout) :: &
       primarySourceTendk, &
       primarySinkTendk

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
   integer (KIND=benthos_i4) :: &
        iReactions, &
        iType, &
        iTracers

   do iReactions  = 1,nPrimaryReactions

      primarySourceTendk(poca_ind) = &
           primarySourceTendk(poca_ind) + primarySourceStoich(poca_ind,iReactions)*&
           primaryRates(1,iReactions)*dt

      primarySourceTendk(pocb_ind) = &
           primarySourceTendk(pocb_ind) + primarySourceStoich(pocb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySourceTendk(pocc_ind) = &
           primarySourceTendk(pocc_ind) + primarySourceStoich(pocc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySourceTendk(pona_ind) = &
           primarySourceTendk(pona_ind) + primarySourceStoich(pona_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySourceTendk(ponb_ind) = &
           primarySourceTendk(ponb_ind) + primarySourceStoich(ponb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySourceTendk(ponc_ind) = &
           primarySourceTendk(ponc_ind) + primarySourceStoich(ponc_ind,iReactions)*&
           primaryRates(3,iReactions)*dt

      primarySourceTendk(popa_ind) = &
           primarySourceTendk(popa_ind) + primarySourceStoich(popa_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySourceTendk(popb_ind) = &
           primarySourceTendk(popb_ind) + primarySourceStoich(popb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySourceTendk(popc_ind) = &
           primarySourceTendk(popc_ind) + primarySourceStoich(popc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySinkTendk(poca_ind) = &
           primarySinkTendk(poca_ind) + primarySinkStoich(poca_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySinkTendk(pocb_ind) = &
           primarySinkTendk(pocb_ind) +    primarySinkStoich(pocb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySinkTendk(pocc_ind) = &
           primarySinkTendk(pocc_ind) +    primarySinkStoich(pocc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySinkTendk(pona_ind) = &
           primarySinkTendk(pona_ind) +    primarySinkStoich(pona_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySinkTendk(ponb_ind) = &
           primarySinkTendk(ponb_ind) +    primarySinkStoich(ponb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySinkTendk(ponc_ind) = &
           primarySinkTendk(ponc_ind) +    primarySinkStoich(ponc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySinkTendk(popa_ind) = &
           primarySinkTendk(popa_ind) +    primarySinkStoich(popa_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySinkTendk(popb_ind) = &
           primarySinkTendk(popb_ind) +    primarySinkStoich(popb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySinkTendk(popc_ind) = &
           primarySinkTendk(popc_ind) +    primarySinkStoich(popc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt


      do iTracers = o2_ind,nBenthicTracers
         do iType = 1,3
            primarySourceTendk(iTracers) = primarySourceTendk(iTracers) + &
                 primarySourceStoich(iTracers,iReactions)*primaryRates(iType,iReactions)*dt

            primarySinkTendk(iTracers) = primarySinkTendk(iTracers) + &
                 primarySinkStoich(iTracers,iReactions)*primaryRates(iType,iReactions)*dt
         end do ! iType
      end do ! iTracers
   end do ! iReactions

 end subroutine primaryBenthosReactions

!*****************************************************************************
!BOP
! !IROUTINE: comp_CO3terms
! !INTERFACE:

subroutine comp_CO3terms_benthos (CO3,HCO3,H2CO3,ph,kw,kb, &
                       ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,dic,ta,pt,sit, temp, salt,depth, &
                       lcomp_co3_coeffs, sit_in, phlo, phhi, dic_in,ta_in,pt_in,k,&
                       oceanBottomDensity,benthosPorosityLevel)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
!   after SUBROUTINE comp_CO3terms
!
!   PURPOSE : Calculate H2CO3, HCO3, CO3 from
!             total alkalinity, total CO2, temp, salinity (s), etc.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

  ! !INPUT PARAMETERS:
  integer (KIND=benthos_i4), intent(in) :: k

  real (KIND=benthos_r8), intent(in) :: &
       temp, &
       salt, &
       depth, &
       sit_in, &
       dic_in, &
       ta_in, &
       pt_in, &
       oceanBottomDensity, &
       benthosPorosityLevel

   logical (KIND=benthos_log), intent(in) :: &
       lcomp_co3_coeffs

! !INPUT/OUTPUT PARAMETERS:
  real(KIND=benthos_r8), intent(inout):: &
       phlo, &
       phhi

! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: &
       CO3, &
       HCO3, &
       H2CO3, &
       ph, &
       dic, &
       ta, &
       pt, &
       sit, &
       ! do we need these in the main program???
       bt, &
       st, &
       ft, &
       kw, &
       kb, &
       ks, &
       kf, &
       k1p, &
       k2p, &
       k3p, &
       ksi

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

   logical (KIND=benthos_log) :: k1_k2_pH_tot

   integer(KIND=benthos_i4) :: &
        i

   real(KIND=benthos_r8) :: &
        htotal2, denom, &
        massToVol, &
        volToMass

   real(KIND=benthos_r8) :: &
        htotal,       & ! free concentration of H ion
        k0,k1,k2,     & ! equilibrium constants for CO2 species
        ff              ! fugacity of CO2

  !  Initialize

   CO3 = c0_benthos
   HCO3 = c0_benthos
   H2CO3 = c0_benthos
   !phlo_out = phlo
   !phhi_out = phhi
   ph = phlo
   kw = c0_benthos
   kb = c0_benthos
   ks = c0_benthos
   kf = c0_benthos
   k1p = c0_benthos
   k2p = c0_benthos
   k3p = c0_benthos
   ksi = c0_benthos
   bt = c0_benthos
   st = c0_benthos
   ft = c0_benthos
   dic = dic_in
   ta = ta_in
   pt = pt_in
   sit = sit_in

    !------------------------------------------------------------------
    !   Convert units of output arguments
    !   benthosPorosityLevel = volumeSolute/totalVolume
    !------------------------------------------------------------------
    massToVol = oceanBottomDensity * 1.0e3_benthos_r8 * benthosPorosityLevel
    volToMass = c1_benthos/massToVol

    !------------------------------------------------------------------------
    !   compute thermodynamic CO3 coefficients
    !------------------------------------------------------------------------

    !    IF (lcomp_co3_coeffs) THEN

    k1_k2_pH_tot=.true.
    ! confirmed!
    call comp_co3_coeffs_benthos(k0,k1,k2,ff,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,k, depth, temp, salt, k1_k2_pH_tot)

    !    END IF

    !------------------------------------------------------------------------
    !   compute htotal
    !------------------------------------------------------------------------
    ! confirmed!
    call  comp_htotal_benthos(htotal,dic,ta,pt,sit, &
         k, temp, dic_in, ta_in, pt_in, sit_in, k1, &
         k2,phlo,phhi,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft, &
         volToMass)

    !------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49-51)
    !------------------------------------------------------------------------

    htotal2  = htotal**2
    denom    = 1 / (htotal2 + k1 * htotal + k1 * k2)

    !  this insures H2CO3 + HCO3 + CO3 = 1

    H2CO3 = dic* htotal2 * denom
    HCO3  = dic * k1 * htotal * denom
    CO3   = dic * k1 * k2 * denom

    ph    = -LOG10(htotal)

    H2CO3 = H2CO3 * massToVol
    HCO3  = HCO3 * massToVol
    CO3   = CO3 * massToVol
    dic = dic * massToVol
    ta = ta * massToVol
    pt = pt * massToVol
    sit = sit * massToVol

   end subroutine comp_CO3terms_benthos

 !*****************************************************************************

   SUBROUTINE comp_co3_coeffs_benthos(sk0,sk1,sk2,sff,kw,kb,ks,kf,k1p,k2p,k3p,ksi,&
        bt,st,ft,k,depth,temp,salt,k1_k2_pH_tot)

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    integer(KIND=benthos_i4), intent(in) :: k
    real(KIND=benthos_r8), intent(in) :: &
         depth,    & ! depth (meters)
         temp,     & ! temperature (degrees C)
         salt        ! salinity (PSU)
    LOGICAL(KIND=benthos_log), intent(in) :: k1_k2_pH_tot

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

!maltrud these are scalar versions--need to copy from array(1) due to shr_vmath
    real(KIND=benthos_r8), intent(out) :: &
         sk0,sk1,sk2,     & ! equilibrium constants for CO2 species
         sff, &              ! fugacity of CO2
         kw, kb, ks, kf, &
         k1p, k2p, k3p, &
         ksi, bt, st, ft
    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    real(KIND=benthos_r8), dimension(1) :: &  ! need to be arrays for shr_vmath
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    integer(KIND=benthos_i4) :: i

    real(KIND=benthos_r8) :: &
         press_bar       ! pressure at level k [bars]

    real(KIND=benthos_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
         salt_lim,     & ! bounded salt
         tk,           & ! temperature (K)
         is,           & ! ionic strength
         scl,          & ! chlorinity
         tk100, tk1002, invtk, dlogtk, is2, sqrtis, &
         s2, sqrts, s15, invRtk, arg, &
         deltaV,Kappa,lnKfac,Kfac, & ! pressure correction terms
         log_1_m_1p005em3_s, &
         log_1_p_tot_sulfate_div_ks

    !---------------------------------------------------------------------------
    !initialize output

    sk0 = c0_benthos
    sk1 = c0_benthos
    sk2 = c0_benthos
    sff = c0_benthos
    kw = c0_benthos
    kb = c0_benthos
    ks = c0_benthos
    kf = c0_benthos
    k1p = c0_benthos
    k2p = c0_benthos
    k3p = c0_benthos
    ksi = c0_benthos
    bt = c0_benthos
    st = c0_benthos
    ft = c0_benthos

!   press_bar = ref_pressure(k)
!  below is from POP ref_pressure
    press_bar = 0.059808_benthos_r8*(exp(-0.025_benthos_r8*depth) - c1_benthos)     &
            + 0.100766_benthos_r8*depth + 2.28405e-7_benthos_r8*depth**2

    !---------------------------------------------------------------------------
    !   Calculate all constants needed to convert between various
    !   measured carbon species. References for each equation are
    !   noted in the code.  Once calculated, the constants are stored
    !   and passed in the common block "const". The original version
    !   of this code was based on the code by Dickson in Version 2 of
    !   "Handbook of Methods for the Analysis of the Various Parameters
    !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
    !   p25-26).
    !   Derive simple terms used more than once
    !---------------------------------------------------------------------------

    salt_lim = max(salt,benthos_salt_min)
    tk       = T0_Kelvin_benthos + temp
    tk100    = tk * 1e-2_benthos_r8
    tk1002   = tk100 * tk100
    invtk    = c1_benthos / tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(tk, dlogtk, 1)
#else
    dlogtk   = LOG(tk)
#endif
    invRtk   = (c1_benthos / 83.1451_benthos_r8) * invtk

    is       = 19.924_benthos_r8 * salt_lim / (1000.0_benthos_r8 - 1.005_benthos_r8 * salt_lim)
    is2      = is * is
#ifdef CCSMCOUPLED
    CALL shr_vmath_sqrt(is, sqrtis, 1)
    CALL shr_vmath_sqrt(salt_lim, sqrts, 1)
#else
    sqrtis   = SQRT(is)
    sqrts    = SQRT(salt_lim)
#endif
    s2       = salt_lim * salt_lim
    scl      = salt_lim / 1.80655_benthos_r8

    arg = c1_benthos - 0.001005_benthos_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(arg, log_1_m_1p005em3_s, 1)
#else
    log_1_m_1p005em3_s = LOG(arg)
#endif

    !---------------------------------------------------------------------------
    !   f = k0(1-pH2O)*correction term for non-ideality
    !   Weiss & Price (1980, Mar. Chem., 8, 347-359
    !                 Eq 13 with table 6 values)
    !---------------------------------------------------------------------------

    arg = -162.8301_benthos_r8 + 218.2968_benthos_r8 / tk100 + &
          90.9241_benthos_r8 * (dlogtk + LOG(1e-2_benthos_r8)) - 1.47696_benthos_r8 * tk1002 + &
          salt_lim * (.025695_benthos_r8 - .025225_benthos_r8 * tk100 + 0.0049867_benthos_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ff, 1)
#else
    ff = EXP(arg)
#endif
    sff = ff(1)

    !---------------------------------------------------------------------------
    !   K0 from Weiss 1974
    !---------------------------------------------------------------------------

    arg = 93.4517_benthos_r8 / tk100 - 60.2409_benthos_r8 + 23.3585_benthos_r8 * (dlogtk + LOG(1e-2_benthos_r8)) + &
          salt_lim * (.023517_benthos_r8 - 0.023656_benthos_r8 * tk100 + 0.0047036_benthos_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k0, 1)
#else
    k0 = EXP(arg)
#endif
    sk0 = k0(1)

    !---------------------------------------------------------------------------
    !   k1 = [H][HCO3]/[H2CO3]
    !   k2 = [H][CO3]/[HCO3]
    !   if k1_k2_pH_tot == .true., then use
    !      Lueker, Dickson, Keeling (2000) using Mehrbach et al. data on total scale
    !   otherwise, use
    !      Millero p.664 (1995) using Mehrbach et al. data on seawater scale
    !      this is only present to be consistent w/ OCMIP2 code
    !      it should not be used for new runs
    !      the only reason to use it is to be compatible with prior
    !      long spun up runs that had used it
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 3633.86_benthos_r8 * invtk - 61.2172_benthos_r8 + &
             9.67770_benthos_r8 * dlogtk - 0.011555_benthos_r8 * salt_lim + &
             0.0001152_benthos_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 3670.7_benthos_r8 * invtk - 62.008_benthos_r8 + &
             9.7944_benthos_r8 * dlogtk - 0.0118_benthos_r8 * salt_lim + &
             0.000116_benthos_r8 * s2
    END IF
    arg = -LOG(c10_benthos) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1, 1)
#else
    k1 = EXP(arg)
#endif
    sk1 = k1(1)

    IF (k > 1) THEN
       deltaV = -25.5_benthos_r8 + 0.1271_benthos_r8 * temp
       Kappa  = (-3.08_benthos_r8 + 0.0877_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k1 = k1 * Kfac
    END IF

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 471.78_benthos_r8 * invtk + 25.9290_benthos_r8 - &
             3.16967_benthos_r8 * dlogtk - 0.01781_benthos_r8 * salt_lim + 0.0001122_benthos_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 1394.7_benthos_r8 * invtk + 4.777_benthos_r8 - &
             0.0184_benthos_r8 * salt_lim + 0.000118_benthos_r8 * s2
    END IF
    arg = -LOG(c10_benthos) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2, 1)
#else
    k2 = EXP(arg)
#endif
    sk2 = k2(1)

    IF (k > 1) THEN
       deltaV = -15.82_benthos_r8 - 0.0219_benthos_r8 * temp
       Kappa  = (1.13_benthos_r8 - 0.1475_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k2 = k2 * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   kb = [H][BO2]/[HBO2]
    !   Millero p.669 (1995) using data from Dickson (1990)
    !   CO2SYS states that this in on total pH scale
    !   pressure correction from Millero 1979, p. 1657
    !      omitting salinity contribution
    !---------------------------------------------------------------------------

    arg = (-8966.90_benthos_r8 - 2890.53_benthos_r8 * sqrts - &
           77.942_benthos_r8 * salt_lim + 1.728_benthos_r8 * salt_lim * sqrts - &
           0.0996_benthos_r8 * s2) * invtk + &
          (148.0248_benthos_r8 + 137.1942_benthos_r8 * sqrts + 1.62142_benthos_r8 * salt_lim) + &
          (-24.4344_benthos_r8 - 25.085_benthos_r8 * sqrts - 0.2474_benthos_r8 * salt_lim) * dlogtk + &
          0.053105_benthos_r8 * sqrts * tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kb, 1)
#else
    kb = exp(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -29.48_benthos_r8 + (0.1622_benthos_r8 - 0.002608_benthos_r8 * temp) * temp
       Kappa  = -2.84_benthos_r8 * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kb = kb * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   k1p = [H][H2PO4]/[H3PO4]
    !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4576.752_benthos_r8 * invtk + 115.525_benthos_r8 - &
          18.453_benthos_r8 * dlogtk + &
          (-106.736_benthos_r8 * invtk + 0.69171_benthos_r8) * sqrts + &
          (-0.65643_benthos_r8 * invtk - 0.01844_benthos_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1p, 1)
#else
    k1p = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -14.51_benthos_r8 + (0.1211_benthos_r8 - 0.000321_benthos_r8 * temp) * temp
       Kappa  = (-2.67_benthos_r8 + 0.0427_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k1p = k1p * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   k2p = [H][HPO4]/[H2PO4]
    !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -8814.715_benthos_r8 * invtk + 172.0883_benthos_r8 - &
          27.927_benthos_r8 * dlogtk + &
          (-160.340_benthos_r8 * invtk + 1.3566_benthos_r8) * sqrts + &
          (0.37335_benthos_r8 * invtk - 0.05778_benthos_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2p, 1)
#else
    k2p = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -23.12_benthos_r8 + (0.1758_benthos_r8 - 0.002647_benthos_r8 * temp) * temp
       Kappa  = (-5.15_benthos_r8 + 0.09_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k2p = k2p * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   k3p = [H][PO4]/[HPO4]
    !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -3070.75_benthos_r8 * invtk - 18.141_benthos_r8 + &
          (17.27039_benthos_r8 * invtk + 2.81197_benthos_r8) * sqrts + &
          (-44.99486_benthos_r8 * invtk - 0.09984_benthos_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k3p, 1)
#else
    k3p = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -26.57_benthos_r8 + (0.202_benthos_r8 - 0.003042_benthos_r8 * temp) * temp
       Kappa  = (-4.08_benthos_r8 + 0.0714_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k3p = k3p * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   ksi = [H][SiO(OH)3]/[Si(OH)4]
    !   Millero p.671 (1995) using data from Yao and Millero (1995)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !      apply boric acid values
    !---------------------------------------------------------------------------

    arg = -8904.2_benthos_r8 * invtk + 117.385_benthos_r8 - &
          19.334_benthos_r8 * dlogtk + &
          (-458.79_benthos_r8 * invtk + 3.5913_benthos_r8) * sqrtis + &
          (188.74_benthos_r8 * invtk - 1.5998_benthos_r8) * is + &
          (-12.1652_benthos_r8 * invtk + 0.07871_benthos_r8) * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ksi, 1)
#else
    ksi = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -29.48_benthos_r8 + (0.1622_benthos_r8 - 0.002608_benthos_r8 * temp) * temp
       Kappa  = -2.84_benthos_r8 * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       ksi = ksi * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   kw = [H][OH]
    !   Millero p.670 (1995) using composite data
    !   following DOE Handbook, 0.015 substracted from constant to
    !   approximately convert from SWS pH scale to total pH scale
    !   pressure correction from Millero 1983
    !      note that deltaV coeffs in Millero 1995 are those actually
    !      freshwater deltaV coeffs from Millero 1983
    !---------------------------------------------------------------------------

    arg = -13847.26_benthos_r8 * invtk + 148.9652_benthos_r8 - 23.6521_benthos_r8 * dlogtk + &
          (118.67_benthos_r8 * invtk - 5.977_benthos_r8 + 1.0495_benthos_r8 * dlogtk) * sqrts - &
          0.01615_benthos_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kw, 1)
#else
    kw = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -20.02_benthos_r8 + (0.1119_benthos_r8 - 0.001409_benthos_r8 * temp) * temp
       Kappa  = (-5.13_benthos_r8 + 0.0794_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kw = kw * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   ks = [H][SO4]/[HSO4], free pH scale
    !   Dickson (1990, J. chem. Thermodynamics 22, 113)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4276.1_benthos_r8 * invtk + 141.328_benthos_r8 - 23.093_benthos_r8 * dlogtk + &
          (-13856.0_benthos_r8 * invtk + 324.57_benthos_r8 - 47.986_benthos_r8 * dlogtk) * sqrtis + &
          (35474.0_benthos_r8 * invtk - 771.54_benthos_r8 + 114.723_benthos_r8 * dlogtk) * is - &
          2698.0_benthos_r8 * invtk * is * sqrtis + &
          1776.0_benthos_r8 * invtk * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ks, 1)
#else
    ks = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -18.03_benthos_r8 + (0.0466_benthos_r8 + 0.000316_benthos_r8 * temp) * temp
       Kappa  = (-4.53_benthos_r8 + 0.09_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       ks = ks * Kfac(1)
    END IF

    !---------------------------------------------------------------------
    !   kf = [H][F]/[HF]
    !   Dickson and Riley (1979) -- change pH scale to total
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------

    arg = c1_benthos + (0.1400_benthos_r8 / 96.062_benthos_r8) * (scl) / ks
#ifdef CCSMCOUPLED
       CALL shr_vmath_log(arg, log_1_p_tot_sulfate_div_ks, 1)
#else
    log_1_p_tot_sulfate_div_ks = LOG(arg)
#endif
    arg = 1590.2_benthos_r8 * invtk - 12.641_benthos_r8 + 1.525_benthos_r8 * sqrtis + &
          log_1_m_1p005em3_s + log_1_p_tot_sulfate_div_ks
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kf, 1)
#else
    kf = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -9.78_benthos_r8 - (0.009_benthos_r8 + 0.000942_benthos_r8 * temp) * temp
       Kappa  = (-3.91_benthos_r8 + 0.054_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kf = kf * Kfac(1)
    END IF

    !---------------------------------------------------------------------
    !   Calculate concentrations for borate, sulfate, and fluoride
    !   bt : Uppstrom (1974)
    !   st : Morris & Riley (1966)
    !   ft : Riley (1965)
    !---------------------------------------------------------------------

    bt = 0.000232_benthos_r8 / 10.811_benthos_r8 * scl(1)
    st = 0.14_benthos_r8 / 96.062_benthos_r8 * scl(1)
    ft = 0.000067_benthos_r8 / 18.9984_benthos_r8 * scl(1)

  END SUBROUTINE comp_co3_coeffs_benthos

!*****************************************************************************
!BOP
! !IROUTINE: comp_htotal_benthos
! !INTERFACE:

  subroutine comp_htotal_benthos(htotal,dic,ta,pt,sit, &
       k, temp, dic_in, ta_in, pt_in, sit_in, &
       k1, k2,phlo,phhi,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,volToMass)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
!   SUBROUTINE comp_htotal
!
!   PURPOSE : Calculate htotal from total alkalinity, total CO2,
!             temp, salinity (s), etc.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

    integer(KIND=benthos_i4), intent(in) :: k
    real(KIND=benthos_r8), intent(in) :: &
         temp,     & ! temperature (degrees C)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         k1,k2       ! equilibrium constants for CO2 species

    real(KIND=benthos_r8), intent(in) :: &
         kw, kb, ks, &
         kf, k1p, k2p, &
         k3p, ksi, bt, st, ft, &
         volToMass

! !INPUT/OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(inout) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

! !OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(out) :: &
         htotal, &      ! free concentration of H ion
         dic, &
         ta, &
         pt, &
         sit

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

    integer(KIND=benthos_i4) :: i

    real(KIND=benthos_r8) :: &
         x1, x2, &         ! bounds on htotal for solver
         xacc

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass
    !  (switch vol_to_mass to volToMass using prognostic density)
    !---------------------------------------------------------------------------

    dic  = max(dic_in,benthos_dic_min) * volToMass
    ta   = max(ta_in,benthos_alk_min)  * volToMass
    pt   = max(pt_in,c0_benthos)       * volToMass
    sit  = max(sit_in,c0_benthos)      * volToMass

    x1 = c10_benthos**(-phhi)
    x2 = c10_benthos**(-phlo)

    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------
!      x1_in = x1
!      x2_in  = x2
      xacc = c10_benthos **(-10)

      call  drtsafe_row_benthos(htotal, k, k1, k2, x1, x2, xacc, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

    end subroutine comp_htotal_benthos

!*****************************************************************************
!BOP
! !IROUTINE: drtsafe_row_benthos
! !INTERFACE:

   subroutine drtsafe_row_benthos(soln, k, k1, k2, x1, x2, xacc, dic,kw,kb,ks,kf,k1p,&
        k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
! from drtsafe_row:
!   Vectorized version of drtsafe, which was a modified version of
!      Numerical Recipes algorithm.
!   Keith Lindsay, Oct 1999
!
!   Algorithm comment :
!      Iteration from Newtons method is used unless it leaves
!      bracketing interval or the dx is > 0.5 the previous dx.
!      In that case, bisection method is used.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

    integer(KIND=benthos_i4), intent(in) :: k
    real(KIND=benthos_r8), intent(in) :: &
         k1, k2,&
         xacc, &
         dic, &
         kw, kb, ks, kf, &
         k1p, k2p, k3p, &
         ksi, bt, st, ft, &
         ta, pt, sit

!  !INPUT/OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(INOUT) :: x1, x2

!  !OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(OUT) :: soln

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

    logical(KIND=benthos_log) :: leave_bracket, dx_decrease, mask
    integer(KIND=benthos_i4) ::  i, it, maxit, max_bracket_grow_it
    real(KIND=benthos_r8) :: temp
    real(KIND=benthos_r8) :: xlo, xhi, flo, fhi, f, df, dxold, dx

!---------------------------------------------------------------------------
!   bracket root at each location and set up first iteration
!---------------------------------------------------------------------------
    max_bracket_grow_it = 3
    maxit = 100
    it = 0

    DO

       call talk_row_benthos(flo,df,k1, k2, x1, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

       call talk_row_benthos(fhi,df,k1, k2, x2, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

       mask = (flo .gt.  c0_benthos .and.  fhi .gt. c0_benthos) .or. &
            (flo .lt. c0_benthos .and. fhi .lt. c0_benthos)

       if (.not. mask) EXIT

       it = it + 1

    !    if (it .gt. max_bracket_grow_it) then
    !         'bounding bracket for pH solution not found'
    !        'drtsafe_row'
    !         return
    !         CALL shr_sys_abort('bounding bracket for pH solution not found')
    !       end if

       dx = sqrt(x2 / x1)
       x2 = x2 * dx
       x1 = x1 / dx
    end do

    if (flo .lt. c0_benthos) then
       xlo = x1
       xhi = x2
    else
       xlo = x2
       xhi = x1
       temp = flo
       flo = fhi
       fhi = temp
    end if
    soln = p5_benthos* (xlo+xhi)

    dxold = abs(xlo - xhi)
    dx = dxold

    call talk_row_benthos(f,df,k1, k2, soln, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)
    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

    mask = .true.
    do it = 1,maxit
       leave_bracket = ((soln - xhi) * df - f) * &
       ((soln - xlo) * df - f) .ge. c0_benthos
       dx_decrease = abs(c2_benthos * f) .le. abs(dxold * df)
       if (leave_bracket .or. .not. dx_decrease) then
          dxold = dx
          dx = p5_benthos * (xhi - xlo)
          soln = xlo + dx
          if (xlo .eq. soln) mask = .false.
       else
          dxold = dx
          dx = -f / df
          temp = soln
          soln = soln + dx
          if (temp .eq. soln) mask = .false.
       end if
       if (abs(dx) .lt. xacc) mask = .false.

       if (.not. mask) return

       call talk_row_benthos(f,df,k1, k2, soln, dic,kw,kb,ks,kf,&
            k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

       if (f .lt. c0_benthos) then
          xlo = soln
          flo = f
       else
          xhi = soln
          fhi = f
       end  if
    end do

#ifdef CCSMCOUPLED
!   CALL shr_sys_abort('lack of convergence in drtsafe_row')
#endif

  end subroutine drtsafe_row_benthos

!*****************************************************************************
!BOP
! !IROUTINE: talk_row_benthos
! !INTERFACE:

   subroutine talk_row_benthos(fn,df, &
        k1, k2, x,dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

! !DESCRIPTION:
! from talk_row:
!---------------------------------------------------------------------------
!   This routine computes TA as a function of DIC, htotal and constants.
!   It also calculates the derivative of this function with respect to
!   htotal. It is used in the iterative solution for htotal. In the call
!   "x" is the input value for htotal, "fn" is the calculated value for
!   TA and "df" is the value for dTA/dhtotal.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
     real(KIND=benthos_r8), intent(in) :: &
          k1, k2, x, &
          dic, kw, kb, ks, kf, &
          k1p, k2p, k3p, ksi, &
          bt, st, ft, ta, pt, sit

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

     real(KIND=benthos_r8), intent(out) :: &
          fn, &
          df

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

     integer(KIND=benthos_i4) :: i
     real(KIND=benthos_r8) :: &
         x1, x1_r, x2, x2_r, x3, k12, k12p, k123p, &
         a, a_r, a2_r, da, b, b_r, b2_r, db, c, c_r, &
         kb_p_x1_r, ksi_p_x1_r, c1_p_c_ks_x1_r_r, c1_p_kf_x1_r_r

     x1 = x
     x1_r = c1_benthos / x1
     x2 = x1 * x1
     x2_r = x1_r * x1_r
     x3 = x2 * x1
     k12 = k1 * k2
     k12p = k1p * k2p
     k123p = k12p * k3p
     a = x3 + k1p * x2 + k12p * x1 + k123p
     a_r = c1_benthos/ a
     a2_r = a_r * a_r
     da = 3.0_benthos_r8 * x2 + c2_benthos * k1p * x1 + k12p
     b = x2 + k1 * x1 + k12
     b_r = c1_benthos / b
     b2_r = b_r * b_r
     db = c2_benthos * x1 + k1
     c = c1_benthos + st/ks
     c_r = c1_benthos/ c
     kb_p_x1_r = c1_benthos / (kb + x1)
     ksi_p_x1_r = c1_benthos/ (ksi + x1)
     c1_p_c_ks_x1_r_r = c1_benthos / (c1_benthos + c * ks * x1_r)
     c1_p_kf_x1_r_r = c1_benthos/ (c1_benthos+ kf * x1_r)

    !---------------------------------------------------------------------
    !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
    !---------------------------------------------------------------------

     fn = k1 * dic * x1 * b_r &
          + c2_benthos * dic * k12 * b_r &
          + bt * kb * kb_p_x1_r &
          + kw * x1_r &
          + pt * k12p * x1 * a_r &
          + c2_benthos * pt * k123p * a_r &
          + sit * ksi * ksi_p_x1_r &
          - x1 * c_r &
          - st * c1_p_c_ks_x1_r_r &
          - ft * c1_p_kf_x1_r_r &
          - pt * x3 * a_r &
          - ta

     !---------------------------------------------------------------------
     !   df = d(fn)/dx
     !---------------------------------------------------------------------

     df = k1 * dic * (b - x1 * db) * b2_r &
          - c2_benthos * dic * k12 * db * b2_r &
          - bt * kb * kb_p_x1_r * kb_p_x1_r &
          - kw * x2_r &
          + (pt * k12p * (a - x1 * da)) * a2_r &
          - c2_benthos * pt * k123p * da * a2_r &
          - sit * ksi * ksi_p_x1_r * ksi_p_x1_r &
          - c1_benthos * c_r &
          - st * c1_p_c_ks_x1_r_r * c1_p_c_ks_x1_r_r * (c * ks * x2_r) &
          - ft * c1_p_kf_x1_r_r * c1_p_kf_x1_r_r * kf * x2_r &
          - pt * x2 * (3.0_benthos_r8 * a - x1 * da) * a2_r

   end subroutine talk_row_benthos

!*****************************************************************************
!BOP
! !IROUTINE: comp_co3_sat_benthos
! !INTERFACE:

subroutine comp_co3_sat_benthos(sat_calc, sat_arag, sat_mgcalc,&
                       K_arag,K_calc,K_mgcalc,temp,salt, &
                       depth,mg_ions,co3_mol,oceanBottomDensity,benthosPorosityLevel)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
!   Based on SUBROUTINE after comp_co3_sat_vals
!
!   PURPOSE : Calculate co3 concentration at calcite,  aragonite and 15! Mg-calcite saturation
!             from temp, salinity (s), press
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

    real(KIND=benthos_r8), intent(in) :: &
         depth,      & ! depth (meters)
         temp,       & ! temperature (degrees C)
         salt,       & ! salinity (PSU)
         mg_ions,    & ! magnesium ion concentration (mmol/m3)
         co3_mol,    & ! bottom ocean density
         oceanBottomDensity, &  ! mass per volume of liquid
         benthosPorosityLevel  ! volume solute/total volume

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:
    real(KIND=benthos_r8), intent(out) :: &
         sat_calc,  & ! co3 concentration at calcite saturation : change to sat_calc  for saturation state
         sat_arag,  & ! co3 concentration at aragonite saturation
         sat_mgcalc ! co3 concentration at 15% magnesium calcite

    real(KIND=benthos_r8), intent(out):: &
         K_calc,       & ! thermodynamic constant for calcite
         K_arag,       & ! thermodynamic constant for aragonite
         K_mgcalc       ! thermodynamic constant for mg calcite

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

    integer(KIND=benthos_i4) :: i, j

    real(KIND=benthos_r8) :: &
         press_bar, &          ! pressure at level k [bars]
         volToMass

    real(KIND=benthos_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
         salt_lim,     & ! bounded salt
         tk,             & ! temperature (K)
         log10tk,invtk,sqrts,s15,invRtk,arg,&
         deltaV,Kappa, & ! pressure correction terms
         lnKfac,Kfac,    & ! pressure correction terms
         inv_K_calc,     & ! inverse K_calc
         inv_K_arag,     & ! inverse K_arag
         inv_K_mgcalc,   & ! inverse K_mgcalc
         inv_Ca, &         ! inverse of Calcium concentration
         ca_ions, &
         ca85, &
         mg15

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !  mmol/kg to mmol/m^3
    !---------------------------------------------------------------------------

    volToMass = c1_benthos/1.0e3_benthos_r8/oceanBottomDensity/benthosPorosityLevel

    press_bar = 0.059808 * (exp(-0.025_benthos_r8 * depth) - c1_benthos) + &
         0.100766_benthos_r8 * depth + 2.28405_benthos_r8*c10_benthos**(-7) * depth**2


    salt_lim = max(salt,benthos_salt_min)
    tk       = T0_Kelvin_benthos + temp

#ifdef CCSMCOUPLED
       CALL shr_vmath_log(tk, log10tk, 1)
#else
       log10tk  = log(tk)
#endif

       !maltrud dividing by non-shared log(10) here.  i guess since it is a scalar?

       log10tk  = log10tk/log(c10_benthos)
       invtk    = c1_benthos/tk
       invRtk   = (c1_benthos/83.1451_benthos_r8)*invtk

#ifdef CCSMCOUPLED
       CALL shr_vmath_sqrt(salt_lim, sqrts, 1)
#else
       sqrts    =sqrt(salt_lim)
#endif
       s15      = sqrts * salt_lim

       !------------------------------------------------------------------------
       !   Compute K_calc, K_arag, apply pressure factor
       !   Mucci, Amer. J. of Science 283:781-799, 1983 & Millero 1979
       !------------------------------------------------------------------------

       arg = -171.9065_benthos_r8 - 0.077993_benthos_r8 * tk + 2839.319_benthos_r8 * invtk + &
            71.595_benthos_r8 * log10tk + &
             (-0.77712_benthos_r8 + 0.0028426_benthos_r8 * tk + 178.34_benthos_r8 * invtk) * sqrts - &
             0.07711_benthos_r8 * salt_lim + 0.0041249_benthos_r8 * s15

       arg = log(c10_benthos) * arg
#ifdef CCSMCOUPLED
!       CALL shr_vmath_exp(arg, K_calc, 1)
#else
       K_calc = exp(arg(1))
#endif

!       if (k .ge. 1) then
          deltaV = -48.76_benthos_r8 + 0.5304_benthos_r8 * temp
          Kappa  = (-11.76_benthos_r8 + 0.3692_benthos_r8 * temp) * 0.001_benthos_r8
          lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk

#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
          Kfac = exp(lnKfac)
#endif

          K_calc = K_calc * Kfac(1)    ! mol/kg
!       end if

       arg = -171.945_benthos_r8 - 0.077993_benthos_r8 * tk + 2903.293_benthos_r8 * invtk + &
            71.595_benthos_r8 * log10tk +  &
           (-0.068393_benthos_r8 + 0.0017276_benthos_r8 * tk + 88.135_benthos_r8 * invtk) * sqrts -  &
            0.10018_benthos_r8 * salt_lim + 0.0059415_benthos_r8 * s15

       arg = log(c10_benthos) * arg
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(arg, K_arag, 1)
#else
      K_arag = exp(arg(1))
#endif

!       if  (k .ge. 1) then
          deltaV = deltaV + 2.8_benthos_r8
          lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk

#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else

          Kfac = exp(lnKfac) !EXP(lnKfac)
#endif
          K_arag = K_arag * Kfac(1)    !mol/kg
!       end if

       !
       ! Mg -Calcite
       !  The solubility factor is assumed to be 21! greater than that of  aragonite
       !  Morse et al (2006).

       K_mgcalc = K_arag * 1.21_benthos_r8    !mol/kg

    !------------------------------------------------------------------
    !   Compute CO3 concentration at calcite & aragonite saturation
    !  but, we are now keeping track of calcium ions...
    !------------------------------------------------------------------

    !   inv_Ca(1) = (35.0_benthos_r8 / 0.01028_benthos_r8) / salt_lim(1)
   !   co3_sat_calc = K_calc(1) * inv_Ca(1)   ! only need to be arrays if using shared math
   !   co3_sat_arag = K_arag(1) * inv_Ca(1)
   !   co3_sat_mgcalc = K_mgcalc(1) * inv_Ca(1)

     inv_K_calc = c1_benthos/K_calc  !(kg/mol)^2
     ca_ions = 0.01028_benthos_r8/35.0_benthos_r8*salt_lim  !(mol/kg)
     sat_calc = co3_mol * volToMass * ca_ions(1) * inv_K_calc(1)  ! (co3_mol from mmol/m3 to mol/kg)  now unitless

     inv_K_arag = c1_benthos/K_arag
     sat_arag = co3_mol* volToMass * ca_ions(1) * inv_K_arag(1)

     inv_K_mgcalc = c1_benthos/K_mgcalc
     ca85 = (ca_ions)**(0.85_benthos_r8)
     mg15 = (mg_ions)**(0.15_benthos_r8)
     sat_mgcalc = ca85(1) * mg15(1) * inv_K_mgcalc(1)

   end subroutine comp_co3_sat_benthos

!*****************************************************************************
!BOP
! !IROUTINE: carbonateRateConstants
! !INTERFACE:

subroutine carbonateRateConstants(carbonateRates,benthosTracerBulkLevels,&
     dt,sat_calc, sat_arag,sat_mgcalc, oceanBottomDensity)

! !DESCRIPTION:
!  compute the rate constants for calcium carbonate dissociation
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       dt, &
       sat_calc, &
       sat_arag, &
       sat_mgcalc, &
       oceanBottomDensity

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulkLevels

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nCarbonateReactions), intent(out) :: &
       carbonateRates

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  integer (KIND=benthos_i4) :: iReactions

  real (KIND=benthos_r8) :: &
       k_p1, k_p2, k_p3  ! mmol/mmol of calcite/s

   !  Constants from Krumins et al.  Something odd about the units

   ! k_p1 = 40.0_benthos_r8/sec_per_year/mM_umolperL  ! /(mmol/m3)/s   from  40.0 per mM per yr
   ! k_p2 = 110.0_benthos_r8/sec_per_year/mM_umolperL ! /(mmol/m3)/s from 110.0 per mM per yr
   ! k_p3 = 58.0_benthos_r8/sec_per_year/mM_umolperL ! /(mmol/m3)/s from 58 per mM per yr

   ! Constants from Walter and Morse (1985) reduced by inhibition factor of 0.1 as in Krumins et al

   k_p1 = kp1 * inhibition / sec_per_hour * mm_calcite * sediment_density! mmol/mmol of calcite/s*kg/m^3 !* sediment_density  ! barnacle Balanus
   k_p2 = kp2 * inhibition / sec_per_hour * mm_aragonite * sediment_density! mmol/mmol of aragonite/s*kg/m^3 !* sediment_density ! green algae Halimeda
   k_p3 = kp3 * inhibition / sec_per_hour * mm_mgcalcite * sediment_density ! !mmol/mmmol  of mg-calcite /s*kg/m^3 !* sediment_density ! foraminifier Peneroplis

   do iReactions = 1,nCarbonateReactions
      carbonateRates(iReactions) = c0_benthos
   end do

   carbonateRates(1) = k_p1 * (c1_benthos - min(c1_benthos,sat_calc))**n1*benthosTracerBulkLevels(caco3a_ind)   ! mmol/m^3/s
   carbonateRates(2) = k_p2 * (c1_benthos - min(c1_benthos,sat_arag))**n2*benthosTracerBulkLevels(caco3b_ind)   ! mmol/m^3/s
   carbonateRates(3) = k_p3 * (c1_benthos - min(c1_benthos,sat_mgcalc))**n3*benthosTracerBulkLevels(camgco3_ind) ! mmol/m^3/s

 end subroutine carbonateRateConstants

!*****************************************************************************
!BOP
! !IROUTINE: carbonateBenthosReactions
! !INTERFACE:

subroutine carbonateBenthosReactions (carbonateSourceTendk,carbonateSinkTendk,&
                    carbonateSourceStoich,carbonateSinkStoich,carbonateRates,dt)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       dt

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       carbonateRates  ! 1/s

  real(KIND=benthos_r8), dimension(:,:), intent(in) :: &
       carbonateSourceStoich, &
       carbonateSinkStoich

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nBenthicTracers), intent(inout) :: &
       carbonateSourceTendk, &
       carbonateSinkTendk

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  integer (KIND=benthos_i4) :: &
       iTracers, &
       iReactions

  do iTracers = 1,nBenthicTracers

     carbonateSourceTendk(iTracers) = c0_benthos
     carbonateSinkTendk(iTracers) = c0_benthos

      do iReactions = 1,nCarbonateReactions

         if (carbonateSourceStoich(iTracers,iReactions) .gt. c0_benthos .and. &
              carbonateRates(iReactions) .gt. puny) then

            carbonateSourceTendk(iTracers) = carbonateSourceTendk(iTracers) +&
                 carbonateSourceStoich(iTracers,iReactions)*carbonateRates(iReactions)*dt

         end if

         if (carbonateSinkStoich(iTracers,iReactions) .gt. c0_benthos .and. &
              carbonateRates(iReactions) .gt. puny) then
            carbonateSinkTendk(iTracers) = carbonateSinkTendk(iTracers) +&
                 carbonateSinkStoich(iTracers,iReactions)*carbonateRates(iReactions)*dt

         end if
      end do ! iReactions
   end do    ! iTracers

 end subroutine carbonateBenthosReactions

!*****************************************************************************
!BOP
! !IROUTINE: secondaryRateConstants
! !INTERFACE:

 subroutine secondaryRateConstants(secondaryRates, benthosTracerBulkLevels)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real(KIND=benthos_r8), dimension(:), intent(in) ::&
        benthosTracerBulkLevels

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

   real(KIND=benthos_r8), dimension(nSecondaryReactions), intent(out) :: &
        secondaryRates
!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  integer (KIND=benthos_i4) :: &
       iReactions

   !  Constants from Reed et al 2011 in benthos_parms.F90
   do iReactions = 1,nSecondaryReactions
      secondaryRates(iReactions) = c0_benthos
   end do

   ! all rates changed to mmol/m3/s

   secondaryRates(1) = max(c0_benthos,k_s1*benthosTracerBulkLevels(nh4_ind)  * benthosTracerBulkLevels(o2_ind))  !mmol/m3/s

   secondaryRates(2) = max(c0_benthos,k_s2 * benthosTracerBulkLevels(mn_ind)* benthosTracerBulkLevels(o2_ind))  ! mmol/m3/s

   secondaryRates(3) =  max(c0_benthos,k_s3*benthosTracerBulkLevels(fe_ind) * benthosTracerBulkLevels(o2_ind) )!mmol/m3/s

   secondaryRates(4) = max(c0_benthos,k_s4 * benthosTracerBulkLevels(fes_ind)* benthosTracerBulkLevels(o2_ind)*sediment_density)  ! mmol/m3/s

   secondaryRates(5) = max(c0_benthos,k_s5 * benthosTracerBulkLevels(fes2_ind)* benthosTracerBulkLevels(o2_ind)*sediment_density)  !  umol/g/s

   secondaryRates(6) = max(c0_benthos,k_s6*benthosTracerBulkLevels(h2s_ind)  * benthosTracerBulkLevels(o2_ind) )   !mmol/m3/s

   secondaryRates(7) =max(c0_benthos,k_s7*benthosTracerBulkLevels(ch4_ind)  * benthosTracerBulkLevels(o2_ind) )   !mmol/m3/s

   secondaryRates(8) = max(c0_benthos,k_s8 * benthosTracerBulkLevels(mno2a_ind)* benthosTracerBulkLevels(fe_ind)*sediment_density )  ! umol/g/s

   secondaryRates(9) =max(c0_benthos,k_s8 * benthosTracerBulkLevels(mno2b_ind)* benthosTracerBulkLevels(fe_ind)*sediment_density)  !  umol/g/s

   secondaryRates(10) = max(c0_benthos,k_s9 * benthosTracerBulkLevels(mno2a_ind)* benthosTracerBulkLevels(h2s_ind)*sediment_density)  ! umol/g/s

   secondaryRates(11) = max(c0_benthos,k_s9 * benthosTracerBulkLevels(mno2b_ind)*benthosTracerBulkLevels(h2s_ind)*sediment_density)  !  umol/g/s

   secondaryRates(12) =max(c0_benthos,k_s10 * benthosTracerBulkLevels(feoh3a_ind)*benthosTracerBulkLevels(h2s_ind)*sediment_density) ! umol/g/s

   secondaryRates(13) = max(c0_benthos,k_s10 * benthosTracerBulkLevels(feoh3b_ind)* benthosTracerBulkLevels(h2s_ind)*sediment_density)   ! umol/g/s

   secondaryRates(14) = max(c0_benthos,k_s11*benthosTracerBulkLevels(h2s_ind)* benthosTracerBulkLevels(fe_ind))    !mmol/m3/s

   secondaryRates(15) = max(c0_benthos,k_s12*benthosTracerBulkLevels(ch4_ind)* benthosTracerBulkLevels(so4_ind)) ! mmol/m3/s

   secondaryRates(16) = max(c0_benthos,k_s13 * benthosTracerBulkLevels(s_ind)*sediment_density)   ! mmol/m3/s

   secondaryRates(17) = max(c0_benthos,k_s14*benthosTracerBulkLevels(s_ind) * benthosTracerBulkLevels(fes_ind)*sediment_density) !  umol/g/s

   secondaryRates(18) = max(c0_benthos,k_s15 * benthosTracerBulkLevels(feoh3a_ind)) !umol/g/s

   secondaryRates(19) = max(c0_benthos,k_s16 * benthosTracerBulkLevels(mno2a_ind)) !umol/g/s

 end subroutine secondaryRateConstants

!*****************************************************************************
!BOP
! !IROUTINE: secondaryBenthosReactions
! !INTERFACE:

subroutine secondaryBenthosReactions(secondarySourceTendk,secondarySinkTendk, &
     secondarySourceStoich,secondarySinkStoich,secondaryRates,dt)

! !DESCRIPTION:
!
!  compute tendencies for secondary reactions
!
! s1. O2 + NH4/2 + HCO3 --> NO3/2 + CO2 + 3/2H2O
! ** s2. O2 + 2Mn + 4HCO3 --> 2MnO2^a + 4CO2 + 2H2O
! ** s3. O2 + 4Fe + 8HCO3 + 2 H20 + 4 gamma H2PO4 + 24 gamma H^+ --> 4Fe(OH)3^a + 4 gamma [Fe-P]^a + 8CO2 + 16 gamma H2O
! s4. O2 + FeS/2 --> SO4/2 + Fe/2
! s5. O2+ 2/7FeS2 +  2/7H2O --> 4/7SO4 +  2/7Fe+ 4/7H
! s6. O2 + H2S/2 + HCO3 -->  SO4/2 + CO2 + H2O
! s7. O2 + CH4 --> CO2 + 4H^+
! s8. MnO2a + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gamma H --> 2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s9)s8b. MnO2b + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gammaH --> 2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s10)s9. Mn02a + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s11)s9b. Mn02b + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s12)s10. Fe(OH)3a+   gamma[Fe-P]a + H2S/2+ 2CO2 + 4 gamma H2O--> Fe +   gamma H2PO4 + S/2 + 2HCO3 + (2+6gamma)H
! (s13)s10b. Fe(OH)3b+   gamma[Fe-P]b + H2S/2+ 2CO2 + 4 gamma H2O--> Fe +  gamma H2PO4+ S/2 + 2HCO3 + (2+6 gamma)H
! (s14)s11. Fe + H2S -->  FeS + 2H
! (s15)s12. SO4 +  CH4 +  CO2 -->  2HCO3 + H2S
! (s16)s13. S +  H2O -->  3/4H2S +  SO4/4 +  H/2
! (s17)s14. FeS +  S --> FeS2
! (s18)s15. Fe(OH)3a + gamma [Fe-P]a --> Fe(OH)3b + gamma [Fe-P]b
! (s19)s16 MnO2a --> MnO2b
!
! 19 total secondary equations with 3 repetitions
! gamma is the fraction of  co-precipitated P with iron

!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       secondaryRates

  real(KIND=benthos_r8), dimension(:,:),intent(in) :: &
       secondarySourceStoich, &
       secondarySinkStoich

  real(KIND=benthos_r8), intent(in) :: &
       dt

  ! !INPUT/OUTPUT PARAMETERS:
  real(KIND=benthos_r8), dimension(nBenthicTracers), intent(inout) :: &
       secondarySourceTendk, &
       secondarySinkTendk

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  integer (KIND=benthos_i4) :: &
       iTracers, &
       iReactions

  do iTracers = 1,nBenthicTracers

     secondarySourceTendk(iTracers) = c0_benthos
     secondarySinkTendk(iTracers) = c0_benthos

      do iReactions = 1,nSecondaryReactions

         if (secondarySourceStoich(iTracers,iReactions) .gt. c0_benthos) then
            secondarySourceTendk(iTracers) = secondarySourceTendk(iTracers) + &

                 secondarySourceStoich(iTracers,iReactions)*secondaryRates(iReactions)*dt

         end if

         if (secondarySinkStoich(iTracers,iReactions) .gt. c0_benthos) then

            secondarySinkTendk(iTracers) = secondarySinkTendk(iTracers) + &
                 secondarySinkStoich(iTracers,iReactions)*secondaryRates(iReactions)*dt
         end if

      end do ! iReactions
   end do ! iTracers

   end subroutine secondaryBenthosReactions

!***********************************************************************
!BOP
! !IROUTINE: benthos_SurfaceFluxes
! !INTERFACE:

 subroutine benthos_SurfaceFluxes(benthos_input, benthos_forcing,   &
                              benthos_flux_diagnostic_fields, &
                              numColumnsMax, numColumns)

! !DESCRIPTION:
!  Compute surface fluxes for tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(benthos_input_type),   intent(in )   :: benthos_input

  integer (KIND=benthos_i4), intent(in)  :: numColumnsMax, numColumns

! !INPUT/OUTPUT PARAMETERS:

  type(benthos_forcing_type), intent(inout) :: benthos_forcing
  type(benthos_flux_diagnostics_type), intent(inout) :: benthos_flux_diagnostic_fields

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (KIND=benthos_r8) :: &
      seaSurfaceTemp
!-----------------------------------------------------------------------
!EOC

 end subroutine benthos_SurfaceFluxes

  !*****************************************************************************

end module ocn_benthos_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
