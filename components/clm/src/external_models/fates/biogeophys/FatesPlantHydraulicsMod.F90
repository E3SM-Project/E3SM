module FatesPlantHydraulicsMod


   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!99
   ! (TODO: THE ROW WIDTH ON THIS MODULE ARE TOO LARGE. NAG COMPILERS
   !  WILL FREAK IF LINES ARE TOO LONG.  BEFORE SUBMITTING THIS TO 
   !  MASTER WE NEED TO GO THROUGH AND GET THESE LINES BELOW
   !  100 spaces (for readability), or 130 (for NAG)
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!99

   use FatesGlobals, only      : endrun => fates_endrun
   use FatesGlobals, only      : fates_log

   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : fates_huge
   use FatesConstantsMod, only : denh2o => dens_fresh_liquid_water
   use FatesConstantsMod, only : grav => grav_earth

   use EDTypesMod        , only : ed_site_type
   use EDTypesMod        , only : ed_patch_type
   use EDTypesMod        , only : ed_cohort_type

   use FatesInterfaceMod  , only : bc_in_type
   use FatesInterfaceMod  , only : bc_out_type
   use FatesInterfaceMod  , only : hlm_numlevsoil
   
   use FatesHydraulicsMemMod, only: ed_site_hydr_type
   use FatesHydraulicsMemMod, only: ed_patch_hydr_type
   use FatesHydraulicsMemMod, only: ed_cohort_hydr_type
   use FatesHydraulicsMemMod, only: npool_leaf
   use FatesHydraulicsMemMod, only: npool_tot
   use FatesHydraulicsMemMod, only: npool_stem
   use FatesHydraulicsMemMod, only: numLWPmem
   use FatesHydraulicsMemMod, only: npool_troot
   use FatesHydraulicsMemMod, only: npool_aroot
   use FatesHydraulicsMemMod, only: n_porous_media
   use FatesHydraulicsMemMod, only: nshell
   use FatesHydraulicsMemMod, only: npool_ag
   use FatesHydraulicsMemMod, only: npool_bg
   use FatesHydraulicsMemMod, only: porous_media
   use FatesHydraulicsMemMod, only: nlevsoi_hyd

   use EDPftvarcon, only : EDPftvarcon_inst

   ! CIME Globals
   use shr_log_mod , only      : errMsg => shr_log_errMsg
   use shr_infnan_mod   , only : isnan => shr_infnan_isnan
   

   implicit none

   private
   integer, parameter :: van_genuchten = 1
   integer, parameter :: campbell      = 2
   integer :: iswc = campbell
   
   ! 1=leaf, 2=stem, 3=troot, 4=aroot
   ! Several of these may be better transferred to the parameter file in due time (RGK)

   integer, public :: use_ed_planthydraulics    =  1      ! 0 => use vanilla btran
                                                          ! 1 => use BC hydraulics; 
                                                          ! 2 => use CX hydraulics
   logical, public :: do_dqtopdth_leaf          = .false. ! should a nonzero dqtopdth_leaf
                                                          ! term be applied to the plant 
                                                          ! hydraulics numerical solution?
   logical, public :: do_dyn_xylemrefill        = .true.  ! should the dynamics of xylem refilling 
                                                          ! (i.e., non-instantaneous) be considered 
                                                          ! within plant hydraulics?
   logical, public :: do_kbound_upstream        = .false. ! should the hydraulic conductance at the 
                                                          ! boundary between nodes be taken to be a
                                                          ! function of the upstream loss of 
                                                          ! conductivity (flc)?
   logical, public :: do_growthrecruiteffects   = .true. ! should size- or root length-dependent 
                                                          ! hydraulic properties and states be 
                                                          ! updated every day when trees grow or 
                                                          ! when recruitment happens?
   logical, public :: do_static_ed              = .true.  ! should growth, mortality, and patch
                                                          ! dynamics be turned off?

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   public :: hydraulics_drive
   public :: InitHydrSites
   public :: HydrSiteColdStart
   public :: BTranForHLMDiagnosticsFromCohortHydr
   public :: InitHydrCohort
   public :: DeallocateHydrCohort
   public :: UpdateH2OVeg
   public :: CopyCohortHydraulics
   public :: FuseCohortHydraulics
   public :: updateSizeDepTreeHydProps
   public :: updateSizeDepTreeHydStates
   public :: initTreeHydStates
   public :: updateSizeDepRhizHydProps

   !------------------------------------------------------------------------------
   ! 01/18/16: Created by Brad Christoffersen
   ! 02/xx/17: Refactoring by Ryan Knox and Brad Christoffersen
   !------------------------------------------------------------------------------
   
contains 

   !------------------------------------------------------------------------------
   subroutine hydraulics_drive( nsites, sites, bc_in,bc_out,dtime )
    
      ! ARGUMENTS:
      ! -----------------------------------------------------------------------------------
      integer,intent(in)                      :: nsites
      type(ed_site_type),intent(inout),target :: sites(nsites)
      type(bc_in_type),intent(in)             :: bc_in(nsites)
      type(bc_out_type),intent(inout)         :: bc_out(nsites)
      real(r8),intent(in)                     :: dtime

      write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
      call endrun(msg=errMsg(sourcefile, __LINE__))
	     
 end subroutine Hydraulics_Drive
  
  ! =====================================================================================

  subroutine initTreeHydStates(cc_p, bc_in)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:

    ! !ARGUMENTS:
    type(ed_cohort_type), intent(inout), target  :: cc_p ! current cohort pointer
    type(bc_in_type)    , intent(in)             :: bc_in 
    
    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))
    
  end subroutine initTreeHydStates

  ! =====================================================================================

  subroutine updateSizeDepTreeHydProps(cc_p,bc_in)
    !
    ! !DESCRIPTION: Updates absorbing root length (total and its vertical distribution)
    !   as well as the consequential change in the size of the 'representative' rhizosphere
    !   shell radii, volumes
    !
    ! !USES:
    use FatesConstantsMod  , only : pi_const
    use shr_sys_mod        , only : shr_sys_abort
    !
    ! !ARGUMENTS:
    type(ed_cohort_type)   , intent(inout), target  :: cc_p    ! current cohort pointer
    type(bc_in_type)       , intent(in)             :: bc_in   ! Boundary Conditions

    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))

  end subroutine updateSizeDepTreeHydProps

  ! =====================================================================================

  subroutine updateSizeDepTreeHydStates(cc_p)
    !
    ! !DESCRIPTION: 
    !
    ! !USES:
    ! !ARGUMENTS:
    type(ed_cohort_type)   , intent(inout), target  :: cc_p ! current cohort pointer
    ! 
    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))
   
  end subroutine updateSizeDepTreeHydStates

  ! =====================================================================================
  
  subroutine CopyCohortHydraulics(newCohort, oldCohort)

     ! Arguments
     type(ed_cohort_type), intent(inout), target :: newCohort
     type(ed_cohort_type), intent(inout), target :: oldCohort
     
     write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
     call endrun(msg=errMsg(sourcefile, __LINE__))
     
  end subroutine CopyCohortHydraulics
  
  ! =====================================================================================

  subroutine FuseCohortHydraulics(currentCohort, nextCohort, bc_in, newn)

     
     type(ed_cohort_type), intent(inout), target :: currentCohort ! current cohort
     type(ed_cohort_type), intent(inout), target :: nextCohort    ! next (donor) cohort
     type(bc_in_type), intent(in)                :: bc_in
     real(r8), intent(in)                        :: newn
     
     write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
     call endrun(msg=errMsg(sourcefile, __LINE__))
     
  end subroutine FuseCohortHydraulics



  ! =====================================================================================
  ! Initialization Routines
  ! =====================================================================================

  subroutine InitHydrCohort(currentCohort)

    ! Arguments
    type(ed_cohort_type), target :: currentCohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))
    
  end subroutine InitHydrCohort

  ! =====================================================================================

  subroutine DeallocateHydrCohort(currentCohort)

    ! Arguments
    type(ed_cohort_type), target :: currentCohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr

    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))

  end subroutine DeallocateHydrCohort


  ! =====================================================================================

   subroutine InitHydrSites(sites)

       ! Arguments
       type(ed_site_type),intent(inout),target :: sites(:)

       ! Locals
       integer :: nsites
       integer :: s
       type(ed_site_hydr_type),pointer :: csite_hydr
       
       write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end subroutine InitHydrSites

    ! ===================================================================================

  subroutine HydrSiteColdStart(sites, bc_in)
       

     ! Arguments
     type(ed_site_type),intent(inout),target :: sites(:)
     type(bc_in_type),intent(inout)      :: bc_in(:)
     
     write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
     call endrun(msg=errMsg(sourcefile, __LINE__))

  end subroutine HydrSiteColdStart

  ! =====================================================================================

  subroutine UpdateH2OVeg(nsites,sites,bc_out)

     ! ----------------------------------------------------------------------------------
     ! This subroutine is called following dynamics. After growth has been updated
     ! there needs to be a re-assesment of the how much liquid water is bound in the
     ! plants.  This value is necessary for water balancing in the HLM.
     ! ----------------------------------------------------------------------------------

     use EDTypesMod, only : AREA

     ! Arguments
     integer, intent(in)                       :: nsites
     type(ed_site_type), intent(inout), target :: sites(nsites)
     type(bc_out_type), intent(inout)          :: bc_out(nsites)

     write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
     call endrun(msg=errMsg(sourcefile, __LINE__))

  end subroutine UpdateH2OVeg

  ! =====================================================================================

  subroutine updateSizeDepRhizHydProps(currentSite, bc_in )
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:
    use FatesConstantsMod    , only : pi_const
    use EDTypesMod           , only : AREA
    
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(bc_in_type)       , intent(in) :: bc_in

    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))

 end subroutine updateSizeDepRhizHydProps
  
 ! =================================================================================
 
 subroutine updateSizeDepRhizHydStates(currentSite, bc_in)
    !
    ! !DESCRIPTION: Updates size of 'representative' rhizosphere -- node radii, volumes.
    ! As fine root biomass (and thus absorbing root length) increases, this characteristic
    ! rhizosphere shrinks even though the total volume of soil tapped by fine roots remains
    ! the same.  
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type), intent(in)              :: bc_in
    !
    write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
    call endrun(msg=errMsg(sourcefile, __LINE__))

end subroutine updateSizeDepRhizHydStates

   ! ====================================================================================

   subroutine BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)

     ! Arguments
     integer,intent(in)                      :: nsites
     type(ed_site_type),intent(inout),target :: sites(nsites)
     type(bc_out_type),intent(inout)         :: bc_out(nsites)

     write(fates_log(),*) 'FATES Plant Hydraulics is still under development, ending run.'
     call endrun(msg=errMsg(sourcefile, __LINE__))

     
   end subroutine BTranForHLMDiagnosticsFromCohortHydr


   ! ====================================================================================


 

end module FatesPlantHydraulicsMod
