module ComputeSeedMod
  !-----------------------------------------------------------------------
  ! Module to compute seed amounts for new patch areas
  !
  ! !USES:
#include "shr_assert.h"

  use shr_kind_mod             , only : r8 => shr_kind_r8
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use decompMod                , only : bounds_type
  use elm_varcon               , only : c3_r2, c4_r2, c14ratio
  use elm_varctl               , only : iulog
  use abortutils               , only : endrun
  use SpeciesMod               , only : CN_SPECIES_C12, CN_SPECIES_C13, CN_SPECIES_C14
  use SpeciesMod               , only : CN_SPECIES_N, CN_SPECIES_P
  use VegetationPropertiesType , only : veg_vp
  use VegetationType           , only : veg_pp
  use LandunitType             , only : lun_pp
  !
  ! !PUBLIC ROUTINES:
  implicit none
  private

  public :: ComputeSeedAmounts

  ! !PRIVATE ROUTINES:

  private :: SpeciesTypeMultiplier
  private :: LeafProportions  ! compute leaf proportions (leaf, storage and xfer)

  ! !PRIVATE DATA:

  integer  , parameter :: COMPONENT_LEAF       = 1
  integer  , parameter :: COMPONENT_DEADWOOD   = 2
  integer  , parameter :: COMPONENT_SEED       = 3
  real(r8) , parameter :: leafc_seed_param     = 1._r8
  real(r8) , parameter :: deadstemc_seed_param = 0.1_r8
  !$acc declare copyin (COMPONENT_LEAF      )
  !$acc declare copyin (COMPONENT_DEADWOOD  )
  !$acc declare copyin (COMPONENT_SEED      )
  !$acc declare copyin (leafc_seed_param    )
  !$acc declare copyin (deadstemc_seed_param)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains


  !-----------------------------------------------------------------------
  subroutine ComputeSeedAmounts(p,                          &
       species,                                                        &
       leaf_patch, leaf_storage_patch, leaf_xfer_patch,                &
       compute_here_patch, ignore_current_state_patch,                 &
       seed_leaf_patch, seed_leaf_storage_patch, seed_leaf_xfer_patch, &
       seed_deadstem_patch, pool_seed_param, pool_seed_patch)
    !
    ! !DESCRIPTION:
    ! Compute seed amounts for patches that increase in area, for various variables, for
    ! the given species (c12, c13, c14 or n).
    !
    ! The output variables are only set for patches inside the filter, where
    ! compute_here_patch is true; for other patches, they remain at their original values.
    !
    ! Note that, regardless of the species, leafc_seed and deadstemc_seed are specified
    ! in terms of gC/m2; these amounts are converted to the amount of the given species
    ! here.
    !
    ! !USES:
    !$acc routine seq
    use pftvarcon       , only : noveg
    !
    ! !ARGUMENTS:
    integer, value    , intent(in)     :: p
    integer           , intent(in)     :: species                ! which C/N species we're operating on; should be one of the values in SpeciesMod
    real(r8)          , intent(in)     :: leaf_patch             ! current leaf C or N content (g/m2)
    real(r8)          , intent(in)     :: leaf_storage_patch     ! current leaf C or N storage content (g/m2)
    real(r8)          , intent(in)     :: leaf_xfer_patch        ! current leaf C or N xfer content (g/m2)
                                                                                  ! whether to compute outputs for each patch
    logical           , intent(in)     :: compute_here_patch
                                                                                  ! If ignore_current_state is true, then use default leaf proportions rather than
                                                                                  ! proportions based on current state.
    logical           , intent(in)     :: ignore_current_state_patch
    real(r8)          , intent(inout)  :: seed_leaf_patch          ! seed amount for leaf itself for this species (g/m2)
    real(r8)          , intent(inout)  :: seed_leaf_storage_patch  ! seed amount for leaf storage for this species (g/m2)
    real(r8)          , intent(inout)  :: seed_leaf_xfer_patch     ! seed amount for leaf xfer for this species (g/m2)
    real(r8)          , intent(inout)  :: seed_deadstem_patch      ! seed amount for deadstem for this species (g/m2)
    real(r8), optional, intent(in)     :: pool_seed_param
    real(r8), optional, intent(inout)  :: pool_seed_patch
    !
    ! !LOCAL VARIABLES:
    real(r8) :: my_leaf_seed
    real(r8) :: my_deadstem_seed
    real(r8) :: my_pool_seed
    integer  :: pft_type
    real(r8) :: pleaf
    real(r8) :: pstor
    real(r8) :: pxfer
    !-----------------------------------------------------------------------
#ifndef _OPENACC
    if (present(pool_seed_patch)) then
       if (.not. present(pool_seed_param)) then
          call endrun( ': pool_seed_patch can only be provided with pool_seed_param')
       end if
    end if
#endif

    if (compute_here_patch) then
          my_leaf_seed = 0._r8
          my_deadstem_seed = 0._r8
          my_pool_seed = 0._r8
          pft_type = veg_pp%itype(p)

          call LeafProportions( &
               ignore_current_state = ignore_current_state_patch , &
               pft_type = pft_type, &
               leaf = leaf_patch , &
               leaf_storage = leaf_storage_patch , &
               leaf_xfer = leaf_xfer_patch , &
               pleaf = pleaf, &
               pstorage = pstor, &
               pxfer = pxfer)

          if (pft_type /= noveg) then
             my_leaf_seed = leafc_seed_param * &
                  SpeciesTypeMultiplier(species, pft_type, COMPONENT_LEAF)
             if (veg_vp%woody(pft_type) >= 1.0_r8) then
                my_deadstem_seed = deadstemc_seed_param * &
                     SpeciesTypeMultiplier(species, pft_type, COMPONENT_DEADWOOD)
             end if
             if (present(pool_seed_param)) then
                my_pool_seed = pool_seed_param     * &
                     SpeciesTypeMultiplier(species, pft_type, COMPONENT_SEED)
             end if
          end if

          seed_leaf_patch          = my_leaf_seed * pleaf
          seed_leaf_storage_patch  = my_leaf_seed * pstor
          seed_leaf_xfer_patch     = my_leaf_seed * pxfer
          seed_deadstem_patch      = my_deadstem_seed
          if (present(pool_seed_param)) then
             pool_seed_patch  = my_pool_seed
          end if
    else
          seed_leaf_patch          = 0._r8
          seed_leaf_storage_patch  = 0._r8
          seed_leaf_xfer_patch     = 0._r8
          seed_deadstem_patch      = 0._r8
    end if

  end subroutine ComputeSeedAmounts

  !-----------------------------------------------------------------------
  subroutine LeafProportions(pft_type, ignore_current_state, &
       leaf, leaf_storage, leaf_xfer, &
       pleaf, pstorage, pxfer)
    !
    ! !DESCRIPTION:
    ! Compute leaf proportions (leaf, storage and xfer)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
      !$acc routine seq
    integer , intent(in)  :: pft_type
    logical , intent(in)  :: ignore_current_state
    real(r8), intent(in)  :: leaf         ! g/m2 leaf C or N
    real(r8), intent(in)  :: leaf_storage ! g/m2 leaf C or N storage
    real(r8), intent(in)  :: leaf_xfer    ! g/m2 leaf C or N transfer

    real(r8), intent(out) :: pleaf        ! proportion in leaf itself
    real(r8), intent(out) :: pstorage     ! proportion in leaf storage
    real(r8), intent(out) :: pxfer        ! proportion in leaf xfer
    !
    ! !LOCAL VARIABLES:
    real(r8) :: tot_leaf

    character(len=*), parameter :: subname = 'LeafProportions'
    !-----------------------------------------------------------------------

    tot_leaf = leaf + leaf_storage + leaf_xfer
    pleaf    = 0._r8
    pstorage = 0._r8
    pxfer    = 0._r8

    if (tot_leaf == 0._r8 .or. ignore_current_state) then
       if (veg_vp%evergreen(pft_type) == 1._r8) then
          pleaf    = 1._r8
       else
          pstorage = 1._r8
       end if
    else
       pleaf    = leaf        /tot_leaf
       pstorage = leaf_storage/tot_leaf
       pxfer    = leaf_xfer   /tot_leaf
    end if

  end subroutine LeafProportions

  !-----------------------------------------------------------------------
  function SpeciesTypeMultiplier(species, pft_type, component) result(multiplier)
    !
    ! !DESCRIPTION:
    ! Returns a multiplier based on the species type. This multiplier is
    ! meant to be applied to some state variable expressed in terms of g C, translating
    ! this value into an appropriate value for c13, c14, n or p.
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    real(r8)            :: multiplier ! function result
    integer, intent(in) :: species    ! which C/N species we're operating on; should be one of the values in SpeciesMod
    integer, intent(in) :: pft_type
    integer, intent(in) :: component  ! which plant component; should be one of the COMPONENT_* parameters defined in this module
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'SpeciesTypeMultiplier' 
    !-----------------------------------------------------------------------
    select case (species)
    case (CN_SPECIES_C12)
       multiplier = 1._r8

    case (CN_SPECIES_C13)
       if (veg_vp%c3psn(pft_type) == 1._r8) then
          multiplier = c3_r2
       else
          multiplier = c4_r2
       end if

    case (CN_SPECIES_C14)
       ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
       multiplier = c14ratio

    case (CN_SPECIES_N)
       select case (component)
       case (COMPONENT_LEAF)
          multiplier = 1._r8 / veg_vp%leafcn(pft_type)
       case (COMPONENT_DEADWOOD)
          multiplier = 1._r8 / veg_vp%deadwdcn(pft_type)
       case (COMPONENT_SEED)
          if (pft_type /= 0) then
             multiplier = 1._r8
          else
             multiplier = 0._r8
          end if
       case default
#ifndef _OPENACC
          write(iulog,*) subname//' ERROR: unknown component: ', component 
          call endrun(subname//' ERROR: unknown component')
#endif 
       end select

    case (CN_SPECIES_P)
       select case (component)
       case (COMPONENT_LEAF)
          multiplier = 1._r8 / veg_vp%leafcp(pft_type)
       case (COMPONENT_DEADWOOD)
          multiplier = 1._r8 / veg_vp%deadwdcp(pft_type)
       case (COMPONENT_SEED)
          if (pft_type /= 0) then
             multiplier = 1._r8
          else
             multiplier = 0._r8
          end if
       case default
#ifndef _OPENACC
          write(iulog,*) subname//' ERROR: unknown component: ', component 
          call endrun(subname//' ERROR: unknown component')
#endif 
       end select

    case default
#ifndef _OPENACC
          write(iulog,*) subname//' ERROR: unknown species: ', species 
          call endrun(subname//' ERROR: unknown species')
#endif 
    end select

  end function SpeciesTypeMultiplier

end module ComputeSeedMod
