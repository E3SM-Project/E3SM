module map_lnd2rof_irrig_mod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module contains routines for mapping the irrigation field from the LND grid onto
  ! the ROF grid.
  !
  ! These routines could go in prep_rof_mod, but are separated into their own module for
  ! the sake of (1) testability: this module has fewer dependencies than prep_rof_mod;
  ! and (2) symmetry with the lnd2glc and glc2lnd custom mapping routines, which also
  ! have their own modules.

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mct_mod
  use seq_map_type_mod, only : seq_map
  use seq_map_mod, only : seq_map_map
  use shr_log_mod, only : errMsg => shr_log_errMsg

  implicit none
  private

  ! ------------------------------------------------------------------------
  ! Public interfaces
  ! ------------------------------------------------------------------------

  public :: map_lnd2rof_irrig  ! map irrigation from lnd -> rof grid

  ! ------------------------------------------------------------------------
  ! Private interfaces
  ! ------------------------------------------------------------------------

  private :: map_rof2lnd_volr  ! map volr from rof -> lnd grid

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  subroutine map_lnd2rof_irrig(l2r_l, r2x_r, irrig_flux_field, &
       avwts_s, avwtsfld_s, mapper_Fl2r, mapper_Fr2l, l2r_r)
    !---------------------------------------------------------------
    ! Description
    ! Do custom mapping for the irrigation flux, from land -> rof.
    !
    ! The basic idea is that we want to pull irrigation out of ROF cells proportionally to
    ! the river volume (volr) in each cell. This is important in cases where the various
    ! ROF cells overlapping a CLM cell have very different volr: If we didn't do this
    ! volr-normalized remapping, we'd try to extract the same amount of water from each
    ! of the ROF cells, which would be more likely to have withdrawals exceeding
    ! available volr.
    !
    ! (Both RTM and MOSART have code to handle excess withdrawals, by pulling the excess
    ! directly out of the ocean, but we'd like to avoid resorting to this as much as
    ! possible.)
    !
    ! This mapping works by:
    !
    ! (1) Normalizing the land's irrigation flux by volr
    !
    ! (2) Mapping this volr-normalized flux to the rof grid
    !
    ! (3) Converting the mapped, volr-normalized flux back to a normal
    !     (non-volr-normalized) flux on the rof grid.
    !
    ! This assumes that the following fields are contained in the attribute vector
    ! arguments:
    !
    ! - l2r_l: field given by irrig_flux_field (read)
    ! - l2r_r: field given by irrig_flux_field (set)
    ! - r2x_r: 'Flrr_volrmch' (read)
    !
    ! Arguments
    type(mct_aVect)  , intent(in)    :: l2r_l            ! lnd -> rof fields on the land grid
    type(mct_aVect)  , intent(in)    :: r2x_r            ! rof -> cpl fields on the rof grid
    character(len=*) , intent(in)    :: irrig_flux_field ! name of irrigation field to remap
    type(mct_aVect)  , intent(in)    :: avwts_s          ! attr vect for source weighting
    character(len=*) , intent(in)    :: avwtsfld_s       ! field in avwts_s to use
    type(seq_map)    , intent(inout) :: mapper_Fl2r      ! flux mapper for mapping lnd -> rof
    type(seq_map)    , intent(inout) :: mapper_Fr2l      ! flux mapper for mapping rof -> lnd
    type(mct_aVect)  , intent(inout) :: l2r_r            ! lnd -> rof fields on the rof grid
    !
    ! Local variables
    integer :: r, l
    integer :: lsize_l  ! number of land points
    integer :: lsize_r  ! number of rof points
    type(mct_avect) :: irrig_l_av  ! temporary attribute vector holding irrigation fluxes on the land grid
    type(mct_avect) :: irrig_r_av  ! temporary attribute vector holding irrigation fluxes on the rof grid

    ! The following need to be pointers to satisfy the MCT interface:
    real(r8), pointer :: volr_r(:)             ! river volume on the rof grid
    real(r8), pointer :: volr_l(:)             ! river volume on the land grid
    real(r8), pointer :: irrig_flux_l(:)       ! irrigation flux on the land grid [kg m-2 s-1]
    real(r8), pointer :: irrig_flux_r(:)       ! irrigation flux on the rof grid [kg m-2 s-1]
    real(r8), pointer :: irrig_normalized_l(:) ! irrigation normalized by volr, land grid
    real(r8), pointer :: irrig_normalized_r(:) ! irrigation normalized by volr, rof grid
    real(r8), pointer :: irrig_volr0_l(:)      ! irrigation where volr <= 0, land grid
    real(r8), pointer :: irrig_volr0_r(:)      ! irrigation where volr <= 0, rof grid

    character(len=*), parameter :: volr_field             = 'Flrr_volrmch'
    character(len=*), parameter :: irrig_normalized_field = 'Flrl_irrig_normalized'
    character(len=*), parameter :: irrig_volr0_field      = 'Flrl_irrig_volr0'
    character(len=*), parameter :: fields_to_remap = &
         irrig_normalized_field // ':' // irrig_volr0_field
    !---------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Determine attribute vector sizes
    ! ------------------------------------------------------------------------

    lsize_l = mct_aVect_lsize(l2r_l)
    lsize_r = mct_aVect_lsize(l2r_r)

    ! ------------------------------------------------------------------------
    ! Extract the necessary fields from attribute vectors
    ! ------------------------------------------------------------------------

    allocate(irrig_flux_l(lsize_l))
    call mct_aVect_exportRattr(l2r_l, irrig_flux_field, irrig_flux_l)

    allocate(volr_r(lsize_r))
    call mct_aVect_exportRattr(r2x_r, volr_field, volr_r)

    ! ------------------------------------------------------------------------
    ! Adjust volr_r, and map it to the land grid
    ! ------------------------------------------------------------------------

    ! Treat any rof point with volr < 0 as if it had volr = 0. Negative volr values can
    ! arise in RTM. This fix is needed to avoid mapping negative irrigation to those
    ! cells: while conservative, this would be unphysical (it would mean that irrigation
    ! actually adds water to those cells).
    do r = 1, lsize_r
       if (volr_r(r) < 0._r8) then
          volr_r(r) = 0._r8
       end if
    end do

    allocate(volr_l(lsize_l))
    call map_rof2lnd_volr(volr_r, mapper_Fr2l, volr_l)

    ! ------------------------------------------------------------------------
    ! Determine irrigation normalized by volr
    !
    ! In order to avoid possible divide by 0, as well as to handle non-sensical negative
    ! volr on the land grid, we divide the land's irrigation flux into two separate flux
    ! components: a component where we have positive volr on the land grid (put in
    ! irrig_normalized_l, which is mapped using volr-normalization) and a component where
    ! we have zero or negative volr on the land grid (put in irrig_volr0_l, which is
    ! mapped as a standard flux). We then remap both of these components to the rof grid,
    ! and then finally add the two components to determine the total irrigation flux on
    ! the rof grid.
    ! ------------------------------------------------------------------------

    allocate(irrig_normalized_l(lsize_l))
    allocate(irrig_volr0_l(lsize_l))
    do l = 1, lsize_l
       if (volr_l(l) > 0._r8) then
          irrig_normalized_l(l) = irrig_flux_l(l) / volr_l(l)
          irrig_volr0_l(l)      = 0._r8
       else
          irrig_normalized_l(l) = 0._r8
          irrig_volr0_l(l)      = irrig_flux_l(l)
       end if
    end do

    ! ------------------------------------------------------------------------
    ! Map irrigation
    ! ------------------------------------------------------------------------

    call mct_aVect_init(irrig_l_av, rList = fields_to_remap, lsize = lsize_l)
    call mct_aVect_importRattr(irrig_l_av, irrig_normalized_field, irrig_normalized_l)
    call mct_aVect_importRattr(irrig_l_av, irrig_volr0_field, irrig_volr0_l)
    call mct_aVect_init(irrig_r_av, rList = fields_to_remap, lsize = lsize_r)

    ! This mapping uses the same options (such as avwts) as is used for mapping all other
    ! fields in prep_rof_calc_l2r_rx
    call seq_map_map(mapper = mapper_Fl2r, &
         av_s = irrig_l_av, &
         av_d = irrig_r_av, &
         fldlist = fields_to_remap, &
         norm = .true., &
         avwts_s = avwts_s, &
         avwtsfld_s = avwtsfld_s)

    allocate(irrig_normalized_r(lsize_r))
    allocate(irrig_volr0_r(lsize_r))
    call mct_aVect_exportRattr(irrig_r_av, irrig_normalized_field, irrig_normalized_r)
    call mct_aVect_exportRattr(irrig_r_av, irrig_volr0_field, irrig_volr0_r)

    ! ------------------------------------------------------------------------
    ! Convert to a total irrigation flux on the ROF grid, and put this in the l2r_rx
    ! attribute vector
    ! ------------------------------------------------------------------------

    allocate(irrig_flux_r(lsize_r))
    do r = 1, lsize_r
       irrig_flux_r(r) = (irrig_normalized_r(r) * volr_r(r)) + irrig_volr0_r(r)
    end do

    call mct_aVect_importRattr(l2r_r, irrig_flux_field, irrig_flux_r)

    ! ------------------------------------------------------------------------
    ! Clean up
    ! ------------------------------------------------------------------------

    deallocate(volr_r)
    deallocate(volr_l)
    deallocate(irrig_flux_l)
    deallocate(irrig_flux_r)
    deallocate(irrig_normalized_l)
    deallocate(irrig_normalized_r)
    deallocate(irrig_volr0_l)
    deallocate(irrig_volr0_r)
    call mct_aVect_clean(irrig_l_av)
    call mct_aVect_clean(irrig_r_av)

  end subroutine map_lnd2rof_irrig

  subroutine map_rof2lnd_volr(volr_r, mapper_Fr2l, volr_l)
    !---------------------------------------------------------------
    ! Description
    ! Map volr from the rof grid to the lnd grid.
    !
    ! This is needed for the volr-normalization that is done in map_lnd2rof_irrig.
    !
    ! Note that this mapping is also done in the course of mapping all rof -> lnd fields
    ! in prep_lnd_calc_r2x_lx. However, we do this mapping ourselves here for two reasons:
    !
    ! (1) For the sake of this normalization, we change all volr < 0 values to 0; this is
    !     not done for the standard rof -> lnd mapping.
    !
    ! (2) It's possible that the driver sequencing would be changed such that this rof ->
    !     lnd mapping happens before the lnd -> rof mapping. If that happened, then volr_l
    !     (i.e., volr that has been mapped to the land grid by prep_lnd_calc_r2x_lx) would
    !     be inconsistent with volr_r, which would be a Bad Thing for the
    !     volr-normalizated mapping (this mapping would no longer be conservative). So we
    !     do the rof -> lnd remapping here to ensure we have a volr_l that is consistent
    !     with volr_r.
    !
    ! The pointer arguments to this routine should already be allocated to be the
    ! appropriate size.
    !
    ! Arguments
    real(r8), pointer, intent(in)    :: volr_r(:)   ! river volume on the rof grid (input)
    type(seq_map)    , intent(inout) :: mapper_Fr2l ! flux mapper for mapping rof -> lnd
    real(r8), pointer, intent(inout) :: volr_l(:)   ! river volume on the lnd grid (output) (technically intent(in) since intent gives the association status of a pointer, but given as intent(inout) to avoid confusion, since its data are modified)
    !
    ! Local variables
    integer :: lsize_r  ! number of rof points
    integer :: lsize_l  ! number of lnd points
    type(mct_avect) :: volr_r_av   ! temporary attribute vector holding volr on the rof grid
    type(mct_avect) :: volr_l_av   ! temporary attribute vector holding volr on the land grid

    ! This volr field name does not need to agree with the volr field name used in the
    ! 'real' attribute vectors
    character(len=*), parameter :: volr_field = 'volr'
    !---------------------------------------------------------------

    SHR_ASSERT_FL(associated(volr_r), sourcefile, __LINE__)
    SHR_ASSERT_FL(associated(volr_l), sourcefile, __LINE__)

    lsize_r = size(volr_r)
    lsize_l = size(volr_l)

    call mct_aVect_init(volr_r_av, rList = volr_field, lsize = lsize_r)
    call mct_aVect_importRattr(volr_r_av, volr_field, volr_r)
    call mct_aVect_init(volr_l_av, rList = volr_field, lsize = lsize_l)

    ! This mapping uses the same options as the standard rof -> lnd mapping done in
    ! prep_lnd_calc_r2x_lx. If that mapping ever changed (e.g., introducing an avwts_s
    ! argument), then it's *possible* that we'd want this mapping to change, too.
    call seq_map_map(mapper = mapper_Fr2l, &
         av_s = volr_r_av, &
         av_d = volr_l_av, &
         fldlist = volr_field, &
         norm = .true.)

    call mct_aVect_exportRattr(volr_l_av, volr_field, volr_l)

    call mct_aVect_clean(volr_r_av)
    call mct_aVect_clean(volr_l_av)

  end subroutine map_rof2lnd_volr

end module map_lnd2rof_irrig_mod
