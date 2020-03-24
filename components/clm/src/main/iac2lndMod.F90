module iac2lndMod

  !------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle coupled data from iac for use in clm
  ! Iac is on the same grid as clm.
  ! !USES:
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : numpft
  use clm_varctl     , only : iulog
  use abortutils     , only : endrun
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp 
  use VegetationType , only : veg_pp

  ! !PUBLIC TPES:
  implicit none
  private
  save

  ! iac -> land structure
  ! Dimensioned by (ngrid,npft)
  type, public :: iac2lnd_type
     real(r8), pointer :: landuse(:,:) => null()

     contains

     procedure, public :: Init
     procedure, public :: Restart
     !procedure, public :: update_glc2lnd

  end type iac2lnd_type
   
  contains

  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION
    ! Allocate and initialize iac variables used by land.  Currently
    ! only coupling landuse percentage.
    ! We need a different landuse field for each pft, since we can only
    ! couple on the grid, landuse is actually an array of grid
    ! variables, dinemsioned by (grid,pft)
    !
    ! !ARGUMENTS:
    class(iac2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg
    integer :: p

    begg = bounds%begg; endg = bounds%endg

    ! Add in 0 as bare ground pft
    !allocate(this%landuse(begg:endg,0:npft))
  end subroutine Init

  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write iac2lnd information to/from restart file.
    ! Not yet implemented
    !
    use ncdio_pio , only : ncd_double 
    use pio       , only : file_desc_t

    class(iac2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*),  intent(in)    :: flag

  end subroutine Restart

  !------------------------------------------------------
  subroutine update_iac2lnd(this, bounds)
    !
    ! !DESCRIPTION:
    ! Extract into clm variables from iac coupled inputs
    ! 
    ! !USES:
    !use gcam_var_mod, only: gcam_active
    !
    ! !ARGUMENTS:
    class(iac2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    ! !LOCAL VARIABLES:
    integer :: g, c, p, pft  ! grid and pft indices
    integer :: begg, endg
    integer :: begp, endp

    character(len=*), parameter :: subname = 'update_iac2lnd'

    begg = bounds%begg; endg = bounds%endg
    begp = bounds%begp; endg = bounds%endp

#if 0
    if (gcam_active) then 
       ! The idea is to convert (ngrid,pft) to a patch-dimensioned 1D array
       ! So, loop over patches, and extract pft and g, and copy over.
       do p = begp,endp
          g=veg_pp%gridcell(p)
          pft=veg_pp%itype(p)
          
          !landuse_patch(p) = iac2lnd_vars%landuse(g,pft)

          ! Is it as simple as this?  Do we need to scale the other
          ! weights, too?
          veg_pp%wtgcell(p) = iac2lnd_vars%landuse(g,pft)

          ! Assuming we do, follow compute_higher_order_weights, line
          ! 247 in subgridWeightsMod.F90: compute wtcol first using
          ! col_pp%wtgcell, then use that to set other weights via col_pp.
          c = veg_pp%column(p)

          ! If column weights are zero, then all bets are off
          if (col_pp%wtgcell(c) /= 0) then 
             veg_pp%wtcol(p) = veg_pp%wtgcell(p)/col_pp%wtgcell(c)

             veg_pp%wttopounit(p) = veg_pp%wtcol(p) * col_pp%wttopounit(c)
             veg_pp%wtlunit(p)    = veg_pp%wtcol(p) * col_pp%wtlunit(c)             
          endif
       end do

       ! Now, send this where it needs to go.
    endif
#endif

  end subroutine update_iac2lnd

end module iac2lndMod
