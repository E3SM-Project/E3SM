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
  use clm_varpar     , only : numpft, numharvest, maxpatch_pft
  use clm_varctl     , only : iulog
  use abortutils     , only : endrun
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp 
  use VegetationType , only : veg_pp
  use dynHarvestMod  , only : harvest, do_cn_harvest

  ! !PUBLIC TPES:
  implicit none
  private
  save

  ! iac -> land structure
  ! Dimensioned by (ngrid,numpft)
  type, public :: iac2lnd_type
     real(r8), pointer :: pct_pft(:,:) => null()
     real(r8), pointer :: pct_pft_prev(:,:) => null()
     real(r8), pointer :: harvest_frac(:,:) => null()

   contains

     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: update_iac2lnd

  end type iac2lnd_type

contains

  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION
    ! Allocate and initialize iac variables used by land.
    ! We need a different landuse field for each pft, since we can only
    ! couple on the grid, landuse is actually an array of grid
    ! variables, dimensioned by (grid,pft)
    ! harvest dimensioned by (grid,harvest)
    !
    ! !ARGUMENTS:
    class(iac2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg
    integer :: p
    integer :: ier        ! error code
    character(len=100) :: errstr

    begg = bounds%begg; endg = bounds%endg

    ! set the cn harveset flag
    do_cn_harvest = .true.

    ! Add in 0 as bare ground pft
    ! note that indices here start at 0
    ! numpft does not include bare ground
    ! nharvest is the total harvest fields
    allocate(this%pct_pft(begg:endg,0:numpft))
    allocate(this%pct_pft_prev(begg:endg,0:numpft))
    allocate(this%harvest_frac(begg:endg,0:(numharvest-1)))

    ! allocate the harvest array; 1d cuz they are added
    allocate(harvest(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for harvest'// &
                   errMsg(__FILE__, __LINE__))
    end if

    ! check that number of pfts is correct
    if ( maxpatch_pft /= numpft+1 )then
       errstr = 'maxpatch_pft does NOT equal numpft+1 - invalid for dyn pft' 
       call endrun(msg=errstr//errMsg(__FILE__, __LINE__) )
    end if

  end subroutine Init

  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write iac2lnd information to/from restart file.
    ! Not yet implemented

    ! !USES:
    use ncdio_pio , only : ncd_double, file_desc_t
    use decompMod , only : bounds_type
    !
    ! !ARGUMENTS:
    class(iac2lnd_type) , intent(inout) :: this
    type(bounds_type)   , intent(in)    :: bounds 
    type(file_desc_t)   , intent(inout) :: ncid ! netcdf id
    character(len=*)    , intent(in)    :: flag ! 'read' or 'write'

  end subroutine Restart

  !------------------------------------------------------
  subroutine update_iac2lnd(this, bounds)
    !
    ! !DESCRIPTION:
    ! Extract into clm variables from iac coupled inputs
    ! 
    ! !USES:
    use clm_time_manager, only : get_curr_yearfrac
    use landunit_varcon , only : istsoil
    use clm_varctl, only: iac_active
    !
    ! !ARGUMENTS:
    class(iac2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    ! !LOCAL VARIABLES:
    integer :: g, p, c, h, l, pft  ! grid, other  indices
    integer :: begg, endg
    integer :: begp, endp
    real(r8) :: wt1    ! weight of time1 (prev time point)

    character(len=*), parameter :: subname = 'update_iac2lnd'

    begg = bounds%begg; endg = bounds%endg
    begp = bounds%begp; endp = bounds%endp

    if (iac_active) then  ! this is also checked in order to this function
       ! The idea is to convert (ngrid,pft) to a patch-dimensioned 1D array
       ! So, loop over patches, and extract pft and g, and copy over.
       ! currently, each patch is a pft
       ! this is from dynpft_interp
       
       ! the the weight of prev time point
       wt1 = 1.0_r8 - get_curr_yearfrac()

       do p = begp,endp
          g=veg_pp%gridcell(p)
          pft=veg_pp%itype(p)
          l = veg_pp%landunit(p)
         
          ! Note that we only deal with the istsoil landunit here, NOT the
          ! istcrop landunit
          ! (if there is one)
          ! (However, currently [as of 5-9-13] the code won't let you run with
          ! transient
          ! Patches combined with create_crop_landunit anyway, so it's a moot
          ! point.)
          if (lun_pp%itype(l) == istsoil) then
             ! interpolate between the yearly data; from dynvartimeinterp
             ! Note that the following assignment assumes that all Patches share a
             ! single column
             ! iac2lnd indices start at 0 to match pft ids
             veg_pp%wtcol(p) = this%pct_pft(g,pft) + &
                             wt1*(this%pct_pft_prev(g,pft) - this%pct_pft(g,pft))
          end if

       end do

       ! sum the harvest data into one field
       harvest(bounds%begg:bounds%endg) = 0._r8
       do h=0,(numharvest-1)
          harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + &
                               this%harvest_frac(bounds%begg:bounds%endg,h)
       end do

    endif
  end subroutine update_iac2lnd
end module iac2lndMod

