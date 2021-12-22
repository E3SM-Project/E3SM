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

    ! initialize the arrays
    this%pct_pft(:,:) = 1.0_r8
    this%pct_pft_prev(:,:) = 1.0_r8
    this%harvest_frac(:,:) = 1.0_r8
    harvest(:) = 0.0

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
    use clm_varctl, only: iac_active, iulog
    use netcdf   !avd - this is for a diagnostic file
    !
    ! !ARGUMENTS:
    class(iac2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    ! !LOCAL VARIABLES:
    integer :: g, p, c, h, l, pft  ! grid, other  indices
    integer :: begg, endg
    integer :: begp, endp
    real(r8) :: wt1    ! weight of time1 (prev time point)

! avd - diagnostic file
character(len=128) :: hfile
integer :: ierr, nmode, ncid, ngcells, pft_varid, pft_prev_varid
integer :: harv_varid, cid_varid
integer :: pdimid(2), hdimid(2)
integer, allocatable :: cell_ids(:)

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

!avd
!write(iulog,*) subname, 'begg,endg=', begg,endg
!write(iulog,*) subname, 'begp,endp=', begp,endp

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

             ! avd - do not actually set the land values because the mapping is
             !       incorrect such that the land model fails
             !veg_pp%wtcol(p) = this%pct_pft(g,pft) + &
             !                wt1*(this%pct_pft_prev(g,pft) - this%pct_pft(g,pft))
          end if

! avd write these to log
!if (this%pct_pft_prev(g,pft) /= 0.) then
!if (g == 1628) then 
!write(iulog,*) 'g=',g,' p=',p,' pft=',pft
!write(iulog,*) 'prev_val=',this%pct_pft_prev(g,pft)
!write(iulog,*) 'val=',this%pct_pft(g,pft)
!write(iulog,*) 'interp=', this%pct_pft(g,pft) + &
!                          wt1*(this%pct_pft_prev(g,pft) - this%pct_pft(g,pft))
!end if

       end do

       ! sum the harvest data into one field
       harvest(begg:endg) = 0._r8
       do h=0,(numharvest-1)
          harvest(begg:endg) = harvest(begg:endg) + &
                               this%harvest_frac(begg:endg,h)
       end do

       ! avd - write these values to a file
       write(iulog,*) subname, 'Writing diagnostic iac2land file'

       ngcells = endg-begg+1
       allocate(cell_ids(ngcells))
       p = 1
       do c = begg,endg
          cell_ids(p) = c
          p = p+1
       end do
       write(hfile,'(a)') './iac2lnd_update.nc'
       ierr = nf90_create(trim(hfile),nf90_clobber,ncid)
       ierr = nf90_def_dim(ncid,'ngcells',ngcells,pdimid(1))
       ierr = nf90_def_dim(ncid,'npfts',17,pdimid(2))
       ierr = nf90_def_dim(ncid,'ngcells',ngcells,hdimid(1))
       ierr = nf90_def_dim(ncid,'nharv',5,hdimid(2))

       ierr = nf90_def_var(ncid,'pftdata',NF90_DOUBLE,pdimid,pft_varid)
       ierr = nf90_def_var(ncid,'prev_pftdata',NF90_DOUBLE,pdimid,pft_prev_varid)
       ierr = nf90_def_var(ncid,'harvdata',NF90_DOUBLE,hdimid,harv_varid)
       ierr = nf90_def_var(ncid,'cell_ids',NF90_INT,pdimid(1),cid_varid)

       ierr = nf90_enddef(ncid)

       ierr = nf90_put_var(ncid,pft_varid,this%pct_pft(begg:endg,:))

       ierr = nf90_put_var(ncid,pft_prev_varid,this%pct_pft_prev(begg:endg,:))

       ierr = nf90_put_var(ncid,harv_varid,this%harvest_frac(begg:endg,:))

       ierr = nf90_put_var(ncid,cid_varid,cell_ids)

       ierr = nf90_close(ncid)

    endif
  end subroutine update_iac2lnd
end module iac2lndMod

