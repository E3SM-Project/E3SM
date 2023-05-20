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
  use elm_varpar     , only : numpft, numharvest, maxpatch_pft
  use elm_varctl     , only : iulog
  use abortutils     , only : endrun
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp 
  use VegetationType , only : veg_pp
  use dynHarvestMod  , only : harvest_rates, do_cn_harvest

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS
  public :: iac_rpointer_write

  ! iac -> land structure
  ! Dimensioned by (ngrid,numpft)
  ! these values are fraction of actual grid cell (not fraction of land)
  type, public :: iac2lnd_type
     real(r8), pointer :: frac_pft(:,:) => null()
     real(r8), pointer :: frac_pft_prev(:,:) => null()
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

    ! Add in 0 as bare ground pft
    ! note that indices here start at 0
    ! numpft does not include bare ground
    ! nharvest is the total harvest fields
    allocate(this%frac_pft(begg:endg,0:numpft))
    allocate(this%frac_pft_prev(begg:endg,0:numpft))
    allocate(this%harvest_frac(begg:endg,0:(numharvest-1)))

    ! initialize the arrays
    this%frac_pft(:,:) = 1.0_r8
    this%frac_pft_prev(:,:) = 1.0_r8
    this%harvest_frac(:,:) = 1.0_r8
    harvest_rates(:,:) = 0.0_r8

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
    use elm_varctl, only : iac_active, iulog
    use netcdf   !avd - this is for a diagnostic file
    use domainMod,   only : ldomain
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
    real(r8) :: temp, temp_prev

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
          c = veg_pp%column(p)        
 
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

             ! these values are fraction of actual grid cell (not fraction of
             !    land), so need to be converted
             ! use the domain landfrac to get veg_pp%wtgcell_iac, which is the
             !    pft fraction of land (lnd model assumes all land cells have
             !    landfrac==1)
             ! since this happens before update_landunit_weights and
             !    compute_higher_order_weights,
             !    store this in a special variable and then do the final calcs
             !    after the comp_h_o_w function and before the checks
             ! use a new function to calculate 
             !    the veg_pp fractions, including veg_pp%wtcol, which is the
             !    pft fraction of the column (i.e., of the veg land unit)
             !    it will renormalize veg_pp%wtcol sum to 1 because landfracs are
             !    consistent across components and grids, and then recalc the
             !    other fractions

             ! ldomain%mask matches ldomain%frac - these are the only land
             !    cells processed, so land mask should be one if frac != 0
             ! so can just set pfts to zero outside of the active cells,
             !    actually, this zero condition shouldn't exist 

             if (ldomain%frac(g) .eq. 0._r8) then
                   veg_pp%wtgcell_iac(p) = 0._r8
             else
                temp = this%frac_pft(g,pft) / (ldomain%frac(g) * ldomain%mask(g))
                temp_prev = this%frac_pft_prev(g,pft) / &
                             (ldomain%frac(g) * ldomain%mask(g)) 
                veg_pp%wtgcell_iac(p) = temp + wt1*(temp_prev - temp)
             end if
          end if
       end do

       ! these are also fractions of actual grid cell
       !    but are for this year so can use the current/prior col_pp%wtgcell for veg in
       !    this cell to get the final fractions of col or veg land unit
       ! land unit type 1 is vegetated land unit
       ! just need to make sure that the values are <= 1

       harvest_rates(:,begg:endg) = 0._r8
       do g = begg,endg
         do c = bounds%begc, bounds%endc
            if (col_pp%itype(c) .eq. 1 .and. col_pp%gridcell(c) .eq. g) then
                ! sum the harvest data into one field 
                do h=0,(numharvest-1)
                   if (.not.(col_pp%wtgcell(c) .eq. 0._r8)) then
                      ! TRS - now using harvest-distinct harvest_rates var
                      !harvest(g) = harvest(g) + &
                      !           this%harvest_frac(g,h) / col_pp%wtgcell(c)
                      ! Note harvest_rates is one-offset in h 
                      ! or not?
                      harvest_rates(h+1,g) = this%harvest_frac(g,h) / col_pp%wtgcell(c)
                   end if
                end do
             end if
          end do
       end do

! avd - don't do this to see if it speeds up; besides, it isn't useful in
! parallel
if (.false.) then
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

       ierr = nf90_put_var(ncid,pft_varid,this%frac_pft(begg:endg,:))

       ierr = nf90_put_var(ncid,pft_prev_varid,this%frac_pft_prev(begg:endg,:))

       ierr = nf90_put_var(ncid,harv_varid,this%harvest_frac(begg:endg,:))

       ierr = nf90_put_var(ncid,cid_varid,cell_ids)

       ierr = nf90_close(ncid)
end if

    endif
  end subroutine update_iac2lnd


  !---------------------------------------------------------------------
  subroutine iac_rpointer_write(rdate)
  !
  ! !DESCRIPTION:
  !
  ! Write the gcam2glm and glm rpointer files so that they point to the correct
  !    year restart files that correspond with the other E3SM restart files
  !
  ! !USES:
  use shr_file_mod, only: shr_file_getunit, shr_file_freeunit
  !
  ! !ARGUMENTS:
  implicit none
  character(len=*),intent(in) :: rdate       ! restart file time stamp for name
  !
  ! !LOCAL VARIABLES:
  integer       :: year         ! 4 digit year extracted from the front of rdate
  character(256) :: filename
  integer       :: iun
  character(len=*),parameter :: gcam2glm_restfile = 'gcam2glm_restart.'
  character(len=*),parameter :: gcam2glm_rpointer = 'rpointer.gcam2glm'
  character(len=*),parameter :: glm_restfile = 'output.glm.restart.state.'
  character(len=*),parameter :: glm_rpointer = 'rpointer.glm'

  !-----------------------------------------------------------------------
 
  ! set the names here cuz the iac mods are not available
  ! use the year from rdate
 
  ! write the gcam2glm rpointer file  

  write(filename,'(a)') trim(gcam2glm_restfile)//'r.'//rdate(1:4)//'.nc'

  iun = shr_file_getunit()
  open(iun,file=trim(gcam2glm_rpointer),form='formatted')
  write(iun,'(a)') trim(filename)
  close(iun)
  call shr_file_freeunit(iun)

 ! write the glm rpointer file  

  write(filename,'(a)') trim(glm_restfile)//rdate(1:4)//'.nc'

  iun = shr_file_getunit()
  open(iun,file=trim(glm_rpointer),form='formatted')
  write(iun,'(a)') trim(filename)
  close(iun)
  call shr_file_freeunit(iun)

  end subroutine iac_rpointer_write


end module iac2lndMod

