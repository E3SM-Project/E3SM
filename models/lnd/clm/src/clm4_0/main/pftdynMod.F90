module pftdynMod

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: pftdynMod
!
! !USES:
  use spmdMod
  use clmtype
  use decompMod   , only : get_proc_bounds
  use clm_varsur  , only : pctspec
  use clm_varpar  , only : max_pft_per_col
  use clm_varctl  , only : iulog, use_c13, use_cn, use_cndv
  use shr_sys_mod , only : shr_sys_flush
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use ncdio_pio
!
! !DESCRIPTION:
! Determine pft weights at current time using dynamic landuse datasets.
! ASSUMES that only have one dynamic landuse dataset.
!
! !PUBLIC TYPES:
  implicit none
  private
  save
  public :: pftdyn_init
  public :: pftdyn_interp
  public :: pftdyn_wbal_init
  public :: pftdyn_wbal
  public :: pftdyn_cnbal
  public :: pftwt_init
  public :: pftwt_interp
  public :: CNHarvest
  public :: CNHarvestPftToColumn
!
! !REVISION HISTORY:
! Created by Peter Thornton
! slevis modified to handle CNDV and crop model
! 19 May 2009: PET - modified to handle harvest fluxes
!
!EOP
!
! ! PRIVATE TYPES
  integer , pointer   :: yearspft(:)
  real(r8), pointer   :: wtpft1(:,:)   
  real(r8), pointer   :: wtpft2(:,:)
  real(r8), pointer   :: harvest(:)   
  real(r8), pointer   :: wtcol_old(:)
  integer :: nt1
  integer :: nt2
  integer :: ntimes
  logical :: do_harvest
  type(file_desc_t)  :: ncid   ! netcdf id
!---------------------------------------------------------------------------

contains
  
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_init
!
! !INTERFACE:
  subroutine pftdyn_init()
!
! !DESCRIPTION:
! Initialize dynamic landuse dataset (position it to the right time samples
! that bound the initial model date)
!
! !USES:
    use clm_time_manager, only : get_curr_date
    use clm_varctl  , only : fpftdyn
    use clm_varpar  , only : numpft, maxpatch_pft
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i,j,m,n,g                       ! indices
    real(r8) :: sumpct                          ! sum for error check
    integer  :: varid                           ! netcdf ids
    integer  :: year                            ! year (0, ...) for nstep+1
    integer  :: mon                             ! month (1, ..., 12) for nstep+1
    integer  :: day                             ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                             ! seconds into current date for nstep+1
    integer  :: ier, ret                        ! error status
    logical  :: found                           ! true => input dataset bounding dates found
    logical  :: readvar	                        ! true => variable is on input dataset
    integer  :: begg,endg                       ! beg/end indices for land gridcells
    integer  :: begl,endl                       ! beg/end indices for land landunits
    integer  :: begc,endc                       ! beg/end indices for land columns
    integer  :: begp,endp                       ! beg/end indices for land pfts
    real(r8), pointer :: pctgla(:)          ! percent of gcell is glacier
    real(r8), pointer :: pctlak(:)          ! percent of gcell is lake
    real(r8), pointer :: pctwet(:)          ! percent of gcell is wetland
    real(r8), pointer :: pcturb(:)          ! percent of gcell is urbanized
    type(gridcell_type), pointer :: gptr        ! pointer to gridcell derived subtype
    character(len=256) :: locfn                 ! local file name
    character(len= 32) :: subname='pftdyn_init' ! subroutine name
 !-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Error check

    if ( maxpatch_pft /= numpft+1 )then
       call endrun( subname//' maxpatch_pft does NOT equal numpft+1 -- this is invalid for dynamic PFT case' )
    end if

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg))

    ! Set pointers into derived type

    gptr => grc

    ! pctspec must be saved between time samples
    ! position to first time sample - assume that first time sample must match starting date
    ! check consistency -  special landunits, grid, frac and mask
    ! only do this once

    ! read data PCT_PFT corresponding to correct year

    allocate(wtpft1(begg:endg,0:numpft), wtpft2(begg:endg,0:numpft), stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for wtpft1, wtpft2' )
    end if
    
    allocate(harvest(begg:endg),stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for harvest')
    end if

    allocate(wtcol_old(begp:endp),stat=ier)
    if (ier /= 0) then
       call endrun( subname//' allocation error for wtcol_old' )
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read pft dynamic landuse data .....'
    end if

    ! Obtain file
    call getfil (fpftdyn, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Obtain pft years from dynamic landuse file
    
    call ncd_inqdid(ncid, 'time', varid)
    call ncd_inqdlen(ncid, varid, ntimes)

    ! Consistency check
    
    call check_dim(ncid, 'lsmpft', numpft+1)

    allocate (yearspft(ntimes), stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname//' allocation error for yearspft'; call endrun()
    end if

    call ncd_io(ncid=ncid, varname='YEAR', flag='read', data=yearspft)

    call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_WETLAND NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname= 'PCT_LAKE', flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname= 'PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on pftdyn file' )

    call ncd_io(ncid=ncid, varname= 'PCT_URBAN'  , flag='read', data=pcturb, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on pftdyn file' )

    ! Consistency check
    do g = begg,endg
    !   this was causing a fail, even though values are the same to within 1e-15
    !   if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
       if (abs((pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g))-pctspec(g)) > 1e-13_r8) then 
          write(iulog,*) subname//'mismatch between input pctspec = ',&
                     pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g),&
                    ' and that obtained from surface dataset ', pctspec(g),' at g= ',g
           call endrun()
       end if
    end do

    ! Determine if current date spans the years
    ! If current year is less than first dynamic PFT timeseries year,
    ! then use the first year from dynamic pft file for both nt1 and nt2,
    ! forcing constant weights until the model year enters the dynamic
    ! pft dataset timeseries range.
    ! If current year is equal to or greater than the last dynamic pft
    ! timeseries year, then use the last year for both nt1 and nt2, 
    ! forcing constant weights for the remainder of the simulation.
    ! This mechanism permits the introduction of a dynamic pft period in the middle
    ! of a simulation, with constant weights before and after the dynamic period.
    ! PET: harvest - since harvest is specified as a rate for each year, this
    ! approach will not work. Instead, need to seta flag that indicates harvest is
    ! zero for the period before the beginning and after the end of the dynpft timeseries.

    call get_curr_date(year, mon, day, sec)

    if (year < yearspft(1)) then
       nt1 = 1
       nt2 = 1
       do_harvest = .false.
    else if (year >= yearspft(ntimes)) then
       nt1 = ntimes
       nt2 = ntimes
       do_harvest = .false.
    else
       found = .false.
       do n = 1,ntimes-1 
          if (year == yearspft(n)) then
             nt1 = n
             nt2 = nt1 + 1
             found = .true.
             do_harvest = .true.
          end if   
       end do
       if (.not. found) then
          write(iulog,*) subname//' error: model year not found in pftdyn timeseries'
          write(iulog,*)'model year = ',year
          call endrun()
       end if
    end if

    ! Get pctpft time samples bracketing the current time

    if (masterproc) then
       write(iulog,*) 'Get PFTDYN data for year: ', yearspft(nt1)
    end if
    call pftdyn_getdata(nt1, wtpft1, begg,endg,0,numpft)
    if (masterproc) then
       write(iulog,*) 'Get PFTDYN data for year: ', yearspft(nt2)
    end if
    call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)
    
    if (use_cn) then
       ! Get harvest rate at the nt1 time
       call pftdyn_getharvest(nt1,begg,endg)
    end if

    ! convert weights from percent to proportion
    do m = 0,numpft
       do g = begg,endg
          wtpft1(g,m) = wtpft1(g,m)/100._r8
          wtpft2(g,m) = wtpft2(g,m)/100._r8
       end do
    end do
       
    deallocate(pctgla,pctlak,pctwet,pcturb)

  end subroutine pftdyn_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_interp
!
! !INTERFACE:
  subroutine pftdyn_interp()
!
! !DESCRIPTION:
! Time interpolate dynamic landuse data to get pft weights for model time
! Note that harvest data are stored as rates (not weights) and so time interpolation is 
! not necessary - the harvest rate is held constant through the year.  This is consistent with
! the treatment of changing PFT weights, where interpolation of the annual endpoint weights leads to 
! a constant rate of change in PFT weight through the year, with abrupt changes in the rate at
! annual boundaries. This routine is still used to get the next harvest time slice, when needed.
! This routine is also used to turn off the harvest switch when the model year runs past the end of
! the dynpft time series.
!
! !USES:
    use clm_time_manager, only : get_curr_date, get_curr_calday, &
                                 get_days_per_year
    use clm_varcon      , only : istsoil
    use clm_varcon      , only : istcrop
    use clm_varpar      , only : numpft
    implicit none
!
!
! !LOCAL VARIABLES:
!EOP
!
! !ARGUMENTS:
    integer  :: begg,endg        ! beg/end indices for land gridcells
    integer  :: begl,endl        ! beg/end indices for land landunits
    integer  :: begc,endc        ! beg/end indices for land columns
    integer  :: begp,endp        ! beg/end indices for land pfts
    integer  :: i,j,m,p,l,g,c    ! indices
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    real(r8) :: cday             ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: days_per_year    ! days per year
    integer  :: ier              ! error status
    integer  :: lbc,ubc
    real(r8) :: wt1              ! time interpolation weights
    real(r8), pointer :: wtpfttot1(:)            ! summation of pft weights for renormalization
    real(r8), pointer :: wtpfttot2(:)            ! summation of pft weights for renormalization
    real(r8), parameter :: wtpfttol = 1.e-10     ! tolerance for pft weight renormalization
    type(gridcell_type), pointer :: gptr         ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_interp' ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Set pointers into derived type

    gptr => grc
    lptr => lun
    pptr => pft

    allocate(wtpfttot1(begc:endc),wtpfttot2(begc:endc))
    wtpfttot1(:) = 0._r8
    wtpfttot2(:) = 0._r8

    ! Interpolat pctpft to current time step - output in pctpft_intp
    ! Map interpolated pctpft to subgrid weights
    ! assumes that maxpatch_pft = numpft + 1, that each landunit has only 1 column, 
    ! SCAM and CNDV have not been defined, and create_croplandunit = .false.

    ! If necessary, obtain new time sample

    ! Get current date

    call get_curr_date(year, mon, day, sec)

    ! Obtain new time sample if necessary.
    ! The first condition is the regular crossing of a year boundary
    ! when within the dynpft timeseries range. The second condition is
    ! the case of the first entry into the dynpft timeseries range from
    ! an earlier period of constant weights.

    if (year > yearspft(nt1) .or. (nt1 == 1 .and. nt2 == 1 .and. year == yearspft(1))) then

       if (year >= yearspft(ntimes)) then
          nt1 = ntimes
          nt2 = ntimes
       else
          nt1        = nt2
          nt2        = nt2 + 1
          do_harvest = .true.
       end if
       
       if (year > yearspft(ntimes)) then
          do_harvest = .false.
       endif
       
       if (nt2 > ntimes .and. masterproc) then
          write(iulog,*)subname,' error - current year is past input data boundary'
       end if
       
       do m = 0,numpft
          do g = begg,endg
             wtpft1(g,m) = wtpft2(g,m)
          end do
       end do

       if (masterproc) then
          write(iulog,*) 'Get PFTDYN data for year: ', yearspft(nt2)
       end if
       call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)
       if (use_cn) then
          call pftdyn_getharvest(nt1,begg,endg)
       end if

       do m = 0,numpft
          do g = begg,endg
             wtpft2(g,m) = wtpft2(g,m)/100._r8
          end do
       end do
    
    end if  ! end of need new data if-block 

    ! Interpolate pft weight to current time

    cday          = get_curr_calday() 
    days_per_year = get_days_per_year()

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year

    do p = begp,endp
       c = pptr%column(p)
       g = pptr%gridcell(p)
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
          m = pptr%itype(p)
          wtcol_old(p)      = pptr%wtcol(p)
!         --- recoded for roundoff performance, tcraig 3/07 from k.lindsay
!         pptr%wtgcell(p)   = wtpft1(g,m)*wt1 + wtpft2(g,m)*wt2
          wtpfttot1(c) = wtpfttot1(c)+pptr%wtgcell(p)    
          pptr%wtgcell(p)   = wtpft2(g,m) + wt1*(wtpft1(g,m)-wtpft2(g,m))
          pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
          pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
          wtpfttot2(c) = wtpfttot2(c)+pptr%wtgcell(p)
       end if

    end do

!   Renormalize pft weights so that sum of pft weights relative to grid cell 
!   remain constant even as land cover changes.  Doing this eliminates 
!   soil balance error warnings.  (DML, 4/8/2009)
    do p = begp,endp
       c = pptr%column(p)
       g = pptr%gridcell(p)
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
          if (wtpfttot2(c) /= 0 .and. &
              abs(wtpfttot1(c)-wtpfttot2(c)) > wtpfttol) then
             pptr%wtgcell(p)   = (wtpfttot1(c)/wtpfttot2(c))*pptr%wtgcell(p)
             pptr%wtlunit(p)   = pptr%wtgcell(p) / lptr%wtgcell(l)
             pptr%wtcol(p)     = pptr%wtgcell(p) / lptr%wtgcell(l)
          end if
       end if

    end do
   
    deallocate(wtpfttot1,wtpfttot2) 
    
  end subroutine pftdyn_interp

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_getdata
!
! !INTERFACE:
  subroutine pftdyn_getdata(ntime, pctpft, begg, endg, pft0, maxpft)
!
! !DESCRIPTION:
! Obtain dynamic landuse data (pctpft) and make sure that
! percentage of PFTs sum to 100% cover for vegetated landunit
!
! !USES:
    use clm_varpar  , only : numpft
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: ntime
    integer , intent(in)  :: begg,endg,pft0,maxpft
    real(r8), intent(out) :: pctpft(begg:endg,pft0:maxpft)
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i,j,m,n
    integer  :: err, ierr, ret
    real(r8) :: sumpct,sumerr                     ! temporary
    real(r8), pointer :: arrayl(:,:)              ! temporary array
    logical  :: readvar
   
    character(len=32) :: subname='pftdyn_getdata' ! subroutine name
!-----------------------------------------------------------------------
    
    allocate(arrayl(begg:endg,pft0:maxpft))	
    call ncd_io(ncid=ncid, varname= 'PCT_PFT', flag='read', data=arrayl, &
         dim1name=grlnd, nt=ntime, readvar=readvar)
    pctpft(begg:endg,pft0:maxpft) = arrayl(begg:endg,pft0:maxpft)
    deallocate(arrayl)		
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_PFT NOT on pftdyn file' )

    err = 0
    do n = begg,endg
       if (pctspec(n) < 100._r8) then
          sumpct = 0._r8
          do m = 0, numpft
             sumpct = sumpct + pctpft(n,m) * 100._r8/(100._r8-pctspec(n))
          end do
          if (abs(sumpct - 100._r8) > 0.1_r8) then
             err = 1; ierr = n; sumerr = sumpct
          end if
          if (sumpct < -0.000001_r8) then
             err = 2; ierr = n; sumerr = sumpct
          end if
       end if
    end do
    if (err == 1) then
       write(iulog,*) subname,' error: sum(pct) over numpft+1 is not = 100.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
       call endrun()
    else if (err == 2) then
       write(iulog,*)subname,' error: sum(pct) over numpft+1 is < 0.',sumerr,ierr,pctspec(ierr),pctpft(ierr,:)
       call endrun()
    end if
    
  end subroutine pftdyn_getdata

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_getharvest
!
! !INTERFACE:
  subroutine pftdyn_getharvest(ntime, begg, endg)
!
! !DESCRIPTION:
! Obtain harvest data 
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: ntime
    integer , intent(IN)  :: begg     ! beg indices for land gridcells
    integer , intent(IN)  :: endg     ! end indices for land gridcells
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8), pointer :: arrayl(:)                   ! temporary array
    logical :: readvar 
    character(len=32) :: subname='pftdyn_getharvest' ! subroutine name
!-----------------------------------------------------------------------
    
    allocate(arrayl(begg:endg))
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_VH1', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_VH1 not on pftdyn file' )
    harvest(begg:endg) = arrayl(begg:endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_VH2', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_VH2 not on pftdyn file' )
    harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_SH1', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH1 not on pftdyn file' )
    harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_SH2', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH2 not on pftdyn file' )
    harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)
    
    call ncd_io(ncid=ncid, varname= 'HARVEST_SH3', flag='read', data=arrayl, dim1name=grlnd, &
         nt=ntime, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: HARVEST_SH3 not on pftdyn file' )
    harvest(begg:endg) = harvest(begg:endg) + arrayl(begg:endg)

    deallocate(arrayl)

  end subroutine pftdyn_getharvest

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal_init
!
! !INTERFACE:
  subroutine pftdyn_wbal_init( begc, endc )
!
! !DESCRIPTION:
! initialize the column-level mass-balance correction term.
! Called in every timestep.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(IN)  :: begc, endc    ! proc beginning and ending column indices
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: c             ! indices
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    cptr => col

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every timestep
    
    do c = begc,endc
       cwf%h2ocan_loss(c) = 0._r8
    end do
    
  end subroutine pftdyn_wbal_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_wbal
!
! !INTERFACE:
  subroutine pftdyn_wbal( begg, endg, begc, endc, begp, endp )
!
! !DESCRIPTION:
! modify pft-level state and flux variables to maintain water balance with
! dynamic pft-weights.
! Canopy water balance does not need to consider harvest fluxes, since pft weights are
! not affected by harvest.
!
! !USES:
    use clm_varcon  , only : istsoil
    use clm_varcon  , only : istcrop
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer, intent(IN)  :: begg     ! beg indices for land gridcells
    integer, intent(IN)  :: endg     ! end indices for land gridcells
    integer, intent(IN)  :: begc     ! beg indices for land columns
    integer, intent(IN)  :: endc     ! end indices for land columns
    integer, intent(IN)  :: begp     ! beg indices for land plant function types
    integer, intent(IN)  :: endp     ! end indices for land plant function types
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: pi,p,c,l,g    ! indices
    integer  :: ier           ! error code
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: loss_h2ocan(:) ! canopy water mass loss due to weight shift
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
    character(len=32) :: subname='pftdyn_wbal' ! subroutine name
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    lptr => lun
    cptr => col
    pptr => pft

    ! Allocate loss_h2ocan
    allocate(loss_h2ocan(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for loss_h2ocan'; call endrun()
    end if

    ! Get time step

    dtime = get_step_size()

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every weight-shifting timestep

    do c = begc,endc
       cwf%h2ocan_loss(c) = 0._r8 ! is this OR pftdyn_wbal_init redundant?
    end do

    do p = begp,endp
       l = pptr%landunit(p)
       loss_h2ocan(p) = 0._r8

       if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then

          ! calculate the change in weight for the timestep
          dwt = pptr%wtcol(p)-wtcol_old(p)
  
          if (dwt > 0._r8) then
          
             ! if the pft gained weight, then the 
             ! initial canopy water state is redistributed over the
             ! new (larger) area, conserving mass.

             pws%h2ocan(p) = pws%h2ocan(p) * (wtcol_old(p)/pptr%wtcol(p))
          
          else if (dwt < 0._r8) then
          
             ! if the pft lost weight on the timestep, then the canopy water
             ! mass associated with the lost weight is directed to a 
             ! column-level flux term that gets added to the precip flux
             ! for every pft calculation in Hydrology1()
             
             init_h2ocan = pws%h2ocan(p) * wtcol_old(p)
             loss_h2ocan(p) = pws%h2ocan(p) * (-dwt)
             new_h2ocan = init_h2ocan - loss_h2ocan(p)
             if (abs(new_h2ocan) < 1e-8_r8) then
                new_h2ocan = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if
             if (pptr%wtcol(p) /= 0._r8) then  
                pws%h2ocan(p) = new_h2ocan/pptr%wtcol(p)
             else
                pws%h2ocan(p) = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if 
       

          end if

       end if
    end do

    do pi = 1,max_pft_per_col
       do c = begc,endc
          if (pi <= cptr%npfts(c)) then
             p = cptr%pfti(c) + pi - 1
             cwf%h2ocan_loss(c) = cwf%h2ocan_loss(c) + loss_h2ocan(p)/dtime
          end if
       end do
    end do

    ! Deallocate loss_h2ocan
    deallocate(loss_h2ocan)
    
  end subroutine pftdyn_wbal
  
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftdyn_cnbal
!
! !INTERFACE:
  subroutine pftdyn_cnbal( begc, endc, begp, endp )
!
! !DESCRIPTION:
! modify pft-level state and flux variables to maintain carbon and nitrogen balance with
! dynamic pft-weights.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_const_mod,only : SHR_CONST_PDB
    use decompMod   , only : get_proc_bounds
    use clm_varcon  , only : istsoil
    use clm_varcon  , only : istcrop
    use clm_varpar  , only : numveg, numpft
    use pftvarcon   , only : pconv, pprod10, pprod100
    use clm_varcon  , only : c13ratio
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer, intent(IN)  :: begp, endp    ! proc beginning and ending pft indices
    integer, intent(IN)  :: begc, endc    ! proc beginning and ending column indices
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: pi,p,c,l,g    ! indices
    integer  :: ier           ! error code
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: dt            ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: dwt_leafc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: conv_cflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod10_cflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod100_cflux(:)      ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)    ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)      ! pft-level mass loss due to weight shift
    real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1         ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1         ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: t1,t2,wt_new,wt_old
    real(r8) :: init_state, change_state, new_state
    real(r8) :: tot_leaf, pleaf, pstor, pxfer
    real(r8) :: leafc_seed, leafn_seed
    real(r8) :: deadstemc_seed, deadstemn_seed
    real(r8) :: leafc13_seed, deadstemc13_seed
    real(r8), pointer :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
    integer          , pointer   :: pcolumn(:)   ! column of corresponding pft
    character(len=32) :: subname='pftdyn_cbal' ! subroutine name
!-----------------------------------------------------------------------
    
    ! Set pointers into derived type

    lptr    => lun
    cptr    => col
    pptr    => pft
    pcolumn => pptr%column

    ! Allocate pft-level mass loss arrays
    allocate(dwt_leafc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc_seed'; call endrun()
    end if
    allocate(dwt_leafn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafn_seed'; call endrun()
    end if
    if (use_c13) then
       allocate(dwt_leafc13_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc13_seed'; call endrun()
       end if
    endif
    allocate(dwt_deadstemc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'; call endrun()
    end if
    allocate(dwt_deadstemn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'; call endrun()
    end if
    if (use_c13) then
       allocate(dwt_deadstemc13_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'; call endrun()
       end if
    endif
    allocate(dwt_frootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
       write(iulog,*)subname,' allocation error for dwt_deadcrootc_to_litter'; call endrun()
    end if
    if (use_c13) then
       allocate(dwt_frootc13_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc13_to_litter'; call endrun()
       end if
       allocate(dwt_livecrootc13_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc13_to_litter'; call endrun()
       end if
       allocate(dwt_deadcrootc13_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc13_to_litter'; call endrun()
       end if
    endif
    allocate(dwt_frootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootn_to_litter'; call endrun()
    end if
    allocate(dwt_livecrootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootn_to_litter'; call endrun()
    end if
    allocate(dwt_deadcrootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootn_to_litter'; call endrun()
    end if
    allocate(conv_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_cflux'; call endrun()
    end if
    allocate(prod10_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_cflux'; call endrun()
    end if
    allocate(prod100_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_cflux'; call endrun()
    end if
    if (use_c13) then
       allocate(conv_c13flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c13flux'; call endrun()
       end if
       allocate(prod10_c13flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c13flux'; call endrun()
       end if
       allocate(prod100_c13flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c13flux'; call endrun()
       end if
    endif
    allocate(conv_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_nflux'; call endrun()
    end if
    allocate(prod10_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_nflux'; call endrun()
    end if
    allocate(prod100_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_nflux'; call endrun()
    end if

    ! Get time step
    dt = real( get_step_size(), r8 )

	do p = begp,endp
                c = pcolumn(p)
		! initialize all the pft-level local flux arrays
		dwt_leafc_seed(p) = 0._r8
		dwt_leafn_seed(p) = 0._r8
                if (use_c13) then
                   dwt_leafc13_seed(p) = 0._r8
                endif
		dwt_deadstemc_seed(p) = 0._r8
		dwt_deadstemn_seed(p) = 0._r8
                if (use_c13) then
                   dwt_deadstemc13_seed(p) = 0._r8
                endif
		dwt_frootc_to_litter(p) = 0._r8
		dwt_livecrootc_to_litter(p) = 0._r8
		dwt_deadcrootc_to_litter(p) = 0._r8
                if (use_c13) then
                   dwt_frootc13_to_litter(p) = 0._r8
                   dwt_livecrootc13_to_litter(p) = 0._r8
                   dwt_deadcrootc13_to_litter(p) = 0._r8
                endif
		dwt_frootn_to_litter(p) = 0._r8
		dwt_livecrootn_to_litter(p) = 0._r8
		dwt_deadcrootn_to_litter(p) = 0._r8
		conv_cflux(p) = 0._r8
		prod10_cflux(p) = 0._r8
		prod100_cflux(p) = 0._r8
                if (use_c13) then
                   conv_c13flux(p) = 0._r8
                   prod10_c13flux(p) = 0._r8
                   prod100_c13flux(p) = 0._r8
                endif
		conv_nflux(p) = 0._r8
		prod10_nflux(p) = 0._r8
		prod100_nflux(p) = 0._r8
       
		l = pptr%landunit(p)
		if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then

			! calculate the change in weight for the timestep
			dwt = pptr%wtcol(p)-wtcol_old(p)

			! PFTs for which weight increases on this timestep
			if (dwt > 0._r8) then

				! first identify PFTs that are initiating on this timestep
				! and set all the necessary state and flux variables
				if (wtcol_old(p) == 0._r8) then

					! set initial conditions for PFT that is being initiated
					! in this time step.  Based on the settings in cnIniTimeVar.

					! pft-level carbon state variables
					pcs%leafc(p)              = 0._r8
					pcs%leafc_storage(p)      = 0._r8
					pcs%leafc_xfer(p)         = 0._r8
					pcs%frootc(p)             = 0._r8
					pcs%frootc_storage(p)     = 0._r8
					pcs%frootc_xfer(p)        = 0._r8
					pcs%livestemc(p)          = 0._r8
					pcs%livestemc_storage(p)  = 0._r8
					pcs%livestemc_xfer(p)     = 0._r8
					pcs%deadstemc(p)          = 0._r8
					pcs%deadstemc_storage(p)  = 0._r8
					pcs%deadstemc_xfer(p)     = 0._r8
					pcs%livecrootc(p)         = 0._r8
					pcs%livecrootc_storage(p) = 0._r8
					pcs%livecrootc_xfer(p)    = 0._r8
					pcs%deadcrootc(p)         = 0._r8
					pcs%deadcrootc_storage(p) = 0._r8
					pcs%deadcrootc_xfer(p)    = 0._r8
					pcs%gresp_storage(p)      = 0._r8
					pcs%gresp_xfer(p)         = 0._r8
					pcs%cpool(p)              = 0._r8
					pcs%xsmrpool(p)           = 0._r8
					pcs%pft_ctrunc(p)         = 0._r8
					pcs%dispvegc(p)           = 0._r8
					pcs%storvegc(p)           = 0._r8
					pcs%totvegc(p)            = 0._r8
					pcs%totpftc(p)            = 0._r8

					! pft-level carbon-13 state variables
                                        if (use_c13) then
                                           pc13s%leafc(p)              = 0._r8
                                           pc13s%leafc_storage(p)      = 0._r8
                                           pc13s%leafc_xfer(p)         = 0._r8
                                           pc13s%frootc(p)             = 0._r8
                                           pc13s%frootc_storage(p)     = 0._r8
                                           pc13s%frootc_xfer(p)        = 0._r8
                                           pc13s%livestemc(p)          = 0._r8
                                           pc13s%livestemc_storage(p)  = 0._r8
                                           pc13s%livestemc_xfer(p)     = 0._r8
                                           pc13s%deadstemc(p)          = 0._r8
                                           pc13s%deadstemc_storage(p)  = 0._r8
                                           pc13s%deadstemc_xfer(p)     = 0._r8
                                           pc13s%livecrootc(p)         = 0._r8
                                           pc13s%livecrootc_storage(p) = 0._r8
                                           pc13s%livecrootc_xfer(p)    = 0._r8
                                           pc13s%deadcrootc(p)         = 0._r8
                                           pc13s%deadcrootc_storage(p) = 0._r8
                                           pc13s%deadcrootc_xfer(p)    = 0._r8
                                           pc13s%gresp_storage(p)      = 0._r8
                                           pc13s%gresp_xfer(p)         = 0._r8
                                           pc13s%cpool(p)              = 0._r8
                                           pc13s%xsmrpool(p)           = 0._r8
                                           pc13s%pft_ctrunc(p)         = 0._r8
                                           pc13s%dispvegc(p)           = 0._r8
                                           pc13s%storvegc(p)           = 0._r8
                                           pc13s%totvegc(p)            = 0._r8
                                           pc13s%totpftc(p)            = 0._r8
                                        endif

					! pft-level nitrogen state variables
					pns%leafn(p)	           = 0._r8
					pns%leafn_storage(p)      = 0._r8
					pns%leafn_xfer(p)         = 0._r8
					pns%frootn(p)	           = 0._r8
					pns%frootn_storage(p)     = 0._r8
					pns%frootn_xfer(p)        = 0._r8
					pns%livestemn(p)	       = 0._r8
					pns%livestemn_storage(p)  = 0._r8
					pns%livestemn_xfer(p)     = 0._r8
					pns%deadstemn(p)	       = 0._r8
					pns%deadstemn_storage(p)  = 0._r8
					pns%deadstemn_xfer(p)     = 0._r8
					pns%livecrootn(p)         = 0._r8
					pns%livecrootn_storage(p) = 0._r8
					pns%livecrootn_xfer(p)    = 0._r8
					pns%deadcrootn(p)         = 0._r8
					pns%deadcrootn_storage(p) = 0._r8
					pns%deadcrootn_xfer(p)    = 0._r8
					pns%retransn(p)	       = 0._r8
					pns%npool(p)	           = 0._r8
					pns%pft_ntrunc(p)         = 0._r8
					pns%dispvegn(p)           = 0._r8
					pns%storvegn(p)           = 0._r8
					pns%totvegn(p)            = 0._r8
					pns%totpftn (p)           = 0._r8

					! initialize same flux and epv variables that are set
					! in CNiniTimeVar
					pcf%psnsun(p) = 0._r8
					pcf%psnsha(p) = 0._r8
                                        if (use_c13) then
                                           pc13f%psnsun(p) = 0._r8
                                           pc13f%psnsha(p) = 0._r8
                                        endif
					pps%laisun(p) = 0._r8
					pps%laisha(p) = 0._r8
					pps%lncsun(p) = 0._r8
					pps%lncsha(p) = 0._r8
					pps%vcmxsun(p) = 0._r8
					pps%vcmxsha(p) = 0._r8
                                        if (use_c13) then
                                           pps%alphapsnsun(p) = 0._r8
                                           pps%alphapsnsha(p) = 0._r8
                                        endif

					pepv%dormant_flag(p) = 1._r8
					pepv%days_active(p) = 0._r8
					pepv%onset_flag(p) = 0._r8
					pepv%onset_counter(p) = 0._r8
					pepv%onset_gddflag(p) = 0._r8
					pepv%onset_fdd(p) = 0._r8
					pepv%onset_gdd(p) = 0._r8
					pepv%onset_swi(p) = 0.0_r8
					pepv%offset_flag(p) = 0._r8
					pepv%offset_counter(p) = 0._r8
					pepv%offset_fdd(p) = 0._r8
					pepv%offset_swi(p) = 0._r8
					pepv%lgsf(p) = 0._r8
					pepv%bglfr(p) = 0._r8
					pepv%bgtr(p) = 0._r8
					! difference from CNiniTimeVar: using column-level
					! information to initialize annavg_t2m.
					pepv%annavg_t2m(p) = cps%cannavg_t2m(c)
					pepv%tempavg_t2m(p) = 0._r8
					pepv%gpp(p) = 0._r8
					pepv%availc(p) = 0._r8
					pepv%xsmrpool_recover(p) = 0._r8
                                        if (use_c13) then
                                           pepv%xsmrpool_c13ratio(p) = c13ratio
                                        endif
					pepv%alloc_pnow(p) = 1._r8
					pepv%c_allometry(p) = 0._r8
					pepv%n_allometry(p) = 0._r8
					pepv%plant_ndemand(p) = 0._r8
					pepv%tempsum_potential_gpp(p) = 0._r8
					pepv%annsum_potential_gpp(p) = 0._r8
					pepv%tempmax_retransn(p) = 0._r8
					pepv%annmax_retransn(p) = 0._r8
					pepv%avail_retransn(p) = 0._r8
					pepv%plant_nalloc(p) = 0._r8
					pepv%plant_calloc(p) = 0._r8
					pepv%excess_cflux(p) = 0._r8
					pepv%downreg(p) = 0._r8
					pepv%prev_leafc_to_litter(p) = 0._r8
					pepv%prev_frootc_to_litter(p) = 0._r8
					pepv%tempsum_npp(p) = 0._r8
					pepv%annsum_npp(p) = 0._r8
                                        if (use_c13) then
                                           pepv%rc13_canair(p) = 0._r8
                                           pepv%rc13_psnsun(p) = 0._r8
                                           pepv%rc13_psnsha(p) = 0._r8
                                        endif

				end if  ! end initialization of new pft

				! (still in dwt > 0 block)

				! set the seed sources for leaf and deadstem
				! leaf source is split later between leaf, leaf_storage, leaf_xfer
				leafc_seed   = 0._r8
				leafn_seed   = 0._r8
                                if (use_c13) then
                                   leafc13_seed = 0._r8
                                endif
				deadstemc_seed   = 0._r8
				deadstemn_seed   = 0._r8
                                if (use_c13) then
                                   deadstemc13_seed = 0._r8
                                endif
				if (pptr%itype(p) /= 0) then
					leafc_seed = 1._r8
					leafn_seed  = leafc_seed / pftcon%leafcn(pptr%itype(p))
					if (pftcon%woody(pptr%itype(p)) == 1._r8) then
						deadstemc_seed = 0.1_r8
						deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(pptr%itype(p))
					end if

                                        if (use_c13) then
                                           ! 13c state is initialized assuming del13c = -28 permil for C3, and -13 permil for C4.
                                           ! That translates to ratios of (13c/(12c+13c)) of 0.01080455 for C3, and 0.01096945 for C4
                                           ! based on the following formulae: 
                                           ! r1 (13/12) = PDB + (del13c * PDB)/1000.0
                                           ! r2 (13/(13+12)) = r1/(1+r1)
                                           ! PDB = 0.0112372_R8  (ratio of 13C/12C in Pee Dee Belemnite, C isotope standard)
                                           c3_del13c = -28._r8
                                           c4_del13c = -13._r8
                                           c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
                                           c3_r2 = c3_r1/(1._r8 + c3_r1)
                                           c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
                                           c4_r2 = c4_r1/(1._r8 + c4_r1)
                                           
                                           if (pftcon%c3psn(pptr%itype(p)) == 1._r8) then
                                              leafc13_seed     = leafc_seed     * c3_r2
                                              deadstemc13_seed = deadstemc_seed * c3_r2
                                           else
                                              leafc13_seed     = leafc_seed     * c4_r2
                                              deadstemc13_seed = deadstemc_seed * c4_r2
                                           end if
                                        endif
                                     end if

				! When PFT area expands (dwt > 0), the pft-level mass density 
				! is modified to conserve the original pft mass distributed
				! over the new (larger) area, plus a term to account for the 
				! introduction of new seed source for leaf and deadstem
				t1 = wtcol_old(p)/pptr%wtcol(p)
				t2 = dwt/pptr%wtcol(p)

				tot_leaf = pcs%leafc(p) + pcs%leafc_storage(p) + pcs%leafc_xfer(p)
				pleaf = 0._r8
				pstor = 0._r8
				pxfer = 0._r8
				if (tot_leaf /= 0._r8) then
					! when adding seed source to non-zero leaf state, use current proportions
					pleaf = pcs%leafc(p)/tot_leaf
					pstor = pcs%leafc_storage(p)/tot_leaf
					pxfer = pcs%leafc_xfer(p)/tot_leaf
				else
					! when initiating from zero leaf state, use evergreen flag to set proportions
					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
						pleaf = 1._r8
					else
						pstor = 1._r8
					end if
				end if 
				pcs%leafc(p)         = pcs%leafc(p)*t1         + leafc_seed*pleaf*t2
				pcs%leafc_storage(p) = pcs%leafc_storage(p)*t1 + leafc_seed*pstor*t2
				pcs%leafc_xfer(p)    = pcs%leafc_xfer(p)*t1    + leafc_seed*pxfer*t2
				pcs%frootc(p)  		   = pcs%frootc(p) 			* t1
				pcs%frootc_storage(p)     = pcs%frootc_storage(p) 	* t1
				pcs%frootc_xfer(p) 	   = pcs%frootc_xfer(p)		* t1
				pcs%livestemc(p)		   = pcs%livestemc(p)  		* t1
				pcs%livestemc_storage(p)  = pcs%livestemc_storage(p)  * t1
				pcs%livestemc_xfer(p)     = pcs%livestemc_xfer(p) 	* t1
				pcs%deadstemc(p)     = pcs%deadstemc(p)*t1     + deadstemc_seed*t2
				pcs%deadstemc_storage(p)  = pcs%deadstemc_storage(p)  * t1
				pcs%deadstemc_xfer(p)     = pcs%deadstemc_xfer(p) 	* t1
				pcs%livecrootc(p)  	   = pcs%livecrootc(p) 		* t1
				pcs%livecrootc_storage(p) = pcs%livecrootc_storage(p) * t1
				pcs%livecrootc_xfer(p)    = pcs%livecrootc_xfer(p)	* t1
				pcs%deadcrootc(p)  	   = pcs%deadcrootc(p) 		* t1
				pcs%deadcrootc_storage(p) = pcs%deadcrootc_storage(p) * t1
				pcs%deadcrootc_xfer(p)    = pcs%deadcrootc_xfer(p)	* t1
				pcs%gresp_storage(p)	   = pcs%gresp_storage(p)  	* t1
				pcs%gresp_xfer(p)  	   = pcs%gresp_xfer(p) 		* t1
				pcs%cpool(p)			   = pcs%cpool(p)  			* t1
				pcs%xsmrpool(p)		   = pcs%xsmrpool(p)			* t1
				pcs%pft_ctrunc(p)  	   = pcs%pft_ctrunc(p) 		* t1
				pcs%dispvegc(p)		   = pcs%dispvegc(p)			* t1
				pcs%storvegc(p)		   = pcs%storvegc(p)			* t1
				pcs%totvegc(p) 		   = pcs%totvegc(p)			* t1
				pcs%totpftc(p) 		   = pcs%totpftc(p)			* t1

				! pft-level carbon-13 state variables 
                                if (use_c13) then
                                   tot_leaf = pc13s%leafc(p) + pc13s%leafc_storage(p) + pc13s%leafc_xfer(p)
                                   pleaf = 0._r8
                                   pstor = 0._r8
                                   pxfer = 0._r8
                                   if (tot_leaf /= 0._r8) then
                                      pleaf = pc13s%leafc(p)/tot_leaf
                                      pstor = pc13s%leafc_storage(p)/tot_leaf
                                      pxfer = pc13s%leafc_xfer(p)/tot_leaf
                                   else
                                      ! when initiating from zero leaf state, use evergreen flag to set proportions
                                      if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
                                         pleaf = 1._r8
                                      else
                                         pstor = 1._r8
                                      end if
                                   end if
                                   pc13s%leafc(p)         = pc13s%leafc(p)*t1         + leafc13_seed*pleaf*t2
                                   pc13s%leafc_storage(p) = pc13s%leafc_storage(p)*t1 + leafc13_seed*pstor*t2
                                   pc13s%leafc_xfer(p)    = pc13s%leafc_xfer(p)*t1    + leafc13_seed*pxfer*t2
                                   pc13s%frootc(p)			 = pc13s%frootc(p) 		* t1
                                   pc13s%frootc_storage(p)	         = pc13s%frootc_storage(p) 	* t1
                                   pc13s%frootc_xfer(p)		 = pc13s%frootc_xfer(p)		* t1
                                   pc13s%livestemc(p) 		 = pc13s%livestemc(p)  		* t1
                                   pc13s%livestemc_storage(p)       = pc13s%livestemc_storage(p)      * t1
                                   pc13s%livestemc_xfer(p)	         = pc13s%livestemc_xfer(p) 	* t1
                                   pc13s%deadstemc(p)               = pc13s%deadstemc(p)*t1     + deadstemc13_seed*t2
                                   pc13s%deadstemc_storage(p)       = pc13s%deadstemc_storage(p)      * t1
                                   pc13s%deadstemc_xfer(p)	         = pc13s%deadstemc_xfer(p) 	* t1
                                   pc13s%livecrootc(p)		 = pc13s%livecrootc(p) 		* t1
                                   pc13s%livecrootc_storage(p)      = pc13s%livecrootc_storage(p)     * t1
                                   pc13s%livecrootc_xfer(p)         = pc13s%livecrootc_xfer(p)	* t1
                                   pc13s%deadcrootc(p)		 = pc13s%deadcrootc(p) 		* t1
                                   pc13s%deadcrootc_storage(p)      = pc13s%deadcrootc_storage(p)     * t1
                                   pc13s%deadcrootc_xfer(p)         = pc13s%deadcrootc_xfer(p)	* t1
                                   pc13s%gresp_storage(p) 	         = pc13s%gresp_storage(p)  	* t1
                                   pc13s%gresp_xfer(p)		 = pc13s%gresp_xfer(p) 		* t1
                                   pc13s%cpool(p) 			 = pc13s%cpool(p)  		* t1
                                   pc13s%xsmrpool(p)  		 = pc13s%xsmrpool(p)		* t1
                                   pc13s%pft_ctrunc(p)		 = pc13s%pft_ctrunc(p) 		* t1
                                   pc13s%dispvegc(p)  		 = pc13s%dispvegc(p)		* t1
                                   pc13s%storvegc(p)  		 = pc13s%storvegc(p)		* t1
                                   pc13s%totvegc(p)		 = pc13s%totvegc(p)		* t1
                                   pc13s%totpftc(p)		 = pc13s%totpftc(p)		* t1
                                endif

				tot_leaf = pns%leafn(p) + pns%leafn_storage(p) + pns%leafn_xfer(p)
				pleaf = 0._r8
				pstor = 0._r8
				pxfer = 0._r8
				if (tot_leaf /= 0._r8) then
					pleaf = pns%leafn(p)/tot_leaf
					pstor = pns%leafn_storage(p)/tot_leaf
					pxfer = pns%leafn_xfer(p)/tot_leaf
				else
					! when initiating from zero leaf state, use evergreen flag to set proportions
					if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
						pleaf = 1._r8
					else
						pstor = 1._r8
					end if
				end if 
				! pft-level nitrogen state variables
				pns%leafn(p)         = pns%leafn(p)*t1         + leafn_seed*pleaf*t2
				pns%leafn_storage(p) = pns%leafn_storage(p)*t1 + leafn_seed*pstor*t2
				pns%leafn_xfer(p)    = pns%leafn_xfer(p)*t1    + leafn_seed*pxfer*t2
				pns%frootn(p)  		   = pns%frootn(p) 		* t1
				pns%frootn_storage(p)         = pns%frootn_storage(p) 	* t1
				pns%frootn_xfer(p) 	   = pns%frootn_xfer(p)		* t1
				pns%livestemn(p)		   = pns%livestemn(p)  		* t1
				pns%livestemn_storage(p)      = pns%livestemn_storage(p)      * t1
				pns%livestemn_xfer(p)         = pns%livestemn_xfer(p) 	* t1
				pns%deadstemn(p)              = pns%deadstemn(p)*t1     + deadstemn_seed*t2
				pns%deadstemn_storage(p)      = pns%deadstemn_storage(p)      * t1
				pns%deadstemn_xfer(p)         = pns%deadstemn_xfer(p) 	* t1
				pns%livecrootn(p)  	   = pns%livecrootn(p) 		* t1
				pns%livecrootn_storage(p)     = pns%livecrootn_storage(p)     * t1
				pns%livecrootn_xfer(p)        = pns%livecrootn_xfer(p)	* t1
				pns%deadcrootn(p)  	   = pns%deadcrootn(p) 		* t1
				pns%deadcrootn_storage(p)     = pns%deadcrootn_storage(p)     * t1
				pns%deadcrootn_xfer(p)        = pns%deadcrootn_xfer(p)        * t1
				pns%retransn(p)		   = pns%retransn(p)		* t1
				pns%npool(p)		   = pns%npool(p)  		* t1
				pns%pft_ntrunc(p)  	   = pns%pft_ntrunc(p)        	* t1
				pns%dispvegn(p)		   = pns%dispvegn(p)		* t1
				pns%storvegn(p)		   = pns%storvegn(p)		* t1
				pns%totvegn(p) 		   = pns%totvegn(p)		* t1
				pns%totpftn(p) 		   = pns%totpftn(p)		* t1

				! update temporary seed source arrays
				! These are calculated in terms of the required contributions from
				! column-level seed source
				dwt_leafc_seed(p)   = leafc_seed   * dwt
                                if (use_c13) then
                                   dwt_leafc13_seed(p) = leafc13_seed * dwt
                                endif
				dwt_leafn_seed(p)   = leafn_seed   * dwt
				dwt_deadstemc_seed(p)   = deadstemc_seed   * dwt
                                if (use_c13) then
                                   dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
                                endif
				dwt_deadstemn_seed(p)   = deadstemn_seed   * dwt

			else if (dwt < 0._r8) then

				! if the pft lost weight on the timestep, then the carbon and nitrogen state
				! variables are directed to litter, CWD, and wood product pools.

				! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
				! as negative values, but the fluxes for pft-to-litter are accumulated as 
				! positive values

				! set local weight variables for this pft
				wt_new = pptr%wtcol(p)
				wt_old = wtcol_old(p)

				!---------------
				! C state update
				!---------------

				! leafc 
				ptr => pcs%leafc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! leafc_storage 
				ptr => pcs%leafc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! leafc_xfer 
				ptr => pcs%leafc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! frootc 
				ptr => pcs%frootc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) - change_state
				else
					ptr = 0._r8
					dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) + init_state
				end if

				! frootc_storage 
				ptr => pcs%frootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! frootc_xfer 
				ptr => pcs%frootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livestemc 
				ptr => pcs%livestemc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livestemc_storage 
				ptr => pcs%livestemc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livestemc_xfer 
				ptr => pcs%livestemc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadstemc 
				ptr => pcs%deadstemc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state*pconv(pptr%itype(p))
					prod10_cflux(p) = prod10_cflux(p) + change_state*pprod10(pptr%itype(p))
					prod100_cflux(p) = prod100_cflux(p) + change_state*pprod100(pptr%itype(p))
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state*pconv(pptr%itype(p))
					prod10_cflux(p) = prod10_cflux(p) - init_state*pprod10(pptr%itype(p))
					prod100_cflux(p) = prod100_cflux(p) - init_state*pprod100(pptr%itype(p))
				end if

				! deadstemc_storage 
				ptr => pcs%deadstemc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadstemc_xfer 
				ptr => pcs%deadstemc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livecrootc 
				ptr => pcs%livecrootc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) - change_state
				else
					ptr = 0._r8
					dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) + init_state
				end if

				! livecrootc_storage 
				ptr => pcs%livecrootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! livecrootc_xfer 
				ptr => pcs%livecrootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadcrootc 
				ptr => pcs%deadcrootc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) - change_state
				else
					ptr = 0._r8
					dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) + init_state
				end if

				! deadcrootc_storage 
				ptr => pcs%deadcrootc_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! deadcrootc_xfer 
				ptr => pcs%deadcrootc_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! gresp_storage 
				ptr => pcs%gresp_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! gresp_xfer 
				ptr => pcs%gresp_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! cpool 
				ptr => pcs%cpool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! xsmrpool 
				ptr => pcs%xsmrpool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

				! pft_ctrunc 
				ptr => pcs%pft_ctrunc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					conv_cflux(p) = conv_cflux(p) + change_state
				else
					ptr = 0._r8
					conv_cflux(p) = conv_cflux(p) - init_state
				end if

                                if (use_c13) then
                                   !-----------------
                                   ! C13 state update
                                   !-----------------
                                   
                                   ! set pointers to the conversion and product pool fluxes for this pft
                                   ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                                   dwt_ptr1 => conv_c13flux(p)
                                   dwt_ptr2 => prod10_c13flux(p)
                                   dwt_ptr3 => prod100_c13flux(p)
                                   
                                   ! leafc 
                                   ptr => pc13s%leafc(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! leafc_storage 
                                   ptr => pc13s%leafc_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! leafc_xfer 
                                   ptr => pc13s%leafc_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! frootc 
                                   ptr => pc13s%frootc(p)
                                   dwt_ptr0 => dwt_frootc13_to_litter(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr0 = dwt_ptr0 - change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr0 = dwt_ptr0 + init_state
                                   end if
                                   
                                   ! frootc_storage 
                                   ptr => pc13s%frootc_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if

                                   ! frootc_xfer 
                                   ptr => pc13s%frootc_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! livestemc 
                                   ptr => pc13s%livestemc(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! livestemc_storage 
                                   ptr => pc13s%livestemc_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! livestemc_xfer 
                                   ptr => pc13s%livestemc_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if

                                   ! deadstemc 
                                   ptr => pc13s%deadstemc(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state*pconv(pptr%itype(p))
                                      dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pptr%itype(p))
                                      dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pptr%itype(p))
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state*pconv(pptr%itype(p))
                                      dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pptr%itype(p))
                                      dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pptr%itype(p))
                                   end if
                                   
                                   ! deadstemc_storage 
                                   ptr => pc13s%deadstemc_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! deadstemc_xfer 
                                   ptr => pc13s%deadstemc_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! livecrootc 
                                   ptr => pc13s%livecrootc(p)
                                   dwt_ptr0 => dwt_livecrootc13_to_litter(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr0 = dwt_ptr0 - change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr0 = dwt_ptr0 + init_state
                                   end if
                                   
                                   ! livecrootc_storage 
                                   ptr => pc13s%livecrootc_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! livecrootc_xfer 
                                   ptr => pc13s%livecrootc_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! deadcrootc 
                                   ptr => pc13s%deadcrootc(p)
                                   dwt_ptr0 => dwt_deadcrootc13_to_litter(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr0 = dwt_ptr0 - change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr0 = dwt_ptr0 + init_state
                                   end if
                                   
                                   ! deadcrootc_storage 
                                   ptr => pc13s%deadcrootc_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                   
                                   ! deadcrootc_xfer 
                                   ptr => pc13s%deadcrootc_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if

                                   ! gresp_storage 
                                   ptr => pc13s%gresp_storage(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if

                                   ! gresp_xfer 
                                   ptr => pc13s%gresp_xfer(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if

                                   ! cpool 
                                   ptr => pc13s%cpool(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if

                                   ! pft_ctrunc 
                                   ptr => pc13s%pft_ctrunc(p)
                                   init_state = ptr*wt_old
                                   change_state = ptr*dwt
                                   new_state = init_state+change_state
                                   if (wt_new /= 0._r8) then
                                      ptr = new_state/wt_new
                                      dwt_ptr1 = dwt_ptr1 + change_state
                                   else
                                      ptr = 0._r8
                                      dwt_ptr1 = dwt_ptr1 - init_state
                                   end if
                                endif

				!---------------
				! N state update
				!---------------

				! set pointers to the conversion and product pool fluxes for this pft
				! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
				dwt_ptr1 => conv_nflux(p)
				dwt_ptr2 => prod10_nflux(p)
				dwt_ptr3 => prod100_nflux(p)

				! leafn 
				ptr => pns%leafn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! leafn_storage  
				ptr => pns%leafn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! leafn_xfer  
				ptr => pns%leafn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! frootn 
				ptr => pns%frootn(p)
				dwt_ptr0 => dwt_frootn_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! frootn_storage 
				ptr => pns%frootn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! frootn_xfer  
				ptr => pns%frootn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemn  
				ptr => pns%livestemn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemn_storage 
				ptr => pns%livestemn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livestemn_xfer 
				ptr => pns%livestemn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadstemn 
				ptr => pns%deadstemn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state*pconv(pptr%itype(p))
					dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pptr%itype(p))
					dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pptr%itype(p))
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state*pconv(pptr%itype(p))
					dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pptr%itype(p))
					dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pptr%itype(p))
				end if

				! deadstemn_storage 
				ptr => pns%deadstemn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadstemn_xfer 
				ptr => pns%deadstemn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livecrootn 
				ptr => pns%livecrootn(p)
				dwt_ptr0 => dwt_livecrootn_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! livecrootn_storage  
				ptr => pns%livecrootn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! livecrootn_xfer  
				ptr => pns%livecrootn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadcrootn 
				ptr => pns%deadcrootn(p)
				dwt_ptr0 => dwt_deadcrootn_to_litter(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr0 = dwt_ptr0 - change_state
				else
					ptr = 0._r8
					dwt_ptr0 = dwt_ptr0 + init_state
				end if

				! deadcrootn_storage  
				ptr => pns%deadcrootn_storage(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! deadcrootn_xfer  
				ptr => pns%deadcrootn_xfer(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! retransn  
				ptr => pns%retransn(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if

				! npool  
				ptr => pns%npool(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if
				
				! pft_ntrunc  
				ptr => pns%pft_ntrunc(p)
				init_state = ptr*wt_old
				change_state = ptr*dwt
				new_state = init_state+change_state
				if (wt_new /= 0._r8) then
					ptr = new_state/wt_new
					dwt_ptr1 = dwt_ptr1 + change_state
				else
					ptr = 0._r8
					dwt_ptr1 = dwt_ptr1 - init_state
				end if
				
			end if       ! weight decreasing
		end if           ! is soil
	end do               ! pft loop
    
	! calculate column-level seeding fluxes
	do pi = 1,max_pft_per_col
		do c = begc, endc
			if ( pi <=  cptr%npfts(c) ) then
				p = cptr%pfti(c) + pi - 1
				
				! C fluxes
				ccf%dwt_seedc_to_leaf(c) = ccf%dwt_seedc_to_leaf(c) + dwt_leafc_seed(p)/dt
				ccf%dwt_seedc_to_deadstem(c) = ccf%dwt_seedc_to_deadstem(c) &
                                                                    + dwt_deadstemc_seed(p)/dt
				
                                ! C13 fluxes
                                if (use_c13) then
                                   cc13f%dwt_seedc_to_leaf(c) = cc13f%dwt_seedc_to_leaf(c) + dwt_leafc13_seed(p)/dt
                                   cc13f%dwt_seedc_to_deadstem(c) = cc13f%dwt_seedc_to_deadstem(c) &
                                                                         + dwt_deadstemc13_seed(p)/dt
                                endif
				
				! N fluxes
				cnf%dwt_seedn_to_leaf(c) = cnf%dwt_seedn_to_leaf(c) + dwt_leafn_seed(p)/dt
				cnf%dwt_seedn_to_deadstem(c) = cnf%dwt_seedn_to_deadstem(c) &
                                                                    + dwt_deadstemn_seed(p)/dt
			end if
		end do
	end do


	! calculate pft-to-column for fluxes into litter and CWD pools
	do pi = 1,max_pft_per_col
		do c = begc, endc
			if ( pi <=  cptr%npfts(c) ) then
				p = cptr%pfti(c) + pi - 1

				! fine root litter carbon fluxes
				ccf%dwt_frootc_to_litr1c(c) = ccf%dwt_frootc_to_litr1c(c) + &
                                                            (dwt_frootc_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
				ccf%dwt_frootc_to_litr2c(c) = ccf%dwt_frootc_to_litr2c(c) + &
                                                            (dwt_frootc_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
				ccf%dwt_frootc_to_litr3c(c) = ccf%dwt_frootc_to_litr3c(c) + &
                                                            (dwt_frootc_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt

				! fine root litter C13 fluxes
                                if (use_c13) then
                                   cc13f%dwt_frootc_to_litr1c(c) = cc13f%dwt_frootc_to_litr1c(c) + &
                                                               (dwt_frootc13_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
                                   cc13f%dwt_frootc_to_litr2c(c) = cc13f%dwt_frootc_to_litr2c(c) + &
                                                               (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
                                   cc13f%dwt_frootc_to_litr3c(c) = cc13f%dwt_frootc_to_litr3c(c) + &
                                                               (dwt_frootc13_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt
                                endif

				! fine root litter nitrogen fluxes
				cnf%dwt_frootn_to_litr1n(c) = cnf%dwt_frootn_to_litr1n(c) + &
                                                            (dwt_frootn_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt
				cnf%dwt_frootn_to_litr2n(c) = cnf%dwt_frootn_to_litr2n(c) + &
                                                            (dwt_frootn_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt
				cnf%dwt_frootn_to_litr3n(c) = cnf%dwt_frootn_to_litr3n(c) + &
                                                            (dwt_frootn_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt

				! livecroot fluxes to cwd
				ccf%dwt_livecrootc_to_cwdc(c) = ccf%dwt_livecrootc_to_cwdc(c) + &
                                                            (dwt_livecrootc_to_litter(p))/dt
                                if (use_c13) then
                                   cc13f%dwt_livecrootc_to_cwdc(c) = cc13f%dwt_livecrootc_to_cwdc(c) + &
                                                               (dwt_livecrootc13_to_litter(p))/dt
                                endif
				cnf%dwt_livecrootn_to_cwdn(c) = cnf%dwt_livecrootn_to_cwdn(c) + &
                                                            (dwt_livecrootn_to_litter(p))/dt

				! deadcroot fluxes to cwd
				ccf%dwt_deadcrootc_to_cwdc(c) = ccf%dwt_deadcrootc_to_cwdc(c) + &
                                                            (dwt_deadcrootc_to_litter(p))/dt
                                if (use_c13) then
                                   cc13f%dwt_deadcrootc_to_cwdc(c) = cc13f%dwt_deadcrootc_to_cwdc(c) + &
                                                               (dwt_deadcrootc13_to_litter(p))/dt
                                endif
				cnf%dwt_deadcrootn_to_cwdn(c) = cnf%dwt_deadcrootn_to_cwdn(c) + &
                                                            (dwt_deadcrootn_to_litter(p))/dt
			end if
		end do
	end do

	! calculate pft-to-column for fluxes into product pools and conversion flux
	do pi = 1,max_pft_per_col
		do c = begc,endc
			if (pi <= cptr%npfts(c)) then
				p = cptr%pfti(c) + pi - 1

				! column-level fluxes are accumulated as positive fluxes.
				! column-level C flux updates
				ccf%dwt_conv_cflux(c) = ccf%dwt_conv_cflux(c) - conv_cflux(p)/dt
				ccf%dwt_prod10c_gain(c) = ccf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
				ccf%dwt_prod100c_gain(c) = ccf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt

				! column-level C13 flux updates
                                if (use_c13) then
                                   cc13f%dwt_conv_cflux(c) = cc13f%dwt_conv_cflux(c) - conv_c13flux(p)/dt
                                   cc13f%dwt_prod10c_gain(c) = cc13f%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
                                   cc13f%dwt_prod100c_gain(c) = cc13f%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt
                                endif

				! column-level N flux updates
				cnf%dwt_conv_nflux(c) = cnf%dwt_conv_nflux(c) - conv_nflux(p)/dt
				cnf%dwt_prod10n_gain(c) = cnf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
				cnf%dwt_prod100n_gain(c) = cnf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt

			end if
		end do
	end do

	! Deallocate pft-level flux arrays
        deallocate(dwt_leafc_seed)
        deallocate(dwt_leafn_seed)
        if (use_c13) then
           deallocate(dwt_leafc13_seed)
        endif
        deallocate(dwt_deadstemc_seed)
        deallocate(dwt_deadstemn_seed)
        if (use_c13) then
           deallocate(dwt_deadstemc13_seed)
        endif
	deallocate(dwt_frootc_to_litter)
	deallocate(dwt_livecrootc_to_litter)
	deallocate(dwt_deadcrootc_to_litter)
        if (use_c13) then
           deallocate(dwt_frootc13_to_litter)
           deallocate(dwt_livecrootc13_to_litter)
           deallocate(dwt_deadcrootc13_to_litter)
        endif
	deallocate(dwt_frootn_to_litter)
	deallocate(dwt_livecrootn_to_litter)
	deallocate(dwt_deadcrootn_to_litter)
	deallocate(conv_cflux)
	deallocate(prod10_cflux)
	deallocate(prod100_cflux)
        if (use_c13) then
           deallocate(conv_c13flux)
           deallocate(prod10_c13flux)
           deallocate(prod100_c13flux)
        endif
	deallocate(conv_nflux)
	deallocate(prod10_nflux)
	deallocate(prod100_nflux)
    
end subroutine pftdyn_cnbal

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftwt_init
!
! !INTERFACE:
  subroutine pftwt_init()
!
! !DESCRIPTION:
! Initialize time interpolation of cndv pft weights from annual to time step
!
! !USES:
  use clm_varctl, only : nsrest, nsrStartup
!
! !ARGUMENTS:
    implicit none
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ier, p                        ! error status, do-loop index
    integer  :: begp,endp                     ! beg/end indices for land pfts
    character(len=32) :: subname='pftwt_init' ! subroutine name
    type(pft_type), pointer :: pptr           ! ponter to pft derived subtype
!-----------------------------------------------------------------------

    pptr => pft

    call get_proc_bounds(begp=begp,endp=endp)

    allocate(wtcol_old(begp:endp),stat=ier)
    if (ier /= 0) then
       call endrun( subname//'::ERROR: pftwt_init allocation error for wtcol_old')
    end if

    if (nsrest == nsrStartup) then
       do p = begp,endp
          pdgvs%fpcgrid(p) = pptr%wtcol(p)
          pdgvs%fpcgridold(p) = pptr%wtcol(p)
          wtcol_old(p) = pptr%wtcol(p)
       end do
     else
       do p = begp,endp
          wtcol_old(p) = pptr%wtcol(p)
       end do
    end if

  end subroutine pftwt_init

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pftwt_interp
!
! !INTERFACE:
  subroutine pftwt_interp( begp, endp )
!
! !DESCRIPTION:
! Time interpolate cndv pft weights from annual to time step
!
! !USES:
    use clm_time_manager, only : get_curr_calday, get_curr_date, &
                                 get_days_per_year
    use clm_time_manager, only : get_step_size, get_nstep
    use clm_varcon      , only : istsoil ! CNDV incompatible with dynLU
    use clm_varctl      , only : finidat
!
! !ARGUMENTS:
    implicit none
    integer, intent(IN)  :: begp,endp                ! beg/end indices for land pfts
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: c,g,l,p            ! indices
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: wt1                ! time interpolation weights
    real(r8) :: dtime              ! model time step
    real(r8) :: days_per_year      ! days per year
    integer  :: nstep              ! time step number
    integer  :: year               ! year (0, ...) at nstep + 1
    integer  :: mon                ! month (1, ..., 12) at nstep + 1
    integer  :: day                ! day of month (1, ..., 31) at nstep + 1
    integer  :: sec                ! seconds into current date at nstep + 1
    type(landunit_type), pointer :: lptr ! pointer to landunit derived subtype
    type(pft_type)     , pointer :: pptr ! ...     to pft derived subtype
    character(len=32) :: subname='pftwt_interp' ! subroutine name

! !CALLED FROM:
!  subr. driver
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    lptr => lun
    pptr => pft

    ! Interpolate pft weight to current time step
    ! Map interpolated pctpft to subgrid weights
    ! assumes maxpatch_pft = numpft + 1, each landunit has 1 column, 
    ! SCAM not defined and create_croplandunit = .false.

    nstep         = get_nstep()
    dtime         = get_step_size()
    cday          = get_curr_calday(offset=-int(dtime))
    days_per_year = get_days_per_year()

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year

    call get_curr_date(year, mon, day, sec, offset=int(dtime))

    do p = begp,endp
       g = pptr%gridcell(p)
       l = pptr%landunit(p)

       if (lptr%itype(l) == istsoil .and. lptr%wtgcell(l) > 0._r8) then ! CNDV incompatible with dynLU
          wtcol_old(p)    = pptr%wtcol(p)
          pptr%wtcol(p)   = pdgvs%fpcgrid(p) + &
                     wt1 * (pdgvs%fpcgridold(p) - pdgvs%fpcgrid(p))
          pptr%wtlunit(p) = pptr%wtcol(p)
          pptr%wtgcell(p) = pptr%wtcol(p) * lptr%wtgcell(l)

          if (mon==1 .and. day==1 .and. sec==dtime .and. nstep>0) then
             pdgvs%fpcgridold(p) = pdgvs%fpcgrid(p)
          end if
       end if
    end do

  end subroutine pftwt_interp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvest
!
! !INTERFACE:
subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Harvest mortality routine for coupled carbon-nitrogen code (CN)
!
! !USES:
   use clmtype
   use pftvarcon       , only : noveg, nbrdlf_evr_shrub, pprodharv10
   use clm_varcon      , only : secspday
   use clm_time_manager, only : get_days_per_year
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! column filter for soil points
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   integer , pointer :: pgridcell(:)   ! pft-level index into gridcell-level quantities
   integer , pointer :: ivt(:)         ! pft vegetation type

   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real(r8), pointer :: hrv_leafc_to_litter(:)
   real(r8), pointer :: hrv_frootc_to_litter(:)
   real(r8), pointer :: hrv_livestemc_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: hrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)
   real(r8), pointer :: hrv_xsmrpool_to_atm(:)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
   real(r8), pointer :: hrv_leafn_to_litter(:)
   real(r8), pointer :: hrv_frootn_to_litter(:)
   real(r8), pointer :: hrv_livestemn_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: hrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_livecrootn_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_to_litter(:)
   real(r8), pointer :: hrv_retransn_to_litter(:)
   real(r8), pointer :: hrv_leafn_storage_to_litter(:)
   real(r8), pointer :: hrv_frootn_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: p                         ! pft index
   integer :: g                         ! gridcell index
   integer :: fp                        ! pft filter index
   real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
   real(r8):: m                         ! rate for fractional harvest mortality (1/s)
   real(r8):: days_per_year             ! days per year
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to pft-level arrays
   pgridcell                      => pft%gridcell
   
   ivt                            => pft%itype
   leafc                          => pcs%leafc
   frootc                         => pcs%frootc
   livestemc                      => pcs%livestemc
   deadstemc                      => pcs%deadstemc
   livecrootc                     => pcs%livecrootc
   deadcrootc                     => pcs%deadcrootc
   xsmrpool                       => pcs%xsmrpool
   leafc_storage                  => pcs%leafc_storage
   frootc_storage                 => pcs%frootc_storage
   livestemc_storage              => pcs%livestemc_storage
   deadstemc_storage              => pcs%deadstemc_storage
   livecrootc_storage             => pcs%livecrootc_storage
   deadcrootc_storage             => pcs%deadcrootc_storage
   gresp_storage                  => pcs%gresp_storage
   leafc_xfer                     => pcs%leafc_xfer
   frootc_xfer                    => pcs%frootc_xfer
   livestemc_xfer                 => pcs%livestemc_xfer
   deadstemc_xfer                 => pcs%deadstemc_xfer
   livecrootc_xfer                => pcs%livecrootc_xfer
   deadcrootc_xfer                => pcs%deadcrootc_xfer
   gresp_xfer                     => pcs%gresp_xfer
   leafn                          => pns%leafn
   frootn                         => pns%frootn
   livestemn                      => pns%livestemn
   deadstemn                      => pns%deadstemn
   livecrootn                     => pns%livecrootn
   deadcrootn                     => pns%deadcrootn
   retransn                       => pns%retransn
   leafn_storage                  => pns%leafn_storage
   frootn_storage                 => pns%frootn_storage
   livestemn_storage              => pns%livestemn_storage
   deadstemn_storage              => pns%deadstemn_storage
   livecrootn_storage             => pns%livecrootn_storage
   deadcrootn_storage             => pns%deadcrootn_storage
   leafn_xfer                     => pns%leafn_xfer
   frootn_xfer                    => pns%frootn_xfer
   livestemn_xfer                 => pns%livestemn_xfer
   deadstemn_xfer                 => pns%deadstemn_xfer
   livecrootn_xfer                => pns%livecrootn_xfer
   deadcrootn_xfer                => pns%deadcrootn_xfer
   hrv_leafc_to_litter              => pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => pcf%hrv_livestemc_to_litter
   hrv_deadstemc_to_prod10c         => pcf%hrv_deadstemc_to_prod10c
   hrv_deadstemc_to_prod100c        => pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => pcf%hrv_deadcrootc_to_litter
   hrv_xsmrpool_to_atm              => pcf%hrv_xsmrpool_to_atm
   hrv_leafc_storage_to_litter      => pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => pnf%hrv_livestemn_to_litter
   hrv_deadstemn_to_prod10n         => pnf%hrv_deadstemn_to_prod10n
   hrv_deadstemn_to_prod100n        => pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => pnf%hrv_deadcrootn_xfer_to_litter


   days_per_year = get_days_per_year()

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      
      ! If this is a tree pft, then
      ! get the annual harvest "mortality" rate (am) from harvest array
      ! and convert to rate per second
      if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

         if (do_harvest) then
            am = harvest(g)
            m  = am/(days_per_year * secspday)
         else
            m = 0._r8
         end if   

         ! pft-level harvest carbon fluxes
         ! displayed pools
         hrv_leafc_to_litter(p)               = leafc(p)               * m
         hrv_frootc_to_litter(p)              = frootc(p)              * m
         hrv_livestemc_to_litter(p)           = livestemc(p)           * m
         hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
         hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
         hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

         ! storage pools
         hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
         hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
         hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
         hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
         hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
         hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
         hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

         ! transfer pools
         hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
         hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
         hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
         hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
         hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
         hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
         hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

         ! pft-level harvest mortality nitrogen fluxes
         ! displayed pools
         hrv_leafn_to_litter(p)               = leafn(p)               * m
         hrv_frootn_to_litter(p)              = frootn(p)              * m
         hrv_livestemn_to_litter(p)           = livestemn(p)           * m
         hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * &
                                                pprodharv10(ivt(p))
         hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * &
                                                (1.0_r8 - pprodharv10(ivt(p)))
         hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
         hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
         hrv_retransn_to_litter(p)            = retransn(p)            * m

         ! storage pools
         hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
         hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
         hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
         hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
         hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
         hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

         ! transfer pools
         hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
         hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
         hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
         hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
         hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
         hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
         
      end if  ! end tree block

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes from harvest to the column
   ! for litter C and N inputs

   call CNHarvestPftToColumn(num_soilc, filter_soilc)

end subroutine CNHarvest
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvestPftToColumn
!
! !INTERFACE:
subroutine CNHarvestPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real(r8), pointer :: hrv_leafc_to_litter(:)
   real(r8), pointer :: hrv_frootc_to_litter(:)
   real(r8), pointer :: hrv_livestemc_to_litter(:)
   real(r8), pointer :: phrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: phrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
   real(r8), pointer :: hrv_leafn_to_litter(:)
   real(r8), pointer :: hrv_frootn_to_litter(:)
   real(r8), pointer :: hrv_livestemn_to_litter(:)
   real(r8), pointer :: phrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: phrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_livecrootn_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_to_litter(:)
   real(r8), pointer :: hrv_retransn_to_litter(:)
   real(r8), pointer :: hrv_leafn_storage_to_litter(:)
   real(r8), pointer :: hrv_frootn_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: hrv_leafc_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_to_litr2c(:)
   real(r8), pointer :: hrv_leafc_to_litr3c(:)
   real(r8), pointer :: hrv_frootc_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_to_litr2c(:)
   real(r8), pointer :: hrv_frootc_to_litr3c(:)
   real(r8), pointer :: hrv_livestemc_to_cwdc(:)
   real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_cwdc(:)
   real(r8), pointer :: hrv_deadcrootc_to_cwdc(:)
   real(r8), pointer :: hrv_leafc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_storage_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_leafn_to_litr1n(:)
   real(r8), pointer :: hrv_leafn_to_litr2n(:)
   real(r8), pointer :: hrv_leafn_to_litr3n(:)
   real(r8), pointer :: hrv_frootn_to_litr1n(:)
   real(r8), pointer :: hrv_frootn_to_litr2n(:)
   real(r8), pointer :: hrv_frootn_to_litr3n(:)
   real(r8), pointer :: hrv_livestemn_to_cwdn(:)
   real(r8), pointer :: chrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: chrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_livecrootn_to_cwdn(:)
   real(r8), pointer :: hrv_deadcrootn_to_cwdn(:)
   real(r8), pointer :: hrv_retransn_to_litr1n(:)
   real(r8), pointer :: hrv_leafn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_frootn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litr1n(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litr1n(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litr1n(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => col%npfts
   pfti                           => col%pfti
   hrv_leafc_to_litr1c              => ccf%hrv_leafc_to_litr1c
   hrv_leafc_to_litr2c              => ccf%hrv_leafc_to_litr2c
   hrv_leafc_to_litr3c              => ccf%hrv_leafc_to_litr3c
   hrv_frootc_to_litr1c             => ccf%hrv_frootc_to_litr1c
   hrv_frootc_to_litr2c             => ccf%hrv_frootc_to_litr2c
   hrv_frootc_to_litr3c             => ccf%hrv_frootc_to_litr3c
   hrv_livestemc_to_cwdc            => ccf%hrv_livestemc_to_cwdc
   chrv_deadstemc_to_prod10c        => ccf%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => ccf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_cwdc           => ccf%hrv_livecrootc_to_cwdc
   hrv_deadcrootc_to_cwdc           => ccf%hrv_deadcrootc_to_cwdc
   hrv_leafc_storage_to_litr1c      => ccf%hrv_leafc_storage_to_litr1c
   hrv_frootc_storage_to_litr1c     => ccf%hrv_frootc_storage_to_litr1c
   hrv_livestemc_storage_to_litr1c  => ccf%hrv_livestemc_storage_to_litr1c
   hrv_deadstemc_storage_to_litr1c  => ccf%hrv_deadstemc_storage_to_litr1c
   hrv_livecrootc_storage_to_litr1c => ccf%hrv_livecrootc_storage_to_litr1c
   hrv_deadcrootc_storage_to_litr1c => ccf%hrv_deadcrootc_storage_to_litr1c
   hrv_gresp_storage_to_litr1c      => ccf%hrv_gresp_storage_to_litr1c
   hrv_leafc_xfer_to_litr1c         => ccf%hrv_leafc_xfer_to_litr1c
   hrv_frootc_xfer_to_litr1c        => ccf%hrv_frootc_xfer_to_litr1c
   hrv_livestemc_xfer_to_litr1c     => ccf%hrv_livestemc_xfer_to_litr1c
   hrv_deadstemc_xfer_to_litr1c     => ccf%hrv_deadstemc_xfer_to_litr1c
   hrv_livecrootc_xfer_to_litr1c    => ccf%hrv_livecrootc_xfer_to_litr1c
   hrv_deadcrootc_xfer_to_litr1c    => ccf%hrv_deadcrootc_xfer_to_litr1c
   hrv_gresp_xfer_to_litr1c         => ccf%hrv_gresp_xfer_to_litr1c
   hrv_leafn_to_litr1n              => cnf%hrv_leafn_to_litr1n
   hrv_leafn_to_litr2n              => cnf%hrv_leafn_to_litr2n
   hrv_leafn_to_litr3n              => cnf%hrv_leafn_to_litr3n
   hrv_frootn_to_litr1n             => cnf%hrv_frootn_to_litr1n
   hrv_frootn_to_litr2n             => cnf%hrv_frootn_to_litr2n
   hrv_frootn_to_litr3n             => cnf%hrv_frootn_to_litr3n
   hrv_livestemn_to_cwdn            => cnf%hrv_livestemn_to_cwdn
   chrv_deadstemn_to_prod10n        => cnf%hrv_deadstemn_to_prod10n
   chrv_deadstemn_to_prod100n       => cnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_cwdn           => cnf%hrv_livecrootn_to_cwdn
   hrv_deadcrootn_to_cwdn           => cnf%hrv_deadcrootn_to_cwdn
   hrv_retransn_to_litr1n           => cnf%hrv_retransn_to_litr1n
   hrv_leafn_storage_to_litr1n      => cnf%hrv_leafn_storage_to_litr1n
   hrv_frootn_storage_to_litr1n     => cnf%hrv_frootn_storage_to_litr1n
   hrv_livestemn_storage_to_litr1n  => cnf%hrv_livestemn_storage_to_litr1n
   hrv_deadstemn_storage_to_litr1n  => cnf%hrv_deadstemn_storage_to_litr1n
   hrv_livecrootn_storage_to_litr1n => cnf%hrv_livecrootn_storage_to_litr1n
   hrv_deadcrootn_storage_to_litr1n => cnf%hrv_deadcrootn_storage_to_litr1n
   hrv_leafn_xfer_to_litr1n         => cnf%hrv_leafn_xfer_to_litr1n
   hrv_frootn_xfer_to_litr1n        => cnf%hrv_frootn_xfer_to_litr1n
   hrv_livestemn_xfer_to_litr1n     => cnf%hrv_livestemn_xfer_to_litr1n
   hrv_deadstemn_xfer_to_litr1n     => cnf%hrv_deadstemn_xfer_to_litr1n
   hrv_livecrootn_xfer_to_litr1n    => cnf%hrv_livecrootn_xfer_to_litr1n
   hrv_deadcrootn_xfer_to_litr1n    => cnf%hrv_deadcrootn_xfer_to_litr1n

   ! assign local pointers to pft-level arrays
   ivt                            => pft%itype
   wtcol                          => pft%wtcol
   pwtgcell                       => pft%wtgcell  
   hrv_leafc_to_litter              => pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => pcf%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => pcf%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => pcf%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => pnf%hrv_livestemn_to_litter
   phrv_deadstemn_to_prod10n        => pnf%hrv_deadstemn_to_prod10n
   phrv_deadstemn_to_prod100n       => pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => pnf%hrv_deadcrootn_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1

            if (pwtgcell(p)>0._r8) then

               ! leaf harvest mortality carbon fluxes
               hrv_leafc_to_litr1c(c) = hrv_leafc_to_litr1c(c) + &
                  hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr2c(c) = hrv_leafc_to_litr2c(c) + &
                  hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr3c(c) = hrv_leafc_to_litr3c(c) + &
                  hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root harvest mortality carbon fluxes
               hrv_frootc_to_litr1c(c) = hrv_frootc_to_litr1c(c) + &
                  hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr2c(c) = hrv_frootc_to_litr2c(c) + &
                  hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr3c(c) = hrv_frootc_to_litr3c(c) + &
                  hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality carbon fluxes
               hrv_livestemc_to_cwdc(c)  = hrv_livestemc_to_cwdc(c)  + &
                  hrv_livestemc_to_litter(p)  * wtcol(p)
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                  phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                  phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               hrv_livecrootc_to_cwdc(c) = hrv_livecrootc_to_cwdc(c) + &
                  hrv_livecrootc_to_litter(p) * wtcol(p)
               hrv_deadcrootc_to_cwdc(c) = hrv_deadcrootc_to_cwdc(c) + &
                  hrv_deadcrootc_to_litter(p) * wtcol(p)

               ! storage harvest mortality carbon fluxes
               hrv_leafc_storage_to_litr1c(c)      = hrv_leafc_storage_to_litr1c(c)      + &
                  hrv_leafc_storage_to_litter(p)      * wtcol(p)
               hrv_frootc_storage_to_litr1c(c)     = hrv_frootc_storage_to_litr1c(c)     + &
                  hrv_frootc_storage_to_litter(p)     * wtcol(p)
               hrv_livestemc_storage_to_litr1c(c)  = hrv_livestemc_storage_to_litr1c(c)  + &
                  hrv_livestemc_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemc_storage_to_litr1c(c)  = hrv_deadstemc_storage_to_litr1c(c)  + &
                  hrv_deadstemc_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootc_storage_to_litr1c(c) = hrv_livecrootc_storage_to_litr1c(c) + &
                  hrv_livecrootc_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootc_storage_to_litr1c(c) = hrv_deadcrootc_storage_to_litr1c(c) + &
                  hrv_deadcrootc_storage_to_litter(p) * wtcol(p)
               hrv_gresp_storage_to_litr1c(c)      = hrv_gresp_storage_to_litr1c(c)      + &
                  hrv_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer harvest mortality carbon fluxes
               hrv_leafc_xfer_to_litr1c(c)      = hrv_leafc_xfer_to_litr1c(c)      + &
                  hrv_leafc_xfer_to_litter(p)      * wtcol(p)
               hrv_frootc_xfer_to_litr1c(c)     = hrv_frootc_xfer_to_litr1c(c)     + &
                  hrv_frootc_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemc_xfer_to_litr1c(c)  = hrv_livestemc_xfer_to_litr1c(c)  + &
                  hrv_livestemc_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemc_xfer_to_litr1c(c)  = hrv_deadstemc_xfer_to_litr1c(c)  + &
                  hrv_deadstemc_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootc_xfer_to_litr1c(c) = hrv_livecrootc_xfer_to_litr1c(c) + &
                  hrv_livecrootc_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootc_xfer_to_litr1c(c) = hrv_deadcrootc_xfer_to_litr1c(c) + &
                  hrv_deadcrootc_xfer_to_litter(p) * wtcol(p)
               hrv_gresp_xfer_to_litr1c(c)      = hrv_gresp_xfer_to_litr1c(c)      + &
                  hrv_gresp_xfer_to_litter(p)      * wtcol(p)

               ! leaf harvest mortality nitrogen fluxes
               hrv_leafn_to_litr1n(c) = hrv_leafn_to_litr1n(c) + &
                  hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafn_to_litr2n(c) = hrv_leafn_to_litr2n(c) + &
                  hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafn_to_litr3n(c) = hrv_leafn_to_litr3n(c) + &
                  hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               hrv_frootn_to_litr1n(c) = hrv_frootn_to_litr1n(c) + &
                  hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootn_to_litr2n(c) = hrv_frootn_to_litr2n(c) + &
                  hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootn_to_litr3n(c) = hrv_frootn_to_litr3n(c) + &
                  hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality nitrogen fluxes
               hrv_livestemn_to_cwdn(c)  = hrv_livestemn_to_cwdn(c)  + &
                  hrv_livestemn_to_litter(p)  * wtcol(p)
               chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
                  phrv_deadstemn_to_prod10n(p)  * wtcol(p)
               chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
                  phrv_deadstemn_to_prod100n(p)  * wtcol(p)
               hrv_livecrootn_to_cwdn(c) = hrv_livecrootn_to_cwdn(c) + &
                  hrv_livecrootn_to_litter(p) * wtcol(p)
               hrv_deadcrootn_to_cwdn(c) = hrv_deadcrootn_to_cwdn(c) + &
                  hrv_deadcrootn_to_litter(p) * wtcol(p)

               ! retranslocated N pool harvest mortality fluxes
               hrv_retransn_to_litr1n(c) = hrv_retransn_to_litr1n(c) + &
                  hrv_retransn_to_litter(p) * wtcol(p)

               ! storage harvest mortality nitrogen fluxes
               hrv_leafn_storage_to_litr1n(c)      = hrv_leafn_storage_to_litr1n(c)      + &
                  hrv_leafn_storage_to_litter(p)      * wtcol(p)
               hrv_frootn_storage_to_litr1n(c)     = hrv_frootn_storage_to_litr1n(c)     + &
                  hrv_frootn_storage_to_litter(p)     * wtcol(p)
               hrv_livestemn_storage_to_litr1n(c)  = hrv_livestemn_storage_to_litr1n(c)  + &
                  hrv_livestemn_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemn_storage_to_litr1n(c)  = hrv_deadstemn_storage_to_litr1n(c)  + &
                  hrv_deadstemn_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootn_storage_to_litr1n(c) = hrv_livecrootn_storage_to_litr1n(c) + &
                  hrv_livecrootn_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootn_storage_to_litr1n(c) = hrv_deadcrootn_storage_to_litr1n(c) + &
                  hrv_deadcrootn_storage_to_litter(p) * wtcol(p)

               ! transfer harvest mortality nitrogen fluxes
               hrv_leafn_xfer_to_litr1n(c)      = hrv_leafn_xfer_to_litr1n(c)      + &
                  hrv_leafn_xfer_to_litter(p)      * wtcol(p)
               hrv_frootn_xfer_to_litr1n(c)     = hrv_frootn_xfer_to_litr1n(c)     + &
                  hrv_frootn_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemn_xfer_to_litr1n(c)  = hrv_livestemn_xfer_to_litr1n(c)  + &
                  hrv_livestemn_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemn_xfer_to_litr1n(c)  = hrv_deadstemn_xfer_to_litr1n(c)  + &
                  hrv_deadstemn_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootn_xfer_to_litr1n(c) = hrv_livecrootn_xfer_to_litr1n(c) + &
                  hrv_livecrootn_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootn_xfer_to_litr1n(c) = hrv_deadcrootn_xfer_to_litr1n(c) + &
                  hrv_deadcrootn_xfer_to_litter(p) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNHarvestPftToColumn
!-----------------------------------------------------------------------

end module pftdynMod
