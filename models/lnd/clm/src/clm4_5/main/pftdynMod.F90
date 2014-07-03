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
  use clm_varctl  , only : iulog, use_c13, use_c14
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
#ifdef CN
  public :: pftdyn_cnbal
#ifdef CNDV
  public :: pftwt_init
  public :: pftwt_interp
#endif
  public :: CNHarvest
  public :: CNHarvestPftToColumn
#endif
!
! !REVISION HISTORY:
! Created by Peter Thornton
! slevis modified to handle CNDV and crop model
! 19 May 2009: PET - modified to handle harvest fluxes
! F. Li and S. Levis (11/06/12)

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
  ! default multiplication factor for epsilon for error checks
  real(r8), private, parameter :: eps_fact = 2._r8
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
    use clm_varpar  , only : numpft, maxpatch_pft, numurbl
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i,j,m,n,g,nl                    ! indices
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
    real(r8), pointer :: pcturb(:,:)        ! percent of gcell is urbanized
    real(r8), pointer :: pcturb_tot(:)      ! percent of grid cell is urban (sum over density classes)
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
    allocate(pctwet(begg:endg),pcturb(begg:endg,numurbl),pcturb_tot(begg:endg))

    ! Set pointers into derived type

    gptr => clm3%g

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
    pcturb_tot(:) = 0._r8
    do n = 1, numurbl
       do nl = begg,endg
          pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
       enddo
    enddo

    ! Consistency check
    do g = begg,endg
    !   this was causing a fail, even though values are the same to within 1e-15
    !   if (pctlak(g)+pctwet(g)+pcturb(g)+pctgla(g) /= pctspec(g)) then 
       if (abs((pctlak(g)+pctwet(g)+pcturb_tot(g)+pctgla(g))-pctspec(g)) > 1e-13_r8) then 
          write(iulog,*) subname//'mismatch between input pctspec = ',&
                     pctlak(g)+pctwet(g)+pcturb_tot(g)+pctgla(g),&
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

    call pftdyn_getdata(nt1, wtpft1, begg,endg,0,numpft)
    call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)
    
#ifdef CN
    ! Get harvest rate at the nt1 time
    call pftdyn_getharvest(nt1,begg,endg)
#endif

    ! convert weights from percent to proportion
    do m = 0,numpft
       do g = begg,endg
          wtpft1(g,m) = wtpft1(g,m)/100._r8
          wtpft2(g,m) = wtpft2(g,m)/100._r8
       end do
    end do
       
    deallocate(pctgla,pctlak,pctwet,pcturb,pcturb_tot)

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
    type(column_type),   pointer :: cptr         ! F. Li and S. Levis
    character(len=32) :: subname='pftdyn_interp' ! subroutine name
    logical  :: readvar    ! F. Li and S. Levis
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Set pointers into derived type
    gptr => clm3%g
    lptr => clm3%g%l
    pptr => clm3%g%l%c%p

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

       call pftdyn_getdata(nt2, wtpft2, begg,endg,0,numpft)

#ifdef CN
       call pftdyn_getharvest(nt1,begg,endg)
#endif

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
       ! THESE CHECKS NEEDS TO BE THE SAME AS IN surfrdMod.F90!
       if (pctspec(n) < 100._r8 * (1._r8 - eps_fact*epsilon(1._r8))) then  ! pctspec not within eps_fact*epsilon of 100
          sumpct = 0._r8
          do m = 0, numpft
             sumpct = sumpct + pctpft(n,m) * 100._r8/(100._r8-pctspec(n))
          end do
          if (abs(sumpct - 100._r8) > 0.1e-4_r8) then
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

#ifdef CN
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
#endif

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

    cptr => clm3%g%l%c

    ! set column-level canopy water mass balance correction flux
    ! term to 0 at the beginning of every timestep
    
    do c = begc,endc
       cptr%cwf%h2ocan_loss(c) = 0._r8
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

    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

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
       cptr%cwf%h2ocan_loss(c) = 0._r8 ! is this OR pftdyn_wbal_init redundant?
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

             pptr%pws%h2ocan(p) = pptr%pws%h2ocan(p) * (wtcol_old(p)/pptr%wtcol(p))
          
          else if (dwt < 0._r8) then
          
             ! if the pft lost weight on the timestep, then the canopy water
             ! mass associated with the lost weight is directed to a 
             ! column-level flux term that gets added to the precip flux
             ! for every pft calculation in Hydrology1()
             
             init_h2ocan = pptr%pws%h2ocan(p) * wtcol_old(p)
             loss_h2ocan(p) = pptr%pws%h2ocan(p) * (-dwt)
             new_h2ocan = init_h2ocan - loss_h2ocan(p)
             if (abs(new_h2ocan) < 1e-8_r8) then
                new_h2ocan = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if
             if (pptr%wtcol(p) /= 0._r8) then  
                pptr%pws%h2ocan(p) = new_h2ocan/pptr%wtcol(p)
             else
                pptr%pws%h2ocan(p) = 0._r8
                loss_h2ocan(p) = init_h2ocan
             end if 
       

          end if

       end if
    end do

    do pi = 1,max_pft_per_col
       do c = begc,endc
          if (pi <= cptr%npfts(c)) then
             p = cptr%pfti(c) + pi - 1
             cptr%cwf%h2ocan_loss(c) = cptr%cwf%h2ocan_loss(c) + loss_h2ocan(p)/dtime
          end if
       end do
    end do

    ! Deallocate loss_h2ocan
    deallocate(loss_h2ocan)
    
  end subroutine pftdyn_wbal
  
#ifdef CN
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
    use clm_varpar  , only : numveg, numpft, nlevdecomp
    use clm_varcon  , only : istcrop
    use pftvarcon   , only : pconv, pprod10, pprod100
    use clm_varcon  , only : c13ratio, c14ratio
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
    integer  :: pi,p,c,l,g,j    ! indices
    integer  :: ier           ! error code
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: dt            ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: dwt_leafc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: conv_cflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod10_cflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod100_cflux(:)      ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)      ! pft-level mass loss due to weight shift
    real(r8) :: t1,t2,wt_new,wt_old
    real(r8) :: init_state, change_state, new_state
    real(r8) :: tot_leaf, pleaf, pstor, pxfer
    real(r8) :: leafc_seed, leafn_seed
    real(r8) :: deadstemc_seed, deadstemn_seed
    real(r8), pointer :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    integer , pointer :: ivt(:)    ! pft vegetation type added by F. Li and S. Levis
    real(r8),   pointer :: lfpftd(:)       !F. Li and S. Levis
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type),   pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)   ,   pointer :: pptr         ! pointer to pft derived subtype
    integer          , pointer   :: pcolumn(:)   ! column of corresponding pft
    character(len=32) :: subname='pftdyn_cbal' ! subroutine name
    !!! C13
    real(r8), allocatable :: dwt_leafc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)    ! pft-level mass loss due to weight shift
    real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1_c13         ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1_c13         ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2_c13         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2_c13         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: leafc13_seed, deadstemc13_seed
    !!! C14
    real(r8), allocatable :: dwt_leafc14_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc14_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c14flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c14flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c14flux(:)    ! pft-level mass loss due to weight shift
    real(r8) :: c3_del14c     ! typical del14C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del14c     ! typical del14C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1_c14         ! isotope ratio (14c/12c) for C3 photosynthesis
    real(r8) :: c4_r1_c14         ! isotope ratio (14c/12c) for C4 photosynthesis
    real(r8) :: c3_r2_c14         ! isotope ratio (14c/[12c+14c]) for C3 photosynthesis
    real(r8) :: c4_r2_c14         ! isotope ratio (14c/[12c+14c]) for C4 photosynthesis
    real(r8) :: leafc14_seed, deadstemc14_seed
   
!-----------------------------------------------------------------------
    
    ! Set pointers into derived type

    lptr    => clm3%g%l
    cptr    => clm3%g%l%c
    pptr    => clm3%g%l%c%p
    pcolumn => pptr%column
    lfpftd  => clm3%g%l%c%p%pps%lfpftd     ! F. Li and S. Levis
    ivt     => clm3%g%l%c%p%itype           ! F. Li and S. Levis

    ! Allocate pft-level mass loss arrays
    allocate(dwt_leafc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc_seed'; call endrun()
    end if
    allocate(dwt_leafn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafn_seed'; call endrun()
    end if
    allocate(dwt_deadstemc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'; call endrun()
    end if
    allocate(dwt_deadstemn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'; call endrun()
    end if
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

    if ( use_c13 ) then
       allocate(dwt_leafc13_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc13_seed'; call endrun()
       end if
       allocate(dwt_deadstemc13_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'; call endrun()
       end if
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
    if ( use_c14 ) then
       allocate(dwt_leafc14_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc14_seed'; call endrun()
       end if
       allocate(dwt_deadstemc14_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc14_seed'; call endrun()
       end if
       allocate(dwt_frootc14_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc14_to_litter'; call endrun()
       end if
       allocate(dwt_livecrootc14_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc14_to_litter'; call endrun()
       end if
       allocate(dwt_deadcrootc14_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc14_to_litter'; call endrun()
       end if
       allocate(conv_c14flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c14flux'; call endrun()
       end if
       allocate(prod10_c14flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c14flux'; call endrun()
       end if
       allocate(prod100_c14flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c14flux'; call endrun()
       end if
    endif
    
    ! Get time step
    dt = real( get_step_size(), r8 )
    
    do p = begp,endp
       c = pcolumn(p)
       ! initialize all the pft-level local flux arrays
       dwt_leafc_seed(p) = 0._r8
       dwt_leafn_seed(p) = 0._r8
       dwt_deadstemc_seed(p) = 0._r8
       dwt_deadstemn_seed(p) = 0._r8
       dwt_frootc_to_litter(p) = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       dwt_frootn_to_litter(p) = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_cflux(p) = 0._r8
       prod10_cflux(p) = 0._r8
       prod100_cflux(p) = 0._r8
       conv_nflux(p) = 0._r8
       prod10_nflux(p) = 0._r8
       prod100_nflux(p) = 0._r8
       
       if ( use_c13 ) then
          dwt_leafc13_seed(p) = 0._r8
          dwt_deadstemc13_seed(p) = 0._r8
          dwt_frootc13_to_litter(p) = 0._r8
          dwt_livecrootc13_to_litter(p) = 0._r8
          dwt_deadcrootc13_to_litter(p) = 0._r8
          conv_c13flux(p) = 0._r8
          prod10_c13flux(p) = 0._r8
          prod100_c13flux(p) = 0._r8
       endif
       
       if ( use_c14 ) then
          dwt_leafc14_seed(p) = 0._r8
          dwt_deadstemc14_seed(p) = 0._r8
          dwt_frootc14_to_litter(p) = 0._r8
          dwt_livecrootc14_to_litter(p) = 0._r8
          dwt_deadcrootc14_to_litter(p) = 0._r8
          conv_c14flux(p) = 0._r8
          prod10_c14flux(p) = 0._r8
          prod100_c14flux(p) = 0._r8
       endif
       
       l = pptr%landunit(p)
       if (lptr%itype(l) == istsoil .or. lptr%itype(l) == istcrop) then
          
          ! calculate the change in weight for the timestep
          dwt = pptr%wtcol(p)-wtcol_old(p)
            lfpftd(p)=-dwt
          ! PFTs for which weight increases on this timestep
          if (dwt > 0._r8) then
             
             ! first identify PFTs that are initiating on this timestep
             ! and set all the necessary state and flux variables
             if (wtcol_old(p) == 0._r8) then
                
                ! set initial conditions for PFT that is being initiated
                ! in this time step.  Based on the settings in cnIniTimeVar.
                
                ! pft-level carbon state variables
                pptr%pcs%leafc(p)              = 0._r8
                pptr%pcs%leafc_storage(p)      = 0._r8
                pptr%pcs%leafc_xfer(p)         = 0._r8
                pptr%pcs%frootc(p)             = 0._r8
                pptr%pcs%frootc_storage(p)     = 0._r8
                pptr%pcs%frootc_xfer(p)        = 0._r8
                pptr%pcs%livestemc(p)          = 0._r8
                pptr%pcs%livestemc_storage(p)  = 0._r8
                pptr%pcs%livestemc_xfer(p)     = 0._r8
                pptr%pcs%deadstemc(p)          = 0._r8
                pptr%pcs%deadstemc_storage(p)  = 0._r8
                pptr%pcs%deadstemc_xfer(p)     = 0._r8
                pptr%pcs%livecrootc(p)         = 0._r8
                pptr%pcs%livecrootc_storage(p) = 0._r8
                pptr%pcs%livecrootc_xfer(p)    = 0._r8
                pptr%pcs%deadcrootc(p)         = 0._r8
                pptr%pcs%deadcrootc_storage(p) = 0._r8
                pptr%pcs%deadcrootc_xfer(p)    = 0._r8
                pptr%pcs%gresp_storage(p)      = 0._r8
                pptr%pcs%gresp_xfer(p)         = 0._r8
                pptr%pcs%cpool(p)              = 0._r8
                pptr%pcs%xsmrpool(p)           = 0._r8
                pptr%pcs%pft_ctrunc(p)         = 0._r8
                pptr%pcs%dispvegc(p)           = 0._r8
                pptr%pcs%storvegc(p)           = 0._r8
                pptr%pcs%totvegc(p)            = 0._r8
                pptr%pcs%totpftc(p)            = 0._r8
                
                if ( use_c13 ) then
                   ! pft-level carbon-13 state variables
                   pptr%pc13s%leafc(p)              = 0._r8
                   pptr%pc13s%leafc_storage(p)      = 0._r8
                   pptr%pc13s%leafc_xfer(p)         = 0._r8
                   pptr%pc13s%frootc(p)             = 0._r8
                   pptr%pc13s%frootc_storage(p)     = 0._r8
                   pptr%pc13s%frootc_xfer(p)        = 0._r8
                   pptr%pc13s%livestemc(p)          = 0._r8
                   pptr%pc13s%livestemc_storage(p)  = 0._r8
                   pptr%pc13s%livestemc_xfer(p)     = 0._r8
                   pptr%pc13s%deadstemc(p)          = 0._r8
                   pptr%pc13s%deadstemc_storage(p)  = 0._r8
                   pptr%pc13s%deadstemc_xfer(p)     = 0._r8
                   pptr%pc13s%livecrootc(p)         = 0._r8
                   pptr%pc13s%livecrootc_storage(p) = 0._r8
                   pptr%pc13s%livecrootc_xfer(p)    = 0._r8
                   pptr%pc13s%deadcrootc(p)         = 0._r8
                   pptr%pc13s%deadcrootc_storage(p) = 0._r8
                   pptr%pc13s%deadcrootc_xfer(p)    = 0._r8
                   pptr%pc13s%gresp_storage(p)      = 0._r8
                   pptr%pc13s%gresp_xfer(p)         = 0._r8
                   pptr%pc13s%cpool(p)              = 0._r8
                   pptr%pc13s%xsmrpool(p)           = 0._r8
                   pptr%pc13s%pft_ctrunc(p)         = 0._r8
                   pptr%pc13s%dispvegc(p)           = 0._r8
                   pptr%pc13s%storvegc(p)           = 0._r8
                   pptr%pc13s%totvegc(p)            = 0._r8
                   pptr%pc13s%totpftc(p)            = 0._r8
                endif
                
                if ( use_c14 ) then
                   ! pft-level carbon-14 state variables
                   pptr%pc14s%leafc(p)              = 0._r8
                   pptr%pc14s%leafc_storage(p)      = 0._r8
                   pptr%pc14s%leafc_xfer(p)         = 0._r8
                   pptr%pc14s%frootc(p)             = 0._r8
                   pptr%pc14s%frootc_storage(p)     = 0._r8
                   pptr%pc14s%frootc_xfer(p)        = 0._r8
                   pptr%pc14s%livestemc(p)          = 0._r8
                   pptr%pc14s%livestemc_storage(p)  = 0._r8
                   pptr%pc14s%livestemc_xfer(p)     = 0._r8
                   pptr%pc14s%deadstemc(p)          = 0._r8
                   pptr%pc14s%deadstemc_storage(p)  = 0._r8
                   pptr%pc14s%deadstemc_xfer(p)     = 0._r8
                   pptr%pc14s%livecrootc(p)         = 0._r8
                   pptr%pc14s%livecrootc_storage(p) = 0._r8
                   pptr%pc14s%livecrootc_xfer(p)    = 0._r8
                   pptr%pc14s%deadcrootc(p)         = 0._r8
                   pptr%pc14s%deadcrootc_storage(p) = 0._r8
                   pptr%pc14s%deadcrootc_xfer(p)    = 0._r8
                   pptr%pc14s%gresp_storage(p)      = 0._r8
                   pptr%pc14s%gresp_xfer(p)         = 0._r8
                   pptr%pc14s%cpool(p)              = 0._r8
                   pptr%pc14s%xsmrpool(p)           = 0._r8
                   pptr%pc14s%pft_ctrunc(p)         = 0._r8
                   pptr%pc14s%dispvegc(p)           = 0._r8
                   pptr%pc14s%storvegc(p)           = 0._r8
                   pptr%pc14s%totvegc(p)            = 0._r8
                   pptr%pc14s%totpftc(p)            = 0._r8
                endif
                
                ! pft-level nitrogen state variables
                pptr%pns%leafn(p)	           = 0._r8
                pptr%pns%leafn_storage(p)      = 0._r8
                pptr%pns%leafn_xfer(p)         = 0._r8
                pptr%pns%frootn(p)	           = 0._r8
                pptr%pns%frootn_storage(p)     = 0._r8
                pptr%pns%frootn_xfer(p)        = 0._r8
                pptr%pns%livestemn(p)	       = 0._r8
                pptr%pns%livestemn_storage(p)  = 0._r8
                pptr%pns%livestemn_xfer(p)     = 0._r8
                pptr%pns%deadstemn(p)	       = 0._r8
                pptr%pns%deadstemn_storage(p)  = 0._r8
                pptr%pns%deadstemn_xfer(p)     = 0._r8
                pptr%pns%livecrootn(p)         = 0._r8
                pptr%pns%livecrootn_storage(p) = 0._r8
                pptr%pns%livecrootn_xfer(p)    = 0._r8
                pptr%pns%deadcrootn(p)         = 0._r8
                pptr%pns%deadcrootn_storage(p) = 0._r8
                pptr%pns%deadcrootn_xfer(p)    = 0._r8
                pptr%pns%retransn(p)	       = 0._r8
                pptr%pns%npool(p)	           = 0._r8
                pptr%pns%pft_ntrunc(p)         = 0._r8
                pptr%pns%dispvegn(p)           = 0._r8
                pptr%pns%storvegn(p)           = 0._r8
                pptr%pns%totvegn(p)            = 0._r8
                pptr%pns%totpftn (p)           = 0._r8
                
                ! initialize same flux and epv variables that are set
                ! in CNiniTimeVar
                pptr%pcf%psnsun(p) = 0._r8
                pptr%pcf%psnsha(p) = 0._r8
                pptr%pps%laisun(p) = 0._r8
                pptr%pps%laisha(p) = 0._r8
                
                pptr%pepv%dormant_flag(p) = 1._r8
                pptr%pepv%days_active(p) = 0._r8
                pptr%pepv%onset_flag(p) = 0._r8
                pptr%pepv%onset_counter(p) = 0._r8
                pptr%pepv%onset_gddflag(p) = 0._r8
                pptr%pepv%onset_fdd(p) = 0._r8
                pptr%pepv%onset_gdd(p) = 0._r8
                pptr%pepv%onset_swi(p) = 0.0_r8
                pptr%pepv%offset_flag(p) = 0._r8
                pptr%pepv%offset_counter(p) = 0._r8
                pptr%pepv%offset_fdd(p) = 0._r8
                pptr%pepv%offset_swi(p) = 0._r8
                pptr%pepv%lgsf(p) = 0._r8
                pptr%pepv%bglfr(p) = 0._r8
                pptr%pepv%bgtr(p) = 0._r8
                ! difference from CNiniTimeVar: using column-level
                ! information to initialize annavg_t2m.
                pptr%pepv%annavg_t2m(p) = cptr%cps%cannavg_t2m(c)
                pptr%pepv%tempavg_t2m(p) = 0._r8
                pptr%pepv%gpp(p) = 0._r8
                pptr%pepv%availc(p) = 0._r8
                pptr%pepv%xsmrpool_recover(p) = 0._r8
                pptr%pepv%alloc_pnow(p) = 1._r8
                pptr%pepv%c_allometry(p) = 0._r8
                pptr%pepv%n_allometry(p) = 0._r8
                pptr%pepv%plant_ndemand(p) = 0._r8
                pptr%pepv%tempsum_potential_gpp(p) = 0._r8
                pptr%pepv%annsum_potential_gpp(p) = 0._r8
                pptr%pepv%tempmax_retransn(p) = 0._r8
                pptr%pepv%annmax_retransn(p) = 0._r8
                pptr%pepv%avail_retransn(p) = 0._r8
                pptr%pepv%plant_nalloc(p) = 0._r8
                pptr%pepv%plant_calloc(p) = 0._r8
                pptr%pepv%excess_cflux(p) = 0._r8
                pptr%pepv%downreg(p) = 0._r8
                pptr%pepv%prev_leafc_to_litter(p) = 0._r8
                pptr%pepv%prev_frootc_to_litter(p) = 0._r8
                pptr%pepv%tempsum_npp(p) = 0._r8
                pptr%pepv%annsum_npp(p) = 0._r8
                
                if ( use_c13 ) then
                   pptr%pc13f%psnsun(p) = 0._r8
                   pptr%pc13f%psnsha(p) = 0._r8
                   
                   pptr%pps%alphapsnsun(p) = 0._r8
                   pptr%pps%alphapsnsha(p) = 0._r8
                   
                   pptr%pepv%xsmrpool_c13ratio(p) = c13ratio
                   
                   pptr%pepv%rc13_canair(p) = 0._r8
                   pptr%pepv%rc13_psnsun(p) = 0._r8
                   pptr%pepv%rc13_psnsha(p) = 0._r8
                endif
                
                if ( use_c14 ) then
                   pptr%pc14f%psnsun(p) = 0._r8
                   pptr%pc14f%psnsha(p) = 0._r8
                   
                   pptr%pepv%rc14_atm(p) = c14ratio
                   
                   ! pptr%pps%alphapsnsun(p) = 0._r8
                   ! pptr%pps%alphapsnsha(p) = 0._r8
                   
                   ! pptr%pepv%xsmrpool_c14ratio(p) = c14ratio
                   
                   pptr%pepv%rc14_atm(p) = 0._r8
                   
                   ! pptr%pepv%rc14_canair(p) = 0._r8
                   ! pptr%pepv%rc14_psnsun(p) = 0._r8
                   ! pptr%pepv%rc14_psnsha(p) = 0._r8
                endif
                
             end if  ! end initialization of new pft
             
             ! (still in dwt > 0 block)
             
             ! set the seed sources for leaf and deadstem
             ! leaf source is split later between leaf, leaf_storage, leaf_xfer
             leafc_seed   = 0._r8
             leafn_seed   = 0._r8
             deadstemc_seed   = 0._r8
             deadstemn_seed   = 0._r8
             if ( use_c13 ) then
                leafc13_seed = 0._r8
                deadstemc13_seed = 0._r8
             endif
             if ( use_c14 ) then
                leafc14_seed = 0._r8
                deadstemc14_seed = 0._r8
             endif
             if (pptr%itype(p) /= 0) then
                leafc_seed = 1._r8
                leafn_seed  = leafc_seed / pftcon%leafcn(pptr%itype(p))
                if (pftcon%woody(pptr%itype(p)) == 1._r8) then
                   deadstemc_seed = 0.1_r8
                   deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(pptr%itype(p))
                end if
                
                if ( use_c13 ) then
                   ! 13c state is initialized assuming del13c = -28 permil for C3, and -13 permil for C4.
                   ! That translates to ratios of (13c/(12c+13c)) of 0.01080455 for C3, and 0.01096945 for C4
                   ! based on the following formulae: 
                   ! r1 (13/12) = PDB + (del13c * PDB)/1000.0
                   ! r2 (13/(13+12)) = r1/(1+r1)
                   ! PDB = 0.0112372_R8  (ratio of 13C/12C in Pee Dee Belemnite, C isotope standard)
                   c3_del13c = -28._r8
                   c4_del13c = -13._r8
                   c3_r1_c13 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
                   c3_r2_c13 = c3_r1_c13/(1._r8 + c3_r1_c13)
                   c4_r1_c13 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
                   c4_r2_c13 = c4_r1_c13/(1._r8 + c4_r1_c13)
                   
                   if (pftcon%c3psn(pptr%itype(p)) == 1._r8) then
                      leafc13_seed     = leafc_seed     * c3_r2_c13
                      deadstemc13_seed = deadstemc_seed * c3_r2_c13
                   else
                      leafc13_seed     = leafc_seed     * c4_r2_c13
                      deadstemc13_seed = deadstemc_seed * c4_r2_c13
                   end if
                endif
                
                if ( use_c14 ) then
                   ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
                   if (pftcon%c3psn(pptr%itype(p)) == 1._r8) then
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   else
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   end if
                endif
             end if
             
             ! When PFT area expands (dwt > 0), the pft-level mass density 
             ! is modified to conserve the original pft mass distributed
             ! over the new (larger) area, plus a term to account for the 
             ! introduction of new seed source for leaf and deadstem
             t1 = wtcol_old(p)/pptr%wtcol(p)
             t2 = dwt/pptr%wtcol(p)
             
             tot_leaf = pptr%pcs%leafc(p) + pptr%pcs%leafc_storage(p) + pptr%pcs%leafc_xfer(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                ! when adding seed source to non-zero leaf state, use current proportions
                pleaf = pptr%pcs%leafc(p)/tot_leaf
                pstor = pptr%pcs%leafc_storage(p)/tot_leaf
                pxfer = pptr%pcs%leafc_xfer(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             pptr%pcs%leafc(p)         = pptr%pcs%leafc(p)*t1         + leafc_seed*pleaf*t2
             pptr%pcs%leafc_storage(p) = pptr%pcs%leafc_storage(p)*t1 + leafc_seed*pstor*t2
             pptr%pcs%leafc_xfer(p)    = pptr%pcs%leafc_xfer(p)*t1    + leafc_seed*pxfer*t2
             pptr%pcs%frootc(p)  		   = pptr%pcs%frootc(p) 			* t1
             pptr%pcs%frootc_storage(p)     = pptr%pcs%frootc_storage(p) 	* t1
             pptr%pcs%frootc_xfer(p) 	   = pptr%pcs%frootc_xfer(p)		* t1
             pptr%pcs%livestemc(p)		   = pptr%pcs%livestemc(p)  		* t1
             pptr%pcs%livestemc_storage(p)  = pptr%pcs%livestemc_storage(p)  * t1
             pptr%pcs%livestemc_xfer(p)     = pptr%pcs%livestemc_xfer(p) 	* t1
             pptr%pcs%deadstemc(p)     = pptr%pcs%deadstemc(p)*t1     + deadstemc_seed*t2
             pptr%pcs%deadstemc_storage(p)  = pptr%pcs%deadstemc_storage(p)  * t1
             pptr%pcs%deadstemc_xfer(p)     = pptr%pcs%deadstemc_xfer(p) 	* t1
             pptr%pcs%livecrootc(p)  	   = pptr%pcs%livecrootc(p) 		* t1
             pptr%pcs%livecrootc_storage(p) = pptr%pcs%livecrootc_storage(p) * t1
             pptr%pcs%livecrootc_xfer(p)    = pptr%pcs%livecrootc_xfer(p)	* t1
             pptr%pcs%deadcrootc(p)  	   = pptr%pcs%deadcrootc(p) 		* t1
             pptr%pcs%deadcrootc_storage(p) = pptr%pcs%deadcrootc_storage(p) * t1
             pptr%pcs%deadcrootc_xfer(p)    = pptr%pcs%deadcrootc_xfer(p)	* t1
             pptr%pcs%gresp_storage(p)	   = pptr%pcs%gresp_storage(p)  	* t1
             pptr%pcs%gresp_xfer(p)  	   = pptr%pcs%gresp_xfer(p) 		* t1
             pptr%pcs%cpool(p)			   = pptr%pcs%cpool(p)  			* t1
             pptr%pcs%xsmrpool(p)		   = pptr%pcs%xsmrpool(p)			* t1
             pptr%pcs%pft_ctrunc(p)  	   = pptr%pcs%pft_ctrunc(p) 		* t1
             pptr%pcs%dispvegc(p)		   = pptr%pcs%dispvegc(p)			* t1
             pptr%pcs%storvegc(p)		   = pptr%pcs%storvegc(p)			* t1
             pptr%pcs%totvegc(p) 		   = pptr%pcs%totvegc(p)			* t1
             pptr%pcs%totpftc(p) 		   = pptr%pcs%totpftc(p)			* t1
             
             if ( use_c13 ) then
                ! pft-level carbon-13 state variables 
                tot_leaf = pptr%pc13s%leafc(p) + pptr%pc13s%leafc_storage(p) + pptr%pc13s%leafc_xfer(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = pptr%pc13s%leafc(p)/tot_leaf
                   pstor = pptr%pc13s%leafc_storage(p)/tot_leaf
                   pxfer = pptr%pc13s%leafc_xfer(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                pptr%pc13s%leafc(p)         = pptr%pc13s%leafc(p)*t1         + leafc13_seed*pleaf*t2
                pptr%pc13s%leafc_storage(p) = pptr%pc13s%leafc_storage(p)*t1 + leafc13_seed*pstor*t2
                pptr%pc13s%leafc_xfer(p)    = pptr%pc13s%leafc_xfer(p)*t1    + leafc13_seed*pxfer*t2
                pptr%pc13s%frootc(p)			 = pptr%pc13s%frootc(p) 		* t1
                pptr%pc13s%frootc_storage(p)	         = pptr%pc13s%frootc_storage(p) 	* t1
                pptr%pc13s%frootc_xfer(p)		 = pptr%pc13s%frootc_xfer(p)		* t1
                pptr%pc13s%livestemc(p) 		 = pptr%pc13s%livestemc(p)  		* t1
                pptr%pc13s%livestemc_storage(p)          = pptr%pc13s%livestemc_storage(p)      * t1
                pptr%pc13s%livestemc_xfer(p)	         = pptr%pc13s%livestemc_xfer(p) 	* t1
                pptr%pc13s%deadstemc(p)                  = pptr%pc13s%deadstemc(p)*t1     + deadstemc13_seed*t2
                pptr%pc13s%deadstemc_storage(p)          = pptr%pc13s%deadstemc_storage(p)      * t1
                pptr%pc13s%deadstemc_xfer(p)	         = pptr%pc13s%deadstemc_xfer(p) 	* t1
                pptr%pc13s%livecrootc(p)		 = pptr%pc13s%livecrootc(p) 		* t1
                pptr%pc13s%livecrootc_storage(p)         = pptr%pc13s%livecrootc_storage(p)     * t1
                pptr%pc13s%livecrootc_xfer(p)	         = pptr%pc13s%livecrootc_xfer(p)	* t1
                pptr%pc13s%deadcrootc(p)		 = pptr%pc13s%deadcrootc(p) 		* t1
                pptr%pc13s%deadcrootc_storage(p)         = pptr%pc13s%deadcrootc_storage(p)     * t1
                pptr%pc13s%deadcrootc_xfer(p)	         = pptr%pc13s%deadcrootc_xfer(p)	* t1
                pptr%pc13s%gresp_storage(p) 	         = pptr%pc13s%gresp_storage(p)  	* t1
                pptr%pc13s%gresp_xfer(p)		 = pptr%pc13s%gresp_xfer(p) 		* t1
                pptr%pc13s%cpool(p) 			 = pptr%pc13s%cpool(p)  		* t1
                pptr%pc13s%xsmrpool(p)  		 = pptr%pc13s%xsmrpool(p)		* t1
                pptr%pc13s%pft_ctrunc(p)		 = pptr%pc13s%pft_ctrunc(p) 		* t1
                pptr%pc13s%dispvegc(p)  		 = pptr%pc13s%dispvegc(p)		* t1
                pptr%pc13s%storvegc(p)  		 = pptr%pc13s%storvegc(p)		* t1
                pptr%pc13s%totvegc(p)			 = pptr%pc13s%totvegc(p)		* t1
                pptr%pc13s%totpftc(p)			 = pptr%pc13s%totpftc(p)		* t1
                
             endif
             
             if ( use_c14 ) then
                ! pft-level carbon-14 state variables 
                tot_leaf = pptr%pc14s%leafc(p) + pptr%pc14s%leafc_storage(p) + pptr%pc14s%leafc_xfer(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = pptr%pc14s%leafc(p)/tot_leaf
                   pstor = pptr%pc14s%leafc_storage(p)/tot_leaf
                   pxfer = pptr%pc14s%leafc_xfer(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                pptr%pc14s%leafc(p)         = pptr%pc14s%leafc(p)*t1         + leafc14_seed*pleaf*t2
                pptr%pc14s%leafc_storage(p) = pptr%pc14s%leafc_storage(p)*t1 + leafc14_seed*pstor*t2
                pptr%pc14s%leafc_xfer(p)    = pptr%pc14s%leafc_xfer(p)*t1    + leafc14_seed*pxfer*t2
                pptr%pc14s%frootc(p)			 = pptr%pc14s%frootc(p) 		* t1
                pptr%pc14s%frootc_storage(p)	         = pptr%pc14s%frootc_storage(p) 	* t1
                pptr%pc14s%frootc_xfer(p)		 = pptr%pc14s%frootc_xfer(p)		* t1
                pptr%pc14s%livestemc(p) 		 = pptr%pc14s%livestemc(p)  		* t1
                pptr%pc14s%livestemc_storage(p)          = pptr%pc14s%livestemc_storage(p)      * t1
                pptr%pc14s%livestemc_xfer(p)	         = pptr%pc14s%livestemc_xfer(p) 	* t1
                pptr%pc14s%deadstemc(p)                  = pptr%pc14s%deadstemc(p)*t1     + deadstemc14_seed*t2
                pptr%pc14s%deadstemc_storage(p)          = pptr%pc14s%deadstemc_storage(p)      * t1
                pptr%pc14s%deadstemc_xfer(p)	         = pptr%pc14s%deadstemc_xfer(p) 	* t1
                pptr%pc14s%livecrootc(p)		 = pptr%pc14s%livecrootc(p) 		* t1
                pptr%pc14s%livecrootc_storage(p)         = pptr%pc14s%livecrootc_storage(p)     * t1
                pptr%pc14s%livecrootc_xfer(p)	         = pptr%pc14s%livecrootc_xfer(p)	* t1
                pptr%pc14s%deadcrootc(p)		 = pptr%pc14s%deadcrootc(p) 		* t1
                pptr%pc14s%deadcrootc_storage(p)         = pptr%pc14s%deadcrootc_storage(p)     * t1
                pptr%pc14s%deadcrootc_xfer(p)	         = pptr%pc14s%deadcrootc_xfer(p)	* t1
                pptr%pc14s%gresp_storage(p) 	         = pptr%pc14s%gresp_storage(p)  	* t1
                pptr%pc14s%gresp_xfer(p)		 = pptr%pc14s%gresp_xfer(p) 		* t1
                pptr%pc14s%cpool(p) 			 = pptr%pc14s%cpool(p)  		* t1
                pptr%pc14s%xsmrpool(p)  		 = pptr%pc14s%xsmrpool(p)		* t1
                pptr%pc14s%pft_ctrunc(p)		 = pptr%pc14s%pft_ctrunc(p) 		* t1
                pptr%pc14s%dispvegc(p)  		 = pptr%pc14s%dispvegc(p)		* t1
                pptr%pc14s%storvegc(p)  		 = pptr%pc14s%storvegc(p)		* t1
                pptr%pc14s%totvegc(p)			 = pptr%pc14s%totvegc(p)		* t1
                pptr%pc14s%totpftc(p)			 = pptr%pc14s%totpftc(p)		* t1
             endif
             
             
             tot_leaf = pptr%pns%leafn(p) + pptr%pns%leafn_storage(p) + pptr%pns%leafn_xfer(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                pleaf = pptr%pns%leafn(p)/tot_leaf
                pstor = pptr%pns%leafn_storage(p)/tot_leaf
                pxfer = pptr%pns%leafn_xfer(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(pptr%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             ! pft-level nitrogen state variables
             pptr%pns%leafn(p)         = pptr%pns%leafn(p)*t1         + leafn_seed*pleaf*t2
             pptr%pns%leafn_storage(p) = pptr%pns%leafn_storage(p)*t1 + leafn_seed*pstor*t2
             pptr%pns%leafn_xfer(p)    = pptr%pns%leafn_xfer(p)*t1    + leafn_seed*pxfer*t2
             pptr%pns%frootn(p)  		   = pptr%pns%frootn(p) 		* t1
             pptr%pns%frootn_storage(p)         = pptr%pns%frootn_storage(p) 	* t1
             pptr%pns%frootn_xfer(p) 	   = pptr%pns%frootn_xfer(p)		* t1
             pptr%pns%livestemn(p)		   = pptr%pns%livestemn(p)  		* t1
             pptr%pns%livestemn_storage(p)      = pptr%pns%livestemn_storage(p)      * t1
             pptr%pns%livestemn_xfer(p)         = pptr%pns%livestemn_xfer(p) 	* t1
             pptr%pns%deadstemn(p)              = pptr%pns%deadstemn(p)*t1     + deadstemn_seed*t2
             pptr%pns%deadstemn_storage(p)      = pptr%pns%deadstemn_storage(p)      * t1
             pptr%pns%deadstemn_xfer(p)         = pptr%pns%deadstemn_xfer(p) 	* t1
             pptr%pns%livecrootn(p)  	   = pptr%pns%livecrootn(p) 		* t1
             pptr%pns%livecrootn_storage(p)     = pptr%pns%livecrootn_storage(p)     * t1
             pptr%pns%livecrootn_xfer(p)        = pptr%pns%livecrootn_xfer(p)	* t1
             pptr%pns%deadcrootn(p)  	   = pptr%pns%deadcrootn(p) 		* t1
             pptr%pns%deadcrootn_storage(p)     = pptr%pns%deadcrootn_storage(p)     * t1
             pptr%pns%deadcrootn_xfer(p)        = pptr%pns%deadcrootn_xfer(p)        * t1
             pptr%pns%retransn(p)		   = pptr%pns%retransn(p)		* t1
             pptr%pns%npool(p)		   = pptr%pns%npool(p)  		* t1
             pptr%pns%pft_ntrunc(p)  	   = pptr%pns%pft_ntrunc(p)        	* t1
             pptr%pns%dispvegn(p)		   = pptr%pns%dispvegn(p)		* t1
             pptr%pns%storvegn(p)		   = pptr%pns%storvegn(p)		* t1
             pptr%pns%totvegn(p) 		   = pptr%pns%totvegn(p)		* t1
             pptr%pns%totpftn(p) 		   = pptr%pns%totpftn(p)		* t1
             
             ! update temporary seed source arrays
             ! These are calculated in terms of the required contributions from
             ! column-level seed source
             dwt_leafc_seed(p)   = leafc_seed   * dwt
             if ( use_c13 ) then
                dwt_leafc13_seed(p) = leafc13_seed * dwt
                dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
             endif
             if ( use_c14 ) then
                dwt_leafc14_seed(p) = leafc14_seed * dwt
                dwt_deadstemc14_seed(p) = deadstemc14_seed * dwt
             endif
             dwt_leafn_seed(p)   = leafn_seed   * dwt
             dwt_deadstemc_seed(p)   = deadstemc_seed   * dwt
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
             ptr => pptr%pcs%leafc(p)
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
             ptr => pptr%pcs%leafc_storage(p)
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
             ptr => pptr%pcs%leafc_xfer(p)
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
             ptr => pptr%pcs%frootc(p)
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
             ptr => pptr%pcs%frootc_storage(p)
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
             ptr => pptr%pcs%frootc_xfer(p)
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
             ptr => pptr%pcs%livestemc(p)
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
             ptr => pptr%pcs%livestemc_storage(p)
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
             ptr => pptr%pcs%livestemc_xfer(p)
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
             ptr => pptr%pcs%deadstemc(p)
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
             ptr => pptr%pcs%deadstemc_storage(p)
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
             ptr => pptr%pcs%deadstemc_xfer(p)
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
             ptr => pptr%pcs%livecrootc(p)
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
             ptr => pptr%pcs%livecrootc_storage(p)
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
             ptr => pptr%pcs%livecrootc_xfer(p)
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
             ptr => pptr%pcs%deadcrootc(p)
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
             ptr => pptr%pcs%deadcrootc_storage(p)
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
             ptr => pptr%pcs%deadcrootc_xfer(p)
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
             ptr => pptr%pcs%gresp_storage(p)
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
             ptr => pptr%pcs%gresp_xfer(p)
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
             ptr => pptr%pcs%cpool(p)
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
             ptr => pptr%pcs%xsmrpool(p)
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
             ptr => pptr%pcs%pft_ctrunc(p)
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

             if ( use_c13 ) then
                !-------------------
                ! C13 state update
                !-------------------
                
                ! set pointers to the conversion and product pool fluxes for this pft
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c13flux(p)
                dwt_ptr2 => prod10_c13flux(p)
                dwt_ptr3 => prod100_c13flux(p)
                
                ! leafc 
                ptr => pptr%pc13s%leafc(p)
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
                ptr => pptr%pc13s%leafc_storage(p)
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
                ptr => pptr%pc13s%leafc_xfer(p)
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
                ptr => pptr%pc13s%frootc(p)
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
                ptr => pptr%pc13s%frootc_storage(p)
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
                ptr => pptr%pc13s%frootc_xfer(p)
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
                ptr => pptr%pc13s%livestemc(p)
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
                ptr => pptr%pc13s%livestemc_storage(p)
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
                ptr => pptr%pc13s%livestemc_xfer(p)
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
                ptr => pptr%pc13s%deadstemc(p)
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
                ptr => pptr%pc13s%deadstemc_storage(p)
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
                ptr => pptr%pc13s%deadstemc_xfer(p)
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
                ptr => pptr%pc13s%livecrootc(p)
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
                ptr => pptr%pc13s%livecrootc_storage(p)
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
                ptr => pptr%pc13s%livecrootc_xfer(p)
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
                ptr => pptr%pc13s%deadcrootc(p)
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
                ptr => pptr%pc13s%deadcrootc_storage(p)
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
                ptr => pptr%pc13s%deadcrootc_xfer(p)
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
                ptr => pptr%pc13s%gresp_storage(p)
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
                ptr => pptr%pc13s%gresp_xfer(p)
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
                ptr => pptr%pc13s%cpool(p)
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
                ptr => pptr%pc13s%pft_ctrunc(p)
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

             if ( use_c14 ) then
                !-------------------
                ! C14 state update
                !-------------------
                
                ! set pointers to the conversion and product pool fluxes for this pft
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c14flux(p)
                dwt_ptr2 => prod10_c14flux(p)
                dwt_ptr3 => prod100_c14flux(p)
                
                ! leafc 
                ptr => pptr%pc14s%leafc(p)
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
                ptr => pptr%pc14s%leafc_storage(p)
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
                ptr => pptr%pc14s%leafc_xfer(p)
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
                ptr => pptr%pc14s%frootc(p)
                dwt_ptr0 => dwt_frootc14_to_litter(p)
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
                ptr => pptr%pc14s%frootc_storage(p)
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
                ptr => pptr%pc14s%frootc_xfer(p)
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
                ptr => pptr%pc14s%livestemc(p)
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
                ptr => pptr%pc14s%livestemc_storage(p)
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
                ptr => pptr%pc14s%livestemc_xfer(p)
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
                ptr => pptr%pc14s%deadstemc(p)
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
                ptr => pptr%pc14s%deadstemc_storage(p)
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
                ptr => pptr%pc14s%deadstemc_xfer(p)
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
                ptr => pptr%pc14s%livecrootc(p)
                dwt_ptr0 => dwt_livecrootc14_to_litter(p)
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
                ptr => pptr%pc14s%livecrootc_storage(p)
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
                ptr => pptr%pc14s%livecrootc_xfer(p)
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
                ptr => pptr%pc14s%deadcrootc(p)
                dwt_ptr0 => dwt_deadcrootc14_to_litter(p)
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
                ptr => pptr%pc14s%deadcrootc_storage(p)
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
                ptr => pptr%pc14s%deadcrootc_xfer(p)
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
                ptr => pptr%pc14s%gresp_storage(p)
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
                ptr => pptr%pc14s%gresp_xfer(p)
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
                ptr => pptr%pc14s%cpool(p)
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
                ptr => pptr%pc14s%pft_ctrunc(p)
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
             ptr => pptr%pns%leafn(p)
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
             ptr => pptr%pns%leafn_storage(p)
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
             ptr => pptr%pns%leafn_xfer(p)
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
             ptr => pptr%pns%frootn(p)
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
             ptr => pptr%pns%frootn_storage(p)
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
             ptr => pptr%pns%frootn_xfer(p)
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
             ptr => pptr%pns%livestemn(p)
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
             ptr => pptr%pns%livestemn_storage(p)
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
             ptr => pptr%pns%livestemn_xfer(p)
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
             ptr => pptr%pns%deadstemn(p)
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
             ptr => pptr%pns%deadstemn_storage(p)
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
             ptr => pptr%pns%deadstemn_xfer(p)
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
             ptr => pptr%pns%livecrootn(p)
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
             ptr => pptr%pns%livecrootn_storage(p)
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
             ptr => pptr%pns%livecrootn_xfer(p)
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
             ptr => pptr%pns%deadcrootn(p)
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
             ptr => pptr%pns%deadcrootn_storage(p)
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
             ptr => pptr%pns%deadcrootn_xfer(p)
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
             ptr => pptr%pns%retransn(p)
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
             ptr => pptr%pns%npool(p)
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
             ptr => pptr%pns%pft_ntrunc(p)
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
             cptr%ccf%dwt_seedc_to_leaf(c) = cptr%ccf%dwt_seedc_to_leaf(c) + dwt_leafc_seed(p)/dt
             cptr%ccf%dwt_seedc_to_deadstem(c) = cptr%ccf%dwt_seedc_to_deadstem(c) &
                  + dwt_deadstemc_seed(p)/dt
             
             if ( use_c13 ) then
                cptr%cc13f%dwt_seedc_to_leaf(c) = cptr%cc13f%dwt_seedc_to_leaf(c) + dwt_leafc13_seed(p)/dt
                cptr%cc13f%dwt_seedc_to_deadstem(c) = cptr%cc13f%dwt_seedc_to_deadstem(c) &
                     + dwt_deadstemc13_seed(p)/dt
             endif
             
             if ( use_c14 ) then	
                cptr%cc14f%dwt_seedc_to_leaf(c) = cptr%cc14f%dwt_seedc_to_leaf(c) + dwt_leafc14_seed(p)/dt
                cptr%cc14f%dwt_seedc_to_deadstem(c) = cptr%cc14f%dwt_seedc_to_deadstem(c) &
                     + dwt_deadstemc14_seed(p)/dt
             endif
             
             ! N fluxes
             cptr%cnf%dwt_seedn_to_leaf(c) = cptr%cnf%dwt_seedn_to_leaf(c) + dwt_leafn_seed(p)/dt
             cptr%cnf%dwt_seedn_to_deadstem(c) = cptr%cnf%dwt_seedn_to_deadstem(c) &
                  + dwt_deadstemn_seed(p)/dt
          end if
       end do
    end do
    
    
    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do pi = 1,max_pft_per_col
          do c = begc, endc
             if ( pi <=  cptr%npfts(c) ) then
                p = cptr%pfti(c) + pi - 1
                
                ! fine root litter carbon fluxes
                cptr%ccf%dwt_frootc_to_litr_met_c(c,j) = cptr%ccf%dwt_frootc_to_litr_met_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                cptr%ccf%dwt_frootc_to_litr_cel_c(c,j) = cptr%ccf%dwt_frootc_to_litr_cel_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                cptr%ccf%dwt_frootc_to_litr_lig_c(c,j) = cptr%ccf%dwt_frootc_to_litr_lig_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                
                
                ! fine root litter nitrogen fluxes
                cptr%cnf%dwt_frootn_to_litr_met_n(c,j) = cptr%cnf%dwt_frootn_to_litr_met_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                cptr%cnf%dwt_frootn_to_litr_cel_n(c,j) = cptr%cnf%dwt_frootn_to_litr_cel_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                cptr%cnf%dwt_frootn_to_litr_lig_n(c,j) = cptr%cnf%dwt_frootn_to_litr_lig_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                
                ! livecroot fluxes to cwd
                cptr%ccf%dwt_livecrootc_to_cwdc(c,j) = cptr%ccf%dwt_livecrootc_to_cwdc(c,j) + &
                     (dwt_livecrootc_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                cptr%cnf%dwt_livecrootn_to_cwdn(c,j) = cptr%cnf%dwt_livecrootn_to_cwdn(c,j) + &
                     (dwt_livecrootn_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                
                ! deadcroot fluxes to cwd
                cptr%ccf%dwt_deadcrootc_to_cwdc(c,j) = cptr%ccf%dwt_deadcrootc_to_cwdc(c,j) + &
                     (dwt_deadcrootc_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                cptr%cnf%dwt_deadcrootn_to_cwdn(c,j) = cptr%cnf%dwt_deadcrootn_to_cwdn(c,j) + &
                     (dwt_deadcrootn_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
             
                if ( use_c13 ) then
                   ! C13 fine root litter fluxes
                   cptr%cc13f%dwt_frootc_to_litr_met_c(c,j) = cptr%cc13f%dwt_frootc_to_litr_met_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                   cptr%cc13f%dwt_frootc_to_litr_cel_c(c,j) = cptr%cc13f%dwt_frootc_to_litr_cel_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                   cptr%cc13f%dwt_frootc_to_litr_lig_c(c,j) = cptr%cc13f%dwt_frootc_to_litr_lig_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                   ! livecroot fluxes to cwd
                   cptr%cc13f%dwt_livecrootc_to_cwdc(c,j) = cptr%cc13f%dwt_livecrootc_to_cwdc(c,j) + &
                        (dwt_livecrootc13_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                   ! deadcroot fluxes to cwd
                   cptr%cc13f%dwt_deadcrootc_to_cwdc(c,j) = cptr%cc13f%dwt_deadcrootc_to_cwdc(c,j) + &
                        (dwt_deadcrootc13_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                   
                endif
                
                if ( use_c14 ) then                   
                   ! C14 fine root litter fluxes
                   cptr%cc14f%dwt_frootc_to_litr_met_c(c,j) = cptr%cc14f%dwt_frootc_to_litr_met_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flab(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                   cptr%cc14f%dwt_frootc_to_litr_cel_c(c,j) = cptr%cc14f%dwt_frootc_to_litr_cel_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_fcel(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                   cptr%cc14f%dwt_frootc_to_litr_lig_c(c,j) = cptr%cc14f%dwt_frootc_to_litr_lig_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flig(pptr%itype(p)))/dt * pptr%pps%froot_prof(p,j)
                   ! livecroot fluxes to cwd
                   cptr%cc14f%dwt_livecrootc_to_cwdc(c,j) = cptr%cc14f%dwt_livecrootc_to_cwdc(c,j) + &
                        (dwt_livecrootc14_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                   ! deadcroot fluxes to cwd
                   cptr%cc14f%dwt_deadcrootc_to_cwdc(c,j) = cptr%cc14f%dwt_deadcrootc_to_cwdc(c,j) + &
                        (dwt_deadcrootc14_to_litter(p))/dt * pptr%pps%croot_prof(p,j)
                endif
                
             end if
          end do
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do pi = 1,max_pft_per_col
       do c = begc,endc
          if (pi <= cptr%npfts(c)) then
             p = cptr%pfti(c) + pi - 1
             
             ! column-level fluxes are accumulated as positive fluxes.
             ! column-level C flux updates
             cptr%ccf%dwt_conv_cflux(c) = cptr%ccf%dwt_conv_cflux(c) - conv_cflux(p)/dt
             cptr%ccf%dwt_prod10c_gain(c) = cptr%ccf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
             cptr%ccf%dwt_prod100c_gain(c) = cptr%ccf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt

             ! These magic constants should be replaced with: nbrdlf_evr_trp_tree and nbrdlf_dcd_trp_tree
             if(ivt(p)==4.or.ivt(p)==6)then
                cptr%ccf%lf_conv_cflux(c) = cptr%ccf%lf_conv_cflux(c) - conv_cflux(p)/dt
             end if
             
             if ( use_c13 ) then
                ! C13 column-level flux updates
                cptr%cc13f%dwt_conv_cflux(c) = cptr%cc13f%dwt_conv_cflux(c) - conv_c13flux(p)/dt
                cptr%cc13f%dwt_prod10c_gain(c) = cptr%cc13f%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
                cptr%cc13f%dwt_prod100c_gain(c) = cptr%cc13f%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt
             endif
             
             if ( use_c14 ) then
                ! C14 column-level flux updates
                cptr%cc14f%dwt_conv_cflux(c) = cptr%cc14f%dwt_conv_cflux(c) - conv_c14flux(p)/dt
                cptr%cc14f%dwt_prod10c_gain(c) = cptr%cc14f%dwt_prod10c_gain(c) - prod10_c14flux(p)/dt
                cptr%cc14f%dwt_prod100c_gain(c) = cptr%cc14f%dwt_prod100c_gain(c) - prod100_c14flux(p)/dt
             endif
             
             ! column-level N flux updates
             cptr%cnf%dwt_conv_nflux(c) = cptr%cnf%dwt_conv_nflux(c) - conv_nflux(p)/dt
             cptr%cnf%dwt_prod10n_gain(c) = cptr%cnf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
             cptr%cnf%dwt_prod100n_gain(c) = cptr%cnf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt
             
          end if
       end do
    end do
    
    ! Deallocate pft-level flux arrays
    deallocate(dwt_leafc_seed)
    deallocate(dwt_leafn_seed)
    deallocate(dwt_deadstemc_seed)
    deallocate(dwt_deadstemn_seed)
    deallocate(dwt_frootc_to_litter)
    deallocate(dwt_livecrootc_to_litter)
    deallocate(dwt_deadcrootc_to_litter)
    deallocate(dwt_frootn_to_litter)
    deallocate(dwt_livecrootn_to_litter)
    deallocate(dwt_deadcrootn_to_litter)
    deallocate(conv_cflux)
    deallocate(prod10_cflux)
    deallocate(prod100_cflux)
    deallocate(conv_nflux)
    deallocate(prod10_nflux)
    deallocate(prod100_nflux)
             
    if ( use_c13 ) then
       deallocate(dwt_leafc13_seed)
       deallocate(dwt_deadstemc13_seed)
       deallocate(dwt_frootc13_to_litter)
       deallocate(dwt_livecrootc13_to_litter)
       deallocate(dwt_deadcrootc13_to_litter)
       deallocate(conv_c13flux)
       deallocate(prod10_c13flux)
       deallocate(prod100_c13flux)
    endif
             
    if ( use_c14 ) then
       deallocate(dwt_leafc14_seed)
       deallocate(dwt_deadstemc14_seed)
       deallocate(dwt_frootc14_to_litter)
       deallocate(dwt_livecrootc14_to_litter)
       deallocate(dwt_deadcrootc14_to_litter)
       deallocate(conv_c14flux)
       deallocate(prod10_c14flux)
       deallocate(prod100_c14flux)
    endif
    
  end subroutine pftdyn_cnbal
#endif

#if (defined CNDV)
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

    pptr => clm3%g%l%c%p

    call get_proc_bounds(begp=begp,endp=endp)

    allocate(wtcol_old(begp:endp),stat=ier)
    if (ier /= 0) then
       call endrun( subname//'::ERROR: pftwt_init allocation error for wtcol_old')
    end if

    if (nsrest == nsrStartup) then
       do p = begp,endp
          pptr%pdgvs%fpcgrid(p) = pptr%wtcol(p)
          pptr%pdgvs%fpcgridold(p) = pptr%wtcol(p)
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

    lptr => clm3%g%l
    pptr => clm3%g%l%c%p

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
          pptr%wtcol(p)   = pptr%pdgvs%fpcgrid(p) + &
                     wt1 * (pptr%pdgvs%fpcgridold(p) - pptr%pdgvs%fpcgrid(p))
          pptr%wtlunit(p) = pptr%wtcol(p)
          pptr%wtgcell(p) = pptr%wtcol(p) * lptr%wtgcell(l)

          if (mon==1 .and. day==1 .and. sec==dtime .and. nstep>0) then
             pptr%pdgvs%fpcgridold(p) = pptr%pdgvs%fpcgrid(p)
          end if
       end if
    end do

  end subroutine pftwt_interp
#endif

#ifdef CN
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
   pgridcell                      => clm3%g%l%c%p%gridcell
   
   ivt                            => clm3%g%l%c%p%itype
   leafc                          => clm3%g%l%c%p%pcs%leafc
   frootc                         => clm3%g%l%c%p%pcs%frootc
   livestemc                      => clm3%g%l%c%p%pcs%livestemc
   deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
   livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
   deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
   xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
   leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
   frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
   livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
   deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
   livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
   deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
   gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
   leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
   frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
   livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
   deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
   livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
   deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
   gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
   leafn                          => clm3%g%l%c%p%pns%leafn
   frootn                         => clm3%g%l%c%p%pns%frootn
   livestemn                      => clm3%g%l%c%p%pns%livestemn
   deadstemn                      => clm3%g%l%c%p%pns%deadstemn
   livecrootn                     => clm3%g%l%c%p%pns%livecrootn
   deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
   retransn                       => clm3%g%l%c%p%pns%retransn
   leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
   frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
   livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
   deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
   livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
   deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
   leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
   frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
   livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
   deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
   livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
   deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
   hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
   hrv_deadstemc_to_prod10c         => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
   hrv_deadstemc_to_prod100c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
   hrv_xsmrpool_to_atm              => clm3%g%l%c%p%pcf%hrv_xsmrpool_to_atm
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
   hrv_deadstemn_to_prod10n         => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
   hrv_deadstemn_to_prod100n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter


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
  use clm_varpar, only : max_pft_per_col, maxpatch_pft, nlevdecomp
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
   logical , pointer :: pactive(:)  ! true=>do computations on this pft (see reweightMod for details)
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
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
   real(r8), pointer :: harvest_c_to_litr_met_c(:,:)               ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_litr_cel_c(:,:)               ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_litr_lig_c(:,:)               ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_cwdc(:,:)                     ! C fluxes associated with harvest to CWD pool (gC/m3/s)
   real(r8), pointer :: harvest_n_to_litr_met_n(:,:)               ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_litr_cel_n(:,:)               ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_litr_lig_n(:,:)               ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_cwdn(:,:)                     ! N fluxes associated with harvest to CWD pool (gN/m3/s)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: chrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: chrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
   real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
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
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   chrv_deadstemc_to_prod10c        => clm3%g%l%c%ccf%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => clm3%g%l%c%ccf%hrv_deadstemc_to_prod100c
   chrv_deadstemn_to_prod10n        => clm3%g%l%c%cnf%hrv_deadstemn_to_prod10n
   chrv_deadstemn_to_prod100n       => clm3%g%l%c%cnf%hrv_deadstemn_to_prod100n
   harvest_c_to_litr_met_c          => clm3%g%l%c%ccf%harvest_c_to_litr_met_c
   harvest_c_to_litr_cel_c          => clm3%g%l%c%ccf%harvest_c_to_litr_cel_c
   harvest_c_to_litr_lig_c          => clm3%g%l%c%ccf%harvest_c_to_litr_lig_c
   harvest_c_to_cwdc                => clm3%g%l%c%ccf%harvest_c_to_cwdc
   harvest_n_to_litr_met_n          => clm3%g%l%c%cnf%harvest_n_to_litr_met_n
   harvest_n_to_litr_cel_n          => clm3%g%l%c%cnf%harvest_n_to_litr_cel_n
   harvest_n_to_litr_lig_n          => clm3%g%l%c%cnf%harvest_n_to_litr_lig_n
   harvest_n_to_cwdn                => clm3%g%l%c%cnf%harvest_n_to_cwdn


   ! assign local pointers to pft-level arrays
   pactive                        => clm3%g%l%c%p%active
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
   phrv_deadstemn_to_prod10n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
   phrv_deadstemn_to_prod100n       => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter
   leaf_prof                      => clm3%g%l%c%p%pps%leaf_prof
   froot_prof                     => clm3%g%l%c%p%pps%froot_prof
   croot_prof                     => clm3%g%l%c%p%pps%croot_prof
   stem_prof                      => clm3%g%l%c%p%pps%stem_prof

   do j = 1, nlevdecomp
      do pi = 1,maxpatch_pft
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (pi <=  npfts(c)) then
               p = pfti(c) + pi - 1
               
               if (pactive(p)) then
                  
                  
                  ! leaf harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood harvest mortality carbon fluxes
                  harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                       hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j) 
                  harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                       hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                       hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 
                  
                  ! storage harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                       hrv_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
                  ! transfer harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                       hrv_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
                  ! leaf harvest mortality nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                       hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                       hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root litter nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_cel_n(c,j) = harvest_n_to_litr_cel_n(c,j) + &
                       hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_lig_n(c,j) = harvest_n_to_litr_lig_n(c,j) + &
                       hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood harvest mortality nitrogen fluxes
                  harvest_n_to_cwdn(c,j)  = harvest_n_to_cwdn(c,j)  + &
                       hrv_livestemn_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                       hrv_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                       hrv_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
                  ! retranslocated N pool harvest mortality fluxes
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                  
                  ! storage harvest mortality nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                       hrv_leafn_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                       hrv_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_livestemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_deadstemn_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_livecrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_deadcrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
                  ! transfer harvest mortality nitrogen fluxes
                  harvest_n_to_litr_met_n(c,j)      = harvest_n_to_litr_met_n(c,j)      + &
                       hrv_leafn_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)     = harvest_n_to_litr_met_n(c,j)     + &
                       hrv_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_livestemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j)  = harvest_n_to_litr_met_n(c,j)  + &
                       hrv_deadstemn_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_livecrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_n_to_litr_met_n(c,j) = harvest_n_to_litr_met_n(c,j) + &
                       hrv_deadcrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
               end if
            end if
            
         end do
         
      end do
   end do
   
   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1
            
            if (pactive(p)) then
               
               
               ! wood harvest mortality carbon fluxes to product pools
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                    phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                    phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               
               
               ! wood harvest mortality nitrogen fluxes to product pools
               chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
                    phrv_deadstemn_to_prod10n(p)  * wtcol(p)
               chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
                    phrv_deadstemn_to_prod100n(p)  * wtcol(p)
            end if
         end if
         
      end do
      
   end do
   
 end subroutine CNHarvestPftToColumn
!-----------------------------------------------------------------------
#endif

end module pftdynMod
