module mo_fstrat
  !---------------------------------------------------------------
  !	... variables for the upper boundary values
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, pver,pverp, begchunk, endchunk
  use chem_mods,    only : gas_pcnst
  use abortutils,   only : endrun
  use cam_pio_utils, only : cam_pio_openfile
  use pio
  use cam_logfile,  only : iulog

  implicit none

  private
  public  :: fstrat_inti
  public  :: set_fstrat_vals
  public  :: set_fstrat_h2o
  public  :: has_fstrat

  save

  real(r8), parameter :: taurelax = 864000._r8        ! 10 days
  integer  :: gndx = 0
  integer  :: table_nox_ndx = -1
  integer  :: table_h2o_ndx = -1
  integer  :: table_ox_ndx  = -1
  integer  :: table_synoz_ndx  = -1
  integer  :: no_ndx
  integer  :: no2_ndx
  integer  :: h2o_ndx
  integer  :: ox_ndx
  integer  :: o3s_ndx
  integer  :: o3inert_ndx
  integer  :: o3a_ndx
  integer  :: xno_ndx
  integer  :: xno2_ndx
  integer  :: synoz_ndx
  integer  :: o3rad_ndx
  real(r8) :: facrelax
  real(r8) :: days(12)
  real(r8), allocatable   :: ub_plevs(:)         ! table midpoint pressure levels (Pa)
  real(r8), allocatable   :: ub_plevse(:)        ! table edge pressure levels (Pa)
  integer                 :: ub_nlevs            ! # of levs in ubc file
  integer                 :: ub_nlat             ! # of lats in ubc file
  integer                 :: ub_nspecies         ! # of species in ubc file
  integer                 :: ub_nmonth           ! # of months in ubc file
  real(r8), allocatable   :: mr_ub(:,:,:,:,:)      ! vmr
  integer,  allocatable   :: map(:)              ! species indices for ubc species
  logical :: sim_has_nox
  integer :: dtime               ! model time step (s)
  logical :: has_fstrat(gas_pcnst)

contains

  subroutine fstrat_inti( fstrat_file, fstrat_list )
    !------------------------------------------------------------------
    !	... initialize upper boundary values
    !------------------------------------------------------------------

    use mo_constants,  only : d2r
    use phys_grid,     only : get_ncols_p, get_rlat_all_p

    use time_manager,  only : get_step_size
    use time_manager,  only : get_calday
    use ioFileMod,     only : getfil
    use spmd_utils,    only : masterproc
#ifdef SPMD
    use mpishorthand,  only : mpicom,mpiint,mpir8
#endif
    use mo_tracname,  only : solsym
    use chem_mods,    only : gas_pcnst
    use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
    use constituents, only : pcnst
    use dyn_grid,     only : get_dyn_grid_parm
    use interpolate_data, only : interp_type, lininterp_init, lininterp
    implicit none

    character(len=*), intent(in) :: fstrat_file
    character(len=*), intent(in) :: fstrat_list(:)

    !------------------------------------------------------------------
    !	... local variables
    !------------------------------------------------------------------
    real(r8), parameter :: mb2pa = 100._r8

    integer :: i, j, nchar
    integer :: spcno, lev, month, ierr
    integer :: vid, ndims
    type(file_desc_t) :: ncid
    integer :: dimid_lat, dimid_lev, dimid_species, dimid_month
    integer :: dimid(4)
    integer :: start(4)
    integer :: count(4)
    integer, parameter :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)
    integer :: ncols, c
    real(r8), allocatable :: mr_ub_in(:,:,:,:)
    real(r8), allocatable :: lat(:)
    real(r8) :: to_lats(pcols)
    character(len=80) :: attribute
    character(len=8)  :: wrk_name
    character(len=25), allocatable :: ub_species_names(:)
    character(len=256) :: locfn
    type(interp_type) :: lat_wgts


    !-----------------------------------------------------------------------
    !       ... get species indicies
    !-----------------------------------------------------------------------    
    no_ndx      = get_spc_ndx( 'NO' )
    no2_ndx     = get_spc_ndx( 'NO2' )
    sim_has_nox = no_ndx > 0 .or. no2_ndx > 0
    ox_ndx      = get_spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
    end if
    o3s_ndx     = get_spc_ndx( 'O3S' )
    o3inert_ndx = get_spc_ndx( 'O3INERT' )
    o3rad_ndx   = get_spc_ndx( 'O3RAD' )
    synoz_ndx   = get_spc_ndx( 'SYNOZ' )
    o3a_ndx     = get_spc_ndx( 'O3A' )
    xno_ndx      = get_spc_ndx( 'XNO' )
    xno2_ndx      = get_spc_ndx( 'XNO2' )

    if (.not.((o3rad_ndx > 0) .eqv. (synoz_ndx > 0))) then
       call endrun('fstrat_inti: Both SYNOZ and O3RAD are required in chemical mechanism.')
    endif
    if(masterproc) then
       if (synoz_ndx > 0) then
          if ( .not. any(fstrat_list == 'O3RAD') ) then
             write(iulog,*) 'fstrat_inti: ***WARNING*** O3RAD is not include in fstrat_list namelist variable'
          endif
       else if (ox_ndx > 0) then
          if ( .not. any(fstrat_list == 'O3') ) then
             write(iulog,*) 'fstrat_inti: ***WARNING*** O3 is not include in fstrat_list namelist variable'
          endif
       endif
    end if
    h2o_ndx     = get_spc_ndx( 'H2O' )
    if( h2o_ndx < 0 ) then
       h2o_ndx  = get_inv_ndx( 'H2O' )
    end if

    has_fstrat(:) = .false.

    do i = 1,pcnst

       if ( len_trim(fstrat_list(i))==0 ) exit

       j = get_spc_ndx(fstrat_list(i))

       if ( j > 0 ) then
          has_fstrat(j) = .true.
       else
          write(iulog,*) 'fstrat_inti: '//trim(fstrat_list(i))//' is not included in species set'
          call endrun('fstrat_inti: invalid fixed stratosphere species')
       endif

    enddo

    if (.not. any(has_fstrat(:))) return

    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (fstrat_file, locfn, 0)
    call cam_pio_openfile (ncid, trim(locfn), PIO_NOWRITE)
    !-----------------------------------------------------------------------
    !       ... get latitude
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'lat', dimid_lat )
    ierr = pio_inq_dimlen( ncid, dimid_lat, ub_nlat )
    allocate( lat(ub_nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, lat )
    lat(:ub_nlat) = lat(:ub_nlat) * d2r

    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: failed to deallocate lat; ierr = ',ierr
       call endrun
    end if

    !-----------------------------------------------------------------------
    !       ... get vertical coordinate (if necessary, convert units to pa)
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'lev', dimid_lev )
    ierr = pio_inq_dimlen( ncid, dimid_lev, ub_nlevs )
    allocate( ub_plevs(ub_nlevs), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: ub_plevs allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( ncid, 'lev', vid )
    ierr = pio_get_var( ncid, vid, ub_plevs )
    attribute(:) = ' '
    call pio_seterrorhandling(ncid, pio_BCAST_error)
    ierr = pio_get_att( ncid, vid, 'units', attribute )
    call pio_seterrorhandling(ncid, pio_INTERNAL_error)
    if( ierr == PIO_noerr )then
       if( trim(attribute) == 'mb' .or. trim(attribute) == 'hpa' )then
          if (masterproc) write(iulog,*) 'fstrat_inti: units for lev = ',trim(attribute),'... converting to pa'
          ub_plevs(:) = mb2pa * ub_plevs(:)
       else if( trim(attribute) /= 'pa' .and. trim(attribute) /= 'pa' )then
          write(iulog,*) 'fstrat_inti: unknown units for lev, units=*',trim(attribute),'*'
          write(iulog,*) 'fstrat_inti: ',attribute=='mb',trim(attribute)=='mb',attribute(1:2)=='mb'
          call endrun
       end if
    else
       write(iulog,*) 'fstrat_inti: warning! units attribute for lev missing, assuming mb'
       ub_plevs(:) = mb2pa * ub_plevs(:)
    end if
    !-----------------------------------------------------------------------
    !       ... get time and species dimensions
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'month', dimid_month )
    ierr = pio_inq_dimlen( ncid, dimid_month, ub_nmonth )
    if( ub_nmonth /= 12 )then
       write(iulog,*) 'fstrat_inti: error! number of months = ',ub_nmonth,', expecting 12'
       call endrun
    end if
    ierr = pio_inq_dimid( ncid, 'species', dimid_species )
    ierr = pio_inq_dimlen( ncid, dimid_species, ub_nspecies )

    !------------------------------------------------------------------
    !	... allocate arrays
    !------------------------------------------------------------------
    allocate( mr_ub_in(ub_nlat,ub_nspecies,ub_nmonth,ub_nlevs), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: mr_ub_in allocation error = ',ierr
       call endrun
    end if
    allocate( mr_ub(pcols,ub_nspecies,ub_nmonth,ub_nlevs,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: mr_ub allocation error = ',ierr
       call endrun
    end if

    !------------------------------------------------------------------
    !	... read in the species names
    !------------------------------------------------------------------

    ierr = pio_inq_varid( ncid, 'specname', vid )
    ierr = pio_inq_vardimid( ncid, vid, dimid )
    ierr = pio_inq_dimlen( ncid, dimid(1), nchar )
    allocate( ub_species_names(ub_nspecies), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: ub_species_names allocation error = ',ierr
       call endrun
    end if
    allocate( map(ub_nspecies), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'fstrat_inti: map allocation error = ',ierr
       call endrun
    end if

    table_loop :  do i = 1,ub_nspecies
       start(:2) = (/ 1, i /)
       count(:2) = (/ nchar, 1 /)
       ub_species_names(i)(:) = ' '
       ierr = pio_get_var( ncid, vid, start(:2), count(:2), ub_species_names(i:i))
       if( trim(ub_species_names(i)) == 'NOX' ) then
          table_nox_ndx = i
       else if( trim(ub_species_names(i)) == 'H2O' ) then
          table_h2o_ndx = i
       else if( trim(ub_species_names(i)) == 'OX' ) then
          table_ox_ndx = i
       else if( trim(ub_species_names(i)) == 'SYNOZ' ) then
          table_synoz_ndx = i
       end if
       map(i) = 0
       do j = 1,gas_pcnst 
          if( trim(ub_species_names(i)) == trim(solsym(j)) ) then
             if( has_fstrat(j) ) then 
                map(i) = j
                if( masterproc ) write(iulog,*) 'fstrat_inti: '//trim(solsym(j))//' is fixed in stratosphere'
                exit
             end if
          endif
       end do
       if( map(i) == 0 ) then
          if( trim(ub_species_names(i)) == 'OX' ) then
             if( o3rad_ndx > 0 ) then
                wrk_name = 'O3RAD'
             else
                wrk_name = 'O3'
             end if
             do j = 1,gas_pcnst
                if( trim(wrk_name) == trim(solsym(j)) ) then
                   if( has_fstrat(j) ) then 
                      if( masterproc ) write(iulog,*) 'fstrat_inti: '//trim(solsym(j))//' is fixed in stratosphere'
                      map(i) = j
                      exit
                   end if
                end if
             end do
             if( map(i) == 0 ) then
                write(iulog,*) 'fstrat_inti: ubc table species ',trim(ub_species_names(i)), ' not used'
             end if
          else if( (trim(ub_species_names(i)) /= 'NOX') .and. (trim(ub_species_names(i)) /= 'H2O') ) then
             write(iulog,*) 'fstrat_inti: ubc table species ',trim(ub_species_names(i)), ' not used'
          end if
       end if
    end do table_loop

    if( table_nox_ndx > 0 ) then
       if ( any(fstrat_list(:) == 'NO') .or. any(fstrat_list(:) == 'NO2') ) then
          map(table_nox_ndx) = gas_pcnst + 1
       else
          write(iulog,*) 'fstrat_inti: ubc table species ',trim(ub_species_names(table_nox_ndx)), ' not used'
       endif
    end if
    if( table_h2o_ndx > 0 ) then
       if ( h2o_ndx > 0 ) then
          map(table_h2o_ndx) = gas_pcnst + 2
       else
          write(iulog,*) 'fstrat_inti: ubc table species ',trim(ub_species_names(table_h2o_ndx)), ' not used'
       endif
    end if

    if (masterproc) write(iulog,*) 'fstrat_inti: h2o_ndx, table_h2o_ndx  = ', h2o_ndx, table_h2o_ndx

    !------------------------------------------------------------------
    !	... check dimensions for vmr variable
    !------------------------------------------------------------------
    ierr = pio_inq_varid( ncid, 'vmr', vid )
    ierr = pio_inq_varndims( ncid, vid, ndims )
    if( ndims /= 4 ) then
       write(iulog,*) 'fstrat_inti: error! variable vmr has ndims = ',ndims,', expecting 4'
       call endrun
    end if
    ierr = pio_inq_vardimid( ncid, vid, dimid )
    if( dimid(1) /= dimid_lat .or. dimid(2) /= dimid_species .or. &
         dimid(3) /= dimid_month .or. dimid(4) /= dimid_lev )then
       write(iulog,*) 'fstrat_inti: error! dimensions in wrong order for variable vmr,'// &
            'expecting (lat,species,month,lev)'
       call endrun
    end if

    !------------------------------------------------------------------
    !	... read in the ub mixing ratio values
    !------------------------------------------------------------------
    start = (/ 1, 1, 1, 1 /)
    count = (/ ub_nlat, ub_nspecies, ub_nmonth, ub_nlevs /)

    ierr = pio_get_var( ncid, vid, start, count, mr_ub_in )
    call pio_closefile (ncid)
    !--------------------------------------------------------------------
    !	... regrid
    !--------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call lininterp_init(lat,ub_nlat, to_lats, ncols, 1, lat_wgts)

       do lev = 1,ub_nlevs
          do month = 1,ub_nmonth
             do spcno = 1,ub_nspecies
                if( map(spcno) > 0 ) then

                   call lininterp(mr_ub_in(:,spcno, month, lev), ub_nlat, mr_ub(:,spcno,month,lev,c), &
                        ncols, lat_wgts)

                   !                call regrid_1d( mr_ub_in(:,spcno,month,lev), mr_ub(:,spcno,month,lev), gndx, &
                   !                     do_lat=.true.,to_lat_min=1, to_lat_max=plat )
#ifdef DEBUG
                   if( lev == 25 .and. month == 1 .and. spcno == 1 ) then
                      write(iulog,*) 'mr_ub_in='
                      write(iulog,'(10f7.1)') mr_ub_in(:,spcno,month,lev)*1.e9_r8
                      write(iulog,*) 'mr_ub='
                      write(iulog,'(10f7.1)') mr_ub(:,spcno,month,lev)*1.e9_r8
                   end if
#endif
                   !                mr_ub(1,spcno,month,lev) = mr_ub(2,spcno,month,lev)
                   !                mr_ub(plat,spcno,month,lev) = mr_ub(plat-1,spcno,month,lev)

                end if

             end do
          end do
       end do
    end do
    !--------------------------------------------------------
    !	... initialize the monthly day of year times
    !--------------------------------------------------------
    do month = 1,12
       days(month) = get_calday( dates(month), 0 )
    end do

    !--------------------------------------------------------
    !   	... set up the relaxation for lower stratosphere
    !--------------------------------------------------------
    ! 	... taurelax = relaxation timescale (in sec)
    !           facrelax = fractional relaxation towards ubc
    !            1 => use ubc
    !            0 => ignore ubc, use model concentrations
    !--------------------------------------------------------
    dtime = get_step_size()
    facrelax = 1._r8 - exp( -real(dtime)/taurelax )

    !--------------------------------------------------------
    ! 	... setup conserving interp for OX
    !--------------------------------------------------------
    if( table_ox_ndx > 0 ) then
       allocate( ub_plevse(ub_nlevs-1), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'fstrat_inti: ub_plevse allocation error = ',ierr
          call endrun
       end if
       ub_plevse(1:ub_nlevs-1) = .5_r8*(ub_plevs(1:ub_nlevs-1) + ub_plevs(2:ub_nlevs))
    end if

  end subroutine fstrat_inti

  subroutine set_fstrat_vals( vmr, pmid, pint, ltrop, calday, ncol,lchnk )

    !--------------------------------------------------------------------
    !	... set the upper boundary values for :
    !           ox, nox, hno3, ch4, co, n2o, n2o5 & stratospheric o3
    !--------------------------------------------------------------------

    use mo_synoz, only : synoz_region => po3

    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    integer,  intent(in)    :: lchnk             ! chunk number
    integer,  intent(in)    :: ncol              ! columns in chunk
    integer,  intent(in)    :: ltrop(pcols)      ! tropopause vertical index
    real(r8), intent(in)    :: calday            ! day of year including fraction
    real(r8), intent(in)    :: pmid(pcols,pver)  ! midpoint pressure (Pa)
    real(r8), intent(in)    :: pint(pcols,pverp) ! interface pressure (Pa)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst) ! species concentrations (mol/mol)

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer, parameter :: zlower = pver
    real(r8), parameter    :: synoz_thres = 100.e-9_r8      ! synoz threshold
    real(r8), parameter    :: o3rad_relax = 0.5_r8*86400._r8 ! 1/2 day relaxation constant
    real(r8), parameter    :: synoz_relax = 2._r8*86400._r8 ! 2 day relaxation constant
    real(r8), parameter    :: synoz_strat_relax = 5._r8*86400._r8 ! 5 day relaxation constant

    integer  ::  m, last, next, i, k, k1, km
    integer  ::  astat
    integer  ::  kmax(ncol)
    integer  ::  levrelax
    integer  ::  kl(ncol,zlower)
    integer  ::  ku(ncol,zlower)
    real(r8) ::  vmrrelax
    real(r8) ::  fac_relax
    real(r8) ::  pinterp
    real(r8) ::  nox_ubc, xno, xno2, rno
    real(r8) ::  dels
    real(r8) ::  delp(ncol,zlower)
    real(r8) ::  pint_vals(2)
    real(r8), allocatable :: table_ox(:)
    logical  ::  found_trop
    integer  :: lat

    if (.not. any(has_fstrat(:))) return

    !--------------------------------------------------------
    !	... setup the time interpolation
    !--------------------------------------------------------
    if( calday < days(1) ) then
       next = 1
       last = 12
       dels = (365._r8 + calday - days(12)) / (365._r8 + days(1) - days(12))
    else if( calday >= days(12) ) then
       next = 1
       last = 12
       dels = (calday - days(12)) / (365._r8 + days(1) - days(12))
    else
       do m = 11,1,-1
          if( calday >= days(m) ) then
             exit
          end if
       end do
       last = m
       next = m + 1
       dels = (calday - days(m)) / (days(m+1) - days(m))
    end if
    dels = max( min( 1._r8,dels ),0._r8 )

    !--------------------------------------------------------
    !	... setup the pressure interpolation
    !--------------------------------------------------------
    do k = 1,zlower
       do i = 1,ncol
          if( pmid(i,k) <= ub_plevs(1) ) then
             kl(i,k) = 1
             ku(i,k) = 1
             delp(i,k) = 0._r8
          else if( pmid(i,k) >= ub_plevs(ub_nlevs) ) then
             kl(i,k) = ub_nlevs
             ku(i,k) = ub_nlevs
             delp(i,k) = 0._r8
          else
             pinterp = pmid(i,k)
             do k1 = 2,ub_nlevs
                if( pinterp <= ub_plevs(k1) ) then
                   ku(i,k) = k1
                   kl(i,k) = k1 - 1
                   delp(i,k) = log( pinterp/ub_plevs(kl(i,k)) ) &
                        / log( ub_plevs(ku(i,k))/ub_plevs(kl(i,k)) )
                   exit
                end if
             end do
          end if
       end do
    end do

    !--------------------------------------------------------
    !	... find max level less than 50 mb
    !           fix UB vals from top of model to this level
    !--------------------------------------------------------
    do i = 1,ncol
       do k = 2,pver
          if( pmid(i,k) > 50.e2_r8 ) then
             kmax(i) = k
             exit
          end if
       end do
    end do

    !--------------------------------------------------------
    !	... setup for ox conserving interp
    !--------------------------------------------------------
    if( table_ox_ndx > 0 ) then
       if( map(table_ox_ndx) > 0 ) then
          allocate( table_ox(ub_nlevs-2),stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'set_fstrat_vals: table_ox allocation error = ',astat
             call endrun
          end if
#ifdef UB_DEBUG
          write(iulog,*) ' '
          write(iulog,*) 'set_fstrat_vals: ub_nlevs = ',ub_nlevs
          write(iulog,*) 'set_fstrat_vals: ub_plevse'
          write(iulog,'(1p5g15.7)') ub_plevse(:)
          write(iulog,*) ' '
#endif
       end if
    end if

    !--------------------------------------------------------
    !	... set the mixing ratios at upper boundary
    !--------------------------------------------------------
    species_loop : do m = 1,ub_nspecies
       ub_overwrite : if( map(m) > 0 ) then
          if( m == table_ox_ndx ) then
             do i = 1,ncol

                table_ox(1:ub_nlevs-2) = mr_ub(i,m,last,2:ub_nlevs-1,lchnk) &
                     + dels*(mr_ub(i,m,next,2:ub_nlevs-1,lchnk) &
                     - mr_ub(i,m,last,2:ub_nlevs-1,lchnk))
#ifdef UB_DEBUG
                write(iulog,*) 'set_fstrat_vals: table_ox @ lat = ',lat
                write(iulog,'(1p5g15.7)') table_ox(:)
                write(iulog,*) ' '
#endif

                km = kmax(i)
#ifdef UB_DEBUG
                write(iulog,*) 'set_fstrat_vals: pint with km = ',km
                write(iulog,'(1p5g15.7)') pint(i,:km+1)
                write(iulog,*) ' '
                write(iulog,*) 'set_fstrat_vals: pmid with km = ',km
                write(iulog,'(1p5g15.7)') pmid(i,:km)
                write(iulog,*) ' '
#endif
                call rebin( ub_nlevs-2, km, ub_plevse, pint(i,:km+1), table_ox, vmr(i,:km,map(m)) )
#ifdef UB_DEBUG
                write(iulog,*) 'set_fstrat_vals: ub o3 @ lat = ',lat
                write(iulog,'(1p5g15.7)') vmr(i,:km,map(m))
#endif
             end do
             deallocate( table_ox )
             cycle species_loop
          end if
          do i = 1,ncol
             do k = 1,kmax(i)
                pint_vals(1) = mr_ub(i,m,last,kl(i,k),lchnk) &
                     + delp(i,k) &
                     * (mr_ub(i,m,last,ku(i,k),lchnk) &
                     - mr_ub(i,m,last,kl(i,k),lchnk))
                pint_vals(2) = mr_ub(i,m,next,kl(i,k),lchnk) &
                     + delp(i,k) &
                     * (mr_ub(i,m,next,ku(i,k),lchnk) &
                     - mr_ub(i,m,next,kl(i,k),lchnk))
                if( m /= table_nox_ndx .and. m /= table_h2o_ndx .and. m /= table_synoz_ndx ) then
                   vmr(i,k,map(m)) = pint_vals(1) &
                        + dels * (pint_vals(2) - pint_vals(1))
                else if( m == table_nox_ndx .and. sim_has_nox ) then

                   nox_ubc = pint_vals(1) + dels * (pint_vals(2) - pint_vals(1))
                   if( no_ndx > 0 ) then
                      xno  = vmr(i,k,no_ndx)
                   else
                      xno  = 0._r8
                   end if
                   if( no2_ndx > 0 ) then
                      xno2 = vmr(i,k,no2_ndx)
                   else
                      xno2 = 0._r8
                   end if
                   rno  = xno / (xno + xno2)
                   if( no_ndx > 0 ) then
                      vmr(i,k,no_ndx)  = rno * nox_ubc
                   end if
                   if( no2_ndx > 0 ) then
                      vmr(i,k,no2_ndx) = (1._r8 - rno) * nox_ubc
                   end if
                end if
             end do
          end do
       end if ub_overwrite
    end do species_loop

    col_loop2 : do i = 1,ncol
       !--------------------------------------------------------
       ! 	... relax lower stratosphere to extended ubc
       !           check to make sure ubc is not being imposed too low
       !           levrelax = lowest model level (highest pressure)
       !                      in which to relax to ubc
       !--------------------------------------------------------
       levrelax = ltrop(i)
       do while( pmid(i,levrelax) > ub_plevs(ub_nlevs) )
          levrelax = levrelax - 1
       end do
#ifdef DEBUG
       if( levrelax /= ltrop(i) ) then
          write(iulog,*) 'warning -- raised ubc: ',lat,i,
          ltrop(i)-1,nint(pmid(i,ltrop(i)-1)/mb2pa),'mb -->',
          levrelax,nint(pmid(i,levrelax)/mb2pa),'mb'
       end if
#endif
       level_loop2 : do k = kmax(i)+1,levrelax
          if( sim_has_nox ) then
             if( no_ndx > 0 ) then
                xno  = vmr(i,k,no_ndx)
             else
                xno  = 0._r8
             end if
             if( no2_ndx > 0 ) then
                xno2 = vmr(i,k,no2_ndx)
             else
                xno2 = 0._r8
             end if
             rno     = xno / (xno + xno2)
             nox_ubc = xno + xno2
          end if
          do m = 1,ub_nspecies
             if( map(m) < 1 ) then
                cycle
             end if
             pint_vals(1) = mr_ub(i,m,last,kl(i,k),lchnk) &
                  + delp(i,k) &
                  * (mr_ub(i,m,last,ku(i,k),lchnk) &
                  - mr_ub(i,m,last,kl(i,k),lchnk))
             pint_vals(2) = mr_ub(i,m,next,kl(i,k),lchnk) &
                  + delp(i,k) &
                  * (mr_ub(i,m,next,ku(i,k),lchnk) &
                  - mr_ub(i,m,next,kl(i,k),lchnk))
             vmrrelax = pint_vals(1) &
                  + dels * (pint_vals(2) - pint_vals(1))
             if( m /= table_nox_ndx .and. m /= table_h2o_ndx  .and. m /= table_synoz_ndx ) then
                vmr(i,k,map(m)) = vmr(i,k,map(m)) &
                     + (vmrrelax - vmr(i,k,map(m))) * facrelax
             else if( m == table_nox_ndx .and. sim_has_nox) then

                nox_ubc = nox_ubc + (vmrrelax - nox_ubc) * facrelax
             end if
          end do
          if( sim_has_nox ) then
             if( no_ndx > 0 ) then
                vmr(i,k,no_ndx)  = rno * nox_ubc
             end if
             if( no2_ndx > 0 ) then
                vmr(i,k,no2_ndx) = (1._r8 - rno) * nox_ubc
             end if
          end if
       end do level_loop2

       has_synoz : if( synoz_ndx > 0 ) then
          if ( synoz_ndx > 0 .and. table_synoz_ndx > 0 ) then
             fac_relax = 1._r8 - exp( -real(dtime) / synoz_strat_relax )
             do k = 1,pver
                m = table_synoz_ndx
                if ( synoz_region(i,k,lchnk) > 0._r8 ) then
                   pint_vals(1) = mr_ub(i,m,last,kl(i,k),lchnk) &
                        + delp(i,k) &
                        * (mr_ub(i,m,last,ku(i,k),lchnk) &
                        - mr_ub(i,m,last,kl(i,k),lchnk))
                   pint_vals(2) = mr_ub(i,m,next,kl(i,k),lchnk) &
                        + delp(i,k) &
                        * (mr_ub(i,m,next,ku(i,k),lchnk) &
                        - mr_ub(i,m,next,kl(i,k),lchnk))
                   vmrrelax = pint_vals(1) &
                        + dels * (pint_vals(2) - pint_vals(1))
                   vmr(i,k,map(m)) = vmr(i,k,map(m)) &
                        + (vmrrelax - vmr(i,k,map(m))) * fac_relax
                endif
             enddo
          endif
 
          !--------------------------------------------------------
          ! 	... special assignments if synoz is present
          !           update ox, o3s, o3inert in the stratosphere
          !--------------------------------------------------------
          if( ox_ndx > 0 ) then
             do k = 1,levrelax
                if( vmr(i,k,synoz_ndx) >= synoz_thres ) then
                   vmr(i,k,ox_ndx) = vmr(i,k,synoz_ndx)
                end if
             end do
          end if
          if( o3s_ndx > 0 ) then
             do k = 1,levrelax
                if( vmr(i,k,synoz_ndx) >= synoz_thres ) then
                   vmr(i,k,o3s_ndx) = vmr(i,k,synoz_ndx)
                end if
             end do
          end if
          if( o3rad_ndx > 0 .and. o3inert_ndx > 0 ) then
             vmr(i,:ltrop(i),o3inert_ndx) = vmr(i,:ltrop(i),o3rad_ndx)
          end if
          !--------------------------------------------------------
          ! 	... O3RAD is relaxed to climatology in the stratosphere
          !           (done above) and OX in the troposphere
          !--------------------------------------------------------
          if( o3rad_ndx > 0 .and. ox_ndx > 0 ) then
             fac_relax = 1._r8 - exp( -real(dtime) / o3rad_relax )
             do k = levrelax+1,pver
                vmr(i,k,o3rad_ndx) = vmr(i,k,o3rad_ndx) &
                     + (vmr(i,k,ox_ndx) - vmr(i,k,o3rad_ndx)) * fac_relax
             end do
          end if
          !--------------------------------------------------------
          ! 	... relax synoz to 25 ppbv in lower troposphere
          !           (p > 500 hPa) with an e-fold time of 2 days
          !           (Mc Linden et al., JGR, p14,660, 2000)
          !--------------------------------------------------------
          fac_relax = 1._r8 - exp( -real(dtime) / synoz_relax )
          vmrrelax = 25.e-9_r8
          do k = levrelax+2,pver
             if( pmid(i,k) >= 50000._r8 ) then
                vmr(i,k,synoz_ndx) = vmr(i,k,synoz_ndx) &
                     + (vmrrelax - vmr(i,k,synoz_ndx)) * fac_relax
             end if
          end do
       else has_synoz

          !--------------------------------------------------------
          !       ... set O3S and O3INERT to OX when no synoz
          !--------------------------------------------------------
          if( ox_ndx > 0 ) then
             if( o3s_ndx > 0 ) then
                vmr(i,:ltrop(i),o3s_ndx)     = vmr(i,:ltrop(i),ox_ndx)
             end if
             if( o3inert_ndx > 0 ) then
                vmr(i,:ltrop(i),o3inert_ndx) = vmr(i,:ltrop(i),ox_ndx)
             end if
          end if

       end if has_synoz

       if ( o3a_ndx > 0 ) then
          vmr(i,:ltrop(i),o3a_ndx) = (1._r8 - facrelax ) * vmr(i,:ltrop(i),o3a_ndx)
       endif
       if ( xno_ndx > 0 ) then
          vmr(i,:ltrop(i),xno_ndx) = (1._r8 - facrelax ) * vmr(i,:ltrop(i),xno_ndx)
       endif
       if ( xno2_ndx > 0 ) then
          vmr(i,:ltrop(i),xno2_ndx) = (1._r8 - facrelax ) * vmr(i,:ltrop(i),xno2_ndx)
       endif

    end do col_loop2

  end subroutine set_fstrat_vals

  subroutine set_fstrat_h2o( h2o, pmid, ltrop, calday, ncol, lchnk )
    !--------------------------------------------------------------------
    !	... set the h2o upper boundary values
    !--------------------------------------------------------------------


    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    integer,  intent(in)    :: ncol              ! columns in chunk
    integer,  intent(in)    :: ltrop(pcols)      ! tropopause vertical index
    integer,  intent(in)    :: lchnk
    real(r8), intent(in)    :: calday            ! day of year including fraction
    real(r8), intent(in)    :: pmid(pcols,pver)  ! midpoint pressure (Pa)
    real(r8), intent(inout) :: h2o(ncol,pver)    ! h2o concentration (mol/mol)

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer, parameter :: zlower = pver

    integer  ::  m, last, next, i, k, k1
    integer  ::  kmax(ncol)
    integer  ::  levrelax
    integer  ::  kl(ncol,zlower)
    integer  ::  ku(ncol,zlower)
    real(r8) ::  vmrrelax
    real(r8) ::  fac_relax
    real(r8) ::  pinterp
    real(r8) ::  dels
    real(r8) ::  delp(ncol,zlower)
    real(r8) ::  pint_vals(2)
    logical  ::  found_trop
    integer  ::  lat

    h2o_overwrite : if( h2o_ndx > 0 .and. table_h2o_ndx > 0 ) then
       !--------------------------------------------------------
       !	... setup the time interpolation
       !--------------------------------------------------------
       if( calday < days(1) ) then
          next = 1
          last = 12
          dels = (365._r8 + calday - days(12)) / (365._r8 + days(1) - days(12))
       else if( calday >= days(12) ) then
          next = 1
          last = 12
          dels = (calday - days(12)) / (365._r8 + days(1) - days(12))
       else
          do m = 11,1,-1
             if( calday >= days(m) ) then
                exit
             end if
          end do
          last = m
          next = m + 1
          dels = (calday - days(m)) / (days(m+1) - days(m))
       end if
       dels = max( min( 1._r8,dels ),0._r8 )

       !--------------------------------------------------------
       !	... setup the pressure interpolation
       !--------------------------------------------------------
       do k = 1,zlower
          do i = 1,ncol
             if( pmid(i,k) <= ub_plevs(1) ) then
                kl(i,k) = 1
                ku(i,k) = 1
                delp(i,k) = 0._r8
             else if( pmid(i,k) >= ub_plevs(ub_nlevs) ) then
                kl(i,k) = ub_nlevs
                ku(i,k) = ub_nlevs
                delp(i,k) = 0._r8
             else
                pinterp = pmid(i,k)
                do k1 = 2,ub_nlevs
                   if( pinterp <= ub_plevs(k1) ) then
                      ku(i,k) = k1
                      kl(i,k) = k1 - 1
                      delp(i,k) = log( pinterp/ub_plevs(kl(i,k)) ) &
                           / log( ub_plevs(ku(i,k))/ub_plevs(kl(i,k)) )
                      exit
                   end if
                end do
             end if
          end do
       end do

       !--------------------------------------------------------
       !	... Find max level less than 50 mb
       !           fix UB vals from top of model to this level
       !--------------------------------------------------------
       do i = 1,ncol
          do k = 2,pver
             if( pmid(i,k) > 50.e2_r8 ) then
                kmax(i) = k
                exit
             end if
          end do
       end do
       !--------------------------------------------------------
       !	... set the mixing ratio at upper boundary
       !--------------------------------------------------------
       m = table_h2o_ndx
       do i = 1,ncol
          do k = 1,kmax(i)
             pint_vals(1) = mr_ub(i,m,last,kl(i,k),lchnk) &
                  + delp(i,k) &
                  * (mr_ub(i,m,last,ku(i,k),lchnk) &
                  - mr_ub(i,m,last,kl(i,k),lchnk))
             pint_vals(2) = mr_ub(i,m,next,kl(i,k),lchnk) &
                  + delp(i,k) &
                  * (mr_ub(i,m,next,ku(i,k),lchnk) &
                  - mr_ub(i,m,next,kl(i,k),lchnk))
             h2o(i,k) = pint_vals(1) &
                  + dels * (pint_vals(2) - pint_vals(1))
          end do
       end do

       col_loop2 : do i = 1,ncol
          !--------------------------------------------------------
          ! 	... relax lower stratosphere to extended ubc
          !           check to make sure ubc is not being imposed too low
          !           levrelax = lowest model level (highest pressure)
          !                      in which to relax to ubc
          !--------------------------------------------------------
          levrelax = ltrop(i)
          do while( pmid(i,levrelax) > ub_plevs(ub_nlevs) )
             levrelax = levrelax - 1
          end do
#ifdef DEBUG
          if( levrelax /= ltrop(i) ) then
             write(iulog,*) 'warning -- raised ubc: ',lat,i,
             ltrop(i)-1,nint(pmid(i,ltrop(i)-1)/100._r8),'mb -->',
             levrelax,nint(pmid(i,levrelax)/100._r8),'mb'
          end if
#endif
          do k = kmax(i)+1,levrelax
             pint_vals(1) = mr_ub(i,m,last,kl(i,k),lchnk) &
                  + delp(i,k) &
                  * (mr_ub(i,m,last,ku(i,k),lchnk) &
                  - mr_ub(i,m,last,kl(i,k),lchnk))
             pint_vals(2) = mr_ub(i,m,next,kl(i,k),lchnk) &
                  + delp(i,k) &
                  * (mr_ub(i,m,next,ku(i,k),lchnk) &
                  - mr_ub(i,m,next,kl(i,k),lchnk))
             vmrrelax = pint_vals(1) &
                  + dels * (pint_vals(2) - pint_vals(1))
             h2o(i,k) = h2o(i,k) + (vmrrelax - h2o(i,k)) * facrelax
          end do
       end do col_loop2
    end if h2o_overwrite

  end subroutine set_fstrat_h2o

  subroutine rebin( nsrc, ntrg, src_x, trg_x, src, trg )
    !---------------------------------------------------------------
    !	... rebin src to trg
    !---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer, intent(in)   :: nsrc                  ! dimension source array
    integer, intent(in)   :: ntrg                  ! dimension target array
    real(r8), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)  :: trg_x(ntrg+1)         ! target coordinates
    real(r8), intent(in)  :: src(nsrc)             ! source array
    real(r8), intent(out) :: trg(ntrg)             ! target array

    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    integer  :: i, l
    integer  :: si, si1
    integer  :: sil, siu
    real(r8) :: y
    real(r8) :: sl, su
    real(r8) :: tl, tu

    !---------------------------------------------------------------
    !	... check interval overlap
    !---------------------------------------------------------------
    !     if( trg_x(1) < src_x(1) .or. trg_x(ntrg+1) > src_x(nsrc+1) ) then
    !        write(iulog,*) 'rebin: target grid is outside source grid'
    !        write(iulog,*) '       target grid from ',trg_x(1),' to ',trg_x(ntrg+1)
    !        write(iulog,*) '       source grid from ',src_x(1),' to ',src_x(nsrc+1)
    !        call endrun
    !     end if

    do i = 1,ntrg
       tl = trg_x(i)
       if( tl < src_x(nsrc+1) ) then
          do sil = 1,nsrc+1
             if( tl <= src_x(sil) ) then
                exit
             end if
          end do
          tu = trg_x(i+1)
          do siu = 1,nsrc+1
             if( tu <= src_x(siu) ) then
                exit
             end if
          end do
          y   = 0._r8
          sil = max( sil,2 )
          siu = min( siu,nsrc+1 )
          do si = sil,siu
             si1 = si - 1
             sl  = max( tl,src_x(si1) )
             su  = min( tu,src_x(si) )
             y   = y + (su - sl)*src(si1)
          end do
          trg(i) = y/(trg_x(i+1) - trg_x(i))
       else
          trg(i) = 0._r8
       end if
    end do

  end subroutine rebin

end module mo_fstrat
