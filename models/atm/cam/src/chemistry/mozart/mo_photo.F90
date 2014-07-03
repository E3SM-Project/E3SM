module mo_photo
  !----------------------------------------------------------------------
  !	... photolysis interp table and related arrays
  !----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, pver, pverp, begchunk, endchunk
  use abortutils,   only : endrun
  use mo_constants, only : pi,r2d,boltz,d2r
  use ref_pres,     only : num_pr_lev, ptop_ref 
  use pio
  use cam_pio_utils,only : cam_pio_openfile
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog
  use phys_control, only : waccmx_is

  implicit none

  private

  public :: photo_inti, table_photo, xactive_photo
  public :: set_ub_col
  public :: setcol 
  public :: photo_timestep_init
  public :: photo_register

  save

  real(r8), parameter :: kg2g = 1.e3_r8
  integer, parameter  :: pverm = pver - 1

  integer ::  jno_ndx
  integer ::  jonitr_ndx
  integer ::  jho2no2_ndx
  integer ::  jch3cho_a_ndx, jch3cho_b_ndx, jch3cho_c_ndx
  integer ::  jo2_a_ndx, jo2_b_ndx
  integer ::  ox_ndx, o3_ndx, o3_inv_ndx, o3rad_ndx
  integer ::  oc1_ndx, oc2_ndx
  integer ::  cb1_ndx, cb2_ndx
  integer ::  soa_ndx
  integer ::  ant_ndx
  integer ::  so4_ndx
  integer ::  sa1_ndx, sa2_ndx, sa3_ndx, sa4_ndx
  integer ::  n2_ndx, no_ndx, o2_ndx, o_ndx
  integer, allocatable :: lng_indexer(:)
  integer, allocatable :: sht_indexer(:)
  integer, allocatable :: euv_indexer(:)

  integer              :: ki
  integer              :: last
  integer              :: next
  integer              :: n_exo_levs
  real(r8)                 :: delp
  real(r8)                 :: dels
  real(r8), allocatable    :: days(:)
  real(r8), allocatable    :: levs(:)
  real(r8), allocatable    :: o2_exo_coldens(:,:,:,:)
  real(r8), allocatable    :: o3_exo_coldens(:,:,:,:)
  logical              :: o_is_inv
  logical              :: o2_is_inv
  logical              :: n2_is_inv
  logical              :: o3_is_inv
  logical              :: no_is_inv
  logical              :: has_o2_col
  logical              :: has_o3_col
  logical              :: has_fixed_press
  logical  :: waccm
  real(r8) :: max_zen_angle       ! degrees

  integer :: jo1d_ndx, jo3p_ndx, jno2_ndx, jn2o5_ndx
  integer :: jhno3_ndx, jno3_ndx, jpan_ndx, jmpan_ndx

  integer :: jo1da_ndx, jo3pa_ndx, jno2a_ndx, jn2o5a_ndx, jn2o5b_ndx
  integer :: jhno3a_ndx, jno3a_ndx, jpana_ndx, jmpana_ndx, jho2no2a_ndx 
  integer :: jonitra_ndx

  logical :: do_jshort = .false.
  logical :: do_diag = .false.

  integer :: ion_rates_idx = -1

contains

  
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine photo_register
    use mo_jeuv,      only : nIonRates
    use physics_buffer,only : pbuf_add_field, dtype_r8

    ! add photo-ionization rates to phys buffer for waccmx ionosphere module

    call pbuf_add_field('IonRates' , 'physpkg', dtype_r8, (/pcols,pver,nIonRates/), ion_rates_idx) ! Ionization rates for O+,O2+,N+,N2+,NO+ 

  endsubroutine photo_register

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine photo_inti( xs_coef_file, xs_short_file, xs_long_file, rsf_file, &
       euvacdat_file, photon_file, electron_file, &
       exo_coldens_file, tuv_xsect_file, o2_xsect_file, xactive_prates, is_waccm )
    !----------------------------------------------------------------------
    !	... initialize photolysis module
    !----------------------------------------------------------------------

    use mo_photoin,    only : photoin_inti
    use mo_tuv_inti,   only : tuv_inti
    use mo_tuv_inti,   only : nlng
    use mo_seto2,      only : o2_xsect_inti      
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use chem_mods,     only : phtcnt
    use chem_mods,     only : ncol_abs => nabscol
    use chem_mods,     only : rxt_tag_lst, pht_alias_lst, pht_alias_mult
    use time_manager,  only : get_calday
    use ioFileMod,     only : getfil
    use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx, get_inv_ndx
    use mo_jlong,      only : jlong_init
    use seasalt_model, only : sslt_names=>seasalt_names, sslt_ncnst=>seasalt_nbin
    use mo_jshort,     only : jshort_init
    use mo_jeuv,       only : jeuv_init
    use dyn_grid,      only : get_dyn_grid_parm
    use phys_grid,     only : get_ncols_p, get_rlat_all_p    
    use solar_data,    only : has_spectrum

    implicit none

    !----------------------------------------------------------------------
    !	... dummy arguments
    !----------------------------------------------------------------------
    character(len=*), intent(in) :: xs_long_file, rsf_file
    character(len=*), intent(in) :: exo_coldens_file
    character(len=*), intent(in) :: tuv_xsect_file
    character(len=*), intent(in) :: o2_xsect_file
    logical, intent(in)          :: xactive_prates
    ! waccm 
    character(len=*), intent(in) :: xs_coef_file
    character(len=*), intent(in) :: xs_short_file
    character(len=*), intent(in) :: euvacdat_file
    character(len=*), intent(in) :: photon_file
    character(len=*), intent(in) :: electron_file

    logical, optional, intent(in) :: is_waccm

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    real(r8), parameter   :: hPa2Pa = 100._r8
    integer           :: k, n
    type(file_desc_t) :: ncid
    type(var_desc_t)  :: vidO2, vidO3, vid
    type(interp_type) :: lat_wgts
    integer           :: dimid
    integer           :: nlat
    integer           :: ntimes
    integer           :: astat
    integer           :: gndx
    integer           :: ndx
    integer           :: spc_ndx
    integer           :: ierr
    integer           :: c, ncols
    integer           :: plev, plevp
    integer, allocatable :: dates(:)
    real(r8)              :: pinterp
    real(r8), allocatable :: lats(:)
    real(r8), allocatable :: coldens(:,:,:)
    character(len=256)    :: locfn
    character(len=256)    :: filespec
    real(r8), parameter :: trop_thrshld = 1._r8 ! Pa
    real(r8) :: to_lats(pcols)


    if( phtcnt < 1 ) then
       return
    end if

    if (present(is_waccm)) then
       waccm = is_waccm
    else
       waccm = .false.
    endif
    
    !----------------------------------------------------------------------------
    !  Need a larger maximum zenith angle for WACCM-X extended to high altitudes
    !----------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
       max_zen_angle = 116._r8
    else if ( waccm ) then
       max_zen_angle = 97.01_r8 ! degrees
    else
       max_zen_angle = 88.85_r8 ! degrees
    endif

    ! jeuv_1,,, jeuv_25 --> need euv calculations --> must be waccm 
    ! how to determine if shrt calc is needed ?? -- use top level pressure => waccm = true ? false

    if ( waccm .and. xactive_prates ) then
       call endrun( " photo_inti: xactive_prates option is not available in WACCM" )
    endif

    if ( .not. has_spectrum ) then
       write(iulog,*) 'photo_inti: solar_data file needs to contain irradiance spectrum'
       call endrun('photo_inti: ERROR -- solar irradiance spectrum is missing')
    endif
    
    plev = get_dyn_grid_parm('plev')
    plevp = get_dyn_grid_parm('plevp')

    !----------------------------------------------------------------------
    !	... allocate indexers
    !----------------------------------------------------------------------
    allocate( lng_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: lng_indexer allocation error = ',astat
       call endrun
    end if
    lng_indexer(:) = 0
    allocate( sht_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: Failed to allocate sht_indexer; error = ',astat
       call endrun
    end if
    sht_indexer(:) = 0
    allocate( euv_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: Failed to allocate euv_indexer; error = ',astat
       call endrun
    end if
    euv_indexer(:) = 0

    jno_ndx     = get_rxt_ndx( 'jno' )
    jo2_a_ndx   = get_rxt_ndx( 'jo2_a' )
    jo2_b_ndx   = get_rxt_ndx( 'jo2_b' )

    jo1da_ndx = get_rxt_ndx( 'jo1da' )
    jo3pa_ndx = get_rxt_ndx( 'jo3pa' )
    jno2a_ndx = get_rxt_ndx( 'jno2a' )
    jn2o5a_ndx = get_rxt_ndx( 'jn2o5a' )
    jn2o5b_ndx = get_rxt_ndx( 'jn2o5b' )
    jhno3a_ndx = get_rxt_ndx( 'jhno3a' )
    jno3a_ndx = get_rxt_ndx( 'jno3a' )
    jpana_ndx = get_rxt_ndx( 'jpana' )
    jmpana_ndx = get_rxt_ndx( 'jmpana' )
    jho2no2a_ndx  = get_rxt_ndx( 'jho2no2a' )
    jonitra_ndx = get_rxt_ndx( 'jonitra' )

    jo1d_ndx = get_rxt_ndx( 'jo1d' )
    jo3p_ndx = get_rxt_ndx( 'jo3p' )
    jno2_ndx = get_rxt_ndx( 'jno2' )
    jn2o5_ndx = get_rxt_ndx( 'jn2o5' )
    jn2o5_ndx = get_rxt_ndx( 'jn2o5' )
    jhno3_ndx = get_rxt_ndx( 'jhno3' )
    jno3_ndx = get_rxt_ndx( 'jno3' )
    jpan_ndx = get_rxt_ndx( 'jpan' )
    jmpan_ndx = get_rxt_ndx( 'jmpan' )
    jho2no2_ndx  = get_rxt_ndx( 'jho2no2' )
    jonitr_ndx = get_rxt_ndx( 'jonitr' )


    ox_ndx     = get_spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx  = get_spc_ndx( 'O3' )
    end if
    o3_ndx     = get_spc_ndx( 'O3' )
    o3rad_ndx  = get_spc_ndx( 'O3RAD' )
    o3_inv_ndx = get_inv_ndx( 'O3' )

    n2_ndx     = get_inv_ndx( 'N2' )
    n2_is_inv  = n2_ndx > 0
    if( .not. n2_is_inv ) then
       n2_ndx = get_spc_ndx( 'N2' )
    end if
    o2_ndx     = get_inv_ndx( 'O2' )
    o2_is_inv  = o2_ndx > 0
    if( .not. o2_is_inv ) then
       o2_ndx = get_spc_ndx( 'O2' )
    end if
    no_ndx     = get_spc_ndx( 'NO' )
    no_is_inv  = no_ndx < 1
    if( no_is_inv ) then
       no_ndx = get_inv_ndx( 'NO' )
    end if
    o3_is_inv  = o3_ndx < 1

    o_ndx     = get_spc_ndx( 'O' )
    o_is_inv  = o_ndx < 1
    if( o_is_inv ) then
       o_ndx = get_inv_ndx( 'O' )
    end if

    do_jshort = o_ndx>0 .and. o2_ndx>0 .and. (o3_ndx>0.or.o3_inv_ndx>0) .and. n2_ndx>0 .and. no_ndx>0

    !----------------------------------------------------------------------
    !	... call module initializers
    !----------------------------------------------------------------------
    is_xactive : if( xactive_prates ) then
       do_jshort = .false.
       jch3cho_a_ndx = get_rxt_ndx( 'jch3cho_a' )
       jch3cho_b_ndx = get_rxt_ndx( 'jch3cho_b' )
       jch3cho_c_ndx = get_rxt_ndx( 'jch3cho_c' )
       jonitr_ndx    = get_rxt_ndx( 'jonitr' )
       jho2no2_ndx   = get_rxt_ndx( 'jho2no2' )
       call tuv_inti( pverp, tuv_xsect_file, lng_indexer )
    else is_xactive
       call jlong_init( xs_long_file, rsf_file, lng_indexer )
       if (waccm) then
          call jeuv_init( euvacdat_file, photon_file, electron_file, euv_indexer )
       endif
       if (do_jshort) then
          call jshort_init( xs_coef_file, xs_short_file, sht_indexer )
       endif
       jho2no2_ndx = get_rxt_ndx( 'jho2no2_b' )
    end if is_xactive

    !----------------------------------------------------------------------
    !        ... check that each photorate is in short or long datasets
    !----------------------------------------------------------------------
    if( any( (abs(sht_indexer(:)) + abs(lng_indexer(:)) + abs(euv_indexer(:))) == 0 ) ) then
       write(iulog,*) ' '
       write(iulog,*) 'photo_inti: the following photorate(s) are not in'
       write(iulog,*) '            either the short or long datasets'
       write(iulog,*) ' '
       do ndx = 1,phtcnt
          if( abs(sht_indexer(ndx)) + abs(lng_indexer(ndx)) == 0 ) then
             write(iulog,*) '           ',trim( rxt_tag_lst(ndx) )
          end if
       end do
       call endrun
    end if

    !----------------------------------------------------------------------
    !        ... output any aliased photorates
    !----------------------------------------------------------------------
    if( masterproc ) then
       if( any( pht_alias_lst(:,1) /= ' ' ) ) then
          write(iulog,*) ' '
          write(iulog,*) 'photo_inti: the following short photorate(s) are aliased'
          write(iulog,*) ' '
          do ndx = 1,phtcnt
             if( pht_alias_lst(ndx,1) /= ' ' ) then
                if( pht_alias_mult(ndx,1) == 1._r8 ) then
                   write(iulog,*) '           ',trim(rxt_tag_lst(ndx)),' -> ',trim(pht_alias_lst(ndx,1))
                else
                   write(iulog,*) '           ',trim(rxt_tag_lst(ndx)),' -> ',pht_alias_mult(ndx,1),'*',trim(pht_alias_lst(ndx,1))
                end if
             end if
          end do
       end if
       if( any( pht_alias_lst(:,2) /= ' ' ) ) then
          write(iulog,*) ' '
          write(iulog,*) 'photo_inti: the following long photorate(s) are aliased'
          write(iulog,*) ' '
          do ndx = 1,phtcnt
             if( pht_alias_lst(ndx,2) /= ' ' ) then
                if( pht_alias_mult(ndx,2) == 1._r8 ) then
                   write(iulog,*) '           ',trim(rxt_tag_lst(ndx)),' -> ',trim(pht_alias_lst(ndx,2))
                else
                   write(iulog,*) '           ',trim(rxt_tag_lst(ndx)),' -> ',pht_alias_mult(ndx,2),'*',trim(pht_alias_lst(ndx,2))
                end if
             end if
          end do
       end if

       write(iulog,*) ' '
       write(iulog,*) '*********************************************'
       write(iulog,*) 'photo_inti: euv_indexer'
       write(iulog,'(10i6)') euv_indexer(:)
       write(iulog,*) 'photo_inti: sht_indexer'
       write(iulog,'(10i6)') sht_indexer(:)
       write(iulog,*) 'photo_inti: lng_indexer'
       write(iulog,'(10i6)') lng_indexer(:)
       write(iulog,*) '*********************************************'
       write(iulog,*) ' '
    endif

    if( xactive_prates ) then
       call o2_xsect_inti( o2_xsect_file )
       call photoin_inti( nlng, lng_indexer )
    end if

    !----------------------------------------------------------------------
    !	... check for o2, o3 absorber columns
    !----------------------------------------------------------------------
    if( ncol_abs > 0 ) then
       spc_ndx = ox_ndx
       if( spc_ndx < 1 ) then
          spc_ndx = o3_ndx
       end if
       if( spc_ndx > 0 ) then
          has_o3_col = .true.
       else
          has_o3_col = .false.
       end if
       if( ncol_abs > 1 ) then
          if( o2_ndx > 1 ) then
             has_o2_col = .true.
          else
             has_o2_col = .false.
          end if
       else
          has_o2_col = .false.
       end if
    else
       has_o2_col = .false.
       has_o3_col = .false.
    end if

    if ( len_trim(exo_coldens_file) == 0 ) then
       has_o2_col = .false.
       has_o3_col = .false.
    endif

    oc1_ndx = get_spc_ndx( 'OC1' )
    oc2_ndx = get_spc_ndx( 'OC2' )
    cb1_ndx = get_spc_ndx( 'CB1' )
    cb2_ndx = get_spc_ndx( 'CB2' )
    soa_ndx = get_spc_ndx( 'SOA' )
    ant_ndx = get_spc_ndx( 'NH4NO3' )
    so4_ndx = get_spc_ndx( 'SO4' )
    if (sslt_ncnst == 4) then
       sa1_ndx = get_spc_ndx( sslt_names(1) )
       sa2_ndx = get_spc_ndx( sslt_names(2) )
       sa3_ndx = get_spc_ndx( sslt_names(3) )
       sa4_ndx = get_spc_ndx( sslt_names(4) )
    endif

    has_abs_columns : if( has_o2_col .or. has_o3_col ) then
       !-----------------------------------------------------------------------
       !	... open exo coldens file
       !-----------------------------------------------------------------------
       filespec = trim( exo_coldens_file )
       call getfil( filespec, locfn, 0 )
       call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE )

       !-----------------------------------------------------------------------
       !       ... get grid dimensions from file
       !-----------------------------------------------------------------------
       !       ... timing
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'month', dimid )
       ierr = pio_inq_dimlen( ncid, dimid, ntimes )

       if( ntimes /= 12 ) then
          call endrun('photo_inti: exo coldens is not annual period')
       end if
       allocate( dates(ntimes),days(ntimes),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: dates,days allocation error = ',astat
          call endrun
       end if
       dates(:) = (/ 116, 214, 316, 415,  516,  615, &
            716, 816, 915, 1016, 1115, 1216 /)
       !-----------------------------------------------------------------------
       !	... initialize the monthly day of year times
       !-----------------------------------------------------------------------
       do n = 1,ntimes
          days(n) = get_calday( dates(n), 0 )
       end do
       deallocate( dates )
       !-----------------------------------------------------------------------
       !       ... latitudes
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'lat', dimid )
       ierr = pio_inq_dimlen( ncid, dimid, nlat )
       allocate( lats(nlat), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: lats allocation error = ',astat
          call endrun
       end if
       ierr = pio_inq_varid( ncid, 'lat', vid )
       ierr = pio_get_var( ncid, vid, lats )
       lats(:nlat) = lats(:nlat) * d2r
       !-----------------------------------------------------------------------
       !       ... levels
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'lev', dimid )
       ierr = pio_inq_dimlen( ncid, dimid, n_exo_levs )
       allocate( levs(n_exo_levs), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: levs allocation error = ',astat
          call endrun
       end if
       ierr = pio_inq_varid( ncid, 'lev', vid )
       ierr = pio_get_var( ncid, vid, levs )
       levs(:n_exo_levs) = levs(:n_exo_levs) * hPa2Pa
       !-----------------------------------------------------------------------
       !       ... set up regridding
       !-----------------------------------------------------------------------

       allocate( coldens(nlat,n_exo_levs,ntimes),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: coldens allocation error = ',astat
          call endrun
       end if
       if( has_o2_col ) then
          allocate( o2_exo_coldens(n_exo_levs,pcols,begchunk:endchunk,ntimes),stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'photo_inti: o2_exo_coldens allocation error = ',astat
             call endrun
          end if
          ierr = pio_inq_varid( ncid, 'O2_column_density', vid )
          ierr = pio_get_var( ncid, vid,coldens )

          do c=begchunk,endchunk
             ncols = get_ncols_p(c)
             call get_rlat_all_p(c, pcols, to_lats)
             call lininterp_init(lats, nlat, to_lats, ncols, 1, lat_wgts)
             do n=1,ntimes
                do k = 1,n_exo_levs
                   call lininterp(coldens(:,k,n), nlat, o2_exo_coldens(k,:,c,n), ncols, lat_wgts)
                end do
             end do
             call lininterp_finish(lat_wgts)
          enddo


       end if
       if( has_o3_col ) then
          allocate( o3_exo_coldens(n_exo_levs,pcols,begchunk:endchunk,ntimes),stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'photo_inti: o3_exo_coldens allocation error = ',astat
             call endrun
          end if
          ierr = pio_inq_varid( ncid, 'O3_column_density', vid )
          ierr = pio_get_var( ncid, vid,coldens )

          do c=begchunk,endchunk
             ncols = get_ncols_p(c)
             call get_rlat_all_p(c, pcols, to_lats)
             call lininterp_init(lats, nlat, to_lats, ncols, 1, lat_wgts)
             do n=1,ntimes
                do k = 1,n_exo_levs
                   call lininterp(coldens(:,k,n), nlat, o3_exo_coldens(k,:,c,n), ncols, lat_wgts)
                end do
             end do
             call lininterp_finish(lat_wgts)
          enddo
       end if
       call pio_closefile (ncid)
       deallocate( coldens,stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo_inti: failed to deallocate coldens; error = ',astat
          call endrun
       end if
       has_fixed_press = (num_pr_lev .ne. 0)
       !-----------------------------------------------------------------------
       !	... setup the pressure interpolation
       !-----------------------------------------------------------------------
       if( has_fixed_press ) then
          pinterp =  ptop_ref
          if( pinterp <= levs(1) ) then
             ki   = 1
             delp = 0._r8
          else
             do ki = 2,n_exo_levs
                if( pinterp <= levs(ki) ) then
                   delp = log( pinterp/levs(ki-1) )/log( levs(ki)/levs(ki-1) )
                   exit
                end if
             end do
          end if
#ifdef DEBUG
          write(iulog,*) '-----------------------------------'
          write(iulog,*) 'photo_inti: diagnostics'
          write(iulog,*) 'ki, delp = ',ki,delp
          write(iulog,*) 'pinterp,levs(ki-1:ki) = ',pinterp,levs(ki-1:ki)
          write(iulog,*) '-----------------------------------'
#endif
       end if
    end if has_abs_columns

  end subroutine photo_inti

  subroutine table_photo( photos, pmid, pdel, temper, zmid, zint, &
                          col_dens, zen_angle, srf_alb, lwc, clouds, &
                          esfact, vmr, invariants, ncol, lchnk, pbuf )
!-----------------------------------------------------------------
!   	... table photorates for wavelengths > 200nm
!-----------------------------------------------------------------

    use chem_mods,   only : ncol_abs => nabscol, phtcnt, gas_pcnst, nfs
    use chem_mods,   only : pht_alias_mult, indexm
    use mo_jshort,   only : nsht => nj, jshort
    use mo_jlong,    only : nlng => numj, jlong
    use mo_jeuv,     only : neuv, jeuv, nIonRates
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none

!-----------------------------------------------------------------
!   	... dummy arguments
!-----------------------------------------------------------------
    integer,  intent(in)    :: lchnk
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: esfact                       ! earth sun distance factor
    real(r8), intent(in)    :: vmr(ncol,pver,max(1,gas_pcnst)) ! vmr
    real(r8), intent(in)    :: col_dens(ncol,pver,ncol_abs) ! column densities (molecules/cm^2)
    real(r8), intent(in)    :: zen_angle(ncol)              ! solar zenith angle (radians)
    real(r8), intent(in)    :: srf_alb(pcols)               ! surface albedo
    real(r8), intent(in)    :: lwc(ncol,pver)               ! liquid water content (kg/kg)
    real(r8), intent(in)    :: clouds(ncol,pver)            ! cloud fraction
    real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
    real(r8), intent(in)    :: pdel(pcols,pver)             ! pressure delta about midpoint (Pa)
    real(r8), intent(in)    :: temper(pcols,pver)           ! midpoint temperature (K)
    real(r8), intent(in)    :: zmid(ncol,pver)              ! midpoint height (km)
    real(r8), intent(in)    :: zint(ncol,pver)              ! interface height (km)
    real(r8), intent(in)    :: invariants(ncol,pver,max(1,nfs)) ! invariant densities (molecules/cm^3)
    real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)
    type(physics_buffer_desc),pointer :: pbuf(:)

!-----------------------------------------------------------------
!    	... local variables
!-----------------------------------------------------------------
    real(r8), parameter :: Pa2mb         = 1.e-2_r8       ! pascals to mb

    integer ::  i, k, m, n                 ! indicies
    integer ::  astat
    integer ::  indxIR                     ! pbuf index for ionization rates
    real(r8) ::  sza
    real(r8) ::  alias_factor
    real(r8) ::  fac1(pver)                ! work space for j(no) calc
    real(r8) ::  fac2(pver)                ! work space for j(no) calc
    real(r8) ::  colo3(pver)               ! vertical o3 column density
    real(r8) ::  parg(pver)                ! vertical pressure array (hPa)

    real(r8) ::  cld_line(pver)            ! vertical cloud array
    real(r8) ::  lwc_line(pver)            ! vertical lwc array
    real(r8) ::  eff_alb(pver)             ! effective albedo from cloud modifications
    real(r8) ::  cld_mult(pver)            ! clould multiplier
    real(r8) ::  tmp(ncol,pver)            ! wrk array
    real(r8), allocatable ::  lng_prates(:,:) ! photorates matrix (1/s)
    real(r8), allocatable ::  sht_prates(:,:) ! photorates matrix (1/s)
    real(r8), allocatable ::  euv_prates(:,:) ! photorates matrix (1/s)
    

    real(r8), allocatable :: zarg(:)
    real(r8), allocatable :: tline(:)               ! vertical temperature array
    real(r8), allocatable :: o_den(:)               ! o density (molecules/cm^3)
    real(r8), allocatable :: o2_den(:)              ! o2 density (molecules/cm^3)
    real(r8), allocatable :: o3_den(:)              ! o3 density (molecules/cm^3)
    real(r8), allocatable :: no_den(:)              ! no density (molecules/cm^3)
    real(r8), allocatable :: n2_den(:)              ! n2 density (molecules/cm^3)
    real(r8), allocatable :: jno_sht(:)             ! no short photorate
    real(r8), allocatable :: jo2_sht(:,:)           ! o2 short photorate

    integer :: n_jshrt_levs, p1, p2
    real(r8) :: ideltaZkm, factor

    if( phtcnt < 1 ) then
       return
    end if

    if ((.not.do_jshort) .or. (waccm)) then
       n_jshrt_levs = pver
       p1 = 1 
       p2 = pver
    else
       n_jshrt_levs = pver+1
       p1 = 2
       p2 = pver+1
    endif

    allocate( zarg(n_jshrt_levs) )
    allocate( tline(n_jshrt_levs) )
    if (do_jshort) then
       allocate( o_den(n_jshrt_levs) )
       allocate( o2_den(n_jshrt_levs) )
       allocate( o3_den(n_jshrt_levs) )
       allocate( no_den(n_jshrt_levs) )
       allocate( n2_den(n_jshrt_levs) )
       allocate( jno_sht(n_jshrt_levs) )
       allocate( jo2_sht(n_jshrt_levs,2) )
    endif

!-----------------------------------------------------------------
!	... allocate short, long temp arrays
!-----------------------------------------------------------------
    if ( waccm ) then
       if (neuv>0) then
          allocate( euv_prates(pver,neuv),stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'photo: Failed to allocate euv_prates; error = ',astat
             call endrun
          end if
       endif
    endif
    
    if (nsht>0) then
       allocate( sht_prates(n_jshrt_levs,nsht),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo: Failed to allocate sht_prates; error = ',astat
          call endrun
       end if
    endif

    if (nlng>0) then
       allocate( lng_prates(nlng,pver),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photo: Failed to allocate lng_prates; error = ',astat
          call endrun
       end if
    endif

!------------------------------------------------------------------------------------------------------------
!  Point to production rates array in physics buffer where rates will be stored for ionosphere module
!  access.  Also, initialize variables to zero before column loop
!------------------------------------------------------------------------------------------------------------
    if (ion_rates_idx>0) then
       call pbuf_set_field(pbuf, ion_rates_idx, 0._r8)
    endif

!-----------------------------------------------------------------
!	... zero all photorates
!-----------------------------------------------------------------
    do m = 1,max(1,phtcnt)
       do k = 1,pver
          photos(:,k,m) = 0._r8
       end do
    end do

    col_loop : do i = 1,ncol
       sza = zen_angle(i)*r2d
       daylight : if( sza >= 0._r8 .and. sza < max_zen_angle ) then
          parg(:)     = Pa2mb*pmid(i,:)
          colo3(:)    = col_dens(i,:,1)
          fac1(:)     = pdel(i,:)
          lwc_line(:) = lwc(i,:)
          cld_line(:) = clouds(i,:)

          
          tline(p1:p2) = temper(i,:pver)

          zarg(p1:p2) = zmid(i,:pver)

          if (do_jshort) then
             if ( o_is_inv ) then
                o_den(p1:p2) = invariants(i,:pver,o_ndx)
             else
                o_den(p1:p2) = vmr(i,:pver,o_ndx) * invariants(i,:pver,indexm)
             endif
             if ( o2_is_inv ) then
                o2_den(p1:p2) = invariants(i,:pver,o2_ndx)
             else
                o2_den(p1:p2) = vmr(i,:pver,o2_ndx) * invariants(i,:pver,indexm)
             endif
             if ( o3_is_inv ) then
                o3_den(p1:p2) = invariants(i,:pver,o3_inv_ndx)
             else
                o3_den(p1:p2) = vmr(i,:,o3_ndx) * invariants(i,:pver,indexm)
             endif
             if ( n2_is_inv ) then
                n2_den(p1:p2) = invariants(i,:,n2_ndx)
             else
                n2_den(p1:p2) = vmr(i,:pver,n2_ndx) * invariants(i,:pver,indexm)
             endif
             if ( no_is_inv ) then
                no_den(p1:p2) = invariants(i,:pver,no_ndx)
             else
                no_den(p1:p2) = vmr(i,:pver,no_ndx) * invariants(i,:pver,indexm)
             endif

             if ( .not. waccm ) then
                !-----------------------------------------------------------------
                ! Only for lower lid versions of CAM (i.e., not for WACCM)
                ! Column O3 and O2 above the top of the model
                ! DEK 20110224
                !-----------------------------------------------------------------
                ideltaZkm = 1._r8/(zint(i,1) - zint(i,2))

                !-----------------------------------------------------------------
                !... set density (units: molecules cm-3)
                !... used for jshort
                !....... Assuming a scale height of 7km for ozone
                !....... Assuming a scale height of 7km for O2
                !-----------------------------------------------------------------
                o3_den(1)        = o3_den(2)*7.0_r8 * ideltaZkm

                o2_den(1)        = o2_den(2)*7.0_r8 * ideltaZkm

                no_den(1)        = no_den(2)*0.9_r8

                n2_den(1)        = n2_den(2)*0.9_r8

                tline(1)         = tline(2) + 5.0_r8

                zarg(1)          = zarg(2) + (zint(i,1) - zint(i,2))

             endif

             !-----------------------------------------------------------------
             !	... short wave length component
             !-----------------------------------------------------------------
             call jshort( n_jshrt_levs, sza, n2_den, o2_den, o3_den, &
                  no_den, tline, zarg, jo2_sht, jno_sht, sht_prates )  

             do m = 1,phtcnt
                if( sht_indexer(m) > 0 ) then
                   alias_factor = pht_alias_mult(m,1)
                   if( alias_factor == 1._r8 ) then
                      photos(i,pver:1:-1,m) = sht_prates(1:pver,sht_indexer(m))
                   else
                      photos(i,pver:1:-1,m) = alias_factor * sht_prates(1:pver,sht_indexer(m))
                   end if
                end if
             end do

             if( jno_ndx > 0 )   photos(i,pver:1:-1,jno_ndx)   = jno_sht(1:pver)
             if( jo2_a_ndx > 0 ) photos(i,pver:1:-1,jo2_a_ndx) = jo2_sht(1:pver,2)
             if( jo2_b_ndx > 0 ) photos(i,pver:1:-1,jo2_b_ndx) = jo2_sht(1:pver,1)
          endif

          if ( waccm ) then
             !-----------------------------------------------------------------
             !	... euv photorates do not include cloud effects ??
             !-----------------------------------------------------------------
             call jeuv( pver, sza, o_den, o2_den, n2_den,  zarg, euv_prates )
             do m = 1,phtcnt
                if( euv_indexer(m) > 0 ) then
                   photos(i,:,m) = esfact * euv_prates(:,euv_indexer(m))
                endif
             enddo
          endif

          !-----------------------------------------------------------------
          !     ... compute eff_alb and cld_mult -- needs to be before jlong
          !-----------------------------------------------------------------
          call cloud_mod( zen_angle(i), cld_line, lwc_line, fac1, srf_alb(i), &
                          eff_alb, cld_mult )
          cld_mult(:) = esfact * cld_mult(:)

          !-----------------------------------------------------------------
          !	... long wave length component
          !-----------------------------------------------------------------
          call jlong( pver, sza, eff_alb, parg, tline, colo3, lng_prates )          
          do m = 1,phtcnt
             if( lng_indexer(m) > 0 ) then
                alias_factor = pht_alias_mult(m,2)
                if( alias_factor == 1._r8 ) then
                   photos(i,:,m) = (photos(i,:,m) + lng_prates(lng_indexer(m),:))*cld_mult(:)
                else
                   photos(i,:,m) = (photos(i,:,m) + alias_factor * lng_prates(lng_indexer(m),:))*cld_mult(:)
                end if
             end if
          end do

          !-----------------------------------------------------------------
          !	... calculate j(no) from formula
          !-----------------------------------------------------------------
          if( (jno_ndx > 0) .and. (.not.do_jshort)) then
             if( has_o2_col .and. has_o3_col ) then
                fac1(:) = 1.e-8_r8 * (abs(col_dens(i,:,2)/cos(zen_angle(i))))**.38_r8
                fac2(:) = 5.e-19_r8 * abs(col_dens(i,:,1)/cos(zen_angle(i)))
                photos(i,:,jno_ndx) = photos(i,:,jno_ndx) + 4.5e-6_r8 * exp( -(fac1(:) + fac2(:)) )
             end if
          end if

          !-----------------------------------------------------------------
          ! 	... add near IR correction to ho2no2
          !-----------------------------------------------------------------
          if( jho2no2_ndx > 0 ) then
             photos(i,:,jho2no2_ndx) = photos(i,:,jho2no2_ndx) + 1.e-5_r8*cld_mult(:)
          endif

          !  Save photo-ionization rates to physics buffer accessed in ionosphere module for WACCMX
          if (ion_rates_idx>0) then
                call pbuf_set_field(pbuf, ion_rates_idx, esfact * euv_prates(pver:1:-1,1:nIonRates), &
                                                         start=(/i,1,1/), kount=(/1,pver,nIonRates/))
          endif

       end if daylight
    end do col_loop

    if ( allocated(lng_prates) ) deallocate( lng_prates )
    if ( allocated(sht_prates) ) deallocate( sht_prates )
    if ( allocated(euv_prates) ) deallocate( euv_prates )

    if ( allocated(zarg) )    deallocate( zarg )
    if ( allocated(tline) )   deallocate( tline )
    if ( allocated(o_den) )   deallocate( o_den )
    if ( allocated(o2_den) )  deallocate( o2_den )
    if ( allocated(o3_den) )  deallocate( o3_den )
    if ( allocated(no_den) )  deallocate( no_den )
    if ( allocated(n2_den) )  deallocate( n2_den )
    if ( allocated(jno_sht) ) deallocate( jno_sht )
    if ( allocated(jo2_sht) ) deallocate( jo2_sht )

    call set_xnox_photo( photos, ncol  )

  end subroutine table_photo

  subroutine xactive_photo( photos, vmr, temper, cwat, cldfr, &
                            pmid, zmid, col_dens, zen_angle, srf_alb, &
                            tdens, ps, ts, esfact, relhum, dust_vmr, &
                            dt_diag, fracday, &
                            ncol, lchnk )
    !-----------------------------------------------------------------
    !   	... fast online photo rates
    !-----------------------------------------------------------------

    use ppgrid,       only : pver, pverp
    use chem_mods,    only : ncol_abs => nabscol, phtcnt, pcnstm1 => gas_pcnst
    use chem_mods,    only : pht_alias_mult
    use mo_params,    only : kz, kw
    use mo_wavelen,   only : nw
    use mo_photoin,   only : photoin
    use mo_tuv_inti,  only : nlng
    use time_manager, only : get_curr_date
    use dust_model,   only : ndst => dust_nbin
    use phys_grid,    only : get_rlat_all_p, get_rlon_all_p

    implicit none

    !----------------------------------------------------------------
    !   	... dummy arguments
    !-----------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: esfact                       ! earth sun distance factor
    real(r8), intent(in)    :: ps(pcols)                    ! surface pressure (Pa)
    real(r8), intent(in)    :: ts(ncol)                     ! surface temperature (K)
    real(r8), intent(in)    :: col_dens(ncol,pver,ncol_abs) ! column densities (molecules/cm^2)
    real(r8), intent(in)    :: zen_angle(ncol)              ! solar zenith angle (radians)
    real(r8), intent(in)    :: srf_alb(pcols)               ! surface albedo
    real(r8), intent(in)    :: tdens(ncol,pver)             ! total atms density (molecules/cm^3)
    real(r8), intent(in)    :: vmr(ncol,pver,pcnstm1)       ! species concentration (mol/mol)
    real(r8), intent(in)    :: pmid(pcols,pver)             ! midpoint pressure (Pa)
    real(r8), intent(in)    :: zmid(ncol,pver)              ! midpoint height (m)
    real(r8), intent(in)    :: temper(pcols,pver)           ! midpoint temperature (K)
    real(r8), intent(in)    :: relhum(ncol,pver)            ! relative humidity
    real(r8), intent(in)    :: cwat(ncol,pver)              ! cloud water (kg/kg)
    real(r8), intent(in)    :: cldfr(ncol,pver)             ! cloud fraction
    real(r8), intent(in)    :: dust_vmr(ncol,pver,ndst)     ! dust concentration (mol/mol)
    real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)
    real(r8), intent(out)   :: dt_diag(pcols,8)              ! od diagnostics
    real(r8), intent(out)   :: fracday(pcols)                ! fraction of day
    !-----------------------------------------------------------------
    !    	... local variables
    !-----------------------------------------------------------------
    integer, parameter  ::  k_diag = 3
    real(r8), parameter :: secant_limit = 50._r8

    integer  ::  astat
    integer  ::  i                      ! index
    integer  ::  k                      ! index
    integer  ::  m                      ! index
    integer  ::  ndx                    ! index
    integer  ::  spc_ndx                ! index
    integer  ::  file                   ! index
    integer  ::  yr, mon, day, tod      ! time of day (seconds past 0Z)
    integer  ::  ncdate                 ! current date(yyyymmdd)

    real(r8)    ::   sza
    real(r8)    ::   secant
    real(r8)    ::   alias_factor
    real(r8)    ::   alat
    real(r8)    ::   along
    real(r8)    ::   ut
    real(r8)    ::   fac1(pver)                    ! work space for j(no) calc
    real(r8)    ::   fac2(pver)                    ! work space for j(no) calc
    real(r8)    ::   tlay(pver)                    ! vertical temperature array at layer midpoint
    real(r8)    ::   tline(pverp)                  ! vertical temperature array
    real(r8)    ::   xlwc(pverp)                   ! cloud water (kg/kg)
    real(r8)    ::   xfrc(pverp)                   ! cloud fraction      xuexi
    real(r8)    ::   airdens(pverp)                ! atmospheric density
    real(r8)    ::   o3line(pverp)                 ! vertical o3 vmr
    real(r8)    ::   aerocs1(pverp)   
    real(r8)    ::   aerocs2(pverp)   
    real(r8)    ::   aercbs1(pverp)   
    real(r8)    ::   aercbs2(pverp)   
    real(r8)    ::   aersoa(pverp)   
    real(r8)    ::   aerant(pverp)   
    real(r8)    ::   aerso4(pverp)   
    real(r8)    ::   aerds(4,pverp)
    real(r8)    ::   rh(pverp)   
    real(r8)    ::   zarg(pverp)                   ! vertical height array
    real(r8)    ::   aersal(pverp,4)
    real(r8)    ::   albedo(kw)                    ! wavelenght dependent albedo
    real(r8)    ::   dt_xdiag(8)                   ! wrk array
    real(r8), allocatable :: prates(:,:)           ! photorates matrix

    logical  ::  zagtz(ncol)                       ! zenith angle > 0 flag array
    real(r8), dimension(ncol)  :: rlats, rlons     ! chunk latitudes and longitudes (radians)

    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )

    !-----------------------------------------------------------------
    !	... any photorates ?
    !-----------------------------------------------------------------
    if( phtcnt < 1 ) then
       return
    end if

    !-----------------------------------------------------------------
    !	... zero all photorates
    !-----------------------------------------------------------------
    do m = 1,phtcnt
       do k = 1,pver
          photos(:,k,m) = 0._r8
       end do
    end do

!-----------------------------------------------------------------
!	... allocate "working" rate array
!-----------------------------------------------------------------
      allocate( prates(pverp,nlng), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'xactive_photo: failed to allocate prates; error = ',astat
         call endrun
      end if

    zagtz(:) = zen_angle(:) < .99_r8*pi/2._r8 .and. zen_angle(:) > 0._r8 !! daylight
    fracday(:) = 0._r8
    dt_diag(:,:) = 0._r8

    call get_curr_date(yr, mon, day, tod, 0)
    ncdate = yr*10000 + mon*100 + day
    ut   = real(tod)/3600._r8
#ifdef DEBUG
    write(iulog,*) 'photo: nj = ',nlng
    write(iulog,*) 'photo: esfact = ',esfact
#endif
    col_loop : do i = 1,ncol
daylight : &
       if( zagtz(i) ) then
          sza    = zen_angle(i)*r2d
          secant = 1._r8 / cos( zen_angle(i) )
secant_in_bounds : &
          if( secant <= secant_limit ) then
             fracday(i) = 1._r8
             zarg(pverp:2:-1)     = zmid(i,:)
             zarg(1)              = 0._r8
             airdens(pverp:2:-1)  = tdens(i,:)
             airdens(1)           = 10._r8 * ps(i) / (boltz*ts(i))
             if( o3rad_ndx > 0 ) then
                spc_ndx = o3rad_ndx
             else
                spc_ndx = ox_ndx
             end if
             if( spc_ndx < 1 ) then
                spc_ndx = o3_ndx
             end if
             if( spc_ndx > 0 ) then
                o3line(pverp:2:-1) = vmr(i,:,spc_ndx)
             else
                o3line(pverp:2:-1) = 0._r8
             end if
             o3line(1)            = o3line(2)
             tline(pverp:2:-1)    = temper(i,:)
             tline(1)             = tline(2)
             rh(pverp:2:-1)       = relhum(i,:)
             rh(1)                = rh(2)
             xlwc(pverp:2:-1)     = cwat(i,:) * pmid(i,:)/(temper(i,:)*287._r8) * kg2g  !! TIE
             xlwc(1)              = xlwc(2)
             xfrc(pverp:2:-1)     = cldfr(i,:)                      ! cloud fraction
             xfrc(1)              = xfrc(2)
             tlay(1:pver)         = .5_r8*(tline(1:pver) + tline(2:pverp))
             albedo(1:nw)       = srf_alb(i)

             alat = rlats(i)
             along= rlons(i)

             if( oc1_ndx > 0 ) then
                aerocs1(pverp:2:-1) = vmr(i,:,oc1_ndx)
             else
                aerocs1(pverp:2:-1) = 0._r8
             end if
             aerocs1(1)            = aerocs1(2)
             if( oc2_ndx > 0 ) then
                aerocs2(pverp:2:-1) = vmr(i,:,oc2_ndx)
             else
                aerocs2(pverp:2:-1) = 0._r8
             end if
             aerocs2(1)          = aerocs2(2)
             if( cb1_ndx > 0 ) then
                aercbs1(pverp:2:-1) = vmr(i,:,cb1_ndx)
             else
                aercbs1(pverp:2:-1) = 0._r8
             end if
             aercbs1(1)          = aercbs1(2)
             if( cb2_ndx > 0 ) then
                aercbs2(pverp:2:-1) = vmr(i,:,cb2_ndx)
             else
                aercbs2(pverp:2:-1) = 0._r8
             end if
             aercbs2(1)          = aercbs2(2)
             if( soa_ndx > 0 ) then
                aersoa(pverp:2:-1) = vmr(i,:,soa_ndx)
             else
                aersoa(pverp:2:-1) = 0._r8
             end if
             aersoa(1)          = aersoa(2)
             if( ant_ndx > 0 ) then
                aerant(pverp:2:-1) = vmr(i,:,ant_ndx)
             else
                aerant(pverp:2:-1) = 0._r8
             end if
             aerant(1)            = aerant(2)
             if( so4_ndx > 0 ) then
                aerso4(pverp:2:-1) = vmr(i,:,so4_ndx)
             else
                aerso4(pverp:2:-1) = 0._r8
             end if
             aerso4(1)            = aerso4(2)
             if ( ndst == 4 ) then
                do ndx = 1,4
                   aerds(ndx,pverp:2:-1) = dust_vmr(i,:,ndx)
                end do
             else 
                do ndx = 1,4
                   aerds(ndx,pverp:2:-1) = 0._r8
                end do
             endif
             aerds(1,1)          = aerds(1,2)
             aerds(2,1)          = aerds(2,2)
             aerds(3,1)          = aerds(3,2)
             aerds(4,1)          = aerds(4,2)
             if( sa1_ndx > 0 ) then
                aersal(pverp:2:-1,1) = vmr(i,:,sa1_ndx)
             else
                aersal(pverp:2:-1,1) = 0._r8
             end if
             if( sa2_ndx > 0 ) then
                aersal(pverp:2:-1,2) = vmr(i,:,sa2_ndx)
             else
                aersal(pverp:2:-1,2) = 0._r8
             end if
             if( sa3_ndx > 0 ) then
                aersal(pverp:2:-1,3) = vmr(i,:,sa3_ndx)
             else
                aersal(pverp:2:-1,3) = 0._r8
             end if
             if( sa4_ndx > 0 ) then
                aersal(pverp:2:-1,4) = vmr(i,:,sa4_ndx)
             else
                aersal(pverp:2:-1,4) = 0._r8
             end if
             aersal(1,:) = aersal(2,:)
             call photoin( ncdate, alat, along, &
                           ut, esfact, col_dens(i,1,1), col_dens(i,1,2), albedo, &
                           zarg, tline, tlay, xlwc, xfrc, &
                           airdens, aerocs1, aerocs2, aercbs1, aercbs2, &
                           aersoa, aerant, aerso4, aersal, aerds, o3line, rh, &
                           prates, sza, nw, dt_xdiag )
             dt_diag(i,:) = dt_xdiag(:) 
             
             do m = 1,phtcnt
                if( lng_indexer(m) > 0 ) then
                   alias_factor = pht_alias_mult(m,2)
                   if( alias_factor == 1._r8 ) then
                      photos(i,:,m) = prates(1:pver,lng_indexer(m))
                   else
                      photos(i,:,m) = alias_factor * prates(1:pver,lng_indexer(m))
                   end if
                end if

#ifdef DEBUG
                if( do_diag ) then
                   write(iulog,'(''xactive_photo: prates('',i2,'',.)'')') m
                   write(iulog,'(1p,5e21.13)') photos(i,:pver,m)
                   write(iulog,*) ' '
                end if
#endif

             end do
!-----------------------------------------------------------------
!	... set j(onitr)
!-----------------------------------------------------------------
               if( jonitr_ndx > 0 ) then
                  if( jch3cho_a_ndx > 0 ) then
                     photos(i,1:pver,jonitr_ndx) = photos(i,1:pver,jch3cho_a_ndx)
                  end if
                  if( jch3cho_b_ndx > 0 ) then
                     photos(i,1:pver,jonitr_ndx) = photos(i,1:pver,jonitr_ndx) + photos(i,1:pver,jch3cho_b_ndx)
                  end if
                  if( jch3cho_c_ndx > 0 ) then
                     photos(i,1:pver,jonitr_ndx) = photos(i,1:pver,jonitr_ndx) + photos(i,1:pver,jch3cho_c_ndx)
                  end if
               end if
!-----------------------------------------------------------------
!	... calculate j(no) from formula
!-----------------------------------------------------------------
               if( jno_ndx > 0 ) then
                  if( has_o2_col .and. has_o3_col ) then
                     fac1(:) = 1.e-8_r8 * (col_dens(i,:,2)/cos(zen_angle(i)))**.38_r8
                     fac2(:) = 5.e-19_r8 * col_dens(i,:,1) / cos(zen_angle(i))
                     photos(i,:,jno_ndx) = 4.5e-6_r8 * exp( -(fac1(:) + fac2(:)) )
                  end if
               end if
!-----------------------------------------------------------------
! 	... add near IR correction to j(ho2no2)
!-----------------------------------------------------------------
               if( jho2no2_ndx > 0 ) then
                  photos(i,:,jho2no2_ndx) = photos(i,:,jho2no2_ndx) + 1.e-5_r8
               endif
          end if secant_in_bounds
       end if daylight
    end do col_loop

    call set_xnox_photo( photos, ncol  )

    deallocate( prates )

  end subroutine xactive_photo

  subroutine cloud_mod( zen_angle, clouds, lwc, delp, srf_alb, &
                        eff_alb, cld_mult )
    !-----------------------------------------------------------------------
    ! 	... cloud alteration factors for photorates and albedo
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    real(r8), intent(in)    ::  zen_angle         ! zenith angle
    real(r8), intent(in)    ::  srf_alb           ! surface albedo
    real(r8), intent(in)    ::  clouds(pver)       ! cloud fraction
    real(r8), intent(in)    ::  lwc(pver)          ! liquid water content (mass mr)
    real(r8), intent(in)    ::  delp(pver)         ! del press about midpoint in pascals
    real(r8), intent(out)   ::  eff_alb(pver)      ! effective albedo
    real(r8), intent(out)   ::  cld_mult(pver)     ! photolysis mult factor

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: k
    real(r8)    :: coschi
    real(r8)    :: del_lwp(pver)
    real(r8)    :: del_tau(pver)
    real(r8)    :: above_tau(pver)
    real(r8)    :: below_tau(pver)
    real(r8)    :: above_cld(pver)
    real(r8)    :: below_cld(pver)
    real(r8)    :: above_tra(pver)
    real(r8)    :: below_tra(pver)
    real(r8)    :: fac1(pver)
    real(r8)    :: fac2(pver)

    real(r8), parameter :: rgrav = 1._r8/9.80616_r8

    !---------------------------------------------------------
    !	... modify lwc for cloud fraction and form
    !	    liquid water path for each layer
    !---------------------------------------------------------
    where( clouds(:) /= 0._r8 )
       del_lwp(:) = rgrav * lwc(:) * delp(:) * 1.e3_r8 / clouds(:)
    elsewhere
       del_lwp(:) = 0._r8
    endwhere
    !---------------------------------------------------------
    !    	... form tau for each model layer
    !---------------------------------------------------------
    where( clouds(:) /= 0._r8 )
       del_tau(:) = del_lwp(:) *.155_r8 * clouds(:)**1.5_r8
    elsewhere
       del_tau(:) = 0._r8
    end where
    !---------------------------------------------------------
    !    	... form integrated tau from top down
    !---------------------------------------------------------
    above_tau(1) = 0._r8
    do k = 1,pverm
       above_tau(k+1) = del_tau(k) + above_tau(k)
    end do
    !---------------------------------------------------------
    !    	... form integrated tau from bottom up
    !---------------------------------------------------------
    below_tau(pver) = 0._r8
    do k = pverm,1,-1
       below_tau(k) = del_tau(k+1) + below_tau(k+1)
    end do
    !---------------------------------------------------------
    !	... form vertically averaged cloud cover above and below
    !---------------------------------------------------------
    above_cld(1) = 0._r8
    do k = 1,pverm
       above_cld(k+1) = clouds(k) * del_tau(k) + above_cld(k)
    end do
    do k = 2,pver
       if( above_tau(k) /= 0._r8 ) then
          above_cld(k) = above_cld(k) / above_tau(k)
       else
          above_cld(k) = above_cld(k-1)
       end if
    end do
    below_cld(pver) = 0._r8
    do k = pverm,1,-1
       below_cld(k) = clouds(k+1) * del_tau(k+1) + below_cld(k+1)
    end do
    do k = pverm,1,-1
       if( below_tau(k) /= 0._r8 ) then
          below_cld(k) = below_cld(k) / below_tau(k)
       else
          below_cld(k) = below_cld(k+1)
       end if
    end do
    !---------------------------------------------------------
    !	... modify above_tau and below_tau via jfm
    !---------------------------------------------------------
    where( above_cld(2:pver) /= 0._r8 )
       above_tau(2:pver) = above_tau(2:pver) / above_cld(2:pver)
    end where
    where( below_cld(:pverm) /= 0._r8 )
       below_tau(:pverm) = below_tau(:pverm) / below_cld(:pverm)
    end where
    where( above_tau(2:pver) < 5._r8 )
       above_cld(2:pver) = 0._r8
    end where
    where( below_tau(:pverm) < 5._r8 )
       below_cld(:pverm) = 0._r8
    end where
    !---------------------------------------------------------
    !	... form transmission factors
    !---------------------------------------------------------
    above_tra(:) = 11.905_r8 / (9.524_r8 + above_tau(:))
    below_tra(:) = 11.905_r8 / (9.524_r8 + below_tau(:))
    !---------------------------------------------------------
    !	... form effective albedo
    !---------------------------------------------------------
    where( below_cld(:) /= 0._r8 )
       eff_alb(:) = srf_alb + below_cld(:) * (1._r8 - below_tra(:)) &
                  * (1._r8 - srf_alb)
    elsewhere
       eff_alb(:) = srf_alb
    end where
    coschi = max( cos( zen_angle ),.5_r8 )
    where( del_lwp(:)*.155_r8 < 5._r8 )
       fac1(:) = 0._r8
    elsewhere
       fac1(:) = 1.4_r8 * coschi - 1._r8
    end where
    fac2(:)     = min( 0._r8,1.6_r8*coschi*above_tra(:) - 1._r8 )
    cld_mult(:) = 1._r8 + fac1(:) * clouds(:) + fac2(:) * above_cld(:)
    cld_mult(:) = max( .05_r8,cld_mult(:) )

  end subroutine cloud_mod

  subroutine set_ub_col( col_delta, vmr, invariants, ptop, pdel, ncol, lchnk )
    !---------------------------------------------------------------
    !        ... set the column densities at the upper boundary
    !---------------------------------------------------------------

    use chem_mods, only : nfs, ncol_abs=>nabscol, indexm
    use chem_mods, only : nabscol, gas_pcnst, indexm, nfs
    use chem_mods, only : gas_pcnst

    implicit none

    !---------------------------------------------------------------
    !        ... dummy args
    !---------------------------------------------------------------
    real(r8), intent(in)    ::  ptop(pcols)                            ! top pressure (Pa)
    integer,  intent(in)    ::  ncol                                   ! number of columns in current chunk
    integer,  intent(in)    ::  lchnk                                  ! latitude indicies in chunk
    real(r8), intent(in)    ::  vmr(ncol,pver,gas_pcnst)               ! xported species vmr
    real(r8), intent(in)    ::  pdel(pcols,pver)                       ! pressure delta about midpoints (Pa)
    real(r8), intent(in)    ::  invariants(ncol,pver,nfs)
    real(r8), intent(out)   ::  col_delta(ncol,0:pver,max(1,nabscol))  ! /cm**2 o2,o3 col dens above model

    !---------------------------------------------------------------
    !        ... local variables
    !---------------------------------------------------------------
    !---------------------------------------------------------------
    !        note: xfactor = 10.*r/(k*g) in cgs units.
    !              the factor 10. is to convert pdel
    !              from pascals to dyne/cm**2.
    !---------------------------------------------------------------
    real(r8), parameter :: xfactor = 2.8704e21_r8/(9.80616_r8*1.38044_r8)
    integer :: k, kl, spc_ndx
    integer :: ku(ncol)
    real(r8)    :: dp(ncol)
    real(r8)    :: tint_vals(2)
    real(r8)    :: o2_exo_col(ncol)
    real(r8)    :: o3_exo_col(ncol)
    integer :: lat, i
 
    !---------------------------------------------------------------
    !        ... assign column density at the upper boundary
    !            the first column is o3 and the second is o2.
    !            add 10 du o3 column above top of model.
    !---------------------------------------------------------------
    !---------------------------------------------------------------
    !	... set exo absorber columns
    !---------------------------------------------------------------
    has_abs_cols : if( has_o2_col .and. has_o3_col ) then
       if( has_fixed_press ) then
          kl = ki - 1
          if( has_o2_col ) then
             do i = 1,ncol
                if ( kl > 0 ) then
                   tint_vals(1) = o2_exo_coldens(kl,i,lchnk,last) &
                        + delp * (o2_exo_coldens(ki,i,lchnk,last) &
                        - o2_exo_coldens(kl,i,lchnk,last))
                   tint_vals(2) = o2_exo_coldens(kl,i,lchnk,next) &
                        + delp * (o2_exo_coldens(ki,i,lchnk,next) &
                        - o2_exo_coldens(kl,i,lchnk,next))
                else
                   tint_vals(1) = o2_exo_coldens( 1,i,lchnk,last) 
                   tint_vals(2) = o2_exo_coldens( 1,i,lchnk,next) 
                endif
                o2_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
             end do
          else
             o2_exo_col(:) = 0._r8
          end if
          if( has_o3_col ) then
             do i = 1,ncol
                if ( kl > 0 ) then
                   tint_vals(1) = o3_exo_coldens(kl,i,lchnk,last) &
                        + delp * (o3_exo_coldens(ki,i,lchnk,last) &
                        - o3_exo_coldens(kl,i,lchnk,last))
                   tint_vals(2) = o3_exo_coldens(kl,i,lchnk,next) &
                        + delp * (o3_exo_coldens(ki,i,lchnk,next) &
                        - o3_exo_coldens(kl,i,lchnk,next))
                else
                   tint_vals(1) = o3_exo_coldens( 1,i,lchnk,last) 
                   tint_vals(2) = o3_exo_coldens( 1,i,lchnk,next) 
                endif
                o3_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
             end do
          else
             o3_exo_col(:) = 0._r8
          end if
#ifdef DEBUG
          write(iulog,*) '-----------------------------------'
          write(iulog,*) 'set_ub_col: diagnostics @ lat = ',lat
          write(iulog,*) 'o2_exo_col'
          write(iulog,'(1p,5g15.7)') o2_exo_col(:)
          write(iulog,*) 'o3_exo_col'
          write(iulog,'(1p,5g15.7)') o3_exo_col(:)
          write(iulog,*) '-----------------------------------'
#endif
       else
          !---------------------------------------------------------------
          !	... do pressure interpolation
          !---------------------------------------------------------------
          call p_interp( lchnk, ncol, ptop, o2_exo_col, o3_exo_col )
       end if
    else
       o2_exo_col(:) = 0._r8
       o3_exo_col(:) = 0._r8
    end if has_abs_cols

    if( o3rad_ndx > 0 ) then
       spc_ndx = o3rad_ndx
    else
       spc_ndx = ox_ndx
    end if
    if( spc_ndx < 1 ) then
       spc_ndx = o3_ndx
    end if
    if( spc_ndx > 0 ) then
       col_delta(:,0,1) = o3_exo_col(:)
       do k = 1,pver
          col_delta(:ncol,k,1) = xfactor * pdel(:ncol,k) * vmr(:ncol,k,spc_ndx)
       end do
    else if( o3_inv_ndx > 0 ) then
       col_delta(:,0,1) = o3_exo_col(:)
       do k = 1,pver
          col_delta(:ncol,k,1) = xfactor * pdel(:ncol,k) * invariants(:ncol,k,o3_inv_ndx)/invariants(:ncol,k,indexm)
       end do
    else
       col_delta(:,:,1) = 0._r8
    end if
    if( ncol_abs > 1 ) then
       if( o2_ndx > 1 ) then
          col_delta(:,0,2) = o2_exo_col(:)
          if( o2_is_inv ) then
             do k = 1,pver
                col_delta(:ncol,k,2) = xfactor * pdel(:ncol,k) * invariants(:ncol,k,o2_ndx)/invariants(:ncol,k,indexm)
             end do
          else
             do k = 1,pver
                col_delta(:ncol,k,2) = xfactor * pdel(:ncol,k) * vmr(:ncol,k,o2_ndx)
             end do
          endif
       else
          col_delta(:,:,2) = 0._r8
       end if
    end if

  end subroutine set_ub_col

  subroutine p_interp( lchnk, ncol, ptop, o2_exo_col, o3_exo_col )
    !---------------------------------------------------------------
    !     	... pressure interpolation for exo col density
    !---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !     	... dummy arguments
    !---------------------------------------------------------------
    integer, intent(in)     :: ncol                      ! no. of columns
    real(r8), intent(out)   :: o2_exo_col(ncol)          ! exo model o2 column density (molecules/cm^2)
    real(r8), intent(out)   :: o3_exo_col(ncol)          ! exo model o3 column density (molecules/cm^2)
    integer, intent(in)     :: lchnk                     ! latitude  indicies in chunk
    real(r8)                :: ptop(pcols)               ! top pressure (Pa)

    !---------------------------------------------------------------
    !     	... local variables
    !---------------------------------------------------------------
    integer  :: i, k, ki, kl
    integer  :: ku(ncol)                               ! interpolation index
    real(r8) :: pinterp
    real(r8) :: delp
    real(r8) :: tint_vals(2)
    real(r8) :: dp(ncol)                               ! pressure interpolation factor

    do i = 1,ncol
       pinterp = ptop(i)
       if( pinterp < levs(1) ) then
          ki   = 0
          delp = 0._r8
       else
          do ki = 2,n_exo_levs
             if( pinterp <= levs(ki) ) then
                delp = log( pinterp/levs(ki-1) )/log( levs(ki)/levs(ki-1) )
                exit
             end if
          end do
       end if
       kl = ki - 1
       if( has_o2_col ) then
          tint_vals(1) = o2_exo_coldens(kl,i,lchnk,last) &
                         + delp * (o2_exo_coldens(ki,i,lchnk,last) &
                                 - o2_exo_coldens(kl,i,lchnk,last))
          tint_vals(2) = o2_exo_coldens(kl,i,lchnk,next) &
                         + delp * (o2_exo_coldens(ki,i,lchnk,next) &
                         - o2_exo_coldens(kl,i,lchnk,next))
          o2_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
       else
          o2_exo_col(i) = 0._r8
       end if
       if( has_o3_col ) then
          tint_vals(1) = o3_exo_coldens(kl,i,lchnk,last) &
                         + delp * (o3_exo_coldens(ki,i,lchnk,last) &
                         - o3_exo_coldens(kl,i,lchnk,last))
          tint_vals(2) = o3_exo_coldens(kl,i,lchnk,next) &
                         + delp * (o3_exo_coldens(ki,i,lchnk,next) &
                         - o3_exo_coldens(kl,i,lchnk,next))
          o3_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
       else
          o3_exo_col(i) = 0._r8
       end if
    end do

  end subroutine p_interp

  subroutine setcol( col_delta, col_dens, vmr, pdel,  ncol )
    !---------------------------------------------------------------
    !     	... set the column densities
    !---------------------------------------------------------------

    use chem_mods, only : ncol_abs=>nabscol, gas_pcnst

    implicit none

    !---------------------------------------------------------------
    !     	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)    :: ncol                              ! no. of columns in current chunk
    real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst)          ! xported species vmr
    real(r8), intent(in)    :: pdel(pcols,pver)                  ! delta about midpoints
    real(r8), intent(in)    :: col_delta(:,0:,:)                 ! layer column densities (molecules/cm^2)
    real(r8), intent(out)   :: col_dens(:,:,:)                   ! column densities ( /cm**2 )

    !---------------------------------------------------------------
    !        the local variables
    !---------------------------------------------------------------
    integer  ::   i, k, km1, m      ! long, alt indicies

    !---------------------------------------------------------------
    !        note: xfactor = 10.*r/(k*g) in cgs units.
    !              the factor 10. is to convert pdel
    !              from pascals to dyne/cm**2.
    !---------------------------------------------------------------
    real(r8), parameter :: xfactor = 2.8704e21_r8/(9.80616_r8*1.38044_r8)

    !---------------------------------------------------------------
    !   	... compute column densities down to the
    !           current eta index in the calling routine.
    !           the first column is o3 and the second is o2.
    !---------------------------------------------------------------
    do m = 1,ncol_abs
       col_dens(:,1,m) = col_delta(:,0,m) + .5_r8 * col_delta(:,1,m)
       do k = 2,pver
          km1 = k - 1
          col_dens(:,k,m) = col_dens(:,km1,m) + .5_r8 * (col_delta(:,km1,m) + col_delta(:,k,m))
       end do
    enddo

  end subroutine setcol

  subroutine photo_timestep_init( calday )
    use time_manager,   only : is_end_curr_day
    use euvac,          only : euvac_set_etf
    use mo_solar_parms, only : get_solar_parms
    use mo_jshort,      only : jshort_timestep_init
    use mo_jlong,       only : jlong_timestep_init

    !-----------------------------------------------------------------------------
    !	... setup the time interpolation
    !-----------------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------------
    !	... dummy arguments
    !-----------------------------------------------------------------------------
    real(r8), intent(in)    ::  calday                                   ! day of year at end of present time step

    !-----------------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------------
    integer :: m
    real(r8) :: f107
    real(r8) :: f107a

    if (waccm) then
       if( is_end_curr_day() ) then
          call get_solar_parms( f107_s = f107, f107a_s = f107a )
          call euvac_set_etf( f107, f107a )
       end if
    endif

    if( has_o2_col .or. has_o3_col ) then
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
#ifdef DEBUG
       write(iulog,*) '-----------------------------------'
       write(iulog,*) 'photo_timestep_init: diagnostics'
       write(iulog,*) 'calday, last, next, dels = ',calday,last,next,dels
       write(iulog,*) '-----------------------------------'
#endif
    end if

    !-----------------------------------------------------------------------
    ! Set jlong etf
    !-----------------------------------------------------------------------
    call jlong_timestep_init
    
    if ( do_jshort ) then
       !-----------------------------------------------------------------------
       ! Set jshort etf
       !-----------------------------------------------------------------------
       call jshort_timestep_init
    endif

  end subroutine photo_timestep_init

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  subroutine set_xnox_photo( photos, ncol )
    use chem_mods,    only : ncol_abs => nabscol, phtcnt, pcnstm1 => gas_pcnst, nfs
    implicit none
    integer, intent(in)     :: ncol
    real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)

    if( jno2a_ndx > 0 .and. jno2_ndx > 0 ) then
       photos(:,:,jno2a_ndx) = photos(:,:,jno2_ndx)
    end if
    if( jn2o5a_ndx > 0 .and. jn2o5_ndx > 0 ) then
       photos(:,:,jn2o5a_ndx) = photos(:,:,jn2o5_ndx)
    end if
    if( jn2o5b_ndx > 0 .and. jn2o5_ndx > 0 ) then
       photos(:,:,jn2o5b_ndx) = photos(:,:,jn2o5_ndx)
    end if
    if( jhno3a_ndx > 0 .and. jhno3_ndx > 0 ) then
       photos(:,:,jhno3a_ndx) = photos(:,:,jhno3_ndx)
    end if

    if( jno3a_ndx > 0 .and. jno3_ndx > 0 ) then
       photos(:,:,jno3a_ndx) = photos(:,:,jno3_ndx)
    end if
    if( jho2no2a_ndx > 0 .and. jho2no2_ndx > 0 ) then
       photos(:,:,jho2no2a_ndx) = photos(:,:,jho2no2_ndx)
    end if
    if( jmpana_ndx > 0 .and. jmpan_ndx > 0 ) then
       photos(:,:,jmpana_ndx) = photos(:,:,jmpan_ndx)
    end if
    if( jpana_ndx > 0 .and. jpan_ndx > 0 ) then
       photos(:,:,jpana_ndx) = photos(:,:,jpan_ndx)
    end if
    if( jonitra_ndx > 0 .and. jonitr_ndx > 0 ) then
       photos(:,:,jonitra_ndx) = photos(:,:,jonitr_ndx)
    end if
    if( jo1da_ndx > 0 .and. jo1d_ndx > 0 ) then
       photos(:,:,jo1da_ndx) = photos(:,:,jo1d_ndx)
    end if
    if( jo3pa_ndx > 0 .and. jo3p_ndx > 0 ) then
       photos(:,:,jo3pa_ndx) = photos(:,:,jo3p_ndx)
    end if

  endsubroutine set_xnox_photo

end module mo_photo
