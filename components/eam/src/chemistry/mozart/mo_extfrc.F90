module mo_extfrc
  !---------------------------------------------------------------
  ! 	... insitu forcing module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, begchunk, endchunk, pver, pverp
  use chem_mods,    only : gas_pcnst, extcnt
  use spmd_utils,   only : masterproc,iam
  use cam_abortutils,   only : endrun
  use cam_history,  only : addfld, horiz_only, outfld, add_default
  use cam_logfile,  only : iulog
  use tracer_data,  only : trfld,trfile
  use phys_grid,    only : get_rlat_all_p, get_rlon_all_p,get_area_all_p
  use time_manager,  only: get_curr_date
  use mo_constants, only : pi, rgrav, rearth, avogadro
  use physconst,         only: rhoh2o, rga, rair
  implicit none

  type :: forcing
     integer           :: frc_ndx
     real(r8)              :: mw
     character(len=265) :: filename
     real(r8), pointer     :: times(:)
     real(r8), pointer     :: levi(:)
     character(len=8)  :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type forcing

  private
  public  :: extfrc_inti
  public  :: extfrc_set
  public  :: extfrc_timestep_init

  save

  integer, parameter :: time_span = 1

  character(len=256) ::   filename

  logical :: has_extfrc(gas_pcnst)
  type(forcing), allocatable  :: forcings(:)
  integer :: extfrc_cnt = 0
  integer, parameter :: nfire = 10 !!  fire emission
  integer :: nfire_count, nba_count
  integer :: PH_emis_m(nfire), PH_emis_n(nfire) ! fire emission indices
  integer :: BA_emis_m(nfire), BA_emis_n(nfire) ! BA indices
  logical :: plumerise = .false.
  logical :: emis_constrained_frp = .false.
  logical :: diag_run_plumerise = .false.
  real(r8) :: ef_bc_a4 = 0.55_r8*1.0e-03_r8*(1.0_r8/0.50_r8)
  real(r8) :: ef_oc_a4 = 10.9_r8*1.0e-03_r8*(1.0_r8/0.50_r8)
  real(r8) :: ef_h2o_a4 = 350.0_r8*1.0e-03_r8*(1.0_r8/0.50_r8)*(12.0_r8/18.0_r8) ! adjust moleculer weight from carbon to h2o
  !real(r8) :: ef_h2o_a4 = (1.0_r8/0.45_r8) ! test the sensitivity
  real(r8) :: fix_plume_height = 0.0_r8
  real(r8) :: detrainment_para = 0.0_r8
  real(r8) :: ba_ratio = 0.0_r8
contains

  subroutine extfrc_inti( extfrc_specifier, extfrc_type, extfrc_cycle_yr, extfrc_fixed_ymd, extfrc_fixed_tod)

    !-----------------------------------------------------------------------
    ! 	... initialize the surface forcings
    !-----------------------------------------------------------------------
    use cam_pio_utils, only : cam_pio_openfile
    use pio, only : pio_inq_dimid, pio_inquire, pio_inq_varndims, pio_closefile, &
         pio_inq_varname, pio_nowrite, file_desc_t
    use pio,              only : pio_inq_vardimid !(zhang73)
    use mo_tracname,   only : solsym
    use mo_chem_utls,  only : get_extfrc_ndx, get_spc_ndx
    use chem_mods,     only : frc_from_dataset
    use tracer_data,   only : trcdata_init
    use phys_control,  only : phys_getopts
    use physics_buffer, only : physics_buffer_desc

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), dimension(:), intent(in) :: extfrc_specifier
    character(len=*), intent(in) :: extfrc_type
    integer  , intent(in)        :: extfrc_cycle_yr
    integer  , intent(in)        :: extfrc_fixed_ymd
    integer  , intent(in)        :: extfrc_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: j, l, m, n, i,mm                          ! Indices
    character(len=16)  :: species
    character(len=16)  :: spc_name
    character(len=256) :: locfn
    character(len=256) :: spc_fnames(gas_pcnst)

    integer ::  vid, ndims, nvars, isec, ierr
    integer :: dimids(8), did, dimid,ncol_dimid,lat_dimid,time_dimid !(zhang73)
    integer, allocatable :: finddim_time(:), finddim_lat_ncol(:) !(zhang73) 
    type(file_desc_t) :: ncid
    character(len=32)  :: varname

    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.
    logical  :: history_aerosol      ! Output the MAM aerosol tendencies
    logical  :: history_verbose      ! produce verbose history output

    !-----------------------------------------------------------------------
 
    call phys_getopts( history_aerosol_out        = history_aerosol, &
                       history_verbose_out        = history_verbose, &
                       emis_constrained_frp_out   = emis_constrained_frp, &
                       diag_run_plumerise_out   = diag_run_plumerise, &
                       plumerise_out              = plumerise, &
                       fix_plume_height_out = fix_plume_height, &
                       detrainment_para_out = detrainment_para   )

    do i = 1, gas_pcnst
       has_extfrc(i) = .false.
       spc_fnames(i) = ''
    enddo

    !-----------------------------------------------------------------------
    ! 	... species has insitu forcing ?
    !-----------------------------------------------------------------------

    !write(iulog,*) 'Species with insitu forcings'

    count_emis: do n=1,gas_pcnst

       if ( len_trim(extfrc_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(extfrc_specifier(n),'->')
       spc_name = trim(adjustl(extfrc_specifier(n)(:i-1)))
       filename = trim(adjustl(extfrc_specifier(n)(i+2:)))

       m = get_extfrc_ndx( spc_name )

       if ( m < 1 ) then
          call endrun('extfrc_inti: '//trim(spc_name)// ' does not have an external source')
       endif

       if ( .not. frc_from_dataset(m) ) then
          call endrun('extfrc_inti: '//trim(spc_name)//' cannot have external forcing from additional dataset')
       endif

       mm = get_spc_ndx(spc_name)
       spc_fnames(mm) = filename

       has_extfrc(mm) = .true.
       !write(iulog,*) '   ',  spc_name ,' : filename = ',trim(spc_fnames(mm)),' spc ndx = ',mm

    enddo count_emis

    extfrc_cnt = count( has_extfrc(:) )

    if( extfrc_cnt < 1 ) then
       if (masterproc) write(iulog,*) 'There are no species with insitu forcings'
       return
    end if

    if (masterproc) write(iulog,*) ' '

    !-----------------------------------------------------------------------
    ! 	... allocate forcings type array
    !-----------------------------------------------------------------------
    allocate( forcings(extfrc_cnt), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'extfrc_inti: failed to allocate forcings array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the forcing type array
    !-----------------------------------------------------------------------
    n = 0
    species_loop : do m = 1,gas_pcnst
       has_forcing : if( has_extfrc(m) ) then
          spc_name = trim( solsym(m) )
          n        = n + 1
          !-----------------------------------------------------------------------
          ! 	... default settings
          !-----------------------------------------------------------------------
          forcings(n)%frc_ndx          = get_extfrc_ndx( spc_name )
          forcings(n)%species          = spc_name
          forcings(n)%filename         = spc_fnames(m)
          call addfld( trim(spc_name)//'_XFRC', (/ 'lev' /), 'A',  'molec/cm3/s', &
                       'external forcing for '//trim(spc_name) )
          call addfld( trim(spc_name)//'_CLXF', horiz_only, 'A',  'molec/cm2/s', &
                       'vertically intergrated external forcing for '//trim(spc_name) )
          if ( history_aerosol ) then 
             if (history_verbose) call add_default( trim(spc_name)//'_XFRC', 1, ' ' )
             call add_default( trim(spc_name)//'_CLXF', 1, ' ' )
          endif
       end if has_forcing
    end do species_loop
    !plume-rise diagnose
    if (plumerise) then
       !call addfld( 'bc_a4_EM_XFRC', (/ 'lev' /), 'A',  'molec/cm3/s', &
       !                'external forcing for UCI emis before plumerise' )
       !call addfld( 'bc_a4_EM_PH_XFRC', (/ 'lev' /), 'A',  'molec/cm3/s', &
       !                'external forcing for UCI emis after plumerise' )
       call addfld( 'plume_height_EM', horiz_only, 'I',  'meter', &
                       'plumerise height caused by EM fires' ) 
       call addfld( 'heat_flux_plume', horiz_only, 'I',  'kw/m2', &
                       'heat flux used in plume-rise model' )
       call addfld( 'burned_area', horiz_only, 'I',  'm2', &
                       'burned_area used in plume-rise model' )
       call addfld( 'zmidr_ph', (/ 'lev' /), 'I',  'km', &
                       'midpoint geopotential in km realitive to surf' )
       call addfld( 'pmid_ph', (/ 'lev' /), 'I',  'Pa', &
                       'midpoint pressures (Pa)' )
       call addfld( 'tfld_ph', (/ 'lev' /), 'I',  'K', &
                       'midpoint temperature (K)' )
       call addfld( 'relhum_ph', (/ 'lev' /), 'I',  'unitless', &
                       'relative humidity' )
       call addfld( 'qh2o_ph', (/ 'lev' /), 'I',  'kg/kg', &
                       'specific humidity' ) 
       call addfld( 'ufld_ph', (/ 'lev' /), 'I',  'm/s', &
                       'zonal velocity (m/s)' )
       call addfld( 'vfld_ph', (/ 'lev' /), 'I',  'm/s', &
                       'meridional velocity (m/s)' )
       
       !wt_ini_e3sm_out, wt_end_e3sm_out, rbuoy_ini_e3sm_out, rbuoy_end_e3sm_out
       call addfld( 'wt_ini_e3sm_out', (/ 'lev' /), 'I',  'm/s', &
                       'verticle velocity (m/s)' )
       call addfld( 'wt_end_e3sm_out', (/ 'lev' /), 'I',  'm/s', &
                       'verticle velocity (m/s)' )
       call addfld( 'rbuoy_ini_e3sm_out', (/ 'lev' /), 'I',  'm/s', &
                       'verticle velocity (m/s)' )
       call addfld( 'rbuoy_end_e3sm_out', (/ 'lev' /), 'I',  'm/s', &
                       'verticle velocity (m/s)' )      
       call addfld( 't_ini_e3sm_out', (/ 'lev' /), 'I',  'temperature', &
                       'verticle velocity (m/s)' )
       call addfld( 't_end_e3sm_out', (/ 'lev' /), 'I',  'temperature', &
                       'verticle velocity (m/s)' )
       call addfld( 'qv_ini_e3sm_out', (/ 'lev' /), 'I',  'kg/kg', &
                       'verticle velocity (m/s)' )
       call addfld( 'qv_end_e3sm_out', (/ 'lev' /), 'I',  'kg/kg', &
                       'verticle velocity (m/s)' )
       call addfld( 'r_ini_e3sm_out', (/ 'lev' /), 'I',  'm', &
                       'verticle velocity (m/s)' )
       call addfld( 'r_end_e3sm_out', (/ 'lev' /), 'I',  'm', &
                       'verticle velocity (m/s)' )
       call addfld( 'rho_end_e3sm_out', (/ 'lev' /), 'I',  'm', &
                       'verticle velocity (m/s)' )
    endif
    !---------------------------------------------------------------------
    if (masterproc) then
       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'extfrc_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'extfrc timing specs'
       write(iulog,*) 'type = ',extfrc_type
       if( extfrc_type == 'FIXED' ) then
          write(iulog,*) ' fixed date = ', extfrc_fixed_ymd
          write(iulog,*) ' fixed time = ', extfrc_fixed_tod
       else if( extfrc_type == 'CYCLICAL' ) then
          write(iulog,*) ' cycle year = ',extfrc_cycle_yr
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',extfrc_cnt,' species with external forcing files'
       do m = 1,extfrc_cnt
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',m
          write(iulog,*) 'species = ',trim(forcings(m)%species)
          write(iulog,*) 'frc ndx = ',forcings(m)%frc_ndx
          write(iulog,*) 'filename= ',trim(forcings(m)%filename)
       end do
       write(iulog,*) ' '
    endif

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    PH_emis_m(:) = -1 ! fire emission type
    PH_emis_n(:) = -1 ! fire emission sector 
    BA_emis_m(:) = -1 ! fire BA type
    BA_emis_n(:) = -1 ! fire BA sector
    nfire_count = 0 ! fire emission counted
    nba_count   = 0 ! fire area count
    frcing_loop: do m = 1, extfrc_cnt

       forcings(m)%nsectors = 0

       call cam_pio_openfile ( ncid, trim(forcings(m)%filename), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nVariables=nvars)

       allocate(finddim_time(nvars))
       allocate(finddim_lat_ncol(nvars))
       finddim_time=0
       finddim_lat_ncol=0
       time_dimid=-9999
       lat_dimid=-9999
       ncol_dimid=-9999

       ierr = pio_inq_dimid(ncid, 'time', dimid)
       if(ierr==0) time_dimid = dimid
       ierr = pio_inq_dimid(ncid, 'lat', dimid)
       if(ierr==0) lat_dimid = dimid
       ierr = pio_inq_dimid(ncid, 'ncol', dimid)
       forcings(m)%file%is_ncol = (ierr==0)
       if(ierr==0) ncol_dimid = dimid
       if(masterproc) write(iulog,*) '(zhang73 extfrc_inti) time_dimid, lat_dimid, ncol_dimid=',time_dimid, lat_dimid, ncol_dimid

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, ndims)

          ierr = pio_inq_vardimid (ncid, vid, dimids(1:ndims)) !(zhang73)
          do did=1,ndims
             if( dimids(did) == time_dimid ) finddim_time(vid)=1
             if(  dimids(did) == lat_dimid ) finddim_lat_ncol(vid)=1
             if( dimids(did) == ncol_dimid ) finddim_lat_ncol(vid)=1
          enddo

          ierr = pio_inq_varname (ncid, vid, varname)
          if( finddim_time(vid)==1 .and. finddim_lat_ncol(vid)==1)then !(zhang73)
             !write(iulog,*) '(zhang73 extfrc_inti) valid var: finddim_time(vid), finddim_lat_ncol(vid)=',trim(varname),finddim_time(vid), finddim_lat_ncol(vid)
             forcings(m)%nsectors = forcings(m)%nsectors+1
             ! kzm note: here assumes fire emission are in ncol emission files 
             if (plumerise)then
                !if (trim(forcings(m)%species) == 'bc_a4' .and. trim(varname) == 'EM')then
                if (trim(varname) == 'EM' &
                    .or. trim(varname) == 'num_a1_BC_ELEV_EM' & 
                    .or. trim(varname) == 'num_a1_POM_ELEV_EM')then
                   nfire_count = nfire_count +1
                   if(masterproc) write(iulog,*) forcings(m)%species
                   if(masterproc) write(iulog,*) 'sector number = ', forcings(m)%nsectors
                   if(masterproc) write(iulog,*) 'UCI wildfire emission in model type ', nfire_count
                   PH_emis_m(nfire_count) = m
                   PH_emis_n(nfire_count) = forcings(m)%nsectors
                   if(masterproc) write(iulog,*) 'PH_emis_m', PH_emis_m(nfire_count)
                   if(masterproc) write(iulog,*) 'PH_emis_n', PH_emis_n(nfire_count)

                elseif (trim(varname) == 'BA')then
                   nba_count = nba_count +1
                   if(masterproc) write(iulog,*) forcings(m)%species
                   if(masterproc) write(iulog,*) 'sector number = ', forcings(m)%nsectors
                   if(masterproc) write(iulog,*) 'UCI wildfire BA in model type ', nba_count
                   BA_emis_m(nba_count) = m
                   BA_emis_n(nba_count) = forcings(m)%nsectors
                   if(masterproc) write(iulog,*) 'BA_emis_m', BA_emis_m(nba_count)
                   if(masterproc) write(iulog,*) 'BA_emis_n', BA_emis_n(nba_count)
                endif
             endif
          else
             write(iulog,*) 'extfrc_inti: Skipping variable ', trim(varname),', ndims = ',ndims,' , species=',trim(forcings(m)%species)
          end if
       enddo

       allocate( forcings(m)%sectors(forcings(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'extfrc_inti: failed to allocate forcings(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1
       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, ndims)
!          if( ndims == dim_thres ) then !(zhang73) check ndims from 4 -> 3 to activate bc_a4_ncol
          if( finddim_time(vid)==1 .and. finddim_lat_ncol(vid)==1)then !(zhang73)
             ierr = pio_inq_varname(ncid, vid, forcings(m)%sectors(isec))
             isec = isec+1
          endif

       enddo
       deallocate(finddim_time)
       deallocate(finddim_lat_ncol)

       call pio_closefile (ncid)

       allocate(forcings(m)%file%in_pbuf(size(forcings(m)%sectors)))
       forcings(m)%file%in_pbuf(:) = .false.
       call trcdata_init( forcings(m)%sectors, &
                          forcings(m)%filename, filelist, datapath, &
                          forcings(m)%fields,  &
                          forcings(m)%file, &
                          rmv_file, extfrc_cycle_yr, extfrc_fixed_ymd, extfrc_fixed_tod, extfrc_type)

    enddo frcing_loop
    if(masterproc) write(iulog,*) 'PH_emis_m', PH_emis_m(:)
    if(masterproc) write(iulog,*) 'PH_emis_n', PH_emis_n(:)

  end subroutine extfrc_inti

  subroutine extfrc_timestep_init( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,extfrc_cnt
       call advance_trcdata( forcings(m)%fields, forcings(m)%file, state, pbuf2d  )
    end do

  end subroutine extfrc_timestep_init

 ! subroutine extfrc_set( lchnk, zint, frcing, ncol )
  subroutine extfrc_set( lchnk, zint, frcing, ncol, &
                         zmidr, pmid, tfld, relhum, qh2o, ufld, vfld )
   
    !--------------------------------------------------------
    !	... form the external forcing
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    :: ncol                  ! columns in chunk
    integer,  intent(in)    :: lchnk                 ! chunk index
    real(r8), intent(in)    :: zint(ncol, pverp)                  ! interface geopot above surface (km)
    real(r8), intent(inout) :: frcing(ncol,pver,extcnt)   ! insitu forcings (molec/cm^3/s)
    ! plume-rise parameters
    real(r8), intent(in)  ::   zmidr(ncol,pver)             ! midpoint geopot height - elevation ( km )
    real(r8), intent(in)  ::   pmid(pcols,pver)            ! midpoint pressure (Pa)
    real(r8), intent(in)  ::   tfld(pcols,pver)            ! midpoint temperature (K)
    real(r8), intent(in)  ::   relhum(ncol,pver)           ! relative humidity (0~1)
    real(r8), intent(in)  ::   qh2o(pcols,pver)            ! specific humidity (kg/kg)
    real(r8), intent(in)  ::   ufld(pcols,pver)            ! zonal velocity (m/s)
    real(r8), intent(in)  ::   vfld(pcols,pver)            ! meridional velocity (m/s)
    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n
    character(len=16) :: xfcname
    real(r8) :: frcing_col(1:ncol)
    integer  :: k, isec
    real(r8),parameter :: km_to_cm = 1.e5_r8
    integer :: icol ! plumerise
    !------------------------------------------------------
    !    ... plume_height variables
    !-----------------------------------------------------
    real(r8) :: plume_height,emis_col
    real(r8) :: plume_height_EM(ncol), pt_v(pver), heat_flux_plume(ncol)
    integer :: ph_z(ncol)   ! index of levels of max plume height
    real(r8) :: frcing_col_plume,frcing_vertical_plume_old(pver),frcing_vertical_plume_new(pver)
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)
    real(r8) :: burnedarea, heatflux,coef_bc,coef_pom
    real(r8) :: burnedarea_plume(pcols)
    real(r8) :: tl, frp, frp_memory(pcols), frp4plume,burnedarea_memory(pcols),heatflux_memory(pcols)! 
    real(r8) :: area(ncol),frcing_col_plume_bc_memory(pcols),frcing_col_plume_pom_memory(pcols)
    integer :: iyear,imo,iday_m,tod,tod_saved
    real(r8) :: wt_ini_e3sm_out(pcols,pver), wt_end_e3sm_out(pcols,pver), &
                rbuoy_ini_e3sm_out(pcols,pver), rbuoy_end_e3sm_out(pcols,pver), &
                t_ini_e3sm_out(pcols,pver), t_end_e3sm_out(pcols,pver),&
                qv_ini_e3sm_out(pcols,pver), qv_end_e3sm_out(pcols,pver),&
                r_ini_e3sm_out(pcols,pver), r_end_e3sm_out(pcols,pver), &
                rho_ini_e3sm_out(pcols,pver), rho_end_e3sm_out(pcols,pver)
    real(r8) :: air_density(pcols,pver) ! (kg/m3)
    real(r8) :: qmass(pver),qratio(pver),qmass_1
    real(r8) :: v_air, m_vapor, f_vapor,wdqdz,qdwdz,phase_out
    logical :: surface_flux_flag = .false.
    !real(r8) :: detrainment_para = 0.5_r8 ! tunable parameter for detrainment fraction from surface plume
    !kzm ++ surface air flux
    !mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
    if (detrainment_para > 0.0_r8) then
       surface_flux_flag = .true.
    else
       surface_flux_flag = .false.
    endif
    if (surface_flux_flag) then
       air_density(:ncol,:) = pmid(:ncol,:)/(rair*tfld(:ncol,:))
    endif
    !kzm --
    if( extfrc_cnt < 1 .or. extcnt < 1 ) then
       return
    end if
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call get_curr_date (iyear,imo,iday_m,tod)
    !--------------------------------------------------------
    !	... set non-zero forcings
    !--------------------------------------------------------
    frp_memory(:) = 0.0_r8
    heatflux_memory(:) = 0.0_r8
    frcing_col_plume_bc_memory(:) = 0.0_r8
    burnedarea_memory(:) = 0.0_r8
    plume_height_EM(:) = 0._r8
    heat_flux_plume(:) = 0._r8
    burnedarea_plume(:) = 0.0_r8
    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2
    coef_bc = 1.7847635e-19
    coef_pom = 1.04986099e-19

    ! loop for calculate fire area, heat flux, emission
    src_loop1 : do m = 1,extfrc_cnt

      n = forcings(m)%frc_ndx

       frcing(:ncol,:,n) = 0._r8
       do isec = 1,forcings(m)%nsectors
          !check wildfire emission EM
          if ((plumerise) .and. (forcings(m)%file%alt_data)) then
             !plume-rise calculation for certain species
             if ( forcings(m)%species == 'bc_a4' .and. forcings(m)%sectors(isec) == 'EM'    )then
                frcing_vertical_plume_old(:pver) = 0._r8
                frcing_vertical_plume_new(:pver) = 0._r8
                ! loop for each col
                do icol=1,pcols ! calculate if EM emitted
                   frcing_col_plume = 0._r8 !set initial value for each column
                   emis_col = sum(forcings(m)%fields(isec)%data(icol,:,lchnk))
                   if ( emis_col  > 0.0_r8) then
                         ! get initial emission from forcing data
                         frcing_vertical_plume_old(1:pver) = forcings(m)%fields(isec)%data(icol,pver:1:-1,lchnk) ! reverse
                         ! calculate the total emission in this column
                         do k = 1,pver
                            ! unit is molecular/m2
                            frcing_col_plume = frcing_col_plume +  &
                                                     frcing_vertical_plume_old(k)*(zint(icol,k)-zint(icol,k+1))*km_to_cm
                           ! write(iulog,*)'kzm_level ',k, 'old emis ', frcing_vertical_plume_old(k)
                         enddo
                         ! heatflux related to total carbon burned
                         frcing_col_plume_bc_memory(icol) = frcing_col_plume
                         heatflux = (frcing_col_plume*1.0e4)/6.022e23_r8*12.0_r8 ! carbon g/m2/s
                         heatflux = heatflux/0.45_r8 ! carbon to fuel, Rowell et al., 2013; unit g/m2/s
                         heatflux = heatflux*1.0E-3_r8 ! convert to kg/m2/s
                         heatflux = heatflux*19.0e6_r8 ! forest heat from burning (J/s/m2), Freitas et al., 2011
                         heatflux = heatflux*1.0e-3_r8 ! KW/m2
                         heatflux_memory(icol) = heatflux
                         if (diag_run_plumerise) then
                            write(iulog,*)'kzm_heat_flux ', heatflux, icol
                         endif
                         tod_saved = tod
                   endif
                 enddo
              elseif ( trim(forcings(m)%species) == 'bc_a4' .and.  &
                       trim(forcings(m)%sectors(isec)) == 'BA'     )then
                 frcing_vertical_plume_old(:pver) = 0._r8
                 frcing_vertical_plume_new(:pver) = 0._r8
                 ! loop for each col 
                 do icol=1,pcols ! calculate if EM emitted
                   frcing_col_plume = 0._r8 !set initial value for each column
                   emis_col = sum(forcings(m)%fields(isec)%data(icol,:,lchnk))
                   if ( emis_col  > 0.0_r8) then
                         ! get initial emission from forcing data
                         frcing_vertical_plume_old(1:pver) = forcings(m)%fields(isec)%data(icol,pver:1:-1,lchnk) ! reverse
                         ! calculate the total emission in this column
                         do k = 1,pver
                            ! unit is molecular/m2
                            frcing_col_plume = frcing_col_plume +  &
                                                     frcing_vertical_plume_old(k)*(zint(icol,k)-zint(icol,k+1))*km_to_cm
                           ! write(iulog,*)'kzm_level ',k, 'old emis ', frcing_vertical_plume_old(k)
                         enddo
                         ! heatflux related to total carbon burned
                         burnedarea = (frcing_col_plume) ! burned  m2/m2
                         burnedarea_memory(icol) = burnedarea*area(icol) ! m2 burned
                   endif
                 enddo
                
              endif ! species and sector ifs
           endif
       enddo  !isec loop

       !xfcname = trim(forcings(m)%species)//'_XFRC'
       !call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )

       !frcing_col(:ncol) = 0._r8
       !do k = 1,pver
       !   frcing_col(:ncol) = frcing_col(:ncol) + frcing(:ncol,k,n)*(zint(:ncol,k)-zint(:ncol,k+1))*km_to_cm
       !enddo
       !xfcname = trim(forcings(m)%species)//'_CLXF'
       !call outfld( xfcname, frcing_col(:ncol), ncol, lchnk )

    end do src_loop1



    src_loop2 : do m = 1,extfrc_cnt

      n = forcings(m)%frc_ndx

       frcing(:ncol,:,n) = 0._r8
       do isec = 1,forcings(m)%nsectors
          ! move this part to the bottom of the loop ---------------
          !if (forcings(m)%file%alt_data) then
          !   frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,pver:1:-1,lchnk)
          !else
          !   frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,:,lchnk)
          !endif
          !----------------------------------------------------------

          !check wildfire emission EM
          if ((plumerise) .and. (forcings(m)%file%alt_data)) then
             !plume-rise calculation for certain species
             if ((m == PH_emis_m(1) .and. isec == PH_emis_n(1)) .or. & 
                 (m == PH_emis_m(2) .and. isec == PH_emis_n(2)) .or. &
                 (m == PH_emis_m(3) .and. isec == PH_emis_n(3)) .or. &
                 (m == PH_emis_m(4) .and. isec == PH_emis_n(4)) .or. &
                 (m == PH_emis_m(5) .and. isec == PH_emis_n(5)) .or. &
                 (m == PH_emis_m(6) .and. isec == PH_emis_n(6)) )then
                !plume_height_EM(:ncol) = 0._r8
                !heat_flux_plume(:ncol) = 0._r8
                !frcing_col_plume = 0._r8
                frcing_vertical_plume_old(:pver) = 0._r8
                frcing_vertical_plume_new(:pver) = 0._r8
                ! loop for each col
                do icol=1,pcols ! calculate if EM emitted
                   frcing_col_plume = 0._r8 !set initial value for each column
                   emis_col = sum(forcings(m)%fields(isec)%data(icol,:,lchnk))
                   if ( emis_col  > 0.0_r8) then
                         ! get initial emission from forcing data
                         frcing_vertical_plume_old(1:pver) = forcings(m)%fields(isec)%data(icol,pver:1:-1,lchnk) ! reverse
                         ! calculate the total emission in this column
                         do k = 1,pver
                            ! unit is molecular/m2 
                            frcing_col_plume = frcing_col_plume +  &
                                                     frcing_vertical_plume_old(k)*(zint(icol,k)-zint(icol,k+1))*km_to_cm
                           ! write(iulog,*)'kzm_level ',k, 'old emis ', frcing_vertical_plume_old(k)
                         enddo
                         if (forcings(m)%species == 'bc_a4' .and. forcings(m)%sectors(isec) == 'EM' &
                             .and. heatflux_memory(icol) > 0.0)then
                         ! convert molecular/cm2/s to kw/m2/s, based on Wooster et al., 2005, equation 14                           
                         ! mass (kg/s) = emis*1.0E4/Avogadr_cst*12/1000
                         ! FRP (kW/m2) = mass/0.368*1000
                            frp =  heatflux_memory(icol)*area(icol)/burnedarea_memory(icol) ! KW/m2
                            burnedarea = burnedarea_memory(icol)
                            !frp = (frcing_col_plume*1.0e4/0.03_r8)/6.022e23_r8*12.0/0.368
                            !frp = frp*2 ! 
                            !frp = (frcing_col_plume*1.0e4/0.03_r8)/6.022e23_r8*12.0_r8 ! carbon g/m2/s
                            !frp = frp*2.0_r8 ! carbon to fuel, Rowell et al., 2013; unit g/m2/s
                            !frp = frp*1.0E-3_r8 ! convert to kg/m2/s
                            !frp = frp*19.6e6_r8 ! forest heat from burning (J/s/m2), Freitas et al., 2011  
                            !frp = frp*1.0e-3_r8 ! KW/m2
                            !frp_memory(icol) = frp 
                            tod_saved = tod
                            call cal_plume_height(plume_height,zmidr(icol,:), pmid(icol,:), &
                                 tfld(icol,:), relhum(icol,:), qh2o(icol,:), ufld(icol,:), &
                                 vfld(icol,:), clat(icol)/(3.1415_r8)*180.0_r8, &
                                 clon(icol)/(3.1415_r8)*180.0_r8, tl, pt_v, frp, burnedarea,frp4plume, &
                                 wt_ini_e3sm_out(icol,:), wt_end_e3sm_out(icol,:), &
                                 rbuoy_ini_e3sm_out(icol,:), rbuoy_end_e3sm_out(icol,:), &
                                 t_ini_e3sm_out(icol,:), t_end_e3sm_out(icol,:),area(icol), &
                                 qv_ini_e3sm_out(icol,:), qv_end_e3sm_out(icol,:),&
                                 r_ini_e3sm_out(icol,:), r_end_e3sm_out(icol,:), &
                                 rho_ini_e3sm_out(icol,:), rho_end_e3sm_out(icol,:) )
                            plume_height = plume_height + zmidr(icol,pver)*1000.0_r8 ! plume height at middle of bottom layer
                            plume_height_EM(icol) = plume_height ! in meter
                            heat_flux_plume(icol) = frp4plume ! in kw/m2 
                            if (emis_constrained_frp) then
                               burnedarea_plume(icol) = burnedarea 
                            else
                               burnedarea_plume(icol) = area(icol)
                            endif
                         elseif (forcings(m)%species == 'pom_a4' .or. forcings(m)%species == 'brc_a4')then
                            !frp = (frcing_col_plume*1.0e4/0.97_r8)/6.022e23_r8*12.0/0.368
                            !frp = frp*10.0_r8 !
                            !call cal_plume_height(plume_height,zmidr(icol,:), pmid(icol,:), &
                            !     tfld(icol,:), relhum(icol,:), qh2o(icol,:), ufld(icol,:), &
                            !     vfld(icol,:), clat(icol)/(3.1415_r8)*180.0_r8, clon(icol)/(3.1415_r8)*180.0_r8, tl, pt_v,frp)
                            plume_height = plume_height_EM(icol) 
                         elseif (forcings(m)%species == 'num_a4' .or. forcings(m)%species == 'H2OFIRE')then
                            !if at same timestep and other fire species        
                             !frp = frp_memory(icol)
                             !call cal_plume_height(plume_height,zmidr(icol,:), pmid(icol,:), &
                             !    tfld(icol,:), relhum(icol,:), qh2o(icol,:), ufld(icol,:), &
                             !    vfld(icol,:), clat(icol)/(3.1415_r8)*180.0_r8, clon(icol)/(3.1415_r8)*180.0_r8, tl, pt_v,frp)    
                            plume_height = plume_height_EM(icol) 
                         endif 
                         
                         ! match plume height to model vertical grid
                         ph_z(icol) = pver
                         do k = 2, pver
                            if (fix_plume_height > 0.0) then ! here use fixed plume height from namelist
                                plume_height = fix_plume_height ! unit: meter
                                plume_height_EM(icol) = fix_plume_height ! unit: meter
                            endif
                            if ((plume_height - zmidr(icol,k)*1000_r8) > 0._r8 .and. (plume_height - zmidr(icol,k-1)*1000_r8 < 0._r8)) then
                               ph_z(icol) = k
                            endif 
                         enddo
                         ! calculate mass flux
                         qmass(:) = 0.0_r8 
                         qratio(:) = 0.0_r8
                         qmass_1 = (3.1415_r8*(wt_end_e3sm_out(icol,pver)*(r_end_e3sm_out(icol,pver)**2) &
                                       *air_density(icol,pver)*qv_end_e3sm_out(icol,pver)))
                         do k = ph_z(icol), pver-1
                            !qmass(k) = 3.1415_r8*(wt_end_e3sm_out(icol,k+1)*(r_end_e3sm_out(icol,k+1)**2) & 
                            !           *air_density(icol,k+1)*qv_end_e3sm_out(icol,k+1) &
                            !           - wt_end_e3sm_out(icol,k-1)*(r_end_e3sm_out(icol,k-1)**2) &
                            !           *air_density(icol,k-1)*qv_end_e3sm_out(icol,k-1)) 
                            wdqdz = -wt_end_e3sm_out(icol,k)*(qv_end_e3sm_out(icol,k+1)-qv_end_e3sm_out(icol,k-1)) &
                                    /(zmidr(icol,k+1)-zmidr(icol,k-1)) 
                            qdwdz = -qv_end_e3sm_out(icol,k)*(wt_end_e3sm_out(icol,k+1)-wt_end_e3sm_out(icol,k-1)) &
                                    /(zmidr(icol,k+1)-zmidr(icol,k-1))
                            qmass(k) = 3.1415_r8*(r_end_e3sm_out(icol,k)**2)*rho_end_e3sm_out(icol,k)*0.001_r8 &
                                       * (wdqdz+qdwdz) * (zint(icol,k)-zint(icol,k+1))
                            !qmass(k) = (3.1415_r8*(r_end_e3sm_out(icol,k+1)**2)*rho_end_e3sm_out(icol,k+1)*0.001_r8 &
                            !           *wt_end_e3sm_out(icol,k+1)*qv_end_e3sm_out(icol,k+1) &
                            !          - 3.1415_r8*(r_end_e3sm_out(icol,k-1)**2)*rho_end_e3sm_out(icol,k-1)*0.001_r8 &
                            !           *wt_end_e3sm_out(icol,k-1)*qv_end_e3sm_out(icol,k-1)) &
                            !            * (zint(icol,k)-zint(icol,k+1))/(zmidr(icol,k+1)-zmidr(icol,k-1))  
                            !qmass(k) = 3.1415_r8*(r_end_e3sm_out(icol,k)**2)*rho_end_e3sm_out(icol,k)*0.001_r8 & ! rho is g/m3
                            !           *qv_end_e3sm_out(icol,k)*(zint(icol,k)-zint(icol,k+1))*1000.0_r8 & ! zint in km
                            !           - 3.1415_r8*(r_end_e3sm_out(icol,k)**2)*air_density(icol,k) & !
                            !           *qh2o(icol,k)*(zint(icol,k)-zint(icol,k+1))*1000.0_r8
                            qmass(k) = max(0.0_r8, qmass(k)) 
                            !qratio(k) = qmass(k) / qmass_1
                         enddo 
                         if (sum((qmass)) > 1.0E-10_r8 .and. qmass_1 > 1.0E-10_r8 ) then
                            qratio(:) = qmass(:) / sum(qmass(:)) ! normalize qratio
                         else
                            qratio(:) = 0.0_r8 ! no water emission if too small
                         endif
                         !v_air = burnedarea_memory(icol)*wt_end_e3sm_out(icol,pver) !(m3/s) The volume at the bottom per second
                         !m_vapor = v_air*air_density(icol,pver)*qh2o(icol,pver) ! The mass of water vapor (kg/s) 
                         !qmass(pver) = 3.1415_r8*(wt_end_e3sm_out(icol,pver)*(r_end_e3sm_out(icol,pver)**2) &
                         !              *air_density(icol,pver)*qv_end_e3sm_out(icol,pver))
                         if (forcings(m)%species == 'bc_a4' .and. diag_run_plumerise )then
                         write(iulog,*) 'fix_plume_height ', fix_plume_height 
                         write(iulog,*) 'kzm_fire_species ', forcings(m)%species, isec
                         write(iulog,*) 'kzm_heatflux ', heatflux_memory(icol)
                         write(iulog,*) 'kzm_FRP ',frp,frp4plume
                         write(iulog,*) 'kzm_ratio_2_area ',burnedarea_memory(icol)/area(icol), area(icol)
                         write(iulog,*)'kzm_plume_rise_calculation_running'
                         write(iulog,*)'kzm_plume_rise_calculation_lat ', clat(icol)/(3.1415_r8)*180.0_r8
                         write(iulog,*)'kzm_plume_rise_calculation_lon ', clon(icol)/(3.1415_r8)*180.0_r8
                         write(iulog,*)'kzm_plume_time ', iyear,imo,iday_m,tod
                         write(iulog,*)'kzm_plume_height ', plume_height, 'local time ',tl
                         write(iulog,*)'kzm_omega ', wt_ini_e3sm_out(icol,pver),wt_end_e3sm_out(icol,pver)
                         write(iulog,*)'kzm_q_env ', qh2o(icol,pver), qh2o(icol,ph_z(icol))
                         write(iulog,*)'kzm_qv ', qv_end_e3sm_out(icol,pver),qv_end_e3sm_out(icol,ph_z(icol))
                         !write(iulog,*)'kzm_r_ini_end ', r_ini_e3sm_out(icol,ph_z(icol)),r_end_e3sm_out(icol,ph_z(icol))
                         !write(iulog,*)'kzm_qmass ', qmass(:)
                         write(iulog,*)'kzm_qratio ', qratio(:)
                         !write(iulog,*)'kzm_qratio_sum ', sum(qratio(:))
                         write(iulog,*)'kzm_qmass_1 ', qmass_1, m_vapor
                        
                         !write(iulog,*)'kzm_plume_FRP ', frp
                         !write(iulog,*)'kzm_plume environment data begin: zmidr pmid tfld relhum qh2o ufld vfld '
                         !do k = pver,50,-1
                         !   write(iulog,*)k, wt_ini_e3sm_out(icol,k),wt_end_e3sm_out(icol,k), &
                         !                   rbuoy_ini_e3sm_out(icol,k),rbuoy_end_e3sm_out(icol,k), &
                         !                    t_ini_e3sm_out(icol,k),t_end_e3sm_out(icol,k)
                         !enddo
                         write(iulog,*)'kzm_col_total_emis ', frcing_col_plume
                         write(iulog,*)'kzm_plume_layer ', ph_z(icol)
                         endif
                         ! redistrict the emission
                         
                         do k = 1,pver
                            ! option: release all emission into plume top layer
                            if (k == ph_z(icol)) then
                               frcing_col_plume = frcing_col_plume_bc_memory(icol)  
                               if (forcings(m)%species == 'bc_a4' .and. forcings(m)%sectors(isec) == 'EM')then
                                  frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm)
                                  frcing_vertical_plume_new(k) =  frcing_vertical_plume_new(k)*ef_bc_a4
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_bc_a4_emis_at_layer ', k, frcing_vertical_plume_new(k)
                                  endif
                               elseif (forcings(m)%species == 'pom_a4' .and. forcings(m)%sectors(isec) == 'EM')then
                                  frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm)
                                  frcing_vertical_plume_new(k) =  frcing_vertical_plume_new(k)*ef_oc_a4
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_oc_a4_emis_at_layer ', k, frcing_vertical_plume_new(k)
                                  endif
                               elseif (forcings(m)%species == 'brc_a4' .and. forcings(m)%sectors(isec) == 'EM')then
                                  frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm)
                                  frcing_vertical_plume_new(k) =  frcing_vertical_plume_new(k)*ef_oc_a4
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_oc_a4_emis_at_layer ', k, frcing_vertical_plume_new(k)
                                  endif
                               elseif (forcings(m)%species == 'H2OFIRE-NON'.and. forcings(m)%sectors(isec) == 'kEM')then
                                  frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm)
                                  frcing_vertical_plume_new(k) =  frcing_vertical_plume_new(k)*ef_h2o_a4
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_H2O_fuel_emis_at_layer ', k, frcing_vertical_plume_new(k)
                                  endif
                                  !if (surface_flux_flag ) then ! water flux for plume-height 5 layers height 
                                  !   v_air = burnedarea_memory(icol)*wt_end_e3sm_out(icol,pver) !(m3/s) The volume at the bottom per second
                                  !   m_vapor = v_air*air_density(icol,pver)*qh2o(icol,pver) ! The mass of water vapor (kg/s) 
                                  !   f_vapor = m_vapor/area(icol)*1000.0_r8/10000.0_r8 ! g/cm2/s 
                                  !   f_vapor = f_vapor/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm) !g/cm3/s
                                  !   f_vapor = f_vapor/18.0_r8*avogadro !moleculer/cm3/s (H2O)  
                                  !   frcing_vertical_plume_new(k) = frcing_vertical_plume_new(k) + f_vapor*detrainment_para  
                                      
                                  !endif
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_H2O_flux_emis_at_layer ', k, detrainment_para, frcing_vertical_plume_new(k)
                                  endif
                               elseif (forcings(m)%species == 'num_a4' .and. forcings(m)%sectors(isec) == 'num_a1_BC_ELEV_EM')then
                                  frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm) 
                                  frcing_vertical_plume_new(k) =  frcing_vertical_plume_new(k)*ef_bc_a4/coef_bc
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_num_a4_emis_at_layer ', k, frcing_vertical_plume_new(k)
                                  endif
                               elseif (forcings(m)%species == 'num_a4' .and. forcings(m)%sectors(isec) == 'num_a1_POM_ELEV_EM')then
                                  frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm)
                                  frcing_vertical_plume_new(k) =  frcing_vertical_plume_new(k)*ef_oc_a4/coef_pom 
                                  if (diag_run_plumerise) then
                                     write(iulog,*) 'kzm_num_a4_emis_at_layer ', k, frcing_vertical_plume_new(k)
                                  endif
                               endif
                            else
                               frcing_vertical_plume_new(k) =  0.0_r8 
                            endif 
                            if (forcings(m)%species == 'H2OFIRE'.and. forcings(m)%sectors(isec) == 'EM')then
                               frcing_vertical_plume_new(k) =  0.0_r8
                               m_vapor = max(qmass_1,0.0_r8)*qratio(k)*detrainment_para
                               ! phasing out h2o transport as burned area increase
                               phase_out = cosd(burnedarea_memory(icol)/area(icol)*90.0_r8) 
                               phase_out = max(phase_out,0.0_r8)
                               m_vapor = m_vapor*phase_out
                               f_vapor = m_vapor*1000.0_r8/10000.0_r8/area(icol) ! g/cm2/s
                               f_vapor = f_vapor/(abs(zint(icol,k)-zint(icol,k+1))*km_to_cm) !g/cm3/s
                               f_vapor = f_vapor/18.0_r8*avogadro !moleculer/cm3/s (H2O)
                               frcing_vertical_plume_new(k) = frcing_vertical_plume_new(k) + f_vapor
                            endif                              

                            ! option: release emission evenly from plume top to surface 
                            !if (k >= ph_z(icol)) then
                            !   frcing_vertical_plume_new(k) =  frcing_col_plume/(abs(zint(icol,pver+1)-zint(icol,ph_z(icol)))*km_to_cm)
                            !else
                            !   frcing_vertical_plume_new(k) = 0.0_r8
                            !endif
                            !write(iulog,*)'kzm_level ',k, 'new emis ', frcing_vertical_plume_new(k) 
                         enddo 
                         forcings(m)%fields(isec)%data(icol,:,lchnk) = frcing_vertical_plume_new(pver:1:-1) ! reverse back    
                     end if
                  enddo 
                endif ! if fire emission released 
            endif !plumerise flag   

            ! back to no fire plume calculation 
            ! add emission from different sectors together
            if (forcings(m)%file%alt_data) then
               ! avoid non EM sector from H2OFIRE and brc_a4
               if (forcings(m)%species == 'H2OFIRE'.and. forcings(m)%sectors(isec) /= 'EM') then
                    frcing(:ncol,:,n) = frcing(:ncol,:,n)
               elseif (forcings(m)%species == 'brc_a4'.and. forcings(m)%sectors(isec) /= 'EM') then    
                    frcing(:ncol,:,n) = frcing(:ncol,:,n)
               else
                    frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,pver:1:-1,lchnk)
               endif
            else
               frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings(m)%fields(isec)%data(:ncol,:,lchnk)
            endif
       enddo

       xfcname = trim(forcings(m)%species)//'_XFRC'
       call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )

       frcing_col(:ncol) = 0._r8
       do k = 1,pver
          frcing_col(:ncol) = frcing_col(:ncol) + frcing(:ncol,k,n)*(zint(:ncol,k)-zint(:ncol,k+1))*km_to_cm
       enddo
       xfcname = trim(forcings(m)%species)//'_CLXF'
       call outfld( xfcname, frcing_col(:ncol), ncol, lchnk )
       
    end do src_loop2
    ! output plume heights and environment
       if (plumerise)then
         call outfld( 'plume_height_EM', plume_height_EM, ncol, lchnk )
         call outfld( 'heat_flux_plume', heat_flux_plume, ncol, lchnk )
         call outfld( 'burned_area', burnedarea_plume, ncol, lchnk )
         call outfld('zmidr_ph', zmidr(:ncol,:), ncol, lchnk)
         call outfld('pmid_ph', pmid(:ncol,:), ncol, lchnk)
         call outfld('tfld_ph', tfld(:ncol,:), ncol, lchnk)
         call outfld('relhum_ph', relhum(:ncol,:), ncol, lchnk)
         call outfld('qh2o_ph', qh2o(:ncol,:), ncol, lchnk)
         call outfld('ufld_ph', ufld(:ncol,:), ncol, lchnk)
         call outfld('vfld_ph', vfld(:ncol,:), ncol, lchnk)
         call outfld('wt_ini_e3sm_out', wt_ini_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('t_ini_e3sm_out', t_ini_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('wt_end_e3sm_out', wt_end_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('t_end_e3sm_out', t_end_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('rbuoy_ini_e3sm_out', rbuoy_ini_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('rbuoy_end_e3sm_out', rbuoy_end_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('qv_ini_e3sm_out', qv_ini_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('qv_end_e3sm_out', qv_end_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('r_ini_e3sm_out', r_ini_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('r_end_e3sm_out', r_end_e3sm_out(:ncol,:), ncol, lchnk)
         call outfld('rho_end_e3sm_out', rho_end_e3sm_out(:ncol,:), ncol, lchnk)
       endif

  end subroutine extfrc_set

! subroutines for plumerise
  subroutine cal_plume_height( plume_height,zmidr_v, pmid_v, tfld_v, relhum_v, qh2o_v, &
                               ufld_v, vfld_v,lat,lon,tl,pt_v,frp,burnedarea,frp4plume,&
                               wt_ini_e3sm_col, wt_end_e3sm_col, rbuoy_ini_e3sm_col, rbuoy_end_e3sm_col,&
                               t_ini_e3sm_col,t_end_e3sm_col,gridarea, &
                               qv_ini_e3sm_col,qv_end_e3sm_col, &
                               r_ini_e3sm_col,r_end_e3sm_col,  &
                               rho_ini_e3sm_col,rho_end_e3sm_col   )
    use smk_plumerise, only : smk_pr_driver  
    !use time_manager,  only: get_curr_date
    implicit none
    ! plume-rise parameters
    real(r8), intent(out)  ::   plume_height
    real(r8), intent(in)  ::   zmidr_v(pver)             ! midpoint geopot height - elevation ( km )
    real(r8), intent(in)  ::   pmid_v(pver)            ! midpoint pressure (Pa)
    real(r8), intent(in)  ::   tfld_v(pver)            ! midpoint temperature (K)
    real(r8), intent(in)  ::   relhum_v(pver)           ! relative humidity (0~1)
    real(r8), intent(in)  ::   qh2o_v(pver)            ! specific humidity (kg/kg)
    real(r8), intent(in)  ::   ufld_v(pver)            ! zonal velocity (m/s)
    real(r8), intent(in)  ::   vfld_v(pver)            ! meridional velocity (m/s)
    real(r8), intent(in)  ::   frp
    real(r8), intent(in)  ::   burnedarea,gridarea
    real(r8), intent(out)  ::   frp4plume, tl
    real(r8), intent(out) ::  wt_ini_e3sm_col(pver), wt_end_e3sm_col(pver), & 
                              rbuoy_ini_e3sm_col(pver), rbuoy_end_e3sm_col(pver), &
                              t_ini_e3sm_col(pver), t_end_e3sm_col(pver), &
                              qv_ini_e3sm_col(pver), qv_end_e3sm_col(pver), &
                              r_ini_e3sm_col(pver), r_end_e3sm_col(pver),& 
                              rho_ini_e3sm_col(pver), rho_end_e3sm_col(pver) 
    ! local variables
    real(r8)  :: env(8, pver) ! meterology profiles for this column 
    real(r8)  :: gfed_area  ! fire parameters
    real(r8)  :: lat,lon !
    real(r8)  :: pt_v(pver)  ! potential temperature
    integer :: i,ihr,imn,iyear,imo,iday_m,tod 
    real(r8) :: frp_peak,frp_h,frp_b,frp_sigma
    real(r8) :: wt_ini_e3sm(pver), wt_end_e3sm(pver), rbuoy_ini_e3sm(pver), rbuoy_end_e3sm(pver)
    real(r8) :: t_ini_e3sm(pver), t_end_e3sm(pver)
    real(r8) :: qv_ini_e3sm(pver), qv_end_e3sm(pver)
    real(r8) :: r_ini_e3sm(pver), r_end_e3sm(pver)
    real(r8) :: rho_ini_e3sm(pver), rho_end_e3sm(pver)
   ! get plume height
   ! env: geopotential height, pressure, temp(state%t), relative humidity(state%),
   ! env: potential T, specific humidity, U, V

    !call smk_pr_driver(plume_height , env, gfed_area, frp, lat )
    ! for fire at each column
    
       call cal_theta(tfld_v(:), pmid_v(:)/100.0_r8, pt_v(:))
       env(1,:) = zmidr_v(pver:1:-1)*1000.0_r8 ! meter
       env(1,:) = env(1,:) - env(1,1) ! set first layer at zero
       env(2,:) = pmid_v(pver:1:-1) 
       env(3,:) = tfld_v(pver:1:-1)
       env(4,:) = relhum_v(pver:1:-1)
       env(5,:) = pt_v(pver:1:-1)
       env(6,:) = qh2o_v(pver:1:-1)
       env(7,:) = ufld_v(pver:1:-1)
       env(8,:) = vfld_v(pver:1:-1)
       !plume_height = 1000.0
       !gfed_area = burnedarea/10000.0_r8
       !lat = 100.0
       !frp diurnal cycle based on Ke et al., 2021 (WTNA region diurnal cycle)
       !this cycle also similar to WRF-chem fire diurnal cycle used in Kumar et al., 2022
       frp_peak = 155.25_r8 !km/m2 ! 15.525 time 10 as FRP *10 rule
       frp_h = 14.002_r8
       frp_b = 0.023
       frp_sigma = 2.826
       ! local time
       call get_curr_date (iyear,imo,iday_m,tod)  ! year, time of day [sec]
       ihr  = tod/3600
       imn  = mod( tod,3600 )/3600.0_r8
       if (lon > 180.0_r8) then
          tl = ihr*1.0_r8 + imn*1.0_r8 + (lon-360.0)/15.0_r8 ! behind of UTC
       else
          tl = ihr*1.0_r8 + imn*1.0_r8 + (360.0-lon)/15.0_r8 ! ahead of UTC
       endif 
       if (tl > 24.0) then
          tl = mod( tl,24.0 )
       elseif (tl < 0.0) then
          tl = tl + 24.0
       endif 
       if (emis_constrained_frp) then
          frp4plume = frp
          !frp4plume= frp*36.0_r8!80._r8
          gfed_area = burnedarea/10000.0_r8
       else 
          frp4plume=frp_peak*(exp(-0.5*(tl-frp_h)*(tl-frp_h)/frp_sigma/frp_sigma)+frp_b) 
          !gfed_area=gridarea
          gfed_area = burnedarea/10000.0_r8 !still use burned area from emission file  
       endif
       call smk_pr_driver(plume_height , env, gfed_area, frp4plume, lat, &
            wt_ini_e3sm, wt_end_e3sm, rbuoy_ini_e3sm, rbuoy_end_e3sm,t_ini_e3sm,t_end_e3sm &
            , qv_ini_e3sm, qv_end_e3sm, r_ini_e3sm, r_end_e3sm, rho_ini_e3sm,rho_end_e3sm )
       wt_ini_e3sm_col(:) = wt_ini_e3sm(pver:1:-1)
       t_ini_e3sm_col(:) = t_ini_e3sm(pver:1:-1)
       wt_end_e3sm_col(:) = wt_end_e3sm(pver:1:-1)
       t_end_e3sm_col(:) = t_end_e3sm(pver:1:-1)
       rbuoy_ini_e3sm_col(:) = rbuoy_ini_e3sm(pver:1:-1)
       rbuoy_end_e3sm_col(:) = rbuoy_end_e3sm(pver:1:-1)
       qv_ini_e3sm_col(:) = qv_ini_e3sm(pver:1:-1)
       qv_end_e3sm_col(:) = qv_end_e3sm(pver:1:-1)
       r_ini_e3sm_col(:) = r_ini_e3sm(pver:1:-1)
       r_end_e3sm_col(:) = r_end_e3sm(pver:1:-1)
       rho_ini_e3sm_col(:) = rho_ini_e3sm(pver:1:-1)
       rho_end_e3sm_col(:) = rho_end_e3sm(pver:1:-1)
       !if(masterproc) write(iulog,*) 'plume_height ', plume_height
  end subroutine cal_plume_height
!----------------------------------------------------------------------------------------
 subroutine cal_theta(temp,pe,theta)
 ! this code is used to calculate the potential temperature
 ! function [theta]=cal_theta(temp, pe)
 ! temp is temperature in K
 ! pe is pressure in mb
 ! temp and pe is array
 implicit none
 integer::i
 real(r8)::P,R,cp,ca1,ca2
 real(r8),intent(in)::temp(pver),pe(pver)
 real(r8),intent(out)::theta(pver)
 P=1000.0 ! reference level pressure 1000 mb
 R=286.9  ! J/(kg*k)
 cp = 1004. ! J/(kg*k)
 ca1 = R/cp
 do i=1,pver
    ca2=(P/pe(i))**ca1
    theta(i)=temp(i)*ca2
 enddo
 end subroutine cal_theta 

end module mo_extfrc
