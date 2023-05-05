module mo_srf_emissions
  !---------------------------------------------------------------
  ! 	... surface emissions module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : gas_pcnst
  use spmd_utils,   only : masterproc,iam
  use mo_tracname,  only : solsym
  use cam_abortutils,   only : endrun
  use ioFileMod,    only : getfil
  use ppgrid,       only : pcols, begchunk, endchunk
  use cam_logfile,  only : iulog
  use tracer_data,  only : trfld,trfile

  implicit none

  type :: emission
     integer           :: spc_ndx
     real(r8)          :: mw
     character(len=256):: filename
     character(len=16) :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type emission

  private

  public  :: srf_emissions_inti, set_srf_emissions, set_srf_emissions_time 

  save

  real(r8), parameter :: amufac = 1.65979e-23_r8         ! 1.e4* kg / amu
  logical :: has_emis(gas_pcnst)
  type(emission), allocatable :: emissions(:)
  integer                     :: n_emis_species 
  integer :: c10h16_ndx, isop_ndx

contains

  subroutine srf_emissions_inti( srf_emis_specifier, emis_type, emis_cycle_yr, emis_fixed_ymd, emis_fixed_tod )

    !-----------------------------------------------------------------------
    ! 	... initialize the surface emissions
    !-----------------------------------------------------------------------

    use chem_mods,        only : adv_mass
    use mo_constants,     only : d2r, pi, rearth
    use string_utils,     only : to_upper
    use mo_chem_utls,     only : get_spc_ndx 
    use tracer_data,      only : trcdata_init
    use cam_pio_utils,    only : cam_pio_openfile
    use pio,              only : pio_inquire, pio_nowrite, pio_closefile, pio_inq_varndims
    use pio,              only : pio_inq_varname, file_desc_t
    use pio,              only : pio_inq_vardimid, pio_inq_dimid !(zhang73)
    use chem_surfvals,    only : flbc_list

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: srf_emis_specifier(:)
    character(len=*), intent(in) :: emis_type
    integer,          intent(in) :: emis_cycle_yr
    integer,          intent(in) :: emis_fixed_ymd
    integer,          intent(in) :: emis_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer  :: astat
    integer  :: j, l, m, n, i, nn                     ! Indices
    character(len=16)  :: spc_name
    character(len=256) :: filename

    character(len=16)  ::    emis_species(gas_pcnst)
    character(len=256) ::    emis_filenam(gas_pcnst)
    integer ::    emis_indexes(gas_pcnst)

    integer :: vid, nvars, isec
    integer :: dimids(8), did, dimid,ncol_dimid,lat_dimid,time_dimid !(zhang73)
    integer, allocatable :: vndims(:)
    integer, allocatable :: finddim_time(:), finddim_lat_ncol(:) !(zhang73) 
    type(file_desc_t) :: ncid
    character(len=32)  :: varname
    character(len=256) :: locfn
    integer :: ierr
    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.

    has_emis(:) = .false.
    nn = 0

    count_emis: do n=1,gas_pcnst
       if ( len_trim(srf_emis_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(srf_emis_specifier(n),'->')
       spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))
       filename = trim(adjustl(srf_emis_specifier(n)(i+2:)))

       m = get_spc_ndx(spc_name)

       if (m > 0) then
          has_emis(m) = .true.
          has_emis(m) = has_emis(m) .and. ( .not. any( flbc_list == spc_name ) )
       else 
          write(iulog,*) 'srf_emis_inti: spc_name ',spc_name,' is not included in the simulation'
          call endrun('srf_emis_inti: invalid surface emission specification')
       endif

       if ( has_emis(m) ) then
          nn = nn+1
          emis_species(nn) = spc_name
          emis_filenam(nn) = filename
          emis_indexes(nn) = m
       endif
    enddo count_emis

    n_emis_species = count(has_emis(:))

    if (masterproc) write(iulog,*) 'srf_emis_inti: n_emis_species = ',n_emis_species

    allocate( emissions(n_emis_species), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'srf_emis_inti: failed to allocate emissions array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the emission type array
    !-----------------------------------------------------------------------
    do m=1,n_emis_species 
       emissions(m)%spc_ndx          = emis_indexes(m)
       emissions(m)%units            = 'Tg/y'
       emissions(m)%species          = emis_species(m)
       emissions(m)%mw               = adv_mass(emis_indexes(m))                     ! g / mole
       emissions(m)%filename         = emis_filenam(m)
    enddo

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    spc_loop: do m = 1, n_emis_species

       emissions(m)%nsectors = 0
       
       call getfil (emissions(m)%filename, locfn, 0)
       call cam_pio_openfile ( ncid, trim(locfn), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nvariables=nvars)

       allocate(vndims(nvars))
       allocate(finddim_time(nvars))
       allocate(finddim_lat_ncol(nvars))
       finddim_time=0
       finddim_lat_ncol=0
       time_dimid=-9999
       lat_dimid=-9999
       ncol_dimid=-9999

       !(zhang73) Qi: get file%is_ncol using their dimension -> pio_inq_dimid(ncid, 'ncol', dimid) 
       !(zhang73) determine vaild sectors(fields) from dim name "ncol/lat" && "time" (finddim_lat_ncol, finddim_time)
       ierr = pio_inq_dimid(ncid, 'time', dimid)
       if(ierr==0) time_dimid = dimid
       ierr = pio_inq_dimid(ncid, 'lat', dimid)
       if(ierr==0) lat_dimid = dimid
       ierr = pio_inq_dimid(ncid, 'ncol', dimid)
       emissions(m)%file%is_ncol = (ierr==0)
       if(ierr==0) ncol_dimid = dimid
       if(masterproc) write(iulog,*) '(zhang73 srf_emis_inti) time_dimid, lat_dimid, ncol_dimid=',time_dimid, lat_dimid, ncol_dimid

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, vndims(vid))

          ierr = pio_inq_vardimid (ncid, vid, dimids(1:vndims(vid))) !(zhang73)
          do did=1,vndims(vid)
             if( dimids(did) == time_dimid ) finddim_time(vid)=1
             if(  dimids(did) == lat_dimid ) finddim_lat_ncol(vid)=1
             if( dimids(did) == ncol_dimid ) finddim_lat_ncol(vid)=1
          enddo
   
          ierr = pio_inq_varname (ncid, vid, varname)
          if( finddim_time(vid)==1 .and. finddim_lat_ncol(vid)==1)then !(zhang73)
             !write(iulog,*) '(zhang73 srf_emis_inti) valid var: finddim_time(vid), finddim_lat_ncol(vid)=',trim(varname),finddim_time(vid), finddim_lat_ncol(vid)
             emissions(m)%nsectors = emissions(m)%nsectors+1
          else
             !write(iulog,*) 'srf_emis_inti: Skipping variable ', trim(varname),', ndims = ',vndims(vid),' , species=',trim(emissions(m)%species)
          end if

       enddo

       allocate( emissions(m)%sectors(emissions(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'srf_emis_inti: failed to allocate emissions(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1

       do vid = 1,nvars
          !if( vndims(vid) == dim_thres ) then !(zhang73) check vndims from 3 -> 2 to activate bc_a4_ncol
          if( finddim_time(vid)==1 .and. finddim_lat_ncol(vid)==1)then !(zhang73)
             ierr = pio_inq_varname(ncid, vid, emissions(m)%sectors(isec))
             isec = isec+1
          endif

       enddo
       deallocate(vndims)
       deallocate(finddim_time)
       deallocate(finddim_lat_ncol)
       call pio_closefile (ncid)

       allocate(emissions(m)%file%in_pbuf(size(emissions(m)%sectors)))
       emissions(m)%file%in_pbuf(:) = .false.
       if(masterproc) write(iulog,*) '(zhang73 srf_emis_inti) nvars, emissions(m)%species, emissions(m)%nsectors = ',nvars, emissions(m)%species,emissions(m)%nsectors
       call trcdata_init( emissions(m)%sectors, &
                          emissions(m)%filename, filelist, datapath, &
                          emissions(m)%fields,  &
                          emissions(m)%file, &
                          rmv_file, emis_cycle_yr, emis_fixed_ymd, emis_fixed_tod, emis_type )

    enddo spc_loop

    c10h16_ndx = get_spc_ndx('C10H16')
    isop_ndx = get_spc_ndx('ISOP')

  end subroutine srf_emissions_inti

  subroutine set_srf_emissions_time( pbuf2d, state )
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

    do m = 1,n_emis_species
       if(masterproc.and.emissions(m)%file%is_ncol) write(iulog,*) '(zhang73 set_srf_emissions_time) emissions(m)%species: ', emissions(m)%species 
       call advance_trcdata( emissions(m)%fields, emissions(m)%file, state, pbuf2d  )
    end do

  end subroutine set_srf_emissions_time

  ! adds surf flux specified in file to sflx
  subroutine set_srf_emissions( lchnk, ncol, sflx )
    !--------------------------------------------------------
    !	... form the surface fluxes for this latitude slice
    !--------------------------------------------------------

    use mo_constants, only : pi
    use time_manager, only : get_curr_calday
    use string_utils, only : to_lower, GLC
    use phys_grid,    only : get_rlat_all_p, get_rlon_all_p

    implicit none

    !--------------------------------------------------------
    !	... Dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)  :: ncol                  ! columns in chunk
    integer,  intent(in)  :: lchnk                 ! chunk index
    real(r8), intent(out) :: sflx(:,:) ! surface emissions ( kg/m^2/s )

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n
    real(r8) ::  factor
    real(r8) ::  dayfrac            ! fration of day in light
    real(r8) ::  iso_off            ! time iso flux turns off
    real(r8) ::  iso_on             ! time iso flux turns on

    logical  :: polar_day,polar_night
    real(r8) :: doy_loc
    real(r8) :: sunon,sunoff
    real(r8) :: loc_angle
    real(r8) :: latitude
    real(r8) :: declination
    real(r8) :: tod
    real(r8) :: calday

    real(r8), parameter :: dayspy = 365._r8
    real(r8), parameter :: twopi = 2.0_r8 * pi
    real(r8), parameter :: pid2  = 0.5_r8 * pi
    real(r8), parameter :: dec_max = 23.45_r8 * pi/180._r8

    real(r8) :: flux(ncol)
    real(r8) :: mfactor
    integer  :: isec

    character(len=12),parameter :: mks_units(4) = (/ "kg/m2/s     ", &
                                                     "kg/m2/sec   ", &
                                                     "kg/m^2/s    ", &
                                                     "kg/m^2/sec  " /)
    character(len=12) :: units

    real(r8), dimension(ncol) :: rlats, rlons 

    sflx(:,:) = 0._r8

    !--------------------------------------------------------
    !	... set non-zero emissions
    !--------------------------------------------------------
    emis_loop : do m = 1,n_emis_species

       n = emissions(m)%spc_ndx

       flux(:) = 0._r8
       do isec = 1,emissions(m)%nsectors
          flux(:ncol) = flux(:ncol) + emissions(m)%fields(isec)%data(:ncol,1,lchnk)
       enddo

       units = to_lower(trim(emissions(m)%fields(1)%units(:GLC(emissions(m)%fields(1)%units))))
       
       if ( any( mks_units(:) == units ) ) then
          sflx(:ncol,n) = flux(:ncol)
       else
          mfactor = amufac * emissions(m)%mw
          sflx(:ncol,n) = flux(:ncol) * mfactor
       endif

    end do emis_loop

    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )

    calday = get_curr_calday()
    doy_loc     = aint( calday )
    declination = dec_max * cos((doy_loc - 172._r8)*twopi/dayspy)
    tod = (calday - doy_loc) + .5_r8

    do i = 1,ncol
       !
       polar_day   = .false.
       polar_night = .false.
       !
       loc_angle = tod * twopi + rlons(i)
       loc_angle = mod( loc_angle,twopi )
       latitude =  rlats(i)
       !
       !------------------------------------------------------------------
       !        determine if in polar day or night
       !        if not in polar day or night then
       !        calculate terminator longitudes
       !------------------------------------------------------------------
       if( abs(latitude) >= (pid2 - abs(declination)) ) then
          if( sign(1._r8,declination) == sign(1._r8,latitude) ) then
             polar_day = .true.
             sunoff = 2._r8*twopi
             sunon  = -twopi
          else
             polar_night = .true.
          end if
       else
          sunoff = acos( -tan(declination)*tan(latitude) )
          sunon  = twopi - sunoff
       end if

       !--------------------------------------------------------
       !	... adjust alpha-pinene for diurnal variation
       !--------------------------------------------------------
       if( c10h16_ndx > 0 ) then
          if( has_emis(c10h16_ndx) ) then
             if( .not. polar_night .and. .not. polar_day ) then
                dayfrac = sunoff / pi
                sflx(i,c10h16_ndx) = sflx(i,c10h16_ndx) / (.7_r8 + .3_r8*dayfrac)
                if( loc_angle >= sunoff .and. loc_angle <= sunon ) then
                   sflx(i,c10h16_ndx) = sflx(i,c10h16_ndx) * .7_r8
                endif
             end if
          end if
       end if

       !--------------------------------------------------------
       !	... adjust isoprene for diurnal variation
       !--------------------------------------------------------
       if( isop_ndx > 0 ) then
          if( has_emis(isop_ndx) ) then
             if( .not. polar_night ) then
                if( polar_day ) then
                   iso_off = .8_r8 * pi
                   iso_on  = 1.2_r8 * pi
                else
                   iso_off = .8_r8 * sunoff
                   iso_on  = 2._r8 * pi - iso_off
                end if
                if( loc_angle >= iso_off .and. loc_angle <= iso_on ) then
                   sflx(i,isop_ndx) = 0._r8
                else
                   factor = loc_angle - iso_on
                   if( factor <= 0._r8 ) then
                      factor = factor + 2._r8*pi
                   end if
                   factor = factor / (2._r8*iso_off + 1.e-6_r8)
                   sflx(i,isop_ndx) = sflx(i,isop_ndx) * 2._r8 / iso_off * pi * (sin(pi*factor))**2
                end if
             else
                sflx(i,isop_ndx) = 0._r8
             end if
          end if
       end if

    end do

  end subroutine set_srf_emissions

end module mo_srf_emissions
