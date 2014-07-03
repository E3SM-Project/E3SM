module dust_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! Original version: PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
use spmd_utils,       only: masterproc
use camsrfexch,       only: cam_in_t, cam_out_t    
use ppgrid,           only: pcols, pver,pverp
use physconst,        only: mwdry, mwh2o, gravit, rair
use constituents,     only: pcnst, cnst_add, cnst_name, cnst_get_ind
use aerodep_flx,      only: aerodep_flx_prescribed
use abortutils,       only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

#if (defined MODAL_AERO)
  integer, parameter:: dust_number = 2
  integer, parameter:: ndst =2
#else
  integer, parameter:: dust_number = 4
  integer, parameter:: ndst =4
#endif
  integer, parameter:: dst_src_nbr =3
  integer, parameter:: sz_nbr =200

  integer :: ncyear
  integer :: ix_dust = -1

#if (defined MODAL_AERO)
  integer, parameter :: ncnst=ndst+2                 ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
#if (defined MODAL_AERO_7MODE)
       cnst_names = (/'dst_a5', 'dst_a7', 'num_a5', 'num_a7'/)
#elif (defined MODAL_AERO_3MODE)
       cnst_names = (/'dst_a1', 'dst_a3', 'num_a1', 'num_a3'/)
#endif
#else
  integer, parameter :: ncnst=ndst                   ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
       cnst_names = (/'DST01', 'DST02', 'DST03', 'DST04'/)
#endif

  !
  ! Public interfaces
  !
  public dust_register_cnst                        ! register consituents
  public dust_implements_cnst                      ! returns true if consituent is implemented by this package
  public dust_init_cnst                            ! initialize mixing ratios if not read from initial file
  public dust_initialize                           ! initialize (history) variables
  public dust_wet_intr                             ! interface to wet deposition
  public dust_emis_intr                            ! interface to emission
  public dust_drydep_intr                          ! interface to tendency computation
  public dust_idx1                                 ! allow other parts of code to know where dust is
  public dust_names  
  public dust_number
  public dust_set_idx
  public dust_has_wet_dep


  character(len=8), dimension(ncnst) :: dust_names = cnst_names
#if (defined MODAL_AERO)
  logical :: dust_has_wet_dep(ncnst) = .true.
#else
  logical :: dust_has_wet_dep(ndst) = .true.
#endif

  real(r8) stk_crc(ndst) ![frc] Correction to Stokes settling velocity
  real(r8) dns_aer       ![kg m-3] Aerosol density
  real(r8) tmp1          !Factor in saltation computation (named as in Charlie's code)
  real(r8) ovr_src_snk_mss(dst_src_nbr,ndst)  
  real(r8) dmt_vwr(ndst) ![m] Mass-weighted mean diameter resolved
  real(r8), allocatable ::  soil_erodibility(:,:)     ! soil erodibility factor
  real(r8), allocatable ::  soil_erodibility_in(:,:)  ! temporary input array

#if (defined MODAL_AERO)
  integer, target   :: spc_ndx(ndst)
  integer, target   :: spc_num_ndx(ndst)
#if (defined MODAL_AERO_7MODE)
  integer, pointer  :: dst_a5_ndx, dst_a7_ndx
  integer, pointer  :: num_a5_ndx, num_a7_ndx
#elif (defined MODAL_AERO_3MODE)
  integer, pointer  :: dst_a1_ndx, dst_a3_ndx
  integer, pointer  :: num_a1_ndx, num_a3_ndx
#endif
#endif

! Namelist variables
real(r8)      :: dust_emis_fact = -1.e36_r8   ! tuning parameter for dust emissions
character(cl) :: soil_erod = 'soil_erod'   ! full pathname for soil erodibility dataset


!===============================================================================
contains
!===============================================================================

  subroutine dust_register_cnst( )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for all aerosols
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: P. J. Rasch
    ! 
    !-----------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    real(r8), parameter :: one  = 1._r8
    real(r8), parameter :: zero  = 0._r8
    real(r8) :: mwght = one

    !-----------------------------------------------------------------------

    ! Set names of variables undergoing evolution
    ! returns m as current index for tracer numbers
    call cnst_add(cnst_names(1), mwght, one, zero, m ) 

    ! and store the start index of dust species used elsewhere in model retrieved by dust_idx1
    call dust_set_idx(m)

    call cnst_add(cnst_names(2), mwght, one, zero, m )
    call cnst_add(cnst_names(3), mwght, one, zero, m )
    call cnst_add(cnst_names(4), mwght, one, zero, m )

    return

  end subroutine dust_register_cnst



  !=======================================================================
  function dust_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this 
    !          package
    ! 
    ! Author: T. Henderson
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: dust_implements_cnst        ! return value
    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    dust_implements_cnst = .false.
    if (ix_dust<1) return

    do m = 1, ncnst
       if (name == cnst_names(m)) then
          dust_implements_cnst = .true.
          return
       end if
    end do
  end function dust_implements_cnst


  !=======================================================================
  subroutine dust_init_cnst(name, q, gcid)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Set initial mass mixing ratios.
    !
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name         ! constituent name

    real(r8), intent(out) :: q(:,:)            !  mass mixing ratio
    integer, intent(in)   :: gcid(:)           ! global column id
    !-----------------------------------------------------------------------

    if ( name == cnst_names(1) ) then
       q = 0._r8
    else if ( name == cnst_names(2) ) then
       q = 0._r8
    else if ( name == cnst_names(3) ) then
       q = 0._r8
    else if ( name == cnst_names(4) ) then
       q = 0._r8
    end if

  end subroutine dust_init_cnst

!===============================================================================

function dust_idx1()
   integer dust_idx1
   dust_idx1 = ix_dust
end function dust_idx1

!===============================================================================

subroutine dust_set_idx(m)
   integer, intent(in) :: m
   ix_dust = m
end subroutine dust_set_idx

!===============================================================================

subroutine dust_initialize(dust_emis_fact_in, soil_erod_in)

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: initialize parameterization of dust chemistry
   !          (declare history variables)
   ! 
   !-----------------------------------------------------------------------

   use cam_history,      only: addfld, add_default, phys_decomp
   use phys_control,     only: phys_getopts
   use ioFileMod,        only: getfil
#if (defined MODAL_AERO)
   use mo_chem_utls,     only: get_spc_ndx
#endif
   use pio
   use cam_pio_utils,    only : cam_pio_openfile
   use phys_grid,        only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
   use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, interp_type
   use ppgrid,           only : begchunk, endchunk
   use mo_constants,     only : pi, d2r

   ! arguments
   real(r8),      intent(in) :: dust_emis_fact_in
   character(cl), intent(in) :: soil_erod_in


   ! local variables
   integer :: m
   integer :: ix             ! start index for dust species
   character(len=32)  :: dummy
   character(cl)      :: soil_erod_file

   integer :: did, vid, nlat, nlon
   type(file_desc_t) :: ncid

   type(interp_type) :: lon_wgts, lat_wgts
   real(r8) :: to_lats(pcols), to_lons(pcols)
   real(r8), allocatable :: dst_lons(:)
   real(r8), allocatable :: dst_lats(:)
   integer :: c, ncols, ierr
   real(r8), parameter :: zero=0._r8, twopi=2._r8*pi
   logical  :: history_aerosol      ! Output the MAM aerosol tendencies
   logical  :: history_amwg         ! Output variables for AMWG diag

   !-----------------------------------------------------------------------

   call phys_getopts( history_aerosol_out     = history_aerosol, &
                      history_amwg_out        = history_amwg  )

   ! set module data from namelist vars read in aerosol_intr module
   dust_emis_fact = dust_emis_fact_in
   soil_erod      = soil_erod_in

#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
   dst_a5_ndx => spc_ndx(1)
   dst_a7_ndx => spc_ndx(2)
   num_a5_ndx => spc_num_ndx(1)
   num_a7_ndx => spc_num_ndx(2)

   dst_a5_ndx = get_spc_ndx( 'dst_a5' )
   dst_a7_ndx = get_spc_ndx( 'dst_a7' )
   num_a5_ndx = get_spc_ndx( 'num_a5' )
   num_a7_ndx = get_spc_ndx( 'num_a7' )
#elif (defined MODAL_AERO_3MODE)
   dst_a1_ndx => spc_ndx(1)
   dst_a3_ndx => spc_ndx(2)
   num_a1_ndx => spc_num_ndx(1)
   num_a3_ndx => spc_num_ndx(2)

   dst_a1_ndx = get_spc_ndx( 'dst_a1' )
   dst_a3_ndx = get_spc_ndx( 'dst_a3' )
   num_a1_ndx = get_spc_ndx( 'num_a1' )
   num_a3_ndx = get_spc_ndx( 'num_a3' )
#endif
#endif

   ! use Sam's dust initialize subroutine:  call equivalent here:
   call Dustini()

   ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
   ! read in soil erodibility factors, similar to Zender's boundary conditions

   ! Get file name.  
   call getfil(soil_erod, soil_erod_file, 0)
   call cam_pio_openfile (ncid, trim(soil_erod_file), PIO_NOWRITE)

   ! Get input data resolution.
   ierr = pio_inq_dimid( ncid, 'lon', did )
   ierr = pio_inq_dimlen( ncid, did, nlon )

   ierr = pio_inq_dimid( ncid, 'lat', did )
   ierr = pio_inq_dimlen( ncid, did, nlat )

   allocate(dst_lons(nlon))
   allocate(dst_lats(nlat))
   allocate(soil_erodibility_in(nlon,nlat))

   ierr = pio_inq_varid( ncid, 'lon', vid )
   ierr = pio_get_var( ncid, vid, dst_lons  )

   ierr = pio_inq_varid( ncid, 'lat', vid )
   ierr = pio_get_var( ncid, vid, dst_lats  )

   ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
   ierr = pio_get_var( ncid, vid, soil_erodibility_in )

   ! start index for dust species in constituent array
   ix = dust_idx1()

   ! Set names of variable tendencies and declare them as history variables
   call addfld ('LND_MBL','frac ',1, 'A','Soil erodibility factor',phys_decomp)
   call addfld ('RAM1','frac ',1, 'A','RAM1',phys_decomp)
   call addfld ('airFV','frac ',1, 'A','FV',phys_decomp)
   call addfld ('ORO','frac ',1, 'A','ORO',phys_decomp)

   if ( history_aerosol ) then  
      call add_default ('LND_MBL', 1, ' ')
      call add_default ('RAM1', 1, ' ')
      call add_default ('airFV', 1, ' ')
   endif
   if ( history_amwg ) then
     call add_default ('ORO', 1, ' ')
   endif 

#if (defined MODAL_AERO)

   do m = 1, ndst
      dummy = trim(cnst_name(ix-spc_ndx(1)+spc_ndx(m))) // 'SF'
      call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-spc_ndx(1)+spc_ndx(m)))//' dust surface emission',phys_decomp)
      if ( history_amwg ) then
         call add_default (dummy, 1, ' ')
      endif
   end do

#else

   do m = 1, 4
      dummy = trim(cnst_name(ix-1+m))
      call addfld (dummy,'kg/kg ',pver, 'A',trim(cnst_name(ix-1+m))//' dust mixing ratio',phys_decomp)
      call add_default(dummy, 1, ' ' )
      dummy = trim(cnst_name(ix-1+m)) // 'SF'
      call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' dust surface emission',phys_decomp)
      call add_default (dummy, 1, ' ')
      dummy = trim(cnst_name(ix-1+m)) // 'TB'
      call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' turbulent dry deposition flux',phys_decomp)
      call add_default (dummy, 1, ' ')
      dummy = trim(cnst_name(ix-1+m)) // 'GV'
      call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' gravitational dry deposition flux',phys_decomp)
      call add_default (dummy, 1, ' ')
      dummy = trim(cnst_name(ix-1+m)) // 'DD'
      call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' dry deposition flux at bottom (grav + turb)',phys_decomp)
      call add_default (dummy, 1, ' ')
      dummy = trim(cnst_name(ix-1+m)) // 'DT'
      call addfld (dummy,'kg/kg/s ',pver, 'A',trim(cnst_name(ix-1+m))//' dry deposition',phys_decomp)
      call add_default (dummy, 1, ' ')
      dummy = trim(cnst_name(ix-1+m)) // 'PP'
      call addfld (dummy,'kg/kg/s ',pver, 'A',trim(cnst_name(ix-1+m))//' wet deposition',phys_decomp)
      call add_default (dummy, 1, ' ')
      dummy = trim(cnst_name(ix-1+m)) // 'DV'
      call addfld (dummy,'m/s ',pver, 'A',trim(cnst_name(ix-1+m))//' deposition velocity',phys_decomp)
      call add_default (dummy, 1, ' ')
   end do

#endif

   call addfld ('DSTSFDRY','kg/m2/s',1, 'A','Dry deposition flux at surface',phys_decomp)
   call addfld ('DSTSFMBL','kg/m2/s',1, 'A','Mobilization flux at surface',phys_decomp)
   call addfld ('DSTSFWET','kg/m2/s',1, 'A','Wet deposition flux at surface',phys_decomp)
   call addfld ('DSTODXC','Tau ',1, 'A','Optical depth for diagnostics',phys_decomp)

   ! Determine default variables
   if (history_aerosol) then
      call add_default ('DSTSFDRY', 1, ' ')
      call add_default ('DSTSFMBL', 1, ' ')
      call add_default ('DSTSFWET', 1, ' ')
      call add_default ('DSTODXC', 1, ' ')
   endif

   ! Summary to log file
   if (masterproc) then
      write(iulog,*) 'Initializing prognostic dust (subroutine dust_initialize): '
      write(iulog,*) ' start index for dust is ', ix
      write(iulog,*) ' dust_emis_fact = ', dust_emis_fact
      write(iulog,*) ' soil erodibility dataset: ', trim(soil_erod)
   end if

!-----------------------------------------------------------------------
!     	... convert to radians and setup regridding
!-----------------------------------------------------------------------
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

   allocate( soil_erodibility(pcols,begchunk:endchunk) )

!-----------------------------------------------------------------------
!     	... regrid ..
!-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       
       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon,nlat , soil_erodibility(:,c), ncols, lon_wgts,lat_wgts)

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( soil_erodibility_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate soil_erod_in, ierr = ',ierr
       call endrun
    end if

end subroutine dust_initialize

!===============================================================================

  subroutine dust_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr, cam_out)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to wet processing of aerosols (source and sinks).
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: B.A. Boville
    ! 
    !-----------------------------------------------------------------------
    use cam_history,   only: outfld
    use physics_types, only: physics_state, physics_ptend
    use wetdep,        only: wetdepa_v1

    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: nstep
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                    ! latitude 
    real(r8), intent(in) :: cme(pcols,pver)            ! local condensation of cloud water
    real(r8), intent(in) :: prain(pcols,pver)            ! production of rain
    real(r8), intent(in) :: evapr(pcols,pver)            ! evaporation of rain
    real(r8), intent(in) :: cldn(pcols,pver)            ! cloud fraction
    real(r8), intent(in) :: cldc(pcols,pver)            ! convective cloud fraction
    real(r8), intent(in) :: cldv(pcols,pver)            ! cloudy volume undergoing scavenging
    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: rainmr(pcols, pver) ! rain mixing ratio

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    type(cam_out_t),     intent(inout) :: cam_out        ! export state

    real(r8), intent(inout) :: fracis(pcols,pver,pcnst)         ! fraction of transported species that are insoluble

    !
    ! Local variables
    !
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tstate(pcols, pver, 4)            ! temporary state vector
    real(r8) :: sflx(pcols)                       ! deposition flux per bin
    real(r8) :: sflx_tot(pcols)                   ! deposition flux (accumulates all bins)
    real(r8) :: obuf(1)
    real(r8), intent(in) :: calday        ! current calendar day
    real(r8) :: iscavt(pcols, pver)
    real(r8) :: scavt(pcols, pver)
    real(r8) :: scavcoef(pcols, pver)
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    integer :: mm,i,k
    integer :: m                                  ! tracer index
    integer :: ixcldliq
    integer :: ixcldice
    real(r8) totcond(pcols, pver) ! total condensate
    real(r8) :: sol_fact

    !-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    sflx_tot(:)=0._r8

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)
    scavcoef(:ncol,:) = 0.1_r8

    do m = 1, ndst

       sflx(:) = 0._r8
       if (dust_has_wet_dep(m) ) then
          mm = ix_dust + m - 1
          scavt=0._r8
          !       write(iulog,*) 'wet dep removed for debugging'
          ptend%lq(mm) = .TRUE.
          ptend%name  = trim(ptend%name)//',dust'
          sol_fact = 0.15_r8
          call wetdepa_v1( state%t, state%pmid, state%q, state%pdel,  &
               cldn, cldc, cmfdqr, conicw, prain, cme,                     &
               evapr, totcond, state%q(:,:,mm), dt,            &
               scavt, iscavt, cldv, fracis(:,:,mm), sol_fact, ncol, scavcoef )
          ptend%q(:ncol,:,mm)=scavt(:ncol,:)

          call outfld( trim(cnst_name(mm))//'PP', ptend%q(:,:,mm), pcols, lchnk)
          !      write(iulog,*) ' range of ptend ix ', minval(ptend%q(:ncol,:,mm)),maxval(ptend%q(:ncol,:,mm))
          do k = 1, pver
             do i = 1, ncol
                sflx(i) = sflx(i) + ptend%q(i,k,mm)*state%pdel(i,k)/gravit
             end do
          end do
          sflx_tot(:ncol) = sflx_tot(:ncol) + sflx(:ncol)

          ! if the user has specified prescribed aerosol dep fluxes then 
          ! do not set cam_out dep fluxes according to the prognostic aerosols
          if (.not. aerodep_flx_prescribed()) then
             ! set deposition in export state
             select case (m)
             case(1)
                do i = 1, ncol
                   cam_out%dstwet1(i) = max(-sflx(i), 0._r8)
                end do
             case(2)
                do i = 1, ncol
                   cam_out%dstwet2(i) = max(-sflx(i), 0._r8)
                end do
             case(3)
                do i = 1, ncol
                   cam_out%dstwet3(i) = max(-sflx(i), 0._r8)
                end do
             case(4)
                do i = 1, ncol
                   cam_out%dstwet4(i) = max(-sflx(i), 0._r8)
                end do
             end select
          endif

       endif
    end do
    call outfld( 'DSTSFWET', sflx_tot, pcols, lchnk)

    !   write(iulog,*) ' dust_wet_intr: pcols, pcols ', pcols, pcols

    return

  end subroutine dust_wet_intr

  subroutine dust_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac,fvin,ram1in, cam_out)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to dry deposition and sedimentation of dust
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: Natalie Mahowald and Phil Rasch
    ! 
    !-----------------------------------------------------------------------
    use cam_history,       only: outfld
    use physics_types,     only: physics_state, physics_ptend
    use phys_grid,         only: get_lat_all_p
    use constituents,      only: cnst_name
    use drydep_mod,        only: setdvel,  d3ddflux, calcram
    use dust_sediment_mod, only: dust_sediment_tend, dust_sediment_vel
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                 ! latitude 
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! obklen
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel--used over oceans and sea ice.
    real(r8), intent(in) :: ts(pcols)                     ! sfc temp
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: hflx(pcols)                  ! sensible heat flux
    real(r8), intent(in) :: prect(pcols)                     ! prect
    real(r8), intent(in) :: snowh(pcols)                     ! snow depth
    real(r8), intent(in) :: pblh(pcols)                     ! pbl height
    integer, intent(in)  :: month
    real(r8), intent(in) :: wvflx(pcols)       ! water vapor flux
    real(r8), intent(in) :: fvin(pcols)        ! for dry dep velocities from land model for dust
    real(r8), intent(in) :: ram1in(pcols)       ! for dry dep velocities from land model for dust

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    type(cam_out_t),     intent(inout) :: cam_out        ! export state

    ! Local variables
    !
    integer :: m                                  ! tracer index
    integer :: mm                                  ! tracer index
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tvs(pcols,pver)
    real(r8) :: dvel(pcols)            ! deposition velocity
    real(r8) :: sflx(pcols)            ! deposition flux
    real(r8) :: vlc_dry(pcols,pver,ndst)            ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,ndst)            ! dep velocity
    real(r8)::  vlc_trb(pcols,ndst)            ! dep velocity
    real(r8)::  dep_trb(pcols)       !kg/m2/s
    real(r8)::  dep_dry(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_dry_tend(pcols,pver)       !kg/kg/s (total of grav and trb)
    real(r8) :: obuf(1)
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real(r8)  pvdust(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: tsflx(pcols)
    real(r8) :: fv(pcols)         ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)       ! for dry dep velocities, from land modified over ocean & ice

    integer :: i,k
    real(r8) :: oro(pcols)
    !
    !-----------------------------------------------------------------------

    ix = dust_idx1()
    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol
    tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
    rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))
    tsflx(:)=0._r8
    ! calculate oro--need to run match

    do i=1,ncol
       oro(i)=1._r8
       if(icefrac(i)>0.5_r8) oro(i)=2._r8
       if(ocnfrac(i)>0.5_r8) oro(i)=0._r8
    enddo
    call outfld( 'ORO', oro, pcols, lchnk )

    !   write(iulog,*) ' dust drydep invoked '

    !   Dry deposition of Dust Aerosols
    !   #################################
    !    call setdvel( ncol, landfrac, icefrac, ocnfrac, .001_r8, .001_r8, .001_r8, dvel )
    !  we get the ram1,fv from the land model as ram1in,fvin, but need to calculate it over oceans and ice.  
    !  better if we got thse from the ocean and ice model
    !  for friction velocity, we use ustar (from taux and tauy), except over land, where we use fv from the land model.

    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
         ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
         state%pdel(:,pver),fvin,fv)
    call outfld( 'airFV', fv(:), pcols, lchnk )
    call outfld( 'RAM1', ram1(:), pcols, lchnk )
    !       call outfld( 'icefrc', icefrac(:), pcols, lchnk )

    call DustDryDep(ncol,state%t(:,:),state%pmid(:,:),ram1,fv,vlc_dry,vlc_trb,vlc_grv,landfrac)

    do m=1,ndst

       mm = dust_idx1() + m - 1

       ! use pvdust instead (means making the top level 0)
       pvdust(:ncol,1)=0._r8
       pvdust(:ncol,2:pverp) = vlc_dry(:ncol,:,m)


       call outfld( trim(cnst_name(mm))//'DV', pvdust(:,2:pverp), pcols, lchnk )
       if(.true.) then ! use phil's method
          !      convert from meters/sec to pascals/sec
          !      pvdust(:,1) is assumed zero, use density from layer above in conversion
          pvdust(:ncol,2:pverp) = pvdust(:ncol,2:pverp) * rho(:ncol,:)*gravit        

          !      calculate the tendencies and sfc fluxes from the above velocities
          call dust_sediment_tend( &
               ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
               state%q(:,:,mm) , pvdust  , ptend%q(:,:,mm), sflx  )
       else   !use charlie's method
          call d3ddflux(ncol, vlc_dry(:,:,m), state%q(:,:,mm),state%pmid,state%pdel, tvs,sflx,ptend%q(:,:,mm),dt)
       endif

       ! apportion dry deposition into turb and gravitational settling for tapes
       do i=1,ncol
          dep_trb(i)=sflx(i)*vlc_trb(i,m)/vlc_dry(i,pver,m)
          dep_grv(i)=sflx(i)*vlc_grv(i,pver,m)/vlc_dry(i,pver,m)
       enddo
       tsflx(:ncol)=tsflx(:ncol)+sflx(:ncol)

       ! if the user has specified prescribed aerosol dep fluxes then 
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not. aerodep_flx_prescribed()) then
          ! set deposition in export state
          select case (m)
          case(1)
             do i = 1, ncol
                cam_out%dstdry1(i) = max(sflx(i), 0._r8)
             end do
          case(2)
             do i = 1, ncol
                cam_out%dstdry2(i) = max(sflx(i), 0._r8)
             end do
          case(3)
             do i = 1, ncol
                cam_out%dstdry3(i) = max(sflx(i), 0._r8)
             end do
          case(4)
             do i = 1, ncol
                cam_out%dstdry4(i) = max(sflx(i), 0._r8)
             end do
          end select
       endif

       call outfld( trim(cnst_name(mm))//'DD',sflx, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'TB', dep_trb, pcols, lchnk )
       call outfld( trim(cnst_name(mm))//'GV', dep_grv, pcols, lchnk )

       call outfld( trim(cnst_name(mm))//'DT',ptend%q(:,:,mm), pcols, lchnk)
       !       write(iulog,*) ' range of tends for dust ', mm, minval(ptend%q(:ncol,:,mm)), maxval(ptend%q(:ncol,:,mm))
       !      call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk)


    end do
    ! output the total dry deposition
    call outfld( 'DSTSFDRY', tsflx, pcols, lchnk)


    ! set flags for tendencies (water and 4 ghg's)
    ptend%name  = ptend%name//'+dust_drydep'
    ptend%lq(ioff+1:ioff+4) = .TRUE.

    return
  end subroutine dust_drydep_intr

subroutine dust_emis_intr(state, cam_in)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model (need #define BGC) and
    ! the soil erodibility factor is applied here.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: Phil Rasch and Natalie Mahowald
    !
    ! 
    !-----------------------------------------------------------------------
    use cam_history,   only: outfld
    use physics_types, only: physics_state

#if (defined MODAL_AERO)
    use physconst,     only: pi,rair
    use modal_aero_data
#endif

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t), target, intent(inout) :: cam_in  ! import state

    ! Local variables

    integer :: lchnk
    integer :: ncol
    integer :: i
    integer :: m
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output
    real(r8) :: soil_erod_tmp(pcols)
    real(r8), pointer :: cflx(:,:)
    !
    real(r8) :: fact  ! tuning factor for dust emissions
       

#if (defined MODAL_AERO)
    integer :: mm,l
    real(r8) :: x_mton   ! (aero # emit)/(aero mass emit) ratio  (#/kg)
#endif 


    lchnk = state%lchnk
    ncol = state%ncol

    sflx(:)=0._r8

    cflx => cam_in%cflx

#if (defined MODAL_AERO)
! zero out the portion of cflx used for modal aerosol number because number is added up from sea salt and dst
! number and mass species
     do m = 1, ntot_amode
         l = numptr_amode(m)
         if (l > 0) cflx(:,l) = 0.0_r8
     end do
#endif

    !
    !         write(40,*) cflx(1,1)
    !         write(41,*) soil_erodibility(1,lat(1))
    ! we use charlie's sandblasting size distribution (instead of what Sam used
    ! in clm (and overwrite them here).
#if (defined MODAL_AERO)
    do i = 1, ncol
       soil_erod_tmp(i) = soil_erodibility( i,lchnk )
       if(soil_erod_tmp(i) .lt. 0.1_r8) soil_erod_tmp(i)=0._r8
    end do
    cflx(:ncol,dust_idx1())=cflx(:ncol,dust_idx1())
    cflx(:ncol,dust_idx1()+spc_ndx(2)-spc_ndx(1))=cflx(:ncol,dust_idx1()+spc_ndx(2)-spc_ndx(1))
    if(spc_num_ndx(1) > 0) then
      cflx(:ncol,dust_idx1()+spc_num_ndx(1)-spc_ndx(1))=cflx(:ncol,dust_idx1())
    endif
    if(spc_num_ndx(2) > 0) then
      cflx(:ncol,dust_idx1()+spc_num_ndx(2)-spc_ndx(1))=cflx(:ncol,dust_idx1()+spc_ndx(2)-spc_ndx(1))
    endif
#else
    do i = 1, ncol

       soil_erod_tmp(i) = soil_erodibility( i,lchnk )

       ! jfl
       ! change test to from 0.1 to 0.001 after
       ! discussion with N. Mahowald April 8 2009
       if(soil_erod_tmp(i) .lt. 0.001_r8) soil_erod_tmp(i)=0._r8
    end do
    cflx(:ncol,dust_idx1())=cflx(:ncol,dust_idx1())*0.038_r8/0.032456_r8
    cflx(:ncol,dust_idx1()+1)=cflx(:ncol,dust_idx1()+1)*0.11_r8/0.174216_r8
    cflx(:ncol,dust_idx1()+2)=cflx(:ncol,dust_idx1()+2)*0.17_r8/0.4085517_r8
    cflx(:ncol,dust_idx1()+3)=cflx(:ncol,dust_idx1()+3)*0.67_r8/0.384811_r8
#endif

    ! tuning factor (fact) for dust emissions
    fact = dust_emis_fact

#if (defined MODAL_AERO)
    do mm=1,ndst
       ! multiply by soil erodibility factor
       if(spc_num_ndx(mm) > 0) then
         m=dust_idx1()+spc_num_ndx(mm)-spc_ndx(1)
         x_mton = 6.0_r8 / (pi*dns_aer*(dmt_vwr(mm)**3.0_r8))
         cflx(:ncol,m)=cflx(:ncol,m)*x_mton*soil_erod_tmp(:ncol)/fact*1.15_r8
       endif
       m=dust_idx1()+spc_ndx(mm)-spc_ndx(1)
       cflx(:ncol,m)=cflx(:ncol,m)*soil_erod_tmp(:ncol)/fact*1.15_r8
       sflx(:ncol)=sflx(:ncol)+cflx(:ncol,m)
#else
    do m=dust_idx1(),dust_idx1()+3
       ! multiply by soil erodibility factor
       cflx(:ncol,m)=cflx(:ncol,m)*soil_erod_tmp(:ncol)/fact*1.15_r8
       sflx(:ncol)=sflx(:ncol)+cflx(:ncol,m)
#endif

       call outfld(trim(cnst_name(m)) // 'SF',cflx(:,m),pcols,lchnk)
       ! this is being done inside of the vertical diffusion automatically
       !         ptend%lq(m) = .true. ! tendencies for all dust on
       !         ptend%q(:ncol,pver,m) = cflx(:ncol,m)*gravit/state%pdel(:ncol,pver)
    enddo

    call outfld('LND_MBL',soil_erod_tmp,pcols,lchnk)
    call outfld('DSTSFMBL',sflx(:),pcols,lchnk)


    !      write(42,*) cflx(1,1)
    return
  end subroutine dust_emis_intr

  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: subroutine Dustini()
  !
  ! !INTERFACE:
  !
  subroutine Dustini()
    !
    ! !DESCRIPTION: 
    !
    ! Compute source efficiency factor from topography
    ! Initialize other variables used in subroutine Dust:
    ! ovr_src_snk_mss(m,n) and tmp1.
    ! Define particle diameter and density needed by atm model
    ! as well as by dry dep model
    ! Source: Paul Ginoux (for source efficiency factor)
    ! Modifications by C. Zender and later by S. Levis
    ! Rest of subroutine from C. Zender's dust model
    !
    ! !USES
    !
    use physconst,    only: pi,rair
    use shr_spfn_mod, only: erf => shr_spfn_erf
    !
    ! !ARGUMENTS:
    !
    implicit none
    !
    ! !REVISION HISTORY
    ! Created by Samual Levis
    ! Revised for CAM by Natalie Mahowald
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !Local Variables
    integer  :: ci,m,n                  !indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ![frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ![frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ![frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 !Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 !Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ![m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ![m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ![m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ![m] Width of size bin
    real(r8) :: slp_crc(ndst)           ![frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ![m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ![m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ![m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ![frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ![frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     !temporary 
    real(r8) :: ln_gsd                  ![frc] ln(gsd)
    real(r8) :: gsd_anl                 ![frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ![m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ![m] Number median particle diameter
    real(r8) :: lgn_dst                 !Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ![frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ![frc] Current relative accuracy
    real(r8) :: itr_idx                 ![idx] Counting index
    real(r8) :: dns_mdp                 ![kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ![m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ![kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ![m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            !Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      !Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ![m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ![m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ![m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ![m] Size Bin widths

    ! constants
    real(r8) :: dmt_vma_src(dst_src_nbr) =    &     ![m] Mass median diameter
         (/ 0.832e-6_r8 , 4.82e-6_r8 , 19.38e-6_r8 /)        !BSM96 p. 73 Table 2
    real(r8) :: gsd_anl_src(dst_src_nbr) =    &     ![frc] Geometric std deviation
         (/ 2.10_r8     ,  1.90_r8   , 1.60_r8     /)        !BSM96 p. 73 Table 2
    real(r8) :: mss_frc_src(dst_src_nbr) =    &     ![frc] Mass fraction 
         (/ 0.036_r8, 0.957_r8, 0.007_r8 /)                  !BSM96 p. 73 Table 2
#if (defined MODAL_AERO)
    real(r8) :: dmt_grd(3) =                  &     ![m] Particle diameter grid
#if (defined MODAL_AERO_7MODE)
         (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
#elif (defined MODAL_AERO_3MODE)
         (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
#endif
#else
    real(r8) :: dmt_grd(5) =                  &     ![m] Particle diameter grid
         (/ 0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.0e-6_r8 /)
#endif
    real(r8), parameter :: dmt_slt_opt = 75.0e-6_r8    ![m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0_r8         ![kg m-3] Density of optimal saltation particles


    ! declare erf intrinsic function
    real(r8) :: dum     !dummy variable for erf test
    !------------------------------------------------------------------------

    ! Sanity check on erf: erf() in SGI /usr/lib64/mips4/libftn.so is bogus

    dum = 1.0_r8
    if (abs(0.8427_r8-erf(dum))/0.8427_r8>0.001_r8) then
       write(iulog,*) 'erf(1.0) = ',erf(dum)
       call endrun ('Dustini: Error function error')
    endif
    dum = 0.0_r8
    if (erf(dum) /= 0.0_r8) then
       write(iulog,*) 'erf(0.0) = ',erf(dum)
       call endrun ('Dustini: Error function error')
    endif

    ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
    !                      and (2) dstszdst.F subroutine dst_szdst_ini
    ! purpose(1): given one set (the "source") of lognormal distributions,
    !             and one set of bin boundaries (the "sink"), compute and return
    !             the overlap factors between the source and sink distributions
    ! purpose(2): set important statistics of size distributions

    do m = 1, dst_src_nbr
       sqrt2lngsdi = sqrt(2.0_r8) * log(gsd_anl_src(m))
       do n = 1, ndst
          lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
          lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
          ovr_src_snk_frc = 0.5_r8 * (erf(lndmaxjovrdmdni/sqrt2lngsdi) - &
               erf(lndminjovrdmdni/sqrt2lngsdi))
          ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
       enddo
    enddo

    !delete portions that have to do only with the sink

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)

    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if dust_number ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here dust_number = dst_nbr

    ! Override automatic grid with preset grid if available
#if (defined MODAL_AERO)
    if (dust_number == 2) then
#else
    if (dust_number == 4) then
#endif
       do n = 1, dust_number
          dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
          dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
          dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
          dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
       end do
    else
#if (defined MODAL_AERO)
       call endrun('Dustini error: dust_number must equal to 2 with current code')
#else
       call endrun('Dustini error: dust_number must equal to 4 with current code')  
#endif
       !see more comments above)
    endif

    ! Bin physical properties
    gsd_anl = 2.0_r8      ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)
    dns_aer = 2.5e+3_r8   ! [kg m-3] Aerosol density
    ! Set a fundamental statistic for each bin
    dmt_vma = 2.524e-6_r8 ! [m] Mass median diameter analytic She84 p.75 Table1
    dmt_vma = 3.5e-6_r8
    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)
    dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]
    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl
    do n = 1, dust_number
       series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_r8/sz_nbr)
       sz_min(1) = dmt_min(n)
       do m = 2, sz_nbr                            ! Loop starts at 2
          sz_min(m) = sz_min(m-1) * series_ratio
       end do

       ! Derived grid values
       do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
          sz_max(m) = sz_min(m+1)                  ! [m]
       end do
       sz_max(sz_nbr) = dmt_max(n)                 ! [m]

       ! Final derived grid values
       do m = 1, sz_nbr
          sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
          sz_dlt(m) = sz_max(m)-sz_min(m)
       end do
       lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*pi))
       dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
       vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved
       do m = 1, sz_nbr
          ! Evaluate lognormal distribution for these sizes (call lgn_evl)
          tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
          lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)
          ! Integrate moments of size distribution
          dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               pi / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
          vlm_rsl(n) = vlm_rsl(n) +                                &
               pi / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
       end do
       dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved
    end do
    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)
    eps_max = 1.0e-4_r8
    dns_mdp = 100000._r8 / (295.0_r8*rair) ![kg m-3] const prs_mdp & tpt_vrt
    ! Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5_r8 * ((295.0_r8/273.0_r8)**1.5_r8) * 393.0_r8 / &
         (295.0_r8+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000._r8*sqrt(8.0_r8/(pi*rair*295.0_r8)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, dust_number
       slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
            (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
            dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
       vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
            gravit * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
    do m = 1, dust_number

       ! Initialize accuracy and counter
       eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
       itr_idx = 0                 ![idx] Counting index

       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(m) = vlc_stk(m)     ![m s-1]
       do while(eps_crr > eps_max)

          ! Save terminal velocity for convergence test
          vlc_grv_old = vlc_grv(m) ![m s-1]
          ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

          ! Update drag coefficient based on new Reynolds number
          if (ryn_nbr_grv(m) < 0.1_r8) then
             cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                  (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                  9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                  log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 500.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                  (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0e5_r8) then
             cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
          else
             write(iulog,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
             call endrun ('Dustini error: Reynolds number too large in stk_crc_get()')
          endif

          ! Update terminal velocity based on new Reynolds number and drag coeff
          ! [m s-1] Terminal veloc SeP97 p.467 (8.44)
          vlc_grv(m) = sqrt(4.0_r8 * gravit * dmt_vwr(m) * slp_crc(m) * dns_aer / &
               (3.0_r8*cff_drg_grv(m)*dns_mdp))   
          eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
             ! due to discontinuities in derivative of drag coefficient
             vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
          endif
          if (itr_idx > 20) then
             write(iulog,*) 'Dustini error: Terminal velocity not converging ',&
                  ' in stk_crc_get(), breaking loop...'
             goto 100                                        !to next iteration
          endif
          itr_idx = itr_idx + 1

       end do                                                !end while
100    continue   !Label to jump to when iteration does not converge
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities
    do m = 1, dust_number
       stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

    return
  end subroutine Dustini

  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: subroutine DustDryDep(c)
  !
  ! !INTERFACE:
  !
  subroutine DustDryDep(ncol,t,pmid,ram1,fv,vlc_dry,vlc_trb,vlc_grv,landfrac)
    !
    ! !DESCRIPTION: 
    !
    ! Determine Turbulent dry deposition for dust. Calculate the turbulent 
    ! component of dust dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the dust dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CCSM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, dustini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
    !
    ! !USES
    !
    use physconst,     only: rair,pi,boltz
    !
    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8) :: t(pcols,pver)       !atm temperature (K)
    real(r8) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8),intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8),intent(in) :: ram1(pcols)           
    !    real(r8) :: ramxo1(pcols)         !aerodynamical resistance (s/m)
    real(r8) :: vlc_trb(pcols,ndst)  !Turbulent deposn velocity (m/s)
    real(r8) :: vlc_grv(pcols,pver,ndst)  !grav deposn velocity (m/s)
    real(r8) :: vlc_dry(pcols,pver,ndst)  !dry deposn velocity (m/s)
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    integer, intent(in) :: ncol
    !
    ! !REVISION HISTORY
    ! Created by Sam Levis
    ! Modified for CAM by Natalie Mahowald
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k          !indices
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(pcols,pver,ndst) ![frc] Slip correction factor
    real(r8) :: rss_lmn(ndst) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp           !temporary 

    ! constants
    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over land
    ! needs fv and ram1 passed in from lnd model

    !------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol
          rho=pmid(i,k)/rair/t(i,k)
          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

          do m = 1, dust_number
             slp_crc(i,k,m) = 1.0_r8 + 2.0_r8 * mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm(i,k)))) / &
                  dmt_vwr(m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
                  gravit * slp_crc(i,k,m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
             vlc_dry(i,k,m)=vlc_grv(i,k,m)
          end do
       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       do m = 1, dust_number
          stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = boltz * t(i,k) * slp_crc(i,k,m) / &    ![m2 s-1]
               (3.0_r8*pi*vsc_dyn_atm(i,k)*dmt_vwr(m)) !SeP97 p.474
          shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
          !          shm_nbr_xpn=shm_nbr_xpn_ocn
          !          if(landfrac(i) .gt. 0.5_r8 ) shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          shm_nbr_xpn = shm_nbr_xpn_lnd 
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965
       end do

       ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)
       do m = 1, dust_number
          rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
          vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
          vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
       end do

    end do !end ncols loop

    return
  end subroutine DustDryDep

end module dust_intr
