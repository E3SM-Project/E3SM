module progseasalts_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! Original version: PJR (extensively modified from chemistry module)
! prognostic sea salt module taken from dust module by nmm 11/12/03
!---------------------------------------------------------------------------------

use shr_kind_mod,      only: r8 => shr_kind_r8
use spmd_utils,        only: masterproc
use ppgrid,            only: pcols, pver, pverp
use physconst,         only: pi, mwdry, mwh2o, gravit, rair
use camsrfexch,        only: cam_in_t
use constituents,      only: pcnst, cnst_add, cnst_name, cnst_get_ind
use drydep_mod,        only: calcram
use dust_sediment_mod, only: dust_sediment_tend, dust_sediment_vel
use cam_logfile,       only: iulog

implicit none
private          ! Make default type private to the module
save

public nsst  ! number of sea salt constituents
public ncnst  ! number of sea salt constituents

#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
  integer, parameter:: nsst =4
  real(r8), parameter :: scalefactor = 1.62_r8
#elif (defined MODAL_AERO_3MODE)
  integer, parameter:: nsst =3
  real(r8), parameter :: scalefactor = 1.35_r8 
#endif
#else
  integer, parameter:: nsst =4
#endif

  integer :: ncyear
  integer :: ix_progseasalts = -1
  integer :: nx_progseasalts

#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
  integer, parameter :: ncnst=nsst+4                 ! number of constituents
#elif (defined MODAL_AERO_3MODE)
  integer, parameter :: ncnst=nsst+3                 ! number of constituents
#endif
#else
  integer, parameter :: ncnst=nsst                   ! number of constituents
#endif

  character(len=8), dimension(ncnst), parameter :: & ! constituent names
#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
       cnst_names = (/'ncl_a1', 'ncl_a2', 'ncl_a4', 'ncl_a6', 'num_a1', 'num_a2', 'num_a4', 'num_a6'/)
#elif (defined MODAL_AERO_3MODE)
       cnst_names = (/'ncl_a1', 'ncl_a2', 'ncl_a3', 'num_a1', 'num_a2', 'num_a3'/)
#endif
#else
       cnst_names = (/'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04'/)
#endif

  !
  ! Public interfaces
  !
  public progseasalts_register_cnst                        ! register consituents
  public progseasalts_implements_cnst                      ! returns true if consituent is implemented by this package
  public progseasalts_init_cnst                            ! initialize mixing ratios if not read from initial file
  public progseasalts_initialize                           ! initialize (history) variables
  public progseasalts_wet_intr                             ! interface to wet deposition
  public progseasalts_emis_intr                            ! interface to emission
  public progseasalts_drydep_intr                          ! interface to tendency computation
!!$  public progseasalts_time_interp                          ! interpolate oxidants and fluxes to current time
  public progseasalts_idx1                                 ! allow other parts of code to know where progseasalts is
  public progseasalts_set_idx
  public progseasalts_names
  public progseasalts_has_wet_dep

  character(len=8), dimension(ncnst) :: progseasalts_names = cnst_names
  logical :: progseasalts_has_wet_dep(ncnst) = .true.

#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
  real(r8) :: stk_crc(nsst)=(/ 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8 /)
#elif (defined MODAL_AERO_3MODE)
  real(r8) :: stk_crc(nsst)=(/ 1.0_r8, 1.0_r8, 1.0_r8 /)
#endif
#else
  real(r8) :: stk_crc(nsst)=(/ 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8 /) 
#endif
  ![frc] Correction to Stokes settling velocity--we are assuming they are close enough 1.0 to be
  !       set to one--unlikely to cause considerable error
  real(r8),parameter:: dns_aer_sst = &
                       2200.0_r8     ![kg m-3] Aerosol density

!
!  real(r8)::  soil_erodibility(plon,plat)            ! soil erodibility factor
#if (defined MODAL_AERO)
  real(r8) :: sst_sz_range_lo(nsst) = &
#if (defined MODAL_AERO_7MODE)
              (/0.08e-6_r8, 0.02e-6_r8, 0.3e-6_r8, 1.0e-6_r8/)  ! accu, aitken, fine, coarse
#elif (defined MODAL_AERO_3MODE)
              (/0.08e-6_r8, 0.02e-6_r8, 1.0e-6_r8/)  ! accu, aitken, coarse
#endif
  real(r8) :: sst_sz_range_hi(nsst) = &
#if (defined MODAL_AERO_7MODE)
              (/0.3e-6_r8, 0.08e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
#elif (defined MODAL_AERO_3MODE)
              (/1.0e-6_r8, 0.08e-6_r8, 10.0e-6_r8/)
#endif
#else
  real(r8) :: sst_sz_range(nsst+1) = &
              (/0.2e-6_r8,1.0e-6_r8,3.0e-6_r8,10.0e-6_r8,20.0e-6_r8/)  ! bottom and tops of size bins for sea salts in m 
#endif
  !(diameter not radius as in tie)
  real(r8) :: sst_source(nsst)  ! source factors for each bin 
#if (defined MODAL_AERO)
  real(r8) :: sst_source_num(nsst)
  integer, target   :: spc_ndx(nsst)
  integer, target   :: spc_num_ndx(nsst)
#if (defined MODAL_AERO_7MODE)
  integer, pointer  :: ncl_a1_ndx, ncl_a2_ndx, ncl_a4_ndx, ncl_a6_ndx
  integer, pointer  :: num_a1_ndx, num_a2_ndx, num_a4_ndx, num_a6_ndx
#elif (defined MODAL_AERO_3MODE)
  integer, pointer  :: ncl_a1_ndx, ncl_a2_ndx, ncl_a3_ndx
  integer, pointer  :: num_a1_ndx, num_a2_ndx, num_a3_ndx
#endif
#endif
  real(r8) :: smt_vwr(nsst) & 
#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
              =(/0.051e-6_r8, 0.15e-6_r8, 0.56e-6_r8, 3.6e-6_r8/) ![m] Mass-weighted mean diameter resolved
#elif (defined MODAL_AERO_3MODE)
              =(/0.051e-6_r8, 0.256e-6_r8, 3.6e-6_r8/) ![m] Mass-weighted mean diameter resolved
#endif
! assume stddev=2.0 for seasalt, Dp is calculated from Dpg which is based on observations
! (Seinfeld and Pandis, Smirnov et al. 2003 etc). Dpg is assumed to be 0.05,0.2,3.0,10.0
#else
              =(/0.52e-6_r8,2.38e-6_r8,4.86e-6_r8,15.14e-6_r8/) ![m] Mass-weighted mean diameter resolved
#endif
  !initializaed in module and used in emissions

#if (defined MODAL_AERO)
! here is used for Ekman's ss
  integer :: sections ! number of sections used to calculate fluxes
  parameter (sections=(31))

! only use up to ~20um
  real(r8), dimension(sections) :: Dg = (/  &
                                0.0020e-5_r8, 0.0025e-5_r8, 0.0032e-5_r8,  &
                                0.0040e-5_r8, 0.0051e-5_r8, 0.0065e-5_r8,  &
                                0.0082e-5_r8, 0.0104e-5_r8, 0.0132e-5_r8,  &
                                0.0167e-5_r8, 0.0211e-5_r8, 0.0267e-5_r8,  &
                                0.0338e-5_r8, 0.0428e-5_r8, 0.0541e-5_r8,  &
                                0.0685e-5_r8, 0.0867e-5_r8, 0.1098e-5_r8,  &
                                0.1389e-5_r8, 0.1759e-5_r8, 0.2226e-5_r8,  &
                                0.2818e-5_r8, 0.3571e-5_r8, 0.4526e-5_r8,  &
                                0.5735e-5_r8, 0.7267e-5_r8, 0.9208e-5_r8,  &
                                1.1668e-5_r8, 1.4786e-5_r8, 1.8736e-5_r8,  &
                                2.3742e-5_r8 /)

  real(r8), dimension(sections) :: bm, rdry, rm
  real(r8), dimension(4,sections) :: consta, constb  !constants for calculating emission polynomial
#endif
  ! size ranges for sea salts
contains


  !===============================================================================
  subroutine progseasalts_register_cnst
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
    !-----------------------------------------------------------------------

    ! Set names of variables undergoing evolution
    ! returns m as current index for tracer numbers
    call cnst_add(cnst_names(1), one, one, zero, m) 

    ! and store the start index of progseasalts species used elsewhere in model retrieved by progseasalts_idx1
    call progseasalts_set_idx(m)

    call cnst_add(cnst_names(2), one, one, zero, m)
    call cnst_add(cnst_names(3), one, one, zero, m)
    call cnst_add(cnst_names(4), one, one, zero, m)
    return
  end subroutine progseasalts_register_cnst



  !=======================================================================
  function progseasalts_implements_cnst(name)
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
    logical :: progseasalts_implements_cnst        ! return value
    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    progseasalts_implements_cnst = .false.
    if ( ix_progseasalts < 1 ) return

    do m = 1, ncnst
       if (name == cnst_names(m)) then
          progseasalts_implements_cnst = .true.
          return
       end if
    end do

  end function progseasalts_implements_cnst


  !=======================================================================
  subroutine progseasalts_init_cnst(name, q, gcid)
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
    integer,  intent(in)  :: gcid(:)           ! global column id
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

  end subroutine progseasalts_init_cnst



  function progseasalts_idx1()
    implicit none
    integer progseasalts_idx1
    progseasalts_idx1 = ix_progseasalts
  end function progseasalts_idx1

  subroutine progseasalts_set_idx(m)
    implicit none
    integer m

    ix_progseasalts = m
  end subroutine progseasalts_set_idx

  !===============================================================================
  subroutine progseasalts_initialize 
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize parameterization of progseasalts chemistry
    !          (declare history variables)
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: NCAR CMS
    !  NMM, november 14,2003
    ! 
    !-----------------------------------------------------------------------
    use cam_history,   only: addfld, add_default, phys_decomp
    use physconst,     only: pi,rair,boltz
    use phys_control,  only: phys_getopts
#if (defined MODAL_AERO)
    use mo_chem_utls,  only: get_spc_ndx
#endif
    implicit none

    !---------------------------Local workspace-----------------------------
    integer :: m,n

    integer :: numsplit,ix,mm
    parameter (numsplit=50)
    real (r8):: b,rmid,dr  ! work variables for entrainment parameters
    character dummy*20
    logical  :: history_aerosol      ! Output the MAM aerosol tendencies
#if (defined MODAL_AERO)
    real(r8)    :: d0, d1, drydens, seasalt_emfac_num, seasalt_emfac_mas
#endif

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol   )

#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
    ncl_a1_ndx => spc_ndx(1)
    ncl_a2_ndx => spc_ndx(2)
    ncl_a4_ndx => spc_ndx(3)
    ncl_a6_ndx => spc_ndx(4)
    num_a1_ndx => spc_num_ndx(1)
    num_a2_ndx => spc_num_ndx(2)
    num_a4_ndx => spc_num_ndx(3)
    num_a6_ndx => spc_num_ndx(4)

    ncl_a1_ndx = get_spc_ndx( 'ncl_a1' )
    ncl_a2_ndx = get_spc_ndx( 'ncl_a2' )
    ncl_a4_ndx = get_spc_ndx( 'ncl_a4' )
    ncl_a6_ndx = get_spc_ndx( 'ncl_a6' )
    num_a1_ndx = get_spc_ndx( 'num_a1' )
    num_a2_ndx = get_spc_ndx( 'num_a2' )
    num_a4_ndx = get_spc_ndx( 'num_a4' )
    num_a6_ndx = get_spc_ndx( 'num_a6' )
#elif (defined MODAL_AERO_3MODE)
    ncl_a1_ndx => spc_ndx(1)
    ncl_a2_ndx => spc_ndx(2)
    ncl_a3_ndx => spc_ndx(3)
    num_a1_ndx => spc_num_ndx(1)
    num_a2_ndx => spc_num_ndx(2)
    num_a3_ndx => spc_num_ndx(3)

    ncl_a1_ndx = get_spc_ndx( 'ncl_a1' )
    ncl_a2_ndx = get_spc_ndx( 'ncl_a2' )
    ncl_a3_ndx = get_spc_ndx( 'ncl_a3' )
    num_a1_ndx = get_spc_ndx( 'num_a1' )
    num_a2_ndx = get_spc_ndx( 'num_a2' )
    num_a3_ndx = get_spc_ndx( 'num_a3' )
#endif
#endif

    ! initialize the emission bins:
    sst_source(:)=0._r8
#if (defined MODAL_AERO)
    sst_source_num(:)=0._r8
#endif
    ! sea salt source apprortionment calculated elsewhere.  Based on
    ! gong et al., with andreas 1998 refinement--see Mahowald et al., 2005 seasalt paper for descriptions

#if (defined MODAL_AERO)
!   do n = 1,nsst
!      d0 = sst_sz_range(n)   * 1.e+2_r8              ! diameter in cm
!      d1 = sst_sz_range(n+1) * 1.e+2_r8              ! diameter in cm
!      drydens = dns_aer_sst * 1.e-3_r8   ! g/cm3
!      call seasalt_emitfactors_1bin( 1, d0, d1, drydens, seasalt_emfac_num, seasalt_emfac_mas )
!      sst_source_num(n) = seasalt_emfac_num                ! #/m2/s
!      sst_source(n)     = seasalt_emfac_mas * 1.0e-3_r8    ! convert from g/m2/s to kg/m2/s
!   end do
! use Ekman's ss
    rdry(:)=Dg(:)/2._r8   ! meter
! multiply rm with 1.814 because it should be RH=80% and not dry particles
! for the parameterization
    rm(:)=1.814_r8*rdry(:)*1.e6_r8   ! um
    bm(:)=(0.380_r8-log10(rm(:)))/0.65_r8  ! use in Manahan

! calculate constants form emission polynomials
    do m=1,sections
       if ((m).le.9)then
         consta(1,m) = (-2.576_r8)*10._r8**35*Dg(m)**4+5.932_r8*10._r8**28  &
                   * Dg(m)**3+(-2.867_r8)*10._r8**21*Dg(m)**2+(-3.003_r8)  &
                   * 10._r8**13*Dg(m) + (-2.881_r8)*10._r8**6
         constb(1,m) = 7.188_r8*10._r8**37  &
                   * Dg(m)**4+(-1.616_r8)*10._r8**31*Dg(m)**3+6.791_r8*10._r8**23  &
                   * Dg(m)**2+1.829_r8*10._r8**16*Dg(m)+7.609_r8*10._r8**8
       elseif ((m).ge.10.and.(m).le.13)then
         consta(2,m) = (-2.452_r8)*10._r8**33*Dg(m)**4+2.404_r8*10._r8**27  &
                   * Dg(m)**3+(-8.148_r8)*10._r8**20*Dg(m)**2+(1.183_r8)*10._r8**14  &
                   * Dg(m)+(-6.743_r8)*10._r8**6
         constb(2,m) = 7.368_r8*10._r8**35  &
                   * Dg(m)**4+(-7.310_r8)*10._r8**29*Dg(m)**3+ 2.528_r8*10._r8**23  &
                   * Dg(m)**2+(-3.787_r8)*10._r8**16*Dg(m)+ 2.279_r8*10._r8**9
       elseif ((m).ge.14.and.(m).lt.22)then
         consta(3,m) = (1.085_r8)*10._r8**29*Dg(m)**4+(-9.841_r8)*10._r8**23  &
                   * Dg(m)**3+(3.132_r8)*10._r8**18*Dg(m)**2+(-4.165_r8)*10._r8**12  &
                   * Dg(m)+(2.181_r8)*10._r8**6
         constb(3,m) = (-2.859_r8)*10._r8**31  &
                   * Dg(m)**4+(2.601_r8)*10._r8**26*Dg(m)**3+(-8.297_r8)*10._r8**20  &
                   * Dg(m)**2+(1.105_r8)*10._r8**15*Dg(m)+(-5.800_r8)*10._r8**8
       elseif (m.ge.22.and.m.le.40)then
! use monahan
         consta(4,m) = (1.373_r8*rm(m)**(-3)*(1+0.057_r8*rm(m)**1.05_r8)  &
                   * 10**(1.19_r8*exp(-bm(m)**2)))  &
                   * (rm(m)-rm(m-1))
       endif
    enddo
#else
    ! gives source in units of kg/m2/s /(m/s)**3.41
    sst_source(1)=4.77e-15_r8
    sst_source(2)=5.19e-14_r8
    sst_source(3)=1.22e-13_r8
    sst_source(4)=6.91e-14_r8
#endif


    ix=ix_progseasalts

    do m = 1, nsst
#if (defined MODAL_AERO)
       dummy = trim(cnst_name(ix-spc_ndx(1)+spc_ndx(m))) // 'SF'
       call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-spc_ndx(1)+spc_ndx(m)))//' progseasalts surface emission',phys_decomp)
       if ( history_aerosol ) then
          call add_default (dummy, 1, ' ')
       endif   
#else
       dummy = trim(cnst_name(ix-1+m))
       call addfld (dummy,'kg/kg ',pver, 'A',trim(cnst_name(ix-1+m))//' progseasalts mixing ratio',phys_decomp)
       call add_default (dummy, 1, ' ')
       dummy = trim(cnst_name(ix-1+m)) // 'SF'
       call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' progseasalts surface emission',phys_decomp)
       call add_default (dummy, 1, ' ')
       dummy = trim(cnst_name(ix-1+m)) // 'OD'
       call addfld (dummy,'Tau ',1, 'A',trim(cnst_name(ix-1+m))//' optical depth',phys_decomp)
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
#if defined (  debugseasalts)
       dummy = trim(cnst_name(ix-1+m)) // 'DI'
       call addfld (dummy,'m/s ',pver, 'A',trim(cnst_name(ix-1+m))//' deposition diameter',phys_decomp)
       call add_default (dummy, 1, ' ')
#endif
#endif
    end do

    call addfld ('SSTSFDRY','kg/m2/s',1, 'A','Dry deposition flux at surface',phys_decomp)
    call addfld ('SSTSFMBL','kg/m2/s',1, 'A','Mobilization flux at surface',phys_decomp)
    call addfld ('SSTSFWET','kg/m2/s',1, 'A','Wet deposition flux at surface',phys_decomp)
    call addfld ('SSTODXC','Tau ',1, 'A','Optical depth for diagnostics',phys_decomp)
    if (history_aerosol) then
        call add_default ('SSTSFDRY', 1, ' ')
        call add_default ('SSTSFMBL', 1, ' ')
        call add_default ('SSTSFWET', 1, ' ')
        call add_default ('SSTODXC', 1, ' ')
    endif

#if defined ( debugseasalts)
    dummy = 'RH'
    call addfld (dummy,'frac',pver, 'A','RH in dry dep calc',phys_decomp)
    call add_default (dummy, 1, ' ')
#endif

    if (masterproc) then
       write(iulog,*) 'Initializing Prog seasalts'
       write(iulog,*) '  start index for progseasalts is ', ix_progseasalts
       write(iulog,*) '  sst_source',sst_source
    end if

  end subroutine progseasalts_initialize


  !===============================================================================
  subroutine progseasalts_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr )
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

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies


    real(r8), intent(inout) :: fracis(pcols,pver,pcnst)         ! fraction of transported species that are insoluble

    !
    ! Local variables
    !
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tstate(pcols, pver, 4)            ! temporary state vector
    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: rainmr(pcols, pver) ! rain mixing ratio
    real(r8) :: sflx(pcols)            ! deposition flux
    real(r8) :: obuf(1)
    real(r8) :: calday        ! current calendar day
    real(r8) :: iscavt(pcols, pver)
    real(r8) :: scavt(pcols, pver)
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    integer :: mm,i,k
    integer :: m                                  ! tracer index
    integer :: ixcldliq
    integer :: ixcldice
    real(r8) totcond(pcols, pver) ! total condensate
    real(r8) :: sol_fact
    real(r8) scavcoef(pcols,pver) ! ! Dana and Hales coefficient (/mm) (0.1)

    !-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    sflx(:)=0._r8
    scavcoef(:ncol,:) = 0.1_r8

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)
    
    do m = 1, nsst

       if (progseasalts_has_wet_dep(m) ) then
          mm = ix_progseasalts + m - 1
          scavt=0._r8
          !       write(iulog,*) 'wet dep removed for debugging'
          ptend%lq(mm) = .TRUE.
          ptend%name  = trim(ptend%name)//',sslt'
          sol_fact = 0.3_r8
          call wetdepa_v1( state%t, state%pmid, state%q, state%pdel,  &
               cldn, cldc, cmfdqr, conicw, prain, cme,                     &
               evapr, totcond, state%q(:,:,mm), dt,            &
               scavt, iscavt, cldv, fracis(:,:,mm), sol_fact, ncol, &
               scavcoef )
          ptend%q(:ncol,:,mm)=scavt(:ncol,:)

          call outfld( trim(cnst_name(mm))//'PP', ptend%q(:,:,mm), pcols, lchnk)
          !      write(iulog,*) ' range of ptend ix ', minval(ptend%q(:ncol,:,mm)),maxval(ptend%q(:ncol,:,mm))
          do k=1,pver
             do i=1,ncol
                sflx(i)=sflx(i)+ptend%q(i,k,mm)*state%pdel(i,k)/gravit
             enddo
          enddo
       endif
    end do
    call outfld( 'SSTSFWET', sflx, pcols, lchnk)

    !   write(iulog,*) ' progseasalts_wet_intr: pcols, pcols ', pcols, pcols

    return

  end subroutine progseasalts_wet_intr

  subroutine progseasalts_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac,fvin,ram1in)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to dry deposition and sedimentation of progseasalts
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
    use drydep_mod,        only: d3ddflux  
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
    real(r8), intent(in) :: fvin(pcols)        ! for dry dep velocities from land model for progseasalts
    real(r8), intent(in) :: ram1in(pcols)       ! for dry dep velocities from land model for progseasalts

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    !
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
    real(r8) :: vlc_dry(pcols,pver,nsst)            ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,nsst)            ! dep velocity
    real(r8)::  vlc_trb(pcols,nsst)            ! dep velocity
    real(r8)::  dep_trb(pcols)       !kg/m2/s
    real(r8)::  dep_dry(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_dry_tend(pcols,pver)       !kg/kg/s (total of grav and trb)
    real(r8) :: obuf(1)
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real(r8)  pvprogseasalts(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: tsflx(pcols)
    real(r8) :: fv(pcols)         ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)       ! for dry dep velocities, from land modified over ocean & ice

    integer :: i,k
    real(r8) :: oro(pcols)
    !
    !-----------------------------------------------------------------------

    ix = progseasalts_idx1()

    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol
    tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
    rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))
    ! calculate oro--need to run match

!!$    do i=1,ncol
!!$       oro(i)=1._r8
!!$       if(icefrac(i)>0.5) oro(i)=2._r8
!!$       if(ocnfrac(i)>0.5) oro(i)=0._r8
!!$    enddo
!!$    call outfld( 'ORO', oro, pcols, lchnk )

    !   write(iulog,*) ' progseasalts drydep invoked '

    !   Dry deposition of Progseasalts Aerosols
    !   #################################
    !    call setdvel( ncol, landfrac, icefrac, ocnfrac, .001_r8, .001_r8, .001_r8, dvel )
    !  we get the ram1,fv from the land model as ram1in,fvin,, but need to calculate it over oceans and ice.  
    !  better if we got thse from the ocean and ice model
    !  for friction velocity, we use ustar (from taux and tauy), except over land, where we use fv from the land model.

    ! copy fv,ram1 values from land model to variables to modify in calcram
    !    fv = fvin    
    !    ram1 = ram1in
    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
         ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
         state%pdel(:,pver),fvin,fv)
    ! this is the same as in the dust model--perhaps there is someway to 
    !calculate them up higher in the model
!!$    call outfld( 'airFV', fv(:), pcols, lchnk )
!!$    call outfld( 'RAM1', ram1(:), pcols, lchnk )
    !       call outfld( 'icefrc', icefrac(:), pcols, lchnk )

    call ProgseasaltsDryDep(ncol,state%t(:,:),state%pmid(:,:),ram1,fv,vlc_dry,vlc_trb,vlc_grv,state)
    tsflx(:)=0._r8
    do m=1,nsst

       mm = progseasalts_idx1() + m - 1
       ! use pvprogseasalts instead (means making the top level 0)
       pvprogseasalts(:ncol,1)=0._r8
       pvprogseasalts(:ncol,2:pverp) = vlc_dry(:ncol,:,m)


       call outfld( trim(cnst_name(mm))//'DV', pvprogseasalts(:,2:pverp), pcols, lchnk )
       if(.true.) then ! use phil's method
          !      convert from meters/sec to pascals/sec
          !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
          pvprogseasalts(:ncol,2:pverp) = pvprogseasalts(:ncol,2:pverp) * rho(:ncol,:)*gravit        

          !      calculate the tendencies and sfc fluxes from the above velocities
          call dust_sediment_tend( &
               ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
               state%q(:,:,mm) , pvprogseasalts  , ptend%q(:,:,mm), sflx  )
       else   !use charlie's method
          call d3ddflux(ncol, vlc_dry(:,:,m), state%q(:,:,mm),state%pmid,state%pdel, tvs,sflx,ptend%q(:,:,mm),dt)
       endif
       ! apportion dry deposition into turb and gravitational settling for tapes
       do i=1,ncol
          dep_trb(i)=sflx(i)*vlc_trb(i,m)/vlc_dry(i,pver,m)
          dep_grv(i)=sflx(i)*vlc_grv(i,pver,m)/vlc_dry(i,pver,m)
       enddo
       tsflx(:ncol)=tsflx(:ncol)+sflx(:ncol)

       call outfld( trim(cnst_name(mm))//'DD',sflx, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'TB', dep_trb, pcols, lchnk )
       call outfld( trim(cnst_name(mm))//'GV', dep_grv, pcols, lchnk )

       call outfld( trim(cnst_name(mm))//'DT',ptend%q(:,:,mm), pcols, lchnk)
       !       write(iulog,*) ' range of tends for progseasalts ', mm, minval(ptend%q(:ncol,:,mm)), maxval(ptend%q(:ncol,:,mm))
       !      call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk)


    end do
    ! output the total dry deposition
    call outfld( 'SSTSFDRY', tsflx, pcols, lchnk)


    ! set flags for tendencies (water and 4 ghg's)
    ptend%name  = ptend%name//'+progseasalts_drydep'
    ptend%lq(ioff+1:ioff+4) = .TRUE.

    return
  end subroutine progseasalts_drydep_intr

!!$  subroutine progseasalts_time_interp
!!$
!!$    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
!!$                            is_perpetual
!!$    !   use progseasaltsbnd, only: progseasaltsbndint
!!$
!!$    implicit none
!!$
!!$    real(r8) calday
!!$    integer :: yr, mon, day, ncsec
!!$
!!$    if ( ix_progseasalts < 1 ) return
!!$
!!$    calday = get_curr_calday()
!!$    if ( is_perpetual() ) then
!!$       call get_perp_date(yr, mon, day, ncsec)
!!$    else
!!$       call get_curr_date(yr, mon, day, ncsec)
!!$    end if
!!$
!!$    write(iulog,*) ' progseasalts_time_interp: interpolating progseasalts emissions ', calday
!!$    !    call progseasaltsbndint(calday)      ! interpolate oxidants
!!$
!!$
!!$  end subroutine progseasalts_time_interp

subroutine progseasalts_emis_intr(state, cam_in)

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Interface to emission of all progseasaltss.
   ! Derives from Xue Xie Tie's seasalts in MOZART, which derives from 
   !   Chin et al., 2001 (JAS) and Gong et al., 1997 (JGR-102,3805-3818).
   ! 
   ! Author: Phil Rasch and Natalie Mahowald
   !
   ! Derives from Martensson et al. (2003) (JGR-108, 4297,doi:10.1029/2002JD002263)
   ! valid from 20nm to ~2500nm dry diameter (based on lab experiment with artificial sea water)
   !
   ! currently we recommend that it is combined with
   ! the parameterisation by Monahan et al. (1986) for
   ! dry diameters > 2-3 um even if it lacks
   ! temperature dependence (despite that Bowyer et
   ! al. (1990) found a similar dependency in the tank
   ! from Monahan et al. (1986))
   !
   !-----------------------------------------------------------------------

   use cam_history,   only: outfld
   use physics_types, only: physics_state
   use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
   use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                      	     is_perpetual
   ! Arguments:
   type(physics_state),    intent(in )   :: state   ! Physics state variables
   type(cam_in_t), target, intent(inout) :: cam_in  ! import state

   ! Local variables
   real(r8), pointer :: cflx(:,:)
   real(r8), pointer :: ocnfrc(:)
   real(r8), pointer :: sst(:)    ! sea surface temp for Ekman's ss

   integer  lat(pcols)                  ! latitude index 
   integer  lon(pcols)                  ! longitude index
   integer lchnk
   integer ncol
   integer i
   integer m,mm
   real(r8):: sflx(pcols)   ! accumulate over all bins for output
   real(r8) :: soil_erod_tmp(pcols)
   real(r8) :: u10cubed(pcols)
   !    real (r8), parameter :: z0=0.5  ! m roughness length over oceans--from Tie.
   real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model
   !
#if (defined MODAL_AERO)
   ! use Ekman's ss
   real (r8) :: fi(pcols,sections)
   real (r8) :: W(pcols)
#endif

   real(r8) :: calday        ! current calendar day
   integer :: yr, mon, day, ncsec
   integer :: ncdate
   !-----------------------------------------------------------------------

   cflx   => cam_in%cflx
   ocnfrc => cam_in%ocnfrac
   sst    => cam_in%sst

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

    lchnk = state%lchnk
    ncol = state%ncol

    call get_lat_all_p(lchnk, ncol, lat)
    call get_lon_all_p(lchnk, ncol, lon)
    sflx(:)=0._r8
    u10cubed(:ncol)=sqrt(state%u(:ncol,pver)**2+state%v(:ncol,pver)**2)
    ! move the winds to 10m high from the midpoint of the gridbox:
    ! follows Tie and Seinfeld and Pandis, p.859 with math.

    u10cubed(:ncol)=u10cubed(:ncol)*log(10._r8/z0)/log(state%zm(:ncol,pver)/z0)

    ! we need them to the 3.41 power, according to Gong et al., 1997:
    u10cubed(:ncol)=u10cubed(:ncol)**3.41_r8

#if (defined MODAL_AERO)
! Calculations of source strength and size distribution
! NB the 0.1 is the dlogDp we have to multiplie with to get the flux, but the value dependence
! of course on what dlogDp you have. You will also have to change the sections of Dg if you use
! a different number of size bins with different intervals.

    W(:ncol)=3.84e-6_r8*u10cubed(:ncol)*0.1_r8 ! whitecap area

! calculate number flux fi (#/m2/s)
    fi(:,:)=0._r8
    do m=1,sections
       if (m.le.9)then
         fi(:ncol,m)=W(:ncol)*((sst(:ncol))*consta(1,m)+constb(1,m))
       elseif (m.ge.10.and.m.le.13)then
         fi(:ncol,m)=W(:ncol)*((sst(:ncol))*consta(2,m)+constb(2,m))
       elseif (m.ge.14.and.m.lt.22)then
         fi(:ncol,m)=W(:ncol)*((sst(:ncol))*consta(3,m)+constb(3,m))
       elseif (m.ge.22.and.m.le.40)then
! use Monahan
         fi(:ncol,m)=consta(4,m)*u10cubed(:ncol)
       endif
    enddo
#endif

    do m=1,nsst
#if (defined MODAL_AERO)
       if(spc_num_ndx(m) > 0) then
         mm=progseasalts_idx1()+spc_num_ndx(m)-spc_ndx(1)
!        cflx(:ncol,mm)=0.0_r8   ! we need to include dst_a1, dst_a3 to num_a1 and num_a3
         do i=1, sections
           if (Dg(i).ge.sst_sz_range_lo(m) .and. Dg(i).lt.sst_sz_range_hi(m)) then
              cflx(:ncol,mm)=cflx(:ncol,mm)+fi(:ncol,i)*ocnfrc(:ncol)*scalefactor
           endif
         enddo
       endif
       mm=progseasalts_idx1()+spc_ndx(m)-spc_ndx(1)
       cflx(:ncol,mm)=0.0_r8
       do i=1, sections
         if (Dg(i).ge.sst_sz_range_lo(m) .and. Dg(i).lt.sst_sz_range_hi(m)) then
           cflx(:ncol,mm)=cflx(:ncol,mm)+fi(:ncol,i)*ocnfrc(:ncol)*scalefactor & 
                          *4._r8/3._r8*pi*rdry(i)**3*dns_aer_sst  ! should use dry size, convert from number to mass flux (kg/m2/s)
         endif
       enddo
       sflx(:ncol)=sflx(:ncol)+cflx(:ncol,mm)

!      if(spc_num_ndx(m) > 0) then
!         mm=progseasalts_idx1()+spc_num_ndx(m)-spc_ndx(1)
!         cflx(:ncol,mm)=sst_source_num(m)* u10cubed(:ncol)*ocnfrc(:ncol)
!      endif
!      mm=progseasalts_idx1()+spc_ndx(m)-spc_ndx(1)
!      cflx(:ncol,mm)=sst_source(m)* u10cubed(:ncol)*ocnfrc(:ncol)
!      sflx(:ncol)=sflx(:ncol)+cflx(:ncol,mm)
#else
       mm=progseasalts_idx1()+m-1
       cflx(:ncol,mm)=sst_source(m)* u10cubed(:ncol)*ocnfrc(:ncol)       
       sflx(:ncol)=sflx(:ncol)+cflx(:ncol,mm)
#endif
       call outfld(trim(cnst_name(mm)) // 'SF',cflx(:,mm),pcols,lchnk)

       ! this is being done inside of the vertical diffusion automatically
       !         ptend%lq(m) = .true. ! tendencies for all progseasalts on
       !         ptend%q(:ncol,pver,mm) = cflx(:ncol,m)*gravit/state%pdel(:ncol,pver)
    enddo

    call outfld('SSTSFMBL',sflx(:),pcols,lchnk)

    !      write(42,*) cflx(1,1)
    return
  end subroutine progseasalts_emis_intr


  !------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: subroutine ProgseasaltsDryDep(c)
  !
  ! !INTERFACE:
  !
  subroutine ProgseasaltsDryDep(ncol,t,pmid,ram1,fv,vlc_dry,vlc_trb,vlc_grv,state)
    !
    ! !DESCRIPTION: 
    !
    ! Dry deposition for seasalts: modified from dust dry deposition following
    ! Xue Xie Tie's MOZART seasalts by NMM.  Sam Levis did hte first dust dry 
    ! deposition in CLM/CCSM2 from Charlie Zender's codes, modified by NMM for CAM.
    ! Sam's Notes in CLM:
    ! Determine Turbulent dry deposition for progseasalts. Calculate the turbulent 
    ! component of progseasalts dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the progseasalts dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CCSM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, progseasaltsini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
    ! Note that because sea salts' radius changes with humidity we cannot 
    ! precalculate slip coefficients, etc. in the initialization subroutine, but
    ! we have to calculate them at the time--how slow???
    !
    ! !USES
    !
    use physconst,     only: rair,pi,boltz
    use wv_saturation, only: qsat
    use physics_types, only: physics_state, physics_ptend
    use cam_history,   only: outfld

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8) :: t(pcols,pver)       !atm temperature (K)
    real(r8) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8) :: rho(pcols,pver)     !atm density (kg/m**3)
    real(r8) :: fv(pcols)           !friction velocity (m/s)
    real(r8) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8) :: vlc_trb(pcols,nsst)  !Turbulent deposn velocity (m/s)
    real(r8) :: vlc_grv(pcols,pver,nsst)  !grav deposn velocity (m/s)
    real(r8) :: vlc_dry(pcols,pver,nsst)  !dry deposn velocity (m/s)
    integer, intent(in) :: ncol
    type(physics_state), intent(in ) :: state          ! Physics state variables
    !
    ! !REVISION HISTORY
    ! Created by Sam Levis
    ! Modified for CAM by Natalie Mahowald
    ! Modified for Seasalts by NMM
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k,ix,lchnk          !indices
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(pcols,pver,nsst) ![frc] Slip correction factor
    real(r8) :: rss_lmn(nsst) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp ,r          !temporary 
    real(r8) :: wetdia(pcols,pver,nsst)        ! wet diameter of seasalts
    real(r8) :: RH(pcols,pver),es(pcols,pver),qs(pcols,pver)  ! for wet radius calculation

    ! constants
    real(r8),parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
    real(r8),parameter::shm_nbr_xpn_ocn=-1._r8/2._r8 ![frc] shm_nbr_xpn over ccean
    real(r8),parameter:: c1=0.7674_r8, c2=3.0790_r8, c3=2.57e-11_r8,c4=-1.424_r8  ! wet radius calculation constants


    ! needs fv and ram1 passed in from lnd model

    !------------------------------------------------------------------------

    call qsat(state%t(:ncol,:), state%pmid(:ncol,:), &
         es(:ncol,:),qs(:ncol,:))
    RH(:ncol,:)=state%q(:ncol,:,1)/qs(:ncol,:)
    RH(:ncol,:)=max(0.01_r8,min(0.99_r8,RH(:ncol,:)))
    ! set stokes correction to 1.0 for now not a bad assumption for our size range)
    do m=1,nsst
       stk_crc(m)=1._r8
    enddo
    do k=1,pver
       do i=1,ncol
          rho(i,k)=pmid(i,k)/rair/t(i,k)
          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho(i,k) ![m2 s-1] Kinematic viscosity of air

          do m = 1, nsst
             r=smt_vwr(m)/2.0_r8
             wetdia(i,k,m)=((r**3+c1*r**c2/(c3*r**c4-log(RH(i,k))))**(1._r8/3._r8))*2.0_r8
             slp_crc(i,k,m) = 1.0_r8 + 2.0_r8 * mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*wetdia(i,k,m)/(2.0_r8*mfp_atm(i,k)))) / &
                  wetdia(i,k,m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * wetdia(i,k,m) * wetdia(i,k,m) * dns_aer_sst * &
                  gravit * slp_crc(i,k,m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
             vlc_dry(i,k,m)=vlc_grv(i,k,m)
          end do

       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do m = 1, nsst
       do i=1,ncol
          r=smt_vwr(m)/2.0_r8
          wetdia(i,k,m)=((r**3+c1*r**c2/(c3*r**c4-log(RH(i,k))))**(1._r8/3._r8))*2.0_r8

          stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = boltz * t(i,k) * slp_crc(i,k,m) / &    ![m2 s-1]
               (3.0_r8*pi*vsc_dyn_atm(i,k)*wetdia(i,k,m)) !SeP97 p.474
          shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
          shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          !           if(ocnfrac.gt.0.5) shm_nbr_xpn=shm_nbr_xpn_ocn
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965

          rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
          vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
          vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
       end do !ncol

    end do

#if defined ( debugseasalts)
    ! debugging options for wetdiameter depedence on relative humidity
    lchnk = state%lchnk
    ix=ix_progseasalts
    do m=1,4
       call outfld( trim(cnst_name(m+ix-1))//'DI',wetdia(:,:,m), pcols, lchnk)
    enddo
    call outfld( 'RH',RH(:,:), pcols, lchnk)

#endif
    return
  end subroutine ProgseasaltsDryDep

#if (defined MODAL_AERO)
  subroutine seasalt_emitfactors_1bin( ireduce_smallr_emit,       &
             dpdrylo_cm, dpdryhi_cm, drydens, emitfact_numb, emitfact_mass )
!c
!c   computes seasalt emissions factors for a specifed
!c   dry particle size range
!c      dpdrylo_cm  = lower dry diameter (cm)
!c      dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c      number   emissions (#/m2/s) == emitfact_numb * (u10*3.41)
!c      dry-mass emissions (g/m2/s) == emitfact_mass * (u10*3.41)
!c
!c   where u10 = 10 m windspeed in m/s
!c !c   uses bubble emissions formula (eqn 5a) from
!c      Gong et al. [JGR, 1997, p 3805-3818]
!c
!c   *** for rdry < rdry_star, this formula overpredicts emissions.
!c      A strictly ad hoc correction is applied to the formula,
!c      based on sea-salt size measurements of
!c      O'Dowd et al. [Atmos Environ, 1997, p 73-80] !c
!c   *** the correction is only applied when ireduce_smallr_emit > 0
!c
        use mo_constants, only : pi

        implicit none

!c   subr arguments
        integer ireduce_smallr_emit
        real(r8), intent(in)  ::  dpdrylo_cm, dpdryhi_cm, drydens
        real(r8), intent(out) ::  emitfact_numb, emitfact_mass

!c   local variables
        integer isub_bin, nsub_bin

        real(r8) ::  alnrdrylo
        real(r8) ::  drydens_43pi_em12
        real(r8) ::  dum, dumadjust, dumb, dumexpb
        real(r8) ::  dumsum_na, dumsum_ma
        real(r8) ::  drwet, dlnrdry
        real(r8) ::  df0drwet, df0dlnrdry, df0dlnrdry_star
        real(r8) ::  relhum
        real(r8) ::  rdry, rdrylo, rdryhi, rdryaa, rdrybb
        real(r8) ::  rdrylowermost, rdryuppermost, rdry_star
        real(r8) ::  rwet, rwetaa, rwetbb
        real(r8) ::  rdry_cm, rwet_cm
        real(r8) ::  sigmag_star
        real(r8) ::  xmdry

        real(r8),parameter:: c1=0.7674_r8, c2=3.0790_r8, c3=2.57e-11_r8,c4=-1.424_r8  ! wet radius calculation constants

!c   factor for radius (micrometers) to mass (g)
        drydens_43pi_em12 = drydens*(4.0_r8/3.0_r8)*pi*1.0e-12_r8

!c   bubble emissions formula assume 80% RH
        relhum = 0.80_r8

!c   rdry_star = dry radius (micrometers) below which the
!c   dF0/dr emission formula is adjusted downwards
        rdry_star = 0.1_r8
        if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20_r8

!c   sigmag_star = geometric standard deviation used for
!c   rdry < rdry_star
        sigmag_star = 1.9_r8

!c   initialize sums
        dumsum_na = 0.0_r8
        dumsum_ma = 0.0_r8

!c   rdrylowermost, rdryuppermost = lower and upper
!c   dry radii (micrometers) for overall integration
        rdrylowermost = dpdrylo_cm*0.5e4_r8
        rdryuppermost = dpdryhi_cm*0.5e4_r8

!c
!c   "section 1"
!c   integrate over rdry > rdry_star, where the dF0/dr emissions formula is applicable
!c   (when ireduce_smallr_emit <= 0, rdry_star = -1.0e20,
!c    and the entire integration is done here)
!c
        if (rdryuppermost .le. rdry_star) goto 2000

!c   rdrylo, rdryhi = lower and upper dry radii (micrometers)
!c   for this part of the integration
        rdrylo = max( rdrylowermost, rdry_star )
        rdryhi = rdryuppermost

        nsub_bin = 1000

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

!c   compute rdry, rwet (micrometers) at lowest size
        rdrybb = exp( alnrdrylo )
        rdry_cm = rdrybb*1.0e-4_r8

        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**0.333_r8
        rwetbb = rwet_cm*1.0e4_r8

        do 1900 isub_bin = 1, nsub_bin

!c   rdry, rwet at sub_bin lower boundary are those
!c   at upper boundary of previous sub_bin
        rdryaa = rdrybb
        rwetaa = rwetbb

!c   compute rdry, rwet (micrometers) at sub_bin upper boundary
        dum = alnrdrylo + isub_bin*dlnrdry
        rdrybb = exp( dum )
        rdry_cm = rdrybb*1.0e-4_r8

        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**0.333_r8
        rwetbb = rwet_cm*1.0e4_r8

!c   geometric mean rdry, rwet (micrometers) for sub_bin
        rdry = sqrt(rdryaa * rdrybb)
        rwet = sqrt(rwetaa * rwetbb)
        drwet = rwetbb - rwetaa

!c   xmdry is dry mass in g
        xmdry = drydens_43pi_em12 * (rdry**3.0_r8)

!c   dumb is "B" in Gong's Eqn 5a
!c   df0drwet is "dF0/dr" in Gong's Eqn 5a
        dumb = ( 0.380_r8 - log10(rwet) ) / 0.650_r8
        dumexpb = exp( -dumb*dumb)
        df0drwet = 1.373_r8 * (rwet**(-3.0_r8)) *                        &
                (1.0_r8 + 0.057_r8*(rwet**1.05_r8)) *                    &
                (10._r8**(1.19_r8*dumexpb))

        dumsum_na = dumsum_na + drwet*df0drwet
        dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry

1900    continue

!c
!c   "section 2"
!c   integrate over rdry < rdry_star, where the dF0/dr emissions
!c   formula is just an extrapolation and predicts too many emissions
!c
!c   1.  compute dF0/dln(rdry) = (dF0/drwet)*(drwet/dlnrdry)
!c      at rdry_star
!c   2.  for rdry < rdry_star, assume dF0/dln(rdry) is lognormal,
!c      with the same lognormal parameters observed in
!c      O'Dowd et al. [1997]
!c

2000    if (rdrylowermost .ge. rdry_star) goto 3000

!c   compute dF0/dln(rdry) at rdry_star
        rdryaa = 0.99_r8*rdry_star
        rdry_cm = rdryaa*1.0e-4_r8
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**0.333_r8
        rwetaa = rwet_cm*1.0e4_r8

        rdrybb = 1.01_r8*rdry_star
        rdry_cm = rdrybb*1.0e-4_r8
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**0.333_r8
        rwetbb = rwet_cm*1.0e4_r8

        rwet = 0.5_r8*(rwetaa + rwetbb)
        dumb = ( 0.380_r8 - log10(rwet) ) / 0.650_r8
        dumexpb = exp( -dumb*dumb)
        df0drwet = 1.373_r8 * (rwet**(-3.0_r8)) *                        &
                (1.0_r8 + 0.057_r8*(rwet**1.05_r8)) *                    &
                (10.0_r8**(1.19_r8*dumexpb))

        drwet = rwetbb - rwetaa
        dlnrdry = log( rdrybb/rdryaa )
        df0dlnrdry_star = df0drwet * (drwet/dlnrdry)

!c   rdrylo, rdryhi = lower and upper dry radii (micrometers)
!c   for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = min( rdryuppermost, rdry_star )

        nsub_bin = 1000

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

        do 2900 isub_bin = 1, nsub_bin

!c   geometric mean rdry (micrometers) for sub_bin
        dum = alnrdrylo + (isub_bin-0.5_r8)*dlnrdry
        rdry = exp( dum )

!c   xmdry is dry mass in g
        xmdry = drydens_43pi_em12 * (rdry**3.0_r8)

!c   dumadjust is adjustment factor to reduce dF0/dr
        dum = log( rdry/rdry_star ) / log( sigmag_star )
        dumadjust = exp( -0.5_r8*dum*dum )

        df0dlnrdry = df0dlnrdry_star * dumadjust

        dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
        dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry

2900    continue
!c
!c  all done
!c
3000    emitfact_numb = dumsum_na*scalefactor !++ ag: scale sea-salt
        emitfact_mass = dumsum_ma*scalefactor !++ ag: scale sea-salt

        return
   end subroutine seasalt_emitfactors_1bin
#endif


end module progseasalts_intr
