

module iondrag
  !-------------------------------------------------------------------------------
  ! Purpose:
  !   Calculate ion drag tendency and apply to horizontal velocities.
  !   Also calculate joule heating tendency and apply to neutral temperature. 
  ! 
  ! Subroutines:
  !   iondrag_init (initialize module)
  !   iondrag_calc (calculate ion drag tensors)
  !   iondrag_tend (ion drag tendency)   
  !   qjoule_tend (joule heating)
  !
  ! Calling sequence:
  !   inti
  !     iondrag_init
  !   tphysac
  !     iondrag_calc
  !       iondrag_tend
  !       qjoule_tend
  !
  ! Dependencies:
  !   Magnetic field from apex module
  !   ExB ion drifts from exbdrift module
  !
  ! Author:
  !   B. Foster Feb, 2004
  !
  !-------------------------------------------------------------------------------

  use shr_kind_mod ,only: r8 => shr_kind_r8
  use ppgrid       ,only: pcols, pver
  use cam_history  ,only: addfld, add_default, outfld
  use physics_types,only: physics_state, physics_ptend, physics_ptend_init
  
  use physics_buffer, only : pbuf_get_index, physics_buffer_desc, pbuf_get_field
  use perf_mod     ,only: t_startf, t_stopf
  use cam_logfile  ,only: iulog

  use interpolate_data, only: lininterp
  use spmd_utils,       only: masterproc

  implicit none

  save

  private                         ! Make default type private to the module


  !-------------------------------------------------------------------------------
  ! Public interfaces:
  !-------------------------------------------------------------------------------
  public :: iondrag_register         ! Register variables in pbuf physics buffer
  public :: iondrag_init             ! Initialization
  public :: iondrag_calc             ! ion drag tensors lxx,lyy,lxy,lyx
  public :: iondrag_readnl
  public :: do_waccm_ions

  interface iondrag_calc
     module procedure iondrag_calc_ions
     module procedure iondrag_calc_ghg
  end interface

  !-------------------------------------------------------------------------------
  ! Private data
  !-------------------------------------------------------------------------------

  ! Namelist variables
  character(len=256) :: efield_lflux_file = 'coeff_lflux.dat'
  character(len=256) :: efield_hflux_file = 'coeff_hflux.dat'
  character(len=256) :: efield_wei96_file = 'wei96.cofcnts'

  real(r8),parameter :: amu    = 1.6605387e-27_r8  ! atomic mass unit (kg)

  integer  :: ntop_lev = 1
  integer  :: nbot_lev = 0
  integer  :: id_xo2, id_xo1             ! indices to tn and major sp
  integer  :: id_o2p, id_op, id_nop      ! indices to ions
  integer  :: id_elec, id_n

  !Physics buffer indices
  integer :: ue_idx      = 0
  integer :: ve_idx      = 0
  integer :: PedConduct_idx  = 0
  integer :: HallConduct_idx = 0


  logical  :: xo2_slvd, xo1_slvd, o2p_slvd, op_slvd, nop_slvd, elec_slvd, n_slvd

  real(r8) :: rmass_elec                 ! mass of electron (g/mole)
  real(r8) :: rmass_op                   ! mass of O+ (g/mole)
  real(r8) :: rmass_o2p                  ! mass of O2+ (g/mole)
  real(r8) :: rmass_nop                  ! mass of NO+ (g/mole)
  real(r8) :: rmass_o1                   ! mass of O (g/mole)
  real(r8) :: rmass_o2                   ! mass of O2 (g/mole)
  real(r8) :: rmass_n2                   ! mass of N2 (g/mole)

  !-------------------------------------------------------------------------------
  ! Inverted masses (for multiply in loops rather than divide):
  !-------------------------------------------------------------------------------
  real(r8) :: rmi_o1
  real(r8) :: rmi_o2
  real(r8) :: rmi_n2
  real(r8) :: rmi_op
  real(r8) :: rmi_o2p
  real(r8) :: rmi_nop
  real(r8) :: rmi_elec
  real(r8) :: rmi_op_kg
  real(r8) :: rmi_o2p_kg
  real(r8) :: rmi_nop_kg

  ! GHG
  !-------------------------------------------------------------------------

  ! Private data
  integer, parameter :: plevtiod = 97   

  real(r8) alamxx(plevtiod)
  real(r8) alamxy(plevtiod)

  real(r8) pshtiod(plevtiod)           ! TIME pressure scale height
  real(r8) pshccm(pver)                ! CCM pressure scale height

  real(r8) alamxxi(pver)              ! alamxx interpolated to waccm grid
  real(r8) alamxyi(pver)              ! alamxy interpoalted to waccm grid

  logical doiodrg
  logical :: do_waccm_ions

  !
  ! Data statement for ALAMXX
  !
  data alamxx /                                                     &
       0.13902E-17_r8, 0.22222E-17_r8, 0.34700E-17_r8, 0.53680E-17_r8, 0.83647E-17_r8, &
       0.13035E-16_r8, 0.20254E-16_r8, 0.31415E-16_r8, 0.48944E-16_r8, 0.75871E-16_r8, &
       0.11584E-15_r8, 0.17389E-15_r8, 0.25786E-15_r8, 0.37994E-15_r8, 0.58088E-15_r8, &
       0.95179E-15_r8, 0.19052E-14_r8, 0.47869E-14_r8, 0.14284E-13_r8, 0.45584E-13_r8, &
       0.14756E-12_r8, 0.48154E-12_r8, 0.14844E-11_r8, 0.39209E-11_r8, 0.83886E-11_r8, &
       0.14213E-10_r8, 0.20304E-10_r8, 0.27449E-10_r8, 0.39276E-10_r8, 0.59044E-10_r8, &
       0.83683E-10_r8, 0.11377E-09_r8, 0.14655E-09_r8, 0.19059E-09_r8, 0.28338E-09_r8, &
       0.46326E-09_r8, 0.73966E-09_r8, 0.11785E-08_r8, 0.18789E-08_r8, 0.31037E-08_r8, &
       0.53919E-08_r8, 0.97251E-08_r8, 0.17868E-07_r8, 0.33041E-07_r8, 0.61265E-07_r8, &
       0.11406E-06_r8, 0.20912E-06_r8, 0.39426E-06_r8, 0.76691E-06_r8, 0.15113E-05_r8, &
       0.29545E-05_r8, 0.55644E-05_r8, 0.97208E-05_r8, 0.16733E-04_r8, 0.28101E-04_r8, &
       0.36946E-04_r8, 0.44277E-04_r8, 0.50982E-04_r8, 0.57526E-04_r8, 0.64190E-04_r8, &
       0.71471E-04_r8, 0.80311E-04_r8, 0.96121E-04_r8, 0.11356E-03_r8, 0.14131E-03_r8, &
       0.18695E-03_r8, 0.26058E-03_r8, 0.36900E-03_r8, 0.50812E-03_r8, 0.66171E-03_r8, &
       0.80763E-03_r8, 0.92583E-03_r8, 0.10038E-02_r8, 0.10382E-02_r8, 0.10333E-02_r8, &
       0.99732E-03_r8, 0.93994E-03_r8, 0.86984E-03_r8, 0.79384E-03_r8, 0.71691E-03_r8, &
       0.64237E-03_r8, 0.57224E-03_r8, 0.50761E-03_r8, 0.44894E-03_r8, 0.39624E-03_r8, &
       0.34929E-03_r8, 0.30767E-03_r8, 0.27089E-03_r8, 0.23845E-03_r8, 0.20985E-03_r8, &
       0.18462E-03_r8, 0.16233E-03_r8, 0.14260E-03_r8, 0.12510E-03_r8, 0.10955E-03_r8, &
       0.95699E-04_r8, 0.83347E-04_r8/

  !
  ! Data statement for ALAMXY
  !
  data alamxy /                                                     &
       0.74471E-24_r8, 0.22662E-23_r8, 0.69004E-23_r8, 0.20345E-22_r8, 0.58465E-22_r8, &
       0.16542E-21_r8, 0.46240E-21_r8, 0.12795E-20_r8, 0.35226E-20_r8, 0.96664E-20_r8, &
       0.26650E-19_r8, 0.76791E-19_r8, 0.25710E-18_r8, 0.10897E-17_r8, 0.56593E-17_r8, &
       0.30990E-16_r8, 0.16792E-15_r8, 0.85438E-15_r8, 0.40830E-14_r8, 0.18350E-13_r8, &
       0.79062E-13_r8, 0.33578E-12_r8, 0.13348E-11_r8, 0.45311E-11_r8, 0.12443E-10_r8, &
       0.27052E-10_r8, 0.49598E-10_r8, 0.86072E-10_r8, 0.15807E-09_r8, 0.30480E-09_r8, &
       0.55333E-09_r8, 0.96125E-09_r8, 0.15757E-08_r8, 0.25896E-08_r8, 0.48209E-08_r8, &
       0.96504E-08_r8, 0.18494E-07_r8, 0.34296E-07_r8, 0.61112E-07_r8, 0.10738E-06_r8, &
       0.18747E-06_r8, 0.32054E-06_r8, 0.52872E-06_r8, 0.83634E-06_r8, 0.12723E-05_r8, &
       0.18748E-05_r8, 0.26362E-05_r8, 0.36986E-05_r8, 0.52079E-05_r8, 0.72579E-05_r8, &
       0.98614E-05_r8, 0.12775E-04_r8, 0.15295E-04_r8, 0.18072E-04_r8, 0.20959E-04_r8, &
       0.19208E-04_r8, 0.16285E-04_r8, 0.13628E-04_r8, 0.11784E-04_r8, 0.11085E-04_r8, &
       0.11916E-04_r8, 0.14771E-04_r8, 0.20471E-04_r8, 0.29426E-04_r8, 0.42992E-04_r8, &
       0.62609E-04_r8, 0.90224E-04_r8, 0.12870E-03_r8, 0.18281E-03_r8, 0.26029E-03_r8, &
       0.37224E-03_r8, 0.53254E-03_r8, 0.75697E-03_r8, 0.10623E-02_r8, 0.14660E-02_r8, &
       0.19856E-02_r8, 0.26393E-02_r8, 0.34473E-02_r8, 0.44327E-02_r8, 0.56254E-02_r8, &
       0.70672E-02_r8, 0.88174E-02_r8, 0.10960E-01_r8, 0.13613E-01_r8, 0.16934E-01_r8, &
       0.21137E-01_r8, 0.26501E-01_r8, 0.33388E-01_r8, 0.42263E-01_r8, 0.53716E-01_r8, &
       0.68491E-01_r8, 0.87521E-01_r8, 0.11196E+00_r8, 0.14320E+00_r8, 0.18295E+00_r8, &
       0.23321E+00_r8, 0.29631E+00_r8/

contains

!==============================================================================     

  subroutine iondrag_register
!-----------------------------------------------------------------------
! Register E and B fields.
!
! Register iondrag variables with physics buffer:
!
! Hall and Pedersen conductivities
! 
! pcols dimension and lchnk assumed here
!
!-----------------------------------------------------------------------
    use exbdrift,        only: exbdrift_register
    use physics_buffer,  only: pbuf_add_field, dtype_r8
    use phys_control,    only: waccmx_is

    ! E and B fields
    call exbdrift_register()

    if ( waccmx_is("ionosphere") .or. waccmx_is("neutral") ) then
       ! Pedersen Conductivity and Hall Conductivity
       call pbuf_add_field('PedConduct',  'physpkg', dtype_r8, (/pcols,pver/), PedConduct_idx )
       call pbuf_add_field('HallConduct', 'physpkg', dtype_r8, (/pcols,pver/), HallConduct_idx)
    end if

  end subroutine iondrag_register

!================================================================================================

  subroutine iondrag_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,  only: endrun

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'iondrag_readnl'

    namelist /iondrag_nl/ efield_lflux_file, efield_hflux_file, efield_wei96_file

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'iondrag_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, iondrag_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)

    end if

#ifdef SPMD
   call mpibcast (efield_lflux_file, len(efield_lflux_file), mpichar, 0, mpicom)
   call mpibcast (efield_hflux_file, len(efield_hflux_file), mpichar, 0, mpicom)
   call mpibcast (efield_wei96_file, len(efield_wei96_file), mpichar, 0, mpicom)
#endif

  end subroutine iondrag_readnl

  !================================================================================================

  subroutine iondrag_init( pref_mid )
    use constituents, only: cnst_get_ind
    use short_lived_species, only: slvd_index

    !-------------------------------------------------------------------------------
    ! Iondrag initialization, called from inti.F90.
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    real(r8), intent(in) :: pref_mid(pver)

    integer :: k

    integer :: cnst_ids(7)

    doiodrg = .false.
    do_waccm_ions = .false.

    !-------------------------------------------------------------------------------
    ! find lower bnd for iondrag
    !-------------------------------------------------------------------------------
    if( pref_mid(1) < 0.1_r8 ) then
       do k = 1, pver
          if (pref_mid(k) < 50._r8) nbot_lev = k
       end do
    end if

    if (nbot_lev > 0) then
       doiodrg = .true.
    endif

    if ( .not. doiodrg .and. masterproc ) then
       write(iulog,*) ' '
       write(iulog,*) 'iondrag_init: Does not have waccm level. Ion drag does not apply. '
       write(iulog,*) ' '
       return
    endif

    if( masterproc ) then
       write(iulog,*) ' '
       write(iulog,*) 'iondrag_init: nbot_lev,press = ',nbot_lev,pref_mid(nbot_lev)
       write(iulog,*) ' '
    end if

    call cnst_get_ind( 'e', id_elec, abort=.false. )
    if (id_elec < 0) then
       id_elec = slvd_index( 'e' )
       if (id_elec > 0) then
          elec_slvd = .true.
       endif
    else
       elec_slvd = .false.
    endif
    call cnst_get_ind( 'Op', id_op, abort=.false. )
    if (id_op < 0) then
       id_op = slvd_index( 'Op' )
       if (id_op > 0) then
          op_slvd = .true.
       endif
    else
       op_slvd = .false.
    endif
    call cnst_get_ind( 'O2p', id_o2p, abort=.false. )
    if (id_o2p < 0) then
       id_o2p = slvd_index( 'O2p' )
       if (id_o2p > 0) then
          o2p_slvd = .true.
       endif
    else
       o2p_slvd = .false.
    endif
    call cnst_get_ind( 'NOp', id_nop, abort=.false. )
    if (id_nop < 0) then
       id_nop = slvd_index( 'NOp' )
       if (id_nop > 0) then
          nop_slvd = .true.
       endif
    else
       nop_slvd = .false.
    endif
    call cnst_get_ind( 'O', id_xo1, abort=.false. )
    if (id_xo1 < 0) then
       id_xo1 = slvd_index( 'O' )
       if (id_xo1 > 0) then
          xo1_slvd = .true.
       endif
    else
       xo1_slvd = .false.
    endif
    call cnst_get_ind( 'O2', id_xo2, abort=.false. )
    if (id_xo2 < 0) then
       id_xo2 = slvd_index( 'O2' )
       if (id_xo2 > 0) then
          xo2_slvd = .true.
       endif
    else
       xo2_slvd = .false.
    endif
    call cnst_get_ind( 'N', id_n, abort=.false. )
    if (id_n < 0) then
       id_n = slvd_index( 'N' )
       if (id_n > 0) then
          n_slvd = .true.
       endif
    else
       n_slvd = .false.
    endif
 
    cnst_ids = (/ id_elec, id_op, id_o2p, id_nop, id_xo1, id_xo2, id_n /)

    if ( all( cnst_ids > 0 ) ) then
       do_waccm_ions = .true.
    endif

    if ( do_waccm_ions ) then
       call ions_init
    else
       call ghg_init(pref_mid)
    endif

    if (.not.doiodrg) return

    ! Add to masterfield list
    call addfld('UIONTEND',(/ 'lev' /),'A','M/S2','u-tendency due to ion drag')
    call addfld('VIONTEND',(/ 'lev' /),'A','M/S2','v-tendency due to ion drag')

    ue_idx = pbuf_get_index('UE')
    ve_idx = pbuf_get_index('VE')

  end subroutine iondrag_init

  !================================================================================================
  subroutine ions_init

    use constituents, only: cnst_mw
    use efield,       only: efield_init
    use exbdrift,     only: exbdrift_init
    use mo_apex,      only: apexmag
    use mo_chem_utls, only: get_spc_ndx
    use chem_mods,    only: adv_mass
    use mo_chem_utls, only: get_inv_ndx
    use chem_mods,    only: fix_mass

    !-------------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------------
    real(r8) :: rpsh
    integer  :: id

    !-------------------------------------------------------------------------------
    ! initialize related packages: electric field
    !-------------------------------------------------------------------------------
    call apexmag

    call efield_init (efield_lflux_file, efield_hflux_file, efield_wei96_file)
    call exbdrift_init

    id = get_spc_ndx('e')
    rmass_elec = adv_mass(id)
    id = get_spc_ndx('Op')
    rmass_op  = adv_mass(id)
    id = get_spc_ndx('O2p')
    rmass_o2p = adv_mass(id)
    id = get_spc_ndx('NOp')
    rmass_nop = adv_mass(id)
    id = get_spc_ndx('O')
    rmass_o1  = adv_mass(id)
    id = get_spc_ndx('O2')
    rmass_o2 = adv_mass(id)
    id = get_inv_ndx('N2')
    rmass_n2 = fix_mass(id)

    rmi_elec   = 1._r8/rmass_elec
    rmi_o1     = 1._r8/rmass_o1
    rmi_o2     = 1._r8/rmass_o2
    rmi_n2     = 1._r8/rmass_n2
    rmi_op     = 1._r8/rmass_op
    rmi_o2p    = 1._r8/rmass_o2p
    rmi_nop    = 1._r8/rmass_nop
    rmi_op_kg  = 1._r8/(rmass_op *amu)
    rmi_o2p_kg = 1._r8/(rmass_o2p*amu)  
    rmi_nop_kg = 1._r8/(rmass_nop*amu)

    !-------------------------------------------------------------------------------
    ! Set up fields to history files.
    !-------------------------------------------------------------------------------

    call addfld('QIONSUM',(/ 'lev' /),'I','S-1' ,'Ion prod sum')
    call addfld('ELECDEN',(/ 'lev' /),'I','CM-3','NE (ion sum)')
    call addfld( 'SIGMAPED', (/ 'lev' /), 'I', 'siemens/m', 'Pederson conductivity' )
    call addfld( 'SIGMAHAL', (/ 'lev' /), 'I', 'siemens/m', 'Hall conductivity' )
    call addfld('LAMDA1'  ,(/ 'lev' /),'I' ,'S-1','LAMDA PED')
    call addfld('LAMDA2'  ,(/ 'lev' /),'I' ,'S-1','LAMDA HALL')

    call addfld('LXX',(/ 'lev' /),'I','S-1','LXX')
    call addfld('LYY',(/ 'lev' /),'I','S-1','LYY')
    call addfld('LXY',(/ 'lev' /),'I','S-1','LXY')
    call addfld('LYX',(/ 'lev' /),'I','S-1','LYX')
    !
    ! Joule heating, and tn before and after joule heating tendencies are applied:
    !
    call addfld( 'QJOULE', (/ 'lev' /), 'I', 'K/s' , 'Joule Heat' )  ! joule heating
    call add_default( 'QJOULE  ', 1, ' ' )                                   ! joule heating (K/s)

  end subroutine ions_init

  !========================================================================

  subroutine ghg_init (pref_mid)

    !
    ! initialization for ion drag calculation
    !

    use ppgrid
    use cam_history, only: add_default, addfld
    use cam_abortutils,  only: endrun

    !------------------Input arguments---------------------------------------

    real(r8), intent(in) :: pref_mid(pver)           ! model ref pressure at midpoint   

    !-----------------local workspace---------------------------------------
    integer k
    integer kinv

    real(r8) rpsh                                    ! ref pressure scale height

    real(r8), parameter :: preftgcm = 5.e-5_r8       ! TIME GCM reference pressure (Pa)

    !------------------------------------------------------------------------

    ! With the defualt values of nbot_lev and ntop_lev, ion drag calcualtion are NOT carried out
    nbot_lev=0
    ntop_lev=1

    do k = 1, pver
       rpsh=log(1e5_r8/pref_mid(k))
       if (rpsh .gt. 14._r8) nbot_lev  = k
    end do
    if (nbot_lev .gt. ntop_lev) doiodrg=.true.
    if (masterproc) then
       write(iulog,fmt='(a15)') 'From IONDRAGI:'
       write(iulog,fmt='(1a12,1i10)') 'NTOP_LEV  =',ntop_lev
       write(iulog,fmt='(1a12,1i10)') 'NBOT_LEV  =',nbot_lev
       write(iulog,*) 'IONDRAG flag is',doiodrg
    endif
    if (.not.doiodrg) return

    !     obtain TIME/GCM pressure scale height
    pshtiod(1)=-17._r8
    do k=2,plevtiod
       pshtiod(k)=pshtiod(k-1)+0.25_r8
    enddo

    !     map TIME-psh into CCM-psh 
    pshtiod=pshtiod-log(preftgcm/1E5_r8)

    !     CCM psh
    !     note that vertical indexing is inverted with respect to CCM standard
    do k=1,pver
       kinv=pver-k+1
       pshccm(kinv)=log(1e5_r8/pref_mid(k))
    enddo

    !     vertical interpolation 
    write(iulog,*) ' '
    write(iulog,*) 'iondragi: before lininterp for alamxx'
    write(iulog,*) '          nlatin,nlatout =',plevtiod,pver
    write(iulog,*) '          yin'
    write(iulog,'(1p,5g15.8)') pshtiod
    write(iulog,*) '          yout'
    write(iulog,'(1p,5g15.8)') pshccm
    write(iulog,*) ' '

    call lininterp (alamxx  ,pshtiod,plevtiod, alamxxi   ,pshccm,pver)

    call lininterp (alamxy  ,pshtiod,plevtiod, alamxyi   ,pshccm,pver)

    !     invert indeces back to CCM convention
    alamxxi(1:pver)=alamxxi(pver:1:-1)
    alamxyi(1:pver)=alamxyi(pver:1:-1)

    return

  end subroutine ghg_init


  !================================================================================================
  subroutine iondrag_calc_ions( lchnk, ncol, state, ptend, pbuf,  delt )
    !-------------------------------------------------------------------------------
    ! Calculate ion drag tensors lxx,lyy,lxy,lyx.
    ! Also calculate Pedersen and Hall conductivities.
    ! This is called from tphysac.
    !-------------------------------------------------------------------------------

    use constituents ,only: cnst_get_ind      ! get constituent index function
    use mo_apex,only: & ! (pcols,begchunk:endchunk)
         bnorth,      & ! northward component of magnetic field (nT)
         beast,       & ! eastward component of magnetic field (nT)
         bdown,       & ! downward component of magnetic field (nT)
         bmag           ! magnetic field magnitude (nT)
    use physconst,     only: avogad, boltz
    use chemistry,     only: imozart
    use mo_mean_mass,  only: set_mean_mass
    use exbdrift,      only: get_exbdrift
    use short_lived_species, only: slvd_pbf_ndx => pbf_idx
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_set_field
    use phys_control,  only: waccmx_is

    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    integer,intent(in)   :: lchnk               ! current chunk index
    integer,intent(in)   :: ncol                ! number of atmospheric columns
    real(r8), intent(in) :: delt                ! time step (s)
    type(physics_state), intent(in), target    :: state ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend ! Physics tendencies
    
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-------------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------------
    integer :: i,k ! loop indices
    real(r8) :: sqrt_te                              ! sqrt(te)
    real(r8) :: sqrt_tnti                            ! sqrt(tnti)
    real(r8) :: wrk
    real(r8),parameter :: dipmin = 0.17_r8           ! minimum dip angle (tuneable)
    real(r8),parameter :: emass  = 9.1093819e-31_r8  ! electron mass (kg)
    real(r8),parameter :: qe     = 1.6021765e-19_r8  ! electronic charge (coulombs)
    real(r8),parameter :: colfac = 1.5_r8            ! collision factor (tuneable)
    real(r8),parameter :: boltzmann  = 1.e7_r8 * boltz
    real(r8),parameter :: avo    = avogad*1.e-3_r8 ! (molecules/mole)
    ! real(r8),parameter :: rmass_op  = 15.9989_r8     ! mass of O+
    ! real(r8),parameter :: rmass_o2p = 31.9983_r8     ! mass of O2+
    ! real(r8),parameter :: rmass_nop = 30.0056_r8     ! mass of NO+
    ! real(r8),parameter :: rmass_o1  = 16._r8         ! mass of O
    ! real(r8),parameter :: rmass_o2  = 32._r8         ! mass of O2
    ! real(r8),parameter :: rmass_n2  = 28._r8         ! mass of N2

    !-------------------------------------------------------------------------------
    ! Inverted masses (for multiply in loops rather than divide):
    !-------------------------------------------------------------------------------
    ! real(r8),parameter :: rmi_o1     = 1._r8/rmass_o1
    ! real(r8),parameter :: rmi_o2     = 1._r8/rmass_o2
    ! real(r8),parameter :: rmi_n2     = 1._r8/rmass_n2
    ! real(r8),parameter :: rmi_op     = 1._r8/rmass_op
    ! real(r8),parameter :: rmi_o2p    = 1._r8/rmass_o2p
    ! real(r8),parameter :: rmi_nop    = 1._r8/rmass_nop
    ! real(r8),parameter :: rmi_op_kg  = 1._r8/(rmass_op *amu)
    ! real(r8),parameter :: rmi_o2p_kg = 1._r8/(rmass_o2p*amu)  
    ! real(r8),parameter :: rmi_nop_kg = 1._r8/(rmass_nop*amu)

    real(r8) :: tn      (pcols,pver)   ! neutral gas temperature (deg K)
    real(r8) :: ti      (pcols,pver)   ! ion temperature (deg K)
    real(r8) :: te      (pcols,pver)   ! electron temperature (deg K)
    real(r8) :: xo2     (pcols,pver)   ! O2  (mmr)
    real(r8) :: xo1     (pcols,pver)   ! O   (mmr)
    real(r8) :: xn2     (pcols,pver)   ! N2  (mmr)
    real(r8) :: o2p     (pcols,pver)   ! O2+ (mmr)
    real(r8) :: op      (pcols,pver)   ! O+  (mmr)
    real(r8) :: nop     (pcols,pver)   ! NO+ (mmr)
    real(r8) :: barm    (pcols,pver)   ! mean molecular weight (g/mole)
    real(r8) :: xnmbar  (pcols,pver)   ! for unit conversion to volume density
    real(r8) :: tnti    (pcols,pver)   ! average of tn and ti
    real(r8) :: o2_cm3  (pcols,pver)   ! o2 volume density (cm-3)
    real(r8) :: o1_cm3  (pcols,pver)   ! o  volume density (cm-3)
    real(r8) :: n2_cm3  (pcols,pver)   ! n2 volume density (cm-3)
    real(r8) :: o2p_cm3 (pcols,pver)   ! O2+ (cm-3)
    real(r8) :: op_cm3  (pcols,pver)   ! O+  (cm-3)
    real(r8) :: nop_cm3 (pcols,pver)   ! NO+ (cm-3)
    real(r8) :: ne      (pcols,pver)   ! electron density (assume o2p+op+nop)
    real(r8) :: lamda1  (pcols,pver)   ! sigped*b**2/rho
    real(r8) :: lamda2  (pcols,pver)   ! sighal*b**2/rho
    real(r8) :: lxxnorot(pcols,pver)   ! XX before rotation
    real(r8) :: lyynorot(pcols,pver)   ! YY before rotation
    real(r8) :: lxynorot(pcols,pver)   ! XY before rotation
    real(r8) :: lyxnorot(pcols,pver)   ! YX before rotation

    !-------------------------------------------------------------------------------
    ! Ion-neutral momentum transfer collision frequencies.
    ! rnu_xxx_xx = ratio of collision to gyro-frequences for O2+, O+, NO+.
    !-------------------------------------------------------------------------------
    real(r8) :: rnu_o2p_o2(pcols,pver)   ! O2+ ~ O2 collision freq (resonant, T dependent)
    real(r8) :: rnu_op_o2 (pcols,pver)   ! O+  ~ O2 collision freq (non-resonant)
    real(r8) :: rnu_nop_o2(pcols,pver)   ! NO+ ~ O2 collision freq (non-resonant)
    real(r8) :: rnu_o2p_o (pcols,pver)   ! O2+ ~ O  collision freq (non-resonant)
    real(r8) :: rnu_op_o  (pcols,pver)   ! O+  ~ O  collision freq (resonant, T dependent)
    real(r8) :: rnu_nop_o (pcols,pver)   ! NO+ ~ O  collision freq (non-resonant)
    real(r8) :: rnu_o2p_n2(pcols,pver)   ! O2+ ~ N2 collision freq (non-resonant)
    real(r8) :: rnu_op_n2 (pcols,pver)   ! O+  ~ N2 collision freq (non-resonant)
    real(r8) :: rnu_nop_n2(pcols,pver)   ! NO+ ~ N2 collision freq (non-resonant)
    real(r8) :: rnu_o2p   (pcols,pver)   ! [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
    real(r8) :: rnu_op    (pcols,pver)   ! [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
    real(r8) :: rnu_nop   (pcols,pver)   ! [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
    real(r8) :: rnu_ne    (pcols,pver)   ! electron ~ neutral collision frequency (s-1)

    real(r8) :: press        (pcols)     ! pressure at interface levels (dyne/cm^2)
    real(r8) :: qe_fac       (pcols)     ! unit conversion factor for conductivities
    real(r8) :: dipmag       (pcols)     ! magnetic dip angle 
    real(r8) :: decmag       (pcols)     ! magnetic declination
    real(r8) :: btesla       (pcols)     ! magnetic field (teslas)
    real(r8) :: sindip       (pcols)     ! sin(dipmag)
    real(r8) :: cosdip       (pcols)     ! cos(dipmag)
    real(r8) :: sin2dip      (pcols)     ! sindip^2
    real(r8) :: cos2dip      (pcols)     ! cosdip^2
    real(r8) :: sindec       (pcols)     ! sin(decmag)
    real(r8) :: cosdec       (pcols)     ! cos(decmag)
    real(r8) :: sin2dec      (pcols)     ! sindec^2
    real(r8) :: cos2dec      (pcols)     ! cosdec^2
    real(r8) :: omega_o2p    (pcols)     ! angular gyrofrequency for o2+ (s-1)
    real(r8) :: omega_op     (pcols)     ! angular gyrofrequency for o+ (s-1)
    real(r8) :: omega_nop    (pcols)     ! angular gyrofrequency for no+ (s-1)
    real(r8) :: omega_e      (pcols)     ! electron angular gyrofrequency (s-1)
    real(r8) :: omega_o2p_inv(pcols)     ! inverse of o2+ gyrofrequency
    real(r8) :: omega_op_inv (pcols)     ! inverse of o+  gyrofrequency
    real(r8) :: omega_nop_inv(pcols)     ! inverse of no+ gyrofrequency
    real(r8) :: omega_e_inv  (pcols)     ! inverse of electron gyrofrequency

    !-------------------------------------------------------------------------------
    ! Ion drag coefficients output:
    !-------------------------------------------------------------------------------
    real(r8) :: lxx(pcols,pver)  ! lambda XX coefficients (s-1)
    real(r8) :: lyy(pcols,pver)  ! lambda YY coefficients (s-1)
    real(r8) :: lxy(pcols,pver)  ! lambda XY coefficients (s-1)
    real(r8) :: lyx(pcols,pver)  ! lambda YX coefficients (s-1)

    !-------------------------------------------------------------------------------
    ! Conductivities output:
    !-------------------------------------------------------------------------------
    real(r8) :: sigma_ped (pcols,pver)   ! pedersen conductivity (siemens/m)
    real(r8) :: sigma_hall(pcols,pver)   ! hall conductivity (siemens/m)

    real(r8) :: qout(pcols,pver)         ! temp for outfld

    real(r8), dimension(:,:), pointer :: q_elec, q_xo1,  q_xo2, q_o2p, q_op, q_nop

    if (.not.doiodrg) return

    if ( elec_slvd ) then
      call pbuf_get_field(pbuf, slvd_pbf_ndx, q_elec, start=(/1,1,id_elec/), kount=(/pcols,pver,1/) )
    else
      q_elec => state%q(:,:,id_elec)
    endif
    if ( xo1_slvd ) then
      call pbuf_get_field(pbuf, slvd_pbf_ndx, q_xo1, start=(/1,1,id_xo1/), kount=(/pcols,pver,1/) )
    else
      q_xo1 => state%q(:,:,id_xo1)
    endif
    if ( xo2_slvd ) then
      call pbuf_get_field(pbuf, slvd_pbf_ndx, q_xo2, start=(/1,1,id_xo2/), kount=(/pcols,pver,1/) )
    else
      q_xo2 => state%q(:,:,id_xo2)
    endif
    if ( o2p_slvd ) then
      call pbuf_get_field(pbuf, slvd_pbf_ndx, q_o2p, start=(/1,1,id_o2p/), kount=(/pcols,pver,1/) )
    else
      q_o2p => state%q(:,:,id_o2p)
    endif
    if ( op_slvd ) then
      call pbuf_get_field(pbuf, slvd_pbf_ndx, q_op, start=(/1,1,id_op/), kount=(/pcols,pver,1/) )
    else
      q_op => state%q(:,:,id_op)
    endif
    if ( nop_slvd ) then
      call pbuf_get_field(pbuf, slvd_pbf_ndx, q_nop, start=(/1,1,id_nop/), kount=(/pcols,pver,1/) )
    else
      q_nop => state%q(:,:,id_nop)
    endif

    !-------------------------------------------------------------------------------
    ! Define local tn and major species from state (mmr): 
    !-------------------------------------------------------------------------------
    do k = 1,pver
       do i = 1,ncol
          tn (i,k) = state%t(i,k)
          ti (i,k) = tn(i,k)                 ! assume ti==tn for now (waccm2, 1/04)
          te (i,k) = tn(i,k)                 ! assume te==tn for now (waccm2, 1/04)
          xo2(i,k) = q_xo2(i,k)              ! o2  (mmr)
          xo1(i,k) = q_xo1(i,k)              ! o   (mmr)
          xn2(i,k) = 1._r8 - (xo2(i,k) + xo1(i,k)) ! n2  (mmr)
          xn2(i,k) = max( 1.e-20_r8,xn2(i,k) )
          o2p(i,k) = q_o2p(i,k)              ! o2+ (mmr)
          op (i,k) = q_op(i,k)               ! o+  (mmr)
          nop(i,k) = q_nop(i,k)              ! no+ (mmr)
       end do
    end do

    ! call outfld ('TN      ',tn  ,pcols,lchnk)
    ! call outfld ('XO2     ',xo2 ,pcols,lchnk)
    ! call outfld ('XO      ',xo1 ,pcols,lchnk)
    ! call outfld ('XN2     ',xn2 ,pcols,lchnk)
    ! call outfld ('O2P     ',o2p ,pcols,lchnk)
    ! call outfld ('OP      ',op  ,pcols,lchnk)
    ! call outfld ('NOP     ',nop ,pcols,lchnk)
    qout(:ncol,:) = o2p(:ncol,:) + op(:ncol,:) + nop(:ncol,:)
    call outfld ('QIONSUM ', qout, pcols, lchnk)

    !-------------------------------------------------------------------------------
    ! calculate ExB drift velocities
    !-------------------------------------------------------------------------------
    call t_startf ( 'exbdrift' )
    call get_exbdrift( lchnk, ncol, pbuf)
    call t_stopf  ( 'exbdrift' )
    !-------------------------------------------------------------------------------

    do i = 1,ncol
       btesla(i) = bmag(i,lchnk)*1.e-9_r8 ! nT to teslas (see bmag in apex module)
       !-------------------------------------------------------------------------------
       ! Angular gyrofrequency of O+, O2+ and NO+ (s-1):
       !-------------------------------------------------------------------------------
       omega_op (i) = qe*btesla(i)*rmi_op_kg
       omega_o2p(i) = qe*btesla(i)*rmi_o2p_kg
       omega_nop(i) = qe*btesla(i)*rmi_nop_kg
       !-------------------------------------------------------------------------------
       ! Electron angular gyrofrequency (s-1):
       !-------------------------------------------------------------------------------
       omega_e(i) = qe*btesla(i)/emass 
       !-------------------------------------------------------------------------------
       ! Invert now, so we can multiply rather than divide in loops below:
       !-------------------------------------------------------------------------------
       omega_op_inv (i) = 1._r8/omega_op(i)
       omega_o2p_inv(i) = 1._r8/omega_o2p(i)
       omega_nop_inv(i) = 1._r8/omega_nop(i)
       omega_e_inv(i)   = 1._r8/omega_e(i)
       !-------------------------------------------------------------------------------
       ! Magnetic field geometry (used below in rotation of lambdas):
       !-------------------------------------------------------------------------------
       dipmag(i) = atan( bdown(i,lchnk)/sqrt(bnorth(i,lchnk)**2+beast(i,lchnk)**2) )
       decmag(i) = -atan2( beast(i,lchnk),bnorth(i,lchnk) )
       cosdec(i) = cos( decmag(i) )
       sindec(i) = sin( decmag(i) )
       if( abs(dipmag(i)) >= dipmin ) then
          sindip(i) = sin(dipmag(i))
          cosdip(i) = cos(dipmag(i))
       else
          if( dipmag(i) >= 0._r8 ) then
             sindip(i) = sin( dipmin )
             cosdip(i) = cos( dipmin )
          else
             sindip(i) = sin( -dipmin )
             cosdip(i) = cos( -dipmin )
          end if
       end if
       sin2dip(i) = sindip(i)**2
       cos2dip(i) = cosdip(i)**2
       sin2dec(i) = sindec(i)**2
       cos2dec(i) = cosdec(i)**2
    end do

    ! write(iulog,"('iondrag: btesla=',   /,(6e12.4))") btesla
    ! write(iulog,"('iondrag: bdown=',    /,(6e12.4))") bdown(:,lchnk)
    ! write(iulog,"('iondrag: beast=',    /,(6e12.4))") beast(:,lchnk)
    ! write(iulog,"('iondrag: bnorth=',   /,(6e12.4))") bnorth(:,lchnk)
    ! write(iulog,"('iondrag: omega_o2p=',/,(6e12.4))") omega_o2p
    ! write(iulog,"('iondrag: omega_op=' ,/,(6e12.4))") omega_op
    ! write(iulog,"('iondrag: omega_nop=',/,(6e12.4))") omega_nop

    !-------------------------------------------------------------------------------
    ! Ion-neutral momentum transfer collision frequency coefficients:
    !-------------------------------------------------------------------------------
    do k = 1,pver
       do i = 1,ncol
          tnti(i,k) = 0.5_r8*(ti(i,k) + tn(i,k))           ! ave of tn & ti
          sqrt_tnti = sqrt( tnti(i,k) )
          wrk       = log10( tnti(i,k) )
          !-------------------------------------------------------------------------------
          ! Collision frequency coefficients with O2 (cm3/s):
          !-------------------------------------------------------------------------------
          rnu_o2p_o2(i,k) = 2.59e-11_r8*sqrt_tnti &        ! O2+ ~ O2 (resonant)
               *(1._r8 - .073_r8*wrk)**2
          rnu_op_o2 (i,k) = 6.64e-10_r8                    ! O+  ~ O2
          rnu_nop_o2(i,k) = 4.27e-10_r8                    ! NO+ ~ O2
          !-------------------------------------------------------------------------------
          ! Collision frequency coefficients with O (cm3/s):
          !-------------------------------------------------------------------------------
          rnu_o2p_o(i,k) = 2.31e-10_r8                     ! O2+ ~ O
          rnu_op_o (i,k) = 3.67e-11_r8*sqrt_tnti  &        ! O+  ~ O (resonant)
               *(1._r8 - .064_r8*wrk)**2*colfac     
          rnu_nop_o(i,k) = 2.44e-10_r8                     ! NO+ ~ O
          !-------------------------------------------------------------------------------
          ! Collision frequency coefficients with N2 (cm3/s):
          !-------------------------------------------------------------------------------
          rnu_o2p_n2(i,k) = 4.13e-10_r8                    ! O2+ ~ N2
          rnu_op_n2 (i,k) = 6.82e-10_r8                    ! O+  ~ N2
          rnu_nop_n2(i,k) = 4.34e-10_r8                    ! NO+ ~ N2
       end do
    end do

    ! call outfld ('N_O2P_O2',rnu_o2p_o2,pcols,lchnk)
    ! call outfld ('N_OP_O2 ',rnu_op_o2 ,pcols,lchnk)
    ! call outfld ('N_NOP_O2',rnu_nop_o2,pcols,lchnk)

    ! call outfld ('N_O2P_O ',rnu_o2p_o ,pcols,lchnk)
    ! call outfld ('N_OP_O  ',rnu_op_o  ,pcols,lchnk)
    ! call outfld ('N_NOP_O ',rnu_nop_o ,pcols,lchnk)

    ! call outfld ('N_O2P_N2',rnu_o2p_n2,pcols,lchnk)
    ! call outfld ('N_OP_N2 ',rnu_op_n2 ,pcols,lchnk)
    ! call outfld ('N_NOP_N2',rnu_nop_n2,pcols,lchnk)
    !
    !-------------------------------------------------------------------------------
    ! Sub set_mean_mass (mo_mean_mass.F90) returns barm(ncol,pver) in g/mole,
    !   however, set_mean_mass sometimes returns zero in top(?) four values 
    !   of the column, so barm is calculated here, see below.
    !
    ! call set_mean_mass(ncol, state%q(1,1,imozart), barm)
    !
    ! Major species and ion number densities (mmr to cm-3):
    !-------------------------------------------------------------------------------

    call set_mean_mass( ncol, state%q(:,:,imozart:), barm )

    do k = 1,pver
       do i = 1,ncol
          press(i)     = 10._r8*state%pmid(i,k) ! from Pa to dyne/cm^2
          !     barm(i,k)   = 1._r8 / (xo2(i,k)*rmi_o2 + xo1(i,k)*rmi_o1 + xn2(i,k)*rmi_n2)
          xnmbar(i,k)  = press(i)*barm(i,k)/(boltzmann*tn(i,k)) 
          o2_cm3(i,k)  = xo2(i,k)*xnmbar(i,k)*rmi_o2    ! o2 (cm-3)
          o1_cm3(i,k)  = xo1(i,k)*xnmbar(i,k)*rmi_o1    ! o  (cm-3)
          n2_cm3(i,k)  = xn2(i,k)*xnmbar(i,k)*rmi_n2    ! n2 (cm-3)
          o2p_cm3(i,k) = o2p(i,k)*xnmbar(i,k)*rmi_o2p  ! o2+ (cm-3)
          op_cm3 (i,k) = op (i,k)*xnmbar(i,k)*rmi_op   ! o+  (cm-3)
          nop_cm3(i,k) = nop(i,k)*xnmbar(i,k)*rmi_nop  ! no+ (cm-3)
          ne(i,k)      = q_elec(i,k)*xnmbar(i,k)*rmi_elec  ! e (cm-3)
       end do
    end do

    ! call outfld ('O2_CM3  ',o2_cm3 ,pcols,lchnk)
    ! call outfld ('O1_CM3  ',o1_cm3 ,pcols,lchnk)
    ! call outfld ('N2_CM3  ',n2_cm3 ,pcols,lchnk)
    ! call outfld ('O2P_CM3 ',o2p_cm3,pcols,lchnk)
    ! call outfld ('OP_CM3  ',op_cm3 ,pcols,lchnk)
    ! call outfld ('NOP_CM3 ',nop_cm3,pcols,lchnk)

    !-------------------------------------------------------------------------------
    ! Multiply collision freq by neutral number density and sum for each ion:
    !
    ! rnu_o2p = [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
    ! rnu_op  = [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
    ! rnu_nop = [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
    !-------------------------------------------------------------------------------
    do k = 1,pver
       do i = 1,ncol
          rnu_o2p(i,k) = rnu_o2p_o2(i,k)*o2_cm3(i,k) &
               + rnu_o2p_o (i,k)*o1_cm3(i,k) &
               + rnu_o2p_n2(i,k)*n2_cm3(i,k)
          rnu_op (i,k) = rnu_op_o2 (i,k)*o2_cm3(i,k) &
               + rnu_op_o  (i,k)*o1_cm3(i,k) &
               + rnu_op_n2 (i,k)*n2_cm3(i,k)
          rnu_nop(i,k) = rnu_nop_o2(i,k)*o2_cm3(i,k) &
               + rnu_nop_o (i,k)*o1_cm3(i,k) &
               + rnu_nop_n2(i,k)*n2_cm3(i,k)
          !-------------------------------------------------------------------------------
          ! Electron collision frequency (s-1):
          !-------------------------------------------------------------------------------
          sqrt_te = sqrt(te(i,k))
          rnu_ne(i,k) = & 
               2.33e-11_r8*n2_cm3(i,k)*te(i,k)*(1._r8 - 1.21e-4_r8*te(i,k)) &
               + 1.82e-10_r8*o2_cm3(i,k)*sqrt_te*(1._r8 + 3.60e-2_r8*sqrt_te) &
               + 8.90e-11_r8*o1_cm3(i,k)*sqrt_te*(1._r8 + 5.70e-4_r8*te(i,k))
       end do
    end do

    ! call outfld ('COLL_O2P',rnu_o2p,pcols,lchnk)
    ! call outfld ('COLL_OP ',rnu_op ,pcols,lchnk)
    ! call outfld ('COLL_NOP',rnu_nop,pcols,lchnk)
    ! call outfld ('COLL_NE ',rnu_ne ,pcols,lchnk)

    !-------------------------------------------------------------------------------
    ! Ratio of collision to gyro frequencies for o2+, o+, no+, ne:
    !-------------------------------------------------------------------------------
    do k = 1,pver
       do i = 1,ncol
          rnu_o2p(i,k) = rnu_o2p(i,k)*omega_o2p_inv(i)
          rnu_op (i,k) = rnu_op (i,k)*omega_op_inv (i)
          rnu_nop(i,k) = rnu_nop(i,k)*omega_nop_inv(i)
          rnu_ne (i,k) = rnu_ne (i,k)*omega_e_inv  (i)
       end do
    end do

    ! call outfld ('RNU_O2P ',rnu_o2p,pcols,lchnk)
    ! call outfld ('RNU_OP  ',rnu_op ,pcols,lchnk)
    ! call outfld ('RNU_NOP ',rnu_nop,pcols,lchnk)
    ! call outfld ('RNU_NE  ',rnu_ne ,pcols,lchnk)

    !-------------------------------------------------------------------------------
    ! Calculate pedersen and Hall conductivities (siemens/m):
    !
    ! Qe_fac: 1.e6 to convert number densities from cm-3 to m-3:
    !-------------------------------------------------------------------------------
    qe_fac(:ncol) = qe*1.e6_r8/btesla(:ncol)

    do k = 1,pver
       do i = 1,ncol
          !-------------------------------------------------------------------------------
          ! ne = electron density assuming charge neutrality (cm-3):
          !-------------------------------------------------------------------------------
          !     ne(i,k) = op_cm3(i,k) + o2p_cm3(i,k) + nop_cm3(i,k)
          !-------------------------------------------------------------------------------
          ! Pedersen conductivity (siemens/m):
          !-------------------------------------------------------------------------------
          sigma_ped(i,k) = qe_fac(i)                           &
               *((op_cm3   (i,k)*rnu_op (i,k)/(1._r8 + rnu_op (i,k)**2)) &
               + (o2p_cm3(i,k)*rnu_o2p(i,k)/(1._r8 + rnu_o2p(i,k)**2)) &
               + (nop_cm3(i,k)*rnu_nop(i,k)/(1._r8 + rnu_nop(i,k)**2)) &
               + (ne     (i,k)*rnu_ne (i,k)/(1._r8 + rnu_ne (i,k)**2)))
          !-------------------------------------------------------------------------------
          ! Hall conductivity (siemens/m):
          !-------------------------------------------------------------------------------
          sigma_hall(i,k) = qe_fac(i)           &
               *(ne       (i,k)/(1._r8 + rnu_ne (i,k)**2) &
               - op_cm3 (i,k)/(1._r8 + rnu_op (i,k)**2) &
               - o2p_cm3(i,k)/(1._r8 + rnu_o2p(i,k)**2) &
               - nop_cm3(i,k)/(1._r8 + rnu_nop(i,k)**2))
       end do
    end do

    call outfld ('ELECDEN ',ne        ,pcols,lchnk)
    call outfld ('SIGMAPED',sigma_ped ,pcols,lchnk)
    call outfld ('SIGMAHAL',sigma_hall,pcols,lchnk)

    !--------------------------------------------------------------------------------------------
    !  Save conductivities in physics buffer using pointer for access in ionosphere module
    !--------------------------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') ) then 
      call pbuf_set_field(pbuf, PedConduct_idx,  sigma_ped(1:ncol,1:pver),  start=(/1,1/), kount=(/ncol,pver/) )
      call pbuf_set_field(pbuf, HallConduct_idx, sigma_hall(1:ncol,1:pver), start=(/1,1/), kount=(/ncol,pver/) )
    endif

    do k = 1,pver
       do i = 1,ncol
          wrk         = btesla(i)**2*avo*1.e-3_r8/xnmbar(i,k)
          lamda1(i,k) = sigma_ped(i,k)*wrk
          lamda2(i,k) = sigma_hall(i,k)*wrk
       end do
    end do

    call outfld ('LAMDA1',lamda1,pcols,lchnk)
    call outfld ('LAMDA2',lamda2,pcols,lchnk)

    do k = 1,pver
       do i = 1,ncol
          lxxnorot(i,k) = lamda1(i,k)
          lyynorot(i,k) = lamda1(i,k)*sin2dip(i)
          lxynorot(i,k) = lamda2(i,k)*sindip(i)
          lyxnorot(i,k) = lxynorot(i,k)
       end do
    end do

    !-------------------------------------------------------------------------------
    ! Rotate lambdas from local magnetic to geographic coordinates:
    !-------------------------------------------------------------------------------
    do k = 1,pver
       do i = 1,ncol
          lxx(i,k) = lxxnorot(i,k)*cos2dec(i) + lyynorot(i,k)*sin2dec(i)
          lyy(i,k) = lyynorot(i,k)*cos2dec(i) + lxxnorot(i,k)*sin2dec(i)
          wrk      = (lyynorot(i,k) - lxxnorot(i,k))*sindec(i)*cosdec(i)
          lyx(i,k) = lxynorot(i,k) - wrk
          lxy(i,k) = lxynorot(i,k) + wrk
       end do
    end do

    call outfld ('LXX     ',lxx,pcols,lchnk)
    call outfld ('LYY     ',lyy,pcols,lchnk)
    call outfld ('LXY     ',lxy,pcols,lchnk)
    call outfld ('LYX     ',lyx,pcols,lchnk)

    call physics_ptend_init(ptend, state%psetcols, "ion drag", lu=.true., lv=.true., ls=.true.)

    !-------------------------------------------------------------------------------
    ! Calculate ion drag tendencies and apply to neutral velocities:
    !-------------------------------------------------------------------------------
    call iondrag_tend( lchnk, ncol, state, ptend, pbuf,  &
         lxx, lyy, lxy, lyx, delt )

    !-------------------------------------------------------------------------------
    ! Calculate joule heating tendency and apply to temperature:
    !-------------------------------------------------------------------------------
    call jouleheat_tend( lchnk, ncol, state, ptend, pbuf,  &
         lxx, lyy, lxy, lyx )

  end subroutine iondrag_calc_ions

  !=========================================================================

  subroutine iondrag_calc_ghg (lchnk,ncol,state,ptend)

    use phys_grid,      only: get_rlat_all_p 
    use cam_history,    only: outfld
    use physics_types,  only: physics_ptend_init


    !
    !     This subroutine calculates ion drag using globally uniform
    !     ion drag tensor:
    !
    !                |alamxx       alamxy   | 
    !                |                      |
    !         lambda=|                      |
    !                |                      |
    !                |alamyx       alamyy   |
    !
    !     alamxx and alamxy are provided is data statements
    !     alamyy is obtaine from alamxx:
    !
    !
    !       alamyy = alamxx (sin(DIP_ANGLE))**2
    !
    !     where 
    !
    !       DIP_ANGLE = arctan(2.*tan(clat))
    !

    !--------------------Input arguments------------------------------------

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns

    type(physics_state), intent(in) :: state
    type(physics_ptend ), intent(out) :: ptend
    


    !---------------------Local workspace-------------------------------------

    real(r8) :: clat(pcols)             ! latitudes(radians) for columns

    real(r8) alamyyi                                 ! ALAMYY
    real(r8) dipan                                   ! dip angle
    real(r8) dut(pcols,pver)
    real(r8) dvt(pcols,pver)

    integer i
    integer k
    integer kinv

    !-------------------------------------------------------------------------

    if (.not.doiodrg) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    call physics_ptend_init(ptend, state%psetcols, "ion drag", lu=.true., lv=.true.)

    call get_rlat_all_p(lchnk, pcols, clat)

    !     calculate zonal wind drag
    dut(:,:)=0.0_r8
    do i=1,ncol
       do k=ntop_lev,nbot_lev
          dut(i,k)=-alamxyi(k)*state%v(i,k)-alamxxi(k)*state%u(i,k)
       enddo
    enddo

    !     calculate meridional wind drag
    dvt(:,:)=0.0_r8
    do i=1,ncol
       dipan=atan(2._r8*tan(clat(i)))
       do k=ntop_lev,nbot_lev
          alamyyi=alamxxi(k)*(sin(dipan))**2._r8
          dvt(i,k)=+alamxyi(k)*state%u(i,k)-alamyyi*state%v(i,k)
       enddo
    enddo

    do i=1,ncol
       do k=ntop_lev,nbot_lev
          ptend%u(i,k)=dut(i,k)
          ptend%v(i,k)=dvt(i,k)
       enddo
    enddo

    !  Write out tendencies
    call outfld('UIONTEND ',dut   ,pcols   ,lchnk  )
    call outfld('VIONTEND ',dvt   ,pcols   ,lchnk  )

    return
  end subroutine iondrag_calc_ghg

  !===================================================================================

  subroutine iondrag_tend( lchnk, ncol, state, ptend, pbuf,  &
       lxx, lyy, lxy, lyx, delt )

    !-------------------------------------------------------------------------------
    ! Calculate tendencies in U and V from ion drag tensors, which were
    !   calculated by sub iondrag_calc (module data lxx,lyy,lxy,lyx).
    ! This is called from sub iondrag_calc.
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    integer,intent(in)    :: lchnk                 ! current chunk index
    integer,intent(in)    :: ncol                  ! number of atmospheric columns
    real(r8), intent(in)  :: delt                  ! time step (s)
    real(r8), intent(in)  :: lxx(pcols,pver)       ! ion drag tensor
    real(r8), intent(in)  :: lyy(pcols,pver)       ! ion drag tensor
    real(r8), intent(in)  :: lxy(pcols,pver)       ! ion drag tensor
    real(r8), intent(in)  :: lyx(pcols,pver)       ! ion drag tensor
    type(physics_state), intent(in)    :: state ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend ! Physics tendencies
    
    type(physics_buffer_desc), pointer :: pbuf(:)


    !-------------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------------
    integer  :: i, k
    real(r8) :: dti
    real(r8) :: detr
    real(r8) :: us, vs
    real(r8) :: l11, l12, l21, l22
    real(r8) :: du(pcols,pver)             ! zonal ion drag tendency
    real(r8) :: dv(pcols,pver)             ! meridional ion drag tendency
    real(r8) :: dui(pcols,pver)            ! zonal ion drag tendency
    real(r8) :: dvi(pcols,pver)            ! meridional ion drag tendency
    real(r8), pointer :: ue(:)             ! pointer to pbuf
    real(r8), pointer :: ve(:)             ! pointer to pbuf

    !-------------------------------------------------------------------------------
    ! Get ion ExB drift from physics buffer (they were defined by exbdrift module)
    ! Ion drifts are 2d arrays, i.e., no vertical dimension.
    !-------------------------------------------------------------------------------

    call pbuf_get_field(pbuf, ue_idx, ue     )
    call pbuf_get_field(pbuf, ve_idx, ve     )

    dti = 1._r8/delt
    !-------------------------------------------------------------------------------
    ! Zonal (du) and meridional (dv) wind drag, using ExB drift velocities
    ! from exbdrift module (pbuf):
    !-------------------------------------------------------------------------------
    do k = ntop_lev,nbot_lev
       do i = 1,ncol
          !-------------------------------------------------------------------------------
          ! 2/28/04 btf:
          ! Full ion-drag, using lambdas and ExB drifts. 
          !   This should succeed with bz = 0 (efield module)
          ! Runs:
          !   bz=-5, nstep=24 min, nsplit=4 (6 min dynamics): crashed after 2 days.
          !   bz=-5, nstep=24 min, nsplit=6 (4 min dynamics): 5 day run succeeded.
          ! See comments in efield module re bz < 0 (efield.F90).
          !-------------------------------------------------------------------------------
          us      = ue(i) - state%u(i,k)
          vs      = ve(i) - state%v(i,k)
          !-------------------------------------------------------------------------------
          ! Exclude ue,ve drift momentum source to avoid crashes when bz < 0 and 
          ! full 30 min timestep (partial ion-drag):
          !-------------------------------------------------------------------------------
          l11     = dti + lxx(i,k)
          l12     = lxy(i,k)
          l21     = -lyx(i,k)
          l22     = dti + lyy(i,k)
          detr    = dti/(l11*l22 - l12*l21)
          dui(i,k) = dti*(detr*(l12*vs - l22*us) + us)
          dvi(i,k) = dti*(detr*(l21*us - l11*vs) + vs)
       end do
    end do

    !-------------------------------------------------------------------------------
    ! Apply to model tendencies:
    !-------------------------------------------------------------------------------
    do k = ntop_lev,nbot_lev
       !-------------------------------------------------------------------------------
       ! Ion drag tendencies:
       !-------------------------------------------------------------------------------
       ptend%u(:ncol,k) = dui(:ncol,k)
       ptend%v(:ncol,k) = dvi(:ncol,k)
       !-------------------------------------------------------------------------------
       ! Turn off ion drag tendency:
       !-------------------------------------------------------------------------------
       !     ptend%u(:ncol,k) = 0._r8
       !     ptend%v(:ncol,k) = 0._r8
    end do
    do k = nbot_lev+1,pver
       dui(:ncol,k)      = 0._r8
       dvi(:ncol,k)      = 0._r8
       ptend%u(:ncol,k) = 0._r8
       ptend%v(:ncol,k) = 0._r8
    end do

    call outfld ( 'UIONTEND', dui, pcols, lchnk )      ! u ion drag tendency
    call outfld ( 'VIONTEND', dvi, pcols, lchnk )      ! v ion drag tendency

  end subroutine iondrag_tend

  !================================================================================================
  subroutine jouleheat_tend( lchnk, ncol, state, ptend, pbuf,  &
       lxx, lyy, lxy, lyx )
    !-------------------------------------------------------------------------------
    ! Calculate tendencies in T due to joule heating.
    ! This is called from sub iondrag_calc.
    !-------------------------------------------------------------------------------

    use physconst, only: cpair,pi,cpairv
    use phys_grid, only: get_rlon_p, get_rlat_p

    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    integer,intent(in)    :: lchnk                    ! current chunk index
    integer,intent(in)    :: ncol                     ! number of atmospheric columns
    real(r8), intent(in)  :: lxx(pcols,pver)          ! ion drag tensor
    real(r8), intent(in)  :: lyy(pcols,pver)          ! ion drag tensor
    real(r8), intent(in)  :: lxy(pcols,pver)          ! ion drag tensor
    real(r8), intent(in)  :: lyx(pcols,pver)          ! ion drag tensor
    type(physics_state), intent(in)    :: state       ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend       ! Physics tendencies (inout)
    
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-------------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------------
    integer  :: k, i
    integer  :: ndx_ue, ndx_ve             ! indices to ion velocities in pbuf
    integer  :: max_ind(2)
    real(r8) :: us, vs
    real(r8) :: max_q
    real(r8) :: qjoule(pcols,pver)         ! joule heating
    real(r8) :: qout(pcols,pver)           ! temp for outfld
    real(r8), pointer :: ue(:)             ! pointer to pbuf
    real(r8), pointer :: ve(:)             ! pointer to pbuf

    !-------------------------------------------------------------------------------
    ! Get ion velocities from physics buffer (they were defined by exbdrift module)
    ! Ion velocities are 2d arrays, i.e., no vertical dimension.
    !-------------------------------------------------------------------------------
    call pbuf_get_field(pbuf, ue_idx, ue     )
    call pbuf_get_field(pbuf, ve_idx, ve     )

    ! write(iulog,"(/,'qjoule: lchnk=',i5,' pcols=',i4,' ncol=',i4)") lchnk,pcols,ncol
    ! write(iulog,"('qjoule: ue=',/,(6e12.4))") ue
    ! write(iulog,"('qjoule: ve=',/,(6e12.4))") ve

    do k = ntop_lev,nbot_lev
       !   write(iulog,"('qjoule: k=',i3,' u=',/,(6e12.4))") k,state%u(:,k)
       !   write(iulog,"('qjoule: k=',i3,' v=',/,(6e12.4))") k,state%v(:,k)
       !   write(iulog,"('qjoule: k=',i3,' lxx=',/,(6e12.4))") k,lxx(:,k)
       !   write(iulog,"('qjoule: k=',i3,' lxy=',/,(6e12.4))") k,lxy(:,k)
       !   write(iulog,"('qjoule: k=',i3,' lyx=',/,(6e12.4))") k,lyx(:,k)
       !   write(iulog,"('qjoule: k=',i3,' lyy=',/,(6e12.4))") k,lyy(:,k)
       do i = 1,ncol
          us           = ue(i) - state%u(i,k)
          vs           = ve(i) - state%v(i,k)
          qjoule(i,k)  = us*us*lxx(i,k) + us*vs*(lxy(i,k) - lyx(i,k)) + vs*vs*lyy(i,k)
          ptend%s(i,k) = qjoule(i,k)        ! joule heating tendency
          !     ptend%s(i,k) = 0._r8              ! no joule heating tendency
       end do
       !   write(iulog,"('qjoule: k=',i3,' qjoule(:,k)=',/,(6e12.4))") k,qjoule(:,k)
    end do
    do k = nbot_lev+1,pver
       qjoule(:ncol,k)  = 0._r8
       ptend%s(:ncol,k) = 0._r8        ! no joule heating tendency
    end do

#ifdef SW_DEBUG
    max_q      = 100._r8*maxval( abs( qjoule(:ncol,ntop_lev:nbot_lev) )/state%t(:ncol,ntop_lev:nbot_lev) )
    max_ind(:) = maxloc( abs( qjoule(:ncol,ntop_lev:nbot_lev) )/state%t(:ncol,ntop_lev:nbot_lev) )
    i = max_ind(1)
    k = max_ind(2)
    if( lchnk == 25 ) then
       i = 14
       k = 3
       write(iulog,*) ' '
       write(iulog,*) '-------------------------------------------------------'
       write(iulog,*) 'jouleheat_tend: lon,lat = ',get_rlon_p(lchnk,14)*180._r8/pi, get_rlat_p(lchnk,14)*180._r8/pi
       write(iulog,*) 'jouleheat_tend: dt,t,max% dt/t = ',qjoule(i,k)/cpairv(i,k,lchnk),state%t(i,k),max_q, &
            ' @ lchnk,i,k = ',lchnk,max_ind(:)
       write(iulog,*) 'jouleheat_tend: lxx,xy,yx,yy   = ',lxx(i,k),lxy(i,k),lyx(i,k),lyy(i,k)
       write(iulog,*) 'jouleheat_tend: u,ue,v,ve      = ',state%u(i,k),ue(i),state%v(i,k),ve(i)
       write(iulog,*) 'jouleheat_tend: us,vs          = ',ue(i) - state%u(i,k),ve(i) - state%v(i,k)
       write(iulog,*) 'jouleheat_tend: du,dv          = ',ptend%u(i,k),ptend%v(i,k)
       write(iulog,*) 'jouleheat_tend: dt'
       write(iulog,'(1p,5g15.7)') qjoule(max_ind(1),ntop_lev:nbot_lev)/cpairv(max_ind(1),ntop_lev:nbot_lev,lchnk)
       write(iulog,*) '-------------------------------------------------------'
       write(iulog,*) ' '
       ! stop 'diagnostics'
    end if
#endif

    qout(:ncol,:) = qjoule(:ncol,:)/cpairv(:ncol,:,lchnk)

    call outfld ( 'QJOULE', qout, pcols, lchnk )

  end subroutine jouleheat_tend

end module iondrag
