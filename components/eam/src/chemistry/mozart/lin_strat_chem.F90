module lin_strat_chem
!     24 Oct 2008 -- Francis Vitt
!      9 Dec 2008 -- Philip Cameron-Smith, LLNL, -- added ltrop
!      4 Jul 2019 -- Qi Tang (LLNL), Juno Hsu (UCI), -- added sfcsink
!      20 Jan 2021 -- Juno Hsu(UCI) added Linoz v3 
!        --new linoz v3 subroutine linv3_strat_chem_solve (linoz_v3 if O3LNZ, N2OLNZ, NOYLNZ and CH4LNZ are defined)
!        --modified lin_strat_chem_solve is now linv2_strat_chem_solve (linoz_v2 if only O3LNZ is defined; use only O3LNZ part of Linoz v3 netcdf file)
!        --modified sfcsink working for v2 or v3 depending on species
!      Spring 2021 -- Juno Hsu
!       --added H2OLNZ, saves H2O water vapor from CH4 oxidation
!       --added prescribed O3LBS in netcdf file (prescribed CMIP6 historical surface ozone below 925 mb). And Linoz surface ozone are relaxed to this profile in 2-days within the last 9 surface layers
!      added 30-day e-fold decay subroutine (lin_strat_efold_decay) for ChemUCI tropospheric species assigned to the namelist, fstrat_efold_list   
!     20 Sep 2021 -- Qi Tang (LLNL), -- added 3D tropopause
!     04 Nov 2021 -- Hsiang-He Lee (LLNL) -- Added flags for two sets of Linoz chemistry to run either O3/CH4x/N2Ox/NOYx set or 
!                    O3LNZ/N2OLNZ/NOYLNZ/CH4LNZ to get LINOZ tendecies 


  use shr_kind_mod , only : r8 => shr_kind_r8
  use ppgrid       , only : begchunk, endchunk
  use physics_types, only : physics_state
  use cam_logfile  , only : iulog
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use chem_mods,        only: gas_pcnst
  !
  implicit none
  !
  private  ! all unless made public

  save
  !
  ! define public components of module
  !
  public :: lin_strat_chem_inti, linv2_strat_chem_solve, linv3_strat_chem_solve
  public :: lin_strat_sfcsink
  public :: do_lin_strat_chem, linoz_v2, linoz_v3
  public :: linoz_readnl   ! read linoz_nl namelist
  public :: has_fstrat_efold, fstrat_efold_inti, fstrat_efold_decay

  integer :: o3lnz_ndx, n2olnz_ndx, noylnz_ndx, ch4lnz_ndx, h2olnz_ndx
  integer :: o3_ndx, n2o_ndx, ch4_ndx, no_ndx, no2_ndx, hno3_ndx 
  integer :: uci1_ndx
 
  logical :: do_lin_strat_chem, linoz_v2, linoz_v3

  real(r8), parameter :: unset_r8   = huge(1.0_r8)
  integer , parameter :: unset_int  = huge(1)
  integer  :: linoz_lbl = unset_int ! number of layers with ozone decay from the surface
  real(r8) :: linoz_sfc = unset_r8  ! boundary layer concentration (ppb) to which Linoz ozone e-fold
  real(r8) :: linoz_tau = unset_r8  ! Linoz e-fold time scale (in seconds) in the boundary layer
  real(r8) :: linoz_psc_T = unset_r8  ! PSC ozone loss T (K) threshold

  integer  :: o3_lbl ! set from namelist input linoz_lbl
  real(r8) :: o3_sfc ! set from namelist input linoz_sfc
  real(r8) :: o3_tau ! set from namelist input linoz_tau
  real(r8) :: psc_T  ! set from namelist input linoz_psc_T
  
  logical :: has_fstrat_efold(gas_pcnst) 

contains

subroutine linoz_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'linoz_readnl'

   namelist /linoz_nl/ linoz_lbl, linoz_sfc, linoz_tau, linoz_psc_T
   !-----------------------------------------------------------------------------

   ! Set default values
   linoz_lbl    = 4
   linoz_sfc    = 30.0e-9_r8
   linoz_tau    = 172800.0_r8
   linoz_psc_T  = 197.5_r8

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'linoz_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, linoz_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      o3_lbl = linoz_lbl
      o3_sfc = linoz_sfc
      o3_tau = linoz_tau
      psc_T  = linoz_psc_T

      ! check
      write(iulog,*) subname // ', linoz_lbl:',   o3_lbl
      write(iulog,*) subname // ', linoz_sfc:',   o3_sfc
      write(iulog,*) subname // ', linoz_tau:',   o3_tau
      write(iulog,*) subname // ', linoz_psc_T:', psc_T

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(o3_lbl,            1, mpiint, 0, mpicom)
   call mpibcast(o3_sfc,            1, mpir8,  0, mpicom)
   call mpibcast(o3_tau,            1, mpir8,  0, mpicom)
   call mpibcast(psc_T,             1, mpir8,  0, mpicom)
#endif

end subroutine linoz_readnl

!--------------------------------------------------------------------
!--------------------------------------------------------------------
    subroutine lin_strat_chem_inti(phys_state)
    !
    ! initialize linearized stratospheric chemistry by reading
    ! input parameters from netcdf file and interpolate to
    ! present model grid
    !
    use linoz_data,   only : linoz_data_init, has_linozv3_data, has_linoz_data
    use ppgrid,       only : pver
    use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
    use cam_history,  only : addfld, horiz_only, add_default
    use physics_buffer, only : physics_buffer_desc

    implicit none


    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    if (.not. has_linoz_data) return
    
    do_lin_strat_chem = has_linoz_data
    
    ! check for synoz

    if( get_spc_ndx( 'SYNOZ' ) > 0 .and. has_linozv3_data) then
       call endrun('lin_strat_chem_inti: cannot have both synoz and linoz')
    endif

    ! linoz_v3 species
 !if not exist, get_spc_ndx return -1 in mo_chem_utls.F90
    linoz_v3 = .false.
    linoz_v2 = .false.

     o3lnz_ndx  =   get_spc_ndx('O3LNZ')
    n2olnz_ndx  =   get_spc_ndx('N2OLNZ')
    noylnz_ndx  =   get_spc_ndx('NOYLNZ')
    ch4lnz_ndx  =   get_spc_ndx('CH4LNZ')
    h2olnz_ndx  =   get_spc_ndx('H2OLNZ')

    uci1_ndx    = get_rxt_ndx('uci1')

    if (uci1_ndx <=0 ) then
       !write(iulog,*) 'skip Linoz, need to have tracer O3LNZ at least '
       write(iulog,*) 'skip Linoz, temporally change '
       do_lin_strat_chem = .false.
       return
    end if

    !linoz_v3= (o3lnz_ndx > 0 .and. n2olnz_ndx >0  .and. noylnz_ndx >0  .and. ch4lnz_ndx > 0)
    !linoz_v2= (o3lnz_ndx > 0 .and. n2olnz_ndx <0  .and. noylnz_ndx <0  .and. ch4lnz_ndx < 0)
    linoz_v3= (n2olnz_ndx >0  .and. noylnz_ndx >0 .and. ch4lnz_ndx > 0)
    linoz_v2= (n2olnz_ndx <0  .and. noylnz_ndx <0 .and. ch4lnz_ndx < 0)
!    write(iulog,*)'linoz_v3=',linoz_v3,'linoz_v2=',linoz_v2
! real o3, ch4, n2o tracers
    o3_ndx   =   get_spc_ndx('O3')
    ch4_ndx  =   get_spc_ndx('CH4')
    n2o_ndx  =   get_spc_ndx('N2O')
    no_ndx   =   get_spc_ndx('NO')
    no2_ndx  =   get_spc_ndx('NO2')
    hno3_ndx =   get_spc_ndx('HNO3')

   if(masterproc)then
   ! write(iulog,*)'in subroutine lin_strat_chem_inti'
   ! write(iulog,*)'O3LNZ=',o3lnz_ndx,'N2OLNZ=',n2olnz_ndx,'NOYLNZ=',noylnz_ndx,'H2OLNZ=',h2olnz_ndx
   ! write(iulog,*)'O3=',o3_ndx,'CH4=',ch4_ndx,'N2O=',n2o_ndx,'NO=',no_ndx,'NO2=',no2_ndx,'HNO3=',hno3_ndx
   !  write(iulog,*)'inside lin_strat_solve for ndx o3, o3lnz, n2o, noylnz, ch4', o3_ndx, o3lnz_ndx, n2o_ndx, noylnz_ndx, ch4_ndx
   end if
    
!     
!   initialize the linoz data

    call linoz_data_init()

    ! define additional output

    if (o3lnz_ndx >0) call addfld( 'LINOZ_DO3LNZ'    , (/ 'lev' /), 'A', '/s'     , 'O3LNZ vmr tendency by linearized ozone chemistry'   )
    if (o3lnz_ndx >0) call addfld( 'LINOZ_DO3LNZ_PSC', (/ 'lev' /), 'A', '/s'     , 'O3LNZ vmr loss by PSCs using Carille et al. (1990)' )
    if (o3lnz_ndx >0) call addfld( 'LINOZ_2DDO3LNZ'    , horiz_only, 'A', 'DU/s'     , 'O3LNZ 2D DU tendency by linearized ozone chemistry'   )
    if (o3lnz_ndx >0) call addfld( 'LINOZ_2DDO3LNZ_PSC', horiz_only, 'A', 'DU/s'     , 'O3LNZ 2D DU loss by PSCs using Carille et al. (1990)' )
    call addfld( 'LINOZ_SSO3'   , (/ 'lev' /), 'A', 'kg'     , 'steady state O3 in LINOZ'                        )
    call addfld( 'LINOZ_O3COL'  , (/ 'lev' /), 'A', 'DU'     , 'ozone column above'                                 )
    call addfld( 'LINOZ_O3CLIM' , (/ 'lev' /), 'A', 'mol/mol', 'climatology of ozone in LINOZ'                      )
    call addfld( 'LINOZ_SZA'    ,    horiz_only, 'A', 'degrees', 'solar zenith angle in LINOZ'                      )
    call addfld( 'LINOZ_O3SFCSINK',  horiz_only, 'A', 'Tg/yr/m2'   ,   'surface ozone sink in LINOZ with an e-fold to a fixed concentration' )
    call addfld( 'LINOZ_DO3'    , (/ 'lev' /), 'A', '/s'     , 'O3 vmr tendency by linearized ozone chemistry'   )
    call addfld( 'LINOZ_DO3_PSC', (/ 'lev' /), 'A', '/s'     , 'O3 vmr loss by PSCs using Carille et al. (1990)' )
    call addfld( 'LINOZ_2DDO3'  , horiz_only, 'A', 'DU/s'     , 'O3 2D DU tendency by linearized ozone chemistry'   )
    call addfld( 'LINOZ_2DDO3_PSC', horiz_only, 'A', 'DU/s'     , 'O3 2D DU loss by PSCs using Carille et al. (1990)' )
    if (noylnz_ndx >0) call addfld( 'LINOZ_NOYSFCSINK', horiz_only, 'A', 'Tg/yr/m2'   , 'surface noylnz sink in LINOZ v3 with an e-fold to a fixed concentration' )
    if (n2olnz_ndx >0) call addfld( 'LINOZ_N2OSFCSRC',  horiz_only, 'A', 'Tg/yr/m2'   , 'surface n2o source in LINOZ v3 with an e-fold to a fixed concentration' )
    if (ch4lnz_ndx >0) call addfld( 'LINOZ_CH4SFCSRC',  horiz_only, 'A', 'Tg/yr/m2'   , 'surface n2o source in LINOZ v3 with an e-fold to a fixed concentration' )

    if (o3lnz_ndx >0) call add_default( 'LINOZ_DO3LNZ'    , 1, ' ' )
    if (o3lnz_ndx >0) call add_default( 'LINOZ_DO3LNZ_PSC', 1, ' ' )
    if (o3lnz_ndx >0) call add_default( 'LINOZ_2DDO3LNZ'    , 1, ' ' )
    if (o3lnz_ndx >0) call add_default( 'LINOZ_2DDO3LNZ_PSC', 1, ' ' )
    call add_default( 'LINOZ_SSO3'   , 1, ' ' )
    call add_default( 'LINOZ_O3COL'  , 1, ' ' )
    call add_default( 'LINOZ_O3CLIM' , 1, ' ' )
    call add_default( 'LINOZ_SZA'    , 1, ' ' )
    call add_default( 'LINOZ_O3SFCSINK', 1, ' ' )
    call add_default( 'LINOZ_2DDO3'    , 1, ' ' )
    call add_default( 'LINOZ_2DDO3_PSC', 1, ' ' )
    if (noylnz_ndx >0) call add_default( 'LINOZ_NOYSFCSINK', 1, ' ' )
    if (n2olnz_ndx >0) call add_default( 'LINOZ_N2OSFCSRC', 1, ' ' )
    if (ch4lnz_ndx >0) call add_default( 'LINOZ_CH4SFCSRC', 1, ' ' )
    return
  end subroutine lin_strat_chem_inti


!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine linv3_strat_chem_solve(ncol, lchnk, xvmr, h2ovmr, xsfc, o3col,temp, sza, pmid, delta_t, rlats, ltrop, pdeldry, tropFlag)

!--------------------------------------------------------------------
! linearized ozone chemistry linoz-v3 (o3-n2o-noy-ch4 prognostic equations, h2o diagnosed from ch4) 
! from Hsu and Prather, grl, 2009   (https://doi.org/10.1029/2009GL042243)
! 
!written by Juno Hsu (junoh@uci.edu), 09/2020 
 
    use chlorine_loading_data, only: chlorine_loading
    use chem_mods,             only: gas_pcnst
    use mo_chem_utls,          only: get_spc_ndx
    use mo_constants,          only: pi, rgrav, rearth, avogadro
    !
    ! this subroutine updates the ozone mixing ratio in the stratosphere
    ! using linearized chemistry 
    !
    use ppgrid,        only : pcols, pver
    use physconst,     only : pi, &
                              grav => gravit, &
                              mw_air => mwdry
    use cam_history,   only : outfld
    use linoz_data,  only : fields, o3_clim_ndx, n2o_clim_ndx, noy_clim_ndx, ch4_clim_ndx, h2o_clim_ndx, t_clim_ndx,o3col_clim_ndx,&
                              no3_PmL_clim_ndx,  no3_dPmL_dO3_ndx,   no3_dPmL_dN2O_ndx, no3_dPmL_dNOy_ndx,   &
                              no3_dPmL_dCH4_ndx, no3_dPmL_dH2O_ndx,  no3_dPmL_dT_ndx,   no3_dPmL_dO3col_ndx, &
                              pn2o_PmL_clim_ndx, pn2o_dPmL_dO3_ndx,  pn2o_dPmL_dN2O_ndx,pn2o_dPmL_dNOy_ndx,  &
                              pn2o_dPmL_dCH4_ndx,pn2o_dPmL_dH2O_ndx, pn2o_dPmL_dT_ndx,  pn2o_dPmL_dO3col_ndx,&
                              ln2o_PmL_clim_ndx, ln2o_dPmL_dO3_ndx,  ln2o_dPmL_dN2O_ndx,ln2o_dPmL_dNOy_ndx,  &
                              ln2o_dPmL_dCH4_ndx,ln2o_dPmL_dH2O_ndx, ln2o_dPmL_dT_ndx,  ln2o_dPmL_dO3col_ndx,&
                              pnoy_PmL_clim_ndx, pnoy_dPmL_dO3_ndx,  pnoy_dPmL_dN2O_ndx,pnoy_dPmL_dNOy_ndx,  &
                              pnoy_dPmL_dCH4_ndx,pnoy_dPmL_dH2O_ndx, pnoy_dPmL_dT_ndx,  pnoy_dPmL_dO3col_ndx,&
                              lnoy_PmL_clim_ndx, lnoy_dPmL_dO3_ndx,  lnoy_dPmL_dN2O_ndx,lnoy_dPmL_dNOy_ndx,  &
                              lnoy_dPmL_dCH4_ndx,lnoy_dPmL_dH2O_ndx, lnoy_dPmL_dT_ndx,  lnoy_dPmL_dO3col_ndx,&
                              nch4_PmL_clim_ndx, nch4_dPmL_dO3_ndx,  nch4_dPmL_dN2O_ndx,nch4_dPmL_dNOy_ndx,  &
                              nch4_dPmL_dCH4_ndx,nch4_dPmL_dH2O_ndx, nch4_dPmL_dT_ndx,  nch4_dPmL_dO3col_ndx,&
                              cariolle_pscs_ndx, o3lbs_ndx
    !
    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index
    real(r8), intent(inout), dimension(ncol ,pver,gas_pcnst) :: xvmr      ! volume mixing ratio for all
    real(r8), intent(in)   , dimension(ncol ,pver) :: h2ovmr              ! h2o vapor volume mixing ratio 
    real(r8), intent(inout)                        :: xsfc(4,ncol)        ! surface o3,n2o,noy,ch4 from linoz table
    real(r8), intent(in)   , dimension(ncol ,pver) :: o3col               ! ozone column above box, can be either "O3" column or "O3LNZ" column (mol/cm^2)
    real(r8), intent(in)   , dimension(pcols,pver) :: temp                ! temperature (K)
    real(r8), intent(in)   , dimension(ncol )      :: sza                 ! local solar zenith angle
    real(r8), intent(in)   , dimension(pcols,pver) :: pmid                ! midpoint pressure (Pa)
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)
    real(r8), intent(in)                           :: rlats(ncol)         ! column latitudes (radians)
    integer,  intent(in)   , dimension(pcols)      :: ltrop               ! chunk index    
    real(r8), intent(in)   , dimension(ncol ,pver) :: pdeldry             !  dry pressure delta about midpoints (Pa) 
    logical, optional, intent(in)                  :: tropFlag(pcols,pver)! 3D tropospheric level flag
    !
    integer  :: i,k,n,ll,lt0,lt, n_dl !,index_lat,index_month
    integer :: kmax
    real(r8) :: o3col_du,delta_temp,delta_o3col    
    real(r8) ::o3_old,   n2o_old,  noy_old,  ch4_old,  h2o_old
    real(r8) ::o3_new,   n2o_new,  noy_new,  ch4_new,  h2o_new
    real(r8) ::o3_clim, ss_x
    real(r8) :: dn2op, dn2ol, dnoyp, dnoyl, delo3, delo3_psc, delch4, lfreq
    real(r8) :: max_sza, psc_loss, ch4max, pw
    real(r8), dimension(ncol) :: lats
    real(r8), dimension(ncol,pver) :: do3_linoz, do3_linoz_psc, ss_o3, o3col_du_diag, o3clim_linoz_diag
    ! local
    real(r8), dimension(ncol,pver) :: o3_vmr, n2o_vmr, noy_vmr, ch4_vmr, h2o_vmr
    real(r8), dimension(ncol,pver) :: dO3, dN2O, dNOY, dCH4, dH2O, dTemp, dCOL
    real(r8), dimension(ncol,pver) :: PL_O3, P_n2o, Lfreq_n2o, Pfreq_noy, Lfreq_noy, Lfreq_ch4
    real(r8), dimension(:,:), pointer :: linoz_o3_clim
    real(r8), dimension(:,:), pointer :: linoz_n2o_clim
    real(r8), dimension(:,:), pointer :: linoz_noy_clim
    real(r8), dimension(:,:), pointer :: linoz_ch4_clim
    real(r8), dimension(:,:), pointer :: linoz_h2o_clim
    real(r8), dimension(:,:), pointer :: linoz_t_clim
    real(r8), dimension(:,:), pointer :: linoz_o3col_clim

    real(r8), dimension(:,:), pointer :: linoz_PmL_clim
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3X
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dN2O
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dNOY
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dCH4
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dH2O
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dT
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3col
    real(r8), dimension(:,:), pointer :: linoz_cariolle_psc
    real(r8), dimension(:,:), pointer :: linoz_o3lbs
    ! real O3 variables
    real(r8), dimension(ncol,pver) :: do3_linoz_du, do3_linoz_psc_du
    real(r8), dimension(ncol) :: twod_do3_linoz
    real(r8), dimension(ncol) :: twod_do3_linoz_psc
    !
    ! parameters
    !
    real(r8), parameter :: convert_to_du = 1._r8/(2.687e16_r8)      ! convert ozone column from mol/cm^2 to DU
    real(r8), parameter :: degrees_to_radians = pi/180._r8          ! conversion factors
    real(r8), parameter :: radians_to_degrees = 180._r8/pi
    real(r8), parameter :: chlorine_loading_1987    = 2.5977_r8     ! EESC value (ppbv)
    real(r8), parameter :: chlorine_loading_bgnd    = 0.0000_r8     ! EESC value (ppbv) for background conditions
    real(r8), parameter :: small    = 1.e-18_r8     ! prevent dividing by zero
    real(r8), parameter :: pressure_threshold       = 237.1375e+2_r8    ! last vaild linoz pressure layer
    !
    if ( .not. do_lin_strat_chem ) return
    if ( .not. linoz_v3) return

!    write(iulog,*)'inside lin_strat_solve for ndx o3, o3lnz, n2o, noylnz, ch4', o3_ndx, o3lnz_ndx, n2o_ndx, noylnz_ndx, ch4_ndx
    !o3_vmr =  xvmr(:,:, o3lnz_ndx)
    o3_vmr =  xvmr(:,:, o3_ndx)
    n2o_vmr =  xvmr(:,:,n2olnz_ndx)
    noy_vmr =  xvmr(:,:,noylnz_ndx)
    ch4_vmr =  xvmr(:,:,ch4lnz_ndx)
    h2o_vmr =  xvmr(:,:,h2olnz_ndx)
 
    ! associate the field pointers
    !
    !Linoz climatological data 
       linoz_o3_clim      => fields(o3_clim_ndx)      %data(:,:,lchnk )
       linoz_n2o_clim     => fields(n2o_clim_ndx)     %data(:,:,lchnk )
       linoz_noy_clim     => fields(noy_clim_ndx)     %data(:,:,lchnk )
       linoz_ch4_clim     => fields(ch4_clim_ndx)     %data(:,:,lchnk )
       linoz_h2o_clim     => fields(h2o_clim_ndx)     %data(:,:,lchnk )
       linoz_t_clim       => fields(t_clim_ndx)       %data(:,:,lchnk )
       linoz_o3col_clim   => fields(o3col_clim_ndx)   %data(:,:,lchnk )
       linoz_o3lbs        => fields(o3lbs_ndx)        %data(:,:,lchnk )

       dO3(:,:)     =   o3_vmr(:,:)  - linoz_o3_clim(:,:)
       dN2O(:,:)    =  n2o_vmr(:,:)  - linoz_n2o_clim(:,:)
       dNOY(:,:)    =  noy_vmr(:,:)  - linoz_noy_clim(:,:) 
       dCH4(:,:)    =  ch4_vmr(:,:)  - linoz_ch4_clim(:,:)
       dH2O(:,:)    =  h2o_vmr(:,:)  - linoz_h2o_clim(:,:) 
       dTemp(:,:)   =  temp(:,:)     - linoz_t_clim(:,:)
       dCOL(:,:)     =  o3col(:,:)*convert_to_du - linoz_o3col_clim(:,:)

! potential water 2*ch4 +h2o !it might be better to get maxch4 3-4 years back in time but this will do for now 
! upated yearly o3,n2o,noy and ch4 value are stored at xx_clim in linoz file at the bottom padding layer
! output surface constant concentration to control surface sink/source
!
!       write(iulog,*)'pver=',pver

       xsfc(1,:)=   linoz_o3lbs(:,pver) !  ozone surface constant (varying in lat)
       xsfc(2,:)=   linoz_n2o_clim(:,pver) !  n2o (constant throughout latitude)
       xsfc(3,:)=   linoz_o3lbs(:,pver)*3.e-3_r8  !noylnz
       xsfc(4,:)=   linoz_ch4_clim(:,pver) ! ch4 (constant throughout latitude)
       ch4max =     maxval(linoz_ch4_clim(1:ncol,pver)) 
       pw= 2.0_r8 * ch4max + 3.65e-6_r8 
 
! OZONE P-L terms !unit vmr/sec
       linoz_PmL_clim     => fields(no3_PmL_clim_ndx)     %data(:,:,lchnk )    !unit vmr/sec
       linoz_dPmL_dO3X    => fields(no3_dPmL_dO3_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dn2o    => fields(no3_dPmL_dN2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dnoy    => fields(no3_dPmL_dNOY_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dch4    => fields(no3_dPmL_dCH4_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dh2o    => fields(no3_dPmL_dH2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dT      => fields(no3_dPmL_dT_ndx)      %data(:,:,lchnk )
       linoz_dPmL_dO3col  => fields(no3_dPmL_dO3col_ndx)  %data(:,:,lchnk )

!for steady-state sol. plug-in, no need for dPmL_dO3 term 
       PL_O3 =  linoz_PmL_clim             &       
            +  linoz_dPmL_dn2o    * dN2O  &
            +  linoz_dPmL_dnoy    * dNOY  &
            +  linoz_dPmL_dch4    * dCH4  &
            +  linoz_dPmL_dh2o    * dH2O  &
            +  linoz_dPmL_dT      * dTemp &
            +  linoz_dPmL_dO3col  * dCOL
! 
!pointer to pn2o
       linoz_PmL_clim     => fields(pn2o_PmL_clim_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dO3     => fields(pn2o_dPmL_dO3_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dn2o    => fields(pn2o_dPmL_dN2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dnoy    => fields(pn2o_dPmL_dNOY_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dch4    => fields(pn2o_dPmL_dCH4_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dh2o    => fields(pn2o_dPmL_dH2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dT      => fields(pn2o_dPmL_dT_ndx)      %data(:,:,lchnk )
       linoz_dPmL_dO3col  => fields(pn2o_dPmL_dO3col_ndx)  %data(:,:,lchnk )
!Taylor-expansion of P (mr/sec)
       P_n2o =  linoz_PmL_clim           &   
            +  linoz_dPmL_dO3    * dO3      &
            +  linoz_dPmL_dn2o   * dN2O     &
            +  linoz_dPmL_dnoy   * dNOY     &
            +  linoz_dPmL_dch4   * dCH4     &
            +  linoz_dPmL_dh2o   * dH2O     &
            +  linoz_dPmL_dT     * dTemp    &
            +  linoz_dPmL_dO3col * dCOL
! 
!pointer to ln2o
       linoz_PmL_clim     => fields(ln2o_PmL_clim_ndx)     %data(:,:,lchnk ) !1/sec not vrm/sec already scaled by n2o vmr
       linoz_dPmL_dO3     => fields(ln2o_dPmL_dO3_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dn2o    => fields(ln2o_dPmL_dN2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dnoy    => fields(ln2o_dPmL_dNOY_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dch4    => fields(ln2o_dPmL_dCH4_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dh2o    => fields(ln2o_dPmL_dH2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dT      => fields(ln2o_dPmL_dT_ndx)      %data(:,:,lchnk )
       linoz_dPmL_dO3col  => fields(ln2o_dPmL_dO3col_ndx)  %data(:,:,lchnk )
! Taylor expanding loss freq (1/sec) of N2O
       Lfreq_n2o =  linoz_PmL_clim           &
            +  linoz_dPmL_dO3    * dO3    &
            +  linoz_dPmL_dn2o   * dN2O   &
            +  linoz_dPmL_dnoy   * dNOY   &
            +  linoz_dPmL_dch4   * dCH4   &
            +  linoz_dPmL_dh2o   * dH2O   &
            +  linoz_dPmL_dT     * dTemp  &
            +  linoz_dPmL_dO3col * dCOL
!
!pointer to production of NOy
       linoz_PmL_clim     => fields(pnoy_PmL_clim_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dO3     => fields(pnoy_dPmL_dO3_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dn2o    => fields(pnoy_dPmL_dN2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dnoy    => fields(pnoy_dPmL_dNOY_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dch4    => fields(pnoy_dPmL_dCH4_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dh2o    => fields(pnoy_dPmL_dH2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dT      => fields(pnoy_dPmL_dT_ndx)      %data(:,:,lchnk )
       linoz_dPmL_dO3col  => fields(pnoy_dPmL_dO3col_ndx)  %data(:,:,lchnk )
!Taylar expanding production of noy divided by fn2o (so unit 1/sec)      
       Pfreq_noy =  linoz_PmL_clim            & 
            +  linoz_dPmL_dO3   * dO3    &
            +  linoz_dPmL_dn2o  * dN2O   &
            +  linoz_dPmL_dnoy  * dNOY   &
            +  linoz_dPmL_dch4  * dCH4   &
            +  linoz_dPmL_dh2o  * dH2O   &
            +  linoz_dPmL_dT    * dTemp  &
            +  linoz_dPmL_dO3col* dCOL
!
!pointer to loss of NOY
       linoz_PmL_clim     => fields(lnoy_PmL_clim_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dO3     => fields(lnoy_dPmL_dO3_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dn2o    => fields(lnoy_dPmL_dN2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dnoy    => fields(lnoy_dPmL_dNOY_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dch4    => fields(lnoy_dPmL_dCH4_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dh2o    => fields(lnoy_dPmL_dH2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dT      => fields(lnoy_dPmL_dT_ndx)      %data(:,:,lchnk )
       linoz_dPmL_dO3col  => fields(lnoy_dPmL_dO3col_ndx)  %data(:,:,lchnk )
!Taylar expanding loss of noy divided by fn2o (so unit 1/sec)      
       Lfreq_noy =  linoz_PmL_clim    & 
            +  linoz_dPmL_dO3   * dO3  &
            +  linoz_dPmL_dn2o  * dN2O &
            +  linoz_dPmL_dnoy  * dNOY &
            +  linoz_dPmL_dch4  * dCH4 &
            +  linoz_dPmL_dh2o  * dH2O &
            +  linoz_dPmL_dT    * dTemp&
            +  linoz_dPmL_dO3col* dCOL

!pointer to loss of CH4
       linoz_PmL_clim     => fields(nch4_PmL_clim_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dO3     => fields(nch4_dPmL_dO3_ndx)     %data(:,:,lchnk )
       linoz_dPmL_dn2o    => fields(nch4_dPmL_dN2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dnoy    => fields(nch4_dPmL_dNOY_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dch4    => fields(nch4_dPmL_dCH4_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dh2o    => fields(nch4_dPmL_dH2O_ndx)    %data(:,:,lchnk )
       linoz_dPmL_dT      => fields(nch4_dPmL_dT_ndx)      %data(:,:,lchnk )
       linoz_dPmL_dO3col  => fields(nch4_dPmL_dO3col_ndx)  %data(:,:,lchnk )
! Taylor expanding loss frequency of CH4
       Lfreq_ch4 =  linoz_PmL_clim            & 
            +  linoz_dPmL_dO3   * dO3    &
            +  linoz_dPmL_dn2o  * dN2O   &
            +  linoz_dPmL_dnoy  * dNOY   &
            +  linoz_dPmL_dch4  * dCH4   &
            +  linoz_dPmL_dh2o  * dH2O   &
            +  linoz_dPmL_dT    * dTemp  &
            +  linoz_dPmL_dO3col* dCOL
 
       linoz_cariolle_psc => fields(cariolle_pscs_ndx)%data(:,:,lchnk )
    
       ! initialize output arrays
       !
       do3_linoz         = 0._r8
       do3_linoz_psc     = 0._r8
       o3col_du_diag     = 0._r8
       o3clim_linoz_diag = 0._r8
       ss_o3             = 0._r8
    !
    ! convert lats from radians to degrees
    !
       lats = rlats * radians_to_degrees

       LOOP_COL: do i=1,ncol
          if (.not. present(tropFlag)) then
             lt0= ltrop(i)
             lt=  ltrop(i)
             n_dl=0
             do ll= lt0,1,-1
                if ( pmid(i,ll) > pressure_threshold ) THEN            
                   lt= ll-1
                endif
             enddo
             n_dl= lt -lt0          
             kmax = ltrop(i)+ n_dl
          else
             kmax = pver
          endif
          LOOP_LEV: do k=1, kmax
             if (present(tropFlag)) then
                 if (tropFlag(i,k) .or. pmid(i,k)>pressure_threshold) cycle
             endif
          ! current mixing ratio
             o3_old =   o3_vmr(i,k)
             n2o_old=   n2o_vmr(i,k)
             noy_old=   noy_vmr(i,k)
             ch4_old=   ch4_vmr(i,k)
             h2o_old=   h2o_vmr(i,k)
    
          ! climatological ozone
             o3_clim = linoz_o3_clim(i,k)
          ! diagnostic for output
          
! convert o3col from mol/cm2
             o3col_du           = o3col(i,k) * convert_to_du !for diagnostics only
! o3          
             ss_x= o3_clim - PL_O3(i,k)/linoz_dPml_dO3X(i,k)             
             delo3 = (ss_x -o3_old)* (1.0_r8 - exp(linoz_dPmL_dO3X(i,k)*delta_t))
             o3_new = o3_old + delo3
!diagnostics
             do3_linoz(i,k) = (o3_new - o3_old)/delta_t
             o3clim_linoz_diag(i,k) = o3_clim
             ss_o3(i,k) = ss_x    
 
!             o3col_du_diag(i,k) = o3col_du      
!n2o  
             dn2op   = max(P_n2o(i,k)* delta_t, 0.0_r8)  !(dp/dt *dt in vmr)
             Lfreq   = max(Lfreq_n2o(i,k), 0.0_r8) ! n2o loss frequency          
             dn2ol   = n2o_old*(exp(-Lfreq * delta_t) - 1.0_r8)
             n2o_new = n2o_old + dn2op + dn2ol
! noy
             dnoyp =   max(Pfreq_noy(i,k)* n2o_old* delta_t, 0.0_r8)  !scaled by fn2o in Pfreq so multiply it back to get mr/sec
             Lfreq =   max(Lfreq_noy(i,k), 0.0_r8) ! n2o loss frequency          
             dnoyl =   noy_old*(exp(-Lfreq*delta_t) - 1.0_r8)
             noy_new = noy_old + dnoyp + dnoyl
! ch4   
             Lfreq   =    max(Lfreq_ch4(i,k), 0.0_r8)                    
             delch4  =    ch4_old*(exp(-Lfreq*delta_t) - 1.0_r8)
             ch4_new =    ch4_old + delch4
             h2o_new =    pw - 2._r8 * ch4_new
!alternative h2o, -2*delch4 gain as the loss of delch4             
!             h2o_new  =    h2o_old - 2._r8* del1ch4          !
          ! PSC activation (follows Cariolle et al 1990.)
          ! use only if abs(latitude) > 40.
          !
             delo3_psc = 0._r8
             if ( abs(lats(i)) > 40._r8  ) then   
                if ( (chlorine_loading-chlorine_loading_bgnd) > 0._r8 ) then
                   if ( temp(i,k) <= psc_T ) then
!
                   ! define maximum SZA for PSC loss (= tangent height at sunset)
                   !
                      max_sza = (90._r8 + sqrt( max( 16._r8*log10(100000._r8/pmid(i,k)),0._r8)))
#ifdef DEBUG
                      write(iulog,'(A, 2f8.1)')'sza/max_saz', sza(i),max_sza
#endif
                      if ( (sza(i)*radians_to_degrees) <= max_sza ) then
                  
!                      write(iulog,*)'linoz_cariolle_psc(i,k)=',linoz_cariolle_psc(i,k)

                         psc_loss = exp(-linoz_cariolle_psc(i,k) &
                              * (chlorine_loading/chlorine_loading_1987)**2 &
                              * delta_t )
!update from last o3_new due to gas chemisry change
                    ! use o3_new to prevent negative value                     
                         delo3_psc= o3_new*(psc_loss -1._r8)
                         o3_new =o3_new + delo3_psc !o3_new:update from gas chem loss

                      !
                      ! output diagnostic
                      !
                         do3_linoz_psc(i,k) = delo3_psc/delta_t
                      !
                      end if
                   end if
                end if
             end if
          !
          ! update vmr
          ! as a defense the assignments are only performed when the species are active. Otherwise the  
          ! index would be an invalid value (-1)
           if (o3lnz_ndx > 0) xvmr(i,k,  o3lnz_ndx) =   o3_new
           if (n2olnz_ndx > 0) xvmr(i,k, n2olnz_ndx)   = n2o_new
           if (noylnz_ndx > 0) xvmr(i,k, noylnz_ndx)   = noy_new
           if (ch4lnz_ndx > 0) xvmr(i,k, ch4lnz_ndx)   = ch4_new
           if (h2olnz_ndx > 0) xvmr(i,k, h2olnz_ndx)   = h2o_new

          !update real o3, ch4, n2o      
           if(o3_ndx  > 0) xvmr(i,k, o3_ndx ) =  delo3   + delo3_psc +  xvmr(i,k, o3_ndx )
           if(ch4_ndx > 0) xvmr(i,k, ch4_ndx) =  delch4  +  xvmr(i,k, ch4_ndx)
           if(n2o_ndx > 0) xvmr(i,k, n2o_ndx) =  (dn2op + dn2ol)  +  xvmr(i,k, n2o_ndx)
           if(no_ndx >0)  xvmr(i,k, no_ndx)   =  0.05 *(dnoyp + dnoyl) + xvmr(i,k, no_ndx)
           if(no2_ndx>0)  xvmr(i,k, no2_ndx)  =  0.05 *(dnoyp + dnoyl) + xvmr(i,k, no2_ndx)
           if(hno3_ndx>0) xvmr(i,k,hno3_ndx)  =  0.90 *(dnoyp + dnoyl) + xvmr(i,k, hno3_ndx)
           
        end do LOOP_LEV

       !use tropospheric specific humidity
        if (.not. present(tropFlag)) then
          LOOP_TROPLEV: do k= ltrop(i)+ n_dl, pver
            xvmr(i,k, h2olnz_ndx)   = h2ovmr(i,k)
          end do LOOP_TROPLEV
        else
          do k= 1, pver
            if (tropFlag(i,k)) then
              xvmr(i,k, h2olnz_ndx)   = h2ovmr(i,k)
            endif
          end do
        endif
        
     end do LOOP_COL
!    save the passed-in o3col for all layers rather than just the stratosphere
     o3col_du_diag(:ncol,:pver) = o3col(:ncol,:pver) * convert_to_du
     do3_linoz_du(:ncol,:) = pdeldry(:ncol,:)*do3_linoz(:ncol,:)*avogadro*rgrav/mw_air*convert_to_du*1.e3_r8
     do3_linoz_psc_du(:ncol,:) = pdeldry(:ncol,:)*do3_linoz_psc(:ncol,:)*avogadro*rgrav/mw_air*convert_to_du*1.e3_r8

     twod_do3_linoz = 0._r8
     twod_do3_linoz_psc = 0._r8
     do k = 1, pver
        twod_do3_linoz(:)     = twod_do3_linoz(:)     + do3_linoz_du(:,k) 
        twod_do3_linoz_psc(:) = twod_do3_linoz_psc(:) + do3_linoz_psc_du(:,k) 
     end do 
    ! output
    !

 
    !call outfld( 'LINOZ_DO3LNZ'      , do3_linoz              , ncol, lchnk )
    !call outfld( 'LINOZ_DO3LNZ_PSC'  , do3_linoz_psc          , ncol, lchnk )
    !call outfld( 'LINOZ_2DDO3LNZ'    , twod_do3_linoz             , ncol, lchnk )
    !call outfld( 'LINOZ_2DDO3LNZ_PSC', twod_do3_linoz_psc         , ncol, lchnk )
    call outfld( 'LINOZ_SSO3'   , ss_o3                  , ncol, lchnk )
    call outfld( 'LINOZ_O3COL'  , o3col_du_diag          , ncol, lchnk )
    call outfld( 'LINOZ_O3CLIM' , o3clim_linoz_diag      , ncol, lchnk )
    call outfld( 'LINOZ_SZA'    ,(sza*radians_to_degrees), ncol, lchnk )

    call outfld( 'LINOZ_DO3'         , do3_linoz              , ncol, lchnk )
    call outfld( 'LINOZ_DO3_PSC'     , do3_linoz_psc          , ncol, lchnk )
    call outfld( 'LINOZ_2DDO3'       , twod_do3_linoz              , ncol, lchnk )
    call outfld( 'LINOZ_2DDO3_PSC'   , twod_do3_linoz_psc          , ncol, lchnk )
    
    return
  end subroutine linv3_strat_chem_solve

 !--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine linv2_strat_chem_solve( ncol, lchnk, xvmr, o3col, temp, sza, pmid, delta_t, rlats, ltrop, pdeldry, tropFlag )
!  modified from Linoz v2 written by   
    use chlorine_loading_data, only: chlorine_loading
    use chem_mods,             only: gas_pcnst
    use mo_constants,          only: pi, rgrav, rearth, avogadro
    !
    ! this subroutine updates the ozone mixing ratio in the stratosphere
    ! using linearized chemistry 
    !
    use ppgrid,        only : pcols, pver
    use physconst,     only : pi, &
                              grav => gravit, &
                              mw_air => mwdry
    use cam_history,   only : outfld

    use linoz_data,  only : fields, o3_clim_ndx, n2o_clim_ndx, noy_clim_ndx, ch4_clim_ndx, h2o_clim_ndx, t_clim_ndx,o3col_clim_ndx,&
         no3_PmL_clim_ndx,  no3_dPmL_dO3_ndx,   no3_dPmL_dN2O_ndx, no3_dPmL_dNOy_ndx,   &
         no3_dPmL_dCH4_ndx, no3_dPmL_dH2O_ndx,  no3_dPmL_dT_ndx,   no3_dPmL_dO3col_ndx, &
         cariolle_pscs_ndx
    !
    ! dummy arguments
    !
    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index
    real(r8), intent(inout), dimension(ncol ,pver ,gas_pcnst) :: xvmr     ! volume mixing ratio for all
    real(r8), intent(in)   , dimension(ncol ,pver) :: o3col               ! ozone column above box (mol/cm^2)
    real(r8), intent(in)   , dimension(pcols,pver) :: temp                ! temperature (K)
    real(r8), intent(in)   , dimension(ncol )      :: sza                 ! local solar zenith angle
    real(r8), intent(in)   , dimension(pcols,pver) :: pmid                ! midpoint pressure (Pa)
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)
    real(r8), intent(in)                           :: rlats(ncol)         ! column latitudes (radians)
    integer,  intent(in)   , dimension(pcols)      :: ltrop               ! chunk index
    real(r8), intent(in)   , dimension(ncol ,pver) :: pdeldry             !  dry pressure delta about midpoints (Pa)
    logical, optional, intent(in)                  :: tropFlag(pcols,pver)! 3D tropospheric level flag
    !
    ! local
    !
    integer :: i,k,n, ll, lt0, lt, n_dl !,index_lat,index_month
    integer :: kmax
    real(r8) :: o3col_du,delta_temp,delta_o3col,o3_old,o3_new,delta_o3,delta_o3_psc
    real(r8) :: max_sza, psc_loss
    real(r8) :: o3_clim
    real(r8), dimension(ncol) :: lats
    real(r8), dimension(ncol,pver) :: do3_linoz,do3_linoz_psc,ss_o3,o3col_du_diag,o3clim_linoz_diag

    real(r8), dimension(:,:), pointer :: linoz_o3_clim
    real(r8), dimension(:,:), pointer :: linoz_t_clim
    real(r8), dimension(:,:), pointer :: linoz_o3col_clim
    real(r8), dimension(:,:), pointer :: linoz_PmL_clim
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dT
    real(r8), dimension(:,:), pointer :: linoz_dPmL_dO3col
    real(r8), dimension(:,:), pointer :: linoz_cariolle_psc
    ! real O3 variables
    real(r8), dimension(ncol,pver) :: do3_linoz_du, do3_linoz_psc_du
    real(r8), dimension(ncol) :: twod_do3_linoz
    real(r8), dimension(ncol) :: twod_do3_linoz_psc
    !
    ! parameters
    !
    real(r8), parameter :: convert_to_du = 1._r8/(2.687e16_r8)      ! convert ozone column from mol/cm^2 to DU
    real(r8), parameter :: degrees_to_radians = pi/180._r8          ! conversion factors
    real(r8), parameter :: radians_to_degrees = 180._r8/pi
    real(r8), parameter :: chlorine_loading_1987    = 2.5977_r8     ! EESC value (ppbv)
    real(r8), parameter :: chlorine_loading_bgnd    = 0.0000_r8     ! EESC value (ppbv) for background conditions
    real(r8), parameter :: pressure_threshold       = 237.1375e+2_r8    ! last vaild linoz pressure layer

    !
    ! skip if no ozone field available
    !
    if ( .not. do_lin_strat_chem ) return
    if ( .not. linoz_v2) return

    ! 
    ! associate the field pointers
    !  old version 
 !   linoz_o3_clim      => fields(o3_clim_ndx)      %data(:,:,lchnk )
 !   linoz_t_clim       => fields(t_clim_ndx)       %data(:,:,lchnk )
 !   linoz_o3col_clim   => fields(o3col_clim_ndx)   %data(:,:,lchnk )
 !   linoz_PmL_clim     => fields(PmL_clim_ndx)     %data(:,:,lchnk )
 !   linoz_dPmL_dO3     => fields(dPmL_dO3_ndx)     %data(:,:,lchnk )
 !   linoz_dPmL_dT      => fields(dPmL_dT_ndx)      %data(:,:,lchnk )
 !   linoz_dPmL_dO3col  => fields(dPmL_dO3col_ndx)  %data(:,:,lchnk )


    !Linoz climatological data 
    linoz_o3_clim      => fields(o3_clim_ndx)         %data(:,:,lchnk )
    linoz_t_clim       => fields(t_clim_ndx)          %data(:,:,lchnk )
    linoz_o3col_clim   => fields(o3col_clim_ndx)      %data(:,:,lchnk )
    linoz_PmL_clim     => fields(no3_PmL_clim_ndx)    %data(:,:,lchnk )    !unit vmr/sec
    linoz_dPmL_dO3     => fields(no3_dPmL_dO3_ndx)    %data(:,:,lchnk )
    linoz_dPmL_dT      => fields(no3_dPmL_dT_ndx)     %data(:,:,lchnk )
    linoz_dPmL_dO3col  => fields(no3_dPmL_dO3col_ndx) %data(:,:,lchnk )
    linoz_cariolle_psc => fields(cariolle_pscs_ndx)   %data(:,:,lchnk )
    !
    ! initialize output arrays
    !
    do3_linoz         = 0._r8
    do3_linoz_psc     = 0._r8
    o3col_du_diag     = 0._r8
    o3clim_linoz_diag = 0._r8
    ss_o3             = 0._r8
    !
    ! convert lats from radians to degrees
    !
    lats = rlats * radians_to_degrees

    LOOP_COL: do i=1,ncol
       if (.not. present(tropFlag)) then
          lt0= ltrop(i)
          lt = ltrop(i)
          n_dl=0
          do ll= lt0,1,-1
             if ( pmid(i,ll) > pressure_threshold ) THEN            
                lt= ll-1
             endif
          enddo
          n_dl = lt - lt0
          kmax = ltrop(i)+ n_dl
       else
          kmax = pver
       endif
       LOOP_LEV: do k=1, kmax
          if (present(tropFlag)) then
            ! skip Linoz if in the troposphere or pressure exceeds threshold
            if (tropFlag(i,k) .or. pmid(i,k)>pressure_threshold) cycle
          endif
          !
          ! climatological ozone
          !
          o3_clim = linoz_o3_clim(i,k)
          !
          ! skip if not in the stratosphere
          !
          !if ( pmid(i,k) > pressure_threshold ) THEN   ! PJC diagnostic
          !   WRITE(iulog,*)'LINOZ WARNING: Exceeded PRESSURE threshold (i,k,p_threshold,pmid,o3)=',&
          !     i,k,nint(pressure_threshold/100._r8),'mb',nint(pmid(i,k)/100._r8),'mb',nint(o3_vmr(i,k)*1e9_r8),'ppb'   !PJC
!         !    cycle LOOP_LEV
          !endif
          !
          ! diagnostic for output
          !
          o3clim_linoz_diag(i,k) = o3_clim
          !
          ! old ozone mixing ratio
          !
          !o3_old = xvmr(i,k, o3lnz_ndx)
          o3_old = xvmr(i,k, o3_ndx)
          !
          ! convert o3col from mol/cm2
          !
          o3col_du = o3col(i,k) * convert_to_du
!          o3col_du_diag(i,k) = o3col_du
          !
          ! compute differences from climatology
          !
          delta_temp  = temp(i,k) - linoz_t_clim    (i,k)
          delta_o3col = o3col_du  - linoz_o3col_clim(i,k)


          !
          ! steady state ozone
          !
          ss_o3(i,k) = o3_clim - (               linoz_PmL_clim   (i,k)   &
                                 + delta_o3col * linoz_dPmL_dO3col(i,k)   &
                                 + delta_temp  * linoz_dPmL_dT    (i,k)   &
                                             ) / linoz_dPmL_dO3   (i,k)


          !
          ! ozone change
          !
          delta_o3 = (ss_o3(i,k)-o3_old) * (1._r8 - exp(linoz_dPmL_dO3(i,k)*delta_t))
          !
          ! define new ozone mixing ratio
          !
          o3_new = o3_old + delta_o3
          !
          ! output diagnostic
          !
          do3_linoz(i,k) = delta_o3/delta_t
          !
          ! PSC activation (follows Cariolle et al 1990.)
          !
          ! use only if abs(latitude) > 40.
          !
          delta_o3_psc = 0._r8
          if ( abs(lats(i)) > 40._r8 ) then   
             if ( (chlorine_loading-chlorine_loading_bgnd) > 0._r8 ) then
                if ( temp(i,k) <= psc_T ) then
                   !
                   ! define maximum SZA for PSC loss (= tangent height at sunset)
                   !
                   max_sza = (90._r8 + sqrt( max( 16._r8*log10(100000._r8/pmid(i,k)),0._r8)))
#ifdef DEBUG
                   write(iulog,'(A, 2f8.1)')'sza/max_saz', sza(i),max_sza
#endif
                   if ( (sza(i)*radians_to_degrees) <= max_sza ) then

                      psc_loss = exp(-linoz_cariolle_psc(i,k) &
                           * (chlorine_loading/chlorine_loading_1987)**2 &
                           * delta_t )
                      ! use o3_new to prevent negative value
                      delta_o3_psc= o3_new*(psc_loss -1._r8)
                      o3_new= delta_o3_psc + o3_new
                      !
                      ! output diagnostic
                      !
                      do3_linoz_psc(i,k) = delta_o3_psc/delta_t
                      !
                   end if
                end if
             end if
          end if
          !
          ! update ozone vmr
          !
          !xvmr(i,k, o3lnz_ndx) = o3_new
          if(o3_ndx >0) xvmr(i,k, o3_ndx) = delta_o3 + delta_o3_psc + xvmr(i,k, o3_ndx)
          
       end do LOOP_LEV
    end do LOOP_COL
    o3col_du_diag(:ncol,:pver) = o3col(:ncol,:pver) * convert_to_du
    do3_linoz_du(:ncol,:) = pdeldry(:ncol,:)*do3_linoz(:ncol,:)*avogadro*rgrav/mw_air*convert_to_du*1.e3_r8
    do3_linoz_psc_du(:ncol,:) = pdeldry(:ncol,:)*do3_linoz_psc(:ncol,:)*avogadro*rgrav/mw_air*convert_to_du*1.e3_r8

    twod_do3_linoz = 0._r8
    twod_do3_linoz_psc = 0._r8
    do k = 1, pver
        twod_do3_linoz(:)     = twod_do3_linoz(:)     + do3_linoz(:,k) 
        twod_do3_linoz_psc(:) = twod_do3_linoz_psc(:) + do3_linoz_psc(:,k) 
    end do 
    ! output
    !
     !call outfld( 'LINOZ_DO3LNZ'      , do3_linoz              , ncol, lchnk )
     !call outfld( 'LINOZ_DO3LNZ_PSC'  , do3_linoz_psc          , ncol, lchnk )
     !call outfld( 'LINOZ_2DDO3LNZ'    , twod_do3_linoz             , ncol, lchnk )
     !call outfld( 'LINOZ_2DDO3LNZ_PSC', twod_do3_linoz_psc         , ncol, lchnk )
     call outfld( 'LINOZ_SSO3'   , ss_o3                  , ncol, lchnk )
     call outfld( 'LINOZ_O3COL'  , o3col_du_diag          , ncol, lchnk )
     call outfld( 'LINOZ_O3CLIM' , o3clim_linoz_diag      , ncol, lchnk )
     call outfld( 'LINOZ_SZA'    ,(sza*radians_to_degrees), ncol, lchnk )

     call outfld( 'LINOZ_DO3'         , do3_linoz              , ncol, lchnk )
     call outfld( 'LINOZ_DO3_PSC'     , do3_linoz_psc          , ncol, lchnk )
     call outfld( 'LINOZ_2DDO3'       , twod_do3_linoz              , ncol, lchnk )
     call outfld( 'LINOZ_2DDO3_PSC'   , twod_do3_linoz_psc          , ncol, lchnk )
    
    return
  end subroutine linv2_strat_chem_solve

 subroutine lin_strat_sfcsink( ncol, lchnk, x_vmr, x_sfc, delta_t, pdel)

    use ppgrid,        only : pcols, pver

    use physconst,     only : mw_air => mwdry
   
    use mo_constants, only : pi, rgrav, rearth

    use phys_grid,    only : get_area_all_p

    use cam_history,   only : outfld

    use chem_mods,             only: gas_pcnst, adv_mass

    implicit none 

    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index
    real(r8), intent(inout), dimension(ncol ,pver, gas_pcnst) :: x_vmr    ! volume mixing ratio of ndx
    real(r8), intent(in)                           :: x_sfc(4,ncol)       ! sfc concentration 
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)    
    real(r8), intent(in)                           :: pdel(ncol,pver)     ! pressure delta about midnpoints (Pa)  
    real(r8), parameter :: mw_noylnz                = 14.0_r8
    real(r8), parameter :: KgtoTg                   = 1.0e-9_r8
    real(r8), parameter :: peryear                  = 86400._r8* 365.0_r8 ! to multiply to convert per second to per year
    
    real(r8) :: area(ncol), mass(ncol,pver)
    real(r8) :: x_old, x_new, efactor, dx, mw_x
    real(r8), dimension(ncol)  :: dx_mass, x_sfcsink
    real(r8), dimension(4,ncol)  :: sfc_const
    real(r8), dimension(4)  :: mw
    integer , dimension(4)  :: nx
    integer i, j, k, n, ms

    if ( .not. do_lin_strat_chem )    return
    if(.not. linoz_v2 .and. .not. linoz_v3)return
 !
 !   initializing array
 !  
    ms =0   
    if (linoz_v2) then
     if (o3lnz_ndx >0) then
       ms=1
       nx(1) =  o3lnz_ndx
       mw(1) =  adv_mass(o3lnz_ndx)
       sfc_const(1,:)= o3_sfc
       o3_lbl=4

     elseif (o3_ndx > 0) then
       ms=1
       nx(1) =  o3_ndx
       mw(1) =  adv_mass(o3_ndx)
       sfc_const(1,:)= o3_sfc
       o3_lbl=4
       !write(iulog,*)'warning: linoz_v2 is on and surface o3 loss is active'
     else
       write(iulog,*)'warning: linoz_v2 is on but no surface o3 loss'      
       return
     endif
    endif
!       
    if (linoz_v3)then
     if (o3lnz_ndx >0) then
       ms=4
       nx(1) = o3lnz_ndx
       nx(2) = n2olnz_ndx
       nx(3) = noylnz_ndx
       nx(4) = ch4lnz_ndx 
       mw(1) =  adv_mass(o3lnz_ndx)       
       mw(2) =  adv_mass(n2olnz_ndx)  
       mw(3) =  adv_mass(noylnz_ndx)
       mw(4) =  adv_mass(ch4lnz_ndx)   
       sfc_const(1:4,:ncol) = x_sfc(1:4,:ncol)
       o3_lbl =9
     else
       ms=3
       nx(1) = n2olnz_ndx
       nx(2) = noylnz_ndx
       nx(3) = ch4lnz_ndx
       mw(1) =  adv_mass(n2olnz_ndx)  
       mw(2) =  adv_mass(noylnz_ndx)
       mw(3) =  adv_mass(ch4lnz_ndx)
       sfc_const(1,:ncol) = x_sfc(2,:ncol)
       sfc_const(2,:ncol) = x_sfc(3,:ncol)
       sfc_const(3,:ncol) = x_sfc(4,:ncol)
       o3_lbl =9
     endif    
    endif    
! 
    do k = 1,pver
       mass(:ncol,k) = pdel(:ncol,k) * rgrav  ! air mass in kg/m2
    enddo

!o3_tau is 2-days applied to all species    
    efactor  = 1.d0 - exp(-delta_t/o3_tau)

    LOOP_N:  do n=1, ms

       x_sfcsink(:ncol) = 0._r8
          
       LOOP_COL: do i=1,ncol

          dx_mass(i) =0._r8

          LOOP_SFC: do k= pver, pver-o3_lbl+1, -1
!
             j= nx(n)
             x_old = x_vmr(i,k,j)  !vmr
             
             dx =  (sfc_const(n,i) - x_old)* efactor !vmr
             x_new  = x_old + dx
! loss in kg/m2 summed over boundary layers within one time step   
             dx_mass(i) = dx_mass(i) + dx* mass(i,k) * mw(n)/mw_air    
             x_vmr(i,k,j) = x_new

          end do  LOOP_SFC !loop-k

       End do  Loop_COL    !loop-col
       
       x_sfcsink(:ncol) = dx_mass(:ncol)/delta_t * KgtoTg * peryear ! saved in Tg/yr/m2 unit
    
       if(o3lnz_ndx > 0 .and. j.eq.  o3lnz_ndx)  call outfld('LINOZ_O3SFCSINK',    x_sfcsink, ncol, lchnk)
       if(n2olnz_ndx >0 .and. j.eq. n2olnz_ndx)  call outfld('LINOZ_N2OSFCSRC',    x_sfcsink, ncol, lchnk)
       if(noylnz_ndx >0 .and. j.eq. noylnz_ndx)  call outfld('LINOZ_NOYSFCSINK',   x_sfcsink, ncol, lchnk)       
       if(ch4lnz_ndx >0 .and. j.eq. ch4lnz_ndx)  call outfld('LINOZ_CH4SFCSRC',    x_sfcsink, ncol, lchnk)

    End do Loop_N

    return     

  end subroutine lin_strat_sfcsink


   subroutine fstrat_efold_inti(fstrat_efold_list)
    use cam_abortutils,       only : endrun
    use chem_mods,            only : gas_pcnst
    use mo_chem_utls,         only : get_spc_ndx   
    implicit none

    character(len=*), intent(in) :: fstrat_efold_list(:)    
    integer i,j
    has_fstrat_efold(:) = .false.

    do i = 1, gas_pcnst

       if ( len_trim(fstrat_efold_list(i))==0 ) exit

       j = get_spc_ndx(fstrat_efold_list(i))

       if ( j > 0 ) then
          has_fstrat_efold(j) = .true.
       else
          write(iulog,*) 'fstrat_efold_inti: '//trim(fstrat_efold_list(i))//' is not included in species set'
          call endrun('fstrat_inti: invalid stratosphere decaying species')
       endif

    enddo

    return

   end subroutine fstrat_efold_inti


   subroutine fstrat_efold_decay(ncol, xvmr, delta_t, ltrop, tropFlag)

 !   use namelist_utils,       only : fstrat_efold_list
    use mo_chem_utls,         only : get_spc_ndx
    use cam_abortutils,       only : endrun
    use chem_mods,            only : gas_pcnst
    use ppgrid,               only : pcols, pver

   
    implicit none

    
    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    real(r8), intent(inout), dimension(ncol ,pver ,gas_pcnst) :: xvmr     ! volume mixing ratio for all
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)
    integer,  intent(in)   , dimension(pcols)      :: ltrop               ! chunk index
    logical, optional, intent(in)                  :: tropFlag(pcols,pver)! 3D tropospheric level flag

    integer :: ispc, i, k, kmax
    real(r8):: efactor, dx

    !real(r8), parameter :: tau_30d = 1._r8/(30._r8*86400._r8)      ! inverse of (30 days*86400 sec/day)
    ! HHLEE 20220112 changing decay time from 30 days to 90 days
    real(r8), parameter :: tau_30d = 90._r8*86400._r8      ! inverse of (30 days*86400 sec/day)
    integer , parameter :: ms = 18      ! 
    
     efactor  = 1._r8 - exp(-delta_t/tau_30d)
   
     LOOP_SPC: do ispc = 1, gas_pcnst
          
        if (has_fstrat_efold(ispc)) then
 !         write(iulog,*)'ispc= ', ispc
           LOOP_COL: do i=1, ncol
              if (.not. present(tropFlag)) then
                kmax = ltrop(i)
              else
                kmax = pver
              endif

              LOOP_STRAT: do k= 1, kmax
                 if (present(tropFlag)) then
                    ! skip if in the troposphere
                    if (tropFlag(i,k)) cycle
                 endif
                   
                 dx   = xvmr(i,k,ispc) *efactor
                 xvmr(i,k,ispc) = xvmr(i,k,ispc) - dx

              end do LOOP_STRAT ! loop from top to tropopause

           End do  LOOP_COL    !loop-col
        endif
     END DO  LOOP_SPC

  return
    
  end subroutine fstrat_efold_decay



end module lin_strat_chem
