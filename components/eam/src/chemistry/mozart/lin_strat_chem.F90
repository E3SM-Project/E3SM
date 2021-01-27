!--------------------------------------------------------------------
module lin_strat_chem
!     24 Oct 2008 -- Francis Vitt
!      9 Dec 2008 -- Philip Cameron-Smith, LLNL, -- added ltrop
!      4 Jul 2019 -- Qi Tang (LLNL), Juno Hsu (UCI), -- added sfcsink
!     20 Jan 2021 -- Juno Hsu(UCI) added Linoz v3 
!        --new linoz v3 subroutine linv3_strat_chem_solve (linoz_v3 if O3LNZ, N2OLNZ, NOYLNZ and CH4LNZ are defined)
!        --modified lin_strat_chem_solve is now linv2_strat_chem_solve (linoz_v2 if only O3LNZ is defined; use only O3LNZ part of Linoz v3 netcdf file)
!        --modified sfcsink working for v2 or v3 depending on species
  use shr_kind_mod , only : r8 => shr_kind_r8
  use ppgrid       , only : begchunk, endchunk
  use physics_types, only : physics_state
  use cam_logfile  , only : iulog
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  
  !
  implicit none
  !
  private  ! all unless made public

  save
  !
  ! define public components of module
  !
  public :: lin_strat_chem_inti, linv2_strat_chem_solve, linv3_strat_chem_solve, lin_strat_sfcsink
  public :: do_lin_strat_chem, linoz_v2, linoz_v3
  public :: linoz_readnl   ! read linoz_nl namelist

  integer :: o3lnz_ndx, n2olnz_ndx, noylnz_ndx, ch4lnz_ndx 
  integer :: o3_ndx, n2o_ndx, ch4_ndx, no_ndx, no2_ndx, hno3_ndx 
  
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
    use mo_chem_utls, only : get_spc_ndx
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

    if (o3lnz_ndx <=0 ) then
       write(iulog,*) 'need to have tracer O3LNZ at least '
       do_lin_strat_chem = .false.
       return
    end if

    linoz_v3= (o3lnz_ndx > 0 .and. n2olnz_ndx >0  .and. noylnz_ndx >0  .and. ch4lnz_ndx > 0)
    linoz_v2= (o3lnz_ndx > 0 .and. n2olnz_ndx <0  .and. noylnz_ndx <0  .and. ch4lnz_ndx < 0)
!    write(iulog,*)'linoz_v3=',linoz_v3,'linoz_v2=',linoz_v2
! real o3, ch4, n2o tracers
    o3_ndx   =   get_spc_ndx('O3')
    ch4_ndx  =   get_spc_ndx('CH4')
    n2o_ndx  =   get_spc_ndx('N2O')
    no_ndx   =   get_spc_ndx('NO')
    no2_ndx  =   get_spc_ndx('NO2')
    hno3_ndx =   get_spc_ndx('HNO3')
!     
!   initialize the linoz data

    call linoz_data_init()

    ! define additional output

    call addfld( 'LINOZ_DO3'    , (/ 'lev' /), 'A', '1/s'     , 'ozone vmr tendency by linearized ozone chemistry'   )
    call addfld( 'LINOZ_DO3_PSC', (/ 'lev' /), 'A', '1/s'     , 'ozone vmr loss by PSCs using Carille et al. (1990)' )
    call addfld( 'LINOZ_SSO3'   , (/ 'lev' /), 'A', 'kg'     , 'steady state ozone in LINOZ'                        )
    call addfld( 'LINOZ_O3COL'  , (/ 'lev' /), 'A', 'DU'     , 'ozone column above'                                 )
    call addfld( 'LINOZ_O3CLIM' , (/ 'lev' /), 'A', 'mol/mol', 'climatology of ozone in LINOZ'                      )
    call addfld( 'LINOZ_SZA'    ,    horiz_only, 'A', 'degrees', 'solar zenith angle in LINOZ'                      )
    call addfld( 'LINOZ_O3SFCSINK',  horiz_only, 'A', 'Tg/yr/m2'   ,   'surface o3lnz sink in LINOZ with an e-fold to a fixed concentration' )
    if(linoz_v3)then
    call addfld( 'LINOZ_NOYSFCSINK', horiz_only, 'A', 'Tg/yr/m2'   , 'surface noylnz sink in LINOZ v3 with an e-fold to a fixed concentration' )
    call addfld( 'LINOZ_N2OSFCSRC',  horiz_only, 'A', 'Tg/yr/m2'   , 'surface n2o source in LINOZ v3 with an e-fold to a fixed concentration' )
    call addfld( 'LINOZ_CH4SFCSRC',  horiz_only, 'A', 'Tg/yr/m2'   , 'surface n2o source in LINOZ v3 with an e-fold to a fixed concentration' )
    endif

    call add_default( 'LINOZ_DO3'    , 1, ' ' )
    call add_default( 'LINOZ_DO3_PSC', 1, ' ' )
    call add_default( 'LINOZ_SSO3'   , 1, ' ' )
    call add_default( 'LINOZ_O3COL'  , 1, ' ' )
    call add_default( 'LINOZ_O3CLIM' , 1, ' ' )
    call add_default( 'LINOZ_SZA'    , 1, ' ' )
    call add_default( 'LINOZ_O3SFCSINK', 1, ' ' )
    if(linoz_v3)then
       call add_default( 'LINOZ_NOYSFCSINK', 1, ' ' )
       call add_default( 'LINOZ_N2OSFCSRC', 1, ' ' )
       call add_default( 'LINOZ_CH4SFCSRC', 1, ' ' )
    endif
    return
  end subroutine lin_strat_chem_inti


!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine linv3_strat_chem_solve(ncol, lchnk, xvmr, h2ovmr, xsfc, o3col, temp, sza, pmid, delta_t, rlats, ltrop)

!--------------------------------------------------------------------
! linearized ozone chemistry linoz-v3 (o3-n2o-noy-ch4 prognostic equations, h2o diagnosed from ch4) 
! from Hsu and Prather, grl, 2009   (https://doi.org/10.1029/2009GL042243)
! 
!written by Juno Hsu (junoh@uci.edu), 09/2020 
 
    use chlorine_loading_data, only: chlorine_loading
    use chem_mods,             only: gas_pcnst
    use mo_chem_utls,          only: get_spc_ndx
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
                              cariolle_pscs_ndx
    !
    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index
    real(r8), intent(inout), dimension(ncol ,pver,gas_pcnst) :: xvmr      ! volume mixing ratio for all
    real(r8), intent(inout), dimension(ncol ,pver) :: h2ovmr              ! vh2o olume mixing ratio for all
    real(r8), intent(inout)                        :: xsfc(4)             ! surface ch4 and n2o from linoz table
    real(r8), intent(in)   , dimension(ncol ,pver) :: o3col               ! ozone column above box (mol/cm^2)
    real(r8), intent(in)   , dimension(pcols,pver) :: temp                ! temperature (K)
    real(r8), intent(in)   , dimension(ncol )      :: sza                 ! local solar zenith angle
    real(r8), intent(in)   , dimension(pcols,pver) :: pmid                ! midpoint pressure (Pa)
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)
    real(r8), intent(in)                           :: rlats(ncol)         ! column latitudes (radians)
    integer,  intent(in)   , dimension(pcols)      :: ltrop               ! chunk index    
    !
    integer  :: i,k,n !,index_lat,index_month
    real(r8) :: o3col_du,delta_temp,delta_o3col    
    real(r8) ::o3_old,   n2o_old,  noy_old,  ch4_old,  h2o_old
    real(r8) ::o3_new,   n2o_new,  noy_new,  ch4_new,  h2o_new
    real(r8) ::o3_clim, ss_x
    real(r8) :: dn2op, dn2ol, dnoyp, dnoyl, delo3, delch4, lfreq
    real(r8) :: max_sza, psc_loss, pw
    real(r8), dimension(ncol) :: lats
    real(r8), dimension(ncol,pver) :: do3_linoz, do3_linoz_psc, ss_o3, o3col_du_diag, o3clim_linoz_diag
    ! local
    real(r8), dimension(ncol,pver) :: o3_vmr, o3lnz_vmr, n2o_vmr, noy_vmr, ch4_vmr, h2o_vmr
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
    !
    ! parameters
    !
    real(r8), parameter :: convert_to_du = 1._r8/(2.687e16_r8)      ! convert ozone column from mol/cm^2 to DU
    real(r8), parameter :: degrees_to_radians = pi/180._r8          ! conversion factors
    real(r8), parameter :: radians_to_degrees = 180._r8/pi
    real(r8), parameter :: chlorine_loading_1987    = 2.5977_r8     ! EESC value (ppbv)
    real(r8), parameter :: chlorine_loading_bgnd    = 0.0000_r8     ! EESC value (ppbv) for background conditions
    real(r8), parameter :: small    = 1.e-18_r8     ! prevent dividing by zero

    !
    if ( .not. do_lin_strat_chem ) return
    if ( .not. linoz_v3) return

!    write(iulog,*)'inside lin_strat_solve for ndx o3, o3lnz, n2o, noylnz, ch4', o3_ndx, o3lnz_ndx, n2o_ndx, noylnz_ndx, ch4_ndx

        o3_vmr =  xvmr(:,:, o3lnz_ndx)
       n2o_vmr =  xvmr(:,:,n2olnz_ndx)
       noy_vmr =  xvmr(:,:,noylnz_ndx)
       ch4_vmr =  xvmr(:,:,ch4lnz_ndx)
       h2o_vmr =  h2ovmr(:,:)


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

       dO3 (:,:)   =  o3_vmr(:,:)  - linoz_o3_clim(:,:)
       dN2O(:,:)  =  n2o_vmr(:,:)  - linoz_n2o_clim(:,:)
       dNOY(:,:)  =  noy_vmr(:,:)  - linoz_noy_clim(:,:) 
       dCH4(:,:)  =  ch4_vmr(:,:)  - linoz_ch4_clim(:,:)
       dH2O(:,:)  =  h2o_vmr(:,:)  - linoz_h2o_clim(:,:) 
       dTemp(:,:) =  temp(:,:)     - linoz_t_clim(:,:)
       dCOL(:,:)  =  o3col(:,:)*convert_to_du - linoz_o3col_clim(:,:)  

! potential water 2*ch4 +h2o !it might be better to get maxch4 3-4 years back in time but this will do for now 
! upated yearly ch4 value is stored at ch4_clim in linoz file at the surface layer
!output surface constant concentration to control surface sink/source
       xsfc(1)=       o3_sfc ! ozone surface constant
       xsfc(2)=       maxval(linoz_n2o_clim(1:ncol,pver)) !n2o
       xsfc(3)=       o3_sfc*3.e-3_r8             !noylnz
       xsfc(4)=       maxval(linoz_ch4_clim(1:ncol,pver)) !ch4
       pw= 2.0_r8 * xsfc(4) + 3.65e-6_r8 

!       write(iulog,*)'xsfc, o3', xsfc(1)
!       write(iulog,*)'xsfc, n2o',xsfc(2)
!       do I= 1, ncol
!          write(iulog,*),'i=',i
!          do k=1, pver
!           write(iulog,*)'k=',k, 'n2oval=', linoz_n2o_clim(i,k)
!          enddo
!       enddo
!       write(iulog,*)'xsfc, noy',xsfc(3)
!       write(iulog,*)'xsfc, ch4',xsfc(4)
 
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
          LOOP_LEV: do k=1,ltrop(i)
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
             o3col_du_diag(i,k) = o3col_du      
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
          !
          ! PSC activation (follows Cariolle et al 1990.)
          ! use only if abs(latitude) > 40.
          !
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
!
                         o3_new = o3_old * psc_loss
                      !
                      ! output diagnostic
                      !
                         do3_linoz_psc(i,k) = (o3_new - o3_old)/delta_t
                      !
                      end if
                   end if
                end if
             end if
          !
          ! update vmr
           
           xvmr(i,k,  o3lnz_ndx) =   o3_new
           xvmr(i,k, n2olnz_ndx)   = n2o_new
           xvmr(i,k, noylnz_ndx)   = noy_new
           xvmr(i,k, ch4lnz_ndx)   = ch4_new
!update real o3, ch4, n2o      
           if(o3_ndx  > 0) xvmr(i,k, o3_ndx ) =  delo3   +  xvmr(i,k, o3_ndx )
           if(ch4_ndx > 0) xvmr(i,k, ch4_ndx) =  delch4  +  xvmr(i,k, ch4_ndx)
           if(n2o_ndx > 0) xvmr(i,k, n2o_ndx) =  (dn2op + dn2ol)  +  xvmr(i,k, n2o_ndx)
!           h2ovmr(i,k) =  pw - 2._r8* xvmr(i, k, ch4lnz_ndx)
!update h2ovmr (this is tied to H2O_gas and used in strat chem if has_strato_chem)     
             h2ovmr(i,k) = -2._r8 * delch4 + h2ovmr(i,k) 
           if(no_ndx >0)  xvmr(i,k, no_ndx)   =  0.05 *(dnoyp + dnoyl) + xvmr(i,k, no_ndx)
           if(no2_ndx>0)  xvmr(i,k, no2_ndx)  =  0.05 *(dnoyp + dnoyl) + xvmr(i,k, no2_ndx)
           if(hno3_ndx>0) xvmr(i,k,hno3_ndx)  =  0.90 *(dnoyp + dnoyl) + xvmr(i,k, hno3_ndx)

        end do LOOP_LEV
     end do LOOP_COL
    !
    ! output
    !
     call outfld( 'LINOZ_DO3'    , do3_linoz              , ncol, lchnk )
     call outfld( 'LINOZ_DO3_PSC', do3_linoz_psc          , ncol, lchnk )
     call outfld( 'LINOZ_SSO3'   , ss_o3                  , ncol, lchnk )
     call outfld( 'LINOZ_O3COL'  , o3col_du_diag          , ncol, lchnk )
     call outfld( 'LINOZ_O3CLIM' , o3clim_linoz_diag      , ncol, lchnk )
     call outfld( 'LINOZ_SZA'    ,(sza*radians_to_degrees), ncol, lchnk )

    return
  end subroutine linv3_strat_chem_solve

 !--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine linv2_strat_chem_solve( ncol, lchnk, xvmr, o3col, temp, sza, pmid, delta_t, rlats, ltrop )
!  modified from Linoz v2 written by   
    use chlorine_loading_data, only: chlorine_loading
    use chem_mods,             only: gas_pcnst
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
    !
    ! local
    !
    integer :: i,k,n !,index_lat,index_month
    real(r8) :: o3col_du,delta_temp,delta_o3col,o3_old,o3_new,delta_o3
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
    !
    ! parameters
    !
    real(r8), parameter :: convert_to_du = 1._r8/(2.687e16_r8)      ! convert ozone column from mol/cm^2 to DU
    real(r8), parameter :: degrees_to_radians = pi/180._r8          ! conversion factors
    real(r8), parameter :: radians_to_degrees = 180._r8/pi
    real(r8), parameter :: chlorine_loading_1987    = 2.5977_r8     ! EESC value (ppbv)
    real(r8), parameter :: chlorine_loading_bgnd    = 0.0000_r8     ! EESC value (ppbv) for background conditions
    real(r8), parameter :: pressure_threshold       = 210.e+2_r8    ! {PJC} for diagnostics only

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
       LOOP_LEV: do k=1,ltrop(i)
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
          o3_old = xvmr(i,k, o3lnz_ndx)
          !
          ! convert o3col from mol/cm2
          !
          o3col_du = o3col(i,k) * convert_to_du
          o3col_du_diag(i,k) = o3col_du
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

                      o3_new = o3_old * psc_loss
                      !
                      ! output diagnostic
                      !
                      do3_linoz_psc(i,k) = (o3_new - o3_old)/delta_t
                      !
                   end if
                end if
             end if
          end if
          !
          ! update ozone vmr
          !
          xvmr(i,k, o3lnz_ndx) = o3_new
          if(o3_ndx >0) xvmr(i,k, o3_ndx) = (o3_new - o3_old) + xvmr(i,k, o3_ndx)
 
         
       end do LOOP_LEV
    end do LOOP_COL
    !
    ! output
    !
    call outfld( 'LINOZ_DO3'    , do3_linoz              , ncol, lchnk )
    call outfld( 'LINOZ_DO3_PSC', do3_linoz_psc          , ncol, lchnk )
    call outfld( 'LINOZ_SSO3'   , ss_o3                  , ncol, lchnk )
    call outfld( 'LINOZ_O3COL'  , o3col_du_diag          , ncol, lchnk )
    call outfld( 'LINOZ_O3CLIM' , o3clim_linoz_diag      , ncol, lchnk )
    call outfld( 'LINOZ_SZA'    ,(sza*radians_to_degrees), ncol, lchnk )
    
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
    real(r8), intent(in)                           :: x_sfc(4)            ! sfc concentration 
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)    
    real(r8), intent(in)                           :: pdel(ncol,pver)     ! pressure delta about midnpoints (Pa)  
    real(r8), parameter :: mw_noylnz                = 14.0_r8
    real(r8), parameter :: KgtoTg                   = 1.0e-9_r8
    real(r8), parameter :: peryear                  = 86400._r8* 365.0_r8 ! to multiply to convert per second to per year
    
    real(r8) :: area(ncol), mass(ncol,pver)
    real(r8) :: x_old, x_new, efactor, dx, mw_x
    real(r8), dimension(ncol)  :: dx_mass, x_sfcsink
    real(r8), dimension(4)  :: mw, sfc_const
    integer , dimension(4)  :: nx
    integer i, j, k, n, ms

    if ( .not. do_lin_strat_chem )    return
    if(.not. linoz_v2 .and. .not. linoz_v3)return
 !
 !   initializing array
 !  
    ms =0   
    if(linoz_v2)then
       ms=1
       nx(1) =  o3lnz_ndx
       mw(1) =  adv_mass(o3lnz_ndx)
       sfc_const(1)= o3_sfc
    endif
!       
    if(linoz_v3)then
       ms=4
       nx(1) = o3lnz_ndx
       nx(2) = n2olnz_ndx
       nx(3) = noylnz_ndx
       nx(4) = ch4lnz_ndx ! turn off if trop ch4 is on
       mw(1) =  adv_mass(o3lnz_ndx)       
       mw(2) =  adv_mass(n2olnz_ndx)  
       mw(3) =  adv_mass(noylnz_ndx)
       mw(4) =  adv_mass(ch4lnz_ndx)   
       sfc_const(1:4) = x_sfc(1:4)
 !      write(iulog,*) 'sfc species ndx=', nx(1), nx(2), nx(3), nx(4)
 !      write(iulog,*) 'sfc species mw=', mw(1), mw(2), mw(3), mw(4)
 !      write(iulog,*) 'sfc species x_sfc=', x_sfc(1),x_sfc(2), x_sfc(3), x_sfc(4)
 !      write(iulog,*) 'sfc cocentration=', sfc_const(1), sfc_const(2), sfc_const(3), sfc_const(4)
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
             dx =  (sfc_const(n) - x_old)* efactor !vmr
             x_new  = x_old + dx
! loss in kg/m2 summed over boundary layers within one time step   
             dx_mass(i) = dx_mass(i) + dx* mass(i,k) * mw(n)/mw_air    
             x_vmr(i,k,j) = x_new

          end do  LOOP_SFC !loop-k

       End do  Loop_COL    !loop-col
       
       x_sfcsink(:ncol) = dx_mass(:ncol)/delta_t * KgtoTg * peryear ! saved in Tg/yr/m2 unit
    
       if(j.eq.  o3lnz_ndx)  call outfld('LINOZ_O3SFCSINK',    x_sfcsink, ncol, lchnk)
       if(j.eq. n2olnz_ndx)  call outfld('LINOZ_N2OSFCSRC',    x_sfcsink, ncol, lchnk)
       if(j.eq. noylnz_ndx)  call outfld('LINOZ_NOYSFCSINK',   x_sfcsink, ncol, lchnk)       
       if(j.eq. ch4lnz_ndx)  call outfld('LINOZ_CH4SFCSRC',    x_sfcsink, ncol, lchnk)

    End do Loop_N

    return     

  end subroutine lin_strat_sfcsink



end module lin_strat_chem
