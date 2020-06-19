
!--------------------------------------------------------------------
! linearized ozone chemistry LINOZ
! from Hsu and Prather, JGR, 2008
!
! written by Jean-Francois Lamarque (September 2008)
! modified by
!     24 Oct 2008 -- Francis Vitt
!      9 Dec 2008 -- Philip Cameron-Smith, LLNL, -- added ltrop
!      4 Jul 2019 -- Qi Tang (LLNL), Juno Hsu (UCI), -- added sfcsink
!--------------------------------------------------------------------
module lin_strat_chem

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : begchunk, endchunk
  use physics_types,only : physics_state
  use cam_logfile,  only : iulog
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  !
  implicit none
  !
  private  ! all unless made public

  save
  !
  ! define public components of module
  !
  public :: lin_strat_chem_inti, lin_strat_chem_solve, lin_strat_sfcsink
  public :: do_lin_strat_chem
  public :: linoz_readnl   ! read linoz_nl namelist

  integer :: index_o3
!  integer :: index_o3l

  logical :: do_lin_strat_chem

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
   linoz_lbl = 4
   linoz_sfc = 30.0e-9_r8
   linoz_tau = 172800.0_r8
   linoz_psc_T = 193.0_r8

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
      write(iulog,*) subname // ', linoz_lbl:', o3_lbl
      write(iulog,*) subname // ', linoz_sfc:', o3_sfc
      write(iulog,*) subname // ', linoz_tau:', o3_tau
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
    use linoz_data,   only : linoz_data_init, has_linoz_data
    use ppgrid,       only : pver
    use mo_chem_utls, only : get_spc_ndx
    use cam_history,  only : addfld, horiz_only, add_default
    use physics_buffer, only : physics_buffer_desc

    implicit none


    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    if (.not.has_linoz_data) return

    !
    ! find index of ozone
    !
    index_o3  =  get_spc_ndx('O3')
!    index_o3l =  get_spc_ndx('O3L') ! introduce O3L that couples to parameterized surface sink

    do_lin_strat_chem = has_linoz_data
!    if ( index_o3 <= 0 .and. index_o3l <=0 ) then !added index_o3l
     if ( index_o3 <= 0 ) then !added index_o3l
       write(iulog,*) ' No ozone in the chemical mechanism, skipping lin_strat_chem'
       do_lin_strat_chem = .false.
       return
    end if

    ! check for synoz

    if( get_spc_ndx( 'SYNOZ' ) > 0 .and. has_linoz_data) then
       call endrun('lin_strat_chem_inti: cannot have both synoz and linoz')
    endif

    ! initialize the linoz data

    call linoz_data_init()

    ! define additional output

    call addfld( 'LINOZ_DO3'    , (/ 'lev' /), 'A', '/s'     , 'ozone vmr tendency by linearized ozone chemistry'   )
    call addfld( 'LINOZ_DO3_PSC', (/ 'lev' /), 'A', '/s'     , 'ozone vmr loss by PSCs using Carille et al. (1990)' )
    call addfld( 'LINOZ_SSO3'   , (/ 'lev' /), 'A', 'kg'     , 'steady state ozone in LINOZ'                        )
    call addfld( 'LINOZ_O3COL'  , (/ 'lev' /), 'A', 'DU'     , 'ozone column above'                                 )
    call addfld( 'LINOZ_O3CLIM' , (/ 'lev' /), 'A', 'mol/mol', 'climatology of ozone in LINOZ'                      )
    call addfld( 'LINOZ_SZA'    ,    horiz_only, 'A', 'degrees', 'solar zenith angle in LINOZ'                      )
    call addfld( 'LINOZ_SFCSINK', horiz_only, 'A', 'Tg/yr/m2'   , 'surface o3 sink in LINOZ with an e-fold to a fixed concentration  ' )

    call add_default( 'LINOZ_DO3'    , 1, ' ' )
    call add_default( 'LINOZ_DO3_PSC', 1, ' ' )
    call add_default( 'LINOZ_SSO3'   , 1, ' ' )
    call add_default( 'LINOZ_O3COL'  , 1, ' ' )
    call add_default( 'LINOZ_O3CLIM' , 1, ' ' )
    call add_default( 'LINOZ_SZA'    , 1, ' ' )
    call add_default( 'LINOZ_SFCSINK', 1, ' ' )
    return
  end subroutine lin_strat_chem_inti


!--------------------------------------------------------------------
!--------------------------------------------------------------------
  subroutine lin_strat_chem_solve( ncol, lchnk, o3_vmr, o3col, temp, sza, pmid, delta_t, rlats, ltrop )
 
    use chlorine_loading_data, only: chlorine_loading

    !
    ! this subroutine updates the ozone mixing ratio in the stratosphere
    ! using linearized chemistry 
    !
    use ppgrid,        only : pcols, pver
    use physconst,     only : pi, &
                              grav => gravit, &
                              mw_air => mwdry
    use cam_history,   only : outfld
    use linoz_data,    only : fields, o3_clim_ndx,t_clim_ndx,o3col_clim_ndx,PmL_clim_ndx,dPmL_dO3_ndx,&
                                      dPmL_dT_ndx,dPmL_dO3col_ndx,cariolle_pscs_ndx
    !
    ! dummy arguments
    !
    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index
    real(r8), intent(inout), dimension(ncol ,pver) :: o3_vmr              ! ozone volume mixing ratio
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

    ! 
    ! associate the field pointers
    !
    linoz_o3_clim      => fields(o3_clim_ndx)      %data(:,:,lchnk )
    linoz_t_clim       => fields(t_clim_ndx)       %data(:,:,lchnk )
    linoz_o3col_clim   => fields(o3col_clim_ndx)   %data(:,:,lchnk )
    linoz_PmL_clim     => fields(PmL_clim_ndx)     %data(:,:,lchnk )
    linoz_dPmL_dO3     => fields(dPmL_dO3_ndx)     %data(:,:,lchnk )
    linoz_dPmL_dT      => fields(dPmL_dT_ndx)      %data(:,:,lchnk )
    linoz_dPmL_dO3col  => fields(dPmL_dO3col_ndx)  %data(:,:,lchnk )
    linoz_cariolle_psc => fields(cariolle_pscs_ndx)%data(:,:,lchnk )

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
          if ( pmid(i,k) > pressure_threshold ) THEN   ! PJC diagnostic
             WRITE(iulog,*)'LINOZ WARNING: Exceeded PRESSURE threshold (i,k,p_threshold,pmid,o3)=',&
               i,k,nint(pressure_threshold/100._r8),'mb',nint(pmid(i,k)/100._r8),'mb',nint(o3_vmr(i,k)*1e9_r8),'ppb'   !PJC
!             cycle LOOP_LEV
          endif
          !
          ! diagnostic for output
          !
          o3clim_linoz_diag(i,k) = o3_clim
          !
          ! old ozone mixing ratio
          !
          o3_old = o3_vmr(i,k)
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
                      do3_linoz_psc(i,k) = (o3_new-o3_old)/delta_t
                      !
                   end if
                end if
             end if
          end if
          !
          ! update ozone vmr
          !
          o3_vmr(i,k) = o3_new

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
  end subroutine lin_strat_chem_solve

  subroutine lin_strat_sfcsink( ncol, lchnk, o3l_vmr, delta_t, pdel)

    use ppgrid,        only : pcols, pver

    use physconst,     only : mw_air => mwdry, &
                              mw_o3  => mwo3   !in grams

    use mo_constants, only : pi, rgrav, rearth

    use phys_grid,    only : get_area_all_p

    use cam_history,   only : outfld

    implicit none 

    integer,  intent(in)                           :: ncol                ! number of columns in chunk
    integer,  intent(in)                           :: lchnk               ! chunk index    
    real(r8), intent(inout), dimension(ncol ,pver) :: o3l_vmr             ! ozone volume mixing ratio
    real(r8), intent(in)                           :: delta_t             ! timestep size (secs)    
    real(r8), intent(in)                           :: pdel(ncol,pver)     ! pressure delta about midpoints (Pa)  
    
! three parameters are applied to Linoz O3 for surface sink, O3l is not coupled to real O3
    real(r8), parameter :: KgtoTg                   = 1.0e-9_r8
    real(r8), parameter :: peryear                  = 86400._r8* 365.0_r8 ! to multiply to convert per second to per year
    
    real(r8) :: area(ncol), mass(ncol,pver)
    real(r8) :: o3l_old, o3l_new, efactor, do3
    real(r8), dimension(ncol)  :: do3mass, o3l_sfcsink
    integer i, k

    if ( .not. do_lin_strat_chem ) return
 !
 !   initializing array
 !

    o3l_sfcsink =0._r8

    do k = 1,pver
       mass(:ncol,k) = pdel(:ncol,k) * rgrav  ! air mass in kg/m2
    enddo
    
    efactor  = (1.d0 - exp(-delta_t/o3_tau))
    LOOP_COL: do i=1,ncol

       do3mass(i) =0._r8

       LOOP_SFC: do k= pver, pver-o3_lbl+1, -1
  
          o3l_old = o3l_vmr(i,k)  !vmr
   
          do3 =  (o3_sfc - o3l_old)* efactor !vmr

          o3l_new  = o3l_old + do3     

          do3mass(i) = do3mass(i) + do3* mass(i,k) * mw_o3/mw_air  ! loss in kg/m2 summed over boundary layers within one time step   

          o3l_vmr(i,k) = o3l_new

       end do  LOOP_SFC
    End do  Loop_COL
       
    o3l_sfcsink(:ncol) = do3mass(:ncol)/delta_t * KgtoTg * peryear ! saved in Tg/yr/m2 unit
    
    call outfld('LINOZ_SFCSINK', o3l_sfcsink, ncol, lchnk)    

    return     

  end subroutine lin_strat_sfcsink



end module lin_strat_chem
