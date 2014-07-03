#ifdef AIX
#define USE_ESSL
#endif
#define USE_BDE

      module mo_jshort

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use physconst,     only : pi
      use mo_constants,  only : d2r
      use abortutils,    only : endrun
      use cam_logfile,   only : iulog
      use spmd_utils,    only : masterproc
      use ppgrid,        only : pver
      use phys_control,  only : waccmx_is

      implicit none

      interface jshort
         module procedure jshort_photo
         module procedure jshort_hrates
      end interface

      private
      public :: jshort_init
      public :: jshort_timestep_init
      public :: jshort
      public :: sphers
      public :: slant_col
      public :: nj

      save

!------------------------------------------------------------------------------
!     ... define altitude and wavelength parameters
!------------------------------------------------------------------------------
      integer, parameter  :: num_ms93tuv = 4       ! wavelength bins for MS, 93
      integer, parameter  :: nw_ms93     = 4       ! wavelength bins for MS, 93
      integer, parameter  :: nsrc_tot    = 19      ! total bins for SRC
      integer, parameter  :: nsrb_tot    = 14      ! total bins <200nm for SRB
      integer, parameter  :: nsrbtuv     = 17      ! total SRB bins in TUV
      real(r8), parameter :: hc          = 6.62608e-34_r8 * 2.9979e8_r8 / 1.e-9_r8
      real(r8), parameter :: wc_o2_a     = 175.05_r8   ! (nm)
      real(r8), parameter :: wc_o2_b     = 242.37_r8   ! (nm)
      real(r8), parameter :: wc_o3_a     = 310.32_r8   ! (nm)
      real(r8), parameter :: wc_o3_b     = 1179.87_r8  ! (nm)
      real(r8), parameter :: we_ms(nw_ms93+1) = (/ 181.6_r8, 183.1_r8, 184.6_r8, 190.2_r8, 192.5_r8 /)

      integer  :: nw                                 ! Number of wavelength bins <200nm
      integer  :: nj                                 ! Number of photorates
      real(r8) :: wtno50(6,2)
      real(r8) :: wtno90(6,2)
      real(r8) :: wtno100(6,2)
      real(r8) :: csno50(6,2)
      real(r8) :: csno90(6,2)
      real(r8) :: csno100(6,2)
      real(r8) :: ac(20,nsrbtuv)
      real(r8) :: bc(20,nsrbtuv)     ! chebyshev polynomial coeffs
      real(r8) :: wave_num(nsrbtuv)
      real(r8), allocatable :: wc(:)
      real(r8), allocatable :: we(:)
      real(r8), allocatable :: wlintv(:)
      real(r8), allocatable :: bde_o2_a(:)
      real(r8), allocatable :: bde_o2_b(:)
      real(r8), allocatable :: bde_o3_a(:)
      real(r8), allocatable :: bde_o3_b(:)
      real(r8), allocatable :: etfphot(:)
      real(r8), allocatable :: etfphot_ms93(:)
      real(r8), allocatable :: xs_o2src(:)
      real(r8), allocatable :: xs_o3a(:)
      real(r8), allocatable :: xs_o3b(:)
      real(r8), allocatable :: xs_wl(:,:)

      contains

      subroutine jshort_init( xs_coef_file, xs_short_file, sht_indexer )
!------------------------------------------------------------------------------
!    ... initialize photorates for 120nm <= lambda <= 200nm
!------------------------------------------------------------------------------            

      use mo_util,        only : rebin
      use solar_data,  only : data_nbins=>nbins, data_we => we, data_etf => sol_etf

      implicit none

!------------------------------------------------------------------------------
!    ... dummy arguments
!------------------------------------------------------------------------------            
      character(len=*), intent(in) :: xs_coef_file, xs_short_file
      integer, intent(inout)       :: sht_indexer(:)

!------------------------------------------------------------------------------
!     ... set the wavelength grid, exoatmospheric flux,
!         and cross sections (for <200nm) - contained in
!         a NetCDF file
!------------------------------------------------------------------------------
      call get_crs( xs_short_file, sht_indexer )
      if(masterproc) then
         write(iulog,*) ' '
         write(iulog,*) '============================================'
         write(iulog,*) 'jshort_init: finished get_crs'
         write(iulog,*) 'jshort_init: nj, nw = ',nj,nw
         write(iulog,*) 'jshort_init: wc'
         write(iulog,*) wc(:)
         write(iulog,*) '============================================'
         write(iulog,*) ' '
      end if
      we(:nw)  = wc(:nw) - .5_r8*wlintv(:nw)
      we(nw+1) = wc(nw) + .5_r8*wlintv(nw)
      if(masterproc) then
         write(iulog,*) ' '
         write(iulog,*) '-------------------------------------------'
         write(iulog,*) 'jshort_init: diagnostics before rebin'
         write(iulog,*) 'jshort_init: data_nbins, nw = ',data_nbins,nw
         write(iulog,*) 'jshort_init: data_we range'
         write(iulog,'(1p,5g15.7)') minval(data_we(:)),maxval(data_we(:))
         write(iulog,*) 'jshort_init: we range'
         write(iulog,'(1p,5g15.7)') minval(we(:)),maxval(we(:))
      end if
      call rebin( data_nbins, nw, data_we, we, data_etf, etfphot )
      if(masterproc) then
         write(iulog,*) 'jshort_init: etfphot'
         write(iulog,'(1p,5g15.7)') etfphot(:)
         write(iulog,*) '-------------------------------------------'
         write(iulog,*) ' '
         write(iulog,*) 'jshort_init: diagnostics for ms93'
         call rebin( data_nbins, nw_ms93, data_we, we_ms, data_etf, etfphot_ms93 )
         write(iulog,'(1p,5g15.7)') etfphot_ms93(:)
         write(iulog,*) '-------------------------------------------'
         write(iulog,*) ' '
      end if
!------------------------------------------------------------------------------
!     ... loads Chebyshev polynomial Coeff
!------------------------------------------------------------------------------
      call xs_init(xs_coef_file)

!------------------------------------------------------------------------------
!     ... initialize no photolysis
!------------------------------------------------------------------------------
      call jno_init

      end subroutine jshort_init

      subroutine get_crs( xs_short_file, sht_indexer )
!=============================================================================
!   PURPOSE:
!   Reads a NetCDF file that contains:
!     Cross_sections*quantum yield data for species <200nm.
!     Exoatmospheric flux on the wavelength grid of the crs
!=============================================================================
!   PARAMETERS:
!     Input:
!      xs_short_file.... NetCDF file that contains the crs*QY for wavenum species
!
!     Output:
!      xs_species.. Cross Sections * QY data for each species
!      etfphot..... Exoatmospheric flux in photons cm-2 s-1 nm-1
!      etfphot_ms93.Minshwanner and Siskind JNO SRB etf (on MS93 grid)
!      wc.......... wavelength center (nm)
!      numwl ...... Number of wavelengths < 200nm in NetCDF input file
!      wlintv ..... Wavelength inteval for grid, nm
!=============================================================================
!   EDIT HISTORY:
!   Created by Doug Kinnison, 1/14/2002
!   Modified by S. Walters, 4/2/2003
!=============================================================================

      use chem_mods,      only : phtcnt, pht_alias_lst, rxt_tag_lst
      use ioFileMod,      only : getfil
      use error_messages, only : alloc_err
      use pio,            only : file_desc_t, pio_get_var, pio_closefile, pio_noerr, &
           pio_seterrorhandling, pio_bcast_error, pio_internal_error, pio_inq_varid, &
           pio_inq_dimlen, pio_nowrite, pio_inq_dimid
      use cam_pio_utils,  only : cam_pio_openfile
      implicit none

!------------------------------------------------------------------------------
!       ... dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout)       :: sht_indexer(:)
      character(len=*), intent(in) :: xs_short_file

!------------------------------------------------------------------------------
!       ... local variables
!------------------------------------------------------------------------------
      integer :: wn
      type(file_desc_t) :: ncid
      integer :: ierr
      integer :: i, m, ndx
      integer :: varid, dimid
      integer :: wrk_ndx(phtcnt)
      real(r8), allocatable :: xs_species(:)
      character(len=256) :: locfn

!------------------------------------------------------------------------------
!       ... open NetCDF File
!------------------------------------------------------------------------------
      call getfil(xs_short_file, locfn, 0)
      call cam_pio_openfile(ncid, trim(locfn), PIO_NOWRITE)

!------------------------------------------------------------------------------
!       ... get dimensions
!------------------------------------------------------------------------------
      ierr = pio_inq_dimid( ncid, 'numwl', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, nw )

!------------------------------------------------------------------------------
!       ... check for cross section in dataset
!------------------------------------------------------------------------------
      call pio_seterrorhandling(ncid, pio_bcast_error)
      do m = 1,phtcnt
         if( pht_alias_lst(m,1) == ' ' ) then
            ierr = pio_inq_varid( ncid, rxt_tag_lst(m), varid )
            if( ierr == PIO_noerr ) then 
               sht_indexer(m) = varid
            end if
         else if( pht_alias_lst(m,1) == 'userdefined' ) then
            sht_indexer(m) = -1
         else
            ierr = pio_inq_varid( ncid, pht_alias_lst(m,1), varid )
            if( ierr == PIO_noerr ) then 
               sht_indexer(m) = varid
            else
               write(iulog,*) 'get_crs : ',rxt_tag_lst(m)(:len_trim(rxt_tag_lst(m))),' alias ', &
                    pht_alias_lst(m,1)(:len_trim(pht_alias_lst(m,1))),' not in dataset'            
               call endrun
            end if
         end if
      end do
      call pio_seterrorhandling(ncid, pio_internal_error)

      if (masterproc) then
         write(iulog,*) ' '
         write(iulog,*) '###############################################'
         write(iulog,*) 'get_crs: sht_indexer'
         write(iulog,'(10i6)') sht_indexer(:)
         write(iulog,*) '###############################################'
         write(iulog,*) ' '
      endif

      nj = 0
      do m = 1,phtcnt
         if( sht_indexer(m) > 0 ) then
            if( any( sht_indexer(:m-1) == sht_indexer(m) ) ) then
               cycle
            end if
            nj = nj + 1
         end if
      end do

!------------------------------------------------------------------------------
!       ... allocate arrays
!------------------------------------------------------------------------------
      allocate( wc(nw),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'wc', nw )
      end if
      allocate( we(nw+1),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'we', nw+1 )
      end if
      allocate( wlintv(nw),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'wlintv', nw )
      end if
      allocate( etfphot(nw),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'etfphot', nw )
      end if
      allocate( bde_o2_a(nw),bde_o2_b(nw),bde_o3_a(nw),bde_o3_b(nw),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'bde_o2_a ... bde_o3_b', nw )
      end if
      allocate( etfphot_ms93(nw_ms93),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'etfphot_ms93', nw_ms93 )
      end if
      allocate( xs_o2src(nw),stat=ierr )
      if( ierr /= 0 ) then 
	 call alloc_err( ierr, 'get_crs', 'xs_o2src', nw )
      end if
      allocate( xs_o3a(nw),xs_o3b(nw),stat=ierr )
      if( ierr /= 0 ) then 
         call alloc_err( ierr, 'get_crs', 'xs_o3a,xs_o3b', nw )
      end if
      allocate( xs_species(nw),xs_wl(nw,nj),stat=ierr )
      if( ierr /= 0 ) then 
         call alloc_err( ierr, 'get_crs', 'xs_species,xs_wl', nw*nj )
      end if

!------------------------------------------------------------------------------
!       ... read arrays
!------------------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'wc', varid )
      ierr = pio_get_var( ncid, varid, wc )
      ierr = pio_inq_varid( ncid, 'wlintv', varid )
      ierr = pio_get_var( ncid, varid, wlintv )
      ierr = pio_inq_varid( ncid, 'xs_o2src', varid )
      ierr = pio_get_var( ncid, varid, xs_o2src )
      ndx = 0
      do m = 1,phtcnt
         if( sht_indexer(m) > 0 ) then
            if( any( sht_indexer(:m-1) == sht_indexer(m) ) ) then
               cycle
            end if
            ierr = pio_get_var( ncid, sht_indexer(m), xs_species )
            ndx = ndx + 1
            xs_wl(:,ndx) = xs_species(:)*wlintv(:)
         end if
      end do
      deallocate( xs_species )
      if( ndx /= nj ) then
         write(iulog,*) 'get_crs : ndx count /= cross section count'
         call endrun
      end if
      !------------------------------------------------------------------------------
      !       ... get jo3 cross sections
      !------------------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'jo3_a', varid )
      ierr = pio_get_var( ncid, varid, xs_o3a )
      ierr = pio_inq_varid( ncid, 'jo3_b', varid )
      ierr = pio_get_var( ncid, varid, xs_o3b )
      !------------------------------------------------------------------------------
      !       ... setup final sht_indexer
      !------------------------------------------------------------------------------
      ndx = 0
      wrk_ndx(:) = sht_indexer(:)
      do m = 1,phtcnt
         if( wrk_ndx(m) > 0 ) then
            ndx = ndx + 1
            i = wrk_ndx(m)
            where( wrk_ndx(:) == i )
               sht_indexer(:) = ndx
               wrk_ndx(:)     = -100000
            end where
         end if
      end do

      call pio_closefile( ncid )

#ifdef USE_BDE
      if (masterproc) write(iulog,*) 'Jshort using bdes'
#else
      if (masterproc) write(iulog,*) 'Jshort not using bdes'
#endif
      do wn = 1,nw
#ifdef USE_BDE
         bde_o2_a(wn) = max(0._r8, hc*(wc_o2_a - wc(wn))/(wc_o2_a*wc(wn)) )
         bde_o2_b(wn) = max(0._r8, hc*(wc_o2_b - wc(wn))/(wc_o2_b*wc(wn)) )
         bde_o3_a(wn) = max(0._r8, hc*(wc_o3_a - wc(wn))/(wc_o3_a*wc(wn)) )
         bde_o3_b(wn) = max(0._r8, hc*(wc_o3_b - wc(wn))/(wc_o3_b*wc(wn)) )
#else
         bde_o2_a(wn) = hc/wc(wn)
         bde_o2_b(wn) = hc/wc(wn)
         bde_o3_a(wn) = hc/wc(wn)
         bde_o3_b(wn) = hc/wc(wn)
#endif
      end do

      end subroutine get_crs

      subroutine xs_init(xs_coef_file)
!-------------------------------------------------------------
!    	... Loads XS_COEFFS containing the Chebyshev
!           polynomial coeffs necessary to calculate O2 effective
!           cross-sections
!-------------------------------------------------------------

      use ioFileMod,     only : getfil
      use units,         only : getunit, freeunit

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------            
      character(len=*), intent(in) :: xs_coef_file

!-------------------------------------------------------------
!     ... Local variables
!-------------------------------------------------------------
      integer :: unit	! file unit number
      integer :: istat	! i/o status
      integer :: i, j
      character(len=256) :: locfn

!----------------------------------------------------------------------
!	... Get first strato photo rate file
!----------------------------------------------------------------------
         call getfil(xs_coef_file, locfn, 0)
!----------------------------------------------------------------------
!	... open file
!----------------------------------------------------------------------
         unit = getunit()
         open( unit   = unit, &
               file   = trim(locfn), &
               status = 'old', &
               form   = 'formatted', &
               iostat = istat )
         if( istat /= 0 ) then
!----------------------------------------------------------------------
!	... Open error exit
!----------------------------------------------------------------------
            write(iulog,*) 'xs_init: error ',istat,' opening file ',trim(locfn)
            call endrun
         end if
!----------------------------------------------------------------------
!	... read file
!----------------------------------------------------------------------
         read(unit,901)
         do i = 1,20
            read(unit,903,iostat=istat) ac(i,:)
            if( istat /= 0 ) then
               write(iulog,*) 'xs_init: error ',istat,' reading ac'
               call endrun
            end if
         end do

         read(unit,901)
         do i = 1,20
            read(unit,903,iostat=istat) bc(i,:)
            if( istat /= 0 ) then
               write(iulog,*) 'xs_init: error ',istat,' reading bc'
               call endrun
            end if
         end do
         close( unit )
         call freeunit( unit )

      wave_num(17:1:-1) = (/ (48250._r8 + (500._r8*i),i=1,17) /)

 901  format( / )
 903  format( 17(e23.14,1x))

      end subroutine xs_init

      subroutine jno_init
!-------------------------------------------------------------
!    	... Initialization for no photolysis
!-------------------------------------------------------------

      implicit none

!-------------------------------------------------------------
!     ... Local variables
!-------------------------------------------------------------
      real(r8), dimension(24) :: a, b, c

!------------------------------------------------------------------------------
!   	... 6 sub-intervals for O2 5-0 at 265K,
!	    2 sub-sub-intervals for NO 0-0 at 250K
!------------------------------------------------------------------------------
      a(:) = (/    0._r8,       0._r8,       0._r8,       0._r8, &
                   5.12e-02_r8, 5.68e-03_r8, 1.32e-18_r8, 4.41e-17_r8, &
                   1.36e-01_r8, 1.52e-02_r8, 6.35e-19_r8, 4.45e-17_r8, &
                   1.65e-01_r8, 1.83e-02_r8, 7.09e-19_r8, 4.50e-17_r8, &
                   1.41e-01_r8, 1.57e-02_r8, 2.18e-19_r8, 2.94e-17_r8, &
                   4.50e-02_r8, 5.00e-03_r8, 4.67e-19_r8, 4.35e-17_r8 /)

!------------------------------------------------------------------------------
!   	... sub-intervals for o2 9-0 band,
!	    2 sub-sub-intervals for no 1-0 at 250 k
!------------------------------------------------------------------------------
      b(:) = (/        0._r8,       0._r8,       0._r8,       0._r8, &
                       0._r8,       0._r8,       0._r8,       0._r8, &
                 1.93e-03_r8, 2.14e-04_r8, 3.05e-21_r8, 3.20e-21_r8, &
                 9.73e-02_r8, 1.08e-02_r8, 5.76e-19_r8, 5.71e-17_r8, &
                 9.75e-02_r8, 1.08e-02_r8, 2.29e-18_r8, 9.09e-17_r8, &
                 3.48e-02_r8, 3.86e-03_r8, 2.21e-18_r8, 6.00e-17_r8 /)

!------------------------------------------------------------------------------
! 	... sub-intervals for o2 10-0 band,
!	    2 sub-sub-intervals for no 1-0 at 250 k
!------------------------------------------------------------------------------
      c(:) = (/  4.50e-02_r8, 5.00e-03_r8, 1.80e-18_r8, 1.40e-16_r8, &
                 1.80e-01_r8, 2.00e-02_r8, 1.50e-18_r8, 1.52e-16_r8, &
                 2.25e-01_r8, 2.50e-02_r8, 5.01e-19_r8, 7.00e-17_r8, &
                 2.25e-01_r8, 2.50e-02_r8, 7.20e-20_r8, 2.83e-17_r8, &
                 1.80e-01_r8, 2.00e-02_r8, 6.72e-20_r8, 2.73e-17_r8, &
                 4.50e-02_r8, 5.00e-03_r8, 1.49e-21_r8, 6.57e-18_r8 /)

      wtno50 (1:6,1) = a(1:24:4)
      wtno50 (1:6,2) = a(2:24:4)
      csno50 (1:6,1) = a(3:24:4)
      csno50 (1:6,2) = a(4:24:4)
      wtno90 (1:6,1) = b(1:24:4)
      wtno90 (1:6,2) = b(2:24:4)
      csno90 (1:6,1) = b(3:24:4)
      csno90 (1:6,2) = b(4:24:4)
      wtno100(1:6,1) = c(1:24:4)
      wtno100(1:6,2) = c(2:24:4)
      csno100(1:6,1) = c(3:24:4)
      csno100(1:6,2) = c(4:24:4)

      end subroutine jno_init

      subroutine jshort_timestep_init
!---------------------------------------------------------------
!	... set etfphot if required
!---------------------------------------------------------------

      use time_manager,   only : is_end_curr_day
      use mo_util,        only : rebin
      use solar_data,     only : data_nbins=>nbins, data_we => we, data_etf => sol_etf

      implicit none

      call rebin( data_nbins, nw,      data_we, we,    data_etf, etfphot )
      call rebin( data_nbins, nw_ms93, data_we, we_ms, data_etf, etfphot_ms93 )

      end subroutine jshort_timestep_init

      subroutine jshort_hrates( nlev, zen, o2_vmr, o3_vmr, o2cc, &
                                o3cc, tlev, zkm, mw, qrs, cparg, &
                                lchnk, long, co2cc, scco2, do_diag )
!==============================================================================!
!   Subroutine Jshort                                                          !
!==============================================================================!
!   Purpose:                                                                   !
!     To calculate thermal heating rates for lamba < 200 nm
!==============================================================================!
!   This routine uses JO2 parameterizations based on:                          !
!        Lyman alpha... Chabrillat and Kockarts, GRL, 25, 2659, 1998           !
!        SRC .......... Brasseur and Solomon, 1986 (from TUV)                  !
!        SRB .......... Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996    !
!                        (supplied by Dan Marsh, NCAR ACD                      !
!   and JNO:                                                                   !
!        SRB .......... Minschwanner and Siskind, JGR< 98, 20401, 1993.        !
!==============================================================================!
!   Input:                                                                     !
!       o2cc....... O2 concentration, molecule cm-3			       !
!	o3cc....... O3 concentration, molecule cm-3			       !
!     zen........ zenith angle, units = degrees                                !
!     tlev....... Temperature Profile (K)                                      !
!     zkm ....... Altitude, km                                                 !
!                                                                              !
!   Output:                                                                    !
!     qrs ... short wavelength thermal heating rates
!==============================================================================!
!                                                                              !
!   Approach:                                                                  !
!                                                                              !
!    1) Call sphers (taken from TUV)                                           !
!         -> derives dsdh and nid used in slant column routines                !
!	  -> zenith angle dependent                                            !
!                                                                              !
!    2) Call  slant_col (taken from TUV)                                       !
!		-> derives the slant column for each species                   !
!                                                                              !
!    3) Calls get_crs                                                          !
!		-> read a NetCDF file                                          !
!		-> returns cross sections*quantum yields for all species that  !
!		   have absorption below 200nm.                                !
!                                                                              !
!    4) Derives transmission and photolysis rates for selective species        !
!                                                                              !
!==============================================================================!
!   EDIT HISTORY:                                                              !
!   Created by Doug Kinnison, 3/14/2002                                        !
!==============================================================================!

       use physconst,       only : avogad
       use error_messages, only : alloc_err

       implicit none

       integer, parameter  :: branch = 2           ! two photolysis pathways for JO2
       real(r8), parameter :: km2cm  = 1.e5_r8

!------------------------------------------------------------------------------
!     ... dummy arguments
!------------------------------------------------------------------------------
       integer, intent(in)     :: nlev                 ! model vertical levels
       integer, intent(in)     :: lchnk                ! chunk index
       integer, intent(in)     :: long                 ! chunk index
       real(r8), intent(in)    :: zen                  ! Zenith angle (degrees)
       real(r8), intent(in)    :: o2_vmr(nlev)         ! o2 conc (mol/mol)
       real(r8), intent(in)    :: o3_vmr(nlev)         ! o3 conc (mol/mol)
       real(r8), intent(in)    :: o2cc(nlev)           ! o2 conc (mol/cm^3)
       real(r8), intent(in)    :: co2cc(nlev)          ! co2 conc (mol/cm^3)
       real(r8), intent(in)    :: o3cc(nlev)           ! o3 conc (mol/cm^3)
       real(r8), intent(in)    :: tlev(nlev)           ! Temperature profile
       real(r8), intent(in)    :: zkm(nlev)            ! Altitude, km
       real(r8), intent(in)    :: mw(nlev)             ! atms molecular weight
       real(r8), intent(in)    :: cparg(nlev)          ! column specific heat capacity
       real(r8), intent(inout) :: qrs(:,:)             ! sw heating rates
       real(r8), intent(out)   :: scco2(nlev)          ! co2 column concentration (molec/cm^2)
       logical,  intent(in)    :: do_diag

!------------------------------------------------------------------------------
!     ... local variables
!------------------------------------------------------------------------------
       integer :: k, k1                     ! Altitude indicies
       integer :: wn                        ! Wavelength index
       integer :: astat
       integer :: nid(0:nlev)               ! Number of layers crossed by the direct
                                            ! beam when travelling from the top of the
                                            ! atmosphere to layer i; NID(i), i = 0..NZ-1
       real(r8) :: hfactor
       real(r8) :: dsdh(0:nlev,nlev)        ! Slant path of direct beam through each
                                            ! layer crossed  when travelling from the top of
                                            ! the atmosphere to layer i; DSDH(i,j), i = 0.
                                            ! NZ-1, j = 1..NZ-1   (see sphers.f)
       real(r8), allocatable :: fnorm(:,:)      ! Normalized ETF
       real(r8), allocatable :: trans_o2(:,:)   ! Transmission o2 (total)
       real(r8), allocatable :: trans_o3(:,:)   ! Transmission, ozone
       real(r8), allocatable :: wrk(:)      ! wrk array
       real(r8) :: jo2_lya(nlev)            ! Total photolytic rate constant for Ly alpha
       real(r8) :: jo2_srb(nlev)            ! Total JO2  for SRB
       real(r8) :: jo2_src(nlev)            ! Total JO2 for SRC
       real(r8) :: delz(nlev)               ! layer thickness (cm)
       real(r8) :: o2scol(nlev)             ! O2 Slant Column
       real(r8) :: o3scol(nlev)             ! O3 Slant Column
       real(r8) :: rmla(nlev)               ! Transmission, Lyman Alpha (other species)
       real(r8) :: ro2la(nlev)              ! Transmission, Lyman Alpha (for JO2)
       real(r8) :: tlevmin(nlev)
       real(r8) :: abs_col(nlev)
       real(r8) :: tsrb(nlev,nsrbtuv)       ! Transmission in the SRB
       real(r8) :: xs_o2srb(nlev,nsrbtuv)   ! Cross section * QY for O2 in SRB

      allocate( fnorm(nlev,nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jshort_hrates', 'fnorm', nw*nlev )
      end if
      allocate( trans_o2(nlev,nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jshort_hrates', 'trans_o2', nw*nlev )
      end if
      allocate( trans_o3(nlev,nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jshort_hrates', 'trans_o3', nw*nlev )
      end if
      allocate( wrk(nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jshort_hrates', 'wrk', nw )
      end if

!------------------------------------------------------------------------------
!     ... derive Slant Path for Spherical Atmosphere
!------------------------------------------------------------------------------
      call sphers( nlev, zkm, zen, dsdh, nid )

!------------------------------------------------------------------------------
!  	... derive O2, O3, and CO2 slant columns
!------------------------------------------------------------------------------
      delz(1:nlev-1) = km2cm*(zkm(1:nlev-1) - zkm(2:nlev))
      call slant_col( nlev, delz, dsdh, nid, o2cc, o2scol )
      call slant_col( nlev, delz, dsdh, nid, o3cc, o3scol )
      call slant_col( nlev, delz, dsdh, nid, co2cc, scco2 )

!------------------------------------------------------------------------------
!  	... transmission due to ozone
!------------------------------------------------------------------------------
      do wn = 1,nw
         abs_col(:)     = min( (xs_o3a(wn) + xs_o3b(wn))*o3scol(:),100._r8 )
         trans_o3(:,wn) = exp( -abs_col(:) )
      end do

!------------------------------------------------------------------------------
!	... derive the cross section and transmission for
!           molecular oxygen Lya, SRC, and SRB's
!------------------------------------------------------------------------------
!  	... transmission due to molecular oxygen in the SRC
!------------------------------------------------------------------------------
      do wn = 1,nsrc_tot
         abs_col(:) = min( xs_o2src(wn)*o2scol(:),100._r8 )
         trans_o2(:,wn) = exp( -abs_col(:) )
      end do

!------------------------------------------------------------------------------
!     	... transmission and cross sections due to O2 at lyman alpha
!------------------------------------------------------------------------------
      call lymana( nlev, o2scol, rmla, ro2la )

!------------------------------------------------------------------------------
!    	... place lya reduction faction in transmission array
!           this must follow the SRC placement (above)
!------------------------------------------------------------------------------
      trans_o2(:,1) = rmla(:)

!------------------------------------------------------------------------------
!     ... Molecular Oxygen, SRB
!------------------------------------------------------------------------------
!     ... Koppers Grid (see Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996)
!        #    wl(i)       wl(i+1)
!        1   174.4        177.0
!        2   177.0        178.6
!        3   178.6        180.2
!        4   180.2        181.8
!        5   181.8        183.5
!        6   183.5        185.2
!        7   185.2        186.9
!        8   186.9        188.7
!        9   188.7        190.5
!        10  190.5        192.3
!        11  192.3        194.2
!        12  194.2        196.1
!        13  196.1        198.0
!        14  198.0        200.0  <<last wl bin for <200nm
!        ----------------------
!        15  200.0        202.0
!        16  202.0        204.1
!        17  204.1        205.8
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     ... Kopper SRB parameterization is used here
!------------------------------------------------------------------------------
      tlevmin(nlev:1:-1) = min( max( tlev(:),150._r8 ),350._r8 )
      call calc_o2srb( nlev, nid, o2scol, tlevmin, tsrb, xs_o2srb )

!------------------------------------------------------------------------------
!     ... Place Koppers SRB transmission (1-14) on user
!         grid #20 (wc=175.7nm)
!------------------------------------------------------------------------------
      do wn = 1,nsrb_tot
         trans_o2(:,wn+nsrc_tot) = tsrb(:,wn)
      end do

!------------------------------------------------------------------------------
!     ... Derive the normalize flux at each altitude,
!         corrected for O2 and O3 absorption
!------------------------------------------------------------------------------
      do wn = 1,nw                                  ! nw = 33 (nsrb_tot+nsrc_tot)
         fnorm(:,wn) = etfphot(wn)*trans_o2(:,wn)*trans_o3(:,wn)
      end do

!------------------------------------------------------------------------------
!     ... Derive the O2 rate constant and apply branching ratio (QY)
!------------------------------------------------------------------------------
!     ... SRC and SRB QY
!         Longward  of 174.65 the product is O2 + hv => O(3P) + O(3P)
!         Shortward of 174.65 the product is O2 + hv => O(3P) + O(1D)
!         The QY is assumed to be unity in both wavelength ranges.
!
!     ... Lyman Alpha QY
!         O2 + hv -> O(3P) + O(3P) at Lyman Alpha has a QY = 0.47
!         O2 + hv -> O(3P) + O(1D) at Lyman Alpha has a QY = 0.53
!         Lacoursiere et al., J. Chem. Phys. 110., 1949-1958, 1999.
!------------------------------------------------------------------------------
!     ... lyman alpha
!------------------------------------------------------------------------------
      jo2_lya(:) = etfphot(1)*ro2la(:)*wlintv(1)

      wrk(1:nsrc_tot) = xs_o2src(1:nsrc_tot)*wlintv(1:nsrc_tot) &
					    *bde_o2_a(1:nsrc_tot)
      wrk(1)          = 0._r8
!------------------------------------------------------------------------------
!     ... o2 src heating
!------------------------------------------------------------------------------
      if( do_diag ) then
      write(iulog,*) '-------------------------------------------------'
      write(iulog,*) 'jshort_hrates: fnorm,wrk at long,lchnk = ',long,lchnk
      write(iulog,'(1p,5g12.5)') fnorm(nlev,1:nsrc_tot)
      write(iulog,*) ' '
      write(iulog,'(1p,5g12.5)') wrk(1:nsrc_tot)
      write(iulog,*) '-------------------------------------------------'
      end if
#ifdef USE_ESSL
      call dgemm( 'N', 'N', nlev, 1, nsrc_tot, &
		  1._r8, fnorm, nlev, wrk, nw, &
		  0._r8, qrs, nlev )
#else
      qrs(:,1) = matmul( fnorm(:,1:nsrc_tot),wrk(1:nsrc_tot) )
#endif
!------------------------------------------------------------------------------
!     ... o2 srb heating
!------------------------------------------------------------------------------
      do k = 1,nlev
         wrk(1:nsrb_tot) = xs_o2srb(k,1:nsrb_tot)*wlintv(nsrc_tot+1:nsrc_tot+nsrb_tot) &
                                                 *bde_o2_b(nsrc_tot+1:nsrc_tot+nsrb_tot)
	 qrs(k,2)        = dot_product( fnorm(k,nsrc_tot+1:nsrc_tot+nsrb_tot),wrk(1:nsrb_tot) )
      end do

      if( do_diag ) then
      write(iulog,*) '-------------------------------------------------'
      write(iulog,*) 'jshort_hrates: lya,bde_o2_a,qrs(nlev) at long,lchnk = ',long,lchnk
      write(iulog,'(1p,5g12.5)') jo2_lya(nlev),bde_o2_a(2),qrs(nlev,1)
      write(iulog,*) '-------------------------------------------------'
      end if
!------------------------------------------------------------------------------
!     ... total o2 heating
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     ... Branch 1, O2 + hv => O(3P) + O(3P); wavelengths >175nm
!------------------------------------------------------------------------------
      qrs(:,2)  = qrs(:,2) + jo2_lya(:)*.47_r8*bde_o2_b(2)
!------------------------------------------------------------------------------
!     ... Branch 2, O2 + hv => O(3P) + O(1D);  wavelengths <175nm
!------------------------------------------------------------------------------
      qrs(:,1)  = qrs(:,1) + jo2_lya(:)*.53_r8*bde_o2_a(2)
      if( do_diag ) then
      write(iulog,*) '-------------------------------------------------'
      write(iulog,*) 'jshort_hrates: o2(1),qrs(nlev) at long,lchnk = ',long,lchnk
      write(iulog,'(1p,5g12.5)') o2_vmr(1),qrs(nlev,1)
      write(iulog,*) '-------------------------------------------------'
      end if

!------------------------------------------------------------------------------
!	... the o3 heating rates
!------------------------------------------------------------------------------
      wrk(:) = xs_o3a(:)*wlintv(:)*bde_o3_a(:)
#ifdef USE_ESSL
      call dgemm( 'N', 'N', nlev, 1, nw, &
		  1._r8, fnorm, nlev, wrk, nw, &
		  0._r8, qrs(1,3), nlev )
#else
      qrs(:,3) = matmul( fnorm,wrk )
#endif
      wrk(:) = xs_o3b(:)*wlintv(:)*bde_o3_b(:)
#ifdef USE_ESSL
      call dgemm( 'N', 'N', nlev, 1, nw, &
		  1._r8, fnorm, nlev, wrk, nw, &
		  0._r8, qrs(1,4), nlev )
#else
      qrs(:,4) = matmul( fnorm,wrk )
#endif

!------------------------------------------------------------------------------
!	... form actual heating rates (k/s)
!------------------------------------------------------------------------------
      do k = 1,nlev
	 k1       = nlev - k + 1
         hfactor  = avogad/(cparg(k1)*mw(k1))
         qrs(k,1) = qrs(k,1)*hfactor*o2_vmr(k1)
         qrs(k,2) = qrs(k,2)*hfactor*o2_vmr(k1)
         qrs(k,3) = qrs(k,3)*hfactor*o3_vmr(k1)
         qrs(k,4) = qrs(k,4)*hfactor*o3_vmr(k1)
      end do

      deallocate( fnorm, trans_o2, trans_o3, wrk )

      end subroutine jshort_hrates

      subroutine jshort_photo( nlev, zen, n2cc, o2cc, o3cc, &
                               nocc, tlev, zkm, jo2_sht, jno_sht, jsht )
!==============================================================================!
!   Subroutine Jshort                                                          !
!                                                                              !
!==============================================================================!
!   Purpose:                                                                   !
!     To calculate the total J for JO2, JNO, and selective species below 200nm.!
!                                                                              !
!==============================================================================!
!   This routine uses JO2 parameterizations based on:                          !
!        Lyman alpha... Chabrillat and Kockarts, GRL, 25, 2659, 1998           !
!        SRC .......... Brasseur and Solomon, 1986 (from TUV)                  !
!        SRB .......... Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996    !
!                        (supplied by Dan Marsh, NCAR ACD                      !
!   and JNO:                                                                   !
!        SRB .......... Minschwanner and Siskind, JGR< 98, 20401, 1993.        !
!                                                                              !
!==============================================================================!
!   Input:                                                                     !
!	n2cc....... N2 concentration, molecule cm-3			       !
!       o2cc....... O2 concentration, molecule cm-3			       !
!	o3cc....... O3 concentration, molecule cm-3			       !
!	nocc....... NO concentration, molecule cm-3			       !
!	n2cc....... N2 concentration, molecule cm-3			       !
!     zen........ zenith angle, units = degrees                                !
!     tlev....... Temperature Profile (K)                                      !
!     zkm ....... Altitude, km                                                 !
!                                                                              !
!   Output:                                                                    !
!     jo2_sht ... O2 photolytic rate constant, sec-1, <200nm                   !
!     jno_sht ... NO photolytic rate constant, sec-1, SRB                      !
!     jsht. Photolytic rate constant for other species below 200nm             !
!                                                                              !
!==============================================================================!
!                                                                              !
!   Approach:                                                                  !
!                                                                              !
!    1) Call sphers (taken from TUV)                                           !
!         -> derives dsdh and nid used in slant column routines                !
!	  -> zenith angle dependent                                            !
!                                                                              !
!    2) Call  slant_col (taken from TUV)                                       !
!		-> derives the slant column for each species                   !
!                                                                              !
!    3) Calls get_crs                                                          !
!		-> read a NetCDF file                                          !
!		-> returns cross sections*quantum yields for all species that  !
!		   have absorption below 200nm.                                !
!                                                                              !
!    4) Derives transmission and photolysis rates for selective species        !
!                                                                              !
!==============================================================================!
!   EDIT HISTORY:                                                              !
!   Created by Doug Kinnison, 3/14/2002                                        !
!==============================================================================!

        use error_messages, only : alloc_err

	implicit none

	integer, parameter :: branch	  = 2       ! two photolysis pathways for JO2
	real(r8), parameter    :: km2cm	  = 1.e5_r8

!------------------------------------------------------------------------------
!     ... dummy arguments
!------------------------------------------------------------------------------
	integer, intent(in)     :: nlev                 ! model vertical levels
	real(r8), intent(in)    :: zen	                ! Zenith angle (degrees)
        real(r8), intent(in)    :: n2cc(nlev)           ! Molecular Nitrogen conc (mol/cm^3)
	real(r8), intent(in)    :: o2cc(nlev)		! Molecular Oxygen conc (mol/cm^3)
	real(r8), intent(in)    :: o3cc(nlev)		! Ozone concentration (mol/cm^3)
        real(r8), intent(in)    :: nocc(nlev)		! Nitric Oxide conc (mol/cm^3)
	real(r8), intent(in)    :: tlev(nlev)		! Temperature profile
	real(r8), intent(in)    :: zkm(nlev)		! Altitude, km
	real(r8), intent(out)   :: jo2_sht(nlev,branch) ! JO2, sec-1, <200nm
	real(r8), intent(out)   :: jno_sht(nlev)        ! JNO, sec-1, SRB
	real(r8), intent(out)   :: jsht(:,:)	        ! Additional J's

!------------------------------------------------------------------------------
!     ... local variables
!------------------------------------------------------------------------------
	integer :: k 			     ! Altitude index
	integer :: wn			     ! Wavelength index
	integer :: astat
        integer :: nid(0:nlev)	             ! Number of layers crossed by the direct
                                             ! beam when travelling from the top of the
			                     ! atmosphere to layer i; NID(i), i = 0..NZ-1
        real(r8) :: hfactor
	real(r8) :: dsdh(0:nlev,nlev)        ! Slant path of direct beam through each
                                             ! layer crossed  when travelling from the top of
				             ! the atmosphere to layer i; DSDH(i,j), i = 0.
					     ! NZ-1, j = 1..NZ-1   (see sphers.f)
        real(r8), allocatable :: fnorm(:,:)      ! Normalized ETF
	real(r8), allocatable :: trans_o2(:,:)   ! Transmission o2 (total)
        real(r8), allocatable :: trans_o3(:,:)   ! Transmission, ozone
	real(r8), allocatable :: wrk(:)	     ! wrk array
	real(r8) :: jo2_lya(nlev)	     ! Total photolytic rate constant for Ly alpha
	real(r8) :: jo2_srb(nlev)	     ! Total JO2  for SRB
	real(r8) :: jo2_src(nlev)            ! Total JO2 for SRC
	real(r8) :: delz(nlev)               ! layer thickness (cm)
	real(r8) :: o2scol(nlev)             ! O2 Slant Column
        real(r8) :: o3scol(nlev)	     ! O3 Slant Column
        real(r8) :: noscol(nlev)	     ! NO Slant Column
        real(r8) :: rmla(nlev)		     ! Transmission, Lyman Alpha (other species)
	real(r8) :: ro2la(nlev)		     ! Transmission, Lyman Alpha (for JO2)
	real(r8) :: tlevmin(nlev)
	real(r8) :: abs_col(nlev)
        real(r8) :: tsrb(nlev,nsrbtuv)       ! Transmission in the SRB
	real(r8) :: xs_o2srb(nlev,nsrbtuv)   ! Cross section * QY for O2 in SRB

      allocate( fnorm(nlev,nw),stat=astat )
      if( astat /= 0 ) then
	 call alloc_err( astat, 'jshort_photo', 'fnorm', nw*nlev )
      end if
      allocate( trans_o2(nlev,nw),stat=astat )
      if( astat /= 0 ) then
	 call alloc_err( astat, 'jshort_photo', 'trans_o2', nw*nlev )
      end if
      allocate( trans_o3(nlev,nw),stat=astat )
      if( astat /= 0 ) then
	 call alloc_err( astat, 'jshort_photo', 'trans_o3', nw*nlev )
      end if
      allocate( wrk(nw),stat=astat )
      if( astat /= 0 ) then
	 call alloc_err( astat, 'jshort_photo', 'wrk', nw )
      end if

!------------------------------------------------------------------------------
!     ... Derive Slant Path for Spherical Atmosphere
!------------------------------------------------------------------------------
      call sphers( nlev, zkm, zen, dsdh, nid )

!------------------------------------------------------------------------------
!     ... Derive O2, O3, and NO Slant Column
!------------------------------------------------------------------------------
      delz(1:nlev-1) = km2cm*(zkm(1:nlev-1) - zkm(2:nlev))
      call slant_col( nlev, delz, dsdh, nid, o2cc, o2scol )
      call slant_col( nlev, delz, dsdh, nid, o3cc, o3scol )
      call slant_col( nlev, delz, dsdh, nid, nocc, noscol )

!------------------------------------------------------------------------------
!     ... Transmission due to ozone
!------------------------------------------------------------------------------
      do wn = 1,nw
         abs_col(:)     = min( (xs_o3a(wn) + xs_o3b(wn))*o3scol(:),100._r8 )
	 trans_o3(:,wn) = exp( -abs_col(:) )
      end do

!------------------------------------------------------------------------------
!    ... Derive the cross section and transmission for
!        molecular oxygen Lya, SRC, and SRB's
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     ... Transmission due to molecular oxygen in the SRC
!------------------------------------------------------------------------------
      do wn = 1,nsrc_tot
         abs_col(:) = min( xs_o2src(wn)*o2scol(:),100._r8 )
	 trans_o2(:,wn) = exp( -abs_col(:) )
      end do

!------------------------------------------------------------------------------
!     ... Transmission and cross sections due to O2 at lyman alpha
!------------------------------------------------------------------------------
      call lymana( nlev, o2scol, rmla, ro2la )

!------------------------------------------------------------------------------
!    ... Place lya reduction faction in transmission array
!    ... This must follow the SRC placement (above)
!------------------------------------------------------------------------------
      trans_o2(:,1) = rmla(:)

!------------------------------------------------------------------------------
!     ... Molecular Oxygen, SRB
!------------------------------------------------------------------------------
!     ... Koppers Grid (see Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996)
!        #    wl(i)       wl(i+1)
!        1   174.4        177.0
!        2   177.0        178.6
!        3   178.6        180.2
!        4   180.2        181.8
!        5   181.8        183.5
!        6   183.5        185.2
!        7   185.2        186.9
!        8   186.9        188.7
!        9   188.7        190.5
!        10  190.5        192.3
!        11  192.3        194.2
!        12  194.2        196.1
!        13  196.1        198.0
!        14  198.0        200.0  <<last wl bin for <200nm
!        ----------------------
!        15  200.0        202.0
!        16  202.0        204.1
!        17  204.1        205.8
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     ... Kopper SRB parameterization is used here
!------------------------------------------------------------------------------
      tlevmin(nlev:1:-1) = min( max( tlev(:),150._r8 ),350._r8 )
      call calc_o2srb( nlev, nid, o2scol, tlevmin, tsrb, xs_o2srb )

!------------------------------------------------------------------------------
!     ... Place Koppers SRB transmission (1-14) on user
!         grid #20 (wc=175.7nm)
!------------------------------------------------------------------------------
      do wn = 1,nsrb_tot
         trans_o2(:,wn+nsrc_tot) = tsrb(:,wn)
      end do

!------------------------------------------------------------------------------
!     ... Derive the normalize flux at each altitude,
!         corrected for O2 and O3 absorption
!------------------------------------------------------------------------------
      do wn = 1,nw                               ! nw = 33 (nsrb_tot+nsrc_tot)
         fnorm(:,wn) = etfphot(wn)*trans_o2(:,wn)*trans_o3(:,wn)
      end do

!------------------------------------------------------------------------------
!     ... Derive the O2 rate constant and apply branching ratio (QY)
!------------------------------------------------------------------------------
!     ... SRC and SRB QY
!         Longward  of 174.65 the product is O2 + hv => O(3P) + O(3P)
!         Shortward of 174.65 the product is O2 + hv => O(3P) + O(1D)
!         The QY is assumed to be unity in both wavelength ranges.
!
!     ... Lyman Alpha QY
!         O2 + hv -> O(3P) + O(3P) at Lyman Alpha has a QY = 0.47
!         O2 + hv -> O(3P) + O(1D) at Lyman Alpha has a QY = 0.53
!         Lacoursiere et al., J. Chem. Phys. 110., 1949-1958, 1999.
!------------------------------------------------------------------------------
!     ... Lyman Alpha
!------------------------------------------------------------------------------
      jo2_lya(:) = etfphot(1)*ro2la(:)*wlintv(1)

      wrk(1:nsrc_tot) = xs_o2src(1:nsrc_tot)*wlintv(1:nsrc_tot)
      wrk(1)          = 0._r8
!------------------------------------------------------------------------------
!     ... o2 src photolysis
!------------------------------------------------------------------------------
#ifdef USE_ESSL
      call dgemm( 'N', 'N', nlev, 1, nsrc_tot, &
		  1._r8, fnorm, nlev, wrk, nw, &
		  0._r8, jo2_src, nlev )
#else
      jo2_src(:) = matmul( fnorm(:,1:nsrc_tot),wrk(1:nsrc_tot) )
#endif
!------------------------------------------------------------------------------
!     ... o2 srb photolysis
!------------------------------------------------------------------------------
      do k = 1,nlev
         wrk(1:nsrb_tot) = xs_o2srb(k,1:nsrb_tot)*wlintv(nsrc_tot+1:nsrc_tot+nsrb_tot)
	 jo2_srb(k)      = dot_product( fnorm(k,nsrc_tot+1:nsrc_tot+nsrb_tot),wrk(1:nsrb_tot) )
      end do

!------------------------------------------------------------------------------
!     ... total o2 photolysis
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     ... Branch 1, O2 + hv => O(3P) + O(3P); wavelengths >175nm
!------------------------------------------------------------------------------
      jo2_sht(:,1) = jo2_lya(:)*.47_r8 + jo2_srb(:)
!------------------------------------------------------------------------------
!     ... Branch 2, O2 + hv => O(3P) + O(1D);  wavelengths <175nm
!------------------------------------------------------------------------------
      jo2_sht(:,2) = jo2_lya(:)*.53_r8 + jo2_src(:)

!------------------------------------------------------------------------------
!     ... Derive the NO rate constant Minsch. and Siskind, JGR, 98, 20401, 1993
!------------------------------------------------------------------------------
      call calc_jno( nlev, etfphot_ms93, n2cc, o2scol, o3scol, &
                     noscol, jno_sht )

!------------------------------------------------------------------------------
!    ... Derive addtional rate constants for species with wl < 200 nm.!
!        Temperature dependence of the cross sections are not included in this
!        version.
!------------------------------------------------------------------------------
#if defined USE_ESSL
      call dgemm( 'N', 'N', nlev, nj, nw, &
 	          1._r8, fnorm, nlev, xs_wl, nw, &
 	          0._r8, jsht, nlev )
#else
      jsht(:,:) = matmul( fnorm,xs_wl )
#endif

      deallocate( fnorm, trans_o2, trans_o3, wrk )

      end subroutine jshort_photo

      subroutine sphers( nlev, z, zenith_angle, dsdh, nid )
!=============================================================================!
!   Subroutine sphers                                                         !
!=============================================================================!
!   PURPOSE:                                                                  !
!   Calculate slant path over vertical depth ds/dh in spherical geometry.     !
!   Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model  !
!   for computing the radiation field available for photolysis and heating    !
!   at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)   !
!=============================================================================!
!   PARAMETERS:                                                               !
!   NZ      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   Z       - REAL, specified altitude working grid (km)                  (I) !
!   ZEN     - REAL, solar zenith angle (degrees)                          (I) !
!   DSDH    - REAL, slant path of direct beam through each layer crossed  (O) !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID     - INTEGER, number of layers crossed by the direct beam when   (O) !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Original: Taken By Doug Kinnison from Sasha Madronich, TUV Code, V4.1a,   !
!             on 1/1/02                                                       !
!=============================================================================!

        use physconst,    only : rearth

      implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in)  :: nlev              ! number model vertical levels
      integer, intent(out) :: nid(0:nlev)       ! see above
      real(r8), intent (in) :: zenith_angle		! zenith_angle
      real(r8), intent (in) :: z(nlev)		! geometric altitude (km)
      real(r8), intent (out) :: dsdh(0:nlev,nlev)   ! see above


!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      real(r8) :: radius
      real(r8) :: re
      real(r8) :: zenrad
      real(r8) :: rpsinz
      real(r8) :: const0
      real(r8) :: rj
      real(r8) :: rjp1
      real(r8) :: dsj
      real(r8) :: dhj
      real(r8) :: ga
      real(r8) :: gb
      real(r8) :: sm
      real(r8) :: zd(0:nlev-1)

      integer :: i
      integer :: j
      integer :: k
      integer :: id
      integer :: nlayer


      radius = rearth*1.e-3_r8   ! radius earth (km)

!------------------------------------------------------------------------------
!       ... set zenith angle in radians
!------------------------------------------------------------------------------
      zenrad = zenith_angle*d2r
      const0 = sin( zenrad )

!------------------------------------------------------------------------------
!       ... set number of layers:
!------------------------------------------------------------------------------
      nlayer = nlev - 1

!------------------------------------------------------------------------------
!       ... include the elevation above sea level to the radius of the earth:
!------------------------------------------------------------------------------
      re = radius + z(nlev)

!------------------------------------------------------------------------------
!       ... inverse coordinate of z
!------------------------------------------------------------------------------
      do k = 0,nlayer
        zd(k) = z(k+1) - z(nlev)
      end do

!------------------------------------------------------------------------------
!       ... initialize dsdh(i,j), nid(i)
!------------------------------------------------------------------------------
      nid(:) = 0
      do j = 1,nlev
        dsdh(:,j) = 0._r8
      end do

!------------------------------------------------------------------------------
!       ... calculate ds/dh of every layer
!------------------------------------------------------------------------------
      do i = 0,nlayer
        rpsinz = (re + zd(i)) * const0
        if( zenith_angle <= 90._r8 .or. rpsinz >= re ) then
!------------------------------------------------------------------------------
! Find index of layer in which the screening height lies
!------------------------------------------------------------------------------
           id = i 
           if( zenith_angle > 90._r8 ) then
              do j = 1,nlayer
                 if( rpsinz < (zd(j-1) + re) .and.  rpsinz >= (zd(j) + re) ) then
		    id = j
		    exit
		 end if
              end do
           end if
 
           do j = 1,id
             sm = 1._r8
             if( j == id .and. id == i .and. zenith_angle > 90._r8 ) then
                sm = -1._r8
             end if
             rj   = re + zd(j-1)
             rjp1 = re + zd(j)
             dhj  = zd(j-1) - zd(j)
             ga   = max( rj*rj - rpsinz*rpsinz,0._r8 )
             gb   = max( rjp1*rjp1 - rpsinz*rpsinz,0._r8 )
             if( id > i .and. j == id ) then
                dsj = sqrt( ga )
             else
                dsj = sqrt( ga ) - sm*sqrt( gb )
             end if
             dsdh(i,j) = dsj / dhj
           end do
           nid(i) = id
        else
           nid(i) = -1
        end if
      end do

      end subroutine sphers

      subroutine slant_col( nlev, delz, dsdh, nid, absden, scol )
!=============================================================================!
!   PURPOSE:                                                                  !
!   Derive Column
!=============================================================================!
!   PARAMETERS:                                                               !
!   NLEV   - INTEGER, number of specified altitude levels in the working  (I) !
!            grid                                                             !
!   DELZ   - REAL, specified altitude working grid (km)                   (I) !
!   DSDH   - REAL, slant path of direct beam through each layer crossed  (O)  !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID    - INTEGER, number of layers crossed by the direct beam when   (O)  !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!            specified altitude at each specified wavelength                  !
!   absden - REAL, absorber concentration, molecules cm-3                     !
!   SCOL   - REAL, absorber Slant Column, molecules cm-2                      !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   09/01  Read in profile from an input file, DEK                            !
!   01/02  Taken from Sasha Madronich's TUV code                              !
!=============================================================================!

      implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in) :: nlev
      integer, intent(in) :: nid(0:nlev)                ! see above
      real(r8), intent(in)    :: delz(nlev)	        ! layer thickness (cm)
      real(r8), intent(in)    :: dsdh(0:nlev,nlev)	! see above
      real(r8), intent(in)    :: absden(nlev)           ! absorber concentration (molec. cm-3)
      real(r8), intent(out)   :: scol(nlev)		! absorber Slant Column (molec. cm-2)

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      real(r8), parameter :: largest = 1.e+36_r8

      real(r8) :: sum
      real(r8) :: hscale
      real(r8) :: numer, denom
      real(r8) :: cz(nlev)

      integer :: id
      integer :: j
      integer :: k

!------------------------------------------------------------------------------
!     ... compute column increments (logarithmic integrals)
!------------------------------------------------------------------------------
      do k = 1,nlev-1
	if( absden(k) /= 0._r8 .and. absden(k+1) /= 0._r8 ) then
           cz(nlev-k) = (absden(k) - absden(k+1))/log( absden(k)/absden(k+1) ) * delz(k)
	else
           cz(nlev-k) = .5_r8*(absden(k) + absden(k+1)) * delz(k)
	end if
      end do

!------------------------------------------------------------------------------
!     ... Include exponential tail integral from infinity to model top
!         specify scale height near top of data.For WACCM-X model, scale
!         height needs to be increased for higher model top
!------------------------------------------------------------------------------
      if (nlev==pver) then
         if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
           hscale     = 20.e5_r8
         else
           hscale     = 10.e5_r8
         endif
         cz(nlev-1) = cz(nlev-1) + hscale * absden(1)
      endif

!------------------------------------------------------------------------------
!       ...  Calculate vertical and slant column from each level:
!            work downward
!------------------------------------------------------------------------------
      do id = 0,nlev-1
         sum = 0._r8
         if( nid(id) >= 0 ) then
!------------------------------------------------------------------------------
!       ...  Single pass layers:
!------------------------------------------------------------------------------
            do j = 1, min(nid(id), id)
               sum = sum + cz(nlev-j)*dsdh(id,j)
            end do
!------------------------------------------------------------------------------
!       ...  Double pass layers:
!------------------------------------------------------------------------------
            do j = min(nid(id),id)+1, nid(id)
               sum = sum + 2._r8*cz(nlev-j)*dsdh(id,j)
            end do
         else
            sum = largest
         end if
         scol(nlev-id) = sum
      end do
      scol(nlev) = .95_r8*scol(nlev-1)

      end subroutine slant_col

      subroutine lymana( nlev, o2scol, rm, ro2 )
!-----------------------------------------------------------------------------!
!   PURPOSE:                                                                  !
!   Calculate the effective absorption cross section of O2 in the Lyman-Alpha !
!   bands and an effective O2 optical depth at all altitudes.  Parameterized  !
!   after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the   !
!   absorption of the solar Lyman-Alpha line, Geophysical Research Letters,   !
!   Vol.24, No.21, pp 2659-2662, 1997.                                        !
!-----------------------------------------------------------------------------!
!   PARAMETERS:                                                               !
!   nz      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   o2scol  - REAL, slant overhead O2 column (molec/cc) at each specified (I) !
!             altitude                                                        !
!   dto2la  - REAL, optical depth due to O2 absorption at each specified  (O) !
!             vertical layer                                                  !
!   xso2la  - REAL, molecular absorption cross section in LA bands        (O) !
!-----------------------------------------------------------------------------!
!   EDIT HISTORY:                                                             !
!   01/15/2002 Taken from Sasha Madronich's TUV Version 4.1a, Doug Kinnison   !                  !
!   01/15/2002 Upgraded to F90, DK                                            !
!-----------------------------------------------------------------------------!

      implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in) :: nlev
      real(r8), intent(in)    :: o2scol(nlev)
      real(r8), intent(out)   :: ro2(nlev)
      real(r8), intent(out)   :: rm(nlev)

!------------------------------------------------------------------------------
!     ... Local variables
!------------------------------------------------------------------------------
      real(r8),save :: b(3)
      real(r8),save :: c(3)
      real(r8),save :: d(3)
      real(r8),save :: e(3)

      data b / 6.8431e-01_r8, 2.29841e-01_r8,  8.65412e-02_r8 /, &
           c / 8.22114e-21_r8, 1.77556e-20_r8,  8.22112e-21_r8 /, &
           d / 6.0073e-21_r8, 4.28569e-21_r8,  1.28059e-20_r8 /, &
           e / 8.21666e-21_r8, 1.63296e-20_r8,  4.85121e-17_r8 /

      integer  :: i, k
      real(r8) :: wrk, term

!------------------------------------------------------------------------------
!     ... Calculate reduction factors at every altitude
!------------------------------------------------------------------------------
      do k = 1,nlev
	 wrk = 0._r8
         do i = 1,2 ! pc Dan Marsh
	    term = e(i)*o2scol(k)
	    if( term < 100._r8 ) then
               wrk = wrk + d(i) * exp( -term )
	    end if
         end do
	 ro2(k) = wrk
	 wrk = 0._r8
	 do i = 1,3
	    term = c(i)*o2scol(k)
	    if( term < 100._r8 ) then
               wrk = wrk + b(i) * exp( -term )
	    end if
         end do
	 rm(k) = wrk
      end do

      end subroutine lymana

      subroutine calc_o2srb( nlev, nid, o2col, tlev, tsrb, xscho2 )
!-----------------------------------------------------------------------------!
!   PURPOSE:                                                                  !
!   Calculate the equivalent absorption cross section of O2 in the SR bands.  !
!   The algorithm is based on parameterization of G.A. Koppers, and           !
!   D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                          !
!   Final values do include effects from the Herzberg continuum.              !
!-----------------------------------------------------------------------------!
!   PARAMETERS:                                                               !
!   NZ      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I) !
!             altitude                                                        !
!   TLEV    - tmeperature at each level                                   (I) !
!   TSRB    - REAL, transmission for the SRB                                  !
!   XSCHO2  - REAL, molecular absorption cross section in SR bands at     (O) !
!             each specified wavelength.  Includes Herzberg continuum         !
!-----------------------------------------------------------------------------!
!   EDIT HISTORY: Taken from TUV, 1/17/2002                                   !
!   This code was supplied to TUV by Dan Marsh.                               !
!-----------------------------------------------------------------------------!

      implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in)     :: nlev
      integer, intent(in)     :: nid(0:nlev)
      real(r8), intent (in)   :: o2col(nlev)
      real(r8), intent (in)   :: tlev(nlev)
      real(r8), intent (out)  :: tsrb(nlev,nsrbtuv)
      real(r8), intent (out)  :: xscho2(nlev,nsrbtuv)

!------------------------------------------------------------------------------
!     ... Local variables
!------------------------------------------------------------------------------
      integer     :: i, k, ktop, ktop1, kbot
      real(r8)    :: x, dto2
      real(r8)    :: den, num
      real(r8)    :: term1, term2
      real(r8)    :: dtsrb(nlev)
      real(r8)    :: tsrb_rev(nlev,nsrbtuv)
      real(r8)    :: xs(nsrbtuv)

!------------------------------------------------------------------------------
!     ... Calculate cross sections
!------------------------------------------------------------------------------
      ktop = nlev
      kbot = 0

      do k = 1,nlev
	x  = log( o2col(k) )
	if( x >= 38._r8 .and. x <= 56._r8 ) then
          call effxs( x, tlev(k), xs )
          xscho2(k,:) = xs(:)
	else if( x < 38._r8 ) then
	   ktop1 = k-1
           ktop  = min( ktop1,ktop )
	else if( x > 56._r8 ) then
	   kbot = k
	end if
      end do
      
      if ( kbot == nlev ) then
         tsrb(:,:) = 0._r8
         xscho2(:,:) = 0._r8
         return
      endif
!------------------------------------------------------
!     ... Fill in cross section where X is out of range
!         by repeating edge table values
!-------------------------------------------------------
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 

         ! Need to be careful with nlev values for kbot and ktop. 
         ! This was handled by Hanli Liu fix.
         if ( kbot < nlev ) then
            do k = 1,kbot
               xscho2(k,:) = xscho2(kbot+1,:)
            end do
            if (ktop < nlev) then
               do k = ktop+1,nlev
                  xscho2(k,:) = xscho2(ktop,:)
               end do
            else
               xscho2(nlev,:) = 2.0e-19_r8
            endif
         else
            do k = 1,kbot
               xscho2(k,:) = 2.0e-19_r8
            enddo
         endif

      else
         do k = 1,kbot
            xscho2(k,:) = xscho2(kbot+1,:)
         end do
         do k = ktop+1,nlev
            xscho2(k,:) = xscho2(ktop,:)
         end do
      endif

!-------------------------------------------------------
!     ... Calculate incremental optical depths
!-------------------------------------------------------
      do i = 1,nsrbtuv 
        do k = 1,nlev-1
	   if( nid(nlev-k) /= -1 ) then
!-------------------------------------------------------
!     ... Calculate an optical depth weighted by density
!-------------------------------------------------------
               num   = xscho2(k+1,i)*o2col(k+1) - xscho2(k,i)*o2col(k)
	       if( num == 0._r8 ) then
                  write(iulog,*) 'calc_o2srb : o2col(k:k+1),xscho2(k:k+1,i) = ',o2col(k:k+1),xscho2(k:k+1,i),' @ i,k = ',i,k
	       end if
               term1 = log( xscho2(k+1,i)/xscho2(k,i) )
               term2 = log( o2col(k+1)/o2col(k) )
	       if( term2 == 0._r8 ) then
                  write(iulog,*) 'calc_o2srb : o2col(k:k+1),xscho2(k:k+1,i) = ',o2col(k:k+1),xscho2(k:k+1,i),' @ i,k = ',i,k
	          call endrun
	       end if
               den = 1._r8 + log( xscho2(k+1,i)/xscho2(k,i) )/log( o2col(k+1)/o2col(k) )
               dto2 = abs(num/den)
	       if( dto2  < 100._r8 ) then
                  dtsrb(k) = exp( -dto2 )
	       else
                  dtsrb(k) = 0._r8
	       end if
          else
            dtsrb(k) = 0._r8
	  end if
        end do
!-----------------------------------------------
!     ... Calculate Transmission for SRB
!-----------------------------------------------
        if (nlev==pver) then  ! waccm
           tsrb(nlev,i) = 1._r8
           do k = nlev-1,1,-1
              tsrb(k,i) = tsrb(k+1,i)*dtsrb(k)
           end do
        else   ! cam-chem
           tsrb(nlev,i) = exp(-xscho2(nlev,i)*o2col(nlev))
           do k = nlev-1,1,-1
              tsrb(k,i) = tsrb(k+1,i)*dtsrb(k)
           end do
        endif

      end do

      end subroutine calc_o2srb

      subroutine effxs( x, t, xs )
!-------------------------------------------------------------
!     Subroutine for evaluating the effective cross section
!     of O2 in the Schumann-Runge bands using parameterization
!     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
!     68-79, 1996]
!
!     method:
!     ln(xs) = A(X)[T-220]+B(X)
!     X = log of slant column of O2
!     A,B calculated from chebyshev polynomial coeffs
!     AC and BC using NR routine chebev.  Assume interval
!     is 38<ln(NO2)<56.
!
!     Revision History:
!
!     drm 2/97  initial coding
!-------------------------------------------------------------

      implicit none

!-------------------------------------------------------------
!	... Dummy arguments
!-------------------------------------------------------------
      real(r8), intent(in)  :: x
      real(r8), intent(in)  :: t
      real(r8), intent(out) :: xs(nsrbtuv)

!-------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------
      real(r8)  :: a(nsrbtuv), b(nsrbtuv)

      call calc_params( x, a, b )

      xs(:) = exp( a(:)*(t - 220._r8) + b(:) )

      end subroutine effxs


      subroutine calc_params( x, a, b )
!-------------------------------------------------------------
!       calculates coefficients (A,B), used in calculating the
!       effective cross section, for 17 wavelength intervals
!       as a function of log O2 column density (X)
!       Wavelength intervals are defined in WMO1985
!-------------------------------------------------------------

      implicit none

!-------------------------------------------------------------
!	... Dummy arguments
!-------------------------------------------------------------
      real(r8), intent(in)  :: x
      real(r8), intent(out) :: a(nsrbtuv), b(nsrbtuv)

!-------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------
      integer :: i

!-------------------------------------------------------------
!     ... call chebyshev evaluation routine to calc a and b from
!	    set of 20 coeficients for each wavelength
!-------------------------------------------------------------
      do i = 1,nsrbtuv
        a(i) = jchebev( 38._r8, 56._r8, ac(1,i), 20, x )
        b(i) = jchebev( 38._r8, 56._r8, bc(1,i), 20, x )
      end do

      contains

      function jchebev( a, b, c, m, x )
!-------------------------------------------------------------
!     Chebyshev evaluation algorithm
!     See Numerical recipes p193
!-------------------------------------------------------------

!-------------------------------------------------------------
!	... Dummy arguments
!-------------------------------------------------------------
      integer, intent(in)     :: m
      real(r8), intent(in)    :: a, b, x
      real(r8), intent(in)    :: c(m)

      real(r8) :: jchebev
!-------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------
      integer  :: j
      real(r8) :: d, dd, sv, y, y2

      if( (x - a)*(x - b) > 0._r8 ) then
	  write(iulog,*) 'x not in range in chebev', x
	  jchebev = 0._r8
	  return
      end if

      d  = 0._r8
      dd = 0._r8
      y  = (2._r8*x - a - b)/(b - a)
      y2 = 2._r8*y
      do j = m,2,-1
        sv = d
        d  = y2*d - dd + c(j)
        dd = sv
      end do

      jchebev = y*d - dd + .5_r8*c(1)

      end function jchebev

      end subroutine calc_params

      subroutine calc_jno( nlev, etfphot_ms93, n2cc, o2scol, o3scol, &
                           noscol, jno )
!-----------------------------------------------------------------------------!
!   PURPOSE:                                                                  !
!   Compute the total photolytic rate constant for NO in the SR bands         !
!     - following the approach of Minshwanner and Siskind, JGR,               !
!       98, D11, 20401-20412, 1993.                                           !
!                                                                             !
!-----------------------------------------------------------------------------!
!   PARAMETERS:                                                               !
!   NZ           - INTEGER, number of specified altitude levels               !
!                                                                             !
!   etfphot_ms93 - Extraterrestrial Flux, within the MS 1993 Grid             !
!                  units of photons cm-2 sec-1 nm-1                           !
!   n2cc         - N2 conc (molecules cm-3)                                   !
!   o3scol       - Ozone Slant Column (molecules cm-2)                        !
!   o2scol       - Oxygen Slant Column (molecules cm-2)                       !
!   noscol       - Nitric Oxide Slant Column(molecules cm-2)                  !
!                                                                             !
!   LOCAL VARIABLES:                                                          !
!   tauo3        - Transmission factor in the Hartley Band of O3              !
!   etfphot_ms93 - Solar Irr. on Minschwaner and Siskind 1993 (MS93) Grid     !
!   xs_o3ms93    - O3 cross section on the MS93 Grid                          !
!                                                                             !
!   OUTPUT VARIABLES:                                                         !
!   jno          - photolytic rate constant                                   !
!                  each specified altitude                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
!   EDIT HISTORY:                                                             !
!   08/01  Created, Doug Kinnison, NCAR, ACD                                  !
!-----------------------------------------------------------------------------!

      implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in) :: nlev
      real(r8), intent(in)    :: etfphot_ms93(num_ms93tuv)
      real(r8), intent(in)    :: n2cc(nlev)
      real(r8), intent(in)    :: o3scol(nlev)
      real(r8), intent(in)    :: o2scol(nlev)
      real(r8), intent(in)    :: noscol(nlev)
      real(r8), intent(out)   :: jno(nlev)

!------------------------------------------------------------------------------
!	... Local variables
!------------------------------------------------------------------------------
      integer     :: i, iw, lev
      real(r8)    :: jno50
      real(r8)    :: jno90
      real(r8)    :: jno100
      real(r8)    :: tauo3(nlev,num_ms93tuv)

!------------------------------------------------------------------------------
!   	... O3 SRB Cross Sections from WMO 1985, interpolated onto MS, 1993 grid
!------------------------------------------------------------------------------
      real(r8), save :: xso3_ms93(num_ms93tuv) = (/ 7.3307600e-19_r8, 6.9660105E-19_r8, 5.9257699E-19_r8, 4.8372219E-19_r8 /)

!------------------------------------------------------------------------------
!   	... delta wavelength of the MS, 1993 grid
!------------------------------------------------------------------------------
      real(r8), save :: wlintv_ms93(num_ms93tuv) = (/ 1.50_r8, 1.50_r8, 5.6_r8, 2.3_r8 /)

!------------------------------------------------------------------------------
!   	... O2 SRB Cross Sections for the six ODF regions, MS, 1993
!------------------------------------------------------------------------------
      real(r8), save :: cs250(6)  = (/ 1.117e-23_r8, 2.447e-23_r8, 7.188e-23_r8, 3.042e-22_r8, 1.748e-21_r8, 1.112e-20_r8 /)
      real(r8), save :: cs290(6)  = (/ 1.350e-22_r8, 2.991e-22_r8, 7.334e-22_r8, 3.074e-21_r8, 1.689e-20_r8, 1.658e-19_r8 /)
      real(r8), save :: cs2100(6) = (/ 2.968e-22_r8, 5.831e-22_r8, 2.053e-21_r8, 8.192e-21_r8, 4.802e-20_r8, 2.655e-19_r8 /)

!------------------------------------------------------------------------------
!     ... derive tauo3 for the three o2 srb
!     ... iw = 1,2, and 4 are used below for jno
!------------------------------------------------------------------------------
      do iw = 1,num_ms93tuv
         tauo3(:,iw) = exp( -xso3_ms93(iw)*o3scol(:) )
      end do

!------------------------------------------------------------------------------
!   	... Call PJNO Function to derive SR Band JNO contributions
!         Called in order of wavelength interval (shortest firs)
!------------------------------------------------------------------------------
      do lev = 1,nlev
	 jno100   = pjno( 1, cs2100, wtno100, csno100 )
         jno90    = pjno( 2, cs290,  wtno90,  csno90 )
         jno50    = pjno( 4, cs250,  wtno50,  csno50 )
         jno(lev) = jno50 + jno90 + jno100
      end do

      contains

      function pjno( w, cso2, wtno, csno )
!------------------------------------------------------------------------------
!   	... uses xsec at center of g subinterval for o2
!           uses mean values for no
!------------------------------------------------------------------------------
       implicit none

!------------------------------------------------------------------------------
!	... parameters
!------------------------------------------------------------------------------
      integer, parameter :: ngint = 6
      integer, parameter :: nno = 2

!----------------------------------------------------------------
!	... Dummy arguments
!----------------------------------------------------------------
      integer, intent(in)     :: w
      real(r8),    intent(in) :: cso2(ngint)
      real(r8),    intent(in) :: csno(ngint,nno)
      real(r8),    intent(in) :: wtno(ngint,nno)

!----------------------------------------------------------------
!	... Function declarations
!----------------------------------------------------------------
      real(r8) :: pjno

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer  ::  jj, i, k
      real(r8) :: tauno
      real(r8) :: transno
      real(r8) :: transo2
      real(r8) :: tauo2
      real(r8) :: jno
      real(r8) :: jno1

!----------------------------------------------------------------
!	... derive the photolysis frequency for no within a given
!         srb (i.e., 5-0, 9-0, 10-0)
!----------------------------------------------------------------
      jno = 0._r8
      do k = 1,ngint
	 tauo2 = o2scol(lev) * cso2(k)
	 if( tauo2 < 50._r8 ) then
	    transo2 = exp( -tauo2 )
	 else
	    transo2 = 0._r8
	 end if
         jno1 = 0._r8
         do jj = 1,nno
            tauno = noscol(lev)*csno(k,jj)
	    if( tauno < 50._r8 ) then
	       transno = exp( -tauno )
	    else
	       transno = 0._r8
	    end if
            jno1 = jno1 + csno(k,jj) * wtno(k,jj) * transno
         end do
         jno = jno + jno1*transo2
      end do

      pjno = wlintv_ms93(w)*etfphot_ms93(w)*tauo3(lev,w)*jno

!----------------------------------------------------------------
!	... correct for the predissociation of the deltq 1-0
!         transition in the srb (5-0)
!----------------------------------------------------------------
      if( w == 4 ) then
        pjno = 1.65e9_r8/(5.1e7_r8 + 1.65e9_r8 + (1.5e-9_r8*n2cc(nlev-lev+1)))*pjno
      end if

      end function pjno

      end subroutine calc_jno

      end module mo_jshort
