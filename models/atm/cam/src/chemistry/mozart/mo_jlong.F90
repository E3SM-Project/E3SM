#ifdef AIX
#define USE_ESSL
#endif
#define USE_BDE

      module mo_jlong

      use shr_kind_mod, only : r4 => shr_kind_r4
      use shr_kind_mod, only : r8 => shr_kind_r8
      use cam_logfile,  only : iulog
      use abortutils,   only : endrun
#ifdef SPMD
      use mpishorthand, only : mpicom,mpiint,mpir8, mpilog, mpir4
#endif
      use spmd_utils,   only : masterproc

      implicit none

      interface jlong
         module procedure jlong_photo
         module procedure jlong_hrates
      end interface

      private
      public :: jlong_init
      public :: jlong_timestep_init
      public :: jlong
      public :: numj

      save

      real(r8), parameter :: hc      = 6.62608e-34_r8 * 2.9979e8_r8 / 1.e-9_r8
      real(r8), parameter :: wc_o2_b = 242.37_r8   ! (nm)
      real(r8), parameter :: wc_o3_a = 310.32_r8   ! (nm)
      real(r8), parameter :: wc_o3_b = 1179.87_r8  ! (nm)

      integer               :: nw      		! wavelengths >200nm
      integer               :: nt      		! number of temperatures in xsection table
      integer               :: np_xs   		! number of pressure levels in xsection table
      integer               :: numj    		! number of photorates in xsqy, rsf
      integer               :: nump    		! number of altitudes in rsf
      integer               :: numsza  		! number of zen angles in rsf
      integer               :: numalb  		! number of albedos in rsf
      integer               :: numcolo3		! number of o3 columns in rsf
      real(r4), allocatable :: xsqy(:,:,:,:)
      real(r8), allocatable :: wc(:)
      real(r8), allocatable :: we(:)
      real(r8), allocatable :: wlintv(:)
      real(r8), allocatable :: etfphot(:)
      real(r8), allocatable :: bde_o2_b(:)
      real(r8), allocatable :: bde_o3_a(:)
      real(r8), allocatable :: bde_o3_b(:)
      real(r8), allocatable :: xs_o2b(:,:,:)
      real(r8), allocatable :: xs_o3a(:,:,:)
      real(r8), allocatable :: xs_o3b(:,:,:)
      real(r8), allocatable :: p(:)
      real(r8), allocatable :: del_p(:)
      real(r8), allocatable :: prs(:)
      real(r8), allocatable :: dprs(:)
      real(r8), allocatable :: sza(:)
      real(r8), allocatable :: del_sza(:)
      real(r8), allocatable :: alb(:)
      real(r8), allocatable :: del_alb(:)
      real(r8), allocatable :: o3rat(:)
      real(r8), allocatable :: del_o3rat(:)
      real(r8), allocatable :: colo3(:)
      real(r4), allocatable :: rsf_tab(:,:,:,:,:)
      logical :: jlong_used = .false.

      contains

      subroutine jlong_init( xs_long_file, rsf_file, lng_indexer )

      use ppgrid,         only : pver
      use time_manager,   only : is_end_curr_day
      use mo_util,        only : rebin
      use solar_data,  only : data_nw => nbins, data_we => we, data_etf => sol_etf

      implicit none

!------------------------------------------------------------------------------
!    ... dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout)       :: lng_indexer(:)
      character(len=*), intent(in) :: xs_long_file, rsf_file

!------------------------------------------------------------------------------
!     ... read Cross Section * QY NetCDF file
!         find temperature index for given altitude
!         derive cross*QY results, returns xsqy(nj,nz,nw)
!------------------------------------------------------------------------------
      call get_xsqy( xs_long_file, lng_indexer )

!------------------------------------------------------------------------------
!     ... read radiative source function NetCDF file
!------------------------------------------------------------------------------
      if(masterproc) write(iulog,*) 'jlong_init: before get_rsf'
      call get_rsf(rsf_file)
      if(masterproc) write(iulog,*) 'jlong_init: after  get_rsf'

      we(:nw)  = wc(:nw) - .5_r8*wlintv(:nw)
      we(nw+1) = wc(nw) + .5_r8*wlintv(nw)

      if (masterproc) then
         write(iulog,*) ' '
         write(iulog,*) '--------------------------------------------------'
      endif
      call rebin( data_nw, nw, data_we, we, data_etf, etfphot )
      if (masterproc) then
         write(iulog,*) 'jlong_init: etfphot after data rebin'
         write(iulog,'(1p,5g15.7)') etfphot(:)
         write(iulog,*) '--------------------------------------------------'
         write(iulog,*) ' '
      endif

      jlong_used = .true.
 
      end subroutine jlong_init

      subroutine get_xsqy( xs_long_file, lng_indexer )
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:                                        !
!     cross section * QY temperature dependence, >200nm                       !
!                                                                             !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      filepath.... NetCDF filepath that contains the "cross sections"        !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 3/14/2002                                       !
!=============================================================================!

      use ioFileMod,      only : getfil
      use error_messages, only : alloc_err
      use chem_mods,      only : phtcnt, pht_alias_lst, rxt_tag_lst
      use netcdf, only: &
           nf90_open, &
           nf90_nowrite, &
           nf90_inq_dimid, &
           nf90_inquire_dimension, &
           nf90_inq_varid, &
           nf90_noerr, &
           nf90_get_var, &
           nf90_close

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout) :: lng_indexer(:)
      character(len=*) :: xs_long_file

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      integer :: varid, dimid, ndx
      integer :: ncid
      integer :: iret
      integer :: i, k, m, n
      integer :: wrk_ndx(phtcnt)
      character(len=256) :: locfn

      Masterproc_only : if( masterproc ) then
         !------------------------------------------------------------------------------
         !       ... open NetCDF File
         !------------------------------------------------------------------------------
         call getfil(xs_long_file, locfn, 0)
         iret = nf90_open(trim(locfn), NF90_NOWRITE, ncid)

         !------------------------------------------------------------------------------
         !       ... get dimensions
         !------------------------------------------------------------------------------
         iret = nf90_inq_dimid( ncid, 'numtemp', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= nt )
         iret = nf90_inq_dimid( ncid, 'numwl', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= nw )
         iret = nf90_inq_dimid( ncid, 'numprs', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= np_xs )
         !------------------------------------------------------------------------------
         !       ... check for cross section in dataset
         !------------------------------------------------------------------------------
         do m = 1,phtcnt
            if( pht_alias_lst(m,2) == ' ' ) then
               iret = nf90_inq_varid( ncid, rxt_tag_lst(m), varid )
               if( iret == nf90_noerr ) then 
                  lng_indexer(m) = varid
               end if
            else if( pht_alias_lst(m,2) == 'userdefined' ) then
               lng_indexer(m) = -1
            else
               iret = nf90_inq_varid( ncid, pht_alias_lst(m,2), varid )
               if( iret == nf90_noerr ) then 
                  lng_indexer(m) = varid
               else
                  write(iulog,*) 'get_xsqy : ',rxt_tag_lst(m)(:len_trim(rxt_tag_lst(m))),' alias ', &
                       pht_alias_lst(m,2)(:len_trim(pht_alias_lst(m,2))),' not in dataset'            
                  call endrun
               end if
            end if
         end do
         numj = 0
         do m = 1,phtcnt
            if( lng_indexer(m) > 0 ) then
               if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
                  cycle
               end if
               numj = numj + 1
            end if
         end do

         !------------------------------------------------------------------------------
         !       ... allocate arrays
         !------------------------------------------------------------------------------

         allocate( xsqy(numj,nw,nt,np_xs),stat=iret )
         if( iret /= 0 ) then 
            call alloc_err( iret, 'get_xsqy', 'xsqy', numj*nt*np_xs*nw )
         end if
         allocate( prs(np_xs),dprs(np_xs-1),stat=iret )
         if( iret /= 0 ) then 
            call alloc_err( iret, 'get_xsqy', 'prs,dprs', np_xs )
         end if
         allocate( xs_o2b(nw,nt,np_xs),xs_o3a(nw,nt,np_xs),xs_o3b(nw,nt,np_xs),stat=iret )
         if( iret /= 0 ) then 
            call alloc_err( iret, 'get_xsqy', 'xs_o2b ... xs_o3b', np_xs )
         end if
         !------------------------------------------------------------------------------
         !       ... read cross sections
         !------------------------------------------------------------------------------
         ndx = 0
         do m = 1,phtcnt
            if( lng_indexer(m) > 0 ) then
               if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
                  cycle
               end if
               ndx = ndx + 1
               iret = nf90_get_var( ncid, lng_indexer(m), xsqy(ndx,:,:,:) )
            end if
         end do
         if( ndx /= numj ) then
            write(iulog,*) 'get_xsqy : ndx count /= cross section count'
            call endrun
         end if
         iret = nf90_inq_varid( ncid, 'jo2_b', varid )
         iret = nf90_get_var( ncid, varid, xs_o2b )
         iret = nf90_inq_varid( ncid, 'jo3_a', varid )
         iret = nf90_get_var( ncid, varid, xs_o3a )
         iret = nf90_inq_varid( ncid, 'jo3_b', varid )
         iret = nf90_get_var( ncid, varid, xs_o3b )
         !------------------------------------------------------------------------------
         !       ... setup final lng_indexer
         !------------------------------------------------------------------------------
         ndx = 0
         wrk_ndx(:) = lng_indexer(:)
         do m = 1,phtcnt
            if( wrk_ndx(m) > 0 ) then
               ndx = ndx + 1
               i = wrk_ndx(m)
               where( wrk_ndx(:) == i )
                  lng_indexer(:) = ndx
                  wrk_ndx(:)     = -100000
               end where
            end if
         end do

         iret = nf90_inq_varid( ncid, 'pressure', varid )
         iret = nf90_get_var( ncid, varid, prs )
         iret = nf90_close( ncid )
      end if Masterproc_only

#ifdef SPMD
!        call mpibarrier( mpicom )
      call mpibcast( numj,  1, mpiint, 0, mpicom )
      call mpibcast( nt,    1, mpiint, 0, mpicom )
      call mpibcast( nw,    1, mpiint, 0, mpicom )
      call mpibcast( np_xs, 1, mpiint, 0, mpicom )
      call mpibcast( lng_indexer, phtcnt, mpiint, 0, mpicom )
#endif
      if( .not. masterproc ) then
         !------------------------------------------------------------------------------
         !       ... allocate arrays
         !------------------------------------------------------------------------------
         allocate( xsqy(numj,nw,nt,np_xs),stat=iret )
         if( iret /= nf90_noerr) then 
            write(iulog,*) 'get_xsqy : failed to allocate xsqy ; error = ',iret
            call endrun
         end if
         allocate( prs(np_xs),dprs(np_xs-1),stat=iret )
         if( iret /= nf90_noerr) then 
            write(iulog,*) 'get_xsqy : failed to allocate prs,dprs ; error = ',iret
            call endrun
         end if
         allocate( xs_o2b(nw,nt,np_xs),xs_o3a(nw,nt,np_xs),xs_o3b(nw,nt,np_xs),stat=iret )
         if( iret /= 0 ) then 
            call alloc_err( iret, 'get_xsqy', 'xs_o2b ... xs_o3b', np_xs )
         end if
      end if
#ifdef SPMD
!        call mpibarrier( mpicom )
      call mpibcast( prs, np_xs, mpir8, 0, mpicom )
      call mpibcast( xsqy, numj*nt*np_xs*nw, mpir4, 0, mpicom )
      call mpibcast( xs_o2b, nw*nt*np_xs, mpir8, 0, mpicom )
      call mpibcast( xs_o3a, nw*nt*np_xs, mpir8, 0, mpicom )
      call mpibcast( xs_o3b, nw*nt*np_xs, mpir8, 0, mpicom )
#endif
      dprs(:np_xs-1) = 1._r8/(prs(1:np_xs-1) - prs(2:np_xs))

      end subroutine get_xsqy

      subroutine get_rsf(rsf_file)
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:
!     Radiative Souce function                                                !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      filepath.... NetCDF file that contains the RSF                         !
!                                                                             !
!     Output:                                                                 !
!      rsf ........ Radiative Source Function (quanta cm-2 sec-1              !
!                                                                             !
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 3/14/2002                                       !
!=============================================================================!

      use ioFileMod,      only : getfil
      use error_messages, only : alloc_err
      use netcdf, only: &
           nf90_open, &
           nf90_nowrite, &
           nf90_inq_dimid, &
           nf90_inquire_dimension, &
           nf90_inq_varid, &
           nf90_noerr, &
           nf90_get_var, &
           nf90_close

!------------------------------------------------------------------------------
!    ... dummy arguments
!------------------------------------------------------------------------------
      character(len=*) :: rsf_file

!------------------------------------------------------------------------------
!       ... local variables
!------------------------------------------------------------------------------
      integer :: varid, dimid
      integer :: ncid
      integer :: i, j, k, l, w
      integer :: iret
      integer :: count(5)
      integer :: start(5)
      real(r8) :: wrk
      character(len=256) :: locfn

      Masterproc_only : if( masterproc ) then
         !------------------------------------------------------------------------------
         !       ... open NetCDF File
         !------------------------------------------------------------------------------
         call getfil(rsf_file, locfn, 0)
         iret = nf90_open(trim(locfn), NF90_NOWRITE, ncid)

         !------------------------------------------------------------------------------
         !       ... get dimensions
         !------------------------------------------------------------------------------
         iret = nf90_inq_dimid( ncid, 'numz', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= nump )
         iret = nf90_inq_dimid( ncid, 'numsza', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= numsza )
         iret = nf90_inq_dimid( ncid, 'numalb', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= numalb )
         iret = nf90_inq_dimid( ncid, 'numcolo3fact', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= numcolo3 )

      end if Masterproc_only
#ifdef SPMD
!        call mpibarrier( mpicom )
      call mpibcast( nump,     1, mpiint, 0, mpicom )
      call mpibcast( numsza,   1, mpiint, 0, mpicom )
      call mpibcast( numalb,   1, mpiint, 0, mpicom )
      call mpibcast( numcolo3, 1, mpiint, 0, mpicom )
#endif
!------------------------------------------------------------------------------
!       ... allocate arrays
!------------------------------------------------------------------------------
      allocate( wc(nw),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'wc', nw )
      end if
      allocate( wlintv(nw),we(nw+1),etfphot(nw),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'wlintv,etfphot', nw )
      end if
      allocate( bde_o2_b(nw),bde_o3_a(nw),bde_o3_b(nw),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'bde', nw )
      end if
      allocate( p(nump),del_p(nump-1),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'p,delp', nump )
      end if
      allocate( sza(numsza),del_sza(numsza-1),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'sza,del_sza', numsza )
      end if
      allocate( alb(numalb),del_alb(numalb-1),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'alb,del_alb', numalb )
      end if
      allocate( o3rat(numcolo3),del_o3rat(numcolo3-1),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'o3rat,del_o3rat', numcolo3 )
      end if
      allocate( colo3(nump),stat=iret )
      if( iret /= 0 ) then 
         call alloc_err( iret, 'get_rsf', 'colo3', nump )
      end if
      allocate( rsf_tab(nw,nump,numsza,numcolo3,numalb),stat=iret )
      if( iret /= 0 ) then 
         write(iulog,*) 'get_rsf : dimensions = ',nw,nump,numsza,numcolo3,numalb
         call alloc_err( iret, 'get_rsf', 'rsf_tab', numalb*numcolo3*numsza*nump )
      end if

      Masterproc_only2 : if( masterproc ) then
         !------------------------------------------------------------------------------
         !       ... read variables
         !------------------------------------------------------------------------------
         iret = nf90_inq_varid( ncid, 'wc', varid )
         iret = nf90_get_var( ncid, varid, wc )
         iret = nf90_inq_varid( ncid, 'wlintv', varid )
         iret = nf90_get_var( ncid, varid, wlintv )
         iret = nf90_inq_varid( ncid, 'pm', varid )
         iret = nf90_get_var( ncid, varid, p )
         iret = nf90_inq_varid( ncid, 'sza', varid )
         iret = nf90_get_var( ncid, varid, sza )
         iret = nf90_inq_varid( ncid, 'alb', varid )
         iret = nf90_get_var( ncid, varid, alb )
         iret = nf90_inq_varid( ncid, 'colo3fact', varid )
         iret = nf90_get_var( ncid, varid, o3rat )
         iret = nf90_inq_varid( ncid, 'colo3', varid )
         iret = nf90_get_var( ncid, varid, colo3 )

         iret = nf90_inq_varid( ncid, 'RSF', varid )
         
         if (masterproc) then
            write(iulog,*) ' '
            write(iulog,*) '----------------------------------------------'
            write(iulog,*) 'get_rsf: numalb, numcolo3, numsza, nump = ', &
                 numalb, numcolo3, numsza, nump
            write(iulog,*) 'get_rsf: size of rsf_tab = ',size(rsf_tab,dim=1),size(rsf_tab,dim=2), &
                 size(rsf_tab,dim=3),size(rsf_tab,dim=4),size(rsf_tab,dim=5)
            write(iulog,*) '----------------------------------------------'
            write(iulog,*) ' '
         endif

         iret = nf90_get_var( ncid, varid, rsf_tab )
         iret = nf90_close( ncid )

         do w = 1,nw
            wrk = wlintv(w)
            rsf_tab(w,:,:,:,:) = wrk*rsf_tab(w,:,:,:,:)
         enddo

      end if Masterproc_only2
#ifdef SPMD
      call mpibcast( wc,      nw,       mpir8, 0, mpicom )
      call mpibcast( wlintv,  nw,       mpir8, 0, mpicom )
      call mpibcast( p,       nump,     mpir8, 0, mpicom )
      call mpibcast( sza,     numsza,   mpir8, 0, mpicom )
      call mpibcast( alb,     numalb,   mpir8, 0, mpicom )
      call mpibcast( o3rat,   numcolo3, mpir8, 0, mpicom )
      call mpibcast( colo3,   nump,     mpir8, 0, mpicom )
      do w = 1,nw
         call mpibcast( rsf_tab(w,:,:,:,:), numalb*numcolo3*numsza*nump, mpir4, 0, mpicom )
      enddo
#endif
#ifdef USE_BDE
      if (masterproc) write(iulog,*) 'Jlong using bdes'
      bde_o2_b(:) = max( 0._r8, hc*(wc_o2_b - wc(:))/(wc_o2_b*wc(:)) )
      bde_o3_a(:) = max( 0._r8, hc*(wc_o3_a - wc(:))/(wc_o3_a*wc(:)) )
      bde_o3_b(:) = max( 0._r8, hc*(wc_o3_b - wc(:))/(wc_o3_b*wc(:)) )
#else
      if (masterproc) write(iulog,*) 'Jlong not using bdes'
      bde_o2_b(:) = hc/wc(:)
      bde_o3_a(:) = hc/wc(:)
      bde_o3_b(:) = hc/wc(:)
#endif

      del_p(:nump-1)         = 1._r8/abs(p(1:nump-1)- p(2:nump))
      del_sza(:numsza-1)     = 1._r8/(sza(2:numsza) - sza(1:numsza-1))
      del_alb(:numalb-1)     = 1._r8/(alb(2:numalb) - alb(1:numalb-1))
      del_o3rat(:numcolo3-1) = 1._r8/(o3rat(2:numcolo3) - o3rat(1:numcolo3-1))

      end subroutine get_rsf

      subroutine jlong_timestep_init
!---------------------------------------------------------------
!	... set etfphot if required
!---------------------------------------------------------------

      use time_manager,   only : is_end_curr_day
      use mo_util,        only : rebin

      use solar_data,  only : data_nw => nbins, data_we => we, data_etf => sol_etf

      implicit none

      if (.not.jlong_used) return

      call rebin( data_nw, nw, data_we, we, data_etf, etfphot )

      end subroutine jlong_timestep_init

      subroutine jlong_hrates( nlev, sza_in, alb_in, p_in, t_in, &
                               mw, o2_vmr, o3_vmr, colo3_in, qrl_col, &
                               cparg, kbot )
!==============================================================================
!   Purpose:                                                                   
!     To calculate the thermal heating rates longward of 200nm.        
!==============================================================================
!   Approach:
!     1) Reads the Cross Section*QY NetCDF file
!     2) Given a temperature profile, derives the appropriate XS*QY
!
!     3) Reads the Radiative Source function (RSF) NetCDF file
!        Units = quanta cm-2 sec-1
!
!     4) Indices are supplied to select a RSF that is consistent with
!        the reference atmosphere in TUV (for direct comparision of J's).
!        This approach will be replaced in the global model. Here colo3, zenith
!        angle, and altitude will be inputed and the correct entry in the table
!        will be derived.
!==============================================================================

	use physconst,       only : avogad
        use error_messages, only : alloc_err

	implicit none

!------------------------------------------------------------------------------
!    	... dummy arguments
!------------------------------------------------------------------------------
      integer, intent (in)     :: nlev               ! number vertical levels
      integer, intent (in)     :: kbot               ! heating levels
      real(r8), intent(in)     :: o2_vmr(nlev)       ! o2 conc (mol/mol)
      real(r8), intent(in)     :: o3_vmr(nlev)       ! o3 conc (mol/mol)
      real(r8), intent(in)     :: sza_in             ! solar zenith angle (degrees)
      real(r8), intent(in)     :: alb_in(nlev)       ! albedo
      real(r8), intent(in)     :: p_in(nlev)         ! midpoint pressure (hPa)
      real(r8), intent(in)     :: t_in(nlev)         ! Temperature profile (K)
      real(r8), intent(in)     :: colo3_in(nlev)     ! o3 column density (molecules/cm^3)
      real(r8), intent(in)     :: mw(nlev)           ! atms molecular weight
      real(r8), intent(in)     :: cparg(nlev)        ! specific heat capacity
      real(r8), intent(inout)  :: qrl_col(:,:)	     ! heating rates

!----------------------------------------------------------------------
!  	... local variables
!----------------------------------------------------------------------
      integer  ::  astat
      integer  ::  k, km, m
      integer  ::  t_index					! Temperature index
      integer  ::  pndx
      real(r8) ::  ptarget
      real(r8) ::  delp
      real(r8) ::  hfactor
      real(r8), allocatable                :: rsf(:,:)	        ! Radiative source function
      real(r8), allocatable                :: xswk(:,:)	        ! working xsection array
      real(r8), allocatable                :: wrk(:)	        ! work vector

!----------------------------------------------------------------------
!        ... allocate variables
!----------------------------------------------------------------------
      allocate( rsf(nw,kbot),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jlong_hrates', 'rsf', nw*nlev )
      end if
      allocate( xswk(nw,3),wrk(nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jlong_hrates', 'xswk,wrk', 3*nw )
      end if

!----------------------------------------------------------------------
!        ... interpolate table rsf to model variables
!----------------------------------------------------------------------
      call interpolate_rsf( alb_in, sza_in, p_in, colo3_in, kbot, rsf )

!------------------------------------------------------------------------------
!     ... calculate thermal heating rates for wavelengths >200nm
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     LLNL LUT approach to finding temperature index...
!     Calculate the temperature index into the cross section
!     data which lists coss sections for temperatures from
!     150 to 350 degrees K.  Make sure the index is a value
!     between 1 and 201.
!------------------------------------------------------------------------------
      qrl_col(kbot+1:nlev,:) = 0._r8
level_loop_1 : &
      do k = 1,kbot
!----------------------------------------------------------------------
!   	... get index into xsqy
!----------------------------------------------------------------------
        t_index = t_in(k) - 148.5_r8
        t_index = min( 201,max( t_index,1) )
!----------------------------------------------------------------------
!   	... find pressure level
!----------------------------------------------------------------------
	ptarget = p_in(k)
	if( ptarget >= prs(1) ) then
	   xswk(:,1) = xs_o2b(:,t_index,1)
	   xswk(:,2) = xs_o3a(:,t_index,1)
	   xswk(:,3) = xs_o3b(:,t_index,1)
	else if( ptarget <= prs(np_xs) ) then
	   xswk(:,1) = xs_o2b(:,t_index,np_xs)
	   xswk(:,2) = xs_o3a(:,t_index,np_xs)
	   xswk(:,3) = xs_o3b(:,t_index,np_xs)
	else
	   do km = 2,np_xs
	      if( ptarget >= prs(km) ) then
		 pndx = km - 1
		 delp = (prs(pndx) - ptarget)*dprs(pndx)
		 exit
	      end if
	   end do
	   xswk(:,1) = xs_o2b(:,t_index,pndx) &
                       + delp*(xs_o2b(:,t_index,pndx+1) - xs_o2b(:,t_index,pndx))
	   xswk(:,2) = xs_o3a(:,t_index,pndx) &
                       + delp*(xs_o3a(:,t_index,pndx+1) - xs_o3a(:,t_index,pndx))
	   xswk(:,3) = xs_o3b(:,t_index,pndx) &
                       + delp*(xs_o3b(:,t_index,pndx+1) - xs_o3b(:,t_index,pndx))
	end if
	hfactor      = avogad/(cparg(k)*mw(k))
	wrk(:)       = xswk(:,1)*bde_o2_b(:)
	qrl_col(k,2) = dot_product( wrk(:),rsf(:,k) ) * o2_vmr(k) * hfactor
	wrk(:)       = xswk(:,2)*bde_o3_a(:)
	qrl_col(k,3) = dot_product( wrk(:),rsf(:,k) ) * o3_vmr(k) * hfactor
	wrk(:)       = xswk(:,3)*bde_o3_b(:)
	qrl_col(k,4) = dot_product( wrk(:),rsf(:,k) ) * o3_vmr(k) * hfactor
      end do level_loop_1

      deallocate( rsf, xswk, wrk )

      end subroutine jlong_hrates

       subroutine jlong_photo( nlev, sza_in, alb_in, p_in, t_in, &
                              colo3_in, j_long )
!==============================================================================
!   Purpose:                                                                   
!     To calculate the total J for selective species longward of 200nm.        
!==============================================================================
!   Approach:
!     1) Reads the Cross Section*QY NetCDF file
!     2) Given a temperature profile, derives the appropriate XS*QY
!
!     3) Reads the Radiative Source function (RSF) NetCDF file
!        Units = quanta cm-2 sec-1
!
!     4) Indices are supplied to select a RSF that is consistent with
!        the reference atmosphere in TUV (for direct comparision of J's).
!        This approach will be replaced in the global model. Here colo3, zenith
!        angle, and altitude will be inputed and the correct entry in the table
!        will be derived.
!==============================================================================

 use spmd_utils,   only : masterproc
        use error_messages, only : alloc_err

	implicit none

!------------------------------------------------------------------------------
!    	... dummy arguments
!------------------------------------------------------------------------------
      integer, intent (in)     :: nlev               ! number vertical levels
      real(r8), intent(in)     :: sza_in             ! solar zenith angle (degrees)
      real(r8), intent(in)     :: alb_in(nlev)       ! albedo
      real(r8), intent(in)     :: p_in(nlev)         ! midpoint pressure (hPa)
      real(r8), intent(in)     :: t_in(nlev)         ! Temperature profile (K)
      real(r8), intent(in)     :: colo3_in(nlev)     ! o3 column density (molecules/cm^3)
      real(r8), intent(out)    :: j_long(:,:)	     ! photo rates (1/s)

!----------------------------------------------------------------------
!  	... local variables
!----------------------------------------------------------------------
      integer  ::  astat
      integer  ::  k, km, m
      integer  ::  wn
      integer  ::  t_index					! Temperature index
      integer  ::  pndx
      real(r8) ::  ptarget
      real(r8) ::  delp
      real(r8) ::  hfactor
      real(r8), allocatable :: rsf(:,:)	        ! Radiative source function
      real(r8), allocatable :: xswk(:,:)	! working xsection array

!----------------------------------------------------------------------
!        ... allocate variables
!----------------------------------------------------------------------
      allocate( rsf(nw,nlev),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jlong_photo', 'rsf', nw*nlev )
      end if
      allocate( xswk(numj,nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jlong_photo', 'xswk', numj*nw )
      end if

!----------------------------------------------------------------------
!        ... interpolate table rsf to model variables
!----------------------------------------------------------------------
      call interpolate_rsf( alb_in, sza_in, p_in, colo3_in, nlev, rsf )

!------------------------------------------------------------------------------
!     ... calculate total Jlong for wavelengths >200nm
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     LLNL LUT approach to finding temperature index...
!     Calculate the temperature index into the cross section
!     data which lists coss sections for temperatures from
!     150 to 350 degrees K.  Make sure the index is a value
!     between 1 and 201.
!------------------------------------------------------------------------------
level_loop_1 : &
      do k = 1,nlev
!----------------------------------------------------------------------
!   	... get index into xsqy
!----------------------------------------------------------------------
        t_index = t_in(k) - 148.5_r8
        t_index = min( 201,max( t_index,1) )
!----------------------------------------------------------------------
!   	... find pressure level
!----------------------------------------------------------------------
	ptarget = p_in(k)
	if( ptarget >= prs(1) ) then
	   do wn = 1,nw
	      xswk(:,wn) = xsqy(:,wn,t_index,1)
	   end do
	else if( ptarget <= prs(np_xs) ) then
	   do wn = 1,nw
	      xswk(:,wn) = xsqy(:,wn,t_index,np_xs)
	   end do
	else
	   do km = 2,np_xs
	      if( ptarget >= prs(km) ) then
		 pndx = km - 1
		 delp = (prs(pndx) - ptarget)*dprs(pndx)
		 exit
	      end if
	   end do
	   do wn = 1,nw
	      xswk(:,wn) = xsqy(:,wn,t_index,pndx) &
                           + delp*(xsqy(:,wn,t_index,pndx+1) - xsqy(:,wn,t_index,pndx))
	   end do
	end if
#ifdef USE_ESSL
        call dgemm( 'N', 'N', numj, 1, nw, &
		    1._r8, xswk, numj, rsf(1,k), nw, &
		    0._r8, j_long(1,k), numj )
#else
        j_long(:,k) = matmul( xswk(:,:),rsf(:,k) )
#endif
      end do level_loop_1

      deallocate( rsf, xswk )

      end subroutine jlong_photo

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!        ... interpolate table rsf to model variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine interpolate_rsf( alb_in, sza_in, p_in, colo3_in, kbot, rsf )

        use error_messages, only : alloc_err

      implicit none

!------------------------------------------------------------------------------
!    	... dummy arguments
!------------------------------------------------------------------------------
      real(r8), intent(in)  :: alb_in(:)       ! albedo
      real(r8), intent(in)  :: sza_in          ! solar zenith angle (degrees)
      integer,  intent(in)  :: kbot            ! heating levels
      real(r8), intent(in)  :: p_in(:)         ! midpoint pressure (hPa)
      real(r8), intent(in)  :: colo3_in(:)     ! o3 column density (molecules/cm^3)
      real(r8), intent(out) :: rsf(:,:)

!----------------------------------------------------------------------
!  	... local variables
!----------------------------------------------------------------------
      integer  ::  astat
      integer  ::  is, iv, ial
      integer  ::  isp1, ivp1, ialp1
      real(r8), dimension(3)               :: dels
      real(r8), dimension(0:1,0:1,0:1)     :: wghtl, wghtu
      real(r8) ::  psum_u
      real(r8), allocatable                :: psum_l(:)
      real(r8) ::  v3ratl, v3ratu
      integer  ::  pind, albind
      real(r8) ::  wrk0, wrk1, wght1
      integer  ::  iz, k, m
      integer  ::  izl
      integer  ::  ratindl, ratindu
      integer  ::  wn

!----------------------------------------------------------------------
!        ... allocate variables
!----------------------------------------------------------------------
      allocate( psum_l(nw),stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'jlong_hrates', 'psum_l', nw )
      end if

!----------------------------------------------------------------------
!        ... find the zenith angle index ( same for all levels )
!----------------------------------------------------------------------
      do is = 1,numsza
         if( sza(is) > sza_in ) then
            exit
         end if
      end do
      is   = max( min( is,numsza ) - 1,1 )
      isp1 = is + 1
      dels(1)  = max( 0._r8,min( 1._r8,(sza_in - sza(is)) * del_sza(is) ) )
      wrk0     = 1._r8 - dels(1)

      izl = 2
Level_loop : &
      do k = kbot,1,-1
!----------------------------------------------------------------------
!        ... find albedo indicies
!----------------------------------------------------------------------
         do ial = 1,numalb
	    if( alb(ial) > alb_in(k) ) then
	       exit
	    end if
	 end do
	 albind = max( min( ial,numalb ) - 1,1 )
!----------------------------------------------------------------------
!        ... find pressure level indicies
!----------------------------------------------------------------------
         if( p_in(k) > p(1) ) then
            pind  = 2
            wght1 = 1._r8
         else if( p_in(k) <= p(nump) ) then
            pind  = nump
            wght1 = 0._r8
         else
            do iz = izl,nump
               if( p(iz) < p_in(k) ) then
	          izl = iz
	          exit
               end if
            end do
            pind  = max( min( iz,nump ),2 )
            wght1 = max( 0._r8,min( 1._r8,(p_in(k) - p(pind)) * del_p(pind-1) ) )
         end if
!----------------------------------------------------------------------
!        ... find "o3 ratios"
!----------------------------------------------------------------------
         v3ratu  = colo3_in(k) / colo3(pind-1)
         do iv = 1,numcolo3
            if( o3rat(iv) > v3ratu ) then
               exit
            end if
         end do
         ratindu = max( min( iv,numcolo3 ) - 1,1 )

         if( colo3(pind) /= 0._r8 ) then
            v3ratl = colo3_in(k) / colo3(pind)
            do iv = 1,numcolo3
               if( o3rat(iv) > v3ratl ) then
                  exit
               end if
            end do
            ratindl = max( min( iv,numcolo3 ) - 1,1 )
	 else
            ratindl = ratindu
            v3ratl  = o3rat(ratindu)
	 end if

!----------------------------------------------------------------------
!        ... compute the weigths
!----------------------------------------------------------------------
	 ial   = albind
	 ialp1 = ial + 1
	 iv    = ratindl

         dels(2)  = max( 0._r8,min( 1._r8,(v3ratl - o3rat(iv)) * del_o3rat(iv) ) )
         dels(3)  = max( 0._r8,min( 1._r8,(alb_in(k) - alb(ial)) * del_alb(ial) ) )

	 wrk1         = (1._r8 - dels(2))*(1._r8 - dels(3))
	 wghtl(0,0,0) = wrk0*wrk1
	 wghtl(1,0,0) = dels(1)*wrk1
	 wrk1         = (1._r8 - dels(2))*dels(3)
	 wghtl(0,0,1) = wrk0*wrk1
	 wghtl(1,0,1) = dels(1)*wrk1
	 wrk1         = dels(2)*(1._r8 - dels(3))
	 wghtl(0,1,0) = wrk0*wrk1
	 wghtl(1,1,0) = dels(1)*wrk1
	 wrk1         = dels(2)*dels(3)
	 wghtl(0,1,1) = wrk0*wrk1
	 wghtl(1,1,1) = dels(1)*wrk1

	 iv  = ratindu
         dels(2)  = max( 0._r8,min( 1._r8,(v3ratu - o3rat(iv)) * del_o3rat(iv) ) )

	 wrk1         = (1._r8 - dels(2))*(1._r8 - dels(3))
	 wghtu(0,0,0) = wrk0*wrk1
	 wghtu(1,0,0) = dels(1)*wrk1
	 wrk1         = (1._r8 - dels(2))*dels(3)
	 wghtu(0,0,1) = wrk0*wrk1
	 wghtu(1,0,1) = dels(1)*wrk1
	 wrk1         = dels(2)*(1._r8 - dels(3))
	 wghtu(0,1,0) = wrk0*wrk1
	 wghtu(1,1,0) = dels(1)*wrk1
	 wrk1         = dels(2)*dels(3)
	 wghtu(0,1,1) = wrk0*wrk1
	 wghtu(1,1,1) = dels(1)*wrk1

	 iz   = pind
	 iv   = ratindl
	 ivp1 = iv + 1
         do wn = 1,nw
            psum_l(wn) = wghtl(0,0,0) * rsf_tab(wn,iz,is,iv,ial) &
                         + wghtl(0,0,1) * rsf_tab(wn,iz,is,iv,ialp1) &
                         + wghtl(0,1,0) * rsf_tab(wn,iz,is,ivp1,ial) &
                         + wghtl(0,1,1) * rsf_tab(wn,iz,is,ivp1,ialp1) &
                         + wghtl(1,0,0) * rsf_tab(wn,iz,isp1,iv,ial) &
                         + wghtl(1,0,1) * rsf_tab(wn,iz,isp1,iv,ialp1) &
                         + wghtl(1,1,0) * rsf_tab(wn,iz,isp1,ivp1,ial) &
                         + wghtl(1,1,1) * rsf_tab(wn,iz,isp1,ivp1,ialp1)
         end do

	 iz   = iz - 1
	 iv   = ratindu
	 ivp1 = iv + 1
         do wn = 1,nw
            psum_u = wghtu(0,0,0) * rsf_tab(wn,iz,is,iv,ial) &
                     + wghtu(0,0,1) * rsf_tab(wn,iz,is,iv,ialp1) &
                     + wghtu(0,1,0) * rsf_tab(wn,iz,is,ivp1,ial) &
                     + wghtu(0,1,1) * rsf_tab(wn,iz,is,ivp1,ialp1) &
                     + wghtu(1,0,0) * rsf_tab(wn,iz,isp1,iv,ial) &
                     + wghtu(1,0,1) * rsf_tab(wn,iz,isp1,iv,ialp1) &
                     + wghtu(1,1,0) * rsf_tab(wn,iz,isp1,ivp1,ial) &
                     + wghtu(1,1,1) * rsf_tab(wn,iz,isp1,ivp1,ialp1)
            rsf(wn,k) = (psum_l(wn) + wght1*(psum_u - psum_l(wn)))
         end do
!------------------------------------------------------------------------------
!      etfphot comes in as photons/cm^2/sec/nm  (rsf includes the wlintv factor -- nm)
!     ... --> convert to photons/cm^2/s 
!------------------------------------------------------------------------------
         rsf(:,k) = etfphot(:) * rsf(:,k)

      end do Level_loop

      deallocate( psum_l )

      end subroutine interpolate_rsf


      end module mo_jlong
