
      module mo_seto2

      use shr_kind_mod, only : r8 => shr_kind_r8
      use abortutils,   only : endrun
      use spmd_utils,   only : masterproc
      use cam_logfile,  only : iulog

      implicit none

      private
      public :: o2_xsect_inti
      public :: set_o2_xsect

      save

      integer :: nsrc
      integer :: ngast
      integer :: nla
      integer :: nwint
      real(r8), allocatable    :: wlint(:)
      real(r8), allocatable    :: xso2int(:)
      real(r8), allocatable    :: wlla(:)
      real(r8), allocatable    :: wlgast(:)

      contains

      subroutine o2_xsect_inti( o2_xsect_file )
!-----------------------------------------------------------------------------
!   purpose:
!   compute equivalent optical depths for o2 absorption, parameterized in
!   the sr bands and the lyman-alpha line.
!-----------------------------------------------------------------------------
!   parameters:
!   nz      - integer, number of specified altitude levels in the working (i)
!             grid
!   z       - real(r8), specified altitude working grid (km)                  (i)
!   nw      - integer, number of specified intervals + 1 in working       (i)
!             wavelength grid
!   wl      - real(r8), vector of lower limits of wavelength intervals in     (i)
!             working wavelength grid
!   cz      - real(r8), number of air molecules per cm^2 at each specified    (i)
!             altitude layer
!   zen     - real(r8), solar zenith angle                                    (i)
!   dto2    - real(r8), optical depth due to o2 absorption at each specified  (o)
!             vertical layer at each specified wavelength
!   xso2    - real(r8), molecular absorption cross section in sr bands at     (o)
!             each specified altitude and wavelength.  includes herzberg
!             continuum.
!-----------------------------------------------------------------------------

      use mo_params,     only : deltax
      use mo_inter,      only : inter2
      use mo_inter,      only : inter_inti
      use mo_wavelen,    only : nw, wl
      use ioFileMod,     only : getfil
      use cam_pio_utils, only : cam_pio_openfile
      use pio, only : file_desc_t, pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
           pio_get_var, pio_closefile, pio_nowrite

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      character(len=*), intent(in) :: o2_xsect_file

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      type(file_desc_t) :: ncid
      integer :: dimid
      integer :: vid
      integer :: astat
      integer :: iret
      integer :: i, wn, n
      integer :: wrk_ind(4)
      real(r8), allocatable :: x1(:)
      real(r8), allocatable :: y1(:)
      character(len=256) :: filespec
      character(len=256) :: locfn
      integer :: ierr
!-----------------------------------------------------------------------------
! 	... cross section data for use outside the sr-bands (combined from
!           brasseur and solomon and the jpl 1994 recommendation)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! 	... read o2 cross section data outside sr-bands
!-----------------------------------------------------------------------------
!	... o2 absorption cross sections:
!           from 116 nm to 245 nm, including schumann-runge continumm
!           from brasseur and solomon 1986.
!-----------------------------------------------------------------------------
      filespec = trim( o2_xsect_file )
      call getfil( filespec, locfn, 0 )
      call cam_pio_openfile( ncid, trim( locfn ), PIO_NOWRITE )
!---------------------------------------------------------------------------
! 	... get the dimensions
!---------------------------------------------------------------------------
      ierr = pio_inq_dimid( ncid, 'nosr', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, nsrc )
      ierr = pio_inq_dimid( ncid, 'ngast', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, ngast )
      ierr = pio_inq_dimid( ncid, 'nla', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, nla )
!---------------------------------------------------------------------------
! 	... allocate arrays
!---------------------------------------------------------------------------
      allocate( wlint(nsrc), xso2int(nsrc), x1(nsrc), y1(nsrc), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'o2_xsect_inti: failed to allocate wlint ... y1; error = ',astat
         call endrun
      end if
      allocate( wlgast(ngast), wlla(nla), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'o2_xsect_inti: failed to allocate wlgast, wlla; error = ',astat
         call endrun
      end if
!---------------------------------------------------------------------------
! 	... read the wave bin coordinates
!---------------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'wl_src', vid )
      ierr = pio_get_var( ncid, vid, x1 )
      ierr = pio_inq_varid( ncid, 'xs_src', vid )
      ierr = pio_get_var( ncid, vid, y1 )
      ierr = pio_inq_varid( ncid, 'wl_gast', vid )
      ierr = pio_get_var( ncid, vid, wlgast )
      ierr = pio_inq_varid( ncid, 'wl_lym', vid )
      ierr = pio_get_var( ncid, vid, wlla )
      call pio_closefile(  ncid )
!-----------------------------------------------------------------------------
! 	... put together the internal grid by "pasting" the lyman-alpha grid and 
!           kockarts grid into the combination of brasseur/solomon and jpl grid
!-----------------------------------------------------------------------------
      wlint(1:9) = x1(1:9)
      nwint = 9
      wlint(nwint+1:nwint+2) = wlla(1:2)
      nwint = 11
      wlint(nwint+1:nwint+36) = x1(12:47)
      nwint = 47
      wlint(nwint+1:nwint+ngast) = wlgast(1:ngast)
      nwint = nwint + ngast
      wlint(nwint+1:nwint+41) = x1(65:105)
      nwint = nwint + 41
      wrk_ind(1:4) = (/ nsrc, ngast, nla, nwint /)
!-----------------------------------------------------------------------------
! 	... initialize interpolation module
!-----------------------------------------------------------------------------
      call inter_inti( nw+1, wl, nsrc, wlint )
!-----------------------------------------------------------------------------
! 	... interpolate brasseur/solomon and jpl data onto internal grid
!-----------------------------------------------------------------------------
      call inter2( nsrc, wlint, xso2int, nsrc, x1, y1, iret )
      deallocate( x1, y1 )


      end subroutine o2_xsect_inti

      subroutine set_o2_xsect( z, nw, wl, cz, &
                               vcol, scol, dto2, xso2 )
!-----------------------------------------------------------------------------
!   purpose:
!   compute equivalent optical depths for o2 absorption, parameterized in
!   the sr bands and the lyman-alpha line.
!-----------------------------------------------------------------------------
!   parameters:
!   nz      - integer, number of specified altitude levels in the working (i)
!             grid
!   z       - real(r8), specified altitude working grid (km)
!   nw      - integer, number of specified intervals + 1 in working
!             wavelength grid
!   wl      - real(r8), vector of lower limits of wavelength intervals in
!             working wavelength grid
!   cz      - real(r8), number of air molecules per cm^2 at each specified
!             altitude layer
!   zen     - real(r8), solar zenith angle
!   dto2    - real(r8), optical depth due to o2 absorption at each specified  (o)
!             vertical layer at each specified wavelength
!   xso2    - real(r8), molecular absorption cross section in sr bands at     (o)
!             each specified altitude and wavelength.  includes herzberg
!             continuum.
!-----------------------------------------------------------------------------
!   edit history:
!   02/98  included lyman-alpha parameterization
!   03/97  fix dto2 problem at top level (nz)
!   02/97  changed offset for grid-end interpolation to relative number
!          (x * (1 +- deltax))
!   08/96  modified for early exit, no redundant read of data and smaller
!          internal grid if possible;  internal grid uses user grid points
!          whenever possible
!   07/96  modified to work on internal grid and interpolate final values
!          onto the user-defined grid
!-----------------------------------------------------------------------------

      use mo_params,  only : kw
      use mo_wavelen, only : delw_bin
      use mo_inter,   only : inter3
      use mo_schu,    only : schu
      use mo_lymana,  only : lymana
      use ppgrid,     only : pver, pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)    :: nw
      real(r8), intent(in)   :: wl(kw)
      real(r8), intent(in)   :: cz(pverp)
      real(r8), intent(in)   :: z(pverp)
      real(r8), intent(in)   :: vcol(pverp)
      real(r8), intent(in)   :: scol(pverp)
      real(r8), intent(out)  :: dto2(pver,nw)
      real(r8), intent(out)  :: xso2(nw,pverp)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: wn, k, igast
      integer  :: astat
      real(r8) :: secchi(pverp)

!-----------------------------------------------------------------------------
! 	... o2 optical depth and equivalent cross section on kockarts grid
!-----------------------------------------------------------------------------
      real(r8), allocatable :: dto2k(:,:)
      real(r8), allocatable :: xso2k(:,:)
!-----------------------------------------------------------------------------
! 	... o2 optical depth and equivalent cross section in the lyman-alpha region
!-----------------------------------------------------------------------------
      real(r8), allocatable :: dto2la(:,:)
      real(r8), allocatable :: xso2la(:,:)
!-----------------------------------------------------------------------------
! 	... temporary one-dimensional storage for optical depth and cross section values
!           xxtmp  - on internal grid
!           xxuser - on user defined grid
!-----------------------------------------------------------------------------
      real(r8), dimension(2*kw) :: dttmp, xstmp
      real(r8) :: dtuser(kw)
      real(r8) :: xsuser(kw)
      real(r8) :: o2col(pverp)

      real(r8) :: x, y
      real(r8) :: delo2

!-----------------------------------------------------------------------------
!	... allocate local variables
!-----------------------------------------------------------------------------
      allocate( dto2k(pver,ngast-1), xso2k(pverp,ngast-1), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'set_o2_xsect: failed to allocate dto2k,xso2k; error = ',astat
         call endrun
      end if
      allocate( dto2la(pver,nla-1), xso2la(pverp,nla-1), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'set_o2_xsect: failed to allocate dto2k,xso2k; error = ',astat
         call endrun
      end if
!-----------------------------------------------------------------------------
! 	... check, whether user grid is in the o2 absorption band at all...
!           if not, set cross section and optical depth values to zero and return
!-----------------------------------------------------------------------------
      dto2(:pver,:nw)  = 0._r8
      xso2(:nw,:pverp) = 0._r8
      if( wl(1) > 243._r8 ) then
         return
      end if

!-----------------------------------------------------------------------------
! 	... sec xhi or chapman calculation
!           for zen > 95 degrees, use zen = 95.  (this is only to compute effective o2
!           cross sections. still, better than setting dto2 = 0. as was done up to 
!           version 4.0) sm 1/2000
!           in future could replace with mu2(iz) (but mu2 is also wavelength-depenedent)
!           or imporved chapman function 
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! 	... slant o2 column 
!-----------------------------------------------------------------------------
      o2col(1:pverp) = 0.2095_r8 * scol(1:pverp)

!-----------------------------------------------------------------------------
! 	... effective secant of solar zenith angle.  use 2.0 if no direct sun. 
!           for nz, use value at nz-1
!-----------------------------------------------------------------------------
      secchi(1:pver) = scol(1:pver)/vcol(1:pver)
      where( secchi(1:pver) == 0._r8 )
         secchi(1:pver) = 2._r8
      endwhere
      secchi(pverp) = secchi(pver)

!-----------------------------------------------------------------------------
! 	... if necessary:
!           kockarts parameterization of the sr bands, output values of o2
!           optical depth and o2 equivalent cross section are on his grid
!-----------------------------------------------------------------------------
      if( wl(1) < wlgast(ngast) .and. wl(nw+1) > wlgast(1) ) then
           call schu( o2col, secchi, dto2k, xso2k )
      else
         dto2k(:,:) = 0._r8
         xso2k(:,:) = 0._r8
      end if

!-----------------------------------------------------------------------------
! 	... lyman-alpha parameterization, output values of o2 opticaldepth
!           and o2 effective (equivalent) cross section
!-----------------------------------------------------------------------------
      if( wl(1) <= wlla(nla) .and. wl(nw+1) >= wlla(1) ) then
         call lymana( o2col, secchi, dto2la, xso2la )
      else
         dto2la(:,:) = 0._r8
         xso2la(:,:) = 0._r8
      end if

!-----------------------------------------------------------------------------
! 	... loop through the altitude levels
!-----------------------------------------------------------------------------
level_loop : &
      do k = 1,pverp
         igast = 0
!-----------------------------------------------------------------------------
! 	... loop through the internal wavelength grid
!-----------------------------------------------------------------------------
         do wn = 1,nwint-1
!-----------------------------------------------------------------------------
! 	... if outside kockarts grid and outside lyman-alpha, use the 
!           jpl/brasseur+solomon data, if inside
!           kockarts grid, use the parameterized values from the call to schu,
!           if inside lyman-alpha, use the paraemterized values from call to lymana
!-----------------------------------------------------------------------------
            if( wlint(wn+1) <= wlgast(1) .or. wlint(wn) >= wlgast(ngast) ) then
              if( wlint(wn+1) <= wlla(1) .or. wlint(wn) >= wlla(nla) ) then
                 xstmp(wn) = xso2int(wn)
              else
                 xstmp(wn) = xso2la(k,1)
              end if
            else
               igast = igast + 1
               xstmp(wn) = xso2k(k,igast)
            end if
!-----------------------------------------------------------------------------
! 	... compute the area in each bin (for correct interpolation purposes only!)
!-----------------------------------------------------------------------------
            xstmp(wn) = xstmp(wn) * (wlint(wn+1) - wlint(wn))
         end do
!-----------------------------------------------------------------------------
! 	... interpolate o2 cross section from the internal grid onto the user grid
!-----------------------------------------------------------------------------
         call inter3( nw+1, wl, xsuser, nwint, wlint, xstmp )
         xso2(:nw,k) = xsuser(:nw) * delw_bin(:nw)
      end do level_loop

      do k = 1,pver
         igast = 0
         delo2 = .2095_r8 * cz(k)    ! vertical o2 column
!-----------------------------------------------------------------------------
! 	... loop through the internal wavelength grid
!-----------------------------------------------------------------------------
         do wn = 1,nwint-1
!-----------------------------------------------------------------------------
! 	... if outside kockarts grid and outside lyman-alpha, use the 
!           jpl/brasseur+solomon data, if inside
!           kockarts grid, use the parameterized values from the call to schu,
!           if inside lyman-alpha, use the paraemterized values from call to lymana
!-----------------------------------------------------------------------------
            if( wlint(wn+1) <= wlgast(1) .or. wlint(wn) >= wlgast(ngast) ) then
              if( wlint(wn+1) <= wlla(1) .or. wlint(wn) >= wlla(nla) ) then
                 dttmp(wn) = xso2int(wn) * delo2
              else
                 dttmp(wn) = dto2la(k,1)
              end if
            else
               igast = igast + 1
               dttmp(wn) = dto2k(k,igast)
            end if
!-----------------------------------------------------------------------------
! 	... compute the area in each bin (for correct interpolation purposes only!)
!-----------------------------------------------------------------------------
            dttmp(wn) = dttmp(wn) * (wlint(wn+1) - wlint(wn))
         end do
!-----------------------------------------------------------------------------
! 	... interpolate o2 optical depth from the internal grid onto the user grid
!-----------------------------------------------------------------------------
         call inter3( nw+1, wl, dtuser, nwint, wlint, dttmp )
         dto2(k,:nw) = dtuser(:nw) * delw_bin(:nw)
      end do

      deallocate( dto2k, xso2k, dto2la, xso2la )

      end subroutine set_o2_xsect

      end module mo_seto2
