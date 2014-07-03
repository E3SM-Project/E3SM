module mo_apex

!-------------------------------------------------------------------------------
! Purpose:
!
!   Calculate apex coordinates and magnetic field magnitudes
!   at global geographic grid for year of current model run. 
!
! Method: 
!
!   The magnetic field parameters output by this module are time and height
!     independent. They are chunked for waccm physics, i.e., allocated as 
!     (pcols,begchunk:endchunk)
!   Interface sub apexmag is called once per run from sub inti.
!     Sub apexmag may be called for years 1900 through 2005.
!   This module is dependent on routines in apex_subs.F (modified IGRF model).
!   Apex_subs has several authors, but has been modified and maintained
!     in recent years by Roy Barnes (bozo@ucar.edu).
!   Subs apxmka and apxmall are called with the current lat x lon grid 
!     resolution.
!
! Author: Ben Foster, foster@ucar.edu (Nov, 2003)
!-------------------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, begchunk, endchunk          ! physics grid
   use phys_grid,       only: get_lat_p, get_lon_p, get_ncols_p
   use cam_history,     only: addfld, phys_decomp, add_default   ! for history saves
   use abortutils,      only: endrun
   use cam_control_mod, only: magfield_fix_year
   use cam_logfile,     only: iulog
   use spmd_utils,      only: masterproc

   implicit none

   private
   public :: apexmag
   public :: alatm, alonm, glonm, bnorth, beast, bdown, bmag
   public :: d1vec, d2vec, colatp, elonp

!-------------------------------------------------------------------------------
! Magnetic field output arrays, chunked for physics:
! (these are allocated (pcols,begchunk:endchunk) by sub apex_alloc)
!-------------------------------------------------------------------------------
   real(r8), allocatable, dimension(:,:) :: & ! (pcols,begchunk:endchunk)
     alatm,  & ! apex mag latitude at each geographic grid point (radians)
     alonm,  & ! apex mag longitude at each geographic grid point (radians)
     glonm,  & ! apex mag longitude at each geographic grid point (radians)
     bnorth, & ! northward component of magnetic field
     beast,  & ! eastward component of magnetic field
     bdown,  & ! downward component of magnetic field
     bmag      ! magnitude of magnetic field
   real(r8), allocatable, dimension(:,:,:) :: & ! (3,pcols,begchunk:endchunk)
     d1vec,    & ! base vectors more-or-less magnetic eastward direction
     d2vec       ! base vectors more-or-less magnetic downward/equatorward direction
   real(r8) :: &
     colatp,   & ! geocentric colatitude of geomagnetic dipole north pole (deg)
     elonp	 ! East longitude of geomagnetic dipole north pole (deg)

contains

subroutine apexmag
!-------------------------------------------------------------------------------
! Driver for apex code to calculate apex magnetic coordinates at 
!   current geographic spatial resolution for given year. This calls
!   routines in apex_subs.F.
!
! This is called once per run from sub inti.
!-------------------------------------------------------------------------------

   use physconst,      only: pi
   use dyn_grid, only : get_dyn_grid_parm, get_horiz_grid_d
!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------
   integer, parameter  :: nalt  = 2             ! number of altitudes
   real(r8), parameter :: re    = 6.378165e8_r8 ! earth radius (cm)
   real(r8), parameter :: h0    = 9.0e6_r8      ! base height (90 km)
   real(r8), parameter :: hs    = 1.3e7_r8
   real(r8), parameter :: eps   = 1.e-6_r8      ! epsilon
   real(r8), parameter :: cm2km = 1.e-5_r8

   integer  :: i, j, ist, lchnk              ! indices
   integer  :: ncol, ic, jc

   real(r8), allocatable :: dlat_mid(:)              ! midpoint latitudes to assist interp
   real(r8), allocatable :: dlon_mid(:)              ! midpoint longitudes to assist interp
   real(r8), allocatable :: dlon_shft(:)             ! shifted midpoint longitudes

   real(r8) :: rtd, dtr                      ! deg<->radians conversions
   real(r8) :: rekm, h0km, alt, hr, ror03, alon, alat, & ! apxmall args
               vmp, w, d, be3, sim, xlatqd, f, si, date1(1), collat, collon
   real(r8) :: gpalt(nalt)                   ! altitudes passed to apxmka
   real(r8), allocatable :: wk(:), clat(:), clon(:)
   integer :: plon, plat, plonp1, platp1, lwk
!-------------------------------------------------------------------------------
! Non-scalar arguments returned by APXMALL:
!-------------------------------------------------------------------------------
   real(r8) :: bhat(3)
   real(r8) :: d3(3)
   real(r8) :: e1(3), e2(3), e3(3)
   real(r8) :: f1(2), f2(2)

   real(r8), allocatable :: glatm(:,:)
   real(r8), allocatable :: bglobal(:,:,:)
   real(r8), allocatable :: d1global(:,:,:)
   real(r8), allocatable :: d2global(:,:,:)
   real(r8), allocatable :: bmglobal(:,:)   

   real(r8) :: rdum

   plat = get_dyn_grid_parm('plat')
   plon = get_dyn_grid_parm('plon')
   platp1 = plat+1
   plonp1 = plon+1
   lwk = platp1*plonp1*nalt*5+platp1 + plonp1 + nalt

   allocate(dlat_mid(platp1), dlon_mid(plonp1), dlon_shft(plonp1), wk(lwk))
   allocate(glatm(plon,plat), bglobal(3,plon,plat), d1global(3,plon,plat), d2global(3,plon,plat))
   allocate(bmglobal(plon,plat))
   allocate(clat(plat), clon(plon))
!-------------------------------------------------------------------------------
! Allocate output arrays
!-------------------------------------------------------------------------------
   call apex_alloc

   rtd = 180._r8/pi     ! radians to degrees
   dtr = pi/180._r8     ! degrees to radians

!-------------------------------------------------------------------------------
! Extend staggered latitudes (midpoints) to assist apxmka build 
! interpolation tables (degrees).
!-------------------------------------------------------------------------------
   call get_horiz_grid_d(plat, clat_d_out=clat)
   call get_horiz_grid_d(plon, clon_d_out=clon)

   dlat_mid(1) = clat(1)*rtd
   do j = 2,plat
      dlat_mid(j) = rtd*(clat(j-1) + clat(j))/2.0_r8
   end do
   dlat_mid(platp1) = rtd*clat(plat)

!-------------------------------------------------------------------------------
! Set longitude midpoints (assume non-reducing grid, so use londeg(:,1))
!-------------------------------------------------------------------------------
   dlon_mid(1) = rtd*(clon(1) - .5_r8*(clon(2) - clon(1)))
   do i = 2,plon
     dlon_mid(i) = .5_r8*rtd*(clon(i) + clon(i-1))
   end do
   dlon_mid(plon+1) = dlon_mid(plon) + rtd*(clon(plon) - clon(plon-1))


#ifdef AURORA_DEBUG
   write(iulog,*) '---------------------------------------'
   write(iulog,*) 'apex: dlon_mid'
   write(iulog,'(1p,5g15.7)') dlon_mid(:)
   write(iulog,*) 'apex: dlat_mid'
   write(iulog,'(1p,5g15.7)') dlat_mid(:)
   write(iulog,*) '---------------------------------------'
#endif

   date1(1) = magfield_fix_year

!-------------------------------------------------------------------------------
! Center min, max altitudes about 130 km
!-------------------------------------------------------------------------------
   gpalt(1) =  90._r8
   gpalt(2) = 170._r8

!-------------------------------------------------------------------------------
! Set up interpolation.
! (Note apxmka expects longitudes in -180 -> +180)
!-------------------------------------------------------------------------------
   dlon_shft(:) = dlon_mid(:) - 180._r8
   call apxmka( iulog, date1, dlat_mid, dlon_shft, gpalt, &
                platp1, plonp1, nalt, wk, lwk, ist )
   if( ist /= 0 ) then
     write(iulog,"(/,'>>> apexmag: Error from apxmka: ist=',i5)") ist
     call endrun
   end if

   rekm  = re*cm2km    ! earth radius (km)
   h0km  = h0*cm2km    ! base height (km)
   alt   = hs*cm2km    ! altitude for apxmall (km)
   hr    = alt         ! reference altitude (km)
   ror03 = ((rekm + alt)/(rekm + h0km))**3

!------------------------------------------------------------------------------
! Apex coords alon, alat are returned for each geographic grid point:
! first form global arrays
!------------------------------------------------------------------------------
   do j = 1,plat
      collat = clat(j)*rtd                       ! latitude of current column (deg)
      do i = 1,plon
         collon = clon(i)*rtd ! longitude of current column (deg)
         call apxmall(                                  &
           collat, collon, alt, hr, wk, lwk,            & ! Inputs
           bglobal(1,i,j), bhat, bmglobal(i,j), si,     & ! Mag Fld
           alon, alat,                                  & ! Apex lon,lat output
           vmp, w, d, be3, sim, d1global(1,i,j), d2global(1,i,j), d3, e1, e2, e3, & ! Mod Apex
           xlatqd, f, f1, f2, ist )                       ! Qsi-Dpl
         if( ist /= 0 ) then
           write(iulog,"(/,'>>> apexmag: Error from apxmall: ist=',i4)") ist
           call endrun
         end if
         glonm(i,j)   =  alon*dtr                         ! mag lons (radians)
         glatm(i,j)   =  alat*dtr                         ! mag lats (radians)
      end do
   end do

#ifdef AURORA_DEBUG
   write(iulog,*) '---------------------------------------'
   write(iulog,*) 'apex: geo lons'
   write(iulog,'(1p,5g15.7)') rtd*clon(:)
   write(iulog,*) 'apex: geo lats'
   write(iulog,'(1p,5g15.7)') rtd*clat(:)
   write(iulog,*) ' '
   do j = 1,plat
      write(iulog,*) 'apexmag: mag lats @ lat ndx = ',j
      write(iulog,'(1p,5g15.7)') glatm(:,j)/dtr
   end do
   write(iulog,*) ' '
   do j = 1,plat
      write(iulog,*) 'apexmag: mag lons @ lat ndx = ',j
      write(iulog,'(1p,5g15.7)') glonm(:,j)/dtr
   end do
   write(iulog,*) ' '
   write(iulog,*) '---------------------------------------'
#endif

!------------------------------------------------------------------------------
! map globals to chunks
!------------------------------------------------------------------------------
chunk_loop : &
   do lchnk = begchunk,endchunk  ! physics chunks (groups of arbitrary columns)
     ncol = get_ncols_p( lchnk )
column_loop : &
     do i = 1,ncol             ! number of columns in each chunk
!------------------------------------------------------------------------------
! set module output from apxmall output:
!------------------------------------------------------------------------------
       ic                = get_lon_p(lchnk,i)
       jc                = get_lat_p(lchnk,i)
!      write(iulog,*) 'apex: lchnk,i,ic,jc = ',lchnk,i,ic,jc
       alatm (i,lchnk)   =  glatm(ic,jc)                ! mag lats (radians)
       alonm (i,lchnk)   =  glonm(ic,jc)                ! mag lons (radians)
!      bnorth(i,lchnk)   =  bglobal(2,ic,jc)*1.e-5_r8   ! northward nT -> gauss
!      beast (i,lchnk)   =  bglobal(1,ic,jc)*1.e-5_r8   ! eastward  nT -> gauss
!      bdown (i,lchnk)   = -bglobal(3,ic,jc)*1.e-5_r8   ! downward  nT -> gauss
!      bmag  (i,lchnk)   =  bmglobal(ic,jc)*1.e-5_r8    ! magnitude nT -> gauss
       bnorth(i,lchnk)   =  bglobal(2,ic,jc)            ! northward nT -> gauss
       beast (i,lchnk)   =  bglobal(1,ic,jc)            ! eastward  nT -> gauss
       bdown (i,lchnk)   = -bglobal(3,ic,jc)            ! downward  nT -> gauss
       bmag  (i,lchnk)   =  bmglobal(ic,jc)             ! magnitude nT -> gauss
       d1vec (:,i,lchnk) =  d1global(:,ic,jc)
       d2vec (:,i,lchnk) =  d2global(:,ic,jc)
     end do column_loop

#ifdef DEBUG
     write(iulog,"('apexmag: lchnk=',i3,' alatm(:,lchnk)*rtd=',/,(1p,6g12.4))") &
       lchnk,alatm(:,lchnk)*rtd
     write(iulog,"('apexmag: lchnk=',i3,' alonm(:,lchnk)*rtd=',/,(1p,6g12.4))") &
       lchnk,alonm(:,lchnk)*rtd
     write(iulog,"('apexmag: lchnk=',i3,' bnorth(:,lchnk)=',/,(6e12.4))") &
       lchnk,bnorth(:,lchnk)
     write(iulog,"('apexmag: lchnk=',i3,' beast(:,lchnk)=' ,/,(6e12.4))") &
       lchnk,beast(:,lchnk)
     write(iulog,"('apexmag: lchnk=',i3,' bdown(:,lchnk)=' ,/,(6e12.4))") &
       lchnk,bdown(:,lchnk)
     write(iulog,"('apexmag: lchnk=',i3,' bmag(:,lchnk)='  ,/,(6e12.4))") &
       lchnk,bmag(:,lchnk)
#endif
   end do chunk_loop

#ifdef DEBUG
   write(iulog,*) '---------------------------------------'
   write(iulog,*) 'apex: begchunk,endchunk = ', begchunk,endchunk
   write(iulog,*) 'apex: mag lons at j = 81'
   write(iulog,'(1p,5g15.7)') glonm(:,81)/dtr
   write(iulog,*) 'apex: mag lons at j = 82'
   write(iulog,'(1p,5g15.7)') glonm(:,82)/dtr
   write(iulog,*) 'apex: mag lons at j = 83'
   write(iulog,'(1p,5g15.7)') glonm(:,83)/dtr
   write(iulog,*) '---------------------------------------'
#endif

   call dypol( colatp, elonp, rdum )	! get geomagnetic dipole north pole 

   if (masterproc) write(iulog,*) 'apex: colatp,elonp ', colatp, elonp

!------------------------------------------------------------------------------
! Add mag field output to master field list:
!------------------------------------------------------------------------------
  call addfld('ALATM   ','RADIANS ',1,'I',&
    'Magnetic latitude at each geographic coordinate',phys_decomp)
  call addfld('ALONM   ','RADIANS ',1,'I',&
    'Magnetic longitude at each geographic coordinate',phys_decomp)
! call addfld('ALATM   ','RADIANS ',1,'A',&
!   'Magnetic latitude at each geographic coordinate',phys_decomp)
! call addfld('ALONM   ','RADIANS ',1,'A',&
!   'Magnetic longitude at each geographic coordinate',phys_decomp)
! call addfld('BNORTH  ','GAUSS',1,'I',&
!   'Northward component of magnetic field',phys_decomp)
! call addfld('BEAST   ','GAUSS',1,'I',&
!   'Eastward component of magnetic field',phys_decomp)
! call addfld('BDOWN   ','GAUSS',1,'I',&
!   'Downward component of magnetic field',phys_decomp)
! call addfld('BMAG    ','GAUSS',1,'I',&
!   'Magnetic field magnitude',phys_decomp)

!------------------------------------------------------------------------------
! Write these fields to history by default:
!------------------------------------------------------------------------------
! call add_default ('ALATM' , 1, ' ')
! call add_default ('ALONM' , 1, ' ')
! call add_default ('ALATM' , 2, ' ')
! call add_default ('ALONM' , 2, ' ')
! call add_default ('BNORTH', 1, ' ')
! call add_default ('BEAST' , 1, ' ')
! call add_default ('BDOWN' , 1, ' ')
! call add_default ('BMAG'  , 1, ' ')

  deallocate(dlat_mid, dlon_mid, dlon_shft, wk)
  deallocate(glatm, bglobal, d1global, d2global)
  deallocate(bmglobal)
  deallocate(clat, clon)
  if (masterproc) write(iulog,"(' APEXMAG: Calculated apex magnetic coordinates for year AD ',i4)") &
                  magfield_fix_year

end subroutine apexmag

subroutine apex_alloc
  use dyn_grid, only : get_dyn_grid_parm
!------------------------------------------------------------------------------
! Allocate module output arrays for chunked physics grid.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! local variables
!------------------------------------------------------------------------------
  integer :: istat  ! status of allocate statements
  integer :: plat, plon

  plat = get_dyn_grid_parm('plat')
  plon = get_dyn_grid_parm('plon')

  if (.not.allocated(alatm)) then
    allocate(alatm(pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of alatm failed: istat=',i5)") istat
      call endrun
    end if
!   write(iulog,"('apex_alloc: allocated alatm: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  end if

  if (.not.allocated(glonm)) then
    allocate(glonm(plon,plat),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of glonm failed: istat=',i5)") istat
      call endrun
    end if
  end if

  if (.not.allocated(alonm)) then
    allocate(alonm(pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of alonm failed: istat=',i5)") istat
      call endrun
    end if
!   write(iulog,"('apex_alloc: allocated alonm: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  end if

  if (.not.allocated(bnorth)) then
    allocate(bnorth(pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of bnorth failed: istat=',i5)") istat
      call endrun
    end if
!   write(iulog,"('apex_alloc: allocated bnorth: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  end if

  if (.not.allocated(beast)) then
    allocate(beast(pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of beast failed: istat=',i5)") istat
      call endrun
    end if
!   write(iulog,"('apex_alloc: allocated beast: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  end if

  if (.not.allocated(bdown)) then
    allocate(bdown(pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of bdown failed: istat=',i5)") istat
      call endrun
    end if
!   write(iulog,"('apex_alloc: allocated bdown: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  end if

  if (.not.allocated(bmag)) then
    allocate(bmag(pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of bmag failed: istat=',i5)") istat
      call endrun
    end if
!   write(iulog,"('apex_alloc: allocated bmag: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  end if
  if (.not.allocated(d1vec)) then
    allocate(d1vec(3,pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of d1vec failed: istat=',i5)") istat
      call endrun
    endif
!   write(iulog,"('apex_alloc: allocated d1vec: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  endif

  if (.not.allocated(d2vec)) then
    allocate(d2vec(3,pcols,begchunk:endchunk),stat=istat)
    if (istat /= 0) then
      write(iulog,"('>>> apex_alloc: allocate of d2vec failed: istat=',i5)") istat
      call endrun
    endif
!   write(iulog,"('apex_alloc: allocated d2vec: pcols=',i3,' beg,endchunk=',2i4)") &
!     pcols,begchunk,endchunk
  endif

end subroutine apex_alloc

subroutine apex_dealloc
!-------------------------------------------------------------------------------
! Deallocate module output arrays:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------
  integer :: istat(8)

  if (allocated(alatm )) deallocate(alatm ,stat=istat(1))
  if (allocated(alonm )) deallocate(alatm ,stat=istat(2))
  if (allocated(bnorth)) deallocate(bnorth,stat=istat(3))
  if (allocated(beast )) deallocate(beast ,stat=istat(4))
  if (allocated(bdown )) deallocate(bdown ,stat=istat(5))
  if (allocated(bmag  )) deallocate(bmag  ,stat=istat(6))
  if (allocated(d1vec )) deallocate(d1vec ,stat=istat(7))
  if (allocated(d2vec )) deallocate(d2vec ,stat=istat(8))
  if (sum(istat) /= 0) then
    write(iulog,"('>>> apex_dealloc: error deallocating apex arrays: istat=',i4)") istat
    call endrun('apex_dealloc: error deallocating apex arrays')
  endif

end subroutine apex_dealloc

end module mo_apex
