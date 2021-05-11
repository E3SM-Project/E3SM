module mkurbanparCommonMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkurbanparCommon
!
! !DESCRIPTION:
! Common routines for making urban parameter data, independent of the method used for
! making the urban parameters (e.g., averages, dominant type, etc.)
!
! (WJS 4-18-12: In the past, this contained routines shared between mkurbanparDomMod and
! mkurbanparAvgMod; now there is just a single module, mkurbanparMod, but I am keeping the
! separate mkurbanparCommonMod in case a similar split comes back in the future. However,
! if such a split seems unlikely in the future, these routines could be moved back into
! mkurbanparMod.)
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!-----------------------------------------------------------------------
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
   use shr_sys_mod , only : shr_sys_flush
   implicit none
   
   private
   
! !PUBLIC MEMBER FUNCTIONS:
   public :: mkurban_pct             ! Make output urban %, given input urban %
   public :: mkurban_pct_diagnostics ! print diagnostics related to pct urban
   public :: mkelev                  ! Get elevation to reduce urban for high elevation areas
!
! !PUBLIC DATA MEMBERS:
!
   real(r8), parameter :: MIN_DENS = 0.1_r8       ! minimum urban density (% of grid cell) - below this value, urban % is set to 0

   public :: MIN_DENS
!
!EOP

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban_pct
!
! !INTERFACE:
subroutine mkurban_pct(ldomain, tdomain, tgridmap, urbn_i, urbn_o)
!
! !DESCRIPTION:
! make percent urban on output grid, given percent urban on input grid
!
! This assumes that we're neither using all_urban or zero_out
!
!
! !USES:
   use mkdomainMod , only : domain_type, domain_checksame
   use mkgridmapMod
   use mkvarctl    , only : mksrf_gridtype
!
! !ARGUMENTS:
   implicit none
   type(domain_type) , intent(in) :: ldomain
   type(domain_type) , intent(in) :: tdomain    ! local domain
   type(gridmap_type), intent(in) :: tgridmap   ! local gridmap
   real(r8)          , intent(in) :: urbn_i(:)  ! input grid: percent urban
   real(r8)          , intent(out):: urbn_o(:)  ! output grid: percent urban
!
! !REVISION HISTORY:
! Author: Bill Sacks
! (Moved from mkurbanparMod Feb, 2012)
!
!
! !LOCAL VARIABLES:
!EOP
   real(r8) :: sum_fldi                        ! global sum of dummy input fld
   real(r8) :: sum_fldo                        ! global sum of dummy output fld
   integer  :: ni,no                           ! indices
   real(r8) :: relerr = 0.00001_r8             ! max error: sum overlap wts ne 1
   character(len=*), parameter :: subname = 'mkurban_pct'
!-----------------------------------------------------------------------

   ! Error checks for array size consistencies

   if (size(urbn_i) /= tdomain%ns .or. &
       size(urbn_o) /= ldomain%ns) then
      write(6,*) subname//' ERROR: array size inconsistencies'
      write(6,*) 'size(urbn_i) = ', size(urbn_i)
      write(6,*) 'tdomain%ns   = ', tdomain%ns
      write(6,*) 'size(urbn_o) = ', size(urbn_o)
      write(6,*) 'ldomain%ns   = ', ldomain%ns
      stop
   end if
   
   ! Error checks for domain and map consistencies
   
   call domain_checksame( tdomain, ldomain, tgridmap )

   ! Determine urbn_o on ouput grid:
   ! Area-average percent cover on input grid to output grid 
   ! and correct according to land landmask
   ! Note that percent cover is in terms of total grid area.   

   call gridmap_areaave(tgridmap, urbn_i, urbn_o, nodata=0._r8)

   ! Check for conservation

   do no = 1, ldomain%ns
      if ((urbn_o(no)) > 100.000001_r8) then
         write (6,*) 'MKURBAN error: urban = ',urbn_o(no), &
              ' greater than 100.000001 for column, row = ',no
         stop
      end if
   enddo

   ! Global sum of output field -- must multiply by fraction of
   ! output grid that is land as determined by input grid

   sum_fldi = 0.0_r8
   do ni = 1,tdomain%ns
      sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
   enddo

   sum_fldo = 0._r8
   do no = 1, ldomain%ns
      sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
   end do

   ! -----------------------------------------------------------------
   ! Error check1
   ! Compare global sum fld_o to global sum fld_i.
   ! -----------------------------------------------------------------

   if (trim(mksrf_gridtype) == 'global') then
      if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
         write (6,*) 'MKURBAN error: input field not conserved'
         write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
         write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
         stop
      end if
   end if

   ! (Error check2 in mkurban_pct_diagnostics, which should be called separately)

end subroutine mkurban_pct
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban_pct_diagnostics
!
! !INTERFACE:
subroutine mkurban_pct_diagnostics(ldomain, tdomain, tgridmap, urbn_i, urbn_o, ndiag, dens_class)
!
! !DESCRIPTION:
! print diagnostics related to pct urban
!
! This is intended to be called after mkurban_pct, but is split out into a separate
! routine so that modifications to urbn_o can be made in between the two calls (e.g.,
! setting urbn_o to 0 wherever it is less than a certain threshold; the rules for doing
! this can't always be applied inline in mkurban_pct).
!
! !USES:
   use mkdomainMod , only : domain_type
   use mkgridmapMod, only : gridmap_type
   use mkvarpar
!
! !ARGUMENTS:
   implicit none
   type(domain_type) , intent(in) :: ldomain
   type(domain_type) , intent(in) :: tdomain    ! local domain
   type(gridmap_type), intent(in) :: tgridmap   ! local gridmap
   real(r8)          , intent(in) :: urbn_i(:)  ! input grid: percent urban
   real(r8)          , intent(in) :: urbn_o(:)  ! output grid: percent urban
   integer           , intent(in) :: ndiag      ! unit number for diag out

   integer , intent(in), optional :: dens_class ! density class
!
! !REVISION HISTORY:
! Author: Bill Sacks
! (Moved from mkurbanparMod Feb, 2012)
!
!
! !LOCAL VARIABLES:
!EOP
   real(r8) :: gurbn_i                         ! input  grid: global urbn
   real(r8) :: garea_i                         ! input  grid: global area
   real(r8) :: gurbn_o                         ! output grid: global urbn
   real(r8) :: garea_o                         ! output grid: global area
   integer  :: ni,no,k                         ! indices
!-----------------------------------------------------------------------

   ! -----------------------------------------------------------------
   ! Error check2
   ! Compare global areas on input and output grids
   ! -----------------------------------------------------------------

   ! Input grid

   gurbn_i = 0._r8
   garea_i = 0._r8

   do ni = 1, tdomain%ns
      garea_i = garea_i + tgridmap%area_src(ni)*re**2
      gurbn_i = gurbn_i + urbn_i(ni)*(tgridmap%area_src(ni)/100._r8)*&
           tgridmap%frac_src(ni)*re**2
   end do

   ! Output grid

   gurbn_o = 0._r8
   garea_o = 0._r8

   do no = 1, ldomain%ns
      garea_o = garea_o + tgridmap%area_dst(no)*re**2
      gurbn_o = gurbn_o + urbn_o(no)* (tgridmap%area_dst(no)/100._r8)*&
           tgridmap%frac_dst(no)*re**2
   end do

   ! Diagnostic output

   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('=',k=1,70)
   if (present(dens_class)) then
      write (ndiag,'(1x,a,i0)') 'Urban Output -- class ', dens_class
   else
      write (ndiag,'(1x,a)') 'Urban Output'
   end if
   write (ndiag,'(1x,70a1)') ('=',k=1,70)

   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
        1x,'                 10**6 km**2      10**6 km**2   ')
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,*)
   write (ndiag,2003) gurbn_i*1.e-06,gurbn_o*1.e-06
   write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'urban       ',f14.3,f17.3)
2003 format (1x,'urban       ',f14.3,f22.8)
2004 format (1x,'all surface ',f14.3,f17.3)

end subroutine mkurban_pct_diagnostics
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkelev
!
! !INTERFACE:
subroutine mkelev(ldomain, mapfname, datfname, varname, ndiag, elev_o)
!
! !DESCRIPTION:
! Make elevation data
!
! !USES:
  use mkdomainMod  , only : domain_type, domain_clean, domain_read, domain_checksame
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  character(len=*)  , intent(in) :: varname   ! topo variable name
  real(r8)          , intent(out):: elev_o(:) ! output elevation data
!
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Keith Oleson
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: elev_i(:)          ! canyon_height to width ratio in
  real(r8), allocatable :: mask_i(:)          ! input grid: mask (0, 1)
  integer  :: ns_i,ns_o                       ! indices
  integer  :: k,l,n,m,ni                      ! indices
  integer  :: ncidi,dimid,varid               ! input netCDF id's
  integer  :: ier                             ! error status
  character(len=256) :: name                  ! name of attribute
  character(len=256) :: unit                  ! units of attribute
  character(len= 32) :: subname = 'mkelev'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make elevation .....'
  call shr_sys_flush(6)

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)

  ns_i = tdomain%ns
  allocate(elev_i(ns_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkelev allocation error'; call abort()
  end if

  write (6,*) 'Open elevation file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)
  call check_ret(nf_inq_varid (ncidi, trim(varname), varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, elev_i), subname)
  call check_ret(nf_close(ncidi), subname)

  ! Read topo elev dataset with unit mask everywhere

  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies
  ! Note that the topo dataset has no landmask - so a unit landmask is assumed

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Determine elev_o on output grid

  elev_o(:) = 0.

  call gridmap_areaave(tgridmap, elev_i, elev_o, nodata=0._r8)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (elev_i)

  write (6,*) 'Successfully made elevation' 
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkelev

!-----------------------------------------------------------------------

end module mkurbanparCommonMod
