module mkpftMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkpft
!
! !DESCRIPTION:
! Make PFT data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkvarctl    , only : numpft
  use mkdomainMod , only : domain_checksame

  implicit none

  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mkpftInit         ! Initialization
  public mkpft             ! Set PFT
  public mkpft_parse_oride ! Parse the string with PFT fraction/index info to override
  public mkirrig           ! Set irrigation
  public mkpftAtt          ! Write out attributes to output file on pft
!
! !PUBLIC DATA MEMBERS: 
!

  !
  ! When pft_idx and pft_frc are set, they must be set together, and they will cause the
  ! entire area to be covered with vegetation and zero out other landunits.
  ! The sum of pft_frc must = 100%, and each pft_idx point in the array corresponds to
  ! the fraction in pft_frc. Only the first few points are used until pft_frc = 0.0.
  !
  integer            :: m                     ! index
  integer, parameter :: maxpft = 20           ! maximum # of PFT
  integer, public    :: pft_idx(0:maxpft) = & ! PFT vegetation index to override with
                             (/ ( -1,  m = 0, maxpft )   /)
  real(r8), public   :: pft_frc(0:maxpft) = & ! PFT vegetation fraction to override with
                             (/ ( 0.0, m = 0, maxpft ) /)
  integer, public :: nonIrrigIdx = 15
  integer, public :: IrrigIdx    = 16
!
! !PRIVATE DATA MEMBERS:
!
  logical, private :: zero_out      = .false. ! Flag to zero out PFT
  logical, private :: use_input_pft = .false. ! Flag to override PFT with input values
  integer, private :: nzero                   ! index of first zero fraction
!
! !PRIVATE MEMBER FUNCTIONS:
!
  private :: mkpft_check_oride  ! Check the pft_frc and pft_idx values for correctness
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpftInit
!
! !INTERFACE:
subroutine mkpftInit( zero_out_l, all_veg )
!
! !DESCRIPTION:
! Initialize of Make PFT data
! !USES:
!
! !ARGUMENTS:
  implicit none
  logical, intent(IN)  :: zero_out_l ! If veg should be zero'ed out
  logical, intent(OUT) :: all_veg    ! If should zero out other fractions so that
                                     ! all land-cover is vegetation
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8), parameter :: hndrd = 100.0_r8  ! A hundred percent
  character(len=32) :: subname = 'mkpftInit:: '
!-----------------------------------------------------------------------

  call mkpft_check_oride( )
  if ( use_input_pft ) then
     if ( maxpft < numpft ) then
        write(6,*) subname//'number PFT is > max allowed!'
        call abort()
     end if
     write(6,*) 'Set PFT fraction to : ', pft_frc(0:nzero-1)
     write(6,*) 'With PFT index      : ', pft_idx(0:nzero-1)
  end if

  all_veg = use_input_pft

  if ( zero_out_l .and. all_veg )then
     write(6,*) subname//'zeroing out vegetation and setting vegetation to 100% is a contradiction!'
     call abort()
  end if

  ! Copy local zero out to module data version
  zero_out = zero_out_l

end subroutine mkpftInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft
!
! !INTERFACE:
subroutine mkpft(ldomain, mapfname, fpft, firrig, ndiag, &
     pctlnd_o, pctirr_o, pctpft_o)
!
! !DESCRIPTION:
! Make PFT data
! This dataset consists of the %cover of the [numpft]+1 PFTs used by
! the model. The input %cover pertains to the "vegetated" portion of the
! grid cell and sums to 100. The real portion of each grid cell
! covered by each PFT is the PFT cover times the fraction of the
! grid cell that is land. This is the quantity preserved when
! area-averaging from the input (1/2 degree) grid to the models grid.
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(inout) :: ldomain
  character(len=*)  , intent(in) :: mapfname              ! input mapping file name
  character(len=*)  , intent(in) :: fpft                  ! input pft dataset file name
  character(len=*)  , intent(in) :: firrig                ! input irrigation dataset file name
  integer           , intent(in) :: ndiag                 ! unit number for diag out
  real(r8)          , intent(in) :: pctirr_o(:)           ! % irrigated area (output grid)
  real(r8)          , intent(out):: pctlnd_o(:)           ! output grid:%land/gridcell
  real(r8)          , pointer    :: pctpft_o(:,:)         ! PFT cover (% of vegetated area)
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)    :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap
  real(r8), allocatable :: pctpft_i(:,:)      ! input grid: PFT percent
  integer  :: numpft_i                        ! num of plant types input data
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: wst(0:numpft)                   ! as pft_o at specific no
  real(r8) :: wst_sum                         ! sum of %pft
  real(r8) :: gpft_o(0:numpft)                ! output grid: global area pfts
  real(r8) :: garea_o                         ! output grid: global area
  real(r8) :: gpft_i(0:numpft)                ! input grid: global area pfts
  real(r8) :: garea_i                         ! input grid: global area
  integer  :: k,n,m,ni,no,ns_i,ns_o           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  integer  :: indxcrop                        ! input grid index for crop
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1

  character(len=35)  veg(0:numpft)            ! vegetation types
  character(len=32) :: subname = 'mkpft'
!-----------------------------------------------------------------------

  write (6,*)
  write (6,*) 'Attempting to make PFTs .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Set the vegetation types
  ! -----------------------------------------------------------------
  if ( numpft >= numstdpft )then
     veg(0:numstdpft) = (/                                &
                   'not vegetated                      ', &
                   'needleleaf evergreen temperate tree', &
                   'needleleaf evergreen boreal tree   ', &
                   'needleleaf deciduous boreal tree   ', &
                   'broadleaf evergreen tropical tree  ', &
                   'broadleaf evergreen temperate tree ', &
                   'broadleaf deciduous tropical tree  ', &
                   'broadleaf deciduous temperate tree ', &
                   'broadleaf deciduous boreal tree    ', &
                   'broadleaf evergreen shrub          ', &
                   'broadleaf deciduous temperate shrub', &
                   'broadleaf deciduous boreal shrub   ', &
                   'c3 arctic grass                    ', &
                   'c3 non-arctic grass                ', &
                   'c4 grass                           ', &
                   'c3_crop                            ', &
                   'c4_crop                            ' /)
     indxcrop = 15   ! c3_crop is active generic crop type
  end if
  if (      numpft == numstdpft )then
     write(6,*)'Creating surface datasets with the standard # of PFTs =', numpft
  else if ( numpft > numstdpft )then
     write(6,*)'Creating surface datasets with extra types for crops; total pfts =', numpft
  else
     write(6,*) subname//': parameter numpft is NOT set to a known value (should be 16 or more) =',numpft
     call abort()
  end if

  ! -----------------------------------------------------------------
  ! Read input PFT file
  ! -----------------------------------------------------------------

  ns_o = ldomain%ns

  if ( .not. use_input_pft ) then
     ! Obtain input grid info, read PCT_PFT

     call domain_read(tdomain,fpft)
     ns_i = tdomain%ns

     write (6,*) 'Open PFT file: ', trim(fpft)
     call check_ret(nf_open(fpft, 0, ncid), subname)

     call check_ret(nf_inq_dimid  (ncid, 'pft', dimid), subname)
     call check_ret(nf_inq_dimlen (ncid, dimid, numpft_i), subname)

     if (numpft_i .ne. numpft+1) then
        write(6,*) subname//': parameter numpft+1= ',numpft+1, &
             'does not equal input dataset numpft= ',numpft_i
        call abort()
     endif

     allocate(pctpft_i(ns_i,0:numpft), stat=ier)
     if (ier/=0) call abort()

     call check_ret(nf_inq_varid (ncid, 'PCT_PFT', varid), subname)
     call check_ret(nf_get_var_double (ncid, varid, pctpft_i), subname)

     call check_ret(nf_close(ncid), subname)

  else
     ns_i = 1
  end if

  ! Determine pctpft_o on output grid

  if ( zero_out ) then

     pctpft_o(:,:) = 0._r8
     pctlnd_o(:)   = 100._r8

  else if ( use_input_pft ) then

     call mkpft_check_oride( )

     ! set PFT based on input pft_frc and pft_idx
     pctpft_o(:,:) = 0._r8
     pctlnd_o(:)   = 100._r8
     do m = 0, numpft
        ! Once reach a PFT where fraction goes to zero -- exit
        if ( pft_frc(m) .eq. 0.0_r8 ) exit
        do no = 1,ns_o
           pctpft_o(no,pft_idx(m)) = pft_frc(m)
        end do
     end do

  else

     ! Compute pctlnd_o, pctpft_o

     call gridmap_mapread(tgridmap, mapfname)

     ! Error checks for domain and map consistencies

     call domain_checksame( tdomain, ldomain, tgridmap )
     ! Area-average percent cover on input grid [pctpft_i] to output grid 
     ! [pctpft_o] and correct [pctpft_o] according to land landmask
     ! Note that percent cover is in terms of total grid area.

     do no = 1,ns_o
        pctlnd_o(no)     = tgridmap%frac_dst(no) * 100._r8
        ldomain%frac(no) = tgridmap%frac_dst(no) 
     end do

     do m = 0,numpft
        call gridmap_areaave(tgridmap, pctpft_i(:,m), pctpft_o(:,m))
        do no = 1,ns_o
           if (pctlnd_o(no) < 1.0e-6) then
              if (m == 0) then
                 pctpft_o(no,m) = 100._r8
              else
                 pctpft_o(no,m) = 0._r8
              endif
           end if
        enddo
     enddo

     ! if irrigation dataset present, split into irrigated and 
     ! non-irrigated crop area
     if (firrig /= ' ') then
        write(6,*) 'Irrigation dataset present; splitting crop PFT into irrigated ',&
                   'and non-irrigated fractions'
        do no = 1,ns_o
           pctpft_o(no,IrrigIdx)    = min(pctpft_o(no,nonIrrigIdx),pctirr_o(no))
           pctpft_o(no,nonIrrigIdx) = pctpft_o(no,nonIrrigIdx) - pctpft_o(no,IrrigIdx)
        enddo
     endif

  end if

  ! Error check: percents should sum to 100 for land grid cells

  if ( .not. zero_out) then
     do no = 1,ns_o
        wst_sum = 0.
        do m = 0,numpft
           wst_sum = wst_sum + pctpft_o(no,m)
        enddo
        if (abs(wst_sum-100._r8) > 0.00001_r8) then
           write (6,*) subname//'error: pft = ', &
                (pctpft_o(no,m), m = 0, numpft), &
                ' do not sum to 100. at no = ',no,' but to ', wst_sum
           stop
        end if
     end do
  end if

  ! -----------------------------------------------------------------
  ! Error check
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  if ( .not. (zero_out .or. use_input_pft) ) then

     ! input grid

     gpft_i(:) = 0.
     garea_i   = 0.
     do ni = 1,ns_i
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        do m = 0, numpft
           gpft_i(m) = gpft_i(m) + pctpft_i(ni,m)*tgridmap%area_src(ni)*&
                                                  tgridmap%frac_src(ni)*re**2
        end do
     end do
     if ( allocated(pctpft_i) ) deallocate (pctpft_i)

     ! output grid

     gpft_o(:) = 0.
     garea_o   = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        do m = 0, numpft
           gpft_o(m) = gpft_o(m) + pctpft_o(no,m)*tgridmap%area_dst(no)*&
                                                  tgridmap%frac_dst(no)*re**2
        end do
     end do

     ! comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'PFTs Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'plant type     ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     do m = 0, numpft
        write (ndiag,1002) veg(m), gpft_i(m)*1.e-06/100.,gpft_o(m)*1.e-06/100.
     end do
1002 format (1x,a35,f16.3,f17.3)
     call shr_sys_flush(ndiag)

  end if

  ! Deallocate dynamic memory

  call domain_clean(tdomain) 
  if ( .not. zero_out .and. .not. use_input_pft ) then
     call gridmap_clean(tgridmap)
  end if

  write (6,*) 'Successfully made PFTs'
  write (6,*)


end subroutine mkpft

!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft_parse_oride
!
! !INTERFACE:
subroutine mkpft_parse_oride( string )
!
! !DESCRIPTION:
! Parse the string with pft fraction and index information on it, to override
! the file with this information rather than reading from a file.
!
! !USES:
   use shr_string_mod, only: shr_string_betweenTags, shr_string_countChar
! !ARGUMENTS:
   character(len=256), intent(IN) :: string  ! String to parse with PFT fraction 
                                             ! and index data
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: rc                         ! error return code
  integer :: num_elms                   ! number of elements
  character(len=256) :: substring       ! string between tags
  character(len=*), parameter :: frc_start = "<pft_f>"
  character(len=*), parameter :: frc_end   = "</pft_f>"
  character(len=*), parameter :: idx_start = "<pft_i>"
  character(len=*), parameter :: idx_end   = "</pft_i>"
  character(len=*), parameter :: subname = 'mkpft_parse_oride'
!-----------------------------------------------------------------------
  call shr_string_betweenTags( string, frc_start, frc_end, substring, rc )
  if ( rc /= 0 )then
     write(6,*) subname//'Trouble finding pft_frac start end tags'
     call abort()
  end if
  num_elms = shr_string_countChar( substring, ",", rc )
  read(substring,*) pft_frc(0:num_elms)
  call shr_string_betweenTags( string, idx_start, idx_end, substring, rc )
  if ( rc /= 0 )then
     write(6,*) subname//'Trouble finding pft_index start end tags'
     call abort()
  end if
  if ( num_elms /= shr_string_countChar( substring, ",", rc ) )then
     write(6,*) subname//'number of elements different between frc and idx fields'
     call abort()
  end if
  read(substring,*) pft_idx(0:num_elms)
!-----------------------------------------------------------------------

end subroutine mkpft_parse_oride

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft_check_oride
!
! !INTERFACE:
subroutine mkpft_check_oride( )
!
! !DESCRIPTION:
! Check that the pft override values are valid
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: i, j                         ! indices
  real(r8) :: sumpft                       ! Sum of pft_frc
  real(r8), parameter :: hndrd = 100.0_r8  ! A hundred percent
  character(len=32) :: subname = 'mkpft_check_oride:: '
!-----------------------------------------------------------------------

  sumpft = sum(pft_frc)
  if (          sumpft == 0.0 )then
    ! PFT fraction is NOT used
    use_input_pft = .false.
  else if ( abs(sumpft - hndrd) > 1.e-6 )then
    write(6,*) subname//'Sum of PFT fraction is NOT equal to 100% =', sumpft
    call abort()
  else
    use_input_pft = .true.
    nzero = 0
    do i = 0, numpft
       if ( pft_frc(i) == 0.0_r8 )then
          nzero = i
          exit
       end if
    end do
    ! PFT fraction IS used, and sum is OK, now check details
    do i = 0, nzero -1
      if ( pft_frc(i) < 0.0_r8 .or. pft_frc(i) > hndrd )then
         write(6,*) subname//'PFT fraction is out of range: pft_frc=', pft_frc(i)
         call abort()
      else if ( pft_frc(i) > 0.0_r8 .and. pft_idx(i) == -1 )then
         write(6,*) subname//'PFT fraction > zero, but index NOT set: pft_idx=', pft_idx(i)
         call abort()
      end if
      ! PFT index out of range
      if ( pft_idx(i) < 0 .or. pft_idx(i) > numpft )then
         write(6,*) subname//'PFT index is out of range: ', pft_idx(i)
         call abort()
      end if
      ! Make sure index values NOT used twice
      do j = 0, i-1
         if ( pft_idx(i) == pft_idx(j) )then
            write(6,*) subname//'Same PFT index is used twice: ', pft_idx(i)
            call abort()
         end if
      end do
    end do
    ! Make sure the rest of the fraction is zero and index are not set as well
    do i = nzero, numpft
      if ( pft_frc(i) /= 0.0_r8 .or. pft_idx(i) /= -1 )then
         write(6,*) subname//'After PFT fraction is zeroed out, fraction is non zero, or index set'
         call abort()
      end if
    end do
  end if

end subroutine mkpft_check_oride

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkirrig
!
! !INTERFACE:
subroutine mkirrig(ldomain, mapfname, datfname, ndiag, irrig_o)
!
! !DESCRIPTION:
! make percent irrigated area
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname       ! input mapping file name
  character(len=*)  , intent(in) :: datfname       ! input dataset file name
  integer           , intent(in) :: ndiag          ! unit number for diag out
  real(r8)          , intent(out):: irrig_o(:)     ! output grid: %irrigated area
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: David Lawrence
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)    :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap
  real(r8), allocatable :: irrig_i(:)         ! input grid: percent irrig
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: girrig_i                        ! input  grid: global irrig
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: girrig_o                        ! output grid: global irrig
  real(r8) :: garea_o                         ! output grid: global area
  integer  :: k,n,m,ni,no,ns_o,ns_i           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkirrig'
!-----------------------------------------------------------------------

  write (6,*)
  write (6,*) 'Attempting to make %irrigated area .....'
  call shr_sys_flush(6)

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain, datfname)
  ns_i = tdomain%ns

  write (6,*) 'Open irrigation file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(irrig_i(ns_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_IRRIG', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, irrig_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Do mapping from input to output grid

  call gridmap_areaave(tgridmap, irrig_i, irrig_o)

  do no = 1, ns_o
     if (irrig_o(no) < 1.) irrig_o(no) = 0.
  enddo

  ! Check for conservation

  do no = 1, ns_o
     if ((irrig_o(no)) > 100.000001_r8) then
        write (6,*) 'MKIRRIG error: irrigated area = ',irrig_o(no), &
                ' greater than 100.000001 for column, row = ',no
        call shr_sys_flush(6)
        stop
     end if
  enddo

  ! Global sum of output field -- must multiply by fraction of
  ! output grid that is land as determined by input grid

  sum_fldi = 0.
  do ni = 1, ns_i
     sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
  enddo

  sum_fldo = 0.
  do no = 1, ns_o
     sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
  end do

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKIRRIG error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  ! Input grid

  girrig_i = 0.
  garea_i = 0.

  do ni = 1,ns_i
     garea_i = garea_i + tgridmap%area_src(ni)
     girrig_i = girrig_i + irrig_i(ni)*(tgridmap%area_src(ni)/100.)*&
                                        tgridmap%frac_src(ni)*re**2
  end do

  ! Output grid

  girrig_o = 0.
  garea_o = 0.

  do no = 1,ns_o
     garea_o = garea_o + tgridmap%area_dst(no)
     girrig_o = girrig_o + irrig_o(no)*(tgridmap%area_dst(no)/100.)*&
                                        tgridmap%frac_dst(no)*re**2
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Irrigated area Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2002) girrig_i*1.e-06,girrig_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'irrigated area    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)
  call shr_sys_flush(ndiag)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (irrig_i)

  write (6,*) 'Successfully made %irrigated area'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkirrig

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpftAtt
!
! !INTERFACE:
subroutine mkpftAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! make PFT attributes on the output file
!
  use mkncdio    , only : check_ret, ncd_defvar
  use fileutils  , only : get_filename
  use mkvarctl   , only : mksrf_fvegtyp, mksrf_firrig, mksrf_flai, outnc_1d
  use mkvarpar   

! !ARGUMENTS:
  implicit none
  include 'netcdf.inc'
  integer, intent(in) :: ncid         ! NetCDF file ID to write out to
  logical, intent(in) :: dynlanduse   ! if dynamic land-use file
  integer, intent(in) :: xtype        ! external type to output real data as
!
! !CALLED FROM:
! subroutine mkfile in module mkfileMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: pftsize              ! size of lsmpft dimension
  integer :: dimid                ! input netCDF id's
  character(len=256) :: str       ! global attribute string
  character(len=32) :: subname = 'mkpftAtt'

  ! Define dimensions
  call check_ret(nf_def_dim (ncid, 'time'   , nf_unlimited, dimid), subname)

  pftsize = numpft + 1
  call check_ret(nf_def_dim (ncid, 'lsmpft' , pftsize     , dimid), subname)

  ! Add global attributes

  str = get_filename(mksrf_firrig)
  call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
       'Irrig_raw_data_file_name', len_trim(str), trim(str)), subname)

  if (.not. dynlanduse) then
     str = get_filename(mksrf_flai)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Lai_raw_data_file_name', len_trim(str), trim(str)), subname)
  end if

  if ( use_input_pft ) then
     str = 'TRUE'
     call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'pft_override', len_trim(str), trim(str)), subname)
  else if ( zero_out )then
     str = 'TRUE'
     call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'zero_out_pft_override', len_trim(str), trim(str)), subname)
  else
     str = get_filename(mksrf_fvegtyp)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
       'Vegetation_type_raw_data_filename', len_trim(str), trim(str)), subname)
  end if

  ! Define variables

  ! PCT_IRRIG
  if (mksrf_firrig /= ' ') then
     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='PCT_IRRIG', xtype=xtype, &
             dim1name='gridcell',&
             long_name='percent irrigated area', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='PCT_IRRIG', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='percent irrigated area', units='unitless')
     end if
  endif

  ! LANDFRAC_PFT
  if (outnc_1d) then
     call ncd_defvar(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
          dim1name='gridcell',&
          long_name='land fraction from pft dataset', units='unitless')
  else
     call ncd_defvar(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
          dim1name='lsmlon', dim2name='lsmlat', &
          long_name='land fraction from pft dataset', units='unitless')
  end if
  
  ! PFTDATA_MASK
  if (outnc_1d) then
     call ncd_defvar(ncid=ncid, varname='PFTDATA_MASK', xtype=nf_int, &
          dim1name='gridcell',&
          long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')
  else
     call ncd_defvar(ncid=ncid, varname='PFTDATA_MASK', xtype=nf_int, &
          dim1name='lsmlon', dim2name='lsmlat', &
          long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')
  end if
  
  ! PCT_PFT
  if (.not. dynlanduse) then
     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
             dim1name='gridcell', dim2name='lsmpft', &
             long_name='percent plant functional type of gridcell', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', &
             long_name='percent plant functional type of gridcell', units='unitless')
     end if
  else
     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
             dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
             long_name='percent plant functional type of gridcell', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
             long_name='percent plant functional type of gridcell', units='unitless')
     end if
  end if

  ! LAI,SAI,HTOP,HBOT
  if (.not. dynlanduse) then
     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='MONTHLY_LAI', xtype=xtype,  &
             dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
             long_name='monthly leaf area index', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='MONTHLY_LAI', xtype=xtype,  &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
             long_name='monthly leaf area index', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='MONTHLY_SAI', xtype=xtype,  &
             dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
             long_name='monthly stem area index', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='MONTHLY_SAI', xtype=xtype,  &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
             long_name='monthly stem area index', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
             dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
             long_name='monthly height top', units='meters')
     else
        call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
             long_name='monthly height top', units='meters')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
             dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
             long_name='monthly height bottom', units='meters')
     else
        call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
             long_name='monthly height bottom', units='meters')
     end if

  end if

  ! OTHER
  if (dynlanduse) then
     call ncd_defvar(ncid=ncid, varname='YEAR', xtype=nf_int,  &
            dim1name='time', &
            long_name='Year of PFT data', units='unitless')
     call ncd_defvar(ncid=ncid, varname='time', xtype=nf_int,  &
            dim1name='time', &
            long_name='year', units='unitless')
     call ncd_defvar(ncid=ncid, varname='input_pftdata_filename', xtype=nf_char,  &
            dim1name='nchar', &
            dim2name='time',  &
            long_name='Input filepath for PFT values for this year', units='unitless')
  else
     call ncd_defvar(ncid=ncid, varname='time', xtype=nf_int,  &
            dim1name='time', &
            long_name='Calendar month', units='month')
  end if

end subroutine mkpftAtt

!-----------------------------------------------------------------------

end module mkpftMod
