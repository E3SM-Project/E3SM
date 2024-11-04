module mkFertMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkFertMod
!
! !DESCRIPTION:
! Make fertilizer data (nitrogen and phosphorus)
!
! !REVISION HISTORY:
! Author: Beth Drewniak
!
!-----------------------------------------------------------------------
! USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  implicit none

  SAVE
  private
!
! PUBLIC MEMBER FUNCTIONS:
  public  :: mkfert

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkfert
!
! !INTERFACE:
subroutine mkfert(ldomain, mapfname, datfname, ndiag, nfert_o, pfert_o)
!
! !DESCRIPTION:
! Make fertilizer data
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkncdio
  use mkvarctl    , only : numpft
  use mkpftMod    , only : num_cft, cft_lb, cft_ub
  use mkchecksMod , only : min_bad
!
! !ARGUMENTS:
  implicit none
  type(domain_type) , intent(in)  :: ldomain
  character(len=*)  , intent(in)  :: mapfname      ! input mapping file name
  character(len=*)  , intent(in)  :: datfname      ! input data file name
  integer           , intent(in)  :: ndiag         ! unit number for diag out
  real(r8)          , intent(out) :: nfert_o(:,:) ! output netcdf file id
  real(r8)          , intent(out) :: pfert_o(:,:) ! output netcdf file id
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Beth Drewniak
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain          ! local domain
  integer  :: l                             ! index
  integer  :: num_cft_i                      ! number of plant types on input
  real(r8), allocatable :: mfert_i(:,:)     ! data on input grid
  integer  :: ncidi, dimid, varid           ! input netCDF id's
  integer  :: ier                           ! error status

  character(len= 32) :: subname = 'mkfert'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make fertilizer .....'
  call shr_sys_flush(6)

  ! --------------------------------------------------------------------
  ! Read domain file, mapping weights and map, check for consistancies
  ! --------------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)

  call gridmap_mapread(tgridmap, mapfname)
  call gridmap_check( tgridmap, subname)
  call domain_checksame( tdomain, ldomain, tgridmap )

  ! --------------------------------------------------------------------
  ! Open input file, read dimensions, check for errors  and allocate memory 
  ! --------------------------------------------------------------------

  write (6,*) 'Open fert file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)
  call check_ret(nf_inq_dimid(ncidi, 'cft', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, num_cft_i), subname)

  if (num_cft_i /= num_cft) then
     write(6,*)'MKFERT: parameter numcft= ',num_cft, &
          'does not equal input dataset numcft= ',num_cft_i
     stop
  endif

  allocate(mfert_i(tdomain%ns,1:num_cft), stat=ier )
  if (ier /= 0) then
     write(6,*)'mkfert allocation error'; call abort()
  end if

  ! --------------------------------------------------------------------
  ! Read in nitrogen fertilizer
  ! --------------------------------------------------------------------

  call check_ret(nf_inq_varid (ncidi, 'NFERT', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, mfert_i), subname)
  
     ! Loop over pft types to do mapping

     do l = 1,num_cft
        call gridmap_areaave(tgridmap, mfert_i(:,l) , nfert_o(:,l) , nodata=0._r8)
        if (min_bad(nfert_o(:,l), 0.0_r8, 'NFERT')) then
            stop
        end if

     enddo

  ! --------------------------------------------------------------------
  ! Read in phosphorus fertilizer
  ! --------------------------------------------------------------------

  call check_ret(nf_inq_varid (ncidi, 'PFERT', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, mfert_i), subname)

     ! Loop over pft types to do mapping

     do l = 1,num_cft
        call gridmap_areaave(tgridmap, mfert_i(:,l) , pfert_o(:,l) , nodata=0._r8)
        if (min_bad(pfert_o(:,l), 0.0_r8, 'PFERT')) then
           stop
        end if

     enddo

  ! --------------------------------------------------------------------
  ! close file and deallocate memory
  ! --------------------------------------------------------------------

  call check_ret(nf_close(ncidi), subname)

  deallocate(mfert_i)
  call gridmap_clean(tgridmap)
  call domain_clean(tdomain) 

  write (6,*) 'Successfully created Fertilizer data'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkfert

  
end module mkFertMod
