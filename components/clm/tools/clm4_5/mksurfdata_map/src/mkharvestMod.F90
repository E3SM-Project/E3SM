module mkharvestMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkharvest
!
! !DESCRIPTION:
! Make harvest and grazing data to add to the dynamic PFT file.
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_sys_mod  , only : shr_sys_flush
  use mkdomainMod  , only : domain_checksame

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mkharvest_init        ! Initialization
  public mkharvest             ! Calculate the harvest values on output grid
  public mkharvest_fieldname   ! Field name
  public mkharvest_longname    ! Long name
  public mkharvest_numtypes    ! Number of harvest types
  public mkharvest_parse_oride ! Parse the over-ride string

! !PRIVATE DATA MEMBERS:

  integer, parameter :: numharv = 6   ! number of harvest and grazing fields
  integer, parameter :: harlen  = 12  ! length of strings for harvest fieldnames
  character(len=harlen), parameter  :: harvest_fieldnames(numharv) = (/ &
                                                        'HARVEST_VH1',  &
                                                        'HARVEST_VH2',  &
                                                        'HARVEST_SH1',  &
                                                        'HARVEST_SH2',  &
                                                        'HARVEST_SH3',  &
                                                        'GRAZING    '   &
                                                      /)
  character(len=CL), parameter :: string_undef = 'STRING_UNDEFINED'
  real(r8),          parameter :: real_undef   = -999.99
  character(len=CL), save :: harvest_longnames(numharv) = string_undef
  real(r8),  pointer     :: oride_harv(:)        ! array that can override harvesting


!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_init
!
! !INTERFACE:
  subroutine mkharvest_init( ns_o, init_val, harvest, fharvest )
!
! !DESCRIPTION:
!              Initialization of mkharvest module.
!
! !USES:
    use mkncdio 
    implicit none
!
! !ARGUMENTS:
    integer         , intent(in) :: ns_o          ! clm output grid resolution
    real(r8)        , intent(in) :: init_val      ! initial value to set to
    real(r8)        , pointer    :: harvest(:,:)  ! output grid: normalized harvesting
    character(len=*), intent(in) :: fharvest      ! input harvest dataset file name
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_init'
    integer  :: ncid,varid                      ! input netCDF id's
    integer  :: ifld                            ! indices
!EOP
!-----------------------------------------------------------------------

    allocate(harvest(ns_o,numharv))
    harvest(:,:) = init_val

    call check_ret(nf_open(fharvest, 0, ncid), subname)
    do ifld = 1, numharv
       call check_ret(nf_inq_varid (   ncid, mkharvest_fieldname(ifld), varid),            subname)
       call check_ret(nf_get_att_text( ncid, varid, 'long_name', harvest_longnames(ifld)), subname)
    end do

    call check_ret(nf_close(ncid), subname)

    allocate( oride_harv(numharv) )
    oride_harv(:) = real_undef

  end subroutine mkharvest_init

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_fieldname
!
! !INTERFACE:
  character(len=harlen) function mkharvest_fieldname( nfield )
!
! !DESCRIPTION:
!              Return harvest fieldname of input field number.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_fieldname'
!EOP
!-----------------------------------------------------------------------

      if (      nfield < 1       )then
         write(6,*) subname, ' error nfield < 1'
         call abort()
      else if ( nfield > numharv )then
         write(6,*) subname, ' error nfield > max fields'
         call abort()
      else
         mkharvest_fieldname = harvest_fieldnames(nfield)
      end if

  end function mkharvest_fieldname

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_longname
!
! !INTERFACE:
  character(len=CL) function mkharvest_longname( nfield )
!
! !DESCRIPTION:
!              Return longname description of given input field number.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_longname'
!EOP
!-----------------------------------------------------------------------

      if (      nfield < 1       )then
         write(6,*) subname, ' error nfield < 1'
         call abort()
      else if ( nfield > numharv )then
         write(6,*) subname, ' error nfield > max fields'
         call abort()
      else
         if ( trim(harvest_longnames(nfield)) .eq. trim(string_undef) )then
            write(6,*) subname, ' error harvest_longnames not set yet'
            call abort()
         end if
         mkharvest_longname = harvest_longnames(nfield)
      end if

  end function mkharvest_longname

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_numtypes
!
! !INTERFACE:
  integer function mkharvest_numtypes( )
!
! !DESCRIPTION:
!              Return number of different harvest field types.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    character(len=*), parameter :: subname = 'mkharvest_numtypes'
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
      mkharvest_numtypes = numharv

  end function mkharvest_numtypes

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest
!
! !INTERFACE:
subroutine mkharvest(ldomain, mapfname, datfname, ndiag, harv_o)
!
! !DESCRIPTION:
! Make harvest data for the dynamic PFT dataset.
! This dataset consists of the normalized harvest or grazing fraction (0-1) of
! the model.
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
  character(len=*)  , intent(in) :: mapfname    ! input mapping file name
  character(len=*)  , intent(in) :: datfname    ! input data file name
  integer           , intent(in) :: ndiag       ! unit number for diag out
  real(r8)          , intent(out):: harv_o(:,:) ! output grid: normalized harvesting
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
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8), allocatable :: harv_i(:,:)        ! input grid: harvest/grazing percent
  real(r8), allocatable :: pctlnd_o(:)        ! output grid: percent land 
  real(r8) :: gharv_o(numharv)                ! output grid: global area harvesting
  real(r8) :: garea_o                         ! output grid: global area
  real(r8) :: gharv_i(numharv)                ! input grid: global area harvesting
  real(r8) :: garea_i                         ! input grid: global area
  integer  :: ifld                            ! indices
  integer  :: k,n,m,ni,no,ns_i,ns_o           ! indices
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status

  character(len=*), parameter :: unit = '10**6 km**2' ! Output units
  real(r8), parameter :: fac = 1.e-06_r8              ! Output factor
  real(r8), parameter :: rat = fac/100._r8            ! Output factor divided by 100%
  character(len=*), parameter :: subname = 'mkharvest'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make harvest fields .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Normally read in the harvesting file, and then regrid to output grid
  ! -----------------------------------------------------------------
  
  if ( all(oride_harv == real_undef ) )then

     ! -----------------------------------------------------------------
     ! Read input harvesting file
     ! -----------------------------------------------------------------

     ! Obtain input grid info, read HARVEST_VH1, HARVEST_VH2, ... GRAZING etc.

     call domain_read(tdomain,datfname)
     ns_i = tdomain%ns
     allocate(harv_i(ns_i,1:numharv), stat=ier)
     if (ier/=0) call abort()
     ns_o = ldomain%ns

     write (6,*) 'Open harvest file: ', trim(datfname)
     call check_ret(nf_open(datfname, 0, ncid), subname)
     do ifld = 1, numharv
        call check_ret(nf_inq_varid(ncid, mkharvest_fieldname(ifld), varid), subname)
        call check_ret(nf_get_var_double (ncid, varid, harv_i(:,ifld)), subname)
     end do
     call check_ret(nf_close(ncid), subname)

     ! Area-average normalized harvest on input grid [harv_i] to output grid [harv_o]

     call gridmap_mapread(tgridmap, mapfname )

     ! Error checks for domain and map consistencies
     
     call domain_checksame( tdomain, ldomain, tgridmap )

     ! Determine harv_o on output grid

     do ifld = 1,numharv
        call gridmap_areaave(tgridmap, harv_i(:,ifld), harv_o(:,ifld), nodata=0._r8)
     end do

     ! -----------------------------------------------------------------
     ! Error check
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     gharv_i(:) = 0.
     garea_i = 0.
     do ni = 1, ns_i
        garea_i = garea_i + tgridmap%area_src(ni)*re**2
        do m = 1, numharv
           gharv_i(m) = gharv_i(m) + harv_i(ni,m)*tgridmap%area_src(ni)* &
                                                  tgridmap%frac_src(ni)*re**2
        end do
     end do

     gharv_o(:) = 0.
     garea_o = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)*re**2
        do m = 1, numharv
           gharv_o(m) = gharv_o(m) + harv_o(no,m)*tgridmap%area_dst(no)* &
                                                  tgridmap%frac_dst(no)*re**2
        end do
     end do

     ! Write out to diagnostic output file
     !

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Harvesting Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001) unit, unit
1001 format (1x,'harvest type   ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     ',A,'      ',A)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     do m = 1, numharv
        write (ndiag,1002) mkharvest_fieldname(m), gharv_i(m)*rat,gharv_o(m)*rat
     end do
1002 format (1x,a35,f16.3,f17.3)

     ! Deallocate dynamic memory

     call domain_clean(tdomain) 
     call gridmap_clean(tgridmap)
     deallocate (harv_i)

  else

     ! -----------------------------------------------------------------
     ! Otherwise override the harvesting with the input harvest values
     ! -----------------------------------------------------------------
  
     if ( any(oride_harv == real_undef ) )then
         write(6,*) subname, ' error some override harvesting fields set ', &
                    'and others are not = ', oride_harv
         call abort()
     end if
     do m = 1, numharv
        if ( oride_harv(m) < 0.0_r8 .or. oride_harv(m) > 100.0_r8 )then
            write(6,*) subname, ' error override harvesting field out of range', &
                       oride_harv(m), ' field = ', mkharvest_fieldname(m)
            call abort()
        end if
     end do
     do no = 1,ns_o
        do m = 1, numharv
           harv_o(no,m) = oride_harv(m)
        end do
     end do

  end if

  write (6,*) 'Successfully made harvest and grazing'
  write (6,*)

end subroutine mkharvest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_parse_oride
!
! !INTERFACE:
subroutine mkharvest_parse_oride( string )
!
! !DESCRIPTION:
! Parse the string with harvest and grazing information on it, to override
! the file with this information rather than reading from a file.
!
! !USES:
   use shr_string_mod, only: shr_string_betweenTags
! !ARGUMENTS:
   character(len=256), intent(IN) :: string  ! String to parse with harvest and grazing data
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
  character(len=256) :: substring       ! substring between tags
  character(len=*), parameter :: harv_start = "<harv>"
  character(len=*), parameter :: harv_end   = "</harv>"
  character(len=*), parameter :: graz_start = "<graz>"
  character(len=*), parameter :: graz_end   = "</graz>"
  character(len=*), parameter :: subname = 'mkharvest_parse_oride'
!-----------------------------------------------------------------------
  call shr_string_betweenTags( string, harv_start, harv_end, substring, rc )
  if ( rc /= 0 )then
     write(6,*) subname//'Trouble finding harvest start end tags'
     call abort()
  end if
  read(substring,*) oride_harv(1:numharv-1)
  call shr_string_betweenTags( string, graz_start, graz_end, substring, rc )
  if ( rc /= 0 )then
     write(6,*) subname//'Trouble finding grazing start end tags'
     call abort()
  end if
  read(substring,*) oride_harv(numharv)
  if ( harvest_fieldnames(numharv) /= 'GRAZING' )then
     write(6,*) subname, ' grazing is NOT last field as was expected'
     call abort()
  end if

!-----------------------------------------------------------------------

end subroutine mkharvest_parse_oride

!-----------------------------------------------------------------------

end module mkharvestMod
