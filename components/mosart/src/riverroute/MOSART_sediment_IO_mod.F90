!
MODULE MOSART_sediment_IO_mod
! Description: core code of MOSART-sediment. Can be incoporated within any land model via a interface module
! 
! Developed by Hongyi Li, 03/2015.
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
  use RunoffMod     , only : Tctl, TUnit, rtmCTL
  use RtmSpmd       , only : masterproc
  use RtmVar        , only : iulog
  use RtmIO         , only : pio_subsystem, ncd_pio_openfile, ncd_pio_closefile
  use rof_cpl_indices, only : nt_rtm
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_flush, shr_sys_abort
  use MOSART_BGC_type  , only : TSedi, TSedi_para
  use netcdf
  use pio

  implicit none
  private

  public SediInput_annual


! !PUBLIC MEMBER FUNCTIONS:
    contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read soil erosion annual inputs
!
! !INTERFACE:
  subroutine SediInput_annual(iYear)
!
! !DESCRIPTION:
! read the input time series for MOSART-sediment

!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li 
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer, intent(in)  :: iYear
    character(len=4) :: strYear
    character(len=1000) :: fname
  type(file_desc_t)  :: ncid       ! pio file desc
  type(var_desc_t)   :: vardesc    ! pio variable desc 
  type(io_desc_t)    :: iodesc_dbl ! pio io desc
  type(io_desc_t)    :: iodesc_int ! pio io desc
  integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
  integer :: dids(2)               ! variable dimension ids 
  integer :: dsizes(2)             ! variable dimension lengths
  integer :: ier                  ! error code
  integer :: begr, endr, iunit, nn, n, cnt, nr, nt
  integer :: numDT_r, numDT_t
  integer :: lsize, gsize
  integer :: igrow, igcol, iwgt
  integer :: tcnt
  character(len=16384) :: rList             ! list of fields for SM multiply
  character(len=*),parameter :: FORMI = '(2A,2i10)'
  character(len=*),parameter :: FORMR = '(2A,2g15.7)'
    character(len=*),parameter :: subname='(MOSART_sediment_annual_inputs)'
     write(strYear,'(I4.4)') iYear
     fname = trim(TSedi_para%inputPath) // 'Soil_erosion_para_annual_' //strYear//'.nc'
     write(iulog,*) subname, ' reading ',trim(fname)

    begr = rtmCTL%begr
    endr = rtmCTL%endr
    if(endr >= begr) then
        ! reservoir parameters
        call ncd_pio_openfile (ncid, fname, 0)
        call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
        allocate(compdof(rtmCTL%lnumr))
        cnt = 0
        do n = rtmCTL%begr,rtmCTL%endr
           cnt = cnt + 1
           compDOF(cnt) = rtmCTL%gindex(n)
        enddo

        ! setup iodesc based on annualQsur dids
        ier = pio_inq_varid(ncid, name='peakQtot', vardesc=vardesc)
        ier = pio_inq_vardimid(ncid, vardesc, dids)
        ier = pio_inq_dimlen(ncid, dids(1),dsizes(1))
        ier = pio_inq_dimlen(ncid, dids(2),dsizes(2))
        call pio_initdecomp(pio_subsystem, pio_double, dsizes, compDOF, iodesc_dbl)
        call pio_initdecomp(pio_subsystem, pio_int   , dsizes, compDOF, iodesc_int)
        deallocate(compdof)
        call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)

        call ncd_pio_openfile(ncid, trim(fname), 0)
        ier = pio_inq_varid (ncid, name='peakQtot', vardesc=vardesc)
        call pio_read_darray(ncid, vardesc, iodesc_dbl, TSedi_para%peak_Qtot, ier)
        call ncd_pio_closefile(ncid)

        do nr=rtmCTL%begr,rtmCTL%endr
           if (TSedi_para%peak_Qtot(nr).lt.0._r8) then
              TSedi_para%peak_Qtot(nr) = 0._r8
           end if
        end do

    end if 
     write(iulog,FORMR) trim(subname),' read annualQtot', minval(TSedi_para%peak_Qtot),maxval(TSedi_para%peak_Qtot)
     call shr_sys_flush(iulog)

  end subroutine SediInput_annual

!-----------------------------------------------------------------------


end MODULE MOSART_sediment_IO_mod 