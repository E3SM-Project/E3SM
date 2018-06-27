

Module gcam_comp_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: gcam_comp_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:

  use iac_fields_mod

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: gcam_init_mod               ! gcam initialization
  public :: gcam_run_mod                ! gcam run phase
  public :: gcam_setdensity_mod         ! gcam set density phase
  public :: gcam_final_mod              ! gcam finalization/cleanup

! !PUBLIC DATA MEMBERS: None


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication


! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_init_mod

! !INTERFACE:
  subroutine gcam_init_mod( EClock, cdata, gcami, gcamo, gcamoemis)

! !DESCRIPTION:
! Initialize interface for gcam

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: gcami(:,:)
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(gcam_init_mod)'
    character(len=128) :: casename
    integer,save :: iulog
    integer      :: i

! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------
  iulog = cdata%i(iac_cdatai_logunit)
  casename = trim(cdata%c(iac_cdatac_casename))


  cdata%l(iac_cdatal_gcam_present) = .true.
  cdata%l(iac_cdatal_gcam_prognostic) = .true.
  cdata%i(iac_cdatai_gcam_nreg) = iac_gcam_nreg
  cdata%i(iac_cdatai_gcam_naez) = iac_gcam_naez
  cdata%i(iac_cdatai_gcam_timestep) = iac_gcam_timestep
  cdata%i(iac_cdatai_gcamo_ntime) = iac_gcamo_ntime
  cdata%i(iac_cdatai_gcamo_nflds) = iac_gcamo_nflds
  cdata%i(iac_cdatai_gcamo_size) = iac_gcam_nreg*iac_gcam_naez * iac_gcamo_ntime
  cdata%i(iac_cdatai_gcamoemis_size) = iac_gcam_nreg*iac_gcam_nsector
  cdata%i(iac_cdatai_gcami_nflds) = iac_gcami_nflds
  cdata%i(iac_cdatai_gcami_size) = iac_gcam_nreg*iac_gcam_naez*iac_gcam_ncrops

  allocate(gcami(iac_gcami_nflds,cdata%i(iac_cdatai_gcami_size)))
  allocate(gcamo(iac_gcamo_nflds,cdata%i(iac_cdatai_gcamo_size)))
  allocate(gcamoemis(iac_gcamoemis_nemis,cdata%i(iac_cdatai_gcamoemis_size)))

  gcami = iac_spval
  gcamo = iac_spval
  gcamoemis = iac_spval

  ! create CCSM_GCAM_interface Object 
  call initCCSMInterface()

  ! Call initcGCAM method of CCSM/GCAM Interface 
  call initcGCAM()

  end subroutine gcam_init_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_run_mod

! !INTERFACE:
  subroutine gcam_run_mod( EClock, cdata, gcami, gcamo, gcamoemis)

! !DESCRIPTION:
! Run interface for gcam

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: gcami(:,:)
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: iu
    integer :: i,j
    character(len=*),parameter :: subname='(gcam_run_mod)'


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------

  restart_now = cdata%l(iac_cdatal_rest)
  iu  = cdata%i(iac_cdatai_logunit)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iu,*) trim(subname),' date= ',ymd,tod

  !  Call runcGCAM method of CCSM Interface 
  call runcGCAM(ymd,tod,gcami,size(gcami,dim=1),size(gcami,dim=2),gcamo,size(gcamo,dim=1),size(gcamo,dim=2),gcamoemis,size(gcamoemis,dim=1),size(gcamoemis,dim=2), cdata%i(iac_cdatai_gcam_yr1),cdata%i(iac_cdatai_gcam_yr2),cdata%l(iac_cdatal_sneakermode),cdata%l(iac_cdatal_write_rest))

  end subroutine gcam_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_setdensity_mod

! !INTERFACE:
  subroutine gcam_setdensity_mod( EClock, cdata, gcami)

! !DESCRIPTION:
! Setdensity interface for gcam

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: gcami(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: iu
    integer :: i,j
    character(len=*),parameter :: subname='(gcam_setdensity_mod)'


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------

  restart_now = cdata%l(iac_cdatal_rest)
  iu  = cdata%i(iac_cdatai_logunit)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iu,*) trim(subname),' date= ',ymd,tod

  !  Call setdensity method of CCSM Interface 
  call setdensitycGCAM(ymd,tod,gcami,size(gcami,dim=1),size(gcami,dim=2))

  end subroutine gcam_setdensity_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_final_mod

! !INTERFACE:
  subroutine gcam_final_mod( )

! !DESCRIPTION:
! Finalize gcam model
!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:

! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!---------------------------------------------------------------------------

  !  Cleanup GCAM 
  call finalizecGCAM()

  !  Cleanup CCSM Interface Object 
  call deleteCCSMInterface()

  end subroutine gcam_final_mod

!====================================================================================

end module gcam_comp_mod
