module RDycoreRestFile

#include <petsc/finclude/petsc.h>

  use shr_kind_mod   , only : r8 => shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use shr_sys_mod    , only : shr_sys_flush, shr_sys_abort
  use rdycore_varctl , only : iulog, caseid,inst_suffix
  use rdycoreSpmdMod , only : masterproc
  use rdycore
  use petsc

  implicit none
  private

  public :: RDycoreRestFileNameBase
  public :: RDycoreRestFileName
  public :: RDycoreWriteRestFile

  private :: RDycoreWritePfile

contains

  !------------------------------------------------------------------------
  character(len=1024) function RDycoreRestFileNameBase( rdate )
    use rdycore_varctl, only : inst_suffix

    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name

    RDycoreRestFileNameBase = "./"//trim(caseid)//".rdycore"//trim(inst_suffix)//".r."//trim(rdate)!//".h5"

  end function RDycoreRestFileNameBase

  !------------------------------------------------------------------------
  character(len=1024) function RDycoreRestFileName( rdate )
    use rdycore_varctl, only : inst_suffix

    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name

    RDycoreRestFileName = "./"//trim(caseid)//".rdycore"//trim(inst_suffix)//".r."//trim(rdate)//".h5"

  end function RDycoreRestFileName

  !-----------------------------------------------------------------------
  subroutine RDycoreWriteRestFile( rdy_, filename_base )
    !
    implicit none
    !
    type(RDy)                     :: rdy_
    character(len=*), intent (in) :: filename_base
    !
    ! !LCOAL VARIABLES:
    PetscErrorCode :: ierr

    PetscCallA(RDyWriteHDF5CheckpointFile(rdy_, filename_base, ierr))
    call RDycoreWritePfile(trim(filename_base)//'.h5')

  end subroutine RDycoreWriteRestFile
  !-----------------------------------------------------------------------
  subroutine RDycoreWritePfile( fnamer )
    !
    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !USES:
    use rdycore_varctl, only : rpntdir, rpntfil, inst_suffix
    use RDycoreFileUtils , only : relavu
    use RDycoreFileUtils , only : getavu, opnfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    !-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )

       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine RDycoreWritePfile

end module RDycoreRestFile
