!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_files

!BOP
! !MODULE: glc_files

! !DESCRIPTION:
!  Manage file names of some special files: namelists, restart pointer file, etc.
!
! !REVISION HISTORY:
!
! !USES:

   use shr_kind_mod,        only: CL=>SHR_KIND_CL

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: set_filenames


!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

  character(CL), public ::  &
      nml_filename  , & ! namelist input file name
      ionml_filename, & ! model IO namelist file name
      ptr_filename      ! restart pointer file name

!EOP
!BOC
!EOC
!***********************************************************************
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: set_filenames
! !INTERFACE:
   subroutine set_filenames
!
! !DESCRIPTION:
! Set module variables that give various file names
!
! This should be done in model initialization, after the ensemble-related variables are
! initialized.
!
! !USES:
     use glc_ensemble, only : get_inst_suffix
!
! !ARGUMENTS:
!
! !LOCAL VARIABLES:
     character(len=16) :: inst_suffix
!EOP
!-----------------------------------------------------------------------
     
     call get_inst_suffix(inst_suffix)

     nml_filename = 'cism_in'//trim(inst_suffix)
     ionml_filename = 'glc_modelio.nml'//trim(inst_suffix)
     ptr_filename = 'rpointer.glc'//trim(inst_suffix)

   end subroutine set_filenames

end module glc_files
