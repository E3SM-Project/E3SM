MODULE moby_parms

  !-----------------------------------------------------------------------------
  !   This module manages the parameter variables for the module moby_mod.
  !
  !   It is based upon ecosys_mod written by Keith Lindsay.
  !
  !   Most of the variables are not parameters in the Fortran sense. In the
  !   the Fortran sense, they are vanilla module variables.
  !
  !   This modules handles initializing the variables to default values and
  !   reading them from the namelist moby_parms. The values used are echoed
  !   to stdout for record keeping purposes.
  !
  !-----------------------------------------------------------------------------

  USE exit_mod, ONLY : sigAbort, exit_POP
  USE communicate, ONLY : my_task, master_task
  USE constants, ONLY : c1
  USE kinds_mod
  USE io_tools, ONLY : document

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !   all module variables are public and should have their values preserved
  !-----------------------------------------------------------------------------

  PUBLIC
  SAVE

  !-----------------------------------------------------------------------------
  !   floating point constants used across moby module
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  !   moby parameters accessible via input file
  !----------------------------------------------------------------------------

  REAL(KIND=r8) :: &
       parm_Fe_bioavail       ! fraction of Fe flux that is bioavailable

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  !     attenuation coefficients for PAR and related parameters
  !---------------------------------------------------------------------
  real(kind=r8), parameter :: &
       k_chl = 0.03e-2_r8, & ! Chl atten. coeff. (1/cm/(mg Chl/m^3))
       k_h2o = 0.04e-2_r8, & ! water atten. coeff (1/cm)
       f_qsw_par = 0.45_r8   ! PAR fraction

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  subroutine moby_parms_init

    USE io_types, ONLY: stdout, nml_in, nml_filename
    USE broadcast, ONLY : broadcast_scalar

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: subname = 'moby_parms:moby_parms_init'

    LOGICAL(KIND=log_kind) :: &
         lnml_found             ! Was moby_parms_nml found ?

    NAMELIST /moby_parms_nml/ &
         parm_Fe_bioavail

    !---------------------------------------------------------------------------
    !   default namelist settings
    !---------------------------------------------------------------------------

    parm_Fe_bioavail    = 0.02_r8

    !---------------------------------------------------------------------------
    !   read in namelist
    !---------------------------------------------------------------------------

    IF (my_task == master_task) THEN
       lnml_found = .FALSE.
       OPEN(UNIT=nml_in, FILE=nml_filename, STATUS='OLD')
10     CONTINUE
       READ(UNIT=nml_in, NML=moby_parms_nml, ERR=10, END=20)
       CLOSE(UNIT=nml_in)
       lnml_found = .TRUE.
20     CONTINUE
    END IF

    CALL broadcast_scalar(lnml_found, master_task)
    IF (.NOT. lnml_found) THEN
       CALL document(subname, 'moby_parms_nml not found')
       CALL exit_POP(sigAbort, 'ERROR : stopping in ' // subname)
    END IF

    !---------------------------------------------------------------------------
    !   broadcast all namelist variables
    !---------------------------------------------------------------------------

    CALL broadcast_scalar(parm_Fe_bioavail, master_task)

    !---------------------------------------------------------------------------
    !   echo all namelist variables to stdout
    !---------------------------------------------------------------------------

    IF (my_task == master_task) THEN
       WRITE (stdout,*) '----------------------------------------'
       WRITE (stdout,*) '----- moby_parms_nml namelist values -----'
       WRITE (stdout,*) 'parm_Fe_bioavail    = ', parm_Fe_bioavail
       WRITE (stdout,*) '----------------------------------------'
    END IF

  end subroutine moby_parms_init

  !*****************************************************************************

END MODULE moby_parms
