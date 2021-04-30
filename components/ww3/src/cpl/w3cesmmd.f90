!/ ------------------------------------------------------------------- /
      MODULE W3CESMMD

!/ ------------------------------------------------------------------- /
!/
      private

      ! runtype is used by W3SRCE (values are startup, branch, continue)
      character(len=16),public :: runtype

      ! if a run is a startup or branch run, then initfile is used
      ! to construct the initial file and used in W3IORSMD
      character(len=256), public :: initfile

      ! if a run is a continue run, then casename is used to construct
      ! the restart filename in W3IORSMD
      character(len=256), public :: casename

      logical, public :: rstwr   ! true => write restart at end of day

      integer, public :: stdout  ! output log file

      integer, public                  :: inst_index            ! number of current instance (ie. 1)
      character(len=16), public :: inst_name   ! fullname of current instance (ie. "wav_0001")
      character(len=16), public :: inst_suffix ! char string associated with instance
!/
!/ End of module W3CESMMD -------------------------------------------- /
!/
      END MODULE W3CESMMD
