

   ! mpas_log parameters

   integer, parameter :: MPAS_LOG_OUT = 1  !< code for message type "output"
   integer, parameter :: MPAS_LOG_WARN = 2  !< code for message type "warning"
   integer, parameter :: MPAS_LOG_ERR = 3  !< code for message type "error"
   integer, parameter :: MPAS_LOG_CRIT = 4  !< code for message type "critical error"

   type mpas_log_file_type
      integer :: unitNum  !< file unit number for the standard log file
      character(len=StrKIND) :: fileName  !< Filename for the log file
      logical :: isOpen !< logical indicating if the log file has already been opened
      logical :: isActive !< logical indicating if the log file is active on this processor
      logical :: openedByLogModule !< logical indicating if the log file was opened by the log manager;
         !< false if a higher level driver opened it OR file was never opened
   end type mpas_log_file_type


   ! Derived type for storing information related to an instance of the logging manager
   type mpas_log_type
      type(mpas_log_file_type), pointer :: outputLog
      type(mpas_log_file_type), pointer :: errorLog

      integer :: taskID !< task number for this instance of the log manager
                   !< (stored here to eliminate the need for dminfo later)
      integer :: nTasks !<  number of total tasks associated with this instance
                   !< (stored here to eliminate the need for dminfo later)
      character(len=StrKIND) :: coreName  !< name of the core to which this log manager instance belongs

      integer :: outputMessageCount  !< counter for number of output messages printed during the run
      integer :: warningMessageCount  !< counter for number of warning messages printed during the run
      integer :: errorMessageCount  !< counter for number of error messages printed during the run
      integer :: criticalErrorMessageCount  !< counter for number of critical error messages printed during the run


   end type mpas_log_type


