

   ! mpas_log parameters

   integer, parameter :: MPAS_LOG_STD = 1  !< code for message type "standard"
   integer, parameter :: MPAS_LOG_WARN = 2  !< code for message type "warning"
   integer, parameter :: MPAS_LOG_ERR = 3  !< code for message type "error"
   integer, parameter :: MPAS_LOG_CRIT = 4  !< code for message type "critical error"

   ! Derived type for storing information related to an instance of the logging manager
   type mpas_log_type


      integer :: logUnit  !< file unit number for the standard log file
      logical :: logActive !< logical indicating if the log file is active on this processor
      integer :: errUnit  !< file unit number for the error/abort log file
      logical :: errUnitOpen !< logical indicating if the err file has been already opened

      integer :: taskID !< task number for this instance of the log manager
                   !< (stored here to eliminate the need for dminfo later)
      integer :: nTasks !<  number of total tasks associated with this instance
                   !< (stored here to eliminate the need for dminfo later)
      character(len=StrKIND) :: coreName  !< name of the core to which this log manager instance belongs

      character(len=StrKIND) :: logFileName  !< Filename for the log file
      character(len=StrKIND) :: errFileName  !< Filename for the abort file

      integer :: standardMessageCount  !< counter for number of standard messages printed during the run
      integer :: warningMessageCount  !< counter for number of warning messages printed during the run
      integer :: errorMessageCount  !< counter for number of error messages printed during the run
      integer :: criticalErrorMessageCount  !< counter for number of critical error messages printed during the run


   end type mpas_log_type


