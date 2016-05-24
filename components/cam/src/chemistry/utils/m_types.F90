
      module m_types
!---------------------------------------------------------------------
! 	... Derived types definition and related parameters
!---------------------------------------------------------------------

      implicit none

      type filespec
         character(len=168) :: &
            local_path, &                    ! local file path info
            remote_path, &                   ! remote path info (only used if NCAR is defined)
            nl_filename                      ! filename
         character(len=30) :: &
            hor_res                          ! horizontal resolution
      end type filespec

      type time_ramp
	 character(len=8) :: type
	 integer          :: cycle_yr
	 integer          :: fixed_ymd
	 integer          :: fixed_tod
      end type time_ramp

      end module m_types
