module gptl
! GPTL module file for user code. Parameter values match their counterparts
! in gptl.h. This file also contains an interface block for parameter checking.

  implicit none
  public

! User-accessible integers 

  integer, parameter :: GPTLsync_mpi       = 0
  integer, parameter :: GPTLwall           = 1
  integer, parameter :: GPTLcpu            = 2
  integer, parameter :: GPTLabort_on_error = 3
  integer, parameter :: GPTLoverhead       = 4
  integer, parameter :: GPTLdepthlimit     = 5
  integer, parameter :: GPTLverbose        = 6
  integer, parameter :: GPTLnarrowprint    = 7
  integer, parameter :: GPTLpercent        = 9
  integer, parameter :: GPTLpersec         = 10
  integer, parameter :: GPTLmultiplex      = 11
  integer, parameter :: GPTLdopr_preamble  = 12
  integer, parameter :: GPTLdopr_threadsort= 13
  integer, parameter :: GPTLdopr_multparent= 14
  integer, parameter :: GPTLdopr_collision = 15
  integer, parameter :: GPTLdopr_memusage  = 27
  integer, parameter :: GPTLprint_method   = 16
  integer, parameter :: GPTLtablesize      = 50
  integer, parameter :: GPTLmaxthreads     = 51

  integer, parameter :: GPTL_IPC           = 17
  integer, parameter :: GPTL_CI            = 18
  integer, parameter :: GPTL_FPC           = 19
  integer, parameter :: GPTL_FPI           = 20
  integer, parameter :: GPTL_LSTPI         = 21
  integer, parameter :: GPTL_DCMRT         = 22
  integer, parameter :: GPTL_LSTPDCM       = 23
  integer, parameter :: GPTL_L2MRT         = 24
  integer, parameter :: GPTL_LSTPL2M       = 25
  integer, parameter :: GPTL_L3MRT         = 26

  integer, parameter :: GPTLgettimeofday   = 1
  integer, parameter :: GPTLnanotime       = 2
  integer, parameter :: GPTLmpiwtime       = 4
  integer, parameter :: GPTLclockgettime   = 5
  integer, parameter :: GPTLpapitime       = 6
  integer, parameter :: GPTLplacebo        = 7
  integer, parameter :: GPTLread_real_time = 3
					                
  integer, parameter :: GPTLfirst_parent   = 1
  integer, parameter :: GPTLlast_parent    = 2
  integer, parameter :: GPTLmost_frequent  = 3
  integer, parameter :: GPTLfull_tree      = 4

! Function prototypes

  interface
     subroutine gptlprocess_namelist (filename, unitno, outret)
       character(len=*) :: filename
       integer :: unitno
       integer :: outret
     end subroutine gptlprocess_namelist

     integer function gptlinitialize ()
     end function gptlinitialize

     integer function gptlfinalize ()
     end function gptlfinalize

     integer function gptlpr (procid)
       integer :: procid
     end function gptlpr

     integer function gptlpr_file (file)
       character(len=*) :: file
     end function gptlpr_file

#ifdef HAVE_MPI
     integer function gptlpr_summary (fcomm)
       integer :: fcomm
     end function gptlpr_summary

     integer function gptlpr_summary_file (fcomm, name)
       integer :: fcomm
       character(len=*) :: name
     end function gptlpr_summary_file

     integer function gptlbarrier (fcomm, name)
       integer :: fcomm
       character(len=*) :: name
     end function gptlbarrier
#else
     integer function gptlpr_summary ()
     end function gptlpr_summary

     integer function gptlpr_summary_file (name)
       character(len=*) :: name
     end function gptlpr_summary_file

     integer function gptlbarrier ()
     end function gptlbarrier
#endif

     integer function gptlreset ()
     end function gptlreset

     integer function gptlstamp (wall, usr, sys)
       real(8) :: wall, usr, sys
     end function gptlstamp

     integer function gptlstart (name)
       character(len=*) :: name
     end function gptlstart

     integer function gptlinit_handle (name, handle)
       character(len=*) :: name
       integer :: handle
     end function gptlinit_handle

     integer function gptlstart_handle (name, handle)
       character(len=*) :: name
       integer :: handle
     end function gptlstart_handle

     integer function gptlstop (name)
       character(len=*) :: name
     end function gptlstop

     integer function gptlstop_handle (name, handle)
       character(len=*) :: name
       integer :: handle
     end function gptlstop_handle

     integer function gptlsetoption (option, val)
       integer :: option, val
     end function gptlsetoption

     integer function gptlenable ()
     end function gptlenable

     integer function gptldisable ()
     end function gptldisable

     integer function gptlsetutr (option)
       integer :: option
     end function gptlsetutr

     integer function gptlquery (name, t, count, onflg, wallclock, &
                                 usr, sys, papicounters_out, maxcounters)
       character(len=*) :: name
       integer :: t, count
       integer :: onflg
       real(8) :: wallclock, usr, sys
       integer(8) :: papicounters_out
       integer :: maxcounters
     end function gptlquery

     integer function gptlquerycounters (name, t, papicounters_out)
       character(len=*) :: name
       integer :: t
       integer(8) :: papicounters_out
     end function gptlquerycounters

     integer function gptlget_wallclock (name, t, value)
       character(len=*) :: name
       integer :: t
       real(8) :: value
     end function gptlget_wallclock

     integer function gptlget_eventvalue (timername, eventname, t, value)
       character(len=*) :: timername
       character(len=*) :: eventname
       integer :: t
       real(8) :: value
     end function gptlget_eventvalue

     integer function gptlget_nregions (t, nregions)
       integer :: t
       integer :: nregions
     end function gptlget_nregions

     integer function gptlget_regionname (t, region, name)
       integer :: t
       integer :: region
       character(len=*) :: name
     end function gptlget_regionname

     integer function gptlget_memusage (size, rss, share, text, datastack)
       integer :: size, rss, share, text, datastack
     end function gptlget_memusage

     integer function gptlprint_memusage (str)
       character(len=*) :: str
     end function gptlprint_memusage

     integer function gptlprint_rusage (str)
       character(len=*) :: str
     end function gptlprint_rusage

     integer function gptlnum_errors ()
     end function gptlnum_errors

     integer function gptlnum_warn ()
     end function gptlnum_warn

     integer function gptlget_count (name, t, count)
       character(len=*) :: name
       integer :: t
       integer :: count
     end function gptlget_count

#ifdef HAVE_PAPI
     integer function gptl_papilibraryinit ()
     end function gptl_papilibraryinit

     integer function gptlevent_name_to_code (str, code)
       character(len=*) :: str
       integer :: code
     end function gptlevent_name_to_code

     integer function gptlevent_code_to_name (code, str)
       integer :: code
       character(len=*) :: str
     end function gptlevent_code_to_name
#endif

  end interface

  contains
! Do-nothing stub needed because some compilers otherwise generate no symbols
! which can cause ar to barf
    subroutine gptldo_nothing
    end subroutine gptldo_nothing
end module gptl
