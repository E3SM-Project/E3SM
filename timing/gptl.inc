!
! $Id: gptl.inc,v 1.41 2009/12/29 02:03:53 rosinski Exp $
!
! Author: Jim Rosinski
!
! #include or Fortran include this file in application codes which invoke GPTL functions.
! Values defined here match their counterparts in gptl.h See that file or man pages
! or web-based documenation for descriptions of each value
!
      integer GPTLsync_mpi
      integer GPTLwall
      integer GPTLcpu
      integer GPTLabort_on_error
      integer GPTLoverhead
      integer GPTLdepthlimit
      integer GPTLverbose
      integer GPTLnarrowprint
      integer GPTLpercent
      integer GPTLpersec
      integer GPTLmultiplex
      integer GPTLdopr_preamble
      integer GPTLdopr_threadsort
      integer GPTLdopr_multparent
      integer GPTLdopr_collision
      integer GPTLprint_method

      integer GPTL_IPC
      integer GPTL_CI
      integer GPTL_FPC
      integer GPTL_FPI
      integer GPTL_LSTPI
      integer GPTL_DCMRT
      integer GPTL_LSTPDCM
      integer GPTL_L2MRT
      integer GPTL_LSTPL2M
      integer GPTL_L3MRT

      integer GPTLnanotime
      integer GPTLrtc
      integer GPTLmpiwtime
      integer GPTLclockgettime
      integer GPTLgettimeofday
      integer GPTLpapitime

      integer GPTLfirst_parent
      integer GPTLlast_parent
      integer GPTLmost_frequent
      integer GPTLfull_tree

      parameter (GPTLsync_mpi       = 0)
      parameter (GPTLwall           = 1)
      parameter (GPTLcpu            = 2)
      parameter (GPTLabort_on_error = 3)
      parameter (GPTLoverhead       = 4)
      parameter (GPTLdepthlimit     = 5)
      parameter (GPTLverbose        = 6)
      parameter (GPTLnarrowprint    = 7)
      parameter (GPTLpercent        = 9)
      parameter (GPTLpersec         = 10)
      parameter (GPTLmultiplex      = 11)
      parameter (GPTLdopr_preamble  = 12)
      parameter (GPTLdopr_threadsort= 13)
      parameter (GPTLdopr_multparent= 14)
      parameter (GPTLdopr_collision = 15)
      parameter (GPTLprint_method   = 16)

      parameter (GPTL_IPC           = 17)
      parameter (GPTL_CI            = 18)
      parameter (GPTL_FPC           = 19)
      parameter (GPTL_FPI           = 20)
      parameter (GPTL_LSTPI         = 21)
      parameter (GPTL_DCMRT         = 22)
      parameter (GPTL_LSTPDCM       = 23)
      parameter (GPTL_L2MRT         = 24)
      parameter (GPTL_LSTPL2M       = 25)
      parameter (GPTL_L3MRT         = 26)

      parameter (GPTLgettimeofday   = 1)
      parameter (GPTLnanotime       = 2)
      parameter (GPTLrtc            = 3)
      parameter (GPTLmpiwtime       = 4)
      parameter (GPTLclockgettime   = 5)
      parameter (GPTLpapitime       = 6)

      parameter (GPTLfirst_parent   = 1)
      parameter (GPTLlast_parent    = 2)
      parameter (GPTLmost_frequent  = 3)
      parameter (GPTLfull_tree      = 4)

! Externals

      integer gptlsetoption
      integer gptlinitialize
      integer gptlstart
      integer gptlstop
      integer gptlstamp 
      integer gptlpr
      integer gptlpr_file
      integer gptlpr_summary
      integer gptlbarrier
      integer gptlreset 
      integer gptlfinalize
      integer gptlget_memusage
      integer gptlprint_memusage
      integer gptlenable
      integer gptldisable
      integer gptlsetutr
      integer gptlquery
      integer gptlquerycounters
      integer gptlget_wallclock
      integer gptlget_eventvalue
      integer gptlget_nregions
      integer gptlget_regionname
      integer gptl_papilibraryinit
      integer gptlevent_name_to_code
      integer gptlevent_code_to_name

      external gptlsetoption
      external gptlinitialize
      external gptlstart
      external gptlstop
      external gptlstamp 
      external gptlpr
      external gptlpr_file
      external gptlpr_summary
      external gptlbarrier
      external gptlreset 
      external gptlfinalize
      external gptlget_memusage
      external gptlprint_memusage
      external gptlenable
      external gptldisable
      external gptlsetutr
      external gptlquery
      external gptlquerycounters
      external gptlget_wallclock
      external gptlget_eventvalue
      external gptlget_nregions
      external gptlget_regionname
      external gptl_papilibraryinit
      external gptlevent_name_to_code
      external gptlevent_code_to_name
