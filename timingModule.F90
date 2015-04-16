module timingModule

!-----------------------------------------------------------------------
! timing module
!-----------------------------------------------------------------------

implicit none

include 'gptl.inc'

#if defined(HAVE_PAPI)
include 'f77papi.h'
#endif

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

private
public timing_init, timing_on, timing_off, timing_clear, timing_prt

integer :: iret
integer :: nregion
integer :: nevent

character(len=64), allocatable, dimension(:) :: regions
integer  (kind=8), allocatable, dimension(:) :: papicounters

#if defined( use_EFFICIENCY_COUNTERS )

  real :: CPU_CYCLES      
  real :: IA64_INST_RETIRED_THIS 
  real :: NOPS_RETIRED            
  real :: BACK_END_BUBBLE_ALL 
  real :: FP_OPS_RETIRED        
  real :: BUS_MEMORY_EQ_128BYTE_SELF
  real :: BUS_MEMORY_LT_128BYTE_SELF 
  
#endif


#if defined( use_STALL_COUNTERS )

  real :: BACK_END_BUBBLE_ALL
  real :: BE_EXE_BUBBLE_GRALL
  real :: BE_EXE_BUBBLE_GRGR
  real :: BE_L1D_FPU_BUBBLE_L1D
  real :: BE_EXE_BUBBLE_FRALL
  real :: BE_L1D_FPU_BUBBLE_FPU
  real :: BE_FLUSH_BUBBLE_BRU
  real :: FE_BUBBLE_BUBBLE
  real :: FE_BUBBLE_BRANCH
  real :: FE_BUBBLE_IMISS
  real :: BACK_END_BUBBLE_FE
  real :: FE_BUBBLE_ALLBUT_IBFULL
  real :: BE_RSE_BUBBLE_ALL

#endif

#if defined( use_D_CACHE_STALLS_COUNTERS )

  real :: BACK_END_BUBBLE_ALL  
  real :: BE_L1D_FPU_BUBBLE_L1D_L2BPRESS   
  real :: BE_EXE_BUBBLE_GRALL          
  real :: BE_EXE_BUBBLE_GRGR            
  real :: BE_L1D_FPU_BUBBLE_L1D_DCURECIR   
  real :: BE_L1D_FPU_BUBBLE_L1D_STBUFRECIR
  real :: BE_L1D_FPU_BUBBLE_L1D_FULLSTBUF 
  real :: BE_L1D_FPU_BUBBLE_L1D_FILLCONF   
  real :: BE_L1D_FPU_BUBBLE_L1D_TLB       
  real :: BE_L1D_FPU_BUBBLE_L1D_HPW      

#endif

contains

!####################################################################### 
 
subroutine timing_init ()
   
!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

  character(len=32), allocatable, dimension(:) :: eventname

  integer :: n, ncount, event

#if defined(HAVE_PAPI)

   iret = PAPI_VER_CURRENT
   call papif_library_init (iret)
   if (iret .ne. PAPI_VER_CURRENT) call error_handler(' PAPI_VER_CURRENT ')
   
#if defined( use_EFFICIENCY_COUNTERS )

   nevent = 7
   ncount = 0
   allocate (eventname(nevent))

!-----------------------------------------------------------------------
!  cycles
!-----------------------------------------------------------------------
   ncount = ncount + 1; eventname(ncount) = 'CPU_CYCLES'

!-----------------------------------------------------------------------
!  instruction counts
!-----------------------------------------------------------------------
   ncount = ncount + 1; eventname(ncount) = 'IA64_INST_RETIRED_THIS'
   ncount = ncount + 1; eventname(ncount) = 'NOPS_RETIRED'

!----------------------------------------------------------------------- 
!  total stalls
!-----------------------------------------------------------------------
   ncount = ncount + 1; eventname(ncount) = 'BACK_END_BUBBLE_ALL'

!----------------------------------------------------------------------- 
!  fp ops counts
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'FP_OPS_RETIRED'

!----------------------------------------------------------------------- 
!  main memory bandwidth
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BUS_MEMORY_EQ_128BYTE_SELF'
   ncount = ncount + 1; eventname(ncount) = 'BUS_MEMORY_LT_128BYTE_SELF'

#endif


#if defined( use_STALL_COUNTERS )

   nevent = 13
   ncount =  0
   allocate (eventname(nevent))

!----------------------------------------------------------------------- 
!  total stalls
!-----------------------------------------------------------------------
   ncount = ncount + 1; eventname(ncount) = 'BACK_END_BUBBLE_ALL'

!----------------------------------------------------------------------- 
!  d-cache stalls
!-----------------------------------------------------------------------
   ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_GRALL'
   ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_GRGR'
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D'

!----------------------------------------------------------------------- 
!  fpu stalls 
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_FRALL'
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_FPU'
   
!----------------------------------------------------------------------- 
!  branch mispredict stalls
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_FLUSH_BUBBLE_BRU'
   ncount = ncount + 1; eventname(ncount) = 'FE_BUBBLE_BUBBLE'
   ncount = ncount + 1; eventname(ncount) = 'FE_BUBBLE_BRANCH'

!----------------------------------------------------------------------- 
!  i-cache stalls 
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'FE_BUBBLE_IMISS'
   ncount = ncount + 1; eventname(ncount) = 'BACK_END_BUBBLE_FE'
   ncount = ncount + 1; eventname(ncount) = 'FE_BUBBLE_ALLBUT_IBFULL'

!----------------------------------------------------------------------- 
!  rse stalls
!----------------------------------------------------------------------- 
!  ncount = ncount + 1; eventname(ncount) = 'BACK_END_BUBBLE_FE'
   ncount = ncount + 1; eventname(ncount) = 'BE_RSE_BUBBLE_ALL'

!----------------------------------------------------------------------- 
!  support register dependency stalls
!----------------------------------------------------------------------- 
!  ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_ARCR_PR_CANCEL_BANK'

!----------------------------------------------------------------------- 
!  integer register dependency stalls
!----------------------------------------------------------------------- 
!  ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_GRGR'
 
#endif

#if defined( use_D_CACHE_STALLS_COUNTERS )

   nevent = 10
   ncount = 0
   allocate (eventname(nevent))

!----------------------------------------------------------------------- 
!  total stalls
!----------------------------------------------------------------------- 

   ncount = ncount + 1; eventname(ncount) = 'BACK_END_BUBBLE_ALL'

!----------------------------------------------------------------------- 
!  l2 Capacity Stalls
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_L2BPRESS'

!----------------------------------------------------------------------- 
!  integer load latency stalls
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_GRALL'
   ncount = ncount + 1; eventname(ncount) = 'BE_EXE_BUBBLE_GRGR'

!----------------------------------------------------------------------- 
!  l2 recirculation stalls
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_DCURECIR'
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_STBUFRECIR'

!----------------------------------------------------------------------- 
!  store related stalls
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_FULLSTBUF'
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_FILLCONF'

!----------------------------------------------------------------------- 
!  virtual memory stalls
!----------------------------------------------------------------------- 
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_TLB'
   ncount = ncount + 1; eventname(ncount) = 'BE_L1D_FPU_BUBBLE_L1D_HPW'
   
#endif
 
   do n=1,nevent
     call PAPIf_event_name_to_code(eventname(n), event, iret)
     if (iret .ne. 0)                     call error_handler ( ' PAPIf_event_name_to_code ' )
     if (gptlsetoption (event, 1) .lt. 0) call error_handler ( ' gptlsetoption ' )
   enddo

   if (gptlsetoption  (gptloverhead   , 0) < 0) call error_handler ( ' gptloverhead   ' )
   if (gptlsetoption  (gptlnarrowprint, 1) < 0) call error_handler ( ' gptlsetoption  ' )
   if (gptlinitialize ()                   < 0) call error_handler ( ' gptlinitialize ' )

#endif

end subroutine timing_init

!####################################################################### 

subroutine timing_on ( blk_name )

!----------------------------------------------------------------------- 
!  timing on
!----------------------------------------------------------------------- 

   character*(*), intent(in) :: blk_name
 
   iret = gptlstart( blk_name )

end subroutine timing_on

!#######################################################################

subroutine timing_off ( blk_name )
  
!----------------------------------------------------------------------- 
!  timing off
!----------------------------------------------------------------------- 

  character*(*), intent(in) :: blk_name

  iret = gptlstop ( blk_name )

end subroutine timing_off

!#######################################################################

subroutine timing_clear ()
  
!----------------------------------------------------------------------- 
!  timing clear
!----------------------------------------------------------------------- 
     
  iret = gptlreset()

end subroutine timing_clear

!#######################################################################

subroutine timing_prt

!-----------------------------------------------------------------------
! timing print
!-----------------------------------------------------------------------

  integer n

  iret = gptlpr ( 0 )

  nregion = 1
  allocate (regions(nregion))

  regions(1) = 'FYPPM4'

  allocate (papicounters(nevent))

  do n=1,nregion
    iret = gptlquerycounters (TRIM(regions(n)), -1, papicounters)
  
    write(6,*)
    write(6,"(' REGION NAME = ', a32)") regions(n)

#if defined( use_EFFICIENCY_COUNTERS )

    CPU_CYCLES                 = papicounters(1)
    IA64_INST_RETIRED_THIS     = papicounters(2)
    NOPS_RETIRED               = papicounters(3)
    BACK_END_BUBBLE_ALL        = papicounters(4)
    FP_OPS_RETIRED             = papicounters(5)
    BUS_MEMORY_EQ_128BYTE_SELF = papicounters(6)
    BUS_MEMORY_LT_128BYTE_SELF = papicounters(7)

    write(6,*)
    write(6,"(' CPU_CYCLES...........................= ', e12.6)") CPU_CYCLES
    write(6,"(' IA64_INST_RETIRED_THIS...............= ', e12.6)") IA64_INST_RETIRED_THIS
    write(6,"(' NOPS_RETIRED.........................= ', e12.6)") NOPS_RETIRED
    write(6,"(' BACK_END_BUBBLE_ALL..................= ', e12.6)") BACK_END_BUBBLE_ALL
    write(6,"(' FP_OPS_RETIRED.......................= ', e12.6)") FP_OPS_RETIRED
    write(6,"(' BUS_MEMORY_EQ_128BYTE_SELF...........= ', e12.6)") BUS_MEMORY_EQ_128BYTE_SELF
    write(6,"(' BUS_MEMORY_LT_128BYTE_SELF...........= ', e12.6)") BUS_MEMORY_LT_128BYTE_SELF
    write(6,*)
    write(6,"(' Useful Ops/Cycle.........................= ', e9.3)") ( IA64_INST_RETIRED_THIS -  NOPS_RETIRED )/CPU_CYCLES
    write(6,"(' NOPS/Cycle ..............................= ', e9.3)")   NOPS_RETIRED /CPU_CYCLES
    write(6,"(' Total Stalls/Cycle.......................= ', e9.3)")   BACK_END_BUBBLE_ALL / CPU_CYCLES
    write(6,"(' FLOPS/Cycle..............................= ', e9.3)")   FP_OPS_RETIRED / CPU_CYCLES
    write(6,"(' Main Memory Bandwidth Used...............= ', e9.3)") ( BUS_MEMORY_EQ_128BYTE_SELF*128 +BUS_MEMORY_LT_128BYTE_SELF*128)/CPU_CYCLES
#endif

#if defined( use_STALL_COUNTERS )

    BACK_END_BUBBLE_ALL     = papicounters(1)
    BE_EXE_BUBBLE_GRALL     = papicounters(2)
    BE_EXE_BUBBLE_GRGR      = papicounters(3)
    BE_L1D_FPU_BUBBLE_L1D   = papicounters(4)
    BE_EXE_BUBBLE_FRALL     = papicounters(5)
    BE_L1D_FPU_BUBBLE_FPU   = papicounters(6)
    BE_FLUSH_BUBBLE_BRU     = papicounters(7)
    FE_BUBBLE_BUBBLE        = papicounters(8)
    FE_BUBBLE_BRANCH        = papicounters(9)
    FE_BUBBLE_IMISS         = papicounters(10)
    BACK_END_BUBBLE_FE      = papicounters(11)
    FE_BUBBLE_ALLBUT_IBFULL = papicounters(12)
    BE_RSE_BUBBLE_ALL       = papicounters(13)

    write(6,"(' BACK_END_BUBBLE_ALL..................= ', e12.6)") BACK_END_BUBBLE_ALL
    write(6,"(' BE_EXE_BUBBLE_GRALL..................= ', e12.6)") BE_EXE_BUBBLE_GRALL
    write(6,"(' BE_EXE_BUBBLE_GRGR...................= ', e12.6)") BE_EXE_BUBBLE_GRGR
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D................= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D
    write(6,"(' BE_EXE_BUBBLE_FRALL..................= ', e12.6)") BE_EXE_BUBBLE_FRALL
    write(6,"(' BE_L1D_FPU_BUBBLE_FPU................= ', e12.6)") BE_L1D_FPU_BUBBLE_FPU
    write(6,"(' BE_FLUSH_BUBBLE_BRU..................= ', e12.6)") BE_FLUSH_BUBBLE_BRU
    write(6,"(' FE_BUBBLE_BUBBLE.....................= ', e12.6)") FE_BUBBLE_BUBBLE
    write(6,"(' FE_BUBBLE_BRANCH.....................= ', e12.6)") FE_BUBBLE_BUBBLE
    write(6,"(' FE_BUBBLE_IMISS......................= ', e12.6)") FE_BUBBLE_IMISS
    write(6,"(' BACK_END_BUBBLE_FE...................= ', e12.6)") BACK_END_BUBBLE_FE
    write(6,"(' FE_BUBBLE_ALLBUT_IBFULL..............= ', e12.6)") FE_BUBBLE_ALLBUT_IBFULL
    write(6,"(' BE_RSE_BUBBLE_ALL....................= ', e12.6)") BE_RSE_BUBBLE_ALL
   
    write(6,*)
    write(6,"(' D-Cache Stalls...........................= ', e9.3)") ( BE_EXE_BUBBLE_GRALL - BE_EXE_BUBBLE_GRGR + BE_L1D_FPU_BUBBLE_L1D ) / BACK_END_BUBBLE_ALL
    write(6,"(' Branch Misprediction Stalls..............= ', e9.3)") ( BE_FLUSH_BUBBLE_BRU + ( FE_BUBBLE_BUBBLE + FE_BUBBLE_BRANCH ) * &
                                                                      ( BACK_END_BUBBLE_FE / FE_BUBBLE_ALLBUT_IBFULL ) ) / BACK_END_BUBBLE_ALL
    write(6,"(' I-Cache Stalls...........................= ', e9.3)") ( FE_BUBBLE_IMISS ) * ( BACK_END_BUBBLE_FE / FE_BUBBLE_ALLBUT_IBFULL ) / BACK_END_BUBBLE_ALL
    write(6,"(' FPU Stalls...............................= ', e9.3)") ( BE_EXE_BUBBLE_FRALL + BE_L1D_FPU_BUBBLE_FPU ) / BACK_END_BUBBLE_ALL
    write(6,"(' RSE Stalls...............................= ', e9.3)")   BE_RSE_BUBBLE_ALL / BACK_END_BUBBLE_ALL
    write(6,"(' Integer Register Dependency Stalls.......= ', e9.3)")   BE_EXE_BUBBLE_GRGR / BACK_END_BUBBLE_ALL
    write(6,"(' Support Register Dependency Stalls.......= ', a2  )")  'na'
#endif

#if defined( use_D_CACHE_STALLS_COUNTERS )

    BACK_END_BUBBLE_ALL              = papicounters(1)
    BE_L1D_FPU_BUBBLE_L1D_L2BPRESS   = papicounters(2)
    BE_EXE_BUBBLE_GRALL              = papicounters(3)
    BE_EXE_BUBBLE_GRGR               = papicounters(4)
    BE_L1D_FPU_BUBBLE_L1D_DCURECIR   = papicounters(5)
    BE_L1D_FPU_BUBBLE_L1D_STBUFRECIR = papicounters(6)
    BE_L1D_FPU_BUBBLE_L1D_FULLSTBUF  = papicounters(7)
    BE_L1D_FPU_BUBBLE_L1D_FILLCONF   = papicounters(8)
    BE_L1D_FPU_BUBBLE_L1D_TLB        = papicounters(9)
    BE_L1D_FPU_BUBBLE_L1D_HPW        = papicounters(10)

    write(6,"(' BACK_END_BUBBLE_ALL..................= ', e12.6)") BACK_END_BUBBLE_ALL
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_L2BPRESS.......= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_L2BPRESS
    write(6,"(' BE_EXE_BUBBLE_GRALL..................= ', e12.6)") BE_EXE_BUBBLE_GRALL
    write(6,"(' BE_EXE_BUBBLE_GRGR...................= ', e12.6)") BE_EXE_BUBBLE_GRGR
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_DCURECIR.......= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_DCURECIR
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_STBUFRECIR.....= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_STBUFRECIR
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_FULLSTBUF......= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_FULLSTBUF
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_FILLCONF.......= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_FILLCONF
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_TLB............= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_TLB
    write(6,"(' BE_L1D_FPU_BUBBLE_L1D_HPW............= ', e12.6)") BE_L1D_FPU_BUBBLE_L1D_HPW

    write(6,*)   
    write(6,"(' L2 Capacity Stalls.......................= ', e9.3)")  BE_L1D_FPU_BUBBLE_L1D_L2BPRESS / BACK_END_BUBBLE_ALL 
    write(6,"(' Integer Load Latency Stalls..............= ', e9.3)") ( BE_EXE_BUBBLE_GRALL - BE_EXE_BUBBLE_GRGR ) / BACK_END_BUBBLE_ALL
    write(6,"(' L2 Recirculation Stalls..................= ', e9.3)") ( BE_L1D_FPU_BUBBLE_L1D_DCURECIR + BE_L1D_FPU_BUBBLE_L1D_STBUFRECIR ) / BACK_END_BUBBLE_ALL
    write(6,"(' Store Related Stalls.....................= ', e9.3)") ( BE_L1D_FPU_BUBBLE_L1D_FULLSTBUF + BE_L1D_FPU_BUBBLE_L1D_FILLCONF ) / BACK_END_BUBBLE_ALL
    write(6,"(' Virtual Memory Stalls....................= ', e9.3)") ( BE_L1D_FPU_BUBBLE_L1D_TLB + BE_L1D_FPU_BUBBLE_L1D_HPW ) / BACK_END_BUBBLE_ALL
#endif
  
  enddo

  iret = gptlfinalize ()

  return
  
end subroutine timing_prt

!#######################################################################

subroutine error_handler ( message ) 
character(len=*), intent(in) :: message
   
  print *, message
  stop

  return
end subroutine error_handler

!#######################################################################

end module TimingModule
