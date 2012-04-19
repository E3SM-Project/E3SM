#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module perfmodel_mod


	use kinds , only : real_kind, log_kind, int_kind

implicit none
private 

	type,  public :: perf_t
	   sequence
	   real(kind=real_kind)   :: onnode
	   real(kind=real_kind)   :: offnode
        end type perf_t 

	type, public :: commtime_t
	   sequence
	   real(kind=real_kind)  :: latency
	   real(kind=real_kind)  :: bandwidth
	end type commtime_t	

	public :: SetPtoPNetworkParams
        public :: SetReduceNetworkParams
	public :: SetSerialParamsExplicit
        public :: SetSerialParamsPrimExplicit
	public :: SetSerialParamsImplicit

contains

	subroutine SetPtoPNetworkParams(offnode,onnode,network,found)

	type (commtime_t),intent(inout)     :: offnode
	type (commtime_t),intent(inout)     :: onnode
        character(len=80)                   :: network
        logical(kind=log_kind),intent(inout)  :: found
        real (kind=real_kind) :: contention

	! Macine specific characteristics
	! bluesky Power4 with dual rail 8-way LPAR 
        ! bandwidth: microseconds per byte
        ! latency:   microseconds

	found=.FALSE.
        if(network(1:7) == "bluesky") then 
           offnode%latency   = 17.9 
           onnode%latency    = .9   
           offnode%bandwidth = 0.002777
           !offnode%bandwidth = 3.*0.002777
           onnode%bandwidth  = 0.001203
           found=.TRUE.
        elseif (network(1:11) == "blackforest") then 
           offnode%latency   = 30.0 
           onnode%latency    = 2.   
           offnode%bandwidth = 0.0125
           onnode%bandwidth  = 0.00416
           found=.TRUE.
        elseif (network(1:8) == "protoBGL") then 
	   ! This is the from the description of the prototype BGL from Oct 2003 meeting
           offnode%latency   = 6. 
           onnode%latency    = .1   
           !contention = 3.0
           !contention = 1.0
           contention =  9.0
           offnode%bandwidth = contention/420.! Note the value of 420 is 3 times what I get on point-to-point
					      ! messaging... I use the three because they only get 3 receive 
					      ! streams going at once.  The other three is due to the fact 
					      ! that the corner processors only have 3 network connections
					      ! to simulate contention 8 neighbors/6 links * 2.25 = 3
           onnode%bandwidth  =1.0/(3.*420.) 
           found=.TRUE.
        elseif (network(1:8) == "BGL-fast") then 
	   ! BGL with network as described in the Paper "???"
           offnode%latency   = 1.5 
           onnode%latency    = .1   
           !contention = 2.25
           contention = 4.5
           offnode%bandwidth = contention/1050.    !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           onnode%bandwidth  =1.0/(3.*1050.) 
           found=.TRUE.
        elseif (network(1:8) == "BGL-slow") then 
	   ! BGL with network running at 1/2 of described in the Paper "???"
           offnode%latency   = 4.5 
           onnode%latency    = .1   
           !contention = 2.25
           contention = 4.5
           offnode%bandwidth = contention/525.    !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           onnode%bandwidth  =1.0/(3.*525.) 
           found=.TRUE.
        else 
	   print *,"SetPtoPNetworkParams:  Undefined network name"
	   print *,"SetPtoPNetworkParams:  Please select from the following or create your own:"
           print *,"SetPtoPNetworkParams:  bluesky     == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetPtoPNetworkParams:  blackforest == Power3 with TBMX network"
           print *,"SetPtoPNetworkParams:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetPtoPNetworkParams:  BGL-fast    == BGL network as described in Paper ??? "
           print *,"SetPtoPNetworkParams:  BGL-slow    == 1/2 perf. of BGL network as described in Paper ???? "
        endif

	end subroutine SetPtoPNetworkParams

	subroutine SetReduceNetworkParams(latency,bandwidth,network,found)

	type (perf_t),intent(inout)         :: latency
	type (perf_t),intent(inout)         :: bandwidth
        character(len=80)                   :: network
        logical(kind=log_kind),intent(inout)  :: found

	! Macine specific characteristics
	! bluesky Power4 with dual rail 8-way LPAR 
        ! bandwidth: microseconds per byte
        ! latency:   microseconds

	found=.FALSE.
        if(network(1:7) == "bluesky") then 
           latency%offnode   = 17.9 
           latency%onnode    = .9   
           bandwidth%offnode = 0.002777
           bandwidth%onnode  = 0.001203
           found=.TRUE.
        elseif (network(1:8) == "protoBGL") then 
	   ! This is the from the description of the prototype BGL from Oct 2003 meeting
           latency%offnode   = 3. 
           latency%onnode    = .1   
           bandwidth%offnode = 1./350.    ! Add 3 here because the corner processors only 
					  ! have three network connections (just a 3-D mesh) at the moment
           bandwidth%onnode  =1.0/(3.*420.) 
           found=.TRUE.
        elseif (network(1:8) == "BGL-fast") then 
	   ! BGL with network as described in the Paper "???"
           latency%offnode   = 1.5 
           latency%onnode    = .1   
           bandwidth%offnode = 1/350.    !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           bandwidth%onnode  =1.0/(3.*350.) 
           found=.TRUE.
        elseif (network(1:8) == "BGL-slow") then 
	   ! BGL with network running at 1/2 of described in the Paper "???"
           latency%offnode   = 1.5 
           latency%onnode    = .1   
           bandwidth%offnode = 1.0/350.    !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           bandwidth%onnode  =1.0/(3.*350.) 
           found=.TRUE.
        else 
	   print *,"SetReduceNetworkParams:  Undefined network name"
	   print *,"SetReduceNetworkParams:  Please select from the following or create your own:"
           print *,"SetReduceNetworkParams:  bluesky    == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetReduceNetworkParams:  protoBGL   == Prototype BGL machine (November 2003 configuration)"
           print *,"SetReduceNetworkParams:  BGL-fast   == BGL network as described in Paper ??? "
           print *,"SetReduceNetworkParams:  BGL-slow   == 1/2 perf. of BGL network as described in Paper ???? "
        endif

	end subroutine SetReduceNetworkParams

	subroutine SetSerialParamsExplicit(microsec,np,nlev,machine,found)

        real (kind=real_kind),intent(out)   :: microsec
        integer(kind=int_kind),intent(in)   :: np
        integer(kind=int_kind),intent(in)   :: nlev
        character(len=80)                   :: machine
        logical(kind=log_kind),intent(inout)  :: found
        logical(kind=log_kind)                :: m_found,np_found,nlev_found

	m_found    = .FALSE.
	np_found   = .FALSE.
	nlev_found = .FALSE.
        if(machine(1:7) == "bluesky") then 
	   m_found = .TRUE. 
           if(np .eq. 6 ) then
	     np_found=.TRUE.
             select case (nlev)
                case(96)
!                   microsec         = 4243.D0 ! SMP numbers (nproc=8)
                    microsec         = 3299.D0 ! single threaded
		    nlev_found=.TRUE.
                case(60)
!                 microsec         = 2545.D0 ! SMP number (nproc=8)
!                 microsec         = 2194.D0 ! SMP number (nproc=2)
                  microsec         = 2013.D0 ! single threaded
		  nlev_found=.TRUE.
                case(48)
!                 microsec         = 1914.D0 ! SMP number (nproc=8)
                  microsec         = 1702.D0 ! single threaded
		  nlev_found=.TRUE.
                case(30)
!                 microsec         = 1142.D0 ! SMP number (nproc=8)
                  microsec         = 1023.D0 ! single threaded
		  nlev_found=.TRUE.
                case(16)
                  microsec         =  540.D0
		  nlev_found=.TRUE.
                case(2)
                  microsec         =   76.D0
		  nlev_found=.TRUE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select   	
          elseif( np .eq. 8) then  
	     np_found=.TRUE.
             select case (nlev)
                case(16)
                  microsec         =  556.D0
		  nlev_found=.TRUE.
                case(2)
                  microsec         =   76.D0
		  nlev_found=.FALSE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select   	
	  elseif (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  microsec        = 9570.D0 ! single threaded
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	     print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np 
          endif
        elseif(machine(1:11) == "blackforest") then 
	   m_found = .TRUE. 
           if(np .eq. 8 ) then
	     np_found=.TRUE.
             select case (nlev)
                case(16)
                  microsec         =  1652.D0  ! single threaded
		  nlev_found=.TRUE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
	     end select   	
	  elseif (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  !microsec        = 9570.D0 ! Please Fix Me
		  nlev_found = .FALSE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	     print *,"SetSerialParamsExplicit:  machine: ",machine(1:11)," np := ",np 
          endif
        elseif(machine(1:8) == "protoBGL") then 
	   m_found = .TRUE. 
          if (np .eq. 8) then 
	     np_found = .TRUE.
             select case (nlev)
                case(16)
		  microsec        = 2087.D0 ! single threaded
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
	  elseif( np .eq. 10 ) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  microsec        = (900./210)*9570.D0 ! This is an estimate
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "	
	     print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np 
          endif
        elseif(machine(1:8) == "BGL-fast") then 
	   m_found = .TRUE. 
	   if (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  ! Note: Bluesky runs the code at 1140 Mflops in 5761.7
    		  !       Estimate BGL-fast runs at about 750 Mflops
		  !       The 
	 	  microsec   = (1140.0/750.0)* 5761.7 
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParams:  Could not find performance estimate for: "	
	          print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
           else
	     print *,"SetSerialParams:  Could not find performance estimate for: "	
	     print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np 
           endif
        elseif(machine(1:8) == "BGL-slow") then 
	   m_found = .TRUE. 
	   if (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  ! Note: Bluesky runs the code at 1140 Mflops in 5761.7
    		  !       Estimate BGL-slow runs at about 450 Mflops
		  !       The 
	 	  microsec   = (1140.0/450.0)* 5761.7 
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParams:  Could not find performance estimate for: "	
	          print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	      end select
            else
	      print *,"SetSerialParams:  Could not find performance estimate for: "	
	      print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np 
	    endif
         else
	   print *,"SetSerialParamsExplicit:   Undefined machine name"
	   print *,"SetSerialParamsExplicit:  Please select from the following or create your own:"
           print *,"SetSerialParamsExplicit:  bluesky     == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetSerialParamsExplicit:  blackforest == Power3 with TBMX network"
           print *,"SetSerialParamsExplicit:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetSerialParamsExplicit:  BGL-fast    == BGL processor clock scaled based on Power3 "
           print *,"SetSerialParamsExplicit:  BGL-slow    == BGL processor clock scaled based on Power4 "
	   return
         endif
	

	if(m_found) then 
	  if(np_found) then
	    if(nlev_found) found=.TRUE.
	  endif
	endif
	

	end subroutine SetSerialParamsExplicit

	subroutine SetSerialParamsPrimExplicit(microsec,np,nlev,machine,found)

        real (kind=real_kind),intent(out)   :: microsec
        integer(kind=int_kind),intent(in)   :: np
        integer(kind=int_kind),intent(in)   :: nlev
        character(len=80)                   :: machine
        logical(kind=log_kind),intent(inout)  :: found
        logical(kind=log_kind)                :: m_found,np_found,nlev_found

	m_found    = .FALSE.
	np_found   = .FALSE.
	nlev_found = .FALSE.
        if(machine(1:8) == "protoBGL") then 
	   m_found = .TRUE. 
          if (np .eq. 8) then 
	     np_found = .TRUE.
             select case (nlev)
                case(26)
		  microsec        = 10706.8 ! single threaded
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsPrimExplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsPrimExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsPrimExplicit:  Could not find performance estimate for: "	
	     print *,"SetSerialParamsPrimExplicit:  machine: ",machine(1:8)," np := ",np 
          endif
        elseif(machine(1:8) == "BGL-slow") then 
	   m_found = .TRUE. 
	   if (np .eq. 8) then 
	     np_found = .TRUE.
             select case (nlev)
                case(26)
		  ! Assume it scales like clock speed from the prototype BGL 
	 	  microsec   = (500./750)* 10706.8 
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialPrimParams:  Could not find performance estimate for: "	
	          print *,"SetSerialPrimParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	      end select
            else
	      print *,"SetSerialPrimParams:  Could not find performance estimate for: "	
	      print *,"SetSerialPrimParams:  machine: ",machine(1:8)," np := ",np 
	    endif
         else
	   print *,"SetSerialParamsPrimExplicit:   Undefined machine name"
	   print *,"SetSerialParamsPrimExplicit:  Please select from the following or create your own:"
           print *,"SetSerialParamsPrimExplicit:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetSerialParamsPrimExplicit:  BGL-fast    == BGL processor clock scaled based on Power3 "
           print *,"SetSerialParamsPrimExplicit:  BGL-slow    == BGL processor clock scaled based on Power4 "
	   return
         endif
	

	if(m_found) then 
	  if(np_found) then
	    if(nlev_found) found=.TRUE.
	  endif
	endif
	



	end subroutine SetSerialParamsPrimExplicit

	subroutine SetSerialParamsImplicit(microsec,microsec_per_iter,cg_iters,np,nelem,nlev,machine,found)

        real (kind=real_kind),intent(out)   :: microsec
        real (kind=real_kind),intent(out)   :: microsec_per_iter
        real (kind=real_kind),intent(out)   :: cg_iters
        integer(kind=int_kind),intent(in)   :: np
        integer(kind=int_kind),intent(in)   :: nelem
        integer(kind=int_kind),intent(in)   :: nlev
        character(len=80)                   :: machine
        logical(kind=log_kind),intent(inout)  :: found
        logical(kind=log_kind)                :: m_found,np_found,nlev_found,nelem_found

        print *,'SetSerialParamsImplicit:  Note all values of broken... please fix'
	m_found    = .FALSE.
	np_found   = .FALSE.
        nelem_found   = .FALSE.
	nlev_found = .FALSE.
        if(machine(1:7) == "bluesky") then 
	   m_found = .TRUE. 
           if(np .eq. 6 ) then
	     np_found=.TRUE.
             select case (nlev)
                case(96)
!                   microsec         = 4243.D0 ! SMP numbers (nproc=8)
                    microsec         = 3299.D0 ! single threaded
		    nlev_found=.TRUE.
                case(60)
!                 microsec         = 2545.D0 ! SMP number (nproc=8)
!                 microsec         = 2194.D0 ! SMP number (nproc=2)
                  microsec         = 2013.D0 ! single threaded
		  nlev_found=.TRUE.
                case(48)
!                 microsec         = 1914.D0 ! SMP number (nproc=8)
                  microsec         = 1702.D0 ! single threaded
		  nlev_found=.TRUE.
                case(30)
!                 microsec         = 1142.D0 ! SMP number (nproc=8)
                  microsec         = 1023.D0 ! single threaded
		  nlev_found=.TRUE.
                case(16)
                  microsec         =  540.D0
		  nlev_found=.TRUE.
                case(2)
                  microsec         =   76.D0
		  nlev_found=.TRUE.
		case default
	          print *,"SetSerialParamsImplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select   	
	  elseif (np .eq. 8) then 
	     np_found = .TRUE.
             select case (nlev)
                case(16)
		  microsec          = 800.
                  microsec_per_iter = 450.
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
	     end select
	  elseif (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(16)
		  microsec          = 1303.
                  microsec_per_iter =  793.
		  nlev_found = .TRUE.
                case(96)
		  microsec        = 9570.D0 ! single threaded
		  nlev_found = .FALSE.
		case default
	          print *,"SetSerialParamsImplicit:  Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsImplicit:  Could not find performance estimate for: "	
	     print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np 
          endif
        elseif(machine(1:11) == "blackforest") then 
	   m_found = .TRUE. 
          if (np .eq. 8) then 
	     np_found = .TRUE.
             select case (nlev)
                case(16)
		  microsec          = 2077.
                  microsec_per_iter = 1252.
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
	     end select
	  elseif( np .eq. 10 ) then 
	     np_found = .TRUE.
             select case (nlev)
                case(16)
		  microsec          = 3200.
                  microsec_per_iter = 2320.
		  nlev_found = .TRUE.
                case(96)
		  nlev_found = .FALSE.
		case default
	          print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	     print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np 
          endif
        elseif(machine(1:8) == "protoBGL") then 
	   m_found = .TRUE. 
          if (np .eq. 8) then 
	     np_found = .TRUE.
             select case (nlev)
                case(16)
		  microsec        = 9187.4D0 ! single threaded
		  cg_iters = 3.5528
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
	  elseif( np .eq. 10 ) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  microsec        = (900./210)*9570.D0 ! This is an estimate
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	          print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
          else
	     print *,"SetSerialParamsImplicit: Could not find performance estimate for: "	
	     print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np 
          endif
        elseif(machine(1:8) == "BGL-fast") then 
	   m_found = .TRUE. 
	   if (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  ! This is an estimate 
	 	  microsec   = (375./700.) * 24218.D0 ! single threaded
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParams: Could not find performance estimate for: "	
	          print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	     end select
           else
	     print *,"SetSerialParams: Could not find performance estimate for: "	
	     print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np 
           endif
        elseif(machine(1:8) == "BGL-slow") then 
	   m_found = .TRUE. 
	   if (np .eq. 10) then 
	     np_found = .TRUE.
             select case (nlev)
                case(96)
		  ! This is an estimate  
		  microsec        = (1300./700.) * 9570.D0 
		  nlev_found = .TRUE.
		case default
	          print *,"SetSerialParams: Could not find performance estimate for: "	
	          print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
	      end select
            else
	      print *,"SetSerialParams: Could not find performance estimate for: "	
	      print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np 
	    endif
         else
	   print *,"SetSerialParams:  Undefined machine name"
	   print *,"SetSerialParams:  Please select from the following or create your own:"
           print *,"SetSerialParams:  bluesky     == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetSerialParams:  blackforest == Power3 with TBMX network"
           print *,"SetSerialParams:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetSerialParams:  BGL-fast    == BGL processor clock scaled based on Power3 "
           print *,"SetSerialParams:  BGL-slow    == BGL processor clock scaled based on Power4 "
	   return
         endif

         ! =====================================
         !   Set the expected iteration count 
         ! =====================================
         select case(np)
             case(8)
                select case (nelem)
                  case(24) 
		     cg_iters = 2.1042
		     nelem_found = .TRUE.
                  case(54)
                     cg_iters = 2.9722
		     nelem_found = .TRUE.
                  case(96)
		     cg_iters = 3.2778
		     nelem_found = .TRUE.
		  case(384)
		     cg_iters = 3.5556
		     nelem_found = .TRUE.
		  case(864)
		     cg_iters = 3.9708
		     nelem_found = .TRUE.
		  case default 
		     print *,"SetSerialParamsImplicit: Could not iteration count estimate for: "
		     print *,'SetSerialParamsImplicit: np := ',np, 'nelem := ',nelem
		     nelem_found = .FALSE.
                  end select
             case(10)
                select case (nelem)
                  case(24) 
		     cg_iters = 2.9583
		     nelem_found = .TRUE.
                  case(54)
                     cg_iters = 3.4427
		     nelem_found = .TRUE.
                  case(96)
		     cg_iters = 3.6958
		     nelem_found = .TRUE.
		  case(384)
		     cg_iters = 3.7979
		     nelem_found = .TRUE.
                  case(1536) 
		     cg_iters = 3.9708
		     nelem_found = .TRUE.
                  case(55296)
		     cg_iters = 5.8511
		     nelem_found = .TRUE.
		  case default 
		     print *,"SetSerialParamsImplicit: Could not iteration count estimate for: "
		     print *,'SetSerialParamsImplicit: np := ',np, 'nelem := ',nelem
		     nelem_found = .FALSE.
                  end select
             case default  
		  print *,"SetSerialParamsImplicit: Could not iteration count estimate for: "
		  print *,'SetSerialParamsImplicit: np := ',np, 'nelem := ',nelem
         end select
	

	if(m_found) then 
	  if(np_found) then
            if(nelem_found) then 
	      if(nlev_found) found=.TRUE.
            endif
	  endif
	endif


	end subroutine SetSerialParamsImplicit

end module perfmodel_mod
